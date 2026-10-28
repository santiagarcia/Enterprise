// s3: study-3 driver (study3/PROTOCOL.md). Builds against the UNMODIFIED Enterprise/warp_hypad
// headers + libraries at 477141d (same flags as study 2's s2). Shared code paths (node sets, residual,
// gradient, compass polish, validation diagnostics, starts) are copied verbatim from study2/driver/s2.cpp.
//
// Modes
//   run     --form A|B --ws ORDER:SUB --start NAME --lambda L --kmax K --maxfev N --deadline EPOCH
//           --diag-every D --out DIR
//   traj    --ys FILE --its i,j,k --out CSV          (8:1 diagnostics for saved iterates; post hoc)
//   assess  --y FILE|BASELINE --seed S --out CSV --tag TAG [--full 1]
//   floor   (baseline minimum over D \ Ab)
//   axi     --y FILE|BASELINE
//   field   --y FILE|BASELINE --out FILE
//   timing  (per-component timing for the pilot)
//   check   (base-set identity vs build_elementwise_set; E_neg vs evaluate_feasibility; R vs repo)
//
// Per-node values are computed in parallel and summed serially in node order (deterministic).

#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <deque>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "o1_checks.hpp"
#include "warp_hypad/design/feasibility.hpp"
#include "warp_hypad/sensitivity/otilib_extract.hpp"
#include "warp_hypad/sensitivity/otilib_math.hpp"
#include "warp_hypad/sensitivity/otilib_scalar.hpp"
#include "warp_hypad/sensitivity/otilib_seed.hpp"

using namespace warp_hypad;
using namespace warp_hypad::checks;
using namespace warp_hypad::design;
using sensitivity::OTIScalar;

namespace {

using Vec = std::vector<double>;
inline std::size_t U(long i) { return static_cast<std::size_t>(i); }

constexpr double kBs2 = -9.60726526e-4;  // study-2 baseline rho*_D (tail threshold, PROTOCOL 3.4)

double dot(const Vec& a, const Vec& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}
double nrm(const Vec& a) { return std::sqrt(dot(a, a)); }
double now() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}
double epoch() { return static_cast<double>(std::time(nullptr)); }

struct Ctx {
    ShiftFamilyConfig cfg;
    F1wSystem sys;
    AdmFamilySpec f1w, f2a;
    int block = 0;
    double tau = 1e-3;
};

Ctx make_ctx() {
    Ctx c{o0c_config(), build_f1w_system(), f1w_spec(), f2a_spec(), 0, 1e-3};
    c.block = c.sys.active;
    c.tau = o0d_settings().tau;
    return c;
}

template <typename S>
AdmFieldCoefficients<S> fields_of(const Ctx& c, const std::vector<S>& y) {
    const int b = c.block;
    const int pc = c.sys.parent.coefficient_count();
    AdmFieldCoefficients<S> f;
    f.log_lapse = c.sys.family.prolong<S>(std::vector<S>(y.begin(), y.begin() + b), pc);
    f.axial_shift = c.sys.family.prolong<S>(std::vector<S>(y.begin() + b, y.begin() + 2 * b), pc);
    f.radial_shift = c.sys.family.prolong<S>(std::vector<S>(y.begin() + 2 * b, y.begin() + 3 * b), pc);
    return f;
}

double rho_pt(const Ctx& c, const AdmFieldCoefficients<double>& f, double x, double y, double z = 0.0) {
    return family_eulerian_density<double>(c.f2a, c.cfg, c.sys.parent, f, 0.0, x, y, z, 0.0,
                                           SourceSemantics::TotalEffective);
}

// ---------------------------------------------------------------------------
// Element-based node sets (leaf rectangles in (xi, s), n x n Gauss each). Verbatim from s2.
struct Elem { double xa, xb, sa, sb; int depth; };

struct NodeSet {
    int order = 4;
    std::vector<Elem> elems;
    Vec x, y, w;
    std::vector<int> elem_of;
    std::size_t base_nodes = 0;
};

std::vector<double> edges_of(const QuinticBSplineBasis1D& b, double lo, double hi, int subdivision) {
    const std::vector<double>& knots = b.distinct_interior_knots();
    double h = 0.5;
    if (knots.size() >= 2) {
        double narrowest = 1e300;
        for (std::size_t i = 0; i + 1 < knots.size(); ++i) narrowest = std::fmin(narrowest, knots[i + 1] - knots[i]);
        if (narrowest > 1e-9 && narrowest < 1e300) h = narrowest;
    }
    std::vector<double> e;
    e.push_back(lo);
    e.push_back(hi);
    for (double t : knots)
        if (t > lo + 1e-12 && t < hi - 1e-12) e.push_back(t);
    if (b.lower() > lo + 1e-12) {
        e.push_back(b.lower());
        for (double t = b.lower() - h; t > lo + 1e-12; t -= h) e.push_back(t);
    }
    if (b.upper() < hi - 1e-12) {
        e.push_back(b.upper());
        for (double t = b.upper() + h; t < hi - 1e-12; t += h) e.push_back(t);
    }
    std::sort(e.begin(), e.end());
    e.erase(std::unique(e.begin(), e.end(), [](double a, double cc) { return std::fabs(a - cc) < 1e-12; }), e.end());
    std::vector<double> out;
    for (std::size_t i = 0; i + 1 < e.size(); ++i)
        for (int k = 0; k < subdivision; ++k) out.push_back(e[i] + (e[i + 1] - e[i]) * k / subdivision);
    out.push_back(e.back());
    return out;
}

void rebuild_nodes(const Ctx& c, NodeSet& ns) {
    std::vector<double> gn, gw;
    gauss_legendre(ns.order, gn, gw);
    const double R0 = c.sys.parent.scale_R0();
    const double kPiL = 3.14159265358979323846;
    const double vol = kPiL * R0 * R0 * R0;
    ns.x.clear(); ns.y.clear(); ns.w.clear(); ns.elem_of.clear();
    for (std::size_t e = 0; e < ns.elems.size(); ++e) {
        const Elem& E = ns.elems[e];
        for (int i = 0; i < ns.order; ++i)
            for (int j = 0; j < ns.order; ++j) {
                const double xi = 0.5 * (E.xa + E.xb) + 0.5 * (E.xb - E.xa) * gn[U(i)];
                const double s = 0.5 * (E.sa + E.sb) + 0.5 * (E.sb - E.sa) * gn[U(j)];
                const double w = 0.5 * (E.xb - E.xa) * gw[U(i)] * 0.5 * (E.sb - E.sa) * gw[U(j)] * vol;
                ns.x.push_back(xi * R0 + c.cfg.velocity * 0.0);
                ns.y.push_back(R0 * std::sqrt(std::fmax(s, 0.0)));
                ns.w.push_back(w);
                ns.elem_of.push_back(static_cast<int>(e));
            }
    }
}

NodeSet base_set(const Ctx& c, int order, int sub, double xe, double se) {
    NodeSet ns;
    ns.order = order;
    const std::vector<double> ex = edges_of(c.sys.parent.axial(), -xe, xe, sub);
    const std::vector<double> es = edges_of(c.sys.parent.radial(), 0.0, se, sub);
    for (std::size_t a = 0; a + 1 < ex.size(); ++a)
        for (std::size_t b = 0; b + 1 < es.size(); ++b) {
            if (!(ex[a + 1] > ex[a]) || !(es[b + 1] > es[b])) continue;
            ns.elems.push_back({ex[a], ex[a + 1], es[b], es[b + 1], 0});
        }
    rebuild_nodes(c, ns);
    ns.base_nodes = ns.x.size();
    return ns;
}

Vec node_rho(const Ctx& c, const AdmFieldCoefficients<double>& f, const Vec& X, const Vec& Y) {
    const long n = static_cast<long>(X.size());
    Vec r(U(n));
#pragma omp parallel for schedule(dynamic, 64)
    for (long p = 0; p < n; ++p) r[U(p)] = rho_pt(c, f, X[U(p)], Y[U(p)]);
    return r;
}

double residual(const Ctx& c, const Vec& yv, const NodeSet& ns) {
    const Vec r = node_rho(c, fields_of<double>(c, yv), ns.x, ns.y);
    double acc = 0.0;
    for (std::size_t p = 0; p < r.size(); ++p) {
        const double s = softplus(-r[p], c.tau);
        acc = acc + s * s * ns.w[p];
    }
    return acc;
}

Vec gradient(const Ctx& c, const Vec& yv, const NodeSet& ns) {
    const int n = static_cast<int>(yv.size());
    Vec g(U(n), 0.0);
    for (int s0 = 0; s0 < n; s0 += 3) {
        const int count = std::min(3, n - s0);
        std::vector<OTIScalar> yy(U(n));
        for (int j = 0; j < n; ++j) yy[U(j)] = sensitivity::make_constant(yv[U(j)]);
        for (int d = 0; d < count; ++d) yy[U(s0 + d)] = sensitivity::seed_parameter(yv[U(s0 + d)], d);
        const AdmFieldCoefficients<OTIScalar> f = fields_of<OTIScalar>(c, yy);
        const long np = static_cast<long>(ns.x.size());
        std::vector<OTIScalar> rho(U(np));
#pragma omp parallel for schedule(dynamic, 16)
        for (long p = 0; p < np; ++p)
            rho[U(p)] = family_eulerian_density<OTIScalar>(c.f2a, c.cfg, c.sys.parent, f, 0.0, ns.x[U(p)],
                                                           ns.y[U(p)], 0.0, 0.0, SourceSemantics::TotalEffective);
        OTIScalar acc = scalar_cast<OTIScalar>(0.0);
        for (long p = 0; p < np; ++p) {
            const OTIScalar s = softplus(-rho[U(p)], c.tau);
            acc = acc + s * s * scalar_cast<OTIScalar>(ns.w[U(p)]);
        }
        for (int d = 0; d < count; ++d) g[U(s0 + d)] = sensitivity::extract_gradient(acc, d);
    }
    return g;
}

// ---------------------------------------------------------------------------
struct Pt { double x, r, v; };

// Compass search (verbatim from s2).
template <typename F, typename Ok>
Pt polish(F f, Ok ok, double x, double r, double v, double h, double hmin) {
    Pt b{x, r, v};
    while (h > hmin) {
        bool moved = false;
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy) {
                if (!dx && !dy) continue;
                const double xn = b.x + dx * h, rn = b.r + dy * h;
                if (!ok(xn, rn)) continue;
                const double vn = f(xn, rn);
                if (std::isfinite(vn) && vn < b.v) { b = {xn, rn, vn}; moved = true; }
            }
        if (!moved) h *= 0.5;
    }
    return b;
}

std::vector<Pt> pick_distinct(std::vector<Pt> pts, std::size_t k, double sep) {
    std::sort(pts.begin(), pts.end(), [](const Pt& a, const Pt& b) { return a.v < b.v; });
    std::vector<Pt> out;
    for (const Pt& p : pts) {
        bool far = true;
        for (const Pt& q : out)
            if (std::hypot(p.x - q.x, p.r - q.r) < sep) { far = false; break; }
        if (far) out.push_back(p);
        if (out.size() >= k) break;
    }
    return out;
}

// Discrete 8-neighbour local minima of an NX x NR lattice (row-major i*NR+j), optionally only v < vmax.
std::vector<Pt> lattice_minima(const Vec& L, int NX, int NR, double x0, double dx, double r0, double dr,
                               double vmax) {
    std::vector<Pt> out;
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NR; ++j) {
            const double v = L[U(i * NR + j)];
            if (!(v < vmax)) continue;
            bool mn = true;
            for (int a = -1; a <= 1 && mn; ++a)
                for (int b = -1; b <= 1; ++b) {
                    if (!a && !b) continue;
                    const int ii = i + a, jj = j + b;
                    if (ii < 0 || jj < 0 || ii >= NX || jj >= NR) continue;
                    if (L[U(ii * NR + jj)] < v) { mn = false; break; }
                }
            if (mn) out.push_back({x0 + dx * i, r0 + dr * j, v});
        }
    std::sort(out.begin(), out.end(), [](const Pt& a, const Pt& b) { return a.v < b.v; });
    return out;
}

// In-loop search on Dw (PROTOCOL 3.2): 241 x 116 lattice; polish discrete lattice minima (rho < 0) in rank order,
// batches of 8, merge polished results within 1e-3, until 8 distinct or 64 polishes.
struct InLoop { std::vector<Pt> mins; int polished = 0; int duplicates = 0; int lattice_minima = 0; };

InLoop inloop_search(const Ctx& c, const AdmFieldCoefficients<double>& f, double xw, double rw) {
    const int NX = 241, NR = 116;
    Vec lat(U(NX * NR));
#pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NR; ++j) {
            const double x = -xw + 2.0 * xw * i / (NX - 1), r = rw * j / (NR - 1);
            lat[U(i * NR + j)] = rho_pt(c, f, x, r);
        }
    const double dx = 2.0 * xw / (NX - 1), dr = rw / (NR - 1);
    const std::vector<Pt> st = lattice_minima(lat, NX, NR, -xw, dx, 0.0, dr, 0.0);
    InLoop out;
    out.lattice_minima = static_cast<int>(st.size());
    auto fr = [&](double x, double r) { return rho_pt(c, f, x, r); };
    auto ok = [&](double x, double r) { return x >= -xw && x <= xw && r >= 0.0 && r <= rw; };
    std::size_t next = 0;
    while (out.mins.size() < 8 && out.polished < 64 && next < st.size()) {
        const int nb = static_cast<int>(std::min<std::size_t>({8, st.size() - next, U(64 - out.polished)}));
        std::vector<Pt> res(U(nb));
#pragma omp parallel for schedule(dynamic, 1)
        for (int k = 0; k < nb; ++k) {
            const Pt& s = st[next + U(k)];
            res[U(k)] = polish(fr, ok, s.x, s.r, s.v, dx, 1e-9);
        }
        for (const Pt& p : res) {  // serial merge in rank order (deterministic)
            bool dup = false;
            for (Pt& q : out.mins)
                if (std::hypot(p.x - q.x, p.r - q.r) < 1e-3) {
                    dup = true;
                    if (p.v < q.v) q = p;
                    break;
                }
            if (dup) ++out.duplicates; else out.mins.push_back(p);
        }
        out.polished += nb;
        next += U(nb);
    }
    std::sort(out.mins.begin(), out.mins.end(), [](const Pt& a, const Pt& b) { return a.v < b.v; });
    if (out.mins.size() > 8) out.mins.resize(8);
    if (out.mins.empty()) {  // no negative lattice minimum: report the lattice minimum itself
        Pt b{0, 0, 1e300};
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j < NR; ++j)
                if (lat[U(i * NR + j)] < b.v) b = {-xw + dx * i, dr * j, lat[U(i * NR + j)]};
        out.mins.push_back(b);
    }
    return out;
}

// B refinement (PROTOCOL 3.3): study-2 rule, every leaf whose closed rectangle contains the point is tested.
int refine_B(const Ctx& c, const AdmFieldCoefficients<double>& f, NodeSet& ns, const std::vector<Pt>& cand,
             std::ofstream& ev, int it, bool& capped) {
    const double R0 = c.sys.parent.scale_R0();
    int splits = 0;
    for (const Pt& p : cand) {
        if (!(p.v < 0.0)) continue;
        const double xi = p.x / R0, s = (p.r / R0) * (p.r / R0);
        for (int guard = 0; guard < 16; ++guard) {
            std::vector<int> Ls;
            for (std::size_t e = 0; e < ns.elems.size(); ++e) {
                const Elem& E = ns.elems[e];
                if (xi >= E.xa && xi <= E.xb && s >= E.sa && s <= E.sb) Ls.push_back(static_cast<int>(e));
            }
            if (Ls.empty()) break;  // outside Dw
            std::vector<int> todo;
            for (int L : Ls) {
                double rhoL = 1e300;
                for (std::size_t q = 0; q < ns.x.size(); ++q)
                    if (ns.elem_of[q] == L) rhoL = std::fmin(rhoL, rho_pt(c, f, ns.x[q], ns.y[q]));
                const Elem& E = ns.elems[U(L)];
                const bool need = p.v < rhoL - 0.05 * std::fabs(p.v);
                ev << it << "," << p.x << "," << p.r << "," << p.v << "," << rhoL << "," << E.depth << "," << L << ","
                   << Ls.size() << "," << (need ? 1 : 0);
                if (!need) { ev << ",resolved\n"; continue; }
                if (E.depth >= 6) { ev << ",max_depth\n"; continue; }
                if (ns.x.size() + 3u * (todo.size() + 1) * static_cast<std::size_t>(ns.order * ns.order) >
                    4u * ns.base_nodes) {
                    ev << ",node_cap\n"; capped = true; continue;
                }
                ev << ",split\n";
                todo.push_back(L);
            }
            if (todo.empty()) break;
            for (int L : todo) {  // indices of untouched leaves are stable: children are appended
                const Elem E = ns.elems[U(L)];
                const double xm = 0.5 * (E.xa + E.xb), sm = 0.5 * (E.sa + E.sb);
                const int d = E.depth + 1;
                ns.elems[U(L)] = {E.xa, xm, E.sa, sm, d};
                ns.elems.push_back({xm, E.xb, E.sa, sm, d});
                ns.elems.push_back({E.xa, xm, sm, E.sb, d});
                ns.elems.push_back({xm, E.xb, sm, E.sb, d});
                ++splits;
            }
            rebuild_nodes(c, ns);
        }
    }
    ev.flush();
    return splits;
}

// 8:1 validation diagnostics (s2 + Dw-restricted E_neg, tail volume, L4). gamma = I => dV = w.
struct Diag { double Eneg, EnegW, negfrac, vtail, L4, rhomin, rx, ry, refined; };
Diag validation_diag(const Ctx& c, const Vec& yv, const IntegrationSet& vs, double eps_rho) {
    const AdmFieldCoefficients<double> f = fields_of<double>(c, yv);
    const Vec r = node_rho(c, f, vs.x, vs.y);
    Diag d{0, 0, 0, 0, 0, 1e300, 0, 0, 0};
    double vneg = 0, vtot = 0, l4 = 0;
    const double XW = 3.4, RW = 2.0 * std::sqrt(2.6);
    for (std::size_t p = 0; p < r.size(); ++p) {
        vtot += vs.weight[p];
        if (r[p] < 0) {
            d.Eneg += -r[p] * vs.weight[p];
            vneg += vs.weight[p];
            l4 += r[p] * r[p] * r[p] * r[p] * vs.weight[p];
            if (std::fabs(vs.x[p]) <= XW && vs.y[p] <= RW) d.EnegW += -r[p] * vs.weight[p];
        }
        if (r[p] < 1.05 * kBs2) d.vtail += vs.weight[p];
        if (r[p] < d.rhomin) { d.rhomin = r[p]; d.rx = vs.x[p]; d.ry = vs.y[p]; }
    }
    d.negfrac = vneg / vtot;
    d.L4 = std::pow(l4, 0.25);
    FeasibilityReport fr;
    fr.rho_min = d.rhomin; fr.rho_min_x = d.rx; fr.rho_min_y = d.ry; fr.rho_min_z = 0;
    refine_around_minimum(fr, c.f2a, c.cfg, c.sys.parent, f, 0.0, 0.35, 24, 0.0, SourceSemantics::TotalEffective,
                          eps_rho);
    d.refined = fr.refined_rho_min;
    return d;
}

Vec read_vec(const std::string& path) {
    std::ifstream f(path);
    Vec v;
    double x;
    while (f >> x) v.push_back(x);
    return v;
}

std::vector<std::pair<std::string, Vec>> starts_of(int n, int block) {
    std::vector<std::pair<std::string, Vec>> s;
    s.push_back({"baseline_zero", Vec(U(n), 0.0)});
    const char* nm[3] = {"lapse_perturbation", "axial_shift_perturbation", "radial_shift_perturbation"};
    for (int field = 0; field < 3; ++field) {
        Vec p(U(n), 0.0);
        for (int i = 0; i < block; ++i) p[U(field * block + i)] = 0.05 * std::sin(1.3 * i + 0.7 + field);
        s.push_back({nm[field], p});
    }
    Vec comb(U(n), 0.0);
    for (int i = 0; i < n; ++i) comb[U(i)] = 0.04 * std::sin(0.9 * i + 0.2);
    s.push_back({"combined_perturbation", comb});
    return s;
}

struct Val { IntegrationSet vs; double eps_rho; };
Val make_val(const Ctx& c) {
    Val v{build_elementwise_set(c.sys.parent, 8, 1, 3.5, 12.25, 0.0, c.cfg.velocity), 0.0};
    double tv = 0; for (double q : v.vs.weight) tv += q;
    const double eps_E = std::fmax(o0d_state().numerical_uncertainty, 1e-9) * 10.0;
    v.eps_rho = eps_E / tv;
    return v;
}

// ---------------------------------------------------------------------------
int mode_run(const Ctx& c, char form, int order, int sub, const std::string& start, double lam, int kmax, int maxfev,
             double deadline, int diag_every, const std::string& out) {
    const double t_start = now();
    const int n = 3 * c.block;
    Vec y;
    for (auto& s : starts_of(n, c.block))
        if (s.first == start) y = s.second;
    if (y.empty()) { std::fprintf(stderr, "unknown start %s\n", start.c_str()); return 2; }
    NodeSet ns = base_set(c, order, sub, 1.7, 2.6);
    const Val V = make_val(c);
    const double xw = 1.7 * 2.0, rw = 2.0 * std::sqrt(2.6);
    {
        std::ofstream cf(out + "/config.txt");
        cf.precision(10);
        cf << "form " << form << "\norder " << order << "\nsub " << sub << "\nstart " << start << "\nlambda " << lam
           << "\nkmax " << kmax << "\nmaxfev " << maxfev << "\ndeadline_epoch " << deadline << "\ndiag_every "
           << diag_every << "\nbase_nodes " << ns.base_nodes << "\ntau " << c.tau << "\nthreads "
           << omp_get_max_threads()
           << "\nF = R + lambda/2 |y|^2\nC1' |gradF| <= 1e-5 F / max(|y|,1)\nC2' 10-window 1e-8 since last set change\n"
              "STALLED C2' & !C1' for 30 consecutive iterates on one set\n"
              "lbfgs m=8 c1=1e-4 trust0=0.5 trustmax=4 trustmin=1e-8 ls_trials=10; pairs kept across set changes;\n"
              "post-reset step p = -gamma_last g\n";
    }
    std::ofstream tr(out + "/trace.csv");
    tr.precision(12);
    tr << "it,F,R,pen,ynorm,gF,gR,step,trust,gamma,reset,nfev,ngev,nodes,leaves,splits,rhomin_w_pre,rhomin_w_post,"
          "diag,Eneg_v8,Eneg_v8_Dw,negfrac_v8,vtail_v8,L4_v8,rhomin_v8,refined_v8,rhostar_in,xstar,rstar,n_distinct,"
          "n_polished,n_dup,n_latmin,C1,C2,C3,t_grad,t_search,t_diag,secs\n";
    std::ofstream ys(out + "/ys.txt");
    ys.precision(17);
    std::ofstream ev(out + "/refine_events.csv");
    ev.precision(12);
    ev << "it,x,r,rho_star,rho_L,depth,leaf,n_containing,need,action\n";
    std::ofstream ls(out + "/linesearch.csv");
    ls.precision(12);
    ls << "it,trial,alpha,trust,F_trial,F_cur,armijo,ok\n";
    std::ofstream im(out + "/inloop_minima.csv");
    im.precision(12);
    im << "it,rank,x,r,rho\n";

    int nfev = 0, ngev = 0, resets = 0;
    bool capped = false;
    double t_grad = 0, t_search = 0, t_diag = 0;
    auto Fval = [&](const Vec& yy, const NodeSet& s, double& Rout) {
        Rout = residual(c, yy, s);
        return Rout + 0.5 * lam * dot(yy, yy);
    };
    auto Fgrad = [&](const Vec& yy, const NodeSet& s, Vec& gR) {
        const double t0 = now();
        gR = gradient(c, yy, s);
        Vec gF = gR;
        for (std::size_t i = 0; i < gF.size(); ++i) gF[i] += lam * yy[i];
        t_grad = now() - t0;
        return gF;
    };
    auto search = [&]() {
        const double t0 = now();
        InLoop r = inloop_search(c, fields_of<double>(c, y), xw, rw);
        t_search = now() - t0;
        return r;
    };
    auto wmin = [&]() {
        const Vec rr = node_rho(c, fields_of<double>(c, y), ns.x, ns.y);
        double m = 1e300; for (double v : rr) m = std::fmin(m, v);
        return m;
    };
    auto log_iter = [&](int it, double F, double R, const Vec& gF, const Vec& gR, double step, double trust,
                        double gamma, int reset, int splits, double rmw_pre, bool C1, bool C2, bool C3,
                        const InLoop& in, bool force_diag) {
        const double rmw_post = splits > 0 ? wmin() : rmw_pre;
        const bool dg = force_diag || (it % diag_every == 0);
        Diag d{NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN};
        t_diag = 0;
        if (dg) { const double t0 = now(); d = validation_diag(c, y, V.vs, V.eps_rho); t_diag = now() - t0; }
        tr << it << "," << F << "," << R << "," << 0.5 * lam * dot(y, y) << "," << nrm(y) << "," << nrm(gF) << ","
           << nrm(gR) << "," << step << "," << trust << "," << gamma << "," << reset << "," << nfev << "," << ngev
           << "," << ns.x.size() << "," << ns.elems.size() << "," << splits << "," << rmw_pre << "," << rmw_post
           << "," << (dg ? 1 : 0) << "," << d.Eneg << "," << d.EnegW << "," << d.negfrac << "," << d.vtail << ","
           << d.L4 << "," << d.rhomin << "," << d.refined << "," << in.mins[0].v << "," << in.mins[0].x << ","
           << in.mins[0].r << "," << in.mins.size() << "," << in.polished << "," << in.duplicates << ","
           << in.lattice_minima << "," << C1 << "," << C2 << "," << C3 << "," << t_grad << "," << t_search << ","
           << t_diag << "," << (now() - t_start) << "\n";
        tr.flush();
        for (std::size_t k = 0; k < in.mins.size(); ++k)
            im << it << "," << k << "," << in.mins[k].x << "," << in.mins[k].r << "," << in.mins[k].v << "\n";
        im.flush();
        for (double v : y) ys << v << " ";
        ys << "\n";
        ys.flush();
    };

    double R = 0.0;
    double f = Fval(y, ns, R); ++nfev;
    Vec gR;
    Vec g = Fgrad(y, ns, gR); ++ngev;
    InLoop cand = search();
    int last_change = 0;
    std::vector<double> fh;
    int sp0 = 0;
    const double rmw0 = wmin();
    if (form == 'B') {
        sp0 = refine_B(c, fields_of<double>(c, y), ns, cand.mins, ev, 0, capped);
        if (sp0 > 0) { f = Fval(y, ns, R); ++nfev; g = Fgrad(y, ns, gR); ++ngev; }
    }
    log_iter(0, f, R, g, gR, 0.0, 0.5, NAN, 0, sp0, rmw0, false, false, sp0 == 0, cand, true);
    fh.push_back(f);
    double trust = 0.5, gamma_last = NAN;
    std::deque<Vec> S, Y;
    std::string term;
    const double kC1 = 1e-4;
    int it = 0, stall_run = 0;
    bool C1 = false, C2 = false;
    while (true) {
        if (it >= kmax) { term = "ITERATION_LIMIT"; break; }
        if (nfev >= maxfev) { term = "EVAL_LIMIT"; break; }
        if (deadline > 0 && epoch() >= deadline) { term = "WALL_LIMIT"; break; }
        Vec p(U(n));
        int reset = 0;
        for (int i = 0; i < n; ++i) p[U(i)] = -g[U(i)];
        if (!S.empty()) {
            Vec q = g;
            std::vector<double> al(S.size(), 0.0);
            for (int k = static_cast<int>(S.size()) - 1; k >= 0; --k) {
                const std::size_t K = U(k);
                const double r2 = 1.0 / dot(Y[K], S[K]);
                al[K] = r2 * dot(S[K], q);
                for (int i = 0; i < n; ++i) q[U(i)] -= al[K] * Y[K][U(i)];
            }
            const std::size_t L = S.size() - 1;
            const double gam = dot(S[L], Y[L]) / dot(Y[L], Y[L]);
            for (int i = 0; i < n; ++i) q[U(i)] *= gam;
            for (std::size_t k = 0; k < S.size(); ++k) {
                const double r2 = 1.0 / dot(Y[k], S[k]);
                const double b = r2 * dot(Y[k], q);
                for (int i = 0; i < n; ++i) q[U(i)] += (al[k] - b) * S[k][U(i)];
            }
            for (int i = 0; i < n; ++i) p[U(i)] = -q[U(i)];
        } else if (it > 0 && std::isfinite(gamma_last)) {
            for (int i = 0; i < n; ++i) p[U(i)] = -gamma_last * g[U(i)];
            reset = 1;
        }
        double dd = dot(p, g);
        if (!(dd < 0.0)) {
            S.clear(); Y.clear();
            const double sc = std::isfinite(gamma_last) ? gamma_last : 1.0;
            for (int i = 0; i < n; ++i) p[U(i)] = -sc * g[U(i)];
            dd = dot(p, g);
            reset = 2;
            ++resets;
        }
        const double pn = nrm(p);
        double alpha = pn > 0.0 ? std::fmin(1.0, trust / pn) : 0.0;
        bool ok = false;
        Vec yn;
        double fn = 0.0, step = 0.0, Rn = 0.0;
        int att = 0;
        for (; att < 10; ++att) {
            yn = y;
            for (int i = 0; i < n; ++i) yn[U(i)] += alpha * p[U(i)];
            step = 0.0;
            for (int i = 0; i < n; ++i) { const double d = yn[U(i)] - y[U(i)]; step += d * d; }
            step = std::sqrt(step);
            bool finite = true;
            for (double q : yn) if (!std::isfinite(q)) finite = false;
            bool armijo = false;
            if (finite) {
                fn = Fval(yn, ns, Rn); ++nfev;
                const MetricJet<double> probe =
                    adm_metric_jet<double>(c.f2a, c.cfg, c.sys.parent, fields_of<double>(c, yn), 0.0, 0.7, 1.3, 0.8);
                armijo = std::isfinite(fn) && probe.valid && fn <= f + kC1 * alpha * dd;
            }
            ls << it + 1 << "," << att << "," << alpha << "," << trust << "," << fn << "," << f << ","
               << f + kC1 * alpha * dd << "," << armijo << "\n";
            if (armijo) { ok = true; break; }
            alpha *= 0.5;
        }
        ls.flush();
        if (!ok) {
            trust *= 0.5;
            if (trust < 1e-8) { term = "LINESEARCH_FAIL"; break; }
            continue;
        }
        Vec gR2;
        const Vec g2 = Fgrad(yn, ns, gR2); ++ngev;
        bool finite = std::isfinite(fn);
        for (double q : g2) if (!std::isfinite(q)) finite = false;
        if (!finite) { term = "INVALID"; break; }
        Vec sv(U(n)), yv(U(n));
        for (int i = 0; i < n; ++i) { sv[U(i)] = yn[U(i)] - y[U(i)]; yv[U(i)] = g2[U(i)] - g[U(i)]; }
        if (dot(sv, yv) > 1e-14 * nrm(sv) * nrm(yv)) {
            S.push_back(sv); Y.push_back(yv);
            while (S.size() > 8) { S.pop_front(); Y.pop_front(); }
            gamma_last = dot(sv, yv) / dot(yv, yv);
        }
        y = yn; g = g2; gR = gR2; f = fn; R = Rn;
        ++it;
        if (att == 0) trust = std::fmin(trust * 2.0, 4.0);
        fh.push_back(f);
        const std::size_t m = fh.size() - 1;
        C1 = nrm(g) <= 1e-5 * f / std::fmax(nrm(y), 1.0);
        C2 = m >= 10 && (fh[m - 10] - fh[m]) / fh[m] <= 1e-8;
        stall_run = (C2 && !C1) ? stall_run + 1 : 0;
        int splits = 0;
        cand = search();
        const double f_acc = f, R_acc = R;
        const Vec g_acc = g, gR_acc = gR;
        const double rmw_pre = wmin();
        if (form == 'B' && !capped) {
            splits = refine_B(c, fields_of<double>(c, y), ns, cand.mins, ev, it, capped);
            if (splits > 0) {
                f = Fval(y, ns, R); ++nfev;
                g = Fgrad(y, ns, gR); ++ngev;
                fh.clear(); fh.push_back(f);  // L-BFGS pairs are KEPT (PROTOCOL 3.1)
                stall_run = 0;
                last_change = it;
                std::ofstream(out + "/set_changes.csv", std::ios::app)
                    << it << "," << splits << "," << ns.x.size() << "," << f_acc << "," << f << "," << nrm(g) << "\n";
            }
        }
        const bool C3 = (form != 'B') || splits == 0;
        const bool conv = C1 && C2 && C3;
        log_iter(it, f_acc, R_acc, g_acc, gR_acc, step, trust, gamma_last, reset, splits, rmw_pre, C1, C2, C3, cand,
                 conv || it >= kmax);
        std::printf("it %3d F %.10e R %.8e |y| %.3f |gF| %.3e (tol %.2e) step %.3e nodes %zu sp %d C %d%d%d %.0fs\n",
                    it, f_acc, R_acc, nrm(y), nrm(g_acc), 1e-5 * f_acc / std::fmax(nrm(y), 1.0), step, ns.x.size(),
                    splits, C1, C2, C3, now() - t_start);
        std::fflush(stdout);
        if (conv) { term = "CONVERGED"; break; }
        if (stall_run >= 30) { term = "STALLED"; break; }
    }
    if (term == "LINESEARCH_FAIL" && C1 && C2) term = "CONVERGED_LS_LIMITED";
    if (term != "CONVERGED") {  // final diagnostics row for the last accepted iterate if it was skipped
        std::ofstream fd(out + "/final_diag.csv");
        fd.precision(12);
        const Diag d = validation_diag(c, y, V.vs, V.eps_rho);
        fd << "it,Eneg_v8,Eneg_v8_Dw,negfrac_v8,vtail_v8,L4_v8,rhomin_v8,refined_v8,rhomin_w\n"
           << it << "," << d.Eneg << "," << d.EnegW << "," << d.negfrac << "," << d.vtail << "," << d.L4 << ","
           << d.rhomin << "," << d.refined << "," << wmin() << "\n";
    }
    {
        std::ofstream fy(out + "/y_final.txt");
        fy.precision(17);
        for (double v : y) fy << v << "\n";
        std::ofstream ft(out + "/termination.txt");
        ft.precision(10);
        ft << term << " it=" << it << " nfev=" << nfev << " ngev=" << ngev << " nodes=" << ns.x.size()
           << " capped=" << capped << " last_set_change=" << last_change << " resets=" << resets
           << " gF=" << nrm(g) << " C1tol=" << 1e-5 * f / std::fmax(nrm(y), 1.0) << " F=" << f << " R=" << R
           << " ynorm=" << nrm(y) << " wall_s=" << (now() - t_start) << "\n";
    }
    std::printf("TERMINATION %s it=%d nfev=%d wall=%.0fs\n", term.c_str(), it, nfev, now() - t_start);
    return 0;
}

// Post-hoc 8:1 diagnostics for saved iterates (row k of ys.txt = iterate k).
int mode_traj(const Ctx& c, const std::string& ysfile, const std::string& its, const std::string& out) {
    std::vector<int> want;
    {
        std::stringstream ss(its);
        std::string t;
        while (std::getline(ss, t, ',')) if (!t.empty()) want.push_back(std::stoi(t));
    }
    std::ifstream f(ysfile);
    std::string line;
    int k = 0;
    const Val V = make_val(c);
    std::ofstream o(out, std::ios::app);
    o.precision(12);
    while (std::getline(f, line)) {
        if (std::find(want.begin(), want.end(), k) != want.end()) {
            std::stringstream ss(line);
            Vec y; double v;
            while (ss >> v) y.push_back(v);
            const Diag d = validation_diag(c, y, V.vs, V.eps_rho);
            o << ysfile << "," << k << "," << d.Eneg << "," << d.EnegW << "," << d.negfrac << "," << d.vtail << ","
              << d.L4 << "," << d.rhomin << "," << d.refined << "\n";
            o.flush();
        }
        ++k;
    }
    return 0;
}

// ---------------------------------------------------------------------------
double halton(unsigned i, unsigned b) {
    double f = 1.0, r = 0.0;
    while (i > 0) { f /= b; r += f * (i % b); i /= b; }
    return r;
}

// 2-D Nelder-Mead (standard coefficients 1, 2, 0.5, 0.5); returns best vertex.
template <typename F>
Pt nelder_mead(F f, double x, double r, double v, double d0, double tol, int maxit) {
    Pt s[3] = {{x, r, v}, {x + d0, r, 0}, {x, r + d0, 0}};
    s[1].v = f(s[1].x, s[1].r);
    s[2].v = f(s[2].x, s[2].r);
    for (int k = 0; k < maxit; ++k) {
        std::sort(s, s + 3, [](const Pt& a, const Pt& b) { return a.v < b.v; });
        const double diam = std::fmax(std::hypot(s[1].x - s[0].x, s[1].r - s[0].r),
                                      std::hypot(s[2].x - s[0].x, s[2].r - s[0].r));
        if (diam < tol) break;
        const double cx = 0.5 * (s[0].x + s[1].x), cr = 0.5 * (s[0].r + s[1].r);
        const Pt xr{cx + (cx - s[2].x), cr + (cr - s[2].r), 0};
        const double vr = f(xr.x, xr.r);
        if (vr < s[0].v) {
            const Pt xe{cx + 2 * (cx - s[2].x), cr + 2 * (cr - s[2].r), 0};
            const double ve = f(xe.x, xe.r);
            s[2] = ve < vr ? Pt{xe.x, xe.r, ve} : Pt{xr.x, xr.r, vr};
        } else if (vr < s[1].v) {
            s[2] = {xr.x, xr.r, vr};
        } else {
            const bool outside = vr < s[2].v;
            const Pt xc = outside ? Pt{cx + 0.5 * (xr.x - cx), cr + 0.5 * (xr.r - cr), 0}
                                  : Pt{cx + 0.5 * (s[2].x - cx), cr + 0.5 * (s[2].r - cr), 0};
            const double vc = f(xc.x, xc.r);
            if (vc < (outside ? vr : s[2].v)) {
                s[2] = {xc.x, xc.r, vc};
            } else {
                for (int j = 1; j < 3; ++j) {
                    s[j].x = s[0].x + 0.5 * (s[j].x - s[0].x);
                    s[j].r = s[0].r + 0.5 * (s[j].r - s[0].r);
                    s[j].v = f(s[j].x, s[j].r);
                }
            }
        }
    }
    std::sort(s, s + 3, [](const Pt& a, const Pt& b) { return a.v < b.v; });
    return s[0];
}

struct BoxSearch {
    Pt best;
    std::vector<Pt> minima;  // distinct polished minima (merge 1e-3)
    double sample_min = 0, min_est = 0, prepolish_best = 0, sub_min = 0;
    std::size_t flagged = 0, starts = 0;
    bool d7 = false;
};

// Assessment search on a box (PROTOCOL 7). ok(x,r) confines the polish (box membership, optional exclusion).
template <typename Ok>
BoxSearch box_search(const Ctx& c, const AdmFieldCoefficients<double>& f, double x0, double x1, double r0, double r1,
                     double h, unsigned seed, int nhalton, Ok ok) {
    auto fr = [&](double x, double r) { return ok(x, r) ? rho_pt(c, f, x, r) : INFINITY; };
    const int NX = static_cast<int>(std::lround((x1 - x0) / h)) + 1, NR = static_cast<int>(std::lround((r1 - r0) / h)) + 1;
    Vec L(U(NX * NR));
#pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NR; ++j) L[U(i * NR + j)] = fr(x0 + h * i, r0 + h * j);
    BoxSearch B;
    // Halton
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    const double sh1 = u01(rng), sh2 = u01(rng);
    std::vector<Pt> H(U(nhalton));
#pragma omp parallel for schedule(dynamic, 64)
    for (int k = 0; k < nhalton; ++k) {
        const double a = std::fmod(halton(static_cast<unsigned>(k + 1), 2) + sh1, 1.0);
        const double b = std::fmod(halton(static_cast<unsigned>(k + 1), 3) + sh2, 1.0);
        const double x = x0 + (x1 - x0) * a, r = r0 + (r1 - r0) * b;
        H[U(k)] = {x, r, fr(x, r)};
    }
    double best = 1e300;
    for (double v : L) best = std::fmin(best, v);
    for (const Pt& p : H) best = std::fmin(best, p.v);
    // cell screening: central-difference |grad| at lattice points
    Vec G(U(NX * NR), 0.0);
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NR; ++j) {
            const int ia = std::max(i - 1, 0), ib = std::min(i + 1, NX - 1);
            const int ja = std::max(j - 1, 0), jb = std::min(j + 1, NR - 1);
            const double a = L[U(ib * NR + j)], b = L[U(ia * NR + j)], cc = L[U(i * NR + jb)], d = L[U(i * NR + ja)];
            if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(cc) || !std::isfinite(d)) { G[U(i * NR + j)] = 0; continue; }
            const double gx = (a - b) / (h * (ib - ia)), gr = (cc - d) / (h * (jb - ja));
            G[U(i * NR + j)] = std::hypot(gx, gr);
        }
    struct Cell { int i, j; double est; };
    std::vector<Cell> flagged;
    double min_est = 1e300;
    for (int i = 0; i + 1 < NX; ++i)
        for (int j = 0; j + 1 < NR; ++j) {
            const double m = std::fmin(std::fmin(L[U(i * NR + j)], L[U((i + 1) * NR + j)]),
                                       std::fmin(L[U(i * NR + j + 1)], L[U((i + 1) * NR + j + 1)]));
            if (!std::isfinite(m)) continue;
            const double Lc = std::fmax(std::fmax(G[U(i * NR + j)], G[U((i + 1) * NR + j)]),
                                        std::fmax(G[U(i * NR + j + 1)], G[U((i + 1) * NR + j + 1)]));
            const double est = m - 2.0 * Lc * h / std::sqrt(2.0);
            min_est = std::fmin(min_est, est);
            if (est < best - 0.1 * std::fabs(best)) flagged.push_back({i, j, est});
        }
    std::sort(flagged.begin(), flagged.end(), [](const Cell& a, const Cell& b) { return a.est < b.est; });
    B.flagged = flagged.size();
    if (flagged.size() > 4096) flagged.resize(4096);
    std::vector<Pt> subbest(flagged.size());
    const long nfl = static_cast<long>(flagged.size());
#pragma omp parallel for schedule(dynamic, 4)
    for (long k = 0; k < nfl; ++k) {
        Pt b{0, 0, INFINITY};
        for (int a = 0; a < 4; ++a)
            for (int q = 0; q < 4; ++q) {
                const double x = x0 + h * (flagged[U(k)].i + (a + 0.5) / 4.0);
                const double r = r0 + h * (flagged[U(k)].j + (q + 0.5) / 4.0);
                const double v = fr(x, r);
                if (v < b.v) b = {x, r, v};
            }
        subbest[U(k)] = b;
    }
    B.sub_min = INFINITY;
    for (const Pt& p : subbest) B.sub_min = std::fmin(B.sub_min, p.v);
    B.prepolish_best = best;
    B.d7 = B.sub_min < best;
    B.min_est = min_est;
    // starts
    std::vector<Pt> pool = lattice_minima(L, NX, NR, x0, h, r0, h, INFINITY);
    for (const Pt& p : subbest) pool.push_back(p);
    {
        std::vector<Pt> hs = H;
        std::sort(hs.begin(), hs.end(), [](const Pt& a, const Pt& b) { return a.v < b.v; });
        for (std::size_t k = 0; k < hs.size() && k < 64; ++k) pool.push_back(hs[k]);
    }
    std::vector<Pt> fin;
    for (const Pt& p : pool) if (std::isfinite(p.v)) fin.push_back(p);
    const std::vector<Pt> st = pick_distinct(fin, 128, 0.01);
    B.starts = st.size();
    B.sample_min = std::fmin(best, B.sub_min);
    std::vector<Pt> res(st.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (int k = 0; k < static_cast<int>(st.size()); ++k) {
        Pt p = nelder_mead(fr, st[U(k)].x, st[U(k)].r, st[U(k)].v, 0.5 * h, 1e-11, 4000);
        for (int rep = 0; rep < 2; ++rep) p = nelder_mead(fr, p.x, p.r, p.v, 0.25 * h, 1e-11, 4000);
        res[U(k)] = p;
    }
    std::sort(res.begin(), res.end(), [](const Pt& a, const Pt& b) { return a.v < b.v; });
    B.minima = pick_distinct(res, 8, 1e-3);
    B.best = B.minima.empty() ? Pt{0, 0, best} : B.minima[0];
    if (B.sample_min < B.best.v) {  // never report a polished value above a sampled one
        for (const Pt& p : fin) if (p.v == B.sample_min) { B.best = p; break; }
    }
    return B;
}

constexpr double kAbX = 2.1, kAbR = 3.1;
bool in_Ab(double x, double r) { return std::fabs(x) <= kAbX && r >= 0 && r <= kAbR; }

int mode_assess(const Ctx& c, const std::string& yfile, unsigned seed, const std::string& out, const std::string& tag,
                bool full) {
    const double t0 = now();
    const int n = 3 * c.block;
    Vec yv = yfile == "BASELINE" ? Vec(U(n), 0.0) : read_vec(yfile);
    const AdmFieldCoefficients<double> f = fields_of<double>(c, yv);
    const AdmFieldCoefficients<double> f0 = fields_of<double>(c, Vec(U(n), 0.0));
    // design-independence outside Ab on a 0.05 lattice over D
    long mism = 0, chk = 0;
    {
        std::vector<char> bad(281 * 141, 0), used(281 * 141, 0);
#pragma omp parallel for schedule(dynamic, 4)
        for (int i = 0; i < 281; ++i)
            for (int j = 0; j < 141; ++j) {
                const double x = -7.0 + 0.05 * i, r = 0.05 * j;
                if (in_Ab(x, r)) continue;
                used[U(i * 141 + j)] = 1;
                bad[U(i * 141 + j)] = rho_pt(c, f, x, r) != rho_pt(c, f0, x, r);
            }
        for (std::size_t k = 0; k < bad.size(); ++k) { mism += bad[k]; chk += used[k]; }
    }
    BoxSearch S;
    bool fallback = false;
    if (mism == 0) {
        S = box_search(c, f, -kAbX, kAbX, 0.0, kAbR, 0.005, seed, 1 << 17,
                       [](double x, double r) { return in_Ab(x, r); });
    } else {
        fallback = true;
        S = box_search(c, f, -7.0, 7.0, 0.0, 7.0, 0.005, seed, 1 << 17,
                       [](double x, double r) { return std::fabs(x) <= 7.0 && r >= 0 && r <= 7.0; });
    }
    // 8:1 metrics
    const Val V = make_val(c);
    const Diag d = validation_diag(c, yv, V.vs, V.eps_rho);
    std::ofstream o(out, std::ios::app);
    o.precision(12);
    o << tag << "," << yfile << "," << seed << "," << S.best.v << "," << S.best.x << "," << S.best.r << ","
      << S.sample_min << "," << S.prepolish_best << "," << S.sub_min << "," << S.min_est << "," << S.flagged << ","
      << S.d7 << "," << S.starts << "," << mism << "," << chk << "," << fallback << "," << d.Eneg << "," << d.EnegW
      << "," << d.negfrac << "," << d.vtail << "," << d.L4 << "," << d.rhomin << "," << d.refined;
    // minima list
    {
        std::ofstream om(out + "." + tag + ".s" + std::to_string(seed) + ".minima.csv");
        om.precision(12);
        om << "x,r,rho\n";
        for (const Pt& p : S.minima) om << p.x << "," << p.r << "," << p.v << "\n";
    }
    if (full) {
        // E16
        const IntegrationSet s16 = build_elementwise_set(c.sys.parent, 16, 1, 3.5, 12.25, 0.0, c.cfg.velocity);
        const Vec r16 = node_rho(c, f, s16.x, s16.y);
        double E16 = 0, m16 = 1e300;
        for (std::size_t p = 0; p < r16.size(); ++p) {
            if (r16[p] < 0) E16 += -r16[p] * s16.weight[p];
            m16 = std::fmin(m16, r16[p]);
        }
        // R on Dw at orders 4, 6, 8
        double Ro[3];
        const int ords[3] = {4, 6, 8};
        for (int k = 0; k < 3; ++k) Ro[k] = residual(c, yv, base_set(c, ords[k], 1, 1.7, 2.6));
        // admissibility on a 0.025 lattice over Ab
        const int AX = 169, AR = 125;
        Vec ll(U(AX * AR)), bm(U(AX * AR));
        std::vector<char> val(U(AX * AR));
#pragma omp parallel for schedule(dynamic, 4)
        for (int i = 0; i < AX; ++i)
            for (int j = 0; j < AR; ++j) {
                const double x = -kAbX + 0.025 * i, r = 0.025 * j;
                const AdmFieldJets<double> J = adm_field_jets<double>(c.f2a, c.cfg, c.sys.parent, f, 0.0, x, r, 0.0);
                ll[U(i * AR + j)] = J.log_lapse.v;
                bm[U(i * AR + j)] = std::sqrt(J.beta[0].v * J.beta[0].v + J.beta[1].v * J.beta[1].v +
                                              J.beta[2].v * J.beta[2].v);
                val[U(i * AR + j)] =
                    adm_metric_jet<double>(c.f2a, c.cfg, c.sys.parent, f, 0.0, x, r, 0.0).valid ? 1 : 0;
            }
        double llmin = 1e300, llmax = -1e300, bmax = 0;
        long nval = 0;
        for (std::size_t k = 0; k < ll.size(); ++k) {
            llmin = std::fmin(llmin, ll[k]); llmax = std::fmax(llmax, ll[k]); bmax = std::fmax(bmax, bm[k]);
            nval += val[k];
        }
        o << "," << E16 << "," << m16 << "," << Ro[0] << "," << Ro[1] << "," << Ro[2] << "," << llmin << "," << llmax
          << "," << bmax << "," << nval << "," << ll.size();
    } else {
        o << ",,,,,,,,,,";
    }
    o << "," << (now() - t0) << "\n";
    std::printf("%s seed %u: rho*_D %.8e at (%.5f,%.5f) | sample %.6e | est_min %.3e flagged %zu d7 %d | "
                "outside-Ab mismatches %ld/%ld | E8 %.6e | %.1fs\n",
                tag.c_str(), seed, S.best.v, S.best.x, S.best.r, S.sample_min, S.min_est, S.flagged, S.d7 ? 1 : 0,
                mism, chk, d.Eneg, now() - t0);
    return 0;
}

int mode_floor(const Ctx& c) {
    const int n = 3 * c.block;
    const AdmFieldCoefficients<double> f0 = fields_of<double>(c, Vec(U(n), 0.0));
    const BoxSearch S = box_search(c, f0, -7.0, 7.0, 0.0, 7.0, 0.005, 1, 1 << 17, [](double x, double r) {
        return std::fabs(x) <= 7.0 && r >= 0 && r <= 7.0 && !in_Ab(x, r);
    });
    std::printf("FLOOR (baseline min over D \\ Ab): %.10e at (%.5f, %.5f); sample %.10e; starts %zu\n", S.best.v,
                S.best.x, S.best.r, S.sample_min, S.starts);
    return 0;
}

int mode_axi(const Ctx& c, const std::string& yfile) {
    const int n = 3 * c.block;
    Vec yv = yfile == "BASELINE" ? Vec(U(n), 0.0) : read_vec(yfile);
    const AdmFieldCoefficients<double> f = fields_of<double>(c, yv);
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> ux(-7, 7), ur(0, 7), up(0, 6.283185307179586);
    const int N = 2000;
    std::vector<double> d(N), m(N), X(N), R(N), P(N);
    for (int k = 0; k < N; ++k) { X[U(k)] = ux(rng); R[U(k)] = ur(rng); P[U(k)] = up(rng); }
#pragma omp parallel for
    for (int k = 0; k < N; ++k) {
        const double x = X[U(k)], r = R[U(k)], ph = P[U(k)];
        const double a = rho_pt(c, f, x, r * std::cos(ph), r * std::sin(ph));
        const double b = rho_pt(c, f, x, r, 0.0);
        d[U(k)] = std::fabs(a - b);
        m[U(k)] = std::fabs(b);
    }
    double dmax = 0, mmax = 0;
    for (int k = 0; k < N; ++k) { dmax = std::fmax(dmax, d[U(k)]); mmax = std::fmax(mmax, m[U(k)]); }
    const bool pass = dmax <= 1e-9 * mmax + 1e-14;
    std::printf("AXI %s: max|rho3D - rho_meridian| %.3e  max|rho| %.3e  ratio %.3e  %s\n", yfile.c_str(), dmax, mmax,
                dmax / mmax, pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}

int mode_field(const Ctx& c, const std::string& yfile, const std::string& out) {
    const int n = 3 * c.block;
    Vec yv = yfile == "BASELINE" ? Vec(U(n), 0.0) : read_vec(yfile);
    const AdmFieldCoefficients<double> f = fields_of<double>(c, yv);
    const int NX = 421, NR = 311;  // Ab at 0.01
    std::vector<double> v(U(NX * NR));
#pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NR; ++j) v[U(i * NR + j)] = rho_pt(c, f, -kAbX + 0.01 * i, 0.01 * j);
    std::ofstream o(out);
    o.precision(9);
    o << "# x in [-2.1,2.1] " << NX << " rows; r in [0,3.1] " << NR << " cols\n";
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NR; ++j) o << v[U(i * NR + j)] << (j + 1 < NR ? " " : "\n");
    return 0;
}

int mode_timing(const Ctx& c) {
    const int n = 3 * c.block;
    const Vec y(U(n), 0.0);
    const AdmFieldCoefficients<double> f = fields_of<double>(c, y);
    for (int ord : {4, 6}) {
        const NodeSet ns = base_set(c, ord, 1, 1.7, 2.6);
        double t = now();
        const Vec g = gradient(c, y, ns);
        std::printf("order %d: gradient %.2fs (|g| %.6e)\n", ord, now() - t, nrm(g));
        t = now();
        const double R = residual(c, y, ns);
        std::printf("order %d: residual %.2fs (R %.12e)\n", ord, now() - t, R);
    }
    double t = now();
    const InLoop in = inloop_search(c, f, 3.4, 2.0 * std::sqrt(2.6));
    std::printf("in-loop search %.2fs: %zu distinct, %d polished, %d dup, %d lattice minima, best %.8e at (%.4f,%.4f)\n",
                now() - t, in.mins.size(), in.polished, in.duplicates, in.lattice_minima, in.mins[0].v, in.mins[0].x,
                in.mins[0].r);
    const Val V = make_val(c);
    t = now();
    const Diag d = validation_diag(c, y, V.vs, V.eps_rho);
    std::printf("8:1 diag %.2fs (E8 %.9e, EnegW %.6e, vtail %.6e, L4 %.6e)\n", now() - t, d.Eneg, d.EnegW, d.vtail, d.L4);
    return 0;
}

int mode_check(const Ctx& c) {
    int bad = 0;
    for (int ord : {4, 6, 8}) {
        const NodeSet ns = base_set(c, ord, 1, ord == 8 ? 3.5 : 1.7, ord == 8 ? 12.25 : 2.6);
        const IntegrationSet s =
            build_elementwise_set(c.sys.parent, ord, 1, ord == 8 ? 3.5 : 1.7, ord == 8 ? 12.25 : 2.6, 0.0, c.cfg.velocity);
        bool same = s.x.size() == ns.x.size();
        for (std::size_t p = 0; same && p < s.x.size(); ++p)
            same = s.x[p] == ns.x[p] && s.y[p] == ns.y[p] && s.weight[p] == ns.w[p];
        std::printf("base set order %d: %zu nodes, bit-identical to build_elementwise_set: %s\n", ord, ns.x.size(),
                    same ? "yes" : "NO");
        if (!same) ++bad;
    }
    const int n = 3 * c.block;
    const Val V = make_val(c);
    for (auto& st : starts_of(n, c.block)) {
        if (st.first != "combined_perturbation" && st.first != "baseline_zero") continue;
        const Diag d = validation_diag(c, st.second, V.vs, 4.64e-12);
        const FeasibilityReport r = evaluate_feasibility(c.f2a, c.cfg, c.sys.parent, fields_of<double>(c, st.second),
                                                         V.vs, c.tau, 0.0, SourceSemantics::TotalEffective);
        std::printf("%s: E_neg parallel %.15e vs evaluate_feasibility %.15e; rho_min %.15e vs %.15e\n",
                    st.first.c_str(), d.Eneg, r.E_neg, d.rhomin, r.rho_min);
        if (std::fabs(d.Eneg - r.E_neg) > 1e-12 * r.E_neg || d.rhomin != r.rho_min) ++bad;
        const NodeSet w4 = base_set(c, 4, 1, 1.7, 2.6);
        const IntegrationSet s4 = build_elementwise_set(c.sys.parent, 4, 1, 1.7, 2.6, 0.0, c.cfg.velocity);
        const double r1 = residual(c, st.second, w4);
        const double r2 = feasibility_residual<double>(c.f2a, c.cfg, c.sys.parent, fields_of<double>(c, st.second), s4,
                                                       c.tau, 0.0, SourceSemantics::TotalEffective);
        std::printf("   R(4:1) %.17e vs feasibility_residual %.17e %s\n", r1, r2, r1 == r2 ? "identical" : "DIFF");
        if (r1 != r2) ++bad;
    }
    std::printf("CHECK %s\n", bad ? "FAIL" : "PASS");
    return bad;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: s3 run|traj|assess|floor|axi|field|timing|check ...\n"); return 2; }
    const std::string mode = argv[1];
    std::string out = ".", start = "baseline_zero", yfile = "BASELINE", tag = "t", ysfile, its;
    char form = 'A';
    int order = 4, sub = 1, kmax = 600, maxfev = 3000, diag_every = 1;
    double lam = 0.0, deadline = 0.0;
    unsigned seed = 1;
    bool full = false;
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        auto nx = [&]() { return std::string(argv[++i]); };
        if (a == "--out") out = nx();
        else if (a == "--form") form = nx()[0];
        else if (a == "--ws") { const std::string s = nx(); order = std::stoi(s.substr(0, s.find(':'))); sub = std::stoi(s.substr(s.find(':') + 1)); }
        else if (a == "--start") start = nx();
        else if (a == "--kmax") kmax = std::stoi(nx());
        else if (a == "--maxfev") maxfev = std::stoi(nx());
        else if (a == "--lambda") lam = std::stod(nx());
        else if (a == "--deadline") deadline = std::stod(nx());
        else if (a == "--diag-every") diag_every = std::stoi(nx());
        else if (a == "--y") yfile = nx();
        else if (a == "--ys") ysfile = nx();
        else if (a == "--its") its = nx();
        else if (a == "--seed") seed = static_cast<unsigned>(std::stoul(nx()));
        else if (a == "--tag") tag = nx();
        else if (a == "--full") full = nx() == "1";
    }
    const Ctx c = make_ctx();
    if (mode == "run") return mode_run(c, form, order, sub, start, lam, kmax, maxfev, deadline, diag_every, out);
    if (mode == "traj") return mode_traj(c, ysfile, its, out);
    if (mode == "assess") return mode_assess(c, yfile, seed, out, tag, full);
    if (mode == "floor") return mode_floor(c);
    if (mode == "axi") return mode_axi(c, yfile);
    if (mode == "field") return mode_field(c, yfile, out);
    if (mode == "timing") return mode_timing(c);
    if (mode == "check") return mode_check(c);
    return 2;
}
