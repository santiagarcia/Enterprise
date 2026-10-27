// s2: study-2 driver (PROTOCOL.md). Builds against the UNMODIFIED Enterprise/warp_hypad
// headers + libraries at 477141d (same flags as study 1's wh_probe).
//
// Modes
//   run     --form A|B --ws ORDER:SUB --start NAME --kmax K --out DIR
//   assess  --y FILE|BASELINE --seed S --out CSV --tag TAG
//   axi     --y FILE|BASELINE
//   field   --y FILE|BASELINE --out FILE
//   check   (base-set identity vs build_elementwise_set; E_neg vs evaluate_feasibility)
//
// Per-node values are computed in parallel and summed serially in node order.

#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
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
// Element-based node sets (leaf rectangles in (xi, s), n x n Gauss each).
struct Elem { double xa, xb, sa, sb; int depth; };

struct NodeSet {
    int order = 4;
    std::vector<Elem> elems;
    Vec x, y, w;            // nodes (z = 0), weights
    std::vector<int> elem_of;
    std::size_t base_nodes = 0;
};

// Copy of build_elementwise_set's edge construction (verified bit-identical in `check`).
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
// Local polish: compass search (8 directions), step halving.
struct Pt { double x, r, v; };

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

// In-loop search on Dw (PROTOCOL 3): 241 x 116 lattice, polish 8 best distinct (sep >= 0.1).
std::vector<Pt> inloop_search(const Ctx& c, const AdmFieldCoefficients<double>& f, double xw, double rw) {
    const int NX = 241, NR = 116;
    std::vector<Pt> lat(U(NX * NR));
#pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NR; ++j) {
            const double x = -xw + 2.0 * xw * i / (NX - 1), r = rw * j / (NR - 1);
            lat[U(i * NR + j)] = {x, r, rho_pt(c, f, x, r)};
        }
    std::vector<Pt> st = pick_distinct(lat, 8, 0.1);
    std::vector<Pt> res(st.size());
    const double h0 = 2.0 * xw / (NX - 1);
    auto fr = [&](double x, double r) { return rho_pt(c, f, x, r); };
    auto ok = [&](double x, double r) { return x >= -xw && x <= xw && r >= 0.0 && r <= rw; };
#pragma omp parallel for schedule(dynamic, 1)
    for (int k = 0; k < static_cast<int>(st.size()); ++k)
        res[U(k)] = polish(fr, ok, st[U(k)].x, st[U(k)].r, st[U(k)].v, h0, 1e-9);
    std::sort(res.begin(), res.end(), [](const Pt& a, const Pt& b) { return a.v < b.v; });
    return res;
}

// B refinement (PROTOCOL 3). Returns number of splits.
int refine_B(const Ctx& c, const AdmFieldCoefficients<double>& f, NodeSet& ns, const std::vector<Pt>& cand,
             std::ofstream& ev, int it, bool& capped) {
    const double R0 = c.sys.parent.scale_R0();
    int splits = 0;
    for (const Pt& p : cand) {
        if (!(p.v < 0.0)) continue;
        const double xi = p.x / R0, s = (p.r / R0) * (p.r / R0);
        for (int guard = 0; guard < 16; ++guard) {
            int L = -1;
            for (std::size_t e = 0; e < ns.elems.size(); ++e) {
                const Elem& E = ns.elems[e];
                if (xi >= E.xa && xi <= E.xb && s >= E.sa && s <= E.sb) { L = static_cast<int>(e); break; }
            }
            if (L < 0) break;  // outside Dw
            double rhoL = 1e300;
            for (std::size_t q = 0; q < ns.x.size(); ++q)
                if (ns.elem_of[q] == L) rhoL = std::fmin(rhoL, rho_pt(c, f, ns.x[q], ns.y[q]));
            const Elem E = ns.elems[U(L)];
            const bool need = p.v < rhoL - 0.05 * std::fabs(p.v);
            ev << it << "," << p.x << "," << p.r << "," << p.v << "," << rhoL << "," << E.depth << "," << (need ? 1 : 0);
            if (!need) { ev << ",resolved\n"; break; }
            if (E.depth >= 6) { ev << ",max_depth\n"; break; }
            if (ns.x.size() + 3u * static_cast<std::size_t>(ns.order * ns.order) > 4u * ns.base_nodes) {
                ev << ",node_cap\n"; capped = true; break;
            }
            ev << ",split\n";
            const double xm = 0.5 * (E.xa + E.xb), sm = 0.5 * (E.sa + E.sb);
            const int d = E.depth + 1;
            ns.elems[U(L)] = {E.xa, xm, E.sa, sm, d};
            ns.elems.push_back({xm, E.xb, E.sa, sm, d});
            ns.elems.push_back({E.xa, xm, sm, E.sb, d});
            ns.elems.push_back({xm, E.xb, sm, E.sb, d});
            rebuild_nodes(c, ns);
            ++splits;
        }
    }
    ev.flush();
    return splits;
}

// Validation-set diagnostics, parallel evaluation, serial sums. gamma = I => dV = w.
struct Diag { double Eneg, negfrac, rhomin, rx, ry, refined; };
Diag validation_diag(const Ctx& c, const Vec& yv, const IntegrationSet& vs, double eps_rho) {
    const AdmFieldCoefficients<double> f = fields_of<double>(c, yv);
    const Vec r = node_rho(c, f, vs.x, vs.y);
    Diag d{0, 0, 1e300, 0, 0, 0};
    double vneg = 0, vtot = 0;
    for (std::size_t p = 0; p < r.size(); ++p) {
        vtot += vs.weight[p];
        if (r[p] < 0) { d.Eneg += -r[p] * vs.weight[p]; vneg += vs.weight[p]; }
        if (r[p] < d.rhomin) { d.rhomin = r[p]; d.rx = vs.x[p]; d.ry = vs.y[p]; }
    }
    d.negfrac = vneg / vtot;
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

// ---------------------------------------------------------------------------
int mode_run(const Ctx& c, char form, int order, int sub, const std::string& start, int kmax, int maxfev,
             const std::string& out) {
    const double t_start = now();
    const int n = 3 * c.block;
    Vec y;
    for (auto& s : starts_of(n, c.block))
        if (s.first == start) y = s.second;
    if (y.empty()) { std::fprintf(stderr, "unknown start %s\n", start.c_str()); return 2; }
    NodeSet ns = base_set(c, order, sub, 1.7, 2.6);
    const IntegrationSet vs = build_elementwise_set(c.sys.parent, 8, 1, 3.5, 12.25, 0.0, c.cfg.velocity);
    double tv = 0; for (double q : vs.weight) tv += q;
    const double eps_E = std::fmax(o0d_state().numerical_uncertainty, 1e-9) * 10.0;
    const double eps_rho = eps_E / tv;
    const double xw = 1.7 * 2.0, rw = 2.0 * std::sqrt(2.6);
    {
        std::ofstream cf(out + "/config.txt");
        cf << "form " << form << "\norder " << order << "\nsub " << sub << "\nstart " << start << "\nkmax " << kmax
           << "\nmaxfev " << maxfev << "\nbase_nodes " << ns.base_nodes << "\ntau " << c.tau
           << "\nthreads " << omp_get_max_threads() << "\nC1 1e-3 rel ||g0||\nC2 5-window 1e-6\nstall 20-window 1e-6\n"
           << "lbfgs m=8 c1=1e-4 trust0=0.5 trustmax=4 trustmin=1e-8 ls_trials=10\n";
    }
    std::ofstream tr(out + "/trace.csv");
    tr.precision(12);
    tr << "it,f,gnorm,gnorm_rel,step,trust,nfev,ngev,nodes,leaves,splits,Eneg_v8,negfrac_v8,rhomin_w,rhomin_v8,"
          "refined_v8,rhostar_in,xstar,rstar,C1,C2,C3,secs\n";
    std::ofstream ys(out + "/ys.txt");
    ys.precision(17);
    std::ofstream ev(out + "/refine_events.csv");
    ev.precision(12);
    ev << "it,x,r,rho_star,rho_L,depth,need,action\n";
    std::ofstream ls(out + "/linesearch.csv");
    ls.precision(12);
    ls << "it,trial,alpha,trust,f_trial,f_cur,armijo,ok\n";

    int nfev = 0, ngev = 0;
    bool capped = false;
    auto search = [&]() { return inloop_search(c, fields_of<double>(c, y), xw, rw); };
    auto log_iter = [&](int it, double f, const Vec& g, double g0n, double step, double trust, int splits, bool C1,
                        bool C2, bool C3, const std::vector<Pt>& in) {
        const AdmFieldCoefficients<double> fd = fields_of<double>(c, y);
        const Vec rw_ = node_rho(c, fd, ns.x, ns.y);
        double rmw = 1e300; for (double v : rw_) rmw = std::fmin(rmw, v);
        const Diag d = validation_diag(c, y, vs, eps_rho);
        tr << it << "," << f << "," << nrm(g) << "," << nrm(g) / g0n << "," << step << "," << trust << "," << nfev
           << "," << ngev << "," << ns.x.size() << "," << ns.elems.size() << "," << splits << "," << d.Eneg << ","
           << d.negfrac << "," << rmw << "," << d.rhomin << "," << d.refined << "," << in[0].v << "," << in[0].x
           << "," << in[0].r << "," << C1 << "," << C2 << "," << C3 << "," << (now() - t_start) << "\n";
        tr.flush();
        for (double v : y) ys << v << " ";
        ys << "\n";
        ys.flush();
    };

    // g0 on the BASE set (PROTOCOL 2).
    double f = residual(c, y, ns); ++nfev;
    Vec g = gradient(c, y, ns); ++ngev;
    const double g0n = nrm(g);
    std::vector<Pt> cand = search();
    int last_change = 0;
    std::vector<double> fh;  // f history since last set change
    int sp0 = 0;
    if (form == 'B') {
        sp0 = refine_B(c, fields_of<double>(c, y), ns, cand, ev, 0, capped);
        if (sp0 > 0) { f = residual(c, y, ns); ++nfev; g = gradient(c, y, ns); ++ngev; }
    }
    log_iter(0, f, g, g0n, 0.0, 0.5, sp0, false, false, sp0 == 0, cand);
    fh.push_back(f);
    double trust = 0.5;
    std::deque<Vec> S, Y;
    std::string term;
    const double kC1 = 1e-4;
    int it = 0;
    while (true) {
        if (it >= kmax) { term = "ITERATION_LIMIT"; break; }
        if (nfev >= maxfev) { term = "EVAL_LIMIT"; break; }
        Vec p(U(n));
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
        }
        double dd = dot(p, g);
        if (!(dd < 0.0)) {
            S.clear(); Y.clear();
            for (int i = 0; i < n; ++i) p[U(i)] = -g[U(i)];
            dd = dot(p, g);
        }
        const double pn = nrm(p);
        double alpha = pn > 0.0 ? std::fmin(1.0, trust / pn) : 0.0;
        bool ok = false;
        Vec yn;
        double fn = 0.0, step = 0.0;
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
                fn = residual(c, yn, ns); ++nfev;
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
        const Vec g2 = gradient(c, yn, ns); ++ngev;
        bool finite = std::isfinite(fn);
        for (double q : g2) if (!std::isfinite(q)) finite = false;
        if (!finite) { term = "INVALID"; break; }
        Vec sv(U(n)), yv(U(n));
        for (int i = 0; i < n; ++i) { sv[U(i)] = yn[U(i)] - y[U(i)]; yv[U(i)] = g2[U(i)] - g[U(i)]; }
        if (dot(sv, yv) > 1e-14 * nrm(sv) * nrm(yv)) {
            S.push_back(sv); Y.push_back(yv);
            while (S.size() > 8) { S.pop_front(); Y.pop_front(); }
        }
        y = yn; g = g2; f = fn;
        ++it;
        if (att == 0) trust = std::fmin(trust * 2.0, 4.0);
        fh.push_back(f);
        // criteria on the current set (before any refinement at this iterate)
        const std::size_t m = fh.size() - 1;
        const bool C1 = nrm(g) <= 1e-3 * g0n;
        const bool C2 = m >= 5 && (fh[m - 5] - fh[m]) / fh[m] <= 1e-6;
        const bool stall = m >= 20 && (fh[m - 20] - fh[m]) / fh[m] <= 1e-6 && !C1;
        // in-loop search (diagnostic for A, refinement driver for B), logged with this iterate
        int splits = 0;
        cand = search();
        const double f_acc = f;
        const Vec g_acc = g;
        if (form == 'B' && !capped) {
            splits = refine_B(c, fields_of<double>(c, y), ns, cand, ev, it, capped);
            if (splits > 0) {
                f = residual(c, y, ns); ++nfev;
                g = gradient(c, y, ns); ++ngev;
                S.clear(); Y.clear();
                fh.clear(); fh.push_back(f);
                last_change = it;
                std::ofstream(out + "/set_changes.csv", std::ios::app)
                    << it << "," << splits << "," << ns.x.size() << "," << f << "," << nrm(g) << "\n";
            }
        }
        const bool C3 = (form != 'B') || splits == 0;
        // trace row: f, |g| of the accepted iterate on the set it was accepted on; node columns on the current set
        log_iter(it, f_acc, g_acc, g0n, step, trust, splits, C1, C2, C3, cand);
        std::printf("it %3d f %.8e |g|/g0 %.3e step %.3e trust %.3g nodes %zu splits %d C1 %d C2 %d C3 %d  %.0fs\n", it,
                    f, nrm(g) / g0n, step, trust, ns.x.size(), splits, C1, C2, C3, now() - t_start);
        std::fflush(stdout);
        if (C1 && C2 && C3) { term = "CONVERGED"; break; }
        if (stall && splits == 0) { term = "STALLED"; break; }
    }
    if (term == "LINESEARCH_FAIL") {
        const std::size_t m = fh.size() - 1;
        const bool C1 = nrm(g) <= 1e-3 * g0n;
        const bool C2 = m >= 5 && (fh[m - 5] - fh[m]) / fh[m] <= 1e-6;
        if (C1 && C2) term = "CONVERGED_LINESEARCH_LIMITED";
    }
    {
        std::ofstream fy(out + "/y_final.txt");
        fy.precision(17);
        for (double v : y) fy << v << "\n";
        std::ofstream ft(out + "/termination.txt");
        ft << term << " it=" << it << " nfev=" << nfev << " ngev=" << ngev << " nodes=" << ns.x.size()
           << " capped=" << capped << " last_set_change=" << last_change << " gnorm_rel=" << nrm(g) / g0n
           << " f=" << f << " wall_s=" << (now() - t_start) << "\n";
    }
    std::printf("TERMINATION %s it=%d nfev=%d wall=%.0fs\n", term.c_str(), it, nfev, now() - t_start);
    return 0;
}

// ---------------------------------------------------------------------------
double halton(unsigned i, unsigned b) {
    double f = 1.0, r = 0.0;
    while (i > 0) { f /= b; r += f * (i % b); i /= b; }
    return r;
}

int mode_assess(const Ctx& c, const std::string& yfile, unsigned seed, const std::string& out, const std::string& tag) {
    const double t0 = now();
    const int n = 3 * c.block;
    Vec yv = yfile == "BASELINE" ? Vec(U(n), 0.0) : read_vec(yfile);
    const AdmFieldCoefficients<double> f = fields_of<double>(c, yv);
    const AdmFieldCoefficients<double> f0 = fields_of<double>(c, Vec(U(n), 0.0));
    Vec yr(U(n));
    for (int i = 0; i < n; ++i) yr[U(i)] = 0.3 * std::sin(1.7 * i + 0.3) + 0.1;
    const AdmFieldCoefficients<double> frand = fields_of<double>(c, yr);
    auto fr = [&](double x, double r) { return rho_pt(c, f, x, r); };
    auto active = [&](double x, double r) { return rho_pt(c, frand, x, r) != rho_pt(c, f0, x, r); };
    const double XD = 7.0, RD = 7.0, XW = 3.4, RW = 2.0 * std::sqrt(2.6);
    std::vector<Pt> all;
    // 1a lattice 1401 x 701
    const int NX = 1401, NR = 701;
    const double dx = 2 * XD / (NX - 1), dr = RD / (NR - 1);
    std::vector<double> L(U(NX * NR));
    std::vector<char> act(U(NX * NR));
#pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NR; ++j) {
            const double x = -XD + dx * i, r = dr * j;
            L[U(i * NR + j)] = fr(x, r);
            act[U(i * NR + j)] = active(x, r) ? 1 : 0;
        }
    // 1b Halton 2^16, Cranley-Patterson shift from seed
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    const double sh1 = u01(rng), sh2 = u01(rng);
    const int NH = 1 << 16;
    std::vector<Pt> H(U(NH));
#pragma omp parallel for schedule(dynamic, 64)
    for (int k = 0; k < NH; ++k) {
        const double a = std::fmod(halton(static_cast<unsigned>(k + 1), 2) + sh1, 1.0);
        const double b = std::fmod(halton(static_cast<unsigned>(k + 1), 3) + sh2, 1.0);
        const double x = -XD + 2 * XD * a, r = RD * b;
        H[U(k)] = {x, r, fr(x, r)};
    }
    double gmin = 1e300;
    for (double v : L) gmin = std::fmin(gmin, v);
    for (const Pt& p : H) gmin = std::fmin(gmin, p.v);
    // 2 adaptive subdivision
    struct Cell { double x, r, hx, hr, m; };
    std::vector<Cell> cells;
    for (int i = 0; i + 1 < NX; ++i)
        for (int j = 0; j + 1 < NR; ++j) {
            const double m = std::fmin(std::fmin(L[U(i * NR + j)], L[U((i + 1) * NR + j)]),
                                       std::fmin(L[U(i * NR + j + 1)], L[U((i + 1) * NR + j + 1)]));
            if (m <= gmin + 0.1 * std::fabs(gmin)) cells.push_back({-XD + dx * i, dr * j, dx, dr, m});
        }
    std::sort(cells.begin(), cells.end(), [](const Cell& a, const Cell& b) { return a.m < b.m; });
    const std::size_t ncells0 = cells.size();
    if (cells.size() > 4096) cells.resize(4096);
    std::vector<Pt> sub;
    for (int level = 0; level < 3; ++level) {
        std::vector<Cell> next;
        std::vector<Pt> pts(cells.size() * 64);
        const long nc = static_cast<long>(cells.size());
#pragma omp parallel for schedule(dynamic, 4)
        for (long k = 0; k < nc; ++k)
            for (int a = 0; a < 8; ++a)
                for (int b = 0; b < 8; ++b) {
                    const Cell& C = cells[U(k)];
                    const double x = C.x + C.hx * (a + 0.5) / 8, r = C.r + C.hr * (b + 0.5) / 8;
                    pts[U(k * 64 + a * 8 + b)] = {x, r, fr(x, r)};
                }
        for (long k = 0; k < nc; ++k)
            for (int q = 0; q < 64; ++q) {
                const Pt& p = pts[U(k * 64 + q)];
                next.push_back({p.x - cells[U(k)].hx / 16, p.r - cells[U(k)].hr / 16, cells[U(k)].hx / 8,
                                cells[U(k)].hr / 8, p.v});
                sub.push_back(p);
            }
        std::sort(next.begin(), next.end(), [](const Cell& a, const Cell& b) { return a.m < b.m; });
        if (next.size() > 32) next.resize(32);
        cells = next;
    }
    // candidate pool
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NR; ++j) all.push_back({-XD + dx * i, dr * j, L[U(i * NR + j)]});
    for (const Pt& p : H) all.push_back(p);
    for (const Pt& p : sub) all.push_back(p);
    double sample_min = 1e300;
    for (const Pt& p : all) sample_min = std::fmin(sample_min, p.v);
    // 3 multi-start polish on D
    auto okD = [&](double x, double r) { return x >= -XD && x <= XD && r >= 0 && r <= RD; };
    auto okW = [&](double x, double r) { return x >= -XW && x <= XW && r >= 0 && r <= RW; };
    auto run_polish = [&](const std::vector<Pt>& st, auto ok, auto fun) {
        std::vector<Pt> res(st.size());
#pragma omp parallel for schedule(dynamic, 1)
        for (int k = 0; k < static_cast<int>(st.size()); ++k)
            res[U(k)] = polish(fun, ok, st[U(k)].x, st[U(k)].r, st[U(k)].v, 0.01, 1e-10);
        std::sort(res.begin(), res.end(), [](const Pt& a, const Pt& b) { return a.v < b.v; });
        return res;
    };
    const std::vector<Pt> stD = pick_distinct(all, 64, 0.05);
    const std::vector<Pt> pD = run_polish(stD, okD, fr);
    std::vector<Pt> inW;
    for (const Pt& p : all) if (okW(p.x, p.r)) inW.push_back(p);
    const std::vector<Pt> pW = run_polish(pick_distinct(inW, 64, 0.05), okW, fr);
    // active region / design-independent floor (lattice-based candidates)
    std::vector<Pt> la, li;
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NR; ++j) {
            const Pt p{-XD + dx * i, dr * j, L[U(i * NR + j)]};
            (act[U(i * NR + j)] ? la : li).push_back(p);
        }
    auto okA = [&](double x, double r) { return okD(x, r) && active(x, r); };
    auto okI = [&](double x, double r) { return okD(x, r) && !active(x, r); };
    auto fb = [&](double x, double r) { return rho_pt(c, f0, x, r); };
    const std::vector<Pt> pA = la.empty() ? std::vector<Pt>{} : run_polish(pick_distinct(la, 16, 0.05), okA, fr);
    std::vector<Pt> liB;
    for (const Pt& p : li) liB.push_back({p.x, p.r, fb(p.x, p.r)});
    const std::vector<Pt> pI = liB.empty() ? std::vector<Pt>{} : run_polish(pick_distinct(liB, 16, 0.05), okI, fb);
    // active bounding box
    double axl = 1e9, axh = -1e9, arh = 0;
    for (const Pt& p : la) { axl = std::fmin(axl, p.x); axh = std::fmax(axh, p.x); arh = std::fmax(arh, p.r); }
    // E_neg on 8:1 and 16:1 (validation extent)
    auto Eneg = [&](int order, int subd, double& negfrac, double& nodemin) {
        const IntegrationSet s = build_elementwise_set(c.sys.parent, order, subd, 3.5, 12.25, 0.0, c.cfg.velocity);
        const Vec r = node_rho(c, f, s.x, s.y);
        double E = 0, vn = 0, vt = 0;
        nodemin = 1e300;
        for (std::size_t p = 0; p < r.size(); ++p) {
            vt += s.weight[p];
            if (r[p] < 0) { E += -r[p] * s.weight[p]; vn += s.weight[p]; }
            nodemin = std::fmin(nodemin, r[p]);
        }
        negfrac = vn / vt;
        return E;
    };
    double nf8, nm8, nf16, nm16;
    const double E8 = Eneg(8, 1, nf8, nm8);
    const double E16 = Eneg(16, 1, nf16, nm16);
    std::ofstream o(out, std::ios::app);
    o.precision(12);
    o << tag << "," << yfile << "," << seed << "," << pD[0].v << "," << pD[0].x << "," << pD[0].r << "," << pW[0].v
      << "," << pW[0].x << "," << pW[0].r << "," << (pA.empty() ? NAN : pA[0].v) << ","
      << (pA.empty() ? NAN : pA[0].x) << "," << (pA.empty() ? NAN : pA[0].r) << "," << (pI.empty() ? NAN : pI[0].v)
      << "," << (pI.empty() ? NAN : pI[0].x) << "," << (pI.empty() ? NAN : pI[0].r) << "," << sample_min << ","
      << ncells0 << "," << E8 << "," << nf8 << "," << nm8 << "," << E16 << "," << nf16 << "," << nm16 << "," << axl
      << "," << axh << "," << arh << "," << (now() - t0) << "\n";
    // top local minima for the map
    std::ofstream om(out + "." + tag + ".minima.csv");
    om.precision(12);
    om << "x,r,rho\n";
    for (const Pt& p : pick_distinct(pD, 64, 1e-3)) om << p.x << "," << p.r << "," << p.v << "\n";
    std::printf("%s seed %u: rho*_D %.8e at (%.5f,%.5f) | rho*_Dw %.8e | active %.6e | ext floor %.6e | E8 %.6e E16 %.6e"
                " | cells %zu | %.1fs\n",
                tag.c_str(), seed, pD[0].v, pD[0].x, pD[0].r, pW[0].v, pA.empty() ? NAN : pA[0].v,
                pI.empty() ? NAN : pI[0].v, E8, E16, ncells0, now() - t0);
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
    const int NX = 701, NR = 351;
    std::vector<double> v(U(NX * NR));
#pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NR; ++j) v[U(i * NR + j)] = rho_pt(c, f, -7.0 + 14.0 * i / (NX - 1), 7.0 * j / (NR - 1));
    std::ofstream o(out);
    o.precision(9);
    o << "# x in [-7,7] " << NX << " rows; r in [0,7] " << NR << " cols\n";
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NR; ++j) o << v[U(i * NR + j)] << (j + 1 < NR ? " " : "\n");
    }
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
    const IntegrationSet vs = build_elementwise_set(c.sys.parent, 8, 1, 3.5, 12.25, 0.0, c.cfg.velocity);
    for (auto& st : starts_of(n, c.block)) {
        if (st.first != "combined_perturbation" && st.first != "baseline_zero") continue;
        const Diag d = validation_diag(c, st.second, vs, 4.64e-12);
        const FeasibilityReport r = evaluate_feasibility(c.f2a, c.cfg, c.sys.parent, fields_of<double>(c, st.second), vs,
                                                         c.tau, 0.0, SourceSemantics::TotalEffective);
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
    if (argc < 2) { std::fprintf(stderr, "usage: s2 run|assess|axi|field|check ...\n"); return 2; }
    const std::string mode = argv[1];
    std::string out = ".", start = "baseline_zero", yfile = "BASELINE", tag = "t";
    char form = 'A';
    int order = 4, sub = 1, kmax = 150, maxfev = 1500;
    unsigned seed = 1;
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        auto nx = [&]() { return std::string(argv[++i]); };
        if (a == "--out") out = nx();
        else if (a == "--form") form = nx()[0];
        else if (a == "--ws") { const std::string s = nx(); order = std::stoi(s.substr(0, s.find(':'))); sub = std::stoi(s.substr(s.find(':') + 1)); }
        else if (a == "--start") start = nx();
        else if (a == "--kmax") kmax = std::stoi(nx());
        else if (a == "--maxfev") maxfev = std::stoi(nx());
        else if (a == "--y") yfile = nx();
        else if (a == "--seed") seed = static_cast<unsigned>(std::stoul(nx()));
        else if (a == "--tag") tag = nx();
    }
    const Ctx c = make_ctx();
    if (mode == "run") return mode_run(c, form, order, sub, start, kmax, maxfev, out);
    if (mode == "assess") return mode_assess(c, yfile, seed, out, tag);
    if (mode == "axi") return mode_axi(c, yfile);
    if (mode == "field") return mode_field(c, yfile, out);
    if (mode == "check") return mode_check(c);
    return 2;
}
