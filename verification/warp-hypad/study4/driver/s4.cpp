// s4: study-4 driver (study4/PROTOCOL-v3.md). Builds against the UNMODIFIED Enterprise/warp_hypad headers +
// libraries at 477141d with the flags of study3/driver/build.sh.
//
// Assembled from study3/driver/s3.cpp: every block not listed below is copied verbatim (node sets, rho_pt,
// node_rho, residual, compass polish, in-loop search, 8:1 diagnostics, traj, Nelder-Mead, box_search, assess,
// floor, axi, field, timing, check). Changed/new blocks (PROTOCOL-v3 sections 1-4, 9):
//   * gradient(): forward sum unchanged (bit-identical, gate V2); per-node contributions are also summed in
//     reverse order, with Neumaier compensation and in long double (H1b/H1c, v3 N3).
//   * refine_B(): split decisions verbatim; templated on the rho evaluator; returns an outcome record; no
//     persistent `capped` flag (v2 R1).
//   * run: dense BFGS initialised from the measured Hessian, refresh rule, B nested loop (sections 3, 4).
//   * hessian / fd / recon / v2dump / selftest modes (sections 2, 9).
//
// Per-node values are computed in parallel and summed serially in node order (deterministic).

#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <ctime>
#include <deque>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// ==== verbatim s3.cpp lines 32-190 ====
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

// ---------------------------------------------------------------------------
// Summation variants of stored per-node contributions (PROTOCOL-v3 2.1, 2.2 H1b; v3 N3).
// fwd: the value used (OTI operator+ in node order, exactly the study-3 accumulation); rev: reverse node order;
// comp: Neumaier compensated, coefficient by coefficient; ld: long double forward; abs: forward sum of |x|.
constexpr int kNC = 10;
static_assert(sizeof(onumm3n2_t) == kNC * sizeof(double), "onumm3n2_t must be 10 packed doubles");
constexpr double kU = 1.1102230246251565e-16;  // 2^-53

inline void to_arr(const OTIScalar& v, double* a) { std::memcpy(a, &v.raw(), sizeof(onumm3n2_t)); }
inline OTIScalar from_arr(const double* a) {
    OTIScalar v;
    std::memcpy(&v.raw(), a, sizeof(onumm3n2_t));
    return v;
}

struct Neumaier {
    double s = 0.0, cc = 0.0;
    void add(double x) {
        const double t = s + x;
        if (std::fabs(s) >= std::fabs(x)) cc += (s - t) + x;
        else cc += (x - t) + s;
        s = t;
    }
    double value() const { return s + cc; }
};

double neumaier_sum(const Vec& x) {
    Neumaier nm;
    for (double v : x) nm.add(v);
    return nm.value();
}

struct Sums {
    double fwd[kNC], rev[kNC], comp[kNC], ld[kNC], abs[kNC];
};

Sums sum_variants(const std::vector<OTIScalar>& contrib) {
    Sums S{};
    const std::size_t np = contrib.size();
    OTIScalar acc = scalar_cast<OTIScalar>(0.0);
    for (std::size_t p = 0; p < np; ++p) acc = acc + contrib[p];
    to_arr(acc, S.fwd);
    OTIScalar acr = scalar_cast<OTIScalar>(0.0);
    for (std::size_t p = np; p-- > 0;) acr = acr + contrib[p];
    to_arr(acr, S.rev);
    Neumaier nm[kNC];
    long double l[kNC];
    double a[kNC];
    for (int k = 0; k < kNC; ++k) { l[k] = 0.0L; a[k] = 0.0; }
    double x[kNC];
    for (std::size_t p = 0; p < np; ++p) {
        to_arr(contrib[p], x);
        for (int k = 0; k < kNC; ++k) {
            nm[k].add(x[k]);
            l[k] = l[k] + static_cast<long double>(x[k]);
            a[k] = a[k] + std::fabs(x[k]);
        }
    }
    for (int k = 0; k < kNC; ++k) {
        S.comp[k] = nm[k].value();
        S.ld[k] = static_cast<double>(l[k]);
        S.abs[k] = a[k];
    }
    return S;
}

// Per-node contribution evaluator for the warp objective: s(rho)^2 w with s = softplus_tau(-rho), in OTIScalar,
// exactly the expression of the study-3 gradient() loop.
struct WarpEval {
    const Ctx* c;
    const NodeSet* ns;
    void operator()(const std::vector<OTIScalar>& yy, std::vector<OTIScalar>& contrib) const {
        const AdmFieldCoefficients<OTIScalar> f = fields_of<OTIScalar>(*c, yy);
        const long np = static_cast<long>(ns->x.size());
        std::vector<OTIScalar> rho(U(np));
#pragma omp parallel for schedule(dynamic, 16)
        for (long p = 0; p < np; ++p)
            rho[U(p)] = family_eulerian_density<OTIScalar>(c->f2a, c->cfg, c->sys.parent, f, 0.0, ns->x[U(p)],
                                                           ns->y[U(p)], 0.0, 0.0, SourceSemantics::TotalEffective);
        contrib.assign(U(np), OTIScalar());
#pragma omp parallel for schedule(static)
        for (long p = 0; p < np; ++p) {
            const OTIScalar s = softplus(-rho[U(p)], c->tau);
            contrib[U(p)] = s * s * scalar_cast<OTIScalar>(ns->w[U(p)]);
        }
    }
};

// One OTI pass: seeds the listed coefficients (<= 3) and returns the summation variants.
template <typename Eval>
Sums run_pass(const Eval& ev, const Vec& yv, const std::vector<int>& seeds) {
    const std::size_t n = yv.size();
    std::vector<OTIScalar> yy(n);
    for (std::size_t j = 0; j < n; ++j) yy[j] = sensitivity::make_constant(yv[j]);
    for (std::size_t d = 0; d < seeds.size(); ++d)
        yy[U(seeds[d])] = sensitivity::seed_parameter(yv[U(seeds[d])], static_cast<int>(d));
    std::vector<OTIScalar> contrib;
    ev(yy, contrib);
    return sum_variants(contrib);
}

struct GradExt {
    Vec fwd, rev, comp, ld, bound;  // dR/dy variants and a-priori bound (n-1) u sum|contrib|
    long nodes = 0;
};

template <typename Eval>
GradExt gradient_ext_gen(const Eval& ev, const Vec& yv, long nodes) {
    const int n = static_cast<int>(yv.size());
    GradExt G;
    G.nodes = nodes;
    G.fwd.assign(U(n), 0.0); G.rev = G.fwd; G.comp = G.fwd; G.ld = G.fwd; G.bound = G.fwd;
    const double fac = static_cast<double>(nodes - 1) * kU;
    for (int s0 = 0; s0 < n; s0 += 3) {
        const int count = std::min(3, n - s0);
        std::vector<int> seeds;
        for (int d = 0; d < count; ++d) seeds.push_back(s0 + d);
        const Sums S = run_pass(ev, yv, seeds);
        for (int d = 0; d < count; ++d) {
            G.fwd[U(s0 + d)] = sensitivity::extract_gradient(from_arr(S.fwd), d);
            G.rev[U(s0 + d)] = sensitivity::extract_gradient(from_arr(S.rev), d);
            G.comp[U(s0 + d)] = sensitivity::extract_gradient(from_arr(S.comp), d);
            G.ld[U(s0 + d)] = sensitivity::extract_gradient(from_arr(S.ld), d);
            G.bound[U(s0 + d)] = fac * sensitivity::extract_gradient(from_arr(S.abs), d);
        }
    }
    return G;
}

GradExt gradient_ext(const Ctx& c, const Vec& yv, const NodeSet& ns) {
    return gradient_ext_gen(WarpEval{&c, &ns}, yv, static_cast<long>(ns.x.size()));
}

// Forward sum: bit-identical to study-3 gradient() (gate V2).
Vec gradient(const Ctx& c, const Vec& yv, const NodeSet& ns) { return gradient_ext(c, yv, ns).fwd; }

struct ResExt { double fwd, rev, comp, ld, bound; };
ResExt residual_ext(const Ctx& c, const Vec& yv, const NodeSet& ns) {
    const Vec r = node_rho(c, fields_of<double>(c, yv), ns.x, ns.y);
    Vec t(r.size());
    for (std::size_t p = 0; p < r.size(); ++p) {
        const double s = softplus(-r[p], c.tau);
        t[p] = s * s * ns.w[p];
    }
    ResExt o{0, 0, 0, 0, 0};
    double acc = 0.0, acr = 0.0, ab = 0.0;
    long double l = 0.0L;
    for (std::size_t p = 0; p < t.size(); ++p) {
        acc = acc + t[p];
        l += static_cast<long double>(t[p]);
        ab += std::fabs(t[p]);
    }
    for (std::size_t p = t.size(); p-- > 0;) acr = acr + t[p];
    o.fwd = acc; o.rev = acr; o.comp = neumaier_sum(t); o.ld = static_cast<double>(l);
    o.bound = static_cast<double>(t.size() - 1) * kU * ab;
    return o;
}

// ---------------------------------------------------------------------------
// Assembled OTI Hessian (PROTOCOL-v3 2.1): one pass per seed set of the registered schedule.
using Passes = std::vector<std::vector<int>>;

Passes read_passes(const std::string& path) {
    Passes P;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<int> s;
        std::stringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ',')) s.push_back(std::stoi(tok));
        P.push_back(s);
    }
    return P;
}

struct HessExt {
    int n = 0;
    long nodes = 0;
    Vec Hf, Hr, Hc, Hl, Bh;  // H_R variants (row-major n x n) and entrywise a-priori bound
    Vec gf;                  // gradient from the first pass containing i (fwd)
    long diag_mismatch = 0, grad_mismatch = 0, pair_errors = 0;
    double secs = 0.0;
};

template <typename Eval>
HessExt assemble_hessian(const Eval& ev, const Vec& yv, const Passes& P, long nodes) {
    const double t0 = now();
    const int n = static_cast<int>(yv.size());
    HessExt H;
    H.n = n;
    H.nodes = nodes;
    const std::size_t nn = U(n) * U(n);
    H.Hf.assign(nn, 0.0); H.Hr = H.Hf; H.Hc = H.Hf; H.Hl = H.Hf; H.Bh = H.Hf;
    H.gf.assign(U(n), 0.0);
    std::vector<int> cnt(nn, 0), seen(U(n), 0);
    const double fac = static_cast<double>(nodes - 1) * kU;
    for (const std::vector<int>& ps : P) {
        const Sums S = run_pass(ev, yv, ps);
        const OTIScalar F = from_arr(S.fwd), Rv = from_arr(S.rev), Cv = from_arr(S.comp), Lv = from_arr(S.ld),
                        Av = from_arr(S.abs);
        const int m = static_cast<int>(ps.size());
        for (int a = 0; a < m; ++a) {
            const int i = ps[U(a)];
            const double gi = sensitivity::extract_gradient(F, a);
            const double hii = sensitivity::extract_hessian(F, a, a);
            const std::size_t ii = U(i) * U(n) + U(i);
            if (!seen[U(i)]) {
                seen[U(i)] = 1;
                H.gf[U(i)] = gi;
                H.Hf[ii] = hii;
                H.Hr[ii] = sensitivity::extract_hessian(Rv, a, a);
                H.Hc[ii] = sensitivity::extract_hessian(Cv, a, a);
                H.Hl[ii] = sensitivity::extract_hessian(Lv, a, a);
                H.Bh[ii] = fac * sensitivity::extract_hessian(Av, a, a);
            } else {
                if (std::memcmp(&gi, &H.gf[U(i)], sizeof(double)) != 0) ++H.grad_mismatch;
                if (std::memcmp(&hii, &H.Hf[ii], sizeof(double)) != 0) ++H.diag_mismatch;
            }
            for (int b = a + 1; b < m; ++b) {
                const int j = ps[U(b)];
                const std::size_t ij = U(i) * U(n) + U(j), ji = U(j) * U(n) + U(i);
                ++cnt[ij];
                ++cnt[ji];
                H.Hf[ij] = H.Hf[ji] = sensitivity::extract_hessian(F, a, b);
                H.Hr[ij] = H.Hr[ji] = sensitivity::extract_hessian(Rv, a, b);
                H.Hc[ij] = H.Hc[ji] = sensitivity::extract_hessian(Cv, a, b);
                H.Hl[ij] = H.Hl[ji] = sensitivity::extract_hessian(Lv, a, b);
                H.Bh[ij] = H.Bh[ji] = fac * sensitivity::extract_hessian(Av, a, b);
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        if (!seen[U(i)]) ++H.pair_errors;
        for (int j = 0; j < n; ++j)
            if (i != j && cnt[U(i) * U(n) + U(j)] != 1) ++H.pair_errors;
    }
    H.secs = now() - t0;
    return H;
}

// Cyclic Jacobi eigensolver for a symmetric n x n matrix (row-major). Eigenvalues in ev; Q columns = vectors.
void jacobi_eig(Vec A, int n, Vec& ev, Vec& Q) {
    const std::size_t N = U(n);
    Q.assign(N * N, 0.0);
    for (std::size_t i = 0; i < N; ++i) Q[i * N + i] = 1.0;
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0, tot = 0.0;
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) {
                tot += A[i * N + j] * A[i * N + j];
                if (i != j) off += A[i * N + j] * A[i * N + j];
            }
        if (off <= 1e-30 * tot) break;
        for (std::size_t p = 0; p + 1 < N; ++p)
            for (std::size_t q = p + 1; q < N; ++q) {
                const double apq = A[p * N + q];
                if (apq == 0.0) continue;
                const double theta = (A[q * N + q] - A[p * N + p]) / (2.0 * apq);
                const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double cs = 1.0 / std::sqrt(t * t + 1.0), sn = t * cs;
                for (std::size_t k = 0; k < N; ++k) {
                    const double akp = A[k * N + p], akq = A[k * N + q];
                    A[k * N + p] = cs * akp - sn * akq;
                    A[k * N + q] = sn * akp + cs * akq;
                }
                for (std::size_t k = 0; k < N; ++k) {
                    const double apk = A[p * N + k], aqk = A[q * N + k];
                    A[p * N + k] = cs * apk - sn * aqk;
                    A[q * N + k] = sn * apk + cs * aqk;
                }
                for (std::size_t k = 0; k < N; ++k) {
                    const double qkp = Q[k * N + p], qkq = Q[k * N + q];
                    Q[k * N + p] = cs * qkp - sn * qkq;
                    Q[k * N + q] = sn * qkp + cs * qkq;
                }
            }
    }
    ev.assign(N, 0.0);
    for (std::size_t i = 0; i < N; ++i) ev[i] = A[i * N + i];
}

// W_0 = Q diag(1 / max(|mu_i|, lambda)) Q^T of H_F = H_R + lambda I (PROTOCOL-v3 section 3).
Vec modified_inverse(const Vec& HR, int n, double lam, int& n_mod, Vec& mu) {
    const std::size_t N = U(n);
    Vec HF = HR;
    for (std::size_t i = 0; i < N; ++i) HF[i * N + i] += lam;
    Vec Q;
    jacobi_eig(HF, n, mu, Q);
    n_mod = 0;
    Vec d(N);
    for (std::size_t i = 0; i < N; ++i) {
        if (!(mu[i] >= lam)) ++n_mod;  // eigenvalue replaced by max(|mu|, lambda)
        d[i] = 1.0 / std::fmax(std::fabs(mu[i]), lam);
    }
    Vec W(N * N, 0.0);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = i; j < N; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < N; ++k) s += Q[i * N + k] * d[k] * Q[j * N + k];
            W[i * N + j] = W[j * N + i] = s;
        }
    return W;
}

void write_bin(const std::string& path, const Vec& v) {
    std::ofstream f(path, std::ios::binary);
    f.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(double)));
}
Vec read_bin(const std::string& path, std::size_t count) {
    Vec v(count, 0.0);
    std::ifstream f(path, std::ios::binary);
    f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(count * sizeof(double)));
    if (!f) { std::fprintf(stderr, "read_bin: short read %s\n", path.c_str()); std::exit(3); }
    return v;
}

// ==== verbatim s3.cpp lines 216-324 ====
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
// B refinement (PROTOCOL 3.3 of study 3; PROTOCOL-v3 4.3, v2 R1): study-2 rule, every leaf whose closed rectangle
// contains the point is tested. Split decisions (need test, depth-6 test, node-cap arithmetic, edge-leaf testing,
// 16-pass guard) are verbatim from s3.cpp. Declared changes: templated on the rho evaluator rho(x, y); returns an
// outcome record counted over all test rows of the check; no persistent `capped` flag.
struct RefineOutcome {
    int n_need = 0, n_split = 0, n_refused_depth = 0, n_refused_cap = 0, n_outside_Dw = 0, n_rows = 0;
};

template <typename RhoFn>
RefineOutcome refine_B(const Ctx& c, const RhoFn& rho_at, NodeSet& ns, const std::vector<Pt>& cand,
                       std::ostream& ev, int it) {
    const double R0 = c.sys.parent.scale_R0();
    RefineOutcome o;
    for (const Pt& p : cand) {
        if (!(p.v < 0.0)) continue;
        const double xi = p.x / R0, s = (p.r / R0) * (p.r / R0);
        for (int guard = 0; guard < 16; ++guard) {
            std::vector<int> Ls;
            for (std::size_t e = 0; e < ns.elems.size(); ++e) {
                const Elem& E = ns.elems[e];
                if (xi >= E.xa && xi <= E.xb && s >= E.sa && s <= E.sb) Ls.push_back(static_cast<int>(e));
            }
            if (Ls.empty()) {  // outside Dw
                if (guard == 0) ++o.n_outside_Dw;
                break;
            }
            std::vector<int> todo;
            for (int L : Ls) {
                double rhoL = 1e300;
                for (std::size_t q = 0; q < ns.x.size(); ++q)
                    if (ns.elem_of[q] == L) rhoL = std::fmin(rhoL, rho_at(ns.x[q], ns.y[q]));
                const Elem& E = ns.elems[U(L)];
                const bool need = p.v < rhoL - 0.05 * std::fabs(p.v);
                ++o.n_rows;
                ev << it << "," << p.x << "," << p.r << "," << p.v << "," << rhoL << "," << E.depth << "," << L << ","
                   << Ls.size() << "," << (need ? 1 : 0);
                if (!need) { ev << ",resolved\n"; continue; }
                ++o.n_need;
                if (E.depth >= 6) { ev << ",max_depth\n"; ++o.n_refused_depth; continue; }
                if (ns.x.size() + 3u * (todo.size() + 1) * static_cast<std::size_t>(ns.order * ns.order) >
                    4u * ns.base_nodes) {
                    ev << ",node_cap\n"; ++o.n_refused_cap; continue;
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
                ++o.n_split;
            }
            rebuild_nodes(c, ns);
        }
    }
    ev.flush();
    return o;
}

// Outcome classification and the outer-loop decision (PROTOCOL-v3 4.3; v2 R1, v3 N5a). Tested by gate V5d.
enum class Chk { PASS, SPLIT, BLOCKED_DEPTH, BLOCKED_CAP };
const char* chk_name(Chk k) {
    switch (k) {
        case Chk::PASS: return "PASS";
        case Chk::SPLIT: return "SPLIT";
        case Chk::BLOCKED_DEPTH: return "BLOCKED_DEPTH";
        default: return "BLOCKED_CAP";
    }
}
Chk classify(const RefineOutcome& o) {
    if (o.n_need == 0) return Chk::PASS;
    if (o.n_split >= 1) return Chk::SPLIT;
    if (o.n_refused_cap >= 1) return Chk::BLOCKED_CAP;
    return Chk::BLOCKED_DEPTH;
}
struct OuterDecision {
    std::string term;  // empty: continue with the next inner solve
    bool apply_split = false;
    double trust = 0.0;
};
constexpr int kOuterCap = 10;
constexpr double kTrust0 = 0.5;
OuterDecision outer_decide(Chk k, int round, double trust, bool ls_limited) {
    OuterDecision d;
    d.trust = trust;
    switch (k) {
        case Chk::PASS: d.term = ls_limited ? "CONVERGED_LS_LIMITED" : "CONVERGED"; break;
        case Chk::BLOCKED_DEPTH: d.term = "REFINE_BLOCKED_DEPTH"; break;
        case Chk::BLOCKED_CAP: d.term = "REFINE_BLOCKED_CAP"; break;
        case Chk::SPLIT:
            if (round >= kOuterCap) { d.term = "OUTER_LIMIT"; break; }
            d.apply_split = true;
            d.trust = kTrust0;  // v3 N5a: trust reset after every SPLIT
            break;
    }
    return d;
}

std::vector<Elem> replay_splits(std::vector<Elem> elems, const std::string& events_csv, int& nsplit_rows) {
    std::ifstream f(events_csv);
    std::string line;
    std::getline(f, line);  // header
    nsplit_rows = 0;
    while (std::getline(f, line)) {
        std::vector<std::string> t;
        std::stringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ',')) t.push_back(tok);
        if (t.size() < 10 || t[9] != "split") continue;
        const int L = std::stoi(t[6]);
        const Elem E = elems[U(L)];
        const double xm = 0.5 * (E.xa + E.xb), sm = 0.5 * (E.sa + E.sb);
        const int d = E.depth + 1;
        elems[U(L)] = {E.xa, xm, E.sa, sm, d};
        elems.push_back({xm, E.xb, E.sa, sm, d});
        elems.push_back({E.xa, xm, sm, E.sb, d});
        elems.push_back({xm, E.xb, sm, E.sb, d});
        ++nsplit_rows;
    }
    return elems;
}

void write_elems(const std::string& path, const NodeSet& ns) {
    std::ofstream f(path);
    f.precision(17);
    f << "# order " << ns.order << " base_nodes " << ns.base_nodes << " leaves " << ns.elems.size() << "\n";
    for (const Elem& E : ns.elems) f << E.xa << " " << E.xb << " " << E.sa << " " << E.sb << " " << E.depth << "\n";
}
std::vector<Elem> read_elems(const std::string& path) {
    std::vector<Elem> v;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        Elem E{0, 0, 0, 0, 0};
        ss >> E.xa >> E.xb >> E.sa >> E.sb >> E.depth;
        v.push_back(E);
    }
    return v;
}

// Working node set for a mode: base set (order:sub), optionally replaced by a replayed study-3 B set (--recon DIR)
// or an element file (--set FILE). base_nodes always refers to the base set (node cap = 4x base).
NodeSet make_set(const Ctx& c, int order, int sub, const std::string& recon, const std::string& setfile) {
    NodeSet ns = base_set(c, order, sub, 1.7, 2.6);
    if (!recon.empty()) {
        int nr = 0;
        ns.elems = replay_splits(ns.elems, recon + "/refine_events.csv", nr);
        rebuild_nodes(c, ns);
    } else if (!setfile.empty()) {
        ns.elems = read_elems(setfile);
        rebuild_nodes(c, ns);
    }
    return ns;
}

// ==== verbatim s3.cpp lines 377-437 ====
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
// Solver pieces (PROTOCOL-v3 section 3), shared by the run mode and gate V5a.

// Line search verbatim from the s3 run loop (Armijo c1 = 1e-4, x0.5 backtracking, <= 10 trials, validity probe).
struct LsResult {
    bool ok = false;
    Vec yn;
    double fn = 0.0, Rn = 0.0, step = 0.0;
    int att = 0;
};
template <typename FvalFn, typename ProbeFn>
LsResult line_search(const FvalFn& Fval, const ProbeFn& probe_ok, const Vec& y, const Vec& p, double f, double dd,
                     double trust, int it_next, std::ostream& ls, int& nfev) {
    const int n = static_cast<int>(y.size());
    const double kC1 = 1e-4;
    const double pn = nrm(p);
    double alpha = pn > 0.0 ? std::fmin(1.0, trust / pn) : 0.0;
    LsResult r;
    for (; r.att < 10; ++r.att) {
        r.yn = y;
        for (int i = 0; i < n; ++i) r.yn[U(i)] += alpha * p[U(i)];
        r.step = 0.0;
        for (int i = 0; i < n; ++i) { const double d = r.yn[U(i)] - y[U(i)]; r.step += d * d; }
        r.step = std::sqrt(r.step);
        bool finite = true;
        for (double q : r.yn) if (!std::isfinite(q)) finite = false;
        bool armijo = false;
        if (finite) {
            r.fn = Fval(r.yn, r.Rn); ++nfev;
            armijo = std::isfinite(r.fn) && probe_ok(r.yn) && r.fn <= f + kC1 * alpha * dd;
        }
        ls << it_next << "," << r.att << "," << alpha << "," << trust << "," << r.fn << "," << f << ","
           << f + kC1 * alpha * dd << "," << armijo << "\n";
        if (armijo) { r.ok = true; break; }
        alpha *= 0.5;
    }
    ls.flush();
    return r;
}

// Dense BFGS on the inverse Hessian approximation W (row-major n x n).
struct DenseBfgs {
    int n = 0;
    Vec W, W0;
    Vec matvec(const Vec& M, const Vec& x) const {
        Vec r(U(n), 0.0);
        for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int j = 0; j < n; ++j) s += M[U(i) * U(n) + U(j)] * x[U(j)];
            r[U(i)] = s;
        }
        return r;
    }
    // p = -W g; if g.p >= 0, reset W to W0 (reset = 1).
    Vec direction(const Vec& g, int& reset) {
        reset = 0;
        Vec p = matvec(W, g);
        for (double& v : p) v = -v;
        if (!(dot(p, g) < 0.0)) {
            W = W0;
            reset = 1;
            p = matvec(W, g);
            for (double& v : p) v = -v;
        }
        return p;
    }
    // Standard inverse update; skipped (returns false) if s.t <= 1e-14 |s||t|.
    bool update(const Vec& s, const Vec& t) {
        const double st = dot(s, t);
        if (!(st > 1e-14 * nrm(s) * nrm(t))) return false;
        const Vec Wt = matvec(W, t);
        const double tWt = dot(t, Wt);
        const double a = (st + tWt) / (st * st), b = 1.0 / st;
        for (int i = 0; i < n; ++i)
            for (int j = i; j < n; ++j) {
                const std::size_t ij = U(i) * U(n) + U(j);
                const double v = W[ij] + a * s[U(i)] * s[U(j)] - b * (Wt[U(i)] * s[U(j)] + s[U(i)] * Wt[U(j)]);
                W[ij] = v;
                W[U(j) * U(n) + U(i)] = v;
            }
        return true;
    }
};

// Refresh eligibility (PROTOCOL-v3 section 3, conditions 1-4). gn[k] = |grad F| at accepted iterate k on the current
// set, skip[k] = 1 if the BFGS update at accepted iterate k was skipped. Returns "" if not eligible.
struct RefreshState {
    int last_meas = 0, last_change = 0, n_refresh = 0;
    Vec gn;
    std::vector<int> skip;
};
std::string refresh_trigger(const RefreshState& st, int k) {
    if (k - st.last_meas < 20) return "";
    if (k - st.last_change < 20) return "";
    if (st.n_refresh >= 2) return "";
    const bool tg = st.gn[U(k)] > 0.5 * st.gn[U(k - 20)];
    int sk = 0;
    for (int j = k - 19; j <= k; ++j) sk += st.skip[U(j)];
    const bool ts = sk > 10;
    if (tg && ts) return "grad+skip";
    if (tg) return "grad";
    if (ts) return "skip";
    return "";
}
double skip_fraction20(const std::vector<int>& skip, int k) {
    int sk = 0, m = 0;
    for (int j = std::max(1, k - 19); j <= k; ++j) { sk += skip[U(j)]; ++m; }
    return m ? static_cast<double>(sk) / m : 0.0;
}

// ---------------------------------------------------------------------------
// Hessian measurement and archive (PROTOCOL-v3 2.1-2.3).
struct HessReport {
    HessExt H;
    GradExt G;
    ResExt Rx;
};
HessReport measure_hessian(const Ctx& c, const Vec& y, const NodeSet& ns, const Passes& P) {
    HessReport hr;
    hr.H = assemble_hessian(WarpEval{&c, &ns}, y, P, static_cast<long>(ns.x.size()));
    hr.G = gradient_ext(c, y, ns);
    hr.Rx = residual_ext(c, y, ns);
    return hr;
}
void write_hessian_dir(const std::string& dir, const HessReport& hr, const Vec& y, const NodeSet& ns, double lam,
                       const std::string& tag) {
    std::filesystem::create_directories(dir);
    write_bin(dir + "/H_R_fwd.bin", hr.H.Hf);
    write_bin(dir + "/H_R_rev.bin", hr.H.Hr);
    write_bin(dir + "/H_R_comp.bin", hr.H.Hc);
    write_bin(dir + "/H_R_ld.bin", hr.H.Hl);
    write_bin(dir + "/B_H.bin", hr.H.Bh);
    write_bin(dir + "/gR_fwd.bin", hr.G.fwd);
    write_bin(dir + "/gR_rev.bin", hr.G.rev);
    write_bin(dir + "/gR_comp.bin", hr.G.comp);
    write_bin(dir + "/gR_ld.bin", hr.G.ld);
    write_bin(dir + "/gR_bound.bin", hr.G.bound);
    write_bin(dir + "/gR_pass.bin", hr.H.gf);
    write_bin(dir + "/y.bin", y);
    write_elems(dir + "/elems.txt", ns);
    long gvg = 0;
    for (std::size_t i = 0; i < y.size(); ++i)
        if (std::memcmp(&hr.H.gf[i], &hr.G.fwd[i], sizeof(double)) != 0) ++gvg;
    std::ofstream m(dir + "/meta.txt");
    m.precision(17);
    m << "tag " << tag << "\nnodes " << ns.x.size() << "\nleaves " << ns.elems.size() << "\nbase_nodes "
      << ns.base_nodes << "\nlambda " << lam << "\nR_fwd " << hr.Rx.fwd << "\nR_rev " << hr.Rx.rev << "\nR_comp "
      << hr.Rx.comp << "\nR_ld " << hr.Rx.ld << "\nR_bound " << hr.Rx.bound << "\nH1a_diag_mismatch "
      << hr.H.diag_mismatch << "\nH1a_passgrad_mismatch " << hr.H.grad_mismatch << "\nH1a_grad_vs_gradient_mismatch "
      << gvg << "\npair_errors " << hr.H.pair_errors << "\nsecs " << hr.H.secs << "\nthreads " << omp_get_max_threads()
      << "\nld_digits " << std::numeric_limits<long double>::digits << "\n";
}

// ---------------------------------------------------------------------------
struct RunArgs {
    char form = 'A';
    int order = 4, sub = 1, kmax = 200, diag_every = 1;
    std::string y0, recon, setfile, hdir, out, passes;
    double lam = 0.0, deadline = 0.0;
};

int mode_run(const Ctx& c, const RunArgs& A) {
    const double t_start = now();
    const int n = 3 * c.block;
    Vec y = read_vec(A.y0);
    if (static_cast<int>(y.size()) != n) { std::fprintf(stderr, "bad y0 %s\n", A.y0.c_str()); return 2; }
    NodeSet ns = make_set(c, A.order, A.sub, A.recon, A.setfile);
    const Val V = make_val(c);
    const Passes P = read_passes(A.passes);
    const double xw = 1.7 * 2.0, rw = 2.0 * std::sqrt(2.6);
    const double lam = A.lam;
    const char form = A.form;
    const std::string& out = A.out;
    const int kmax = A.kmax, maxfev = 3 * A.kmax;
    // Measured Hessian at the start design (step 0); must match y0 and the node set.
    {
        const Vec yh = read_bin(A.hdir + "/y.bin", U(n));
        if (std::memcmp(yh.data(), y.data(), U(n) * sizeof(double)) != 0) {
            std::fprintf(stderr, "Hessian design %s does not match y0\n", A.hdir.c_str());
            return 4;
        }
    }
    DenseBfgs bf;
    bf.n = n;
    int nmod0 = 0;
    Vec mu0;
    bf.W0 = modified_inverse(read_bin(A.hdir + "/H_R_fwd.bin", U(n) * U(n)), n, lam, nmod0, mu0);
    bf.W = bf.W0;
    {
        std::ofstream cf(out + "/config.txt");
        cf.precision(17);
        cf << "form " << form << "\norder " << A.order << "\nsub " << A.sub << "\ny0 " << A.y0 << "\nrecon " << A.recon
           << "\nset " << A.setfile << "\nhessian " << A.hdir << "\nlambda " << lam << "\nkmax " << kmax
           << "\nmaxfev " << maxfev << "\ndeadline_epoch " << A.deadline << "\ndiag_every " << A.diag_every
           << "\nbase_nodes " << ns.base_nodes << "\nnodes " << ns.x.size() << "\nleaves " << ns.elems.size()
           << "\ntau " << c.tau << "\nthreads " << omp_get_max_threads() << "\nW0_modified_eigs " << nmod0
           << "\nW0_mu_min " << *std::min_element(mu0.begin(), mu0.end()) << "\nW0_mu_max "
           << *std::max_element(mu0.begin(), mu0.end())
           << "\nF = R + lambda/2 |y|^2\nC1' |gradF| <= 1e-5 F / max(|y|,1)\nC2' 10-window 1e-8 since last set change\n"
              "STALLED label: C2' & !C1' for 30 consecutive iterates on one set (not terminal)\n"
              "dense BFGS, W0 = Q diag(1/max(|mu|,lambda)) Q^T; c1=1e-4 trust0=0.5 trustmax=4 trustmin=1e-8 "
              "ls_trials=10; skip if s.t <= 1e-14|s||t|; refresh rule PROTOCOL-v3 s3; B: refine only at inner "
              "convergence, outer cap 10, trust reset 0.5 after SPLIT\n";
    }
    std::ofstream tr(out + "/trace.csv");
    tr.precision(12);
    tr << "it,F,R,pen,ynorm,gF,gR,step,trust,gamma,reset,nfev,ngev,nodes,leaves,splits,rhomin_w_pre,rhomin_w_post,"
          "diag,Eneg_v8,Eneg_v8_Dw,negfrac_v8,vtail_v8,L4_v8,rhomin_v8,refined_v8,rhostar_in,xstar,rstar,n_distinct,"
          "n_polished,n_dup,n_latmin,C1,C2,C3,t_grad,t_search,t_diag,secs,"
          "pnorm,skip,wreset,skipfrac20,refresh,refresh_trigger,e_g,e_g_c,eg_hat,tolC1,prec_limited,stalled,check,"
          "outer_round,inner_it\n";
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
    std::ofstream oc(out + "/outer.csv");
    oc.precision(12);
    oc << "check,it,outcome,n_need,n_split,n_refused_depth,n_refused_cap,n_outside_Dw,n_rows,trust_before,"
          "trust_after,ls_limited,nodes_after\n";
    std::ofstream rf(out + "/refresh.csv");
    rf.precision(12);
    rf << "it,trigger,gn_k,gn_k20,skips20,secs,modified_eigs\n";

    int nfev = 0, ngev = 0, wresets = 0, nskips = 0;
    double t_grad = 0, t_search = 0, t_diag = 0;
    auto Fval = [&](const Vec& yy, double& Rout) {
        Rout = residual(c, yy, ns);
        return Rout + 0.5 * lam * dot(yy, yy);
    };
    auto probe_ok = [&](const Vec& yy) {
        return adm_metric_jet<double>(c.f2a, c.cfg, c.sys.parent, fields_of<double>(c, yy), 0.0, 0.7, 1.3, 0.8).valid;
    };
    struct GInfo { Vec g, gR; double e_g, e_g_c; };
    auto Fgrad = [&](const Vec& yy) {
        const double t0 = now();
        const GradExt G = gradient_ext(c, yy, ns);
        GInfo gi;
        gi.gR = G.fwd;
        gi.g = G.fwd;
        Vec gr = G.rev, gc = G.comp;
        for (std::size_t i = 0; i < gi.g.size(); ++i) {
            gi.g[i] += lam * yy[i];
            gr[i] += lam * yy[i];
            gc[i] += lam * yy[i];
        }
        double a = 0, b = 0;
        for (std::size_t i = 0; i < gi.g.size(); ++i) {
            a += (gi.g[i] - gr[i]) * (gi.g[i] - gr[i]);
            b += (gi.g[i] - gc[i]) * (gi.g[i] - gc[i]);
        }
        gi.e_g = std::sqrt(a);
        gi.e_g_c = std::sqrt(b);
        t_grad = now() - t0;
        return gi;
    };
    auto search = [&]() {
        const double t0 = now();
        InLoop r = inloop_search(c, fields_of<double>(c, y), xw, rw);
        t_search = now() - t0;
        return r;
    };
    auto wmin = [&]() {
        const Vec rr = node_rho(c, fields_of<double>(c, y), ns.x, ns.y);
        double m = 1e300;
        for (double v : rr) m = std::fmin(m, v);
        return m;
    };
    struct Row {
        int it = 0;
        double F = 0, R = 0, step = 0, trust = 0, pnorm = 0, e_g = 0, e_g_c = 0, rmw_pre = 0, skipfrac = 0;
        int wreset = 0, skip = 0, splits = 0, refresh = 0, stalled = 0, outer_round = 0, inner_it = 0;
        bool C1 = false, C2 = false;
        int check = 0;  // 0: no check at this iterate; 1: check made
        int C3 = -1;    // -1: empty (no check)
        std::string trig;
    };
    auto log_iter = [&](const Row& rw_, const Vec& gF, const Vec& gR, const InLoop& in, bool force_diag) {
        const double rmw_post = rw_.splits > 0 ? wmin() : rw_.rmw_pre;
        const bool dg = force_diag || (rw_.it % A.diag_every == 0);
        Diag d{NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN};
        t_diag = 0;
        if (dg) { const double t0 = now(); d = validation_diag(c, y, V.vs, V.eps_rho); t_diag = now() - t0; }
        const double tol = 1e-5 * rw_.F / std::fmax(nrm(y), 1.0);
        const double eg_hat = std::fmax(rw_.e_g, rw_.e_g_c);
        tr << rw_.it << "," << rw_.F << "," << rw_.R << "," << 0.5 * lam * dot(y, y) << "," << nrm(y) << ","
           << nrm(gF) << "," << nrm(gR) << "," << rw_.step << "," << rw_.trust << "," << NAN << "," << rw_.wreset
           << "," << nfev << "," << ngev << "," << ns.x.size() << "," << ns.elems.size() << "," << rw_.splits << ","
           << rw_.rmw_pre << "," << rmw_post << "," << (dg ? 1 : 0) << "," << d.Eneg << "," << d.EnegW << ","
           << d.negfrac << "," << d.vtail << "," << d.L4 << "," << d.rhomin << "," << d.refined << ","
           << in.mins[0].v << "," << in.mins[0].x << "," << in.mins[0].r << "," << in.mins.size() << ","
           << in.polished << "," << in.duplicates << "," << in.lattice_minima << "," << rw_.C1 << "," << rw_.C2
           << ",";
        if (rw_.C3 >= 0) tr << rw_.C3;
        tr << "," << t_grad << "," << t_search << "," << t_diag << "," << (now() - t_start) << "," << rw_.pnorm << ","
           << rw_.skip << "," << rw_.wreset << "," << rw_.skipfrac << "," << rw_.refresh << "," << rw_.trig << ","
           << rw_.e_g << "," << rw_.e_g_c << "," << eg_hat << "," << tol << "," << (eg_hat > 0.1 * tol ? 1 : 0)
           << "," << rw_.stalled << "," << rw_.check << "," << rw_.outer_round << "," << rw_.inner_it << "\n";
        tr.flush();
        for (std::size_t k = 0; k < in.mins.size(); ++k)
            im << rw_.it << "," << k << "," << in.mins[k].x << "," << in.mins[k].r << "," << in.mins[k].v << "\n";
        im.flush();
        for (double v : y) ys << v << " ";
        ys << "\n";
        ys.flush();
    };

    double R = 0.0;
    double f = Fval(y, R); ++nfev;
    GInfo G = Fgrad(y); ++ngev;
    Vec g = G.g, gR = G.gR;
    InLoop cand = search();
    RefreshState rs;
    rs.gn.push_back(nrm(g));
    rs.skip.push_back(0);
    std::vector<double> fh;
    fh.push_back(f);
    {
        Row r0;
        r0.it = 0; r0.F = f; r0.R = R; r0.trust = kTrust0; r0.e_g = G.e_g; r0.e_g_c = G.e_g_c; r0.rmw_pre = wmin();
        log_iter(r0, g, gR, cand, true);
    }
    double trust = kTrust0;
    std::string term;
    int it = 0, stall_run = 0, n_stalled = 0, inner_start = 0;
    int n_check = 0, n_outer = 0, k_chk1 = -1, cum_need = 0, cum_split = 0, cum_rdepth = 0, cum_rcap = 0, cum_out = 0;
    RefineOutcome last_o;
    std::string last_outcome = "none";
    bool C1 = false, C2 = false;
    NodeSet ns_before_check;
    // Refinement check at inner convergence (B). Returns true if the run terminates.
    auto do_check = [&](bool ls_limited, Row* row) {
        ++n_check;
        if (k_chk1 < 0) k_chk1 = it;
        ns_before_check = ns;
        const AdmFieldCoefficients<double> fc = fields_of<double>(c, y);
        auto rho_at = [&](double x, double r) { return rho_pt(c, fc, x, r); };
        const RefineOutcome o = refine_B(c, rho_at, ns, cand.mins, ev, it);
        last_o = o;
        const Chk k = classify(o);
        last_outcome = chk_name(k);
        cum_need += o.n_need; cum_split += o.n_split; cum_rdepth += o.n_refused_depth; cum_rcap += o.n_refused_cap;
        cum_out += o.n_outside_Dw;
        const OuterDecision dec = outer_decide(k, n_check, trust, ls_limited);
        if (!dec.apply_split && k == Chk::SPLIT) ns = ns_before_check;  // OUTER_LIMIT: split not applied
        if (row) { row->check = 1; row->C3 = (o.n_need == 0) ? 1 : 0; row->splits = dec.apply_split ? o.n_split : 0; }
        oc << n_check << "," << it << "," << chk_name(k) << "," << o.n_need << "," << o.n_split << ","
           << o.n_refused_depth << "," << o.n_refused_cap << "," << o.n_outside_Dw << "," << o.n_rows << "," << trust
           << "," << dec.trust << "," << (ls_limited ? 1 : 0) << "," << ns.x.size() << "\n";
        oc.flush();
        if (!dec.term.empty()) { term = dec.term; return true; }
        // SPLIT applied: recompute F and grad on the new set; keep W; trust reset; restart C2 history.
        ++n_outer;
        trust = dec.trust;
        const double f_acc = f;
        f = Fval(y, R); ++nfev;
        G = Fgrad(y); ++ngev;
        g = G.g; gR = G.gR;
        fh.clear(); fh.push_back(f);
        stall_run = 0;
        rs.last_change = it;
        rs.gn[U(it)] = nrm(g);
        inner_start = it;
        C1 = false; C2 = false;
        std::ofstream(out + "/set_changes.csv", std::ios::app)
            << it << "," << o.n_split << "," << ns.x.size() << "," << f_acc << "," << f << "," << nrm(g) << "\n";
        return false;
    };

    while (true) {
        if (it >= kmax) { term = "ITERATION_LIMIT"; break; }
        if (nfev >= maxfev) { term = "EVAL_LIMIT"; break; }
        if (A.deadline > 0 && epoch() >= A.deadline) { term = "WALL_LIMIT"; break; }
        int wreset = 0;
        const Vec p = bf.direction(g, wreset);
        if (wreset) ++wresets;
        const double dd = dot(p, g);
        const LsResult L = line_search([&](const Vec& yy, double& Ro) { return Fval(yy, Ro); }, probe_ok, y, p, f, dd,
                                       trust, it + 1, ls, nfev);
        if (!L.ok) {
            trust *= 0.5;
            if (trust < 1e-8) {
                if (C1 && C2) {
                    if (form != 'B') { term = "CONVERGED_LS_LIMITED"; break; }
                    if (do_check(true, nullptr)) break;
                    continue;
                }
                term = "LINESEARCH_FAIL";
                break;
            }
            continue;
        }
        const GInfo G2 = Fgrad(L.yn); ++ngev;
        bool finite = std::isfinite(L.fn);
        for (double q : G2.g) if (!std::isfinite(q)) finite = false;
        if (!finite) { term = "INVALID"; break; }
        Vec sv(U(n)), tv(U(n));
        for (int i = 0; i < n; ++i) { sv[U(i)] = L.yn[U(i)] - y[U(i)]; tv[U(i)] = G2.g[U(i)] - g[U(i)]; }
        const bool updated = bf.update(sv, tv);
        if (!updated) ++nskips;
        y = L.yn; G = G2; g = G2.g; gR = G2.gR; f = L.fn; R = L.Rn;
        ++it;
        if (L.att == 0) trust = std::fmin(trust * 2.0, 4.0);
        fh.push_back(f);
        rs.gn.push_back(nrm(g));
        rs.skip.push_back(updated ? 0 : 1);
        const std::size_t m = fh.size() - 1;
        C1 = nrm(g) <= 1e-5 * f / std::fmax(nrm(y), 1.0);
        C2 = m >= 10 && (fh[m - 10] - fh[m]) / fh[m] <= 1e-8;
        stall_run = (C2 && !C1) ? stall_run + 1 : 0;
        const int stalled = stall_run >= 30 ? 1 : 0;
        n_stalled += stalled;
        cand = search();
        Row row;
        row.it = it; row.F = f; row.R = R; row.step = L.step; row.trust = trust; row.pnorm = nrm(p);
        row.e_g = G.e_g; row.e_g_c = G.e_g_c; row.rmw_pre = wmin(); row.wreset = wreset; row.skip = updated ? 0 : 1;
        row.skipfrac = skip_fraction20(rs.skip, it); row.stalled = stalled; row.C1 = C1; row.C2 = C2;
        row.inner_it = it - inner_start;
        const Vec g_acc = g, gR_acc = gR;
        bool stop = false;
        if (C1 && C2) {
            if (form != 'B') { term = "CONVERGED"; stop = true; }
            else stop = do_check(false, &row);
        }
        row.outer_round = n_check;
        if (!stop) {
            row.trig = refresh_trigger(rs, it);
            row.refresh = row.trig.empty() ? 0 : 1;
        }
        // Row reports the accepted iterate (pre-split values), as in s3.
        {
            const Vec ysave = y;
            log_iter(row, g_acc, gR_acc, cand, stop || it >= kmax);
            (void)ysave;
        }
        std::printf("it %3d F %.10e |y| %.3f |gF| %.3e (tol %.2e) step %.3e nodes %zu skip %d wr %d C %d%d %s%.0fs\n",
                    it, row.F, nrm(y), nrm(g_acc), 1e-5 * row.F / std::fmax(nrm(y), 1.0), L.step, ns.x.size(),
                    row.skip, wreset, C1, C2, row.check ? (std::string("chk ") + last_outcome + " ").c_str() : "",
                    now() - t_start);
        std::fflush(stdout);
        if (stop) break;
        if (row.refresh) {
            const double t0 = now();
            const HessReport hr = measure_hessian(c, y, ns, P);
            write_hessian_dir(out + "/refresh_it" + std::to_string(it), hr, y, ns, lam, "refresh");
            int nm = 0;
            Vec mu;
            bf.W0 = modified_inverse(hr.H.Hf, n, lam, nm, mu);
            bf.W = bf.W0;
            int sk = 0;
            for (int j = it - 19; j <= it; ++j) sk += rs.skip[U(j)];
            rf << it << "," << row.trig << "," << rs.gn[U(it)] << "," << rs.gn[U(it - 20)] << "," << sk << ","
               << now() - t0 << "," << nm << "\n";
            rf.flush();
            rs.last_meas = it;
            ++rs.n_refresh;
        }
    }
    {
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
        write_elems(out + "/set_final.txt", ns);
        std::ofstream ft(out + "/termination.txt");
        ft.precision(10);
        ft << term << " it=" << it << " nfev=" << nfev << " ngev=" << ngev << " nodes=" << ns.x.size()
           << " leaves=" << ns.elems.size() << " n_outer=" << n_outer << " n_check=" << n_check << " k_chk1="
           << (k_chk1 >= 0 ? std::to_string(k_chk1) : std::string("")) << " last_outcome=" << last_outcome
           << " last_record=" << last_o.n_need << "/" << last_o.n_split << "/" << last_o.n_refused_depth << "/"
           << last_o.n_refused_cap << "/" << last_o.n_outside_Dw << " cum_need=" << cum_need << " cum_split="
           << cum_split << " cum_refused_depth=" << cum_rdepth << " cum_refused_cap=" << cum_rcap
           << " cum_outside_Dw=" << cum_out << " last_set_change=" << rs.last_change << " w_resets=" << wresets
           << " bfgs_skips=" << nskips << " refreshes=" << rs.n_refresh << " stalled_iterates=" << n_stalled
           << " C1=" << C1 << " C2=" << C2 << " gF=" << nrm(g) << " C1tol=" << 1e-5 * f / std::fmax(nrm(y), 1.0)
           << " F=" << f << " R=" << R << " ynorm=" << nrm(y) << " wall_s=" << (now() - t_start) << "\n";
    }
    std::printf("TERMINATION %s it=%d nfev=%d wall=%.0fs\n", term.c_str(), it, nfev, now() - t_start);
    return 0;
}

// ---------------------------------------------------------------------------
Vec y_from_args(const std::string& yfile, const std::string& ysfile, int line) {
    if (ysfile.empty()) return read_vec(yfile);
    std::ifstream f(ysfile);
    std::string s;
    for (int k = 1; k <= line; ++k) std::getline(f, s);
    std::stringstream ss(s);
    Vec y;
    double v;
    while (ss >> v) y.push_back(v);
    return y;
}

int mode_hessian(const Ctx& c, const Vec& y, const NodeSet& ns, const std::string& passes, double lam,
                 const std::string& out, const std::string& tag) {
    const Passes P = read_passes(passes);
    const HessReport hr = measure_hessian(c, y, ns, P);
    write_hessian_dir(out, hr, y, ns, lam, tag);
    std::printf("HESSIAN %s nodes %zu passes %zu secs %.1f H1a diag_mismatch %ld passgrad_mismatch %ld pair_errors %ld\n",
                tag.c_str(), ns.x.size(), P.size(), hr.H.secs, hr.H.diag_mismatch, hr.H.grad_mismatch,
                hr.H.pair_errors);
    return (hr.H.pair_errors == 0) ? 0 : 1;
}

// FD data for H2/H3: each line of DIRS is "label h d_0 ... d_{n-1}"; writes gradients and residuals at y +- h d.
int mode_fd(const Ctx& c, const Vec& y, const NodeSet& ns, const std::string& dirs, const std::string& out) {
    std::ifstream f(dirs);
    std::ofstream o(out);
    o.precision(17);
    const double R0 = residual(c, y, ns);
    o << "R0 " << R0 << "\n";
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string label;
        double h;
        ss >> label >> h;
        Vec d;
        double v;
        while (ss >> v) d.push_back(v);
        Vec yp = y, ym = y;
        for (std::size_t i = 0; i < y.size(); ++i) { yp[i] += h * d[i]; ym[i] -= h * d[i]; }
        const Vec gp = gradient(c, yp, ns), gm = gradient(c, ym, ns);
        const double Rp = residual(c, yp, ns), Rm = residual(c, ym, ns);
        o << label << " " << h << " " << Rp << " " << Rm;
        for (double q : gp) o << " " << q;
        for (double q : gm) o << " " << q;
        o << "\n";
        o.flush();
        std::printf("fd %s h %.3e done\n", label.c_str(), h);
        std::fflush(stdout);
    }
    return 0;
}

int mode_recon(const Ctx& c, const std::string& run, int order, int sub, const std::string& out) {
    NodeSet ns = base_set(c, order, sub, 1.7, 2.6);
    int nr = 0;
    ns.elems = replay_splits(ns.elems, run + "/refine_events.csv", nr);
    rebuild_nodes(c, ns);
    const Vec y = read_vec(run + "/y_final.txt");
    const Vec r = node_rho(c, fields_of<double>(c, y), ns.x, ns.y);
    double m = 1e300;
    for (double v : r) m = std::fmin(m, v);
    write_elems(out, ns);
    std::printf("RECON %s split_rows %d nodes %zu leaves %zu rhomin_w %.17g\n", run.c_str(), nr, ns.x.size(),
                ns.elems.size(), m);
    return 0;
}

int mode_v2dump(const Ctx& c, const Vec& y, const NodeSet& ns, double lam) {
    const double R = residual(c, y, ns);
    const double F = R + 0.5 * lam * dot(y, y);
    Vec g = gradient(c, y, ns);
    std::printf("nodes %zu\nR %a\nF %a\n", ns.x.size(), R, F);
    for (std::size_t i = 0; i < g.size(); ++i) {
        const double gF = g[i] + lam * y[i];
        std::printf("g %zu %a %a\n", i, g[i], gF);
    }
    return 0;
}

// ---------------------------------------------------------------------------
// Gates V5a-V5e (PROTOCOL-v3 section 9). Synthetic problems only; no warp objective evaluation except the node-set
// geometry for V5d.

struct SynthEval {  // f(y) = 1/2 y^T M y + sum_m softplus_1(w_m^T y)^2, as 108 row terms + 5 softplus terms
    const Vec* M;
    const std::vector<Vec>* w;
    int n;
    void operator()(const std::vector<OTIScalar>& yy, std::vector<OTIScalar>& contrib) const {
        const std::size_t nw = w->size();
        contrib.assign(U(n) + nw, OTIScalar());
#pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            OTIScalar my = scalar_cast<OTIScalar>(0.0);
            for (int j = 0; j < n; ++j) my = my + scalar_cast<OTIScalar>((*M)[U(i) * U(n) + U(j)]) * yy[U(j)];
            contrib[U(i)] = scalar_cast<OTIScalar>(0.5) * yy[U(i)] * my;
        }
        for (std::size_t m = 0; m < nw; ++m) {
            OTIScalar z = scalar_cast<OTIScalar>(0.0);
            for (int j = 0; j < n; ++j) z = z + scalar_cast<OTIScalar>((*w)[m][U(j)]) * yy[U(j)];
            const OTIScalar s = softplus(z, 1.0);
            contrib[U(n) + m] = s * s;
        }
    }
};

int selftest(const Ctx& c, const std::string& passes, const std::string& out) {
    std::filesystem::create_directories(out);
    int fails = 0;
    auto report = [&](const char* gate, const char* what, bool ok, const std::string& detail) {
        std::printf("%s %-62s %s  %s\n", gate, what, ok ? "PASS" : "FAIL", detail.c_str());
        std::fflush(stdout);
        if (!ok) ++fails;
    };
    auto fmt = [](const char* f, double a, double b = 0.0) {
        char buf[256];
        std::snprintf(buf, sizeof buf, f, a, b);
        return std::string(buf);
    };
    const int n = 3 * c.block;
    std::mt19937_64 rng(20261027);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::uniform_real_distribution<double> ud(-1.0, 1.0);

    // ---- V5b numeric (supplement to the analytic check): both softplus branches at u = 0
    {
        const double tau = c.tau;
        OTIScalar up = sensitivity::seed_parameter(0.0, 0);
        OTIScalar um = sensitivity::seed_parameter(-1e-300, 0);
        const OTIScalar sp = softplus(up, tau), sm = softplus(um, tau);
        const double d1p = sensitivity::extract_gradient(sp, 0), d1m = sensitivity::extract_gradient(sm, 0);
        const double d2p = sensitivity::extract_hessian(sp, 0, 0), d2m = sensitivity::extract_hessian(sm, 0, 0);
        const double j1 = std::fabs(d1p - d1m) / std::fabs(d1p), j2 = std::fabs(d2p - d2m) / std::fabs(d2p);
        report("V5b", "softplus branch jump f', f'' at u=0 (numeric OTI) <= 1e-12 rel", j1 <= 1e-12 && j2 <= 1e-12,
               fmt("f' %.3e f'' %.3e", j1, j2) + fmt(" (f'=%.17g, f''*4tau=%.17g)", d1p, d2p * 4 * tau));
    }

    // ---- V5e (i) Neumaier unit cases
    {
        const Vec a{1.0, 1e100, 1.0, -1e100};
        const double s = neumaier_sum(a);
        report("V5e", "Neumaier [1,1e100,1,-1e100] = 2", s == 2.0, fmt("got %.17g", s));
        Vec x;
        std::uniform_int_distribution<int> ue(-40, 40);
        for (int k = 0; k < 50000; ++k) {
            const double v = ud(rng) * std::ldexp(1.0, ue(rng));
            x.push_back(v);
            x.push_back(-v * (1.0 + 1e-9 * ud(rng)));
        }
        std::shuffle(x.begin(), x.end(), rng);
        const double ns_ = neumaier_sum(x);
        double sa = 0.0, sf = 0.0;
        for (double v : x) { sa += std::fabs(v); sf = sf + v; }
        std::ofstream f(out + "/v5e_terms.txt");
        f << "neumaier " << std::hexfloat << ns_ << "\nforward " << sf << "\nsumabs " << sa << "\n";
        for (double v : x) f << v << "\n";
        std::printf("V5e 1e5-term sequence written (verified exactly with fractions.Fraction by v5_check.py): "
                    "neumaier %.17g forward %.17g sum|x| %.6g\n", ns_, sf, sa);
    }
    std::printf("V5e long double digits %d\n", std::numeric_limits<long double>::digits);

    // ---- V5c assembly on a known dense Hessian (+ V5e ii)
    {
        Vec M(U(n) * U(n));
        for (int i = 0; i < n; ++i)
            for (int j = i; j < n; ++j) {
                const double v = ud(rng) + 1e-3 * (i * n + j) / (n * n);
                M[U(i) * U(n) + U(j)] = M[U(j) * U(n) + U(i)] = v;
            }
        Vec ent;
        for (int i = 0; i < n; ++i) for (int j = i; j < n; ++j) ent.push_back(M[U(i) * U(n) + U(j)]);
        std::sort(ent.begin(), ent.end());
        bool distinct = std::adjacent_find(ent.begin(), ent.end()) == ent.end();
        report("V5c", "M has 5886 pairwise-distinct independent entries", distinct && ent.size() == 5886,
               fmt("%.0f entries", static_cast<double>(ent.size())));
        std::vector<Vec> w(5, Vec(U(n)));
        for (Vec& wm : w) for (double& v : wm) v = 2.0 * nd(rng) / std::sqrt(static_cast<double>(n));
        Vec y(U(n));
        for (double& v : y) v = nd(rng);
        // analytic (long double)
        Vec Hx(U(n) * U(n)), gx(U(n));
        for (int i = 0; i < n; ++i) {
            long double gi = 0;
            for (int j = 0; j < n; ++j) {
                Hx[U(i) * U(n) + U(j)] = M[U(i) * U(n) + U(j)];
                gi += static_cast<long double>(M[U(i) * U(n) + U(j)]) * y[U(j)];
            }
            gx[U(i)] = static_cast<double>(gi);
        }
        for (const Vec& wm : w) {
            long double z = 0;
            for (int j = 0; j < n; ++j) z += static_cast<long double>(wm[U(j)]) * y[U(j)];
            const long double sp = std::log1p(std::exp(z)), sg = 1.0L / (1.0L + std::exp(-z));
            const long double gco = 2.0L * sp * sg, hco = 2.0L * (sg * sg + sp * sg * (1.0L - sg));
            for (int i = 0; i < n; ++i) {
                gx[U(i)] = static_cast<double>(gx[U(i)] + gco * wm[U(i)]);
                for (int j = 0; j < n; ++j)
                    Hx[U(i) * U(n) + U(j)] = static_cast<double>(Hx[U(i) * U(n) + U(j)] + hco * wm[U(i)] * wm[U(j)]);
            }
        }
        const Passes P = read_passes(passes);
        const SynthEval se{&M, &w, n};
        const long nodes = n + 5;
        const HessExt H = assemble_hessian(se, y, P, nodes);
        const GradExt G = gradient_ext_gen(se, y, nodes);
        double hmax = 0, gmax = 0;
        for (double v : Hx) hmax = std::fmax(hmax, std::fabs(v));
        for (double v : gx) gmax = std::fmax(gmax, std::fabs(v));
        auto maxerr = [&](const Vec& Hv) {
            double e = 0;
            for (int i = 0; i < n; ++i)
                for (int j = i; j < n; ++j)
                    e = std::fmax(e, std::fabs(Hv[U(i) * U(n) + U(j)] - Hx[U(i) * U(n) + U(j)]));
            return e;
        };
        auto gerr = [&](const Vec& gv) {
            double e = 0;
            for (int i = 0; i < n; ++i) e = std::fmax(e, std::fabs(gv[U(i)] - gx[U(i)]));
            return e;
        };
        report("V5c", "pass schedule covers every pair once (pair_errors = 0)", H.pair_errors == 0,
               fmt("pair_errors %.0f, passes %.0f", static_cast<double>(H.pair_errors), static_cast<double>(P.size())));
        report("V5c", "H1a on synthetic: 54 diag / grad values bit-identical", H.diag_mismatch == 0 && H.grad_mismatch == 0,
               fmt("diag %.0f grad %.0f", static_cast<double>(H.diag_mismatch), static_cast<double>(H.grad_mismatch)));
        const double eF = maxerr(H.Hf), eC = maxerr(H.Hc), eR = maxerr(H.Hr), eL = maxerr(H.Hl);
        report("V5c", "every entry of H (fwd) within 1e-12 max|H| of analytic", eF <= 1e-12 * hmax,
               fmt("max err %.3e (max|H| %.3e)", eF, hmax));
        report("V5c", "gradient (pass and gradient_ext fwd) within 1e-12 max|g|",
               gerr(H.gf) <= 1e-12 * gmax && gerr(G.fwd) <= 1e-12 * gmax,
               fmt("pass %.3e gradient_ext %.3e", gerr(H.gf), gerr(G.fwd)));
        report("V5e", "(ii) H_comp, grad_comp within V5c tolerance", eC <= 1e-12 * hmax && gerr(G.comp) <= 1e-12 * gmax,
               fmt("H %.3e g %.3e", eC, gerr(G.comp)));
        std::printf("V5c info: rev max err %.3e, ld max err %.3e\n", eR, eL);
        long gvg = 0;
        for (int i = 0; i < n; ++i) if (std::memcmp(&H.gf[U(i)], &G.fwd[U(i)], sizeof(double)) != 0) ++gvg;
        report("V5c", "pass gradient bit-identical to gradient_ext (H1a path)", gvg == 0,
               fmt("%.0f mismatches", static_cast<double>(gvg)));
        // archive for v5_check.py (spectral-norm checks e_H^c <= |B_H|_F)
        const std::string d = out + "/v5c";
        std::filesystem::create_directories(d);
        write_bin(d + "/H_R_fwd.bin", H.Hf); write_bin(d + "/H_R_rev.bin", H.Hr); write_bin(d + "/H_R_comp.bin", H.Hc);
        write_bin(d + "/H_R_ld.bin", H.Hl); write_bin(d + "/B_H.bin", H.Bh); write_bin(d + "/H_exact.bin", Hx);
        write_bin(d + "/g_exact.bin", gx); write_bin(d + "/gR_fwd.bin", G.fwd); write_bin(d + "/gR_comp.bin", G.comp);
        // R value: forward vs analytic (informative)
        std::ofstream(d + "/meta.txt") << "nodes " << nodes << "\n";
    }

    // ---- V5a BFGS unit checks
    {
        // SPD quadratic f = 1/2 (y-ys)^T A (y-ys), A = Q diag(L) Q^T, L in [1, 1e3]
        Vec S0(U(n) * U(n));
        for (int i = 0; i < n; ++i)
            for (int j = i; j < n; ++j) S0[U(i) * U(n) + U(j)] = S0[U(j) * U(n) + U(i)] = nd(rng);
        Vec ev0, Q;
        jacobi_eig(S0, n, ev0, Q);
        Vec Lm(U(n));
        for (int i = 0; i < n; ++i) Lm[U(i)] = std::pow(10.0, 3.0 * (i + 0.5) / n);
        Vec Am(U(n) * U(n), 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                double s = 0;
                for (int k = 0; k < n; ++k) s += Q[U(i) * U(n) + U(k)] * Lm[U(k)] * Q[U(j) * U(n) + U(k)];
                Am[U(i) * U(n) + U(j)] = s;
            }
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) Am[U(j) * U(n) + U(i)] = Am[U(i) * U(n) + U(j)];
        Vec ystar(U(n)), y(U(n)), dlt(U(n));
        for (double& v : ystar) v = nd(rng);
        for (double& v : dlt) v = nd(rng);
        const double dn = nrm(dlt);
        for (int i = 0; i < n; ++i) y[U(i)] = ystar[U(i)] + 0.1 * dlt[U(i)] / dn;
        DenseBfgs bq;
        bq.n = n;
        auto qf = [&](const Vec& yy, double& Ro) {
            Vec d(U(n));
            for (int i = 0; i < n; ++i) d[U(i)] = yy[U(i)] - ystar[U(i)];
            const Vec Ad = bq.matvec(Am, d);
            Ro = 0.5 * dot(d, Ad);
            return Ro;
        };
        auto qg = [&](const Vec& yy) {
            Vec d(U(n));
            for (int i = 0; i < n; ++i) d[U(i)] = yy[U(i)] - ystar[U(i)];
            return bq.matvec(Am, d);
        };
        int nmod = 0;
        Vec mu;
        bq.W0 = modified_inverse(Am, n, 1e-12, nmod, mu);
        bq.W = bq.W0;
        double Rq = 0;
        const double f0 = qf(y, Rq);
        const Vec g0 = qg(y);
        int wr = 0, nfev = 0;
        const Vec p = bq.direction(g0, wr);
        std::ostringstream lsink;
        const LsResult L = line_search(qf, [](const Vec&) { return true; }, y, p, f0, dot(p, g0), kTrust0, 1, lsink,
                                       nfev);
        const Vec g1 = qg(L.yn);
        report("V5a", "SPD quadratic, W0 = H^-1 (Jacobi path): one step converges", L.ok && L.att == 0 && nrm(g1) <= 1e-10 * nrm(g0),
               fmt("|g1|/|g0| %.3e, trials %.0f", nrm(g1) / nrm(g0), static_cast<double>(L.att + 1)));
        // secant equation
        DenseBfgs bs;
        bs.n = n;
        bs.W = bq.W0;
        bs.W0 = bq.W0;
        Vec s(U(n)), t(U(n));
        for (double& v : s) v = nd(rng);
        t = bq.matvec(Am, s);
        for (int i = 0; i < n; ++i) t[U(i)] += 0.3 * nd(rng);
        if (dot(s, t) <= 0) for (double& v : t) v = -v;
        const bool upd = bs.update(s, t);
        const Vec Wt = bs.matvec(bs.W, t);
        double se = 0;
        for (int i = 0; i < n; ++i) se += (Wt[U(i)] - s[U(i)]) * (Wt[U(i)] - s[U(i)]);
        se = std::sqrt(se) / nrm(s);
        report("V5a", "inverse update satisfies secant W+ t = s to 1e-12", upd && se <= 1e-12, fmt("rel err %.3e", se));
        // skip
        Vec t2 = s;
        for (double& v : t2) v = -v;  // s.t < 0
        const Vec Wb = bs.W;
        const bool upd2 = bs.update(s, t2);
        report("V5a", "skip rule fires (s.t <= 1e-14|s||t|), W unchanged", !upd2 && bs.W == Wb, "");
        Vec t3(U(n), 0.0);
        t3[0] = 1.0;
        Vec s3v(U(n), 0.0);
        s3v[1] = 1.0;  // orthogonal: s.t = 0
        report("V5a", "skip rule fires for s.t = 0", !bs.update(s3v, t3), "");
        // reset
        DenseBfgs br;
        br.n = n;
        br.W0 = bq.W0;
        br.W.assign(U(n) * U(n), 0.0);
        for (int i = 0; i < n; ++i) br.W[U(i) * U(n) + U(i)] = -1.0;  // negative definite: g.p > 0
        int wr2 = 0;
        const Vec p2 = br.direction(g0, wr2);
        report("V5a", "reset rule fires on g.p >= 0 and restores W0", wr2 == 1 && br.W == br.W0 && dot(p2, g0) < 0, "");
        // skip fraction and refresh eligibility
        RefreshState st;
        st.gn.assign(41, 1.0);
        st.skip.assign(41, 0);
        st.last_meas = 0; st.last_change = 0; st.n_refresh = 0;
        for (int k = 30; k <= 40; ++k) st.skip[U(k)] = 1;  // 11 of iterates 21..40 skipped
        st.gn[40] = 0.1; st.gn[20] = 1.0;                  // gradient fell: grad trigger off
        report("V5a", "skip-fraction trigger fires (11/20 skipped)", refresh_trigger(st, 40) == "skip",
               "trigger=" + refresh_trigger(st, 40) + fmt(" skipfrac %.2f", skip_fraction20(st.skip, 40)));
        st.skip[30] = 0;  // 10/20: not > 50 %
        report("V5a", "skip-fraction trigger silent at 10/20", refresh_trigger(st, 40).empty(), "");
        st.gn[40] = 0.6;  // > 0.5 |g_{k-20}|
        report("V5a", "gradient-ratio trigger fires (|g_k| > 0.5 |g_k-20|)", refresh_trigger(st, 40) == "grad", "");
        st.last_meas = 25;
        report("V5a", "not eligible < 20 iterates since last measurement", refresh_trigger(st, 40).empty(), "");
        st.last_meas = 0; st.last_change = 21;
        report("V5a", "not eligible < 20 iterates since last set change (B)", refresh_trigger(st, 40).empty(), "");
        st.last_change = 20;
        report("V5a", "eligible exactly 20 iterates after a set change", refresh_trigger(st, 40) == "grad", "");
        st.n_refresh = 2;
        report("V5a", "not eligible after 2 refreshes", refresh_trigger(st, 40).empty(), "");
    }

    // ---- V5d refinement outcome logic (synthetic rho field on the 4:1 base set)
    {
        const NodeSet base = base_set(c, 4, 1, 1.7, 2.6);
        const std::size_t L0 = base.elems.size() / 2, L1 = base.elems.size() / 4;
        NodeSet tmp = base;
        {
            const Elem E = tmp.elems[L0];
            const double xm = 0.5 * (E.xa + E.xb), sm = 0.5 * (E.sa + E.sb);
            tmp.elems[L0] = {E.xa, xm, E.sa, sm, E.depth + 1};
            tmp.elems.push_back({xm, E.xb, E.sa, sm, E.depth + 1});
            tmp.elems.push_back({E.xa, xm, sm, E.sb, E.depth + 1});
            tmp.elems.push_back({xm, E.xb, sm, E.sb, E.depth + 1});
            rebuild_nodes(c, tmp);
        }
        std::size_t q0 = 0;
        while (tmp.elem_of[q0] != static_cast<int>(L0)) ++q0;
        const Pt pc{tmp.x[q0], tmp.y[q0], -2.0};  // a node of the first child of L0
        const double R0 = c.sys.parent.scale_R0();
        const Elem E1 = base.elems[L1];
        const double xi1 = 0.5 * (E1.xa + E1.xb), s1 = 0.5 * (E1.sa + E1.sb);
        const Pt pd{xi1 * R0, R0 * std::sqrt(s1), -2.0};  // centre of L1
        auto rho_syn = [&](double x, double r) {
            const double d2 = (x - pc.x) * (x - pc.x) + (r - pc.r) * (r - pc.r);
            return -1.0 - std::exp(-d2 / 1e-18);
        };
        std::ostringstream evs;
        {
            NodeSet ns = base;
            const RefineOutcome o1 = refine_B(c, rho_syn, ns, {pc}, evs, 1);
            const Chk k1 = classify(o1);
            report("V5d", "needed split below both limits -> SPLIT", k1 == Chk::SPLIT && o1.n_split == 1 &&
                   o1.n_refused_depth == 0 && o1.n_refused_cap == 0, std::string("got ") + chk_name(k1));
            const RefineOutcome o2 = refine_B(c, rho_syn, ns, {pc}, evs, 2);
            report("V5d", "resolved candidate -> PASS (C3 true)", classify(o2) == Chk::PASS && o2.n_need == 0,
                   std::string("got ") + chk_name(classify(o2)));
        }
        {
            NodeSet ns = base;
            ns.elems[L1].depth = 6;
            const RefineOutcome o = refine_B(c, rho_syn, ns, {pd}, evs, 1);
            report("V5d", "split needed at depth 6 -> BLOCKED_DEPTH, C3 false",
                   classify(o) == Chk::BLOCKED_DEPTH && o.n_need >= 1 && o.n_split == 0,
                   std::string("got ") + chk_name(classify(o)));
        }
        {
            NodeSet ns = base;
            ns.base_nodes = ns.x.size() / 4;  // set is at the node cap
            const RefineOutcome o = refine_B(c, rho_syn, ns, {pc}, evs, 1);
            report("V5d", "set at node cap -> BLOCKED_CAP, C3 false",
                   classify(o) == Chk::BLOCKED_CAP && o.n_need >= 1 && o.n_split == 0,
                   std::string("got ") + chk_name(classify(o)));
            const RefineOutcome o2 = refine_B(c, rho_syn, ns, {pc}, evs, 2);
            report("V5d", "check after a cap refusal is still executed (no persistent flag)",
                   o2.n_rows >= 1 && classify(o2) == Chk::BLOCKED_CAP, fmt("rows %.0f", static_cast<double>(o2.n_rows)));
        }
        {
            NodeSet ns = base;
            ns.elems[L1].depth = 6;
            const RefineOutcome o = refine_B(c, rho_syn, ns, {pc, pd}, evs, 1);
            report("V5d", "mixed check (one split + one depth refusal) -> SPLIT",
                   classify(o) == Chk::SPLIT && o.n_split == 1 && o.n_refused_depth == 1,
                   std::string("got ") + chk_name(classify(o)));
            const RefineOutcome o2 = refine_B(c, rho_syn, ns, {pc, pd}, evs, 2);
            report("V5d", "refused row re-tested at the next check -> BLOCKED_DEPTH",
                   classify(o2) == Chk::BLOCKED_DEPTH && o2.n_refused_depth == 1,
                   std::string("got ") + chk_name(classify(o2)));
        }
        {
            NodeSet ns = base;
            const RefineOutcome o = refine_B(c, rho_syn, ns, {Pt{100.0, 1.0, -2.0}}, evs, 1);
            report("V5d", "candidate outside Dw counted, does not affect C3",
                   o.n_outside_Dw == 1 && classify(o) == Chk::PASS, "");
        }
        const OuterDecision d10 = outer_decide(Chk::SPLIT, 10, 0.3, false);
        report("V5d", "SPLIT at outer round 10 -> OUTER_LIMIT", d10.term == "OUTER_LIMIT" && !d10.apply_split, d10.term);
        const OuterDecision d9 = outer_decide(Chk::SPLIT, 9, 0.3, false);
        report("V5d", "SPLIT at outer round 9 -> next inner solve", d9.term.empty() && d9.apply_split, "");
        const OuterDecision dt = outer_decide(Chk::SPLIT, 2, 3e-9, true);
        report("V5d", "after SPLIT with trust < 1e-8 the next inner solve starts at 0.5",
               dt.apply_split && dt.trust == 0.5, fmt("trust %.3g", dt.trust));
        report("V5d", "PASS -> CONVERGED; PASS after LS limit -> CONVERGED_LS_LIMITED",
               outer_decide(Chk::PASS, 1, 0.5, false).term == "CONVERGED" &&
                   outer_decide(Chk::PASS, 1, 1e-9, true).term == "CONVERGED_LS_LIMITED", "");
        std::ofstream(out + "/v5d_events.csv") << "it,x,r,rho_star,rho_L,depth,leaf,n_containing,need,action\n"
                                               << evs.str();
    }
    std::printf("SELFTEST %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails;
}

// ==== verbatim s3.cpp lines 693-1087 ====
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
    if (argc < 2) {
        std::fprintf(stderr, "usage: s4 run|hessian|fd|recon|v2dump|selftest|traj|assess|floor|axi|field|timing|check ...\n");
        return 2;
    }
    const std::string mode = argv[1];
    std::string out = ".", yfile = "BASELINE", tag = "t", ysfile, its, recon, setfile, passes, dirs, run, hdir;
    int order = 4, sub = 1, kmax = 200, diag_every = 1, line = 0;
    double lam = 0.0, deadline = 0.0;
    char form = 'A';
    unsigned seed = 1;
    bool full = false;
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        auto nx = [&]() { return std::string(argv[++i]); };
        if (a == "--out") out = nx();
        else if (a == "--form") form = nx()[0];
        else if (a == "--ws") { const std::string s = nx(); order = std::stoi(s.substr(0, s.find(':'))); sub = std::stoi(s.substr(s.find(':') + 1)); }
        else if (a == "--kmax") kmax = std::stoi(nx());
        else if (a == "--lambda") lam = std::stod(nx());
        else if (a == "--deadline") deadline = std::stod(nx());
        else if (a == "--diag-every") diag_every = std::stoi(nx());
        else if (a == "--y") yfile = nx();
        else if (a == "--ys") ysfile = nx();
        else if (a == "--line") line = std::stoi(nx());
        else if (a == "--its") its = nx();
        else if (a == "--seed") seed = static_cast<unsigned>(std::stoul(nx()));
        else if (a == "--tag") tag = nx();
        else if (a == "--full") full = nx() == "1";
        else if (a == "--recon") recon = nx();
        else if (a == "--set") setfile = nx();
        else if (a == "--passes") passes = nx();
        else if (a == "--dirs") dirs = nx();
        else if (a == "--run") run = nx();
        else if (a == "--hessian") hdir = nx();
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    const Ctx c = make_ctx();
    if (mode == "run") {
        RunArgs A;
        A.form = form; A.order = order; A.sub = sub; A.kmax = kmax; A.diag_every = diag_every; A.y0 = yfile;
        A.recon = recon; A.setfile = setfile; A.hdir = hdir; A.out = out; A.passes = passes; A.lam = lam;
        A.deadline = deadline;
        return mode_run(c, A);
    }
    if (mode == "hessian" || mode == "fd" || mode == "v2dump") {
        const Vec y = y_from_args(yfile, ysfile, line);
        const NodeSet ns = make_set(c, order, sub, recon, setfile);
        if (mode == "hessian") return mode_hessian(c, y, ns, passes, lam, out, tag);
        if (mode == "fd") return mode_fd(c, y, ns, dirs, out);
        return mode_v2dump(c, y, ns, lam);
    }
    if (mode == "recon") return mode_recon(c, run, order, sub, out);
    if (mode == "selftest") return selftest(c, passes, out);
    if (mode == "traj") return mode_traj(c, ysfile, its, out);
    if (mode == "assess") return mode_assess(c, yfile, seed, out, tag, full);
    if (mode == "floor") return mode_floor(c);
    if (mode == "axi") return mode_axi(c, yfile);
    if (mode == "field") return mode_field(c, yfile, out);
    if (mode == "timing") return mode_timing(c);
    if (mode == "check") return mode_check(c);
    return 2;
}
