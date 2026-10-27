// wh_probe: independent-of-the-gate driver for the publication feasibility round.
//
// Builds against the UNMODIFIED Enterprise/warp_hypad headers + libraries at
// commit 477141d. Nothing in the repository is patched. Three modes:
//
//   o2   : replicate the O2 multistart F2A search (the L-BFGS loop is ported
//          verbatim from tests/o2_verify.cpp) with a configurable working
//          (optimizer) integration set and validation (reporting) set, then
//          evaluate every final design on a ladder of validation sets.
//   fd   : full-pipeline derivative verification. OTI (order 2, 3 directions)
//          gradient/Hessian of J(y) vs central finite differences in y with a
//          step sweep, for J in {R (O2 residual), E_neg, E_neg_tau, rho_min}.
//   eval : evaluate stored designs on a ladder of sets.
//
// Per-node values are computed in parallel (OpenMP) but SUMMED SERIALLY in node
// order, so every integral is bit-identical to the repository's serial loop.

#include <omp.h>  // must precede OTIlib headers (they include omp.h inside extern "C")
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <iostream>
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

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}
double nrm(const std::vector<double>& a) { return std::sqrt(dot(a, a)); }
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

// Per-node Eulerian density (and det gamma via double path when S == double).
template <typename S>
std::vector<S> densities(const Ctx& c, const AdmFamilySpec& spec,
                         const AdmFieldCoefficients<S>& f, const IntegrationSet& set) {
    const long n = static_cast<long>(set.x.size());
    std::vector<S> rho(static_cast<std::size_t>(n));
#pragma omp parallel for schedule(dynamic, 16)
    for (long p = 0; p < n; ++p) {
        const std::size_t q = static_cast<std::size_t>(p);
        rho[q] = family_eulerian_density<S>(spec, c.cfg, c.sys.parent, f, set.time,
                                            set.x[q], set.y[q], set.z[q], 0.0,
                                            SourceSemantics::TotalEffective);
    }
    return rho;
}

// All objectives from one density sweep. F2A has gamma = I, so dV = weight
// (checked against evaluate_feasibility's sqrt(det gamma) weighting in `o2`).
template <typename S>
struct Objs {
    S R, Eneg, Eneg_tau, rho_min;
    std::size_t argmin = 0;
    double gap12 = 0.0;   // real gap between smallest and 2nd smallest rho
    std::size_t near_zero = 0;  // nodes with |rho| < 1e-12 (E_neg kink)
    std::vector<std::size_t> ties;   // nodes within 1e-12 relative of the min
    std::vector<S> tie_rho;
};

template <typename S>
Objs<S> objectives(const Ctx& c, const std::vector<S>& y, const IntegrationSet& set) {
    const AdmFieldCoefficients<S> f = fields_of<S>(c, y);
    const std::vector<S> rho = densities<S>(c, c.f2a, f, set);
    Objs<S> o;
    o.R = scalar_cast<S>(0.0);
    o.Eneg = scalar_cast<S>(0.0);
    o.Eneg_tau = scalar_cast<S>(0.0);
    double m1 = 1e300, m2 = 1e300;
    for (std::size_t p = 0; p < rho.size(); ++p) {
        const S w = scalar_cast<S>(set.weight[p]);
        const S s = softplus(-rho[p], c.tau);
        o.R = o.R + s * s * w;
        o.Eneg_tau = o.Eneg_tau + s * w;
        const double r = real_part(rho[p]);
        // max(-rho, 0): branch on the real part -> one-sided derivative at 0.
        if (r < 0.0) o.Eneg = o.Eneg - rho[p] * w;
        if (std::fabs(r) < 1e-12) ++o.near_zero;
        if (r < m1) { m2 = m1; m1 = r; o.argmin = p; }
        else if (r < m2) m2 = r;
    }
    o.rho_min = rho[o.argmin];   // derivative of the ACTIVE node only
    o.gap12 = m2 - m1;
    for (std::size_t p = 0; p < rho.size(); ++p)
        if (real_part(rho[p]) - m1 <= 1e-10 * std::fabs(m1)) { o.ties.push_back(p); o.tie_rho.push_back(rho[p]); }
    return o;
}

struct SetSpec { int order, sub; double xe, se; };

IntegrationSet make_set(const Ctx& c, const SetSpec& s) {
    return build_elementwise_set(c.sys.parent, s.order, s.sub, s.xe, s.se, 0.0,
                                 c.cfg.velocity);
}

Ctx make_ctx() {
    Ctx c{o0c_config(), build_f1w_system(), f1w_spec(), f2a_spec(), 0, 1e-3};
    c.f1w = f1w_spec();
    c.f2a = f2a_spec();
    c.block = c.sys.active;
    c.tau = o0d_settings().tau;
    return c;
}

std::vector<double> read_vec(const std::string& path) {
    std::ifstream f(path);
    std::vector<double> v;
    double x;
    while (f >> x) v.push_back(x);
    return v;
}
void write_vec(const std::string& path, const std::vector<double>& v) {
    std::ofstream f(path);
    f.precision(17);
    for (double x : v) f << x << "\n";
}

// Gradient via batched OTI exactly as o2_verify does (3 coordinates per pass).
std::vector<double> gradient_R(const Ctx& c, const std::vector<double>& y,
                               const IntegrationSet& set) {
    const int n = static_cast<int>(y.size());
    std::vector<double> g(static_cast<std::size_t>(n), 0.0);
    for (int s0 = 0; s0 < n; s0 += 3) {
        const int count = std::min(3, n - s0);
        std::vector<OTIScalar> yy(static_cast<std::size_t>(n));
        for (int j = 0; j < n; ++j) yy[static_cast<std::size_t>(j)] = sensitivity::make_constant(y[static_cast<std::size_t>(j)]);
        for (int d = 0; d < count; ++d)
            yy[static_cast<std::size_t>(s0 + d)] = sensitivity::seed_parameter(y[static_cast<std::size_t>(s0 + d)], d);
        const AdmFieldCoefficients<OTIScalar> f = fields_of<OTIScalar>(c, yy);
        const std::vector<OTIScalar> rho = densities<OTIScalar>(c, c.f2a, f, set);
        OTIScalar acc = scalar_cast<OTIScalar>(0.0);
        for (std::size_t p = 0; p < rho.size(); ++p) {
            const OTIScalar s = softplus(-rho[p], c.tau);
            acc = acc + s * s * scalar_cast<OTIScalar>(set.weight[p]);
        }
        for (int d = 0; d < count; ++d)
            g[static_cast<std::size_t>(s0 + d)] = sensitivity::extract_gradient(acc, d);
    }
    return g;
}

double residual_R(const Ctx& c, const std::vector<double>& y, const IntegrationSet& set) {
    return real_part(objectives<double>(c, y, set).R);
}

// ---------------------------------------------------------------------------
int mode_o2(const Ctx& c, SetSpec ws, SetSpec vs, int max_acc, const std::string& out,
            const std::vector<SetSpec>& ladder, const std::string& only) {
    const int block = c.block;
    const int n = 3 * block;
    const IntegrationSet wset = make_set(c, ws);
    const IntegrationSet vset = make_set(c, vs);
    std::printf("working set: order %d sub %d -> %zu nodes; validation: order %d sub %d -> %zu nodes; n=%d\n",
                ws.order, ws.sub, wset.x.size(), vs.order, vs.sub, vset.x.size(), n);
    double total_volume = 0.0;
    for (double q : vset.weight) total_volume += q;
    // o2_verify: o0d_state() is never populated in that process, so
    // numerical_uncertainty == 0 and epsilon_E falls to the 1e-9 floor * 10.
    const double epsilon_E = std::fmax(o0d_state().numerical_uncertainty, 1e-9) * 10.0;
    const double epsilon_rho = epsilon_E / total_volume;

    // Baseline exactly as the gate: F1W spec, zero coefficients, vset.
    AdmFieldCoefficients<double> bf;
    bf.axial_shift = c.sys.family.prolong<double>(std::vector<double>(static_cast<std::size_t>(block), 0.0),
                                                  c.sys.parent.coefficient_count());
    FeasibilityReport base = evaluate_feasibility(c.f1w, c.cfg, c.sys.parent, bf, vset, c.tau, 0.0,
                                                  SourceSemantics::TotalEffective);
    refine_around_minimum(base, c.f1w, c.cfg, c.sys.parent, bf, 0.0, 0.35, 24, 0.0,
                          SourceSemantics::TotalEffective, epsilon_rho);
    std::printf("baseline rho_min %.6e  E_neg %.6e  negfrac %.4f  refined %.6e  eps_E %.3e eps_rho %.3e\n",
                base.rho_min, base.E_neg, base.negative_volume_fraction, base.refined_rho_min,
                epsilon_E, epsilon_rho);

    struct Start { std::string label; std::vector<double> y; };
    std::vector<Start> starts;
    starts.push_back({"baseline_zero", std::vector<double>(static_cast<std::size_t>(n), 0.0)});
    const char* nm[3] = {"lapse_perturbation", "axial_shift_perturbation", "radial_shift_perturbation"};
    for (int field = 0; field < 3; ++field) {
        std::vector<double> p(static_cast<std::size_t>(n), 0.0);
        for (int i = 0; i < block; ++i)
            p[static_cast<std::size_t>(field * block + i)] = 0.05 * std::sin(1.3 * i + 0.7 + field);
        starts.push_back({nm[field], p});
    }
    std::vector<double> comb(static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) comb[static_cast<std::size_t>(i)] = 0.04 * std::sin(0.9 * i + 0.2);
    starts.push_back({"combined_perturbation", comb});

    std::ofstream sum(out + "/o2_starts.csv");
    sum.precision(10);
    sum << "start,accepted,rejected,R_init,R_final,rho_min_init_w,rho_min_final_v,refined_rho_min,"
           "E_neg_init_w,E_neg_final_v,negfrac_final,termination,feasible,seconds\n";
    std::ofstream trace(out + "/o2_trace.csv");
    trace.precision(12);
    trace << "start,iteration,residual,rho_min_w,E_neg_w,negfrac_w,grad_norm,step_norm,y_norm,"
             "rho_min_v,refined_rho_min_v,E_neg_v\n";

    const double kC1 = 1e-4;
    double best_rho = -1e300; std::string best_label; double best_E = 0;
    for (const Start& sp : starts) {
        if (!only.empty() && sp.label != only) continue;
        const double t0 = now();
        std::vector<double> y = sp.y;
        double f = residual_R(c, y, wset);
        FeasibilityReport r0 = evaluate_feasibility(c.f2a, c.cfg, c.sys.parent, fields_of<double>(c, y), wset,
                                                    c.tau, 0.0, SourceSemantics::TotalEffective);
        const double Rinit = f;
        std::vector<double> g = gradient_R(c, y, wset);
        double trust = 0.5;
        std::deque<std::vector<double>> S, Y;
        int accepted = 0, rejected = 0;
        std::string term;
        for (int it = 1; it <= max_acc; ++it) {
            const double gn = nrm(g);
            if (gn < 1e-14) { term = "gradient_norm_below_threshold"; break; }
            std::vector<double> p(static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i) p[static_cast<std::size_t>(i)] = -g[static_cast<std::size_t>(i)];
            if (!S.empty()) {
                std::vector<double> q = g;
                std::vector<double> al(S.size(), 0.0);
                for (int k = static_cast<int>(S.size()) - 1; k >= 0; --k) {
                    const std::size_t K = static_cast<std::size_t>(k);
                    const double rho2 = 1.0 / dot(Y[K], S[K]);
                    al[K] = rho2 * dot(S[K], q);
                    for (int i = 0; i < n; ++i) q[static_cast<std::size_t>(i)] -= al[K] * Y[K][static_cast<std::size_t>(i)];
                }
                const std::size_t L = S.size() - 1;
                const double gam = dot(S[L], Y[L]) / dot(Y[L], Y[L]);
                for (int i = 0; i < n; ++i) q[static_cast<std::size_t>(i)] *= gam;
                for (std::size_t k = 0; k < S.size(); ++k) {
                    const double rho2 = 1.0 / dot(Y[k], S[k]);
                    const double b = rho2 * dot(Y[k], q);
                    for (int i = 0; i < n; ++i) q[static_cast<std::size_t>(i)] += (al[k] - b) * S[k][static_cast<std::size_t>(i)];
                }
                for (int i = 0; i < n; ++i) p[static_cast<std::size_t>(i)] = -q[static_cast<std::size_t>(i)];
            }
            double dd = dot(p, g);
            if (!(dd < 0.0)) {
                S.clear(); Y.clear();
                for (int i = 0; i < n; ++i) p[static_cast<std::size_t>(i)] = -g[static_cast<std::size_t>(i)];
                dd = dot(p, g);
            }
            const double pn = nrm(p);
            double alpha = pn > 0.0 ? std::fmin(1.0, trust / pn) : 0.0;
            bool ok = false;
            std::vector<double> yn;
            double fn = 0.0, step = 0.0;
            int att = 0;
            for (; att < 10; ++att) {
                yn = y;
                for (int i = 0; i < n; ++i) yn[static_cast<std::size_t>(i)] += alpha * p[static_cast<std::size_t>(i)];
                step = 0.0;
                for (int i = 0; i < n; ++i) { const double d = yn[static_cast<std::size_t>(i)] - y[static_cast<std::size_t>(i)]; step += d * d; }
                step = std::sqrt(step);
                bool finite = true;
                for (double q : yn) if (!std::isfinite(q)) finite = false;
                if (finite) {
                    fn = residual_R(c, yn, wset);
                    const MetricJet<double> probe = adm_metric_jet<double>(c.f2a, c.cfg, c.sys.parent,
                                                                          fields_of<double>(c, yn), 0.0, 0.7, 1.3, 0.8);
                    if (std::isfinite(fn) && probe.valid && fn <= f + kC1 * alpha * dd) { ok = true; break; }
                }
                ++rejected;
                alpha *= 0.5;
            }
            if (!ok) {
                trust *= 0.5;
                if (trust < 1e-8) { term = "trust_radius_below_minimum"; break; }
                continue;
            }
            const std::vector<double> gn2 = gradient_R(c, yn, wset);
            std::vector<double> sv(static_cast<std::size_t>(n)), yv(static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i) {
                sv[static_cast<std::size_t>(i)] = yn[static_cast<std::size_t>(i)] - y[static_cast<std::size_t>(i)];
                yv[static_cast<std::size_t>(i)] = gn2[static_cast<std::size_t>(i)] - g[static_cast<std::size_t>(i)];
            }
            if (dot(sv, yv) > 1e-14 * nrm(sv) * nrm(yv)) {
                S.push_back(sv); Y.push_back(yv);
                while (S.size() > 8) { S.pop_front(); Y.pop_front(); }
            }
            y = yn; g = gn2; f = fn;
            ++accepted;
            const FeasibilityReport ri = evaluate_feasibility(c.f2a, c.cfg, c.sys.parent, fields_of<double>(c, y), wset,
                                                              c.tau, 0.0, SourceSemantics::TotalEffective);
            // Added diagnostic (not in o2_verify): the validation-set view of
            // every accepted iterate, plus the gate's own local refinement.
            FeasibilityReport rv = evaluate_feasibility(c.f2a, c.cfg, c.sys.parent, fields_of<double>(c, y), vset,
                                                        c.tau, 0.0, SourceSemantics::TotalEffective);
            refine_around_minimum(rv, c.f2a, c.cfg, c.sys.parent, fields_of<double>(c, y), 0.0, 0.35, 24, 0.0,
                                  SourceSemantics::TotalEffective, epsilon_rho);
            trace << sp.label << "," << it << "," << f << "," << ri.rho_min << "," << ri.E_neg << ","
                  << ri.negative_volume_fraction << "," << nrm(g) << "," << step << "," << nrm(y) << ","
                  << rv.rho_min << "," << rv.refined_rho_min << "," << rv.E_neg << "\n";
            trace.flush();
            if (att == 0) trust = std::fmin(trust * 2.0, 4.0);
        }
        if (term.empty()) term = "accepted_iteration_limit";
        FeasibilityReport rf = evaluate_feasibility(c.f2a, c.cfg, c.sys.parent, fields_of<double>(c, y), vset, c.tau,
                                                    0.0, SourceSemantics::TotalEffective);
        refine_around_minimum(rf, c.f2a, c.cfg, c.sys.parent, fields_of<double>(c, y), 0.0, 0.35, 24, 0.0,
                              SourceSemantics::TotalEffective, epsilon_rho);
        classify_feasibility(rf, epsilon_rho, epsilon_E);
        const double secs = now() - t0;
        sum << sp.label << "," << accepted << "," << rejected << "," << Rinit << "," << f << "," << r0.rho_min << ","
            << rf.rho_min << "," << rf.refined_rho_min << "," << r0.E_neg << "," << rf.E_neg << ","
            << rf.negative_volume_fraction << "," << term << "," << (rf.feasible ? 1 : 0) << "," << secs << "\n";
        sum.flush();
        std::printf("  %-26s R %.4e -> %.4e  rho_min(v) %.6e  refined %.6e  E_neg(v) %.6e  acc %d  %s  %.0fs\n",
                    sp.label.c_str(), Rinit, f, rf.rho_min, rf.refined_rho_min, rf.E_neg, accepted, term.c_str(), secs);
        std::fflush(stdout);
        write_vec(out + "/y_" + sp.label + ".txt", y);
        if (rf.rho_min > best_rho) { best_rho = rf.rho_min; best_label = sp.label; best_E = rf.E_neg; }
    }
    std::printf("BASELINE rho_min %.6e E_neg %.6e | BEST(%s) rho_min %.6e E_neg %.6e | dE %.2f%% drho %.2f%%\n",
                base.rho_min, base.E_neg, best_label.c_str(), best_rho, best_E,
                100.0 * (best_E - base.E_neg) / base.E_neg,
                100.0 * (best_rho - base.rho_min) / std::fabs(base.rho_min));
    std::ofstream bfo(out + "/best.txt");
    bfo.precision(12);
    bfo << best_label << " " << base.rho_min << " " << base.E_neg << " " << best_rho << " " << best_E << "\n";
    (void)ladder;
    return 0;
}

// Evaluate designs on a ladder of sets: E_neg, node rho_min, refined rho_min,
// plus a dense global (x, y) scan at z = 0 for the continuum minimum.
int mode_eval(const Ctx& c, const std::vector<std::string>& files, const std::vector<SetSpec>& ladder,
              const std::string& out, int scan_n) {
    std::ofstream o(out, std::ios::app);
    o.precision(12);
    for (const std::string& fpath : files) {
        std::vector<double> y;
        bool baseline = (fpath == "BASELINE");
        if (!baseline) y = read_vec(fpath);
        else y.assign(static_cast<std::size_t>(3 * c.block), 0.0);
        const AdmFieldCoefficients<double> f = fields_of<double>(c, y);
        for (const SetSpec& s : ladder) {
            const double t0 = now();
            const IntegrationSet set = make_set(c, s);
            const Objs<double> ob = objectives<double>(c, y, set);
            o << fpath << "," << s.order << "," << s.sub << "," << s.xe << "," << s.se << "," << set.x.size() << ","
              << ob.Eneg << "," << ob.rho_min << "," << set.x[ob.argmin] << "," << set.y[ob.argmin] << ","
              << set.z[ob.argmin] << "," << ob.R << "," << ob.Eneg_tau << "," << (now() - t0) << "\n";
            o.flush();
            std::printf("%s o%d s%d: nodes %zu E_neg %.10e rho_min %.8e at (%.3f,%.3f,%.3f)\n", fpath.c_str(),
                        s.order, s.sub, set.x.size(), ob.Eneg, ob.rho_min, set.x[ob.argmin], set.y[ob.argmin],
                        set.z[ob.argmin]);
            std::fflush(stdout);
        }
        if (scan_n > 0) {
            // Dense uniform scan of the meridian half-plane z = 0, y >= 0, then
            // a golden local polish (coordinate search) from the best lattice node.
            const double xl = -5.0, xh = 5.0;
            const double yh = 5.0;
            std::vector<double> vals(static_cast<std::size_t>((scan_n + 1) * (scan_n + 1)));
#pragma omp parallel for schedule(dynamic, 8)
            for (int i = 0; i <= scan_n; ++i)
                for (int j = 0; j <= scan_n; ++j) {
                    const double x = xl + (xh - xl) * i / scan_n;
                    const double yy = yh * j / scan_n;
                    vals[static_cast<std::size_t>(i * (scan_n + 1) + j)] =
                        family_eulerian_density<double>(c.f2a, c.cfg, c.sys.parent, f, 0.0, x, yy, 0.0, 0.0,
                                                        SourceSemantics::TotalEffective);
                }
            std::size_t am = 0;
            for (std::size_t k = 1; k < vals.size(); ++k) if (vals[k] < vals[am]) am = k;
            double bx = xl + (xh - xl) * static_cast<double>(am / static_cast<std::size_t>(scan_n + 1)) / scan_n;
            double by = yh * static_cast<double>(am % static_cast<std::size_t>(scan_n + 1)) / scan_n;
            double bv = vals[am];
            double h = (xh - xl) / scan_n;
            auto rho_at = [&](double x, double yy) {
                return family_eulerian_density<double>(c.f2a, c.cfg, c.sys.parent, f, 0.0, x, std::fmax(yy, 0.0), 0.0,
                                                       0.0, SourceSemantics::TotalEffective);
            };
            while (h > 1e-7) {
                bool moved = false;
                for (int dx = -1; dx <= 1; ++dx)
                    for (int dy = -1; dy <= 1; ++dy) {
                        if (!dx && !dy) continue;
                        const double v = rho_at(bx + dx * h, by + dy * h);
                        if (v < bv) { bv = v; bx += dx * h; by = std::fmax(by + dy * h, 0.0); moved = true; }
                    }
                if (!moved) h *= 0.5;
            }
            o << fpath << ",scan" << scan_n << ",,,,," << ",," << bv << "," << bx << "," << by << ",0,,,\n";
            std::printf("%s continuum scan(%d)+polish: rho_min %.8e at (%.5f, %.5f, 0)\n", fpath.c_str(), scan_n, bv,
                        bx, by);
            std::fflush(stdout);
        }
    }
    return 0;
}

// Full-pipeline derivative verification along three directions V (n x 3).
int mode_fd(const Ctx& c, const std::vector<double>& y0, const std::vector<std::vector<double>>& V,
            const SetSpec& ss, const std::string& out, const std::string& tag) {
    const IntegrationSet set = make_set(c, ss);
    const int n = static_cast<int>(y0.size());
    // OTI: y = y0 + sum_d e_d V_d
    std::vector<OTIScalar> yy(static_cast<std::size_t>(n));
    for (int j = 0; j < n; ++j) {
        OTIScalar v = sensitivity::make_constant(y0[static_cast<std::size_t>(j)]);
        for (int d = 0; d < 3; ++d) {
            const double a = V[static_cast<std::size_t>(d)][static_cast<std::size_t>(j)];
            if (a != 0.0) v = v + sensitivity::seed_parameter(0.0, d) * scalar_cast<OTIScalar>(a);
        }
        yy[static_cast<std::size_t>(j)] = v;
    }
    const double t0 = now();
    const Objs<OTIScalar> oo = objectives<OTIScalar>(c, yy, set);
    const double toti = now() - t0;
    const double t1 = now();
    const Objs<double> od = objectives<double>(c, y0, set);
    const double tdbl = now() - t1;
    std::printf("[%s] set o%d s%d nodes %zu | OTI eval %.2fs, double eval %.2fs | rho_min gap12 %.3e, |rho|<1e-12 nodes %zu\n",
                tag.c_str(), ss.order, ss.sub, set.x.size(), toti, tdbl, od.gap12, od.near_zero);


    auto pick = [](const Objs<double>& o, int k) {
        switch (k) { case 0: return o.R; case 1: return o.Eneg; case 2: return o.Eneg_tau; default: return o.rho_min; }
    };
    auto pickO = [](const Objs<OTIScalar>& o, int k) {
        switch (k) { case 0: return o.R; case 1: return o.Eneg; case 2: return o.Eneg_tau; default: return o.rho_min; }
    };
    const char* names[4] = {"R_feas", "E_neg", "E_neg_tau", "rho_min"};
    std::ofstream f(out, std::ios::app);
    f.precision(17);
    // header: tag,objective,kind,i,j,h,oti,fd,relerr
    for (int k = 0; k < 4; ++k) {
        const OTIScalar Jo = pickO(oo, k);
        if (std::fabs(Jo.real() - pick(od, k)) > 1e-14 * std::fmax(1.0, std::fabs(pick(od, k))))
            std::printf("  WARNING real-part mismatch %s: %.17e vs %.17e\n", names[k], Jo.real(), pick(od, k));
        double g[3], H[3][3];
        for (int d = 0; d < 3; ++d) g[d] = sensitivity::extract_gradient(Jo, d);
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) H[a][b] = sensitivity::extract_hessian(Jo, a, b);
        f << tag << "," << names[k] << ",value,,," << "," << Jo.real() << ",,\n";
        for (int d = 0; d < 3; ++d) f << tag << "," << names[k] << ",oti_grad," << d << ",,," << g[d] << ",,\n";
        for (int a = 0; a < 3; ++a) for (int b = a; b < 3; ++b)
            f << tag << "," << names[k] << ",oti_hess," << a << "," << b << ",," << H[a][b] << ",,\n";
    }
    // rho_min is a min over nodes: with several (near-)tied active nodes it is
    // only directionally differentiable. d+ = min_k g_k.v, d- = max_k g_k.v.
    for (std::size_t t = 0; t < oo.ties.size(); ++t) {
        const std::size_t p = oo.ties[t];
        std::printf("  tie node %zu at (%.5f,%.5f,%.5f) rho %.12e  dir-derivs:", p, set.x[p], set.y[p], set.z[p],
                    oo.tie_rho[t].real());
        for (int d = 0; d < 3; ++d) std::printf(" % .8e", sensitivity::extract_gradient(oo.tie_rho[t], d));
        std::printf("\n");
        for (int d = 0; d < 3; ++d)
            f << tag << ",rho_min,tie_grad," << d << "," << p << ",," << sensitivity::extract_gradient(oo.tie_rho[t], d) << ",,\n";
    }
    // FD sweep. Shared evaluations for all 4 objectives.
    const double hs[] = {1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9};
    auto shifted = [&](double ha, int a, double hb, int b) {
        std::vector<double> y = y0;
        for (int j = 0; j < n; ++j) {
            const std::size_t J2 = static_cast<std::size_t>(j);
            if (a >= 0) y[J2] += ha * V[static_cast<std::size_t>(a)][J2];
            if (b >= 0) y[J2] += hb * V[static_cast<std::size_t>(b)][J2];
        }
        return objectives<double>(c, y, set);
    };
    const double t2 = now();
    for (double h : hs) {
        Objs<double> P[3], M[3];
        for (int a = 0; a < 3; ++a) { P[a] = shifted(h, a, 0, -1); M[a] = shifted(-h, a, 0, -1); }
        // mixed (0,1),(0,2),(1,2): 4-point
        Objs<double> pp[3], pm[3], mp[3], mm[3];
        const int pa[3] = {0, 0, 1}, pb[3] = {1, 2, 2};
        for (int q = 0; q < 3; ++q) {
            pp[q] = shifted(h, pa[q], h, pb[q]); pm[q] = shifted(h, pa[q], -h, pb[q]);
            mp[q] = shifted(-h, pa[q], h, pb[q]); mm[q] = shifted(-h, pa[q], -h, pb[q]);
        }
        for (int k = 0; k < 4; ++k) {
            const OTIScalar Jo = pickO(oo, k);
            const double J0 = pick(od, k);
            for (int a = 0; a < 3; ++a) {
                const double gfd = (pick(P[a], k) - pick(M[a], k)) / (2 * h);
                const double go = sensitivity::extract_gradient(Jo, a);
                const double hfd = (pick(P[a], k) - 2 * J0 + pick(M[a], k)) / (h * h);
                const double ho = sensitivity::extract_hessian(Jo, a, a);
                // argmin identity at +/- h, so rho_min branch switches are visible
                f << tag << "," << names[k] << ",grad," << a << ",," << h << "," << go << "," << gfd << ","
                  << std::fabs(gfd - go) / std::fmax(std::fabs(go), 1e-300) << "," << P[a].argmin << "," << M[a].argmin << "\n";
                f << tag << "," << names[k] << ",fwd," << a << ",," << h << "," << go << "," << (pick(P[a], k) - J0) / h
                  << ",," << P[a].argmin << ",\n";
                f << tag << "," << names[k] << ",bwd," << a << ",," << h << "," << go << "," << (J0 - pick(M[a], k)) / h
                  << ",,," << M[a].argmin << "\n";
                f << tag << "," << names[k] << ",hess," << a << "," << a << "," << h << "," << ho << "," << hfd << ","
                  << std::fabs(hfd - ho) / std::fmax(std::fabs(ho), 1e-300) << "," << P[a].argmin << "," << M[a].argmin << "\n";
            }
            for (int q = 0; q < 3; ++q) {
                const double hfd = (pick(pp[q], k) - pick(pm[q], k) - pick(mp[q], k) + pick(mm[q], k)) / (4 * h * h);
                const double ho = sensitivity::extract_hessian(Jo, pa[q], pb[q]);
                f << tag << "," << names[k] << ",hess," << pa[q] << "," << pb[q] << "," << h << "," << ho << "," << hfd
                  << "," << std::fabs(hfd - ho) / std::fmax(std::fabs(ho), 1e-300) << ",,\n";
            }
        }
        f.flush();
    }
    std::printf("[%s] FD sweep %.1fs; argmin(y0)=%zu at (%.4f,%.4f,%.4f)\n", tag.c_str(), now() - t2, od.argmin,
                set.x[od.argmin], set.y[od.argmin], set.z[od.argmin]);
    return 0;
}

SetSpec parse_set(const std::string& s) {
    // "order:sub:xi_extent:s_extent"
    SetSpec r{8, 1, 3.5, 12.25};
    std::stringstream ss(s);
    std::string t;
    std::vector<std::string> p;
    while (std::getline(ss, t, ':')) p.push_back(t);
    if (p.size() > 0) r.order = std::stoi(p[0]);
    if (p.size() > 1) r.sub = std::stoi(p[1]);
    if (p.size() > 2) r.xe = std::stod(p[2]);
    if (p.size() > 3) r.se = std::stod(p[3]);
    return r;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: wh_probe o2|eval|fd ...\n"); return 2; }
    const std::string mode = argv[1];
    std::string out = ".", tag = "run";
    SetSpec ws{4, 1, 1.7, 2.6}, vs{8, 1, 3.5, 12.25}, fs{4, 1, 1.7, 2.6};
    int max_acc = 12, scan = 0;
    std::vector<SetSpec> ladder;
    std::vector<std::string> files;
    std::string y0file, dirseed = "1", only;
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        auto nx = [&]() { return std::string(argv[++i]); };
        if (a == "--out") out = nx();
        else if (a == "--ws") ws = parse_set(nx());
        else if (a == "--vs") vs = parse_set(nx());
        else if (a == "--fs") fs = parse_set(nx());
        else if (a == "--maxacc") max_acc = std::stoi(nx());
        else if (a == "--ladder") ladder.push_back(parse_set(nx()));
        else if (a == "--y") files.push_back(nx());
        else if (a == "--y0") y0file = nx();
        else if (a == "--tag") tag = nx();
        else if (a == "--dirs") dirseed = nx();
        else if (a == "--scan") scan = std::stoi(nx());
        else if (a == "--only") only = nx();
    }
    const Ctx c = make_ctx();
    std::printf("block (active coeffs/field) = %d, tau = %g, parent coeffs = %d\n", c.block, c.tau,
                c.sys.parent.coefficient_count());
    if (mode == "o2") return mode_o2(c, ws, vs, max_acc, out, ladder, only);
    if (mode == "eval") return mode_eval(c, files, ladder, out, scan);
    if (mode == "fd") {
        const int n = 3 * c.block;
        std::vector<double> y0(static_cast<std::size_t>(n), 0.0);
        if (!y0file.empty() && y0file != "ZERO") y0 = read_vec(y0file);
        // Directions: d0 = one lapse coefficient, d1 = one axial-shift coefficient,
        // d2 = a deterministic dense unit vector over all 3 fields.
        std::vector<std::vector<double>> V(3, std::vector<double>(static_cast<std::size_t>(n), 0.0));
        const int k = std::stoi(dirseed);
        V[0][static_cast<std::size_t>((7 * k) % c.block)] = 1.0;
        V[1][static_cast<std::size_t>(c.block + (11 * k) % c.block)] = 1.0;
        double s2 = 0;
        for (int j = 0; j < n; ++j) { const double v = std::sin(2.3 * j + 0.4 * k + 0.1); V[2][static_cast<std::size_t>(j)] = v; s2 += v * v; }
        for (double& v : V[2]) v /= std::sqrt(s2);
        std::printf("dirs: e_%d (lapse), e_%d (axial), dense unit\n", (7 * k) % c.block, c.block + (11 * k) % c.block);
        return mode_fd(c, y0, V, fs, out, tag);
    }
    return 2;
}
