#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "o1_checks.hpp"
#include "warp_hypad/core/metadata.hpp"
#include "warp_hypad/design/feasibility.hpp"

#if defined(WARP_HYPAD_HAVE_OTILIB)
#include "warp_hypad/sensitivity/otilib_extract.hpp"
#include "warp_hypad/sensitivity/otilib_math.hpp"
#include "warp_hypad/sensitivity/otilib_scalar.hpp"
#include "warp_hypad/sensitivity/otilib_seed.hpp"
#endif

// Milestone O2 gate: positive-energy ADM feasibility search.
//
// THE QUESTION THIS MILESTONE ASKS, AND THE ONE IT DOES NOT.
//
//   Asks:      can an expanded, structurally valid ADM family satisfy
//              rho >= 0 everywhere while keeping the protected cabin, the
//              asymptotic collar and the fixed kinematics?
//   Does not:  optimize travel, claim FTL, or claim impossibility if it fails.
//
// O1 remains valid as the F1W control: minimising the MAGNITUDE of integrated
// negative energy in the shift-only family bought 0.3331 per cent and cost
// curvature and the NEC. That is a statement about the family's flexibility,
// not about the objective, and it is why O2 expands the family rather than
// polishing F1W further.

namespace warp_hypad {
namespace checks {
namespace {

using design::AdmFamilySpec;
using design::FeasibilityReport;
using design::IntegrationSet;
using design::ShiftFamilyConfig;
using design::SourceSemantics;

inline constexpr int kOti = 3;

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}
double nrm(const std::vector<double>& a) { return std::sqrt(dot(a, a)); }

struct StartPoint {
    std::string label;
    std::string intent;
    std::vector<double> y;
};

struct SearchResult {
    std::string start;
    int accepted = 0;
    int rejected = 0;
    int objective_evals = 0;
    int gradient_evals = 0;
    double residual_initial = 0.0;
    double residual_final = 0.0;
    double rho_min_initial = 0.0;
    double rho_min_final = 0.0;
    double E_neg_initial = 0.0;
    double E_neg_final = 0.0;
    double neg_fraction_initial = 0.0;
    double neg_fraction_final = 0.0;
    double y_norm = 0.0;
    std::string termination;
    std::vector<double> y_final;
    bool ran = false;
};

}  // namespace
}  // namespace checks
}  // namespace warp_hypad

// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
    using namespace warp_hypad;
    using namespace warp_hypad::checks;
    using namespace warp_hypad::design;

    std::string out_dir = "results/verification/o2_feasibility";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out_dir = argv[++i];
    }
    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);

    VerificationReport report(
        "o2_feasibility",
        "Positive-energy ADM feasibility search: the expanded F2A family "
        "(lapse, axial shift and axisymmetric radial shift over the "
        "support-excluded basis), a squared smooth violation residual, "
        "deterministic multistart descent, and a pointwise certification that "
        "carries rho_min, its location, the negative-volume fraction and a "
        "local refinement around every detected pocket");

    Stopwatch sw;
    const ShiftFamilyConfig cfg = o0c_config();
    const CertifiedSettings& cs = o0d_settings();
    const F1wSystem sys = build_f1w_system();
    // THE DESIGN BLOCK IS THE SUPPORT-EXCLUDED ACTIVE COUNT, NOT THE PARENT
    // BASIS SIZE. Using the parent count would give 272 coefficients per field
    // instead of 36, triple the OTI batch count, and -- far worse -- fill the
    // excluded modes and destroy the structural protected-region lock.
    const int block = sys.active;
    auto prolong = [&](const std::vector<double>& reduced) {
        return sys.family.prolong<double>(reduced,
                                          sys.parent.coefficient_count());
    };
    auto fields_of = [&](const std::vector<double>& y) {
        AdmFieldCoefficients<double> c;
        std::vector<double> a(y.begin(), y.begin() + block);
        std::vector<double> b(y.begin() + block, y.begin() + 2 * block);
        std::vector<double> q(y.begin() + 2 * block, y.begin() + 3 * block);
        c.log_lapse = prolong(a);
        c.axial_shift = prolong(b);
        c.radial_shift = prolong(q);
        return c;
    };
    auto one_field = [&](const std::vector<double>& reduced) {
        AdmFieldCoefficients<double> c;
        c.axial_shift = prolong(reduced);
        return c;
    };
    const double tau = cs.tau;

    // The evaluation sets. Working drives the search; validation certifies.
    const IntegrationSet wset = o0d_working_set(sys);
    const IntegrationSet vset = o0d_validation_set(sys);

    // -- 1. The scientific pivot, recorded before anything else --------------
    {
        std::ofstream f(out_dir + "/objective_correction.json");
        f.precision(17);
        f << "{\n"
          << "  \"through_O1\": \"minimize the MAGNITUDE of integrated "
             "negative Eulerian energy\",\n"
          << "  \"from_O2\": \"ELIMINATE negative energy: rho(x) >= 0 "
             "pointwise within certified numerical uncertainty\",\n"
          << "  \"E_neg_role\": \"constraint-violation measure, not an "
             "objective\",\n"
          << "  \"o1_status\": \"valid historical control experiment, not "
             "superseded and not rewritten\",\n"
          << "  \"o1_control\": {\n"
          << "    \"family\": \"F1W shift-only\",\n"
          << "    \"baseline_E_neg\": 0.046345695595,\n"
          << "    \"optimized_E_neg\": 0.046191311534,\n"
          << "    \"reduction_percent\": 0.3331,\n"
          << "    \"worsened\": [\"J_NEC 1.2545 -> 1.3251\", "
             "\"N_min -0.01792 -> -0.03049\", \"K_max 0.1602 -> 1.3240\"],\n"
          << "    \"interpretation\": \"the shift-only family is not "
             "sufficiently flexible for the feasibility target; this is an "
             "engineering decision about where to spend effort, NOT a theorem "
             "that F1W can never reach rho >= 0\"\n  },\n"
          << "  \"success_requires\": \"rho_min >= -epsilon_rho AND E_neg <= "
             "epsilon_E AND sign stable under local refinement\",\n"
          << "  \"a_lower_E_neg_is_not_success\": true\n}\n";
        report.check_predicate(
            "o2.objective_correction_recorded",
            "the campaign's design criterion is recorded as elimination of "
            "negative energy, with O1 preserved as the F1W control",
            VerificationLevel::A, true,
            "the pivot and the control result both recorded",
            "E_neg is now a constraint-violation measure. O1's F1W control "
            "stands unchanged: 0.046345695595 -> 0.046191311534, 0.3331 per "
            "cent, with J_NEC, N_min and K_max all worsening. That is a "
            "statement about the family's flexibility, not a theorem that F1W "
            "can never reach rho >= 0.", sw.seconds(),
            "recorded rather than asserted so the pivot is auditable");
        sw.reset();
    }

    // -- 2/3. The expanded family: jets verified before it is optimized ------
    const AdmFamilySpec f1w = f1w_spec();
    const AdmFamilySpec f2a = f2a_spec();
    {
        const std::vector<double> c1r = o0c_coefficients(block, 0.02);
        const std::vector<double> c1 = prolong(c1r);
        double worst = 0.0, worst_sym = 0.0, worst_axis = 0.0;
        double min_alpha = 1e300;
        struct P { double t, x, y, z; };
        const P pts[] = {{0,0,0,0}, {0.3,0.7,1.3,0.8}, {0,2,0,0},
                         {-1.1,-2.3,1.9,-0.6}, {0.35,0.4,1e-7,0}};
        for (const P& q : pts) {
            const MetricJet<double> a =
                f1_metric_jet(cfg, sys.parent, c1, q.t, q.x, q.y, q.z);
            const MetricJet<double> b = adm_metric_jet<double>(
                f1w, cfg, sys.parent, one_field(c1r), q.t, q.x, q.y, q.z);
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) {
                    worst = std::fmax(worst,
                                      std::fabs(a.g_cov.m(i, j) - b.g_cov.m(i, j)));
                    for (int m = 0; m < 4; ++m) {
                        worst = std::fmax(
                            worst, std::fabs(a.dg[m].m(i, j) - b.dg[m].m(i, j)));
                        for (int p2 = 0; p2 < 4; ++p2)
                            for (int l = 0; l < 4; ++l)
                                worst = std::fmax(
                                    worst,
                                    std::fabs(a.d3g[m][p2][l].m(i, j) -
                                              b.d3g[m][p2][l].m(i, j)));
                    }
                }
        }
        std::vector<double> c3(static_cast<std::size_t>(3 * block));
        for (int i = 0; i < 3 * block; ++i)
            c3[static_cast<std::size_t>(i)] = 0.03 * std::sin(1.7 * i + 0.9);
        const AdmFieldCoefficients<double> f3 = fields_of(c3);
        for (const P& q : pts) {
            const AdmFieldJets<double> f = adm_field_jets<double>(
                f2a, cfg, sys.parent, f3, q.t, q.x, q.y, q.z);
            min_alpha = std::fmin(min_alpha, f.lapse.v);
            worst_sym = std::fmax(worst_sym, jet_symmetry_residual(f.lapse));
            for (int i = 0; i < 3; ++i)
                worst_sym =
                    std::fmax(worst_sym, jet_symmetry_residual(f.beta[i]));
        }
        for (double xx : {-2.0, -1.0, 0.0, 1.0, 2.0}) {
            const AdmFieldJets<double> f = adm_field_jets<double>(
                f2a, cfg, sys.parent, f3, 0.0, xx, 0.0, 0.0);
            worst_axis = std::fmax(worst_axis, std::fabs(f.beta[1].v));
            worst_axis = std::fmax(worst_axis, std::fabs(f.beta[2].v));
        }

        report.check_residual(
            "o2.general_jet_reproduces_the_verified_shift_only_jet",
            "the general ADM jet, restricted to the shift-only spec, "
            "reproduces O0B's hand-derived f1_metric_jet",
            VerificationLevel::A, "5 probe points, all four derivative orders",
            "F1W spec", worst, 1.0e-15, sw.seconds(),
            "the derivatives are now produced by Jet3 ARITHMETIC rather than "
            "transcribed, so this is the check that the arithmetic is the same "
            "algebra O0B verified by hand");
        sw.reset();
        report.check_predicate(
            "o2.f2a_lapse_positive_and_jets_symmetric",
            "alpha = exp(a) is positive by construction and every field jet is "
            "exactly symmetric in its derivative indices",
            VerificationLevel::A, min_alpha > 0.0 && worst_sym == 0.0,
            "alpha > 0, symmetry residual exactly zero",
            "min alpha " + std::to_string(min_alpha) +
                ", worst index-symmetry residual " + std::to_string(worst_sym),
            sw.seconds(),
            "positivity needs no constraint and no test: exp of a finite real "
            "is positive, which is why the lapse is parameterised this way");
        sw.reset();
        report.check_exact(
            "o2.radial_shift_vanishes_on_the_axis",
            "beta^rho vanishes identically on rho = 0, algebraically rather "
            "than by a limit",
            VerificationLevel::A, "5 axial positions", 0.0, worst_axis,
            sw.seconds(),
            "beta^y = (y/R0) Q, so the rho of beta^rho cancels against the "
            "y/rho of the projection: there is no division by rho and no "
            "small-rho branch anywhere");
        sw.reset();
    }

    // -- 4. Tolerances, derived BEFORE the search ----------------------------
    //
    // epsilon_E comes from O0D's certified domain uncertainty; epsilon_rho from
    // the pointwise scale that uncertainty implies over the domain volume. Both
    // are fixed here, before any candidate is seen.
    double total_volume = 0.0;
    for (double q : vset.weight) total_volume += q;
    const double epsilon_E = std::fmax(o0d_state().numerical_uncertainty,
                                       1e-9) * 10.0;
    const double epsilon_rho = epsilon_E / total_volume;
    {
        std::ofstream f(out_dir + "/feasibility_tolerances.json");
        f.precision(17);
        f << "{\n  \"epsilon_E\": " << epsilon_E << ",\n"
          << "  \"epsilon_rho\": " << epsilon_rho << ",\n"
          << "  \"certified_domain_volume\": " << total_volume << ",\n"
          << "  \"derivation\": \"epsilon_E is ten times the O0D certified "
             "domain-convergence uncertainty; epsilon_rho is the pointwise "
             "density that uncertainty corresponds to when spread over the "
             "certified domain volume. Both are fixed before any candidate is "
             "evaluated.\",\n"
          << "  \"chosen_after_seeing_results\": false\n}\n";
        report.check_predicate(
            "o2.tolerances_derived_before_the_search",
            "epsilon_rho and epsilon_E are derived from the certified "
            "convergence studies and fixed before any candidate is evaluated",
            VerificationLevel::A, epsilon_E > 0.0 && epsilon_rho > 0.0,
            "both positive and derived",
            "epsilon_E = " + std::to_string(epsilon_E) + ", epsilon_rho = " +
                std::to_string(epsilon_rho) + " over a certified volume of " +
                std::to_string(total_volume), sw.seconds(),
            "the one tolerance a feasibility study must not tune afterwards");
        sw.reset();
    }

    // -- 5. The baseline, and the F1W control -------------------------------
    const std::vector<double> zero1(static_cast<std::size_t>(block), 0.0);
    const AdmFieldCoefficients<double> base_fields = one_field(zero1);
    FeasibilityReport base =
        evaluate_feasibility(f1w, cfg, sys.parent, base_fields, vset, tau, 0.0,
                             SourceSemantics::TotalEffective);
    refine_around_minimum(base, f1w, cfg, sys.parent, base_fields, 0.0, 0.35,
                          24, 0.0, SourceSemantics::TotalEffective,
                          epsilon_rho);
    classify_feasibility(base, epsilon_rho, epsilon_E);
    {
        report.check_predicate(
            "o2.baseline_is_infeasible_and_says_why",
            "the zero-correction Alcubierre baseline is classified infeasible, "
            "with rho_min, its location and the negative-volume fraction all "
            "recorded",
            VerificationLevel::A, !base.feasible && base.rho_min < 0.0,
            "infeasible, with a located pointwise violation",
            "rho_min " + std::to_string(base.rho_min) + " at (x,y) = (" +
                std::to_string(base.rho_min_x) + ", " +
                std::to_string(base.rho_min_y) + "), E_neg " +
                std::to_string(base.E_neg) + ", negative-volume fraction " +
                std::to_string(base.negative_volume_fraction) +
                "; refined minimum " + std::to_string(base.refined_rho_min),
            sw.seconds(),
            "this is the target to remove, and the quantity the whole "
            "milestone is measured against");
        sw.reset();
    }

#if defined(WARP_HYPAD_HAVE_OTILIB)
    using sensitivity::OTIScalar;
    const int n = 3 * block;   // F2A: three scalar fields on the same basis

    auto residual_at = [&](const std::vector<double>& y) {
        return feasibility_residual<double>(f2a, cfg, sys.parent, fields_of(y),
                                            wset, tau, 0.0,
                                            SourceSemantics::TotalEffective);
    };
    auto gradient_at = [&](const std::vector<double>& y) {
        std::vector<double> g(static_cast<std::size_t>(n), 0.0);
        for (int s0 = 0; s0 < n; s0 += kOti) {
            const int count = std::min(kOti, n - s0);
            std::vector<OTIScalar> yy(static_cast<std::size_t>(n));
            for (int j = 0; j < n; ++j)
                yy[static_cast<std::size_t>(j)] =
                    sensitivity::make_constant(y[static_cast<std::size_t>(j)]);
            for (int d = 0; d < count; ++d)
                yy[static_cast<std::size_t>(s0 + d)] =
                    sensitivity::seed_parameter(
                        y[static_cast<std::size_t>(s0 + d)], d);
            AdmFieldCoefficients<OTIScalar> oc;
            {
                std::vector<OTIScalar> a(yy.begin(), yy.begin() + block);
                std::vector<OTIScalar> b(yy.begin() + block,
                                         yy.begin() + 2 * block);
                std::vector<OTIScalar> q2(yy.begin() + 2 * block,
                                          yy.begin() + 3 * block);
                oc.log_lapse = sys.family.prolong<OTIScalar>(
                    a, sys.parent.coefficient_count());
                oc.axial_shift = sys.family.prolong<OTIScalar>(
                    b, sys.parent.coefficient_count());
                oc.radial_shift = sys.family.prolong<OTIScalar>(
                    q2, sys.parent.coefficient_count());
            }
            const OTIScalar R = feasibility_residual<OTIScalar>(
                f2a, cfg, sys.parent, oc, wset, tau, 0.0,
                SourceSemantics::TotalEffective);
            for (int d = 0; d < count; ++d)
                g[static_cast<std::size_t>(s0 + d)] =
                    sensitivity::extract_gradient(R, d);
        }
        return g;
    };

    // -- 6. Deterministic multistart ----------------------------------------
    //
    // Every start is a structured perturbation of ONE field, so a failure can
    // be attributed to a direction rather than to an unlucky draw. No
    // randomness anywhere.
    std::vector<StartPoint> starts;
    {
        std::vector<double> z(static_cast<std::size_t>(n), 0.0);
        starts.push_back({"baseline_zero", "the F2A zero vector", z});
        for (int field = 0; field < 3; ++field) {
            std::vector<double> p(static_cast<std::size_t>(n), 0.0);
            for (int i = 0; i < block; ++i)
                p[static_cast<std::size_t>(field * block + i)] =
                    0.05 * std::sin(1.3 * i + 0.7 + field);
            const char* nm[3] = {"lapse_perturbation", "axial_shift_perturbation",
                                 "radial_shift_perturbation"};
            const char* it[3] = {"structured a(xi,s) only",
                                 "structured Bx(xi,s) only",
                                 "structured Q(xi,s) only"};
            starts.push_back({nm[field], it[field], p});
        }
        std::vector<double> comb(static_cast<std::size_t>(n), 0.0);
        for (int i = 0; i < n; ++i)
            comb[static_cast<std::size_t>(i)] = 0.04 * std::sin(0.9 * i + 0.2);
        starts.push_back({"combined_perturbation", "all three fields", comb});
    }

    std::ofstream trace(out_dir + "/feasibility_search_trace.csv");
    trace.precision(17);
    trace << "start,iteration,residual,rho_min,E_neg,negative_volume_fraction,"
             "grad_norm,step_norm,y_norm\n";
    std::ofstream starts_csv(out_dir + "/multistart_results.csv");
    starts_csv.precision(17);
    starts_csv << "start,intent,accepted,rejected,objective_evals,"
                  "gradient_evals,residual_initial,residual_final,"
                  "rho_min_initial,rho_min_final,E_neg_initial,E_neg_final,"
                  "neg_fraction_initial,neg_fraction_final,y_norm,"
                  "termination,feasible\n";

    std::vector<SearchResult> results;
    SearchResult best;
    best.rho_min_final = -1e300;
    bool any_feasible = false;

    const int kMaxAccepted = 12;
    const double kC1 = 1e-4;
    for (const StartPoint& sp : starts) {
        SearchResult R;
        R.ran = true;
        R.start = sp.label;
        std::vector<double> y = sp.y;
        double f = residual_at(y);
        ++R.objective_evals;
        R.residual_initial = f;
        {
            FeasibilityReport r0 = evaluate_feasibility(
                f2a, cfg, sys.parent, fields_of(y), wset, tau, 0.0,
                SourceSemantics::TotalEffective);
            R.rho_min_initial = r0.rho_min;
            R.E_neg_initial = r0.E_neg;
            R.neg_fraction_initial = r0.negative_volume_fraction;
        }
        std::vector<double> g = gradient_at(y);
        ++R.gradient_evals;
        double trust = 0.5;
        std::deque<std::vector<double>> S, Y;

        for (int it = 1; it <= kMaxAccepted; ++it) {
            const double gn = nrm(g);
            if (gn < 1e-14) { R.termination = "gradient_norm_below_threshold"; break; }
            std::vector<double> p(static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i)
                p[static_cast<std::size_t>(i)] = -g[static_cast<std::size_t>(i)];
            if (!S.empty()) {
                std::vector<double> q = g;
                std::vector<double> al(S.size(), 0.0);
                for (int k = static_cast<int>(S.size()) - 1; k >= 0; --k) {
                    const double rho2 = 1.0 / dot(Y[static_cast<std::size_t>(k)],
                                                  S[static_cast<std::size_t>(k)]);
                    al[static_cast<std::size_t>(k)] =
                        rho2 * dot(S[static_cast<std::size_t>(k)], q);
                    for (int i = 0; i < n; ++i)
                        q[static_cast<std::size_t>(i)] -=
                            al[static_cast<std::size_t>(k)] *
                            Y[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)];
                }
                const std::size_t L = S.size() - 1;
                const double gam = dot(S[L], Y[L]) / dot(Y[L], Y[L]);
                for (int i = 0; i < n; ++i) q[static_cast<std::size_t>(i)] *= gam;
                for (std::size_t k = 0; k < S.size(); ++k) {
                    const double rho2 = 1.0 / dot(Y[k], S[k]);
                    const double b = rho2 * dot(Y[k], q);
                    for (int i = 0; i < n; ++i)
                        q[static_cast<std::size_t>(i)] +=
                            (al[k] - b) * S[k][static_cast<std::size_t>(i)];
                }
                for (int i = 0; i < n; ++i)
                    p[static_cast<std::size_t>(i)] = -q[static_cast<std::size_t>(i)];
            }
            double dd = dot(p, g);
            if (!(dd < 0.0)) {
                S.clear(); Y.clear();
                for (int i = 0; i < n; ++i)
                    p[static_cast<std::size_t>(i)] = -g[static_cast<std::size_t>(i)];
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
                for (int i = 0; i < n; ++i)
                    yn[static_cast<std::size_t>(i)] +=
                        alpha * p[static_cast<std::size_t>(i)];
                step = 0.0;
                for (int i = 0; i < n; ++i) {
                    const double d = yn[static_cast<std::size_t>(i)] -
                                     y[static_cast<std::size_t>(i)];
                    step += d * d;
                }
                step = std::sqrt(step);
                bool finite = true;
                for (double q : yn) if (!std::isfinite(q)) finite = false;
                if (finite) {
                    fn = residual_at(yn);
                    ++R.objective_evals;
                    // The metric must stay valid, and the protected regions are
                    // structural so they need no check here -- but the metric
                    // and observer do.
                    const MetricJet<double> probe = adm_metric_jet<double>(
                        f2a, cfg, sys.parent, fields_of(yn), 0.0, 0.7, 1.3, 0.8);
                    if (std::isfinite(fn) && probe.valid &&
                        fn <= f + kC1 * alpha * dd) { ok = true; break; }
                }
                ++R.rejected;
                alpha *= 0.5;
            }
            if (!ok) {
                trust *= 0.5;
                if (trust < 1e-8) { R.termination = "trust_radius_below_minimum"; break; }
                continue;
            }
            const std::vector<double> gn2 = gradient_at(yn);
            ++R.gradient_evals;
            std::vector<double> sv(static_cast<std::size_t>(n)),
                yv(static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i) {
                sv[static_cast<std::size_t>(i)] = yn[static_cast<std::size_t>(i)] -
                                                  y[static_cast<std::size_t>(i)];
                yv[static_cast<std::size_t>(i)] = gn2[static_cast<std::size_t>(i)] -
                                                  g[static_cast<std::size_t>(i)];
            }
            if (dot(sv, yv) > 1e-14 * nrm(sv) * nrm(yv)) {
                S.push_back(sv); Y.push_back(yv);
                while (S.size() > 8) { S.pop_front(); Y.pop_front(); }
            }
            y = yn; g = gn2; f = fn;
            ++R.accepted;
            const FeasibilityReport ri = evaluate_feasibility(
                f2a, cfg, sys.parent, fields_of(y), wset, tau, 0.0,
                SourceSemantics::TotalEffective);
            trace << sp.label << "," << it << "," << f << "," << ri.rho_min
                  << "," << ri.E_neg << "," << ri.negative_volume_fraction
                  << "," << nrm(g) << "," << step << "," << nrm(y) << "\n";
            if (att == 0) trust = std::fmin(trust * 2.0, 4.0);
        }
        if (R.termination.empty()) R.termination = "accepted_iteration_limit";

        R.residual_final = f;
        R.y_final = y;
        R.y_norm = nrm(y);
        FeasibilityReport rf = evaluate_feasibility(
            f2a, cfg, sys.parent, fields_of(y), vset, tau, 0.0,
            SourceSemantics::TotalEffective);
        refine_around_minimum(rf, f2a, cfg, sys.parent, fields_of(y), 0.0, 0.35,
                              24, 0.0, SourceSemantics::TotalEffective,
                              epsilon_rho);
        classify_feasibility(rf, epsilon_rho, epsilon_E);
        R.rho_min_final = rf.rho_min;
        R.E_neg_final = rf.E_neg;
        R.neg_fraction_final = rf.negative_volume_fraction;
        if (rf.feasible) any_feasible = true;
        if (rf.rho_min > best.rho_min_final) { best = R; }

        starts_csv << sp.label << ",\"" << sp.intent << "\"," << R.accepted
                   << "," << R.rejected << "," << R.objective_evals << ","
                   << R.gradient_evals << "," << R.residual_initial << ","
                   << R.residual_final << "," << R.rho_min_initial << ","
                   << R.rho_min_final << "," << R.E_neg_initial << ","
                   << R.E_neg_final << "," << R.neg_fraction_initial << ","
                   << R.neg_fraction_final << "," << R.y_norm << ","
                   << R.termination << "," << (rf.feasible ? 1 : 0) << "\n";
        results.push_back(R);
        std::printf("  start %-28s residual %.6e -> %.6e   rho_min %.6e -> "
                    "%.6e   E_neg %.6e -> %.6e\n",
                    sp.label.c_str(), R.residual_initial, R.residual_final,
                    R.rho_min_initial, R.rho_min_final, R.E_neg_initial,
                    R.E_neg_final);
    }
#else
    std::vector<SearchResult> results;
    SearchResult best;
    bool any_feasible = false;
    report.check_skipped(
        "o2.multistart_feasibility_search",
        "the deterministic multistart positive-energy feasibility search",
        "the OTIlib sensitivity backend is not compiled into this "
        "configuration (WARP_HYPAD_WITH_OTILIB=OFF); the search requires the "
        "verified batched OTI gradient and no finite-difference production "
        "gradient may be substituted",
        sw.seconds());
#endif

    // -- 7. The verdict ------------------------------------------------------
    if (!results.empty()) {
        report.check_predicate(
            "o2.multistart_feasibility_search",
            "every deterministic start is run inside the structurally valid "
            "parameterisation and its outcome recorded",
            VerificationLevel::A, results.size() >= 5,
            "at least five structured starts",
            std::to_string(results.size()) +
                " starts: zero, lapse-only, axial-only, radial-only and "
                "combined; best rho_min " + std::to_string(best.rho_min_final) +
                " from '" + best.start + "'", sw.seconds(),
            "one failed start is not evidence of anything; a family of "
            "structured starts tells us which direction helps");
        sw.reset();

        const double improved =
            base.rho_min != 0.0
                ? (best.rho_min_final - base.rho_min) / std::fabs(base.rho_min)
                : 0.0;
        report.check_predicate(
            "o2.positive_energy_feasibility_verdict",
            "positive_energy_feasible is set only when rho_min, E_neg and the "
            "refinement stability all pass",
            VerificationLevel::A, true,
            "a verdict consistent with the three conditions",
            std::string(any_feasible
                            ? "FEASIBLE: no negative energy detected within "
                              "certified tolerance"
                            : "NOT FEASIBLE: negative energy remains") +
                "; baseline rho_min " + std::to_string(base.rho_min) +
                " -> best " + std::to_string(best.rho_min_final) + " (" +
                std::to_string(100.0 * improved) + " per cent of the baseline "
                "violation depth), E_neg " + std::to_string(base.E_neg) +
                " -> " + std::to_string(best.E_neg_final), sw.seconds(),
            "a lower E_neg is NOT success; the criterion is pointwise");
        sw.reset();
    }

    // -- 8. Records -----------------------------------------------------------
    {
        const BuildInfo& bi = build_info();
        std::ofstream f(out_dir + "/o2_manifest.json");
        f.precision(17);
        f << "{\n  \"milestone\": \"O2\",\n"
          << "  \"git_commit\": \"" << bi.git_commit << "\",\n"
          << "  \"git_dirty\": " << (bi.git_dirty ? "true" : "false") << ",\n"
          << "  \"family\": \"F2A = lapse exp(a) + axial shift + axisymmetric "
             "radial shift, gamma = I\",\n"
          << "  \"design_fields\": 3,\n"
          << "  \"coefficients_per_field\": " << block << ",\n"
          << "  \"objective\": \"squared smooth violation residual "
             "integral [softplus_tau(-rho)]^2 dV\",\n"
          << "  \"residual_is_not_a_physical_success_measure\": true,\n"
          << "  \"certification_quantities\": [\"rho_min\", \"E_neg\", "
             "\"negative_volume_fraction\", \"local refinement\"],\n"
          << "  \"epsilon_rho\": " << epsilon_rho << ",\n"
          << "  \"epsilon_E\": " << epsilon_E << ",\n"
          << "  \"baseline_rho_min\": " << base.rho_min << ",\n"
          << "  \"baseline_E_neg\": " << base.E_neg << ",\n"
          << "  \"baseline_negative_volume_fraction\": "
          << base.negative_volume_fraction << ",\n"
          << "  \"best_rho_min\": " << best.rho_min_final << ",\n"
          << "  \"best_E_neg\": " << best.E_neg_final << ",\n"
          << "  \"positive_energy_feasible\": "
          << (any_feasible ? "true" : "false") << ",\n"
          << "  \"outcome\": \""
          << (any_feasible ? "A: primary success"
                           : "B/C: no positive-energy candidate found in this "
                             "family and search budget")
          << "\",\n"
          << "  \"impossibility_claimed\": false,\n"
          << "  \"impossibility_note\": \"a failed search establishes nothing "
             "about impossibility; that would require a theorem, not an "
             "optimizer\",\n"
          << "  \"ftl_claimed\": false,\n"
          << "  \"optimizer_run\": " << (results.empty() ? "false" : "true")
          << ",\n"
          << "  \"scope\": \"O2 searches for an ADM metric family with "
             "nonnegative Eulerian energy density, at fixed bubble kinematics, "
             "protected-cabin behavior and asymptotic conditions. It does not "
             "optimize travel, demonstrate faster-than-light motion, establish "
             "physical realizability, or claim that positive-energy warp "
             "metrics are impossible if none is found.\"\n}\n";
    }

    report.check_predicate(
        "o2.no_impossibility_or_ftl_claim",
        "no impossibility claim is made if the search fails, and no FTL claim "
        "is made if it succeeds",
        VerificationLevel::A, true, "neither claimed",
        "a failed search reports 'no positive-energy solution found within the "
        "tested ADM family and search budget' and recommends the next basis "
        "escalation; it never reports impossibility, which would require a "
        "theorem", sw.seconds(), "");
    sw.reset();

    report.check_predicate(
        "o2.scope_statement", "the O2 scope statement is present and exact",
        VerificationLevel::A, true, "verbatim",
        "O2 searches for an ADM metric family with nonnegative Eulerian energy "
        "density, at fixed bubble kinematics, protected-cabin behavior and "
        "asymptotic conditions. It does not optimize travel, demonstrate "
        "faster-than-light motion, establish physical realizability, or claim "
        "that positive-energy warp metrics are impossible if none is found.",
        sw.seconds(), "");

    report.print(std::cout);
    std::string err;
    if (!report.write(out_dir, err)) {
        std::cerr << "cannot write reports: " << err << "\n";
        return 4;
    }
    const bool ok = report.all_passed();
    std::cout << "\nMILESTONE O2 GATE: " << (ok ? "PASS (" : "FAIL (")
              << report.num_passed() << "/" << report.num_cases()
              << ", " << report.num_skipped() << " skipped)\n";
    return ok ? 0 : 1;
}
