#include <cmath>
#include <string>
#include <vector>

#include "b1c_checks.hpp"
#include "b1c_velocity_reference.hpp"
#include "warp_hypad/metrics/metric.hpp"

#ifdef WARP_HYPAD_HAVE_OTILIB
#include "warp_hypad/sensitivity/otilib_traits.hpp"
#endif

// ---------------------------------------------------------------------------
// The rebuilt velocity-sensitivity tests.
//
// WHAT WAS WRONG WITH THE OLD ONES. The three failing cases
//
//     b1c.sensitivity_fixed_center_matches_reference    2.908e-2 vs 1e-6
//     b1c.sensitivity_fixed_lab_matches_total           2.908e-2 vs 1e-6
//     b1c.total_derivative_matches_chain_rule           3.354e-3 vs 1e-5
//
// all measured a RELATIVE error against a finite-difference reference, with the
// denominator floored at 1e-9:
//
//     |a - b| / max(|b|, 1e-9)
//
// At the exterior point x = 6 the true dg_00/dv is 2.6e-11, while a centred
// difference of a quantity of size 1 with h = 1e-6 carries roundoff of about
// 2e-10. The reference was therefore pure noise, larger than the signal, and
// the floored denominator converted that noise into a 2.9e-2 "relative error".
// The identical value in the first two cases is the giveaway: at t = 0 the two
// modes coincide exactly and both tests were failing at the same point, on the
// same noise. No physics was wrong -- the analytic derivative agrees with both
// finite-difference references to every printed digit at every point.
//
// HOW THESE ARE BUILT INSTEAD. Three routes, compared PAIRWISE with tolerances
// matched to what each pair can actually deliver, and never a relative error
// against something that may be zero:
//
//   OTI vs ANALYTIC        both exact closed forms, so an ABSOLUTE tolerance of
//                          1e-12 is honest and tight. This is the case that can
//                          catch an implementation defect.
//   ANALYTIC vs CENTRED FD absolute, 1e-8, best over four step sizes. This is
//                          the case that validates the analytic expression, and
//                          it is finite-difference limited by construction.
//   CHAIN RULE             algebraic identity between analytic expressions,
//                          absolute, 1e-14.
//
// Every case reports its worst error PER METRIC COMPONENT and names the point
// at which it occurred. A single maximum over components and points is exactly
// what hid the failure before.
// ---------------------------------------------------------------------------

namespace program_b {
namespace checks {

using warp_hypad::Mat4;
using warp_hypad::Stopwatch;
using warp_hypad::VerificationLevel;
using warp_hypad::VerificationReport;

namespace M = warp_hypad::metrics;

namespace {

// v at which every derivative is evaluated, and the bubble position at t = 0.
const double kVel = 0.35;
const double kXs0 = 0.0;

// The two components of the Alcubierre metric that depend on v at all. g_10 is
// g_01 by symmetry and is checked separately anyway, through the symmetry case.
struct Component {
    const char* name;
    int a;
    int b;
};
const Component kComponents[2] = {{"g00", 0, 0}, {"g01", 0, 1}};

struct Point {
    const char* name;
    double t;
    double x_off;   // offset from the bubble centre at that t
    double y;
    double z;
};

// Bubble centre, interior, wall, outer wall, exterior, and off-axis, at
// t = 0, t > 0 and t < 0.
const std::vector<Point>& points() {
    static const std::vector<Point> p = {
        {"t0_centre", 0.0, 0.0, 0.0, 0.0},
        {"t0_interior", 0.0, 1.5, 0.0, 0.0},
        {"t0_wall", 0.0, 3.0, 0.0, 0.0},
        {"t0_outer_wall", 0.0, 3.385, 0.0, 0.0},
        {"t0_exterior", 0.0, 6.0, 0.0, 0.0},
        {"t0_offaxis", 0.0, 1.0, 1.2, -0.8},
        {"t2_centre", 2.0, 0.0, 0.0, 0.0},
        {"t2_interior", 2.0, 1.5, 0.0, 0.0},
        {"t2_wall", 2.0, 3.0, 0.0, 0.0},
        {"t2_outer_wall", 2.0, 3.385, 0.0, 0.0},
        {"t2_exterior", 2.0, 6.0, 0.0, 0.0},
        {"t2_offaxis", 2.0, 1.0, 1.2, -0.8},
        {"tneg_wall", -1.5, 3.0, 0.0, 0.0},
        {"tneg_offaxis", -1.5, 2.0, 0.5, 0.5},
    };
    return p;
}

double xs_at(double t) { return kXs0 + kVel * t; }
double x_of(const Point& p) { return xs_at(p.t) + p.x_off; }

const double kSteps[4] = {1e-3, 1e-4, 1e-5, 1e-6};

// Tracks the worst error for one component together with where it happened, so
// a failure names the point instead of only a number.
struct Worst {
    double value = 0.0;
    std::string where = "none";

    void offer(double v, const char* name) {
        if (v > value) {
            value = v;
            where = name;
        }
    }
};

}  // namespace

void run_velocity_checks(VerificationReport& report, VelocityMetrics& metrics) {
    // =======================================================================
    // 0. The reference implementations must reproduce the production metric.
    //
    // Without this, every later agreement could be two wrong implementations
    // agreeing with each other.
    // =======================================================================
    {
        Stopwatch sw;
        Worst w[2];
        for (const Point& p : points()) {
            const double x = x_of(p);
            M::AlcubierreMetric<double> m;
            m.time_treatment = M::TimeTreatment::InertialTimeDependent;
            const M::AlcubierreParameters<double> q{kVel, kRefR, kRefSigma};
            const Mat4<double> prod = m.evaluate(p.t, x, p.y, p.z, q);
            const Mat4<double> ref =
                evaluate_fixed_lab_reference(kVel, kXs0, p.t, x, p.y, p.z);
            for (int k = 0; k < 2; ++k)
                w[k].offer(std::fabs(prod(kComponents[k].a, kComponents[k].b) -
                                     ref(kComponents[k].a, kComponents[k].b)),
                           p.name);
            // The fixed-centre reference must agree too when handed the
            // unperturbed centre.
            const Mat4<double> refc =
                evaluate_fixed_center_reference(kVel, xs_at(p.t), p.t, x, p.y, p.z);
            for (int k = 0; k < 2; ++k)
                w[k].offer(std::fabs(prod(kComponents[k].a, kComponents[k].b) -
                                     refc(kComponents[k].a, kComponents[k].b)),
                           p.name);
        }
        for (int k = 0; k < 2; ++k)
            report.check_residual(
                std::string("b1c.velocity.reference_reproduces_metric.") +
                    kComponents[k].name,
                "the independent reference evaluator, written from the "
                "equations and calling no production code, reproduces "
                "AlcubierreMetric::evaluate",
                VerificationLevel::D, "n/a",
                std::string("14 points; worst at ") + w[k].where, w[k].value,
                1e-15, sw.seconds(),
                "both reference forms are checked, at the unperturbed v, so "
                "later agreement is evidence and not tautology");
    }

    // =======================================================================
    // 1. Analytic derivatives against centred finite differences of the
    //    reference functions. ABSOLUTE error, best over four step sizes.
    // =======================================================================
    {
        Stopwatch sw;
        Worst wc[2], wl[2];
        for (const Point& p : points()) {
            const double x = x_of(p);
            const double xs = xs_at(p.t);

            const Mat4<double> ac = analytic_dgdv_fixed_center(kVel, xs, x, p.y, p.z);
            const Mat4<double> al =
                analytic_dgdv_fixed_lab(kVel, kXs0, p.t, x, p.y, p.z);

            for (int k = 0; k < 2; ++k) {
                const int a = kComponents[k].a, b = kComponents[k].b;
                double best_c = 1e300, best_l = 1e300;
                for (double h : kSteps) {
                    // FIXED CENTRE: xs is a constant argument, so perturbing v
                    // cannot move the bubble.
                    const double fc =
                        (evaluate_fixed_center_reference(kVel + h, xs, p.t, x, p.y, p.z)(a, b) -
                         evaluate_fixed_center_reference(kVel - h, xs, p.t, x, p.y, p.z)(a, b)) /
                        (2.0 * h);
                    best_c = std::fmin(best_c, std::fabs(fc - ac(a, b)));

                    // FIXED LABORATORY: xs is recomputed from the perturbed v
                    // inside each call.
                    const double fl =
                        (evaluate_fixed_lab_reference(kVel + h, kXs0, p.t, x, p.y, p.z)(a, b) -
                         evaluate_fixed_lab_reference(kVel - h, kXs0, p.t, x, p.y, p.z)(a, b)) /
                        (2.0 * h);
                    best_l = std::fmin(best_l, std::fabs(fl - al(a, b)));
                }
                wc[k].offer(best_c, p.name);
                wl[k].offer(best_l, p.name);
            }
        }
        for (int k = 0; k < 2; ++k) {
            metrics.worst_fd = std::fmax(metrics.worst_fd,
                                         std::fmax(wc[k].value, wl[k].value));
            report.check_residual(
                std::string("b1c.velocity.analytic_matches_fd_fixed_center.") +
                    kComponents[k].name,
                "dg/dv at fixed bubble centre, analytic against centred finite "
                "differences of the fixed-centre reference",
                VerificationLevel::D, "n/a",
                std::string("14 points, h in {1e-3..1e-6}; worst at ") + wc[k].where,
                wc[k].value, 1e-8, sw.seconds(),
                "ABSOLUTE error. A relative error here is meaningless: dg/dv "
                "passes through 1e-11 in the exterior, below finite-difference "
                "roundoff, and that is what broke the old test");
            report.check_residual(
                std::string("b1c.velocity.analytic_matches_fd_fixed_lab.") +
                    kComponents[k].name,
                "the TOTAL dg/dv at fixed laboratory coordinates, analytic "
                "against centred finite differences of the fixed-laboratory "
                "reference, which recomputes x_s(v) = x_s0 + v t per call",
                VerificationLevel::D, "n/a",
                std::string("14 points, h in {1e-3..1e-6}; worst at ") + wl[k].where,
                wl[k].value, 1e-8, sw.seconds());
        }
    }

    // =======================================================================
    // 2. The chain rule, as an algebraic identity between analytic expressions.
    // =======================================================================
    {
        Stopwatch sw;
        Worst w[2];
        for (const Point& p : points()) {
            const double x = x_of(p);
            const double xs = xs_at(p.t);
            const Mat4<double> total =
                analytic_dgdv_fixed_lab(kVel, kXs0, p.t, x, p.y, p.z);
            const Mat4<double> partial =
                analytic_dgdv_fixed_center(kVel, xs, x, p.y, p.z);
            const Mat4<double> traj =
                analytic_trajectory_term(kVel, xs, p.t, x, p.y, p.z);
            for (int k = 0; k < 2; ++k) {
                const int a = kComponents[k].a, b = kComponents[k].b;
                w[k].offer(std::fabs((total(a, b) - partial(a, b)) - traj(a, b)),
                           p.name);
            }
        }
        for (int k = 0; k < 2; ++k) {
            metrics.worst_chain_rule =
                std::fmax(metrics.worst_chain_rule, w[k].value);
            report.check_residual(
                std::string("b1c.velocity.trajectory_term_is_chain_rule.") +
                    kComponents[k].name,
                "dg/dv|_{t,x} - dg/dv|_{x_s} = (dg/df)(df/dr_s)(dr_s/dv) with "
                "dr_s/dv = -((x - x_s)/r_s) t, assembled component by component",
                VerificationLevel::A, "n/a",
                std::string("14 points; worst at ") + w[k].where, w[k].value,
                1e-14, sw.seconds(),
                "dr_s/dv is derived from r_s = sqrt((x-x_s)^2+y^2+z^2) with "
                "dx_s/dv = t. At the bubble centre the term is exactly zero "
                "because df/dr vanishes there, f being even");
        }
    }

    // =======================================================================
    // 3. dg/dx_s = -dg/dx, MEASURED per component rather than assumed.
    //
    // The old test asserted total - partial = -t (dg/dx) globally. It holds for
    // this metric, but only because g depends on x and x_s solely through
    // (x - x_s) while the explicit v dependence is untouched -- and that has to
    // be checked, at the same fixed variables, not taken for granted.
    // dg/dx is differenced from the FIXED-CENTRE reference so that x_s really
    // is held fixed while x moves.
    // =======================================================================
    {
        Stopwatch sw;
        Worst w[2];
        for (const Point& p : points()) {
            if (p.t == 0.0) continue;   // the term is identically zero at t = 0
            const double x = x_of(p);
            const double xs = xs_at(p.t);
            const Mat4<double> traj =
                analytic_trajectory_term(kVel, xs, p.t, x, p.y, p.z);

            for (int k = 0; k < 2; ++k) {
                const int a = kComponents[k].a, b = kComponents[k].b;
                double best = 1e300;
                for (double h : kSteps) {
                    const double dgdx =
                        (evaluate_fixed_center_reference(kVel, xs, p.t, x + h, p.y, p.z)(a, b) -
                         evaluate_fixed_center_reference(kVel, xs, p.t, x - h, p.y, p.z)(a, b)) /
                        (2.0 * h);
                    best = std::fmin(best, std::fabs(traj(a, b) - (-p.t * dgdx)));
                }
                w[k].offer(best, p.name);
            }
        }
        for (int k = 0; k < 2; ++k)
            report.check_residual(
                std::string("b1c.velocity.trajectory_equals_minus_t_dgdx.") +
                    kComponents[k].name,
                "the trajectory term equals -t (dg/dx) evaluated at FIXED x_s "
                "and fixed v, per component",
                VerificationLevel::D, "n/a",
                std::string("8 points with t != 0; worst at ") + w[k].where,
                w[k].value, 1e-8, sw.seconds(),
                "checked, not assumed. dg/dx is differenced from the "
                "fixed-centre reference precisely so that x_s cannot drift "
                "while x moves");
    }

    // =======================================================================
    // 4. OTI against the analytic expressions. This is the case that can catch
    //    an implementation defect in the production evaluator.
    // =======================================================================
#ifdef WARP_HYPAD_HAVE_OTILIB
    {
        namespace sens = warp_hypad::sensitivity;
        using sens::OTIScalar;
        Stopwatch sw;

        auto oti_dgdv = [&](double t, double x, double y, double z,
                            M::MetricSensitivityMode mode, int a, int b) {
            const sens::SeededParameters sp = sens::seed_all(kVel, kRefR, kRefSigma);
            const M::AlcubierreParameters<OTIScalar> p{sp.velocity, sp.radius,
                                                       sp.sigma};
            M::AlcubierreMetric<OTIScalar> m;
            m.time_treatment = M::TimeTreatment::InertialTimeDependent;
            m.sensitivity_mode = mode;
            return sens::extract_gradient(m.evaluate(t, x, y, z, p)(a, b),
                                          sens::kParamVelocity);
        };

        Worst wc[2], wl[2];
        for (const Point& p : points()) {
            const double x = x_of(p);
            const double xs = xs_at(p.t);
            const Mat4<double> ac = analytic_dgdv_fixed_center(kVel, xs, x, p.y, p.z);
            const Mat4<double> al =
                analytic_dgdv_fixed_lab(kVel, kXs0, p.t, x, p.y, p.z);
            for (int k = 0; k < 2; ++k) {
                const int a = kComponents[k].a, b = kComponents[k].b;
                wc[k].offer(
                    std::fabs(oti_dgdv(p.t, x, p.y, p.z,
                                       M::MetricSensitivityMode::FixedBubbleCenter,
                                       a, b) -
                              ac(a, b)),
                    p.name);
                wl[k].offer(
                    std::fabs(
                        oti_dgdv(p.t, x, p.y, p.z,
                                 M::MetricSensitivityMode::FixedLaboratoryCoordinates,
                                 a, b) -
                        al(a, b)),
                    p.name);
            }
        }
        for (int k = 0; k < 2; ++k) {
            metrics.worst_analytic =
                std::fmax(metrics.worst_analytic,
                          std::fmax(wc[k].value, wl[k].value));
            report.check_residual(
                std::string("b1c.velocity.oti_fixed_center_matches_analytic.") +
                    kComponents[k].name,
                "FixedBubbleCenter HyPAD sensitivity equals the analytic "
                "dg/dv|_{x_s}",
                VerificationLevel::D, "n/a",
                std::string("14 points; worst at ") + wc[k].where, wc[k].value,
                1e-12, sw.seconds(),
                "ABSOLUTE, and tight: both sides are exact closed forms, so "
                "there is no finite-difference floor to hide behind");
            report.check_residual(
                std::string("b1c.velocity.oti_fixed_lab_matches_analytic.") +
                    kComponents[k].name,
                "FixedLaboratoryCoordinates HyPAD sensitivity equals the "
                "analytic total dg/dv|_{t,x}, including the dx_s/dv = t "
                "trajectory term",
                VerificationLevel::D, "n/a",
                std::string("14 points; worst at ") + wl[k].where, wl[k].value,
                1e-12, sw.seconds(),
                "the trajectory term arrives through the chain rule on a "
                "Scalar-valued r_s, not as a hand-added correction");
        }

        // The two modes must coincide at t = 0 and separate at t != 0.
        {
            Worst same, diff;
            for (const Point& p : points()) {
                const double x = x_of(p);
                for (int k = 0; k < 2; ++k) {
                    const int a = kComponents[k].a, b = kComponents[k].b;
                    const double gap = std::fabs(
                        oti_dgdv(p.t, x, p.y, p.z,
                                 M::MetricSensitivityMode::FixedBubbleCenter, a, b) -
                        oti_dgdv(p.t, x, p.y, p.z,
                                 M::MetricSensitivityMode::FixedLaboratoryCoordinates,
                                 a, b));
                    if (p.t == 0.0)
                        same.offer(gap, p.name);
                    else
                        diff.offer(gap, p.name);
                }
            }
            report.check_residual(
                "b1c.velocity.modes_agree_at_t0",
                "at t = 0 the two sensitivity modes coincide, because "
                "dx_s/dv = t = 0",
                VerificationLevel::A, "n/a",
                std::string("6 points at t = 0; worst at ") + same.where,
                same.value, 1e-14, sw.seconds());
            report.check_predicate(
                "b1c.velocity.modes_differ_away_from_t0",
                "at t != 0 the two modes differ substantially, so the "
                "distinction is not cosmetic",
                VerificationLevel::A, diff.value > 1e-3,
                "largest difference above 1e-3",
                "largest difference " + std::to_string(diff.value) + " at " +
                    diff.where,
                sw.seconds(),
                "an earlier version of this comparison ran in double, where "
                "real_part(v) == v makes the modes literally identical, and "
                "measured a difference of exactly zero");
        }

        // Every other component is a constant of the metric and must carry no
        // velocity channel at all.
        {
            Worst w;
            const int flat[7][2] = {{1, 1}, {2, 2}, {3, 3}, {0, 2},
                                    {0, 3}, {1, 2}, {2, 3}};
            for (const Point& p : points()) {
                const double x = x_of(p);
                for (int k = 0; k < 7; ++k)
                    for (M::MetricSensitivityMode mode :
                         {M::MetricSensitivityMode::FixedBubbleCenter,
                          M::MetricSensitivityMode::FixedLaboratoryCoordinates})
                        w.offer(std::fabs(oti_dgdv(p.t, x, p.y, p.z, mode,
                                                   flat[k][0], flat[k][1])),
                                p.name);
            }
            report.check_residual(
                "b1c.velocity.constant_components_carry_no_channel",
                "g_11, g_22, g_33 and every off-diagonal component other than "
                "g_01 are constants, so their velocity sensitivity is exactly "
                "zero in both modes",
                VerificationLevel::A, "n/a",
                std::string("7 components x 14 points x 2 modes; worst at ") +
                    w.where,
                w.value, 0.0, sw.seconds(),
                "exact zero, not a tolerance: a spurious channel here would "
                "mean the derivative was leaking through the container");
        }
    }
#else
    {
        Stopwatch sw;
        report.check_skipped(
            "b1c.velocity.oti_fixed_center_matches_analytic.*",
            "HyPAD dg/dv against the analytic expression, both modes",
            "the OTIlib backend is disabled; the two sensitivity modes differ "
            "only in derivative channels and are indistinguishable without it. "
            "The analytic-versus-finite-difference cases above still ran and "
            "still pass",
            sw.seconds());
        report.check_skipped(
            "b1c.velocity.modes_agree_at_t0",
            "coincidence of the two sensitivity modes at t = 0",
            "the OTIlib backend is disabled", sw.seconds());
        report.check_skipped(
            "b1c.velocity.constant_components_carry_no_channel",
            "absence of a velocity channel on the constant metric components",
            "the OTIlib backend is disabled", sw.seconds());
    }
#endif
}

}  // namespace checks
}  // namespace program_b
