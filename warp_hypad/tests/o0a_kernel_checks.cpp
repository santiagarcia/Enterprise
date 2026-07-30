#include "o0a_checks.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>

#include "warp_hypad/design/quintic_bspline.hpp"
#include "warp_hypad/design/radial_validation.hpp"

// ---------------------------------------------------------------------------
// O0A cases 1-10: the univariate quintic B-spline kernel.
//
// The continuity cases are the interesting ones, and they need a measure that
// distinguishes "continuous with a kink" from "discontinuous". Comparing a
// derivative at t-d and t+d is not enough: for a CONTINUOUS derivative with a
// bounded one-sided slope, that difference is O(d) and shrinks, while for a
// DISCONTINUOUS one it tends to a nonzero constant. Both look "small" at a
// single d. So each continuity case measures the jump at two step sizes and
// tests the SCALING, which separates the two cases decisively:
//
//     simple interior knot, order 4:  relative jump 1.1e-3 -> 1.1e-4 -> 1.1e-5
//     double interior knot, order 4:  relative jump 1.57 -> 1.57 -> 1.57
//
// The first is C4 and the second is C3, which is exactly p - m for m = 1 and 2.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using design::BasisEvaluation1D;
using design::QuinticBSplineBasis1D;

// The fourth derivative from a one-sided difference of the ANALYTIC third. Exact
// up to rounding for a quintic, whose fourth derivative is piecewise linear.
double d4_left(const QuinticBSplineBasis1D& b, int i, double t, double d) {
    return (b.derivative_of(i, t - d, 3) - b.derivative_of(i, t - 2 * d, 3)) / d;
}
double d4_right(const QuinticBSplineBasis1D& b, int i, double t, double d) {
    return (b.derivative_of(i, t + 2 * d, 3) - b.derivative_of(i, t + d, 3)) / d;
}

// The worst relative jump in the fourth derivative across t.
double relative_jump4(const QuinticBSplineBasis1D& b, double t, double d) {
    double jump = 0.0, scale = 0.0;
    for (int i = 0; i < b.size(); ++i) {
        const double l = d4_left(b, i, t, d), r = d4_right(b, i, t, d);
        jump = std::fmax(jump, std::fabs(r - l));
        scale = std::fmax(scale, std::fmax(std::fabs(l), std::fabs(r)));
    }
    return (scale > 0.0) ? jump / scale : 0.0;
}

// The worst jump in derivatives 0..order across t, absolute.
double jump_upto(const QuinticBSplineBasis1D& b, double t, int order, double d) {
    double j = 0.0;
    for (int i = 0; i < b.size(); ++i)
        for (int k = 0; k <= order; ++k)
            j = std::fmax(j, std::fabs(b.derivative_of(i, t + d, k) -
                                       b.derivative_of(i, t - d, k)));
    return j;
}

// A basis with one deliberately doubled interior knot.
QuinticBSplineBasis1D repeated_knot_basis() {
    std::vector<double> kv;
    for (int i = 0; i <= 5; ++i) kv.push_back(0.0);
    kv.push_back(0.25);
    kv.push_back(0.5);
    kv.push_back(0.5);
    kv.push_back(0.75);
    for (int i = 0; i <= 5; ++i) kv.push_back(1.0);
    return QuinticBSplineBasis1D(std::move(kv));
}

}  // namespace

void run_o0a_kernel_checks(VerificationReport& report,
                           const std::string& out_dir) {
    Stopwatch sw;
    const design::RadialValidationConfig f0cfg;
    QuinticBSplineBasis1D basis(0.0, f0cfg.s_max, f0cfg.n_basis);
    const std::string params = "degree 5, 12 basis functions, clamped, s in [0,4]";

    // --- 1. A valid degree-five knot vector ---------------------------------
    sw.reset();
    const std::size_t expect_knots =
        static_cast<std::size_t>(f0cfg.n_basis + design::kQuinticDegree + 1);
    bool nondecreasing = true;
    for (std::size_t i = 1; i < basis.knots().size(); ++i)
        if (basis.knots()[i] < basis.knots()[i - 1]) nondecreasing = false;
    int end_multiplicity_ok = 1;
    for (int k = 0; k <= design::kQuinticDegree; ++k) {
        if (basis.knots()[static_cast<std::size_t>(k)] != basis.lower())
            end_multiplicity_ok = 0;
        if (basis.knots()[basis.knots().size() - 1 -
                          static_cast<std::size_t>(k)] != basis.upper())
            end_multiplicity_ok = 0;
    }
    report.check_predicate(
        "o0a.valid_quintic_knot_vector",
        "n + p + 1 knots, non-decreasing, clamped with multiplicity p+1 = 6 at "
        "each end, and n - p - 1 = 6 simple interior knots",
        VerificationLevel::A,
        basis.valid() && basis.knots().size() == expect_knots &&
            nondecreasing && end_multiplicity_ok == 1 &&
            basis.distinct_interior_knots().size() == 6,
        "18 knots, clamped, 6 simple interior",
        std::to_string(basis.knots().size()) + " knots, " +
            std::to_string(basis.distinct_interior_knots().size()) +
            " distinct interior, clamped=" +
            std::to_string(end_multiplicity_ok),
        sw.seconds(),
        "the knot vector IS the basis: an off-by-one here changes the function "
        "space silently and every downstream residual with it");

    // --- 2. The minimum basis count is enforced -----------------------------
    sw.reset();
    QuinticBSplineBasis1D too_small(0.0, 1.0, 4);
    QuinticBSplineBasis1D exactly_min(0.0, 1.0, design::kMinQuinticBasisCount);
    report.check_predicate(
        "o0a.minimum_basis_count_enforced",
        "a degree-5 family with fewer than 6 functions is REJECTED, and 6 is "
        "accepted",
        VerificationLevel::A,
        !too_small.valid() && exactly_min.valid() &&
            !too_small.rejection_reason().empty(),
        "n = 4 rejected with a reason, n = 6 accepted",
        std::string("n=4 valid=") + (too_small.valid() ? "yes" : "no") +
            " (\"" + too_small.rejection_reason() + "\"), n=6 valid=" +
            (exactly_min.valid() ? "yes" : "no"),
        sw.seconds(),
        "the 6 x 4 lattice proposed earlier is inadmissible in the second "
        "direction: 4 < 6 leaves no interior knot at all, so the family is a "
        "single Bernstein patch with no local support");

    // --- 3. Partition of unity ----------------------------------------------
    sw.reset();
    double worst_pu = 0.0;
    for (int k = 0; k <= 2000; ++k) {
        const double x = f0cfg.s_max * static_cast<double>(k) / 2000.0;
        const BasisEvaluation1D e = basis.evaluate(x);
        double s = 0.0;
        for (std::size_t i = 0; i < e.size(); ++i) s += e.value[i];
        worst_pu = std::fmax(worst_pu, std::fabs(s - 1.0));
    }
    report.check_residual(
        "o0a.partition_of_unity", "sum_i N_i(x) = 1 on the clamped interval",
        VerificationLevel::A, "2001 samples", params, worst_pu, 1.0e-14,
        sw.seconds(),
        "checked on the CLAMPED validation basis. The correction lattices are "
        "deliberately not partitions of unity: they taper to zero at the edge of "
        "the free band, which is what makes the collar exact");

    // --- 4. Local support ----------------------------------------------------
    sw.reset();
    std::size_t max_active = 0, min_active = 99;
    for (int k = 0; k <= 2000; ++k) {
        const double x = f0cfg.s_max * static_cast<double>(k) / 2000.0;
        const std::size_t n = basis.evaluate(x).size();
        if (n > max_active) max_active = n;
        if (n < min_active) min_active = n;
    }
    report.check_predicate(
        "o0a.local_support",
        "at most p+1 = 6 functions are active at any point, and exactly 6 "
        "inside the interval",
        VerificationLevel::A, max_active == 6 && min_active == 6,
        "6 active everywhere",
        "min " + std::to_string(min_active) + ", max " +
            std::to_string(max_active),
        sw.seconds(),
        "compact support is what makes the sparse stencil sparse and the collar "
        "exclusion possible");

    // --- 5. Exactly zero outside support ------------------------------------
    sw.reset();
    double leak = 0.0;
    std::size_t outside_samples = 0;
    for (int i = 0; i < basis.size(); ++i) {
        const double lo = basis.support_lower(i), hi = basis.support_upper(i);
        for (int k = 0; k <= 2000; ++k) {
            const double x = f0cfg.s_max * static_cast<double>(k) / 2000.0;
            if (x >= lo - 1e-12 && x <= hi + 1e-12) continue;
            ++outside_samples;
            for (int ord = 0; ord <= 3; ++ord)
                leak = std::fmax(leak, std::fabs(basis.derivative_of(i, x, ord)));
        }
    }
    report.check_exact(
        "o0a.exact_zero_outside_support",
        "N_i and all its derivatives are EXACTLY zero outside "
        "[knot_i, knot_{i+6}]",
        VerificationLevel::A, std::to_string(outside_samples) + " samples", 0.0,
        leak, sw.seconds(),
        "exact rather than small: the flat-collar guarantee of section 6 is a "
        "statement about supports, and it only holds if a function outside its "
        "support contributes bit-zero rather than a rounding-level value");

    // --- 6-8. Derivative orders against an independent finite difference ----
    //
    // The analytic recurrence is the production path; centred differences are the
    // INDEPENDENT check, and the thing being tested is the convergence ORDER. A
    // single step size cannot distinguish a correct derivative from a wrong one
    // with a small coefficient, so each order is measured at three step sizes and
    // the observed rate is reported.
    const char* names[3] = {"o0a.first_derivative_analytic",
                            "o0a.second_derivative_analytic",
                            "o0a.third_derivative_analytic"};
    const char* refs[3] = {
        "N_i' from the recurrence, against a centred difference of N_i",
        "N_i'' from the recurrence, against a centred second difference",
        "N_i''' from the recurrence, against a centred third difference"};
    // ONE LADDER PER ORDER, because the rounding floor of a centred
    // difference of order n grows like eps/h^n while its truncation error falls
    // like h^2. The two cross at h ~ eps^(1/(n+2)): around 1e-5 for a first
    // difference but around 1e-3 for a third. A single ladder ending at h = 1e-4
    // put the third difference below its own noise floor -- the measured error
    // ROSE from 9.5e-4 to 2.9e-3 and the observed rate collapsed to 0.3, which is
    // the signature of rounding, not of a wrong derivative.
    double err_at[3][3];
    const double ladder[3][3] = {{1e-2, 1e-3, 1e-4},   // first derivative
                                 {1e-2, 1e-3, 1e-4},   // second derivative
                                 {4e-2, 2e-2, 1e-2}};  // third derivative
    for (int si = 0; si < 3; ++si) {
        double e1 = 0.0, e2 = 0.0, e3 = 0.0;
        const std::vector<double> interior = basis.distinct_interior_knots();
        // Stay clear of the ends and of every knot line by more than the largest
        // stencil reach, so no difference straddles a discontinuity in the fifth
        // derivative.
        const double clear = 0.14;
        for (int k = 0; k <= 400; ++k) {
            const double x = f0cfg.s_max * static_cast<double>(k) / 400.0;
            if (x < clear || x > f0cfg.s_max - clear) continue;
            bool near = false;
            for (double t : interior)
                if (std::fabs(x - t) < clear) near = true;
            if (near) continue;
            for (int i = 0; i < basis.size(); ++i) {
                const double h1 = ladder[0][si], h2 = ladder[1][si],
                             h3 = ladder[2][si];
                e1 = std::fmax(
                    e1, std::fabs((basis.value_of(i, x + h1) -
                                   basis.value_of(i, x - h1)) / (2 * h1) -
                                  basis.derivative_of(i, x, 1)));
                e2 = std::fmax(
                    e2, std::fabs((basis.value_of(i, x + h2) -
                                   2 * basis.value_of(i, x) +
                                   basis.value_of(i, x - h2)) / (h2 * h2) -
                                  basis.derivative_of(i, x, 2)));
                e3 = std::fmax(
                    e3, std::fabs((basis.value_of(i, x + 2 * h3) -
                                   2 * basis.value_of(i, x + h3) +
                                   2 * basis.value_of(i, x - h3) -
                                   basis.value_of(i, x - 2 * h3)) /
                                      (2 * h3 * h3 * h3) -
                                  basis.derivative_of(i, x, 3)));
            }
        }
        err_at[0][si] = e1;
        err_at[1][si] = e2;
        err_at[2][si] = e3;
    }
    for (int ord = 0; ord < 3; ++ord) {
        sw.reset();
        const double r1 = err_at[ord][0] / err_at[ord][1];
        const double r2 = err_at[ord][1] / err_at[ord][2];
        char buf[320];
        std::snprintf(buf, sizeof(buf),
                      "errors %.3e, %.3e, %.3e at h = %.0e, %.0e, %.0e; "
                      "observed rates %.1f and %.1f (second order is 100 for the "
                      "first two ladders and 4 for the third, whose steps fall by "
                      "2x to stay above the eps/h^3 rounding floor)",
                      err_at[ord][0], err_at[ord][1], err_at[ord][2],
                      ladder[ord][0], ladder[ord][1], ladder[ord][2], r1, r2);
        // The expected ratio is (h_coarse/h_fine)^2: 100 for the decade
        // ladders, 4 for the halving ladder the third derivative needs.
        const double expect =
            (ladder[ord][0] / ladder[ord][1]) * (ladder[ord][0] / ladder[ord][1]);
        report.check_predicate(
            names[ord], refs[ord], VerificationLevel::B,
            r1 > 0.5 * expect && r1 < 2.0 * expect && r2 > 0.5 * expect &&
                r2 < 2.0 * expect,
            "second-order convergence: ratio near " + std::to_string(expect) +
                " per step",
            buf, sw.seconds(),
            "finite differences are used ONLY as this independent path; "
            "production derivatives never come from a difference");
    }

    // --- 9. C4 continuity at an ordinary interior knot ----------------------
    sw.reset();
    const std::vector<double> interior = basis.distinct_interior_knots();
    const double t_simple = interior[2];
    const double j03_a = jump_upto(basis, t_simple, 3, 1e-5);
    const double j03_b = jump_upto(basis, t_simple, 3, 1e-6);
    const double rel4_a = relative_jump4(basis, t_simple, 1e-5);
    const double rel4_b = relative_jump4(basis, t_simple, 1e-6);
    {
        char buf[400];
        std::snprintf(buf, sizeof(buf),
                      "orders 0-3 jump %.2e -> %.2e as d falls 10x; relative "
                      "order-4 jump %.2e -> %.2e (falls with d, so the 4th "
                      "derivative is continuous too)",
                      j03_a, j03_b, rel4_a, rel4_b);
        report.check_predicate(
            "o0a.c4_continuity_ordinary_knot",
            "at a multiplicity-1 interior knot the spline is C4: derivatives 0 "
            "through 4 are continuous",
            VerificationLevel::A,
            j03_b < 0.2 * j03_a && rel4_b < 0.2 * rel4_a && rel4_b < 1e-3,
            "every jump falls linearly with the probe offset",
            buf, sw.seconds(),
            "degree 5 is required for exactly this: the metric is "
            "differentiated twice for curvature and three times for the "
            "Bianchi diagnostic, and a C2 second derivative would put knot "
            "lines into the curvature");
    }

    // --- 10. A repeated knot reduces continuity correctly -------------------
    sw.reset();
    const QuinticBSplineBasis1D rb = repeated_knot_basis();
    const double dbl_a = relative_jump4(rb, 0.5, 1e-5);
    const double dbl_b = relative_jump4(rb, 0.5, 1e-6);
    const double smp_a = relative_jump4(rb, 0.25, 1e-5);
    const double smp_b = relative_jump4(rb, 0.25, 1e-6);
    double pu_rep = 0.0;
    for (int k = 0; k <= 1000; ++k) {
        const double x = static_cast<double>(k) / 1000.0;
        const BasisEvaluation1D e = rb.evaluate(x);
        double s = 0.0;
        for (std::size_t i = 0; i < e.size(); ++i) s += e.value[i];
        pu_rep = std::fmax(pu_rep, std::fabs(s - 1.0));
    }
    {
        char buf[440];
        std::snprintf(buf, sizeof(buf),
                      "multiplicity 2: relative order-4 jump %.3f then %.3f "
                      "(constant, so discontinuous); multiplicity 1 on the same "
                      "basis: %.2e then %.2e (falls, so continuous). Partition "
                      "of unity still %.2e",
                      dbl_a, dbl_b, smp_a, smp_b, pu_rep);
        report.check_predicate(
            "o0a.repeated_knot_continuity_reduced",
            "a multiplicity-m interior knot leaves C^(p-m): at m = 2 the fourth "
            "derivative jumps, while a simple knot on the same basis does not",
            VerificationLevel::A,
            rb.valid() && rb.interior_knot_multiplicity(0.5) == 2 &&
                dbl_b > 0.5 * dbl_a && dbl_b > 0.1 && smp_b < 0.2 * smp_a &&
                pu_rep < 1e-14,
            "the multiple knot's jump does not shrink; the simple knot's does",
            buf, sw.seconds(),
            "the production lattices use simple knots only; this case exists so "
            "that a future family which needs reduced continuity gets what it "
            "asks for rather than silently keeping C4");
    }

    // --- Records -------------------------------------------------------------
    {
        std::ofstream f(out_dir + "/radial_squared_knot_vector.csv");
        f << "index,knot,is_clamped_end\n";
        for (std::size_t i = 0; i < basis.knots().size(); ++i) {
            const bool clamped =
                (i <= static_cast<std::size_t>(design::kQuinticDegree)) ||
                (i + static_cast<std::size_t>(design::kQuinticDegree) + 1 >=
                 basis.knots().size());
            f << i << "," << basis.knots()[i] << "," << (clamped ? 1 : 0)
              << "\n";
        }
    }
    {
        std::ofstream f(out_dir + "/univariate_basis_profiles.csv");
        f << "s,index,value\n";
        std::ofstream g(out_dir + "/univariate_basis_derivatives.csv");
        g << "s,index,first,second,third\n";
        for (int k = 0; k <= 800; ++k) {
            const double x = f0cfg.s_max * static_cast<double>(k) / 800.0;
            for (int i = 0; i < basis.size(); ++i) {
                f << x << "," << i << "," << basis.value_of(i, x) << "\n";
                g << x << "," << i << "," << basis.derivative_of(i, x, 1) << ","
                  << basis.derivative_of(i, x, 2) << ","
                  << basis.derivative_of(i, x, 3) << "\n";
            }
        }
    }
}

}  // namespace checks
}  // namespace warp_hypad
