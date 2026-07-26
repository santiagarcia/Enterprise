#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "metric_checks.hpp"
#include "warp_hypad/core/test_scalar.hpp"
#include "warp_hypad/geometry/signature.hpp"
#include "warp_hypad/metrics/metric_field.hpp"

namespace warp_hypad {
namespace checks {
namespace {

using metrics::AlcubierreMetric;
using metrics::AlcubierreParameters;
using metrics::MinkowskiMetric;

constexpr double kV = 0.2;
constexpr double kR = 3.0;
constexpr double kSigma = 2.0;

Grid make_grid(int n, double half) {
    GridSpec s;
    s.nx = s.ny = s.nz = n;
    s.x_min = -half; s.x_max = half;
    s.y_min = -half; s.y_max = half;
    s.z_min = -half; s.z_max = half;
    return Grid(s);
}

AlcubierreParameters<double> params(double v = kV, double R = kR,
                                    double s = kSigma) {
    return AlcubierreParameters<double>(v, R, s);
}

}  // namespace

void run_metric_checks(VerificationReport& report, const std::string& out_dir) {
    const AlcubierreMetric<double> warp;
    const MinkowskiMetric<double> flat;

    // -----------------------------------------------------------------------
    // Minkowski: exact values, exact determinant, self-inverse.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        const Grid g = make_grid(21, 5.0);
        double worst_component = 0.0, worst_det = 0.0, worst_sym = 0.0;
        for (std::size_t p = 0; p < g.num_points(); ++p) {
            int i, j, k;
            g.unflatten(p, i, j, k);
            const Mat4<double> m =
                flat.evaluate(0.0, g.x(i), g.y(j), g.z(k), {});
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b) {
                    const double want =
                        (a != b) ? 0.0 : (a == 0 ? -1.0 : 1.0);
                    worst_component =
                        std::fmax(worst_component, std::fabs(m(a, b) - want));
                }
            worst_det = std::fmax(worst_det, std::fabs(det(m) + 1.0));
            worst_sym = std::fmax(worst_sym, symmetry_residual(m));
        }
        report.check_residual(
            "metric.minkowski.exact_values",
            "g_{mu nu} == diag(-1,1,1,1) at every grid point",
            VerificationLevel::A, "21x21x21", "n/a", worst_component, 0.0,
            sw.seconds(), "exact equality: every entry is 0 or +-1");
        report.check_residual(
            "metric.minkowski.determinant",
            "det g == -1", VerificationLevel::A, "21x21x21", "n/a", worst_det,
            0.0, sw.seconds());
        report.check_residual(
            "metric.minkowski.symmetry", "g_{mu nu} == g_{nu mu}",
            VerificationLevel::A, "21x21x21", "n/a", worst_sym, 0.0, sw.seconds());
    }

    // -----------------------------------------------------------------------
    // Shape function: exact interior value, decay, monotonicity, range.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        const AlcubierreParameters<double> p = params();

        // f(0) = [tanh(sR) - tanh(-sR)] / (2 tanh(sR)) = 1, algebraically exact.
        const double f0 = metrics::shape_function(0.0, p);
        report.check_value("metric.shape.interior_value",
                           "f(0) == 1 exactly, by construction of the numerator",
                           VerificationLevel::B, "n/a", "R=3, sigma=2", 1.0, f0,
                           1e-15, sw.seconds());

        const double f_far = metrics::shape_function(60.0, p);
        report.check_residual(
            "metric.shape.exterior_decay",
            "f(r) -> 0 as r -> infinity (asymptotic flatness of the shape "
            "function)",
            VerificationLevel::B, "n/a", "r = 60, R = 3, sigma = 2", f_far,
            1e-12, sw.seconds());

        // Monotone decreasing in r, and bounded in [0,1].
        double worst_increase = 0.0, worst_range = 0.0;
        double prev = f0;
        for (int i = 1; i <= 2000; ++i) {
            const double r = 0.02 * i;
            const double f = metrics::shape_function(r, p);
            worst_increase = std::fmax(worst_increase, f - prev);
            worst_range = std::fmax(worst_range,
                                    std::fmax(f - 1.0, -f));
            prev = f;
        }
        report.check_residual(
            "metric.shape.monotone_decreasing",
            "f is non-increasing in r on [0, 40]", VerificationLevel::A, "n/a",
            "2000 samples", worst_increase, 1e-15, sw.seconds());
        report.check_residual(
            "metric.shape.bounded",
            "0 <= f <= 1 everywhere", VerificationLevel::A, "n/a",
            "2000 samples", worst_range, 1e-15, sw.seconds());
    }

    // -----------------------------------------------------------------------
    // Parameter effects behave as the physics requires.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        // Wall location: the radius where f crosses 1/2 should track R.
        auto wall_radius = [](double R, double s) {
            const AlcubierreParameters<double> p = params(kV, R, s);
            double lo = 0.0, hi = 10.0 * R + 20.0;
            for (int it = 0; it < 200; ++it) {
                const double mid = 0.5 * (lo + hi);
                if (metrics::shape_function(mid, p) > 0.5) lo = mid; else hi = mid;
            }
            return 0.5 * (lo + hi);
        };
        const double w2 = wall_radius(2.0, kSigma);
        const double w4 = wall_radius(4.0, kSigma);
        const double w6 = wall_radius(6.0, kSigma);
        report.check_predicate(
            "metric.shape.radius_shifts_wall",
            "increasing R moves the f = 1/2 wall outward",
            VerificationLevel::A, (w4 > w2) && (w6 > w4),
            "wall radius strictly increasing with R",
            "R=2 -> " + std::to_string(w2) + ", R=4 -> " + std::to_string(w4) +
                ", R=6 -> " + std::to_string(w6),
            sw.seconds());

        // Steepness: |df/dr| at the wall should grow with sigma.
        auto wall_slope = [&](double s) {
            const AlcubierreParameters<double> p = params(kV, kR, s);
            const double r = wall_radius(kR, s);
            const double h = 1e-6;
            return std::fabs((metrics::shape_function(r + h, p) -
                              metrics::shape_function(r - h, p)) /
                             (2.0 * h));
        };
        const double s1 = wall_slope(1.0);
        const double s2 = wall_slope(2.0);
        const double s4 = wall_slope(4.0);
        report.check_predicate(
            "metric.shape.sigma_steepens_wall",
            "increasing sigma increases |df/dr| at the wall",
            VerificationLevel::A, (s2 > s1) && (s4 > s2),
            "wall slope strictly increasing with sigma",
            "sigma=1 -> " + std::to_string(s1) + ", 2 -> " + std::to_string(s2) +
                ", 4 -> " + std::to_string(s4),
            sw.seconds());
    }

    // Velocity must not touch the shape function or the spatial block.
    {
        Stopwatch sw;
        double worst = 0.0;
        for (double r = 0.0; r <= 10.0; r += 0.25)
            worst = std::fmax(worst,
                              std::fabs(metrics::shape_function(r, params(0.1)) -
                                        metrics::shape_function(r, params(0.9))));
        report.check_residual(
            "metric.shape.velocity_independent",
            "f depends on R and sigma only, never on v",
            VerificationLevel::A, "n/a", "v = 0.1 vs 0.9", worst, 0.0,
            sw.seconds(), "exact: v does not appear in the shape function");
    }
    {
        Stopwatch sw;
        double worst = 0.0;
        for (double x = -6.0; x <= 6.0; x += 0.5) {
            const Mat4<double> a = warp.evaluate(0.0, x, 0.3, -0.2, params(0.1));
            const Mat4<double> b = warp.evaluate(0.0, x, 0.3, -0.2, params(0.9));
            for (int i = 1; i < 4; ++i)
                for (int j = 1; j < 4; ++j)
                    worst = std::fmax(worst, std::fabs(a(i, j) - b(i, j)));
        }
        report.check_residual(
            "metric.alcubierre.velocity_affects_expected_components",
            "changing v leaves the purely spatial block g_{ij} unchanged; only "
            "g_00 and g_0x respond",
            VerificationLevel::A, "n/a", "v = 0.1 vs 0.9", worst, 0.0,
            sw.seconds());
    }

    // -----------------------------------------------------------------------
    // Alcubierre metric structure: determinant, symmetry, inverse, signature.
    // -----------------------------------------------------------------------
    {
        const Grid g = make_grid(41, 8.0);
        Stopwatch sw;
        const Rank2Field<double> gf =
            metrics::build_metric_field<double>(g, warp, params());

        const metrics::ValidityReport vr =
            metrics::check_metric_validity(g, gf, -1.0);

        report.check_residual(
            "metric.alcubierre.determinant_exact",
            "det g == -1 identically, for every parameter value and every point "
            "(the (t,x) block has determinant (-1 + v^2 f^2) - (v f)^2 = -1)",
            VerificationLevel::B, "41x41x41", "v=0.2, R=3, sigma=2",
            vr.worst_det_error, 1e-15, sw.seconds(),
            "parameter-independent, so this is a strong pointwise test rather "
            "than a tolerance one");

        report.check_residual(
            "metric.alcubierre.symmetry", "g_{mu nu} == g_{nu mu} everywhere",
            VerificationLevel::A, "41x41x41", "v=0.2, R=3, sigma=2",
            vr.worst_symmetry_residual, 0.0, sw.seconds());

        report.check_residual(
            "metric.alcubierre.inverse_identity",
            "g_{mu alpha} g^{alpha nu} == delta_mu^nu at every point",
            VerificationLevel::A, "41x41x41", "v=0.2, R=3, sigma=2",
            vr.worst_inverse_residual, 1e-12, sw.seconds());

        report.check_predicate(
            "metric.alcubierre.lorentzian_signature",
            "every point has exactly one negative and three positive metric "
            "eigenvalues",
            VerificationLevel::A, vr.num_bad_signature == 0,
            "0 points with wrong signature",
            std::to_string(vr.num_bad_signature) + " of " +
                std::to_string(g.num_points()) + " points",
            sw.seconds(), "computed by Jacobi eigenvalue rotation, not by "
                          "Sylvester's criterion");

        report.check_predicate(
            "metric.alcubierre.all_points_valid",
            "no point is non-finite, asymmetric, singular, ill-conditioned or "
            "of wrong signature",
            VerificationLevel::A, vr.all_valid(), "0 invalid points",
            std::to_string(vr.num_invalid) + " invalid", sw.seconds());
    }

    // Numerical inverse must agree with the closed-form inverse.
    {
        Stopwatch sw;
        const Grid g = make_grid(21, 8.0);
        double worst = 0.0;
        for (std::size_t p = 0; p < g.num_points(); ++p) {
            int i, j, k;
            g.unflatten(p, i, j, k);
            const double x = g.x(i), y = g.y(j), z = g.z(k);
            const Mat4<double> m = warp.evaluate(0.0, x, y, z, params());
            const InverseResult<double> num = inverse(m);
            const Mat4<double> ana = warp.analytic_inverse(0.0, x, y, z, params());
            if (!num.ok) { worst = 1.0; break; }
            worst = std::fmax(worst, max_abs_diff(num.inverse, ana));
        }
        report.check_residual(
            "metric.alcubierre.inverse_vs_analytic",
            "the numerical adjugate inverse agrees with the closed-form inverse "
            "g^00=-1, g^01=-vf, g^11=1-v^2f^2",
            VerificationLevel::B, "21x21x21", "v=0.2, R=3, sigma=2", worst,
            1e-14, sw.seconds(),
            "the pipeline inverts numerically so one code path serves every "
            "metric; the closed form is a reference, not a shortcut");
    }

    // Asymptotic flatness: far from the bubble the metric must become Minkowski.
    {
        Stopwatch sw;
        const Mat4<double> eta = Mat4<double>::minkowski();
        double worst = 0.0;
        const double far = 60.0;
        const double pts[4][3] = {{far, 0, 0}, {0, far, 0}, {0, 0, far},
                                  {far, far, far}};
        for (const auto& q : pts) {
            const Mat4<double> m = warp.evaluate(0.0, q[0], q[1], q[2], params());
            worst = std::fmax(worst, max_abs_diff(m, eta));
        }
        report.check_residual(
            "metric.alcubierre.asymptotic_flatness",
            "g -> eta far from the bubble",
            VerificationLevel::B, "n/a", "|x| = 60, R = 3, sigma = 2", worst,
            1e-12, sw.seconds());
    }

    // At v = 0 the metric must be exactly Minkowski, for any R and sigma.
    {
        Stopwatch sw;
        const Mat4<double> eta = Mat4<double>::minkowski();
        double worst = 0.0;
        for (double x = -5; x <= 5; x += 1.0)
            for (double y = -2; y <= 2; y += 1.0)
                worst = std::fmax(
                    worst, max_abs_diff(warp.evaluate(0.0, x, y, 0.5, params(0.0)),
                                        eta));
        report.check_residual(
            "metric.alcubierre.zero_velocity_is_minkowski",
            "v = 0 gives exactly the Minkowski metric regardless of R, sigma",
            VerificationLevel::A, "n/a", "v = 0", worst, 0.0, sw.seconds());
    }

    // -----------------------------------------------------------------------
    // Signature analyser is itself tested, against matrices of known signature.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        Mat4<double> euclid = Mat4<double>::identity();
        Mat4<double> lorentz = Mat4<double>::minkowski();
        Mat4<double> degenerate = Mat4<double>::identity();
        degenerate(2, 2) = 0.0;
        Mat4<double> two_negative = Mat4<double>::minkowski();
        two_negative(1, 1) = -1.0;

        const auto se = geometry::analyse_signature(euclid);
        const auto sl = geometry::analyse_signature(lorentz);
        const auto sd = geometry::analyse_signature(degenerate);
        const auto s2 = geometry::analyse_signature(two_negative);

        const bool ok = !se.lorentzian && se.positive == 4 &&
                        sl.lorentzian && sl.negative == 1 && sl.positive == 3 &&
                        !sd.lorentzian && sd.zero == 1 &&
                        !s2.lorentzian && s2.negative == 2;

        report.check_predicate(
            "metric.signature_analyser_works",
            "the eigenvalue signature classifier gives (4+,0-) for the "
            "identity, (3+,1-) for Minkowski, detects a zero eigenvalue, and "
            "detects a second negative eigenvalue",
            VerificationLevel::A, ok, "all four classified correctly",
            ok ? "all four correct"
               : ("euclid+" + std::to_string(se.positive) + " lorentz" +
                  std::to_string(sl.lorentzian) + " degen0" +
                  std::to_string(sd.zero) + " two-" + std::to_string(s2.negative)),
            sw.seconds(),
            "a classifier that always said 'Lorentzian' would make the metric "
            "signature gate vacuous");
    }

    // Signature eigenvalues of a rotated Minkowski metric are unchanged --
    // eigenvalues are invariant under similarity, so this checks the solver
    // rather than the metric.
    {
        Stopwatch sw;
        const double c = std::cos(0.7), s = std::sin(0.7);
        Mat4<double> rot = Mat4<double>::identity();
        rot(1, 1) = c;  rot(1, 2) = -s;
        rot(2, 1) = s;  rot(2, 2) = c;
        const Mat4<double> eta = Mat4<double>::minkowski();
        const Mat4<double> rotated = mat_mul(mat_mul(transpose(rot), eta), rot);
        const auto sig = geometry::analyse_signature(rotated);
        double worst = 0.0;
        const double expect[4] = {-1, 1, 1, 1};
        for (int i = 0; i < 4; ++i)
            worst = std::fmax(worst, std::fabs(sig.eigenvalues[i] - expect[i]));
        report.check_residual(
            "metric.signature_rotation_invariant",
            "eigenvalues of R^T eta R equal those of eta",
            VerificationLevel::A, "n/a", "rotation by 0.7 rad in the x-y plane",
            worst, 1e-14, sw.seconds());
    }

    // -----------------------------------------------------------------------
    // Scalar genericity for the metric itself.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        const AlcubierreMetric<TestScalar> warp_ts;
        const AlcubierreMetric<long double> warp_ld;
        // Braces, not parentheses: with parentheses this is a most-vexing-parse
        // function declaration, since TestScalar(kV) reads as a parameter.
        const AlcubierreParameters<TestScalar> pts{TestScalar(kV), TestScalar(kR),
                                                   TestScalar(kSigma)};
        const AlcubierreParameters<long double> pld{kV, kR, kSigma};

        double worst_ts = 0.0, worst_ld = 0.0;
        for (double x = -6; x <= 6; x += 0.75) {
            const Mat4<double> a = warp.evaluate(0.0, x, 0.4, -0.3, params());
            const Mat4<TestScalar> b = warp_ts.evaluate(0.0, x, 0.4, -0.3, pts);
            const Mat4<long double> c2 = warp_ld.evaluate(0.0, x, 0.4, -0.3, pld);
            for (int i = 0; i < 16; ++i) {
                worst_ts = std::fmax(worst_ts,
                                     std::fabs(a.m[i] - real_part(b.m[i])));
                worst_ld = std::fmax(worst_ld,
                                     std::fabs(a.m[i] - static_cast<double>(c2.m[i])));
            }
        }
        report.check_residual(
            "metric.generic.test_scalar_bit_exact",
            "the metric evaluates bit-identically with double and a custom "
            "scalar type",
            VerificationLevel::A, "n/a", "v=0.2, R=3, sigma=2", worst_ts, 0.0,
            sw.seconds());
        report.check_residual(
            "metric.generic.long_double_agrees",
            "long double agrees with double to double precision",
            VerificationLevel::D, "n/a", "v=0.2, R=3, sigma=2", worst_ld, 1e-14,
            sw.seconds());
    }

    // -----------------------------------------------------------------------
    // Frozen-time contract: t must not change anything, and must be recorded.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        double worst = 0.0;
        for (double t : {0.0, 1.0, -5.0, 1000.0}) {
            const Mat4<double> a = warp.evaluate(0.0, 1.3, 0.4, -0.2, params());
            const Mat4<double> b = warp.evaluate(t, 1.3, 0.4, -0.2, params());
            worst = std::fmax(worst, max_abs_diff(a, b));
        }
        report.check_residual(
            "metric.frozen_time_contract",
            "in the frozen-time formulation the metric is independent of t, "
            "exactly",
            VerificationLevel::A, "n/a", "t = 0, 1, -5, 1000", worst, 0.0,
            sw.seconds(),
            "this documents the approximation rather than hiding it: a moving "
            "bubble WOULD depend on t, and this result must change when "
            "TimeTreatment::MovingBubble is implemented");
    }

    // --- Exports for the figures -------------------------------------------
    {
        std::ofstream os(out_dir + "/shape_profiles.csv", std::ios::binary);
        if (os) {
            os << "R,sigma,r,f\n";
            os.setf(std::ios::scientific);
            os.precision(10);
            const double Rs[3] = {2.0, 3.0, 5.0};
            const double Ss[3] = {1.0, 2.0, 5.0};
            for (double R : Rs)
                for (double S : Ss)
                    for (int i = 0; i <= 400; ++i) {
                        const double r = 0.03 * i;
                        os << R << ',' << S << ',' << r << ','
                           << metrics::shape_function(r, params(kV, R, S)) << '\n';
                    }
        }
    }
    {
        const Grid g = make_grid(81, 8.0);
        const int kc = g.center_index_z();
        std::ofstream os(out_dir + "/metric_slice_xy.csv", std::ios::binary);
        if (os && kc >= 0) {
            os << "x,y,g00,g01,g11,det,inverse_residual,signature_ok\n";
            os.setf(std::ios::scientific);
            os.precision(10);
            for (int i = 0; i < g.nx(); ++i)
                for (int j = 0; j < g.ny(); ++j) {
                    const double x = g.x(i), y = g.y(j), z = g.z(kc);
                    const Mat4<double> m = warp.evaluate(0.0, x, y, z, params());
                    const InverseResult<double> inv = inverse(m);
                    const double res = inv.ok ? identity_residual(m, inv.inverse) : -1.0;
                    const auto sig = geometry::analyse_signature(m);
                    os << x << ',' << y << ',' << m(0, 0) << ',' << m(0, 1) << ','
                       << m(1, 1) << ',' << det(m) << ',' << res << ','
                       << (sig.lorentzian ? 1 : 0) << '\n';
                }
        }
    }
    {
        // Centreline profile through the bubble along x.
        std::ofstream os(out_dir + "/centerline_x.csv", std::ios::binary);
        if (os) {
            os << "x,f,g00,g01\n";
            os.setf(std::ios::scientific);
            os.precision(10);
            for (int i = -400; i <= 400; ++i) {
                const double x = 0.02 * i;
                const Mat4<double> m = warp.evaluate(0.0, x, 0.0, 0.0, params());
                os << x << ',' << metrics::shape_function(std::fabs(x), params())
                   << ',' << m(0, 0) << ',' << m(0, 1) << '\n';
            }
        }
    }
}

}  // namespace checks
}  // namespace warp_hypad
