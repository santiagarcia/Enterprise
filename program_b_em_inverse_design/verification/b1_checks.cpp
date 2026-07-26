#include <cmath>
#include <fstream>
#include <sstream>

#include "b1_checks.hpp"
#include "program_b/effective_metric/metric_medium_map.hpp"
#include "warp_hypad/core/test_scalar.hpp"
#include "warp_hypad/metrics/metric.hpp"

#ifdef WARP_HYPAD_HAVE_OTILIB
#include "warp_hypad/sensitivity/otilib_traits.hpp"
#endif

namespace program_b {
namespace checks {

using warp_hypad::Mat4;
using warp_hypad::Stopwatch;
using warp_hypad::VerificationLevel;
using warp_hypad::VerificationReport;

namespace {

Mat4<double> diag4(double a, double b, double c, double d) {
    Mat4<double> m;
    m(0, 0) = a; m(1, 1) = b; m(2, 2) = c; m(3, 3) = d;
    return m;
}

}  // namespace

void run_b1_checks(VerificationReport& report, const std::string& out_dir) {
    // -----------------------------------------------------------------------
    // CASE 1 -- Vacuum Minkowski.
    //
    // Fixes the overall normalisation of the map. In Heaviside-Lorentz units
    // with c = 1, Minkowski must give exactly eps = mu = I and w = 0.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        const Mat4<double> eta = Mat4<double>::minkowski();
        const ForwardResult<double> f = metric_to_medium(eta);

        double worst = f.ok ? 0.0 : 1.0;
        if (f.ok) {
            worst = std::fmax(worst,
                              mat3_max_abs_diff(f.response.epsilon,
                                                Mat3<double>::identity()));
            worst = std::fmax(worst, mat3_max_abs_diff(f.response.mu,
                                                       Mat3<double>::identity()));
            worst = std::fmax(worst, mat3_max_abs(f.response.xi));
            worst = std::fmax(worst, mat3_max_abs(f.response.zeta));
        }

        report.check_residual(
            "b1.minkowski_is_vacuum",
            "eta_{mu nu} maps to eps = mu = I and zero magnetoelectric coupling",
            VerificationLevel::B, "n/a", "Heaviside-Lorentz, c = 1", worst, 1e-15,
            sw.seconds(),
            "fixes the overall normalisation of the Plebanski map in these "
            "units");
    }

    // -----------------------------------------------------------------------
    // CASE 2 -- Static isotropic impedance-matched medium.
    //
    // g = diag(-1, a, a, a) has eps = mu = a^{1/2} I. Check against the closed
    // form and check the implied wave speed.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        double worst = 0.0;
        for (double a : {0.25, 1.0, 4.0, 9.0}) {
            const Mat4<double> g = diag4(-1.0, a, a, a);
            const ForwardResult<double> f = metric_to_medium(g);
            if (!f.ok) { worst = 1.0; continue; }
            // sqrt(-g) = sqrt(a^3), g^ij = I/a, g_00 = -1
            //   eps = sqrt(a^3)/a = sqrt(a)
            const double expect = std::sqrt(a);
            for (int i = 0; i < 3; ++i)
                worst = std::fmax(worst,
                                  std::fabs(f.response.epsilon(i, i) - expect));
        }
        report.check_residual(
            "b1.isotropic_forward",
            "g = diag(-1,a,a,a) gives eps = mu = sqrt(a) I",
            VerificationLevel::B, "n/a", "a = 0.25, 1, 4, 9", worst, 1e-14,
            sw.seconds());
    }

    // Refractive index relation, via the Gordon metric.
    {
        Stopwatch sw;
        double worst = 0.0;
        for (double n : {1.0, 1.5, 2.0, 3.0}) {
            const Mat4<double> gbar = gordon_metric_static_isotropic(n);
            // Null condition g_00 dt^2 + g_11 dx^2 = 0 gives
            // (dx/dt)^2 = -g_00/g_11, NOT -g_11/g_00. For
            // diag(-1/n^2, 1, 1, 1) that is 1/n.
            const double speed = std::sqrt(-gbar(0, 0) / gbar(1, 1));
            worst = std::fmax(worst, std::fabs(speed - 1.0 / n));
        }
        report.check_residual(
            "b1.gordon_optical_metric_speed",
            "the static isotropic Gordon metric diag(-1/n^2,1,1,1) has null "
            "coordinate speed 1/n",
            VerificationLevel::B, "n/a", "n = 1, 1.5, 2, 3", worst, 1e-15,
            sw.seconds(),
            "this is a RAY statement for a medium that need not be impedance "
            "matched, and is a different claim from the Plebanski map");
    }

    // The two notions must be kept apart: a plain dielectric is NOT
    // metric-representable in the Plebanski sense even though it has a Gordon
    // metric. This is the distinction most often conflated.
    {
        Stopwatch sw;
        ConstitutiveResponse<double> plain;
        plain.epsilon = Mat3<double>::diagonal(4.0, 4.0, 4.0);  // n = 2
        plain.mu = Mat3<double>::identity();                    // non-magnetic
        const ResponseClassification c = classify(plain);

        report.check_predicate(
            "b1.plain_dielectric_not_metric_representable",
            "a non-magnetic dielectric (eps = 4I, mu = I) is rejected as "
            "metric-representable, despite having a well-defined Gordon "
            "optical metric",
            VerificationLevel::A, !c.metric_representable && !c.impedance_matched,
            "rejected, impedance mismatch reported",
            c.metric_representable ? "ACCEPTED (wrong)" : c.obstruction,
            sw.seconds(),
            "separates full-Maxwell metric equivalence from ray-level optical "
            "metrics");
    }

    // -----------------------------------------------------------------------
    // CASE 3 -- Diagonal anisotropic medium, manufactured.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        const Mat4<double> g = diag4(-1.0, 2.0, 3.0, 5.0);
        const ForwardResult<double> f = metric_to_medium(g);
        // sqrt(-g) = sqrt(2*3*5) = sqrt(30); g^ij = diag(1/2,1/3,1/5)
        // eps = sqrt(30) * diag(1/2,1/3,1/5)
        const double root = std::sqrt(30.0);
        const double expect[3] = {root / 2.0, root / 3.0, root / 5.0};
        double worst = f.ok ? 0.0 : 1.0;
        for (int i = 0; i < 3 && f.ok; ++i)
            worst = std::fmax(worst, std::fabs(f.response.epsilon(i, i) - expect[i]));

        report.check_residual(
            "b1.anisotropic_forward",
            "g = diag(-1,2,3,5) gives eps^ii = sqrt(30)/g_ii",
            VerificationLevel::B, "n/a", "manufactured diagonal metric", worst,
            1e-14, sw.seconds());
    }

    // -----------------------------------------------------------------------
    // CASE 4 -- Conformal invariance of the forward map.
    //
    // This is the structural fact that forces the inverse to be
    // conformal-class-valued, so it is verified directly rather than assumed.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        const Mat4<double> g = diag4(-1.0, 2.0, 3.0, 5.0);
        const ForwardResult<double> base = metric_to_medium(g);
        double worst = 0.0;
        for (double omega2 : {0.25, 0.5, 2.0, 7.0}) {
            Mat4<double> scaled;
            for (int i = 0; i < 16; ++i) scaled.m[i] = omega2 * g.m[i];
            const ForwardResult<double> f = metric_to_medium(scaled);
            if (!f.ok) { worst = 1.0; continue; }
            worst = std::fmax(worst,
                              mat3_max_abs_diff(f.response.epsilon,
                                                base.response.epsilon));
            for (int i = 0; i < 3; ++i)
                worst = std::fmax(worst, std::fabs(f.w[i] - base.w[i]));
        }
        report.check_residual(
            "b1.conformal_invariance",
            "the medium is unchanged under g -> Omega^2 g, for Omega^2 = 0.25, "
            "0.5, 2, 7",
            VerificationLevel::A, "n/a", "diag(-1,2,3,5)", worst, 1e-13,
            sw.seconds(),
            "sqrt(-g) contributes Omega^4, g^{ij} contributes Omega^-2 and "
            "1/g_00 contributes Omega^-2, so the three cancel exactly. This is "
            "WHY the inverse can only return a conformal class");
    }

    // -----------------------------------------------------------------------
    // CASE 5 -- Round trip g -> medium -> ghat, up to conformal factor.
    //
    // Componentwise equality is the wrong test and would fail; conformal
    // equivalence is what is physically determined.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        double worst_conformal = 0.0;
        double biggest_omega_deviation = 0.0;
        int failures = 0;

        const Mat4<double> cases[4] = {
            Mat4<double>::minkowski(),
            diag4(-1.0, 2.0, 3.0, 5.0),
            diag4(-4.0, 4.0, 4.0, 4.0),
            diag4(-0.25, 1.5, 0.75, 2.25),
        };

        for (const Mat4<double>& g : cases) {
            const ForwardResult<double> f = metric_to_medium(g);
            if (!f.ok) { ++failures; continue; }
            const InverseResultMetric<double> inv = medium_to_metric(f.response);
            if (!inv.ok()) { ++failures; continue; }

            const ConformalComparison<double> cmp =
                compare_up_to_conformal(g, inv.metric);
            worst_conformal = std::fmax(worst_conformal, cmp.residual);
            biggest_omega_deviation =
                std::fmax(biggest_omega_deviation, std::fabs(cmp.omega_squared - 1.0));
        }

        report.check_residual(
            "b1.roundtrip_conformal",
            "g -> (eps,mu) -> ghat recovers the CONFORMAL CLASS of g",
            VerificationLevel::A, "n/a", "4 diagonal metrics",
            worst_conformal + failures, 1e-13, sw.seconds(),
            "componentwise equality is deliberately NOT demanded: the forward "
            "map discards the conformal factor, so requiring it would be "
            "demanding information the physics does not determine");

        // And the conformal factor really is not 1 in general -- if it were,
        // the previous check would be silently testing componentwise equality.
        report.check_predicate(
            "b1.roundtrip_factor_is_nontrivial",
            "at least one round trip returns a representative with Omega^2 != 1, "
            "confirming the conformal test is not a disguised equality test",
            VerificationLevel::A, biggest_omega_deviation > 1e-6,
            "some Omega^2 differs from 1",
            "largest |Omega^2 - 1| = " + std::to_string(biggest_omega_deviation),
            sw.seconds());
    }

    // Null-cone equivalence: the conformal representative must have the same
    // null directions as the original. This is the physically meaningful
    // invariant, and it is checked independently of the conformal-factor test.
    {
        Stopwatch sw;
        const Mat4<double> g = diag4(-1.0, 2.0, 3.0, 5.0);
        const ForwardResult<double> f = metric_to_medium(g);
        const InverseResultMetric<double> inv = medium_to_metric(f.response);

        double worst = inv.ok() ? 0.0 : 1.0;
        if (inv.ok()) {
            // For each spatial direction, the null speed |dx/dt| from
            // g_00 dt^2 + g_ii dx^2 = 0 must agree between g and ghat.
            for (int i = 1; i < 4 && inv.ok(); ++i) {
                const double v_orig = std::sqrt(-g(0, 0) / g(i, i));
                const double v_hat =
                    std::sqrt(-inv.metric(0, 0) / inv.metric(i, i));
                worst = std::fmax(worst, std::fabs(v_orig - v_hat) /
                                             std::fmax(v_orig, 1e-300));
            }
        }
        report.check_residual(
            "b1.null_cone_equivalence",
            "the reconstructed metric has the same null coordinate speeds along "
            "each axis as the original",
            VerificationLevel::A, "n/a", "diag(-1,2,3,5)", worst, 1e-13,
            sw.seconds(),
            "a conformally invariant observable, so this is the right thing to "
            "demand of a reconstruction");
    }

    // -----------------------------------------------------------------------
    // CASE 6 -- Non-metric constitutive tensors must be DETECTED.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        int detected = 0, total = 0;
        std::ostringstream detail;

        // (a) Birefringent: eps != mu, uniaxial.
        {
            ++total;
            ConstitutiveResponse<double> r;
            r.epsilon = Mat3<double>::diagonal(2.0, 2.0, 3.0);
            r.mu = Mat3<double>::identity();
            const InverseResultMetric<double> inv = medium_to_metric(r);
            if (inv.status == InverseStatus::NotMetricRepresentable) ++detected;
            detail << "uniaxial birefringent: " << to_string(inv.status) << "; ";
        }
        // (b) Axion (Tellegen) part: xi = zeta = alpha I, symmetric not skew.
        {
            ++total;
            ConstitutiveResponse<double> r;
            r.epsilon = Mat3<double>::identity();
            r.mu = Mat3<double>::identity();
            r.xi = mat3_scale(Mat3<double>::identity(), 0.3);
            r.zeta = r.xi;
            const InverseResultMetric<double> inv = medium_to_metric(r);
            if (inv.status == InverseStatus::NotMetricRepresentable) ++detected;
            detail << "axion/Tellegen: " << to_string(inv.status) << "; ";
        }
        // (c) Skewon-like: eps not symmetric.
        {
            ++total;
            ConstitutiveResponse<double> r;
            r.epsilon = Mat3<double>::identity();
            r.epsilon(0, 1) = 0.4;   // asymmetric
            r.mu = r.epsilon;
            const InverseResultMetric<double> inv = medium_to_metric(r);
            if (inv.status == InverseStatus::NotMetricRepresentable) ++detected;
            detail << "asymmetric eps: " << to_string(inv.status) << "; ";
        }

        report.check_predicate(
            "b1.nonmetric_detected",
            "constitutive tensors that no single pseudo-Riemannian metric can "
            "produce -- birefringent, axion/Tellegen, and skewon-like -- are "
            "each REPORTED rather than fitted",
            VerificationLevel::A, detected == total,
            std::to_string(total) + " of " + std::to_string(total) + " detected",
            std::to_string(detected) + " detected. " + detail.str(), sw.seconds(),
            "the requirement is that the code refuse, not that it produce a "
            "best-fit metric with no warning");
    }

    // A valid medium must still be ACCEPTED, or the detector above would be
    // passing by rejecting everything.
    {
        Stopwatch sw;
        ConstitutiveResponse<double> r;
        r.epsilon = Mat3<double>::diagonal(2.0, 3.0, 4.0);
        r.mu = r.epsilon;
        const InverseResultMetric<double> inv = medium_to_metric(r);
        report.check_predicate(
            "b1.valid_medium_accepted",
            "an impedance-matched symmetric medium with no magnetoelectric "
            "coupling is accepted and reconstructed",
            VerificationLevel::A, inv.ok(), "accepted",
            std::string(to_string(inv.status)), sw.seconds(),
            "guards against a classifier that rejects everything");
    }

    // -----------------------------------------------------------------------
    // CASE 7 -- The warp metric as a medium, and where it stops existing.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        using warp_hypad::metrics::AlcubierreMetric;
        using warp_hypad::metrics::AlcubierreParameters;

        const AlcubierreParameters<double> p{0.2, 3.0, 2.0};
        AlcubierreMetric<double> warp;

        // Subluminal: v f <= 0.2 < 1, so g_00 < 0 everywhere and the medium
        // exists at every point.
        int ok_count = 0, total = 0;
        for (double x = -8; x <= 8; x += 0.25) {
            ++total;
            const Mat4<double> g = warp.evaluate(0.0, x, 0.0, 0.0, p);
            if (metric_to_medium(g).ok) ++ok_count;
        }
        report.check_predicate(
            "b1.warp_subluminal_medium_exists",
            "for v_s = 0.2 the Alcubierre metric has an equivalent medium at "
            "every sampled point",
            VerificationLevel::B, ok_count == total,
            std::to_string(total) + " of " + std::to_string(total),
            std::to_string(ok_count) + " of " + std::to_string(total),
            sw.seconds());

        // The magnetoelectric coupling is nonzero exactly where the shift is.
        const Mat4<double> gw = warp.evaluate(0.0, 3.0, 0.0, 0.0, p);
        const ForwardResult<double> fw = metric_to_medium(gw);
        report.check_predicate(
            "b1.warp_is_bianisotropic",
            "the warp medium has a nonzero magnetoelectric vector w, i.e. it is "
            "bianisotropic, not a simple dielectric",
            VerificationLevel::B, fw.ok && std::fabs(fw.w[0]) > 1e-6,
            "w_x != 0 at the bubble wall",
            fw.ok ? ("w_x = " + std::to_string(fw.w[0])) : "forward map failed",
            sw.seconds(),
            "a realizability fact: no arrangement of ordinary dielectrics can "
            "reproduce this without magnetoelectric coupling");

        // And its inverse is correctly reported as out of B1 scope rather than
        // silently approximated.
        const InverseResultMetric<double> iw = medium_to_metric(fw.response);
        report.check_predicate(
            "b1.warp_inverse_reports_unsupported",
            "reconstructing a metric from a magnetoelectric medium is reported "
            "as unsupported in B1 rather than approximated",
            VerificationLevel::A,
            iw.status == InverseStatus::MagnetoelectricUnsupported,
            "MagnetoelectricUnsupported", std::string(to_string(iw.status)),
            sw.seconds());
    }

    // SUPERLUMINAL: the map must fail where v_s f > 1, and it must fail by
    // reporting, not by returning nonsense.
    {
        Stopwatch sw;
        using warp_hypad::metrics::AlcubierreMetric;
        using warp_hypad::metrics::AlcubierreParameters;

        const AlcubierreParameters<double> p{2.0, 3.0, 2.0};  // v_s = 2 > 1
        AlcubierreMetric<double> warp;

        int failed = 0, total = 0;
        double first_fail_x = 0.0;
        bool found = false;
        for (double x = -8; x <= 8; x += 0.05) {
            ++total;
            const Mat4<double> g = warp.evaluate(0.0, x, 0.0, 0.0, p);
            const ForwardResult<double> f = metric_to_medium(g);
            if (!f.ok) {
                ++failed;
                if (!found) { first_fail_x = x; found = true; }
            }
        }

        std::ostringstream got;
        got << failed << " of " << total << " sample points have no equivalent "
            << "medium; first at x = " << first_fail_x;

        report.check_predicate(
            "b1.superluminal_medium_fails_reported",
            "for v_s = 2 there are points where g_00 >= 0 and NO equivalent "
            "medium exists, and the map reports this instead of returning a "
            "value",
            VerificationLevel::B, failed > 0, "some points rejected", got.str(),
            sw.seconds(),
            "this is the central realizability limitation of the whole "
            "programme: the superluminal regime is exactly the regime in which "
            "the metric-to-medium correspondence breaks down in this frame. It "
            "is the same v_s f > 1 condition that costs Program A its "
            "coordinate observers");
    }

    // -----------------------------------------------------------------------
    // Scalar genericity.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        using warp_hypad::TestScalar;
        Mat4<TestScalar> g;
        g(0, 0) = TestScalar(-1.0);
        g(1, 1) = TestScalar(2.0);
        g(2, 2) = TestScalar(3.0);
        g(3, 3) = TestScalar(5.0);
        const ForwardResult<TestScalar> ft = metric_to_medium(g);

        const Mat4<double> gd = diag4(-1.0, 2.0, 3.0, 5.0);
        const ForwardResult<double> fd = metric_to_medium(gd);

        double worst = (ft.ok && fd.ok) ? 0.0 : 1.0;
        for (int i = 0; i < 9 && ft.ok && fd.ok; ++i)
            worst = std::fmax(worst,
                              std::fabs(warp_hypad::real_part(ft.response.epsilon.m[i]) -
                                        fd.response.epsilon.m[i]));

        report.check_residual(
            "b1.generic.test_scalar",
            "the constitutive mapping instantiates with a custom scalar having "
            "no ordering, no equality and no implicit double conversion",
            VerificationLevel::A, "n/a", "diag(-1,2,3,5)", worst, 0.0,
            sw.seconds());
    }

    // --- Export for the figures --------------------------------------------
    {
        using warp_hypad::metrics::AlcubierreMetric;
        using warp_hypad::metrics::AlcubierreParameters;
        const AlcubierreParameters<double> p{0.2, 3.0, 2.0};
        AlcubierreMetric<double> warp;

        std::ofstream os(out_dir + "/warp_medium_profile.csv", std::ios::binary);
        if (os) {
            os << "x,eps_xx,eps_yy,eps_zz,w_x,medium_exists\n";
            os.setf(std::ios::scientific);
            os.precision(10);
            for (int i = -320; i <= 320; ++i) {
                const double x = 0.025 * i;
                const Mat4<double> g = warp.evaluate(0.0, x, 0.0, 0.0, p);
                const ForwardResult<double> f = metric_to_medium(g);
                os << x << ',';
                if (f.ok) {
                    os << f.response.epsilon(0, 0) << ','
                       << f.response.epsilon(1, 1) << ','
                       << f.response.epsilon(2, 2) << ',' << f.w[0] << ",1\n";
                } else {
                    os << "nan,nan,nan,nan,0\n";
                }
            }
        }
    }
    {
        // Where the correspondence exists, as a function of v_s.
        using warp_hypad::metrics::AlcubierreMetric;
        using warp_hypad::metrics::AlcubierreParameters;
        AlcubierreMetric<double> warp;
        std::ofstream os(out_dir + "/realizability_map.csv", std::ios::binary);
        if (os) {
            os << "v_s,x,medium_exists\n";
            for (int vi = 1; vi <= 30; ++vi) {
                const double vs = 0.1 * vi;
                const AlcubierreParameters<double> p{vs, 3.0, 2.0};
                for (int i = -160; i <= 160; ++i) {
                    const double x = 0.05 * i;
                    const Mat4<double> g = warp.evaluate(0.0, x, 0.0, 0.0, p);
                    os << vs << ',' << x << ','
                       << (metric_to_medium(g).ok ? 1 : 0) << '\n';
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OTI sensitivities through the constitutive mapping.
// ---------------------------------------------------------------------------

#ifndef WARP_HYPAD_HAVE_OTILIB

void run_b1_sensitivity_checks(VerificationReport& report,
                               const std::string& out_dir) {
    (void)out_dir;
    Stopwatch sw;
    report.check_skipped(
        "b1.sensitivity.*",
        "HyPAD derivatives of the constitutive response with respect to metric "
        "parameters, the refractive index, and the conformal normalisation",
        "the OTIlib backend is disabled; configure with "
        "-DWARP_HYPAD_WITH_OTILIB=ON to run these",
        sw.seconds());
}

#else

void run_b1_sensitivity_checks(VerificationReport& report,
                               const std::string& out_dir) {
    using warp_hypad::sensitivity::extract_gradient;
    using warp_hypad::sensitivity::OTIScalar;
    namespace sens = warp_hypad::sensitivity;

    // (a) METRIC PARAMETERS: d eps / d(v,R,sigma) through the warp metric.
    {
        Stopwatch sw;
        using warp_hypad::metrics::AlcubierreMetric;
        using warp_hypad::metrics::AlcubierreParameters;

        const double v0 = 0.2, R0 = 3.0, s0 = 2.0;
        const double x = 3.0;

        const sens::SeededParameters sp = sens::seed_all(v0, R0, s0);
        const AlcubierreParameters<OTIScalar> p{sp.velocity, sp.radius, sp.sigma};
        AlcubierreMetric<OTIScalar> warp;
        const ForwardResult<OTIScalar> f =
            metric_to_medium(warp.evaluate(0.0, x, 0.0, 0.0, p));

        // Reference: centred differences in each parameter, all in double,
        // swept over steps and taking the best.
        auto eps_yy_double = [&](double v, double rr, double ss) {
            AlcubierreMetric<double> w;
            const AlcubierreParameters<double> q{v, rr, ss};
            return metric_to_medium(w.evaluate(0.0, x, 0.0, 0.0, q))
                .response.epsilon(1, 1);
        };

        const double base[3] = {v0, R0, s0};
        const double steps[4] = {1e-3, 1e-4, 1e-5, 1e-6};
        double worst = f.ok ? 0.0 : 1.0;

        for (int d = 0; d < 3 && f.ok; ++d) {
            const double oti = extract_gradient(f.response.epsilon(1, 1), d);
            double best = 1e300;
            for (double h : steps) {
                double lo[3] = {base[0], base[1], base[2]};
                double hi[3] = {base[0], base[1], base[2]};
                lo[d] -= h; hi[d] += h;
                const double fd = (eps_yy_double(hi[0], hi[1], hi[2]) -
                                   eps_yy_double(lo[0], lo[1], lo[2])) / (2 * h);
                best = std::fmin(best, std::fabs(fd - oti) /
                                           std::fmax(std::fabs(oti), 1e-9));
            }
            worst = std::fmax(worst, best);
        }

        report.check_residual(
            "b1.sensitivity.metric_parameters",
            "d eps_yy / d(v, R, sigma) through the Plebanski map matches "
            "centred finite differences",
            VerificationLevel::D, "n/a", "v=0.2, R=3, sigma=2 at the wall x=3",
            worst, 1e-6, sw.seconds(),
            "the constitutive mapping is templated, so OTI channels flow from "
            "the metric parameters through to the medium with no extra code");
    }

    // (b) CONSTITUTIVE PARAMETER: d(reconstructed metric)/dn for an isotropic
    //     impedance-matched medium of index n, against the closed form.
    //
    //     eps = mu = n^2 I  ->  h = det(eps) eps^{-1} = n^6 / n^2 I = n^4 I
    //     so d h_xx / dn = 4 n^3, exactly.
    {
        Stopwatch sw;
        const double n0 = 1.7;
        const OTIScalar n = sens::seed_parameter(n0, sens::kParamVelocity);

        ConstitutiveResponse<OTIScalar> r;
        const OTIScalar n2 = n * n;
        r.epsilon = Mat3<OTIScalar>::diagonal(n2, n2, n2);
        r.mu = r.epsilon;

        const InverseResultMetric<OTIScalar> inv = medium_to_metric(r);
        double worst = inv.ok() ? 0.0 : 1.0;
        if (inv.ok()) {
            const double got =
                extract_gradient(inv.metric(1, 1), sens::kParamVelocity);
            const double expect = 4.0 * n0 * n0 * n0;
            worst = std::fabs(got - expect) / std::fabs(expect);
            // and the value itself
            worst = std::fmax(worst,
                              std::fabs(warp_hypad::real_part(inv.metric(1, 1)) -
                                        std::pow(n0, 4.0)) / std::pow(n0, 4.0));
        }

        report.check_residual(
            "b1.sensitivity.refractive_index",
            "for eps = mu = n^2 I the reconstruction gives h = n^4 I and "
            "dh/dn = 4n^3, both matched by HyPAD",
            VerificationLevel::B, "n/a", "n = 1.7", worst, 1e-13, sw.seconds(),
            "an analytic reference, not a finite-difference one");
    }

    // (c) CONFORMAL NORMALISATION parameter: the forward map must be
    //     INSENSITIVE to it. d eps / d Omega == 0 is the differential statement
    //     of conformal invariance, and is a sharper test than comparing values.
    {
        Stopwatch sw;
        const double omega0 = 1.4;
        const OTIScalar omega = sens::seed_parameter(omega0, sens::kParamRadius);
        const OTIScalar omega2 = omega * omega;

        Mat4<OTIScalar> g;
        const double base[4] = {-1.0, 2.0, 3.0, 5.0};
        for (int i = 0; i < 4; ++i)
            g(i, i) = omega2 * warp_hypad::scalar_cast<OTIScalar>(base[i]);

        const ForwardResult<OTIScalar> f = metric_to_medium(g);
        double worst = f.ok ? 0.0 : 1.0;
        for (int i = 0; i < 9 && f.ok; ++i)
            worst = std::fmax(worst,
                              std::fabs(extract_gradient(f.response.epsilon.m[i],
                                                         sens::kParamRadius)));

        report.check_residual(
            "b1.sensitivity.conformal_normalisation",
            "d eps / d Omega == 0: the differential form of conformal "
            "invariance of the forward map",
            VerificationLevel::A, "n/a", "Omega = 1.4 on diag(-1,2,3,5)", worst,
            1e-13, sw.seconds(),
            "sharper than comparing values at discrete Omega: it asserts the "
            "derivative vanishes identically, which a near-invariance would "
            "fail");
    }

    // (d) A DERIVED OBSERVABLE: the null coordinate speed along y for an
    //     isotropic impedance-matched medium is 1/n^2, so d(speed)/dn = -2/n^3.
    {
        Stopwatch sw;
        const double n0 = 1.7;
        const OTIScalar n = sens::seed_parameter(n0, sens::kParamSigma);

        ConstitutiveResponse<OTIScalar> r;
        const OTIScalar n2 = n * n;
        r.epsilon = Mat3<OTIScalar>::diagonal(n2, n2, n2);
        r.mu = r.epsilon;

        const InverseResultMetric<OTIScalar> inv = medium_to_metric(r);
        double worst = inv.ok() ? 0.0 : 1.0;
        if (inv.ok()) {
            using std::sqrt;
            const OTIScalar speed =
                sqrt(-inv.metric(0, 0) / inv.metric(1, 1));
            const double got = extract_gradient(speed, sens::kParamSigma);
            const double expect = -2.0 / (n0 * n0 * n0);
            worst = std::fabs(got - expect) / std::fabs(expect);
            worst = std::fmax(
                worst, std::fabs(warp_hypad::real_part(speed) - 1.0 / (n0 * n0)) *
                           (n0 * n0));
        }

        report.check_residual(
            "b1.sensitivity.observable_null_speed",
            "the null coordinate speed 1/n^2 of the reconstructed medium, and "
            "its derivative -2/n^3, both match HyPAD",
            VerificationLevel::B, "n/a", "n = 1.7", worst, 1e-12, sw.seconds(),
            "an observable derived from the effective metric, not a raw "
            "component");
    }

    // --- Export sensitivity profile ----------------------------------------
    {
        using warp_hypad::metrics::AlcubierreMetric;
        using warp_hypad::metrics::AlcubierreParameters;
        const sens::SeededParameters sp = sens::seed_all(0.2, 3.0, 2.0);
        const AlcubierreParameters<OTIScalar> p{sp.velocity, sp.radius, sp.sigma};
        AlcubierreMetric<OTIScalar> warp;

        std::ofstream os(out_dir + "/medium_sensitivity.csv", std::ios::binary);
        if (os) {
            os << "x,deps_yy_dv,deps_yy_dR,deps_yy_dsigma,dw_dv,dw_dR,dw_dsigma\n";
            os.setf(std::ios::scientific);
            os.precision(10);
            for (int i = -320; i <= 320; ++i) {
                const double x = 0.025 * i;
                const ForwardResult<OTIScalar> f =
                    metric_to_medium(warp.evaluate(0.0, x, 0.0, 0.0, p));
                if (!f.ok) continue;
                os << x;
                for (int d = 0; d < 3; ++d)
                    os << ',' << extract_gradient(f.response.epsilon(1, 1), d);
                for (int d = 0; d < 3; ++d)
                    os << ',' << extract_gradient(f.w[0], d);
                os << '\n';
            }
        }
    }
}

#endif  // WARP_HYPAD_HAVE_OTILIB

}  // namespace checks
}  // namespace program_b
