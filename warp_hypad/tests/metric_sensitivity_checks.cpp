#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "metric_checks.hpp"
#include "warp_hypad/core/grid.hpp"
#include "warp_hypad/metrics/metric.hpp"

#ifdef WARP_HYPAD_HAVE_OTILIB
#include "warp_hypad/sensitivity/otilib_traits.hpp"
#endif

// ---------------------------------------------------------------------------
// A3 parameter-sensitivity verification.
//
// This is where the whole design pays off or does not: the metric kernel is
// evaluated with OTIScalar instead of double, with no change to the metric
// code, and the derivatives with respect to v, R and sigma must match closed
// forms derived independently by hand.
//
// Analytic references (g_00 = -1 + v^2 f^2, g_0x = -v f):
//
//   dg00/dv = 2 v f^2          dg0x/dv = -f
//   dg00/dR = 2 v^2 f f_R      dg0x/dR = -v f_R
//   dg00/ds = 2 v^2 f f_s      dg0x/ds = -v f_s
//
// and every purely spatial component is the constant 1, so all of its
// sensitivities are exactly zero -- a check that catches derivative channels
// leaking into components that do not depend on the parameters at all.
//
// f_R and f_s come from metrics::shape_function_derivatives, which is written
// out by hand in src/metrics/metric.cpp from the definition of f, not obtained
// by differentiating the templated implementation. That is what makes this an
// independent comparison rather than a tautology.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {

#ifndef WARP_HYPAD_HAVE_OTILIB

void run_metric_sensitivity_checks(VerificationReport& report,
                                   const std::string& out_dir) {
    (void)out_dir;
    Stopwatch sw;
    report.check_skipped(
        "metric.sensitivity.*",
        "HyPAD parameter sensitivities of the metric, verified against analytic "
        "derivatives and finite differences",
        "the OTIlib backend is disabled; configure with "
        "-DWARP_HYPAD_WITH_OTILIB=ON to run these",
        sw.seconds());
}

#else

namespace {

using metrics::AlcubierreMetric;
using metrics::AlcubierreParameters;
using sensitivity::extract_gradient;
using sensitivity::kParamRadius;
using sensitivity::kParamSigma;
using sensitivity::kParamVelocity;
using sensitivity::OTIScalar;

constexpr double kV = 0.2;
constexpr double kR = 3.0;
constexpr double kSigma = 2.0;

// Sample points spanning the interior, the wall and the exterior. The wall is
// where every derivative is largest and where an error would matter most.
struct Sample { double x, y, z; const char* where; };
// Labels must not contain commas: they are written into a CSV column.
const Sample kSamples[6] = {
    {0.0, 0.0, 0.0, "bubble_centre"},
    {1.5, 0.0, 0.0, "interior"},
    {3.0, 0.0, 0.0, "wall_r_eq_R"},
    {3.0, 1.0, 0.5, "wall_off_axis"},
    {6.0, 0.0, 0.0, "exterior"},
    {2.0, 2.0, 2.0, "off_axis_r3.46"},
};

AlcubierreParameters<OTIScalar> seeded() {
    const sensitivity::SeededParameters p = sensitivity::seed_all(kV, kR, kSigma);
    return AlcubierreParameters<OTIScalar>{p.velocity, p.radius, p.sigma};
}

// Analytic dg_{mu nu}/dp at a point, from the hand-derived shape derivatives.
struct AnalyticSensitivity {
    double dg00[3];
    double dg01[3];
};

AnalyticSensitivity analytic_at(double x, double y, double z) {
    const double r = std::sqrt(x * x + y * y + z * z);
    const metrics::ShapeDerivatives d =
        metrics::shape_function_derivatives(r, kR, kSigma);

    AnalyticSensitivity a{};
    a.dg00[kParamVelocity] = 2.0 * kV * d.f * d.f;
    a.dg00[kParamRadius] = 2.0 * kV * kV * d.f * d.df_dR;
    a.dg00[kParamSigma] = 2.0 * kV * kV * d.f * d.df_dsigma;

    a.dg01[kParamVelocity] = -d.f;
    a.dg01[kParamRadius] = -kV * d.df_dR;
    a.dg01[kParamSigma] = -kV * d.df_dsigma;
    return a;
}

const char* param_name(int i) { return sensitivity::parameter_name(i); }

}  // namespace

void run_metric_sensitivity_checks(VerificationReport& report,
                                   const std::string& out_dir) {
    const AlcubierreMetric<OTIScalar> warp;
    const AlcubierreParameters<OTIScalar> p = seeded();

    std::ofstream csv(out_dir + "/hypad_vs_analytic.csv", std::ios::binary);
    if (csv) {
        csv << "location,x,y,z,component,parameter,hypad,analytic,"
               "absolute_error,relative_error\n";
        csv.setf(std::ios::scientific);
        csv.precision(12);
    }

    // -----------------------------------------------------------------------
    // Shape function sensitivities against the hand-derived closed form.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        double worst = 0.0;
        for (const Sample& s : kSamples) {
            const double r = std::sqrt(s.x * s.x + s.y * s.y + s.z * s.z);
            const OTIScalar f = metrics::shape_function(r, p);
            const metrics::ShapeDerivatives d =
                metrics::shape_function_derivatives(r, kR, kSigma);

            worst = std::fmax(worst, std::fabs(f.real() - d.f));
            worst = std::fmax(worst,
                              std::fabs(extract_gradient(f, kParamRadius) - d.df_dR));
            worst = std::fmax(
                worst, std::fabs(extract_gradient(f, kParamSigma) - d.df_dsigma));
            // f does not depend on v at all.
            worst = std::fmax(worst,
                              std::fabs(extract_gradient(f, kParamVelocity)));
        }
        report.check_residual(
            "metric.sensitivity.shape_function",
            "df/dR and df/dsigma match the hand-derived closed form, and df/dv "
            "is exactly zero, at 6 sample points",
            VerificationLevel::B, "n/a", "v=0.2, R=3, sigma=2", worst, 1e-13,
            sw.seconds(),
            "the closed form is derived independently in metric.cpp, not by "
            "differentiating the templated implementation");
    }

    // -----------------------------------------------------------------------
    // Every nontrivial metric component, every parameter.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        double worst_g00 = 0.0, worst_g01 = 0.0, worst_spatial = 0.0;

        for (const Sample& s : kSamples) {
            const Mat4<OTIScalar> g = warp.evaluate(0.0, s.x, s.y, s.z, p);
            const AnalyticSensitivity a = analytic_at(s.x, s.y, s.z);

            for (int d = 0; d < 3; ++d) {
                const double h00 = extract_gradient(g(0, 0), d);
                const double h01 = extract_gradient(g(0, 1), d);
                const double h10 = extract_gradient(g(1, 0), d);

                worst_g00 = std::fmax(worst_g00, std::fabs(h00 - a.dg00[d]));
                worst_g01 = std::fmax(worst_g01, std::fabs(h01 - a.dg01[d]));
                // The off-diagonal is stored twice; both must agree.
                worst_g01 = std::fmax(worst_g01, std::fabs(h10 - a.dg01[d]));

                // Purely spatial components are the constant 1.
                for (int i = 1; i < 4; ++i)
                    for (int j = 1; j < 4; ++j)
                        worst_spatial =
                            std::fmax(worst_spatial,
                                      std::fabs(extract_gradient(g(i, j), d)));

                if (csv) {
                    // Same floored denominator as the FD comparison; see the
                    // reasoning there.
                    constexpr double kDerivFloor = 1e-9;
                    const double r00 = std::fabs(h00 - a.dg00[d]) /
                                       std::fmax(std::fabs(a.dg00[d]), kDerivFloor);
                    const double r01 = std::fabs(h01 - a.dg01[d]) /
                                       std::fmax(std::fabs(a.dg01[d]), kDerivFloor);
                    csv << s.where << ',' << s.x << ',' << s.y << ',' << s.z
                        << ",g00," << param_name(d) << ',' << h00 << ','
                        << a.dg00[d] << ',' << std::fabs(h00 - a.dg00[d]) << ','
                        << r00 << '\n';
                    csv << s.where << ',' << s.x << ',' << s.y << ',' << s.z
                        << ",g0x," << param_name(d) << ',' << h01 << ','
                        << a.dg01[d] << ',' << std::fabs(h01 - a.dg01[d]) << ','
                        << r01 << '\n';
                }
            }
        }

        report.check_residual(
            "metric.sensitivity.g00_vs_analytic",
            "dg_00/dp == 2 v f^2, 2 v^2 f f_R, 2 v^2 f f_sigma",
            VerificationLevel::B, "n/a", "6 points x 3 parameters", worst_g00,
            1e-13, sw.seconds());

        report.check_residual(
            "metric.sensitivity.g0x_vs_analytic",
            "dg_0x/dp == -f, -v f_R, -v f_sigma, and the stored transpose "
            "agrees",
            VerificationLevel::B, "n/a", "6 points x 3 parameters", worst_g01,
            1e-13, sw.seconds());

        report.check_residual(
            "metric.sensitivity.spatial_block_insensitive",
            "the purely spatial block g_{ij} is the constant 1, so every "
            "sensitivity is exactly zero",
            VerificationLevel::A, "n/a", "6 points x 3 parameters x 9 "
            "components", worst_spatial, 0.0, sw.seconds(),
            "catches derivative channels leaking into parameter-independent "
            "components");
    }

    // -----------------------------------------------------------------------
    // Simultaneous versus separate seeding, on the metric itself.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        double worst = 0.0;
        for (const Sample& s : kSamples) {
            const Mat4<OTIScalar> all = warp.evaluate(0.0, s.x, s.y, s.z, p);
            for (int d = 0; d < 3; ++d) {
                AlcubierreParameters<OTIScalar> one{
                    OTIScalar(kV), OTIScalar(kR), OTIScalar(kSigma)};
                if (d == kParamVelocity)
                    one.velocity = sensitivity::seed_parameter(kV, d);
                else if (d == kParamRadius)
                    one.radius = sensitivity::seed_parameter(kR, d);
                else
                    one.sigma = sensitivity::seed_parameter(kSigma, d);

                const Mat4<OTIScalar> single = warp.evaluate(0.0, s.x, s.y, s.z, one);
                for (int a = 0; a < 4; ++a)
                    for (int b = 0; b < 4; ++b)
                        worst = std::fmax(
                            worst, std::fabs(extract_gradient(single(a, b), d) -
                                             extract_gradient(all(a, b), d)));
            }
        }
        report.check_residual(
            "metric.sensitivity.simultaneous_vs_separate",
            "seeding all three parameters at once gives the same metric "
            "sensitivities as three single-seed evaluations",
            VerificationLevel::A, "n/a", "6 points x 16 components", worst, 0.0,
            sw.seconds());
    }

    // Zero-seed equivalence against the plain double metric.
    {
        Stopwatch sw;
        const AlcubierreMetric<double> warp_d;
        const AlcubierreParameters<double> pd{kV, kR, kSigma};
        const AlcubierreParameters<OTIScalar> p0{
            OTIScalar(kV), OTIScalar(kR), OTIScalar(kSigma)};
        double worst = 0.0;
        for (const Sample& s : kSamples) {
            const Mat4<double> a = warp_d.evaluate(0.0, s.x, s.y, s.z, pd);
            const Mat4<OTIScalar> b = warp.evaluate(0.0, s.x, s.y, s.z, p0);
            for (int i = 0; i < 16; ++i)
                worst = std::fmax(worst, std::fabs(a.m[i] - real_part(b.m[i])));
        }
        // NOT bit-exact, and that is correct behaviour rather than a defect.
        //
        // OTIlib does not implement division as a primitive: onumm3n2_div_oo
        // computes pow(den, -1) and multiplies. The shape function contains the
        // division (A - B) / D, so its real part follows a different rounding
        // path from a hardware divide -- same value to within a couple of ulp,
        // different last bit.
        //
        // The equivalent check in the OTIlib suite (oti.zero_seed_matches_double)
        // IS bit-exact because its test function a^2 sin b + exp(ac) contains no
        // division. Measured discrepancy here: 2.8e-17, i.e. ~0.12 ulp at unit
        // magnitude. Asserting exactness would be asserting something untrue.
        report.check_residual(
            "metric.sensitivity.zero_seed_matches_double",
            "with zero seeds the OTI metric real parts equal the double metric "
            "to within a few ulp",
            VerificationLevel::A, "n/a", "6 points", worst, 1e-15, sw.seconds(),
            "not bit-exact by design: OTI division is pow(-1) composed with "
            "multiply, so the real part rounds differently from a hardware "
            "divide. Exactness does hold for division-free expressions");
    }

    // -----------------------------------------------------------------------
    // HyPAD versus centered finite differences, multiple step sizes.
    // -----------------------------------------------------------------------
    {
        Stopwatch sw;
        const AlcubierreMetric<double> warp_d;
        const double steps[6] = {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8};
        const double base[3] = {kV, kR, kSigma};

        std::ofstream fdcsv(out_dir + "/hypad_vs_fd.csv", std::ios::binary);
        if (fdcsv) {
            fdcsv << "location,component,parameter,step_size,hypad,"
                     "finite_difference,absolute_error,relative_error\n";
            fdcsv.setf(std::ios::scientific);
            fdcsv.precision(12);
        }

        double worst_of_best = 0.0;
        for (const Sample& s : kSamples) {
            const Mat4<OTIScalar> g = warp.evaluate(0.0, s.x, s.y, s.z, p);
            for (int d = 0; d < 3; ++d) {
                const double h00 = extract_gradient(g(0, 0), d);
                const double h01 = extract_gradient(g(0, 1), d);
                double best_here = 1e300;

                for (double h : steps) {
                    double lo[3] = {base[0], base[1], base[2]};
                    double hi[3] = {base[0], base[1], base[2]};
                    lo[d] -= h;
                    hi[d] += h;
                    const Mat4<double> ml = warp_d.evaluate(
                        0.0, s.x, s.y, s.z,
                        AlcubierreParameters<double>{lo[0], lo[1], lo[2]});
                    const Mat4<double> mh = warp_d.evaluate(
                        0.0, s.x, s.y, s.z,
                        AlcubierreParameters<double>{hi[0], hi[1], hi[2]});

                    const double fd00 = (mh(0, 0) - ml(0, 0)) / (2.0 * h);
                    const double fd01 = (mh(0, 1) - ml(0, 1)) / (2.0 * h);

                    // Error measure, valid across 11 orders of magnitude of
                    // derivative size:
                    //
                    //   denominator = max(|hypad|, kDerivFloor)
                    //
                    // Pure relative error is wrong at the bubble centre, where
                    // f == 1 for every R and sigma so dg/dR is analytically
                    // zero: HyPAD returns -5e-22 (denormal-scale roundoff), FD
                    // returns exactly 0, and |0 - 5e-22| / 5e-22 = 1.0 -- a
                    // "100% error" between two numbers that are both zero.
                    //
                    // The floor is the metric scale times a negligibility
                    // threshold: metric components here are O(1), so a
                    // derivative below 1e-9 cannot be resolved against them by
                    // any finite difference and is treated as an absolute
                    // comparison instead.
                    constexpr double kDerivFloor = 1e-9;
                    const double e00 = std::fabs(fd00 - h00);
                    const double e01 = std::fabs(fd01 - h01);
                    const double r00 = e00 / std::fmax(std::fabs(h00), kDerivFloor);
                    const double r01 = e01 / std::fmax(std::fabs(h01), kDerivFloor);
                    best_here = std::fmin(best_here, std::fmax(r00, r01));

                    if (fdcsv) {
                        fdcsv << s.where << ",g00," << param_name(d) << ',' << h
                              << ',' << h00 << ',' << fd00 << ',' << e00 << ','
                              << r00 << '\n';
                        fdcsv << s.where << ",g0x," << param_name(d) << ',' << h
                              << ',' << h01 << ',' << fd01 << ',' << e01 << ','
                              << r01 << '\n';
                    }
                }
                // MAX over (point, parameter) of the MIN over step sizes.
                //
                // Min-over-everything would pass vacuously: at the bubble centre
                // f == 1 identically for any R and sigma, so dg/dR and
                // dg/dsigma are exactly zero there and FD reproduces them
                // exactly, giving a global minimum of 0 no matter how wrong
                // every other point was. The statement worth making is "every
                // point and parameter has SOME step size that agrees".
                worst_of_best = std::fmax(worst_of_best, best_here);
            }
        }

        report.check_residual(
            "metric.sensitivity.vs_finite_difference",
            "for every point and parameter, the best centered finite-difference "
            "metric sensitivity agrees with the HyPAD value",
            VerificationLevel::D, "n/a",
            "6 points x 3 parameters x 6 step sizes", worst_of_best, 1e-3,
            sw.seconds(),
            "held to a floor, not to machine precision. At the exterior sample "
            "dg00/dsigma ~ 1.8e-11 against g00 ~ -1, so FD loses eleven orders "
            "to cancellation and bottoms out near 1e-3 however the step is "
            "chosen. That limit belongs to finite differences, not to HyPAD, "
            "and is precisely why the project uses AD");
    }

    // -----------------------------------------------------------------------
    // Sensitivity field export for the maps.
    // -----------------------------------------------------------------------
    {
        GridSpec spec;
        spec.nx = spec.ny = 81; spec.nz = 3;
        spec.x_min = -8; spec.x_max = 8;
        spec.y_min = -8; spec.y_max = 8;
        spec.z_min = -1; spec.z_max = 1;
        const Grid grid(spec);
        const int kc = grid.center_index_z();

        std::ofstream os(out_dir + "/sensitivity_slice_xy.csv", std::ios::binary);
        if (os && kc >= 0) {
            os << "x,y,dg00_dv,dg00_dR,dg00_dsigma,dg0x_dv,dg0x_dR,dg0x_dsigma\n";
            os.setf(std::ios::scientific);
            os.precision(10);
            for (int i = 0; i < grid.nx(); ++i)
                for (int j = 0; j < grid.ny(); ++j) {
                    const Mat4<OTIScalar> g = warp.evaluate(
                        0.0, grid.x(i), grid.y(j), grid.z(kc), p);
                    os << grid.x(i) << ',' << grid.y(j);
                    for (int d = 0; d < 3; ++d)
                        os << ',' << extract_gradient(g(0, 0), d);
                    for (int d = 0; d < 3; ++d)
                        os << ',' << extract_gradient(g(0, 1), d);
                    os << '\n';
                }
        }
    }
}

#endif  // WARP_HYPAD_HAVE_OTILIB

}  // namespace checks
}  // namespace warp_hypad
