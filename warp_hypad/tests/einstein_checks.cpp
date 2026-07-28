#include "einstein_checks.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "reference/alcubierre_curvature_reference.hpp"
#include "warp_hypad/core/units.hpp"
#include "warp_hypad/geometry/einstein.hpp"
#include "warp_hypad/metrics/analytic_metrics.hpp"
#include "warp_hypad/metrics/metric.hpp"
#include "warp_hypad/metrics/metric_field.hpp"

// ---------------------------------------------------------------------------
// Milestone A5 verification.
//
// The tolerance discipline is inherited from A4 unchanged: a residual is either
// an exactness claim at roundoff, or a truncation claim gated at twice the
// measured value WITH a companion convergence case. G_{mu nu} is a linear
// combination of R_{mu nu}, g_{mu nu} and R, so it inherits their error
// directly -- including the O(h^2) antisymmetry that A4 traced to the discrete
// product rule.
//
// SCOPE. Everything here is geometry. G_{mu nu} is computed, and the tensor
// called `inferred_stress_energy` is what Einstein's equation says would have to
// be present to source the chosen metric. Nothing in this file establishes that
// such a source can exist, that any energy condition holds, or that a
// coordinate component of G is an energy density -- it is not, and case
// einstein.g00_is_not_an_energy_density exists to make that concrete.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using geometry::CosmologicalConstantPolicy;
using geometry::CosmologicalConstantTreatment;
using geometry::EinsteinFields;
using geometry::TemporalPolicy;

constexpr double kPi = 3.14159265358979323846;
constexpr int kShellDepth = 2;

constexpr double kSchwarzschildMass = 1.0;
constexpr double kDeSitterLength = 2.0;
constexpr double kThetaMin = 1.0;
constexpr double kThetaMax = kPi - 1.0;
constexpr int kNPhi = 5;

Grid spherical_grid(int nr, int nth, double r0, double r1) {
    GridSpec s;
    s.nx = nr;  s.x_min = r0;        s.x_max = r1;
    s.ny = nth; s.y_min = kThetaMin; s.y_max = kThetaMax;
    s.nz = kNPhi; s.z_min = 0.0;     s.z_max = 1.0;
    return Grid(s);
}

Grid cartesian_grid(int n, double half) {
    GridSpec s;
    s.nx = s.ny = s.nz = n;
    s.x_min = s.y_min = s.z_min = -half;
    s.x_max = s.y_max = s.z_max = half;
    return Grid(s);
}

std::string grid_label(const Grid& g) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%dx%dx%d", g.nx(), g.ny(), g.nz());
    return std::string(buf);
}

template <typename MetricT>
auto sampler_for(const Grid& grid, const MetricT& m,
                 const typename MetricT::Parameters& p) {
    return [&grid, &m, &p](double t) {
        return metrics::build_metric_field<double>(grid, m, p, t);
    };
}

TemporalPolicy static_policy(const char* which) {
    return TemporalPolicy::exactly_static(
        std::string(which) + " has no t dependence in this chart", 1.0);
}

std::size_t center_point(const Grid& g) {
    return g.flat(g.center_index_x(), g.center_index_y(), g.center_index_z());
}

// Worst interior value of a per-point quantity.
template <typename Fn>
double interior_worst(const Grid& g, Fn&& fn, int depth = kShellDepth) {
    double worst = 0.0;
    for (std::size_t p = 0; p < g.num_points(); ++p) {
        if (geometry::is_shell_point(g, p, depth)) continue;
        const double v = fn(p);
        if (v > worst) worst = v;
    }
    return worst;
}

template <typename Fn>
double boundary_worst(const Grid& g, Fn&& fn, int depth = kShellDepth) {
    double worst = 0.0;
    for (std::size_t p = 0; p < g.num_points(); ++p) {
        if (!geometry::is_shell_point(g, p, depth)) continue;
        const double v = fn(p);
        if (v > worst) worst = v;
    }
    return worst;
}

const char* kScope =
    "A5 computes the geometric Einstein tensor and the stress-energy tensor "
    "inferred through Einstein's equation for the selected metric. It does not "
    "establish that this source can be generated physically, that energy "
    "conditions are satisfied, or that FTL travel is realizable.";

}  // namespace

// ===========================================================================
// Group: reference -- geometries whose Einstein tensor must vanish
// ===========================================================================

void run_einstein_reference_checks(VerificationReport& report,
                                   const std::string& out_dir) {
    Stopwatch sw;
    const auto lp = CosmologicalConstantPolicy::none();

    // --- Minkowski -----------------------------------------------------------
    {
        const Grid g = cartesian_grid(9, 2.0);
        metrics::MinkowskiMetric<double> m;
        metrics::MinkowskiMetric<double>::Parameters p;
        sw.reset();
        auto e = geometry::compute_einstein<double>(
            g, sampler_for(g, m, p), 0.0, static_policy("Minkowski"), lp);

        const double worst = e.ok ? interior_worst(g, [&](std::size_t q) {
            return max_abs_real(e.raw_einstein.matrix_at(q));
        }) : 1.0;

        report.check_exact("einstein.minkowski_vanishes",
                           "G_{mu nu} = 0 in flat Cartesian coordinates",
                           VerificationLevel::A, grid_label(g), 0.0, worst,
                           sw.seconds(),
                           "constant metric, so every stencil sums to exactly "
                           "zero -- an exactness claim");

        std::ofstream f(out_dir + "/einstein_minkowski.csv");
        f << "x,G00,G01,G11,G22,G33,trace\n";
        const int jm = g.center_index_y(), km = g.center_index_z();
        for (int i = 0; i < g.nx(); ++i) {
            const std::size_t q = g.flat(i, jm, km);
            const Mat4<double> G = e.raw_einstein.matrix_at(q);
            f << g.x(i) << "," << G(0, 0) << "," << G(0, 1) << "," << G(1, 1)
              << "," << G(2, 2) << "," << G(3, 3) << "," << e.trace.at(q) << "\n";
        }
    }

    // --- Flat space in spherical coordinates ---------------------------------
    //
    // The case that matters most in this group. A large connection must NOT
    // become a false source: if coordinate terms leaked into G, this is where a
    // vacuum would acquire matter.
    {
        const Grid g = spherical_grid(33, 33, 1.0, 5.0);
        metrics::FlatSphericalMetric<double> m;
        metrics::FlatSphericalMetric<double>::Parameters p;
        sw.reset();
        auto e = geometry::compute_einstein<double>(
            g, sampler_for(g, m, p), 0.0, static_policy("flat spherical"), lp);

        const double win = interior_worst(g, [&](std::size_t q) {
            return max_abs_real(e.raw_einstein.matrix_at(q));
        });
        const double wbd = boundary_worst(g, [&](std::size_t q) {
            return max_abs_real(e.raw_einstein.matrix_at(q));
        });

        report.check_residual(
            "einstein.flat_spherical_vanishes",
            "G_{mu nu} = 0 for flat space in ANY chart, despite Gamma != 0",
            VerificationLevel::B, grid_label(g), "r in [1,5], theta in [1,pi-1]",
            win, 2.0e-2, sw.seconds(),
            "measured 6.7e-3; the case where a coordinate connection would "
            "become false stress-energy. Boundary worst = " +
                std::to_string(wbd));

        std::ofstream f(out_dir + "/einstein_flat_spherical.csv");
        f << "r,max_abs_G,trace,einstein_square\n";
        const int jm = g.center_index_y(), km = g.center_index_z();
        for (int i = 0; i < g.nx(); ++i) {
            const std::size_t q = g.flat(i, jm, km);
            f << g.x(i) << "," << max_abs_real(e.raw_einstein.matrix_at(q)) << ","
              << e.trace.at(q) << "," << e.square.at(q) << "\n";
        }
    }

    // --- Schwarzschild -------------------------------------------------------
    {
        const Grid g = spherical_grid(33, 33, 4.0, 6.0);
        metrics::SchwarzschildMetric<double> m;
        metrics::SchwarzschildParameters<double> p(kSchwarzschildMass);
        sw.reset();
        const auto c = geometry::compute_curvature<double>(
            g, sampler_for(g, m, p), 0.0, static_policy("Schwarzschild"));
        auto e = geometry::assemble_einstein(g, c, lp);

        const double win = interior_worst(g, [&](std::size_t q) {
            return max_abs_real(e.raw_einstein.matrix_at(q));
        });
        report.check_residual("einstein.schwarzschild_vanishes",
                              "G_{mu nu} = 0 (vacuum Einstein equation)",
                              VerificationLevel::B, grid_label(g), "M=1, r in [4,6]",
                              win, 4.0e-3, sw.seconds(),
                              "measured 1.7e-3");

        // ... while the spacetime is emphatically not flat.
        const double k_mid = c.kretschmann.at(center_point(g));
        report.check_value("einstein.schwarzschild_still_curved",
                           "K = 48 M^2 / r^6 remains nonzero where G vanishes",
                           VerificationLevel::B, grid_label(g), "M=1, r=5",
                           metrics::SchwarzschildMetric<double>::kretschmann(
                               5.0, kSchwarzschildMass),
                           k_mid, 1.0e-3, sw.seconds(),
                           "a vanishing Einstein tensor is not a vanishing "
                           "curvature; without this the case above would be "
                           "satisfied by a pipeline that returns zero");

        std::ofstream f(out_dir + "/einstein_schwarzschild.csv");
        f << "r,max_abs_G,kretschmann,trace\n";
        const int jm = g.center_index_y(), km = g.center_index_z();
        for (int i = 0; i < g.nx(); ++i) {
            const std::size_t q = g.flat(i, jm, km);
            f << g.x(i) << "," << max_abs_real(e.raw_einstein.matrix_at(q)) << ","
              << c.kretschmann.at(q) << "," << e.trace.at(q) << "\n";
        }
    }

    // --- Chart invariance of a scalar contraction ----------------------------
    //
    // G_{mu nu} G^{mu nu} is a scalar, so it must agree between charts
    // describing the same geometry. Both are zero for Schwarzschild, so the
    // informative version of this test uses de Sitter, where the square is
    // nonzero -- see the desitter group.
    {
        const Grid gs = spherical_grid(33, 33, 4.0, 6.0);
        metrics::SchwarzschildMetric<double> ms;
        metrics::SchwarzschildParameters<double> ps(kSchwarzschildMass);
        const auto cs = geometry::compute_curvature<double>(
            gs, sampler_for(gs, ms, ps), 0.0, static_policy("Schwarzschild"));
        auto es = geometry::assemble_einstein(gs, cs, lp);

        const double rho0 = 2.0 + std::sqrt(15.0) / 2.0;
        GridSpec is;
        is.nx = 33; is.x_min = rho0 - 1.0;  is.x_max = rho0 + 1.0;
        is.ny = 33; is.y_min = kThetaMin;   is.y_max = kThetaMax;
        is.nz = kNPhi; is.z_min = 0.0;      is.z_max = 1.0;
        const Grid gi(is);
        metrics::SchwarzschildIsotropicMetric<double> mi;
        const auto ci = geometry::compute_curvature<double>(
            gi, sampler_for(gi, mi, ps), 0.0, static_policy("isotropic"));
        auto ei = geometry::assemble_einstein(gi, ci, lp);

        sw.reset();
        const double sq_s = es.square.at(center_point(gs));
        const double sq_i = ei.square.at(center_point(gi));
        report.check_residual(
            "einstein.square_chart_independent",
            "G_{mu nu} G^{mu nu} is a scalar and agrees between charts",
            VerificationLevel::A, grid_label(gs), "M=1, areal radius exactly 5",
            std::fabs(sq_s - sq_i), 1.0e-5, sw.seconds(),
            "both are near zero for a vacuum; the de Sitter group carries the "
            "same test with a nonzero value");
    }
}

// ===========================================================================
// Group: desitter -- where G != 0 and the spacetime is still a vacuum
// ===========================================================================

void run_einstein_desitter_checks(VerificationReport& report,
                                  const std::string& out_dir) {
    Stopwatch sw;
    const double L = kDeSitterLength;
    const double lambda = 3.0 / (L * L);

    const Grid g = spherical_grid(33, 33, 0.5, 1.5);
    metrics::DeSitterStaticMetric<double> m;
    metrics::DeSitterParameters<double> p(L);

    // Geometric Einstein tensor only: G = -(3/L^2) g. NONZERO.
    const auto lp_geom = CosmologicalConstantPolicy::none();
    sw.reset();
    const auto c = geometry::compute_curvature<double>(
        g, sampler_for(g, m, p), 0.0, static_policy("static de Sitter"));
    auto e = geometry::assemble_einstein(g, c, lp_geom);

    const double win = interior_worst(g, [&](std::size_t q) {
        const Mat4<double> G = e.symmetric_einstein.matrix_at(q);
        const Mat4<double> gm = c.g.matrix_at(q);
        double w = 0.0;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                w = std::fmax(w, std::fabs(G(a, b) + lambda * gm(a, b)));
        return w;
    });

    report.check_residual(
        "einstein.desitter_analytic_value",
        "G_{mu nu} = -(3/L^2) g_{mu nu} for static de Sitter",
        VerificationLevel::B, grid_label(g), "L=2, r in [0.5,1.5]", win, 6.0e-2,
        sw.seconds(),
        "measured 2.9e-2. NOTE the sign: with R_{mu nu} = +(3/L^2) g and "
        "R = +12/L^2, G = R_{mu nu} - g R/2 = (3/L^2 - 6/L^2) g = -(3/L^2) g");

    // G + Lambda g = 0 with Lambda = 3/L^2: de Sitter is a VACUUM solution with
    // a cosmological constant. It contains no matter.
    sw.reset();
    const Rank2Field<double> gl =
        geometry::einstein_plus_lambda_g(g, e, c.g, lambda);
    const double wvac = interior_worst(g, [&](std::size_t q) {
        return max_abs_real(gl.matrix_at(q));
    });
    report.check_residual(
        "einstein.desitter_vacuum_with_lambda",
        "G_{mu nu} + Lambda g_{mu nu} = 0 with Lambda = 3/L^2",
        VerificationLevel::B, grid_label(g), "L=2, Lambda=0.75", wvac, 6.0e-2,
        sw.seconds(),
        "the case that separates the GEOMETRIC Einstein tensor from MATTER "
        "stress-energy: G is large here, yet the spacetime contains nothing");

    // The three treatments are numerically different objects.
    sw.reset();
    CosmologicalConstantPolicy lp_eff;
    lp_eff.treatment = CosmologicalConstantTreatment::LambdaAsEffectiveStressEnergy;
    lp_eff.lambda = lambda;
    lp_eff.justification = "de Sitter, Lambda moved to the source side";
    auto e_eff = geometry::assemble_einstein(g, c, lp_eff);

    const std::size_t q0 = center_point(g);
    const double t_geom = max_abs_real(e.inferred_stress_energy.matrix_at(q0));
    const double t_eff = max_abs_real(e_eff.inferred_stress_energy.matrix_at(q0));

    report.check_predicate(
        "einstein.lambda_treatments_distinct",
        "GeometricEinsteinOnly and LambdaAsEffectiveStressEnergy give different "
        "inferred sources",
        VerificationLevel::A, std::fabs(t_geom - t_eff) > 1.0e-3,
        "the two policies differ",
        "geometric |T| = " + std::to_string(t_geom) + ", effective |T| = " +
            std::to_string(t_eff),
        sw.seconds(),
        "under the effective policy de Sitter's inferred source is ~0 because "
        "Lambda cancels G; under the geometric policy it is not. Choosing one "
        "silently is how a vacuum comes to be reported as containing matter");

    report.check_predicate(
        "einstein.lambda_policy_recorded",
        "every result carries the treatment that produced it",
        VerificationLevel::A,
        std::string(geometry::to_string(e.lambda_policy.treatment)) ==
                "geometric_einstein_only" &&
            std::string(geometry::to_string(e_eff.lambda_policy.treatment)) ==
                "lambda_as_effective_stress_energy",
        "policies recorded and distinct",
        std::string(geometry::to_string(e.lambda_policy.treatment)) + " / " +
            geometry::to_string(e_eff.lambda_policy.treatment),
        sw.seconds(), "");

    std::ofstream f(out_dir + "/einstein_desitter.csv");
    f << "r,G00,G11,G22,G33,expected_G00,trace,ricci_scalar,"
         "einstein_plus_lambda_g_max\n";
    const int jm = g.center_index_y(), km = g.center_index_z();
    for (int i = 0; i < g.nx(); ++i) {
        const std::size_t q = g.flat(i, jm, km);
        const Mat4<double> G = e.symmetric_einstein.matrix_at(q);
        const Mat4<double> gm = c.g.matrix_at(q);
        f << g.x(i) << "," << G(0, 0) << "," << G(1, 1) << "," << G(2, 2) << ","
          << G(3, 3) << "," << (-lambda * gm(0, 0)) << "," << e.trace.at(q) << ","
          << c.ricci_scalar.at(q) << "," << max_abs_real(gl.matrix_at(q)) << "\n";
    }
}

// ===========================================================================
// Group: alcubierre -- against the independent symbolic reference
// ===========================================================================

void run_einstein_alcubierre_checks(VerificationReport& report,
                                    const std::string& out_dir) {
    Stopwatch sw;
    const auto lp = CosmologicalConstantPolicy::none();

    const Grid g = cartesian_grid(33, 4.0);
    metrics::AlcubierreMetric<double> m(
        0.0, metrics::TimeTreatment::InertialTimeDependent);
    metrics::AlcubierreParameters<double> p(reference::kAlcubierreRefVelocity,
                                            reference::kAlcubierreRefRadius,
                                            reference::kAlcubierreRefSigma);
    const auto policy = TemporalPolicy::analytic_resample(
        g.dx() / reference::kAlcubierreRefVelocity, "stencil-aligned resampling");

    sw.reset();
    const auto c =
        geometry::compute_curvature<double>(g, sampler_for(g, m, p), 0.0, policy);
    auto e = geometry::assemble_einstein(g, c, lp);

    double worst = 0.0;
    std::string worst_at;
    for (int n = 0; n < reference::kAlcubierreCurvatureCount; ++n) {
        const auto& ref = reference::kAlcubierreCurvature[n];
        int bi = -1, bj = -1, bk = -1;
        for (int i = 0; i < g.nx(); ++i) if (g.x(i) == ref.x) bi = i;
        for (int j = 0; j < g.ny(); ++j) if (g.y(j) == ref.y) bj = j;
        for (int k = 0; k < g.nz(); ++k) if (g.z(k) == ref.z) bk = k;
        if (bi < 0 || bj < 0 || bk < 0) continue;

        const Mat4<double> G = e.symmetric_einstein.matrix_at(g.flat(bi, bj, bk));
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) {
                const double d = std::fabs(G(a, b) - ref.einstein[4 * a + b]);
                if (d > worst) {
                    worst = d;
                    worst_at = "(" + std::to_string(ref.x) + "," +
                               std::to_string(ref.y) + "," +
                               std::to_string(ref.z) + ") component " +
                               std::to_string(a) + std::to_string(b);
                }
            }
    }

    report.check_residual(
        "einstein.alcubierre_vs_symbolic",
        "G_{mu nu} componentwise against exact symbolic differentiation",
        VerificationLevel::B, grid_label(g), "v=0.5 R=2 sigma=1, t=0", worst,
        3.0e-2, sw.seconds(),
        "worst at " + worst_at +
            "; the reference assembles G from ITS OWN symbolic Ricci, never "
            "from this pipeline. Measured 1.2e-2 -- O(h^2), inherited from the "
            "Ricci tensor it is built from");

    // The point the naming discipline exists to protect.
    sw.reset();
    const std::size_t ctr = center_point(g);
    const double G00 = e.symmetric_einstein.matrix_at(ctr)(0, 0);
    report.check_predicate(
        "einstein.g00_is_not_an_energy_density",
        "G_00 is a coordinate component, not an observer-measured scalar",
        VerificationLevel::A, true, "recorded as a coordinate component",
        "G_00 = " + std::to_string(G00) + " at the bubble centre",
        sw.seconds(),
        "G_00 is POSITIVE across the sampled points for this configuration. "
        "That is not a statement that the inferred source has positive energy "
        "density: an energy density is T_{mu nu} u^mu u^nu for a specified "
        "normalised observer, which A6 introduces and A5 does not have. No "
        "output of A5 may be labelled an energy density");

    // --- Records -------------------------------------------------------------
    const int km = g.center_index_z();
    std::ofstream f(out_dir + "/alcubierre_einstein_slice_xy.csv");
    f << "x,y,g00,g01,ricci_scalar,G00_raw,G01_raw,G11_raw,G22_raw,G33_raw,"
         "G00_symmetric,G01_symmetric,G11_symmetric,G22_symmetric,G33_symmetric,"
         "einstein_trace,einstein_square,antisymmetry_residual\n";
    for (int i = 0; i < g.nx(); ++i)
        for (int j = 0; j < g.ny(); ++j) {
            const std::size_t q = g.flat(i, j, km);
            const Mat4<double> Gr = e.raw_einstein.matrix_at(q);
            const Mat4<double> Gs = e.symmetric_einstein.matrix_at(q);
            f << g.x(i) << "," << g.y(j) << "," << c.g.at(q, 0, 0) << ","
              << c.g.at(q, 0, 1) << "," << c.ricci_scalar.at(q) << ","
              << Gr(0, 0) << "," << Gr(0, 1) << "," << Gr(1, 1) << ","
              << Gr(2, 2) << "," << Gr(3, 3) << ","
              << Gs(0, 0) << "," << Gs(0, 1) << "," << Gs(1, 1) << ","
              << Gs(2, 2) << "," << Gs(3, 3) << "," << e.trace.at(q) << ","
              << e.square.at(q) << "," << symmetry_residual(Gr) << "\n";
        }

    const int jm = g.center_index_y();
    std::ofstream fc(out_dir + "/alcubierre_einstein_centerline.csv");
    fc << "x,g00,g01,ricci_scalar,G00_raw,G01_raw,G11_raw,G22_raw,G33_raw,"
          "G00_symmetric,G01_symmetric,G11_symmetric,G22_symmetric,"
          "G33_symmetric,einstein_trace,einstein_square,antisymmetry_residual\n";
    for (int i = 0; i < g.nx(); ++i) {
        const std::size_t q = g.flat(i, jm, km);
        const Mat4<double> Gr = e.raw_einstein.matrix_at(q);
        const Mat4<double> Gs = e.symmetric_einstein.matrix_at(q);
        fc << g.x(i) << "," << c.g.at(q, 0, 0) << "," << c.g.at(q, 0, 1) << ","
           << c.ricci_scalar.at(q) << "," << Gr(0, 0) << "," << Gr(0, 1) << ","
           << Gr(1, 1) << "," << Gr(2, 2) << "," << Gr(3, 3) << ","
           << Gs(0, 0) << "," << Gs(0, 1) << "," << Gs(1, 1) << "," << Gs(2, 2)
           << "," << Gs(3, 3) << "," << e.trace.at(q) << "," << e.square.at(q)
           << "," << symmetry_residual(Gr) << "\n";
    }
}

// ===========================================================================
// Group: identities
// ===========================================================================

void run_einstein_identity_checks(VerificationReport& report,
                                  const std::string& out_dir) {
    (void)out_dir;
    Stopwatch sw;
    const auto lp = CosmologicalConstantPolicy::none();

    const Grid g = cartesian_grid(17, 4.0);
    metrics::AlcubierreMetric<double> m(
        0.0, metrics::TimeTreatment::InertialTimeDependent);
    metrics::AlcubierreParameters<double> p(0.5, 2.0, 1.0);
    const auto policy =
        TemporalPolicy::analytic_resample(g.dx() / 0.5, "stencil aligned");

    sw.reset();
    const auto c =
        geometry::compute_curvature<double>(g, sampler_for(g, m, p), 0.0, policy);
    auto e = geometry::assemble_einstein(g, c, lp);

    // The raw antisymmetry is MEASURED, not projected away.
    const double asym_in = interior_worst(g, [&](std::size_t q) {
        return symmetry_residual(e.raw_einstein.matrix_at(q));
    });
    report.check_residual(
        "einstein.raw_antisymmetry_measured",
        "max |G_{mu nu} - G_{nu mu}| on the directly computed tensor",
        VerificationLevel::A, grid_label(g), "v=0.5 R=2 sigma=1", asym_in,
        2.0e-2, sw.seconds(),
        "measured 4.7e-3, inherited from the O(h^2) Ricci asymmetry A4 traced "
        "to the discrete product rule. Reported rather than hidden: "
        "symmetrising in place would destroy the only direct measure of how "
        "far the discretisation is from the continuum. Convergence is "
        "demonstrated separately");

    // The projection really is exactly symmetric.
    sw.reset();
    double proj = 0.0;
    for (std::size_t q = 0; q < g.num_points(); ++q)
        proj = std::fmax(proj, symmetry_residual(e.symmetric_einstein.matrix_at(q)));
    report.check_exact("einstein.symmetric_projection_is_symmetric",
                       "(1/2)(G + G^T) is exactly symmetric",
                       VerificationLevel::A, grid_label(g), 0.0, proj,
                       sw.seconds(),
                       "exact in floating point: a + b and b + a are the same "
                       "expression, so the halved sum is bit-identical either "
                       "way round");

    // Both remain reachable. A silent in-place symmetrisation would make these
    // the same object.
    sw.reset();
    double diff = 0.0;
    for (std::size_t q = 0; q < g.num_points(); ++q)
        diff = std::fmax(diff, max_abs_diff(e.raw_einstein.matrix_at(q),
                                            e.symmetric_einstein.matrix_at(q)));
    report.check_predicate(
        "einstein.raw_and_symmetric_both_available",
        "the raw tensor is preserved alongside its symmetric projection",
        VerificationLevel::A, diff > 1.0e-6, "the two fields differ",
        "max |raw - symmetric| = " + std::to_string(diff), sw.seconds(),
        "if this ever reads zero, the raw tensor has been silently "
        "symmetrised and the antisymmetry measurement above is vacuous");

    // G^mu_mu = -R in four dimensions.
    sw.reset();
    const double tr_in = interior_worst(g, [&](std::size_t q) {
        return std::fabs(e.trace.at(q) + c.ricci_scalar.at(q));
    });
    report.check_residual(
        "einstein.trace_identity",
        "G^mu_mu = -R in four dimensions",
        VerificationLevel::A, grid_label(g), "v=0.5 R=2 sigma=1", tr_in,
        1.0e-13, sw.seconds(),
        "algebraic: g^{mu nu}(R_{mu nu} - g_{mu nu} R/2) = R - 4R/2 = -R, using "
        "g^{mu nu} g_{mu nu} = 4. Roundoff, not truncation -- it needs no "
        "product rule and holds for whatever R_{mu nu} the discretisation "
        "produced");
}

// ===========================================================================
// Group: bianchi
// ===========================================================================

void run_bianchi_checks(VerificationReport& report, const std::string& out_dir) {
    (void)out_dir;
    Stopwatch sw;
    const auto lp = CosmologicalConstantPolicy::none();

    // Two geometries, so the absolute and relative measures can be compared
    // against each other rather than each asserted in isolation.
    double schw_scale = 0.0, schw_div = 0.0, schw_rel = 0.0;
    {
        const Grid g = spherical_grid(21, 21, 4.0, 6.0);
        metrics::SchwarzschildMetric<double> m;
        metrics::SchwarzschildParameters<double> p(kSchwarzschildMass);
        auto e = geometry::compute_einstein<double>(
            g, sampler_for(g, m, p), 0.0, static_policy("Schwarzschild"), lp);
        schw_scale = e.einstein_scale;
        schw_div = e.divergence_residual;
        schw_rel = e.divergence_relative;
    }

    // Alcubierre: G is genuinely nonzero, so both measures mean something.
    {
        const Grid g = cartesian_grid(21, 4.0);
        metrics::AlcubierreMetric<double> m(
            0.0, metrics::TimeTreatment::InertialTimeDependent);
        metrics::AlcubierreParameters<double> p(0.5, 2.0, 1.0);
        sw.reset();
        auto e = geometry::compute_einstein<double>(
            g, sampler_for(g, m, p), 0.0,
            TemporalPolicy::analytic_resample(g.dx() / 0.5, "stencil aligned"),
            lp);

        // The point of reporting both. Schwarzschild is a vacuum, so its G is
        // pure truncation error and its ABSOLUTE divergence is correspondingly
        // tiny -- which says nothing whatever about how well the identity is
        // satisfied. Normalising by max|G| makes the two geometries
        // comparable, and it is the relative figure that carries meaning.
        report.check_predicate(
            "bianchi.absolute_residual_is_not_enough",
            "the absolute divergence tracks the size of G, not the quality of "
            "the identity",
            VerificationLevel::A,
            schw_div < e.divergence_residual && schw_scale < e.einstein_scale,
            "the geometry with the smaller |G| has the smaller absolute residual",
            "Schwarzschild: max|G| = " + std::to_string(schw_scale) +
                ", |div G| = " + std::to_string(schw_div) + ", relative = " +
                std::to_string(schw_rel) + ";  Alcubierre: max|G| = " +
                std::to_string(e.einstein_scale) + ", |div G| = " +
                std::to_string(e.divergence_residual) + ", relative = " +
                std::to_string(e.divergence_relative),
            sw.seconds(),
            "a vacuum would otherwise appear to satisfy the contracted Bianchi "
            "identity better than a curved spacetime, purely because there is "
            "less of everything to be wrong about");

        report.check_residual(
            "bianchi.contracted_residual_absolute",
            "nabla_mu G^{mu nu} = 0, absolute",
            VerificationLevel::A, grid_label(g), "v=0.5 R=2 sigma=1",
            e.divergence_residual, 5.0e-2, sw.seconds(),
            "NOT an exact finite-difference identity: it needs a third nested "
            "derivative pass and the product rule. max|G| = " +
                std::to_string(e.einstein_scale) + ". Interior only, shell "
                "depth 3 -- the divergence is contaminated three points deep");

        report.check_residual(
            "bianchi.contracted_residual_relative",
            "nabla_mu G^{mu nu} = 0, normalised by max |G_{mu nu}|",
            VerificationLevel::A, grid_label(g), "v=0.5 R=2 sigma=1",
            e.divergence_relative, 2.0e-1, sw.seconds(),
            "the meaningful measure: an absolute residual can be made "
            "arbitrarily small by choosing a geometry with little curvature");
    }
}

// ===========================================================================
// Group: convergence
// ===========================================================================

void run_einstein_convergence_checks(VerificationReport& report,
                                     const std::string& out_dir) {
    Stopwatch sw;
    const auto lp = CosmologicalConstantPolicy::none();
    const int ns[3] = {17, 25, 33};

    // --- Antisymmetry --------------------------------------------------------
    {
        std::ofstream f(out_dir + "/einstein_convergence.csv");
        f << "case,n,h,antisymmetry_at_probe\n";
        double err[3] = {0, 0, 0}, hs[3] = {0, 0, 0};

        sw.reset();
        for (int a = 0; a < 3; ++a) {
            const Grid g = cartesian_grid(ns[a], 4.0);
            metrics::AlcubierreMetric<double> m(
                0.0, metrics::TimeTreatment::InertialTimeDependent);
            metrics::AlcubierreParameters<double> p(0.5, 2.0, 1.0);
            const auto c = geometry::compute_curvature<double>(
                g, sampler_for(g, m, p), 0.0,
                TemporalPolicy::analytic_resample(g.dx() / 0.5, "aligned"));
            auto e = geometry::assemble_einstein(g, c, lp);
            hs[a] = g.dx();

            // Fixed integer probe points, present on all three grids, and NOT
            // the bubble centre, where the antisymmetry vanishes by symmetry
            // and the measured order would be 0/0.
            const double probe[2][3] = {{2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
            double w = 0.0;
            for (int s = 0; s < 2; ++s) {
                int bi = -1, bj = -1, bk = -1;
                for (int i = 0; i < g.nx(); ++i) if (g.x(i) == probe[s][0]) bi = i;
                for (int j = 0; j < g.ny(); ++j) if (g.y(j) == probe[s][1]) bj = j;
                for (int k = 0; k < g.nz(); ++k) if (g.z(k) == probe[s][2]) bk = k;
                if (bi < 0 || bj < 0 || bk < 0) continue;
                w = std::fmax(w, symmetry_residual(
                                     e.raw_einstein.matrix_at(g.flat(bi, bj, bk))));
            }
            err[a] = w;
            f << "einstein_antisymmetry," << ns[a] << "," << hs[a] << "," << w
              << "\n";
        }

        const double order = std::log(err[1] / err[2]) / std::log(hs[1] / hs[2]);
        report.check_value(
            "einstein.antisymmetry_converges",
            "the antisymmetric part of G falls as O(h^2)",
            VerificationLevel::E, "17,25,33", "v=0.5 R=2 sigma=1", 2.0, order,
            0.4, sw.seconds(),
            "order 17->25 = " +
                std::to_string(std::log(err[0] / err[1]) /
                               std::log(hs[0] / hs[1])) +
                ". This is what licenses treating the antisymmetry as "
                "truncation error rather than as a defect, and what makes the "
                "symmetric projection legitimate rather than a cover-up");
    }

    // --- Contracted Bianchi --------------------------------------------------
    {
        std::ofstream f(out_dir + "/bianchi_convergence.csv");
        f << "case,n,h,divergence_absolute,divergence_relative,einstein_scale\n";
        double err[3] = {0, 0, 0}, hs[3] = {0, 0, 0};

        sw.reset();
        for (int a = 0; a < 3; ++a) {
            const Grid g = cartesian_grid(ns[a], 4.0);
            metrics::AlcubierreMetric<double> m(
                0.0, metrics::TimeTreatment::InertialTimeDependent);
            metrics::AlcubierreParameters<double> p(0.5, 2.0, 1.0);
            auto e = geometry::compute_einstein<double>(
                g, sampler_for(g, m, p), 0.0,
                TemporalPolicy::analytic_resample(g.dx() / 0.5, "aligned"), lp);
            hs[a] = g.dx();
            err[a] = e.divergence_relative;
            f << "contracted_bianchi," << ns[a] << "," << hs[a] << ","
              << e.divergence_residual << "," << e.divergence_relative << ","
              << e.einstein_scale << "\n";
        }

        const double order = std::log(err[1] / err[2]) / std::log(hs[1] / hs[2]);
        report.check_value(
            "bianchi.converges",
            "the contracted Bianchi residual falls under refinement",
            VerificationLevel::E, "17,25,33", "v=0.5 R=2 sigma=1, relative", 2.0,
            order, 0.9, sw.seconds(),
            "order 17->25 = " +
                std::to_string(std::log(err[0] / err[1]) /
                               std::log(hs[0] / hs[1])) +
                ". A wide band is asserted deliberately: the divergence is a "
                "THIRD nested derivative, so at these grids it carries more "
                "pre-asymptotic error than any other quantity in the project. "
                "What matters is that it converges at all -- a non-converging "
                "Bianchi residual would mean the Einstein tensor is wrong");
    }
}

// ===========================================================================
// Group: units
// ===========================================================================

void run_units_checks(VerificationReport& report, const std::string& out_dir) {
    Stopwatch sw;

    // T = G / (8 pi), geometrized. Exact by construction.
    sw.reset();
    const double G_sample = 0.0671161386;  // a real Alcubierre G_00 value
    const double t_geom = units::geometrized_stress_energy(G_sample);
    report.check_value("units.geometrized_conversion",
                       "T_{mu nu} = G_{mu nu} / (8 pi) in geometrized units",
                       VerificationLevel::A, "n/a (pointwise)", "G = c = 1",
                       G_sample / (8.0 * kPi), t_geom, 1.0e-15, sw.seconds(),
                       "the conversion the A5 gate verifies; SI is a "
                       "presentation layer on top of it");

    // SI requires a scale, and refuses to invent one.
    sw.reset();
    units::CoordinateScale unscaled;
    const double nan_result = units::si_stress_energy(G_sample, unscaled);
    report.check_predicate(
        "units.si_requires_an_explicit_scale",
        "no SI value is produced without a declared coordinate scale",
        VerificationLevel::A, std::isnan(nan_result), "NaN",
        std::isnan(nan_result) ? "NaN" : std::to_string(nan_result),
        sw.seconds(),
        "a grid coordinate is a pure number, not metres. Multiplying a "
        "dimensionless curvature by c^4/(8 pi G_N) and printing J/m^3 is the "
        "standard way to produce a confidently wrong number");

    // Dimensional consistency: T[J/m^3] = c^4/(8 pi G_N) * G / L0^2.
    sw.reset();
    units::CoordinateScale scale;
    scale.length_metres = 100.0;
    scale.name = "one grid unit = 100 m, so the R = 2 bubble is 200 m across";
    const double t_si = units::si_stress_energy(G_sample, scale);
    const double expected =
        units::einstein_coupling_si() * G_sample / (100.0 * 100.0);
    report.check_value("units.si_conversion_dimensional_consistency",
                       "T^SI = c^4/(8 pi G_N) * G_dimensionless / L0^2",
                       VerificationLevel::A, "n/a (pointwise)",
                       "L0 = 100 m", expected, t_si, 1.0e-6 * expected,
                       sw.seconds(),
                       "coupling c^4/(8 pi G_N) = " +
                           std::to_string(units::einstein_coupling_si()) +
                           " Pa per (1/m^2)");

    // Scaling law: halving the length scale quadruples the SI stress-energy.
    sw.reset();
    units::CoordinateScale half;
    half.length_metres = 50.0;
    const double t_half = units::si_stress_energy(G_sample, half);
    report.check_value("units.si_scales_as_inverse_length_squared",
                       "T^SI scales as 1/L0^2",
                       VerificationLevel::A, "n/a (pointwise)",
                       "L0 = 100 m vs 50 m", 4.0, t_half / t_si,
                       1.0e-9, sw.seconds(),
                       "the physical content of the conversion: the same "
                       "dimensionless geometry means a very different source "
                       "depending on how big you declare it to be");

    std::ofstream f(out_dir + "/units_and_scaling.json");
    const units::UnitsRecord rec;
    f << "{\n";
    f << "  \"system\": \"" << rec.system << "\",\n";
    f << "  \"coordinate_length_units\": \"" << rec.coordinate_length_units
      << "\",\n";
    f << "  \"coordinate_time_units\": \"" << rec.coordinate_time_units << "\",\n";
    f << "  \"curvature_units\": \"" << rec.curvature_units << "\",\n";
    f << "  \"stress_energy_units_geometrized\": \""
      << rec.stress_energy_units_geometrized << "\",\n";
    f << "  \"stress_energy_units_si\": \"" << rec.stress_energy_units_si
      << "\",\n";
    f << "  \"einstein_coupling_si_Pa_per_inverse_m2\": "
      << units::einstein_coupling_si() << ",\n";
    f << "  \"speed_of_light_m_per_s\": " << units::kSpeedOfLight << ",\n";
    f << "  \"newton_G_si\": " << units::kNewtonG << ",\n";
    f << "  \"example_scale_length_metres\": " << scale.length_metres << ",\n";
    f << "  \"example_scale_time_seconds\": " << scale.time_seconds() << ",\n";
    f << "  \"example_scale_name\": \"" << scale.name << "\",\n";
    f << "  \"conversion_note\": \"" << rec.conversion_note << "\",\n";
    f << "  \"scope\": \"" << kScope << "\"\n";
    f << "}\n";
}

// ===========================================================================
// Group: policy
// ===========================================================================

void run_einstein_policy_checks(VerificationReport& report,
                                const std::string& out_dir) {
    (void)out_dir;
    Stopwatch sw;
    const auto lp = CosmologicalConstantPolicy::none();
    const Grid g = cartesian_grid(9, 4.0);
    metrics::AlcubierreParameters<double> p(0.5, 2.0, 1.0);

    // Frozen time is refused all the way through A5, not just at A4.
    sw.reset();
    auto refused = geometry::compute_einstein<double>(
        g,
        [&](double t) {
            metrics::AlcubierreMetric<double> mm(
                0.0, metrics::TimeTreatment::FrozenTimeDiagnostic);
            return metrics::build_metric_field<double>(g, mm, p, t);
        },
        0.0,
        TemporalPolicy::refused("FrozenTimeDiagnostic discards d_0 g"), lp);

    report.check_predicate(
        "einstein.frozen_time_refused",
        "no Einstein tensor is produced from a discarded temporal derivative",
        VerificationLevel::A, !refused.ok, "ok == false",
        refused.ok ? "ok == true" : refused.refusal_reason, sw.seconds(),
        "the refusal propagates: A5 does not re-derive it, it inherits it");

    // Determinism.
    sw.reset();
    metrics::AlcubierreMetric<double> m(
        0.0, metrics::TimeTreatment::InertialTimeDependent);
    const auto pol =
        TemporalPolicy::analytic_resample(g.dx() / 0.5, "stencil aligned");
    auto e1 = geometry::compute_einstein<double>(g, sampler_for(g, m, p), 0.0,
                                                 pol, lp);
    auto e2 = geometry::compute_einstein<double>(g, sampler_for(g, m, p), 0.0,
                                                 pol, lp);
    double drift = 0.0;
    for (std::size_t q = 0; q < g.num_points(); ++q)
        drift = std::fmax(drift, max_abs_diff(e1.raw_einstein.matrix_at(q),
                                              e2.raw_einstein.matrix_at(q)));
    report.check_exact("einstein.deterministic",
                       "two identical runs produce bit-identical tensors",
                       VerificationLevel::A, grid_label(g), 0.0, drift,
                       sw.seconds(),
                       "no ordering, threading or accumulation nondeterminism "
                       "in the assembly");

    // The scope statement travels with the milestone.
    sw.reset();
    report.check_predicate(
        "einstein.scope_statement_recorded",
        "every A5 report carries its scope limitation",
        VerificationLevel::A, std::string(kScope).find("does not establish") !=
                                  std::string::npos,
        "scope statement present", kScope, sw.seconds(), "");
}

}  // namespace checks
}  // namespace warp_hypad
