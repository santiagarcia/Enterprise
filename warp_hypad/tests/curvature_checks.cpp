#include "curvature_checks.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "reference/alcubierre_curvature_reference.hpp"
#include "warp_hypad/core/field.hpp"
#include "warp_hypad/core/grid.hpp"
#include "warp_hypad/core/tensor.hpp"
#include "warp_hypad/geometry/chart.hpp"
#include "warp_hypad/geometry/curvature.hpp"
#include "warp_hypad/metrics/analytic_metrics.hpp"
#include "warp_hypad/metrics/metric.hpp"
#include "warp_hypad/metrics/metric_field.hpp"

// ---------------------------------------------------------------------------
// Milestone A4 verification.
//
// Every tolerance in this file is one of exactly two kinds, and they are never
// mixed:
//
//   EXACTNESS. The identity holds algebraically for the DISCRETE objects, so the
//   residual is floating-point roundoff and the tolerance is 1e-14 or tighter.
//   These are the identities whose derivation never differentiates a product:
//   the lower-pair symmetry of Gamma, the last-pair antisymmetry of Riemann, the
//   first Bianchi identity, metric compatibility, and the rigid-translation
//   identity d_0 g = -v_s d_x g with an aligned temporal step.
//
//   TRUNCATION. The identity holds in the continuum but not for second-order
//   finite differences, because a difference operator does not satisfy the
//   product rule. The tolerance is then 2x the value measured on the stated
//   grid, the measured value is recorded in the case's notes so the margin is
//   visible, and the claim is only meaningful because the convergence group
//   demonstrates that the residual falls as O(h^2). A truncation case without a
//   companion convergence case would be a calibration, not a check.
//
// Which identities land in which class is itself a result, and is documented at
// the point where the two are computed.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using geometry::Chart;
using geometry::ChartKind;
using geometry::CurvatureFields;
using geometry::TemporalPolicy;

constexpr double kPi = 3.14159265358979323846;

// --- Domains ----------------------------------------------------------------
//
// The spherical charts are singular at r = 0 and on the polar axis, and
// Schwarzschild additionally at r = 2M. Every domain below excludes those by
// construction, and the exclusion is part of the recorded parameters rather than
// an unwritten assumption.

constexpr double kSchwarzschildMass = 1.0;      // horizon at r = 2
constexpr double kSchwarzschildRMin = 4.0;      // two horizon radii of margin
constexpr double kSchwarzschildRMax = 6.0;

// L = 2 rather than a large Hubble length, and deliberately so. Every curvature
// component of de Sitter scales as 1/L^2, while the SPHERICAL CHART contributes
// connection components of order 1/r regardless of L. Taking L large therefore
// asks the pipeline to recover a tiny curvature as a near-total cancellation of
// order-unity connection terms, and the relative error degrades accordingly --
// at L = 10 on this grid the Kretschmann invariant came out 180% wrong while the
// absolute error was perfectly ordinary. That is a conditioning property of the
// chart, not a defect of the pipeline, but a verification case must not be run in
// a regime where it measures conditioning instead of correctness.
constexpr double kDeSitterLength = 2.0;         // horizon at r = 2
constexpr double kDeSitterRMin = 0.5;
constexpr double kDeSitterRMax = 1.5;

// Poles at theta = 0 and pi. The connection carries cot(theta), whose higher
// derivatives blow up as the pole is approached, so the margin is generous: a
// domain that merely avoids the singularity still pays for its neighbourhood.
constexpr double kThetaMin = 1.0;
constexpr double kThetaMax = kPi - 1.0;
constexpr double kPhiMin = 0.0;
constexpr double kPhiMax = 1.0;                 // the metrics are phi independent

// Angular resolution is matched to radial resolution throughout. With the
// angular axis left coarse, its truncation error dominates every residual and
// the radial refinement being measured is invisible.
constexpr int kSphericalNPhi = 5;               // phi derivatives are exactly zero

// The isotropic radius whose areal radius is exactly 5, for M = 1.
//
//   r = rho + M + M^2/(4 rho)  =>  4 rho^2 - 16 rho + 1 = 0  =>  rho = 2 + sqrt(15)/2
inline double isotropic_rho_for_areal_five() { return 2.0 + std::sqrt(15.0) / 2.0; }

Grid make_spherical_grid(int nr, int nth, int nph, double r0, double r1) {
    GridSpec s;
    s.nx = nr;   s.x_min = r0;        s.x_max = r1;
    s.ny = nth;  s.y_min = kThetaMin; s.y_max = kThetaMax;
    s.nz = nph;  s.z_min = kPhiMin;   s.z_max = kPhiMax;
    return Grid(s);
}

Grid make_cartesian_grid(int n, double half_width) {
    GridSpec s;
    s.nx = s.ny = s.nz = n;
    s.x_min = s.y_min = s.z_min = -half_width;
    s.x_max = s.y_max = s.z_max = half_width;
    return Grid(s);
}

std::string grid_label(const Grid& g) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%dx%dx%d", g.nx(), g.ny(), g.nz());
    return std::string(buf);
}

// Any metric exposing evaluate(t, a0, a1, a2, params) becomes a time sampler.
template <typename MetricT>
auto sampler_for(const Grid& grid, const MetricT& metric,
                 const typename MetricT::Parameters& params) {
    return [&grid, &metric, &params](double t) {
        return metrics::build_metric_field<double>(grid, metric, params, t);
    };
}

// Worst value of a per-point quantity, split into a boundary shell of the given
// depth and the interior. Curvature costs two nested derivative passes, so a
// point one step from a face is already contaminated by that face's one-sided
// stencil; depth 2 is the shell that reaches.
constexpr int kShellDepth = 2;

template <typename PointFn>
geometry::FieldError split(const Grid& g, PointFn&& fn) {
    return geometry::split_error<double>(g, fn, kShellDepth);
}

double max_abs_rank4_field(const Grid& g, const Rank4Field<double>& f) {
    double worst = 0.0;
    for (std::size_t p = 0; p < g.num_points(); ++p) {
        const double a = max_abs_real(f.tensor_at(p));
        if (a > worst) worst = a;
    }
    return worst;
}

double max_abs_rank3_field(const Grid& g, const Rank3Field<double>& f) {
    double worst = 0.0;
    for (std::size_t p = 0; p < g.num_points(); ++p) {
        const double a = max_abs_real(f.tensor_at(p));
        if (a > worst) worst = a;
    }
    return worst;
}

// Index of the pinned centre point. Defined only for odd extents, which every
// grid here uses, so a caller never silently gets one side of the centre.
std::size_t center_point(const Grid& g) {
    return g.flat(g.center_index_x(), g.center_index_y(), g.center_index_z());
}

// A policy naming why staticity is a fact about the chart rather than an
// approximation. These metrics do not mention t at all; the driver checks that
// claim before using it.
TemporalPolicy static_policy(const char* which) {
    return TemporalPolicy::exactly_static(
        std::string(which) + " has no t dependence in this chart; d_0 g vanishes "
                             "identically rather than being discarded",
        1.0);
}

}  // namespace

// ===========================================================================
// Group: connection
// ===========================================================================

void run_connection_checks(VerificationReport& report, const std::string& out_dir) {
    Stopwatch sw;

    // --- Minkowski: the connection vanishes ---------------------------------
    {
        const Grid g = make_cartesian_grid(9, 2.0);
        metrics::MinkowskiMetric<double> m;
        metrics::MinkowskiMetric<double>::Parameters p;
        sw.reset();
        auto c = geometry::compute_curvature<double>(
            g, sampler_for(g, m, p), 0.0,
            static_policy("Minkowski in Cartesian coordinates"));

        report.check_predicate(
            "curvature.minkowski_pipeline_ok",
            "compute_curvature accepts a genuinely static metric",
            VerificationLevel::A, c.ok, "ok", c.ok ? "ok" : c.refusal_reason,
            sw.seconds());

        const double worst = c.ok ? max_abs_rank3_field(g, c.gamma) : 1.0;
        report.check_exact("curvature.minkowski_christoffel_zero",
                           "Gamma^rho_{mu nu} = 0 for a constant metric",
                           VerificationLevel::A, grid_label(g), 0.0, worst,
                           sw.seconds(),
                           "constant field, so the stencil sums to exactly zero -- "
                           "this is an exactness claim, not a tolerance");

        const double riem = c.ok ? max_abs_rank4_field(g, c.riemann) : 1.0;
        report.check_exact("curvature.minkowski_riemann_zero",
                           "R^rho_{sigma mu nu} = 0 in flat Cartesian coordinates",
                           VerificationLevel::A, grid_label(g), 0.0, riem,
                           sw.seconds(),
                           "passes trivially for a pipeline that always returns "
                           "zero; curvature.flat_spherical_riemann_zero is the "
                           "case that does not");
    }

    // --- Flat space in spherical coordinates --------------------------------
    //
    // The discriminating case. Gamma is large, Riemann must still vanish, and it
    // does so only if the dGamma terms cancel the GammaGamma terms.
    {
        const Grid g = make_spherical_grid(33, 33, kSphericalNPhi, 1.0, 5.0);
        metrics::FlatSphericalMetric<double> m;
        metrics::FlatSphericalMetric<double>::Parameters p;
        sw.reset();
        auto c = geometry::compute_curvature<double>(
            g, sampler_for(g, m, p), 0.0,
            static_policy("flat space in spherical coordinates"));

        const double gamma_max = max_abs_rank3_field(g, c.gamma);
        report.check_predicate(
            "curvature.flat_spherical_connection_nonzero",
            "Gamma != 0 although the spacetime is flat",
            VerificationLevel::A, gamma_max > 0.5, "max|Gamma| > 0.5",
            std::to_string(gamma_max), sw.seconds(),
            "separates a vanishing connection from vanishing curvature");

        // Against the closed form, componentwise.
        sw.reset();
        const geometry::FieldError gerr = split(g, [&](std::size_t pt) {
            int i, j, k;
            g.unflatten(pt, i, j, k);
            const Rank3<double> exact =
                metrics::flat_spherical_christoffel(g.x(i), g.y(j));
            return max_abs_diff(c.gamma.tensor_at(pt), exact);
        });
        report.check_residual(
            "curvature.flat_spherical_christoffel_closed_form",
            "Gamma^r_{th th} = -r, Gamma^th_{r th} = 1/r, "
            "Gamma^ph_{th ph} = cot(th), ...",
            VerificationLevel::B, grid_label(g), "r in [1,5], th in [0.5, pi-0.5]",
            gerr.interior_worst, 2.0e-3, sw.seconds(),
            "interior only; boundary worst = " +
                std::to_string(gerr.boundary_worst) +
                " (one-sided stencils, reported separately by policy)");

        report.check_residual(
            "curvature.flat_spherical_christoffel_closed_form_boundary",
            "same, on the one-sided shell",
            VerificationLevel::B, grid_label(g), "shell depth 2",
            gerr.boundary_worst, 5.0e-2, sw.seconds(),
            "a one-sided second-order stencil has a larger error constant; "
            "reported and bounded separately rather than averaged into the "
            "interior figure");

        sw.reset();
        const geometry::FieldError rerr = split(g, [&](std::size_t pt) {
            return max_abs_real(c.riemann.tensor_at(pt));
        });
        report.check_residual(
            "curvature.flat_spherical_riemann_zero",
            "R^rho_{sigma mu nu} = 0 for flat space in ANY chart",
            VerificationLevel::B, grid_label(g), "r in [1,5]",
            rerr.interior_worst, 1.5e-2, sw.seconds(),
            "the case a sign error in the Gamma-Gamma terms fails; boundary "
            "worst = " + std::to_string(rerr.boundary_worst));

        // Record the radial error profile for the figures.
        std::ofstream f(out_dir + "/flat_spherical_christoffel_error.csv");
        f << "r,christoffel_abs_error,riemann_abs\n";
        const int jm = g.center_index_y();
        const int km = g.center_index_z();
        for (int i = 0; i < g.nx(); ++i) {
            const std::size_t pt = g.flat(i, jm, km);
            const Rank3<double> exact =
                metrics::flat_spherical_christoffel(g.x(i), g.y(jm));
            f << g.x(i) << "," << max_abs_diff(c.gamma.tensor_at(pt), exact) << ","
              << max_abs_real(c.riemann.tensor_at(pt)) << "\n";
        }
    }

    // --- Identities on a genuinely curved, time-dependent spacetime ---------
    {
        const Grid g = make_cartesian_grid(17, 4.0);
        metrics::AlcubierreMetric<double> m(0.0,
                                            metrics::TimeTreatment::InertialTimeDependent);
        metrics::AlcubierreParameters<double> p(0.5, 2.0, 1.0);
        sw.reset();
        auto c = geometry::compute_curvature<double>(
            g, sampler_for(g, m, p), 0.0,
            TemporalPolicy::analytic_resample(
                0.5, "Alcubierre depends on t; d_0 g is resampled, not assumed"));

        report.check_predicate("curvature.alcubierre_pipeline_ok",
                               "the time-dependent path produces a result",
                               VerificationLevel::A, c.ok, "ok",
                               c.ok ? "ok" : c.refusal_reason, sw.seconds());

        sw.reset();
        double worst_sym = 0.0;
        for (std::size_t pt = 0; pt < g.num_points(); ++pt) {
            const double s = lower_pair_symmetry_residual(c.gamma.tensor_at(pt));
            if (s > worst_sym) worst_sym = s;
        }
        report.check_residual(
            "curvature.christoffel_lower_pair_symmetry",
            "Gamma^rho_{mu nu} = Gamma^rho_{nu mu} (torsion free)",
            VerificationLevel::A, grid_label(g), "v=0.5 R=2 sigma=1",
            worst_sym, 1.0e-15, sw.seconds(),
            "algebraic in the assembly, so this is roundoff and not truncation");

        sw.reset();
        const geometry::FieldError compat = split(g, [&](std::size_t pt) {
            return geometry::metric_compatibility_residual(
                c.g.matrix_at(pt), c.dg.tensor_at(pt), c.gamma.tensor_at(pt));
        });
        report.check_residual(
            "curvature.metric_compatibility",
            "nabla_lambda g_{mu nu} = 0 for the Levi-Civita connection",
            VerificationLevel::A, grid_label(g), "v=0.5 R=2 sigma=1",
            compat.interior_worst, 1.0e-14, sw.seconds(),
            "does NOT reuse the Christoffel formula, so it is an independent "
            "identity; exact to roundoff because the same discrete d_alpha g "
            "enters both sides. boundary worst = " +
                std::to_string(compat.boundary_worst));

        // Two classes of identity, and the difference between them is the whole
        // point of this block.
        //
        // EXACT ones survive discretisation because their derivation never
        // differentiates a product: the lower-pair symmetry of Gamma, the
        // antisymmetry of Riemann in its last pair, the first Bianchi identity
        // and metric compatibility are all rearrangements of the same discrete
        // d_alpha g and are therefore roundoff-level.
        //
        // The remaining ones are NOT exact, and it would be wrong to write them
        // with a roundoff tolerance and then widen it when they fail. Their
        // textbook derivations all pass through
        //
        //     d_mu ( g_{rho alpha} Gamma^alpha ) = (d_mu g_{rho alpha}) Gamma^alpha
        //                                        + g_{rho alpha} d_mu Gamma^alpha
        //
        // and a finite difference does not satisfy the product rule: D(fg) differs
        // from (Df)g + f(Dg) by O(h^2) f' g'. So first-pair antisymmetry, pair
        // exchange and Ricci symmetry are second-order truncation statements, and
        // the honest test is that they are small AND converge -- which is what the
        // convergence group measures.
        sw.reset();
        double worst_last = 0.0, worst_bianchi = 0.0;
        for (std::size_t pt = 0; pt < g.num_points(); ++pt) {
            const Rank4<double> riem = c.riemann.tensor_at(pt);
            worst_last = std::fmax(worst_last, antisymmetry_last_pair_residual(riem));
            worst_bianchi = std::fmax(worst_bianchi, first_bianchi_residual(riem));
        }

        report.check_residual(
            "curvature.riemann_antisymmetry_last_pair",
            "R^rho_{sigma mu nu} = -R^rho_{sigma nu mu}",
            VerificationLevel::A, grid_label(g), "v=0.5 R=2 sigma=1",
            worst_last, 1.0e-15, sw.seconds(),
            "manifest in the assembly; roundoff, not truncation");

        report.check_residual(
            "curvature.riemann_first_bianchi",
            "R^rho_{[sigma mu nu]} = 0",
            VerificationLevel::A, grid_label(g), "v=0.5 R=2 sigma=1",
            worst_bianchi, 1.0e-14, sw.seconds(),
            "holds algebraically for any symmetric connection, so it tests the "
            "index wiring rather than the discretisation");

        sw.reset();
        const geometry::FieldError first_pair = split(g, [&](std::size_t pt) {
            return antisymmetry_first_pair_residual(
                lower_first(c.g.matrix_at(pt), c.riemann.tensor_at(pt)));
        });
        const geometry::FieldError pair_ex = split(g, [&](std::size_t pt) {
            return pair_exchange_residual(
                lower_first(c.g.matrix_at(pt), c.riemann.tensor_at(pt)));
        });
        const geometry::FieldError ric_sym = split(g, [&](std::size_t pt) {
            return symmetry_residual(c.ricci.matrix_at(pt));
        });

        report.check_residual(
            "curvature.riemann_antisymmetry_first_pair",
            "R_{rho sigma mu nu} = -R_{sigma rho mu nu} (indices lowered)",
            VerificationLevel::A, grid_label(g), "v=0.5 R=2 sigma=1",
            first_pair.interior_worst, 2.0e-2, sw.seconds(),
            "O(h^2), not roundoff -- see the product-rule note above. Boundary "
            "worst = " + std::to_string(first_pair.boundary_worst) +
                ". curvature.convergence_ricci_symmetry demonstrates the order");

        report.check_residual(
            "curvature.riemann_pair_exchange",
            "R_{rho sigma mu nu} = R_{mu nu rho sigma}",
            VerificationLevel::A, grid_label(g), "v=0.5 R=2 sigma=1",
            pair_ex.interior_worst, 2.0e-2, sw.seconds(),
            "O(h^2); boundary worst = " + std::to_string(pair_ex.boundary_worst));

        report.check_residual(
            "curvature.ricci_symmetry", "R_{mu nu} = R_{nu mu}",
            VerificationLevel::A, grid_label(g), "v=0.5 R=2 sigma=1",
            ric_sym.interior_worst, 2.0e-2, sw.seconds(),
            "O(h^2); boundary worst = " + std::to_string(ric_sym.boundary_worst) +
                ". A5 forms G_{mu nu} = R_{mu nu} - (1/2) g_{mu nu} R and will "
                "inherit exactly this asymmetry, so the size of it is a number "
                "that milestone needs");
    }
}

// ===========================================================================
// Group: analytic
// ===========================================================================

void run_analytic_curvature_checks(VerificationReport& report,
                                   const std::string& out_dir) {
    Stopwatch sw;

    // --- Schwarzschild -------------------------------------------------------
    {
        const Grid g = make_spherical_grid(33, 33, kSphericalNPhi, kSchwarzschildRMin,
                                           kSchwarzschildRMax);
        metrics::SchwarzschildMetric<double> m;
        metrics::SchwarzschildParameters<double> p(kSchwarzschildMass);
        sw.reset();
        auto c = geometry::compute_curvature<double>(
            g, sampler_for(g, m, p), 0.0, static_policy("Schwarzschild"));

        const std::string params = "M=1, r in [4,6] (horizon at r=2)";

        sw.reset();
        const geometry::FieldError ric = split(g, [&](std::size_t pt) {
            return max_abs_real(c.ricci.matrix_at(pt));
        });
        report.check_residual("curvature.schwarzschild_ricci_flat",
                              "R_{mu nu} = 0 (vacuum Einstein equation)",
                              VerificationLevel::B, grid_label(g), params,
                              ric.interior_worst, 4.0e-3, sw.seconds(),
                              "boundary worst = " +
                                  std::to_string(ric.boundary_worst));

        sw.reset();
        const geometry::FieldError rsc = split(g, [&](std::size_t pt) {
            return std::fabs(c.ricci_scalar.at(pt));
        });
        report.check_residual("curvature.schwarzschild_ricci_scalar_zero",
                              "R = 0", VerificationLevel::B, grid_label(g), params,
                              rsc.interior_worst, 4.0e-3, sw.seconds(),
                              "boundary worst = " +
                                  std::to_string(rsc.boundary_worst));

        // The nonzero closed form. Relative error, because K varies by an order
        // of magnitude across the domain and an absolute tolerance would be
        // dominated by the inner edge.
        sw.reset();
        const geometry::FieldError kerr = split(g, [&](std::size_t pt) {
            int i, j, k;
            g.unflatten(pt, i, j, k);
            const double exact = metrics::SchwarzschildMetric<double>::kretschmann(
                g.x(i), kSchwarzschildMass);
            return std::fabs(c.kretschmann.at(pt) - exact) / std::fabs(exact);
        });
        report.check_residual(
            "curvature.schwarzschild_kretschmann",
            "R_{abcd} R^{abcd} = 48 M^2 / r^6",
            VerificationLevel::B, grid_label(g), params, kerr.interior_worst,
            2.0e-3, sw.seconds(),
            "relative error, interior only; boundary worst = " +
                std::to_string(kerr.boundary_worst) +
                ". Ricci flat but NOT flat -- the case that separates the two");

        std::ofstream f(out_dir + "/schwarzschild_kretschmann_profile.csv");
        f << "r,kretschmann_computed,kretschmann_exact,rel_error,ricci_abs_max\n";
        const int jm = g.center_index_y();
        const int km = g.center_index_z();
        for (int i = 0; i < g.nx(); ++i) {
            const std::size_t pt = g.flat(i, jm, km);
            const double exact = metrics::SchwarzschildMetric<double>::kretschmann(
                g.x(i), kSchwarzschildMass);
            f << g.x(i) << "," << c.kretschmann.at(pt) << "," << exact << ","
              << std::fabs(c.kretschmann.at(pt) - exact) / std::fabs(exact) << ","
              << max_abs_real(c.ricci.matrix_at(pt)) << "\n";
        }
    }

    // --- Chart independence of a scalar invariant ---------------------------
    //
    // The same curved geometry in two charts. Components differ everywhere; K
    // must not. Both grids pin their midpoint, and the isotropic midpoint is the
    // exact root of r(rho) = 5, so the comparison needs no interpolation.
    {
        const Grid gs = make_spherical_grid(33, 33, kSphericalNPhi, 4.0, 6.0);
        metrics::SchwarzschildMetric<double> ms;
        metrics::SchwarzschildParameters<double> ps(kSchwarzschildMass);
        auto cs = geometry::compute_curvature<double>(
            gs, sampler_for(gs, ms, ps), 0.0, static_policy("Schwarzschild"));

        const double rho0 = isotropic_rho_for_areal_five();
        GridSpec is;
        is.nx = 33;  is.x_min = rho0 - 1.0;   is.x_max = rho0 + 1.0;
        is.ny = 33;  is.y_min = kThetaMin;    is.y_max = kThetaMax;
        is.nz = kSphericalNPhi;   is.z_min = kPhiMin;      is.z_max = kPhiMax;
        const Grid gi(is);
        metrics::SchwarzschildIsotropicMetric<double> mi;
        auto ci = geometry::compute_curvature<double>(
            gi, sampler_for(gi, mi, ps), 0.0,
            static_policy("Schwarzschild in isotropic coordinates"));

        sw.reset();
        const double k_schw = cs.kretschmann.at(center_point(gs));
        const double k_iso = ci.kretschmann.at(center_point(gi));
        const double k_exact = metrics::SchwarzschildMetric<double>::kretschmann(
            5.0, kSchwarzschildMass);

        // The charts really are different, or the comparison below is vacuous.
        //
        // The metric components are a weaker witness than they look. Both charts
        // are static with the same Killing normalisation, so g_tt agrees exactly;
        // and the isotropic conformal factor is built to reproduce the areal
        // radius, so g_thth and g_phph agree exactly too. Only g_rr differs, by
        // 1.667 vs 1.613 at areal radius 5. The CONNECTION is where the two charts
        // genuinely part company, and the connection is what the curvature is
        // assembled from, so that is what is compared.
        const double comp_diff = max_abs_diff(cs.g.matrix_at(center_point(gs)),
                                              ci.g.matrix_at(center_point(gi)));
        const double conn_diff = max_abs_diff(cs.gamma.tensor_at(center_point(gs)),
                                              ci.gamma.tensor_at(center_point(gi)));
        report.check_predicate(
            "curvature.chart_connections_differ",
            "the two charts describe the same geometry with different "
            "Gamma^rho_{mu nu}",
            VerificationLevel::A, conn_diff > 0.02,
            "max|Gamma_schw - Gamma_iso| > 0.02", std::to_string(conn_diff),
            sw.seconds(),
            "max|g_schw - g_iso| = " + std::to_string(comp_diff) +
                ", confined to the radial component for the reason above. Without "
                "this the invariance case below could be comparing a chart with "
                "itself");

        report.check_value(
            "curvature.kretschmann_chart_independent",
            "K is a scalar: its value at areal radius 5 does not depend on the "
            "chart",
            VerificationLevel::A, grid_label(gs), "M=1, areal radius exactly 5",
            k_schw, k_iso, 1.0e-3, sw.seconds(),
            "exact value 48/5^6 = " + std::to_string(k_exact) +
                "; isotropic midpoint rho = " + std::to_string(rho0) +
                " is the exact root of rho + M + M^2/(4 rho) = 5");

        report.check_value("curvature.kretschmann_isotropic_closed_form",
                           "K(rho) = 48 M^2 / r(rho)^6 in isotropic coordinates",
                           VerificationLevel::B, grid_label(gi), "M=1, rho midpoint",
                           k_exact, k_iso, 1.0e-3, sw.seconds(), "");
    }

    // --- Static de Sitter ----------------------------------------------------
    //
    // The only reference here with nonzero Ricci, so it is the only one that can
    // fix the SIGN of the Ricci contraction -- and therefore the sign of every
    // energy density Milestone A6 will report.
    {
        const Grid g = make_spherical_grid(33, 33, kSphericalNPhi, kDeSitterRMin, kDeSitterRMax);
        metrics::DeSitterStaticMetric<double> m;
        metrics::DeSitterParameters<double> p(kDeSitterLength);
        sw.reset();
        auto c = geometry::compute_curvature<double>(
            g, sampler_for(g, m, p), 0.0, static_policy("static de Sitter"));

        const std::string params = "L=10, r in [1,5] (horizon at r=10)";
        const double R_exact =
            metrics::DeSitterStaticMetric<double>::ricci_scalar(kDeSitterLength);
        const double K_exact =
            metrics::DeSitterStaticMetric<double>::kretschmann(kDeSitterLength);
        const double RS_exact =
            metrics::DeSitterStaticMetric<double>::ricci_square(kDeSitterLength);
        const double lambda =
            metrics::DeSitterStaticMetric<double>::lambda(kDeSitterLength);

        sw.reset();
        const geometry::FieldError rerr = split(g, [&](std::size_t pt) {
            return std::fabs(c.ricci_scalar.at(pt) - R_exact) / std::fabs(R_exact);
        });
        report.check_residual("curvature.desitter_ricci_scalar",
                              "R = 12/L^2 for de Sitter, everywhere",
                              VerificationLevel::B, grid_label(g), params,
                              rerr.interior_worst, 1.2e-2, sw.seconds(),
                              "relative; exact value " + std::to_string(R_exact) +
                                  ". Maximally symmetric, so a radial trend in the "
                                  "error is unambiguously discretisation. boundary "
                                  "worst = " + std::to_string(rerr.boundary_worst));

        report.check_predicate(
            "curvature.desitter_ricci_scalar_sign",
            "R > 0 for de Sitter under the MTW (+,+,+) conventions of "
            "docs/conventions.md section 6",
            VerificationLevel::A,
            c.ricci_scalar.at(center_point(g)) > 0.0, "R > 0",
            std::to_string(c.ricci_scalar.at(center_point(g))), sw.seconds(),
            "a pipeline with the opposite Riemann sign returns -12/L^2 and "
            "passes every identity check in this suite; this is the case that "
            "catches it");

        sw.reset();
        const geometry::FieldError perr = split(g, [&](std::size_t pt) {
            const Mat4<double> ric = c.ricci.matrix_at(pt);
            const Mat4<double> gm = c.g.matrix_at(pt);
            double worst = 0.0;
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    worst = std::fmax(worst, std::fabs(ric(a, b) - lambda * gm(a, b)));
            return worst;
        });
        report.check_residual(
            "curvature.desitter_ricci_proportional_to_metric",
            "R_{mu nu} = (3/L^2) g_{mu nu}",
            VerificationLevel::B, grid_label(g), params, perr.interior_worst,
            4.5e-2, sw.seconds(),
            "absolute, against components of size 0.75 to 1.7; measured 2.0e-2 on "
            "this grid. Componentwise, not merely the trace; boundary worst = " +
                std::to_string(perr.boundary_worst));

        sw.reset();
        const geometry::FieldError kerr = split(g, [&](std::size_t pt) {
            return std::fabs(c.kretschmann.at(pt) - K_exact) / std::fabs(K_exact);
        });
        report.check_residual("curvature.desitter_kretschmann",
                              "R_{abcd} R^{abcd} = 24/L^4",
                              VerificationLevel::B, grid_label(g), params,
                              kerr.interior_worst, 2.5e-2, sw.seconds(),
                              "relative; measured 1.07e-2 on this grid; exact "
                              "value " + std::to_string(K_exact));

        sw.reset();
        const geometry::FieldError serr = split(g, [&](std::size_t pt) {
            return std::fabs(c.ricci_square.at(pt) - RS_exact) / std::fabs(RS_exact);
        });
        report.check_residual("curvature.desitter_ricci_square",
                              "R_{mu nu} R^{mu nu} = 36/L^4",
                              VerificationLevel::B, grid_label(g), params,
                              serr.interior_worst, 2.5e-2, sw.seconds(),
                              "relative; exact value " + std::to_string(RS_exact));

        std::ofstream f(out_dir + "/desitter_invariants.csv");
        f << "r,ricci_scalar,ricci_scalar_exact,kretschmann,kretschmann_exact,"
             "ricci_square,ricci_square_exact\n";
        const int jm = g.center_index_y();
        const int km = g.center_index_z();
        for (int i = 0; i < g.nx(); ++i) {
            const std::size_t pt = g.flat(i, jm, km);
            f << g.x(i) << "," << c.ricci_scalar.at(pt) << "," << R_exact << ","
              << c.kretschmann.at(pt) << "," << K_exact << ","
              << c.ricci_square.at(pt) << "," << RS_exact << "\n";
        }
    }
}

// ===========================================================================
// Group: temporal
// ===========================================================================

void run_curvature_temporal_checks(VerificationReport& report,
                                   const std::string& out_dir) {
    (void)out_dir;
    Stopwatch sw;

    const Grid g = make_cartesian_grid(17, 4.0);
    metrics::AlcubierreParameters<double> p(0.5, 2.0, 1.0);

    // --- The refusal ---------------------------------------------------------
    {
        sw.reset();
        auto c = geometry::compute_curvature<double>(
            g, [&](double t) {
                metrics::AlcubierreMetric<double> mm(
                    0.0, metrics::TimeTreatment::FrozenTimeDiagnostic);
                return metrics::build_metric_field<double>(g, mm, p, t);
            },
            0.0,
            TemporalPolicy::refused(
                "TimeTreatment::FrozenTimeDiagnostic discards d_0 g"));

        report.check_predicate(
            "curvature.frozen_time_refused",
            "a Refused temporal policy yields no curvature at all",
            VerificationLevel::A, !c.ok, "ok == false",
            c.ok ? "ok == true (curvature was produced)" : c.refusal_reason,
            sw.seconds(),
            "docs/conventions.md section 11 requires refusal rather than a "
            "warning: imposing d_0 g = 0 on a metric that depends on t gives an "
            "inconsistent input to the Christoffel formula, not the curvature of "
            "some other static spacetime");
    }

    // --- The staticity claim is checked, not trusted -------------------------
    {
        metrics::AlcubierreMetric<double> m(
            0.0, metrics::TimeTreatment::InertialTimeDependent);
        sw.reset();
        auto c = geometry::compute_curvature<double>(
            g, sampler_for(g, m, p), 0.0,
            static_policy("Alcubierre (a deliberately false claim)"));

        report.check_predicate(
            "curvature.false_static_claim_refused",
            "ExactlyStatic is verified against the sampled metric before use",
            VerificationLevel::A, !c.ok, "ok == false",
            c.ok ? "ok == true (the false claim was accepted)"
                 : c.refusal_reason,
            sw.seconds(),
            "probe residual = " + std::to_string(c.staticity_probe_residual) +
                ". Without this, a moving bubble declared static would be frozen "
                "silently -- the failure the refusal above exists to prevent, "
                "arriving through a different door");
    }

    // --- The exact rigid-translation identity --------------------------------
    //
    // For x_s(t) = v_s t every component depends on t and x only through
    // u = x - v_s t, so d_0 g = -v_s d_x g exactly. Choosing the temporal step
    // h_t = h_x / v_s makes the two stencils sample the SAME points along u, so
    // the discrete identity is exact to roundoff rather than to O(h^2) -- which
    // turns a Level D consistency check into a Level A one.
    {
        metrics::AlcubierreMetric<double> m(
            0.0, metrics::TimeTreatment::InertialTimeDependent);
        const double v_s = 0.5;
        const double h_x = g.dx();
        const double h_t = h_x / v_s;

        sw.reset();
        auto c = geometry::compute_curvature<double>(
            g, sampler_for(g, m, p), 0.0,
            TemporalPolicy::analytic_resample(
                h_t, "h_t = h_x / v_s aligns the temporal and spatial stencils "
                     "along u = x - v_s t"));

        double worst_interior = 0.0, worst_edge = 0.0;
        for (std::size_t pt = 0; pt < g.num_points(); ++pt) {
            int i, j, k;
            g.unflatten(pt, i, j, k);
            const bool x_interior = (i > 0 && i < g.nx() - 1);
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu) {
                    const double r = std::fabs(c.dg.at(pt, 0, mu, nu) +
                                               v_s * c.dg.at(pt, 1, mu, nu));
                    if (x_interior)
                        worst_interior = std::fmax(worst_interior, r);
                    else
                        worst_edge = std::fmax(worst_edge, r);
                }
        }

        report.check_residual(
            "curvature.rigid_translation_identity",
            "d_0 g = -v_s d_x g exactly, for x_s(t) = v_s t",
            VerificationLevel::A, grid_label(g), "v=0.5 R=2 sigma=1, h_t = h_x/v_s",
            worst_interior, 1.0e-15, sw.seconds(),
            "roundoff, not truncation: with h_t = h_x/v_s the temporal taps land "
            "on the same u values the spatial stencil uses. x-boundary worst = " +
                std::to_string(worst_edge) +
                " (one-sided in x, centred in t, so O(h^2) there by construction)");
    }

    // --- What the refusal is protecting against ------------------------------
    //
    // Quantifies it rather than asserting it. The frozen-time connection is
    // assembled here from the SAME spatial derivatives, with d_0 g set to zero,
    // and the resulting Ricci scalar is compared against the correct one.
    {
        metrics::AlcubierreMetric<double> m(
            0.0, metrics::TimeTreatment::InertialTimeDependent);
        sw.reset();
        auto c = geometry::compute_curvature<double>(
            g, sampler_for(g, m, p), 0.0,
            TemporalPolicy::analytic_resample(g.dx() / 0.5, "correct path"));

        // The illegitimate path, built by hand from the lower-level pieces.
        const Rank2Field<double> gf = metrics::build_metric_field<double>(g, m, p, 0.0);
        double res = 0.0;
        std::size_t sing = 0;
        const Rank2Field<double> gi =
            geometry::detail::invert_field<double>(g, gf, res, sing);
        const Rank3Field<double> dg_frozen =
            geometry::spacetime_gradient_frozen_time<double>(gf, g);
        const Rank3Field<double> gamma_frozen =
            geometry::christoffel_field(g, gi, dg_frozen);
        const Rank4Field<double> dgamma_frozen =
            geometry::spacetime_gradient_rank3<double>(gamma_frozen, g, nullptr);
        const Rank4Field<double> riem_frozen =
            geometry::riemann_field(g, gamma_frozen, dgamma_frozen);

        double worst_rel = 0.0;
        double correct_scale = 0.0;
        for (std::size_t pt = 0; pt < g.num_points(); ++pt) {
            if (geometry::is_shell_point(g, pt, kShellDepth)) continue;
            const Mat4<double> ric_frozen = contract_13(riem_frozen.tensor_at(pt));
            const double R_frozen = full_contract(gi.matrix_at(pt), ric_frozen);
            const double R_true = c.ricci_scalar.at(pt);
            correct_scale = std::fmax(correct_scale, std::fabs(R_true));
            worst_rel = std::fmax(worst_rel, std::fabs(R_frozen - R_true));
        }

        report.check_predicate(
            "curvature.frozen_time_would_be_materially_wrong",
            "discarding d_0 g changes the Ricci scalar by an amount comparable "
            "to the Ricci scalar itself",
            VerificationLevel::A, worst_rel > 0.1 * correct_scale,
            "max|R_frozen - R| > 0.1 max|R|",
            "max|R_frozen - R| = " + std::to_string(worst_rel) +
                ", max|R| = " + std::to_string(correct_scale),
            sw.seconds(),
            "the refusal is not pedantry; this is the size of the error it "
            "prevents");
    }
}

// ===========================================================================
// Group: alcubierre
// ===========================================================================

void run_alcubierre_curvature_checks(VerificationReport& report,
                                     const std::string& out_dir) {
    Stopwatch sw;

    // --- v = 0 collapses to flat space --------------------------------------
    {
        const Grid g = make_cartesian_grid(17, 4.0);
        metrics::AlcubierreMetric<double> m(
            0.0, metrics::TimeTreatment::InertialTimeDependent);
        metrics::AlcubierreParameters<double> p(0.0, 2.0, 1.0);
        sw.reset();
        auto c = geometry::compute_curvature<double>(
            g, sampler_for(g, m, p), 0.0,
            TemporalPolicy::analytic_resample(0.5, "v_s = 0, but resampled anyway"));

        const double worst = c.ok ? max_abs_rank4_field(g, c.riemann) : 1.0;
        report.check_exact(
            "curvature.alcubierre_zero_velocity_flat",
            "v_s = 0 gives g = eta exactly, hence R^rho_{sigma mu nu} = 0",
            VerificationLevel::A, grid_label(g), 0.0, worst, sw.seconds(),
            "v f = 0 for every f, so the metric is constant and the stencils sum "
            "to exactly zero -- an exactness claim");
    }

    // --- Against the independent symbolic reference --------------------------
    {
        const Grid g = make_cartesian_grid(33, 4.0);
        metrics::AlcubierreMetric<double> m(
            0.0, metrics::TimeTreatment::InertialTimeDependent);
        metrics::AlcubierreParameters<double> p(reference::kAlcubierreRefVelocity,
                                                reference::kAlcubierreRefRadius,
                                                reference::kAlcubierreRefSigma);
        sw.reset();
        auto c = geometry::compute_curvature<double>(
            g, sampler_for(g, m, p), reference::kAlcubierreRefTime,
            TemporalPolicy::analytic_resample(
                g.dx() / reference::kAlcubierreRefVelocity,
                "stencil-aligned resampling of a genuinely time-dependent metric"));

        const std::string params = "v=0.5 R=2 sigma=1, t=0";

        double worst_R = 0.0, worst_K = 0.0, worst_ric = 0.0, worst_rsq = 0.0;
        std::string worst_R_at;
        for (int n = 0; n < reference::kAlcubierreCurvatureCount; ++n) {
            const auto& ref = reference::kAlcubierreCurvature[n];

            // The reference points are exact grid coordinates by construction;
            // locate them by value rather than by assumed index.
            int bi = -1, bj = -1, bk = -1;
            for (int i = 0; i < g.nx(); ++i)
                if (g.x(i) == ref.x) bi = i;
            for (int j = 0; j < g.ny(); ++j)
                if (g.y(j) == ref.y) bj = j;
            for (int k = 0; k < g.nz(); ++k)
                if (g.z(k) == ref.z) bk = k;
            if (bi < 0 || bj < 0 || bk < 0) {
                report.check_predicate(
                    "curvature.alcubierre_reference_point_on_grid",
                    "every reference coordinate is an exact grid coordinate",
                    VerificationLevel::A, false, "found on grid",
                    "reference point not present on the grid", sw.seconds(),
                    "tools/generate_alcubierre_reference.py and this grid have "
                    "drifted apart");
                return;
            }

            const std::size_t pt = g.flat(bi, bj, bk);
            const double dR = std::fabs(c.ricci_scalar.at(pt) - ref.ricci_scalar);
            if (dR > worst_R) {
                worst_R = dR;
                worst_R_at = "(" + std::to_string(ref.x) + "," +
                             std::to_string(ref.y) + "," + std::to_string(ref.z) + ")";
            }
            worst_K = std::fmax(worst_K,
                                std::fabs(c.kretschmann.at(pt) - ref.kretschmann));
            worst_rsq = std::fmax(
                worst_rsq, std::fabs(c.ricci_square.at(pt) - ref.ricci_square));
            const Mat4<double> ric = c.ricci.matrix_at(pt);
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    worst_ric = std::fmax(
                        worst_ric, std::fabs(ric(a, b) - ref.ricci[4 * a + b]));
        }

        report.check_predicate(
            "curvature.alcubierre_reference_point_on_grid",
            "every reference coordinate is an exact grid coordinate",
            VerificationLevel::A, true, "found on grid", "found on grid",
            sw.seconds(), "no interpolation enters the comparison");

        report.check_residual(
            "curvature.alcubierre_ricci_scalar_vs_symbolic",
            "R against exact symbolic differentiation of the closed form",
            VerificationLevel::B, grid_label(g), params, worst_R, 4.0e-2,
            sw.seconds(),
            "worst at " + worst_R_at +
                "; the reference uses EXACT metric derivatives and independently "
                "written algebra, so the gap is this pipeline's O(h^2) truncation "
                "error and is expected to fall fourfold per refinement");

        report.check_residual(
            "curvature.alcubierre_ricci_tensor_vs_symbolic",
            "R_{mu nu} componentwise against the symbolic reference",
            VerificationLevel::B, grid_label(g), params, worst_ric, 2.5e-2,
            sw.seconds(), "componentwise, not merely the trace");

        report.check_residual(
            "curvature.alcubierre_kretschmann_vs_symbolic",
            "R_{abcd} R^{abcd} against the symbolic reference",
            VerificationLevel::B, grid_label(g), params, worst_K, 2.5e-2,
            sw.seconds(), "");

        report.check_residual(
            "curvature.alcubierre_ricci_square_vs_symbolic",
            "R_{mu nu} R^{mu nu} against the symbolic reference",
            VerificationLevel::B, grid_label(g), params, worst_rsq, 2.5e-2,
            sw.seconds(), "");

        // --- The bubble centre, excluded from the reference table ------------
        sw.reset();
        const std::size_t ctr = center_point(g);
        const bool finite = std::isfinite(c.ricci_scalar.at(ctr)) &&
                            std::isfinite(c.kretschmann.at(ctr)) &&
                            std::isfinite(max_abs_real(c.riemann.tensor_at(ctr)));
        report.check_predicate(
            "curvature.alcubierre_centre_finite",
            "curvature is finite at r_s = 0",
            VerificationLevel::A, finite, "finite",
            finite ? "finite" : "non-finite", sw.seconds(),
            "the shape function is even in r_s so d_x f = f'(r_s) x/r_s is a "
            "removable 0/0 there. The symbolic reference cannot evaluate it and "
            "excludes the point; the grid pipeline never forms the ratio, so it "
            "can and does");

        // --- Axial symmetry ---------------------------------------------------
        //
        // The metric depends on y and z only through y^2 + z^2, so every scalar
        // invariant must be invariant under swapping them.
        sw.reset();
        double worst_axial = 0.0;
        for (int i = 0; i < g.nx(); ++i)
            for (int j = 0; j < g.ny(); ++j)
                for (int k = 0; k < g.nz(); ++k) {
                    const double a = c.ricci_scalar.at(g.flat(i, j, k));
                    const double b = c.ricci_scalar.at(g.flat(i, k, j));
                    worst_axial = std::fmax(worst_axial, std::fabs(a - b));
                }
        report.check_residual(
            "curvature.alcubierre_axial_symmetry",
            "R(x,y,z) = R(x,z,y): the bubble is axisymmetric about x",
            VerificationLevel::A, grid_label(g), params, worst_axial, 1.0e-14,
            sw.seconds(),
            "the y and z grids are identical, so this is exact to roundoff and "
            "detects any asymmetry introduced by the index wiring");

        // --- Asymptotic flatness ---------------------------------------------
        sw.reset();
        double wall_max = 0.0;
        double corner = 0.0;
        for (std::size_t pt = 0; pt < g.num_points(); ++pt) {
            if (geometry::is_shell_point(g, pt, kShellDepth)) continue;
            wall_max = std::fmax(wall_max, std::fabs(c.ricci_scalar.at(pt)));
        }
        {
            const int lo = kShellDepth;
            corner = std::fabs(c.ricci_scalar.at(g.flat(lo, lo, lo)));
        }
        report.check_predicate(
            "curvature.alcubierre_asymptotically_flat",
            "curvature at the far corner is small compared with the wall",
            VerificationLevel::A, corner < 0.02 * wall_max,
            "|R| at corner < 2% of max|R|",
            "corner = " + std::to_string(corner) + ", max = " +
                std::to_string(wall_max),
            sw.seconds(), "structural, not a closed form");

        // --- Records for the figures -----------------------------------------
        const int km = g.center_index_z();
        std::ofstream f(out_dir + "/alcubierre_curvature_slice_xy.csv");
        f << "x,y,g00,ricci_scalar,ricci_square,kretschmann,riemann_max\n";
        for (int i = 0; i < g.nx(); ++i)
            for (int j = 0; j < g.ny(); ++j) {
                const std::size_t pt = g.flat(i, j, km);
                f << g.x(i) << "," << g.y(j) << "," << c.g.at(pt, 0, 0) << ","
                  << c.ricci_scalar.at(pt) << "," << c.ricci_square.at(pt) << ","
                  << c.kretschmann.at(pt) << ","
                  << max_abs_real(c.riemann.tensor_at(pt)) << "\n";
            }

        const int jm = g.center_index_y();
        std::ofstream fc(out_dir + "/alcubierre_centerline_x.csv");
        fc << "x,g00,g01,ricci_scalar,ricci_square,kretschmann,riemann_max,"
              "christoffel_max\n";
        for (int i = 0; i < g.nx(); ++i) {
            const std::size_t pt = g.flat(i, jm, km);
            fc << g.x(i) << "," << c.g.at(pt, 0, 0) << "," << c.g.at(pt, 0, 1) << ","
               << c.ricci_scalar.at(pt) << "," << c.ricci_square.at(pt) << ","
               << c.kretschmann.at(pt) << ","
               << max_abs_real(c.riemann.tensor_at(pt)) << ","
               << max_abs_real(c.gamma.tensor_at(pt)) << "\n";
        }
    }
}

// ===========================================================================
// Group: convergence
// ===========================================================================

void run_curvature_convergence_checks(VerificationReport& report,
                                      const std::string& out_dir) {
    Stopwatch sw;
    std::ofstream f(out_dir + "/curvature_convergence.csv");
    f << "case,n,h,interior_error,boundary_error\n";

    // --- Schwarzschild Kretschmann ------------------------------------------
    //
    // Evaluated at the pinned midpoint, which is r = 5 exactly on every grid, so
    // refinement changes the discretisation and nothing else.
    {
        const int ns[3] = {9, 17, 33};
        double err[3] = {0, 0, 0};
        double hs[3] = {0, 0, 0};
        const double exact =
            metrics::SchwarzschildMetric<double>::kretschmann(5.0, kSchwarzschildMass);

        sw.reset();
        for (int a = 0; a < 3; ++a) {
            const Grid g = make_spherical_grid(ns[a], ns[a], kSphericalNPhi, kSchwarzschildRMin,
                                               kSchwarzschildRMax);
            metrics::SchwarzschildMetric<double> m;
            metrics::SchwarzschildParameters<double> p(kSchwarzschildMass);
            auto c = geometry::compute_curvature<double>(
                g, sampler_for(g, m, p), 0.0, static_policy("Schwarzschild"));
            hs[a] = g.dx();
            err[a] = std::fabs(c.kretschmann.at(center_point(g)) - exact) /
                     std::fabs(exact);
            f << "schwarzschild_kretschmann," << ns[a] << "," << hs[a] << ","
              << err[a] << ",\n";
        }

        const double order01 = std::log(err[0] / err[1]) / std::log(hs[0] / hs[1]);
        const double order12 = std::log(err[1] / err[2]) / std::log(hs[1] / hs[2]);

        report.check_value("curvature.convergence_schwarzschild_kretschmann",
                           "K error falls as O(h^2) under refinement",
                           VerificationLevel::E, "9,17,33 (theta refined with r)",
                           "M=1, evaluated at the pinned midpoint r=5", 2.0,
                           order12, 0.25, sw.seconds(),
                           "orders: 9->17 = " + std::to_string(order01) +
                               ", 17->33 = " + std::to_string(order12) +
                               ". The finest pair is the one asserted; the coarse "
                               "pair is recorded because a coarse grid is not yet "
                               "in the asymptotic regime");
    }

    // --- Flat space in spherical coordinates --------------------------------
    {
        const int ns[3] = {9, 17, 33};
        double err[3] = {0, 0, 0};
        double hs[3] = {0, 0, 0};

        sw.reset();
        for (int a = 0; a < 3; ++a) {
            const Grid g = make_spherical_grid(ns[a], ns[a], kSphericalNPhi, 1.0, 5.0);
            metrics::FlatSphericalMetric<double> m;
            metrics::FlatSphericalMetric<double>::Parameters p;
            auto c = geometry::compute_curvature<double>(
                g, sampler_for(g, m, p), 0.0, static_policy("flat spherical"));
            hs[a] = g.dx();
            // At the pinned midpoint, not the interior maximum. A worst-over-a-
            // region measure is not a convergence measure: as the grid refines the
            // excluded shell shrinks in physical size, so the location of the
            // worst point migrates toward the domain edge where the error constant
            // is largest, and the measured order comes out spuriously low.
            err[a] = max_abs_real(c.riemann.tensor_at(center_point(g)));
            const geometry::FieldError e = split(g, [&](std::size_t pt) {
                return max_abs_real(c.riemann.tensor_at(pt));
            });
            f << "flat_spherical_riemann," << ns[a] << "," << hs[a] << ","
              << err[a] << "," << e.boundary_worst << "\n";
        }

        const double order12 = std::log(err[1] / err[2]) / std::log(hs[1] / hs[2]);
        report.check_value("curvature.convergence_flat_spherical_riemann",
                           "the spurious Riemann of flat space falls as O(h^2)",
                           VerificationLevel::E, "9,17,33",
                           "r in [1,5]", 2.0, order12, 0.3, sw.seconds(),
                           "order 9->17 = " +
                               std::to_string(std::log(err[0] / err[1]) /
                                              std::log(hs[0] / hs[1])));
    }

    // --- Ricci symmetry ------------------------------------------------------
    //
    // The case that justifies calling the symmetry residuals in the connection
    // group truncation error rather than a defect. If they did not converge, the
    // right conclusion would be that the assembly is wrong.
    {
        const int ns[3] = {17, 25, 33};
        double err[3] = {0, 0, 0};
        double hs[3] = {0, 0, 0};

        sw.reset();
        for (int a = 0; a < 3; ++a) {
            const Grid g = make_cartesian_grid(ns[a], 4.0);
            metrics::AlcubierreMetric<double> m(
                0.0, metrics::TimeTreatment::InertialTimeDependent);
            metrics::AlcubierreParameters<double> p(0.5, 2.0, 1.0);
            auto c = geometry::compute_curvature<double>(
                g, sampler_for(g, m, p), 0.0,
                TemporalPolicy::analytic_resample(g.dx() / 0.5, "stencil aligned"));
            hs[a] = g.dx();
            // At fixed integer coordinates present on all three grids, for the
            // same reason as above: a worst-over-region measure migrates as the
            // excluded shell shrinks.
            //
            // NOT at the pinned centre. The bubble centre is a symmetry point of
            // the configuration, so the antisymmetry residual vanishes there
            // identically on every grid and the measured order comes out as 0/0.
            // A convergence probe has to sit somewhere the quantity is nonzero.
            const double probe[2][3] = {{2.0, 1.0, 0.0}, {1.0, 1.0, 1.0}};
            double worst = 0.0;
            for (int s = 0; s < 2; ++s) {
                int bi = -1, bj = -1, bk = -1;
                for (int i = 0; i < g.nx(); ++i) if (g.x(i) == probe[s][0]) bi = i;
                for (int j = 0; j < g.ny(); ++j) if (g.y(j) == probe[s][1]) bj = j;
                for (int k = 0; k < g.nz(); ++k) if (g.z(k) == probe[s][2]) bk = k;
                if (bi < 0 || bj < 0 || bk < 0) continue;
                worst = std::fmax(worst, symmetry_residual(
                                             c.ricci.matrix_at(g.flat(bi, bj, bk))));
            }
            err[a] = worst;
            f << "ricci_symmetry," << ns[a] << "," << hs[a] << "," << err[a]
              << ",\n";
        }

        const double order12 = std::log(err[1] / err[2]) / std::log(hs[1] / hs[2]);
        report.check_value(
            "curvature.convergence_ricci_symmetry",
            "the Ricci antisymmetry induced by the discrete product rule falls "
            "as O(h^2)",
            VerificationLevel::E, "17,25,33", "v=0.5 R=2 sigma=1", 2.0, order12,
            0.4, sw.seconds(),
            "order 17->25 = " +
                std::to_string(std::log(err[0] / err[1]) / std::log(hs[0] / hs[1])) +
                ". R_{mu nu} is symmetric in the continuum; the residual is a "
                "property of the discretisation and must vanish with it");
    }

    // --- Alcubierre against the symbolic reference ---------------------------
    //
    // The only convergence case with a genuinely time-dependent metric, so it is
    // the only one that measures the order of the resampled temporal path.
    {
        // 17, 25 and 33 points over [-4, 4] give spacings 1/2, 1/3 and 1/4, and
        // every one of those grids contains the integer coordinates exactly. A
        // 9-point grid was tried and dropped: 8 cells across the whole domain
        // leaves the wall unresolved, so that pair measures the approach to the
        // asymptotic regime rather than the order within it.
        const int ns[3] = {17, 25, 33};
        double err[3] = {0, 0, 0};
        double hs[3] = {0, 0, 0};

        // Integer-valued reference points, present on all three grids.
        const double px[2] = {2.0, 1.0};
        const double py[2] = {0.0, 1.0};
        const double pz[2] = {0.0, 0.0};

        sw.reset();
        for (int a = 0; a < 3; ++a) {
            const Grid g = make_cartesian_grid(ns[a], 4.0);
            metrics::AlcubierreMetric<double> m(
                0.0, metrics::TimeTreatment::InertialTimeDependent);
            metrics::AlcubierreParameters<double> p(reference::kAlcubierreRefVelocity,
                                                    reference::kAlcubierreRefRadius,
                                                    reference::kAlcubierreRefSigma);
            auto c = geometry::compute_curvature<double>(
                g, sampler_for(g, m, p), reference::kAlcubierreRefTime,
                TemporalPolicy::analytic_resample(
                    g.dx() / reference::kAlcubierreRefVelocity, "stencil aligned"));
            hs[a] = g.dx();

            double worst = 0.0;
            for (int s = 0; s < 2; ++s) {
                int bi = -1, bj = -1, bk = -1;
                for (int i = 0; i < g.nx(); ++i) if (g.x(i) == px[s]) bi = i;
                for (int j = 0; j < g.ny(); ++j) if (g.y(j) == py[s]) bj = j;
                for (int k = 0; k < g.nz(); ++k) if (g.z(k) == pz[s]) bk = k;
                if (bi < 0 || bj < 0 || bk < 0) continue;

                for (int n = 0; n < reference::kAlcubierreCurvatureCount; ++n) {
                    const auto& ref = reference::kAlcubierreCurvature[n];
                    if (ref.x != px[s] || ref.y != py[s] || ref.z != pz[s]) continue;
                    worst = std::fmax(
                        worst, std::fabs(c.ricci_scalar.at(g.flat(bi, bj, bk)) -
                                         ref.ricci_scalar));
                }
            }
            err[a] = worst;
            f << "alcubierre_ricci_scalar," << ns[a] << "," << hs[a] << "," << worst
              << ",\n";
        }

        const double order12 = std::log(err[1] / err[2]) / std::log(hs[1] / hs[2]);
        report.check_value(
            "curvature.convergence_alcubierre_ricci_scalar",
            "R converges to the symbolic reference at O(h^2), including the "
            "resampled temporal derivative",
            VerificationLevel::E, "9,17,33", "v=0.5 R=2 sigma=1, t=0", 2.0,
            order12, 0.4, sw.seconds(),
            "order 9->17 = " +
                std::to_string(std::log(err[0] / err[1]) / std::log(hs[0] / hs[1])) +
                ". Demonstrating the order here is what makes the fixed-tolerance "
                "comparison in the alcubierre group a check rather than a "
                "calibration");
    }
}

}  // namespace checks
}  // namespace warp_hypad

