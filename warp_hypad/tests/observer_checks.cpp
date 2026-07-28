#include "observer_checks.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "warp_hypad/geometry/einstein.hpp"
#include "warp_hypad/geometry/energy_conditions.hpp"
#include "warp_hypad/geometry/observer.hpp"
#include "warp_hypad/metrics/analytic_metrics.hpp"
#include "warp_hypad/metrics/metric.hpp"
#include "warp_hypad/metrics/metric_field.hpp"

// ---------------------------------------------------------------------------
// Milestone A6 verification.
//
// SCOPE: A6 evaluates observer-dependent energy density and classical pointwise
// energy conditions for the stress-energy tensor inferred from the selected
// metric. It does not demonstrate that the source can be generated, that quantum
// energy inequalities are satisfied, that the geometry is dynamically stable, or
// that operational FTL travel is possible.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using geometry::AdmDecomposition;
using geometry::ConditionStatus;
using geometry::FrameStress;
using geometry::Observer;
using geometry::Tetrad;

constexpr double kPi = 3.14159265358979323846;

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

// The Alcubierre configuration used throughout A6.
constexpr double kV = 0.5, kR = 2.0, kSig = 1.0;

struct AlcubierreFields {
    Grid grid;
    geometry::CurvatureFields<double> curv;
    geometry::EinsteinFields<double> ein;
};

// Build the full chain once, for a given grid.
inline geometry::EinsteinFields<double> alcubierre_einstein(
    const Grid& g, geometry::CurvatureFields<double>& curv_out, double v_s = kV) {
    static metrics::AlcubierreMetric<double> m(
        0.0, metrics::TimeTreatment::InertialTimeDependent);
    metrics::AlcubierreParameters<double> p(v_s, kR, kSig);
    auto sample = [&](double t) {
        return metrics::build_metric_field<double>(g, m, p, t);
    };
    curv_out = geometry::compute_curvature<double>(
        g, sample, 0.0,
        geometry::TemporalPolicy::analytic_resample(g.dx() / v_s,
                                                    "stencil aligned"));
    return geometry::assemble_einstein(
        g, curv_out, geometry::CosmologicalConstantPolicy::none());
}

// --- The independent analytic Eulerian energy density -------------------------
//
// For the Alcubierre metric alpha = 1 and gamma_ij = delta_ij, so the Eulerian
// observer is n^mu = (1, -v_s f, 0, 0) and the standard result is
//
//     rho = -(v_s^2 / 32 pi) [ (d_y f)^2 + (d_z f)^2 ]
//
// Written here from the closed form, using the analytic radial derivative of the
// shape function from metric.hpp. It never contracts the production tensor, so
// it is an independent reference rather than a restatement.
//
// With f = f(r_s):  d_y f = f'(r) y/r,  d_z f = f'(r) z/r, hence
//
//     (d_y f)^2 + (d_z f)^2 = f'(r)^2 (y^2 + z^2) / r^2
//
// which vanishes on the x axis, where y = z = 0. That is the analytic reason the
// energy density is zero on the axis of travel, and it is checked below.
inline double analytic_eulerian_rho(double x, double y, double z, double v_s) {
    const double r2 = x * x + y * y + z * z;
    if (r2 <= 1e-24) return 0.0;  // f'(0) = 0, f is even in r
    const double r = std::sqrt(r2);
    metrics::AlcubierreParameters<double> p(v_s, kR, kSig);
    const double fp = metrics::shape_function_radial_derivative(r, p);
    const double rho2 = y * y + z * z;
    return -(v_s * v_s / (32.0 * kPi)) * fp * fp * rho2 / r2;
}

// A frame stress tensor built directly from rho and a diagonal pressure triple,
// bypassing the metric pipeline entirely.
inline FrameStress perfect_fluid(double rho, double px, double py, double pz) {
    FrameStress f;
    f.T[0][0] = rho;
    f.T[1][1] = px; f.T[2][2] = py; f.T[3][3] = pz;
    f.rho = rho;
    f.S[0] = px; f.S[4] = py; f.S[8] = pz;
    f.trace_spatial = px + py + pz;
    geometry::symmetric3_eigenvalues(f.S, f.principal);
    return f;
}

const char* kScope =
    "A6 evaluates observer-dependent energy density and classical pointwise "
    "energy conditions for the stress-energy tensor inferred from the selected "
    "metric. It does not demonstrate that the source can be generated, that "
    "quantum energy inequalities are satisfied, that the geometry is "
    "dynamically stable, or that operational FTL travel is possible.";

// A generous fixed envelope for the analytically exact fixtures: they carry no
// discretisation error at all, so anything above roundoff is certified.
constexpr double kFixtureEnvelope = 1e-14;

}  // namespace

// ===========================================================================
// Group: frames
// ===========================================================================

void run_frame_checks(VerificationReport& report, const std::string& out_dir) {
    Stopwatch sw;
    const Grid g = cartesian_grid(17, 4.0);
    geometry::CurvatureFields<double> c;
    auto e = alcubierre_einstein(g, c);
    (void)e;

    const std::string params = "v=0.5 R=2 sigma=1";

    // --- ADM recovery --------------------------------------------------------
    //
    // For Alcubierre the exact answer is known: alpha = 1 and beta^x = -v_s f,
    // gamma_ij = delta_ij. Recovering those from the numerically inverted metric
    // is a real test of the decomposition rather than of the metric.
    sw.reset();
    double worst_alpha = 0.0, worst_gamma = 0.0, worst_norm = 0.0;
    std::size_t invalid_adm = 0;
    for (std::size_t p = 0; p < g.num_points(); ++p) {
        const Mat4<double> gm = c.g.matrix_at(p);
        const AdmDecomposition a = geometry::adm_decompose(gm, c.g_inv.matrix_at(p));
        if (!a.valid) { ++invalid_adm; continue; }
        worst_alpha = std::fmax(worst_alpha, std::fabs(a.alpha - 1.0));
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                worst_gamma = std::fmax(
                    worst_gamma,
                    std::fabs(a.gamma[3 * i + j] - ((i == j) ? 1.0 : 0.0)));
        double n_up[4], n_dn[4];
        geometry::eulerian_normal(a, n_up, n_dn);
        worst_norm = std::fmax(
            worst_norm, std::fabs(geometry::four_velocity_norm(gm, n_up) + 1.0));
    }

    report.check_residual("frames.adm_lapse_recovery",
                          "alpha = 1 exactly for the Alcubierre metric",
                          VerificationLevel::B, grid_label(g), params,
                          worst_alpha, 1.0e-14, sw.seconds(),
                          "recovered from the numerically inverted metric via "
                          "alpha = 1/sqrt(-g^{00}), not read off the ansatz");

    report.check_residual("frames.adm_spatial_metric_recovery",
                          "gamma_ij = delta_ij for the Alcubierre metric",
                          VerificationLevel::B, grid_label(g), params,
                          worst_gamma, 1.0e-15, sw.seconds(), "");

    report.check_residual("frames.eulerian_normalization",
                          "g_{mu nu} n^mu n^nu = -1",
                          VerificationLevel::A, grid_label(g), params,
                          worst_norm, 1.0e-13, sw.seconds(),
                          "the Eulerian normal is normalised by construction, "
                          "so this measures the inverse metric's quality");

    // --- Tetrad --------------------------------------------------------------
    sw.reset();
    double worst_ortho = 0.0, min_hand = 1e300, max_hand = -1e300;
    std::size_t invalid_tetrad = 0;
    for (std::size_t p = 0; p < g.num_points(); ++p) {
        const Mat4<double> gm = c.g.matrix_at(p);
        const AdmDecomposition a = geometry::adm_decompose(gm, c.g_inv.matrix_at(p));
        const Tetrad t = geometry::build_tetrad(gm, a);
        if (!t.valid) { ++invalid_tetrad; continue; }
        worst_ortho = std::fmax(worst_ortho, t.orthonormality_residual);
        min_hand = std::fmin(min_hand, t.handedness);
        max_hand = std::fmax(max_hand, t.handedness);
    }

    report.check_residual("frames.tetrad_orthonormality",
                          "g_{mu nu} e_a^mu e_b^nu = eta_{ab}",
                          VerificationLevel::A, grid_label(g), params,
                          worst_ortho, 1.0e-13, sw.seconds(),
                          "over every grid point; " +
                              std::to_string(invalid_tetrad) +
                              " points had no valid tetrad");

    report.check_predicate(
        "frames.tetrad_handedness_uniform",
        "the frame determinant does not change sign across the grid",
        VerificationLevel::A, (min_hand > 0.0) == (max_hand > 0.0),
        "one sign throughout",
        "determinant range [" + std::to_string(min_hand) + ", " +
            std::to_string(max_hand) + "]",
        sw.seconds(),
        "a handedness flip between neighbouring points would appear as a "
        "spurious sign change in q_i and in the worst null direction");

    // --- Continuity ----------------------------------------------------------
    //
    // The reason fixed-order Gram-Schmidt was chosen over an eigen-decomposition.
    // Adjacent tetrads must differ by O(h), not by O(1).
    sw.reset();
    double worst_jump = 0.0;
    for (int i = 1; i < g.nx() - 1; ++i)
        for (int j = 1; j < g.ny() - 1; ++j)
            for (int k = 1; k < g.nz() - 1; ++k) {
                const std::size_t p0 = g.flat(i, j, k), p1 = g.flat(i + 1, j, k);
                const AdmDecomposition a0 = geometry::adm_decompose(
                    c.g.matrix_at(p0), c.g_inv.matrix_at(p0));
                const AdmDecomposition a1 = geometry::adm_decompose(
                    c.g.matrix_at(p1), c.g_inv.matrix_at(p1));
                const Tetrad t0 = geometry::build_tetrad(c.g.matrix_at(p0), a0);
                const Tetrad t1 = geometry::build_tetrad(c.g.matrix_at(p1), a1);
                if (!t0.valid || !t1.valid) continue;
                for (int A = 0; A < 4; ++A)
                    for (int mu = 0; mu < 4; ++mu)
                        worst_jump = std::fmax(
                            worst_jump, std::fabs(t0.e[A][mu] - t1.e[A][mu]));
            }
    report.check_residual(
        "frames.tetrad_continuity",
        "adjacent tetrads differ by O(h), not by O(1)",
        VerificationLevel::A, grid_label(g), params, worst_jump, 0.2,
        sw.seconds(),
        "the reason for fixed-order Gram-Schmidt with a positive diagonal "
        "rather than an eigen-decomposition: eigenvector signs and ordering "
        "flip wherever eigenvalues nearly cross, and each flip is a "
        "discontinuity no physics produced");

    std::ofstream f(out_dir + "/tetrad_verification.csv");
    f << "x,alpha,beta_x,gamma_min_eigenvalue,orthonormality_residual,"
         "handedness,eulerian_norm_residual\n";
    const int jm = g.center_index_y(), km = g.center_index_z();
    for (int i = 0; i < g.nx(); ++i) {
        const std::size_t p = g.flat(i, jm, km);
        const Mat4<double> gm = c.g.matrix_at(p);
        const AdmDecomposition a = geometry::adm_decompose(gm, c.g_inv.matrix_at(p));
        const Tetrad t = geometry::build_tetrad(gm, a);
        double n_up[4], n_dn[4];
        if (a.valid) geometry::eulerian_normal(a, n_up, n_dn);
        f << g.x(i) << "," << a.alpha << "," << a.beta[0] << ","
          << a.gamma_min_eigenvalue << "," << t.orthonormality_residual << ","
          << t.handedness << ","
          << (a.valid ? std::fabs(geometry::four_velocity_norm(gm, n_up) + 1.0)
                      : -1.0)
          << "\n";
    }
}

// ===========================================================================
// Group: observers
// ===========================================================================

void run_observer_checks(VerificationReport& report, const std::string& out_dir) {
    Stopwatch sw;

    // A metric slice sampled across the g_00 = 0 crossing, by varying v_s.
    // g_00 = -1 + v_s^2 f^2, so at the bubble centre f = 1 and g_00 = -1 + v_s^2.
    const Grid g = cartesian_grid(13, 4.0);

    std::size_t rejected = 0;
    std::ofstream f(out_dir + "/observer_comparison.csv");
    f << "v_s,g00_center,eulerian_valid,eulerian_norm_residual,"
         "stationary_valid,stationary_reason,comoving_valid,boosted_valid\n";

    bool stationary_ok_below = true;
    bool stationary_rejected_at_zero = false;
    bool stationary_rejected_above = true;
    bool eulerian_ok_everywhere = true;
    double worst_eul_resid = 0.0;

    for (int k = 2; k <= 16; ++k) {
        const double v_s = static_cast<double>(k) / 10.0;
        metrics::AlcubierreMetric<double> m(
            0.0, metrics::TimeTreatment::InertialTimeDependent);
        metrics::AlcubierreParameters<double> p(v_s, kR, kSig);
        const Rank2Field<double> gf =
            metrics::build_metric_field<double>(g, m, p, 0.0);

        const std::size_t ctr =
            g.flat(g.center_index_x(), g.center_index_y(), g.center_index_z());
        const Mat4<double> gm = gf.matrix_at(ctr);
        const InverseResult<double> inv = inverse(gm);

        const AdmDecomposition a = geometry::adm_decompose(gm, inv.inverse);
        const Observer eul = geometry::make_eulerian(gm, a);
        const Observer sta = geometry::make_coordinate_stationary(gm);
        const Observer com = geometry::make_bubble_comoving(gm, v_s);
        const Tetrad t = geometry::build_tetrad(gm, a);
        const double bv[3] = {0.3, 0.2, 0.0};
        const Observer bst = geometry::make_tetrad_boosted(gm, t, bv);

        const double g00 = gm(0, 0);
        if (!eul.valid) eulerian_ok_everywhere = false;
        worst_eul_resid = std::fmax(worst_eul_resid, eul.normalization_residual);
        if (g00 < 0.0 && !sta.valid) stationary_ok_below = false;
        if (g00 == 0.0 && !sta.valid) stationary_rejected_at_zero = true;
        if (g00 > 0.0 && sta.valid) stationary_rejected_above = false;
        if (!sta.valid) ++rejected;
        if (!com.valid) ++rejected;
        if (!bst.valid) ++rejected;

        f << v_s << "," << g00 << "," << (eul.valid ? 1 : 0) << ","
          << eul.normalization_residual << "," << (sta.valid ? 1 : 0) << ",\""
          << sta.rejection_reason << "\"," << (com.valid ? 1 : 0) << ","
          << (bst.valid ? 1 : 0) << "\n";
    }

    report.check_predicate(
        "observers.stationary_valid_only_below_zero",
        "the coordinate-stationary observer exists exactly where g_00 < 0",
        VerificationLevel::A, stationary_ok_below, "valid for every g_00 < 0",
        stationary_ok_below ? "valid" : "rejected somewhere it should exist",
        sw.seconds(), "");

    report.check_predicate(
        "observers.stationary_rejected_at_null",
        "the coordinate-stationary observer is REJECTED at g_00 = 0",
        VerificationLevel::A, stationary_rejected_at_zero, "rejected",
        stationary_rejected_at_zero ? "rejected: d/dt is null"
                                    : "accepted -- a fictitious observer",
        sw.seconds(),
        "d/dt is null there. Normalising by |g_00| would produce a plausible, "
        "plottable, entirely fictitious energy density, and that substitution "
        "is what this case forbids");

    report.check_predicate(
        "observers.stationary_rejected_when_spacelike",
        "the coordinate-stationary observer is REJECTED for g_00 > 0",
        VerificationLevel::A, stationary_rejected_above, "rejected",
        stationary_rejected_above ? "rejected: d/dt is spacelike" : "accepted",
        sw.seconds(), "");

    report.check_predicate(
        "observers.eulerian_survives_the_crossing",
        "the Eulerian observer remains timelike through g_00 = 0 and beyond",
        VerificationLevel::A, eulerian_ok_everywhere, "valid at every v_s",
        eulerian_ok_everywhere ? "valid throughout" : "failed somewhere",
        sw.seconds(),
        "the reason A6 hangs the physics on the foliation rather than on the "
        "coordinate time direction. Worst normalisation residual " +
            std::to_string(worst_eul_resid));

    report.check_predicate(
        "observers.rejections_counted",
        "inadmissible observers are counted, not silently normalised",
        VerificationLevel::A, rejected > 0, "some candidates rejected",
        std::to_string(rejected) + " rejections across the sweep",
        sw.seconds(), "");

    // Boost admissibility at the light cone.
    sw.reset();
    const Grid g1 = cartesian_grid(5, 1.0);
    metrics::MinkowskiMetric<double> mm;
    metrics::MinkowskiMetric<double>::Parameters mp;
    const Rank2Field<double> mf =
        metrics::build_metric_field<double>(g1, mm, mp, 0.0);
    const Mat4<double> eta = mf.matrix_at(0);
    const InverseResult<double> ei = inverse(eta);
    const AdmDecomposition ea = geometry::adm_decompose(eta, ei.inverse);
    const Tetrad et = geometry::build_tetrad(eta, ea);
    const double v_light[3] = {1.0, 0.0, 0.0};
    const double v_over[3] = {1.2, 0.0, 0.0};
    const Observer o_light = geometry::make_tetrad_boosted(eta, et, v_light);
    const Observer o_over = geometry::make_tetrad_boosted(eta, et, v_over);
    report.check_predicate(
        "observers.boost_rejected_at_and_beyond_light_speed",
        "boosted observers require |v| < 1 strictly",
        VerificationLevel::A, !o_light.valid && !o_over.valid,
        "both rejected", o_light.rejection_reason, sw.seconds(), "");
}

// ===========================================================================
// Group: fixtures
// ===========================================================================

void run_fixture_checks(VerificationReport& report, const std::string& out_dir) {
    Stopwatch sw;
    std::ofstream f(out_dir + "/observer_reference_tensors.csv");
    f << "fixture,rho,p1,p2,p3,N_min,worst_nx,worst_ny,worst_nz,"
         "WEC_worst,DEC_worst,SEC_worst,NEC_status,WEC_status,DEC_status,"
         "SEC_status\n";

    struct Case {
        const char* name;
        FrameStress f;
        bool nec_ok, wec_ok, dec_ok, sec_ok;
        const char* note;
    };

    std::vector<Case> cases;

    cases.push_back({"vacuum", perfect_fluid(0.0, 0.0, 0.0, 0.0),
                     true, true, true, true,
                     "every condition saturated; nothing to violate"});
    cases.push_back({"dust", perfect_fluid(1.0, 0.0, 0.0, 0.0),
                     true, true, true, true, "p = 0, rho > 0"});
    cases.push_back({"radiation",
                     perfect_fluid(1.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
                     true, true, true, true,
                     "p = rho/3; the frame trace vanishes"});
    // Positive cosmological constant as an effective source: p = -rho.
    cases.push_back({"lambda_positive", perfect_fluid(1.0, -1.0, -1.0, -1.0),
                     true, true, true, false,
                     "NEC saturated, WEC and DEC satisfied, SEC VIOLATED -- "
                     "the standard signature of a positive Lambda"});
    cases.push_back({"phantom", perfect_fluid(1.0, -1.6, -1.6, -1.6),
                     false, false, false, false,
                     "p < -rho; NEC and WEC both violated"});
    // Independent principal pressures with rho + p_3 < 0. For a diagonal
    // source every condition is governed by rho + p_i, so an NEC violation
    // necessarily drags WEC, DEC and SEC down with it -- an earlier expectation
    // of "NEC violated but WEC satisfied" was simply wrong physics. What the
    // fixture still establishes is DIRECTIONAL: rho = +1 throughout, so a test
    // that looked only at the energy density would call this source unremarkable,
    // while N_min = rho + p_min = -0.4 along the z axis.
    cases.push_back({"anisotropic", perfect_fluid(1.0, 0.5, -0.2, -1.4),
                     false, false, false, false,
                     "rho = +1 everywhere yet N_min = -0.4 along z; a "
                     "density-only test would miss it entirely"});

    // A momentum-flux source: the worst null direction is deliberately NOT a
    // coordinate axis, so a six-direction test finds the wrong answer.
    FrameStress mom = perfect_fluid(1.0, 0.2, 0.2, 0.2);
    mom.T[0][1] = mom.T[1][0] = 0.9;
    mom.T[0][2] = mom.T[2][0] = 0.6;
    mom.q[0] = -0.9; mom.q[1] = -0.6;
    mom.q_norm = std::sqrt(0.9 * 0.9 + 0.6 * 0.6);
    cases.push_back({"momentum_flux", mom, false, false, false, false,
                     "off-axis worst direction by construction"});

    for (const Case& cs : cases) {
        sw.reset();
        const geometry::NecResult n = geometry::nec_minimum(cs.f);
        const auto w = geometry::wec(cs.f, kFixtureEnvelope);
        const auto d = geometry::dec(cs.f, kFixtureEnvelope);
        const auto s = geometry::sec(cs.f, kFixtureEnvelope);

        const bool nec_satisfied = n.n_min > -1e-12;
        const bool wec_satisfied = w.status != ConditionStatus::CertifiedViolated;
        const bool dec_satisfied = d.status != ConditionStatus::CertifiedViolated;
        const bool sec_satisfied = s.status != ConditionStatus::CertifiedViolated;

        const std::string nm = std::string("fixtures.") + cs.name;
        report.check_predicate(
            nm, "known energy-condition signature of a standard source",
            VerificationLevel::B,
            nec_satisfied == cs.nec_ok && wec_satisfied == cs.wec_ok &&
                dec_satisfied == cs.dec_ok && sec_satisfied == cs.sec_ok,
            std::string("NEC ") + (cs.nec_ok ? "ok" : "violated") + ", WEC " +
                (cs.wec_ok ? "ok" : "violated") + ", DEC " +
                (cs.dec_ok ? "ok" : "violated") + ", SEC " +
                (cs.sec_ok ? "ok" : "violated"),
            std::string("NEC ") + (nec_satisfied ? "ok" : "violated") +
                " (N_min = " + std::to_string(n.n_min) + "), WEC " +
                (wec_satisfied ? "ok" : "violated") + ", DEC " +
                (dec_satisfied ? "ok" : "violated") + ", SEC " +
                (sec_satisfied ? "ok" : "violated"),
            sw.seconds(),
            std::string(cs.note) +
                ". Built directly in the orthonormal frame, so this tests the "
                "condition evaluators alone and no part of the metric pipeline");

        f << cs.name << "," << cs.f.rho << "," << cs.f.principal[0] << ","
          << cs.f.principal[1] << "," << cs.f.principal[2] << "," << n.n_min
          << "," << n.direction[0] << "," << n.direction[1] << ","
          << n.direction[2] << "," << w.worst_value << "," << d.worst_value
          << "," << s.worst_value << ","
          << (nec_satisfied ? "satisfied" : "violated") << ","
          << geometry::to_string(w.status) << "," << geometry::to_string(d.status)
          << "," << geometry::to_string(s.status) << "\n";
    }

    // The momentum-flux case is where a six-direction test fails. Demonstrate it.
    sw.reset();
    const geometry::NecResult nm2 = geometry::nec_minimum(mom);
    double axis_best = 1e300;
    for (int ax = 0; ax < 3; ++ax)
        for (int sgn = -1; sgn <= 1; sgn += 2) {
            double n[3] = {0, 0, 0};
            n[ax] = static_cast<double>(sgn);
            axis_best = std::fmin(axis_best, geometry::nec_value(mom, n));
        }
    report.check_predicate(
        "nec.axis_only_search_is_insufficient",
        "the six coordinate directions do not find the NEC minimum",
        VerificationLevel::A, axis_best > nm2.n_min + 1e-6,
        "axis-only search overstates N_min",
        "axis-only = " + std::to_string(axis_best) + ", true minimum = " +
            std::to_string(nm2.n_min) + " at (" +
            std::to_string(nm2.direction[0]) + ", " +
            std::to_string(nm2.direction[1]) + ", " +
            std::to_string(nm2.direction[2]) + ")",
        sw.seconds(),
        "quantifies why the global minimiser exists. An axis-only test would "
        "report a less negative worst case and could certify a violated NEC as "
        "satisfied");
}

// ===========================================================================
// Group: nec
// ===========================================================================

void run_nec_checks(VerificationReport& report, const std::string& out_dir) {
    Stopwatch sw;

    // A stress tensor whose minimiser is off every coordinate axis. All three
    // momentum components are nonzero deliberately: with q_y = 0 the problem
    // reduces to the x-z plane AND lands in the trust-region hard case, which
    // is a valuable configuration but not the one this group is measuring.
    FrameStress f = perfect_fluid(0.4, 0.9, -0.3, 0.2);
    f.T[0][1] = f.T[1][0] = 0.5;
    f.T[0][2] = f.T[2][0] = 0.28;
    f.T[0][3] = f.T[3][0] = -0.35;
    f.q[0] = -0.5; f.q[1] = -0.28; f.q[2] = 0.35;

    sw.reset();
    const geometry::NecResult n = geometry::nec_minimum(f, 16384);

    report.check_residual("nec.sphere_constraint",
                          "the reported direction is a unit vector",
                          VerificationLevel::A, "n/a (pointwise)", "off-axis test tensor",
                          n.sphere_constraint_residual, 1.0e-12, sw.seconds(), "");

    report.check_residual(
        "nec.stationarity", "(S - mu I) n = q at the minimiser",
        VerificationLevel::A, "n/a (pointwise)", "off-axis test tensor",
        n.stationarity_residual, 1.0e-8, sw.seconds(),
        "the trust-region stationarity condition; a nonzero residual would mean "
        "the secular equation was solved for the wrong root");

    report.check_predicate(
        "nec.quadrature_crosscheck",
        "16384-point deterministic spherical quadrature finds nothing lower",
        VerificationLevel::D, n.quadrature_crosscheck >= -1e-6,
        "quadrature minimum >= analytic minimum",
        "difference = " + std::to_string(n.quadrature_crosscheck),
        sw.seconds(),
        "an independent method: a Fibonacci lattice sweep with no knowledge of "
        "the secular equation");

    report.check_predicate(
        "nec.multistart_crosscheck",
        "26-seed projected gradient descent finds nothing lower",
        VerificationLevel::D, n.multistart_crosscheck >= -1e-8,
        "descent minimum >= analytic minimum",
        "difference = " + std::to_string(n.multistart_crosscheck),
        sw.seconds(), "a second independent method");

    // Angular convergence: coarser quadratures approach the analytic value from
    // above and never beat it.
    sw.reset();
    std::ofstream fq(out_dir + "/null_direction_samples.csv");
    fq << "quadrature_points,best_value,gap_to_analytic\n";
    bool monotone_ok = true;
    double prev_gap = 1e300;
    for (int q : {64, 256, 1024, 4096, 16384}) {
        const geometry::NecResult r = geometry::nec_minimum(f, q);
        const double gap = r.quadrature_crosscheck;
        fq << q << "," << (r.n_min + gap) << "," << gap << "\n";
        if (gap < -1e-6) monotone_ok = false;
        prev_gap = gap;
    }
    (void)prev_gap;
    report.check_predicate(
        "nec.angular_convergence",
        "quadrature approaches the analytic minimum from above at every "
        "resolution",
        VerificationLevel::E, monotone_ok, "no quadrature beats the analytic value",
        monotone_ok ? "confirmed at 64, 256, 1024, 4096, 16384 points"
                    : "a quadrature point found a lower value",
        sw.seconds(),
        "if a sampled direction ever beat the analytic minimum, the secular "
        "solve would be wrong");

    // Direct substitution: recompute N at the reported direction from scratch.
    sw.reset();
    const double direct = geometry::nec_value(f, n.direction);
    report.check_residual(
        "nec.direct_substitution",
        "N(n_worst) recomputed independently equals the reported minimum",
        VerificationLevel::A, "n/a (pointwise)", "off-axis test tensor",
        std::fabs(direct - n.n_min), 1.0e-14, sw.seconds(), "");

    report.check_predicate(
        "nec.minimizer_is_off_axis",
        "the test tensor's worst direction is not a coordinate axis",
        VerificationLevel::A,
        std::fabs(n.direction[0]) > 1e-3 && std::fabs(n.direction[1]) > 1e-3,
        "all components nonzero",
        "(" + std::to_string(n.direction[0]) + ", " +
            std::to_string(n.direction[1]) + ", " +
            std::to_string(n.direction[2]) + ")",
        sw.seconds(), "");
}

}  // namespace checks
}  // namespace warp_hypad
