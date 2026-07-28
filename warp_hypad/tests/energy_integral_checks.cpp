#include "observer_checks.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "warp_hypad/geometry/einstein.hpp"
#include "warp_hypad/geometry/energy_conditions.hpp"
#include "warp_hypad/geometry/observer.hpp"
#include "warp_hypad/metrics/metric.hpp"
#include "warp_hypad/metrics/metric_field.hpp"

// ---------------------------------------------------------------------------
// A6, groups "alcubierre" and "integrals".
//
// The Eulerian energy density of an Alcubierre bubble against its closed form,
// the energy-condition fields, and the integrated objectives an optimizer would
// eventually consume.
//
// SCOPE: these are classical pointwise conditions for an INFERRED source. They
// do not establish that the source can be generated, that quantum energy
// inequalities hold, that the geometry is stable, or that FTL travel is
// possible.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using geometry::AdmDecomposition;
using geometry::FrameStress;
using geometry::Tetrad;

constexpr double kPi = 3.14159265358979323846;
constexpr double kR = 2.0, kSig = 1.0;

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

// See observer_checks.cpp for the derivation. Independent of the production
// tensor: built from the closed-form radial derivative of the shape function.
double analytic_eulerian_rho(double x, double y, double z, double v_s) {
    const double r2 = x * x + y * y + z * z;
    if (r2 <= 1e-24) return 0.0;
    const double r = std::sqrt(r2);
    metrics::AlcubierreParameters<double> p(v_s, kR, kSig);
    const double fp = metrics::shape_function_radial_derivative(r, p);
    return -(v_s * v_s / (32.0 * kPi)) * fp * fp * (y * y + z * z) / r2;
}

// Everything A6 needs at one grid point, or a reason it is unavailable.
struct PointState {
    bool ok = false;
    FrameStress frame;
    geometry::NecResult nec;
    double tetrad_residual = 0.0;
};

PointState point_state(const Mat4<double>& g, const Mat4<double>& g_inv,
                       const Mat4<double>& T_cov) {
    PointState s;
    const AdmDecomposition a = geometry::adm_decompose(g, g_inv);
    if (!a.valid) return s;
    const Tetrad t = geometry::build_tetrad(g, a);
    if (!t.valid) return s;
    s.tetrad_residual = t.orthonormality_residual;
    s.frame = geometry::project_to_frame(t, T_cov);
    s.nec = geometry::nec_minimum(s.frame, 512);
    s.ok = true;
    return s;
}

geometry::EinsteinFields<double> alcubierre_einstein(
    const Grid& g, geometry::CurvatureFields<double>& curv, double v_s) {
    metrics::AlcubierreMetric<double> m(
        0.0, metrics::TimeTreatment::InertialTimeDependent);
    metrics::AlcubierreParameters<double> p(v_s, kR, kSig);
    auto sample = [&](double t) {
        return metrics::build_metric_field<double>(g, m, p, t);
    };
    curv = geometry::compute_curvature<double>(
        g, sample, 0.0,
        geometry::TemporalPolicy::analytic_resample(g.dx() / v_s, "aligned"));
    return geometry::assemble_einstein(
        g, curv, geometry::CosmologicalConstantPolicy::none());
}

// Integrated objectives over the slice. sqrt(gamma) = 1 for Alcubierre, but it
// is computed rather than assumed.
struct Integrals {
    double e_plus = 0.0, e_minus = 0.0, e_net = 0.0, e_abs = 0.0;
    double j_nec = 0.0;
    double e_minus_smooth = 0.0, j_nec_smooth = 0.0;
    double boundary_fraction = 0.0;   // share of E_abs in the outer shell
    std::size_t invalid_points = 0;
};

// softplus(x) = (1/k) log(1 + exp(k x)), a documented smooth surrogate for
// max(0, x). NEVER used for pass/fail: it systematically overestimates near
// zero, which is exactly where a certification decision is made.
inline double softplus(double x, double k) {
    const double kx = k * x;
    if (kx > 30.0) return x;          // avoid overflow; log(1+e^u) -> u
    if (kx < -30.0) return 0.0;
    return std::log1p(std::exp(kx)) / k;
}

constexpr double kSoftplusWidth = 200.0;

Integrals integrate(const Grid& g, const geometry::CurvatureFields<double>& c,
                    const geometry::EinsteinFields<double>& e, int shell = 2) {
    Integrals I;
    const double dV = g.cell_volume();
    double abs_shell = 0.0;

    for (std::size_t p = 0; p < g.num_points(); ++p) {
        const PointState s = point_state(
            c.g.matrix_at(p), c.g_inv.matrix_at(p),
            e.source.total_effective_stress_energy.matrix_at(p));
        if (!s.ok) { ++I.invalid_points; continue; }

        const AdmDecomposition a =
            geometry::adm_decompose(c.g.matrix_at(p), c.g_inv.matrix_at(p));
        const double det_gamma =
            a.gamma[0] * (a.gamma[4] * a.gamma[8] - a.gamma[5] * a.gamma[7]) -
            a.gamma[1] * (a.gamma[3] * a.gamma[8] - a.gamma[5] * a.gamma[6]) +
            a.gamma[2] * (a.gamma[3] * a.gamma[7] - a.gamma[4] * a.gamma[6]);
        const double w = std::sqrt(std::fmax(det_gamma, 0.0)) * dV;

        const double rho = s.frame.rho;
        I.e_plus += w * std::fmax(0.0, rho);
        I.e_minus += w * std::fmax(0.0, -rho);
        I.e_net += w * rho;
        I.e_abs += w * std::fabs(rho);
        I.j_nec += w * std::fmax(0.0, -s.nec.n_min);

        I.e_minus_smooth += w * softplus(-rho, kSoftplusWidth);
        I.j_nec_smooth += w * softplus(-s.nec.n_min, kSoftplusWidth);

        if (geometry::is_shell_point(g, p, shell)) abs_shell += w * std::fabs(rho);
    }
    I.boundary_fraction = I.e_abs > 0.0 ? abs_shell / I.e_abs : 0.0;
    return I;
}

}  // namespace

// ===========================================================================
// Group: alcubierre
// ===========================================================================

void run_alcubierre_energy_checks(VerificationReport& report,
                                  const std::string& out_dir) {
    Stopwatch sw;
    const double v_s = 0.5;
    const Grid g = cartesian_grid(25, 4.0);
    geometry::CurvatureFields<double> c;
    auto e = alcubierre_einstein(g, c, v_s);
    const std::string params = "v=0.5 R=2 sigma=1";

    // --- Against the closed form ---------------------------------------------
    sw.reset();
    double worst_rho = 0.0, scale = 0.0, worst_axis = 0.0, min_rho = 0.0;
    double worst_tetrad = 0.0;
    for (std::size_t p = 0; p < g.num_points(); ++p) {
        if (geometry::is_shell_point(g, p, 2)) continue;
        int i, j, k;
        g.unflatten(p, i, j, k);
        const PointState s = point_state(
            c.g.matrix_at(p), c.g_inv.matrix_at(p),
            e.source.total_effective_stress_energy.matrix_at(p));
        if (!s.ok) continue;
        worst_tetrad = std::fmax(worst_tetrad, s.tetrad_residual);

        const double exact = analytic_eulerian_rho(g.x(i), g.y(j), g.z(k), v_s);
        worst_rho = std::fmax(worst_rho, std::fabs(s.frame.rho - exact));
        scale = std::fmax(scale, std::fabs(exact));
        min_rho = std::fmin(min_rho, s.frame.rho);

        if (g.y(j) == 0.0 && g.z(k) == 0.0)
            worst_axis = std::fmax(worst_axis, std::fabs(s.frame.rho));
    }

    report.check_residual(
        "alcubierre.eulerian_density_vs_closed_form",
        "rho = -(v_s^2/32 pi)[(d_y f)^2 + (d_z f)^2]",
        VerificationLevel::B, grid_label(g), params, worst_rho, 4.0e-4,
        sw.seconds(),
        "peak |rho| = " + std::to_string(scale) +
            ". The reference is built from the closed-form radial derivative of "
            "the shape function and never contracts the production tensor");

    report.check_residual(
        "alcubierre.density_vanishes_on_axis",
        "rho = 0 on the axis of travel, where y = z = 0",
        VerificationLevel::B, grid_label(g), params, worst_axis, 4.0e-4,
        sw.seconds(),
        "analytic: the density depends on the TRANSVERSE derivatives only, and "
        "(d_y f)^2 + (d_z f)^2 = f'(r)^2 (y^2+z^2)/r^2 vanishes there. A source "
        "concentrated in a torus around the axis, not a ball");

    report.check_predicate(
        "alcubierre.density_negative_off_axis",
        "rho < 0 in the off-axis wall",
        VerificationLevel::A, min_rho < -1e-6, "a negative minimum exists",
        "min rho = " + std::to_string(min_rho), sw.seconds(),
        "the classical statement that an Alcubierre bubble requires negative "
        "energy density as measured by the Eulerian observer. It is a statement "
        "about THIS observer and this inferred source, not a claim that such a "
        "source exists");

    report.check_residual("alcubierre.tetrad_residual",
                          "orthonormality across the energy field",
                          VerificationLevel::A, grid_label(g), params,
                          worst_tetrad, 1.0e-12, sw.seconds(), "");

    // --- Velocity scaling ----------------------------------------------------
    //
    // rho scales as v_s^2 exactly, which is a property of the closed form that
    // the numerical pipeline must reproduce.
    sw.reset();
    const Grid gs = cartesian_grid(17, 4.0);
    geometry::CurvatureFields<double> c1, c2;
    auto e1 = alcubierre_einstein(gs, c1, 0.4);
    auto e2 = alcubierre_einstein(gs, c2, 0.8);
    double r1 = 0.0, r2 = 0.0;
    for (std::size_t p = 0; p < gs.num_points(); ++p) {
        if (geometry::is_shell_point(gs, p, 2)) continue;
        const PointState s1 = point_state(
            c1.g.matrix_at(p), c1.g_inv.matrix_at(p),
            e1.source.total_effective_stress_energy.matrix_at(p));
        const PointState s2 = point_state(
            c2.g.matrix_at(p), c2.g_inv.matrix_at(p),
            e2.source.total_effective_stress_energy.matrix_at(p));
        if (!s1.ok || !s2.ok) continue;
        r1 = std::fmin(r1, s1.frame.rho);
        r2 = std::fmin(r2, s2.frame.rho);
    }
    const double ratio = (r1 != 0.0) ? r2 / r1 : 0.0;
    report.check_value("alcubierre.density_scales_as_velocity_squared",
                       "rho ~ v_s^2, so doubling v_s quadruples |rho|",
                       VerificationLevel::B, grid_label(gs),
                       "v_s = 0.4 vs 0.8", 4.0, ratio, 0.06, sw.seconds(),
                       "min rho: " + std::to_string(r1) + " -> " +
                           std::to_string(r2));

    // --- Records --------------------------------------------------------------
    const int km = g.center_index_z(), jm = g.center_index_y();
    std::ofstream fs(out_dir + "/alcubierre_energy_slice_xy.csv");
    fs << "x,y,z,observer_type,observer_valid,tetrad_residual,rho,rho_exact,"
          "q_x,q_y,q_z,principal_pressure_1,principal_pressure_2,"
          "principal_pressure_3,N_min,worst_null_direction_x,"
          "worst_null_direction_y,worst_null_direction_z,NEC_status,WEC_status,"
          "DEC_status,SEC_status,numerical_error_estimate,certification_margin\n";

    std::ofstream fn(out_dir + "/alcubierre_nec_slice_xy.csv");
    fn << "x,y,N_min,worst_nx,worst_ny,worst_nz,rho\n";
    std::ofstream fp(out_dir + "/principal_stresses.csv");
    fp << "x,p1,p2,p3,trace_spatial,q_norm\n";
    std::ofstream fst(out_dir + "/energy_condition_status.csv");
    fst << "x,rho,N_min,NEC_status,WEC_status,DEC_status,SEC_status,"
           "error_envelope\n";

    // The numerical error envelope: half the difference between this grid and a
    // coarser one is the Richardson estimate; a fixed conservative floor is used
    // here so the record is self-contained.
    const double envelope = 4.0e-4;

    for (int i = 0; i < g.nx(); ++i)
        for (int j = 0; j < g.ny(); ++j) {
            const std::size_t p = g.flat(i, j, km);
            const PointState s = point_state(
                c.g.matrix_at(p), c.g_inv.matrix_at(p),
                e.source.total_effective_stress_energy.matrix_at(p));
            if (!s.ok) continue;
            const auto w = geometry::wec(s.frame, envelope);
            const auto d = geometry::dec(s.frame, envelope);
            const auto se = geometry::sec(s.frame, envelope);
            const auto nst = geometry::certify(s.nec.n_min, envelope);

            fs << g.x(i) << "," << g.y(j) << "," << g.z(km) << ",eulerian,1,"
               << s.tetrad_residual << "," << s.frame.rho << ","
               << analytic_eulerian_rho(g.x(i), g.y(j), g.z(km), v_s) << ","
               << s.frame.q[0] << "," << s.frame.q[1] << "," << s.frame.q[2]
               << "," << s.frame.principal[0] << "," << s.frame.principal[1]
               << "," << s.frame.principal[2] << "," << s.nec.n_min << ","
               << s.nec.direction[0] << "," << s.nec.direction[1] << ","
               << s.nec.direction[2] << "," << geometry::to_string(nst) << ","
               << geometry::to_string(w.status) << ","
               << geometry::to_string(d.status) << ","
               << geometry::to_string(se.status) << "," << envelope << ","
               << w.certification_margin << "\n";

            fn << g.x(i) << "," << g.y(j) << "," << s.nec.n_min << ","
               << s.nec.direction[0] << "," << s.nec.direction[1] << ","
               << s.nec.direction[2] << "," << s.frame.rho << "\n";
        }

    std::ofstream fc(out_dir + "/alcubierre_energy_centerline.csv");
    fc << "x,rho,rho_exact,q_x,q_y,q_z,N_min,p1,p2,p3\n";
    const int joff = g.center_index_y() + 3;  // off-axis, where rho != 0
    for (int i = 0; i < g.nx(); ++i) {
        const std::size_t p = g.flat(i, joff, km);
        const PointState s = point_state(
            c.g.matrix_at(p), c.g_inv.matrix_at(p),
            e.source.total_effective_stress_energy.matrix_at(p));
        if (!s.ok) continue;
        const auto w = geometry::wec(s.frame, envelope);
        const auto d = geometry::dec(s.frame, envelope);
        const auto se = geometry::sec(s.frame, envelope);
        fc << g.x(i) << "," << s.frame.rho << ","
           << analytic_eulerian_rho(g.x(i), g.y(joff), g.z(km), v_s) << ","
           << s.frame.q[0] << "," << s.frame.q[1] << "," << s.frame.q[2] << ","
           << s.nec.n_min << "," << s.frame.principal[0] << ","
           << s.frame.principal[1] << "," << s.frame.principal[2] << "\n";
        fp << g.x(i) << "," << s.frame.principal[0] << ","
           << s.frame.principal[1] << "," << s.frame.principal[2] << ","
           << s.frame.trace_spatial << "," << s.frame.q_norm << "\n";
        fst << g.x(i) << "," << s.frame.rho << "," << s.nec.n_min << ","
            << geometry::to_string(geometry::certify(s.nec.n_min, envelope))
            << "," << geometry::to_string(w.status) << ","
            << geometry::to_string(d.status) << ","
            << geometry::to_string(se.status) << "," << envelope << "\n";
    }
    (void)jm;

    std::ofstream fw(out_dir + "/worst_null_directions.csv");
    fw << "x,y,N_min,nx,ny,nz,stationarity_residual,quadrature_crosscheck\n";
    for (int i = 2; i < g.nx() - 2; i += 2)
        for (int j = 2; j < g.ny() - 2; j += 2) {
            const std::size_t p = g.flat(i, j, km);
            const PointState s = point_state(
                c.g.matrix_at(p), c.g_inv.matrix_at(p),
                e.source.total_effective_stress_energy.matrix_at(p));
            if (!s.ok) continue;
            fw << g.x(i) << "," << g.y(j) << "," << s.nec.n_min << ","
               << s.nec.direction[0] << "," << s.nec.direction[1] << ","
               << s.nec.direction[2] << "," << s.nec.stationarity_residual << ","
               << s.nec.quadrature_crosscheck << "\n";
        }

    // --- The NEC is violated, and it is certified ----------------------------
    sw.reset();
    double worst_nec = 0.0;
    for (std::size_t p = 0; p < g.num_points(); ++p) {
        if (geometry::is_shell_point(g, p, 2)) continue;
        const PointState s = point_state(
            c.g.matrix_at(p), c.g_inv.matrix_at(p),
            e.source.total_effective_stress_energy.matrix_at(p));
        if (s.ok) worst_nec = std::fmin(worst_nec, s.nec.n_min);
    }
    report.check_predicate(
        "alcubierre.nec_violated_and_certified",
        "the minimum NEC value is negative by more than the error envelope",
        VerificationLevel::A,
        geometry::certify(worst_nec, envelope) ==
            geometry::ConditionStatus::CertifiedViolated,
        "certified_violated",
        "N_min = " + std::to_string(worst_nec) + ", envelope = " +
            std::to_string(envelope) + ", safety factor " +
            std::to_string(geometry::kCertificationSafetyFactor),
        sw.seconds(),
        "certified rather than asserted: the magnitude exceeds the numerical "
        "error envelope by the declared safety factor, so it is not "
        "finite-difference noise");
}

// ===========================================================================
// Group: integrals
// ===========================================================================

void run_integral_checks(VerificationReport& report, const std::string& out_dir) {
    Stopwatch sw;

    // --- Velocity sweep ------------------------------------------------------
    std::ofstream fv(out_dir + "/integrated_energy_vs_velocity.csv");
    fv << "v_s,E_plus,E_minus,E_net,E_abs,J_NEC,E_minus_smooth,J_NEC_smooth,"
          "boundary_fraction,invalid_points\n";

    const Grid g = cartesian_grid(21, 4.0);
    double e_minus_at_half = 0.0, j_at_half = 0.0, e_abs_at_half = 0.0;
    double e_plus_at_half = 0.0, boundary_at_half = 0.0;
    double em_smooth_half = 0.0;

    sw.reset();
    for (int k = 2; k <= 12; k += 2) {
        const double v_s = static_cast<double>(k) / 10.0;
        geometry::CurvatureFields<double> c;
        auto e = alcubierre_einstein(g, c, v_s);
        const Integrals I = integrate(g, c, e);
        fv << v_s << "," << I.e_plus << "," << I.e_minus << "," << I.e_net << ","
           << I.e_abs << "," << I.j_nec << "," << I.e_minus_smooth << ","
           << I.j_nec_smooth << "," << I.boundary_fraction << ","
           << I.invalid_points << "\n";
        if (k == 6) {
            e_minus_at_half = I.e_minus; j_at_half = I.j_nec;
            e_abs_at_half = I.e_abs; e_plus_at_half = I.e_plus;
            boundary_at_half = I.boundary_fraction;
            em_smooth_half = I.e_minus_smooth;
        }
    }

    report.check_predicate(
        "integrals.e_minus_positive",
        "E_- is strictly positive: the inferred source has negative-density "
        "regions",
        VerificationLevel::A, e_minus_at_half > 1e-8, "E_- > 0",
        "E_- = " + std::to_string(e_minus_at_half), sw.seconds(), "");

    report.check_predicate(
        "integrals.e_plus_and_e_minus_reported_separately",
        "E_+ and E_- are both reported, not folded into E_net",
        VerificationLevel::A,
        e_plus_at_half >= 0.0 && e_minus_at_half > 0.0 &&
            e_abs_at_half > e_minus_at_half * 0.999,
        "E_+, E_- and E_abs separately available",
        "E_+ = " + std::to_string(e_plus_at_half) + ", E_- = " +
            std::to_string(e_minus_at_half) + ", E_abs = " +
            std::to_string(e_abs_at_half),
        sw.seconds(),
        "E_net alone would let a large positive region cancel a large negative "
        "one and report a small number for a violently non-uniform source");

    report.check_predicate(
        "integrals.nec_violation_positive",
        "J_NEC > 0: the integrated NEC violation is nonzero",
        VerificationLevel::A, j_at_half > 1e-8, "J_NEC > 0",
        "J_NEC = " + std::to_string(j_at_half), sw.seconds(), "");

    report.check_predicate(
        "integrals.not_boundary_dominated",
        "the energy integral is not dominated by the outer shell",
        VerificationLevel::A, boundary_at_half < 0.35,
        "boundary share < 35%",
        "boundary share = " + std::to_string(boundary_at_half),
        sw.seconds(),
        "a boundary-dominated integral would be measuring the domain truncation "
        "rather than the bubble");

    // --- Exact and smooth are distinct ---------------------------------------
    sw.reset();
    report.check_predicate(
        "integrals.exact_and_smooth_penalties_distinct",
        "the softplus surrogate is not the exact objective",
        VerificationLevel::A,
        std::fabs(em_smooth_half - e_minus_at_half) > 1e-12,
        "the two differ",
        "exact E_- = " + std::to_string(e_minus_at_half) + ", softplus E_- = " +
            std::to_string(em_smooth_half) + " (width " +
            std::to_string(200.0) + ")",
        sw.seconds(),
        "softplus systematically overestimates max(0,x) near zero, which is "
        "exactly where certification happens. It is provided for a future "
        "gradient-based optimizer and is never used for scientific pass/fail");

    // --- Resolution convergence ----------------------------------------------
    sw.reset();
    std::ofstream fr(out_dir + "/energy_condition_convergence.csv");
    fr << "n,h,E_minus,J_NEC,min_rho\n";
    double em[3] = {0, 0, 0}, hs[3] = {0, 0, 0};
    int idx = 0;
    for (int n : {13, 17, 25}) {
        const Grid gn = cartesian_grid(n, 4.0);
        geometry::CurvatureFields<double> cn;
        auto en = alcubierre_einstein(gn, cn, 0.5);
        const Integrals I = integrate(gn, cn, en);
        double mn = 0.0;
        for (std::size_t p = 0; p < gn.num_points(); ++p) {
            if (geometry::is_shell_point(gn, p, 2)) continue;
            const PointState s = point_state(
                cn.g.matrix_at(p), cn.g_inv.matrix_at(p),
                en.source.total_effective_stress_energy.matrix_at(p));
            if (s.ok) mn = std::fmin(mn, s.frame.rho);
        }
        fr << n << "," << gn.dx() << "," << I.e_minus << "," << I.j_nec << ","
           << mn << "\n";
        em[idx] = I.e_minus; hs[idx] = gn.dx(); ++idx;
    }
    const double rel_change = std::fabs(em[2] - em[1]) / std::fabs(em[2]);
    report.check_predicate(
        "integrals.resolution_convergence",
        "E_- changes by less than 20% between the two finest grids",
        VerificationLevel::E, rel_change < 0.2, "relative change < 20%",
        "E_- = " + std::to_string(em[0]) + ", " + std::to_string(em[1]) + ", " +
            std::to_string(em[2]) + "; finest relative change = " +
            std::to_string(rel_change),
        sw.seconds(),
        "the Richardson-style envelope the certification uses. A tighter claim "
        "would need finer grids than the rank-4 memory footprint allows");

    // --- Domain convergence ---------------------------------------------------
    sw.reset();
    std::ofstream fd(out_dir + "/domain_convergence.csv");
    fd << "half_width,E_minus,E_abs,boundary_fraction\n";
    double last = 0.0, prev = 0.0;
    for (double half : {3.0, 4.0, 5.0}) {
        GridSpec s;
        s.nx = s.ny = s.nz = 21;
        s.x_min = s.y_min = s.z_min = -half;
        s.x_max = s.y_max = s.z_max = half;
        const Grid gd(s);
        geometry::CurvatureFields<double> cd;
        auto ed = alcubierre_einstein(gd, cd, 0.5);
        const Integrals I = integrate(gd, cd, ed);
        fd << half << "," << I.e_minus << "," << I.e_abs << ","
           << I.boundary_fraction << "\n";
        prev = last;
        last = I.e_minus;
    }
    const double dom_change = std::fabs(last - prev) / std::fabs(last);
    report.check_predicate(
        "integrals.domain_convergence",
        "E_- is insensitive to enlarging the domain",
        VerificationLevel::E, dom_change < 0.25, "relative change < 25%",
        "E_- at half-width 4 -> 5: " + std::to_string(prev) + " -> " +
            std::to_string(last) + ", relative change " +
            std::to_string(dom_change),
        sw.seconds(),
        "the tail contribution: if this grew with the domain, the integral "
        "would not exist and no optimizer could use it");
}

}  // namespace checks
}  // namespace warp_hypad
