#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "warp_hypad/core/metadata.hpp"
#include "warp_hypad/geometry/einstein.hpp"
#include "warp_hypad/geometry/energy_conditions.hpp"
#include "warp_hypad/geometry/observer.hpp"
#include "warp_hypad/metrics/metric.hpp"
#include "warp_hypad/metrics/metric_field.hpp"

// ---------------------------------------------------------------------------
// A6 observer-energy sweep exporter.
//
// The SMALLEST exporter that lets the A6 visuals be drawn from production
// output. Every number below comes out of observer.hpp, energy_conditions.hpp
// or einstein.hpp; this file contains no physics of its own beyond one explicit
// contraction, which is validated against production in the record itself (see
// contract_uu).
//
// WHY IT EXISTS. The A6 gate writes a single z = 0 slice at v_s = 0.5. The
// visuals need the same fields across the established velocity sequence, plus
// several quantities production computes and then discards: the per-condition
// certification margins, NecResult::direct_substitution, the tetrad handedness
// and positive-diagonal witness, the interior/boundary split of the integrated
// objectives, and the contracted-Bianchi diagnostic (A6 calls assemble_einstein,
// which never computes the divergence; this calls compute_einstein, which does).
// See visualization/program_a6/a6_capability_audit.json for the full audit.
//
// OBSERVER SEMANTICS ARE NOT NEGOTIABLE HERE. The coordinate-stationary observer
// is admissible only where g_00 < 0. At g_00 = 0 the vector d/dt is null and for
// g_00 > 0 it is spacelike; in both cases production REJECTS it, and this
// exporter writes an EMPTY energy-density cell plus the rejection reason rather
// than a number. There is no |g_00| anywhere in this file. The Eulerian observer
// is built from the ADM normal and remains timelike straight through, and it is
// recorded as a separate row so the two can never be conflated downstream.
//
// SCOPE. A6 evaluates observer-dependent energy density and classical pointwise
// energy conditions for the stress-energy tensor INFERRED from the selected
// metric. It does not demonstrate that the source can be generated, that quantum
// energy inequalities are satisfied, that the geometry is dynamically stable, or
// that operational FTL travel is possible.
// ---------------------------------------------------------------------------

namespace {

using namespace warp_hypad;

// The velocity sequence is the one V0.1A established and A4 and A5 reuse:
// v_s = k/10 for k = 2..16. k = 10 is EXACTLY 1.0 in binary floating point, so
// the onset frame is exact rather than nearly exact.
constexpr int kFirstK = 2;
constexpr int kLastK = 16;

constexpr int kGridN = 21;
constexpr double kHalfWidth = 4.0;
constexpr double kTailHalfWidth = 5.0;   // the exterior-tail probe
constexpr int kCoarseN = 13;             // the two-grid error estimate
constexpr double kRadius = 2.0;
constexpr double kSigma = 1.0;
constexpr double kTime = 0.0;
constexpr int kShell = 2;                // boundary-shell depth, as in the gate
constexpr int kNecQuadrature = 512;
constexpr int kSphereSamples = 1024;
constexpr double kReferenceVelocity = 0.5;   // the A6 gate's configuration

// The numerical error envelope. A6 uses a FIXED conservative floor rather than a
// per-point Richardson estimate, and that choice is inherited rather than
// second-guessed: a visual must never certify something the science left
// uncertified. The exporter additionally measures a two-grid difference per
// frame and takes whichever is LARGER, so the envelope can only ever grow.
constexpr double kEnvelopeFloor = 4.0e-4;

// softplus(x) = (1/k) log(1 + exp(k x)). The documented smooth surrogate a
// future gradient-based optimizer would consume. It is exported ALONGSIDE the
// exact objective and never in place of it: softplus systematically
// overestimates max(0, x) near zero, which is exactly where certification
// happens.
constexpr double kSoftplusWidth = 200.0;

Grid cartesian_grid(int n, double half) {
    GridSpec s;
    s.nx = s.ny = s.nz = n;
    s.x_min = s.y_min = s.z_min = -half;
    s.x_max = s.y_max = s.z_max = half;
    return Grid(s);
}

double softplus(double x, double k) {
    const double kx = k * x;
    if (kx > 30.0) return x;
    if (kx < -30.0) return 0.0;
    return std::log1p(std::exp(kx)) / k;
}

// T_{mu nu} u^mu u^nu.
//
// The one contraction this file performs. It is the same one
// geometry::project_to_frame applies to the tetrad legs, written out because
// production exposes it only through a Tetrad and the observer families that
// are not the Eulerian one do not come with one. It is VALIDATED IN THE RECORD:
// for the Eulerian observer the result is compared against
// FrameStress::rho and the residual is written to
// observer_comparison.csv, so a reader can confirm the helper agrees with
// production rather than take it on trust.
double contract_uu(const Mat4<double>& T, const double u[4]) {
    double s = 0.0;
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu) s += T(mu, nu) * u[mu] * u[nu];
    return s;
}

// Everything A6 knows at one grid point, or a reason there is nothing to know.
struct PointState {
    bool ok = false;
    geometry::FrameStress frame;
    geometry::NecResult nec;
    double tetrad_residual = 0.0;
    double tetrad_handedness = 0.0;
    // min over a = 1..3 of e_a^a. Fixed-order Gram-Schmidt with a positive root
    // makes these positive BY CONSTRUCTION, so a negative one is an arbitrary
    // sign flip and nothing else.
    double tetrad_min_diagonal = 0.0;
    double eulerian_norm_residual = 0.0;
    bool eulerian_valid = false;
    double sqrt_gamma = 0.0;
    double nec_eigen_gap = 0.0;   // lambda_2 - lambda_1 of S: the degeneracy witness
};

PointState point_state(const Mat4<double>& g, const Mat4<double>& g_inv,
                       const Mat4<double>& T_cov, bool with_nec) {
    PointState s;
    const geometry::AdmDecomposition a = geometry::adm_decompose(g, g_inv);
    if (!a.valid) return s;
    const geometry::Tetrad t = geometry::build_tetrad(g, a);
    if (!t.valid) return s;

    const geometry::Observer eul = geometry::make_eulerian(g, a);
    s.eulerian_valid = eul.valid;
    s.eulerian_norm_residual = eul.normalization_residual;

    s.tetrad_residual = t.orthonormality_residual;
    s.tetrad_handedness = t.handedness;
    s.tetrad_min_diagonal = t.e[1][1];
    if (t.e[2][2] < s.tetrad_min_diagonal) s.tetrad_min_diagonal = t.e[2][2];
    if (t.e[3][3] < s.tetrad_min_diagonal) s.tetrad_min_diagonal = t.e[3][3];

    const double det_gamma =
        a.gamma[0] * (a.gamma[4] * a.gamma[8] - a.gamma[5] * a.gamma[7]) -
        a.gamma[1] * (a.gamma[3] * a.gamma[8] - a.gamma[5] * a.gamma[6]) +
        a.gamma[2] * (a.gamma[3] * a.gamma[7] - a.gamma[4] * a.gamma[6]);
    s.sqrt_gamma = std::sqrt(std::fmax(det_gamma, 0.0));

    s.frame = geometry::project_to_frame(t, T_cov);
    if (with_nec) {
        s.nec = geometry::nec_minimum(s.frame, kNecQuadrature);
        double lam[3];
        geometry::symmetric3_eigenvalues(s.frame.S, lam);
        s.nec_eigen_gap = lam[1] - lam[0];
    }
    s.ok = true;
    return s;
}

// The integrated objectives, split so that a boundary-dominated result cannot
// hide inside a single number.
struct Integrals {
    double e_plus = 0.0, e_minus = 0.0, e_net = 0.0, e_abs = 0.0, j_nec = 0.0;
    double e_minus_smooth = 0.0, j_nec_smooth = 0.0;
    double e_plus_interior = 0.0, e_minus_interior = 0.0, e_abs_interior = 0.0;
    double e_plus_shell = 0.0, e_minus_shell = 0.0, e_abs_shell = 0.0;
    double j_nec_interior = 0.0, j_nec_shell = 0.0;
    double boundary_fraction = 0.0;
    double min_rho = 0.0;
    double min_rho_x = 0.0, min_rho_y = 0.0, min_rho_z = 0.0;
    double n_min = 0.0;
    double n_dir[3] = {0, 0, 0};
    double n_min_x = 0.0, n_min_y = 0.0, n_min_z = 0.0;
    std::size_t invalid_points = 0;
};

Integrals integrate(const Grid& g, const geometry::CurvatureFields<double>& c,
                    const geometry::EinsteinFields<double>& e, bool with_nec) {
    Integrals I;
    const double dV = g.cell_volume();

    for (std::size_t p = 0; p < g.num_points(); ++p) {
        const PointState s = point_state(
            c.g.matrix_at(p), c.g_inv.matrix_at(p),
            e.source.total_effective_stress_energy.matrix_at(p), with_nec);
        if (!s.ok) { ++I.invalid_points; continue; }

        const double w = s.sqrt_gamma * dV;
        const double rho = s.frame.rho;
        const bool shell = geometry::is_shell_point(g, p, kShell);

        const double ep = w * std::fmax(0.0, rho);
        const double em = w * std::fmax(0.0, -rho);
        const double ea = w * std::fabs(rho);

        I.e_plus += ep;
        I.e_minus += em;
        // Accumulated independently of E_+ and E_-, which is what makes the
        // closure check |E_net - (E_+ - E_-)| < eps a test and not an identity.
        I.e_net += w * rho;
        I.e_abs += ea;
        I.e_minus_smooth += w * softplus(-rho, kSoftplusWidth);

        if (shell) {
            I.e_plus_shell += ep; I.e_minus_shell += em; I.e_abs_shell += ea;
        } else {
            I.e_plus_interior += ep; I.e_minus_interior += em;
            I.e_abs_interior += ea;
        }

        if (with_nec) {
            const double jn = w * std::fmax(0.0, -s.nec.n_min);
            I.j_nec += jn;
            I.j_nec_smooth += w * softplus(-s.nec.n_min, kSoftplusWidth);
            if (shell) I.j_nec_shell += jn; else I.j_nec_interior += jn;
        }

        // Extrema are reported over the INTERIOR only. The shell carries the
        // one-sided derivative treatment, so a minimum found there would be
        // reporting the boundary rather than the bubble.
        if (!shell) {
            int i, j, k;
            g.unflatten(p, i, j, k);
            if (rho < I.min_rho) {
                I.min_rho = rho;
                I.min_rho_x = g.x(i); I.min_rho_y = g.y(j); I.min_rho_z = g.z(k);
            }
            if (with_nec && s.nec.n_min < I.n_min) {
                I.n_min = s.nec.n_min;
                for (int d = 0; d < 3; ++d) I.n_dir[d] = s.nec.direction[d];
                I.n_min_x = g.x(i); I.n_min_y = g.y(j); I.n_min_z = g.z(k);
            }
        }
    }
    I.boundary_fraction = I.e_abs > 0.0 ? I.e_abs_shell / I.e_abs : 0.0;
    return I;
}

// The Alcubierre chain, one velocity, one grid.
geometry::EinsteinFields<double> alcubierre_einstein(
    const Grid& g, geometry::CurvatureFields<double>& curv, double v_s,
    bool with_bianchi) {
    metrics::AlcubierreMetric<double> m(
        0.0, metrics::TimeTreatment::InertialTimeDependent);
    metrics::AlcubierreParameters<double> p(v_s, kRadius, kSigma);
    auto sample = [&](double t) {
        return metrics::build_metric_field<double>(g, m, p, t);
    };
    const auto policy = geometry::TemporalPolicy::analytic_resample(
        g.dx() / v_s,
        "h_t = h_x/v_s aligns the temporal and spatial stencils along "
        "u = x - v_s t");
    const auto lp = geometry::CosmologicalConstantPolicy::none();

    if (with_bianchi) {
        // compute_einstein performs THREE full curvature computations, because
        // the contracted-Bianchi residual needs d_0 G^{mu nu}. A6's own checks
        // call assemble_einstein and therefore leave divergence_relative at
        // zero; the frame overlay is required to show it, so it is paid for
        // here rather than reported as a fabricated zero.
        auto e = geometry::compute_einstein<double>(g, sample, kTime, policy, lp);
        curv = geometry::compute_curvature<double>(g, sample, kTime, policy);
        return e;
    }
    curv = geometry::compute_curvature<double>(g, sample, kTime, policy);
    return geometry::assemble_einstein(g, curv, lp);
}

const char* status_name(geometry::ConditionStatus s) {
    return geometry::to_string(s);
}

// Cells resolving the bubble wall, MEASURED rather than assumed: the count of
// grid samples on the +x half of the axis where the ADM shift magnitude
// |beta^x| = v_s f lies strictly between 5% and 95% of v_s. Production output
// throughout -- beta comes from adm_decompose.
int wall_cells(const Grid& g, const Rank2Field<double>& gf,
               const geometry::CurvatureFields<double>& c, double v_s) {
    const int jm = g.center_index_y(), km = g.center_index_z();
    int n = 0;
    for (int i = 0; i < g.nx(); ++i) {
        if (!(g.x(i) > 0.0)) continue;
        const std::size_t p = g.flat(i, jm, km);
        const geometry::AdmDecomposition a =
            geometry::adm_decompose(c.g.matrix_at(p), c.g_inv.matrix_at(p));
        if (!a.valid) continue;
        const double frac = std::fabs(a.beta[0]) / v_s;
        if (frac > 0.05 && frac < 0.95) ++n;
    }
    (void)gf;
    return n;
}

std::string csv_quote(const std::string& s) { return "\"" + s + "\""; }

// One observer sample site, evaluated for all four families.
struct SiteRow {
    const char* site;
    int i, j, k;
};

}  // namespace

int main(int argc, char** argv) {
    std::string out_dir = "results/verification/v0_program_a6_energy";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: a6_energy_sweep [--out <dir>]\n";
            return 2;
        }
    }

    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) {
        std::cerr << "cannot create " << out_dir << ": " << ec.message() << "\n";
        return 3;
    }

    const Grid grid = cartesian_grid(kGridN, kHalfWidth);
    const int im = grid.center_index_x();
    const int jm = grid.center_index_y();
    const int km = grid.center_index_z();
    if (im < 0 || jm < 0 || km < 0) {
        std::cerr << "grid extents must be odd so the centre is pinned\n";
        return 3;
    }

    // The two fixed observer sample sites. FIXED, so observer identity and site
    // identity are both stable across the animation: a site that moved would be
    // indistinguishable from a physical change.
    //
    //   bubble_center  (0,0,0)   -- where g_00 crosses zero first, since f = 1
    //   wall_probe     (0,R,0)   -- off the axis of travel, where rho != 0
    int j_wall = jm;
    for (int j = 0; j < grid.ny(); ++j)
        if (std::fabs(grid.y(j) - kRadius) < std::fabs(grid.y(j_wall) - kRadius))
            j_wall = j;
    const SiteRow kSites[2] = {{"bubble_center", im, jm, km},
                               {"wall_probe", im, j_wall, km}};

    std::ofstream sweep(out_dir + "/energy_velocity_sweep.csv");
    sweep << std::setprecision(12);
    // The *_valid_sites columns are COUNTS over the two fixed observer sites,
    // not booleans: the coordinate-stationary observer can be admissible at the
    // wall and inadmissible at the centre in the same frame, and a single flag
    // would have to pick one and hide the other.
    sweep << "frame,v_s,g00_center,eulerian_valid_sites,"
             "stationary_valid_sites,comoving_valid_sites,"
             "boosted_valid_sites,observer_sites,rho_min,rho_min_x,rho_min_y,"
             "rho_min_z,rho_max,slice_rho_min,slice_rho_min_x,slice_rho_min_y,"
             "slice_N_min,slice_N_min_x,slice_N_min_y,"
             "N_min,worst_direction_x,worst_direction_y,"
             "worst_direction_z,N_min_x,N_min_y,N_min_z,E_plus,E_minus,E_net,"
             "E_abs,J_NEC,numerical_error_estimate,richardson_difference,"
             "certification_margin,rho_sign_status,NEC_status,WEC_status,"
             "DEC_status,SEC_status,certified_violated_cells,"
             "indeterminate_cells,certified_satisfied_cells,"
             "bianchi_relative_residual,grid_resolution,domain_extent,"
             "wall_cells,boundary_fraction,grid_brackets_g00_zero\n";

    std::ofstream necsweep(out_dir + "/nec_velocity_sweep.csv");
    necsweep << std::setprecision(12);
    necsweep << "frame,v_s,N_min,worst_direction_x,worst_direction_y,"
                "worst_direction_z,nec_stationarity_residual,"
                "nec_crosscheck_residual,nec_multistart_residual,"
                "nec_direct_substitution,nec_eigen_gap,nec_hard_case,J_NEC,"
                "J_NEC_smooth,softplus_width,numerical_error_estimate,"
                "certification_margin,NEC_status,nec_violating_cells,"
                "slice_cells\n";

    std::ofstream wnd(out_dir + "/worst_null_directions.csv");
    wnd << std::setprecision(12);
    wnd << "frame,v_s,x,y,z,N_min,worst_direction_x,worst_direction_y,"
           "worst_direction_z,nec_stationarity_residual,"
           "nec_crosscheck_residual,nec_eigen_gap,rho,NEC_status\n";

    std::ofstream obs(out_dir + "/observer_comparison.csv");
    obs << std::setprecision(12);
    obs << "frame,v_s,site,x,y,z,g00,observer_type,observer_valid,"
           "observer_rejection_reason,observer_normalization_residual,"
           "tetrad_orthonormality_residual,rho_measured,"
           "eulerian_rho_crosscheck_residual,admissibility_rule\n";

    std::ofstream integ(out_dir + "/integrated_objectives.csv");
    integ << std::setprecision(12);
    integ << "frame,v_s,E_plus,E_minus,E_net,E_abs,J_NEC,exact_value,"
             "smooth_surrogate_value,softplus_width,J_NEC_smooth,"
             "net_closure_residual,abs_closure_residual,"
             "interior_contribution,boundary_shell_contribution,"
             "tail_estimate,tail_estimate_e_abs,domain_extent,"
             "tail_domain_extent,grid_resolution,wall_cells,boundary_fraction,"
             "e_plus_interior,e_plus_shell,e_minus_interior,e_minus_shell,"
             "j_nec_interior,j_nec_shell,invalid_points,"
             "optimization_ready\n";

    int exact_onset_frame = -1;
    int first_contour_frame = -1;
    int reference_frame = -1;
    bool all_finite = true;
    double worst_eulerian_residual = 0.0;
    double worst_tetrad_residual = 0.0;
    double worst_net_closure = 0.0, worst_abs_closure = 0.0;
    double handedness_min = 1e300, handedness_max = -1e300;
    double tetrad_diag_min = 1e300;
    int frame_count = 0;
    int total_frame_rows = 0;

    for (int k = kFirstK; k <= kLastK; ++k) {
        const int frame = k - kFirstK;
        const double v_s = static_cast<double>(k) / 10.0;
        if (v_s == kReferenceVelocity) reference_frame = frame;

        geometry::CurvatureFields<double> c;
        const auto e = alcubierre_einstein(grid, c, v_s, /*with_bianchi=*/true);
        if (!e.ok || !c.ok) {
            std::cerr << "frame " << frame << " (v_s=" << v_s
                      << ") refused: " << e.refusal_reason << "\n";
            return 1;
        }

        metrics::AlcubierreMetric<double> metric(
            0.0, metrics::TimeTreatment::InertialTimeDependent);
        metrics::AlcubierreParameters<double> params(v_s, kRadius, kSigma);
        const Rank2Field<double> gf =
            metrics::build_metric_field<double>(grid, metric, params, kTime);

        // --- The integrated objectives, on the production domain -------------
        const Integrals I = integrate(grid, c, e, /*with_nec=*/true);

        // --- The exterior tail, measured ------------------------------------
        //
        // A second integration on an enlarged domain. The point count is held
        // fixed, so h grows with the extent -- the same convention the A6 gate's
        // domain_convergence uses, and stated rather than hidden. The NEC is not
        // re-minimised here: this probe answers "how much energy is outside the
        // box", and paying for the sphere minimisation again would treble the
        // sweep for a number it does not report.
        const Grid gtail = cartesian_grid(kGridN, kTailHalfWidth);
        geometry::CurvatureFields<double> ctail;
        const auto etail =
            alcubierre_einstein(gtail, ctail, v_s, /*with_bianchi=*/false);
        const Integrals Itail = integrate(gtail, ctail, etail, /*with_nec=*/false);
        const double tail_e_minus = std::fabs(Itail.e_minus - I.e_minus);
        const double tail_e_abs = std::fabs(Itail.e_abs - I.e_abs);

        // --- The two-grid error estimate ------------------------------------
        const Grid gcoarse = cartesian_grid(kCoarseN, kHalfWidth);
        geometry::CurvatureFields<double> ccoarse;
        const auto ecoarse =
            alcubierre_einstein(gcoarse, ccoarse, v_s, /*with_bianchi=*/false);
        const Integrals Icoarse =
            integrate(gcoarse, ccoarse, ecoarse, /*with_nec=*/false);
        const double richardson = std::fabs(I.min_rho - Icoarse.min_rho);
        const double envelope = std::fmax(kEnvelopeFloor, richardson);

        // --- Closure -----------------------------------------------------------
        const double net_closure = std::fabs(I.e_net - (I.e_plus - I.e_minus));
        const double abs_closure = std::fabs(I.e_abs - (I.e_plus + I.e_minus));
        worst_net_closure = std::fmax(worst_net_closure, net_closure);
        worst_abs_closure = std::fmax(worst_abs_closure, abs_closure);

        // --- The z = 0 slice record -------------------------------------------
        char name[64];
        std::snprintf(name, sizeof(name), "/energy_frame_%03d.csv", frame);
        std::ofstream f(out_dir + name);
        f << std::setprecision(12);
        f << "frame,v_s,x,y,z,g00,observer_type,observer_valid,"
             "observer_rejection_reason,observer_normalization_residual,"
             "tetrad_orthonormality_residual,tetrad_handedness,"
             "tetrad_min_diagonal,sqrt_gamma,rho,q_x,q_y,q_z,"
             "principal_pressure_1,principal_pressure_2,principal_pressure_3,"
             "N_min,worst_direction_x,worst_direction_y,worst_direction_z,"
             "nec_stationarity_residual,nec_crosscheck_residual,"
             "nec_direct_substitution,nec_eigen_gap,nec_hard_case,"
             "rho_sign_status,rho_certification_margin,NEC_status,WEC_status,"
             "DEC_status,SEC_status,nec_certification_margin,"
             "wec_certification_margin,dec_certification_margin,"
             "sec_certification_margin,numerical_error_estimate,"
             "certification_margin,bianchi_relative_residual\n";

        double rho_max = -1e300, rho_min_slice = 1e300;
        double rho_min_slice_x = 0.0, rho_min_slice_y = 0.0;
        double n_min_slice_x = 0.0, n_min_slice_y = 0.0;
        int cert_violated = 0, cert_satisfied = 0, indeterminate = 0;
        int nec_violating = 0, slice_cells = 0;
        bool saw_pos = false, saw_neg = false;
        double slice_n_min = 0.0;
        geometry::NecResult slice_worst;
        double slice_worst_eigen_gap = 0.0;

        for (int i = 0; i < grid.nx(); ++i)
            for (int j = 0; j < grid.ny(); ++j) {
                const std::size_t p = grid.flat(i, j, km);
                const double g00 = gf.at(p, 0, 0);
                if (g00 > 0.0) saw_pos = true;
                if (g00 < 0.0) saw_neg = true;

                const PointState s = point_state(
                    c.g.matrix_at(p), c.g_inv.matrix_at(p),
                    e.source.total_effective_stress_energy.matrix_at(p), true);
                if (!s.ok) continue;
                ++slice_cells;
                ++total_frame_rows;

                const auto w = geometry::wec(s.frame, envelope);
                const auto d = geometry::dec(s.frame, envelope);
                const auto se = geometry::sec(s.frame, envelope);
                const auto nst = geometry::certify(s.nec.n_min, envelope);
                const auto rst = geometry::certify(s.frame.rho, envelope);

                if (nst == geometry::ConditionStatus::CertifiedViolated)
                    ++cert_violated, ++nec_violating;
                else if (nst == geometry::ConditionStatus::CertifiedSatisfied)
                    ++cert_satisfied;
                else
                    ++indeterminate;

                if (!std::isfinite(s.frame.rho) || !std::isfinite(s.nec.n_min) ||
                    !std::isfinite(g00))
                    all_finite = false;

                rho_max = std::fmax(rho_max, s.frame.rho);
                if (s.frame.rho < rho_min_slice) {
                    rho_min_slice = s.frame.rho;
                    rho_min_slice_x = grid.x(i);
                    rho_min_slice_y = grid.y(j);
                }
                if (s.nec.n_min < slice_n_min) {
                    slice_n_min = s.nec.n_min;
                    slice_worst = s.nec;
                    slice_worst_eigen_gap = s.nec_eigen_gap;
                    n_min_slice_x = grid.x(i);
                    n_min_slice_y = grid.y(j);
                }
                worst_eulerian_residual =
                    std::fmax(worst_eulerian_residual, s.eulerian_norm_residual);
                worst_tetrad_residual =
                    std::fmax(worst_tetrad_residual, s.tetrad_residual);
                handedness_min = std::fmin(handedness_min, s.tetrad_handedness);
                handedness_max = std::fmax(handedness_max, s.tetrad_handedness);
                tetrad_diag_min =
                    std::fmin(tetrad_diag_min, s.tetrad_min_diagonal);

                const double nec_margin = std::fabs(s.nec.n_min) -
                    geometry::kCertificationSafetyFactor * envelope;
                const double rho_margin = std::fabs(s.frame.rho) -
                    geometry::kCertificationSafetyFactor * envelope;

                f << frame << "," << v_s << "," << grid.x(i) << ","
                  << grid.y(j) << "," << grid.z(km) << "," << g00
                  << ",eulerian," << (s.eulerian_valid ? 1 : 0) << ",\"\","
                  << s.eulerian_norm_residual << "," << s.tetrad_residual << ","
                  << s.tetrad_handedness << "," << s.tetrad_min_diagonal << ","
                  << s.sqrt_gamma << "," << s.frame.rho << "," << s.frame.q[0]
                  << "," << s.frame.q[1] << "," << s.frame.q[2] << ","
                  << s.frame.principal[0] << "," << s.frame.principal[1] << ","
                  << s.frame.principal[2] << "," << s.nec.n_min << ","
                  << s.nec.direction[0] << "," << s.nec.direction[1] << ","
                  << s.nec.direction[2] << "," << s.nec.stationarity_residual
                  << "," << s.nec.quadrature_crosscheck << ","
                  << s.nec.direct_substitution << "," << s.nec_eigen_gap << ","
                  << (s.nec.hard_case ? 1 : 0) << "," << status_name(rst) << ","
                  << rho_margin << "," << status_name(nst) << ","
                  << status_name(w.status) << "," << status_name(d.status) << ","
                  << status_name(se.status) << "," << nec_margin << ","
                  << w.certification_margin << "," << d.certification_margin
                  << "," << se.certification_margin << "," << envelope << ","
                  << nec_margin << "," << e.divergence_relative << "\n";

                if (i >= 2 && j >= 2 && i < grid.nx() - 2 && j < grid.ny() - 2 &&
                    (i % 2 == 0) && (j % 2 == 0)) {
                    wnd << frame << "," << v_s << "," << grid.x(i) << ","
                        << grid.y(j) << "," << grid.z(km) << "," << s.nec.n_min
                        << "," << s.nec.direction[0] << ","
                        << s.nec.direction[1] << "," << s.nec.direction[2]
                        << "," << s.nec.stationarity_residual << ","
                        << s.nec.quadrature_crosscheck << "," << s.nec_eigen_gap
                        << "," << s.frame.rho << "," << status_name(nst) << "\n";
                }
            }

        const bool brackets = saw_pos && saw_neg;
        if (brackets && first_contour_frame < 0) first_contour_frame = frame;

        const std::size_t centre = grid.flat(im, jm, km);
        const double g00_centre = gf.at(centre, 0, 0);
        if (v_s == 1.0 && g00_centre == 0.0 && exact_onset_frame < 0)
            exact_onset_frame = frame;

        // --- The four observer families, at the two fixed sites ---------------
        int n_eul_valid = 0, n_sta_valid = 0, n_com_valid = 0, n_bst_valid = 0;
        for (const SiteRow& site : kSites) {
            const std::size_t p = grid.flat(site.i, site.j, site.k);
            const Mat4<double> gm = gf.matrix_at(p);
            const InverseResult<double> inv = inverse(gm);
            const geometry::AdmDecomposition a =
                geometry::adm_decompose(gm, inv.inverse);
            const geometry::Tetrad t = geometry::build_tetrad(gm, a);
            const Mat4<double> T =
                e.source.total_effective_stress_energy.matrix_at(p);
            const double g00 = gm(0, 0);

            const double bv[3] = {0.3, 0.2, 0.0};
            const geometry::Observer fam[4] = {
                geometry::make_eulerian(gm, a),
                geometry::make_coordinate_stationary(gm),
                geometry::make_bubble_comoving(gm, v_s),
                geometry::make_tetrad_boosted(gm, t, bv)};
            const char* rule[4] = {
                "the ADM normal; timelike wherever the foliation exists, "
                "including through g_00 = 0",
                "admissible ONLY where g_00 < 0. At g_00 = 0 d/dt is null and "
                "for g_00 > 0 it is spacelike; neither is an observer, and "
                "|g_00| is never substituted",
                "the bubble-centre coordinate trajectory; timelikeness is "
                "checked, not assumed",
                "gamma(1, v) in the orthonormal frame, |v| < 1 strictly"};

            for (int fi = 0; fi < 4; ++fi) {
                const geometry::Observer& o = fam[fi];
                if (o.valid) {
                    if (fi == 0) ++n_eul_valid;
                    if (fi == 1) ++n_sta_valid;
                    if (fi == 2) ++n_com_valid;
                    if (fi == 3) ++n_bst_valid;
                }
                obs << frame << "," << v_s << "," << site.site << ","
                    << grid.x(site.i) << "," << grid.y(site.j) << ","
                    << grid.z(site.k) << "," << g00 << ","
                    << geometry::to_string(o.kind) << ","
                    << (o.valid ? 1 : 0) << ","
                    << csv_quote(o.rejection_reason) << ","
                    << o.normalization_residual << ","
                    << (t.valid ? t.orthonormality_residual : -1.0) << ",";
                // An INADMISSIBLE observer measures nothing. The cell is left
                // EMPTY rather than filled with a plausible number: a renderer
                // that tried to plot it would have to notice.
                //
                // The condition is o.valid alone, deliberately. T_{mu nu} u^mu
                // u^nu needs the observer and the tensor, not a tetrad; only
                // the Eulerian cross-check below does, and make_tetrad_boosted
                // already refuses to produce an observer without one. Gating
                // the value on the tetrad as well would leave a row marked
                // admissible with an empty measurement, which is a state the
                // reader is right to reject.
                if (o.valid) {
                    const double rho_u = contract_uu(T, o.u);
                    obs << rho_u << ",";
                    if (fi == 0 && t.valid) {
                        const geometry::FrameStress fs =
                            geometry::project_to_frame(t, T);
                        obs << std::fabs(rho_u - fs.rho) << ",";
                    } else {
                        obs << "0,";
                    }
                } else {
                    obs << ",,";
                }
                obs << csv_quote(rule[fi]) << "\n";
            }
        }

        const double nec_margin_frame =
            std::fabs(I.n_min) - geometry::kCertificationSafetyFactor * envelope;
        const auto frame_nec_status = geometry::certify(I.n_min, envelope);
        const auto frame_rho_status = geometry::certify(I.min_rho, envelope);
        const PointState centre_state = point_state(
            c.g.matrix_at(centre), c.g_inv.matrix_at(centre),
            e.source.total_effective_stress_energy.matrix_at(centre), true);
        const auto wf = geometry::wec(centre_state.frame, envelope);
        const auto df = geometry::dec(centre_state.frame, envelope);
        const auto sf = geometry::sec(centre_state.frame, envelope);
        const int cells = wall_cells(grid, gf, c, v_s);

        sweep << frame << "," << v_s << "," << g00_centre << ","
              << n_eul_valid << "," << n_sta_valid << "," << n_com_valid << ","
              << n_bst_valid << "," << 2 << "," << I.min_rho << ","
              << I.min_rho_x << ","
              << I.min_rho_y << "," << I.min_rho_z << "," << rho_max << ","
              << rho_min_slice << "," << rho_min_slice_x << ","
              << rho_min_slice_y << "," << slice_n_min << "," << n_min_slice_x
              << "," << n_min_slice_y << ","
              << I.n_min << "," << I.n_dir[0] << "," << I.n_dir[1] << ","
              << I.n_dir[2] << "," << I.n_min_x << "," << I.n_min_y << ","
              << I.n_min_z << "," << I.e_plus << "," << I.e_minus << ","
              << I.e_net << "," << I.e_abs << "," << I.j_nec << "," << envelope
              << "," << richardson << "," << nec_margin_frame << ","
              << status_name(frame_rho_status) << ","
              << status_name(frame_nec_status) << "," << status_name(wf.status)
              << "," << status_name(df.status) << "," << status_name(sf.status)
              << "," << cert_violated << "," << indeterminate << ","
              << cert_satisfied << "," << e.divergence_relative << ","
              << kGridN << "," << kHalfWidth << "," << cells << ","
              << I.boundary_fraction << "," << (brackets ? 1 : 0) << "\n";

        necsweep << frame << "," << v_s << "," << slice_n_min << ","
                 << slice_worst.direction[0] << "," << slice_worst.direction[1]
                 << "," << slice_worst.direction[2] << ","
                 << slice_worst.stationarity_residual << ","
                 << slice_worst.quadrature_crosscheck << ","
                 << slice_worst.multistart_crosscheck << ","
                 << slice_worst.direct_substitution << ","
                 << slice_worst_eigen_gap << ","
                 << (slice_worst.hard_case ? 1 : 0) << "," << I.j_nec << ","
                 << I.j_nec_smooth << "," << kSoftplusWidth << "," << envelope
                 << "," << nec_margin_frame << ","
                 << status_name(geometry::certify(slice_n_min, envelope)) << ","
                 << nec_violating << "," << slice_cells << "\n";

        // "optimization_ready" is a REFUSABLE claim, not a label: an objective
        // dominated by the numerical boundary is measuring the domain
        // truncation, and an objective whose exterior tail is comparable to its
        // value has not converged. Both disqualify it.
        const bool opt_ready =
            I.boundary_fraction < 0.35 &&
            (I.e_minus > 0.0 && tail_e_minus < 0.25 * I.e_minus);

        integ << frame << "," << v_s << "," << I.e_plus << "," << I.e_minus
              << "," << I.e_net << "," << I.e_abs << "," << I.j_nec << ","
              << I.e_minus << "," << I.e_minus_smooth << "," << kSoftplusWidth
              << "," << I.j_nec_smooth << "," << net_closure << ","
              << abs_closure << "," << I.e_abs_interior << ","
              << I.e_abs_shell << "," << tail_e_minus << "," << tail_e_abs
              << "," << kHalfWidth << "," << kTailHalfWidth << "," << kGridN
              << "," << cells << "," << I.boundary_fraction << ","
              << I.e_plus_interior << "," << I.e_plus_shell << ","
              << I.e_minus_interior << "," << I.e_minus_shell << ","
              << I.j_nec_interior << "," << I.j_nec_shell << ","
              << I.invalid_points << "," << (opt_ready ? 1 : 0) << "\n";

        std::cout << "frame " << frame << "  v_s=" << v_s
                  << "  g00(c)=" << g00_centre << "  min rho=" << I.min_rho
                  << "  N_min=" << I.n_min << "  E_-=" << I.e_minus
                  << "  J_NEC=" << I.j_nec << "  env=" << envelope
                  << "  bianchi_rel=" << e.divergence_relative << "\n";
        ++frame_count;
    }

    if (reference_frame < 0) {
        std::cerr << "the sweep did not contain the reference velocity\n";
        return 1;
    }

    // --- Resolution certification -------------------------------------------
    //
    // At the reference velocity. Reports what each grid CERTIFIES, not merely
    // what it measures: the whole point of the study is that the same physical
    // claim can be certified on one grid and indeterminate on a coarser one.
    {
        std::ofstream fr(out_dir + "/resolution_certification.csv");
        fr << std::setprecision(12);
        fr << "n,h,domain_extent,wall_cells,E_minus,J_NEC,min_rho,N_min,"
              "richardson_difference,numerical_error_estimate,"
              "certification_margin,rho_sign_status,NEC_status,"
              "boundary_fraction,grid_resolution\n";
        double prev_min_rho = 0.0;
        bool have_prev = false;
        for (int n : {13, 17, 21, 25}) {
            const Grid gn = cartesian_grid(n, kHalfWidth);
            geometry::CurvatureFields<double> cn;
            const auto en =
                alcubierre_einstein(gn, cn, kReferenceVelocity, false);
            const Integrals In = integrate(gn, cn, en, true);
            const double diff =
                have_prev ? std::fabs(In.min_rho - prev_min_rho) : 0.0;
            const double env = std::fmax(kEnvelopeFloor, diff);
            metrics::AlcubierreMetric<double> mn(
                0.0, metrics::TimeTreatment::InertialTimeDependent);
            metrics::AlcubierreParameters<double> pn(kReferenceVelocity,
                                                     kRadius, kSigma);
            const Rank2Field<double> gfn =
                metrics::build_metric_field<double>(gn, mn, pn, kTime);
            fr << n << "," << gn.dx() << "," << kHalfWidth << ","
               << wall_cells(gn, gfn, cn, kReferenceVelocity) << ","
               << In.e_minus << "," << In.j_nec << "," << In.min_rho << ","
               << In.n_min << "," << diff << "," << env << ","
               << (std::fabs(In.n_min) -
                   geometry::kCertificationSafetyFactor * env)
               << "," << status_name(geometry::certify(In.min_rho, env)) << ","
               << status_name(geometry::certify(In.n_min, env)) << ","
               << In.boundary_fraction << "," << n << "\n";
            prev_min_rho = In.min_rho;
            have_prev = true;
            std::cout << "resolution n=" << n << "  min rho=" << In.min_rho
                      << "  E_-=" << In.e_minus << "\n";
        }
    }

    // --- Domain convergence --------------------------------------------------
    {
        std::ofstream fd(out_dir + "/domain_convergence.csv");
        fd << std::setprecision(12);
        fd << "half_width,domain_extent,grid_resolution,h,E_plus,E_minus,E_abs,"
              "J_NEC,boundary_fraction,tail_estimate,min_rho,N_min,"
              "wall_cells\n";
        double prev = 0.0;
        bool have_prev = false;
        for (double half : {3.0, 4.0, 5.0}) {
            const Grid gd = cartesian_grid(kGridN, half);
            geometry::CurvatureFields<double> cd;
            const auto ed =
                alcubierre_einstein(gd, cd, kReferenceVelocity, false);
            const Integrals Id = integrate(gd, cd, ed, true);
            const double tail = have_prev ? std::fabs(Id.e_minus - prev) : 0.0;
            metrics::AlcubierreMetric<double> md(
                0.0, metrics::TimeTreatment::InertialTimeDependent);
            metrics::AlcubierreParameters<double> pd(kReferenceVelocity,
                                                     kRadius, kSigma);
            const Rank2Field<double> gfd =
                metrics::build_metric_field<double>(gd, md, pd, kTime);
            fd << half << "," << half << "," << kGridN << "," << gd.dx() << ","
               << Id.e_plus << "," << Id.e_minus << "," << Id.e_abs << ","
               << Id.j_nec << "," << Id.boundary_fraction << "," << tail << ","
               << Id.min_rho << "," << Id.n_min << ","
               << wall_cells(gd, gfd, cd, kReferenceVelocity) << "\n";
            prev = Id.e_minus;
            have_prev = true;
            std::cout << "domain half=" << half << "  E_-=" << Id.e_minus
                      << "  boundary=" << Id.boundary_fraction << "\n";
        }
    }

    // --- N(n) over the whole sphere -----------------------------------------
    //
    // At representative points, so the null-direction figure can show that the
    // minimum is genuinely off-axis rather than assert it. The directions are
    // production's own deterministic Fibonacci lattice and the values are
    // production's nec_value; nothing here searches.
    {
        const Grid gs = cartesian_grid(kGridN, kHalfWidth);
        geometry::CurvatureFields<double> cs;
        const auto es = alcubierre_einstein(gs, cs, kReferenceVelocity, false);

        // Where the negative energy actually is, on the axis of travel where it
        // is not, on the nominal wall, and outside the bubble.
        const Integrals Is = integrate(gs, cs, es, true);
        int mi = im, mj = jm;
        for (int i = 0; i < gs.nx(); ++i)
            if (std::fabs(gs.x(i) - Is.min_rho_x) <
                std::fabs(gs.x(mi) - Is.min_rho_x)) mi = i;
        for (int j = 0; j < gs.ny(); ++j)
            if (std::fabs(gs.y(j) - Is.min_rho_y) <
                std::fabs(gs.y(mj) - Is.min_rho_y)) mj = j;

        struct Probe { const char* label; int i, j; };
        const Probe probes[4] = {
            {"rho_minimum", mi, mj},
            {"axis_of_travel", im, jm},
            {"nominal_wall", im, j_wall},
            {"exterior", 1, 1}};

        std::ofstream fs(out_dir + "/null_direction_sphere.csv");
        fs << std::setprecision(12);
        fs << "point_label,x,y,z,sample_index,nx,ny,nz,theta,phi,N_of_n\n";
        std::ofstream fp(out_dir + "/null_direction_sphere_points.csv");
        fp << std::setprecision(12);
        // The angles of the two marked directions are exported alongside the
        // directions themselves, so the renderer places its markers without
        // performing any trigonometry of its own.
        fp << "point_label,x,y,z,rho,N_min,production_direction_x,"
              "production_direction_y,production_direction_z,"
              "production_theta,production_phi,"
              "quadrature_minimum,quadrature_direction_x,"
              "quadrature_direction_y,quadrature_direction_z,"
              "quadrature_theta,quadrature_phi,"
              "sampled_minimum,nec_stationarity_residual,"
              "nec_crosscheck_residual,nec_direct_substitution,nec_eigen_gap,"
              "nec_hard_case,near_flat_fraction,sample_count\n";

        for (const Probe& pr : probes) {
            const std::size_t p = gs.flat(pr.i, pr.j, km);
            const PointState s = point_state(
                cs.g.matrix_at(p), cs.g_inv.matrix_at(p),
                es.source.total_effective_stress_energy.matrix_at(p), true);
            if (!s.ok) continue;

            double best = 1e300, bestn[3] = {0, 0, 0};
            int near_flat = 0;
            std::vector<double> vals(static_cast<std::size_t>(kSphereSamples));
            std::vector<double> dirs(static_cast<std::size_t>(kSphereSamples) * 3);
            for (int q = 0; q < kSphereSamples; ++q) {
                double n[3];
                geometry::detail::fibonacci_direction(q, kSphereSamples, n);
                const double v = geometry::nec_value(s.frame, n);
                vals[static_cast<std::size_t>(q)] = v;
                for (int d = 0; d < 3; ++d)
                    dirs[static_cast<std::size_t>(q) * 3 +
                         static_cast<std::size_t>(d)] = n[d];
                if (v < best) {
                    best = v;
                    for (int d = 0; d < 3; ++d) bestn[d] = n[d];
                }
            }
            const double span = std::fabs(s.nec.n_min) > 0.0
                                    ? std::fabs(s.nec.n_min)
                                    : 1.0;
            for (int q = 0; q < kSphereSamples; ++q) {
                const double v = vals[static_cast<std::size_t>(q)];
                if (std::fabs(v - s.nec.n_min) < 0.01 * span) ++near_flat;
                const double nx = dirs[static_cast<std::size_t>(q) * 3 + 0];
                const double ny = dirs[static_cast<std::size_t>(q) * 3 + 1];
                const double nz = dirs[static_cast<std::size_t>(q) * 3 + 2];
                // Spherical angles of a direction production already produced.
                // Presentation geometry, exported here rather than in the
                // renderer so no trigonometry at all lives on the Python side.
                fs << pr.label << "," << gs.x(pr.i) << "," << gs.y(pr.j) << ","
                   << gs.z(km) << "," << q << "," << nx << "," << ny << ","
                   << nz << "," << std::acos(nz < -1.0 ? -1.0
                                             : (nz > 1.0 ? 1.0 : nz))
                   << "," << std::atan2(ny, nx) << ","
                   << vals[static_cast<std::size_t>(q)] << "\n";
            }

            auto polar = [](const double n[3], double& th, double& ph) {
                const double c = n[2] < -1.0 ? -1.0 : (n[2] > 1.0 ? 1.0 : n[2]);
                th = std::acos(c);
                ph = std::atan2(n[1], n[0]);
            };
            double pth = 0.0, pph = 0.0, qth = 0.0, qph = 0.0;
            polar(s.nec.direction, pth, pph);
            polar(bestn, qth, qph);

            fp << pr.label << "," << gs.x(pr.i) << "," << gs.y(pr.j) << ","
               << gs.z(km) << "," << s.frame.rho << "," << s.nec.n_min << ","
               << s.nec.direction[0] << "," << s.nec.direction[1] << ","
               << s.nec.direction[2] << "," << pth << "," << pph << ","
               << best << "," << bestn[0] << ","
               << bestn[1] << "," << bestn[2] << "," << qth << "," << qph
               << "," << best << ","
               << s.nec.stationarity_residual << ","
               << s.nec.quadrature_crosscheck << ","
               << s.nec.direct_substitution << "," << s.nec_eigen_gap << ","
               << (s.nec.hard_case ? 1 : 0) << ","
               << static_cast<double>(near_flat) /
                      static_cast<double>(kSphereSamples)
               << "," << kSphereSamples << "\n";
        }
    }

    // --- Manifest -------------------------------------------------------------
    {
        std::ofstream man(out_dir + "/a6_energy_manifest.json");
        man << std::setprecision(12);
        man << "{\n";
        man << "  \"git_commit\": \"" << build_info().git_commit << "\",\n";
        man << "  \"git_dirty\": " << (build_info().git_dirty ? "true" : "false")
            << ",\n";
        man << "  \"generated_utc\": \"" << utc_timestamp() << "\",\n";
        man << "  \"milestone\": \"A6\",\n";
        man << "  \"grid\": \"" << kGridN << "^3 over [-" << kHalfWidth << ","
            << kHalfWidth << "]^3\",\n";
        man << "  \"grid_resolution\": " << kGridN << ",\n";
        man << "  \"domain_extent\": " << kHalfWidth << ",\n";
        man << "  \"tail_domain_extent\": " << kTailHalfWidth << ",\n";
        man << "  \"coarse_grid_resolution\": " << kCoarseN << ",\n";
        man << "  \"slice\": \"z = 0\",\n";
        man << "  \"radius\": " << kRadius << ",\n";
        man << "  \"sigma\": " << kSigma << ",\n";
        man << "  \"reference_velocity\": " << kReferenceVelocity << ",\n";
        man << "  \"reference_frame\": " << reference_frame << ",\n";
        man << "  \"frame_count\": " << frame_count << ",\n";
        man << "  \"frame_rows\": " << total_frame_rows << ",\n";
        man << "  \"velocity_sequence\": \"v_s = k/10 for k = " << kFirstK
            << ".." << kLastK << "\",\n";
        man << "  \"exact_onset_frame\": " << exact_onset_frame << ",\n";
        man << "  \"first_spatial_contour_frame\": " << first_contour_frame
            << ",\n";
        man << "  \"numerical_error_envelope_floor\": " << kEnvelopeFloor
            << ",\n";
        man << "  \"certification_safety_factor\": "
            << geometry::kCertificationSafetyFactor << ",\n";
        man << "  \"softplus_width\": " << kSoftplusWidth << ",\n";
        man << "  \"nec_quadrature_points\": " << kNecQuadrature << ",\n";
        man << "  \"sphere_sample_count\": " << kSphereSamples << ",\n";
        man << "  \"boundary_shell_depth\": " << kShell << ",\n";
        man << "  \"worst_eulerian_normalization_residual\": "
            << worst_eulerian_residual << ",\n";
        man << "  \"worst_tetrad_orthonormality_residual\": "
            << worst_tetrad_residual << ",\n";
        man << "  \"tetrad_handedness_range\": [" << handedness_min << ", "
            << handedness_max << "],\n";
        man << "  \"tetrad_min_diagonal\": " << tetrad_diag_min << ",\n";
        man << "  \"worst_net_closure_residual\": " << worst_net_closure
            << ",\n";
        man << "  \"worst_abs_closure_residual\": " << worst_abs_closure
            << ",\n";
        man << "  \"all_finite\": " << (all_finite ? "true" : "false") << ",\n";
        man << "  \"observer_families\": [\"eulerian\", "
               "\"coordinate_stationary\", \"bubble_comoving\", "
               "\"tetrad_boosted\"],\n";
        man << "  \"observer_admissibility\": \"The coordinate-stationary "
               "observer is admissible ONLY where g_00 < 0. At g_00 = 0 the "
               "coordinate vector d/dt is NULL and for g_00 > 0 it is "
               "SPACELIKE; in both cases it is rejected with a reason and its "
               "energy-density cell is left EMPTY. |g_00| is never substituted. "
               "The Eulerian observer is built from the ADM normal and remains "
               "timelike through the crossing; it is recorded separately.\",\n";
        man << "  \"g00_zero_meaning\": \"g_00 = 0 is where the coordinate "
               "vector d/dt becomes null: an ergosurface-type condition on the "
               "TIME COORDINATE. It is not a metric degeneracy, not a signature "
               "failure, not a curvature singularity and not automatically a "
               "horizon. det(g) and the spatial metric remain regular there and "
               "the ADM foliation persists.\",\n";
        man << "  \"exact_vs_smooth\": \"E_minus and J_NEC are the EXACT "
               "objectives, integrals of max(0, -rho) and max(0, -N_min). The "
               "softplus columns are a smooth surrogate for a future "
               "gradient-based optimizer; softplus systematically overestimates "
               "max(0,x) near zero, which is exactly where certification "
               "happens, so it is never substituted for the exact value.\",\n";
        man << "  \"bianchi_note\": \"bianchi_relative_residual is the A5 "
               "contracted-Bianchi quality metric, carried as a NUMERICAL "
               "DIAGNOSTIC. It is not evidence that local conservation is well "
               "resolved.\",\n";
        man << "  \"quantity\": \"observer-measured energy density, momentum "
               "density, principal stresses and pointwise energy-condition "
               "certification in the orthonormal frame of a named observer\",\n";
        man << "  \"renderer_computes_physics\": false,\n";
        man << "  \"scope\": \"A6 evaluates observer-dependent energy density "
               "and classical pointwise energy conditions for the "
               "stress-energy tensor inferred from the selected metric. It does "
               "not demonstrate that the source can be generated, that quantum "
               "energy inequalities are satisfied, that the geometry is "
               "dynamically stable, or that operational FTL travel is "
               "possible.\"\n";
        man << "}\n";
    }

    std::cout << "\nframes: " << frame_count
              << "\nexact_onset_frame: " << exact_onset_frame
              << "\nfirst_spatial_contour_frame: " << first_contour_frame
              << "\nworst eulerian normalisation residual: "
              << worst_eulerian_residual
              << "\nworst tetrad orthonormality residual: "
              << worst_tetrad_residual
              << "\nworst E_net closure residual: " << worst_net_closure
              << "\nworst E_abs closure residual: " << worst_abs_closure
              << "\nall finite: " << (all_finite ? "yes" : "NO") << "\n";

    if (exact_onset_frame != 8) {
        std::cerr << "exact_onset_frame is " << exact_onset_frame
                  << ", expected 8 (v_s = 1 exactly)\n";
        return 1;
    }
    if (first_contour_frame != 9) {
        std::cerr << "first_spatial_contour_frame is " << first_contour_frame
                  << ", expected 9\n";
        return 1;
    }
    if (!all_finite) {
        std::cerr << "non-finite observer-frame values encountered\n";
        return 1;
    }
    if (handedness_min <= 0.0 || tetrad_diag_min <= 0.0) {
        std::cerr << "tetrad handedness or positive-diagonal convention "
                     "violated\n";
        return 1;
    }
    return 0;
}
