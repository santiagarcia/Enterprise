#include "o0c_checks.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// O0C cases 18-26: observer-frame quantities, energy conditions, the integrated
// objectives and their sign conventions, and the smooth surrogate.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using design::AxisymmetricBasis;
using design::IntegrationSet;
using design::MetricJet;
using design::ShiftFamilyConfig;
using design::SourceSemantics;

const double kPi = 3.14159265358979323846;

struct Site { const char* label; double x, y, z; };

std::vector<Site> observer_sites() {
    return {
        {"cabin_centre", 0.0, 0.0, 0.0},
        {"cabin_offaxis", 0.0, 0.7, 0.0},
        {"axis_front_wall", 2.0, 0.0, 0.0},
        {"axis_back_wall", -2.0, 0.0, 0.0},
        {"wall_transverse", 0.0, 2.0, 0.0},
        {"wall_off_axis", 0.9, 1.6, 0.0},
        {"wall_off_axis_z", 0.9, 0.9, 1.4},
        {"exterior", 0.0, 3.6, 0.0},
    };
}

}  // namespace

void run_o0c_energy_checks(VerificationReport& report,
                           const std::string& out_dir) {
    Stopwatch sw;
    const ShiftFamilyConfig cfg = o0c_config();
    AxisymmetricBasis f1 = o0c_f1_basis();
    const std::vector<double> zero(
        static_cast<std::size_t>(f1.coefficient_count()), 0.0);

    // --- 18. The Eulerian observer, and the scalar-generic rho ---------------
    sw.reset();
    {
        std::ofstream f(out_dir + "/observer_frame_samples.csv");
        f.precision(17);
        f << "site,x,y,z,g00,observer,valid,rejection_reason,norm_residual,"
             "tetrad_orthonormality,sqrt_gamma,rho,q_x,q_y,q_z,"
             "principal_1,principal_2,principal_3,trace_spatial,"
             "scalar_generic_rho,scalar_generic_residual\n";

        double worst_norm = 0.0, worst_ortho = 0.0, worst_generic = 0.0;
        int valid_eulerian = 0;
        for (const Site& s : observer_sites()) {
            const MetricJet<double> jet =
                design::f1_metric_jet(cfg, f1, zero, 0.0, s.x, s.y, s.z);
            const design::PointGeometry pg = design::point_geometry(jet);
            const Mat4<double> T = design::stress_energy_from_einstein(
                pg.einstein, pg.g, 0.0, SourceSemantics::TotalEffective);
            const geometry::AdmDecomposition adm =
                geometry::adm_decompose(jet.g_cov.m, jet.g_contra.m);
            const geometry::Observer eul =
                geometry::make_eulerian(jet.g_cov.m, adm);
            const geometry::Tetrad tet =
                geometry::build_tetrad(jet.g_cov.m, adm);
            const geometry::FrameStress fs =
                geometry::project_to_frame(tet, T);
            const double gen = design::eulerian_density<double>(
                cfg, f1, zero, 0.0, s.x, s.y, s.z, 0.0,
                SourceSemantics::TotalEffective);

            if (eul.valid) ++valid_eulerian;
            worst_norm = std::fmax(worst_norm, eul.normalization_residual);
            worst_ortho = std::fmax(worst_ortho, tet.orthonormality_residual);
            worst_generic = std::fmax(worst_generic, std::fabs(gen - fs.rho));

            f << s.label << "," << s.x << "," << s.y << "," << s.z << ","
              << jet.g_cov.m(0, 0) << ",eulerian," << (eul.valid ? 1 : 0)
              << ",\"" << eul.rejection_reason << "\","
              << eul.normalization_residual << ","
              << tet.orthonormality_residual << ",1," << fs.rho << ","
              << fs.q[0] << "," << fs.q[1] << "," << fs.q[2] << ","
              << fs.principal[0] << "," << fs.principal[1] << ","
              << fs.principal[2] << "," << fs.trace_spatial << "," << gen
              << "," << std::fabs(gen - fs.rho) << "\n";
        }

        report.check_predicate(
            "o0c.eulerian_observer_valid",
            "the ADM Eulerian observer exists and is unit-timelike at every "
            "sample point, including through the g_00 crossing",
            VerificationLevel::A,
            valid_eulerian == static_cast<int>(observer_sites().size()) &&
                worst_norm < 1e-12 && worst_ortho < 1e-12,
            "valid everywhere, |g(n,n)+1| and tetrad residuals at roundoff",
            std::to_string(valid_eulerian) + " of " +
                std::to_string(observer_sites().size()) +
                " valid; worst normalisation " + std::to_string(worst_norm) +
                ", worst tetrad orthonormality " + std::to_string(worst_ortho),
            sw.seconds(),
            "n^mu = (1/alpha, -beta^i/alpha) is timelike wherever the foliation "
            "exists, which is why it and not d/dt is the energy observer");
        sw.reset();

        report.check_residual(
            "o0c.scalar_generic_density_matches_a6",
            "the differentiable rho used by the optimizer path reproduces A6's "
            "make_eulerian + project_to_frame exactly",
            VerificationLevel::A, "8 sample points", "Lambda = 0, total "
            "effective", worst_generic, 1.0e-18, sw.seconds(),
            "this agreement is what licenses differentiating the compact form "
            "instead of the tetrad construction, which branches on eigenvalues");
        sw.reset();
    }

    // --- 19. The coordinate-stationary observer must be REJECTED -------------
    {
        int rejected = 0, admitted = 0, wrongly_normalised = 0;
        std::string sample_reason;
        std::ofstream f(out_dir + "/observer_frame_samples.csv",
                        std::ios::app);
        f.precision(17);
        // Sweep the velocity so g_00 = -1 + v^2 f^2 crosses zero on the axis.
        for (double v : {0.5, 0.9, 1.0, 1.1, 1.6}) {
            ShiftFamilyConfig c = cfg;
            c.velocity = v;
            for (const Site& s : observer_sites()) {
                const MetricJet<double> jet =
                    design::f1_metric_jet(c, f1, zero, 0.0, s.x, s.y, s.z);
                const geometry::Observer cs =
                    geometry::make_coordinate_stationary(jet.g_cov.m);
                const double g00 = jet.g_cov.m(0, 0);
                if (g00 >= 0.0) {
                    if (cs.valid) ++wrongly_normalised; else ++rejected;
                    if (sample_reason.empty() && !cs.valid)
                        sample_reason = cs.rejection_reason;
                } else if (cs.valid) {
                    ++admitted;
                }
                f << s.label << ",;v=" << v << "," << s.y << "," << s.z << ","
                  << g00 << ",coordinate_stationary," << (cs.valid ? 1 : 0)
                  << ",\"" << cs.rejection_reason << "\","
                  << cs.normalization_residual << ",0,1,,,,,,,,,\n";
            }
        }
        report.check_predicate(
            "o0c.invalid_observer_rejected",
            "the coordinate-stationary observer is rejected wherever "
            "g_00 >= 0, and |g_00| is never substituted to keep it alive",
            VerificationLevel::A, wrongly_normalised == 0 && rejected > 0,
            "zero invalid observers normalised, and the rejection actually "
            "exercised",
            std::to_string(rejected) + " rejected where g_00 >= 0, " +
                std::to_string(admitted) + " admitted where g_00 < 0, " +
                std::to_string(wrongly_normalised) +
                " wrongly normalised; example reason: " + sample_reason,
            sw.seconds(),
            "at g_00 = 0 the coordinate vector d/dt is NULL and for g_00 > 0 it "
            "is SPACELIKE; neither is an observer");
        sw.reset();
    }

    // --- 20. The other observers --------------------------------------------
    {
        int comoving_valid = 0, comoving_rejected = 0, boosted_valid = 0;
        std::string comoving_reason;
        for (double v : {0.5, 0.9, 1.4}) {
            ShiftFamilyConfig c = cfg;
            c.velocity = v;
            for (const Site& s : observer_sites()) {
                const MetricJet<double> jet =
                    design::f1_metric_jet(c, f1, zero, 0.0, s.x, s.y, s.z);
                const geometry::Observer bc =
                    geometry::make_bubble_comoving(jet.g_cov.m, v);
                if (bc.valid) ++comoving_valid;
                else { ++comoving_rejected;
                       if (comoving_reason.empty())
                           comoving_reason = bc.rejection_reason; }
                const geometry::AdmDecomposition adm =
                    geometry::adm_decompose(jet.g_cov.m, jet.g_contra.m);
                const geometry::Tetrad tet =
                    geometry::build_tetrad(jet.g_cov.m, adm);
                const double vv[3] = {0.3, 0.2, -0.1};
                const geometry::Observer tb =
                    geometry::make_tetrad_boosted(jet.g_cov.m, tet, vv);
                if (tb.valid) ++boosted_valid;
            }
        }
        report.check_predicate(
            "o0c.additional_observers_evaluated",
            "the bubble-comoving and tetrad-boosted observers are evaluated, "
            "with timelikeness checked rather than assumed",
            VerificationLevel::B, boosted_valid > 0,
            "boosted observers valid; comoving validity reported either way",
            "bubble-comoving valid at " + std::to_string(comoving_valid) +
                " sites, rejected at " + std::to_string(comoving_rejected) +
                (comoving_reason.empty() ? "" : " (" + comoving_reason + ")") +
                "; tetrad-boosted valid at " + std::to_string(boosted_valid),
            sw.seconds(),
            "the bubble-centre coordinate trajectory is not automatically a "
            "worldline, and for a superluminal bubble in the exterior it is not "
            "one -- a rejection here is the correct physics, not a failure");
        sw.reset();
    }

    // --- 21. The analytic Alcubierre density --------------------------------
    {
        std::ofstream f(out_dir + "/alcubierre_density_reference.csv");
        f.precision(17);
        f << "site,x,y,z,rho_computed,rho_analytic,absolute_difference,"
             "relative_difference\n";
        double worst_abs = 0.0, worst_rel = 0.0;
        int nonzero_sites = 0;
        for (const Site& s : observer_sites()) {
            const double rho = design::eulerian_density<double>(
                cfg, f1, zero, 0.0, s.x, s.y, s.z, 0.0,
                SourceSemantics::TotalEffective);
            const design::ShiftJet sj =
                design::shift_jet(cfg, f1, zero, 0.0, s.x, s.y, s.z, false);
            const double fy = -sj.d1[2] / cfg.velocity;
            const double fz = -sj.d1[3] / cfg.velocity;
            const double ref = -(cfg.velocity * cfg.velocity) / (32.0 * kPi) *
                               (fy * fy + fz * fz);
            const double da = std::fabs(rho - ref);
            worst_abs = std::fmax(worst_abs, da);
            if (std::fabs(ref) > 1e-12) {
                ++nonzero_sites;
                worst_rel = std::fmax(worst_rel, da / std::fabs(ref));
            }
            f << s.label << "," << s.x << "," << s.y << "," << s.z << "," << rho
              << "," << ref << "," << da << ","
              << (std::fabs(ref) > 1e-12 ? da / std::fabs(ref) : 0.0) << "\n";
        }
        report.check_residual(
            "o0c.analytic_alcubierre_density",
            "rho_Eulerian = -(v_s^2/32 pi) [ (d_y f)^2 + (d_z f)^2 ] for the "
            "shift-only baseline",
            VerificationLevel::A,
            "8 sites: axis, cabin, wall on three axes, exterior",
            "v_s = 0.6, R = 2, sigma = 1", worst_abs, 1.0e-16, sw.seconds(),
            "worst relative over the " + std::to_string(nonzero_sites) +
                " sites where the reference is nonzero: " +
                std::to_string(worst_rel) +
                "; the axis and cabin sites have an exactly zero reference "
                "because d_y f = d_z f = 0 there, so they are checked "
                "absolutely");
        sw.reset();
    }

    // --- 22. v_s^2 scaling ---------------------------------------------------
    {
        std::ofstream f(out_dir + "/velocity_scaling.csv");
        f.precision(17);
        f << "v_s,site,rho,rho_over_v_squared,ratio_to_first\n";
        double worst_dev = 0.0;
        for (const Site& s : {observer_sites()[5], observer_sites()[4]}) {
            double base = 0.0;
            bool first = true;
            for (double v : {0.2, 0.4, 0.6, 0.8, 1.2, 1.6}) {
                ShiftFamilyConfig c = cfg;
                c.velocity = v;
                const double rho = design::eulerian_density<double>(
                    c, f1, zero, 0.0, s.x, s.y, s.z, 0.0,
                    SourceSemantics::TotalEffective);
                const double scaled = rho / (v * v);
                if (first) { base = scaled; first = false; }
                const double ratio = base != 0.0 ? scaled / base : 1.0;
                worst_dev = std::fmax(worst_dev, std::fabs(ratio - 1.0));
                f << v << "," << s.label << "," << rho << "," << scaled << ","
                  << ratio << "\n";
            }
        }
        report.check_residual(
            "o0c.velocity_squared_scaling",
            "the Eulerian density scales exactly as v_s^2 at a fixed geometric "
            "point",
            VerificationLevel::A, "6 velocities x 2 wall sites",
            "v_s in {0.2, 0.4, 0.6, 0.8, 1.2, 1.6}", worst_dev, 1.0e-14,
            sw.seconds(),
            "the shape function carries no velocity dependence, so rho/v_s^2 "
            "must be velocity-independent; a spurious v-dependence anywhere in "
            "the chain would show here");
        sw.reset();
    }

    // --- 23/24. Energy conditions and the exact null minimum ----------------
    {
        std::ofstream f(out_dir + "/energy_condition_samples.csv");
        f.precision(17);
        f << "site,x,y,z,rho,N_min,direction_x,direction_y,direction_z,"
             "stationarity_residual,quadrature_crosscheck,direct_substitution,"
             "hard_case,solver_ok,NEC_status,WEC_status,DEC_status,SEC_status\n";
        std::ofstream n(out_dir + "/null_minimizer_checks.csv");
        n.precision(17);
        n << "site,n_min,independent_quadrature_value,difference,"
             "multistart_difference,degenerate_or_nonunique,solver_status\n";

        double worst_cross = 0.0;
        int violated = 0, satisfied = 0, indeterminate = 0, invalid = 0;
        int degenerate = 0;
        for (const Site& s : observer_sites()) {
            const MetricJet<double> jet =
                design::f1_metric_jet(cfg, f1, zero, 0.0, s.x, s.y, s.z);
            const design::PointGeometry pg = design::point_geometry(jet);
            const Mat4<double> T = design::stress_energy_from_einstein(
                pg.einstein, pg.g, 0.0, SourceSemantics::TotalEffective);
            const geometry::AdmDecomposition adm =
                geometry::adm_decompose(jet.g_cov.m, jet.g_contra.m);
            const geometry::Tetrad tet =
                geometry::build_tetrad(jet.g_cov.m, adm);
            if (!tet.valid) { ++invalid; continue; }
            const geometry::FrameStress fs =
                geometry::project_to_frame(tet, T);
            const geometry::NecResult nec = geometry::nec_minimum(fs, 4096);

            // The numerical envelope: this route takes no differences, so the
            // envelope is roundoff scaled by the local magnitude, not a
            // truncation estimate.
            const double scale =
                std::fmax(std::fabs(fs.rho),
                          std::fmax(std::fabs(fs.principal[0]),
                                    std::fabs(fs.principal[2])));
            const double envelope = std::fmax(1e-18, 1e-13 * scale);
            const geometry::ConditionStatus nec_status =
                geometry::certify(nec.n_min, envelope);
            const geometry::ConditionResult w = geometry::wec(fs, envelope);
            const geometry::ConditionResult d = geometry::dec(fs, envelope);
            const geometry::ConditionResult sc = geometry::sec(fs, envelope);

            worst_cross = std::fmax(worst_cross,
                                    std::fabs(nec.quadrature_crosscheck));
            if (nec.hard_case) ++degenerate;
            if (nec_status == geometry::ConditionStatus::CertifiedViolated)
                ++violated;
            else if (nec_status == geometry::ConditionStatus::CertifiedSatisfied)
                ++satisfied;
            else if (nec_status ==
                     geometry::ConditionStatus::IndeterminateWithinDiscretization)
                ++indeterminate;

            f << s.label << "," << s.x << "," << s.y << "," << s.z << ","
              << fs.rho << "," << nec.n_min << "," << nec.direction[0] << ","
              << nec.direction[1] << "," << nec.direction[2] << ","
              << nec.stationarity_residual << "," << nec.quadrature_crosscheck
              << "," << nec.direct_substitution << "," << (nec.hard_case ? 1 : 0)
              << "," << (nec.ok ? 1 : 0) << ","
              << geometry::to_string(nec_status) << ","
              << geometry::to_string(w.status) << ","
              << geometry::to_string(d.status) << ","
              << geometry::to_string(sc.status) << "\n";
            n << s.label << "," << nec.n_min << ","
              << (nec.n_min + nec.quadrature_crosscheck) << ","
              << nec.quadrature_crosscheck << "," << nec.multistart_crosscheck
              << "," << (nec.hard_case ? 1 : 0) << ","
              << (nec.ok ? "ok" : "failed") << "\n";
        }

        // THE QUADRATURE GAP IS NOT AN ERROR, AND TESTING IT AGAINST A FIXED
        // TOLERANCE WOULD BE THE WRONG STATEMENT.
        //
        // The angular quadrature evaluates N at a finite set of directions, so
        // its best value is an UPPER bound on the true minimum: the gap must be
        // >= 0 by construction, and it measures the sampling resolution rather
        // than the solver. Near a nondegenerate minimum N is quadratic, so the
        // gap falls like h^2 with h the angular spacing, i.e. like 1/M for M
        // Fibonacci directions. Two things are therefore checked, and neither
        // is a magic number:
        //
        //   the solver's own direction reproduces its reported minimum EXACTLY
        //   when N is recomputed from scratch there (direct_substitution), and
        //
        //   the gap is non-negative and shrinks like 1/M when M quadruples.
        double worst_direct = 0.0;
        double gap_coarse = 0.0, gap_fine = 0.0;
        bool gap_nonnegative = true;
        for (const Site& s : observer_sites()) {
            const MetricJet<double> jet =
                design::f1_metric_jet(cfg, f1, zero, 0.0, s.x, s.y, s.z);
            const design::PointGeometry pg = design::point_geometry(jet);
            const Mat4<double> T = design::stress_energy_from_einstein(
                pg.einstein, pg.g, 0.0, SourceSemantics::TotalEffective);
            const geometry::AdmDecomposition adm =
                geometry::adm_decompose(jet.g_cov.m, jet.g_contra.m);
            const geometry::Tetrad tet =
                geometry::build_tetrad(jet.g_cov.m, adm);
            if (!tet.valid) continue;
            const geometry::FrameStress fs =
                geometry::project_to_frame(tet, T);
            const geometry::NecResult a = geometry::nec_minimum(fs, 1024);
            const geometry::NecResult b = geometry::nec_minimum(fs, 4096);
            if (!a.ok || !b.ok) continue;
            worst_direct = std::fmax(
                worst_direct, std::fabs(b.direct_substitution - b.n_min));
            if (a.quadrature_crosscheck < -1e-15 ||
                b.quadrature_crosscheck < -1e-15)
                gap_nonnegative = false;
            gap_coarse = std::fmax(gap_coarse, a.quadrature_crosscheck);
            gap_fine = std::fmax(gap_fine, b.quadrature_crosscheck);
        }
        const double shrink =
            gap_fine > 0.0 ? gap_coarse / gap_fine : 0.0;

        report.check_residual(
            "o0c.exact_nec_minimisation",
            "N recomputed from scratch at the solver's minimising direction "
            "reproduces the reported minimum exactly",
            VerificationLevel::A, "8 sites", "N(n) = rho - 2 q.n + n^T S n",
            worst_direct, 1.0e-15, sw.seconds(),
            "this is the exactness statement; the angular quadrature is a "
            "separate resolution check below. " + std::to_string(degenerate) +
                " sites reported the degenerate hard case");
        sw.reset();

        report.check_predicate(
            "o0c.nec_quadrature_crosscheck_resolution",
            "the independent angular quadrature bounds the minimum from above "
            "and its gap shrinks like the direction count",
            VerificationLevel::B,
            gap_nonnegative && (gap_fine == 0.0 || shrink > 1.5),
            "gap >= 0 always, and shrinking by more than 1.5x when the "
            "direction count quadruples",
            "worst gap " + std::to_string(gap_coarse) + " at 1024 directions, " +
                std::to_string(gap_fine) + " at 4096, ratio " +
                std::to_string(shrink) + "; non-negative everywhere: " +
                (gap_nonnegative ? "yes" : "NO"), sw.seconds(),
            "a fixed tolerance here would be testing the sphere sampling, not "
            "the minimiser: the gap is an upper-bound artefact that must be "
            "positive, and what matters is that it converges");
        sw.reset();
        (void)worst_cross;

        report.check_predicate(
            "o0c.energy_condition_classification",
            "WEC, DEC, SEC and the NEC are reported with an explicit "
            "indeterminate category, never collapsed into a verdict",
            VerificationLevel::A,
            violated + satisfied + indeterminate + invalid ==
                static_cast<int>(observer_sites().size()),
            "every site lands in exactly one of four categories",
            std::to_string(violated) + " certified violated, " +
                std::to_string(satisfied) + " certified satisfied, " +
                std::to_string(indeterminate) + " indeterminate, " +
                std::to_string(invalid) + " observer/tetrad invalid",
            sw.seconds(),
            "an indeterminate value is never turned into a violation or a "
            "satisfaction; the axis and cabin sites are exactly the ones where "
            "the density is zero to roundoff and the honest answer is "
            "'indeterminate'");
        sw.reset();
    }

    // --- 25. Objective definitions, sign conventions, closures ---------------
    {
        const design::DesignQuadrature q = design::build_quadrature(f1);
        IntegrationSet set;
        set.time = 0.0;
        set.route = "elementwise_analytic_jet";
        set.quadrature_order = 8;
        set.subdivision = 1;
        set.xi_extent = o0c_f1_region().xi_max;
        set.s_extent = o0c_f1_region().s_max;
        for (std::size_t i = 0; i < q.point.size(); ++i) {
            set.x.push_back(q.point[i].xi * cfg.scale_R0);
            set.y.push_back(q.point[i].y);
            set.z.push_back(q.point[i].z);
            set.weight.push_back(q.weight[i]);
        }
        const design::ObjectiveValues v = design::evaluate_objectives(
            cfg, f1, zero, set, 0.0, SourceSemantics::TotalEffective,
            /*want_nec=*/false, /*want_kretschmann=*/true);

        std::ofstream f(out_dir + "/objective_definitions.json");
        f.precision(17);
        f << "{\n"
          << "  \"E_net\": \"integral rho dV_gamma, SIGNED\",\n"
          << "  \"E_plus\": \"integral max(rho, 0) dV_gamma, positive\",\n"
          << "  \"E_neg\": \"integral max(-rho, 0) dV_gamma, the POSITIVE "
             "MAGNITUDE of the negative-energy content\",\n"
          << "  \"E_abs\": \"integral |rho| dV_gamma = E_plus + E_neg\",\n"
          << "  \"J_NEC\": \"integral max(-N_min, 0) dV_gamma\",\n"
          << "  \"identity\": \"E_net = E_plus - E_neg\",\n"
          << "  \"E_minus_is_not_used\": \"A5/A6 wrote E_minus for what is here "
             "E_neg, also as a positive magnitude. O0C does not use the name "
             "E_minus at all, so no reader has to decide whether it is the "
             "signed integral or the magnitude.\",\n"
          << "  \"volume_element\": \"dV_gamma = sqrt(gamma) d^3x; in the "
             "axisymmetric design chart d^3x = pi R0^3 dxi ds, so dV_gamma = "
             "pi R0^3 sqrt(gamma) dxi ds\",\n"
          << "  \"sqrt_gamma_range\": [" << v.sqrt_gamma_min << ", "
          << v.sqrt_gamma_max << "],\n"
          << "  \"source_semantics\": \""
          << design::to_string(v.semantics) << "\",\n"
          << "  \"lambda\": " << v.lambda << "\n}\n";

        std::ofstream ov(out_dir + "/objective_values.csv");
        ov.precision(17);
        ov << "family,route,order,subdivision,xi_extent,s_extent,points,"
              "E_net,E_plus,E_neg,E_abs,rho_min,rho_max,K_max,volume,"
              "net_closure_residual,abs_closure_residual,invalid_points\n";
        ov << "F1_zero_correction," << set.route << ",8,1," << set.xi_extent
           << "," << set.s_extent << "," << v.points << "," << v.E_net << ","
           << v.E_plus << "," << v.E_neg << "," << v.E_abs << "," << v.rho_min
           << "," << v.rho_max << "," << v.K_max << "," << v.volume << ","
           << v.net_closure_residual << "," << v.abs_closure_residual << ","
           << v.invalid_points << "\n";

        report.check_predicate(
            "o0c.objective_sign_conventions",
            "E_neg is a positive magnitude, E_net is signed, and the name "
            "E_minus is not used at all",
            VerificationLevel::A,
            v.E_neg >= 0.0 && v.E_plus >= 0.0 && v.E_net <= 0.0,
            "E_neg >= 0, E_plus >= 0, and for the Alcubierre baseline "
            "E_net <= 0",
            "E_net " + std::to_string(v.E_net) + ", E_plus " +
                std::to_string(v.E_plus) + ", E_neg " + std::to_string(v.E_neg) +
                " (positive magnitude), E_abs " + std::to_string(v.E_abs),
            sw.seconds(),
            "E_plus is exactly zero for the shift-only baseline because "
            "rho = -(v^2/32 pi)[(d_y f)^2 + (d_z f)^2] is nonpositive "
            "everywhere -- a property of the model, not of the integrator");
        sw.reset();

        report.check_residual(
            "o0c.energy_decomposition_identity",
            "E_net = E_plus - E_neg and E_abs = E_plus + E_neg",
            VerificationLevel::A,
            std::to_string(v.points) + " quadrature points",
            "elementwise analytic-jet route",
            std::fmax(v.net_closure_residual, v.abs_closure_residual), 1.0e-14,
            sw.seconds(), "an accumulation identity, so it is exact rather than "
                          "convergent");
        sw.reset();

        // The volume-element identity: integrating 1 must give the exact
        // pi R0^3 * (xi range) * (s range), which is what proves the weights
        // carry the measure the objectives claim.
        const design::DesignRegion reg = o0c_f1_region();
        const double exact_vol = kPi * cfg.scale_R0 * cfg.scale_R0 *
                                 cfg.scale_R0 * (reg.xi_max - reg.xi_min) *
                                 reg.s_max;
        report.check_residual(
            "o0c.volume_element_identity",
            "integrating dV_gamma over the design region reproduces "
            "pi R0^3 (xi range)(s range) exactly, and sqrt(gamma) = 1 for a "
            "shift-only family",
            VerificationLevel::A, "the same weights the objectives use",
            "pi R0^3 = " + std::to_string(kPi * 8.0),
            std::fmax(std::fabs(v.volume - exact_vol),
                      std::fabs(v.sqrt_gamma_max - 1.0)),
            1.0e-9, sw.seconds(),
            "computed " + std::to_string(v.volume) + " against " +
                std::to_string(exact_vol) + "; sqrt(gamma) in [" +
                std::to_string(v.sqrt_gamma_min) + ", " +
                std::to_string(v.sqrt_gamma_max) +
                "] -- measured, because the identity stops holding the moment a "
                "spatial-metric family is added");
        sw.reset();
    }

    // --- 26. The smooth surrogate --------------------------------------------
    {
        double worst_bias_error = 0.0;
        std::ofstream f(out_dir + "/smooth_surrogate_bias.csv");
        f.precision(17);
        f << "tau,softplus_at_zero,tau_log2,bias_error,softplus_at_one,"
             "softplus_at_minus_one,overestimates_max\n";
        bool always_over = true;
        for (double tau : {1e-1, 1e-2, 1e-3, 1e-4, 1e-5}) {
            const double at0 = design::softplus<double>(0.0, tau);
            const double want = design::softplus_bias_at_zero(tau);
            worst_bias_error = std::fmax(worst_bias_error,
                                         std::fabs(at0 - want));
            for (double u : {-3.0, -0.5, 0.0, 0.5, 3.0}) {
                const double s = design::softplus<double>(u, tau);
                if (s < std::fmax(u, 0.0) - 1e-15) always_over = false;
            }
            f << tau << "," << at0 << "," << want << ","
              << std::fabs(at0 - want) << ","
              << design::softplus<double>(1.0, tau) << ","
              << design::softplus<double>(-1.0, tau) << ","
              << (always_over ? 1 : 0) << "\n";
        }
        report.check_predicate(
            "o0c.softplus_surrogate_properties",
            "s_tau(u) = tau log(1 + exp(u/tau)) is evaluated without overflow, "
            "equals tau log 2 at zero, and overestimates max(u, 0) everywhere",
            VerificationLevel::A,
            worst_bias_error < 1e-15 && always_over,
            "s_tau(0) = tau log 2 to roundoff, and s_tau >= max(u,0)",
            "worst |s_tau(0) - tau log 2| = " + std::to_string(worst_bias_error) +
                " over tau in {1e-1 ... 1e-5}; overestimates max everywhere: " +
                (always_over ? "yes" : "NO"),
            sw.seconds(),
            "the overestimate is one-sided and largest exactly where rho "
            "crosses zero, which is where certification happens -- so the "
            "surrogate is never substituted for the exact objective, only "
            "differentiated");
        sw.reset();
    }
}

}  // namespace checks
}  // namespace warp_hypad
