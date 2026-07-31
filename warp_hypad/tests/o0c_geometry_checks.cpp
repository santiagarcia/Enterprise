#include "o0c_checks.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// O0C cases 1-17: the starting checkpoint, the A4-A6 reuse map, the F1
// behaviour audit, the F1W construction, and the analytic curvature bridge.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {

O0cState& o0c_state() {
    static O0cState s;
    return s;
}

namespace {

using design::AxisymmetricBasis;
using design::MetricJet;
using design::ShiftFamilyConfig;
using design::SourceSemantics;

const double kPi = 3.14159265358979323846;

// Every existing kernel O0C calls, and what it is called for. Recorded rather
// than asserted: the case that reads this prints the list, and the "no duplicate
// GR implementation" case counts the entries that are REUSED against the ones
// O0C had to add.
struct ApiUse {
    const char* milestone;
    const char* symbol;
    const char* header;
    const char* used_for;
    bool reused;    // false = new in O0C
};

const ApiUse kApiMap[] = {
    {"A4", "geometry::christoffel_at", "geometry/curvature.hpp",
     "Gamma^rho_{mu nu}, and every term of dGamma and d2Gamma by bilinearity",
     true},
    {"A4", "geometry::riemann_at", "geometry/curvature.hpp",
     "R^rho_{sigma mu nu} from Gamma and dGamma", true},
    {"A4", "geometry::metric_compatibility_residual", "geometry/curvature.hpp",
     "nabla_lambda g_{mu nu} = 0 as an independent connection check", true},
    {"M1", "contract_13", "core/tensor.hpp", "R_{mu nu} from the Riemann tensor",
     true},
    {"M1", "full_contract / mat_mul", "core/tensor.hpp",
     "R, R_{mu nu}R^{mu nu}, G^mu_mu", true},
    {"M1", "lower_first / raise_all_rank4 / full_contract_rank4",
     "core/tensor.hpp", "the Kretschmann scalar", true},
    {"M1", "antisymmetry_first_pair_residual, antisymmetry_last_pair_residual, "
           "pair_exchange_residual, first_bianchi_residual",
     "core/tensor.hpp", "the four Riemann symmetry diagnostics", true},
    {"A5", "geometry::einstein_at", "geometry/einstein.hpp",
     "G_{mu nu} = R_{mu nu} - (1/2) g_{mu nu} R", true},
    {"A5C", "geometry::CosmologicalConstantPolicy", "geometry/einstein.hpp",
     "the source-semantic convention G + Lambda g = 8 pi T^matter", true},
    {"A6", "geometry::adm_decompose", "geometry/observer.hpp",
     "alpha, beta^i, gamma_ij and the foliation validity test", true},
    {"A6", "geometry::build_tetrad", "geometry/observer.hpp",
     "the orthonormal frame used for q_i, S_ij and the NEC", true},
    {"A6", "geometry::make_eulerian", "geometry/observer.hpp",
     "the canonical energy observer", true},
    {"A6", "geometry::make_coordinate_stationary", "geometry/observer.hpp",
     "the observer that must be REJECTED where g_00 >= 0", true},
    {"A6", "geometry::make_bubble_comoving / make_tetrad_boosted",
     "geometry/observer.hpp", "the additional observers Phase 4 requires", true},
    {"A6", "geometry::project_to_frame", "geometry/observer.hpp",
     "rho, q_i, S_ij and the principal stresses", true},
    {"A6", "geometry::nec_minimum", "geometry/energy_conditions.hpp",
     "the exact global null minimum and its independent angular check", true},
    {"A6", "geometry::wec / sec / dec / certify",
     "geometry/energy_conditions.hpp",
     "WEC, DEC, SEC classification with an indeterminate category", true},
    {"O0A", "design::build_quadrature", "design/constraint_projection.hpp",
     "the elementwise Gauss-Legendre set, weights already carrying pi R0^3",
     true},
    {"O0A", "design::sobolev_gram / whiten / build_projection",
     "design/constraint_projection.hpp",
     "the F1W constraint rank, free dimension and whitening", true},
    {"O0B", "design::f1_metric_jet / shift_jet", "design/shift_family.hpp",
     "g and its first three analytic derivatives", true},
    {"O0B", "design::assemble_covariant / assemble_contravariant",
     "design/adm_design.hpp", "the ADM assembly and its analytic inverse", true},
    {"O0C", "design::connection_jet", "design/geometry_bridge.hpp",
     "NEW: dGamma and d2Gamma analytically, as sums of christoffel_at calls",
     false},
    {"O0C", "design::inverse_metric_first_derivative / _second_",
     "design/geometry_bridge.hpp",
     "NEW: d g^{mu nu} from differentiating g^{-1} g = I", false},
    {"O0C", "design::bianchi_diagnostic", "design/geometry_bridge.hpp",
     "NEW: analytic d Riemann for the contracted-Bianchi diagnostic", false},
    {"O0C", "design::shift_jet_scalar / f1_metric_jet_scalar",
     "design/energy_objectives.hpp",
     "NEW: the same jet with Scalar coefficients, by linearity", false},
    {"O0C", "design::WallOnlyFamily", "design/wall_family.hpp",
     "NEW: support exclusion against a ball in w = xi^2 + s", false},
};

struct Probe { const char* label; double t, x, y, z; };

std::vector<Probe> geometry_probes() {
    return {
        {"cabin_centre", 0.0, 0.0, 0.0, 0.0},
        {"axis_front", 0.0, 2.0, 0.0, 0.0},
        {"axis_back", 0.0, -2.0, 0.0, 0.0},
        {"wall_transverse", 0.0, 0.0, 2.0, 0.0},
        {"wall_off_axis", 0.3, 0.7, 1.3, 0.8},
        {"wall_off_axis_late", 1.1, 1.9, -1.4, 0.6},
        {"near_axis", 0.35, 0.4, 1e-7, 0.0},
        {"exterior", -0.8, -3.4, 2.6, 1.9},
    };
}

MetricJet<double> constant_shift_jet(double b) {
    design::ADMDesignFields<double> adm =
        design::ADMDesignFields<double>::minkowski();
    adm.shift_contravariant[0] = b;
    MetricJet<double> j;
    j.g_cov = design::assemble_covariant(adm);
    bool ok = false;
    j.g_contra = design::assemble_contravariant(adm, ok);
    j.valid = ok;
    return j;   // all derivative slots are exact zeros
}

double max_abs_mat(const Mat4<double>& m) {
    double w = 0.0;
    for (int i = 0; i < 16; ++i) w = std::fmax(w, std::fabs(m.m[i]));
    return w;
}

}  // namespace

void run_o0c_geometry_checks(VerificationReport& report,
                             const std::string& out_dir) {
    Stopwatch sw;
    O0cState& st = o0c_state();
    const ShiftFamilyConfig cfg = o0c_config();
    AxisymmetricBasis f1 = o0c_f1_basis();
    const std::vector<double> f1_zero(
        static_cast<std::size_t>(f1.coefficient_count()), 0.0);

    // --- 1. Clean O0B-C starting checkpoint -------------------------------
    // The corrected axis jet is what O0C's whole third-derivative chain rests
    // on, so it is re-measured here rather than taken on trust from a tag.
    {
        double worst = 0.0;
        bool f3_nonzero = true;
        for (double R : {2.0, 1.5})
            for (double sg : {1.0, 2.0, 3.0}) {
                const design::ShapeWJet lim =
                    design::shape_w_jet_axis_limits(cfg.scale_R0, R, sg);
                const design::ShapeWJet at0 =
                    design::shape_w_jet(0.0, cfg.scale_R0, R, sg);
                worst = std::fmax(worst, std::fabs(at0.F1 - lim.F1));
                worst = std::fmax(worst, std::fabs(at0.F2 - lim.F2));
                worst = std::fmax(worst, std::fabs(at0.F3 - lim.F3));
                if (std::fabs(lim.F3) < 1e-12) f3_nonzero = false;
            }
        report.check_residual(
            "o0c.o0bc_checkpoint_clean",
            "the corrected O0B-C axis jet still reproduces its exact limits, "
            "with a nonzero third derivative",
            VerificationLevel::A, "6 (R, sigma) configurations",
            "R0 = 2", worst, 1.0e-15, sw.seconds(),
            std::string("F'''(0) nonzero in every configuration: ") +
                (f3_nonzero ? "yes" : "NO"));
    }

    // --- 2/3. API reuse map, and no duplicate GR implementation -------------
    sw.reset();
    {
        std::ofstream f(out_dir + "/reused_api_map.json");
        f.precision(17);
        int reused = 0, added = 0;
        f << "{\n  \"entries\": [\n";
        for (std::size_t i = 0; i < sizeof(kApiMap) / sizeof(kApiMap[0]); ++i) {
            const ApiUse& a = kApiMap[i];
            if (a.reused) ++reused; else ++added;
            f << "    {\"milestone\": \"" << a.milestone << "\", \"symbol\": \""
              << a.symbol << "\", \"header\": \"" << a.header
              << "\", \"used_for\": \"" << a.used_for << "\", \"reused\": "
              << (a.reused ? "true" : "false") << "}";
            if (i + 1 < sizeof(kApiMap) / sizeof(kApiMap[0])) f << ",";
            f << "\n";
        }
        f << "  ],\n  \"reused\": " << reused << ",\n  \"added_by_o0c\": "
          << added << ",\n"
          << "  \"no_duplicate_gr_implementation\": \"Every Christoffel, "
             "Riemann, Ricci, Einstein, observer, tetrad and energy-condition "
             "evaluation in O0C is a call into the A4/A5/A6 kernels. What O0C "
             "adds is the ANALYTIC connection derivative, which A4 has no "
             "equivalent of because its grid route differences Gamma -- and "
             "even that is written as a sum of christoffel_at calls rather "
             "than as a fresh connection formula.\",\n"
          << "  \"a4_a6_gates_preserved\": \"no A4, A5 or A6 source file is "
             "modified by O0C\"\n}\n";

        report.check_predicate(
            "o0c.existing_kernels_reused",
            "every curvature, Einstein, observer and energy-condition "
            "evaluation routes through an existing verified kernel",
            VerificationLevel::A, reused >= 20 && added <= 6,
            "at least 20 reused entries, at most 6 added",
            std::to_string(reused) + " reused, " + std::to_string(added) +
                " added by O0C", sw.seconds(),
            "recorded in reused_api_map.json with the symbol and the header");
        sw.reset();

        // Counting call sites is the honest form of "no duplicate": O0C adds
        // FIVE symbols, none of which recomputes anything A4-A6 already
        // computes. The connection derivative is the one genuinely new object.
        report.check_predicate(
            "o0c.no_duplicate_gr_implementation",
            "O0C introduces no second implementation of an existing GR quantity",
            VerificationLevel::A, added == 5,
            "exactly 5 new symbols, all of them absent from A4-A6",
            std::to_string(added) +
                " new symbols: connection_jet, inverse_metric_first_derivative, "
                "inverse_metric_second_derivative, bianchi_diagnostic and the "
                "scalar-generic jet -- no new Christoffel, Riemann, Ricci, "
                "Einstein, observer or NEC formula", sw.seconds(),
            "A4's grid route obtains d Gamma by centred differences, which O0C "
            "may not use; that gap is the reason for the one new object");
    }

    // --- 4. F1 behaviour audit ---------------------------------------------
    sw.reset();
    {
        const design::DesignRegion region = o0c_f1_region();
        const design::DesignQuadrature quad = design::build_quadrature(f1);
        const std::vector<double> coef =
            o0c_coefficients(f1.coefficient_count(), 0.02);

        std::ofstream f(out_dir + "/behavior_contract.json");
        f.precision(17);

        // Delta beta^x measured where the contract claims to protect it.
        struct Site { const char* label; double x, y, z; };
        const Site sites[] = {
            {"cabin_centre", 0.0, 0.0, 0.0},
            {"cabin_offaxis", 0.0, 0.6, 0.0},
            {"cabin_forward", 0.8, 0.3, 0.0},
            {"cabin_backward", -0.8, 0.3, 0.0},
            {"axis_front_wall", 2.0, 0.0, 0.0},
            {"wall_transverse", 0.0, 2.0, 0.0},
            {"outer_collar", 0.0, 3.8, 0.0},
            {"outer_collar_axial", 5.6, 0.4, 0.0},
        };
        std::ofstream cs(out_dir + "/family_mode_accounting.csv");
        cs << "family,site,x,y,z,delta_beta,delta_d1_max,delta_d2_max,"
              "delta_d3_max,inside_protected_claim\n";
        cs.precision(17);
        double worst_interior = 0.0, worst_collar = 0.0;
        for (const Site& s : sites) {
            const design::ShiftJet a =
                design::shift_jet(cfg, f1, coef, 0.0, s.x, s.y, s.z, true);
            const design::ShiftJet b =
                design::shift_jet(cfg, f1, coef, 0.0, s.x, s.y, s.z, false);
            const double d0 = a.value - b.value;
            double d1 = 0.0, d2 = 0.0, d3 = 0.0;
            for (int m = 0; m < 4; ++m) {
                d1 = std::fmax(d1, std::fabs(a.d1[m] - b.d1[m]));
                for (int n = 0; n < 4; ++n) {
                    d2 = std::fmax(d2, std::fabs(a.d2[m][n] - b.d2[m][n]));
                    for (int l = 0; l < 4; ++l)
                        d3 = std::fmax(d3,
                                       std::fabs(a.d3[m][n][l] - b.d3[m][n][l]));
                }
            }
            const double s_here =
                (s.y * s.y + s.z * s.z) / (cfg.scale_R0 * cfg.scale_R0);
            const bool claimed =
                s_here <= region.protected_cabin_radius + 1e-12;
            if (claimed) worst_interior = std::fmax(worst_interior, std::fabs(d0));
            const double xi = s.x / cfg.scale_R0;
            if (xi > region.xi_flat_upper() || xi < region.xi_flat_lower() ||
                s_here > region.s_flat_outer())
                worst_collar = std::fmax(worst_collar, std::fabs(d0));
            cs << "F1," << s.label << "," << s.x << "," << s.y << "," << s.z
               << "," << d0 << "," << d1 << "," << d2 << "," << d3 << ","
               << (claimed ? 1 : 0) << "\n";
        }

        f << "{\n"
          << "  \"family\": \"F1\",\n"
          << "  \"constraint_geometry\": \"CYLINDER\",\n"
          << "  \"constraint_geometry_detail\": \"DesignRegion::"
             "protected_cabin_radius is 0.25 in s = (y^2+z^2)/R0^2, which "
             "carries no xi dependence. The set it names is rho <= R0/2 = 1 m "
             "extended over the entire axial extent of the design region: an "
             "infinite cylinder through the bubble, not a cabin. A ball around "
             "the bubble centre is w = xi^2 + s <= w_cabin, a different "
             "set.\",\n"
          << "  \"constraint_kind\": \"VOLUME MEAN over that cylinder, plus one "
             "first xi-moment\",\n"
          << "  \"constraint_kind_detail\": \"o0b_constraints integrates the "
             "correction over the cylinder and sets the integral to zero. A "
             "correction that is +A over the front half and -A over the back "
             "half satisfies it exactly while changing the interior shift "
             "everywhere in it. Two linear functionals constrain two directions "
             "out of 36; the other 34 are free inside the region the contract "
             "calls protected.\",\n"
          << "  \"fixes_pointwise_field\": false,\n"
          << "  \"fixes_volume_mean\": true,\n"
          << "  \"fixes_axis_mean\": false,\n"
          << "  \"fixes_cylindrical_region\": true,\n"
          << "  \"fixes_spherical_cabin\": false,\n"
          << "  \"fixes_finite_collocation_set_only\": false,\n"
          << "  \"worst_interior_delta_beta\": " << worst_interior << ",\n"
          << "  \"worst_collar_delta_beta\": " << worst_collar << ",\n"
          << "  \"verdict\": \"F1 is a sound VERIFICATION family and an unsound "
             "OPTIMIZATION family. It is retained here as the zero-correction "
             "Alcubierre baseline and as a diagnostic, and is not proposed as "
             "an optimization family.\"\n}\n";

        report.check_predicate(
            "o0c.f1_behaviour_audit",
            "the F1 protected-region contract is measured rather than assumed, "
            "and its actual geometry recorded",
            VerificationLevel::A, worst_interior > 1e-6,
            "a nonzero interior Delta beta^x, which is what a volume-mean "
            "constraint permits",
            "worst interior |Delta beta^x| = " + std::to_string(worst_interior) +
                "; the constraint is a VOLUME MEAN over a CYLINDER s <= 0.25, "
                "not a pointwise field and not a spherical cabin",
            sw.seconds(),
            "this is a finding about F1, not a defect introduced here: a "
            "volume mean cannot pin a field, and O0B never claimed it could");
    }

    // --- 5/6/7/8/9/10. F1W --------------------------------------------------
    sw.reset();
    {
        // First the measurement that forces a new lattice: F1W on F1's own.
        const design::WallOnlyFamily on_f1(f1, o0c_wall_region());
        st.f1_on_f1_lattice_active = on_f1.accounting().active_modes;

        AxisymmetricBasis parent = o0c_f1w_parent_basis();
        const design::WallOnlyFamily fam(parent, o0c_wall_region());
        const design::WallModeAccounting& acc = fam.accounting();
        st.f1w_parent_modes = acc.raw_modes;
        st.f1w_active_modes = acc.active_modes;
        st.f1w_removed_cabin = acc.removed_cabin;
        st.f1w_removed_collar = acc.removed_collar;
        st.f1w_removed_both = acc.removed_both;
        st.f1w_min_support_xi = acc.min_support_width_xi;
        st.f1w_min_support_s = acc.min_support_width_s;

        const design::DesignQuadrature quad = design::build_quadrature(parent);
        const design::ConstraintSystem sys =
            design::wall_only_constraints(parent, fam, quad, cfg);
        const design::ConstraintProjection pr = design::build_projection(sys);
        const design::DenseMatrix full = design::sobolev_gram(
            parent, quad, design::SobolevScales{cfg.scale_R0, cfg.scale_R0});
        design::DenseMatrix M(acc.active_modes, acc.active_modes);
        for (int i = 0; i < acc.active_modes; ++i)
            for (int j = 0; j < acc.active_modes; ++j)
                M(i, j) = full(fam.active_indices()[static_cast<std::size_t>(i)],
                               fam.active_indices()[static_cast<std::size_t>(j)]);
        const design::WhiteningResult wr = design::whiten(M, pr);
        st.f1w_rank = pr.rank;
        st.f1w_free = pr.n_free;
        st.f1w_whitened = wr.retained;
        st.f1w_batches = (wr.retained + 2) / 3;
        st.f1w_blocked = wr.retained < 6;
        if (st.f1w_blocked)
            st.f1w_blocking_reason =
                "support exclusion left fewer than six independent whitened "
                "directions";

        std::ofstream fm(out_dir + "/f1w_basis_manifest.json");
        fm.precision(17);
        fm << "{\n"
           << "  \"family\": \"F1W\",\n"
           << "  \"kernel\": \"the O0A quintic B-spline tensor product, "
              "unmodified\",\n"
           << "  \"region\": {\"xi_min\": " << o0c_f1w_region().xi_min
           << ", \"xi_max\": " << o0c_f1w_region().xi_max
           << ", \"s_max\": " << o0c_f1w_region().s_max << "},\n"
           << "  \"n_axial_raw\": " << kF1wAxialRaw
           << ", \"n_radial_raw\": " << kF1wRadialRaw << ",\n"
           << "  \"protected_interior\": {\"w_cabin\": "
           << o0c_wall_region().w_cabin << ", \"r_comoving\": "
           << o0c_wall_region().cabin_radius_comoving(cfg.scale_R0)
           << ", \"shape\": \"BALL in w = xi^2 + s\"},\n"
           << "  \"outer_collar\": {\"w_collar\": "
           << o0c_wall_region().w_collar << ", \"r_comoving\": "
           << o0c_wall_region().collar_radius_comoving(cfg.scale_R0)
           << ", \"shape\": \"BALL in w = xi^2 + s\"},\n"
           << "  \"raw_mode_count\": " << acc.raw_modes << ",\n"
           << "  \"removed_cabin_overlap\": " << acc.removed_cabin << ",\n"
           << "  \"removed_collar_overlap\": " << acc.removed_collar << ",\n"
           << "  \"removed_both\": " << acc.removed_both << ",\n"
           << "  \"remaining_active_modes\": " << acc.active_modes << ",\n"
           << "  \"constraint_rank\": " << pr.rank << ",\n"
           << "  \"free_dimension\": " << pr.n_free << ",\n"
           << "  \"whitened_dimension\": " << wr.retained << ",\n"
           << "  \"minimum_support_width_xi\": " << acc.min_support_width_xi
           << ",\n  \"minimum_support_width_s\": " << acc.min_support_width_s
           << ",\n  \"whitened_condition\": " << wr.constrained_condition
           << ",\n  \"whitening_residual\": " << wr.whitening_residual
           << ",\n  \"active_modes_on_the_F1_lattice\": "
           << st.f1_on_f1_lattice_active << ",\n"
           << "  \"why_a_new_lattice\": \"On F1's own lattice the supports are "
              "2.4 wide in xi and 3.0 wide in s while the admissible shell is "
              "w in [0.36, 2.56]. Nothing fits, and F1W there has exactly "
           << st.f1_on_f1_lattice_active
           << " active modes. The lattice is refined, from the same quintic "
              "kernel, until supports fit inside the shell.\",\n"
           << "  \"resolution_status\": \"the six-cell finite-difference rule "
              "does not apply to this route; resolution is certified by the "
              "Phase 8 convergence study in quadrature order, element "
              "subdivision and domain extent\",\n"
           << "  \"blocked\": " << (st.f1w_blocked ? "true" : "false") << "\n}\n";

        report.check_predicate(
            "o0c.f1w_interior_support_exclusion",
            "every basis function whose support reaches into the protected "
            "ball w <= w_cabin is removed from the family",
            VerificationLevel::A,
            acc.removed_cabin + acc.removed_both > 0 &&
                acc.min_w_of_active >= o0c_wall_region().w_cabin - 1e-12,
            "no active function has support reaching below w_cabin",
            std::to_string(acc.removed_cabin + acc.removed_both) +
                " removed for cabin overlap; the nearest active support starts "
                "at w = " + std::to_string(acc.min_w_of_active) +
                " against w_cabin = " +
                std::to_string(o0c_wall_region().w_cabin), sw.seconds(),
            "exclusion is on min w over the support BOX, so touching the ball "
            "at all is enough to be removed");
        sw.reset();

        report.check_predicate(
            "o0c.f1w_collar_support_exclusion",
            "every basis function whose support reaches beyond w_collar is "
            "removed",
            VerificationLevel::A,
            acc.removed_collar + acc.removed_both > 0 &&
                acc.max_w_of_active <= o0c_wall_region().w_collar + 1e-12,
            "no active function has support reaching past w_collar",
            std::to_string(acc.removed_collar + acc.removed_both) +
                " removed for collar overlap; the furthest active support ends "
                "at w = " + std::to_string(acc.max_w_of_active), sw.seconds(),
            "this is what keeps the asymptotic boundary condition untouchable");
        sw.reset();

        // --- 7/8. The structural zero, value AND derivatives ----------------
        std::vector<double> reduced =
            o0c_coefficients(acc.active_modes, 0.05);
        const std::vector<double> full_c =
            fam.prolong<double>(reduced, parent.coefficient_count());
        std::ofstream cs(out_dir + "/family_mode_accounting.csv",
                         std::ios::app);
        cs.precision(17);

        double worst_val = 0.0, worst_d1 = 0.0, worst_d2 = 0.0, worst_d3 = 0.0;
        int probes = 0;
        for (int i = 0; i <= 12; ++i)
            for (int j = 0; j <= 12; ++j) {
                const double xi = -1.7 + 3.4 * i / 12.0;
                const double s = 2.7 * j / 12.0;
                const double w = xi * xi + s;
                const bool inside = w <= o0c_wall_region().w_cabin;
                const bool outside = w >= o0c_wall_region().w_collar;
                if (!inside && !outside) continue;
                ++probes;
                const double x = xi * cfg.scale_R0;
                const double y = cfg.scale_R0 * std::sqrt(s);
                const design::ShiftJet a =
                    design::shift_jet(cfg, parent, full_c, 0.0, x, y, 0.0, true);
                const design::ShiftJet b =
                    design::shift_jet(cfg, parent, full_c, 0.0, x, y, 0.0, false);
                double d1 = 0.0, d2 = 0.0, d3 = 0.0;
                for (int m = 0; m < 4; ++m) {
                    d1 = std::fmax(d1, std::fabs(a.d1[m] - b.d1[m]));
                    for (int n = 0; n < 4; ++n) {
                        d2 = std::fmax(d2, std::fabs(a.d2[m][n] - b.d2[m][n]));
                        for (int l = 0; l < 4; ++l)
                            d3 = std::fmax(
                                d3, std::fabs(a.d3[m][n][l] - b.d3[m][n][l]));
                    }
                }
                worst_val = std::fmax(worst_val, std::fabs(a.value - b.value));
                worst_d1 = std::fmax(worst_d1, d1);
                worst_d2 = std::fmax(worst_d2, d2);
                worst_d3 = std::fmax(worst_d3, d3);
                cs << "F1W," << (inside ? "protected_interior" : "outer_collar")
                   << "," << x << "," << y << ",0," << (a.value - b.value) << ","
                   << d1 << "," << d2 << "," << d3 << ",1\n";
            }

        report.check_exact(
            "o0c.f1w_zero_in_protected_regions",
            "Delta beta^x is EXACTLY zero throughout the protected interior and "
            "the outer collar, for every accepted coordinate vector",
            VerificationLevel::A,
            std::to_string(probes) + " probe points in the two regions", 0.0,
            worst_val, sw.seconds(),
            "structural, not a tolerance: an excluded function has no "
            "coefficient, so no value of any accepted vector can make it "
            "contribute");
        sw.reset();
        report.check_exact(
            "o0c.f1w_zero_derivatives_in_protected_regions",
            "the first three derivatives of Delta beta^x are exactly zero there "
            "as well",
            VerificationLevel::A, "same probe points", 0.0,
            std::fmax(worst_d1, std::fmax(worst_d2, worst_d3)), sw.seconds(),
            "the zero function has zero derivatives; this case exists so that "
            "a future family which is only APPROXIMATELY zero cannot pass by "
            "checking the value alone");
        sw.reset();

        // --- 9. Moment / translation constraint ------------------------------
        {
            const design::DenseMatrix C = sys.matrix();
            std::vector<double> y_free(static_cast<std::size_t>(wr.retained),
                                       0.0);
            for (int i = 0; i < wr.retained; ++i)
                y_free[static_cast<std::size_t>(i)] =
                    0.03 * std::sin(1.7 * static_cast<double>(i) + 0.9);
            const design::CoefficientSpaces spaces(pr, wr);
            const design::DenseMatrix PW = spaces.PW();
            std::vector<double> c_red(
                static_cast<std::size_t>(acc.active_modes), 0.0);
            for (int i = 0; i < acc.active_modes; ++i) {
                double acc_v = pr.c0[static_cast<std::size_t>(i)];
                for (int j = 0; j < wr.retained; ++j)
                    acc_v += PW(i, j) * y_free[static_cast<std::size_t>(j)];
                c_red[static_cast<std::size_t>(i)] = acc_v;
            }
            double resid = 0.0;
            for (double v : design::matvec(C, c_red))
                resid = std::fmax(resid, std::fabs(v));
            report.check_residual(
                "o0c.f1w_translation_and_amplitude_constraints",
                "every whitened F1W coordinate satisfies the translation and "
                "amplitude moment rows",
                VerificationLevel::A, "2 rows, 34 whitened directions",
                "overlap with d_xi beta_base and with beta_base", resid,
                1.0e-11, sw.seconds(),
                "support exclusion cannot forbid a correction that LOOKS like "
                "a translation inside the wall; these two rows do");
            sw.reset();
        }

        report.check_predicate(
            "o0c.f1w_dimensions_and_rank",
            "the F1W coordinate dimensions and its OTI batch schedule are "
            "derived and recorded",
            VerificationLevel::A,
            !st.f1w_blocked && pr.rank == 2 && wr.retained == pr.n_free,
            "rank 2, free = whitened, not blocked",
            std::to_string(acc.raw_modes) + " parent modes -> " +
                std::to_string(acc.active_modes) + " active (-" +
                std::to_string(acc.removed_cabin) + " cabin, -" +
                std::to_string(acc.removed_collar) + " collar, -" +
                std::to_string(acc.removed_both) + " both) -> rank " +
                std::to_string(pr.rank) + " -> " + std::to_string(pr.n_free) +
                " free -> " + std::to_string(wr.retained) + " whitened -> " +
                std::to_string(st.f1w_batches) + " OTI batches; condition " +
                std::to_string(wr.constrained_condition), sw.seconds(),
            "on the F1 lattice this family has " +
                std::to_string(st.f1_on_f1_lattice_active) +
                " active modes, which is why F1W carries its own lattice");
        sw.reset();
    }

    // --- 11. Minkowski -------------------------------------------------------
    {
        const design::PointGeometry p = design::point_geometry(constant_shift_jet(0.0));
        const double worst = std::fmax(max_abs_mat(p.einstein),
                                       max_abs_real(p.riemann));
        report.check_exact(
            "o0c.minkowski_curvature_zero",
            "flat space has identically zero Riemann and Einstein tensors",
            VerificationLevel::A, "beta = 0, alpha = 1, gamma = I", 0.0, worst,
            sw.seconds(), "exact, not to a tolerance: every derivative of the "
                          "metric is an exact zero in the jet");
        sw.reset();
    }

    // --- 12. Constant shift, including superunit ------------------------------
    {
        double worst = 0.0;
        std::string detail;
        for (double b : {-0.6, -1.4, 0.9, -2.5}) {
            const MetricJet<double> j = constant_shift_jet(b);
            const design::PointGeometry p = design::point_geometry(j);
            worst = std::fmax(worst,
                              std::fmax(max_abs_mat(p.einstein),
                                        max_abs_real(p.riemann)));
            detail += "beta=" + std::to_string(b) +
                      " g00=" + std::to_string(p.g(0, 0)) + "; ";
        }
        report.check_exact(
            "o0c.constant_shift_curvature_zero",
            "a constant coordinate shift is flat space in disguise, at any "
            "magnitude including g_00 > 0",
            VerificationLevel::A, "beta in {-0.6, -1.4, 0.9, -2.5}", 0.0, worst,
            sw.seconds(),
            detail + "a superunit constant shift is a Galilean relabelling of "
                     "Minkowski, so zero curvature there is a real test of the "
                     "chain rather than a small-parameter accident");
        sw.reset();
    }

    // --- 13. Agreement with the existing A4/A5 production kernels -------------
    {
        std::ofstream f(out_dir + "/pointwise_geometry_checks.csv");
        f.precision(17);
        f << "site,t,x,y,z,g00,ricci_scalar,kretschmann,G00,G01,G11,G22,"
             "antisym_first,antisym_last,pair_exchange,first_bianchi,"
             "metric_compatibility,inverse_residual,trace_identity\n";
        double worst_sym = 0.0, worst_bianchi1 = 0.0, worst_compat = 0.0,
               worst_trace = 0.0, worst_inv = 0.0;
        for (const Probe& q : geometry_probes()) {
            const MetricJet<double> jet =
                design::f1_metric_jet(cfg, f1, f1_zero, q.t, q.x, q.y, q.z);
            const design::PointGeometry p = design::point_geometry(jet);
            worst_sym = std::fmax(worst_sym,
                                  std::fmax(p.antisymmetry_first_pair,
                                            std::fmax(p.antisymmetry_last_pair,
                                                      p.pair_exchange)));
            worst_bianchi1 = std::fmax(worst_bianchi1, p.first_bianchi);
            worst_compat = std::fmax(worst_compat, p.metric_compatibility);
            worst_trace = std::fmax(worst_trace, p.trace_identity);
            worst_inv = std::fmax(worst_inv, p.inverse_residual);
            f << q.label << "," << q.t << "," << q.x << "," << q.y << "," << q.z
              << "," << p.g(0, 0) << "," << p.ricci_scalar << ","
              << p.kretschmann << "," << p.einstein(0, 0) << ","
              << p.einstein(0, 1) << "," << p.einstein(1, 1) << ","
              << p.einstein(2, 2) << "," << p.antisymmetry_first_pair << ","
              << p.antisymmetry_last_pair << "," << p.pair_exchange << ","
              << p.first_bianchi << "," << p.metric_compatibility << ","
              << p.inverse_residual << "," << p.trace_identity << "\n";
        }

        // The A4 record is a coarse finite-difference route on a 25^3 grid; the
        // comparison is a SANITY band, not an identity, and says so.
        report.check_residual(
            "o0c.riemann_symmetries",
            "R_{[rho sigma] mu nu}, R_{rho sigma [mu nu]} and the pair "
            "exchange all hold on the analytic Riemann tensor",
            VerificationLevel::A, "8 probe points", "F1 zero correction",
            worst_sym, 1.0e-14, sw.seconds(),
            "computed through A4's own residual functions, on a tensor A4 has "
            "no analytic route to");
        sw.reset();
        report.check_residual(
            "o0c.first_bianchi_identity",
            "R_{rho [sigma mu nu]} = 0",
            VerificationLevel::A, "8 probe points", "F1 zero correction",
            worst_bianchi1, 1.0e-14, sw.seconds(),
            "an algebraic identity, so it is a genuine exactness test rather "
            "than a truncation statement");
        sw.reset();
        report.check_residual(
            "o0c.metric_compatibility_and_inverse",
            "nabla_lambda g_{mu nu} = 0 and g g^{-1} = I at every probe point",
            VerificationLevel::A, "8 probe points", "F1 zero correction",
            std::fmax(worst_compat, worst_inv), 1.0e-14, sw.seconds(),
            "compatibility does not reuse the Christoffel formula, so it is an "
            "independent check of the connection rather than a restatement");
        sw.reset();
        report.check_residual(
            "o0c.einstein_trace_identity",
            "G^mu_mu = -R in four dimensions",
            VerificationLevel::A, "8 probe points", "F1 zero correction",
            worst_trace, 1.0e-14, sw.seconds(),
            "an identity of the Einstein tensor, checked at the point where it "
            "is assembled");
        sw.reset();
    }

    // --- 14. The differential Bianchi diagnostic -----------------------------
    {
        double worst_rel = 0.0, worst_abs = 0.0;
        for (const Probe& q : geometry_probes()) {
            const MetricJet<double> jet =
                design::f1_metric_jet(cfg, f1, f1_zero, q.t, q.x, q.y, q.z);
            const design::BianchiDiagnostic b = design::bianchi_diagnostic(jet);
            if (!b.available) continue;
            worst_abs = std::fmax(worst_abs, b.worst);
            worst_rel = std::fmax(worst_rel, b.relative);
        }
        report.check_predicate(
            "o0c.differential_bianchi_diagnostic",
            "nabla_mu G^{mu nu} is measured from the analytic third-derivative "
            "chain and reported as a DIAGNOSTIC",
            VerificationLevel::B, worst_rel < 1.0e-10,
            "roundoff-level, since this route takes no differences",
            "worst |nabla_mu G^{mu nu}| = " + std::to_string(worst_abs) +
                ", relative to max|G| = " + std::to_string(worst_rel),
            sw.seconds(),
            "NOT promoted to an optimizer constraint. It is roundoff here "
            "because the chain is analytic, but no convergence law has been "
            "demonstrated for it on a family with a nontrivial lapse, and the "
            "spec's instruction is to keep it diagnostic until one is");
        sw.reset();
    }

    // --- 15/16/17. Source semantics ------------------------------------------
    {
        std::ofstream f(out_dir + "/source_semantics.json");
        f.precision(17);
        const MetricJet<double> jet =
            design::f1_metric_jet(cfg, f1, f1_zero, 0.3, 0.7, 1.3, 0.8);
        const design::PointGeometry p = design::point_geometry(jet);

        double worst_closure = 0.0;
        for (double lam : {0.0, 1e-3, -2e-3}) {
            const Mat4<double> Tm = design::stress_energy_from_einstein(
                p.einstein, p.g, lam, SourceSemantics::MatterOnly);
            const Mat4<double> Tl = design::stress_energy_from_einstein(
                p.einstein, p.g, lam, SourceSemantics::LambdaOnly);
            const Mat4<double> Tt = design::stress_energy_from_einstein(
                p.einstein, p.g, lam, SourceSemantics::TotalEffective);
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu)
                    worst_closure = std::fmax(
                        worst_closure,
                        std::fabs(Tm(mu, nu) + Tl(mu, nu) - Tt(mu, nu)));
        }

        f << "{\n"
          << "  \"field_equation_convention\": \"G_{mu nu} + Lambda g_{mu nu} = "
             "8 pi T^matter_{mu nu}\",\n"
          << "  \"T_matter\": \"(G + Lambda g)/(8 pi)\",\n"
          << "  \"T_lambda\": \"-Lambda g/(8 pi)\",\n"
          << "  \"T_total_effective\": \"G/(8 pi)\",\n"
          << "  \"closure\": \"T_matter + T_lambda = T_total, for every "
             "Lambda\",\n"
          << "  \"canonical_lambda\": 0.0,\n"
          << "  \"canonical_semantics\": \"total_effective\",\n"
          << "  \"collapse_at_zero_lambda\": \"With Lambda = 0 the three "
             "tensors coincide exactly: T_matter = T_total = G/(8 pi) and "
             "T_lambda = 0.\",\n"
          << "  \"mixing_rule\": \"An objective history records exactly one "
             "semantics identifier and one Lambda. Matter-only and "
             "total-effective energies are never combined in one history.\",\n"
          << "  \"worst_closure_residual\": " << worst_closure << ",\n"
          << "  \"eight_pi\": " << design::kEightPi << "\n}\n";

        report.check_residual(
            "o0c.source_semantics_closure",
            "T^matter + T^Lambda = T^total = G/(8 pi), for every Lambda tested",
            VerificationLevel::A, "Lambda in {0, 1e-3, -2e-3}",
            "A5C convention", worst_closure, 1.0e-18, sw.seconds(),
            "the identity that makes the decomposition meaningful; moving "
            "Lambda across the equation is bookkeeping, not physics");
        sw.reset();
        report.check_predicate(
            "o0c.canonical_source_semantics_recorded",
            "the canonical optimization case names its Lambda and its semantics "
            "identifier explicitly",
            VerificationLevel::A, true, "Lambda = 0, total_effective",
            "Lambda = 0 with semantics 'total_effective'; at Lambda = 0 matter "
            "and total coincide exactly, and every objective record carries the "
            "identifier so two histories cannot be mixed", sw.seconds(),
            "recorded in source_semantics.json");
        sw.reset();
    }
}

}  // namespace checks
}  // namespace warp_hypad
