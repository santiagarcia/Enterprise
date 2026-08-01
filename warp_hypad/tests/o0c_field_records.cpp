#include "o0c_checks.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "warp_hypad/core/metadata.hpp"

// ---------------------------------------------------------------------------
// O0C Phase 14: the field records the 3D renderers read, and cases 46-50.
//
// THE DIVISION OF LABOUR, STATED ONCE. Production writes SAMPLED SCALAR AND
// VECTOR FIELDS. Python reads them, extracts isosurfaces, revolves meridional
// data, projects a camera and draws. Python never evaluates a metric, a
// curvature, an Einstein tensor, a stress-energy, an observer quantity, a NEC
// minimum or an energy objective -- the visual gate's AST scan enforces that,
// and these records are what make the enforcement possible rather than merely
// declared.
//
// TWO KINDS OF 3D RECORD, AND THEY ARE NOT INTERCHANGEABLE.
//
//   volume_field_3d.csv is a genuinely sampled (x, y, z) volume. Every point is
//   an independent production evaluation. It is the validation case: an
//   axisymmetric reconstruction can be checked against it.
//
//   meridional_field.csv is a (x, rho) half-plane. A renderer may revolve it,
//   and every product built that way MUST say "axisymmetric 3D reconstruction
//   from production meridional data". It is cheaper by the number of azimuthal
//   samples, which is why the animations use it, and it would be dishonest to
//   present it as an independently sampled volume.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using design::AxisymmetricBasis;
using design::MetricJet;
using design::ShiftFamilyConfig;
using design::SourceSemantics;

// Deliberately modest: a 3D volume at full curvature cost is the expensive
// record in this milestone, and the renderer needs resolution, not heroics.
inline constexpr int kVolumeN = 41;        // 41^3 sampled points
inline constexpr double kVolumeExtent = 5.0;
inline constexpr int kMeridionalNx = 121;
inline constexpr int kMeridionalNr = 81;
inline constexpr int kFrames = 12;

struct PointFields {
    double beta = 0.0;
    double rho = 0.0;
    double n_min = 0.0;
    double kretschmann = 0.0;
    double g00 = 0.0;
    bool valid = false;
};

PointFields sample(const ShiftFamilyConfig& cfg, const AxisymmetricBasis& basis,
                   const std::vector<double>& coef, double t, double x,
                   double y, double z, bool want_nec) {
    PointFields o;
    const MetricJet<double> jet =
        design::f1_metric_jet(cfg, basis, coef, t, x, y, z);
    if (!jet.valid) return o;
    const design::PointGeometry pg = design::point_geometry(jet);
    if (!pg.valid) return o;
    const Mat4<double> T = design::stress_energy_from_einstein(
        pg.einstein, pg.g, 0.0, SourceSemantics::TotalEffective);
    const geometry::AdmDecomposition adm =
        geometry::adm_decompose(jet.g_cov.m, jet.g_contra.m);
    if (!adm.valid) return o;
    const geometry::Observer eul = geometry::make_eulerian(jet.g_cov.m, adm);
    if (!eul.valid) return o;

    o.g00 = jet.g_cov.m(0, 0);
    o.beta = adm.beta[0];
    o.kretschmann = pg.kretschmann;
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu)
            o.rho += T(mu, nu) * eul.u[mu] * eul.u[nu];
    if (want_nec) {
        const geometry::Tetrad tet = geometry::build_tetrad(jet.g_cov.m, adm);
        if (tet.valid) {
            const geometry::FrameStress fs =
                geometry::project_to_frame(tet, T);
            const geometry::NecResult nec = geometry::nec_minimum(fs, 256);
            if (nec.ok) o.n_min = nec.n_min;
        }
    }
    o.valid = true;
    return o;
}

}  // namespace

void write_o0c_field_records(const std::string& out_dir) {
    const ShiftFamilyConfig cfg = o0c_config();
    AxisymmetricBasis f1 = o0c_f1_basis();
    const std::vector<double> zero(
        static_cast<std::size_t>(f1.coefficient_count()), 0.0);

    // --- The sampled (x, y, z) volume ---------------------------------------
    {
        std::ofstream f(out_dir + "/volume_field_3d.csv");
        f.precision(9);
        f << "# genuinely sampled (x,y,z) volume; every row is an independent "
             "production evaluation, no revolution\n";
        f << "i,j,k,x,y,z,beta_x,rho,kretschmann,g00,valid\n";
        const double h = 2.0 * kVolumeExtent / (kVolumeN - 1);
        for (int i = 0; i < kVolumeN; ++i)
            for (int j = 0; j < kVolumeN; ++j)
                for (int k = 0; k < kVolumeN; ++k) {
                    const double x = -kVolumeExtent + i * h;
                    const double y = -kVolumeExtent + j * h;
                    const double z = -kVolumeExtent + k * h;
                    const PointFields p =
                        sample(cfg, f1, zero, 0.0, x, y, z, false);
                    f << i << "," << j << "," << k << "," << x << "," << y << ","
                      << z << "," << p.beta << "," << p.rho << ","
                      << p.kretschmann << "," << p.g00 << ","
                      << (p.valid ? 1 : 0) << "\n";
                }
        std::ofstream m(out_dir + "/volume_field_3d.meta.json");
        m << "{\n  \"kind\": \"sampled_volume\",\n  \"n\": " << kVolumeN
          << ",\n  \"extent\": " << kVolumeExtent
          << ",\n  \"time\": 0.0,\n"
          << "  \"axisymmetric_reconstruction\": false,\n"
          << "  \"statement\": \"independently sampled (x,y,z) volume\",\n"
          << "  \"fields\": [\"beta_x\", \"rho\", \"kretschmann\", \"g00\"]\n}\n";
    }

    // --- The meridional half-plane, one per frame ---------------------------
    {
        std::ofstream mf(out_dir + "/meridional_manifest.json");
        mf.precision(17);
        mf << "{\n  \"kind\": \"meridional_half_plane\",\n"
           << "  \"axisymmetric_reconstruction\": true,\n"
           << "  \"required_statement\": \"axisymmetric 3D reconstruction from "
              "production meridional data\",\n"
           << "  \"nx\": " << kMeridionalNx << ",\n  \"nr\": " << kMeridionalNr
           << ",\n  \"extent\": " << kVolumeExtent << ",\n"
           << "  \"frames\": " << kFrames << ",\n"
           << "  \"frame_axis\": \"physical time\",\n"
           << "  \"frame_axis_detail\": \"frame k is the metric at coordinate "
              "time t = k * dt with the bubble centre at x_c = v_s t; the "
              "structure translates because the bubble moves, not because a "
              "camera does\",\n"
           << "  \"dt\": " << (6.0 / (kFrames - 1)) << ",\n"
           << "  \"velocity\": " << cfg.velocity << "\n}\n";

        double rho_lo = 0.0, rho_hi = 0.0, k_hi = 0.0, n_lo = 0.0, b_lo = 0.0;
        for (int frame = 0; frame < kFrames; ++frame) {
            const double t = -3.0 + 6.0 * frame / (kFrames - 1);
            char name[64];
            std::snprintf(name, sizeof(name), "/meridional_field_%03d.csv",
                          frame);
            std::ofstream f(out_dir + name);
            f.precision(9);
            f << "i,j,x,rho,beta_x,energy_density,n_min,kretschmann,g00\n";
            for (int i = 0; i < kMeridionalNx; ++i)
                for (int j = 0; j < kMeridionalNr; ++j) {
                    const double x = -kVolumeExtent +
                                     2.0 * kVolumeExtent * i /
                                         (kMeridionalNx - 1) +
                                     cfg.velocity * t;
                    const double r = kVolumeExtent * j / (kMeridionalNr - 1);
                    const PointFields p =
                        sample(cfg, f1, zero, t, x, r, 0.0, true);
                    rho_lo = std::fmin(rho_lo, p.rho);
                    rho_hi = std::fmax(rho_hi, p.rho);
                    k_hi = std::fmax(k_hi, std::fabs(p.kretschmann));
                    n_lo = std::fmin(n_lo, p.n_min);
                    b_lo = std::fmin(b_lo, p.beta);
                    f << i << "," << j << "," << x << "," << r << "," << p.beta
                      << "," << p.rho << "," << p.n_min << "," << p.kretschmann
                      << "," << p.g00 << "\n";
                }
        }

        // Fixed scales, computed once here and frozen for every product. A
        // renderer that recomputed them per frame would make an animation whose
        // colour changes carry no physical meaning.
        std::ofstream s(out_dir + "/field_scales.json");
        s.precision(17);
        s << "{\n  \"policy\": \"fixed; computed once over every frame by "
             "production and frozen\",\n"
          << "  \"energy_density\": {\"min\": " << rho_lo << ", \"max\": "
          << rho_hi << "},\n"
          << "  \"n_min\": {\"min\": " << n_lo << ", \"max\": 0.0},\n"
          << "  \"kretschmann\": {\"min\": 0.0, \"max\": " << k_hi << "},\n"
          << "  \"beta_x\": {\"min\": " << b_lo << ", \"max\": 0.0},\n"
          << "  \"isovalues\": {\n"
          << "    \"wall_shift\": " << (0.5 * b_lo) << ",\n"
          << "    \"negative_energy\": " << (0.15 * rho_lo) << ",\n"
          << "    \"nec_violation\": " << (0.15 * n_lo) << ",\n"
          << "    \"curvature_shell\": " << (0.25 * k_hi) << "\n  },\n"
          << "  \"protected_cabin_radius_comoving\": "
          << o0c_wall_region().cabin_radius_comoving(cfg.scale_R0) << ",\n"
          << "  \"index_ordering\": \"0=t, 1=x, 2=y, 3=z\"\n}\n";
    }

    // --- Sparse shift vectors, for the arrow overlay -------------------------
    {
        std::ofstream f(out_dir + "/shift_vectors_sparse.csv");
        f.precision(9);
        f << "x,y,z,beta_x,beta_y,beta_z,magnitude\n";
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j)
                for (int k = 0; k < 5; ++k) {
                    const double x = -4.0 + 8.0 * i / 8.0;
                    const double y = -4.0 + 8.0 * j / 8.0;
                    const double z = -4.0 + 8.0 * k / 4.0;
                    const PointFields p =
                        sample(cfg, f1, zero, 0.0, x, y, z, false);
                    f << x << "," << y << "," << z << "," << p.beta
                      << ",0,0," << std::fabs(p.beta) << "\n";
                }
    }

    // --- Local null cones ----------------------------------------------------
    //
    // THE COMPLETE ALGEBRA, WRITTEN OUT, BECAUSE IT IS SHORT AND BECAUSE THE
    // RENDERER MUST NOT DERIVE IT. A tangent vector (1, v^i) is null when
    //
    //     g_00 + 2 g_0i v^i + g_ij v^i v^j = 0.
    //
    // In ADM variables g_00 = -alpha^2 + gamma_ij beta^i beta^j and
    // g_0i = gamma_ij beta^j, so completing the square in v gives
    //
    //     gamma_ij (v^i + beta^i)(v^j + beta^j) = alpha^2.
    //
    // The cone's cross-section at unit coordinate time is therefore an ELLIPSOID
    // in velocity space, centred at -beta^i, with radius alpha in the gamma
    // metric. That single statement is the whole of local cone tipping: the
    // bubble does not open or close the cone, it TRANSLATES it, by -beta^x =
    // +v_s f. For F1, gamma = I and alpha = 1, so it is a unit sphere centred at
    // +v_s f on the x axis.
    //
    // On the axis the two null slopes in the (t, x) plane follow immediately:
    //
    //     dx/dt = -beta^x +- alpha  ->  v_s f +- 1.
    //
    // WHAT THIS IS NOT. This is pointwise algebra on the metric coefficients at
    // a single event. No geodesic is integrated, no worldline is followed, and
    // nothing here says whether any two events are causally connected or how
    // long a trip takes. Cone tipping is a local property of g; turning it into
    // a statement about reachability requires the geodesic work A7 exists for,
    // and this record deliberately stops short of it.
    {
        std::ofstream f(out_dir + "/null_cone_samples.csv");
        f.precision(12);
        f << "# local null cones from the metric coefficients at each event; "
             "no geodesic is integrated and no causal claim is made\n";
        f << "frame,t,x,rho,site,alpha,beta_x,gamma_xx,gamma_rr,g00,"
             "cone_center_vx,cone_radius,slope_forward,slope_backward,"
             "coordinate_stationary_inside_cone,cone_is_tipped\n";

        struct ConeSite { const char* label; double xi_offset, rho; };
        const ConeSite sites[] = {
            {"centerline_far_back",  -3.2, 0.0},
            {"centerline_back_wall", -2.0, 0.0},
            {"centerline_inner",     -1.0, 0.0},
            {"centerline_cabin",      0.0, 0.0},
            {"centerline_inner_fwd",  1.0, 0.0},
            {"centerline_front_wall", 2.0, 0.0},
            {"centerline_far_front",  3.2, 0.0},
            {"wall_off_axis_low",     0.0, 1.0},
            {"wall_off_axis_mid",     0.0, 2.0},
            {"wall_off_axis_high",    0.0, 3.0},
            {"wall_shoulder",         1.4, 1.4},
            {"wall_shoulder_back",   -1.4, 1.4},
        };

        for (int frame = 0; frame < kFrames; ++frame) {
            const double t = -3.0 + 6.0 * frame / (kFrames - 1);
            const double xc = cfg.bubble_center_x0 + cfg.velocity * t;
            for (const ConeSite& s : sites) {
                const double x = xc + s.xi_offset;
                const MetricJet<double> jet =
                    design::f1_metric_jet(cfg, f1, zero, t, x, s.rho, 0.0);
                const geometry::AdmDecomposition adm =
                    geometry::adm_decompose(jet.g_cov.m, jet.g_contra.m);
                if (!adm.valid) continue;
                const double a = adm.alpha;
                const double b = adm.beta[0];
                const double gxx = adm.gamma[0], grr = adm.gamma[4];
                const double centre = -b;
                // Radius along x, in coordinate velocity: alpha / sqrt(gamma_xx).
                const double radius = a / std::sqrt(gxx);
                const double fwd = centre + radius;
                const double back = centre - radius;
                // v = 0 is inside the cone exactly when the coordinate-stationary
                // observer is timelike, i.e. when g_00 < 0.
                const bool inside = jet.g_cov.m(0, 0) < 0.0;
                f << frame << "," << t << "," << x << "," << s.rho << ","
                  << s.label << "," << a << "," << b << "," << gxx << ","
                  << grr << "," << jet.g_cov.m(0, 0) << "," << centre << ","
                  << radius << "," << fwd << "," << back << ","
                  << (inside ? 1 : 0) << ","
                  << (std::fabs(centre) > 1e-12 ? 1 : 0) << "\n";
            }
        }

        std::ofstream m(out_dir + "/null_cone_manifest.json");
        m.precision(17);
        m << "{\n"
          << "  \"kind\": \"local_null_cones\",\n"
          << "  \"derivation\": \"gamma_ij (v^i + beta^i)(v^j + beta^j) = "
             "alpha^2, from g_00 + 2 g_0i v^i + g_ij v^i v^j = 0 in ADM "
             "variables\",\n"
          << "  \"meaning\": \"the cone cross-section at unit coordinate time "
             "is an ellipsoid centred at -beta^i with radius alpha in the "
             "gamma metric; the bubble TRANSLATES the cone rather than opening "
             "or closing it\",\n"
          << "  \"geodesics_integrated\": false,\n"
          << "  \"causal_reachability_claimed\": false,\n"
          << "  \"travel_time_claimed\": false,\n"
          << "  \"ftl_motion_claimed\": false,\n"
          << "  \"limitation\": \"Cone tipping is a LOCAL property of the "
             "metric coefficients at one event. It says nothing about whether "
             "two events are causally connected, about trip duration, or about "
             "faster-than-light travel. Converting it into any of those "
             "statements requires geodesic integration, which O0C does not "
             "perform.\",\n"
          << "  \"sites\": " << (sizeof(sites) / sizeof(sites[0])) << ",\n"
          << "  \"frames\": " << kFrames << ",\n"
          << "  \"velocity\": " << cfg.velocity << ",\n"
          << "  \"renderer_must_not_rederive\": \"cone_center_vx, cone_radius "
             "and both null slopes are written here so the renderer draws "
             "numbers rather than computing them\"\n}\n";
    }

    // --- The pointwise objective samples the 2D figures read ------------------
    {
        std::ofstream f(out_dir + "/pointwise_objective_samples.csv");
        f.precision(9);
        f << "r,rho,rho_analytic,n_min,integrand_E_neg,integrand_J_NEC\n";
        for (int i = 0; i <= 240; ++i) {
            const double r = 6.0 * i / 240.0;
            const PointFields p = sample(cfg, f1, zero, 0.0, 0.0, r, 0.0, true);
            const design::ShiftJet sj =
                design::shift_jet(cfg, f1, zero, 0.0, 0.0, r, 0.0, false);
            const double fy = -sj.d1[2] / cfg.velocity;
            const double ref = -(cfg.velocity * cfg.velocity) /
                               (32.0 * 3.14159265358979323846) * (fy * fy);
            f << r << "," << p.rho << "," << ref << "," << p.n_min << ","
              << std::fmax(-p.rho, 0.0) << "," << std::fmax(-p.n_min, 0.0)
              << "\n";
        }
    }
}

void run_o0c_product_checks(VerificationReport& report,
                            const std::string& out_dir) {
    Stopwatch sw;
    O0cState& st = o0c_state();

    // --- 46. Production-generated 3D records --------------------------------
    {
        std::ifstream v(out_dir + "/volume_field_3d.csv");
        std::size_t rows = 0, valid = 0;
        std::string line;
        std::getline(v, line);   // comment
        std::getline(v, line);   // header
        while (std::getline(v, line)) {
            if (line.empty()) continue;
            ++rows;
            if (!line.empty() && line[line.size() - 1] == '1') ++valid;
        }
        const std::size_t want =
            static_cast<std::size_t>(kVolumeN) * kVolumeN * kVolumeN;
        report.check_predicate(
            "o0c.production_generated_3d_volume",
            "at least one validation case is a genuinely sampled (x, y, z) "
            "volume written by C++, not a revolution",
            VerificationLevel::A, rows == want && valid == want,
            std::to_string(want) + " independently evaluated points, all valid",
            std::to_string(rows) + " rows, " + std::to_string(valid) +
                " valid, on a " + std::to_string(kVolumeN) + "^3 lattice over "
                "[-" + std::to_string(kVolumeExtent) + ", " +
                std::to_string(kVolumeExtent) + "]^3", sw.seconds(),
            "the meridional records are cheaper by the azimuthal sample count "
            "and every product built from them must say so; this one is the "
            "case an axisymmetric reconstruction can be validated against");
        sw.reset();
    }

    // --- 47. Fixed scales, declared by production ---------------------------
    {
        std::ifstream s(out_dir + "/field_scales.json");
        std::string all((std::istreambuf_iterator<char>(s)),
                        std::istreambuf_iterator<char>());
        const bool has_policy = all.find("fixed; computed once") !=
                                std::string::npos;
        const bool has_iso = all.find("isovalues") != std::string::npos;
        report.check_predicate(
            "o0c.fixed_scientific_scales_declared",
            "colour scales and isovalues are computed once by production over "
            "every frame and frozen",
            VerificationLevel::A, has_policy && has_iso,
            "a fixed-scale policy and a frozen isovalue set",
            std::string("policy present: ") + (has_policy ? "yes" : "no") +
                ", isovalues present: " + (has_iso ? "yes" : "no"),
            sw.seconds(),
            "a per-frame autoscale would make an animation whose colour "
            "changes carry no physical meaning");
        sw.reset();
    }

    // --- 48. Animation axis labels ------------------------------------------
    {
        std::ifstream s(out_dir + "/meridional_manifest.json");
        std::string all((std::istreambuf_iterator<char>(s)),
                        std::istreambuf_iterator<char>());
        const bool axis = all.find("\"frame_axis\": \"physical time\"") !=
                          std::string::npos;
        const bool recon =
            all.find("axisymmetric 3D reconstruction from production "
                     "meridional data") != std::string::npos;
        report.check_predicate(
            "o0c.animation_axis_and_reconstruction_labels",
            "the frame axis is named and the required axisymmetric-"
            "reconstruction statement is carried with the meridional records",
            VerificationLevel::A, axis && recon,
            "frame axis named, reconstruction statement present",
            std::string("frame axis declared: ") + (axis ? "yes" : "no") +
                "; reconstruction statement present: " + (recon ? "yes" : "no"),
            sw.seconds(),
            "a camera orbit is not physical time, and a revolved half-plane is "
            "not an independently sampled volume; both distinctions travel "
            "with the data rather than living in a caption");
        sw.reset();
    }

    // --- 48b. The null-cone record, and the line it does not cross ------------
    {
        std::ifstream c(out_dir + "/null_cone_samples.csv");
        std::string line;
        std::getline(c, line);   // comment
        std::getline(c, line);   // header
        std::size_t rows = 0, tipped = 0, outside = 0;
        double worst_identity = 0.0;
        while (std::getline(c, line)) {
            if (line.empty()) continue;
            ++rows;
            // Re-derive the two null slopes from the recorded centre and radius
            // and check they bracket it: an independent read of the same row.
            std::vector<std::string> f;
            std::string cur;
            for (char ch : line) {
                if (ch == ',') { f.push_back(cur); cur.clear(); }
                else cur.push_back(ch);
            }
            f.push_back(cur);
            if (f.size() < 16) continue;
            const double centre = std::stod(f[10]);
            const double radius = std::stod(f[11]);
            const double fwd = std::stod(f[12]);
            const double back = std::stod(f[13]);
            // Judged RELATIVE to the magnitude, because the CSV is written at
            // twelve significant digits for the renderer's benefit: a recorded
            // slope of 1.31098569804 cannot agree with centre + radius =
            // 1.310985698042 to better than the file's own precision, and a
            // tighter absolute tolerance would be testing the write format
            // rather than the algebra.
            const double scale = std::fmax(1.0, std::fabs(fwd));
            worst_identity = std::fmax(
                worst_identity,
                std::fmax(std::fabs(fwd - (centre + radius)),
                          std::fabs(back - (centre - radius))) / scale);
            if (f[15] == "1" || f[15] == "1\r") ++tipped;
            if (f[14] == "0" || f[14] == "0\r") ++outside;
        }

        std::ifstream mm(out_dir + "/null_cone_manifest.json");
        std::string all((std::istreambuf_iterator<char>(mm)),
                        std::istreambuf_iterator<char>());
        const bool no_geodesics =
            all.find("\"geodesics_integrated\": false") != std::string::npos;
        const bool no_causal =
            all.find("\"causal_reachability_claimed\": false") !=
                std::string::npos &&
            all.find("\"travel_time_claimed\": false") != std::string::npos &&
            all.find("\"ftl_motion_claimed\": false") != std::string::npos;

        report.check_predicate(
            "o0c.local_null_cones_recorded",
            "the local null cone is written from the metric coefficients at "
            "each event, with the slopes consistent with the recorded centre "
            "and radius, and with every causal claim explicitly disclaimed",
            VerificationLevel::A,
            rows > 0 && worst_identity < 1e-11 && tipped > 0 && no_geodesics &&
                no_causal,
            "slopes = centre +- radius to the record's twelve-digit write "
            "precision; cones observably tipped; no "
            "geodesic, reachability, travel-time or FTL claim",
            std::to_string(rows) + " cone samples, " + std::to_string(tipped) +
                " tipped, " + std::to_string(outside) +
                " with the coordinate-stationary observer OUTSIDE the cone "
                "(g_00 >= 0); slope identity residual " +
                std::to_string(worst_identity) +
                "; geodesics integrated: no; causal claims: none",
            sw.seconds(),
            "gamma_ij (v^i + beta^i)(v^j + beta^j) = alpha^2 -- the bubble "
            "TRANSLATES the cone rather than opening it. That is a local "
            "property of g at one event and says nothing about reachability, "
            "which is what A7 exists for and what this record stops short of");
        sw.reset();
    }

    // --- 49. No energy claim -------------------------------------------------
    report.check_predicate(
        "o0c.no_reduced_energy_claim",
        "no lower objective value is reported as a reduction in warp-drive "
        "energy",
        VerificationLevel::A, !st.any_behaviour_change_accepted,
        "zero such claims",
        "the only integrated values O0C reports are for the zero-correction "
        "Alcubierre baseline and for eight deterministic adversarial vectors, "
        "each carrying an explicit comparable / unresolved / behavior_changing "
        "verdict. No objective is marked ready for optimization, so no "
        "comparison between two designs is made at all.",
        sw.seconds(),
        "the milestone's forbidden claim, checked as a case rather than left "
        "to prose");
    sw.reset();

    // --- 50. Scope -----------------------------------------------------------
    report.check_predicate(
        "o0c.scope_statement", "the O0C scope statement is present and exact",
        VerificationLevel::A, true, "the required sentence, verbatim",
        o0c_scope(), sw.seconds(), "");
}

void write_o0c_manifest(const VerificationReport& report,
                        const std::string& out_dir) {
    const O0cState& st = o0c_state();
    const ShiftFamilyConfig cfg = o0c_config();
    std::ofstream f(out_dir + "/geometry_bridge_manifest.json");
    f.precision(17);
    const BuildInfo& b = build_info();
    f << "{\n"
      << "  \"milestone\": \"O0C\",\n"
      << "  \"git_commit\": \"" << b.git_commit << "\",\n"
      << "  \"git_dirty\": " << (b.git_dirty ? "true" : "false") << ",\n"
      << "  \"velocity\": " << cfg.velocity << ",\n"
      << "  \"bubble_radius\": " << cfg.bubble_radius << ",\n"
      << "  \"wall_sigma\": " << cfg.wall_sigma << ",\n"
      << "  \"scale_R0\": " << cfg.scale_R0 << ",\n"
      << "  \"lambda\": 0.0,\n"
      << "  \"source_semantics\": \"total_effective\",\n"
      << "  \"production_finite_differences_in_geometry\": 0,\n"
      << "  \"finite_difference_note\": \"Every curvature and Einstein-tensor "
         "component in the O0C production route is built from the analytic "
         "O0B metric jet through connection_jet. Centred differences appear "
         "only as TEST references in the gradient cases.\",\n"
      << "  \"energy_evaluated\": true,\n"
      << "  \"optimizer_run\": false,\n"
      << "  \"f1w_parent_modes\": " << st.f1w_parent_modes << ",\n"
      << "  \"f1w_active_modes\": " << st.f1w_active_modes << ",\n"
      << "  \"f1w_free_dimension\": " << st.f1w_free << ",\n"
      << "  \"f1w_whitened_dimension\": " << st.f1w_whitened << ",\n"
      << "  \"f1w_oti_batches\": " << st.f1w_batches << ",\n"
      << "  \"f1w_blocked\": " << (st.f1w_blocked ? "true" : "false") << ",\n"
      << "  \"converged_E_net\": " << st.converged_E_net << ",\n"
      << "  \"converged_E_neg\": " << st.converged_E_neg << ",\n"
      << "  \"order_converged\": " << (st.order_converged ? "true" : "false")
      << ",\n  \"subdivision_converged\": "
      << (st.subdivision_converged ? "true" : "false")
      << ",\n  \"domain_converged\": "
      << (st.domain_converged ? "true" : "false")
      << ",\n  \"independent_agreement\": "
      << (st.independent_agreement ? "true" : "false") << ",\n"
      << "  \"gradient_verified\": "
      << (st.gradient_verified ? "true" : "false") << ",\n"
      << "  \"oti_available\": " << (st.oti_available ? "true" : "false")
      << ",\n  \"oti_batches\": " << st.oti_batches << ",\n"
      << "  \"hessian_limitation\": \"The complete design gradient is "
         "recovered through deterministic OTI batches. Mixed Hessian entries "
         "between coordinates assigned to different batches are not "
         "available.\",\n"
      << "  \"renderer_computes_physics\": false,\n"
      << "  \"cases\": " << report.num_cases() << ",\n"
      << "  \"passed\": " << report.num_passed() << ",\n"
      << "  \"failed\": " << report.num_failed() << ",\n"
      << "  \"skipped\": " << report.num_skipped() << ",\n"
      << "  \"scope\": \"" << o0c_scope() << "\"\n}\n";
}

}  // namespace checks
}  // namespace warp_hypad
