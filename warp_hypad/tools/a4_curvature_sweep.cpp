#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "warp_hypad/core/metadata.hpp"
#include "warp_hypad/geometry/curvature.hpp"
#include "warp_hypad/metrics/metric.hpp"
#include "warp_hypad/metrics/metric_field.hpp"

// ---------------------------------------------------------------------------
// A4 curvature sweep exporter.
//
// Sweeps the bubble velocity through and past v_s = 1 and records, for every
// frame, the z = 0 slice of the metric and of the curvature invariants. The
// renderer reads these files and nothing else.
//
// WHY THIS EXPORTER EXISTS AND THE RENDERER DOES NOT COMPUTE.
// The project's rule is that a figure is redrawn from saved output without
// rerunning the solver. Curvature is the strongest case for that rule so far:
// it needs the inverse metric, two nested derivative passes and a temporal
// resampling policy, and a Python reimplementation of any of those would drift
// from the verified C++ silently.
//
// WHAT v_s = 1 IS AND IS NOT.
// It is the velocity at which g_00 first touches zero at the bubble centre,
// because g_00 = -1 + v_s^2 f^2 and f(0) = 1 exactly. Following the V0.1A
// convention this milestone distinguishes:
//
//   exact_onset_frame           v_s = 1 exactly; g_00(centre) = 0 exactly
//   first_spatial_contour_frame the first frame whose SAVED GRID brackets zero
//                               and therefore supports an interpolated contour
//
// It is NOT metric degeneracy, NOT signature loss and NOT automatically a
// horizon. det g = -1 identically for every v_s and f, and the curvature stays
// finite through the crossing -- which is a claim this sweep makes checkable
// rather than one it asserts.
//
// Build-local output only. Nothing here writes into the source tree.
// ---------------------------------------------------------------------------

namespace {

using namespace warp_hypad;

constexpr int kGridN = 25;          // 25^3; a Rank4Field at 33^3 is 74 MB a copy
constexpr double kHalfWidth = 4.0;
constexpr double kRadius = 2.0;
constexpr double kSigma = 1.0;
constexpr double kTime = 0.0;

// v_s = k/10 for k = 2..16. k = 10 gives EXACTLY 1.0 in binary floating point,
// which is what makes the onset frame an exact statement rather than a nearby
// one. A sequence built by repeated addition of 0.1 would not contain 1.0.
constexpr int kFirstK = 2;
constexpr int kLastK = 16;

double velocity_for(int k) { return static_cast<double>(k) / 10.0; }

}  // namespace

int main(int argc, char** argv) {
    std::string out_dir = "results/verification/v0_program_a4_curvature";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: a4_curvature_sweep [--out <dir>]\n";
            return 2;
        }
    }

    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec) {
        std::cerr << "cannot create " << out_dir << ": " << ec.message() << "\n";
        return 3;
    }

    GridSpec spec;
    spec.nx = spec.ny = spec.nz = kGridN;
    spec.x_min = spec.y_min = spec.z_min = -kHalfWidth;
    spec.x_max = spec.y_max = spec.z_max = kHalfWidth;
    const Grid grid(spec);

    const int km = grid.center_index_z();
    const int im = grid.center_index_x();
    const int jm = grid.center_index_y();
    if (km < 0 || im < 0 || jm < 0) {
        std::cerr << "grid extents must be odd so the centre is pinned\n";
        return 3;
    }

    std::ofstream summary(out_dir + "/curvature_velocity_sweep.csv");
    summary << "frame,v_s,g00_center,g01_center,ricci_scalar_center,"
               "kretschmann_center,max_abs_ricci_scalar,max_abs_kretschmann,"
               "max_abs_riemann,grid_brackets_g00_zero,curvature_finite\n";

    int exact_onset_frame = -1;
    int first_contour_frame = -1;
    bool all_finite = true;

    int frame = 0;
    for (int k = kFirstK; k <= kLastK; ++k, ++frame) {
        const double v_s = velocity_for(k);

        metrics::AlcubierreMetric<double> metric(
            0.0, metrics::TimeTreatment::InertialTimeDependent);
        metrics::AlcubierreParameters<double> params(v_s, kRadius, kSigma);

        auto sample = [&](double t) {
            return metrics::build_metric_field<double>(grid, metric, params, t);
        };

        // h_t = h_x / v_s aligns the temporal taps with the spatial stencil
        // along u = x - v_s t, which makes d_0 g = -v_s d_x g exact rather than
        // second order. See the A4 suite's rigid_translation_identity case.
        const auto policy = geometry::TemporalPolicy::analytic_resample(
            grid.dx() / v_s,
            "Alcubierre depends on t; h_t = h_x/v_s aligns the stencils");

        auto c = geometry::compute_curvature<double>(grid, sample, kTime, policy);
        if (!c.ok) {
            std::cerr << "frame " << frame << " (v_s=" << v_s
                      << ") refused: " << c.refusal_reason << "\n";
            return 1;
        }

        char name[64];
        std::snprintf(name, sizeof(name), "/frame_%03d.csv", frame);
        std::ofstream f(out_dir + name);
        f << "x,y,g00,g01,ricci_scalar,ricci_square,kretschmann,riemann_max,"
             "christoffel_max\n";

        double max_R = 0.0, max_K = 0.0, max_riem = 0.0;
        bool saw_positive = false, saw_negative = false;

        for (int i = 0; i < grid.nx(); ++i) {
            for (int j = 0; j < grid.ny(); ++j) {
                const std::size_t p = grid.flat(i, j, km);
                const double g00 = c.g.at(p, 0, 0);
                const double R = c.ricci_scalar.at(p);
                const double K = c.kretschmann.at(p);
                const double riem = max_abs_real(c.riemann.tensor_at(p));

                if (g00 > 0.0) saw_positive = true;
                if (g00 < 0.0) saw_negative = true;

                max_R = std::fmax(max_R, std::fabs(R));
                max_K = std::fmax(max_K, std::fabs(K));
                max_riem = std::fmax(max_riem, riem);

                if (!std::isfinite(g00) || !std::isfinite(R) ||
                    !std::isfinite(K) || !std::isfinite(riem)) {
                    all_finite = false;
                }

                f << grid.x(i) << "," << grid.y(j) << "," << g00 << ","
                  << c.g.at(p, 0, 1) << "," << R << ","
                  << c.ricci_square.at(p) << "," << K << "," << riem << ","
                  << max_abs_real(c.gamma.tensor_at(p)) << "\n";
            }
        }

        const std::size_t centre = grid.flat(im, jm, km);
        const double g00_centre = c.g.at(centre, 0, 0);

        // The exact onset: v_s = 1 makes g_00 = -1 + v_s^2 f(0)^2 = 0 at the
        // centre, exactly, because f(0) = 1 identically. Tested with ==
        // deliberately -- this is an exactness claim, not a tolerance.
        if (v_s == 1.0 && g00_centre == 0.0 && exact_onset_frame < 0)
            exact_onset_frame = frame;

        const bool brackets = saw_positive && saw_negative;
        if (brackets && first_contour_frame < 0) first_contour_frame = frame;

        summary << frame << "," << v_s << "," << g00_centre << ","
                << c.g.at(centre, 0, 1) << "," << c.ricci_scalar.at(centre) << ","
                << c.kretschmann.at(centre) << "," << max_R << "," << max_K << ","
                << max_riem << "," << (brackets ? 1 : 0) << ","
                << (all_finite ? 1 : 0) << "\n";

        std::cout << "frame " << frame << "  v_s=" << v_s
                  << "  g00(centre)=" << g00_centre << "  max|R|=" << max_R
                  << "  max|K|=" << max_K << "\n";
    }

    const int n_frames = frame;

    std::ofstream man(out_dir + "/a4_curvature_manifest.json");
    man << "{\n";
    man << "  \"git_commit\": \"" << build_info().git_commit << "\",\n";
    man << "  \"git_dirty\": " << (build_info().git_dirty ? "true" : "false")
        << ",\n";
    man << "  \"generated_utc\": \"" << utc_timestamp() << "\",\n";
    man << "  \"grid\": \"" << kGridN << "x" << kGridN << "x" << kGridN
        << " over [-" << kHalfWidth << "," << kHalfWidth << "]^3\",\n";
    man << "  \"slice\": \"z = 0\",\n";
    man << "  \"radius\": " << kRadius << ",\n";
    man << "  \"sigma\": " << kSigma << ",\n";
    man << "  \"time\": " << kTime << ",\n";
    man << "  \"frame_count\": " << n_frames << ",\n";
    man << "  \"velocity_sequence\": \"v_s = k/10 for k = " << kFirstK << ".."
        << kLastK << "\",\n";
    man << "  \"exact_onset_frame\": " << exact_onset_frame << ",\n";
    man << "  \"first_spatial_contour_frame\": " << first_contour_frame << ",\n";
    man << "  \"onset_meaning\": \"At the exact onset v_s = 1 the lapse term "
           "g_00 = -1 + v_s^2 f^2 touches zero at the bubble centre because "
           "f(0) = 1 exactly. The saved grid need not yet bracket zero, so an "
           "interpolated spatial contour is only supported from the later "
           "frame. This is NOT metric degeneracy, NOT signature loss and NOT "
           "automatically a horizon: det g = -1 identically for every v_s.\",\n";
    man << "  \"curvature_finite_through_crossing\": "
        << (all_finite ? "true" : "false") << ",\n";
    man << "  \"temporal_policy\": \"analytic_resample with h_t = h_x / v_s\",\n";
    man << "  \"scope\": \"These outputs visualize the curvature of a chosen "
           "metric. They do not demonstrate a stress-energy distribution, an "
           "energy condition, or any statement about whether this geometry can "
           "be sourced. No quantity here may be called an energy.\"\n";
    man << "}\n";

    std::cout << "\nframes: " << n_frames
              << "\nexact_onset_frame: " << exact_onset_frame
              << "\nfirst_spatial_contour_frame: " << first_contour_frame
              << "\ncurvature finite through crossing: "
              << (all_finite ? "yes" : "NO") << "\n";

    if (exact_onset_frame < 0) {
        std::cerr << "the sweep did not contain v_s = 1 exactly\n";
        return 1;
    }
    if (!all_finite) {
        std::cerr << "non-finite curvature encountered\n";
        return 1;
    }
    return 0;
}
