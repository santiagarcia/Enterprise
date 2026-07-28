#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "warp_hypad/core/metadata.hpp"
#include "warp_hypad/geometry/einstein.hpp"
#include "warp_hypad/metrics/metric.hpp"
#include "warp_hypad/metrics/metric_field.hpp"

// ---------------------------------------------------------------------------
// A5 Einstein-tensor sweep exporter.
//
// The same velocity sequence as the A4 curvature sweep, so the two animations
// are directly comparable frame for frame, but carrying the Einstein tensor and
// its diagnostics instead of curvature invariants.
//
// SCOPE. These are COORDINATE-BASIS components of the geometric Einstein tensor.
// They are not observer-measured energy densities and no frame of the resulting
// animation may be read as one. An energy density is T_{mu nu} u^mu u^nu for a
// specified normalised observer; A5 has no observers.
//
// COST. Every frame runs compute_einstein, which for a resampling temporal
// policy performs THREE full curvature computations -- one at the frame time and
// one on each neighbouring time slice, because the contracted Bianchi residual
// needs d_0 G^{mu nu}. That is why this sweep uses a 21^3 grid where the A4
// curvature sweep used 25^3.
// ---------------------------------------------------------------------------

namespace {

using namespace warp_hypad;

constexpr int kGridN = 21;
constexpr double kHalfWidth = 4.0;
constexpr double kRadius = 2.0;
constexpr double kSigma = 1.0;
constexpr double kTime = 0.0;

// v_s = k/10 for k = 2..16. k = 10 is EXACTLY 1.0 in binary floating point.
constexpr int kFirstK = 2;
constexpr int kLastK = 16;

}  // namespace

int main(int argc, char** argv) {
    std::string out_dir = "results/verification/v0_program_a5_einstein";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else {
            std::cerr << "usage: a5_einstein_sweep [--out <dir>]\n";
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

    const int im = grid.center_index_x();
    const int jm = grid.center_index_y();
    const int km = grid.center_index_z();
    if (im < 0 || jm < 0 || km < 0) {
        std::cerr << "grid extents must be odd so the centre is pinned\n";
        return 3;
    }

    const auto lp = geometry::CosmologicalConstantPolicy::none();

    std::ofstream summary(out_dir + "/einstein_velocity_sweep.csv");
    summary << "frame,v_s,g00_center,G00_center,G01_center,trace_center,"
               "max_abs_G,antisymmetry_residual,bianchi_absolute,"
               "bianchi_relative,grid_brackets_g00_zero,all_finite\n";

    int exact_onset_frame = -1;
    int first_contour_frame = -1;
    bool all_finite = true;

    int frame = 0;
    for (int k = kFirstK; k <= kLastK; ++k, ++frame) {
        const double v_s = static_cast<double>(k) / 10.0;

        metrics::AlcubierreMetric<double> metric(
            0.0, metrics::TimeTreatment::InertialTimeDependent);
        metrics::AlcubierreParameters<double> params(v_s, kRadius, kSigma);
        auto sample = [&](double t) {
            return metrics::build_metric_field<double>(grid, metric, params, t);
        };

        const auto policy = geometry::TemporalPolicy::analytic_resample(
            grid.dx() / v_s, "h_t = h_x/v_s aligns the temporal and spatial "
                             "stencils along u = x - v_s t");

        auto e = geometry::compute_einstein<double>(grid, sample, kTime, policy, lp);
        if (!e.ok) {
            std::cerr << "frame " << frame << " (v_s=" << v_s
                      << ") refused: " << e.refusal_reason << "\n";
            return 1;
        }

        // The metric on this slice, for g_00 and the contour test.
        const Rank2Field<double> gf =
            metrics::build_metric_field<double>(grid, metric, params, kTime);

        char name[64];
        std::snprintf(name, sizeof(name), "/frame_%03d.csv", frame);
        std::ofstream f(out_dir + name);
        f << "x,y,g00,G00_raw,G01_raw,G11_raw,G00_symmetric,G01_symmetric,"
             "G11_symmetric,G22_symmetric,einstein_trace,einstein_square,"
             "antisymmetry_residual,bianchi_nu0\n";

        double max_G = 0.0;
        bool saw_pos = false, saw_neg = false;

        for (int i = 0; i < grid.nx(); ++i) {
            for (int j = 0; j < grid.ny(); ++j) {
                const std::size_t p = grid.flat(i, j, km);
                const double g00 = gf.at(p, 0, 0);
                if (g00 > 0.0) saw_pos = true;
                if (g00 < 0.0) saw_neg = true;

                const Mat4<double> Gr = e.raw_einstein.matrix_at(p);
                const Mat4<double> Gs = e.symmetric_einstein.matrix_at(p);
                max_G = std::fmax(max_G, max_abs_real(Gs));

                if (!std::isfinite(g00) || !std::isfinite(Gs(0, 0)) ||
                    !std::isfinite(e.trace.at(p)))
                    all_finite = false;

                f << grid.x(i) << "," << grid.y(j) << "," << g00 << ","
                  << Gr(0, 0) << "," << Gr(0, 1) << "," << Gr(1, 1) << ","
                  << Gs(0, 0) << "," << Gs(0, 1) << "," << Gs(1, 1) << ","
                  << Gs(2, 2) << "," << e.trace.at(p) << "," << e.square.at(p)
                  << "," << symmetry_residual(Gr) << ","
                  << e.divergence.at(p, 0) << "\n";
            }
        }

        const std::size_t centre = grid.flat(im, jm, km);
        const double g00_centre = gf.at(centre, 0, 0);
        if (v_s == 1.0 && g00_centre == 0.0 && exact_onset_frame < 0)
            exact_onset_frame = frame;

        const bool brackets = saw_pos && saw_neg;
        if (brackets && first_contour_frame < 0) first_contour_frame = frame;

        const Mat4<double> Gc = e.symmetric_einstein.matrix_at(centre);
        summary << frame << "," << v_s << "," << g00_centre << "," << Gc(0, 0)
                << "," << Gc(0, 1) << "," << e.trace.at(centre) << "," << max_G
                << "," << e.antisymmetry_residual << "," << e.divergence_residual
                << "," << e.divergence_relative << "," << (brackets ? 1 : 0)
                << "," << (all_finite ? 1 : 0) << "\n";

        std::cout << "frame " << frame << "  v_s=" << v_s
                  << "  g00(c)=" << g00_centre << "  G00(c)=" << Gc(0, 0)
                  << "  max|G|=" << max_G
                  << "  bianchi_rel=" << e.divergence_relative << "\n";
    }

    std::ofstream man(out_dir + "/a5_einstein_manifest.json");
    man << "{\n";
    man << "  \"git_commit\": \"" << build_info().git_commit << "\",\n";
    man << "  \"git_dirty\": " << (build_info().git_dirty ? "true" : "false")
        << ",\n";
    man << "  \"generated_utc\": \"" << utc_timestamp() << "\",\n";
    man << "  \"grid\": \"" << kGridN << "^3 over [-" << kHalfWidth << ","
        << kHalfWidth << "]^3\",\n";
    man << "  \"slice\": \"z = 0\",\n";
    man << "  \"radius\": " << kRadius << ",\n";
    man << "  \"sigma\": " << kSigma << ",\n";
    man << "  \"frame_count\": " << frame << ",\n";
    man << "  \"velocity_sequence\": \"v_s = k/10 for k = " << kFirstK << ".."
        << kLastK << "\",\n";
    man << "  \"exact_onset_frame\": " << exact_onset_frame << ",\n";
    man << "  \"first_spatial_contour_frame\": " << first_contour_frame << ",\n";
    man << "  \"cosmological_constant_treatment\": \""
        << geometry::to_string(lp.treatment) << "\",\n";
    man << "  \"lambda\": " << lp.lambda << ",\n";
    man << "  \"all_finite\": " << (all_finite ? "true" : "false") << ",\n";
    man << "  \"quantity\": \"coordinate-basis components of the geometric "
           "Einstein tensor G_{mu nu}\",\n";
    man << "  \"not_an_energy_density\": \"These are coordinate-basis "
           "Einstein-tensor components. They are NOT observer-measured energy "
           "densities. An energy density is T_{mu nu} u^mu u^nu for a specified "
           "normalised observer, which Milestone A6 introduces.\",\n";
    man << "  \"scope\": \"A5 computes the geometric Einstein tensor and the "
           "stress-energy tensor inferred through Einstein's equation for the "
           "selected metric. It does not establish that this source can be "
           "generated physically, that energy conditions are satisfied, or that "
           "FTL travel is realizable.\"\n";
    man << "}\n";

    std::cout << "\nframes: " << frame
              << "\nexact_onset_frame: " << exact_onset_frame
              << "\nfirst_spatial_contour_frame: " << first_contour_frame
              << "\nall finite: " << (all_finite ? "yes" : "NO") << "\n";

    if (exact_onset_frame < 0) {
        std::cerr << "the sweep did not contain v_s = 1 exactly\n";
        return 1;
    }
    if (!all_finite) {
        std::cerr << "non-finite Einstein components encountered\n";
        return 1;
    }
    return 0;
}
