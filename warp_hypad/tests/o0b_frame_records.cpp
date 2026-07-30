#include "o0b_checks.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// O0B animation frame records.
//
// Production computes every frame. A Python renderer could form a metric from
// saved coefficients by multiplying P and W and summing basis values, and that
// would be an unverified second implementation of the design layer -- exactly
// what the visualization contract forbids. So the frames are written here, from
// the same verified path the gate just checked, and the renderer only draws.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using design::AxisymmetricBasis;
using design::MetricJet;
using design::ShiftFamilyConfig;

struct Slice {
    int nx = 61, ny = 25;
    double x_lo = -6.0, x_hi = 6.0, y_lo = 0.0, y_hi = 5.0;
    double x(int i) const { return x_lo + (x_hi - x_lo) * i / (nx - 1.0); }
    double y(int j) const { return y_lo + (y_hi - y_lo) * j / (ny - 1.0); }
};

}  // namespace

void write_o0b_frame_records(const std::string& out_dir) {
    const Slice sl;
    const ShiftFamilyConfig cfg = o0b_config();
    const design::DesignRegion region = o0b_region();
    AxisymmetricBasis basis = o0b_basis();
    const design::DesignQuadrature quad = design::build_quadrature(basis);
    const design::ConstraintProjection pr =
        design::build_projection(o0b_constraints(basis, quad));
    const design::DenseMatrix M = design::sobolev_gram(
        basis, quad, design::SobolevScales{cfg.scale_R0, cfg.scale_R0});
    const design::WhiteningResult wr = design::whiten(M, pr);
    const design::CoefficientSpaces spaces(pr, wr);
    const std::vector<double> zero(
        static_cast<std::size_t>(basis.coefficient_count()), 0.0);

    // --- one frame per active mode ------------------------------------------
    {
        std::ofstream f(out_dir + "/frames_single_shift_mode.csv");
        f.precision(10);
        f << "mode,x,y,beta_x,beta_baseline,beta_delta,g00,g01\n";
        for (int k = 0; k < basis.coefficient_count(); ++k) {
            std::vector<double> c(zero.size(), 0.0);
            c[static_cast<std::size_t>(k)] = 0.08;
            for (int i = 0; i < sl.nx; ++i)
                for (int j = 0; j < sl.ny; ++j) {
                    const double x = sl.x(i), y = sl.y(j);
                    const MetricJet<double> jc = design::f1_metric_jet(
                        cfg, basis, c, 0.0, x, y, 0.0, true);
                    const MetricJet<double> jb = design::f1_metric_jet(
                        cfg, basis, zero, 0.0, x, y, 0.0, false);
                    f << k << "," << x << "," << y << "," << jc.g_cov.m(0, 1)
                      << "," << jb.g_cov.m(0, 1) << ","
                      << (jc.g_cov.m(0, 1) - jb.g_cov.m(0, 1)) << ","
                      << jc.g_cov.m(0, 0) << "," << jc.g_cov.m(0, 1) << "\n";
                }
        }
    }

    // --- amplitude sweep: baseline -> corrected ------------------------------
    {
        std::ofstream f(out_dir + "/frames_amplitude_sweep.csv");
        f.precision(10);
        f << "frame,amplitude,x,y,beta_x,g00,g01,det_tx,inverse_residual\n";
        const std::vector<double> base =
            o0b_coefficients(basis.coefficient_count(), 0.05);
        for (int t = 0; t <= 12; ++t) {
            const double amp = t / 12.0;
            std::vector<double> c(base.size());
            for (std::size_t i = 0; i < base.size(); ++i) c[i] = amp * base[i];
            for (int i = 0; i < sl.nx; ++i)
                for (int j = 0; j < sl.ny; ++j) {
                    const double x = sl.x(i), y = sl.y(j);
                    const MetricJet<double> jj = design::f1_metric_jet(
                        cfg, basis, c, 0.0, x, y, 0.0, true);
                    const double d = jj.g_cov.m(0, 0) * jj.g_cov.m(1, 1) -
                                     jj.g_cov.m(0, 1) * jj.g_cov.m(1, 0);
                    f << t << "," << amp << "," << x << "," << y << ","
                      << jj.g_cov.m(0, 1) << "," << jj.g_cov.m(0, 0) << ","
                      << jj.g_cov.m(0, 1) << "," << d << ","
                      << design::inverse_residual(jj.g_cov, jj.g_contra) << "\n";
                }
        }
    }

    // --- analytic time translation ------------------------------------------
    {
        std::ofstream f(out_dir + "/frames_time_translation.csv");
        f.precision(10);
        f << "frame,t,bubble_center_x,x,y,beta_x,g00,d_t_g01,d_tt_g01\n";
        const std::vector<double> c =
            o0b_coefficients(basis.coefficient_count(), 0.02);
        for (int t = 0; t <= 14; ++t) {
            const double tt = -3.0 + 6.0 * t / 14.0;
            const double xc = cfg.bubble_center_x0 + cfg.velocity * tt;
            for (int i = 0; i < sl.nx; ++i)
                for (int j = 0; j < sl.ny; ++j) {
                    const double x = sl.x(i), y = sl.y(j);
                    const MetricJet<double> jj = design::f1_metric_jet(
                        cfg, basis, c, tt, x, y, 0.0, true);
                    f << t << "," << tt << "," << xc << "," << x << "," << y
                      << "," << jj.g_cov.m(0, 1) << "," << jj.g_cov.m(0, 0)
                      << "," << jj.dg[0].m(0, 1) << "," << jj.d2g[0][0].m(0, 1)
                      << "\n";
                }
        }
    }

    // --- the three coordinate systems ----------------------------------------
    {
        std::ofstream f(out_dir + "/frames_coordinate_systems.csv");
        f.precision(10);
        f << "frame,system,index,x,y,beta_delta,coefficients_touched,"
             "total_coefficients\n";
        const int n = std::min(12, pr.n_free);
        for (int t = 0; t < n; ++t) {
            // The SAME index in three spaces. They are different fields, which
            // is the whole point of keeping the spaces typed and distinct.
            for (int sys = 0; sys < 3; ++sys) {
                std::vector<double> c(zero.size(), 0.0);
                if (sys == 0) {
                    if (t < static_cast<int>(c.size()))
                        c[static_cast<std::size_t>(t)] = 0.08;
                } else if (sys == 1) {
                    design::ConstrainedSplineCoordinates z;
                    z.z.assign(static_cast<std::size_t>(pr.n_free), 0.0);
                    z.z[static_cast<std::size_t>(t)] = 0.08;
                    design::RawSplineCoefficients rc;
                    spaces.to_raw(z, rc);
                    for (std::size_t i = 0; i < c.size(); ++i)
                        c[i] = rc.c[i] - pr.c0[i];
                } else {
                    design::WhitenedDesignCoordinates y;
                    y.y.assign(static_cast<std::size_t>(wr.retained), 0.0);
                    y.y[static_cast<std::size_t>(t)] = 0.08;
                    design::RawSplineCoefficients rc;
                    spaces.to_raw(y, rc);
                    for (std::size_t i = 0; i < c.size(); ++i)
                        c[i] = rc.c[i] - pr.c0[i];
                }
                double mx = 0.0;
                for (double v : c) mx = std::fmax(mx, std::fabs(v));
                int touched = 0;
                for (double v : c)
                    if (mx > 0.0 && std::fabs(v) > 1e-6 * mx) ++touched;
                const char* name = (sys == 0) ? "raw_c"
                                 : (sys == 1) ? "constrained_z" : "whitened_y";
                for (int i = 0; i < sl.nx; ++i)
                    for (int j = 0; j < sl.ny; ++j) {
                        const double x = sl.x(i), y2 = sl.y(j);
                        const MetricJet<double> jc = design::f1_metric_jet(
                            cfg, basis, c, 0.0, x, y2, 0.0, true);
                        const MetricJet<double> jb = design::f1_metric_jet(
                            cfg, basis, zero, 0.0, x, y2, 0.0, false);
                        f << t << "," << name << "," << t << "," << x << ","
                          << y2 << ","
                          << (jc.g_cov.m(0, 1) - jb.g_cov.m(0, 1)) << ","
                          << touched << "," << pr.n_raw << "\n";
                    }
            }
        }
    }

    // --- sampling comparison --------------------------------------------------
    {
        const design::UniformGridCost cost =
            design::uniform_grid_cost(0.3018, 4.0, 0.4);
        const design::QuadratureScaffold sc =
            design::build_quadrature_scaffold(basis, region);
        std::ofstream f(out_dir + "/frames_sampling_comparison.csv");
        f.precision(10);
        f << "frame,route,label,xi,s,weight_or_spacing,cumulative_points\n";
        for (int t = 0; t < 12; ++t) {
            const std::size_t upto =
                (sc.element.size() * static_cast<std::size_t>(t + 1)) / 12;
            for (std::size_t e = 0; e < upto; ++e) {
                const design::QuadratureElement& el = sc.element[e];
                f << t << ",B,element," << 0.5 * (el.xi_lower + el.xi_upper)
                  << "," << 0.5 * (el.s_lower + el.s_upper) << "," << el.volume
                  << ","
                  << (upto * static_cast<std::size_t>(
                                 sc.nodes_per_direction *
                                 sc.nodes_per_direction))
                  << "\n";
            }
            f << t << ",A,spacing," << cost.required_spacing_m << ","
              << cost.total_cells << "," << cost.current_spacing_m << ","
              << cost.total_cells << "\n";
        }
    }
}

}  // namespace checks
}  // namespace warp_hypad
