#include "o0a_checks.hpp"

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "warp_hypad/design/constraint_projection.hpp"
#include "warp_hypad/design/radial_validation.hpp"

// ---------------------------------------------------------------------------
// Animation frame records.
//
// WHY THESE EXIST AT ALL. The O0A animations show four things that are the
// RESULT of a production operation: a raw perturbation projected onto the
// constrained subspace, the field that projection produces, an F0 reconstruction
// built one coefficient at a time, and the coefficient pattern a single whitened
// coordinate excites. A Python renderer could form each of those by multiplying
// saved matrices -- and that would be a second, unverified implementation of the
// projection and of the basis evaluation, which is exactly what the
// visualization contract forbids. So production computes every frame and writes
// it; the renderer only draws.
//
// These writers add no scientific claim. They re-use the canonical basis, the
// canonical constraints and the verified projection, and they are called after
// every group has run so the numbers they contain are the numbers the gate just
// checked.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using design::AxisymmetricBasis;
using design::DesignCoordinates;

// A coarse (xi, s) slice, shared by every frame record so the animations all sit
// on one grid and a reader can compare them frame for frame.
struct SliceGrid {
    int n_xi = 41;
    int n_s = 25;
    double xi_min, xi_max, s_max;
    SliceGrid() {
        const design::DesignRegion r = canonical_region();
        xi_min = r.xi_min;
        xi_max = r.xi_max;
        s_max = r.s_max;
    }
    double xi(int i) const {
        return xi_min + (xi_max - xi_min) * static_cast<double>(i) /
                            static_cast<double>(n_xi - 1);
    }
    double s(int j) const {
        return s_max * static_cast<double>(j) / static_cast<double>(n_s - 1);
    }
};

}  // namespace

void write_o0a_frame_records(const std::string& out_dir) {
    const SliceGrid slice;
    AxisymmetricBasis basis = canonical_basis();
    const design::DesignQuadrature quad = design::build_quadrature(basis);
    const design::ConstraintSystem sys = canonical_constraints(basis, quad);
    const design::ConstraintProjection pr = design::build_projection(sys);
    const design::DenseMatrix M = design::sobolev_gram(
        basis, quad, design::SobolevScales{kO0aR0, kO0aR0});
    const design::WhiteningResult wr = design::whiten(M, pr);
    const design::DenseMatrix C = sys.matrix();

    // --- Every active tensor product, on the slice ---------------------------
    //
    // Replaces the six-mode atlas for the animation: one frame per active mode.
    {
        std::ofstream f(out_dir + "/tensor_mode_frames.csv");
        f.precision(10);
        f << "flat_index,index_xi,index_s,xi,s,value,d_xi,d_s\n";
        for (int k = 0; k < basis.coefficient_count(); ++k) {
            for (int i = 0; i < slice.n_xi; ++i)
                for (int j = 0; j < slice.n_s; ++j) {
                    const DesignCoordinates c = design::from_design_coordinates(
                        slice.xi(i), slice.s(j), kO0aR0);
                    const design::TensorBasisEvaluation e = basis.evaluate(c);
                    for (std::size_t m = 0; m < e.size(); ++m) {
                        if (e.flat_index[m] != k) continue;
                        f << k << "," << e.index_xi[m] << "," << e.index_s[m]
                          << "," << slice.xi(i) << "," << slice.s(j) << ","
                          << e.value[m] << "," << e.d_xi[m] << "," << e.d_s[m]
                          << "\n";
                    }
                }
        }
    }

    // --- Per-mode statistics, for the individual-mode animations ------------
    {
        std::ofstream f(out_dir + "/mode_statistics.csv");
        f.precision(10);
        f << "family,mode_index,support_lower,support_upper,status,"
             "max_abs_value,max_abs_first,max_abs_second,max_abs_third\n";
        struct Which {
            const char* name;
            const design::QuinticBSplineBasis1D* b;
            const design::LatticeClassification* c;
            double lo, hi;
        };
        const Which which[2] = {
            {"axial_xi", &basis.axial(), &basis.axial_class(), slice.xi_min,
             slice.xi_max},
            {"radial_s", &basis.radial(), &basis.radial_class(), 0.0,
             slice.s_max}};
        for (const Which& w : which)
            for (int i = 0; i < w.b->size(); ++i) {
                bool raw = false;
                for (int r : w.c->raw)
                    if (r == i) raw = true;
                if (!raw) continue;
                double v = 0, d1 = 0, d2 = 0, d3 = 0;
                for (int k = 0; k <= 400; ++k) {
                    const double t = w.lo + (w.hi - w.lo) * k / 400.0;
                    v = std::fmax(v, std::fabs(w.b->value_of(i, t)));
                    d1 = std::fmax(d1, std::fabs(w.b->derivative_of(i, t, 1)));
                    d2 = std::fmax(d2, std::fabs(w.b->derivative_of(i, t, 2)));
                    d3 = std::fmax(d3, std::fabs(w.b->derivative_of(i, t, 3)));
                }
                f << w.name << "," << i << "," << w.b->support_lower(i) << ","
                  << w.b->support_upper(i) << ","
                  << design::to_string(w.c->reason[static_cast<std::size_t>(i)])
                  << "," << v << "," << d1 << "," << d2 << "," << d3 << "\n";
            }
        // F0's modes, from the validation family.
        const design::RadialValidationConfig f0cfg;
        const design::RadialFamily f0 = design::make_f0_family(f0cfg);
        for (int i = 0; i < f0.basis.size(); ++i) {
            double v = 0, d1 = 0, d2 = 0, d3 = 0;
            for (int k = 0; k <= 400; ++k) {
                const double t = f0cfg.s_max * k / 400.0;
                v = std::fmax(v, std::fabs(f0.basis.value_of(i, t)));
                d1 = std::fmax(d1, std::fabs(f0.basis.derivative_of(i, t, 1)));
                d2 = std::fmax(d2, std::fabs(f0.basis.derivative_of(i, t, 2)));
                d3 = std::fmax(d3, std::fabs(f0.basis.derivative_of(i, t, 3)));
            }
            f << "f0_radial," << i << "," << f0.basis.support_lower(i) << ","
              << f0.basis.support_upper(i) << ","
              << design::to_string(f0.reason[static_cast<std::size_t>(i)]) << ","
              << v << "," << d1 << "," << d2 << "," << d3 << "\n";
        }
    }

    // --- Raw perturbation, projected -----------------------------------------
    //
    // One frame per raw coefficient direction: the raw unit perturbation, the
    // component the constraints REMOVE, and what survives. Production performs
    // the projection; the renderer never multiplies P by anything.
    {
        std::ofstream f(out_dir + "/raw_to_constrained_frames.csv");
        f.precision(12);
        f << "frame,raw_index,coefficient_index,raw_value,projected_value,"
             "removed_value,constraint_residual_raw,constraint_residual_projected,"
             "removed_norm,projected_norm\n";
        for (int k = 0; k < basis.coefficient_count(); ++k) {
            std::vector<double> dc(static_cast<std::size_t>(pr.n_raw), 0.0);
            dc[static_cast<std::size_t>(k)] = 1.0;
            // The projection onto the constrained subspace: P (P^T dc).
            std::vector<double> zt(static_cast<std::size_t>(pr.n_free), 0.0);
            for (int j = 0; j < pr.n_free; ++j) {
                double acc = 0.0;
                for (int i = 0; i < pr.n_raw; ++i) acc += pr.P(i, j) * dc[
                    static_cast<std::size_t>(i)];
                zt[static_cast<std::size_t>(j)] = acc;
            }
            std::vector<double> proj(static_cast<std::size_t>(pr.n_raw), 0.0);
            for (int i = 0; i < pr.n_raw; ++i) {
                double acc = 0.0;
                for (int j = 0; j < pr.n_free; ++j)
                    acc += pr.P(i, j) * zt[static_cast<std::size_t>(j)];
                proj[static_cast<std::size_t>(i)] = acc;
            }
            const std::vector<double> r_raw = design::matvec(C, dc);
            const std::vector<double> r_proj = design::matvec(C, proj);
            double rr = 0.0, rp = 0.0, rn = 0.0, pn = 0.0;
            for (double v : r_raw) rr = std::fmax(rr, std::fabs(v));
            for (double v : r_proj) rp = std::fmax(rp, std::fabs(v));
            for (int i = 0; i < pr.n_raw; ++i) {
                const double rem = dc[static_cast<std::size_t>(i)] -
                                   proj[static_cast<std::size_t>(i)];
                rn += rem * rem;
                pn += proj[static_cast<std::size_t>(i)] *
                      proj[static_cast<std::size_t>(i)];
            }
            rn = std::sqrt(rn);
            pn = std::sqrt(pn);
            for (int i = 0; i < pr.n_raw; ++i)
                f << k << "," << k << "," << i << ","
                  << dc[static_cast<std::size_t>(i)] << ","
                  << proj[static_cast<std::size_t>(i)] << ","
                  << (dc[static_cast<std::size_t>(i)] -
                      proj[static_cast<std::size_t>(i)])
                  << "," << rr << "," << rp << "," << rn << "," << pn << "\n";
        }
    }

    // --- The same perturbation as a FIELD, before and after projection ------
    {
        std::ofstream f(out_dir + "/constraint_effect_frames.csv");
        f.precision(10);
        f << "frame,xi,s,unconstrained,projected,difference,"
             "constraint_residual_raw,constraint_residual_projected,"
             "in_outer_collar,in_protected_cabin\n";
        const design::DesignRegion region = canonical_region();
        for (int t = 0; t < 12; ++t) {
            // A deterministic raw perturbation, different each frame.
            std::vector<double> dc(static_cast<std::size_t>(pr.n_raw));
            for (int i = 0; i < pr.n_raw; ++i)
                dc[static_cast<std::size_t>(i)] =
                    std::sin(0.7 * static_cast<double>(i) +
                             0.53 * static_cast<double>(t));
            std::vector<double> zt(static_cast<std::size_t>(pr.n_free), 0.0);
            for (int j = 0; j < pr.n_free; ++j) {
                double acc = 0.0;
                for (int i = 0; i < pr.n_raw; ++i)
                    acc += pr.P(i, j) * dc[static_cast<std::size_t>(i)];
                zt[static_cast<std::size_t>(j)] = acc;
            }
            std::vector<double> proj(static_cast<std::size_t>(pr.n_raw), 0.0);
            for (int i = 0; i < pr.n_raw; ++i) {
                double acc = 0.0;
                for (int j = 0; j < pr.n_free; ++j)
                    acc += pr.P(i, j) * zt[static_cast<std::size_t>(j)];
                proj[static_cast<std::size_t>(i)] = acc;
            }
            double rr = 0.0, rp = 0.0;
            for (double v : design::matvec(C, dc)) rr = std::fmax(rr, std::fabs(v));
            for (double v : design::matvec(C, proj)) rp = std::fmax(rp, std::fabs(v));
            for (int i = 0; i < slice.n_xi; ++i)
                for (int j = 0; j < slice.n_s; ++j) {
                    const DesignCoordinates c = design::from_design_coordinates(
                        slice.xi(i), slice.s(j), kO0aR0);
                    const double a = basis.field_value<double>(dc, c);
                    const double b = basis.field_value<double>(proj, c);
                    const bool collar =
                        slice.xi(i) < region.xi_flat_lower() ||
                        slice.xi(i) > region.xi_flat_upper() ||
                        slice.s(j) > region.s_flat_outer();
                    f << t << "," << slice.xi(i) << "," << slice.s(j) << ","
                      << a << "," << b << "," << (a - b) << "," << rr << ","
                      << rp << "," << (collar ? 1 : 0) << ","
                      << (slice.s(j) < region.protected_cabin_radius ? 1 : 0)
                      << "\n";
                }
        }
    }

    // --- F0 reconstruction, one coefficient at a time ------------------------
    {
        const design::RadialValidationConfig f0cfg;
        const design::RadialFamily f0 = design::make_f0_family(f0cfg);
        const design::RadialProjection proj =
            design::project_shape_function(f0);
        std::ofstream f(out_dir + "/f0_reconstruction_frames.csv");
        f.precision(10);
        f << "frame,coefficients_included,s_r,radius,f_reference,f_spline,"
             "f_first_reference,f_first_spline,f_second_reference,"
             "f_second_spline,zone\n";
        for (std::size_t n = 0; n <= f0.active.size(); ++n) {
            std::vector<double> partial(proj.coefficients.size(), 0.0);
            for (std::size_t m = 0; m < n; ++m)
                partial[static_cast<std::size_t>(f0.active[m])] =
                    proj.coefficients[static_cast<std::size_t>(f0.active[m])];
            for (int k = 0; k <= 200; ++k) {
                const double s = f0cfg.s_max * k / 200.0;
                const design::RadialEvaluation q =
                    design::evaluate_radial(f0, partial, s);
                if (!q.inside) continue;
                const double r = f0cfg.scale_R0 * std::sqrt(std::fmax(s, 0.0));
                const design::ShapeReference ref = design::shape_reference(
                    r, f0cfg.bubble_radius, f0cfg.wall_sigma);
                f << n << "," << n << "," << s << "," << r << "," << ref.value
                  << "," << q.value << "," << ref.first << "," << q.d_r << ","
                  << ref.second << "," << q.d_rr << ","
                  << design::to_string(design::zone_of(f0cfg, s)) << "\n";
            }
        }
    }

    // --- One whitened coordinate at a time -----------------------------------
    //
    // c = c0 + P W e_j. The occupancy column is the point: a single whitened
    // coordinate moves MANY local coefficients, which is why whitening must never
    // be described as preserving local support.
    {
        const design::DenseMatrix PW = design::multiply(pr.P, wr.w);
        std::ofstream f(out_dir + "/whitened_perturbation_frames.csv");
        f.precision(10);
        f << "frame,whitened_index,coefficient_index,coefficient_delta,"
             "coefficients_touched,total_coefficients,max_abs_delta\n";
        std::ofstream g(out_dir + "/whitened_field_frames.csv");
        g.precision(10);
        g << "frame,whitened_index,xi,s,field,in_outer_collar\n";
        const design::DesignRegion region = canonical_region();
        const int nframes = wr.retained < 12 ? wr.retained : 12;
        for (int t = 0; t < nframes; ++t) {
            const int j = t * (wr.retained / (nframes > 0 ? nframes : 1));
            std::vector<double> dc(static_cast<std::size_t>(pr.n_raw), 0.0);
            double mx = 0.0;
            for (int i = 0; i < pr.n_raw; ++i) {
                dc[static_cast<std::size_t>(i)] = PW(i, j);
                mx = std::fmax(mx, std::fabs(PW(i, j)));
            }
            int touched = 0;
            for (int i = 0; i < pr.n_raw; ++i)
                if (std::fabs(dc[static_cast<std::size_t>(i)]) > 1e-6 * mx)
                    ++touched;
            for (int i = 0; i < pr.n_raw; ++i)
                f << t << "," << j << "," << i << ","
                  << dc[static_cast<std::size_t>(i)] << "," << touched << ","
                  << pr.n_raw << "," << mx << "\n";
            std::vector<double> c_full(static_cast<std::size_t>(pr.n_raw));
            for (int i = 0; i < pr.n_raw; ++i)
                c_full[static_cast<std::size_t>(i)] =
                    pr.c0[static_cast<std::size_t>(i)] +
                    dc[static_cast<std::size_t>(i)];
            for (int a = 0; a < slice.n_xi; ++a)
                for (int b = 0; b < slice.n_s; ++b) {
                    const DesignCoordinates c = design::from_design_coordinates(
                        slice.xi(a), slice.s(b), kO0aR0);
                    const bool collar =
                        slice.xi(a) < region.xi_flat_lower() ||
                        slice.xi(a) > region.xi_flat_upper() ||
                        slice.s(b) > region.s_flat_outer();
                    g << t << "," << j << "," << slice.xi(a) << ","
                      << slice.s(b) << ","
                      << basis.field_value<double>(c_full, c) << ","
                      << (collar ? 1 : 0) << "\n";
                }
        }
    }
}

}  // namespace checks
}  // namespace warp_hypad
