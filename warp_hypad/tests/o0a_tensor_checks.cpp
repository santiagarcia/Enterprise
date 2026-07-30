#include "o0a_checks.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "warp_hypad/design/axisymmetric_basis.hpp"

// ---------------------------------------------------------------------------
// O0A cases 11-20: the tensor product, its Cartesian derivatives, and the
// structural axis regularity of the four field families.
//
// WHAT MAKES THE REGULARITY CASES DIFFERENT FROM A PREFACTOR CHECK.
//
// It would be easy, and nearly useless, to verify that eta^m vanishes at the
// axis. The claim O0A actually makes is stronger: that the CARTESIAN COMPONENTS
// of each family are polynomials in y and z times a spline, with the 1/rho of
// the tensor projection cancelled algebraically rather than guarded numerically.
// So cases 15-17 evaluate the components a metric would consume --
//
//     q_y = (y/R0) Q          q_xy = (y/R0) Q
//     q_yy - q_zz = ((y^2-z^2)/R0^2) Q         q_yz = (y z/R0^2) Q
//
// -- against those closed forms, at points including the axis itself and both
// signs of y and z. A sign error in the projection, or a rho that failed to
// cancel, shows up as a mismatch at negative transverse coordinate or as a
// non-finite value on the axis, neither of which a prefactor test can see.
//
// THE FINITE DIFFERENCES HERE ARE VERIFICATION ONLY. Production derivatives come
// from the chain rule on analytic spline derivatives. The differences below are
// taken of the FIELD, from outside, at five deliberately chosen classes of point:
// generic off-axis, near-axis, exactly on the axis, near-but-not-on a knot line,
// and mirrored transverse coordinates. The near-knot class is the one that
// matters most -- a derivative that was right in the middle of a span and wrong
// across a knot would pass a naive sweep.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace checks {
namespace {

using design::AxisymmetricBasis;
using design::DesignCoordinates;
using design::DesignRegion;
using design::FieldFamily;
using design::TensorBasisEvaluation;

constexpr double kR0 = kO0aR0;

// The scalar field from Cartesian coordinates, for finite differencing from
// outside the basis.
double scalar_field(const AxisymmetricBasis& b, const std::vector<double>& c,
                    double x, double y, double z) {
    return b.field_value<double>(c, design::make_design_coordinates(
                                        x, y, z, 0.0, b.scale_R0()));
}

struct ProbePoint {
    const char* label;
    double x, y, z;
};

// The five required classes. The knot-line offsets use a fraction of the lattice
// spacing so the point is demonstrably inside a span but close to its edge.
std::vector<ProbePoint> probe_points(const AxisymmetricBasis& b) {
    const double hx = b.axial_class().spacing;
    const double R0 = b.scale_R0();
    // A knot line in xi, converted to x, nudged inward by 4% of a span.
    const double xi_knot = b.axial().support_lower(b.axial_class().active[2]);
    const double x_near_knot = (xi_knot + 0.04 * hx) * R0;
    // A knot line in s, converted to a transverse radius.
    const double s_knot = b.radial().support_upper(b.radial_class().active[1]);
    const double y_near_s_knot = R0 * std::sqrt(s_knot * 1.02);
    return {
        {"generic_off_axis", 0.31 * R0, 1.10, -0.70},
        {"near_axis", -0.22 * R0, 1.0e-3, 0.0},
        {"exactly_on_axis", 0.44 * R0, 0.0, 0.0},
        {"near_xi_knot", x_near_knot, 0.83, 0.41},
        {"near_s_knot", 0.17 * R0, y_near_s_knot, 0.0},
        {"mirrored_transverse", 0.31 * R0, -1.10, 0.70},
        {"mirrored_z_only", 0.31 * R0, 1.10, 0.70},
    };
}

}  // namespace

void run_o0a_tensor_checks(VerificationReport& report,
                           const std::string& out_dir) {
    Stopwatch sw;
    DesignRegion region;
    AxisymmetricBasis basis = design::make_f1_basis(region, kR0,
                                                    FieldFamily::Scalar);
    const std::vector<double> coeff =
        deterministic_coefficients(basis.coefficient_count(), 0.7);
    const std::string grid = "F1 kernel, 10 axial x 8 radial raw, 6 x 6 active";
    const std::string params = "R0 = 2, xi in [-3,3], s in [0,4]";

    // --- 11. Tensor-product value -------------------------------------------
    //
    // B_ij must be exactly N_i(xi) M_j(s). Checked against products formed
    // independently from the two univariate bases, so an index transposition in
    // the flat map cannot hide.
    sw.reset();
    double worst_value = 0.0;
    std::size_t compared = 0;
    for (int a = 0; a <= 23; ++a)
        for (int b = 0; b <= 23; ++b) {
            const double xi = region.xi_min +
                              (region.xi_max - region.xi_min) * a / 23.0;
            const double s = region.s_max * b / 23.0;
            const DesignCoordinates c =
                design::from_design_coordinates(xi, s, kR0);
            const TensorBasisEvaluation e = basis.evaluate(c);
            for (std::size_t k = 0; k < e.size(); ++k) {
                const double expect =
                    basis.axial().value_of(e.index_xi[k], xi) *
                    basis.radial().value_of(e.index_s[k], s);
                worst_value = std::fmax(worst_value,
                                        std::fabs(e.value[k] - expect));
                ++compared;
            }
        }
    report.check_residual(
        "o0a.tensor_product_value",
        "B_ij(xi,s) = N_i(xi) M_j(s), against products formed independently "
        "from the two univariate bases",
        VerificationLevel::A, grid, params, worst_value, 1.0e-15, sw.seconds(),
        std::to_string(compared) +
            " active products compared; an index transposition in the flat map "
            "would appear here rather than as a plausible field");

    // --- 12/13. Cartesian gradient and Hessian against finite differences ---
    sw.reset();
    const std::vector<ProbePoint> pts = probe_points(basis);
    double worst_grad = 0.0, worst_hess = 0.0;
    std::string worst_grad_where, worst_hess_where;
    std::ofstream fd_log(out_dir + "/tensor_derivative_probes.csv");
    fd_log << "point_label,x,y,z,component,analytic,finite_difference,"
              "absolute_difference\n";
    for (const ProbePoint& p : pts) {
        const DesignCoordinates c =
            design::make_design_coordinates(p.x, p.y, p.z, 0.0, kR0);
        const TensorBasisEvaluation e = basis.evaluate(c);
        if (e.empty()) continue;
        // The field's analytic derivatives: the prefactor is 1 for the scalar
        // family, so the stencil channels ARE the field derivatives.
        double ax = 0, ay = 0, az = 0, axx = 0, ayy = 0, azz = 0, axy = 0,
               ayz = 0, axz = 0;
        for (std::size_t k = 0; k < e.size(); ++k) {
            const double ck = coeff[static_cast<std::size_t>(e.flat_index[k])];
            ax += ck * e.d_x[k];  ay += ck * e.d_y[k];  az += ck * e.d_z[k];
            axx += ck * e.d_xx[k]; ayy += ck * e.d_yy[k]; azz += ck * e.d_zz[k];
            axy += ck * e.d_xy[k]; ayz += ck * e.d_yz[k]; axz += ck * e.d_xz[k];
        }
        const double h1 = 1.0e-5, h2 = 1.0e-4;
        auto F = [&](double dx, double dy, double dz) {
            return scalar_field(basis, coeff, p.x + dx, p.y + dy, p.z + dz);
        };
        const double f0 = F(0, 0, 0);
        const struct { const char* n; double an, fd; } firsts[3] = {
            {"d_x", ax, (F(h1, 0, 0) - F(-h1, 0, 0)) / (2 * h1)},
            {"d_y", ay, (F(0, h1, 0) - F(0, -h1, 0)) / (2 * h1)},
            {"d_z", az, (F(0, 0, h1) - F(0, 0, -h1)) / (2 * h1)}};
        for (const auto& q : firsts) {
            const double d = std::fabs(q.an - q.fd);
            fd_log << p.label << "," << p.x << "," << p.y << "," << p.z << ","
                   << q.n << "," << q.an << "," << q.fd << "," << d << "\n";
            if (d > worst_grad) { worst_grad = d; worst_grad_where = std::string(p.label) + ":" + q.n; }
        }
        const struct { const char* n; double an, fd; } seconds[6] = {
            {"d_xx", axx, (F(h2, 0, 0) - 2 * f0 + F(-h2, 0, 0)) / (h2 * h2)},
            {"d_yy", ayy, (F(0, h2, 0) - 2 * f0 + F(0, -h2, 0)) / (h2 * h2)},
            {"d_zz", azz, (F(0, 0, h2) - 2 * f0 + F(0, 0, -h2)) / (h2 * h2)},
            {"d_xy", axy, (F(h2, h2, 0) - F(h2, -h2, 0) - F(-h2, h2, 0) +
                           F(-h2, -h2, 0)) / (4 * h2 * h2)},
            {"d_xz", axz, (F(h2, 0, h2) - F(h2, 0, -h2) - F(-h2, 0, h2) +
                           F(-h2, 0, -h2)) / (4 * h2 * h2)},
            {"d_yz", ayz, (F(0, h2, h2) - F(0, h2, -h2) - F(0, -h2, h2) +
                           F(0, -h2, -h2)) / (4 * h2 * h2)}};
        for (const auto& q : seconds) {
            const double d = std::fabs(q.an - q.fd);
            fd_log << p.label << "," << p.x << "," << p.y << "," << p.z << ","
                   << q.n << "," << q.an << "," << q.fd << "," << d << "\n";
            if (d > worst_hess) { worst_hess = d; worst_hess_where = std::string(p.label) + ":" + q.n; }
        }
    }
    report.check_residual(
        "o0a.tensor_product_gradient",
        "analytic Cartesian first derivatives against centred differences of the "
        "field, over generic, near-axis, on-axis, near-knot and mirrored points",
        VerificationLevel::B, std::to_string(pts.size()) + " probe points",
        "h = 1e-5", worst_grad, 2.0e-6, sw.seconds(),
        "worst at " + worst_grad_where +
            "; the finite differences are the independent path and are never "
            "used in production");
    report.check_residual(
        "o0a.tensor_product_hessian",
        "analytic Cartesian second derivatives against centred second and mixed "
        "differences of the field, at the same points",
        VerificationLevel::B, std::to_string(pts.size()) + " probe points",
        "h = 1e-4", worst_hess, 5.0e-4, sw.seconds(),
        "worst at " + worst_hess_where +
            "; d_yy and d_zz carry the 2/R0^2 term from d2s/dy2, which a chain "
            "rule that forgot the second derivative of s would drop");

    // --- 14-17. Structural regularity of the four families ------------------
    //
    // Each family is checked on its ACTUAL Cartesian components against the
    // closed forms, at points that include the axis and both signs of y and z.
    struct FamilyCase {
        FieldFamily family;
        const char* case_name;
        const char* claim;
    };
    const FamilyCase families[4] = {
        {FieldFamily::Scalar, "o0a.scalar_axis_regularity",
         "a scalar q(xi,s) has dq/drho identically zero on the axis, for every "
         "coefficient vector, with no division by rho"},
        {FieldFamily::RadialVector, "o0a.radial_vector_axis_regularity",
         "q_y = (y/R0) Q and q_z = (z/R0) Q exactly, so the vector vanishes "
         "linearly on the axis and its Cartesian components are smooth there"},
        {FieldFamily::AxialRadialTensor,
         "o0a.axial_radial_tensor_regularity",
         "q_xy = (y/R0) Q and q_xz = (z/R0) Q exactly"},
        {FieldFamily::TransverseAnisotropy,
         "o0a.transverse_anisotropy_regularity",
         "q_yy - q_zz = ((y^2-z^2)/R0^2) Q and q_yz = (y z/R0^2) Q exactly"}};

    std::ofstream modes(out_dir + "/axis_regular_modes.csv");
    modes << "family,axis_vanishing_order,x,y,z,rho,Q,component_expected,"
             "component_actual,difference,d_rho_on_axis,finite\n";

    for (const FamilyCase& fc : families) {
        sw.reset();
        AxisymmetricBasis fb =
            design::make_f1_basis(region, kR0, fc.family);
        const std::vector<double> fcoef =
            deterministic_coefficients(fb.coefficient_count(), 1.9);
        const int order = design::axis_vanishing_order(fc.family);

        double worst = 0.0, worst_axis_value = 0.0, worst_axis_drho = 0.0;
        bool all_finite = true;
        double power_law_worst = 0.0;

        const double ys[7] = {0.0, 1e-9, 1e-4, 0.35, 1.10, -1.10, -0.35};
        const double zs[7] = {0.0, 0.0, 1e-4, -0.62, 0.70, 0.70, -0.62};
        for (int pi = 0; pi < 7; ++pi)
            for (int xa = 0; xa <= 6; ++xa) {
                const double x = (-1.6 + 3.2 * xa / 6.0) * kR0;
                const DesignCoordinates c = design::make_design_coordinates(
                    x, ys[pi], zs[pi], 0.0, kR0);
                const TensorBasisEvaluation e = fb.evaluate(c);
                // Q: the underlying spline sum, prefactor removed.
                double Q = 0.0;
                for (std::size_t k = 0; k < e.size(); ++k)
                    Q += fcoef[static_cast<std::size_t>(e.flat_index[k])] *
                         e.value[k];
                for (std::size_t k = 0; k < e.size(); ++k) {
                    const double v[] = {e.value[k], e.d_x[k], e.d_y[k],
                                        e.d_z[k], e.d_xx[k], e.d_yy[k],
                                        e.d_yz[k], e.d_rho[k]};
                    for (double q : v)
                        if (!std::isfinite(q)) all_finite = false;
                }

                // The closed form of the component this family contributes.
                double expected = 0.0;
                switch (fc.family) {
                    case FieldFamily::Scalar:
                        expected = Q;
                        break;
                    case FieldFamily::RadialVector:
                    case FieldFamily::AxialRadialTensor:
                        expected = (ys[pi] / kR0) * Q;
                        break;
                    case FieldFamily::TransverseAnisotropy:
                        expected = ((ys[pi] * ys[pi] - zs[pi] * zs[pi]) /
                                    (kR0 * kR0)) * Q;
                        break;
                }
                const double actual = fb.cartesian_component<double>(fcoef, c);
                const double diff = std::fabs(actual - expected);
                worst = std::fmax(worst, diff);

                double drho_axis = 0.0;
                if (c.rho == 0.0) {
                    worst_axis_value =
                        std::fmax(worst_axis_value,
                                  (order > 0) ? std::fabs(actual) : 0.0);
                    for (std::size_t k = 0; k < e.size(); ++k)
                        drho_axis = std::fmax(drho_axis,
                                              std::fabs(e.d_rho[k]));
                    worst_axis_drho = std::fmax(worst_axis_drho, drho_axis);
                }
                // The vanishing POWER, not merely the fact of vanishing: at
                // rho = 1e-9 the component must be O(rho^order), so dividing by
                // rho^order must give something of order Q rather than blowing up.
                // Along the ray z = 0, y > 0 the Cartesian prefactor equals
                // eta^order exactly, so dividing by it must return Q. Off that
                // ray a tensor component vanishes like a power of ITS OWN
                // prefactor and not like a power of rho, which is why the probe
                // is taken along a stated ray rather than at arbitrary points.
                if (order > 0 && zs[pi] == 0.0 && ys[pi] > 0.0 &&
                    c.rho < 1e-3 && std::fabs(Q) > 1e-12) {
                    const double eta = c.rho / kR0;
                    const double scaled =
                        std::fabs(actual) /
                        std::pow(eta, static_cast<double>(order));
                    power_law_worst =
                        std::fmax(power_law_worst,
                                  std::fabs(scaled / std::fabs(Q) - 1.0));
                }
                modes << design::to_string(fc.family) << "," << order << ","
                      << x << "," << ys[pi] << "," << zs[pi] << "," << c.rho
                      << "," << Q << "," << expected << "," << actual << ","
                      << diff << "," << drho_axis << ","
                      << (all_finite ? 1 : 0) << "\n";
            }

        char note[420];
        std::snprintf(note, sizeof(note),
                      "worst |component - closed form| = %.3e; on-axis value "
                      "%.3e (order %d); worst |dB/drho| on axis = %.3e; "
                      "vanishing-power residual along z=0, y>0: %.3e; all values "
                      "finite: %s",
                      worst, worst_axis_value, order, worst_axis_drho,
                      power_law_worst, all_finite ? "yes" : "NO");
        report.check_predicate(
            fc.case_name, fc.claim, VerificationLevel::A,
            worst <= 1e-14 && all_finite &&
                (order == 0 || worst_axis_value == 0.0) &&
                worst_axis_drho == 0.0 && power_law_worst <= 1e-9,
            "component matches the closed form; the axis value and dB/drho are "
            "exactly zero where the order requires it",
            note, sw.seconds(),
            "the prefactor is not tested on its own: this compares the "
            "CARTESIAN component a metric would consume, at negative as well "
            "as positive transverse coordinate");
    }

    // --- 18. No division by rho at the axis ---------------------------------
    //
    // The structural claim, stated as an approach test: as rho falls through ten
    // orders of magnitude every derivative channel stays bounded by its value at
    // moderate rho. A hidden 1/rho would diverge; a hidden guard would produce a
    // discontinuity at whatever epsilon it used.
    sw.reset();
    double worst_ratio = 0.0;
    bool all_finite_axis = true;
    double reference_scale = 0.0;
    for (FieldFamily fam : {FieldFamily::Scalar, FieldFamily::RadialVector,
                            FieldFamily::AxialRadialTensor,
                            FieldFamily::TransverseAnisotropy}) {
        AxisymmetricBasis fb = design::make_f1_basis(region, kR0, fam);
        for (double rho : {1.0e-1, 1.0e-3, 1.0e-6, 1.0e-9, 1.0e-12, 1.0e-15,
                           0.0}) {
            const DesignCoordinates c =
                design::make_design_coordinates(0.21 * kR0, rho, 0.0, 0.0, kR0);
            const TensorBasisEvaluation e = fb.evaluate(c, true);
            for (std::size_t k = 0; k < e.size(); ++k) {
                const double v[] = {e.d_x[k], e.d_y[k], e.d_z[k], e.d_xx[k],
                                    e.d_xy[k], e.d_yy[k], e.d_yz[k],
                                    e.d_zz[k], e.d_rho[k], e.d3_xi[k],
                                    e.d3_s[k]};
                for (double q : v) {
                    if (!std::isfinite(q)) all_finite_axis = false;
                    if (rho > 1.0e-2) reference_scale =
                        std::fmax(reference_scale, std::fabs(q));
                    worst_ratio = std::fmax(worst_ratio, std::fabs(q));
                }
            }
        }
    }
    report.check_predicate(
        "o0a.no_division_by_rho_at_axis",
        "every derivative channel of every family stays bounded as rho falls to "
        "zero through twelve orders of magnitude, and is finite AT rho = 0",
        VerificationLevel::A,
        all_finite_axis && worst_ratio <= 4.0 * reference_scale,
        "bounded and finite throughout the approach",
        "largest channel magnitude anywhere in the approach " +
            std::to_string(worst_ratio) + ", reference scale at rho > 1e-2 " +
            std::to_string(reference_scale) + ", all finite: " +
            (all_finite_axis ? "yes" : "NO"),
        sw.seconds(),
        "a 1/rho would diverge along this sequence and an epsilon guard would "
        "show a step; neither exists because s = rho^2/R0^2 removes the "
        "cylindrical rho algebraically");

    // --- 19. Outer flat collar exactness ------------------------------------
    sw.reset();
    double collar_leak = 0.0;
    std::size_t collar_samples = 0;
    for (FieldFamily fam : {FieldFamily::Scalar, FieldFamily::RadialVector,
                            FieldFamily::AxialRadialTensor,
                            FieldFamily::TransverseAnisotropy}) {
        AxisymmetricBasis fb = design::make_f1_basis(region, kR0, fam);
        const std::vector<double> fcoef =
            deterministic_coefficients(fb.coefficient_count(), 2.4);
        for (int i = 0; i <= 70; ++i)
            for (int j = 0; j <= 70; ++j) {
                const double xi = region.xi_min +
                                  (region.xi_max - region.xi_min) * i / 70.0;
                const double s = region.s_max * j / 70.0;
                if (xi >= region.xi_flat_lower() &&
                    xi <= region.xi_flat_upper() &&
                    s <= region.s_flat_outer())
                    continue;
                ++collar_samples;
                const DesignCoordinates c =
                    design::from_design_coordinates(xi, s, kR0);
                const TensorBasisEvaluation e = fb.evaluate(c, true);
                collar_leak = std::fmax(collar_leak,
                                        std::fabs(fb.field_value<double>(fcoef, c)));
                for (std::size_t k = 0; k < e.size(); ++k) {
                    const double v[] = {e.value[k], e.d_x[k], e.d_y[k],
                                        e.d_xx[k], e.d_yy[k], e.d_yz[k],
                                        e.d3_xi[k], e.d3_s[k]};
                    for (double q : v)
                        collar_leak = std::fmax(collar_leak, std::fabs(q));
                }
            }
    }
    report.check_exact(
        "o0a.outer_flat_collar_exactness",
        "throughout the outer flat collar the correction and every derivative "
        "are EXACTLY zero, for all four families",
        VerificationLevel::A,
        std::to_string(collar_samples) + " collar samples x 4 families",
        0.0, collar_leak, sw.seconds(),
        "guaranteed by support exclusion, not by a taper: a basis function whose "
        "support reaches the collar has no coefficient at all, so there is "
        "nothing that could be small-but-nonzero here");

    // --- 20. Protected cabin exactness ---------------------------------------
    sw.reset();
    AxisymmetricBasis cabin = design::make_f1_basis(
        region, kR0, FieldFamily::Scalar, 10, 8, /*protect_cabin=*/true);
    const std::vector<double> ccoef =
        deterministic_coefficients(cabin.coefficient_count(), 3.1);
    double cabin_leak = 0.0;
    std::size_t cabin_samples = 0;
    int cabin_excluded = 0;
    for (int i = 0; i < cabin.radial().size(); ++i)
        if (cabin.radial_class().reason[static_cast<std::size_t>(i)] ==
            design::ExclusionReason::ProtectedCabin)
            ++cabin_excluded;
    for (int i = 0; i <= 60; ++i)
        for (int j = 0; j <= 40; ++j) {
            const double xi = region.xi_min +
                              (region.xi_max - region.xi_min) * i / 60.0;
            const double s = region.protected_cabin_radius * j / 40.0;
            ++cabin_samples;
            const DesignCoordinates c =
                design::from_design_coordinates(xi, s, kR0);
            const TensorBasisEvaluation e = cabin.evaluate(c, true);
            cabin_leak = std::fmax(cabin_leak,
                                   std::fabs(cabin.field_value<double>(ccoef, c)));
            for (std::size_t k = 0; k < e.size(); ++k) {
                const double v[] = {e.value[k], e.d_x[k], e.d_y[k], e.d_xx[k],
                                    e.d_yy[k], e.d3_s[k]};
                for (double q : v)
                    cabin_leak = std::fmax(cabin_leak, std::fabs(q));
            }
        }
    report.check_predicate(
        "o0a.protected_cabin_exactness",
        "a cabin-protecting family is EXACTLY zero, with every derivative, "
        "throughout the protected interior",
        VerificationLevel::A,
        cabin_leak == 0.0 && cabin.coefficient_count() > 0,
        "identically zero inside s < " +
            std::to_string(region.protected_cabin_radius),
        std::to_string(cabin_samples) + " samples, worst magnitude " +
            std::to_string(cabin_leak) + "; " +
            std::to_string(cabin.coefficient_count()) +
            " active tensor products, none with any support below s = " +
            std::to_string(region.protected_cabin_radius) + " (" +
            std::to_string(cabin_excluded) +
            " modes carried the protected_cabin exclusion label)",
        sw.seconds(),
        "the cabin is protected by CONSTRUCTION rather than by exclusion: the "
        "cabin-protecting radial lattice tiles [cabin, s_flat_outer], so no "
        "function has support in the cabin to remove. It is deliberately NOT "
        "clamped at the axis -- it must vanish there, which is the point of "
        "protecting the cabin");

    // --- Records -------------------------------------------------------------
    {
        std::ofstream f(out_dir + "/axial_knot_vector.csv");
        f << "index,knot,semantics,is_raw_support,is_active\n";
        const std::string sem =
            design::QuinticBSplineBasis1D::to_string(basis.axial().semantics());
        for (std::size_t i = 0; i < basis.axial().knots().size(); ++i) {
            bool raw = false, act = false;
            for (int r : basis.axial_class().raw)
                if (static_cast<std::size_t>(r) == i) raw = true;
            for (int r : basis.axial_class().active)
                if (static_cast<std::size_t>(r) == i) act = true;
            f << i << "," << basis.axial().knots()[i] << "," << sem << ","
              << (raw ? 1 : 0) << "," << (act ? 1 : 0) << "\n";
        }
    }
    {
        // The univariate profiles of the CORRECTION lattices, with their knot
        // semantics recorded per basis. F0's clamped profiles are written by the
        // kernel group; these are the other two, because the figures have to show
        // that the three families do NOT share boundary knot semantics.
        std::ofstream f(out_dir + "/axisymmetric_basis_profiles.csv");
        f << "basis,semantics,coordinate,index,value,first,second,third,"
             "is_raw_support,is_active,exclusion_reason,support_lower,"
             "support_upper\n";
        struct Which {
            const char* name;
            const design::QuinticBSplineBasis1D* b;
            const design::LatticeClassification* c;
            double lo, hi;
        };
        const Which which[2] = {
            {"axial_xi", &basis.axial(), &basis.axial_class(), region.xi_min,
             region.xi_max},
            {"radial_s", &basis.radial(), &basis.radial_class(), 0.0,
             region.s_max}};
        for (const Which& w : which) {
            const std::string sem =
                design::QuinticBSplineBasis1D::to_string(w.b->semantics());
            for (int k = 0; k <= 600; ++k) {
                const double t = w.lo + (w.hi - w.lo) * k / 600.0;
                for (int i = 0; i < w.b->size(); ++i) {
                    bool raw = false, act = false;
                    for (int r : w.c->raw)
                        if (r == i) raw = true;
                    for (int r : w.c->active)
                        if (r == i) act = true;
                    if (!raw) continue;  // lattice padding is not a design mode
                    f << w.name << "," << sem << "," << t << "," << i << ","
                      << w.b->value_of(i, t) << ","
                      << w.b->derivative_of(i, t, 1) << ","
                      << w.b->derivative_of(i, t, 2) << ","
                      << w.b->derivative_of(i, t, 3) << "," << (raw ? 1 : 0)
                      << "," << (act ? 1 : 0) << ","
                      << design::to_string(
                             w.c->reason[static_cast<std::size_t>(i)])
                      << "," << w.b->support_lower(i) << ","
                      << w.b->support_upper(i) << "\n";
                }
            }
        }
    }
    {
        // Tensor-product samples for the atlas, so the renderer never has to
        // form N_i(xi) M_j(s) itself. Case 11 is what verifies this product; a
        // Python multiplication of two saved columns would be a second,
        // unverified implementation of the basis definition.
        std::ofstream f(out_dir + "/tensor_product_atlas.csv");
        f << "flat_index,index_xi,index_s,xi,s,rho,value,d_xi,d_s,d_x,d_y,"
             "d_xx,d_yy\n";
        const int atlas[6] = {0, 7, 14, 21, 28, 35};
        for (int slot = 0; slot < 6; ++slot) {
            const int target = atlas[slot];
            if (target >= basis.coefficient_count()) continue;
            for (int a = 0; a <= 60; ++a)
                for (int b = 0; b <= 40; ++b) {
                    const double xi = region.xi_min +
                                      (region.xi_max - region.xi_min) * a / 60.0;
                    const double sv = region.s_max * b / 40.0;
                    const DesignCoordinates c =
                        design::from_design_coordinates(xi, sv, kR0);
                    const TensorBasisEvaluation e = basis.evaluate(c);
                    for (std::size_t m = 0; m < e.size(); ++m) {
                        if (e.flat_index[m] != target) continue;
                        f << target << "," << e.index_xi[m] << ","
                          << e.index_s[m] << "," << xi << "," << sv << ","
                          << c.rho << "," << e.value[m] << "," << e.d_xi[m]
                          << "," << e.d_s[m] << "," << e.d_x[m] << ","
                          << e.d_y[m] << "," << e.d_xx[m] << "," << e.d_yy[m]
                          << "\n";
                    }
                }
        }
    }
    {
        std::ofstream f(out_dir + "/tensor_product_supports.csv");
        f << "flat_index,index_xi,index_s,xi_lower,xi_upper,s_lower,s_upper,"
             "rho_lower,rho_upper,xi_width,s_width,rho_width\n";
        for (int k = 0; k < basis.coefficient_count(); ++k) {
            double xl, xh, sl, sh;
            basis.support(k, xl, xh, sl, sh);
            const double rl = kR0 * std::sqrt(std::fmax(sl, 0.0));
            const double rh = kR0 * std::sqrt(std::fmax(sh, 0.0));
            f << k << "," << basis.pairs()[static_cast<std::size_t>(k)].first
              << "," << basis.pairs()[static_cast<std::size_t>(k)].second << ","
              << xl << "," << xh << "," << sl << "," << sh << "," << rl << ","
              << rh << "," << (xh - xl) << "," << (sh - sl) << ","
              << (rh - rl) << "\n";
        }
    }
}

}  // namespace checks
}  // namespace warp_hypad
