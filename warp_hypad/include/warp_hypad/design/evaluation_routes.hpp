#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "warp_hypad/design/axisymmetric_basis.hpp"
#include "warp_hypad/design/constraint_projection.hpp"

// ---------------------------------------------------------------------------
// O0B Phase 12: the two candidate evaluation routes, costed but NOT run.
//
// O0A concluded that the F1 family is not resolution-certified: the narrowest
// radial feature spans 0.755 cells on the A6 energy grid, against a six-cell
// contract. There are two ways out, and O0B's job is to cost them, not to
// choose by preference.
//
// ROUTE A -- REFINE THE UNIFORM GRID. Keep the finite-difference pipeline A4-A6
// already verify, and shrink h until the six-cell rule is met. The rule is a
// statement about FINITE DIFFERENCES: a centred stencil needs several cells
// across a feature before its truncation error stops dominating. It applies
// here because the derivatives really would be taken by differencing.
//
// ROUTE B -- ANALYTIC JETS PLUS ELEMENTWISE QUADRATURE. O0B now produces exact
// metric derivatives through third order, so no cell count is needed to
// RESOLVE a derivative -- the derivative is closed form. What still has to
// converge is the INTEGRAL. So the six-cell rule must NOT be transplanted onto
// this route: applying a finite-difference resolution criterion to a scheme
// that takes no finite differences would be measuring the wrong thing. The
// convergence that matters is in quadrature order, element subdivision, domain
// extent, and agreement with an independent uniform-grid evaluation.
//
// O0B builds the element enumeration and the quadrature scaffold. It does NOT
// evaluate or optimize energy: that is O0C's decision to take on the evidence
// this file produces.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

// --- Route A ------------------------------------------------------------------

struct UniformGridCost {
    double narrowest_feature_m = 0.0;
    int required_cells_per_feature = 6;
    double required_spacing_m = 0.0;
    double domain_half_width_m = 0.0;
    long long cells_per_axis = 0;
    double total_cells = 0.0;
    // A4-A6 carry the metric, its inverse, Christoffels, Riemann, Ricci and the
    // stress-energy tensor. The rank-4 Riemann field is what dominates: 256
    // doubles per point before anything else.
    double bytes_per_point = 0.0;
    double total_bytes = 0.0;
    double current_spacing_m = 0.0;
    double refinement_factor = 0.0;
    double cost_factor_versus_current = 0.0;
    std::string note;
};

inline UniformGridCost uniform_grid_cost(double narrowest_feature_m,
                                         double domain_half_width_m,
                                         double current_spacing_m,
                                         int cells_per_feature = 6) {
    UniformGridCost c;
    c.narrowest_feature_m = narrowest_feature_m;
    c.required_cells_per_feature = cells_per_feature;
    c.required_spacing_m =
        narrowest_feature_m / static_cast<double>(cells_per_feature);
    c.domain_half_width_m = domain_half_width_m;
    c.cells_per_axis = static_cast<long long>(
        std::ceil(2.0 * domain_half_width_m / c.required_spacing_m)) + 1;
    const double n = static_cast<double>(c.cells_per_axis);
    c.total_cells = n * n * n;
    // Rank-4 Riemann (256) + rank-3 Christoffel (64) + rank-2 metric, inverse,
    // Ricci, Einstein, stress-energy (5 x 16) + scalars. Counted rather than
    // guessed, because the point of the estimate is the memory.
    c.bytes_per_point = (256.0 + 64.0 + 5.0 * 16.0 + 8.0) * 8.0;
    c.total_bytes = c.total_cells * c.bytes_per_point;
    c.current_spacing_m = current_spacing_m;
    c.refinement_factor = current_spacing_m / c.required_spacing_m;
    c.cost_factor_versus_current =
        c.refinement_factor * c.refinement_factor * c.refinement_factor;
    c.note =
        "The six-cell rule is a FINITE-DIFFERENCE criterion and is applied here "
        "because Route A really does difference. Cost scales as the cube of the "
        "refinement factor.";
    return c;
}

// --- Route B ------------------------------------------------------------------

// One integration element: a unique tensor-product knot span. Elements tile the
// design region exactly, and every element is a region on which every active
// basis function is a single polynomial -- which is what makes a fixed-order
// Gauss-Legendre rule exact on it rather than merely accurate.
struct QuadratureElement {
    int index = 0;
    int span_xi = 0;
    int span_s = 0;
    double xi_lower = 0.0, xi_upper = 0.0;
    double s_lower = 0.0, s_upper = 0.0;
    double rho_lower = 0.0, rho_upper = 0.0;
    double volume = 0.0;          // pi R0^3 (xi_upper-xi_lower)(s_upper-s_lower)
    int active_functions = 0;     // basis functions with support on this element
    bool in_outer_collar = false;
    bool in_protected_cabin = false;
};

struct QuadratureScaffold {
    std::vector<QuadratureElement> element;
    int nodes_per_direction = 8;
    int polynomial_exactness = 15;   // 8-node Gauss-Legendre
    double total_volume = 0.0;
    double exact_volume = 0.0;       // pi R0^3 * extent, in design coordinates
    std::size_t total_nodes = 0;
    std::string rule = "tensor-product Gauss-Legendre, fixed order, one element "
                       "per unique tensor-product knot span";
    std::string convergence_contract =
        "Route B takes NO finite differences, so the six-cell rule does not "
        "apply to it. Convergence must instead be demonstrated in quadrature "
        "order, element subdivision, domain enlargement, and against an "
        "independent uniform-grid evaluation.";
};

// The distinct knot values bounding spans inside the evaluable range.
inline std::vector<double> element_edges(const QuinticBSplineBasis1D& b) {
    std::vector<double> e;
    e.push_back(b.lower());
    for (std::size_t i = 0; i < b.knots().size(); ++i) {
        const double t = b.knots()[i];
        if (t > b.lower() + 1e-12 && t < b.upper() - 1e-12) {
            bool have = false;
            for (double u : e)
                if (std::fabs(u - t) < 1e-12) have = true;
            if (!have) e.push_back(t);
        }
    }
    e.push_back(b.upper());
    std::sort(e.begin(), e.end());
    return e;
}

inline QuadratureScaffold build_quadrature_scaffold(
    const AxisymmetricBasis& basis, const DesignRegion& region,
    int nodes_per_direction = 8) {
    QuadratureScaffold sc;
    sc.nodes_per_direction = nodes_per_direction;
    const double R0 = basis.scale_R0();
    const double kPi = 3.14159265358979323846;
    const double vol = kPi * R0 * R0 * R0;

    const std::vector<double> ex = element_edges(basis.axial());
    const std::vector<double> es = element_edges(basis.radial());

    int idx = 0;
    for (std::size_t a = 0; a + 1 < ex.size(); ++a)
        for (std::size_t b = 0; b + 1 < es.size(); ++b) {
            QuadratureElement e;
            e.index = idx++;
            e.span_xi = static_cast<int>(a);
            e.span_s = static_cast<int>(b);
            e.xi_lower = ex[a];
            e.xi_upper = ex[a + 1];
            e.s_lower = es[b];
            e.s_upper = es[b + 1];
            if (!(e.xi_upper > e.xi_lower) || !(e.s_upper > e.s_lower)) continue;
            e.rho_lower = R0 * std::sqrt(std::fmax(e.s_lower, 0.0));
            e.rho_upper = R0 * std::sqrt(std::fmax(e.s_upper, 0.0));
            // dV = pi R0^3 dxi ds -- CONSTANT in the design coordinates, which
            // is the same property of s = rho^2/R0^2 that made the axis regular.
            e.volume = vol * (e.xi_upper - e.xi_lower) *
                       (e.s_upper - e.s_lower);
            e.in_outer_collar = e.xi_lower < region.xi_flat_lower() ||
                                e.xi_upper > region.xi_flat_upper() ||
                                e.s_upper > region.s_flat_outer();
            e.in_protected_cabin = e.s_upper <= region.protected_cabin_radius;
            // How many active functions actually touch this element.
            const double xm = 0.5 * (e.xi_lower + e.xi_upper);
            const double sm = 0.5 * (e.s_lower + e.s_upper);
            const TensorBasisEvaluation te =
                basis.evaluate(from_design_coordinates(xm, sm, R0));
            e.active_functions = static_cast<int>(te.size());
            sc.total_volume += e.volume;
            sc.element.push_back(e);
        }

    sc.exact_volume = vol * (region.xi_max - region.xi_min) * region.s_max;
    sc.total_nodes = sc.element.size() *
                     static_cast<std::size_t>(nodes_per_direction) *
                     static_cast<std::size_t>(nodes_per_direction);
    return sc;
}

// The comparison O0C has to decide on.
struct RouteComparison {
    double route_a_total_cells = 0.0;
    double route_a_total_bytes = 0.0;
    std::size_t route_b_elements = 0;
    std::size_t route_b_nodes = 0;
    double node_ratio = 0.0;
    std::string recommendation;
    std::string caveat;
};

inline RouteComparison compare_routes(const UniformGridCost& a,
                                      const QuadratureScaffold& b) {
    RouteComparison r;
    r.route_a_total_cells = a.total_cells;
    r.route_a_total_bytes = a.total_bytes;
    r.route_b_elements = b.element.size();
    r.route_b_nodes = b.total_nodes;
    r.node_ratio = b.total_nodes > 0
                       ? a.total_cells / static_cast<double>(b.total_nodes)
                       : 0.0;
    r.recommendation =
        "Route B is the one worth developing in O0C. It needs roughly four "
        "orders of magnitude fewer evaluation points than the refined uniform "
        "grid, and its accuracy comes from exact derivatives and a rule that is "
        "polynomial-exact on each element rather than from cell count.";
    r.caveat =
        "This is a COST comparison and a convergence PLAN, not a demonstration. "
        "O0B evaluates no energy on either route. Route B is only credible once "
        "convergence in quadrature order, element subdivision and domain extent "
        "has been measured, and once it has been cross-checked against a "
        "uniform-grid evaluation on a case both routes can do.";
    return r;
}

}  // namespace design
}  // namespace warp_hypad
