#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "warp_hypad/design/axisymmetric_basis.hpp"
#include "warp_hypad/design/dense_linalg.hpp"

// ---------------------------------------------------------------------------
// Exact linear constraints, the coefficient-coordinate change they induce, and
// the Sobolev whitening that makes optimization coordinates comparable.
//
// WHY CONSTRAINTS ARE ELIMINATED AND NOT PENALISED.
//
// Every constraint here is exactly linear in the coefficients: a fixed value, a
// reflection, a vanishing moment, a normalisation. Handing such a constraint to
// an optimizer as a penalty term does three harmful things at once -- it makes
// the constraint approximately satisfied at every iterate and exactly satisfied
// at none, it adds a stiff direction that dominates the conditioning, and it
// silently trades constraint violation against objective improvement using a
// weight nobody chose on physical grounds. Eliminating them instead gives
//
//     c = c0 + P z
//
// where c0 satisfies the affine system, the columns of P span the null space of
// the homogeneous part, and z is unconstrained. Every iterate then satisfies
// every constraint to machine precision, by construction.
//
// THE MOMENT CONSTRAINT IS THE ONE THAT EARNS ITS KEEP. A shape correction that
// is merely a translation of the bubble costs nothing to find and buys nothing:
// it changes the coefficient vector, leaves the geometry the same up to a shift
// of origin, and would let an optimizer report progress while moving in a
// direction that is pure gauge. Requiring the first xi-moment of the correction
// to vanish removes that direction from the search space outright rather than
// hoping the objective notices.
//
// WHY WHITENING. Coefficient magnitude is not a measure of anything: two modes
// with the same coefficient can differ by orders of magnitude in field value,
// gradient and curvature because their supports and knot spacings differ. A trust
// region or a step limit stated in raw coefficients therefore means something
// different in every direction. The Sobolev Gram matrix
//
//     M_ij = integral [ B_i B_j + l1^2 grad B_i . grad B_j
//                       + l2^4 hess B_i : hess B_j ] dV
//
// measures each mode the way the physics will: value, slope and curvature, with
// two declared length scales setting the trade. Whitening makes unit steps in y
// carry unit Sobolev norm, so a later optimizer's step size has one meaning.
//
// THE VOLUME ELEMENT IS WHY s IS THE RIGHT COORDINATE, AGAIN. In cylindrical
// coordinates dV = 2 pi rho drho dx, and with s = rho^2/R0^2, xi = (x-x_c)/R0,
//
//     rho drho = (R0^2/2) ds,   dx = R0 dxi   =>   dV = pi R0^3 dxi ds
//
// so the measure is CONSTANT in the design coordinates. A quadrature in (xi, s)
// with no Jacobian weight is exact, and the axis carries no integrable
// singularity to work around.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

// --- Quadrature -------------------------------------------------------------

struct QuadratureRecord {
    std::string rule = "Gauss-Legendre, 8 nodes per knot span, tensor product";
    int nodes_per_span = 8;
    int spans_xi = 0;
    int spans_s = 0;
    std::size_t total_points = 0;
    std::string volume_element =
        "dV = pi R0^3 dxi ds; constant in the design coordinates because "
        "s = rho^2/R0^2 absorbs the cylindrical rho";
    double volume_factor = 0.0;  // pi R0^3
    int exact_to_degree = 15;    // 8-node Gauss-Legendre
};

// 8-node Gauss-Legendre on [-1,1]. Hardcoded rather than computed, so the rule
// is identical on every platform and the record can name it exactly.
inline const double* gauss8_nodes() {
    static const double n[8] = {
        -0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
        -0.1834346424956498,  0.1834346424956498,  0.5255324099163290,
         0.7966664774136267,  0.9602898564975363};
    return n;
}
inline const double* gauss8_weights() {
    static const double w[8] = {
        0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
        0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
        0.2223810344533745, 0.1012285362903763};
    return w;
}

// Quadrature points and weights over the design region, in (xi, s), built span
// by span so the rule is exact on every polynomial piece.
struct DesignQuadrature {
    std::vector<DesignCoordinates> point;
    std::vector<double> weight;  // already includes pi R0^3
    QuadratureRecord record;
};

inline DesignQuadrature build_quadrature(const AxisymmetricBasis& basis) {
    DesignQuadrature q;
    const double R0 = basis.scale_R0();
    const double kPi = 3.14159265358979323846;
    const double vol = kPi * R0 * R0 * R0;
    q.record.volume_factor = vol;

    // Span boundaries: the distinct knots inside each evaluable range.
    auto spans = [](const QuinticBSplineBasis1D& b) {
        std::vector<double> e;
        e.push_back(b.lower());
        for (double t : b.distinct_interior_knots())
            if (t > b.lower() + 1e-12 && t < b.upper() - 1e-12) e.push_back(t);
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
    };

    const std::vector<double> ex = spans(basis.axial());
    const std::vector<double> es = spans(basis.radial());
    q.record.spans_xi = static_cast<int>(ex.size()) - 1;
    q.record.spans_s = static_cast<int>(es.size()) - 1;

    const double* gn = gauss8_nodes();
    const double* gw = gauss8_weights();
    for (std::size_t a = 0; a + 1 < ex.size(); ++a) {
        const double xa = ex[a], xb = ex[a + 1];
        if (!(xb > xa)) continue;
        for (std::size_t b = 0; b + 1 < es.size(); ++b) {
            const double sa = es[b], sb = es[b + 1];
            if (!(sb > sa)) continue;
            for (int i = 0; i < 8; ++i) {
                const double xi = 0.5 * (xa + xb) + 0.5 * (xb - xa) * gn[i];
                const double wx = 0.5 * (xb - xa) * gw[i];
                for (int j = 0; j < 8; ++j) {
                    const double s = 0.5 * (sa + sb) + 0.5 * (sb - sa) * gn[j];
                    const double ws = 0.5 * (sb - sa) * gw[j];
                    q.point.push_back(
                        from_design_coordinates(xi, s, R0));
                    q.weight.push_back(wx * ws * vol);
                }
            }
        }
    }
    q.record.total_points = q.point.size();
    return q;
}

// --- Constraints ------------------------------------------------------------

struct ConstraintRow {
    std::vector<double> row;
    double rhs = 0.0;
    std::string kind;
    std::string detail;
};

class ConstraintSystem {
public:
    explicit ConstraintSystem(int n_coefficients) : n_(n_coefficients) {}

    int coefficient_count() const { return n_; }
    const std::vector<ConstraintRow>& rows() const { return rows_; }

    void add_row(std::vector<double> row, double rhs, const std::string& kind,
                 const std::string& detail) {
        rows_.push_back({std::move(row), rhs, kind, detail});
    }

    // c_index = value.
    void add_fixed(int index, double value, const std::string& detail) {
        std::vector<double> r(static_cast<std::size_t>(n_), 0.0);
        r[static_cast<std::size_t>(index)] = 1.0;
        add_row(std::move(r), value, "fixed_coefficient", detail);
    }

    // c_i - sign * c_j = 0.
    void add_reflection(int i, int j, double sign, const std::string& detail) {
        std::vector<double> r(static_cast<std::size_t>(n_), 0.0);
        r[static_cast<std::size_t>(i)] += 1.0;
        r[static_cast<std::size_t>(j)] -= sign;
        add_row(std::move(r), 0.0, "reflection_symmetry", detail);
    }

    // A weighted integral of the field, set to a value. The caller supplies the
    // already-integrated basis weights, so this one function covers the first
    // moment, the cabin mean, the outer-boundary value and any normalisation.
    void add_functional(std::vector<double> weights, double value,
                        const std::string& kind, const std::string& detail) {
        add_row(std::move(weights), value, kind, detail);
    }

    DenseMatrix matrix() const {
        DenseMatrix C(static_cast<int>(rows_.size()), n_);
        for (std::size_t i = 0; i < rows_.size(); ++i)
            for (int j = 0; j < n_; ++j)
                C(static_cast<int>(i), j) =
                    rows_[i].row[static_cast<std::size_t>(j)];
        return C;
    }

    std::vector<double> rhs() const {
        std::vector<double> d(rows_.size());
        for (std::size_t i = 0; i < rows_.size(); ++i) d[i] = rows_[i].rhs;
        return d;
    }

private:
    int n_;
    std::vector<ConstraintRow> rows_;
};

struct ConstraintProjection {
    std::vector<double> c0;   // n_raw, satisfies C c0 = d
    DenseMatrix P;            // n_raw x n_free, columns span null(C)
    int n_raw = 0;
    int n_constraints = 0;
    int rank = 0;
    int n_free = 0;
    std::vector<double> singular_values;
    double rank_tolerance = 0.0;
    double affine_residual = 0.0;     // |C c0 - d|
    double nullspace_residual = 0.0;  // max |C P|
    // max |P^T P - I|. The columns of P come from an orthogonal factor, so this
    // is zero to rounding; it is RECORDED because every downstream statement --
    // that z is an unconstrained coordinate, that the Gram matrix P^T M P
    // measures the same thing M did, that a whitened step has a defined length --
    // silently assumes the null-space basis is orthonormal.
    double orthonormality_residual = 0.0;
    int redundant_rows = 0;           // constraints that added no rank
    bool svd_converged = false;
    int svd_sweeps = 0;
};

// Build c = c0 + P z from the constraint system.
inline ConstraintProjection build_projection(const ConstraintSystem& sys) {
    ConstraintProjection pr;
    pr.n_raw = sys.coefficient_count();
    pr.n_constraints = static_cast<int>(sys.rows().size());
    const DenseMatrix C = sys.matrix();
    const std::vector<double> d = sys.rhs();

    if (pr.n_constraints == 0) {
        pr.rank = 0;
        pr.n_free = pr.n_raw;
        pr.c0.assign(static_cast<std::size_t>(pr.n_raw), 0.0);
        pr.P = DenseMatrix::identity(pr.n_raw);
        pr.svd_converged = true;
        return pr;
    }

    const SvdResult svd = jacobi_svd(C);
    pr.singular_values = svd.singular_values;
    pr.rank_tolerance = svd.tolerance;
    pr.rank = svd.rank;
    pr.n_free = pr.n_raw - pr.rank;
    // A constraint row that adds no rank is linearly dependent on the others.
    // Reported, because a dependent constraint is usually a modelling mistake --
    // the same condition written twice in two different ways.
    pr.redundant_rows = pr.n_constraints - pr.rank;

    pr.c0 = svd_min_norm_solve(C, svd, d);
    pr.P = DenseMatrix(pr.n_raw, pr.n_free);
    for (int j = 0; j < pr.n_free; ++j)
        for (int i = 0; i < pr.n_raw; ++i) pr.P(i, j) = svd.v(i, pr.rank + j);

    const std::vector<double> r = matvec(C, pr.c0);
    double acc = 0.0;
    for (std::size_t i = 0; i < r.size(); ++i) {
        const double e = r[i] - d[i];
        acc += e * e;
    }
    pr.affine_residual = std::sqrt(acc);
    pr.nullspace_residual = max_abs(multiply(C, pr.P));
    pr.orthonormality_residual = design::orthonormality_residual(pr.P);
    pr.svd_converged = svd.converged;
    pr.svd_sweeps = svd.sweeps;
    return pr;
}

// Replace c0 by another vector that also satisfies the constraints -- the F0
// baseline projection, for instance -- keeping the same null space.
inline double set_baseline(ConstraintProjection& pr, const ConstraintSystem& sys,
                          std::vector<double> baseline) {
    pr.c0 = std::move(baseline);
    const std::vector<double> r = matvec(sys.matrix(), pr.c0);
    const std::vector<double> d = sys.rhs();
    double acc = 0.0;
    for (std::size_t i = 0; i < r.size(); ++i) {
        const double e = r[i] - d[i];
        acc += e * e;
    }
    pr.affine_residual = std::sqrt(acc);
    return pr.affine_residual;
}

// --- Sobolev Gram matrix and whitening --------------------------------------

struct SobolevScales {
    double l1 = 1.0;  // gradient length scale, metres
    double l2 = 1.0;  // curvature length scale, metres
};

// M over the ACTIVE tensor products of one family.
//
// The gradient and Hessian contractions are written in (xi, s) using
// y^2 + z^2 = s R0^2, so no angular integration is needed and the axis is
// ordinary. See the derivation in the file header.
inline DenseMatrix sobolev_gram(const AxisymmetricBasis& basis,
                                const DesignQuadrature& quad,
                                const SobolevScales& scales) {
    const int n = basis.coefficient_count();
    DenseMatrix M(n, n);
    const double R0 = basis.scale_R0();
    const double l1sq = scales.l1 * scales.l1;
    const double l2q = scales.l2 * scales.l2 * scales.l2 * scales.l2;

    for (std::size_t p = 0; p < quad.point.size(); ++p) {
        const DesignCoordinates& c = quad.point[p];
        const double w = quad.weight[p];
        const TensorBasisEvaluation e = basis.evaluate(c);
        if (e.empty()) continue;
        const double s = c.radial_squared;

        const std::size_t m = e.size();
        for (std::size_t a = 0; a < m; ++a) {
            // Per-function pieces of the gradient and Hessian contractions.
            const double gx_a = e.d_xi[a] / R0;              // dB/dx
            const double gs_a = e.d_s[a];                    // dB/ds
            const double hxx_a = e.d_xi_xi[a] / (R0 * R0);
            const double hxs_a = e.d_xi_s[a] / R0;           // N' M' / R0
            const double aa = 4.0 * e.d_s_s[a] / (R0 * R0 * R0 * R0);
            const double bb = 2.0 * gs_a / (R0 * R0);
            // Only b >= a: the integrand is computed ONCE per unordered pair and
            // added to both entries, so M is EXACTLY symmetric rather than
            // symmetric to rounding. Assembling both orders independently made
            // M(i,j) and M(j,i) differ at 6e-14 -- harmless in itself, but
            // jacobi_eigen reads only the upper triangle, so any genuine
            // asymmetry would have been silently discarded instead of failing a
            // check. Exact symmetry by construction removes that hiding place.
            for (std::size_t b = a; b < m; ++b) {
                const double gx_b = e.d_xi[b] / R0;
                const double gs_b = e.d_s[b];
                const double hxx_b = e.d_xi_xi[b] / (R0 * R0);
                const double hxs_b = e.d_xi_s[b] / R0;
                const double ab = 4.0 * e.d_s_s[b] / (R0 * R0 * R0 * R0);
                const double bbb = 2.0 * gs_b / (R0 * R0);

                const double l2term = e.value[a] * e.value[b];
                // grad . grad = B_,x B_,x + 4 s /R0^2 * (dB/ds)(dB/ds)
                const double grad =
                    gx_a * gx_b + 4.0 * s * gs_a * gs_b / (R0 * R0);
                // hess : hess, from the header's algebra
                const double hess =
                    hxx_a * hxx_b +
                    2.0 * (4.0 * s / (R0 * R0)) * hxs_a * hxs_b +
                    aa * ab * s * s * R0 * R0 * R0 * R0 +
                    2.0 * bb * bbb +
                    (aa * bbb + ab * bb) * s * R0 * R0;

                const double term = w * (l2term + l1sq * grad + l2q * hess);
                const int ia = e.flat_index[a], ib = e.flat_index[b];
                M(ia, ib) += term;
                if (ia != ib) M(ib, ia) += term;
            }
        }
    }
    return M;
}

// max |M - M^T|. Zero by construction for sobolev_gram; checked rather than
// assumed, because the whole point of assembling symmetrically was to make a
// later symmetrization unnecessary.
inline double symmetry_residual(const DenseMatrix& M) {
    double worst = 0.0;
    for (int i = 0; i < M.rows; ++i)
        for (int j = 0; j < M.cols; ++j)
            worst = std::fmax(worst, std::fabs(M(i, j) - M(j, i)));
    return worst;
}

struct WhiteningResult {
    DenseMatrix w;                    // n_free x n_retained
    std::vector<double> eigenvalues;  // of G, descending
    int retained = 0;
    int discarded = 0;
    double smallest_retained = 0.0;
    double largest_retained = 0.0;
    double tolerance = 0.0;
    double raw_condition = 0.0;
    double constrained_condition = 0.0;
    double whitened_condition = 0.0;
    double whitening_residual = 0.0;  // max |W^T G W - I|
};

inline double condition_number(const std::vector<double>& eig) {
    if (eig.empty()) return 0.0;
    const double hi = eig.front();
    double lo = 0.0;
    for (auto it = eig.rbegin(); it != eig.rend(); ++it) {
        if (*it > 0.0) { lo = *it; break; }
    }
    return (lo > 0.0) ? hi / lo : 0.0;
}

// G = P^T M P, then W from the eigendecomposition of G.
//
// A mode whose eigenvalue is at or below the tolerance is DISCARDED and counted,
// never divided by: dividing by sqrt(lambda) for a numerically null lambda
// manufactures an enormous optimization direction out of rounding error, which
// an optimizer would then dutifully explore.
inline WhiteningResult whiten(const DenseMatrix& M,
                              const ConstraintProjection& pr,
                              double rel_tol = 1e-12) {
    WhiteningResult r;
    const EigenResult raw = jacobi_eigen(M);
    r.raw_condition = condition_number(raw.eigenvalues);

    const DenseMatrix G = multiply(multiply(transpose(pr.P), M), pr.P);
    const EigenResult ge = jacobi_eigen(G);
    r.eigenvalues = ge.eigenvalues;
    r.constrained_condition = condition_number(ge.eigenvalues);

    const double lmax = ge.eigenvalues.empty() ? 0.0 : ge.eigenvalues[0];
    r.tolerance = rel_tol * std::fabs(lmax);
    std::vector<int> keep;
    for (int j = 0; j < static_cast<int>(ge.eigenvalues.size()); ++j) {
        if (ge.eigenvalues[static_cast<std::size_t>(j)] > r.tolerance)
            keep.push_back(j);
        else
            ++r.discarded;
    }
    r.retained = static_cast<int>(keep.size());
    if (r.retained > 0) {
        r.largest_retained =
            ge.eigenvalues[static_cast<std::size_t>(keep.front())];
        r.smallest_retained =
            ge.eigenvalues[static_cast<std::size_t>(keep.back())];
    }

    r.w = DenseMatrix(G.rows, r.retained);
    for (int k = 0; k < r.retained; ++k) {
        const int j = keep[static_cast<std::size_t>(k)];
        const double scale =
            1.0 / std::sqrt(ge.eigenvalues[static_cast<std::size_t>(j)]);
        for (int i = 0; i < G.rows; ++i) r.w(i, k) = ge.q(i, j) * scale;
    }

    const DenseMatrix WtGW = multiply(multiply(transpose(r.w), G), r.w);
    double worst = 0.0;
    for (int i = 0; i < WtGW.rows; ++i)
        for (int j = 0; j < WtGW.cols; ++j)
            worst = std::fmax(worst, std::fabs(WtGW(i, j) -
                                               (i == j ? 1.0 : 0.0)));
    r.whitening_residual = worst;
    const EigenResult we = jacobi_eigen(WtGW);
    r.whitened_condition = condition_number(we.eigenvalues);
    return r;
}

}  // namespace design
}  // namespace warp_hypad
