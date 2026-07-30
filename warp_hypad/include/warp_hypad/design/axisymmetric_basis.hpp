#pragma once

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "warp_hypad/core/scalar_traits.hpp"
#include "warp_hypad/design/design_coordinates.hpp"
#include "warp_hypad/design/quintic_bspline.hpp"

// ---------------------------------------------------------------------------
// The axisymmetric tensor-product basis B_ij(xi, s) = N_i(xi) M_j(s), its
// Cartesian derivatives, the structurally regular field families built on it,
// and the sparse stencils a fixed grid precomputes.
//
// HOW THE FLAT COLLAR IS MADE EXACT, AND WHY IT DICTATES THE KNOT VECTORS.
//
// Section 6 of the O0A specification requires that a correction vanish -- with
// every derivative -- throughout the outer collar, and says how: no active basis
// function may OVERLAP that collar. That is a statement about supports, and it
// is not satisfiable by an arbitrary knot vector.
//
// A clamped degree-5 basis of n functions on an interval of width L has knot
// spacing L/(n-5) and support width 6L/(n-5), which exceeds L for every n <= 11.
// Every function then touches both ends, so support-exclusion of any collar at
// all would empty the family. Setting the first and last five coefficients to
// zero instead gives C4 contact at the boundary, but consumes ten coefficients
// and so also empties a ten-function axial space.
//
// The construction that works is a UNIFORM lattice, padded by the degree on each
// side so the recurrence has the knots it needs, whose raw space is exactly the
// set of functions with support inside the design region. Each such function has
// support width 6h and tapers to zero with C4 contact at both ends of its own
// support, so a family assembled from the collar-admissible subset is
// identically zero in the collar AND smooth across its inner edge. Nothing is
// penalised and no limit is taken; a function either has support in the collar
// and is removed, or it does not and cannot contribute there.
//
// With the mandated raw spaces this produces, for the design region below:
//
//     axial   raw 10, collar-excluded 4, admissible 6
//     radial  raw  8, collar-excluded 2, admissible 6
//
// The radial lattice is clamped at s = 0 because the axis is a coordinate
// singularity of the chart, not a boundary of the region: a correction must be
// free to be nonzero at the bubble centre, and s already makes that regular.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

// Why a basis function is not in the active set. Every removed function carries
// one of these, so a shrunken family is always explained.
enum class ExclusionReason {
    Active,
    OutsideDesignRegion,   // lattice padding: support not inside the region
    OuterFlatCollar,
    ProtectedCabin,
    SymmetryConstraint,
    MomentConstraint,
    LinearlyDependent
};

inline const char* to_string(ExclusionReason r) {
    switch (r) {
        case ExclusionReason::Active: return "active";
        case ExclusionReason::OutsideDesignRegion: return "outside_design_region";
        case ExclusionReason::OuterFlatCollar: return "outer_flat_collar";
        case ExclusionReason::ProtectedCabin: return "protected_cabin";
        case ExclusionReason::SymmetryConstraint: return "symmetry_constraint";
        case ExclusionReason::MomentConstraint: return "moment_constraint";
        case ExclusionReason::LinearlyDependent: return "linearly_dependent";
    }
    return "unknown";
}

// --- Lattice construction ---------------------------------------------------

struct LatticeClassification {
    std::vector<int> raw;                       // support inside the region
    std::vector<int> active;                    // raw and collar-admissible
    std::vector<ExclusionReason> reason;        // one per basis index
    double spacing = 0.0;
};

// A uniform quintic lattice on [lo, hi], padded by the degree at each end, whose
// raw space is exactly `n_raw` functions with support inside [lo, hi].
inline QuinticBSplineBasis1D make_uniform_lattice(double lo, double hi,
                                                  int n_raw, double* spacing) {
    const int p = kQuinticDegree;
    const double h = (hi - lo) / static_cast<double>(n_raw + p);
    if (spacing) *spacing = h;
    std::vector<double> knots;
    knots.reserve(static_cast<std::size_t>(n_raw + 2 * p + 6));
    for (int k = 0; k < n_raw + 2 * p + 6; ++k)
        knots.push_back(lo + (static_cast<double>(k) - p) * h);
    return QuinticBSplineBasis1D(std::move(knots));
}

// A lattice clamped at `lo` and uniform beyond, whose raw space is exactly
// `n_raw` functions with support inside [lo, hi]. Used for the radial direction,
// where lo is the axis.
inline QuinticBSplineBasis1D make_axis_clamped_lattice(double lo, double hi,
                                                       int n_raw,
                                                       double* spacing) {
    const int p = kQuinticDegree;
    const int K = n_raw + p;
    const double h = (hi - lo) / static_cast<double>(K - p);
    if (spacing) *spacing = h;
    std::vector<double> knots;
    for (int k = 0; k <= p; ++k) knots.push_back(lo);
    for (int k = 1; k <= K; ++k)
        knots.push_back(lo + static_cast<double>(k) * h);
    return QuinticBSplineBasis1D(std::move(knots));
}

// Classify every function of a lattice: padding, collar-excluded, or active.
inline LatticeClassification classify_lattice(const QuinticBSplineBasis1D& b,
                                              double region_lo,
                                              double region_hi,
                                              double free_lo, double free_hi) {
    LatticeClassification c;
    const double tol = 1e-12;
    c.reason.assign(static_cast<std::size_t>(b.size()),
                    ExclusionReason::OutsideDesignRegion);
    for (int i = 0; i < b.size(); ++i) {
        const double lo = b.support_lower(i), hi = b.support_upper(i);
        const bool in_region = (lo >= region_lo - tol && hi <= region_hi + tol);
        if (!in_region) continue;
        c.raw.push_back(i);
        if (lo >= free_lo - tol && hi <= free_hi + tol) {
            c.reason[static_cast<std::size_t>(i)] = ExclusionReason::Active;
            c.active.push_back(i);
        } else {
            c.reason[static_cast<std::size_t>(i)] =
                ExclusionReason::OuterFlatCollar;
        }
    }
    return c;
}

// --- Tensor-product evaluation ----------------------------------------------

struct TensorBasisEvaluation {
    // One entry per active tensor-product function at this point.
    std::vector<int> flat_index;   // index into the family's coefficient vector
    std::vector<int> index_xi;
    std::vector<int> index_s;

    std::vector<double> value;
    // Derivatives in the DESIGN coordinates (xi, s), dimensionless.
    std::vector<double> d_xi, d_s;
    std::vector<double> d_xi_xi, d_xi_s, d_s_s;
    std::vector<double> d3_xi, d3_s;

    // Derivatives in PHYSICAL Cartesian coordinates, per metre and per metre^2.
    std::vector<double> d_x, d_y, d_z;
    std::vector<double> d_xx, d_xy, d_xz, d_yy, d_yz, d_zz;

    // The radial derivative, written as (2 rho/R0^2) dq/ds. Exactly zero on the
    // axis, with no division by rho anywhere in its construction.
    std::vector<double> d_rho;

    std::size_t size() const { return flat_index.size(); }
    bool empty() const { return flat_index.empty(); }
};

// The tensor-product basis of one field family.
class AxisymmetricBasis {
public:
    AxisymmetricBasis(QuinticBSplineBasis1D axial,
                      QuinticBSplineBasis1D radial,
                      LatticeClassification axial_class,
                      LatticeClassification radial_class,
                      DesignRegion region, double scale_R0,
                      FieldFamily family)
        : axial_(std::move(axial)), radial_(std::move(radial)),
          axial_class_(std::move(axial_class)),
          radial_class_(std::move(radial_class)), region_(region),
          scale_R0_(scale_R0), family_(family) {
        // The coefficient vector is indexed over the ACTIVE tensor products
        // only, in a fixed (axial-major) order. A function that was excluded has
        // no coefficient at all, which is why the collar guarantee cannot be
        // broken by a stray nonzero.
        flat_of_.assign(
            static_cast<std::size_t>(axial_.size()) *
                static_cast<std::size_t>(radial_.size()), -1);
        for (int ia : axial_class_.active)
            for (int jb : radial_class_.active) {
                flat_of_[pair_key(ia, jb)] = static_cast<int>(pairs_.size());
                pairs_.push_back({ia, jb});
            }
    }

    int coefficient_count() const { return static_cast<int>(pairs_.size()); }
    int raw_coefficient_count() const {
        return static_cast<int>(axial_class_.raw.size() *
                                radial_class_.raw.size());
    }
    const std::vector<std::pair<int, int>>& pairs() const { return pairs_; }
    const QuinticBSplineBasis1D& axial() const { return axial_; }
    const QuinticBSplineBasis1D& radial() const { return radial_; }
    const LatticeClassification& axial_class() const { return axial_class_; }
    const LatticeClassification& radial_class() const { return radial_class_; }
    const DesignRegion& region() const { return region_; }
    double scale_R0() const { return scale_R0_; }
    FieldFamily family() const { return family_; }

    // Support of the flat-indexed function, in design coordinates.
    void support(int flat, double& xi_lo, double& xi_hi, double& s_lo,
                 double& s_hi) const {
        const std::pair<int, int>& p = pairs_[static_cast<std::size_t>(flat)];
        xi_lo = axial_.support_lower(p.first);
        xi_hi = axial_.support_upper(p.first);
        s_lo = radial_.support_lower(p.second);
        s_hi = radial_.support_upper(p.second);
    }

    // Every derivative the design needs, at one point.
    //
    // THE CHAIN RULE, WRITTEN OUT. With xi = (x - x_c)/R0 and s = (y^2+z^2)/R0^2,
    //
    //     ds/dy = 2y/R0^2,   ds/dz = 2z/R0^2,   d2s/dy2 = d2s/dz2 = 2/R0^2
    //
    // so for B = N(xi) M(s):
    //
    //     B_,y  = N M' (2y/R0^2)
    //     B_,yy = N M'' (2y/R0^2)^2 + N M' (2/R0^2)
    //     B_,yz = N M'' (2y/R0^2)(2z/R0^2)
    //
    // Every one of these is a polynomial in y and z times a spline derivative.
    // At y = z = 0 they are finite and the first derivatives vanish identically.
    // There is no 1/rho anywhere, so there is no axis limit to take and no
    // small-rho branch that could be reached with the wrong sign.
    TensorBasisEvaluation evaluate(const DesignCoordinates& c,
                                   bool want_third = false) const {
        TensorBasisEvaluation out;
        const BasisEvaluation1D a = axial_.evaluate(c.xi);
        const BasisEvaluation1D r = radial_.evaluate(c.radial_squared);
        if (a.empty() || r.empty()) return out;

        const double R0 = scale_R0_;
        const double inv_R0 = 1.0 / R0;
        const double sy = 2.0 * c.y / (R0 * R0);
        const double sz = 2.0 * c.z / (R0 * R0);
        const double syy = 2.0 / (R0 * R0);

        for (std::size_t ka = 0; ka < a.size(); ++ka) {
            const int ia = a.active_indices[ka];
            if (ia < 0 || ia >= axial_.size()) continue;
            for (std::size_t kr = 0; kr < r.size(); ++kr) {
                const int jb = r.active_indices[kr];
                if (jb < 0 || jb >= radial_.size()) continue;
                const int flat = flat_of_[pair_key(ia, jb)];
                if (flat < 0) continue;  // excluded: contributes nothing, ever

                const double N = a.value[ka], N1 = a.first[ka],
                             N2 = a.second[ka], N3 = a.third[ka];
                const double M = r.value[kr], M1 = r.first[kr],
                             M2 = r.second[kr], M3 = r.third[kr];

                out.flat_index.push_back(flat);
                out.index_xi.push_back(ia);
                out.index_s.push_back(jb);

                out.value.push_back(N * M);
                out.d_xi.push_back(N1 * M);
                out.d_s.push_back(N * M1);
                out.d_xi_xi.push_back(N2 * M);
                out.d_xi_s.push_back(N1 * M1);
                out.d_s_s.push_back(N * M2);
                if (want_third) {
                    out.d3_xi.push_back(N3 * M);
                    out.d3_s.push_back(N * M3);
                }

                out.d_x.push_back(N1 * M * inv_R0);
                out.d_y.push_back(N * M1 * sy);
                out.d_z.push_back(N * M1 * sz);

                out.d_xx.push_back(N2 * M * inv_R0 * inv_R0);
                out.d_xy.push_back(N1 * M1 * sy * inv_R0);
                out.d_xz.push_back(N1 * M1 * sz * inv_R0);
                out.d_yy.push_back(N * M2 * sy * sy + N * M1 * syy);
                out.d_zz.push_back(N * M2 * sz * sz + N * M1 * syy);
                out.d_yz.push_back(N * M2 * sy * sz);

                // dq/drho = (2 rho / R0^2) dq/ds. Zero at rho = 0 by the
                // prefactor, not by cancellation.
                out.d_rho.push_back(N * M1 * 2.0 * c.rho / (R0 * R0));
            }
        }
        return out;
    }

    // The family's field value at a point: the prefactor times the spline sum.
    // Templated on the coefficient scalar; the basis itself is always double.
    template <typename Scalar>
    Scalar field_value(const std::vector<Scalar>& coefficients,
                       const DesignCoordinates& c) const {
        const TensorBasisEvaluation e = evaluate(c);
        Scalar acc = scalar_cast<Scalar>(0.0);
        for (std::size_t k = 0; k < e.size(); ++k)
            acc = acc + coefficients[static_cast<std::size_t>(e.flat_index[k])] *
                            scalar_cast<Scalar>(e.value[k]);
        const AxisPrefactor pre = axis_prefactor(family_, c.rho, scale_R0_);
        return acc * scalar_cast<Scalar>(pre.value);
    }

    // The Cartesian component this family contributes, again with no division by
    // rho: for the vector and tensor families the 1/rho of the projection
    // cancels the rho of the prefactor algebraically.
    template <typename Scalar>
    Scalar cartesian_component(const std::vector<Scalar>& coefficients,
                               const DesignCoordinates& c) const {
        const TensorBasisEvaluation e = evaluate(c);
        Scalar acc = scalar_cast<Scalar>(0.0);
        for (std::size_t k = 0; k < e.size(); ++k)
            acc = acc + coefficients[static_cast<std::size_t>(e.flat_index[k])] *
                            scalar_cast<Scalar>(e.value[k]);
        return acc * scalar_cast<Scalar>(
                         cartesian_prefactor(family_, c.y, c.z, scale_R0_));
    }

private:
    std::size_t pair_key(int ia, int jb) const {
        return static_cast<std::size_t>(ia) *
                   static_cast<std::size_t>(radial_.size()) +
               static_cast<std::size_t>(jb);
    }

    QuinticBSplineBasis1D axial_;
    QuinticBSplineBasis1D radial_;
    LatticeClassification axial_class_;
    LatticeClassification radial_class_;
    DesignRegion region_;
    double scale_R0_ = 1.0;
    FieldFamily family_ = FieldFamily::Scalar;
    std::vector<std::pair<int, int>> pairs_;
    std::vector<int> flat_of_;
};

// The standard F1 kernel: the mandated raw spaces, on the stated region.
inline AxisymmetricBasis make_f1_basis(const DesignRegion& region, double R0,
                                       FieldFamily family,
                                       int n_axial_raw = 10,
                                       int n_radial_raw = 8,
                                       bool protect_cabin = false) {
    double h_xi = 0.0, h_s = 0.0;
    QuinticBSplineBasis1D axial = make_uniform_lattice(
        region.xi_min, region.xi_max, n_axial_raw, &h_xi);
    LatticeClassification ac =
        classify_lattice(axial, region.xi_min, region.xi_max,
                         region.xi_flat_lower(), region.xi_flat_upper());
    ac.spacing = h_xi;

    // A cabin-protecting family excludes the axis band as well, so its radial
    // lattice tiles [cabin, s_flat_outer] and is NOT clamped at the axis: it
    // must vanish there, which is the whole point of protecting the cabin.
    QuinticBSplineBasis1D radial =
        protect_cabin
            ? make_uniform_lattice(region.protected_cabin_radius,
                                   region.s_max, n_radial_raw, &h_s)
            : make_axis_clamped_lattice(0.0, region.s_max, n_radial_raw, &h_s);
    LatticeClassification rc = classify_lattice(
        radial, protect_cabin ? region.protected_cabin_radius : 0.0,
        region.s_max, protect_cabin ? region.protected_cabin_radius : 0.0,
        region.s_flat_outer());
    rc.spacing = h_s;
    if (protect_cabin) {
        // Relabel: a function removed because its support reaches below the
        // cabin radius is not a collar exclusion, and the record must say so.
        for (int i = 0; i < radial.size(); ++i) {
            if (rc.reason[static_cast<std::size_t>(i)] !=
                ExclusionReason::OuterFlatCollar)
                continue;
            if (radial.support_lower(i) < region.protected_cabin_radius - 1e-12)
                rc.reason[static_cast<std::size_t>(i)] =
                    ExclusionReason::ProtectedCabin;
        }
    }

    return AxisymmetricBasis(std::move(axial), std::move(radial), std::move(ac),
                             std::move(rc), region, R0, family);
}

// --- Sparse stencils --------------------------------------------------------

// One grid point's worth of basis data. The layout is deliberately
// point-major-sparse: a dense N_grid x N_basis array would be 97% zeros here,
// because only 36 of the tensor products are active and at most 36 of those
// touch any single point.
struct SparseBasisStencil {
    std::vector<int> index;
    std::vector<double> value;
    std::vector<double> d_x, d_y, d_z;
    std::vector<double> d_xx, d_xy, d_xz, d_yy, d_yz, d_zz;
    std::vector<double> d3_xi, d3_s;
    bool has_third = false;
    double prefactor = 1.0;            // the family's axis prefactor at this point
    double cartesian_prefactor_ = 1.0;

    std::size_t nonzeros() const { return index.size(); }
};

// The required evaluation entry point. Basis values are double; only the
// coefficients carry derivative channels, so an OTI coefficient vector
// propagates sensitivities through a sum of double-weighted terms.
template <typename Scalar>
Scalar evaluate(const std::vector<Scalar>& coefficients,
                const SparseBasisStencil& stencil) {
    Scalar acc = scalar_cast<Scalar>(0.0);
    for (std::size_t k = 0; k < stencil.index.size(); ++k)
        acc = acc + coefficients[static_cast<std::size_t>(stencil.index[k])] *
                        scalar_cast<Scalar>(stencil.value[k]);
    return acc * scalar_cast<Scalar>(stencil.prefactor);
}

// The same sum against a chosen derivative channel of the stencil.
enum class StencilDerivative { Value, Dx, Dy, Dz, Dxx, Dxy, Dxz, Dyy, Dyz, Dzz };

template <typename Scalar>
Scalar evaluate_derivative(const std::vector<Scalar>& coefficients,
                           const SparseBasisStencil& stencil,
                           StencilDerivative which) {
    const std::vector<double>* w = &stencil.value;
    switch (which) {
        case StencilDerivative::Value: w = &stencil.value; break;
        case StencilDerivative::Dx: w = &stencil.d_x; break;
        case StencilDerivative::Dy: w = &stencil.d_y; break;
        case StencilDerivative::Dz: w = &stencil.d_z; break;
        case StencilDerivative::Dxx: w = &stencil.d_xx; break;
        case StencilDerivative::Dxy: w = &stencil.d_xy; break;
        case StencilDerivative::Dxz: w = &stencil.d_xz; break;
        case StencilDerivative::Dyy: w = &stencil.d_yy; break;
        case StencilDerivative::Dyz: w = &stencil.d_yz; break;
        case StencilDerivative::Dzz: w = &stencil.d_zz; break;
    }
    Scalar acc = scalar_cast<Scalar>(0.0);
    for (std::size_t k = 0; k < stencil.index.size(); ++k)
        acc = acc + coefficients[static_cast<std::size_t>(stencil.index[k])] *
                        scalar_cast<Scalar>((*w)[k]);
    return acc;
}

struct StencilTableStatistics {
    std::size_t points = 0;
    std::size_t basis_functions = 0;
    std::size_t nonzeros = 0;
    std::size_t max_active = 0;
    std::size_t min_active = 0;
    double average_active = 0.0;
    std::size_t dense_equivalent_bytes = 0;
    std::size_t sparse_bytes = 0;
    double reduction_factor = 0.0;
    int stored_channels = 0;
};

class SparseStencilTable {
public:
    void build(const AxisymmetricBasis& basis,
               const std::vector<DesignCoordinates>& points,
               bool want_third = false) {
        stencils_.clear();
        stencils_.reserve(points.size());
        for (const DesignCoordinates& c : points) {
            const TensorBasisEvaluation e = basis.evaluate(c, want_third);
            SparseBasisStencil st;
            st.has_third = want_third;
            st.index = e.flat_index;
            st.value = e.value;
            st.d_x = e.d_x; st.d_y = e.d_y; st.d_z = e.d_z;
            st.d_xx = e.d_xx; st.d_xy = e.d_xy; st.d_xz = e.d_xz;
            st.d_yy = e.d_yy; st.d_yz = e.d_yz; st.d_zz = e.d_zz;
            if (want_third) { st.d3_xi = e.d3_xi; st.d3_s = e.d3_s; }
            st.prefactor =
                axis_prefactor(basis.family(), c.rho, basis.scale_R0()).value;
            st.cartesian_prefactor_ =
                cartesian_prefactor(basis.family(), c.y, c.z, basis.scale_R0());
            stencils_.push_back(std::move(st));
        }
        n_basis_ = static_cast<std::size_t>(basis.coefficient_count());
        channels_ = want_third ? 12 : 10;
    }

    const SparseBasisStencil& at(std::size_t p) const { return stencils_[p]; }
    std::size_t size() const { return stencils_.size(); }

    StencilTableStatistics statistics() const {
        StencilTableStatistics s;
        s.points = stencils_.size();
        s.basis_functions = n_basis_;
        s.stored_channels = channels_;
        s.min_active = stencils_.empty() ? 0 : stencils_[0].nonzeros();
        for (const SparseBasisStencil& st : stencils_) {
            const std::size_t n = st.nonzeros();
            s.nonzeros += n;
            if (n > s.max_active) s.max_active = n;
            if (n < s.min_active) s.min_active = n;
        }
        s.average_active =
            s.points ? static_cast<double>(s.nonzeros) /
                           static_cast<double>(s.points)
                     : 0.0;
        // A dense table would store every channel for every (point, function)
        // pair; the sparse one stores an int index plus the same channels only
        // where the function is actually nonzero.
        s.dense_equivalent_bytes = s.points * s.basis_functions *
                                   static_cast<std::size_t>(s.stored_channels) *
                                   sizeof(double);
        s.sparse_bytes =
            s.nonzeros * (static_cast<std::size_t>(s.stored_channels) *
                              sizeof(double) +
                          sizeof(int));
        s.reduction_factor =
            s.sparse_bytes ? static_cast<double>(s.dense_equivalent_bytes) /
                                 static_cast<double>(s.sparse_bytes)
                           : 0.0;
        return s;
    }

private:
    std::vector<SparseBasisStencil> stencils_;
    std::size_t n_basis_ = 0;
    int channels_ = 10;
};

}  // namespace design
}  // namespace warp_hypad
