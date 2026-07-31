#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "warp_hypad/design/axisymmetric_basis.hpp"
#include "warp_hypad/design/constraint_projection.hpp"
#include "warp_hypad/design/shift_family.hpp"

// ---------------------------------------------------------------------------
// O0C: F1W, the wall-only shift correction family.
//
// WHY A NEW FAMILY RATHER THAN A NEW CONSTRAINT ON F1. F1's protected-region
// contract does not survive being read carefully, and O0C is the milestone that
// has to read it, because F1 is about to become an OPTIMIZATION family rather
// than a verification family. Three findings, all measured in the gate:
//
//   THE "PROTECTED CABIN" IS A CYLINDER. DesignRegion::protected_cabin_radius
//   is 0.25 in s, and s = (y^2+z^2)/R0^2 carries no xi dependence at all. The
//   region it names is therefore rho <= R0/2, extended over the WHOLE axial
//   extent of the design region -- an infinite cylinder through the bubble, not
//   a cabin. A sphere around the bubble centre is w = xi^2 + s <= w_cabin, and
//   that is a different set.
//
//   THE CONSTRAINT FIXES A VOLUME MEAN, NOT A FIELD. o0b_constraints integrates
//   the correction over that cylinder and sets the integral to zero. A
//   correction that is +A over the front half and -A over the back half
//   satisfies it exactly while changing the interior shift everywhere in it.
//
//   TWO SCALARS CANNOT PIN A FIELD. The moment row is likewise one number. Two
//   linear functionals constrain two directions out of 36; the remaining 34 are
//   free to do whatever they like inside the "protected" region.
//
// None of that was wrong for O0B, whose job was to verify a metric jet. It is
// wrong for an optimizer, which will find exactly those 34 directions.
//
// SO F1W USES SUPPORT EXCLUSION, WHICH IS STRUCTURAL. A basis function that is
// removed from the family has no coefficient; there is no value of any accepted
// coordinate vector that makes it contribute. Delta beta^x is then zero in the
// protected regions IDENTICALLY -- not to a tolerance, not on a collocation
// set, and with every derivative zero too, because the zero function has zero
// derivatives. That is a guarantee a penalty or a quadrature constraint cannot
// give.
//
// THE GEOMETRY OF THE EXCLUSION. Each tensor-product function is supported on a
// BOX [xi_lo, xi_hi] x [s_lo, s_hi] in design coordinates, and the protected
// interior is the ball w = xi^2 + s <= w_cabin. Over a box,
//
//     min w = ( min |xi| on [xi_lo, xi_hi] )^2 + s_lo
//     max w = max( xi_lo^2, xi_hi^2 )         + s_hi
//
// both exact, because w is separable and monotone in each variable away from
// xi = 0. A function is removed when min w < w_cabin (it reaches into the ball)
// or when max w > w_collar (it reaches into the outer flat region). Removing on
// min/max rather than on centres is what makes the exclusion CONSERVATIVE:
// touching the region at all is enough.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

// The protected geometry, in design coordinates. Both radii are in w = xi^2 + s,
// i.e. (r_comoving/R0)^2, so both name BALLS and neither can be mistaken for a
// cylinder in s.
struct WallOnlyRegion {
    // The protected physical interior: a ball of comoving radius R0 sqrt(w).
    double w_cabin = 0.36;      // r_comoving <= 0.6 R0 = 1.2

    // Outside this the correction must vanish so the metric is exactly the
    // baseline in the collar and the asymptotic condition cannot be touched.
    double w_collar = 2.56;     // r_comoving <= 1.6 R0 = 3.2

    double cabin_radius_comoving(double R0) const { return R0 * std::sqrt(w_cabin); }
    double collar_radius_comoving(double R0) const { return R0 * std::sqrt(w_collar); }
};

// Why one function survived or did not. Recorded per raw tensor product.
enum class WallExclusion {
    Active,
    CabinOverlap,
    CollarOverlap,
    BothOverlap
};

inline const char* to_string(WallExclusion e) {
    switch (e) {
        case WallExclusion::Active: return "active";
        case WallExclusion::CabinOverlap: return "cabin_overlap";
        case WallExclusion::CollarOverlap: return "collar_overlap";
        case WallExclusion::BothOverlap: return "both_overlap";
    }
    return "unknown";
}

// The exact extremes of w = xi^2 + s over a support box.
inline void support_w_range(double xi_lo, double xi_hi, double s_lo,
                            double s_hi, double& w_min, double& w_max) {
    // min xi^2 is 0 if the box straddles the axis plane xi = 0, else the nearer
    // endpoint squared.
    double xi2_min;
    if (xi_lo <= 0.0 && xi_hi >= 0.0) xi2_min = 0.0;
    else xi2_min = std::fmin(xi_lo * xi_lo, xi_hi * xi_hi);
    const double xi2_max = std::fmax(xi_lo * xi_lo, xi_hi * xi_hi);
    w_min = xi2_min + s_lo;
    w_max = xi2_max + s_hi;
}

struct WallModeAccounting {
    int raw_modes = 0;              // active tensor products of the parent F1 basis
    int removed_cabin = 0;
    int removed_collar = 0;
    int removed_both = 0;
    int active_modes = 0;
    double min_support_width_xi = 0.0;
    double min_support_width_s = 0.0;
    double min_w_of_active = 0.0;   // how far the nearest active function stays out
    double max_w_of_active = 0.0;
};

// The F1W family: a mask over the parent F1 basis, not a second basis object.
//
// Keeping the parent basis intact matters. Every verified O0A/O0B property --
// the quintic kernel, the axis regularity, the tensor-product evaluation, the
// third-derivative channels -- belongs to that object, and F1W must not fork
// it. What F1W owns is the SUBSET of flat indices it is allowed to use, and a
// prolongation that writes a short coordinate vector into a full-length
// coefficient vector with structural zeros everywhere else.
class WallOnlyFamily {
public:
    WallOnlyFamily(const AxisymmetricBasis& basis, const WallOnlyRegion& region)
        : region_(region) {
        const int n = basis.coefficient_count();
        reason_.assign(static_cast<std::size_t>(n), WallExclusion::Active);
        w_min_.assign(static_cast<std::size_t>(n), 0.0);
        w_max_.assign(static_cast<std::size_t>(n), 0.0);
        acc_.raw_modes = n;
        acc_.min_support_width_xi = 1e300;
        acc_.min_support_width_s = 1e300;
        acc_.min_w_of_active = 1e300;

        for (int k = 0; k < n; ++k) {
            double xl, xh, sl, sh;
            basis.support(k, xl, xh, sl, sh);
            double wmin = 0.0, wmax = 0.0;
            support_w_range(xl, xh, sl, sh, wmin, wmax);
            w_min_[static_cast<std::size_t>(k)] = wmin;
            w_max_[static_cast<std::size_t>(k)] = wmax;

            const bool hits_cabin = wmin < region_.w_cabin;
            const bool hits_collar = wmax > region_.w_collar;
            WallExclusion r = WallExclusion::Active;
            if (hits_cabin && hits_collar) r = WallExclusion::BothOverlap;
            else if (hits_cabin) r = WallExclusion::CabinOverlap;
            else if (hits_collar) r = WallExclusion::CollarOverlap;
            reason_[static_cast<std::size_t>(k)] = r;

            if (r == WallExclusion::Active) {
                active_.push_back(k);
                acc_.min_support_width_xi =
                    std::fmin(acc_.min_support_width_xi, xh - xl);
                acc_.min_support_width_s =
                    std::fmin(acc_.min_support_width_s, sh - sl);
                acc_.min_w_of_active = std::fmin(acc_.min_w_of_active, wmin);
                acc_.max_w_of_active = std::fmax(acc_.max_w_of_active, wmax);
            } else if (r == WallExclusion::CabinOverlap) {
                ++acc_.removed_cabin;
            } else if (r == WallExclusion::CollarOverlap) {
                ++acc_.removed_collar;
            } else {
                ++acc_.removed_both;
            }
        }
        acc_.active_modes = static_cast<int>(active_.size());
        if (active_.empty()) {
            acc_.min_support_width_xi = 0.0;
            acc_.min_support_width_s = 0.0;
            acc_.min_w_of_active = 0.0;
        }
    }

    int active_count() const { return static_cast<int>(active_.size()); }
    const std::vector<int>& active_indices() const { return active_; }
    WallExclusion reason(int flat) const {
        return reason_[static_cast<std::size_t>(flat)];
    }
    double support_w_min(int flat) const {
        return w_min_[static_cast<std::size_t>(flat)];
    }
    double support_w_max(int flat) const {
        return w_max_[static_cast<std::size_t>(flat)];
    }
    const WallModeAccounting& accounting() const { return acc_; }
    const WallOnlyRegion& region() const { return region_; }

    // Short F1W vector -> full-length parent coefficient vector. Every index
    // this family does not own is written as an EXACT zero, which is what makes
    // the protected-region guarantee structural.
    template <typename Scalar>
    std::vector<Scalar> prolong(const std::vector<Scalar>& reduced,
                                int parent_count) const {
        std::vector<Scalar> full(static_cast<std::size_t>(parent_count),
                                 scalar_cast<Scalar>(0.0));
        for (std::size_t j = 0; j < active_.size() && j < reduced.size(); ++j)
            full[static_cast<std::size_t>(active_[j])] = reduced[j];
        return full;
    }

private:
    WallOnlyRegion region_;
    std::vector<int> active_;
    std::vector<WallExclusion> reason_;
    std::vector<double> w_min_;
    std::vector<double> w_max_;
    WallModeAccounting acc_;
};

// The F1W constraint set, expressed in the REDUCED coordinates.
//
// Support exclusion already forbids a correction from touching the cabin or the
// collar. What it does not forbid is a correction that, inside the wall, looks
// like the derivative of the baseline -- which is a translation of the bubble --
// or like the baseline itself, which is a change of velocity. Both are behaviour
// changes disguised as improvements, so both get a row.
//
//     translation row:  integral of B_k * d_xi beta_base   over the region
//     amplitude row:    integral of B_k * beta_base        over the region
//
// Two rows, both from the same verified quadrature O0A built, both measured
// against the same baseline table the metric jet uses.
inline ConstraintSystem wall_only_constraints(const AxisymmetricBasis& basis,
                                              const WallOnlyFamily& family,
                                              const DesignQuadrature& quad,
                                              const ShiftFamilyConfig& cfg) {
    const int n = family.active_count();
    ConstraintSystem sys(n);
    std::vector<double> translation(static_cast<std::size_t>(n), 0.0);
    std::vector<double> amplitude(static_cast<std::size_t>(n), 0.0);

    // Reverse map: parent flat index -> reduced index, or -1.
    std::vector<int> reduced_of(
        static_cast<std::size_t>(basis.coefficient_count()), -1);
    for (int j = 0; j < n; ++j)
        reduced_of[static_cast<std::size_t>(family.active_indices()[
            static_cast<std::size_t>(j)])] = j;

    for (std::size_t p = 0; p < quad.point.size(); ++p) {
        const DesignCoordinates& q = quad.point[p];
        const DesignDerivativeTable base = baseline_shift_table(
            q.xi, q.radial_squared, cfg.velocity, cfg.scale_R0,
            cfg.bubble_radius, cfg.wall_sigma);
        const double b0 = base.at(0, 0);      // beta_base
        const double b_xi = base.at(1, 0);    // d beta_base / d xi
        const TensorBasisEvaluation e = basis.evaluate(q);
        for (std::size_t k = 0; k < e.size(); ++k) {
            const int r = reduced_of[static_cast<std::size_t>(e.flat_index[k])];
            if (r < 0) continue;
            translation[static_cast<std::size_t>(r)] +=
                quad.weight[p] * b_xi * e.value[k];
            amplitude[static_cast<std::size_t>(r)] +=
                quad.weight[p] * b0 * e.value[k];
        }
    }

    sys.add_functional(
        translation, 0.0, "translation_moment",
        "zero overlap with d_xi beta_base: removes the component of a "
        "correction that merely translates the bubble");
    sys.add_functional(
        amplitude, 0.0, "amplitude_moment",
        "zero overlap with beta_base: removes the component that merely "
        "rescales the requested velocity");
    return sys;
}

// The blocking decision. F1W is only usable if support exclusion leaves enough
// independent directions AND the whitened Gram is not singular. Both are
// recorded; neither is worked around by weakening the protected region.
struct WallFamilyStatus {
    bool blocked = false;
    std::string blocking_reason;
    int active_modes = 0;
    int constraint_rank = 0;
    int free_dimension = 0;
    int whitened_dimension = 0;
    double gram_condition_number = 0.0;
    double smallest_eigenvalue = 0.0;
    bool resolution_certified = false;
    std::string resolution_note;
};

}  // namespace design
}  // namespace warp_hypad
