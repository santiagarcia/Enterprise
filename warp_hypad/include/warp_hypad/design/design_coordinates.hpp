#pragma once

#include <cmath>
#include <string>

// ---------------------------------------------------------------------------
// The design coordinate system, and the units contract that goes with it.
//
// THE FIVE THINGS THAT MUST NEVER BE MIXED, and how this file keeps them apart:
//
//     rho          PHYSICAL cylindrical radius, metres. Field `rho`.
//     eta          rho / R0, DIMENSIONLESS. Accessor `eta()`. No stored field:
//                  storing both eta and s invites one to drift from the other.
//     s            eta^2 = (y^2+z^2)/R0^2, DIMENSIONLESS. Field
//                  `radial_squared`, named for what it is rather than for a
//                  one-letter symbol.
//     xi           (x - x_c)/R0, DIMENSIONLESS. Field `xi`.
//     R0           the length scale itself, metres. Field `scale_R0`.
//
// A double called `r` could be any of the first four. Every one of them appears
// in this milestone, so the struct names them and the constructor is the only
// place the conversions are written.
//
// WHY s AND NOT eta IS THE RADIAL BASIS COORDINATE.
//
// A scalar field on an axisymmetric geometry must be even in rho and have
// vanishing radial derivative on the axis, or the metric it builds is singular
// there. Written in eta, a spline can violate both: nothing stops
// sum_j c_j M_j(eta) from having a nonzero slope at eta = 0, and a penalty term
// added to an optimizer's objective only discourages it.
//
// Written in s, regularity is STRUCTURAL:
//
//     q = q(xi, s),   s = rho^2/R0^2
//     d q / d rho = (dq/ds)(ds/drho) = (2 rho / R0^2) dq/ds
//
// which vanishes at rho = 0 for ANY spline coefficients, exactly, with no
// division by rho anywhere in the expression. That identity is the reason this
// milestone exists before any metric is designed.
//
// The vector and tensor families need more than evenness -- they must vanish to
// a definite ORDER on the axis -- and they get it the same structural way, from
// a polynomial prefactor rather than a penalty. See FieldFamily below.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

struct DesignCoordinates {
    double xi = 0.0;             // (x - x_c)/R0, dimensionless
    double radial_squared = 0.0; // s = (y^2+z^2)/R0^2, dimensionless
    double rho = 0.0;            // sqrt(y^2+z^2), METRES
    double scale_R0 = 1.0;       // R0, METRES

    // rho/R0. Provided as an accessor, not a field, so it cannot disagree with
    // radial_squared.
    double eta() const { return rho / scale_R0; }

    // The two Cartesian transverse coordinates are recovered only up to a
    // rotation about the axis, which is exactly the axisymmetry assumption.
    // Callers that need them keep their own y and z; see cartesian_y/z below.
    double y = 0.0;              // METRES
    double z = 0.0;              // METRES
};

// The single place the coordinate transformation is written.
//
//   x, y, z, x_c, R0 are all in METRES.
inline DesignCoordinates make_design_coordinates(double x, double y, double z,
                                                 double x_c, double R0) {
    DesignCoordinates c;
    c.scale_R0 = R0;
    c.xi = (x - x_c) / R0;
    c.y = y;
    c.z = z;
    const double rho2 = y * y + z * z;
    c.rho = std::sqrt(rho2);
    c.radial_squared = rho2 / (R0 * R0);
    return c;
}

// A point given directly in design coordinates, on the y >= 0 meridian. Used by
// the diagnostics and by quadrature, which integrate in (xi, s) and never need
// an angle.
inline DesignCoordinates from_design_coordinates(double xi, double s,
                                                 double R0) {
    DesignCoordinates c;
    c.scale_R0 = R0;
    c.xi = xi;
    c.radial_squared = s;
    c.rho = R0 * std::sqrt(std::fmax(s, 0.0));
    c.y = c.rho;
    c.z = 0.0;
    return c;
}

// --- The design region and its collars --------------------------------------

struct DesignRegion {
    double xi_min = -3.0;
    double xi_max = 3.0;
    double s_max = 4.0;  // s = 4 is rho = 2 R0

    // A correction basis function whose SUPPORT reaches into this outer band is
    // removed from the family, not merely damped. Measured in xi from xi_max
    // (and, for the radial direction, in s from s_max).
    double outer_flat_collar = 0.6;

    // The interior that must stay flat, in s. A family that protects the cabin
    // excludes every basis function whose support reaches inside it.
    double protected_cabin_radius = 0.25;  // in s

    // The band between the protected cabin and the free interior. Reported so a
    // reader can see how much room the transition has; it is not itself an
    // exclusion zone.
    double cabin_transition_collar = 0.25;  // in s

    double xi_flat_lower() const { return xi_min + outer_flat_collar; }
    double xi_flat_upper() const { return xi_max - outer_flat_collar; }
    double s_flat_outer() const { return s_max - outer_flat_collar; }
};

// --- Structurally regular field families ------------------------------------
//
// Each family is a spline scalar Q(xi, s) multiplied by a fixed prefactor in
// eta = rho/R0 whose power is exactly the vanishing order the axis demands:
//
//   Scalar               q            = Q                 (even in rho, d/drho = 0)
//   RadialVector         q_rho        = eta   * Q         (vanishes linearly)
//   AxialRadialTensor    q_{x rho}    = eta   * Q         (vanishes linearly)
//   TransverseAnisotropy q_rr - q_pp  = eta^2 * Q         (vanishes quadratically)
//
// THE CARTESIAN FORM IS WHY THIS IS STRUCTURAL RATHER THAN COSMETIC. Projecting
// the rho-component onto y gives
//
//     q_y = q_rho * (y/rho) = (rho/R0) Q (y/rho) = (y/R0) Q
//
// and the rho cancels ALGEBRAICALLY. There is no division by rho to guard, no
// small-rho branch, and no limit to take: the Cartesian component is a
// polynomial in y times a spline. The same cancellation gives
// q_{yy} - q_{zz} = ((y^2 - z^2)/R0^2) Q and q_{yz} = (y z / R0^2) Q.
enum class FieldFamily {
    Scalar,
    RadialVector,
    AxialRadialTensor,
    TransverseAnisotropy
};

inline const char* to_string(FieldFamily f) {
    switch (f) {
        case FieldFamily::Scalar: return "scalar";
        case FieldFamily::RadialVector: return "radial_vector";
        case FieldFamily::AxialRadialTensor: return "axial_radial_tensor";
        case FieldFamily::TransverseAnisotropy: return "transverse_anisotropy";
    }
    return "unknown";
}

// The power of eta in the prefactor: the exact order to which the family
// vanishes on the axis.
inline int axis_vanishing_order(FieldFamily f) {
    switch (f) {
        case FieldFamily::Scalar: return 0;
        case FieldFamily::RadialVector: return 1;
        case FieldFamily::AxialRadialTensor: return 1;
        case FieldFamily::TransverseAnisotropy: return 2;
    }
    return 0;
}

// The prefactor eta^m and its first two derivatives with respect to rho, in
// units of 1/metre and 1/metre^2. Written out rather than looped so the m = 0
// case is visibly exact.
struct AxisPrefactor {
    double value = 1.0;
    double d_rho = 0.0;
    double d_rho_rho = 0.0;
};

inline AxisPrefactor axis_prefactor(FieldFamily f, double rho, double R0) {
    AxisPrefactor p;
    const int m = axis_vanishing_order(f);
    const double e = rho / R0;
    if (m == 0) {
        p.value = 1.0;
        p.d_rho = 0.0;
        p.d_rho_rho = 0.0;
    } else if (m == 1) {
        p.value = e;
        p.d_rho = 1.0 / R0;
        p.d_rho_rho = 0.0;
    } else {
        p.value = e * e;
        p.d_rho = 2.0 * rho / (R0 * R0);
        p.d_rho_rho = 2.0 / (R0 * R0);
    }
    return p;
}

// The Cartesian prefactor for the component this family contributes, evaluated
// WITHOUT any division by rho. Which component that is depends on the family:
//
//   Scalar               the field itself
//   RadialVector         the y component of the vector
//   AxialRadialTensor    the xy component of the tensor
//   TransverseAnisotropy the yy - zz combination
inline double cartesian_prefactor(FieldFamily f, double y, double z,
                                  double R0) {
    switch (f) {
        case FieldFamily::Scalar:
            return 1.0;
        case FieldFamily::RadialVector:
        case FieldFamily::AxialRadialTensor:
            return y / R0;
        case FieldFamily::TransverseAnisotropy:
            return (y * y - z * z) / (R0 * R0);
    }
    return 1.0;
}

// The off-diagonal transverse companion of the anisotropy family, also free of
// any division by rho.
inline double cartesian_prefactor_yz(FieldFamily f, double y, double z,
                                     double R0) {
    if (f == FieldFamily::TransverseAnisotropy) return y * z / (R0 * R0);
    if (f == FieldFamily::RadialVector || f == FieldFamily::AxialRadialTensor)
        return z / R0;
    return 0.0;
}

}  // namespace design
}  // namespace warp_hypad
