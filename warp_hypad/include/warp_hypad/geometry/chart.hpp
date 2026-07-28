#pragma once

#include <string>

// ---------------------------------------------------------------------------
// Coordinate charts.
//
// The Grid stores three axes. Up to A3 those axes were always (x, y, z) and the
// identification was implicit. Curvature verification breaks that: the strongest
// available closed-form references -- Schwarzschild and static de Sitter -- are
// written in (r, theta, phi), and the flat-space-in-curvilinear-coordinates test
// exists precisely to separate "the connection vanishes" from "the curvature
// vanishes".
//
// Nothing about the finite-difference operators changes. A stencil differentiates
// with respect to a grid axis, and which coordinate that axis carries is a
// property of the metric being sampled, not of the differencing. What DOES change
// is that a results directory is now ambiguous unless the chart is recorded, so
// the chart is a first-class object that travels into metadata rather than a
// comment.
//
// The zeroth coordinate is t in every chart this project uses. Charts here
// describe the three GRIDDED coordinates only.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace geometry {

enum class ChartKind {
    // (t, x, y, z). Minkowski, Alcubierre.
    CartesianTXYZ,

    // (t, r, theta, phi). Schwarzschild, static de Sitter, flat space in
    // spherical coordinates. The grid axes are (r, theta, phi) in that order.
    //
    // The domain must exclude r = 0 and the polar axis sin(theta) = 0, where
    // the chart itself is singular. That is a coordinate artefact, not a
    // curvature singularity, and a grid that straddles it produces garbage that
    // looks like physics.
    SphericalTRThetaPhi,

    // (t, rho, theta, phi) with rho the isotropic radial coordinate. Used only
    // to check that scalar invariants of a genuinely curved spacetime agree
    // between two charts.
    IsotropicTRhoThetaPhi
};

struct Chart {
    ChartKind kind = ChartKind::CartesianTXYZ;

    const char* axis0() const {
        switch (kind) {
            case ChartKind::CartesianTXYZ: return "x";
            case ChartKind::SphericalTRThetaPhi: return "r";
            case ChartKind::IsotropicTRhoThetaPhi: return "rho";
        }
        return "?";
    }
    const char* axis1() const {
        return kind == ChartKind::CartesianTXYZ ? "y" : "theta";
    }
    const char* axis2() const {
        return kind == ChartKind::CartesianTXYZ ? "z" : "phi";
    }

    const char* name() const {
        switch (kind) {
            case ChartKind::CartesianTXYZ: return "cartesian_t_x_y_z";
            case ChartKind::SphericalTRThetaPhi: return "spherical_t_r_theta_phi";
            case ChartKind::IsotropicTRhoThetaPhi: return "isotropic_t_rho_theta_phi";
        }
        return "unknown";
    }

    // True where the chart has a coordinate singularity that no grid may
    // straddle. Reported so a caller can refuse a domain rather than produce
    // large finite numbers that are pure discretisation of a 1/r pole.
    bool has_axis_singularity() const {
        return kind != ChartKind::CartesianTXYZ;
    }
};

}  // namespace geometry
}  // namespace warp_hypad
