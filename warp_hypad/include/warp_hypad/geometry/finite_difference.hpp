#pragma once

#include <cassert>
#include <cstddef>
#include <string>

#include "warp_hypad/core/field.hpp"
#include "warp_hypad/core/grid.hpp"
#include "warp_hypad/core/scalar_traits.hpp"

// ---------------------------------------------------------------------------
// Second-order finite-difference operators on the Cartesian grid.
//
// Interior points use centred stencils; boundary points use one-sided stencils
// of the SAME formal order, so the scheme is second order everywhere rather
// than second order in the interior and first order at the edges. That matters
// here because warp-bubble curvature plots are routinely spoiled by boundary
// artefacts, and a first-order edge would put an O(h) error ring around every
// figure.
//
// Stencils (h = spacing along the axis):
//
//   first derivative, centred        (f[i+1] - f[i-1]) / (2h)
//   first derivative, forward        (-3f[0] + 4f[1] - f[2]) / (2h)
//   first derivative, backward       (3f[n-1] - 4f[n-2] + f[n-3]) / (2h)
//
//   second derivative, centred       (f[i-1] - 2f[i] + f[i+1]) / h^2
//   second derivative, forward       (2f[0] - 5f[1] + 4f[2] - f[3]) / h^2
//   second derivative, backward      (2f[n-1] - 5f[n-2] + 4f[n-3] - f[n-4]) / h^2
//
// The one-sided second-derivative stencil deliberately uses FOUR points. The
// familiar three-point form (f0 - 2f1 + f2)/h^2 is only O(h) when applied at an
// endpoint, and using it would quietly drop the observed convergence order to 1
// on the boundary while the interior still reported 2.
//
// Mixed derivatives are formed by composing two first-derivative operators.
// Because each operator mixes only its own axis index, and the stencil choice
// at a point depends only on that axis's index, the two orderings commute
// EXACTLY (up to floating-point rounding) at every point including corners --
// not merely to truncation order. The commutativity test therefore checks
// arithmetic, not discretisation.
//
// Everything is templated on Scalar and contains no comparison of, or branch
// on, a Scalar value: stencil selection depends on the grid index only.
// Coefficients are exact small rationals lifted through scalar_cast, so they
// carry no derivative channels.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace geometry {

enum class Axis : int { X = 0, Y = 1, Z = 2 };

// Minimum extent along a differentiated axis. Three points are enough for the
// one-sided first derivative; the one-sided second derivative needs four.
inline constexpr int kMinPointsFirst = 3;
inline constexpr int kMinPointsSecond = 4;

// A stencil is at most four taps. `offsets` are index offsets along the axis
// relative to the evaluation point; `coeffs` already include the 1/h or 1/h^2
// factor.
struct Stencil {
    int offsets[4];
    double coeffs[4];
    int count;
    bool one_sided;
};

// Second-order first-derivative stencil at index `i` of an axis with `n` points
// and spacing `h`.
inline Stencil first_derivative_stencil(int i, int n, double h) {
    assert(n >= kMinPointsFirst);
    assert(i >= 0 && i < n);
    Stencil s{};
    const double inv = 1.0 / (2.0 * h);
    if (i == 0) {
        s.count = 3;
        s.one_sided = true;
        s.offsets[0] = 0; s.coeffs[0] = -3.0 * inv;
        s.offsets[1] = 1; s.coeffs[1] = 4.0 * inv;
        s.offsets[2] = 2; s.coeffs[2] = -1.0 * inv;
    } else if (i == n - 1) {
        s.count = 3;
        s.one_sided = true;
        s.offsets[0] = 0;  s.coeffs[0] = 3.0 * inv;
        s.offsets[1] = -1; s.coeffs[1] = -4.0 * inv;
        s.offsets[2] = -2; s.coeffs[2] = 1.0 * inv;
    } else {
        s.count = 2;
        s.one_sided = false;
        s.offsets[0] = -1; s.coeffs[0] = -1.0 * inv;
        s.offsets[1] = 1;  s.coeffs[1] = 1.0 * inv;
    }
    return s;
}

// Second-order second-derivative stencil at index `i`.
inline Stencil second_derivative_stencil(int i, int n, double h) {
    assert(n >= kMinPointsSecond);
    assert(i >= 0 && i < n);
    Stencil s{};
    const double inv = 1.0 / (h * h);
    if (i == 0) {
        s.count = 4;
        s.one_sided = true;
        s.offsets[0] = 0; s.coeffs[0] = 2.0 * inv;
        s.offsets[1] = 1; s.coeffs[1] = -5.0 * inv;
        s.offsets[2] = 2; s.coeffs[2] = 4.0 * inv;
        s.offsets[3] = 3; s.coeffs[3] = -1.0 * inv;
    } else if (i == n - 1) {
        s.count = 4;
        s.one_sided = true;
        s.offsets[0] = 0;  s.coeffs[0] = 2.0 * inv;
        s.offsets[1] = -1; s.coeffs[1] = -5.0 * inv;
        s.offsets[2] = -2; s.coeffs[2] = 4.0 * inv;
        s.offsets[3] = -3; s.coeffs[3] = -1.0 * inv;
    } else {
        s.count = 3;
        s.one_sided = false;
        s.offsets[0] = -1; s.coeffs[0] = 1.0 * inv;
        s.offsets[1] = 0;  s.coeffs[1] = -2.0 * inv;
        s.offsets[2] = 1;  s.coeffs[2] = 1.0 * inv;
    }
    return s;
}

// --- Axis helpers ----------------------------------------------------------

inline int axis_extent(const Grid& g, Axis a) {
    switch (a) {
        case Axis::X: return g.nx();
        case Axis::Y: return g.ny();
        case Axis::Z: return g.nz();
    }
    return 0;
}

inline double axis_spacing(const Grid& g, Axis a) {
    switch (a) {
        case Axis::X: return g.dx();
        case Axis::Y: return g.dy();
        case Axis::Z: return g.dz();
    }
    return 0.0;
}

inline std::size_t axis_stride(const Grid& g, Axis a) {
    switch (a) {
        case Axis::X: return g.stride_x();
        case Axis::Y: return g.stride_y();
        case Axis::Z: return g.stride_z();
    }
    return 0;
}

inline int axis_index(const Grid& g, std::size_t p, Axis a) {
    int i, j, k;
    g.unflatten(p, i, j, k);
    switch (a) {
        case Axis::X: return i;
        case Axis::Y: return j;
        case Axis::Z: return k;
    }
    return 0;
}

// True when the point lies on a face normal to `a`, i.e. where the one-sided
// stencil is used. Reported separately from interior error, always.
inline bool is_axis_boundary(const Grid& g, std::size_t p, Axis a) {
    const int i = axis_index(g, p, a);
    return i == 0 || i == axis_extent(g, a) - 1;
}

// Returns a human-readable description of why differentiation is impossible,
// or an empty string when the grid is usable.
std::string validate_for_derivatives(const Grid& g, bool need_second);

// --- Pointwise application -------------------------------------------------

template <typename Scalar>
Scalar apply_stencil(const ScalarField<Scalar>& f, const Grid& g, std::size_t p,
                     Axis a, const Stencil& s) {
    const std::size_t stride = axis_stride(g, a);
    Scalar acc = scalar_cast<Scalar>(0.0);
    for (int t = 0; t < s.count; ++t) {
        // Index arithmetic in signed space, then back: offsets can be negative.
        const std::size_t q =
            static_cast<std::size_t>(static_cast<std::ptrdiff_t>(p) +
                                     static_cast<std::ptrdiff_t>(s.offsets[t]) *
                                         static_cast<std::ptrdiff_t>(stride));
        acc = acc + scalar_cast<Scalar>(s.coeffs[t]) * f.at(q);
    }
    return acc;
}

template <typename Scalar>
Scalar first_derivative_at(const ScalarField<Scalar>& f, const Grid& g,
                           std::size_t p, Axis a) {
    const Stencil s = first_derivative_stencil(axis_index(g, p, a),
                                               axis_extent(g, a),
                                               axis_spacing(g, a));
    return apply_stencil(f, g, p, a, s);
}

template <typename Scalar>
Scalar second_derivative_at(const ScalarField<Scalar>& f, const Grid& g,
                            std::size_t p, Axis a) {
    const Stencil s = second_derivative_stencil(axis_index(g, p, a),
                                                axis_extent(g, a),
                                                axis_spacing(g, a));
    return apply_stencil(f, g, p, a, s);
}

// --- Whole-field operators -------------------------------------------------

template <typename Scalar>
ScalarField<Scalar> first_derivative(const ScalarField<Scalar>& f,
                                     const Grid& g, Axis a) {
    ScalarField<Scalar> out(g);
    const std::size_t n = g.num_points();
    for (std::size_t p = 0; p < n; ++p) out.at(p) = first_derivative_at(f, g, p, a);
    return out;
}

template <typename Scalar>
ScalarField<Scalar> second_derivative(const ScalarField<Scalar>& f,
                                      const Grid& g, Axis a) {
    ScalarField<Scalar> out(g);
    const std::size_t n = g.num_points();
    for (std::size_t p = 0; p < n; ++p) out.at(p) = second_derivative_at(f, g, p, a);
    return out;
}

// d^2 f / da db, formed as D_a(D_b f).
//
// When a == b this is NOT the same as second_derivative(): composing two
// centred first derivatives gives the wider (f[i+2] - 2f[i] + f[i-2])/(4h^2)
// stencil, which is also second order but has a four-times larger error
// constant. Callers wanting d^2/dx^2 should use second_derivative(); this
// function exists for genuinely mixed pairs, and the distinction is verified.
template <typename Scalar>
ScalarField<Scalar> mixed_derivative(const ScalarField<Scalar>& f,
                                     const Grid& g, Axis a, Axis b) {
    const ScalarField<Scalar> db = first_derivative(f, g, b);
    return first_derivative(db, g, a);
}

// --- Tensor-field differentiation ------------------------------------------
//
// Produces the spacetime gradient of a rank-2 field:
//
//     out(p, alpha, mu, nu) = d_alpha g_{mu nu}
//
// with alpha in {0=t, 1=x, 2=y, 3=z}. This is exactly the object the
// Christoffel assembly consumes.
//
// The alpha = 0 (temporal) slot is filled by `time_derivative`, supplied by the
// caller. It is a required argument rather than defaulting to zero so that a
// frozen-time approximation has to be stated at the call site rather than
// assumed silently -- see docs/conventions.md section 11.
template <typename Scalar>
Rank3Field<Scalar> spacetime_gradient(const Rank2Field<Scalar>& gfield,
                                      const Grid& g,
                                      const Rank3Field<Scalar>* time_derivative) {
    Rank3Field<Scalar> out(g);
    const std::size_t n = g.num_points();

    // Scratch for one component at a time. Differentiating componentwise keeps
    // the stencil code shared with the scalar path rather than duplicated.
    ScalarField<Scalar> comp(g);

    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            for (std::size_t p = 0; p < n; ++p) comp.at(p) = gfield.at(p, mu, nu);

            for (int ax = 0; ax < 3; ++ax) {
                const Axis a = static_cast<Axis>(ax);
                for (std::size_t p = 0; p < n; ++p) {
                    out.at(p, ax + 1, mu, nu) = first_derivative_at(comp, g, p, a);
                }
            }
        }
    }

    const Scalar zero = scalar_cast<Scalar>(0.0);
    for (std::size_t p = 0; p < n; ++p) {
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                out.at(p, 0, mu, nu) =
                    time_derivative ? time_derivative->at(p, 0, mu, nu) : zero;
            }
        }
    }
    return out;
}

// Frozen-time form: d_0 g_{mu nu} := 0.
//
// A separately named function rather than a defaulted null argument, so the
// approximation is spelled out at every call site. `spacetime_gradient(f, g,
// nullptr)` would read as an absent option; this reads as a physical claim,
// which is what it is. See docs/conventions.md section 11.
template <typename Scalar>
Rank3Field<Scalar> spacetime_gradient_frozen_time(const Rank2Field<Scalar>& gfield,
                                                  const Grid& g) {
    return spacetime_gradient<Scalar>(gfield, g, nullptr);
}

// Spacetime gradient of a rank-3 field:
//
//     out(p, alpha, a, b, c) = d_alpha T(p, a, b, c)
//
// The object the Riemann assembly consumes, with T the Christoffel symbol. Same
// contract as the rank-2 form above: the alpha = 0 slot is supplied by the
// caller and is a REQUIRED argument, so a vanishing temporal derivative is a
// stated claim at the call site rather than a default.
//
// Note the storage cost. The result is a Rank4Field, 256 scalars per grid point:
// 41^3 points is 141 MB in double and proportionally more with derivative
// channels. Curvature verification therefore runs on modest grids by design, and
// the memory scaling is a limitation recorded in docs/limitations.md rather than
// an oversight.
template <typename Scalar>
Rank4Field<Scalar> spacetime_gradient_rank3(const Rank3Field<Scalar>& tfield,
                                            const Grid& g,
                                            const Rank4Field<Scalar>* time_derivative) {
    Rank4Field<Scalar> out(g);
    const std::size_t n = g.num_points();

    ScalarField<Scalar> comp(g);

    for (int a = 0; a < 4; ++a) {
        for (int b = 0; b < 4; ++b) {
            for (int c = 0; c < 4; ++c) {
                for (std::size_t p = 0; p < n; ++p) comp.at(p) = tfield.at(p, a, b, c);

                for (int ax = 0; ax < 3; ++ax) {
                    const Axis axis = static_cast<Axis>(ax);
                    for (std::size_t p = 0; p < n; ++p) {
                        out.at(p, ax + 1, a, b, c) =
                            first_derivative_at(comp, g, p, axis);
                    }
                }
            }
        }
    }

    const Scalar zero = scalar_cast<Scalar>(0.0);
    for (std::size_t p = 0; p < n; ++p)
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int c = 0; c < 4; ++c)
                    out.at(p, 0, a, b, c) =
                        time_derivative ? time_derivative->at(p, 0, a, b, c) : zero;

    return out;
}

}  // namespace geometry
}  // namespace warp_hypad
