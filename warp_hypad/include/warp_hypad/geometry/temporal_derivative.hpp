#pragma once

#include <functional>

#include "warp_hypad/core/scalar_traits.hpp"
#include "warp_hypad/geometry/finite_difference.hpp"

// ---------------------------------------------------------------------------
// Temporal derivative operators.
//
// Time is NOT gridded in this project. The metric is a closed-form function of
// t, so a temporal derivative is obtained by re-evaluating that function at
// shifted times rather than by storing a 4D array. This costs a few metric
// evaluations instead of a factor of n_t in memory, and it removes the question
// of what the "temporal boundary" of the domain even is.
//
// The stencils are the SAME ones used spatially -- literally the same
// coefficient generator from finite_difference.hpp -- so the temporal scheme is
// second order interior and second order one-sided at an end, and the
// convergence study measures the same object. Sharing the generator is
// deliberate: two independently written stencil tables would eventually
// disagree.
//
// "Interior" versus "boundary" still has meaning here even without a grid. A
// caller integrating over a finite time window may not be able to sample
// outside it, so the one-sided forms are provided and verified; a caller
// evaluating an analytic metric can always sample both sides and should use the
// centred form.
//
// Everything is templated on Scalar. The sampling callable returns Scalar, so
// the same operator differentiates a double metric and an OTIScalar metric --
// and in the latter case the OTI parameter channels flow through the temporal
// derivative automatically, because the stencil is a fixed-coefficient linear
// combination.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace geometry {

// Where in a time window the evaluation point sits. Chooses the stencil exactly
// as the grid index does spatially.
enum class TimeStencilPosition {
    Interior,   // both t-h and t+h are available: centred
    WindowStart,  // only t and later: forward
    WindowEnd,    // only t and earlier: backward
};

namespace detail {

inline int stencil_index(TimeStencilPosition p, int n) {
    switch (p) {
        case TimeStencilPosition::WindowStart: return 0;
        case TimeStencilPosition::WindowEnd: return n - 1;
        case TimeStencilPosition::Interior: return n / 2;
    }
    return n / 2;
}

// A synthetic axis long enough that first_derivative_stencil() selects the
// centred form for an interior index. The stencil's offsets are what we use;
// the axis itself never exists.
inline constexpr int kSyntheticExtentFirst = 9;
inline constexpr int kSyntheticExtentSecond = 9;

}  // namespace detail

// d/dt of a callable f(t) -> Scalar, at time t with step h.
template <typename Scalar, typename F>
Scalar time_derivative(F&& f, double t, double h,
                       TimeStencilPosition pos = TimeStencilPosition::Interior) {
    const Stencil s = first_derivative_stencil(
        detail::stencil_index(pos, detail::kSyntheticExtentFirst),
        detail::kSyntheticExtentFirst, h);
    Scalar acc = scalar_cast<Scalar>(0.0);
    for (int i = 0; i < s.count; ++i) {
        acc = acc + scalar_cast<Scalar>(s.coeffs[i]) *
                        f(t + static_cast<double>(s.offsets[i]) * h);
    }
    return acc;
}

// d2/dt2 of a callable f(t) -> Scalar.
template <typename Scalar, typename F>
Scalar time_second_derivative(
    F&& f, double t, double h,
    TimeStencilPosition pos = TimeStencilPosition::Interior) {
    const Stencil s = second_derivative_stencil(
        detail::stencil_index(pos, detail::kSyntheticExtentSecond),
        detail::kSyntheticExtentSecond, h);
    Scalar acc = scalar_cast<Scalar>(0.0);
    for (int i = 0; i < s.count; ++i) {
        acc = acc + scalar_cast<Scalar>(s.coeffs[i]) *
                        f(t + static_cast<double>(s.offsets[i]) * h);
    }
    return acc;
}

// Mixed d/dt of a whole spatial field: differentiates in time first, then hands
// back a field the spatial operators can act on.
//
// Ordering matters for cost, not for the answer: the temporal operator has
// fixed coefficients and the spatial one is linear, so they commute exactly
// (see the commutativity case in the A2T suite). Doing time first means the
// spatial pass runs once instead of once per temporal tap.
template <typename Scalar, typename SampleField>
ScalarField<Scalar> time_derivative_field(
    SampleField&& sample, const Grid& g, double t, double h,
    TimeStencilPosition pos = TimeStencilPosition::Interior) {
    const Stencil s = first_derivative_stencil(
        detail::stencil_index(pos, detail::kSyntheticExtentFirst),
        detail::kSyntheticExtentFirst, h);

    ScalarField<Scalar> out(g);
    for (std::size_t p = 0; p < g.num_points(); ++p)
        out.at(p) = scalar_cast<Scalar>(0.0);

    for (int i = 0; i < s.count; ++i) {
        const ScalarField<Scalar> tap =
            sample(t + static_cast<double>(s.offsets[i]) * h);
        const Scalar c = scalar_cast<Scalar>(s.coeffs[i]);
        for (std::size_t p = 0; p < g.num_points(); ++p)
            out.at(p) = out.at(p) + c * tap.at(p);
    }
    return out;
}

// d_0 g_{mu nu} as a Rank3Field slice, for feeding spacetime_gradient().
//
// `sample` returns the full rank-2 metric field at a given time. The result
// carries data only in the alpha = 0 slot; the spatial slots are left zero and
// are filled by spacetime_gradient itself.
template <typename Scalar, typename SampleMetric>
Rank3Field<Scalar> metric_time_derivative(
    SampleMetric&& sample, const Grid& g, double t, double h,
    TimeStencilPosition pos = TimeStencilPosition::Interior) {
    const Stencil s = first_derivative_stencil(
        detail::stencil_index(pos, detail::kSyntheticExtentFirst),
        detail::kSyntheticExtentFirst, h);

    Rank3Field<Scalar> out(g);
    for (int i = 0; i < s.count; ++i) {
        const Rank2Field<Scalar> tap =
            sample(t + static_cast<double>(s.offsets[i]) * h);
        const Scalar c = scalar_cast<Scalar>(s.coeffs[i]);
        for (std::size_t p = 0; p < g.num_points(); ++p)
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu)
                    out.at(p, 0, mu, nu) =
                        out.at(p, 0, mu, nu) + c * tap.at(p, mu, nu);
    }
    return out;
}

}  // namespace geometry
}  // namespace warp_hypad
