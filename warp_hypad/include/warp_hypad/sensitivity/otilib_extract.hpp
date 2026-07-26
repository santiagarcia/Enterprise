#pragma once

#include <cassert>

#include "warp_hypad/sensitivity/otilib_scalar.hpp"
#include "warp_hypad/sensitivity/otilib_seed.hpp"

// ---------------------------------------------------------------------------
// Derivative extraction.
//
// OTI stores TAYLOR COEFFICIENTS, not derivatives. For a repeated direction the
// stored coefficient is 1/k! times the derivative:
//
//     e_ii  = (1/2!) d2f/dpi2        but      e_ij (i != j) = d2f/dpi dpj
//
// Reading the struct fields directly and calling them "the Hessian" would
// therefore be wrong by a factor of 2 on the diagonal only -- a bug that leaves
// gradients perfect, keeps the Hessian symmetric, and is easy to miss.
//
// OTIlib already handles this: onumm3n2_get_deriv() converts a stored
// coefficient into the true derivative. Verified against the analytic Hessian of
// a^2 sin(b) + exp(ac): all six second derivatives matched to machine precision
// through get_deriv, while get_item disagreed on exactly the three diagonal
// entries. Everything below goes through get_deriv for that reason.
//
// Order-2 direction index (0-based), from the library's own switch:
//
//     idx:  0=(0,0) 1=(0,1) 2=(1,1) 3=(0,2) 4=(1,2) 5=(2,2)
//
// which is the column-major upper triangle, idx = j(j+1)/2 + i for i <= j.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace sensitivity {

inline double extract_real(const OTIScalar& value) { return value.real(); }

// Index of the order-2 coefficient for the unordered pair (i, j).
inline int hessian_index(int i, int j) {
    assert(i >= 0 && i < kOtiDirections && j >= 0 && j < kOtiDirections);
    const int lo = (i < j) ? i : j;
    const int hi = (i < j) ? j : i;
    return hi * (hi + 1) / 2 + lo;
}

// df/dp_direction.
inline double extract_gradient(const OTIScalar& value, int direction) {
    assert(direction >= 0 && direction < kOtiDirections);
    return onumm3n2_get_deriv(static_cast<imdir_t>(direction), 1,
                              const_cast<onumm3n2_t*>(&value.raw()));
}

// d2f/dp_i dp_j, with the 1/k! convention already undone by get_deriv.
inline double extract_hessian(const OTIScalar& value, int i, int j) {
    return onumm3n2_get_deriv(static_cast<imdir_t>(hessian_index(i, j)), 2,
                              const_cast<onumm3n2_t*>(&value.raw()));
}

// General form required by the project interface. `derivative_order` selects
// which tier to read; `direction` indexes within that tier (so for order 2 it is
// the packed pair index from hessian_index, not a parameter index).
inline double extract_derivative(const OTIScalar& value, int direction,
                                 int derivative_order) {
    assert(derivative_order >= 0 && derivative_order <= kOtiOrder);
    if (derivative_order == 0) return value.real();
    return onumm3n2_get_deriv(static_cast<imdir_t>(direction),
                              static_cast<ord_t>(derivative_order),
                              const_cast<onumm3n2_t*>(&value.raw()));
}

struct Gradient {
    double d[kOtiDirections];
    double operator[](int i) const { return d[i]; }
};

struct Hessian {
    double h[kOtiDirections][kOtiDirections];
    double operator()(int i, int j) const { return h[i][j]; }
};

inline Gradient extract_full_gradient(const OTIScalar& value) {
    Gradient g{};
    for (int i = 0; i < kOtiDirections; ++i) g.d[i] = extract_gradient(value, i);
    return g;
}

// Symmetric by construction: (i,j) and (j,i) read the same packed coefficient.
// That means Hessian symmetry is NOT evidence the extraction is correct, and
// the verification suite tests it against analytic values instead.
inline Hessian extract_full_hessian(const OTIScalar& value) {
    Hessian h{};
    for (int i = 0; i < kOtiDirections; ++i)
        for (int j = 0; j < kOtiDirections; ++j) h.h[i][j] = extract_hessian(value, i, j);
    return h;
}

}  // namespace sensitivity
}  // namespace warp_hypad
