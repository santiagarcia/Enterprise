#pragma once

// ---------------------------------------------------------------------------
// C++ scalar wrapper around OTIlib's onumm3n2 static OTI type.
//
// onumm3n2_t is an Order Truncated Imaginary number with 3 imaginary directions
// and truncation order 2:
//
//     x = r + e1 x1 + e2 x2 + e3 x3
//           + e11 x11 + e12 x12 + e22 x22 + e13 x13 + e23 x23 + e33 x33
//
// Three directions is exactly the warp parameter set (velocity, radius, sigma);
// order 2 gives the full gradient AND Hessian of any quantity computed through
// it, in a single evaluation.
//
// Why a wrapper rather than using onumm3n2_t directly:
//
//   OTIlib's C API is out-parameter style -- onumm3n2_mul_oo_to(&a,&b,&res) --
//   which cannot be handed to a kernel written as `a*b`. The wrapper supplies
//   operator overloads so that every kernel templated on Scalar in this project
//   compiles unchanged with OTIScalar. Nothing else in the geometry code learns
//   an OTIlib name.
//
// THE ALIASING RULE (the single most important property of this file):
//
//   OTIlib's *_to functions do NOT support output aliasing. Calling
//   onumm3n2_mul_oo_to(&x, &y, &x) does not crash -- it silently returns WRONG
//   derivatives while leaving the real part correct, which is the worst
//   possible failure mode. This was observed directly during integration: an
//   aliased evaluation of a^2 sin(b) + exp(ac) returned df/db = a^2 cos(b)sin(b)
//   instead of a^2 cos(b).
//
//   Every operator below therefore writes into a fresh local `res` and returns
//   it by value. Compound assignment is expressed as `*this = *this OP rhs`,
//   never as an in-place _to call. No member of this class may ever pass its own
//   storage as both an input and the output of an OTIlib call.
//
// Licence: OTIlib is GPL-3.0 (external/otilib/LICENSE). Linking warp_hypad
// against it makes the combined work GPL-3.0. See
// docs/otilib_discovery_report.md -- this is an unresolved project decision and
// the reason this adapter is an opt-in CMake target rather than part of the core
// library.
// ---------------------------------------------------------------------------

extern "C" {
#include "oti/static/onumm3n2.h"
}

#include "warp_hypad/core/scalar_traits.hpp"

namespace warp_hypad {
namespace sensitivity {

// Number of independent parameters this scalar carries derivatives for.
inline constexpr int kOtiDirections = 3;

// Truncation order: 2 means gradient and Hessian are both available.
inline constexpr int kOtiOrder = 2;

class OTIScalar {
public:
    OTIScalar() { onumm3n2_set_r(0.0, &v_); }

    explicit OTIScalar(double x) { onumm3n2_set_r(x, &v_); }

    // Raw access for the seed/extract helpers. Kernels must not use these.
    const onumm3n2_t& raw() const { return v_; }
    onumm3n2_t& raw() { return v_; }

    double real() const { return v_.r; }

    // --- Arithmetic --------------------------------------------------------
    // Each returns a fresh value. See THE ALIASING RULE above.

    friend OTIScalar operator+(const OTIScalar& a, const OTIScalar& b) {
        OTIScalar r;
        onumm3n2_sum_oo_to(a.ptr(), b.ptr(), &r.v_);
        return r;
    }
    friend OTIScalar operator-(const OTIScalar& a, const OTIScalar& b) {
        OTIScalar r;
        onumm3n2_sub_oo_to(a.ptr(), b.ptr(), &r.v_);
        return r;
    }
    friend OTIScalar operator*(const OTIScalar& a, const OTIScalar& b) {
        OTIScalar r;
        onumm3n2_mul_oo_to(a.ptr(), b.ptr(), &r.v_);
        return r;
    }
    friend OTIScalar operator/(const OTIScalar& a, const OTIScalar& b) {
        OTIScalar r;
        onumm3n2_div_oo_to(a.ptr(), b.ptr(), &r.v_);
        return r;
    }
    friend OTIScalar operator-(const OTIScalar& a) {
        OTIScalar r;
        onumm3n2_neg_to(a.ptr(), &r.v_);
        return r;
    }

    OTIScalar& operator+=(const OTIScalar& o) { *this = *this + o; return *this; }
    OTIScalar& operator-=(const OTIScalar& o) { *this = *this - o; return *this; }
    OTIScalar& operator*=(const OTIScalar& o) { *this = *this * o; return *this; }
    OTIScalar& operator/=(const OTIScalar& o) { *this = *this / o; return *this; }

    // Deliberately absent, exactly as in TestScalar: operator<, operator==, and
    // any conversion to double. A kernel that needs an ordering decision must
    // go through real_part() and say so explicitly.

private:
    // OTIlib takes non-const pointers even for pure inputs. The scalar
    // functions used here do not modify their inputs -- verified against
    // analytic derivatives for every operation in the Milestone-1 equivalent
    // suite -- so const_cast is safe and avoids a copy per operand.
    onumm3n2_t* ptr() const { return const_cast<onumm3n2_t*>(&v_); }

    onumm3n2_t v_;
};

}  // namespace sensitivity

// The trait specialisation lives here so that including this one header is
// enough to make every generic kernel in the project work with OTIScalar.
template <>
struct ScalarTraits<sensitivity::OTIScalar> {
    static double real(const sensitivity::OTIScalar& x) { return x.real(); }
    static sensitivity::OTIScalar from_double(double x) {
        return sensitivity::OTIScalar(x);
    }
    // A constant lifted from a double carries zero derivative in every
    // direction, which is what from_double produces.
    static constexpr bool is_differential = true;
    static constexpr int derivative_channels = sensitivity::kOtiDirections;
};

}  // namespace warp_hypad
