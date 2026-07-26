#pragma once

#include "warp_hypad/sensitivity/otilib_scalar.hpp"

// ---------------------------------------------------------------------------
// Elementary functions for OTIScalar, found by ADL from generic kernels.
//
// The set here mirrors exactly the set TestScalar provides, so a kernel that
// compiles against TestScalar compiles against OTIScalar. Every one is smooth
// on its domain; nothing here branches on a value.
//
// abs() is deliberately NOT provided even though OTIlib exposes
// onumm3n2_abs: |x| is not differentiable at zero, and providing it would let a
// kernel silently introduce a derivative discontinuity. Code that genuinely
// needs a magnitude uses real_part() and says so.
//
// Aliasing: as in otilib_scalar.hpp, every call writes to a fresh result.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace sensitivity {

#define WARP_HYPAD_OTI_UNARY(NAME, CFUNC)                       \
    inline OTIScalar NAME(const OTIScalar& x) {                 \
        OTIScalar r;                                            \
        CFUNC(const_cast<onumm3n2_t*>(&x.raw()), &r.raw());     \
        return r;                                               \
    }

WARP_HYPAD_OTI_UNARY(sqrt, onumm3n2_sqrt_to)
WARP_HYPAD_OTI_UNARY(exp, onumm3n2_exp_to)
WARP_HYPAD_OTI_UNARY(log, onumm3n2_log_to)
WARP_HYPAD_OTI_UNARY(sin, onumm3n2_sin_to)
WARP_HYPAD_OTI_UNARY(cos, onumm3n2_cos_to)
WARP_HYPAD_OTI_UNARY(tan, onumm3n2_tan_to)
WARP_HYPAD_OTI_UNARY(sinh, onumm3n2_sinh_to)
WARP_HYPAD_OTI_UNARY(cosh, onumm3n2_cosh_to)
WARP_HYPAD_OTI_UNARY(tanh, onumm3n2_tanh_to)

#undef WARP_HYPAD_OTI_UNARY

// pow with a real exponent. OTIlib exposes only the by-value form for this one.
inline OTIScalar pow(const OTIScalar& x, double e) {
    OTIScalar r;
    r.raw() = onumm3n2_pow(const_cast<onumm3n2_t*>(&x.raw()), e);
    return r;
}

}  // namespace sensitivity
}  // namespace warp_hypad
