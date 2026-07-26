#pragma once

#include <cassert>

#include "warp_hypad/sensitivity/otilib_scalar.hpp"

// ---------------------------------------------------------------------------
// Parameter seeding: turning a plain value into an OTI number that carries a
// derivative direction.
//
// Seeding is the ONLY place in the project where a parameter acquires
// sensitivity. Everything downstream is ordinary arithmetic on the generic
// kernels. That is the whole point of the design: there are no sensitivity
// equations to derive or maintain.
//
// Direction indices are 0-based and fixed by ParameterIndex below, so that the
// meaning of "direction 1" is defined in one place rather than at each call
// site.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace sensitivity {

// The warp metric parameter set. The order here defines the derivative channel
// layout for every gradient and Hessian the project produces, and is recorded
// in run metadata so a saved sensitivity field can be interpreted later.
enum ParameterIndex : int {
    kParamVelocity = 0,  // bubble velocity v
    kParamRadius = 1,    // bubble radius R
    kParamSigma = 2,     // wall steepness sigma
    kParamCount = 3
};

const char* parameter_name(int direction);

// A constant: real value, no derivative in any direction.
inline OTIScalar make_constant(double value) { return OTIScalar(value); }

// A seeded parameter: real value, unit first-order coefficient in `direction`.
//
// `derivative_order` is accepted for interface symmetry with the project spec
// and must be 1 -- seeding a second-order direction directly is meaningful in
// OTI algebra but is not how gradients/Hessians are obtained here. The Hessian
// falls out of order-2 propagation from first-order seeds; seeding order 2
// would compute something else entirely.
inline OTIScalar seed_parameter(double value, int direction,
                                int derivative_order = 1) {
    assert(direction >= 0 && direction < kOtiDirections);
    assert(derivative_order == 1 &&
           "seed first-order directions only; the Hessian comes from order-2 "
           "propagation, not from seeding order 2");
    (void)derivative_order;

    OTIScalar s(value);
    // set_im_r writes the coefficient of the given imaginary direction at the
    // given order. Using the library accessor rather than touching s.raw().e1
    // keeps this correct if the coefficient layout ever changes.
    //
    // The index is 0-BASED: onumm3n2_set_im_r delegates to set_item, whose
    // order-1 switch maps idx 0,1,2 -> e1,e2,e3 (external/otilib
    // src/c/static/onumm3n2/scalar/base.c). Passing direction+1 here would
    // silently seed the wrong parameter for directions 0 and 1 and write
    // nothing at all for direction 2.
    onumm3n2_set_im_r(1.0, static_cast<imdir_t>(direction), 1, &s.raw());
    return s;
}

// Seed all three warp parameters at once. OTI propagates every direction
// simultaneously in a single evaluation -- this is the property that makes the
// gradient cost independent of the parameter count, unlike finite differences.
struct SeededParameters {
    OTIScalar velocity;
    OTIScalar radius;
    OTIScalar sigma;
};

inline SeededParameters seed_all(double velocity, double radius, double sigma) {
    SeededParameters p;
    p.velocity = seed_parameter(velocity, kParamVelocity);
    p.radius = seed_parameter(radius, kParamRadius);
    p.sigma = seed_parameter(sigma, kParamSigma);
    return p;
}

}  // namespace sensitivity
}  // namespace warp_hypad
