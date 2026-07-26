#pragma once

// ---------------------------------------------------------------------------
// Convenience aggregate header for the OTIlib adapter.
//
// The ScalarTraits<OTIScalar> specialisation itself lives in otilib_scalar.hpp,
// deliberately: a translation unit that constructs an OTIScalar must never be
// able to see the type without also seeing its trait, or a generic kernel could
// silently pick up the primary ScalarTraits template and try
// static_cast<double>(OTIScalar), which would fail to compile in a confusing
// place far from the cause.
//
// Include this header to get the whole adapter:
//
//     #include "warp_hypad/sensitivity/otilib_traits.hpp"
//
//     using warp_hypad::sensitivity::OTIScalar;
//     auto p = seed_all(0.2, 3.0, 2.0);
//     OTIScalar J = evaluate_objective<OTIScalar>(p, grid, settings);
//     auto g = extract_full_gradient(J);
//     auto H = extract_full_hessian(J);
//
// Nothing outside include/warp_hypad/sensitivity/ may include an oti/ header or
// name an onumm3n2 symbol. That is what keeps the geometry code independent of
// which AD library is in use, and what would let a different backend be
// substituted by writing one more directory here.
// ---------------------------------------------------------------------------

#include "warp_hypad/sensitivity/otilib_extract.hpp"
#include "warp_hypad/sensitivity/otilib_math.hpp"
#include "warp_hypad/sensitivity/otilib_scalar.hpp"
#include "warp_hypad/sensitivity/otilib_seed.hpp"
