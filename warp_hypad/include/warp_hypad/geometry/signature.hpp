#pragma once

#include "warp_hypad/core/mat4.hpp"

// ---------------------------------------------------------------------------
// Metric signature classification.
//
// Verifying that a metric is Lorentzian -- one timelike and three spacelike
// directions -- requires the SIGNS of its eigenvalues. Sylvester's criterion is
// not sufficient: it tests positive-definiteness, and it breaks down whenever a
// leading principal minor vanishes, which happens for perfectly valid metrics
// under a coordinate ordering that puts a null direction first.
//
// So this computes actual eigenvalues, by cyclic Jacobi rotation. Jacobi is
// chosen over anything faster because at 4x4 the cost is irrelevant and its
// correctness is obvious by inspection: it drives the off-diagonal norm to zero
// by exact plane rotations, each of which preserves the spectrum.
//
// This operates on DOUBLE, on real parts only, deliberately. Signature is a
// discrete classification -- a count of signs -- and therefore has no useful
// derivative. Computing it in a differentiable type would propagate derivative
// channels through comparisons and branches that cannot represent them.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace geometry {

struct SignatureResult {
    double eigenvalues[4];  // ascending
    int negative = 0;
    int positive = 0;
    int zero = 0;

    // Exactly one negative and three positive eigenvalues, in the (-,+,+,+)
    // convention this project uses.
    bool lorentzian = false;

    // Largest |off-diagonal| remaining after Jacobi swept out; a sanity check
    // that the rotation actually converged.
    double residual_offdiagonal = 0.0;

    bool converged = false;
};

// `relative_zero_tolerance` classifies an eigenvalue as zero when
// |lambda| <= tol * max|lambda|. A metric with a zero eigenvalue is degenerate
// and is neither Lorentzian nor Riemannian.
SignatureResult analyse_signature(const Mat4<double>& g,
                                  double relative_zero_tolerance = 1e-12);

}  // namespace geometry
}  // namespace warp_hypad
