#pragma once

#include <string>

#include "warp_hypad/core/mat4.hpp"

// ---------------------------------------------------------------------------
// Laboratory versus effective metric: two different objects that happen to have
// the same storage type.
//
// The whole of Program B rests on a distinction that a bare Mat4 cannot carry:
//
//   g_lab   the PHYSICAL background metric of the laboratory. It is what
//           raises and lowers every index in Schuster & Visser (their NOTATION
//           section is explicit: "Raising and lowering of indices shall always
//           be done employing the physical background metric g"), what
//           normalises a physical observer's four-velocity, and what the
//           Levi-Civita tensor is built from. For the analogue experiment this
//           project is aiming at, g_lab = eta.
//
//   g_eff   the EFFECTIVE metric whose light cones the medium reproduces. It is
//           an output of the design problem, never a background. Nothing is
//           normalised against it, no physical observer moves in it, and the
//           Levi-Civita tensor is NOT built from it.
//
// Passing one where the other is expected is a silent physics error: the code
// still runs, every residual is still finite, and the answer is wrong. B1's
// retracted "no electromagnetic analogue exists" conclusion came from exactly
// this family of mistake -- a decomposition performed against the wrong
// timelike direction. The wrappers below make the substitution a compile error
// instead.
//
// They are deliberately minimal: a named struct with one member and an explicit
// constructor. Explicit, so a raw Mat4 never converts silently; a struct rather
// than a strong typedef, so `.value` marks in the source exactly where the
// distinction is being deliberately dropped.
// ---------------------------------------------------------------------------

namespace program_b {

using warp_hypad::Mat4;

// The physical background metric of the laboratory.
template <typename Scalar>
struct LaboratoryMetric {
    Mat4<Scalar> value;

    LaboratoryMetric() = default;
    explicit LaboratoryMetric(const Mat4<Scalar>& g) : value(g) {}

    // The initial analogue experiment: a flat laboratory in Cartesian
    // coordinates. Non-Minkowski laboratory backgrounds are permitted and are
    // exercised as general tensor tests, but they are a different physical
    // claim and are labelled as such at every call site.
    static LaboratoryMetric minkowski() {
        return LaboratoryMetric(Mat4<Scalar>::minkowski());
    }
};

// The effective metric a medium is to reproduce. Never a background.
template <typename Scalar>
struct EffectiveMetric {
    Mat4<Scalar> value;

    EffectiveMetric() = default;
    explicit EffectiveMetric(const Mat4<Scalar>& g) : value(g) {}
};

template <typename Scalar>
LaboratoryMetric<Scalar> laboratory_metric(const Mat4<Scalar>& g) {
    return LaboratoryMetric<Scalar>(g);
}

template <typename Scalar>
EffectiveMetric<Scalar> effective_metric(const Mat4<Scalar>& g) {
    return EffectiveMetric<Scalar>(g);
}

}  // namespace program_b
