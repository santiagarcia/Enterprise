#pragma once

#include <vector>

#include "program_b/inverse_design/complex_scalar.hpp"

// ---------------------------------------------------------------------------
// B3: frequency-dependent constitutive models.
//
// ===========================================================================
// THE FOURIER CONVENTION, FIXED HERE AND NOWHERE ELSE
// ===========================================================================
//
//     f(t) = Re[ f-hat(omega) e^{-i omega t} ]                      e^{-i w t}
//
// Every sign below follows from that one choice, and none of them is free once
// it is made. Nothing in B3 may test passivity before this section is read.
//
//   TIME DERIVATIVE          d/dt  ->  -i omega
//
//   CONSTITUTIVE RELATION    the project's own (E,B) form, carried unchanged
//                            from B1C:
//                                D = eps E - zeta B
//                                H = zeta_adj E + [mu^-1] B
//                            with zeta_adj = zeta^T, the UNCONJUGATED transpose.
//
//   FIELD ENERGY / POWER     Poynting gives, for time-averaged dissipation,
//                                <p> = (omega/2) Im[ E* . D + H* . B ]
//                            because <E . dD/dt> = (omega/2) Im[E-hat* . D-hat]
//                            under e^{-i omega t}.
//
//   POYNTING VECTOR          <S> = (1/2) Re[ E-hat x H-hat* ]
//
//   SIGN OF Im eps           PASSIVE media have Im eps >= 0 for omega > 0.
//                            For a scalar, <p> = (omega/2) Im[eps] |E|^2, so
//                            Im eps < 0 would be gain. Under the opposite
//                            convention e^{+i omega t} every one of these signs
//                            flips, which is why the convention is recorded
//                            rather than assumed.
//
//   CAUSALITY                the response must be analytic in the UPPER half
//                            omega-plane, so every pole sits in the LOWER half.
//
//   REALITY                  a real time-domain response requires
//                            eps(-omega) = conj(eps(omega)) for real omega.
//
// The model forms the milestone specifies are consistent with this convention
// and with no other. A Lorentz term F/(w0^2 - w^2 - i gamma w) has
// Im = F gamma w / |denominator|^2 > 0 for F, gamma, w > 0 -- passive. A Drude
// term -wp^2/(w^2 + i gamma w) has Im = +wp^2 gamma w / |denominator|^2 > 0 --
// also passive. Had the intended convention been e^{+i omega t}, both would have
// come out as gain, which is the check that pins it.
//
// ===========================================================================
// CAUSALITY IS STRUCTURAL, NOT NUMERICAL
// ===========================================================================
//
// The poles of F/(w0^2 - w^2 - i gamma w) solve w^2 + i gamma w - w0^2 = 0:
//
//     w = +- sqrt(w0^2 - gamma^2/4) - i gamma/2
//
// so Im w = -gamma/2 < 0 for gamma > 0, in the lower half-plane, for BOTH roots
// and for every parameter value with positive damping. Causality is therefore
// enforced by the PARAMETERIZATION -- gamma is mapped through a strictly
// positive transform -- and not argued from sampled data. Kramers-Kronig appears
// only as a diagnostic over a wide band; it is not the proof.
// ---------------------------------------------------------------------------

namespace program_b {
namespace inverse_design {

// --- Recorded conventions, for the report ------------------------------------

struct FrequencyConventions {
    const char* fourier_kernel = "exp(-i*omega*t)";
    const char* time_derivative = "d/dt -> -i*omega";
    const char* constitutive_relation =
        "D = eps E - zeta B ; H = zeta^T E + [mu^-1] B   (E,B representation)";
    const char* reciprocity = "zeta_adjoint = zeta^T, unconjugated transpose";
    const char* dissipated_power =
        "<p> = (omega/2) Im[ E* . D + H* . B ] = (omega/2) v^dagger L v, "
        "L = (C - C^dagger)/(2i)";
    const char* poynting = "<S> = (1/2) Re[ E x H* ]";
    const char* passive_sign =
        "Im eps >= 0 for omega > 0; more generally L positive semidefinite";
    const char* causal_half_plane =
        "response analytic in the UPPER half omega-plane; poles in the LOWER half";
    const char* reality = "eps(-omega) = conj(eps(omega)) for real omega";
};

// --- Strictly positive and box-bounded parameter transforms ------------------
//
// These are what makes HARD feasibility possible, as opposed to B2's soft
// penalties: a parameter that must be positive is expressed as a strictly
// positive function of an unconstrained one, so no optimizer iterate can ever
// leave the feasible set, however long the step.
//
// Both need only sqrt and contain no comparison on a Scalar.

// pos(u) > 0 for every finite u: the positive root of x^2 - u x - c = 0, which
// exists and is unique because the product of the roots is -c < 0.
//
// WRITTEN IN TWO ALGEBRAICALLY IDENTICAL FORMS, and it has to be. The naive
//     pos = (u + sqrt(u^2 + 4c)) / 2
// loses the answer completely to cancellation for large negative u: at
// u = -1e6 with c = 1e-6, the term 4c = 4e-6 is far below the ULP of
// u^2 = 1e12, so sqrt(u^2 + 4c) rounds to exactly 1e6 and the sum evaluates to
// EXACTLY ZERO. That is not a small error -- it is a total loss of the
// positivity the transform exists to guarantee, and it was caught by the gate
// case that asserts strict positivity at extreme inputs.
//
// The stable pairing uses the root-product identity x+ x- = -c:
//     u >= 0:  pos = (u + s)/2          both terms positive, no cancellation
//     u <  0:  pos = 2c / (s - u)       both terms positive, no cancellation
//
// The branch is on real_part(u) only -- a numerical-conditioning decision, which
// is exactly what the scalar seam permits it for. It introduces NO derivative
// discontinuity: the two expressions are the same analytic function, so the
// value and every derivative agree across the switch, and only the rounding
// differs.
template <typename Scalar>
Scalar positive_transform(const Scalar& u, double c) {
    using std::sqrt;
    const Scalar s = sqrt(u * u + scalar_cast<Scalar>(4.0 * c));
    if (real_part(u) >= 0.0)
        return scalar_cast<Scalar>(0.5) * (u + s);
    return scalar_cast<Scalar>(2.0 * c) / (s - u);
}

// The inverse, so a known physical value can be turned back into a design
// variable: u = pos - c/pos.
inline double positive_transform_inverse(double value, double c) {
    return value - c / value;
}

// Smooth sigmoid onto (0,1), built from sqrt so no exp is required.
template <typename Scalar>
Scalar unit_sigmoid(const Scalar& u) {
    using std::sqrt;
    const Scalar one = scalar_cast<Scalar>(1.0);
    return scalar_cast<Scalar>(0.5) * (one + u / sqrt(u * u + one));
}

template <typename Scalar>
Scalar box_transform(const Scalar& u, double lo, double hi) {
    return scalar_cast<Scalar>(lo) +
           scalar_cast<Scalar>(hi - lo) * unit_sigmoid(u);
}

inline double box_transform_inverse(double value, double lo, double hi) {
    const double s = (value - lo) / (hi - lo);          // in (0,1)
    const double t = 2.0 * s - 1.0;                     // in (-1,1)
    return t / std::sqrt(std::fmax(1.0 - t * t, 1e-300));
}

// --- Model terms --------------------------------------------------------------

// One Lorentz oscillator: F / (w0^2 - w^2 - i gamma w).
template <typename Scalar>
Complex<Scalar> lorentz_term(const Scalar& F, const Scalar& w0,
                             const Scalar& gamma, double omega) {
    const Scalar w = scalar_cast<Scalar>(omega);
    const Complex<Scalar> den(w0 * w0 - w * w, -gamma * w);
    return Complex<Scalar>(F, scalar_cast<Scalar>(0.0)) / den;
}

// One Drude term: -wp^2 / (w^2 + i gamma w).
template <typename Scalar>
Complex<Scalar> drude_term(const Scalar& wp, const Scalar& gamma, double omega) {
    const Scalar w = scalar_cast<Scalar>(omega);
    const Complex<Scalar> den(w * w, gamma * w);
    return Complex<Scalar>(-wp * wp, scalar_cast<Scalar>(0.0)) / den;
}

// Pole locations of a Lorentz term, for the causality diagnostic. Real parts
// only: this is a report about where the poles sit, not a differentiated
// quantity.
struct PolePair {
    double real_part_a = 0.0;
    double imag_part_a = 0.0;
    double real_part_b = 0.0;
    double imag_part_b = 0.0;
    bool both_in_lower_half_plane = false;
    bool overdamped = false;
};

inline PolePair lorentz_poles(double w0, double gamma) {
    PolePair p;
    const double disc = w0 * w0 - 0.25 * gamma * gamma;
    if (disc >= 0.0) {
        const double r = std::sqrt(disc);
        p.real_part_a = r;
        p.real_part_b = -r;
        p.imag_part_a = -0.5 * gamma;
        p.imag_part_b = -0.5 * gamma;
    } else {
        // Overdamped: both roots are purely imaginary, and both must still be
        // in the lower half-plane.
        p.overdamped = true;
        const double r = std::sqrt(-disc);
        p.real_part_a = 0.0;
        p.real_part_b = 0.0;
        p.imag_part_a = -0.5 * gamma + r;
        p.imag_part_b = -0.5 * gamma - r;
    }
    p.both_in_lower_half_plane = (p.imag_part_a < 0.0) && (p.imag_part_b < 0.0);
    return p;
}

// --- The design model ---------------------------------------------------------
//
// eps(w)     = eps_inf + sum_n Lorentz_n(w) [+ optional Drude]
// mu(w)      = mu_inf  + sum_m Lorentz_m(w)
// [mu^-1](w) = 1 / mu(w)
// zeta(w)    = zeta_amp * Lorentz_z(w) * K
//
// WHICH QUANTITY IS PARAMETERIZED DIRECTLY. mu is parameterized and mu^-1 is
// DERIVED by inversion, not the other way round. A Lorentz resonance is a
// resonance of the permeability; parameterizing mu^-1 as a sum of Lorentz terms
// would place poles of mu^-1 -- that is, ZEROS of mu -- at the resonance
// frequencies, which is a different and generally non-passive medium. The
// distinction is invisible for a single weak pole and decisive for a strong one.
//
// K is the fixed real antisymmetric magnetoelectric shape, K_23 = +1,
// K_32 = -1, which is the moving-medium and warp-analogue form established at
// B1C. With zeta = z(w) K, the reciprocity condition zeta_adjoint = zeta^T is
// structural: zeta^T = z(w) K^T = -z(w) K. Independent arbitrary complex
// entries are NOT used, because they would break reciprocity by construction.

inline constexpr int kMaxEpsPoles = 2;
inline constexpr int kMaxMuPoles = 1;
inline constexpr int kPhysicalSlots = 16;

// Fixed slot layout. Pinned by the round-trip case, because a silent reordering
// would leave every other case passing while the reported model meant something
// else.
enum PhysicalSlot {
    kSlotEpsInf = 0,
    kSlotEpsPole0F = 1, kSlotEpsPole0W0 = 2, kSlotEpsPole0Gamma = 3,
    kSlotEpsPole1F = 4, kSlotEpsPole1W0 = 5, kSlotEpsPole1Gamma = 6,
    kSlotDrudeWp = 7, kSlotDrudeGamma = 8,
    kSlotMuInf = 9,
    kSlotMuPole0F = 10, kSlotMuPole0W0 = 11, kSlotMuPole0Gamma = 12,
    kSlotZetaAmp = 13, kSlotZetaW0 = 14, kSlotZetaGamma = 15,
};

enum class ParameterMode {
    Direct,             // q holds physical values; feasibility is soft-penalised
    BoundedTransform,   // q is unconstrained; physical values are strictly
                        // feasible by construction -- HARD feasibility
};

const char* to_string(ParameterMode m);

struct DispersiveModelSpec {
    int eps_poles = 2;
    int mu_poles = 1;
    bool drude = false;
    bool magnetoelectric = true;
    ParameterMode mode = ParameterMode::BoundedTransform;

    // Floors for the strictly positive transforms. Also the stated hard lower
    // bounds when mode == BoundedTransform.
    double positive_floor = 1e-6;

    std::vector<int> active_slots() const {
        std::vector<int> s;
        s.push_back(kSlotEpsInf);
        if (eps_poles > 0) { s.push_back(kSlotEpsPole0F); s.push_back(kSlotEpsPole0W0); s.push_back(kSlotEpsPole0Gamma); }
        if (eps_poles > 1) { s.push_back(kSlotEpsPole1F); s.push_back(kSlotEpsPole1W0); s.push_back(kSlotEpsPole1Gamma); }
        if (drude) { s.push_back(kSlotDrudeWp); s.push_back(kSlotDrudeGamma); }
        s.push_back(kSlotMuInf);
        if (mu_poles > 0) { s.push_back(kSlotMuPole0F); s.push_back(kSlotMuPole0W0); s.push_back(kSlotMuPole0Gamma); }
        if (magnetoelectric) { s.push_back(kSlotZetaAmp); s.push_back(kSlotZetaW0); s.push_back(kSlotZetaGamma); }
        return s;
    }
    int parameter_count() const { return static_cast<int>(active_slots().size()); }
};

// Every slot except the magnetoelectric amplitude must be strictly positive:
// oscillator strengths and the background permittivity for passivity, resonance
// frequencies and dampings for passivity AND causality. The zeta amplitude
// carries the handedness of the magnetoelectric response and must be free to
// take either sign.
inline bool slot_is_strictly_positive(int slot) {
    return slot != kSlotZetaAmp;
}

template <typename Scalar>
struct PhysicalDesignVariables {
    DispersiveModelSpec spec;
    Scalar q[kPhysicalSlots];

    PhysicalDesignVariables() {
        for (int i = 0; i < kPhysicalSlots; ++i) q[i] = scalar_cast<Scalar>(0.0);
    }
    int size() const { return spec.parameter_count(); }
};

// Expand the compact design vector into physical slot values, applying the
// feasibility transform. This is the ONLY place the transform is applied.
template <typename Scalar>
void physical_values(const PhysicalDesignVariables<Scalar>& v, Scalar out[kPhysicalSlots]) {
    for (int i = 0; i < kPhysicalSlots; ++i) out[i] = scalar_cast<Scalar>(0.0);
    const std::vector<int> slots = v.spec.active_slots();
    for (std::size_t k = 0; k < slots.size(); ++k) {
        const int slot = slots[k];
        const Scalar raw = v.q[k];
        if (v.spec.mode == ParameterMode::BoundedTransform &&
            slot_is_strictly_positive(slot)) {
            out[slot] = positive_transform(raw, v.spec.positive_floor);
        } else {
            out[slot] = raw;
        }
    }
    // Inactive magnetic background must still be 1, or mu^-1 divides by zero.
    if (v.spec.mu_poles == 0 && v.spec.mode == ParameterMode::Direct &&
        real_part(out[kSlotMuInf]) == 0.0)
        out[kSlotMuInf] = scalar_cast<Scalar>(1.0);
}

// Seed a design vector from known physical values, inverting the transform so
// that physical_values() returns exactly what was asked for.
template <typename Scalar>
PhysicalDesignVariables<Scalar> design_from_physical(
    const DispersiveModelSpec& spec, const double physical[kPhysicalSlots]) {
    PhysicalDesignVariables<Scalar> v;
    v.spec = spec;
    const std::vector<int> slots = spec.active_slots();
    for (std::size_t k = 0; k < slots.size(); ++k) {
        const int slot = slots[k];
        double raw = physical[slot];
        if (spec.mode == ParameterMode::BoundedTransform &&
            slot_is_strictly_positive(slot))
            raw = positive_transform_inverse(physical[slot], spec.positive_floor);
        v.q[k] = scalar_cast<Scalar>(raw);
    }
    return v;
}

// --- The frequency response --------------------------------------------------

template <typename Scalar>
struct ResponseAtFrequency {
    Complex<Scalar> epsilon;      // scalar multiplier of the identity
    Complex<Scalar> mu;           // parameterized directly
    Complex<Scalar> mu_inverse;   // derived: 1 / mu
    Complex<Scalar> zeta;         // scalar multiplier of K
};

template <typename Scalar>
ResponseAtFrequency<Scalar> evaluate_response(
    const PhysicalDesignVariables<Scalar>& v, double omega) {
    Scalar p[kPhysicalSlots];
    physical_values(v, p);

    ResponseAtFrequency<Scalar> r;
    r.epsilon = Complex<Scalar>(p[kSlotEpsInf], scalar_cast<Scalar>(0.0));
    if (v.spec.eps_poles > 0)
        r.epsilon = r.epsilon + lorentz_term(p[kSlotEpsPole0F], p[kSlotEpsPole0W0],
                                             p[kSlotEpsPole0Gamma], omega);
    if (v.spec.eps_poles > 1)
        r.epsilon = r.epsilon + lorentz_term(p[kSlotEpsPole1F], p[kSlotEpsPole1W0],
                                             p[kSlotEpsPole1Gamma], omega);
    if (v.spec.drude)
        r.epsilon = r.epsilon + drude_term(p[kSlotDrudeWp], p[kSlotDrudeGamma], omega);

    r.mu = Complex<Scalar>(p[kSlotMuInf], scalar_cast<Scalar>(0.0));
    if (v.spec.mu_poles > 0)
        r.mu = r.mu + lorentz_term(p[kSlotMuPole0F], p[kSlotMuPole0W0],
                                   p[kSlotMuPole0Gamma], omega);
    r.mu_inverse = reciprocal(r.mu);

    if (v.spec.magnetoelectric)
        r.zeta = lorentz_term(p[kSlotZetaAmp], p[kSlotZetaW0], p[kSlotZetaGamma],
                              omega);
    return r;
}

// The fixed real magnetoelectric shape K: K_23 = +1, K_32 = -1.
template <typename Scalar>
CMat3<Scalar> magnetoelectric_shape() {
    CMat3<Scalar> K;
    K(1, 2) = complex_from_real<Scalar>(1.0);
    K(2, 1) = complex_from_real<Scalar>(-1.0);
    return K;
}

// The three complex 3x3 blocks at one frequency.
template <typename Scalar>
struct ComplexBlocks {
    CMat3<Scalar> epsilon;
    CMat3<Scalar> mu_inverse;
    CMat3<Scalar> zeta;
    CMat3<Scalar> zeta_adjoint;   // = zeta^T, structural
    CMat3<Scalar> mu;             // retained: the passivity form needs mu itself
};

template <typename Scalar>
ComplexBlocks<Scalar> blocks_at_frequency(const PhysicalDesignVariables<Scalar>& v,
                                          double omega) {
    const ResponseAtFrequency<Scalar> r = evaluate_response(v, omega);
    ComplexBlocks<Scalar> b;
    b.epsilon = CMat3<Scalar>::scaled_identity(r.epsilon);
    b.mu_inverse = CMat3<Scalar>::scaled_identity(r.mu_inverse);
    b.mu = CMat3<Scalar>::scaled_identity(r.mu);
    b.zeta = scale(magnetoelectric_shape<Scalar>(), r.zeta);
    b.zeta_adjoint = transpose(b.zeta);
    return b;
}

}  // namespace inverse_design
}  // namespace program_b
