#pragma once

#include "program_b/inverse_design/maxwell.hpp"

// ---------------------------------------------------------------------------
// B4: a one-dimensional frequency-domain layered-medium solver.
//
// SCATTERING MATRICES, NOT TRANSFER MATRICES. A transfer matrix carries
// e^{+i k_- d} for the backward wave, whose magnitude GROWS as e^{|Im k| d} in a
// thick evanescent layer, and the cascade overflows long before the physics goes
// wrong. The scattering formulation uses only the decaying exponentials
// e^{+i k_+ d} and e^{-i k_- d}, both of modulus <= 1, and layers are combined by
// the Redheffer star product. Nothing overflows however thick the stack.
//
// GEOMETRY. Normal incidence along x through a stack of homogeneous slabs. For
// the model family B3 produces -- eps and mu^-1 isotropic, zeta = z K with K the
// fixed antisymmetric shape -- the transverse fields separate, and the y and z
// polarizations obey the SAME dispersion relation. So the problem reduces
// exactly to a two-wave (forward, backward) scalar problem per polarization,
// with no approximation.
//
// THE LAYER DISPERSION, derived from this project's own (E,B) constitutive form.
// With E = E_y y-hat and propagation e^{i k x}, Faraday and Ampere give
// B_z = (k/w) E_y, B_y = -(k/w) E_z, and
//
//     D_y = eps E_y - z B_z,      H_z = z E_y + mu^-1 B_z
//
// Substituting into D_y = (k/w) H_z with kappa = k/w:
//
//     mu^-1 kappa^2 + 2 z kappa - eps = 0
//
//     kappa_+- = ( -z +- sqrt(z^2 + eps mu^-1) ) / mu^-1
//
// At z = 0 this is kappa = +- sqrt(eps/mu^-1) = +- sqrt(eps mu), the familiar
// result. At z != 0 the two roots are NOT symmetric: the forward and backward
// phase velocities differ. That asymmetry is physical -- zeta = z K is the
// moving-medium form, and this is Fresnel drag -- not a sign error.
//
// The wave admittance follows from H_z = (z + mu^-1 kappa) E_y.
// ---------------------------------------------------------------------------

// SCALAR-GENERIC. Everything below is templated on Scalar, so the same solver
// runs at double and at the HyPAD scalar and a physical parameter can be
// differentiated through the ACTUAL stack assembly rather than through a
// hand-written surrogate. Nothing here branches on a Scalar value; the
// comparisons live in the reporting and branch-selection helpers at the bottom,
// which are deliberately outside the differentiated path.
namespace program_b {
namespace inverse_design {

// exp of a complex, built from the Scalar's own exp/cos/sin. Found by ADL, so it
// works for double and for OTIScalar alike.
template <typename Scalar>
Complex<Scalar> cexp(const Complex<Scalar>& z) {
    using std::cos;
    using std::exp;
    using std::sin;
    const Scalar e = exp(z.re);
    return Complex<Scalar>(e * cos(z.im), e * sin(z.im));
}

template <typename Scalar>
struct LayerMode {
    Complex<Scalar> kappa_forward;      // k/omega for the forward wave
    Complex<Scalar> kappa_backward;
    Complex<Scalar> admittance_forward;   // H_z / E_y
    Complex<Scalar> admittance_backward;
};

// The two modes of a homogeneous layer at one frequency.
template <typename Scalar>
LayerMode<Scalar> layer_modes(const Complex<Scalar>& eps,
                              const Complex<Scalar>& mu_inverse,
                              const Complex<Scalar>& zeta) {
    LayerMode<Scalar> m;
    const Complex<Scalar> disc = sqrt(zeta * zeta + eps * mu_inverse);
    m.kappa_forward = (-zeta + disc) / mu_inverse;
    m.kappa_backward = (-zeta - disc) / mu_inverse;
    m.admittance_forward = zeta + mu_inverse * m.kappa_forward;
    m.admittance_backward = zeta + mu_inverse * m.kappa_backward;
    return m;
}

// A 2x2 scattering matrix with scalar blocks, in the basis
//   in  = (a from the left, b from the right)
//   out = (b to the left,   a to the right)
template <typename Scalar>
struct SMatrix {
    Complex<Scalar> s11;   // reflection from the left
    Complex<Scalar> s12;   // transmission right -> left
    Complex<Scalar> s21;   // transmission left -> right
    Complex<Scalar> s22;   // reflection from the right

    static SMatrix identity() {
        SMatrix s;
        s.s21 = complex_from_real<Scalar>(1.0);
        s.s12 = complex_from_real<Scalar>(1.0);
        return s;
    }
};

// Redheffer star product: cascade A (ports 1,2) with B (ports 2,3).
template <typename Scalar>
SMatrix<Scalar> star(const SMatrix<Scalar>& A, const SMatrix<Scalar>& B) {
    const Complex<Scalar> one = complex_from_real<Scalar>(1.0);
    const Complex<Scalar> d1 = one - B.s11 * A.s22;
    const Complex<Scalar> d2 = one - A.s22 * B.s11;
    SMatrix<Scalar> S;
    S.s11 = A.s11 + A.s12 * (B.s11 / d1) * A.s21;
    S.s12 = A.s12 * (B.s12 / d1);
    S.s21 = B.s21 * (A.s21 / d2);
    S.s22 = B.s22 + B.s21 * (A.s22 / d2) * B.s12;
    return S;
}

// Interface between medium A (left) and medium B (right). Continuity of E_y and
// H_z gives a 2x2 system for the outgoing amplitudes.
template <typename Scalar>
SMatrix<Scalar> interface_smatrix(const LayerMode<Scalar>& A,
                                  const LayerMode<Scalar>& B) {
    using C = Complex<Scalar>;
    // Unknowns b_A (reflected left) and a_B (transmitted right):
    //   b_A - a_B                     = b_B - a_A
    //   Y-_A b_A - Y+_B a_B           = Y-_B b_B - Y+_A a_A
    const C det = A.admittance_backward - B.admittance_forward;
    const C one = complex_from_real<Scalar>(1.0);
    const C minus_one = complex_from_real<Scalar>(-1.0);
    const C minus_Yf_B = -B.admittance_forward;

    // Column for a_A = 1, b_B = 0:  rhs = (-1, -Y+_A)
    const C r1a = minus_one;
    const C r2a = -A.admittance_forward;
    // Cramer on [[1,-1],[Y-_A, -Y+_B]]
    const C bA_from_aA = (r1a * minus_Yf_B - r2a * minus_one) / det;
    const C aB_from_aA = (one * r2a - A.admittance_backward * r1a) / det;

    // Column for a_A = 0, b_B = 1:  rhs = (1, Y-_B)
    const C r1b = one;
    const C r2b = B.admittance_backward;
    const C bA_from_bB = (r1b * minus_Yf_B - r2b * minus_one) / det;
    const C aB_from_bB = (one * r2b - A.admittance_backward * r1b) / det;

    SMatrix<Scalar> S;
    S.s11 = bA_from_aA;
    S.s12 = bA_from_bB;
    S.s21 = aB_from_aA;
    S.s22 = aB_from_bB;
    return S;
}

// Propagation through a slab of thickness d. Only decaying exponentials appear.
template <typename Scalar>
SMatrix<Scalar> propagation_smatrix(const LayerMode<Scalar>& m,
                                    const Scalar& omega, const Scalar& d) {
    const Complex<Scalar> i = complex_i<Scalar>();
    const Complex<Scalar> wd(omega * d, scalar_cast<Scalar>(0.0));
    const Complex<Scalar> ikp = i * m.kappa_forward * wd;
    const Complex<Scalar> imk = -(i * m.kappa_backward * wd);
    SMatrix<Scalar> S;
    S.s21 = cexp(ikp);   // forward across the slab
    S.s12 = cexp(imk);   // backward across the slab
    return S;
}

template <typename Scalar>
struct Layer {
    Complex<Scalar> epsilon;
    Complex<Scalar> mu_inverse;
    Complex<Scalar> zeta;
    Scalar thickness = scalar_cast<Scalar>(0.0);
};

template <typename Scalar>
struct StackResult {
    double reflection = 0.0;
    double transmission = 0.0;
    double absorption_from_balance = 0.0;    // 1 - R - T
    double absorption_from_dissipation = 0.0;  // independently integrated
    double energy_balance_residual = 0.0;    // |R + T + A_dissipation - 1|
    Scalar phase_local = scalar_cast<Scalar>(0.0);
    double phase = 0.0;                      // arg(t), a REPORTING quantity
    Complex<Scalar> r, t;
    bool ok = false;
};

// Solve a stack embedded between vacuum half-spaces.
//
// SCALAR-GENERIC, so a physical parameter can be seeded and differentiated
// through the ACTUAL assembly -- interfaces, propagation, Redheffer cascade --
// rather than through a hand-written surrogate.
//
// ABSORPTION IS COMPUTED TWICE, on purpose. Once as 1 - R - T, which is a
// definition and proves nothing, and once by integrating the local dissipated
// power <p> = (w/2) Im[E* D + H* B] through the internal field profile. Their
// agreement is the genuine energy-balance check.
//
// The dissipation integral and the R/T/A magnitudes are REPORTING quantities and
// are reduced to double. The differentiable output is the transmission amplitude
// t, kept as Complex<Scalar>.
template <typename Scalar>
StackResult<Scalar> solve_stack(const std::vector<Layer<Scalar>>& layers,
                                const Scalar& omega, int samples_per_layer = 400) {
    using C = Complex<Scalar>;
    StackResult<Scalar> out;
    if (real_part(omega) <= 0.0) return out;

    const LayerMode<Scalar> vac =
        layer_modes(complex_from_real<Scalar>(1.0), complex_from_real<Scalar>(1.0),
                    complex_from_real<Scalar>(0.0));

    std::vector<LayerMode<Scalar>> modes;
    modes.reserve(layers.size());
    for (const Layer<Scalar>& L : layers)
        modes.push_back(layer_modes(L.epsilon, L.mu_inverse, L.zeta));

    SMatrix<Scalar> S = SMatrix<Scalar>::identity();
    LayerMode<Scalar> left = vac;
    for (std::size_t i = 0; i < layers.size(); ++i) {
        S = star(S, interface_smatrix(left, modes[i]));
        S = star(S, propagation_smatrix(modes[i], omega, layers[i].thickness));
        left = modes[i];
    }
    S = star(S, interface_smatrix(left, vac));

    out.r = S.s11;
    out.t = S.s21;
    out.reflection = real_part(out.r.re) * real_part(out.r.re) +
                     real_part(out.r.im) * real_part(out.r.im);
    out.transmission = real_part(out.t.re) * real_part(out.t.re) +
                       real_part(out.t.im) * real_part(out.t.im);
    out.absorption_from_balance = 1.0 - out.reflection - out.transmission;
    out.phase = std::atan2(real_part(out.t.im), real_part(out.t.re));

    // --- Internal fields, and the dissipation integral ----------------------
    C a = complex_from_real<Scalar>(1.0), b = out.r;
    LayerMode<Scalar> cur = vac;
    double absorbed = 0.0;

    for (std::size_t i = 0; i < layers.size(); ++i) {
        const LayerMode<Scalar>& nx = modes[i];
        const C E = a + b;
        const C H = cur.admittance_forward * a + cur.admittance_backward * b;
        const C den = nx.admittance_forward - nx.admittance_backward;
        const C ap = (H - nx.admittance_backward * E) / den;
        const C bp = (nx.admittance_forward * E - H) / den;

        const double d = real_part(layers[i].thickness);
        const double w = real_part(omega);
        const int n = (samples_per_layer < 2) ? 2 : samples_per_layer;
        const C i_unit = complex_i<Scalar>();
        const C shift = cexp(i_unit * nx.kappa_backward *
                             complex_from_real<Scalar>(w * d));
        double slab = 0.0;
        for (int sIdx = 0; sIdx <= n; ++sIdx) {
            const double x = d * static_cast<double>(sIdx) / static_cast<double>(n);
            const C pf = cexp(i_unit * nx.kappa_forward *
                              complex_from_real<Scalar>(w * x));
            const C pb = cexp(i_unit * nx.kappa_backward *
                              complex_from_real<Scalar>(w * (x - d)));
            const C Ey = ap * pf + (bp * shift) * pb;
            const C Hz = nx.admittance_forward * (ap * pf) +
                         nx.admittance_backward * ((bp * shift) * pb);
            const C Bz = (Hz - layers[i].zeta * Ey) / layers[i].mu_inverse;
            const C Dy = layers[i].epsilon * Ey - layers[i].zeta * Bz;
            const double eyr = real_part(Ey.re), eyi = real_part(Ey.im);
            const double dyr = real_part(Dy.re), dyi = real_part(Dy.im);
            const double hzr = real_part(Hz.re), hzi = real_part(Hz.im);
            const double bzr = real_part(Bz.re), bzi = real_part(Bz.im);
            const double pw = 0.5 * w * ((eyr * dyi - eyi * dyr) +
                                         (hzr * bzi - hzi * bzr));
            const double wgt = (sIdx == 0 || sIdx == n) ? 0.5 : 1.0;
            slab += wgt * pw * (d / static_cast<double>(n));
        }
        absorbed += slab;

        const C pf_d = cexp(i_unit * nx.kappa_forward *
                            complex_from_real<Scalar>(w * d));
        a = ap * pf_d;
        b = bp;
        cur = nx;
    }

    out.absorption_from_dissipation = absorbed / 0.5;
    out.energy_balance_residual = std::fabs(
        out.reflection + out.transmission + out.absorption_from_dissipation - 1.0);
    out.ok = true;
    return out;
}

// --- The differentiable phase observable --------------------------------------
//
// arg(t) needs atan2, which is branch sensitive and not something to hand an AD
// scalar. The DERIVATIVE of the phase does not need it:
//
//     d(phase) = Im( dt / t )   for t != 0
//
// so the phase sensitivity comes from the channels of t alone. With
// t = (t_re, t_im),
//
//     d(phase)/dq = ( t_re dt_im/dq - t_im dt_re/dq ) / |t|^2
//
// which involves no inverse trigonometry at all.
//
// This is NOT a claim to differentiate through branch sorting, mode switching,
// poles, zeros of transmission, or degeneracies. Where |t| approaches zero the
// derivative diverges, and the caller is told so rather than handed a number.
struct PhaseSensitivity {
    double d_phase = 0.0;
    double transmission_magnitude = 0.0;
    bool reliable = false;
    std::string note;
};

inline PhaseSensitivity phase_sensitivity(double t_re, double t_im, double dt_re,
                                          double dt_im,
                                          double magnitude_floor = 1e-8) {
    PhaseSensitivity p;
    const double m2 = t_re * t_re + t_im * t_im;
    p.transmission_magnitude = std::sqrt(m2);
    if (p.transmission_magnitude <= magnitude_floor) {
        p.note = "|t| below the floor: the phase derivative is undefined near a "
                 "transmission zero, so no value is reported";
        return p;
    }
    p.d_phase = (t_re * dt_im - t_im * dt_re) / m2;
    p.reliable = true;
    p.note = "d(phase) = Im(dt/t), no inverse trigonometry involved";
    return p;
}

}  // namespace inverse_design
}  // namespace program_b
