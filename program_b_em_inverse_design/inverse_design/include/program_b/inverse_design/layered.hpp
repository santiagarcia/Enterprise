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

namespace program_b {
namespace inverse_design {

using CD = Complex<double>;

struct LayerMode {
    CD kappa_forward;      // k/omega for the forward wave
    CD kappa_backward;
    CD admittance_forward;   // H_z / E_y
    CD admittance_backward;
};

// The two modes of a homogeneous layer at one frequency.
inline LayerMode layer_modes(const CD& eps, const CD& mu_inverse, const CD& zeta) {
    LayerMode m;
    const CD disc = sqrt(zeta * zeta + eps * mu_inverse);
    m.kappa_forward = (CD(-zeta.re, -zeta.im) + disc) / mu_inverse;
    m.kappa_backward = (CD(-zeta.re, -zeta.im) - disc) / mu_inverse;
    m.admittance_forward = zeta + mu_inverse * m.kappa_forward;
    m.admittance_backward = zeta + mu_inverse * m.kappa_backward;
    return m;
}

// A 2x2 scattering matrix with scalar blocks, in the basis
//   in  = (a from the left, b from the right)
//   out = (b to the left,   a to the right)
struct SMatrix {
    CD s11;   // reflection from the left
    CD s12;   // transmission right -> left
    CD s21;   // transmission left -> right
    CD s22;   // reflection from the right

    static SMatrix identity() {
        SMatrix s;
        s.s21 = CD(1.0, 0.0);
        s.s12 = CD(1.0, 0.0);
        return s;
    }
};

// Redheffer star product: cascade A (ports 1,2) with B (ports 2,3).
inline SMatrix star(const SMatrix& A, const SMatrix& B) {
    const CD one(1.0, 0.0);
    const CD d1 = one - B.s11 * A.s22;
    const CD d2 = one - A.s22 * B.s11;
    SMatrix S;
    S.s11 = A.s11 + A.s12 * (B.s11 / d1) * A.s21;
    S.s12 = A.s12 * (B.s12 / d1);
    S.s21 = B.s21 * (A.s21 / d2);
    S.s22 = B.s22 + B.s21 * (A.s22 / d2) * B.s12;
    return S;
}

// Interface between medium A (left) and medium B (right). Continuity of E_y and
// H_z gives a 2x2 system for the outgoing amplitudes.
inline SMatrix interface_smatrix(const LayerMode& A, const LayerMode& B) {
    // Unknowns b_A (reflected left) and a_B (transmitted right):
    //   b_A - a_B                     = b_B - a_A
    //   Y-_A b_A - Y+_B a_B           = Y-_B b_B - Y+_A a_A
    const CD det = A.admittance_backward - B.admittance_forward;

    // Column for a_A = 1, b_B = 0:  rhs = (-1, -Y+_A)
    const CD r1a(-1.0, 0.0);
    const CD r2a = CD(-A.admittance_forward.re, -A.admittance_forward.im);
    // Cramer on [[1,-1],[Y-_A, -Y+_B]]
    const CD bA_from_aA = (r1a * CD(-B.admittance_forward.re,
                                    -B.admittance_forward.im) - r2a * CD(-1.0, 0.0)) / det;
    const CD aB_from_aA = (CD(1.0, 0.0) * r2a - A.admittance_backward * r1a) / det;

    // Column for a_A = 0, b_B = 1:  rhs = (1, Y-_B)
    const CD r1b(1.0, 0.0);
    const CD r2b = B.admittance_backward;
    const CD bA_from_bB = (r1b * CD(-B.admittance_forward.re,
                                    -B.admittance_forward.im) - r2b * CD(-1.0, 0.0)) / det;
    const CD aB_from_bB = (CD(1.0, 0.0) * r2b - A.admittance_backward * r1b) / det;

    SMatrix S;
    S.s11 = bA_from_aA;
    S.s12 = bA_from_bB;
    S.s21 = aB_from_aA;
    S.s22 = aB_from_bB;
    return S;
}

// Propagation through a slab of thickness d. Only decaying exponentials appear.
inline SMatrix propagation_smatrix(const LayerMode& m, double omega, double d) {
    const CD ikp = CD(0.0, 1.0) * m.kappa_forward * CD(omega * d, 0.0);
    const CD imk = CD(0.0, -1.0) * m.kappa_backward * CD(omega * d, 0.0);
    auto cexp = [](const CD& z) {
        const double e = std::exp(z.re);
        return CD(e * std::cos(z.im), e * std::sin(z.im));
    };
    SMatrix S;
    S.s21 = cexp(ikp);   // forward across the slab
    S.s12 = cexp(imk);   // backward across the slab
    return S;
}

struct Layer {
    CD epsilon;
    CD mu_inverse;
    CD zeta;
    double thickness = 0.0;
};

struct StackResult {
    double reflection = 0.0;
    double transmission = 0.0;
    double absorption_from_balance = 0.0;    // 1 - R - T
    double absorption_from_dissipation = 0.0;  // independently integrated
    double energy_balance_residual = 0.0;    // |R + T + A_dissipation - 1|
    double phase = 0.0;                      // arg(t)
    CD r, t;
    bool ok = false;
};

// Solve a stack embedded between vacuum half-spaces.
//
// ABSORPTION IS COMPUTED TWICE, on purpose. Once as 1 - R - T, which is a
// definition and proves nothing, and once by integrating the local dissipated
// power <p> = (w/2) Im[E* D + H* B] through the internal field profile. Their
// agreement is the genuine energy-balance check; asserting R + T + (1-R-T) = 1
// would be circular.
inline StackResult solve_stack(const std::vector<Layer>& layers, double omega,
                               int samples_per_layer = 400) {
    StackResult out;
    if (omega <= 0.0) return out;

    const LayerMode vac = layer_modes(CD(1.0, 0.0), CD(1.0, 0.0), CD(0.0, 0.0));

    std::vector<LayerMode> modes;
    modes.reserve(layers.size());
    for (const Layer& L : layers)
        modes.push_back(layer_modes(L.epsilon, L.mu_inverse, L.zeta));

    SMatrix S = SMatrix::identity();
    LayerMode left = vac;
    for (std::size_t i = 0; i < layers.size(); ++i) {
        S = star(S, interface_smatrix(left, modes[i]));
        S = star(S, propagation_smatrix(modes[i], omega, layers[i].thickness));
        left = modes[i];
    }
    S = star(S, interface_smatrix(left, vac));

    out.r = S.s11;
    out.t = S.s21;
    out.reflection = out.r.re * out.r.re + out.r.im * out.r.im;
    out.transmission = out.t.re * out.t.re + out.t.im * out.t.im;
    out.absorption_from_balance = 1.0 - out.reflection - out.transmission;
    out.phase = std::atan2(out.t.im, out.t.re);

    // --- Internal fields, by forward substitution with unit incidence --------
    //
    // Amplitudes are propagated layer by layer using the same interface
    // relations, then the dissipated power is integrated.
    CD a(1.0, 0.0), b = out.r;   // in vacuum on the left
    LayerMode cur = vac;
    double absorbed = 0.0;

    for (std::size_t i = 0; i < layers.size(); ++i) {
        // Match across the interface into layer i.
        const LayerMode& nx = modes[i];
        // E and H continuous:  a + b = a' + b',  Ya a + Yb b = Y'a a' + Y'b b'
        const CD E = a + b;
        const CD H = cur.admittance_forward * a + cur.admittance_backward * b;
        const CD den = nx.admittance_forward - nx.admittance_backward;
        const CD ap = (H - nx.admittance_backward * E) / den;
        const CD bp = (nx.admittance_forward * E - H) / den;

        // Integrate dissipation through the slab.
        const double d = layers[i].thickness;
        const int n = std::max(2, samples_per_layer);
        auto cexp = [](const CD& z) {
            const double e = std::exp(z.re);
            return CD(e * std::cos(z.im), e * std::sin(z.im));
        };
        double slab = 0.0;
        for (int s = 0; s <= n; ++s) {
            const double x = d * static_cast<double>(s) / static_cast<double>(n);
            const CD pf = cexp(CD(0.0, 1.0) * nx.kappa_forward * CD(omega * x, 0.0));
            const CD pb = cexp(CD(0.0, 1.0) * nx.kappa_backward *
                               CD(omega * (x - d), 0.0));
            const CD Ey = ap * pf + (bp * cexp(CD(0.0, 1.0) * nx.kappa_backward *
                                               CD(omega * d, 0.0))) * pb;
            const CD Hz = nx.admittance_forward * (ap * pf) +
                          nx.admittance_backward *
                              ((bp * cexp(CD(0.0, 1.0) * nx.kappa_backward *
                                          CD(omega * d, 0.0))) * pb);
            // B_z = (H_z - z E_y)/mu^-1 ;  D_y = eps E_y - z B_z
            const CD Bz = (Hz - layers[i].zeta * Ey) / layers[i].mu_inverse;
            const CD Dy = layers[i].epsilon * Ey - layers[i].zeta * Bz;
            // <p> = (w/2) Im[ conj(E) D + conj(H) B ]
            const double p = 0.5 * omega *
                             ((Ey.re * Dy.im - Ey.im * Dy.re) +
                              (Hz.re * Bz.im - Hz.im * Bz.re));
            const double wgt = (s == 0 || s == n) ? 0.5 : 1.0;
            slab += wgt * p * (d / static_cast<double>(n));
        }
        absorbed += slab;

        // Advance to the far side of the slab.
        const CD pf_d = cexp(CD(0.0, 1.0) * nx.kappa_forward * CD(omega * d, 0.0));
        a = ap * pf_d;
        b = bp;
        cur = nx;
    }

    // Incident power in vacuum with unit E amplitude: <S> = (1/2)Re(E H*) = 1/2.
    out.absorption_from_dissipation = absorbed / 0.5;
    out.energy_balance_residual = std::fabs(
        out.reflection + out.transmission + out.absorption_from_dissipation - 1.0);
    out.ok = true;
    return out;
}

}  // namespace inverse_design
}  // namespace program_b
