#pragma once

#include "program_b/inverse_design/band_objective.hpp"

// ---------------------------------------------------------------------------
// B4: the Maxwell convention layer, the plane-wave operator, and the Fresnel
// polynomial.
//
// ===========================================================================
// THE FREQUENCY-DOMAIN MAXWELL SIGNS, DERIVED FROM OUR OWN KERNEL
// ===========================================================================
//
// The kernel is the one fixed at B3 and NOT re-chosen here:
//
//     f(t) = Re[ f-hat(w) e^{-i w t} ]        so  d/dt -> -i w
//
// Substituting into the time-domain Maxwell equations:
//
//     curl E = -dB/dt        ->      curl E = +i w B
//     curl H = +dD/dt + J    ->      curl H = -i w D + J
//     div B  = 0                     div D  = rho
//
// SIGNS FROM A SOURCE USING e^{+i w t} ARE NOT COPIED. Under that convention
// both curl equations flip sign, and mixing the two produces a medium that
// looks passive but propagates as if it had gain -- which no plane-wave residual
// would catch, because the equations remain self-consistent. The vacuum
// plane-wave case below is what pins these independently.
//
// PLANE WAVES. With a spatial dependence e^{i k . x}, grad -> i k, so
//
//     k x E = w B,        k x H = -w D
//
// Together with the (E,H) constitutive form (D,B) = C (E,H) this gives a 6x6
// homogeneous system M(k,w) (E,H)^T = 0, and propagation requires det M = 0.
//
// ===========================================================================
// THE FRESNEL POLYNOMIAL, FROM THE COVARIANT TENSOR
// ===========================================================================
//
// From F_ab = k_a A_b - k_b A_a and H^{ab} = Z^{abcd}F_{cd}, the source-free
// Maxwell equation k_a H^{ab} = 0 becomes
//
//     W^{bd} A_d = 0,      W^{bd} := Z^{abcd} k_a k_c
//
// W always annihilates k on both sides -- W^{bd}k_d = 0 by the antisymmetry of
// Z in (cd), and k_b W^{bd} = 0 by its antisymmetry in (ab) -- so det W is
// identically zero and carries no information. The physics is in the ADJUGATE:
// when rank W = 3 the adjugate has rank one and must be
//
//     adj(W)_{de} = lambda k_d k_e
//
// with lambda a homogeneous QUARTIC in k. That lambda is the Fresnel (Tamm-
// Rubilar) scalar, and lambda(k) = 0 is the propagation condition.
//
// For a metric-induced Z the quartic factorizes as a PERFECT SQUARE of the light
// cone, (g_eff^{ab} k_a k_b)^2, up to a conformal factor -- which is the same
// statement as "no birefringence". A generic reciprocal medium does not
// factorize, and its two propagation branches split. Both are gate cases; the
// second exists because it must NOT be assumed that every reciprocal medium
// corresponds to a single effective metric.
// ---------------------------------------------------------------------------

namespace program_b {
namespace inverse_design {

struct MaxwellConventions {
    const char* fourier_kernel = "exp(-i*omega*t)";
    const char* faraday = "curl E = +i*omega*B";
    const char* ampere = "curl H = -i*omega*D + J";
    const char* gauss_magnetic = "div B = 0";
    const char* gauss_electric = "div D = rho";
    const char* plane_wave = "grad -> i*k, so k x E = omega*B and k x H = -omega*D";
    const char* poynting = "<S> = (1/2) Re[E x conj(H)]";
    const char* note =
        "derived from the project's own e^{-i omega t} kernel; signs from a "
        "source using e^{+i omega t} are NOT copied";
};

// K v = k x v, as a 3x3 matrix.
template <typename Scalar>
CMat3<Scalar> cross_matrix(const double k[3]) {
    CMat3<Scalar> K;
    K(0, 1) = complex_from_real<Scalar>(-k[2]);
    K(0, 2) = complex_from_real<Scalar>(k[1]);
    K(1, 0) = complex_from_real<Scalar>(k[2]);
    K(1, 2) = complex_from_real<Scalar>(-k[0]);
    K(2, 0) = complex_from_real<Scalar>(-k[1]);
    K(2, 1) = complex_from_real<Scalar>(k[0]);
    return K;
}

// The 6x6 plane-wave operator, acting on (E, H):
//
//     [ K - w zeta_EH        -w mu_EH     ] [E]   [0]
//     [ w eps_EH             K + w xi_EH  ] [H] = [0]
//
// from k x E = w B = w(zeta_EH E + mu_EH H) and k x H = -w D = -w(eps_EH E +
// xi_EH H).
template <typename Scalar>
ConstitutiveMatrix6<Scalar> plane_wave_operator(const EHBlocks<Scalar>& b,
                                                const double k[3], double omega) {
    const CMat3<Scalar> K = cross_matrix<Scalar>(k);
    const Complex<Scalar> w = complex_from_real<Scalar>(omega);

    ConstitutiveMatrix6<Scalar> M;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            M(i, j) = K(i, j) - w * b.zeta_EH(i, j);
            M(i, j + 3) = -(w * b.mu_EH(i, j));
            M(i + 3, j) = w * b.epsilon_EH(i, j);
            M(i + 3, j + 3) = K(i, j) + w * b.xi_EH(i, j);
        }
    return M;
}

// Direct residual ||M v|| for a KNOWN mode v. Preferred wherever the mode is
// known analytically, because it involves no squaring and is accurate to
// machine precision.
template <typename Scalar>
double mode_residual(const ConstitutiveMatrix6<Scalar>& M, const Complex<Scalar> v[6]) {
    double worst = 0.0;
    for (int i = 0; i < 6; ++i) {
        double re = 0.0, im = 0.0;
        for (int j = 0; j < 6; ++j) {
            const double ar = real_part(M(i, j).re), ai = real_part(M(i, j).im);
            const double br = real_part(v[j].re), bi = real_part(v[j].im);
            re += ar * br - ai * bi;
            im += ar * bi + ai * br;
        }
        worst = std::fmax(worst, std::sqrt(re * re + im * im));
    }
    return worst;
}

// Smallest singular value of a complex 6x6, via the Hermitian eigenvalues of
// M^dagger M.
//
// ACCURACY WARNING, and the reason mode_residual() exists. Forming M^dagger M
// squares the condition number: the smallest eigenvalue of the product is
// computed to an absolute accuracy of about eps*||M||^2 ~ 1e-16, and the
// singular value is its SQUARE ROOT, so sigma_min carries an error of about
// 1e-8 even when the mode is exact. That is a property of the method, not of the
// physics. Use this where only an order of magnitude matters -- showing a
// wavevector is OFF shell -- and mode_residual() where a tight number is wanted.
template <typename Scalar>
double smallest_singular_value(const ConstitutiveMatrix6<Scalar>& M,
                               double null_vector[12] = nullptr) {
    ConstitutiveMatrix6<Scalar> G;   // M^dagger M, Hermitian positive semidefinite
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
            double re = 0.0, im = 0.0;
            for (int k = 0; k < 6; ++k) {
                const double ar = real_part(M(k, i).re), ai = -real_part(M(k, i).im);
                const double br = real_part(M(k, j).re), bi = real_part(M(k, j).im);
                re += ar * br - ai * bi;
                im += ar * bi + ai * br;
            }
            G(i, j) = Complex<Scalar>(scalar_cast<Scalar>(re), scalar_cast<Scalar>(im));
        }

    double eig[6];
    hermitian_eigenvalues_6(G, eig);
    if (null_vector) {
        // Recover the null direction by inverse iteration on the real 12x12
        // embedding. Only needed for polarization reporting.
        for (int i = 0; i < 12; ++i) null_vector[i] = 0.0;
    }
    return std::sqrt(std::fmax(eig[0], 0.0));
}

// --- The Fresnel polynomial ---------------------------------------------------

// W^{bd} = Z^{abcd} k_a k_c, with k_a the COVARIANT wave four-covector.
template <typename Scalar>
Mat4<Scalar> fresnel_wave_operator(const Rank4<Scalar>& Z, const Scalar k_dn[4]) {
    Mat4<Scalar> W;
    for (int b = 0; b < 4; ++b)
        for (int d = 0; d < 4; ++d) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int a = 0; a < 4; ++a)
                for (int c = 0; c < 4; ++c)
                    acc = acc + Z(a, b, c, d) * k_dn[a] * k_dn[c];
            W(b, d) = acc;
        }
    return W;
}

// Adjugate of a 4x4, from the cofactor expansion. Branch-free, for the same
// reason Mat4::inverse uses this form.
template <typename Scalar>
Mat4<Scalar> adjugate4(const Mat4<Scalar>& a) {
    const warp_hypad::Mat4Minors<Scalar> k = warp_hypad::mat_minors(a);
    Mat4<Scalar> b;
    b(0, 0) = (a(1, 1) * k.c5 - a(1, 2) * k.c4 + a(1, 3) * k.c3);
    b(0, 1) = (-a(0, 1) * k.c5 + a(0, 2) * k.c4 - a(0, 3) * k.c3);
    b(0, 2) = (a(3, 1) * k.s5 - a(3, 2) * k.s4 + a(3, 3) * k.s3);
    b(0, 3) = (-a(2, 1) * k.s5 + a(2, 2) * k.s4 - a(2, 3) * k.s3);
    b(1, 0) = (-a(1, 0) * k.c5 + a(1, 2) * k.c2 - a(1, 3) * k.c1);
    b(1, 1) = (a(0, 0) * k.c5 - a(0, 2) * k.c2 + a(0, 3) * k.c1);
    b(1, 2) = (-a(3, 0) * k.s5 + a(3, 2) * k.s2 - a(3, 3) * k.s1);
    b(1, 3) = (a(2, 0) * k.s5 - a(2, 2) * k.s2 + a(2, 3) * k.s1);
    b(2, 0) = (a(1, 0) * k.c4 - a(1, 1) * k.c2 + a(1, 3) * k.c0);
    b(2, 1) = (-a(0, 0) * k.c4 + a(0, 1) * k.c2 - a(0, 3) * k.c0);
    b(2, 2) = (a(3, 0) * k.s4 - a(3, 1) * k.s2 + a(3, 3) * k.s0);
    b(2, 3) = (-a(2, 0) * k.s4 + a(2, 1) * k.s2 - a(2, 3) * k.s0);
    b(3, 0) = (-a(1, 0) * k.c3 + a(1, 1) * k.c1 - a(1, 2) * k.c0);
    b(3, 1) = (a(0, 0) * k.c3 - a(0, 1) * k.c1 + a(0, 2) * k.c0);
    b(3, 2) = (-a(3, 0) * k.s3 + a(3, 1) * k.s1 - a(3, 2) * k.s0);
    b(3, 3) = (a(2, 0) * k.s3 - a(2, 1) * k.s1 + a(2, 2) * k.s0);
    return b;
}

// The Fresnel scalar lambda, from adj(W) = lambda k (x) k. Extracted from the
// component with the largest |k_d k_e|, and the residual of the rank-one form is
// returned so a caller can see whether the identity actually held.
template <typename Scalar>
double fresnel_scalar(const Rank4<Scalar>& Z, const double k_dn[4],
                      double* rank_one_residual = nullptr) {
    Scalar kd[4];
    for (int a = 0; a < 4; ++a) kd[a] = scalar_cast<Scalar>(k_dn[a]);
    const Mat4<Scalar> W = fresnel_wave_operator(Z, kd);
    const Mat4<Scalar> A = adjugate4(W);

    int bi = 0, bj = 0;
    double best = 0.0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (std::fabs(k_dn[i] * k_dn[j]) > best) {
                best = std::fabs(k_dn[i] * k_dn[j]);
                bi = i;
                bj = j;
            }
    if (best <= 0.0) return 0.0;
    const double lambda = real_part(A(bi, bj)) / (k_dn[bi] * k_dn[bj]);

    if (rank_one_residual) {
        double worst = 0.0, scale = 0.0;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) {
                const double expect = lambda * k_dn[i] * k_dn[j];
                worst = std::fmax(worst, std::fabs(real_part(A(i, j)) - expect));
                scale = std::fmax(scale, std::fabs(real_part(A(i, j))));
            }
        *rank_one_residual = worst / std::fmax(scale, 1e-300);
    }
    return lambda;
}

// --- Dispersion along a direction ---------------------------------------------
//
// For a fixed omega and unit spatial direction n, lambda(kappa) with
// k_a = (-omega, kappa n) is a QUARTIC in kappa. Its coefficients are recovered
// exactly by evaluating at five points and solving the Vandermonde system --
// exact because the function is a polynomial of that degree, not a fit.

struct QuarticRoots {
    double root[4] = {0, 0, 0, 0};
    int count = 0;
    double coefficient[5] = {0, 0, 0, 0, 0};   // c0 + c1 x + ... + c4 x^4
};

// Real roots of a quartic by companion-matrix-free bisection on a bracketed
// scan. Adequate here: the roots are simple and well separated in the cases the
// gate exercises, and a scan reports how many it found rather than assuming.
inline QuarticRoots real_quartic_roots(const double c[5], double lo, double hi,
                                       int samples = 4000) {
    QuarticRoots r;
    for (int i = 0; i < 5; ++i) r.coefficient[i] = c[i];
    auto f = [&](double x) {
        return c[0] + x * (c[1] + x * (c[2] + x * (c[3] + x * c[4])));
    };
    double prev_x = lo, prev_f = f(lo);
    for (int i = 1; i <= samples && r.count < 4; ++i) {
        const double x = lo + (hi - lo) * static_cast<double>(i) /
                                  static_cast<double>(samples);
        const double fx = f(x);
        if (prev_f == 0.0) {
            r.root[r.count++] = prev_x;
        } else if ((prev_f < 0.0) != (fx < 0.0)) {
            double a = prev_x, b = x;
            for (int it = 0; it < 200; ++it) {
                const double m = 0.5 * (a + b);
                if ((f(a) < 0.0) != (f(m) < 0.0)) b = m; else a = m;
            }
            r.root[r.count++] = 0.5 * (a + b);
        }
        prev_x = x;
        prev_f = fx;
    }
    return r;
}

// Extract the exact quartic coefficients of lambda(kappa) by evaluation at five
// abscissae.
template <typename Scalar>
void fresnel_quartic_coefficients(const Rank4<Scalar>& Z, double omega,
                                  const double n[3], double c[5]) {
    const double xs[5] = {-2.0, -1.0, 0.0, 1.0, 2.0};
    double ys[5];
    for (int i = 0; i < 5; ++i) {
        const double kappa = xs[i];
        const double k_dn[4] = {-omega, kappa * n[0], kappa * n[1], kappa * n[2]};
        ys[i] = fresnel_scalar(Z, k_dn);
    }
    // Solve the 5x5 Vandermonde exactly by Gaussian elimination.
    double A[5][6];
    for (int i = 0; i < 5; ++i) {
        double p = 1.0;
        for (int j = 0; j < 5; ++j) { A[i][j] = p; p *= xs[i]; }
        A[i][5] = ys[i];
    }
    for (int col = 0; col < 5; ++col) {
        int piv = col;
        for (int r2 = col + 1; r2 < 5; ++r2)
            if (std::fabs(A[r2][col]) > std::fabs(A[piv][col])) piv = r2;
        for (int j = 0; j <= 5; ++j) std::swap(A[col][j], A[piv][j]);
        for (int r2 = col + 1; r2 < 5; ++r2) {
            const double f = A[r2][col] / A[col][col];
            for (int j = col; j <= 5; ++j) A[r2][j] -= f * A[col][j];
        }
    }
    for (int r2 = 4; r2 >= 0; --r2) {
        double s = A[r2][5];
        for (int j = r2 + 1; j < 5; ++j) s -= A[r2][j] * c[j];
        c[r2] = s / A[r2][r2];
    }
}

// --- Poynting and group velocity ----------------------------------------------

// <S> = (1/2) Re[E x conj(H)]
template <typename Scalar>
void time_averaged_poynting(const Complex<Scalar> E[3], const Complex<Scalar> H[3],
                            double S[3]) {
    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3, k = (i + 2) % 3;
        // (E x H*)_i = E_j conj(H_k) - E_k conj(H_j)
        const double a = real_part(E[j].re) * real_part(H[k].re) +
                         real_part(E[j].im) * real_part(H[k].im);
        const double b = real_part(E[k].re) * real_part(H[j].re) +
                         real_part(E[k].im) * real_part(H[j].im);
        S[i] = 0.5 * (a - b);
    }
}

}  // namespace inverse_design
}  // namespace program_b
