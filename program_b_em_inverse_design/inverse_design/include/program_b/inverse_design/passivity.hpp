#pragma once

#include <algorithm>

#include "program_b/inverse_design/dispersion.hpp"

// ---------------------------------------------------------------------------
// Passivity, derived rather than assumed.
//
// ===========================================================================
// THE CONDITION
// ===========================================================================
//
// Under the e^{-i omega t} convention fixed in dispersion.hpp, the time-averaged
// dissipated power density is
//
//     <p> = (omega/2) Im[ E-hat* . D-hat + H-hat* . B-hat ]
//
// because <E . dD/dt> = (1/2) Re[E-hat* . (-i omega D-hat)]
//                     = (omega/2) Im[E-hat* . D-hat].
//
// Collect the fields as v = (E-hat, H-hat) and the response as
//
//     ( D-hat )   ( eps_EH   xi_EH )  ( E-hat )
//     (       ) = (                )  (       )   =:  C v
//     ( B-hat )   ( zeta_EH  mu_EH )  ( H-hat )
//
// Then E* . D + H* . B = v^dagger C v, and
//
//     <p> = (omega/2) Im[ v^dagger C v ] = (omega/2) v^dagger L v,
//
//     L := (C - C^dagger) / (2i)          Hermitian by construction
//
// PASSIVITY IS  L POSITIVE SEMIDEFINITE, for every omega > 0.
//
// NOT "all imaginary parts are positive". For a scalar eps that reduces to
// Im eps >= 0, but for a coupled bianisotropic medium L has off-diagonal blocks
//
//     L_12 = (xi_EH - zeta_EH^dagger) / (2i)
//
// built from the magnetoelectric response, and a medium with Im eps > 0 and
// Im mu > 0 can still be ACTIVE if that coupling is strong enough to push an
// eigenvalue of L negative. The gate exercises exactly that case.
//
// ===========================================================================
// THE (E,B) TO (E,H) CONVERSION -- this project's exact block convention
// ===========================================================================
//
// B1C fixed the constitutive relation in the (E,B) representation:
//
//     D = eps E - zeta B
//     H = zeta^T E + [mu^-1] B
//
// Solving the second for B and substituting:
//
//     B = mu (H - zeta^T E),          mu := [mu^-1]^{-1}
//     D = (eps + zeta mu zeta^T) E - zeta mu H
//
// so
//
//     eps_EH  = eps + zeta mu zeta^T
//     xi_EH   = - zeta mu
//     zeta_EH = - mu zeta^T
//     mu_EH   = mu
//
// RECIPROCITY IN THIS CONVENTION. Substituting the above,
//
//     xi_EH^T = (-zeta mu)^T = -mu^T zeta^T = -mu zeta^T = zeta_EH
//
// so reciprocity reads  xi_EH = zeta_EH^T,  the UNCONJUGATED transpose. Many
// references write xi = -zeta^T instead; that sign traces to the opposite sign
// convention for zeta in the (E,B) relation, and this project's sign was DERIVED
// at B1C and cross-checked against Minkowski's moving-medium relations. The
// relation is verified numerically rather than imported.
// ---------------------------------------------------------------------------

namespace program_b {
namespace inverse_design {

// The 6x6 constitutive matrix in the (E,H) representation.
template <typename Scalar>
struct ConstitutiveMatrix6 {
    Complex<Scalar> c[36];

    Complex<Scalar>& operator()(int r, int col) { return c[6 * r + col]; }
    const Complex<Scalar>& operator()(int r, int col) const { return c[6 * r + col]; }
};

template <typename Scalar>
struct EHBlocks {
    CMat3<Scalar> epsilon_EH;
    CMat3<Scalar> xi_EH;
    CMat3<Scalar> zeta_EH;
    CMat3<Scalar> mu_EH;
};

template <typename Scalar>
EHBlocks<Scalar> to_EH_representation(const ComplexBlocks<Scalar>& b) {
    EHBlocks<Scalar> o;
    const CMat3<Scalar> zeta_T = transpose(b.zeta);
    o.mu_EH = b.mu;
    o.epsilon_EH = b.epsilon + (b.zeta * b.mu) * zeta_T;
    o.xi_EH = CMat3<Scalar>::zero() - (b.zeta * b.mu);
    o.zeta_EH = CMat3<Scalar>::zero() - (b.mu * zeta_T);
    return o;
}

template <typename Scalar>
ConstitutiveMatrix6<Scalar> assemble_C(const EHBlocks<Scalar>& b) {
    ConstitutiveMatrix6<Scalar> C;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            C(i, j) = b.epsilon_EH(i, j);
            C(i, j + 3) = b.xi_EH(i, j);
            C(i + 3, j) = b.zeta_EH(i, j);
            C(i + 3, j + 3) = b.mu_EH(i, j);
        }
    return C;
}

// L = (C - C^dagger) / (2i). Hermitian by construction, which the gate checks
// rather than trusts.
//
//   (C - C^dagger)_{jk} = (C_jk - conj(C_kj))
//   dividing by 2i:  Re L = Im[(C_jk - conj(C_kj))]/2, Im L = -Re[...]/2
template <typename Scalar>
ConstitutiveMatrix6<Scalar> loss_matrix(const ConstitutiveMatrix6<Scalar>& C) {
    ConstitutiveMatrix6<Scalar> L;
    const Scalar half = scalar_cast<Scalar>(0.5);
    for (int j = 0; j < 6; ++j)
        for (int k = 0; k < 6; ++k) {
            // d = C_jk - conj(C_kj)
            const Scalar dre = C(j, k).re - C(k, j).re;
            const Scalar dim = C(j, k).im + C(k, j).im;
            // d / (2i) = (dim - i dre) / 2
            L(j, k) = Complex<Scalar>(half * dim, -half * dre);
        }
    return L;
}

// --- Hermitian eigenvalues, via the real 2n x 2n embedding -------------------
//
// A Hermitian H = A + iB with A symmetric and B antisymmetric is similar to the
// real symmetric matrix
//
//     M = [  A  -B ]
//         [  B   A ]
//
// whose eigenvalues are exactly those of H, each appearing twice. Reusing a real
// symmetric solver this way avoids writing a complex Jacobi rotation, which is
// the step most likely to be wrong.
//
// Real parts only: eigenvalues are a diagnostic and a feasibility verdict, both
// of which are comparisons, and comparisons never belong inside a differentiated
// expression.

// Cyclic Jacobi for a real symmetric n x n, n <= 12. Ascending eigenvalues.
inline void symmetric_eigenvalues_n(const double* a_in, int n, double* out) {
    std::vector<double> a(static_cast<std::size_t>(n) * n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            a[static_cast<std::size_t>(i) * n + j] =
                0.5 * (a_in[static_cast<std::size_t>(i) * n + j] +
                       a_in[static_cast<std::size_t>(j) * n + i]);

    auto A = [&](int i, int j) -> double& {
        return a[static_cast<std::size_t>(i) * n + j];
    };

    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0, scale = 0.0;
        for (int i = 0; i < n; ++i) {
            scale = std::fmax(scale, std::fabs(A(i, i)));
            for (int j = i + 1; j < n; ++j) off = std::fmax(off, std::fabs(A(i, j)));
        }
        if (off <= 1e-16 * std::fmax(1.0, scale)) break;

        for (int p = 0; p < n - 1; ++p)
            for (int q = p + 1; q < n; ++q) {
                if (std::fabs(A(p, q)) <= 1e-300) continue;
                const double th = (A(q, q) - A(p, p)) / (2.0 * A(p, q));
                const double t = (th >= 0.0 ? 1.0 : -1.0) /
                                 (std::fabs(th) + std::sqrt(th * th + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0);
                const double s = t * c;
                for (int k = 0; k < n; ++k) {
                    const double akp = A(k, p), akq = A(k, q);
                    A(k, p) = c * akp - s * akq;
                    A(k, q) = s * akp + c * akq;
                }
                for (int k = 0; k < n; ++k) {
                    const double apk = A(p, k), aqk = A(q, k);
                    A(p, k) = c * apk - s * aqk;
                    A(q, k) = s * apk + c * aqk;
                }
            }
    }
    for (int i = 0; i < n; ++i) out[i] = A(i, i);
    std::sort(out, out + n);
}

// Eigenvalues of a Hermitian 6x6, returned as six ascending real numbers.
template <typename Scalar>
void hermitian_eigenvalues_6(const ConstitutiveMatrix6<Scalar>& H, double out[6]) {
    double M[144];
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
            const double re = real_part(H(i, j).re);
            const double im = real_part(H(i, j).im);
            M[i * 12 + j] = re;              // A
            M[i * 12 + (j + 6)] = -im;       // -B
            M[(i + 6) * 12 + j] = im;        //  B
            M[(i + 6) * 12 + (j + 6)] = re;  //  A
        }
    double eig12[12];
    symmetric_eigenvalues_n(M, 12, eig12);
    // Each eigenvalue appears twice; take every second one.
    for (int i = 0; i < 6; ++i) out[i] = eig12[2 * i];
}

// --- The report ---------------------------------------------------------------

struct PassivityReport {
    double min_loss_eigenvalue = 0.0;     // over the whole band
    double max_violation = 0.0;           // max(0, -lambda_min)
    double worst_frequency = 0.0;
    double worst_hermiticity_residual = 0.0;
    double worst_reciprocity_residual = 0.0;   // |xi_EH - zeta_EH^T|
    bool passive_across_band = false;
    int samples = 0;
};

// Passivity at one frequency: the smallest eigenvalue of L.
template <typename Scalar>
double passivity_margin_at(const PhysicalDesignVariables<Scalar>& v, double omega,
                           double* hermiticity = nullptr,
                           double* reciprocity = nullptr) {
    const ComplexBlocks<Scalar> b = blocks_at_frequency(v, omega);
    const EHBlocks<Scalar> eh = to_EH_representation(b);
    const ConstitutiveMatrix6<Scalar> C = assemble_C(eh);
    const ConstitutiveMatrix6<Scalar> L = loss_matrix(C);

    if (hermiticity) {
        double w = 0.0;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) {
                w = std::fmax(w, std::fabs(real_part(L(i, j).re) -
                                           real_part(L(j, i).re)));
                w = std::fmax(w, std::fabs(real_part(L(i, j).im) +
                                           real_part(L(j, i).im)));
            }
        *hermiticity = w;
    }
    if (reciprocity)
        *reciprocity = cmat3_max_diff(eh.xi_EH, transpose(eh.zeta_EH));

    double eig[6];
    hermitian_eigenvalues_6(L, eig);
    return eig[0];
}

template <typename Scalar>
PassivityReport assess_passivity(const PhysicalDesignVariables<Scalar>& v,
                                const std::vector<double>& omegas) {
    PassivityReport r;
    r.min_loss_eigenvalue = 1e300;
    r.samples = static_cast<int>(omegas.size());
    for (double w : omegas) {
        double herm = 0.0, recip = 0.0;
        const double lam = passivity_margin_at(v, w, &herm, &recip);
        r.worst_hermiticity_residual = std::fmax(r.worst_hermiticity_residual, herm);
        r.worst_reciprocity_residual = std::fmax(r.worst_reciprocity_residual, recip);
        if (lam < r.min_loss_eigenvalue) {
            r.min_loss_eigenvalue = lam;
            r.worst_frequency = w;
        }
    }
    r.max_violation = std::fmax(0.0, -r.min_loss_eigenvalue);
    r.passive_across_band = r.min_loss_eigenvalue >= -1e-12;
    return r;
}

// --- Kramers-Kronig, a DIAGNOSTIC only ----------------------------------------
//
// Not the causality proof. Causality is enforced structurally by the pole
// placement (see dispersion.hpp); this is a wide-band consistency check whose
// residual is dominated by truncation of the principal-value integral at the
// ends of a finite band, and it is reported as such.
//
//     Re eps(w) - eps_inf = (2/pi) P int_0^inf w' Im eps(w') / (w'^2 - w^2) dw'
//
// SINGULARITY SUBTRACTION, which this needs to be worth anything. Writing
// g(w') = w' Im eps(w'),
//
//     g(w')/(w'^2 - w^2) = [g(w') - g(w)]/(w'^2 - w^2)  +  g(w)/(w'^2 - w^2)
//
// the first term is regular at w' = w -- the singularity is removable -- and the
// second has a closed-form principal value:
//
//     P int_a^b dw'/(w'^2 - w^2) = (1/2w) [ ln|(w'-w)/(w'+w)| ]_a^b
//
// A first attempt simply SKIPPED the samples nearest the pole on a log grid.
// That is wrong and stays wrong however wide the band: principal-value
// cancellation needs the grid to be symmetric about the pole, and a logarithmic
// grid is not, so the error plateaued at 1.5e-2 instead of converging. The
// plateau is exactly how it was caught -- a truncation error must fall when the
// band widens, and this one did not.
template <typename Scalar>
double kramers_kronig_residual(const PhysicalDesignVariables<Scalar>& v,
                              double omega_probe, double w_min, double w_max,
                              int samples, double eps_inf) {
    const double w = omega_probe;
    const ResponseAtFrequency<Scalar> at_probe = evaluate_response(v, w);
    const double g_at_w = w * real_part(at_probe.epsilon.im);

    const double lo = std::log(w_min), hi = std::log(w_max);
    const double dl = (hi - lo) / static_cast<double>(samples - 1);

    double regular = 0.0;
    for (int i = 0; i < samples; ++i) {
        const double wp = std::exp(lo + dl * static_cast<double>(i));
        const double denom = wp * wp - w * w;
        if (std::fabs(denom) < 1e-12 * std::fmax(1.0, wp * wp)) continue;
        const ResponseAtFrequency<Scalar> r = evaluate_response(v, wp);
        const double g = wp * real_part(r.epsilon.im);
        const double integrand = (g - g_at_w) / denom;
        const double weight = (i == 0 || i == samples - 1) ? 0.5 : 1.0;
        regular += weight * integrand * wp * dl;   // dw' = w' d(log w')
    }

    const double analytic =
        g_at_w * (0.5 / w) *
        (std::log(std::fabs((w_max - w) / (w_max + w))) -
         std::log(std::fabs((w_min - w) / (w_min + w))));

    const double predicted =
        eps_inf + (2.0 / 3.14159265358979323846) * (regular + analytic);
    return std::fabs(predicted - real_part(at_probe.epsilon.re));
}

}  // namespace inverse_design
}  // namespace program_b
