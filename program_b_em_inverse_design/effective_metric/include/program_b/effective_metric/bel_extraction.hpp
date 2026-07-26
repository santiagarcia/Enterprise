#pragma once

#include "program_b/effective_metric/hodge_dual.hpp"

// ---------------------------------------------------------------------------
// B1C steps 2 and 3: the Bel extraction of Schuster & Visser (A7a)-(A7d), and
// the spatial projector.
//
// This REPLACES the guessed decomposition that lived in
// covariant_constitutive.hpp. That one built mu^-1 and zeta from ad-hoc
// three-dimensional Levi-Civita contractions of frame components; the paper
// obtains them by Hodge-dualising Z on one or both index pairs BEFORE
// contracting with the observer. The old code survives only under
// legacy_unverified, which no translation unit includes any more.
//
//     W_eps^{bd}      = V_a V_c   Z^{abcd}        (A7a)
//     W_mu^{bd}       = V_a V_c (*Z*)^{abcd}      (A7b)
//     [W_zeta^T]^{bd} = V_a V_c (*Z )^{abcd}      (A7c)
//     W_zeta^{bd}     = V_a V_c ( Z*)^{abcd}      (A7d)
//
// Physical blocks, with the normalisation DERIVED in
// docs/b1c_source_equations.md (contraction of A7a into A8, confirmed against
// the arbitrary-epsilon worked example B2 -- not inferred from vacuum):
//
//     eps^{bd}    = -2 W_eps^{bd}                 (D1)
//     [mu^-1]^{bd} = +2 W_mu^{bd}                 (D2)   NOTE the opposite sign
//     zeta^{bd}   = c_zeta W_zeta^{bd}            (D4)
//     [zeta+]^{bd} = c_zeta [W_zeta^T]^{bd}       (D3)
//
// c_zeta = +2, DETERMINED at step 5 against (B17). See kZetaNormalisation.
//
// SPATIALITY IS STRUCTURAL, not accidental. Z is antisymmetric in (ab) and in
// (cd), so V_b W^{bd} = V_a V_b V_c Z^{abcd} = 0 because V_aV_b is symmetric and
// Z^{ab..} antisymmetric. The same argument on (cd) gives V_d W^{bd} = 0. The
// tests therefore expect machine-precision zeros, not small residuals.
//
// zeta and zeta+ are stored SEPARATELY. zeta+ = zeta^T holds only under
// reciprocity, which is a property of a particular medium and convention, not
// an identity, and the generic code does not assume it.
// ---------------------------------------------------------------------------

namespace program_b {

// V_a = g_lab_{ab} V^b
template <typename Scalar>
void lower_vector(const Mat4<Scalar>& g_lab, const Scalar V_up[4],
                  Scalar V_dn[4]) {
    for (int a = 0; a < 4; ++a) {
        Scalar acc = scalar_cast<Scalar>(0.0);
        for (int b = 0; b < 4; ++b) acc = acc + g_lab(a, b) * V_up[b];
        V_dn[a] = acc;
    }
}

// h^a_b = delta^a_b + V^a V_b
template <typename Scalar>
Mat4<Scalar> projector_mixed(const Mat4<Scalar>& g_lab, const Scalar V_up[4]) {
    Scalar V_dn[4];
    lower_vector(g_lab, V_up, V_dn);
    Mat4<Scalar> h;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            h(a, b) = scalar_cast<Scalar>(a == b ? 1.0 : 0.0) + V_up[a] * V_dn[b];
    return h;
}

// h^{ab} = g_lab^{ab} + V^a V^b
template <typename Scalar>
Mat4<Scalar> projector_upper(const Mat4<Scalar>& g_lab, const Scalar V_up[4]) {
    const warp_hypad::InverseResult<Scalar> gi = warp_hypad::inverse(g_lab);
    Mat4<Scalar> h;
    if (!gi.ok) return h;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            h(a, b) = gi.inverse(a, b) + V_up[a] * V_up[b];
    return h;
}

// h_{ab} = g_lab_{ab} + V_a V_b
template <typename Scalar>
Mat4<Scalar> projector_lower(const Mat4<Scalar>& g_lab, const Scalar V_up[4]) {
    Scalar V_dn[4];
    lower_vector(g_lab, V_up, V_dn);
    Mat4<Scalar> h;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) h(a, b) = g_lab(a, b) + V_dn[a] * V_dn[b];
    return h;
}

// Project a rank-2 contravariant tensor: Xhat^{ab} = h^a_c h^b_d X^{cd}.
template <typename Scalar>
Mat4<Scalar> project_tensor(const Mat4<Scalar>& X, const Mat4<Scalar>& g_lab,
                            const Scalar V_up[4]) {
    const Mat4<Scalar> h = projector_mixed(g_lab, V_up);
    Mat4<Scalar> out;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d) acc = acc + h(a, c) * h(b, d) * X(c, d);
            out(a, b) = acc;
        }
    return out;
}

// The generic A7 contraction: W^{bd} = V_a V_c T^{abcd}.
// Indices 0 and 2 of T are contracted; 1 and 3 stay free.
template <typename Scalar>
Mat4<Scalar> contract_a7(const Rank4<Scalar>& T, const Scalar V_dn[4]) {
    Mat4<Scalar> W;
    for (int b = 0; b < 4; ++b)
        for (int d = 0; d < 4; ++d) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int a = 0; a < 4; ++a)
                for (int c = 0; c < 4; ++c)
                    acc = acc + V_dn[a] * V_dn[c] * T(a, b, c, d);
            W(b, d) = acc;
        }
    return W;
}

template <typename Scalar>
struct BelBlocks {
    Mat4<Scalar> W_epsilon;   // (A7a)
    Mat4<Scalar> W_mu;        // (A7b)
    Mat4<Scalar> W_zeta_T;    // (A7c)
    Mat4<Scalar> W_zeta;      // (A7d)
};

template <typename Scalar>
BelBlocks<Scalar> bel_extract(const Rank4<Scalar>& Z, const Mat4<Scalar>& g_lab,
                              const Scalar V_up[4]) {
    Scalar V_dn[4];
    lower_vector(g_lab, V_up, V_dn);

    BelBlocks<Scalar> b;
    b.W_epsilon = contract_a7(Z, V_dn);
    b.W_mu = contract_a7(dual_both(Z, g_lab), V_dn);
    b.W_zeta_T = contract_a7(dual_left(Z, g_lab), V_dn);
    b.W_zeta = contract_a7(dual_right(Z, g_lab), V_dn);
    return b;
}

// DETERMINED at step 5 from the (B17) moving-medium benchmark.
//
//     zeta^{bd} = +2 W_zeta^{bd}
//
// The value and its sign were measured, not chosen: across five boost
// configurations (three directions, four eps/mu pairs) the (A7d) block and the
// (B17) closed form are proportional componentwise with residual exactly 0.0,
// and the single ratio is +2 to within 6.7e-16.
//
// NOTE the sign. The provisional value was -2, taken by symmetry with (D1);
// that was WRONG. The pattern is +2 for both mu^-1 (D2) and zeta (D4), with
// eps (D1) the sole -2 -- which is why guessing from symmetry failed and why a
// benchmark with a definite sign was required.
inline constexpr double kZetaNormalisation = 2.0;

template <typename Scalar>
struct ObserverDecomposition4 {
    Mat4<Scalar> epsilon;
    Mat4<Scalar> mu_inverse;
    Mat4<Scalar> zeta;
    Mat4<Scalar> zeta_adjoint;

    // Pinned at step 5 against (B17). Readers of the magnetoelectric blocks
    // may now rely on their scale and sign.
    static constexpr bool zeta_normalisation_verified = true;
};

template <typename Scalar>
ObserverDecomposition4<Scalar> decompose_A7(const Rank4<Scalar>& Z,
                                            const Mat4<Scalar>& g_lab,
                                            const Scalar V_up[4]) {
    const BelBlocks<Scalar> b = bel_extract(Z, g_lab, V_up);
    ObserverDecomposition4<Scalar> d;

    const Scalar m2 = scalar_cast<Scalar>(-2.0);
    const Scalar p2 = scalar_cast<Scalar>(2.0);
    const Scalar cz = scalar_cast<Scalar>(kZetaNormalisation);

    for (int i = 0; i < 16; ++i) {
        d.epsilon.m[i] = m2 * b.W_epsilon.m[i];       // (D1)
        d.mu_inverse.m[i] = p2 * b.W_mu.m[i];         // (D2)
        d.zeta.m[i] = cz * b.W_zeta.m[i];             // (D4)
        d.zeta_adjoint.m[i] = cz * b.W_zeta_T.m[i];   // (D3)
    }
    return d;
}

// Worst |V_a X^{ab}| and |X^{ab} V_b|, reported separately so a failure names
// which index lost spatiality.
template <typename Scalar>
void spatiality_residuals(const Mat4<Scalar>& X, const Mat4<Scalar>& g_lab,
                          const Scalar V_up[4], double& first_index,
                          double& second_index) {
    Scalar V_dn[4];
    lower_vector(g_lab, V_up, V_dn);
    first_index = 0.0;
    second_index = 0.0;
    for (int b = 0; b < 4; ++b) {
        Scalar s1 = scalar_cast<Scalar>(0.0), s2 = scalar_cast<Scalar>(0.0);
        for (int a = 0; a < 4; ++a) {
            s1 = s1 + V_dn[a] * X(a, b);
            s2 = s2 + X(b, a) * V_dn[a];
        }
        first_index = std::fmax(first_index, std::fabs(real_part(s1)));
        second_index = std::fmax(second_index, std::fabs(real_part(s2)));
    }
}

}  // namespace program_b
