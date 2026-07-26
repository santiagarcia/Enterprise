#pragma once

#include "program_b/effective_metric/bel_extraction.hpp"
#include "program_b/effective_metric/metric_types.hpp"

// ---------------------------------------------------------------------------
// B1C step 4: the full electromagnetic field reconstruction.
//
//     H^{ab} = Z^{abcd} F_{cd}                                       (task)
//            = G^{ab}                                (Schuster & Visser eq 4)
//
// The excitation tensor is called G^{ab} in the paper and H^{ab} in the task
// statement. They are the same object; the three-vector magnetic excitation is
// also called H, so this file uses `excitation` for the rank-2 tensor and
// reserves the bare name H for the three-vector.
//
// TWO INDEPENDENT PATHS, which is the entire point of the step:
//
//   DIRECT   H^{ab} = Z^{abcd} F_{cd}, one contraction, no observer anywhere.
//   BLOCK    split F into (E,B) for the observer V, apply eps, mu^-1, zeta and
//            zeta+ to get (D,H), rebuild the excitation from (D,H).
//
// They must agree exactly. The direct path never mentions V; the block path
// mentions nothing else. If the observer decomposition were wrong in scale,
// sign, transposition or index pair, the two would part company on any medium
// with a nonzero magnetoelectric response.
//
// --- Observer split of a two-form ------------------------------------------
//
// Schuster & Visser (A1): for every four-velocity V^a there exist unique E and
// B with
//
//     F_{ab} = V_a E_b - V_b E_a + eps_{abcd} V^c B^d               (A1)
//
// Raising both free indices puts this in the form used throughout this file,
// which is the same for the field strength and for the excitation:
//
//     X^{ab} = V^a x^b - V^b x^a + eps^{abcd} V_c y_d
//
// with the inverse relations (derived below, and both verified numerically):
//
//     x^a = X^{ab} V_b                 y^a = -(*X)^{ab} V_b
//
// Derivation of the second. Contract (A1) with (1/2) eps^{efab} and use
// eps^{abef} eps_{abcd} = -2 (delta^e_c delta^f_d - delta^e_d delta^f_c):
//
//     (1/2) eps^{ef}_{~~ab} X^{ab} = eps^{efab} V_a x_b - (V^e y^f - V^f y^e)
//
// then contract with V_f. The first term dies because eps^{efab} V_a V_f = 0,
// and V.y = 0, V.V = -1 leave exactly -y^e. Hence y^e = -(*X)^{ef} V_f.
//
// Applied to F this gives (E, B); applied to H^{ab} it gives (D, H). Vacuum
// fixes the naming: with Z^{abcd} = (1/2)(g^{ac}g^{bd} - g^{ad}g^{bc}) the
// excitation is H^{ab} = F^{ab}, so D = E and H = B, as it must be.
//
// --- The constitutive relations in terms of the A7 blocks ------------------
//
// DERIVED here, not transcribed. Substituting (A1) into H^{ab} = Z^{abcd}F_{cd}
// and contracting with V_b:
//
//     D^e = H^{ef} V_f
//         = 2 Z^{efcd} V_f V_c E_d  +  Z^{efcd} V_f eps_{cdmn} V^m B^n
//         = -2 W_eps^{ed} E_d       +  2 (Z*)^{efcd} V_f V_c B_d
//         = -2 W_eps^{ed} E_d       -  2 W_zeta^{ed} B_d
//
// using the antisymmetry of Z and of (Z*) in their FIRST index pair to move
// V_f from slot 1 to slot 0, which costs one minus sign each time. With
// eps = -2 W_eps (D1) and zeta = +2 W_zeta (D4, pinned at step 5 against B17):
//
//     D^a = eps^{ab} E_b - zeta^{ab} B_b                             (S4-D)
//
// The same substitution for H^e = -(*H)^{ef} V_f, where the left dual turns Z
// into (*Z) and (*Z*) and therefore brings in (A7c) and (A7b):
//
//     H^e = -[ 2 (*Z)^{efmn} V_f V_m E_n + (*Z)^{efmn} V_f eps_{mnpq} V^p B^q ]
//         = -[ -2 [W_zeta^T]^{en} E_n  -  2 W_mu^{eq} B_q ]
//
//     H^a = zeta_adjoint^{ab} E_b + [mu^-1]^{ab} B_b                 (S4-H)
//
// NOTE THE ASYMMETRY between (S4-D) and (S4-H). It is not a typo and it is not
// a free convention: D is read off the excitation directly while H is read off
// its dual, and that single extra minus sign is the whole difference. The same
// asymmetry appears in the Lagrangian form. With F_{ab}H^{ab} = -2 E.D + 2 B.H
// the action (9) is S = (1/2) INT (E.D - B.H), and (S4-D)/(S4-H) are exactly
//
//     D = +dL/dE,    H = -dL/dB,    L = (1/2)[ E.eps.E - 2 E.zeta.B - B.mu^-1.B ]
//
// which is symmetric in the only sense that matters: the quadratic form is
// symmetric, i.e. Z^{abcd} = Z^{cdab}.
//
// Relative to the paper's 3+1 equation (15), written D = eps E + zeta B and
// H = zeta+ E + mu^-1 B, the D relation here carries the opposite sign on zeta.
// That is a convention gap between the paper's three-dimensional matrices of
// equation (20) and its own covariant blocks of (A6)/(A7), NOT a choice made
// here: with eps^{0123} = -1/sqrt(-det g) forced by raising four indices in
// (-,+,+,+), and with zeta = +2 W_zeta measured against (B17), (S4-D) follows
// with no freedom left. It is checked against a source outside the paper --
// Minkowski's moving-medium relations, D + v x H = eps(E + v x B), which to
// first order give D = eps E + (eps - mu^-1)(v x B) -- and that independent
// check confirms the sign in (S4-D), including its handedness. See the
// handedness benchmark in the B1C-S456 gate.
// ---------------------------------------------------------------------------

namespace program_b {

// A vector orthogonal to the observer, carried in both index positions so no
// call site has to remember which one it holds.
template <typename Scalar>
struct SpatialVector {
    Scalar up[4];
    Scalar dn[4];

    SpatialVector() {
        const Scalar z = scalar_cast<Scalar>(0.0);
        for (int a = 0; a < 4; ++a) {
            up[a] = z;
            dn[a] = z;
        }
    }
};

template <typename Scalar>
SpatialVector<Scalar> spatial_from_upper(const LaboratoryMetric<Scalar>& g_lab,
                                         const Scalar u[4]) {
    SpatialVector<Scalar> v;
    for (int a = 0; a < 4; ++a) v.up[a] = u[a];
    lower_vector(g_lab.value, u, v.dn);
    return v;
}

// Worst |V_a x^a| -- zero for every one of E, B, D, H, structurally.
template <typename Scalar>
double spatiality_residual(const SpatialVector<Scalar>& x,
                           const LaboratoryMetric<Scalar>& g_lab,
                           const Scalar V_up[4]) {
    Scalar V_dn[4];
    lower_vector(g_lab.value, V_up, V_dn);
    Scalar s = scalar_cast<Scalar>(0.0);
    for (int a = 0; a < 4; ++a) s = s + V_dn[a] * x.up[a];
    return std::fabs(real_part(s));
}

// --- Index gymnastics on two-forms ------------------------------------------

// X^{ab} = g^{ac} g^{bd} X_{cd}
template <typename Scalar>
Mat4<Scalar> raise_two_form(const Mat4<Scalar>& X_dn,
                            const LaboratoryMetric<Scalar>& g_lab) {
    const warp_hypad::InverseResult<Scalar> gi = warp_hypad::inverse(g_lab.value);
    Mat4<Scalar> out;
    if (!gi.ok) return out;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d)
                    acc = acc + gi.inverse(a, c) * gi.inverse(b, d) * X_dn(c, d);
            out(a, b) = acc;
        }
    return out;
}

// X_{ab} = g_{ac} g_{bd} X^{cd}
template <typename Scalar>
Mat4<Scalar> lower_two_form(const Mat4<Scalar>& X_up,
                            const LaboratoryMetric<Scalar>& g_lab) {
    Mat4<Scalar> out;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d)
                    acc = acc + g_lab.value(a, c) * g_lab.value(b, d) * X_up(c, d);
            out(a, b) = acc;
        }
    return out;
}

// --- The observer split -----------------------------------------------------

template <typename Scalar>
struct TwoFormSplit {
    SpatialVector<Scalar> electric;   // x^a = X^{ab} V_b
    SpatialVector<Scalar> magnetic;   // y^a = -(*X)^{ab} V_b
};

template <typename Scalar>
TwoFormSplit<Scalar> split_two_form(const Mat4<Scalar>& X_up,
                                    const LaboratoryMetric<Scalar>& g_lab,
                                    const Scalar V_up[4]) {
    Scalar V_dn[4];
    lower_vector(g_lab.value, V_up, V_dn);
    const Mat4<Scalar> Xd = dual_2form(X_up, g_lab.value);

    Scalar x_up[4], y_up[4];
    for (int a = 0; a < 4; ++a) {
        Scalar sx = scalar_cast<Scalar>(0.0);
        Scalar sy = scalar_cast<Scalar>(0.0);
        for (int b = 0; b < 4; ++b) {
            sx = sx + X_up(a, b) * V_dn[b];
            sy = sy - Xd(a, b) * V_dn[b];
        }
        x_up[a] = sx;
        y_up[a] = sy;
    }

    TwoFormSplit<Scalar> s;
    s.electric = spatial_from_upper(g_lab, x_up);
    s.magnetic = spatial_from_upper(g_lab, y_up);
    return s;
}

// X^{ab} = V^a x^b - V^b x^a + eps^{abcd} V_c y_d
template <typename Scalar>
Mat4<Scalar> two_form_from_split(const SpatialVector<Scalar>& x,
                                 const SpatialVector<Scalar>& y,
                                 const LaboratoryMetric<Scalar>& g_lab,
                                 const Scalar V_up[4]) {
    Scalar V_dn[4];
    lower_vector(g_lab.value, V_up, V_dn);
    const Rank4<Scalar> e_up = levi_civita_upper(g_lab.value);

    Mat4<Scalar> X;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
            Scalar acc = V_up[a] * x.up[b] - V_up[b] * x.up[a];
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d)
                    acc = acc + e_up(a, b, c, d) * V_dn[c] * y.dn[d];
            X(a, b) = acc;
        }
    return X;
}

// --- The two paths ----------------------------------------------------------

// H^{ab} = Z^{abcd} F_{cd}. No observer.
template <typename Scalar>
Mat4<Scalar> excitation_direct(const Rank4<Scalar>& Z, const Mat4<Scalar>& F_dn) {
    Mat4<Scalar> H;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d) acc = acc + Z(a, b, c, d) * F_dn(c, d);
            H(a, b) = acc;
        }
    return H;
}

template <typename Scalar>
struct ObserverFields {
    SpatialVector<Scalar> E;   // E^a, E_a
    SpatialVector<Scalar> B;
    SpatialVector<Scalar> D;
    SpatialVector<Scalar> H;   // magnetic excitation
};

// (D, H) from (E, B) through the four constitutive blocks -- equations (S4-D)
// and (S4-H) above.
template <typename Scalar>
void apply_constitutive_blocks(const ObserverDecomposition4<Scalar>& blocks,
                               const SpatialVector<Scalar>& E,
                               const SpatialVector<Scalar>& B,
                               const LaboratoryMetric<Scalar>& g_lab,
                               SpatialVector<Scalar>& D,
                               SpatialVector<Scalar>& H) {
    Scalar D_up[4], H_up[4];
    for (int a = 0; a < 4; ++a) {
        Scalar d = scalar_cast<Scalar>(0.0);
        Scalar h = scalar_cast<Scalar>(0.0);
        for (int b = 0; b < 4; ++b) {
            d = d + blocks.epsilon(a, b) * E.dn[b] - blocks.zeta(a, b) * B.dn[b];
            h = h + blocks.zeta_adjoint(a, b) * E.dn[b] +
                blocks.mu_inverse(a, b) * B.dn[b];
        }
        D_up[a] = d;
        H_up[a] = h;
    }
    D = spatial_from_upper(g_lab, D_up);
    H = spatial_from_upper(g_lab, H_up);
}

// The full block path: F_{ab} -> (E,B) -> (D,H) -> H^{ab}.
template <typename Scalar>
ObserverFields<Scalar> decompose_fields(const ObserverDecomposition4<Scalar>& blocks,
                                        const Mat4<Scalar>& F_dn,
                                        const LaboratoryMetric<Scalar>& g_lab,
                                        const Scalar V_up[4]) {
    const Mat4<Scalar> F_up = raise_two_form(F_dn, g_lab);
    const TwoFormSplit<Scalar> s = split_two_form(F_up, g_lab, V_up);

    ObserverFields<Scalar> f;
    f.E = s.electric;
    f.B = s.magnetic;
    apply_constitutive_blocks(blocks, f.E, f.B, g_lab, f.D, f.H);
    return f;
}

template <typename Scalar>
Mat4<Scalar> excitation_from_blocks(const ObserverDecomposition4<Scalar>& blocks,
                                    const Mat4<Scalar>& F_dn,
                                    const LaboratoryMetric<Scalar>& g_lab,
                                    const Scalar V_up[4]) {
    const ObserverFields<Scalar> f =
        decompose_fields(blocks, F_dn, g_lab, V_up);
    return two_form_from_split(f.D, f.H, g_lab, V_up);
}

// --- Two-form constructors used by the tests --------------------------------

// F_{ab} from three-vector components in a Minkowski laboratory with
// V = (1,0,0,0): F_{0i} = -E_i, F_{ij} = [ijk] B_k. Provided here rather than
// in the test file because it is the concrete meaning of (A1) and belongs with
// it.
template <typename Scalar>
Mat4<Scalar> field_strength_from_3vectors(const Scalar E3[3], const Scalar B3[3]) {
    Mat4<Scalar> F;
    for (int i = 0; i < 3; ++i) {
        F(0, i + 1) = -E3[i];
        F(i + 1, 0) = E3[i];
    }
    F(1, 2) = B3[2];
    F(2, 1) = -B3[2];
    F(2, 3) = B3[0];
    F(3, 2) = -B3[0];
    F(3, 1) = B3[1];
    F(1, 3) = -B3[1];
    return F;
}

template <typename Scalar>
double antisymmetry_residual(const Mat4<Scalar>& X) {
    double worst = 0.0;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            worst = std::fmax(worst,
                              std::fabs(real_part(X(a, b)) + real_part(X(b, a))));
    return worst;
}

}  // namespace program_b
