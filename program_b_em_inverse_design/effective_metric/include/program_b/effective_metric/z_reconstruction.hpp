#pragma once

#include "program_b/effective_metric/bel_extraction.hpp"
#include "program_b/effective_metric/metric_types.hpp"

// ---------------------------------------------------------------------------
// B1C step 6: rebuild Z^{abcd} from the observer blocks, Schuster & Visser
// equation (A11).
//
//     ( eps, mu^-1, zeta, zeta_adjoint, V, g_lab )  ->  Zhat^{abcd}
//
// Transcribed from the extracted source text, with h^{ab} = g^{ab} + V^aV^b:
//
//   Z^{abcd} = 1/2 [  V^d V^a eps^{bc} - V^c V^a eps^{bd}
//                   + V^c V^b eps^{ad} - V^d V^b eps^{ac}
//                   + h^{ad} [mu^-1]^{cb} - h^{ac} [mu^-1]^{db}
//                   + h^{bc} [mu^-1]^{ad} - h^{bd} [mu^-1]^{ac}
//                   + ( h^{bd} h^{ac} - h^{bc} h^{ad} ) [mu^-1]^e_e
//                   + eps^{fabe} ( V^d Y_e^{~c} - V^c Y_e^{~d} ) V_f
//                   + eps^{fcde} ( V^b Y_e^{~a} - V^a Y_e^{~b} ) V_f ]
//
// NOT A FIT. Every coefficient is read straight off the equation; nothing is
// adjusted to make a round trip close.
//
// --- WHICH magnetoelectric block Y is ---------------------------------------
//
// The source writes the last two terms with `zeta_e^{~c}`. In this project's
// conventions that slot is filled by ZETA_ADJOINT, not zeta:
//
//     Y_e^{~c} = [zeta_adjoint]_e^{~c} = g_{ef} [zeta_adjoint]^{fc}
//
// DETERMINED, not chosen, and determined by a test that no sign convention
// could pass. The candidates were zeta and zeta_adjoint, each with either
// overall sign -- four possibilities. Measured round-trip error max|Zhat - Z|:
//
//     medium                      zeta(+)   zeta(-)   zetaAdj(+)  zetaAdj(-)
//     generic electric block      0         0         0           0
//     generic magnetic block      1.1e-16   1.1e-16   1.1e-16     1.1e-16
//     zeta diagonal               0         1.0       0           1.0
//     zeta purely off-diagonal    5.0e-1    5.0e-1    0           1.0
//     moving isotropic, boost x   1.1       2.8e-17   2.8e-17     1.1
//     random reciprocal 6x6       1.8       2.4       2.2e-16     3.0
//
// Only the zeta_adjoint(+) column is exact everywhere. The two rows that
// discriminate are the ones with a magnetoelectric block that is neither
// symmetric nor antisymmetric: for a diagonal zeta the transpose is the
// identity and for the moving-medium zeta the transpose is a sign flip, so
// EITHER of those alone would have licensed the wrong reading -- a sign
// constant of -1 reproduces the moving medium exactly and is still wrong. The
// purely off-diagonal case and the random reciprocal tensor are what settle it,
// and they are in the gate for that reason (b1c.a11.magnetoelectric_slot).
//
// This is the same convention gap already visible in step 4: the paper's
// three-dimensional zeta of equations (15)/(20) is the transpose of the
// covariant block that (A7d) produces, which is also why the derived
// constitutive relation reads D = eps E - zeta B rather than the (15) form.
// One consistent transposition explains both.
//
// --- SCOPE: reciprocal media only -------------------------------------------
//
// (A11) reconstructs Z from three blocks, which is 6 + 6 + 9 = 21 numbers --
// exactly the degrees of freedom of a Z carrying the paper's symmetries (12)
// and (13). Pair symmetry Z^{abcd} = Z^{cdab} is therefore a PRECONDITION, not
// a bonus, and under it zeta_adjoint = zeta^T identically:
//
//     [W_zeta^T]^{bd} = V_a V_c (*Z)^{abcd}
//                     = V_a V_c (Z*)^{cdab}        [pair symmetry]
//                     = W_zeta^{db}                [relabel a <-> c]
//
// A skewon-carrying Z has zeta_adjoint != zeta^T, falls outside (A11), and must
// be detected rather than silently reconstructed wrongly. Both blocks are
// passed in and the residual below is checked at every call site in the gate.
// ---------------------------------------------------------------------------

namespace program_b {

template <typename Scalar>
Rank4<Scalar> reconstruct_Z_A11(const ObserverDecomposition4<Scalar>& b,
                                const LaboratoryMetric<Scalar>& g_lab,
                                const Scalar V_up[4]) {
    Scalar V_dn[4];
    lower_vector(g_lab.value, V_up, V_dn);

    const Mat4<Scalar> h = projector_upper(g_lab.value, V_up);
    const Rank4<Scalar> e_up = levi_civita_upper(g_lab.value);

    // Y_e^{~c} = g_{ef} [zeta_adjoint]^{fc}, the first index lowered.
    Mat4<Scalar> Y;
    for (int e = 0; e < 4; ++e)
        for (int c = 0; c < 4; ++c) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int f = 0; f < 4; ++f)
                acc = acc + g_lab.value(e, f) * b.zeta_adjoint(f, c);
            Y(e, c) = acc;
        }

    // [mu^-1]^e_e = g_{ef} [mu^-1]^{ef}
    Scalar mu_trace = scalar_cast<Scalar>(0.0);
    for (int e = 0; e < 4; ++e)
        for (int f = 0; f < 4; ++f)
            mu_trace = mu_trace + g_lab.value(e, f) * b.mu_inverse(e, f);

    // eps^{fabe} V_f, contracted once and reused: the same rank-3 object serves
    // both magnetoelectric terms, and the slot ORDER of the Levi-Civita tensor
    // is the easiest thing in (A11) to get wrong, so it is written once.
    Scalar eV[4][4][4];
    for (int a = 0; a < 4; ++a)
        for (int bb = 0; bb < 4; ++bb)
            for (int e = 0; e < 4; ++e) {
                Scalar acc = scalar_cast<Scalar>(0.0);
                for (int f = 0; f < 4; ++f)
                    acc = acc + e_up(f, a, bb, e) * V_dn[f];
                eV[a][bb][e] = acc;
            }

    const Scalar half = scalar_cast<Scalar>(0.5);
    Rank4<Scalar> Z;
    for (int a = 0; a < 4; ++a)
        for (int bb = 0; bb < 4; ++bb)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d) {
                    Scalar acc = V_up[d] * V_up[a] * b.epsilon(bb, c) -
                                 V_up[c] * V_up[a] * b.epsilon(bb, d) +
                                 V_up[c] * V_up[bb] * b.epsilon(a, d) -
                                 V_up[d] * V_up[bb] * b.epsilon(a, c);

                    acc = acc + h(a, d) * b.mu_inverse(c, bb) -
                          h(a, c) * b.mu_inverse(d, bb) +
                          h(bb, c) * b.mu_inverse(a, d) -
                          h(bb, d) * b.mu_inverse(a, c);

                    acc = acc + (h(bb, d) * h(a, c) - h(bb, c) * h(a, d)) * mu_trace;

                    for (int e = 0; e < 4; ++e) {
                        acc = acc + eV[a][bb][e] *
                                        (V_up[d] * Y(e, c) - V_up[c] * Y(e, d));
                        acc = acc + eV[c][d][e] *
                                        (V_up[bb] * Y(e, a) - V_up[a] * Y(e, bb));
                    }

                    Z(a, bb, c, d) = half * acc;
                }
    return Z;
}

// Reciprocity residual: max |Z^{abcd} - Z^{cdab}|. Zero for every Z that comes
// from an effective metric, and the precondition for (A11).
template <typename Scalar>
double pair_symmetry_residual(const Rank4<Scalar>& Z) {
    double worst = 0.0;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d)
                    worst = std::fmax(worst,
                                      std::fabs(real_part(Z(a, b, c, d)) -
                                                real_part(Z(c, d, a, b))));
    return worst;
}

// max |zeta_adjoint - zeta^T|. Zero exactly when the medium is reciprocal.
template <typename Scalar>
double reciprocity_residual(const ObserverDecomposition4<Scalar>& b) {
    double worst = 0.0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            worst = std::fmax(worst, std::fabs(real_part(b.zeta_adjoint(i, j)) -
                                               real_part(b.zeta(j, i))));
    return worst;
}

}  // namespace program_b
