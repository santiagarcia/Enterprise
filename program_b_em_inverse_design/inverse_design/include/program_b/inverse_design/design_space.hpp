#pragma once

#include <string>

#include "program_b/effective_metric/bel_extraction.hpp"
#include "program_b/effective_metric/covariant_constitutive.hpp"
#include "program_b/effective_metric/z_reconstruction.hpp"

// ---------------------------------------------------------------------------
// B2 design space: an unconstrained parameter vector q, mapped into the three
// observer blocks so that spatiality, symmetry and reciprocity are STRUCTURAL.
//
// Nothing is projected or symmetrised after the fact, so no residual can appear
// and no optimizer iterate can wander off the constraint surface.
//
// THE OBSERVER FRAME. A design is expressed relative to one fixed laboratory
// observer. Its blocks must be spatial with respect to V, so the parameters are
// coefficients in an orthonormal spatial triad adapted to V:
//
//     eps^{ab} = sum_ij eps_ij e_(i)^a e_(j)^b,     g_lab(e_(i), e_(j)) = delta_ij
//
// The triad is a CHOICE OF BASIS, built by real-part Gram-Schmidt. That is
// legitimate here in a way it was not in the quarantined legacy decomposition:
// there a basis choice was used to extract physics FROM a tensor; here it only
// supplies coordinates ON the design space. The reachable set of tensors does not
// depend on it, and it is built from the fixed inputs (g_lab, V) alone -- never
// from the design variables -- so no derivative channel flows through it.
//
// For g_lab = eta with V = (1,0,0,0) the triad is the coordinate triad and
// eps_ij is literally eps^{(i+1)(j+1)}.
// ---------------------------------------------------------------------------

namespace program_b {
namespace inverse_design {

using warp_hypad::Mat4;
using warp_hypad::Rank4;
using warp_hypad::real_part;
using warp_hypad::scalar_cast;

// --- Observer frame ----------------------------------------------------------

template <typename Scalar>
struct ObserverFrame {
    Scalar V[4];         // normalised laboratory four-velocity, g_lab(V,V) = -1
    Scalar e[3][4];      // spatial triad e_(i)^a
    Scalar e_dn[3][4];   // e_(i)_a = g_lab_{ab} e_(i)^b
    bool ok = false;
    std::string message;
};

// Gram-Schmidt against g_lab starting from V, on real parts. Deterministic: the
// seed directions are tried in a fixed order, so the same (g_lab, V) always
// yields the same triad and therefore the same parameterization.
template <typename Scalar>
ObserverFrame<Scalar> build_observer_frame(const LaboratoryMetric<Scalar>& g_lab,
                                          const Observer<Scalar>& obs) {
    ObserverFrame<Scalar> f;
    if (!obs.admissible) {
        f.message = "the laboratory observer is not admissible: " + obs.reason;
        return f;
    }
    for (int a = 0; a < 4; ++a) f.V[a] = obs.V[a];

    double g[4][4];
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) g[a][b] = real_part(g_lab.value(a, b));

    auto dot = [&](const double u[4], const double v[4]) {
        double s = 0.0;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) s += g[a][b] * u[a] * v[b];
        return s;
    };

    double Vr[4];
    for (int a = 0; a < 4; ++a) Vr[a] = real_part(obs.V[a]);

    double basis[3][4] = {};
    int filled = 0;
    for (int cand = 0; cand < 4 && filled < 3; ++cand) {
        double v[4] = {0, 0, 0, 0};
        v[cand] = 1.0;

        // Remove the timelike part: V has norm -1, so the coefficient is
        // +dot(v,V) rather than -dot(v,V).
        const double cv = dot(v, Vr);
        for (int a = 0; a < 4; ++a) v[a] += cv * Vr[a];

        for (int k = 0; k < filled; ++k) {
            const double c = dot(v, basis[k]);
            for (int a = 0; a < 4; ++a) v[a] -= c * basis[k][a];
        }
        const double n = dot(v, v);
        if (n <= 1e-10) continue;   // degenerate seed direction
        const double s = 1.0 / std::sqrt(n);
        for (int a = 0; a < 4; ++a) basis[filled][a] = v[a] * s;
        ++filled;
    }
    if (filled < 3) {
        f.message = "could not build a spatial triad orthogonal to V.";
        return f;
    }

    for (int i = 0; i < 3; ++i)
        for (int a = 0; a < 4; ++a) f.e[i][a] = scalar_cast<Scalar>(basis[i][a]);
    for (int i = 0; i < 3; ++i) lower_vector(g_lab.value, f.e[i], f.e_dn[i]);

    f.ok = true;
    f.message = "orthonormal spatial triad adapted to V.";
    return f;
}

// --- Families ----------------------------------------------------------------

enum class DesignFamily {
    UnrestrictedReciprocal,   // 21: the whole reciprocal design space
    PositiveDefinite,         // 21: eps = LL^T + lambda_min I, likewise mu^-1
    ZeroMagnetoelectric,      // 12: zeta == 0
    DiagonalAnisotropic,      // 6:  eps, mu^-1 diagonal; zeta == 0
    Isotropic,                // 2:  eps = e0 I, mu^-1 = n0 I; zeta == 0
};

const char* to_string(DesignFamily f);
const char* describe(DesignFamily f);

inline int parameter_count(DesignFamily f) {
    switch (f) {
        case DesignFamily::UnrestrictedReciprocal: return 21;
        case DesignFamily::PositiveDefinite: return 21;
        case DesignFamily::ZeroMagnetoelectric: return 12;
        case DesignFamily::DiagonalAnisotropic: return 6;
        case DesignFamily::Isotropic: return 2;
    }
    return 0;
}

inline constexpr int kMaxDesignParameters = 21;

// --- Design variables --------------------------------------------------------

template <typename Scalar>
struct ConstitutiveDesignVariables {
    DesignFamily family = DesignFamily::UnrestrictedReciprocal;
    Scalar q[kMaxDesignParameters];

    // Floor added to the diagonal in the PositiveDefinite family. Not a design
    // variable: it is the stated lower bound on the spectrum.
    double positive_definite_floor = 1e-3;

    ConstitutiveDesignVariables() {
        const Scalar z = scalar_cast<Scalar>(0.0);
        for (int i = 0; i < kMaxDesignParameters; ++i) q[i] = z;
    }

    int size() const { return parameter_count(family); }
};

// The six independent entries of a symmetric 3x3, in the fixed order
// (11, 12, 13, 22, 23, 33).
inline void symmetric_index(int slot, int& i, int& j) {
    static const int I[6] = {0, 0, 0, 1, 1, 2};
    static const int J[6] = {0, 1, 2, 1, 2, 2};
    i = I[slot];
    j = J[slot];
}

// The three 3x3 coefficient matrices a design vector denotes. Held as plain
// 3x3 arrays because at this stage they are coordinates, not tensors.
template <typename Scalar>
struct BlockCoefficients {
    Scalar epsilon[3][3];
    Scalar mu_inverse[3][3];
    Scalar zeta[3][3];

    BlockCoefficients() {
        const Scalar z = scalar_cast<Scalar>(0.0);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                epsilon[i][j] = z;
                mu_inverse[i][j] = z;
                zeta[i][j] = z;
            }
    }
};

// q -> coefficients. Symmetry of eps and mu^-1 is structural: both (i,j) and
// (j,i) read the SAME parameter, so no symmetrisation step exists to forget.
template <typename Scalar>
BlockCoefficients<Scalar> coefficients_from_parameters(
    const ConstitutiveDesignVariables<Scalar>& v) {
    BlockCoefficients<Scalar> c;
    const Scalar floor = scalar_cast<Scalar>(v.positive_definite_floor);

    switch (v.family) {
        case DesignFamily::UnrestrictedReciprocal: {
            for (int s = 0; s < 6; ++s) {
                int i, j;
                symmetric_index(s, i, j);
                c.epsilon[i][j] = v.q[s];
                c.epsilon[j][i] = v.q[s];
                c.mu_inverse[i][j] = v.q[6 + s];
                c.mu_inverse[j][i] = v.q[6 + s];
            }
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) c.zeta[i][j] = v.q[12 + 3 * i + j];
            break;
        }
        case DesignFamily::PositiveDefinite: {
            // eps = L L^T + floor * I, with L lower triangular from
            // (q0; q1,q2; q3,q4,q5) read row by row.
            static const int LR[6] = {0, 1, 1, 2, 2, 2};
            static const int LC[6] = {0, 0, 1, 0, 1, 2};
            Scalar Le[3][3] = {}, Lm[3][3] = {};
            const Scalar z = scalar_cast<Scalar>(0.0);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) { Le[i][j] = z; Lm[i][j] = z; }
            for (int s = 0; s < 6; ++s) {
                Le[LR[s]][LC[s]] = v.q[s];
                Lm[LR[s]][LC[s]] = v.q[6 + s];
            }
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    Scalar ae = z, am = z;
                    for (int k = 0; k < 3; ++k) {
                        ae = ae + Le[i][k] * Le[j][k];
                        am = am + Lm[i][k] * Lm[j][k];
                    }
                    c.epsilon[i][j] = ae + (i == j ? floor : z);
                    c.mu_inverse[i][j] = am + (i == j ? floor : z);
                }
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) c.zeta[i][j] = v.q[12 + 3 * i + j];
            break;
        }
        case DesignFamily::ZeroMagnetoelectric: {
            for (int s = 0; s < 6; ++s) {
                int i, j;
                symmetric_index(s, i, j);
                c.epsilon[i][j] = v.q[s];
                c.epsilon[j][i] = v.q[s];
                c.mu_inverse[i][j] = v.q[6 + s];
                c.mu_inverse[j][i] = v.q[6 + s];
            }
            break;
        }
        case DesignFamily::DiagonalAnisotropic: {
            for (int i = 0; i < 3; ++i) {
                c.epsilon[i][i] = v.q[i];
                c.mu_inverse[i][i] = v.q[3 + i];
            }
            break;
        }
        case DesignFamily::Isotropic: {
            for (int i = 0; i < 3; ++i) {
                c.epsilon[i][i] = v.q[0];
                c.mu_inverse[i][i] = v.q[1];
            }
            break;
        }
    }
    return c;
}

// Coefficients -> q, for the families where the map is invertible. Used by the
// round-trip case and to seed an optimizer from a known medium. Not defined for
// PositiveDefinite, whose map is not injective (L and -L agree).
template <typename Scalar>
bool parameters_from_coefficients(const BlockCoefficients<Scalar>& c,
                                 DesignFamily family,
                                 ConstitutiveDesignVariables<Scalar>& out) {
    out = ConstitutiveDesignVariables<Scalar>();
    out.family = family;
    switch (family) {
        case DesignFamily::UnrestrictedReciprocal:
        case DesignFamily::ZeroMagnetoelectric: {
            for (int s = 0; s < 6; ++s) {
                int i, j;
                symmetric_index(s, i, j);
                out.q[s] = c.epsilon[i][j];
                out.q[6 + s] = c.mu_inverse[i][j];
            }
            if (family == DesignFamily::UnrestrictedReciprocal)
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j) out.q[12 + 3 * i + j] = c.zeta[i][j];
            return true;
        }
        case DesignFamily::DiagonalAnisotropic: {
            for (int i = 0; i < 3; ++i) {
                out.q[i] = c.epsilon[i][i];
                out.q[3 + i] = c.mu_inverse[i][i];
            }
            return true;
        }
        case DesignFamily::Isotropic: {
            out.q[0] = c.epsilon[0][0];
            out.q[1] = c.mu_inverse[0][0];
            return true;
        }
        case DesignFamily::PositiveDefinite:
            return false;   // not injective; no inverse is offered
    }
    return false;
}

// --- Coefficients -> four-dimensional blocks --------------------------------

template <typename Scalar>
Mat4<Scalar> block_from_coefficients(const Scalar m[3][3],
                                     const ObserverFrame<Scalar>& frame) {
    Mat4<Scalar> out;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    acc = acc + m[i][j] * frame.e[i][a] * frame.e[j][b];
            out(a, b) = acc;
        }
    return out;
}

// The inverse read-off, using the orthonormality of the triad:
//     m_ij = X^{ab} e_(i)_a e_(j)_b
template <typename Scalar>
void coefficients_from_block(const Mat4<Scalar>& X,
                             const ObserverFrame<Scalar>& frame, Scalar m[3][3]) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    acc = acc + X(a, b) * frame.e_dn[i][a] * frame.e_dn[j][b];
            m[i][j] = acc;
        }
}

// q -> the four observer blocks. zeta_adjoint = zeta^T is imposed here, once,
// structurally: for the reciprocal design space it is an identity, not a choice.
template <typename Scalar>
ObserverDecomposition4<Scalar> blocks_from_parameters(
    const ConstitutiveDesignVariables<Scalar>& v,
    const ObserverFrame<Scalar>& frame) {
    const BlockCoefficients<Scalar> c = coefficients_from_parameters(v);
    ObserverDecomposition4<Scalar> d;
    d.epsilon = block_from_coefficients(c.epsilon, frame);
    d.mu_inverse = block_from_coefficients(c.mu_inverse, frame);
    d.zeta = block_from_coefficients(c.zeta, frame);

    Scalar zt[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) zt[i][j] = c.zeta[j][i];
    d.zeta_adjoint = block_from_coefficients(zt, frame);
    return d;
}

// q -> Z_design, through the A11 implementation VERIFIED AT B1C. There is no
// second production reconstruction formula anywhere in B2.
template <typename Scalar>
Rank4<Scalar> reconstruct_design_tensor(
    const ConstitutiveDesignVariables<Scalar>& v,
    const LaboratoryMetric<Scalar>& g_lab, const ObserverFrame<Scalar>& frame) {
    return reconstruct_Z_A11(blocks_from_parameters(v, frame), g_lab, frame.V);
}

}  // namespace inverse_design
}  // namespace program_b
