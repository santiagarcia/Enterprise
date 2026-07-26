#pragma once

#include "program_b/effective_metric/covariant_constitutive.hpp"

// ---------------------------------------------------------------------------
// QUARANTINE. Nothing in this file is on any gate path, and nothing in this
// file may be called from production code or from a new test.
//
// What lives here is the FIRST attempt at the observer decomposition of
// Z^{abcd}: build an orthonormal tetrad adapted to V, transform Z into that
// frame, and read the blocks off the frame components using three-dimensional
// Levi-Civita contractions:
//
//     eps^{ij}       = -2 Z^{0i0j}
//     [mu^{-1}]_{ij} = (1/2) eps_{ikl} eps_{jmn} Z^{klmn}
//     zeta^{ij}      =       eps_{jkl} Z^{0ikl}
//
// Why it is quarantined rather than trusted:
//
//   * eps and mu^-1 happen to be right, but only because the vacuum case
//     pins them and vacuum has eps = mu^-1 = 1. The agreement is not evidence
//     for the construction.
//   * the zeta normalisation here was GUESSED. It was never derived from a
//     primary source and never benchmarked against a medium with a known
//     nonzero magnetoelectric response. The determination that finally settled
//     it (Schuster & Visser B17, zeta = +2 W_zeta) was carried out against the
//     production path in bel_extraction.hpp, not against this one.
//   * the whole construction goes through a tetrad, which is a CHOICE of basis
//     built by Gram-Schmidt with pivot-like comparisons on real parts. The
//     production path (A7a)-(A7d) contracts the observer directly into Hodge
//     duals of Z and never picks a basis at all.
//
// It is retained, rather than deleted, for one reason only: the old B1C gate
// cases still compare against it, and deleting it would silently drop those
// comparisons rather than record them. When those cases are rewritten against
// the (A7)/(A11) path, this file should go.
//
// The namespace is the enforcement: `legacy_unverified::` has to be typed at
// every call site, so no new code reaches it by accident.
// ---------------------------------------------------------------------------

namespace program_b {
namespace legacy_unverified {

using program_b::Mat3;
using program_b::MetricMediumStatus;
using program_b::Observer;
using program_b::symmetric_eigenvalues_3x3;

// --- Orthonormal frame adapted to an observer --------------------------------
//
// Gram-Schmidt against g_lab, starting from V. Columns of the returned matrix
// are the frame vectors e^a_{(A)} with e_{(0)} = V, so
// g_lab(e_(A), e_(B)) = eta_{(A)(B)}.
//
// Real-part arithmetic: the frame is a choice of basis, a discrete construction
// with comparisons and a pivot selection. Propagating derivative channels
// through a basis CHOICE is meaningless -- what carries sensitivity is the
// tensor, not the ruler used to read it.

template <typename Scalar>
struct Tetrad {
    Mat4<double> e;        // e[a][A]
    Mat4<double> e_inv;    // co-frame theta^{(A)}_a
    bool ok = false;
};

template <typename Scalar>
Tetrad<Scalar> build_tetrad(const Mat4<Scalar>& g_lab, const Observer<Scalar>& obs) {
    Tetrad<Scalar> t;
    if (!obs.admissible) return t;

    double g[4][4];
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) g[a][b] = real_part(g_lab(a, b));

    auto dot = [&](const double u[4], const double v[4]) {
        double s = 0.0;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) s += g[a][b] * u[a] * v[b];
        return s;
    };

    double basis[4][4] = {};
    for (int a = 0; a < 4; ++a) basis[0][a] = real_part(obs.V[a]);

    // Seed the spatial legs with coordinate directions, choosing at each step
    // the candidate least aligned with what is already spanned.
    int filled = 1;
    for (int cand = 1; cand < 4 && filled < 4; ++cand) {
        double v[4] = {0, 0, 0, 0};
        v[cand] = 1.0;
        for (int k = 0; k < filled; ++k) {
            const double nk = dot(basis[k], basis[k]);
            if (std::fabs(nk) < 1e-300) continue;
            const double c = dot(v, basis[k]) / nk;
            for (int a = 0; a < 4; ++a) v[a] -= c * basis[k][a];
        }
        const double n = dot(v, v);
        if (n <= 1e-12) continue;  // degenerate direction, skip
        const double s = 1.0 / std::sqrt(n);
        for (int a = 0; a < 4; ++a) basis[filled][a] = v[a] * s;
        ++filled;
    }
    if (filled < 4) return t;

    for (int a = 0; a < 4; ++a)
        for (int A = 0; A < 4; ++A) t.e(a, A) = basis[A][a];

    const warp_hypad::InverseResult<double> ei = warp_hypad::inverse(t.e);
    if (!ei.ok) return t;
    t.e_inv = ei.inverse;
    t.ok = true;
    return t;
}

// --- Bel decomposition relative to an observer --------------------------------
//
// Transform Z into the observer's orthonormal frame and read off the blocks.
// With upper indices, the transformation contracts with the CO-frame:
//
//     Z^{(A)(B)(C)(D)} = theta^{(A)}_a theta^{(B)}_b theta^{(C)}_c theta^{(D)}_d Z^{abcd}
//
// then, in that frame:
//
//     eps^{ij}       = -2 Z^{0i0j}
//     [mu^{-1}]_{ij} = (1/2) eps_{ikl} eps_{jmn} Z^{klmn}
//     zeta^{ij}      = eps_{jkl} Z^{0ikl}
//
// Vacuum check, which fixes the first two factors and NOT the third:
//   Z^{0i0j} = (1/2)(eta^{00}eta^{ij}) = -delta^{ij}/2  =>  eps = I           OK
//   (1/2)eps_{ikl}eps_{jmn}(1/2)(d^{km}d^{ln}-d^{kn}d^{lm}) = delta_{ij}      OK
//   Z^{0ikl} = 0 since eta^{0k} = 0                     =>  zeta = 0   VACUOUS

struct LeviCivita3 {
    static int e(int i, int j, int k) {
        if (i == j || j == k || i == k) return 0;
        if ((i == 0 && j == 1 && k == 2) || (i == 1 && j == 2 && k == 0) ||
            (i == 2 && j == 0 && k == 1))
            return 1;
        return -1;
    }
};

template <typename Scalar>
struct ObserverDecomposition {
    Mat3<Scalar> epsilon;
    Mat3<Scalar> mu_inverse;
    Mat3<Scalar> zeta;

    MetricMediumStatus status =
        MetricMediumStatus::SelectedObserverDecompositionInvalid;
    std::string observer_name;
    std::string message;

    // Conditioning diagnostics on the real parts.
    double epsilon_cond = 0.0;
    double epsilon_eig_min = 0.0;
    double epsilon_eig_max = 0.0;

    bool ok() const { return status == MetricMediumStatus::CovariantMetricMediumValid; }
};

template <typename Scalar>
ObserverDecomposition<Scalar> decompose(const CovariantConstitutive<Scalar>& cc,
                                        const Mat4<Scalar>& g_lab,
                                        const Observer<Scalar>& obs,
                                        double condition_limit = 1e12) {
    ObserverDecomposition<Scalar> out;
    out.observer_name = obs.name;

    if (!cc.ok) {
        out.status = cc.status;
        out.message = cc.message;
        return out;
    }
    if (!obs.admissible) {
        out.status = MetricMediumStatus::SelectedObserverDecompositionInvalid;
        out.message = obs.reason;
        return out;
    }

    const Tetrad<Scalar> tet = build_tetrad(g_lab, obs);
    if (!tet.ok) {
        out.status = MetricMediumStatus::SelectedObserverDecompositionInvalid;
        out.message = "could not build an orthonormal frame adapted to V.";
        return out;
    }

    // Transform Z into the observer frame.
    warp_hypad::Rank4<Scalar> Zf;
    for (int A = 0; A < 4; ++A)
        for (int B = 0; B < 4; ++B)
            for (int C = 0; C < 4; ++C)
                for (int D = 0; D < 4; ++D) {
                    Scalar acc = scalar_cast<Scalar>(0.0);
                    for (int a = 0; a < 4; ++a)
                        for (int b = 0; b < 4; ++b)
                            for (int c = 0; c < 4; ++c)
                                for (int d = 0; d < 4; ++d)
                                    acc = acc +
                                          scalar_cast<Scalar>(
                                              tet.e_inv(A, a) * tet.e_inv(B, b) *
                                              tet.e_inv(C, c) * tet.e_inv(D, d)) *
                                              cc.Z(a, b, c, d);
                    Zf(A, B, C, D) = acc;
                }

    const Scalar minus_two = scalar_cast<Scalar>(-2.0);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            out.epsilon(i, j) = minus_two * Zf(0, i + 1, 0, j + 1);

    const Scalar half = scalar_cast<Scalar>(0.5);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    for (int m = 0; m < 3; ++m)
                        for (int n = 0; n < 3; ++n) {
                            const int s = LeviCivita3::e(i, k, l) *
                                          LeviCivita3::e(j, m, n);
                            if (s == 0) continue;
                            acc = acc + scalar_cast<Scalar>(static_cast<double>(s)) *
                                            Zf(k + 1, l + 1, m + 1, n + 1);
                        }
            out.mu_inverse(i, j) = half * acc;
        }

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l) {
                    const int s = LeviCivita3::e(j, k, l);
                    if (s == 0) continue;
                    acc = acc + scalar_cast<Scalar>(static_cast<double>(s)) *
                                    Zf(0, i + 1, k + 1, l + 1);
                }
            out.zeta(i, j) = acc;
        }

    // Conditioning of eps.
    double er[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) er[i][j] = real_part(out.epsilon(i, j));
    double eig[3];
    symmetric_eigenvalues_3x3(er, eig);
    out.epsilon_eig_min = std::fmin(std::fabs(eig[0]),
                                    std::fmin(std::fabs(eig[1]), std::fabs(eig[2])));
    out.epsilon_eig_max = std::fmax(std::fabs(eig[0]),
                                    std::fmax(std::fabs(eig[1]), std::fabs(eig[2])));
    out.epsilon_cond = (out.epsilon_eig_min > 0.0)
                           ? out.epsilon_eig_max / out.epsilon_eig_min
                           : std::numeric_limits<double>::infinity();

    bool finite = true;
    for (int i = 0; i < 9; ++i) {
        if (!std::isfinite(real_part(out.epsilon.m[i])) ||
            !std::isfinite(real_part(out.mu_inverse.m[i])) ||
            !std::isfinite(real_part(out.zeta.m[i])))
            finite = false;
    }

    if (!finite || out.epsilon_cond > condition_limit) {
        out.status = MetricMediumStatus::ConstitutiveParametersSingular;
        out.message =
            "the covariant tensor is finite and the observer is admissible, but "
            "the 3+1 material tensors for THIS observer diverge or are too "
            "ill conditioned to be usable. A different observer may still work.";
        return out;
    }

    out.status = MetricMediumStatus::CovariantMetricMediumValid;
    out.message = "covariant tensor valid, observer admissible, 3+1 tensors "
                  "finite and well conditioned.";
    return out;
}

// --- Full status report --------------------------------------------------------
//
// Reports the four questions SEPARATELY, which is the correction the covariant
// module exists to make. Quarantined with decompose(), which it calls.

struct MetricMediumReport {
    bool metric_lorentzian_nondegenerate = false;
    bool covariant_tensor_exists = false;
    bool observer_admissible = false;
    bool constitutive_finite = false;
    MetricMediumStatus status =
        MetricMediumStatus::MetricDegenerateOrNonLorentzian;
    double epsilon_cond = 0.0;
    std::string message;
};

template <typename Scalar>
MetricMediumReport analyse(const Mat4<Scalar>& g_eff, const Mat4<Scalar>& g_lab,
                           const Observer<Scalar>& obs) {
    MetricMediumReport r;
    const CovariantConstitutive<Scalar> cc = covariant_constitutive(g_eff, g_lab);

    r.metric_lorentzian_nondegenerate = cc.ok;
    r.covariant_tensor_exists = cc.ok;
    if (!cc.ok) {
        r.status = MetricMediumStatus::MetricDegenerateOrNonLorentzian;
        r.message = cc.message;
        return r;
    }

    r.observer_admissible = obs.admissible;
    const ObserverDecomposition<Scalar> d = decompose(cc, g_lab, obs);
    r.epsilon_cond = d.epsilon_cond;
    r.constitutive_finite =
        d.status == MetricMediumStatus::CovariantMetricMediumValid;
    r.status = d.status;
    r.message = d.message;
    return r;
}

}  // namespace legacy_unverified
}  // namespace program_b
