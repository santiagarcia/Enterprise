#pragma once

#include <cmath>

#include "warp_hypad/core/scalar_traits.hpp"

// ---------------------------------------------------------------------------
// Third-order spacetime jets, with arithmetic.
//
// WHY THIS EXISTS. O0B's metric jet is written out longhand for ONE family:
// alpha = 1, gamma = I, a single shift component. Its third derivatives of
// g_00 = -1 + (beta^x)^2 are four hand-derived product-rule terms. That was
// tractable because the algebra was small.
//
// F2A is not small. With alpha = exp(a) and a radial shift the covariant
// components become
//
//     g_00 = -exp(2a) + (beta^x)^2 + (beta^y)^2 + (beta^z)^2
//     beta^y = (y/R0) Q(xi, s)
//
// and hand-deriving d^3 of exp(2a) plus the products would be roughly forty
// terms per component, each a chance to lose a sign. Terms lost that way do not
// crash: they produce a curvature that is plausible and wrong, which is the
// failure mode this project exists to prevent.
//
// So the derivatives are propagated by ARITHMETIC instead. A Jet3 carries a
// value and its first three spacetime derivatives; +, * and exp implement the
// Leibniz and chain rules once, here, where they can be tested in isolation.
// Every metric component is then written the way it appears on paper.
//
// SYMMETRY IS STORED, NOT ASSUMED. d2[i][j] and d3[i][j][k] are kept as full
// arrays with all index orders filled. That costs a little memory and buys the
// ability to check symmetry as a measured residual rather than as a promise --
// which is exactly what O0B's permutation cases do.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

template <typename Scalar>
struct Jet3 {
    Scalar v;
    Scalar d1[4];
    Scalar d2[4][4];
    Scalar d3[4][4][4];

    Jet3() { set_zero(); }

    void set_zero() {
        const Scalar z = scalar_cast<Scalar>(0.0);
        v = z;
        for (int i = 0; i < 4; ++i) {
            d1[i] = z;
            for (int j = 0; j < 4; ++j) {
                d2[i][j] = z;
                for (int k = 0; k < 4; ++k) d3[i][j][k] = z;
            }
        }
    }

    static Jet3 constant(const Scalar& c) {
        Jet3 j;
        j.v = c;
        return j;
    }
};

template <typename Scalar>
Jet3<Scalar> operator+(const Jet3<Scalar>& a, const Jet3<Scalar>& b) {
    Jet3<Scalar> r;
    r.v = a.v + b.v;
    for (int i = 0; i < 4; ++i) {
        r.d1[i] = a.d1[i] + b.d1[i];
        for (int j = 0; j < 4; ++j) {
            r.d2[i][j] = a.d2[i][j] + b.d2[i][j];
            for (int k = 0; k < 4; ++k)
                r.d3[i][j][k] = a.d3[i][j][k] + b.d3[i][j][k];
        }
    }
    return r;
}

template <typename Scalar>
Jet3<Scalar> operator-(const Jet3<Scalar>& a) {
    Jet3<Scalar> r;
    r.v = -a.v;
    for (int i = 0; i < 4; ++i) {
        r.d1[i] = -a.d1[i];
        for (int j = 0; j < 4; ++j) {
            r.d2[i][j] = -a.d2[i][j];
            for (int k = 0; k < 4; ++k) r.d3[i][j][k] = -a.d3[i][j][k];
        }
    }
    return r;
}

template <typename Scalar>
Jet3<Scalar> operator-(const Jet3<Scalar>& a, const Jet3<Scalar>& b) {
    return a + (-b);
}

// The Leibniz rule to third order, written once.
//
//     (ab)_i     = a_i b + a b_i
//     (ab)_ij    = a_ij b + a_i b_j + a_j b_i + a b_ij
//     (ab)_ijk   = a_ijk b + a_ij b_k + a_ik b_j + a_jk b_i
//                + a_i b_jk + a_j b_ik + a_k b_ij + a b_ijk
template <typename Scalar>
Jet3<Scalar> operator*(const Jet3<Scalar>& a, const Jet3<Scalar>& b) {
    Jet3<Scalar> r;
    r.v = a.v * b.v;
    for (int i = 0; i < 4; ++i) {
        r.d1[i] = a.d1[i] * b.v + a.v * b.d1[i];
        for (int j = 0; j < 4; ++j) {
            r.d2[i][j] = a.d2[i][j] * b.v + a.d1[i] * b.d1[j] +
                         a.d1[j] * b.d1[i] + a.v * b.d2[i][j];
            for (int k = 0; k < 4; ++k)
                r.d3[i][j][k] =
                    a.d3[i][j][k] * b.v + a.d2[i][j] * b.d1[k] +
                    a.d2[i][k] * b.d1[j] + a.d2[j][k] * b.d1[i] +
                    a.d1[i] * b.d2[j][k] + a.d1[j] * b.d2[i][k] +
                    a.d1[k] * b.d2[i][j] + a.v * b.d3[i][j][k];
        }
    }
    return r;
}

template <typename Scalar>
Jet3<Scalar> operator*(const Jet3<Scalar>& a, double s) {
    Jet3<Scalar> r;
    const Scalar c = scalar_cast<Scalar>(s);
    r.v = a.v * c;
    for (int i = 0; i < 4; ++i) {
        r.d1[i] = a.d1[i] * c;
        for (int j = 0; j < 4; ++j) {
            r.d2[i][j] = a.d2[i][j] * c;
            for (int k = 0; k < 4; ++k) r.d3[i][j][k] = a.d3[i][j][k] * c;
        }
    }
    return r;
}

// exp, by the chain rule with f' = f. Writing E for exp(u):
//
//     E_i    = E u_i
//     E_ij   = E (u_ij + u_i u_j)
//     E_ijk  = E (u_ijk + u_ij u_k + u_ik u_j + u_jk u_i + u_i u_j u_k)
//
// This is the single place the lapse's positivity becomes structural: exp of a
// finite real is positive, so alpha > 0 needs no constraint and no test.
template <typename Scalar>
Jet3<Scalar> exp_jet(const Jet3<Scalar>& u) {
    using std::exp;
    Jet3<Scalar> r;
    const Scalar E = exp(u.v);
    r.v = E;
    for (int i = 0; i < 4; ++i) {
        r.d1[i] = E * u.d1[i];
        for (int j = 0; j < 4; ++j) {
            r.d2[i][j] = E * (u.d2[i][j] + u.d1[i] * u.d1[j]);
            for (int k = 0; k < 4; ++k)
                r.d3[i][j][k] =
                    E * (u.d3[i][j][k] + u.d2[i][j] * u.d1[k] +
                         u.d2[i][k] * u.d1[j] + u.d2[j][k] * u.d1[i] +
                         u.d1[i] * u.d1[j] * u.d1[k]);
        }
    }
    return r;
}

// The coordinate functions themselves, as jets. axis: 0=t, 1=x, 2=y, 3=z.
template <typename Scalar>
Jet3<Scalar> coordinate_jet(int axis, double value) {
    Jet3<Scalar> j;
    j.v = scalar_cast<Scalar>(value);
    j.d1[axis] = scalar_cast<Scalar>(1.0);
    return j;
}

// The worst deviation from full index symmetry, measured rather than assumed.
template <typename Scalar>
double jet_symmetry_residual(const Jet3<Scalar>& a) {
    double worst = 0.0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            worst = std::fmax(worst, real_abs(a.d2[i][j] - a.d2[j][i]));
            for (int k = 0; k < 4; ++k) {
                const Scalar ref = a.d3[i][j][k];
                worst = std::fmax(worst, real_abs(ref - a.d3[i][k][j]));
                worst = std::fmax(worst, real_abs(ref - a.d3[j][i][k]));
                worst = std::fmax(worst, real_abs(ref - a.d3[j][k][i]));
                worst = std::fmax(worst, real_abs(ref - a.d3[k][i][j]));
                worst = std::fmax(worst, real_abs(ref - a.d3[k][j][i]));
            }
        }
    return worst;
}

}  // namespace design
}  // namespace warp_hypad
