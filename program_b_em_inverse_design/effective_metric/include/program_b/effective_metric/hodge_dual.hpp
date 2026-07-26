#pragma once

#include <cmath>

#include "warp_hypad/core/mat4.hpp"
#include "warp_hypad/core/tensor.hpp"

// ---------------------------------------------------------------------------
// Four-dimensional Hodge duals, for the Bel extraction of Schuster & Visser
// equations (A7a)-(A7d).
//
// CONVENTION, fixed here and documented in
// docs/metric_constitutive_conventions.md section 11:
//
//     eps_{0123} = + sqrt(-det g_lab)          (fully covariant)
//     eps^{0123} = - 1 / sqrt(-det g_lab)      (fully contravariant)
//
// The contravariant sign follows from raising four indices in the (-,+,+,+)
// signature: for Minkowski, eps^{0123} = g^{00}g^{11}g^{22}g^{33} eps_{0123}
// = (-1)(+1)(+1)(+1)(+1) = -1. It is NOT a free choice once the covariant
// convention and the signature are fixed.
//
// The Levi-Civita TENSOR is built from the laboratory metric at every point.
// The bare permutation SYMBOL appears only inside these two constructors --
// hard-coding a Minkowski symbol into the dual operations themselves would
// silently give wrong answers for any curved or non-Cartesian laboratory
// background, which is exactly the generality the covariant formulation exists
// to provide.
//
// Duals of a tensor carrying the antisymmetries of Z (antisymmetric in the
// first index pair and in the second):
//
//     ( *Z)^{abcd} = (1/2) eps^{ab}_{~~ef} Z^{efcd}        left dual
//     ( Z*)^{abcd} = (1/2) Z^{abef} eps_{ef}^{~~cd}        right dual
//     (*Z*)^{abcd} = left then right
//
// In four-dimensional Lorentzian signature the dual of a two-form pair squares
// to MINUS the identity:
//
//     **F = -F,      **Z = -Z  (on either index pair)
//
// which follows from eps_{abcd} eps^{cdef} = -2(delta_a^e delta_b^f
// - delta_a^f delta_b^e). That minus sign is the one thing most easily lost,
// so it is asserted directly rather than assumed.
// ---------------------------------------------------------------------------

namespace program_b {

using warp_hypad::Mat4;
using warp_hypad::Rank4;
using warp_hypad::real_part;
using warp_hypad::scalar_cast;

// Permutation symbol with [0123] = +1. A symbol, not a tensor: it has no
// metric dependence and never appears outside the two constructors below.
inline int permutation_symbol(int a, int b, int c, int d) {
    if (a == b || a == c || a == d || b == c || b == d || c == d) return 0;
    int idx[4] = {a, b, c, d};
    int sign = 1;
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j) {
            if (idx[i] > idx[j]) {
                const int t = idx[i];
                idx[i] = idx[j];
                idx[j] = t;
                sign = -sign;
            }
        }
    return sign;
}

// eps_{abcd} = + sqrt(-det g) [abcd]
template <typename Scalar>
Rank4<Scalar> levi_civita_lower(const Mat4<Scalar>& g_lab) {
    using std::sqrt;
    const Scalar root = sqrt(-warp_hypad::det(g_lab));
    Rank4<Scalar> e;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d) {
                    const int s = permutation_symbol(a, b, c, d);
                    if (s != 0)
                        e(a, b, c, d) =
                            scalar_cast<Scalar>(static_cast<double>(s)) * root;
                }
    return e;
}

// eps^{abcd} = - [abcd] / sqrt(-det g)
template <typename Scalar>
Rank4<Scalar> levi_civita_upper(const Mat4<Scalar>& g_lab) {
    using std::sqrt;
    const Scalar inv = scalar_cast<Scalar>(-1.0) / sqrt(-warp_hypad::det(g_lab));
    Rank4<Scalar> e;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d) {
                    const int s = permutation_symbol(a, b, c, d);
                    if (s != 0)
                        e(a, b, c, d) =
                            scalar_cast<Scalar>(static_cast<double>(s)) * inv;
                }
    return e;
}

// eps^{ab}_{~~ef} = eps^{abcd} g_{ce} g_{df}
template <typename Scalar>
Rank4<Scalar> levi_civita_up_up_dn_dn(const Mat4<Scalar>& g_lab) {
    const Rank4<Scalar> eu = levi_civita_upper(g_lab);
    Rank4<Scalar> out;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int e = 0; e < 4; ++e)
                for (int f = 0; f < 4; ++f) {
                    Scalar acc = scalar_cast<Scalar>(0.0);
                    for (int c = 0; c < 4; ++c)
                        for (int d = 0; d < 4; ++d)
                            acc = acc + eu(a, b, c, d) * g_lab(c, e) * g_lab(d, f);
                    out(a, b, e, f) = acc;
                }
    return out;
}

// eps_{ef}^{~~cd} = eps_{efgh} g^{gc} g^{hd}
template <typename Scalar>
Rank4<Scalar> levi_civita_dn_dn_up_up(const Mat4<Scalar>& g_lab) {
    const Rank4<Scalar> el = levi_civita_lower(g_lab);
    const warp_hypad::InverseResult<Scalar> gi = warp_hypad::inverse(g_lab);
    Rank4<Scalar> out;
    if (!gi.ok) return out;
    for (int e = 0; e < 4; ++e)
        for (int f = 0; f < 4; ++f)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d) {
                    Scalar acc = scalar_cast<Scalar>(0.0);
                    for (int g = 0; g < 4; ++g)
                        for (int h = 0; h < 4; ++h)
                            acc = acc + el(e, f, g, h) * gi.inverse(g, c) *
                                            gi.inverse(h, d);
                    out(e, f, c, d) = acc;
                }
    return out;
}

// ( *Z)^{abcd} = (1/2) eps^{ab}_{~~ef} Z^{efcd}
template <typename Scalar>
Rank4<Scalar> dual_left(const Rank4<Scalar>& Z, const Mat4<Scalar>& g_lab) {
    const Rank4<Scalar> e = levi_civita_up_up_dn_dn(g_lab);
    const Scalar half = scalar_cast<Scalar>(0.5);
    Rank4<Scalar> out;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d) {
                    Scalar acc = scalar_cast<Scalar>(0.0);
                    for (int ee = 0; ee < 4; ++ee)
                        for (int f = 0; f < 4; ++f)
                            acc = acc + e(a, b, ee, f) * Z(ee, f, c, d);
                    out(a, b, c, d) = half * acc;
                }
    return out;
}

// ( Z*)^{abcd} = (1/2) Z^{abef} eps_{ef}^{~~cd}
template <typename Scalar>
Rank4<Scalar> dual_right(const Rank4<Scalar>& Z, const Mat4<Scalar>& g_lab) {
    const Rank4<Scalar> e = levi_civita_dn_dn_up_up(g_lab);
    const Scalar half = scalar_cast<Scalar>(0.5);
    Rank4<Scalar> out;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d) {
                    Scalar acc = scalar_cast<Scalar>(0.0);
                    for (int ee = 0; ee < 4; ++ee)
                        for (int f = 0; f < 4; ++f)
                            acc = acc + Z(a, b, ee, f) * e(ee, f, c, d);
                    out(a, b, c, d) = half * acc;
                }
    return out;
}

template <typename Scalar>
Rank4<Scalar> dual_both(const Rank4<Scalar>& Z, const Mat4<Scalar>& g_lab) {
    return dual_right(dual_left(Z, g_lab), g_lab);
}

// Dual of a two-form, for the **F = -F check:
//     (*F)^{ab} = (1/2) eps^{ab}_{~~cd} F^{cd}
template <typename Scalar>
Mat4<Scalar> dual_2form(const Mat4<Scalar>& F_up, const Mat4<Scalar>& g_lab) {
    const Rank4<Scalar> e = levi_civita_up_up_dn_dn(g_lab);
    const warp_hypad::InverseResult<Scalar> gi = warp_hypad::inverse(g_lab);
    Mat4<Scalar> out;
    if (!gi.ok) return out;
    const Scalar half = scalar_cast<Scalar>(0.5);
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d)
                    acc = acc + e(a, b, c, d) * F_up(c, d);
            out(a, b) = half * acc;
        }
    return out;
}

template <typename Scalar>
double rank4_max_abs_diff(const Rank4<Scalar>& a, const Rank4<Scalar>& b) {
    double worst = 0.0;
    for (int i = 0; i < 256; ++i)
        worst = std::fmax(worst, std::fabs(real_part(a.t[i]) - real_part(b.t[i])));
    return worst;
}

template <typename Scalar>
Rank4<Scalar> rank4_negate(const Rank4<Scalar>& a) {
    Rank4<Scalar> o;
    for (int i = 0; i < 256; ++i) o.t[i] = -a.t[i];
    return o;
}

}  // namespace program_b
