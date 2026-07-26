#pragma once

#include "warp_hypad/core/mat4.hpp"
#include "warp_hypad/core/scalar_traits.hpp"

// ---------------------------------------------------------------------------
// Fixed-size rank-3 and rank-4 spacetime tensors, plus the specific
// contractions the GR pipeline needs.
//
// This is deliberately NOT a general tensor library. Only the ranks and the
// contractions actually required by Milestones 1-8 exist here. Adding a
// general einsum-style facility would buy nothing and would obscure the index
// placement, which is the thing most likely to be wrong in a GR code.
//
// Storage: row-major, last index fastest, matching Mat4.
//
//     Rank3 t;  t(a,b,c)   -> data[16*a + 4*b + c]      (64 components)
//     Rank4 t;  t(a,b,c,d) -> data[64*a + 16*b + 4*c + d]  (256 components)
//
// Index-position discipline (see docs/conventions.md): these containers do not
// know which slots are contravariant. The naming convention in the physics code
// is that the leading slot(s) are the raised ones, matching how the objects are
// written on paper:
//
//     Rank3 christoffel;  christoffel(rho, mu, nu)  <->  Gamma^rho_{mu nu}
//     Rank4 riemann;      riemann(rho, sig, mu, nu) <->  R^rho_{sig mu nu}
// ---------------------------------------------------------------------------

namespace warp_hypad {

template <typename Scalar>
struct Rank3 {
    Scalar t[64];

    Rank3() {
        const Scalar z = scalar_cast<Scalar>(0.0);
        for (int i = 0; i < 64; ++i) t[i] = z;
    }

    Scalar& operator()(int a, int b, int c) { return t[16 * a + 4 * b + c]; }
    const Scalar& operator()(int a, int b, int c) const {
        return t[16 * a + 4 * b + c];
    }
};

template <typename Scalar>
struct Rank4 {
    Scalar t[256];

    Rank4() {
        const Scalar z = scalar_cast<Scalar>(0.0);
        for (int i = 0; i < 256; ++i) t[i] = z;
    }

    Scalar& operator()(int a, int b, int c, int d) {
        return t[64 * a + 16 * b + 4 * c + d];
    }
    const Scalar& operator()(int a, int b, int c, int d) const {
        return t[64 * a + 16 * b + 4 * c + d];
    }
};

// --- Contractions ----------------------------------------------------------

// Trace of a mixed rank-2 tensor: M^mu_mu.
template <typename Scalar>
Scalar trace_mixed(const Mat4<Scalar>& m) {
    Scalar s = scalar_cast<Scalar>(0.0);
    for (int i = 0; i < 4; ++i) s = s + m(i, i);
    return s;
}

// Full contraction of two rank-2 tensors with opposite index placement:
//     result = A^{mu nu} B_{mu nu}
// Used for the Ricci scalar R = g^{mu nu} R_{mu nu} and for invariants such as
// R_{mu nu} R^{mu nu}.
template <typename Scalar>
Scalar full_contract(const Mat4<Scalar>& a, const Mat4<Scalar>& b) {
    Scalar s = scalar_cast<Scalar>(0.0);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) s = s + a(i, j) * b(i, j);
    return s;
}

// Raise the first index: result^mu_nu = g^{mu alpha} T_{alpha nu}.
// Identical to mat_mul(g_inv, t); named separately so call sites read as physics.
template <typename Scalar>
Mat4<Scalar> raise_first(const Mat4<Scalar>& g_inv, const Mat4<Scalar>& t) {
    return mat_mul(g_inv, t);
}

// Contract a rank-2 against the FIRST slot of a rank-3:
//     result(r, m, n) = sum_s a(r, s) * d(s, m, n)
// This is exactly the step that turns the Christoffel symbol of the first kind
// into the second kind: Gamma^rho_{mu nu} = (1/2) g^{rho sigma} D_{sigma mu nu}.
template <typename Scalar>
Rank3<Scalar> contract_first(const Mat4<Scalar>& a, const Rank3<Scalar>& d) {
    Rank3<Scalar> out;
    for (int r = 0; r < 4; ++r) {
        for (int m = 0; m < 4; ++m) {
            for (int n = 0; n < 4; ++n) {
                Scalar sum = scalar_cast<Scalar>(0.0);
                for (int s = 0; s < 4; ++s) sum = sum + a(r, s) * d(s, m, n);
                out(r, m, n) = sum;
            }
        }
    }
    return out;
}

// Contract slots 1 and 3 of a rank-4, the Riemann -> Ricci step:
//     result(s, n) = sum_r t(r, s, r, n)
// i.e. R_{sigma nu} = R^rho_{sigma rho nu}.
template <typename Scalar>
Mat4<Scalar> contract_13(const Rank4<Scalar>& t) {
    Mat4<Scalar> out;
    for (int s = 0; s < 4; ++s) {
        for (int n = 0; n < 4; ++n) {
            Scalar sum = scalar_cast<Scalar>(0.0);
            for (int r = 0; r < 4; ++r) sum = sum + t(r, s, r, n);
            out(s, n) = sum;
        }
    }
    return out;
}

// Lower the leading index of a rank-4:
//     result(a, b, c, d) = sum_r g(a, r) * t(r, b, c, d)
// Needed to test the pair symmetry R_{rho sigma mu nu} = R_{mu nu rho sigma},
// which only holds once all four indices are down.
template <typename Scalar>
Rank4<Scalar> lower_first(const Mat4<Scalar>& g, const Rank4<Scalar>& t) {
    Rank4<Scalar> out;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d) {
                    Scalar sum = scalar_cast<Scalar>(0.0);
                    for (int r = 0; r < 4; ++r) sum = sum + g(a, r) * t(r, b, c, d);
                    out(a, b, c, d) = sum;
                }
    return out;
}

// --- Residual helpers ------------------------------------------------------

// Largest |real| violation of symmetry in the trailing index pair:
//     t(a, b, c) == t(a, c, b)
// The Levi-Civita connection satisfies this identically.
template <typename Scalar>
double lower_pair_symmetry_residual(const Rank3<Scalar>& t) {
    double worst = 0.0;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = b + 1; c < 4; ++c) {
                const double d = real_part(t(a, b, c)) - real_part(t(a, c, b));
                const double ad = d < 0.0 ? -d : d;
                if (ad > worst) worst = ad;
            }
    return worst;
}

template <typename Scalar>
double max_abs_real(const Rank3<Scalar>& t) {
    double worst = 0.0;
    for (int i = 0; i < 64; ++i) {
        const double a = real_abs(t.t[i]);
        if (a > worst) worst = a;
    }
    return worst;
}

template <typename Scalar>
double max_abs_real(const Rank4<Scalar>& t) {
    double worst = 0.0;
    for (int i = 0; i < 256; ++i) {
        const double a = real_abs(t.t[i]);
        if (a > worst) worst = a;
    }
    return worst;
}

template <typename Scalar>
double max_abs_real(const Mat4<Scalar>& t) {
    double worst = 0.0;
    for (int i = 0; i < 16; ++i) {
        const double a = real_abs(t.m[i]);
        if (a > worst) worst = a;
    }
    return worst;
}

template <typename Scalar>
double max_abs_diff(const Rank3<Scalar>& a, const Rank3<Scalar>& b) {
    double worst = 0.0;
    for (int i = 0; i < 64; ++i) {
        const double d = real_part(a.t[i]) - real_part(b.t[i]);
        const double ad = d < 0.0 ? -d : d;
        if (ad > worst) worst = ad;
    }
    return worst;
}

}  // namespace warp_hypad
