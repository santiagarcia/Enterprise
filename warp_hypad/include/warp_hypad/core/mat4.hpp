#pragma once

#include <cstddef>

#include "warp_hypad/core/scalar_traits.hpp"

// ---------------------------------------------------------------------------
// Fixed-size 4x4 matrix over a generic scalar type.
//
// Index convention: m(row, col), row-major, so the flat offset is 4*row + col.
// For a spacetime tensor the first index is the first (leftmost) tensor slot:
//
//     Mat4 g;  g(mu, nu)  <->  g_{mu nu}
//     Mat4 h;  h(mu, nu)  <->  h^{mu nu}
//
// The matrix type is index-position agnostic; it stores 16 numbers and nothing
// more. Which slots are up and which are down is tracked by the *caller* and
// documented at every call site. See docs/conventions.md.
//
// Determinant and inverse use the closed-form cofactor (adjugate) expansion
// rather than LU with partial pivoting. This is deliberate: pivoting selects a
// row by comparing magnitudes, which is a branch on a Scalar value, which is
// exactly what a hypercomplex AD type cannot tolerate. The cofactor form is
// branch-free straight-line arithmetic, so it differentiates cleanly and gives
// bit-identical structure for double and HyPADScalar. The cost is ~2x more
// flops than LU for a 4x4 and worse behaviour on ill-conditioned input, which
// is why singularity and conditioning are reported explicitly below.
// ---------------------------------------------------------------------------

namespace warp_hypad {

// Relative floor for the singularity test. |det| is compared against
// kSingularDetRelTol * norm1(A)^4 -- the fourth power makes the test scale
// invariant, since det scales as the fourth power of the matrix entries.
inline constexpr double kSingularDetRelTol = 1e-13;

// Accuracy target for the inverse. An inverse whose identity residual may
// exceed this is reported as ill-conditioned.
inline constexpr double kInverseAccuracyTarget = 1e-8;

// 1-norm condition number above which the inverse is reported as untrustworthy.
//
// Derived, not chosen: the standard bound for the residual of a computed
// inverse is ||A A^-1 - I|| ~ cond(A) * eps. Measurement on deliberately
// near-singular matrices (see algebra.ill_conditioning_flag_fires) gives
// residual ~ cond1 * 6e-17 for this adjugate implementation, i.e. about
// 0.3 * cond1 * eps -- so cond1 * eps is a conservative upper bound, which is
// what a warning flag should use.
//
// Setting the limit to target/eps ~ 4.5e7 therefore means: flag exactly when
// the identity residual could exceed kInverseAccuracyTarget.
//
// An earlier round number of 1e12 was wrong by four orders of magnitude -- it
// corresponded to a tolerated residual of about 6e-5, so matrices whose
// inverses had lost every useful digit were being returned unflagged.
inline constexpr double kConditionWarnLimit =
    kInverseAccuracyTarget / std::numeric_limits<double>::epsilon();

template <typename Scalar>
struct Mat4 {
    Scalar m[16];

    Mat4() {
        const Scalar z = scalar_cast<Scalar>(0.0);
        for (int i = 0; i < 16; ++i) m[i] = z;
    }

    Scalar& operator()(int row, int col) { return m[4 * row + col]; }
    const Scalar& operator()(int row, int col) const { return m[4 * row + col]; }

    static Mat4 zero() { return Mat4(); }

    static Mat4 identity() {
        Mat4 r;
        const Scalar one = scalar_cast<Scalar>(1.0);
        for (int i = 0; i < 4; ++i) r(i, i) = one;
        return r;
    }

    // Minkowski metric in the (-,+,+,+) signature, Cartesian coordinates.
    static Mat4 minkowski() {
        Mat4 r;
        r(0, 0) = scalar_cast<Scalar>(-1.0);
        r(1, 1) = scalar_cast<Scalar>(1.0);
        r(2, 2) = scalar_cast<Scalar>(1.0);
        r(3, 3) = scalar_cast<Scalar>(1.0);
        return r;
    }
};

// C = A * B, i.e. C(i,k) = sum_j A(i,j) B(j,k).
// In index terms this contracts the second slot of A against the first of B.
template <typename Scalar>
Mat4<Scalar> mat_mul(const Mat4<Scalar>& a, const Mat4<Scalar>& b) {
    Mat4<Scalar> c;
    for (int i = 0; i < 4; ++i) {
        for (int k = 0; k < 4; ++k) {
            Scalar sum = scalar_cast<Scalar>(0.0);
            for (int j = 0; j < 4; ++j) sum = sum + a(i, j) * b(j, k);
            c(i, k) = sum;
        }
    }
    return c;
}

template <typename Scalar>
Mat4<Scalar> transpose(const Mat4<Scalar>& a) {
    Mat4<Scalar> r;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) r(i, j) = a(j, i);
    return r;
}

template <typename Scalar>
Mat4<Scalar> mat_sub(const Mat4<Scalar>& a, const Mat4<Scalar>& b) {
    Mat4<Scalar> r;
    for (int i = 0; i < 16; ++i) r.m[i] = a.m[i] - b.m[i];
    return r;
}

template <typename Scalar>
Mat4<Scalar> mat_add(const Mat4<Scalar>& a, const Mat4<Scalar>& b) {
    Mat4<Scalar> r;
    for (int i = 0; i < 16; ++i) r.m[i] = a.m[i] + b.m[i];
    return r;
}

template <typename Scalar>
Mat4<Scalar> mat_scale(const Mat4<Scalar>& a, const Scalar& s) {
    Mat4<Scalar> r;
    for (int i = 0; i < 16; ++i) r.m[i] = a.m[i] * s;
    return r;
}

// The six 2x2 minors drawn from rows {0,1} and the six from rows {2,3}.
// Shared between det() and inverse() so the two can never disagree.
template <typename Scalar>
struct Mat4Minors {
    Scalar s0, s1, s2, s3, s4, s5;  // from rows 0,1
    Scalar c0, c1, c2, c3, c4, c5;  // from rows 2,3
};

template <typename Scalar>
Mat4Minors<Scalar> mat_minors(const Mat4<Scalar>& a) {
    Mat4Minors<Scalar> k;
    k.s0 = a(0, 0) * a(1, 1) - a(1, 0) * a(0, 1);
    k.s1 = a(0, 0) * a(1, 2) - a(1, 0) * a(0, 2);
    k.s2 = a(0, 0) * a(1, 3) - a(1, 0) * a(0, 3);
    k.s3 = a(0, 1) * a(1, 2) - a(1, 1) * a(0, 2);
    k.s4 = a(0, 1) * a(1, 3) - a(1, 1) * a(0, 3);
    k.s5 = a(0, 2) * a(1, 3) - a(1, 2) * a(0, 3);

    k.c5 = a(2, 2) * a(3, 3) - a(3, 2) * a(2, 3);
    k.c4 = a(2, 1) * a(3, 3) - a(3, 1) * a(2, 3);
    k.c3 = a(2, 1) * a(3, 2) - a(3, 1) * a(2, 2);
    k.c2 = a(2, 0) * a(3, 3) - a(3, 0) * a(2, 3);
    k.c1 = a(2, 0) * a(3, 2) - a(3, 0) * a(2, 2);
    k.c0 = a(2, 0) * a(3, 1) - a(3, 0) * a(2, 1);
    return k;
}

template <typename Scalar>
Scalar det(const Mat4<Scalar>& a) {
    const Mat4Minors<Scalar> k = mat_minors(a);
    return k.s0 * k.c5 - k.s1 * k.c4 + k.s2 * k.c3 + k.s3 * k.c2 - k.s4 * k.c1 +
           k.s5 * k.c0;
}

// Maximum absolute column sum of the real parts. Used for conditioning
// diagnostics only, so evaluating on real parts is correct.
template <typename Scalar>
double norm1(const Mat4<Scalar>& a) {
    double best = 0.0;
    for (int col = 0; col < 4; ++col) {
        double sum = 0.0;
        for (int row = 0; row < 4; ++row) sum += real_abs(a(row, col));
        if (sum > best) best = sum;
    }
    return best;
}

template <typename Scalar>
struct InverseResult {
    Mat4<Scalar> inverse;

    // Determinant in full scalar precision -- carries derivative channels.
    Scalar determinant;

    // Real-part diagnostics.
    double det_real = 0.0;
    double norm1_a = 0.0;
    double norm1_inv = 0.0;
    double cond1 = 0.0;

    // ok == false means `inverse` is meaningless and must not be used.
    bool ok = false;

    // ok == true but the result has lost most of its significant digits.
    bool ill_conditioned = false;
};

// Inverse via the adjugate. Never throws; the caller must check .ok.
//
// The singularity test is performed on the real part of the determinant only,
// which is the correct choice for an AD scalar: whether a matrix is invertible
// is a property of the primal value, not of its sensitivities.
template <typename Scalar>
InverseResult<Scalar> inverse(const Mat4<Scalar>& a) {
    InverseResult<Scalar> out;
    const Mat4Minors<Scalar> k = mat_minors(a);

    out.determinant = k.s0 * k.c5 - k.s1 * k.c4 + k.s2 * k.c3 + k.s3 * k.c2 -
                      k.s4 * k.c1 + k.s5 * k.c0;
    out.det_real = real_part(out.determinant);
    out.norm1_a = norm1(a);

    const double scale = out.norm1_a * out.norm1_a * out.norm1_a * out.norm1_a;
    const double det_abs = out.det_real < 0.0 ? -out.det_real : out.det_real;

    const bool finite_det = (out.det_real == out.det_real) &&
                            det_abs != std::numeric_limits<double>::infinity();

    if (!finite_det || det_abs <= kSingularDetRelTol * scale || scale == 0.0) {
        out.ok = false;
        return out;  // `inverse` left zero-initialised
    }

    const Scalar inv_det = scalar_cast<Scalar>(1.0) / out.determinant;
    Mat4<Scalar>& b = out.inverse;

    b(0, 0) = (a(1, 1) * k.c5 - a(1, 2) * k.c4 + a(1, 3) * k.c3) * inv_det;
    b(0, 1) = (-a(0, 1) * k.c5 + a(0, 2) * k.c4 - a(0, 3) * k.c3) * inv_det;
    b(0, 2) = (a(3, 1) * k.s5 - a(3, 2) * k.s4 + a(3, 3) * k.s3) * inv_det;
    b(0, 3) = (-a(2, 1) * k.s5 + a(2, 2) * k.s4 - a(2, 3) * k.s3) * inv_det;

    b(1, 0) = (-a(1, 0) * k.c5 + a(1, 2) * k.c2 - a(1, 3) * k.c1) * inv_det;
    b(1, 1) = (a(0, 0) * k.c5 - a(0, 2) * k.c2 + a(0, 3) * k.c1) * inv_det;
    b(1, 2) = (-a(3, 0) * k.s5 + a(3, 2) * k.s2 - a(3, 3) * k.s1) * inv_det;
    b(1, 3) = (a(2, 0) * k.s5 - a(2, 2) * k.s2 + a(2, 3) * k.s1) * inv_det;

    b(2, 0) = (a(1, 0) * k.c4 - a(1, 1) * k.c2 + a(1, 3) * k.c0) * inv_det;
    b(2, 1) = (-a(0, 0) * k.c4 + a(0, 1) * k.c2 - a(0, 3) * k.c0) * inv_det;
    b(2, 2) = (a(3, 0) * k.s4 - a(3, 1) * k.s2 + a(3, 3) * k.s0) * inv_det;
    b(2, 3) = (-a(2, 0) * k.s4 + a(2, 1) * k.s2 - a(2, 3) * k.s0) * inv_det;

    b(3, 0) = (-a(1, 0) * k.c3 + a(1, 1) * k.c1 - a(1, 2) * k.c0) * inv_det;
    b(3, 1) = (a(0, 0) * k.c3 - a(0, 1) * k.c1 + a(0, 2) * k.c0) * inv_det;
    b(3, 2) = (-a(3, 0) * k.s3 + a(3, 1) * k.s1 - a(3, 2) * k.s0) * inv_det;
    b(3, 3) = (a(2, 0) * k.s3 - a(2, 1) * k.s1 + a(2, 2) * k.s0) * inv_det;

    out.norm1_inv = norm1(b);
    out.cond1 = out.norm1_a * out.norm1_inv;
    out.ok = true;
    out.ill_conditioned = out.cond1 > kConditionWarnLimit;
    return out;
}

// Largest |real| entry of A*B - I. The primary Level A residual for Milestone 1.
template <typename Scalar>
double identity_residual(const Mat4<Scalar>& a, const Mat4<Scalar>& b) {
    const Mat4<Scalar> p = mat_mul(a, b);
    double worst = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            const double target = (i == j) ? 1.0 : 0.0;
            const double d = real_part(p(i, j)) - target;
            const double ad = d < 0.0 ? -d : d;
            if (ad > worst) worst = ad;
        }
    }
    return worst;
}

// Largest |real| entry of A - A^T. Zero for a metric, Ricci or Einstein tensor.
template <typename Scalar>
double symmetry_residual(const Mat4<Scalar>& a) {
    double worst = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            const double d = real_part(a(i, j)) - real_part(a(j, i));
            const double ad = d < 0.0 ? -d : d;
            if (ad > worst) worst = ad;
        }
    }
    return worst;
}

// Largest |real| componentwise difference between two matrices.
template <typename Scalar>
double max_abs_diff(const Mat4<Scalar>& a, const Mat4<Scalar>& b) {
    double worst = 0.0;
    for (int i = 0; i < 16; ++i) {
        const double d = real_part(a.m[i]) - real_part(b.m[i]);
        const double ad = d < 0.0 ? -d : d;
        if (ad > worst) worst = ad;
    }
    return worst;
}

}  // namespace warp_hypad
