#pragma once

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

// ---------------------------------------------------------------------------
// Independent reference linear algebra, for Level D verification.
//
// This header exists to disagree with mat4.hpp. It must therefore share no code
// with it and, more importantly, no ALGORITHM with it:
//
//   * mat4.hpp inverts via the adjugate (closed-form cofactors, branch-free,
//     chosen so a hypercomplex AD scalar can flow through it).
//     Here: Gauss-Jordan elimination with partial pivoting, which branches
//     freely because it only ever runs on double.
//
//   * mat4.hpp computes the determinant by pairing 2x2 minors from rows {0,1}
//     and {2,3}.
//     Here: recursive Laplace expansion along the first row, an entirely
//     different summation order and grouping.
//
//   * mat4.hpp multiplies with an accumulate-into-Scalar inner loop.
//     Here: the same mathematical operation but written independently, in
//     plain double, indexed differently.
//
// If both implementations agree to 1e-12 on thousands of random matrices, the
// probability that they share a bug is small. If they shared an algorithm, the
// agreement would prove nothing at all.
//
// Everything here is plain double on purpose. A reference that was itself
// templated on Scalar could not catch a genericity bug in the code under test.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace reference {

// Row-major 4x4 in a flat array, deliberately a different container type from
// Mat4 so the two cannot be accidentally interchanged.
struct RefMat4 {
    double a[4][4];

    double& operator()(int r, int c) { return a[r][c]; }
    double operator()(int r, int c) const { return a[r][c]; }
};

inline RefMat4 ref_zero() {
    RefMat4 m;
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c) m.a[r][c] = 0.0;
    return m;
}

inline RefMat4 ref_identity() {
    RefMat4 m = ref_zero();
    for (int i = 0; i < 4; ++i) m.a[i][i] = 1.0;
    return m;
}

// Naive triple loop, written independently of mat_mul.
inline RefMat4 ref_mul(const RefMat4& x, const RefMat4& y) {
    RefMat4 out = ref_zero();
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
            double acc = 0.0;
            for (int t = 0; t < 4; ++t) acc += x.a[r][t] * y.a[t][c];
            out.a[r][c] = acc;
        }
    }
    return out;
}

// Determinant of an n x n matrix by recursive Laplace expansion along row 0.
// Exponential in n, which is irrelevant at n = 4 and keeps the code obviously
// correct by inspection -- the point of a reference implementation.
inline double ref_det_recursive(const std::vector<std::vector<double>>& m) {
    const std::size_t n = m.size();
    if (n == 1) return m[0][0];
    if (n == 2) return m[0][0] * m[1][1] - m[0][1] * m[1][0];

    double acc = 0.0;
    double sign = 1.0;
    for (std::size_t col = 0; col < n; ++col) {
        std::vector<std::vector<double>> minor;
        minor.reserve(n - 1);
        for (std::size_t r = 1; r < n; ++r) {
            std::vector<double> row;
            row.reserve(n - 1);
            for (std::size_t c = 0; c < n; ++c) {
                if (c == col) continue;
                row.push_back(m[r][c]);
            }
            minor.push_back(std::move(row));
        }
        acc += sign * m[0][col] * ref_det_recursive(minor);
        sign = -sign;
    }
    return acc;
}

inline double ref_det(const RefMat4& m) {
    std::vector<std::vector<double>> v(4, std::vector<double>(4));
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c) v[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] = m.a[r][c];
    return ref_det_recursive(v);
}

// Gauss-Jordan elimination with partial pivoting. Returns false when a pivot
// column is numerically empty, which is a completely different singularity
// criterion from the determinant-magnitude test in mat4.hpp -- so the two can
// disagree, which is what makes the comparison informative.
inline bool ref_inverse(const RefMat4& in, RefMat4& out) {
    double a[4][8];
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) a[r][c] = in.a[r][c];
        for (int c = 0; c < 4; ++c) a[r][4 + c] = (r == c) ? 1.0 : 0.0;
    }

    for (int col = 0; col < 4; ++col) {
        int pivot = col;
        double best = std::fabs(a[col][col]);
        for (int r = col + 1; r < 4; ++r) {
            const double v = std::fabs(a[r][col]);
            if (v > best) { best = v; pivot = r; }
        }
        if (!(best > 0.0) || !std::isfinite(best)) return false;

        if (pivot != col) {
            for (int c = 0; c < 8; ++c) {
                const double t = a[col][c];
                a[col][c] = a[pivot][c];
                a[pivot][c] = t;
            }
        }

        const double inv_pivot = 1.0 / a[col][col];
        for (int c = 0; c < 8; ++c) a[col][c] *= inv_pivot;

        for (int r = 0; r < 4; ++r) {
            if (r == col) continue;
            const double factor = a[r][col];
            if (factor == 0.0) continue;
            for (int c = 0; c < 8; ++c) a[r][c] -= factor * a[col][c];
        }
    }

    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c) out.a[r][c] = a[r][4 + c];
    return true;
}

// --- Deterministic test matrix generation ---------------------------------
//
// Seeded explicitly so that a failing case can be reproduced exactly from the
// verification report. An unseeded RNG in a verification suite makes failures
// unreproducible, which defeats the purpose.

class MatrixGenerator {
public:
    explicit MatrixGenerator(std::uint64_t seed) : rng_(seed) {}

    RefMat4 general() {
        std::uniform_real_distribution<double> d(-2.0, 2.0);
        RefMat4 m;
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c) m.a[r][c] = d(rng_);
        return m;
    }

    RefMat4 symmetric() {
        RefMat4 m = general();
        for (int r = 0; r < 4; ++r)
            for (int c = r + 1; c < 4; ++c) m.a[c][r] = m.a[r][c];
        return m;
    }

    RefMat4 diagonal() {
        std::uniform_real_distribution<double> d(0.5, 3.0);
        std::bernoulli_distribution flip(0.5);
        RefMat4 m = ref_zero();
        for (int i = 0; i < 4; ++i) m.a[i][i] = flip(rng_) ? -d(rng_) : d(rng_);
        return m;
    }

    // A perturbed Lorentzian metric: signature (-,+,+,+) with small symmetric
    // off-diagonal terms, i.e. the shape the Alcubierre metric actually takes.
    RefMat4 lorentzian_like() {
        std::uniform_real_distribution<double> off(-0.3, 0.3);
        RefMat4 m = ref_zero();
        m.a[0][0] = -1.0 + off(rng_) * 0.5;
        for (int i = 1; i < 4; ++i) m.a[i][i] = 1.0 + off(rng_) * 0.2;
        for (int r = 0; r < 4; ++r) {
            for (int c = r + 1; c < 4; ++c) {
                const double v = off(rng_);
                m.a[r][c] = v;
                m.a[c][r] = v;
            }
        }
        return m;
    }

    // Exactly singular: row 3 is a linear combination of rows 0 and 1.
    RefMat4 singular() {
        RefMat4 m = general();
        for (int c = 0; c < 4; ++c) m.a[3][c] = 2.0 * m.a[0][c] - 3.0 * m.a[1][c];
        return m;
    }

    // Nearly singular: as above, plus a perturbation of size `eps`.
    RefMat4 near_singular(double eps) {
        RefMat4 m = singular();
        std::uniform_real_distribution<double> d(-1.0, 1.0);
        for (int c = 0; c < 4; ++c) m.a[3][c] += eps * d(rng_);
        return m;
    }

private:
    std::mt19937_64 rng_;
};

}  // namespace reference
}  // namespace warp_hypad
