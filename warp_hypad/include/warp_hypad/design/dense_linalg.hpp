#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>

// ---------------------------------------------------------------------------
// Small dense linear algebra: one-sided Jacobi SVD and cyclic Jacobi symmetric
// eigendecomposition.
//
// WHY THESE TWO ALGORITHMS AND NOT A LIBRARY.
//
// The project vendors no linear algebra (see the top of CMakeLists.txt: a
// library that could not be instantiated with the custom scalar type would have
// to be removed again later). The matrices here are small -- a constraint matrix
// with a few dozen columns, a Gram matrix of the same size -- so an O(n^3)
// rotation method costs nothing and buys the two properties that actually
// matter for this milestone:
//
//   RANK REVELATION. The constraint null space IS the answer O0A must produce.
//   Gaussian elimination with partial pivoting can report a full-rank matrix
//   that is numerically singular; one-sided Jacobi computes singular values to
//   high relative accuracy even when they span many orders of magnitude, so a
//   dependent constraint row is visible as a small singular value rather than
//   as a wrong null space.
//
//   NO SILENT LOSS OF SMALL EIGENVALUES. Whitening divides by sqrt(lambda). A
//   Cholesky-based whitening of a nearly singular Gram matrix fails or produces
//   enormous entries; the eigendecomposition lets a numerically null mode be
//   IDENTIFIED and DISCARDED with its eigenvalue recorded, which is what
//   section 9 of the O0A specification requires.
//
// Both routines are deterministic: fixed cyclic sweep order, fixed convergence
// threshold, no pivoting on magnitude ties, no randomisation. Two runs on the
// same input agree bit for bit.
//
// This file is plain double throughout, deliberately. Spline geometry, knots,
// constraint matrices and the whitening transform are all double; only field
// COEFFICIENTS carry derivative channels, and no coefficient enters here.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

struct DenseMatrix {
    int rows = 0;
    int cols = 0;
    std::vector<double> a;  // row-major

    DenseMatrix() = default;
    DenseMatrix(int r, int c)
        : rows(r), cols(c),
          a(static_cast<std::size_t>(r) * static_cast<std::size_t>(c), 0.0) {}

    double& operator()(int i, int j) {
        return a[static_cast<std::size_t>(i) * static_cast<std::size_t>(cols) +
                 static_cast<std::size_t>(j)];
    }
    double operator()(int i, int j) const {
        return a[static_cast<std::size_t>(i) * static_cast<std::size_t>(cols) +
                 static_cast<std::size_t>(j)];
    }

    static DenseMatrix identity(int n) {
        DenseMatrix m(n, n);
        for (int i = 0; i < n; ++i) m(i, i) = 1.0;
        return m;
    }
};

inline DenseMatrix multiply(const DenseMatrix& A, const DenseMatrix& B) {
    DenseMatrix C(A.rows, B.cols);
    for (int i = 0; i < A.rows; ++i)
        for (int k = 0; k < A.cols; ++k) {
            const double aik = A(i, k);
            if (aik == 0.0) continue;
            for (int j = 0; j < B.cols; ++j) C(i, j) += aik * B(k, j);
        }
    return C;
}

inline DenseMatrix transpose(const DenseMatrix& A) {
    DenseMatrix T(A.cols, A.rows);
    for (int i = 0; i < A.rows; ++i)
        for (int j = 0; j < A.cols; ++j) T(j, i) = A(i, j);
    return T;
}

inline std::vector<double> matvec(const DenseMatrix& A,
                                  const std::vector<double>& x) {
    std::vector<double> y(static_cast<std::size_t>(A.rows), 0.0);
    for (int i = 0; i < A.rows; ++i) {
        double s = 0.0;
        for (int j = 0; j < A.cols; ++j)
            s += A(i, j) * x[static_cast<std::size_t>(j)];
        y[static_cast<std::size_t>(i)] = s;
    }
    return y;
}

inline double max_abs(const DenseMatrix& A) {
    double m = 0.0;
    for (double v : A.a) m = std::fmax(m, std::fabs(v));
    return m;
}

inline double norm2(const std::vector<double>& x) {
    double s = 0.0;
    for (double v : x) s += v * v;
    return std::sqrt(s);
}

// --- Deterministic sign convention ------------------------------------------
//
// A singular vector and an eigenvector are each defined only up to sign, and a
// rotation method's output sign depends on the sweep order. Two runs of the same
// code agree, but a figure that plots "mode 7" is otherwise at the mercy of an
// arbitrary choice, and a reader comparing two builds cannot tell a sign flip
// from a change in the physics.
//
// The convention: the largest-magnitude component of each column is made
// POSITIVE, and when two components tie in magnitude the LOWEST INDEX decides.
// Applied to every column of every factor this file returns.
//
// This fixes a sign, not a subspace. Within a degenerate eigenvalue group the
// individual vectors remain arbitrary -- any orthonormal rotation of the group
// is an equally valid answer -- so degeneracy is detected and reported
// separately, and no caller may present a single vector from such a group as
// uniquely determined.
inline void canonicalize_columns(DenseMatrix& m) {
    for (int j = 0; j < m.cols; ++j) {
        int best = 0;
        double best_abs = -1.0;
        for (int i = 0; i < m.rows; ++i) {
            const double a = std::fabs(m(i, j));
            if (a > best_abs) {  // strict: a tie keeps the lower index
                best_abs = a;
                best = i;
            }
        }
        if (m(best, j) < 0.0)
            for (int i = 0; i < m.rows; ++i) m(i, j) = -m(i, j);
    }
}

// Indices of columns that share an eigenvalue to within `rel_tol`, grouped. A
// group of size one is not degenerate; a group of size >= 2 is, and its members
// are individually arbitrary.
inline std::vector<std::vector<int>> degenerate_groups(
    const std::vector<double>& eigenvalues, double rel_tol = 1e-10) {
    std::vector<std::vector<int>> groups;
    if (eigenvalues.empty()) return groups;
    double scale = 0.0;
    for (double v : eigenvalues) scale = std::fmax(scale, std::fabs(v));
    const double tol = rel_tol * (scale > 0.0 ? scale : 1.0);
    std::vector<int> current{0};
    for (std::size_t i = 1; i < eigenvalues.size(); ++i) {
        if (std::fabs(eigenvalues[i] - eigenvalues[i - 1]) <= tol) {
            current.push_back(static_cast<int>(i));
        } else {
            if (current.size() > 1) groups.push_back(current);
            current.assign(1, static_cast<int>(i));
        }
    }
    if (current.size() > 1) groups.push_back(current);
    return groups;
}

// --- One-sided Jacobi SVD ---------------------------------------------------

struct SvdResult {
    // Columns of `v` are the right singular vectors, ordered to match
    // `singular_values` (descending) and sign-canonicalized. The null space of
    // the input is the block of columns from `rank` onward.
    DenseMatrix v;
    std::vector<double> singular_values;
    int rank = 0;
    double tolerance = 0.0;
    int sweeps = 0;
    bool converged = false;

    // U's j-th column, for singular values above the tolerance. Recomputed on
    // demand rather than stored: only the pseudo-inverse needs it.
    // (Left singular vectors: u_j = A v_j / sigma_j.)
};

// A = U S V^T, by orthogonalising the columns of A with plane rotations.
//
// `relative_tolerance` scales the rank threshold; 0 selects the conventional
// max(m,n) * eps * sigma_max.
inline SvdResult jacobi_svd(const DenseMatrix& A_in,
                            double relative_tolerance = 0.0) {
    SvdResult r;
    const int m = A_in.rows, n = A_in.cols;
    r.v = DenseMatrix::identity(n);
    r.singular_values.assign(static_cast<std::size_t>(n), 0.0);
    if (m == 0 || n == 0) {
        r.converged = true;
        return r;
    }

    DenseMatrix B = A_in;
    const double eps = 2.220446049250313e-16;
    constexpr int kMaxSweeps = 60;

    for (int sweep = 0; sweep < kMaxSweeps; ++sweep) {
        // The largest column norm sets the scale below which a column carries no
        // information. WITHOUT THIS FLOOR THE METHOD DOES NOT CONVERGE on a
        // rank-deficient matrix: once a rotation has driven a column down to
        // rounding level, its inner product with a surviving column is pure
        // cancellation noise, and |gamma|/sqrt(alpha*beta) for that pair stays
        // O(1) forever. A rank-1 constraint matrix ran all 60 sweeps and returned
        // converged = false -- correct singular values and null space, but no
        // convergence guarantee behind them. Negligible columns are already
        // orthogonal to everything for every purpose this file has, so they are
        // skipped rather than chased.
        double col_max2 = 0.0;
        for (int j = 0; j < n; ++j) {
            double c2 = 0.0;
            for (int i = 0; i < m; ++i) c2 += B(i, j) * B(i, j);
            col_max2 = std::fmax(col_max2, c2);
        }
        const double negligible = col_max2 * eps * eps;

        double worst = 0.0;
        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                double alpha = 0.0, beta = 0.0, gamma = 0.0;
                for (int i = 0; i < m; ++i) {
                    const double bp = B(i, p), bq = B(i, q);
                    alpha += bp * bp;
                    beta += bq * bq;
                    gamma += bp * bq;
                }
                if (alpha <= negligible || beta <= negligible) continue;
                const double scale = std::sqrt(alpha * beta);
                if (scale == 0.0) continue;
                const double ortho = std::fabs(gamma) / scale;
                if (ortho <= eps) continue;
                worst = std::fmax(worst, ortho);

                // The rotation that annihilates gamma exactly.
                const double zeta = (beta - alpha) / (2.0 * gamma);
                const double sign = (zeta >= 0.0) ? 1.0 : -1.0;
                const double t = sign / (std::fabs(zeta) +
                                         std::sqrt(1.0 + zeta * zeta));
                const double c = 1.0 / std::sqrt(1.0 + t * t);
                const double s = c * t;

                for (int i = 0; i < m; ++i) {
                    const double bp = B(i, p), bq = B(i, q);
                    B(i, p) = c * bp - s * bq;
                    B(i, q) = s * bp + c * bq;
                }
                for (int i = 0; i < n; ++i) {
                    const double vp = r.v(i, p), vq = r.v(i, q);
                    r.v(i, p) = c * vp - s * vq;
                    r.v(i, q) = s * vp + c * vq;
                }
            }
        }
        r.sweeps = sweep + 1;
        if (worst <= eps * static_cast<double>(n)) {
            r.converged = true;
            break;
        }
    }

    std::vector<double> sigma(static_cast<std::size_t>(n), 0.0);
    for (int j = 0; j < n; ++j) {
        double s = 0.0;
        for (int i = 0; i < m; ++i) s += B(i, j) * B(i, j);
        sigma[static_cast<std::size_t>(j)] = std::sqrt(s);
    }

    // Descending order, with a stable tie-break on the original column index so
    // the ordering is reproducible when two singular values coincide.
    std::vector<int> order(static_cast<std::size_t>(n));
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int x, int y) {
        return sigma[static_cast<std::size_t>(x)] >
               sigma[static_cast<std::size_t>(y)];
    });

    DenseMatrix v_sorted(n, n);
    for (int j = 0; j < n; ++j) {
        const int src = order[static_cast<std::size_t>(j)];
        r.singular_values[static_cast<std::size_t>(j)] =
            sigma[static_cast<std::size_t>(src)];
        for (int i = 0; i < n; ++i) v_sorted(i, j) = r.v(i, src);
    }
    r.v = v_sorted;
    canonicalize_columns(r.v);

    const double smax = r.singular_values.empty() ? 0.0 : r.singular_values[0];
    const double factor = (relative_tolerance > 0.0)
                              ? relative_tolerance
                              : static_cast<double>(std::max(m, n)) * eps;
    r.tolerance = factor * smax;
    r.rank = 0;
    for (double s : r.singular_values)
        if (s > r.tolerance) ++r.rank;
    return r;
}

// The minimum-norm solution of A x = b, through the SVD. Components along
// singular directions at or below the tolerance are DROPPED rather than
// divided by: that is the difference between a least-squares answer and a
// division by numerical noise.
inline std::vector<double> svd_min_norm_solve(const DenseMatrix& A,
                                             const SvdResult& svd,
                                             const std::vector<double>& b) {
    const int m = A.rows, n = A.cols;
    std::vector<double> x(static_cast<std::size_t>(n), 0.0);
    for (int j = 0; j < svd.rank; ++j) {
        const double sj = svd.singular_values[static_cast<std::size_t>(j)];
        // u_j = A v_j / sigma_j, formed on the fly.
        double dot = 0.0;
        for (int i = 0; i < m; ++i) {
            double av = 0.0;
            for (int k = 0; k < n; ++k) av += A(i, k) * svd.v(k, j);
            dot += (av / sj) * b[static_cast<std::size_t>(i)];
        }
        const double coeff = dot / sj;
        for (int k = 0; k < n; ++k)
            x[static_cast<std::size_t>(k)] += coeff * svd.v(k, j);
    }
    return x;
}

// --- Cyclic Jacobi symmetric eigendecomposition -----------------------------

struct EigenResult {
    DenseMatrix q;                    // columns are eigenvectors
    std::vector<double> eigenvalues;  // descending
    int sweeps = 0;
    bool converged = false;
};

// S = Q L Q^T for symmetric S. The input's strict lower triangle is ignored;
// the upper triangle defines the matrix, so a caller that assembled a Gram
// matrix asymmetrically cannot have that asymmetry silently averaged away.
inline EigenResult jacobi_eigen(const DenseMatrix& S_in) {
    EigenResult r;
    const int n = S_in.rows;
    r.q = DenseMatrix::identity(n);
    r.eigenvalues.assign(static_cast<std::size_t>(n), 0.0);
    if (n == 0) {
        r.converged = true;
        return r;
    }

    DenseMatrix S(n, n);
    for (int i = 0; i < n; ++i) {
        S(i, i) = S_in(i, i);
        for (int j = i + 1; j < n; ++j) {
            S(i, j) = S_in(i, j);
            S(j, i) = S_in(i, j);
        }
    }

    const double eps = 2.220446049250313e-16;
    constexpr int kMaxSweeps = 100;
    double scale = 0.0;
    for (int i = 0; i < n; ++i) scale = std::fmax(scale, std::fabs(S(i, i)));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) scale = std::fmax(scale, std::fabs(S(i, j)));
    if (scale == 0.0) {
        r.converged = true;
        return r;
    }

    for (int sweep = 0; sweep < kMaxSweeps; ++sweep) {
        double off = 0.0;
        for (int p = 0; p < n - 1; ++p)
            for (int q = p + 1; q < n; ++q) off += S(p, q) * S(p, q);
        r.sweeps = sweep;
        if (std::sqrt(off) <= eps * scale * static_cast<double>(n)) {
            r.converged = true;
            break;
        }
        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                const double apq = S(p, q);
                if (std::fabs(apq) <= eps * scale) continue;
                const double theta = (S(q, q) - S(p, p)) / (2.0 * apq);
                const double sign = (theta >= 0.0) ? 1.0 : -1.0;
                const double t = sign / (std::fabs(theta) +
                                         std::sqrt(1.0 + theta * theta));
                const double c = 1.0 / std::sqrt(1.0 + t * t);
                const double s = c * t;
                for (int k = 0; k < n; ++k) {
                    const double skp = S(k, p), skq = S(k, q);
                    S(k, p) = c * skp - s * skq;
                    S(k, q) = s * skp + c * skq;
                }
                for (int k = 0; k < n; ++k) {
                    const double spk = S(p, k), sqk = S(q, k);
                    S(p, k) = c * spk - s * sqk;
                    S(q, k) = s * spk + c * sqk;
                }
                for (int k = 0; k < n; ++k) {
                    const double qkp = r.q(k, p), qkq = r.q(k, q);
                    r.q(k, p) = c * qkp - s * qkq;
                    r.q(k, q) = s * qkp + c * qkq;
                }
            }
        }
    }

    std::vector<double> lam(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) lam[static_cast<std::size_t>(i)] = S(i, i);
    std::vector<int> order(static_cast<std::size_t>(n));
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int x, int y) {
        return lam[static_cast<std::size_t>(x)] >
               lam[static_cast<std::size_t>(y)];
    });
    DenseMatrix q_sorted(n, n);
    for (int j = 0; j < n; ++j) {
        const int src = order[static_cast<std::size_t>(j)];
        r.eigenvalues[static_cast<std::size_t>(j)] =
            lam[static_cast<std::size_t>(src)];
        for (int i = 0; i < n; ++i) q_sorted(i, j) = r.q(i, src);
    }
    r.q = q_sorted;
    canonicalize_columns(r.q);
    return r;
}

// --- Reconstruction residuals, for the independent fixtures -----------------

// max |A V_null|: the null-space columns must be annihilated by A.
inline double null_space_residual(const DenseMatrix& A, const SvdResult& svd) {
    const int m = A.rows, n = A.cols;
    double worst = 0.0;
    for (int j = svd.rank; j < n; ++j)
        for (int i = 0; i < m; ++i) {
            double acc = 0.0;
            for (int k = 0; k < n; ++k) acc += A(i, k) * svd.v(k, j);
            worst = std::fmax(worst, std::fabs(acc));
        }
    return worst;
}

// max |V^T V - I| over the whole factor, so a loss of orthogonality anywhere is
// visible and not only in the retained block.
inline double orthonormality_residual(const DenseMatrix& V) {
    const DenseMatrix G = multiply(transpose(V), V);
    double worst = 0.0;
    for (int i = 0; i < G.rows; ++i)
        for (int j = 0; j < G.cols; ++j)
            worst = std::fmax(worst, std::fabs(G(i, j) - (i == j ? 1.0 : 0.0)));
    return worst;
}

// max |A A^+ A - A|, the Moore-Penrose identity. This is the check that a
// rank-deficient pseudo-inverse is the right one rather than merely finite.
inline double pseudo_inverse_residual(const DenseMatrix& A,
                                      const SvdResult& svd) {
    const int m = A.rows, n = A.cols;
    // A^+ = sum_{j<rank} v_j u_j^T / sigma_j, u_j = A v_j / sigma_j.
    DenseMatrix Aplus(n, m);
    for (int j = 0; j < svd.rank; ++j) {
        const double sj = svd.singular_values[static_cast<std::size_t>(j)];
        std::vector<double> u(static_cast<std::size_t>(m), 0.0);
        for (int i = 0; i < m; ++i) {
            double acc = 0.0;
            for (int k = 0; k < n; ++k) acc += A(i, k) * svd.v(k, j);
            u[static_cast<std::size_t>(i)] = acc / sj;
        }
        for (int k = 0; k < n; ++k)
            for (int i = 0; i < m; ++i)
                Aplus(k, i) += svd.v(k, j) * u[static_cast<std::size_t>(i)] / sj;
    }
    const DenseMatrix AAA = multiply(multiply(A, Aplus), A);
    double worst = 0.0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            worst = std::fmax(worst, std::fabs(AAA(i, j) - A(i, j)));
    return worst;
}

// max |Q L Q^T - S|, the eigendecomposition identity.
inline double eigen_reconstruction_residual(const DenseMatrix& S,
                                           const EigenResult& e) {
    const int n = S.rows;
    double worst = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double acc = 0.0;
            for (int k = 0; k < n; ++k)
                acc += e.q(i, k) *
                       e.eigenvalues[static_cast<std::size_t>(k)] * e.q(j, k);
            worst = std::fmax(worst, std::fabs(acc - S(i, j)));
        }
    return worst;
}

// Solve S x = b for symmetric positive semidefinite S, dropping directions at
// or below `rel_tol` times the largest eigenvalue. Used for the constrained
// least-squares projection, where the normal matrix can be singular if the
// design region leaves a basis direction unsampled.
inline std::vector<double> spd_solve(const DenseMatrix& S,
                                     const std::vector<double>& b,
                                     double rel_tol = 1e-13,
                                     int* dropped = nullptr) {
    const EigenResult e = jacobi_eigen(S);
    const int n = S.rows;
    const double lmax = e.eigenvalues.empty() ? 0.0 : e.eigenvalues[0];
    const double cut = rel_tol * std::fabs(lmax);
    std::vector<double> x(static_cast<std::size_t>(n), 0.0);
    int drop = 0;
    for (int j = 0; j < n; ++j) {
        const double lj = e.eigenvalues[static_cast<std::size_t>(j)];
        if (!(lj > cut)) { ++drop; continue; }
        double dot = 0.0;
        for (int i = 0; i < n; ++i)
            dot += e.q(i, j) * b[static_cast<std::size_t>(i)];
        const double coeff = dot / lj;
        for (int i = 0; i < n; ++i)
            x[static_cast<std::size_t>(i)] += coeff * e.q(i, j);
    }
    if (dropped) *dropped = drop;
    return x;
}

}  // namespace design
}  // namespace warp_hypad
