#include "warp_hypad/geometry/signature.hpp"

#include <algorithm>
#include <cmath>

namespace warp_hypad {
namespace geometry {
namespace {

double max_offdiagonal(const double a[4][4]) {
    double m = 0.0;
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j) m = std::fmax(m, std::fabs(a[i][j]));
    return m;
}

}  // namespace

SignatureResult analyse_signature(const Mat4<double>& g,
                                  double relative_zero_tolerance) {
    SignatureResult out;

    // Work on the symmetric part. A metric is symmetric by definition; if the
    // caller hands us something that is not, symmetrising is the only sensible
    // reading and the asymmetry is caught by a separate check.
    double a[4][4];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) a[i][j] = 0.5 * (g(i, j) + g(j, i));

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (!std::isfinite(a[i][j])) {
                out.converged = false;
                for (int k = 0; k < 4; ++k) out.eigenvalues[k] = std::nan("");
                return out;
            }

    // Cyclic Jacobi. Each sweep annihilates every off-diagonal entry once;
    // convergence is quadratic, so a handful of sweeps is ample at 4x4.
    constexpr int kMaxSweeps = 60;
    int sweep = 0;
    for (; sweep < kMaxSweeps; ++sweep) {
        const double off = max_offdiagonal(a);
        if (off <= 1e-300) break;

        double scale = 0.0;
        for (int i = 0; i < 4; ++i) scale = std::fmax(scale, std::fabs(a[i][i]));
        if (off <= 1e-18 * std::fmax(1.0, scale)) break;

        for (int p = 0; p < 3; ++p) {
            for (int q = p + 1; q < 4; ++q) {
                if (std::fabs(a[p][q]) <= 1e-300) continue;

                // Standard Jacobi rotation angle, computed in the numerically
                // stable form that avoids cancellation when a[p][p] ~ a[q][q].
                const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                const double t =
                    (theta >= 0.0 ? 1.0 : -1.0) /
                    (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0);
                const double s = t * c;

                for (int k = 0; k < 4; ++k) {
                    const double akp = a[k][p];
                    const double akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 4; ++k) {
                    const double apk = a[p][k];
                    const double aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
        }
    }

    out.converged = sweep < kMaxSweeps;
    out.residual_offdiagonal = max_offdiagonal(a);

    for (int i = 0; i < 4; ++i) out.eigenvalues[i] = a[i][i];
    std::sort(out.eigenvalues, out.eigenvalues + 4);

    double biggest = 0.0;
    for (int i = 0; i < 4; ++i)
        biggest = std::fmax(biggest, std::fabs(out.eigenvalues[i]));
    const double zero_cut = relative_zero_tolerance * std::fmax(biggest, 1e-300);

    for (int i = 0; i < 4; ++i) {
        const double v = out.eigenvalues[i];
        if (std::fabs(v) <= zero_cut) {
            ++out.zero;
        } else if (v < 0.0) {
            ++out.negative;
        } else {
            ++out.positive;
        }
    }

    out.lorentzian = out.converged && out.negative == 1 && out.positive == 3 &&
                     out.zero == 0;
    return out;
}

}  // namespace geometry
}  // namespace warp_hypad
