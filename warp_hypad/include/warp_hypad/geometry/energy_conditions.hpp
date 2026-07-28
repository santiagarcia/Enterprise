#pragma once

#include <cmath>
#include <string>

#include "warp_hypad/geometry/observer.hpp"

// ---------------------------------------------------------------------------
// Milestone A6, part 2: the pointwise classical energy conditions.
//
// --- Why the NEC needs a real minimiser --------------------------------------
//
// For a null vector k^a = (1, n) with |n| = 1,
//
//     N(n) = T_ab k^a k^b = rho + 2 T_{0i} n_i + S_ij n_i n_j
//          = rho - 2 q_i n_i + S_ij n_i n_j
//
// and the NEC asks for min over the whole unit sphere. Testing +-x, +-y, +-z
// samples six points of a continuum and will miss the minimum whenever the
// worst direction is off-axis -- which, for an Alcubierre bubble, it always is:
// the momentum density q is large and points along the wall, so the linear term
// tilts the minimiser away from every coordinate axis.
//
// This is a quadratic on a sphere, so it has an exact solution rather than
// needing a general optimiser. Minimising S_ij n_i n_j - 2 q_i n_i subject to
// |n| = 1 gives the stationarity condition
//
//     (S - mu I) n = q
//
// with mu a Lagrange multiplier, and the minimiser corresponds to the SMALLEST
// root mu of the secular equation
//
//     phi(mu) = | (S - mu I)^{-1} q |^2 - 1 = 0,      mu < lambda_min(S)
//
// which is the classical trust-region subproblem. phi is monotone decreasing on
// (-inf, lambda_min), so a bracketed bisection converges unconditionally. The
// "hard case", q orthogonal to the lowest eigenvector, is handled by falling
// back to the eigenvector direction.
//
// The result is then cross-checked three independent ways -- deterministic
// spherical quadrature, multi-start local descent, and direct substitution of
// the reported direction -- because a wrong minimiser would silently understate
// every NEC violation the project ever reports.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace geometry {

// Every condition result is one of these. "pass" is not available: a condition
// that cannot be certified within the numerical error envelope is reported as
// indeterminate, never as satisfied.
enum class ConditionStatus {
    CertifiedSatisfied,
    CertifiedViolated,
    IndeterminateWithinDiscretization,
    ObserverOrTetradInvalid
};

inline const char* to_string(ConditionStatus s) {
    switch (s) {
        case ConditionStatus::CertifiedSatisfied: return "certified_satisfied";
        case ConditionStatus::CertifiedViolated: return "certified_violated";
        case ConditionStatus::IndeterminateWithinDiscretization:
            return "indeterminate_within_discretization";
        case ConditionStatus::ObserverOrTetradInvalid:
            return "observer_or_tetrad_invalid";
    }
    return "unknown";
}

// A measured value is certified only when its magnitude exceeds the numerical
// error envelope by this factor. Chosen rather than derived: the envelope itself
// is a Richardson estimate with its own uncertainty, and a factor of 3 keeps a
// borderline value out of the certified categories instead of forcing it into
// one.
inline constexpr double kCertificationSafetyFactor = 3.0;

inline ConditionStatus certify(double value, double error_envelope,
                               bool satisfied_when_nonnegative = true) {
    const double margin = kCertificationSafetyFactor * error_envelope;
    if (std::fabs(value) <= margin)
        return ConditionStatus::IndeterminateWithinDiscretization;
    const bool ok = satisfied_when_nonnegative ? (value > 0.0) : (value < 0.0);
    return ok ? ConditionStatus::CertifiedSatisfied
              : ConditionStatus::CertifiedViolated;
}

// --- The NEC functional and its exact minimiser -------------------------------

struct NecResult {
    double n_min = 0.0;              // the global minimum of N over the sphere
    double direction[3] = {0, 0, 0}; // the minimising n
    double stationarity_residual = 0.0;   // |(S - mu I) n - q|
    double sphere_constraint_residual = 0.0;  // ||n| - 1|
    double quadrature_crosscheck = 0.0;   // best quadrature value minus n_min
    double multistart_crosscheck = 0.0;   // best local descent minus n_min
    double direct_substitution = 0.0;     // N(direction) recomputed from scratch
    bool hard_case = false;
    bool ok = false;
};

// N(n) = rho - 2 q.n + n^T S n.
inline double nec_value(const FrameStress& f, const double n[3]) {
    double v = f.rho;
    for (int i = 0; i < 3; ++i) v -= 2.0 * f.q[i] * n[i];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) v += f.S[3 * i + j] * n[i] * n[j];
    return v;
}

namespace detail {

// Solve (S - mu I) x = q for symmetric 3x3, by cofactor inverse. Returns false
// when the matrix is numerically singular.
inline bool solve_shifted(const double S[9], double mu, const double q[3],
                          double x[3]) {
    double m[9];
    for (int i = 0; i < 9; ++i) m[i] = S[i];
    m[0] -= mu; m[4] -= mu; m[8] -= mu;

    const double c00 = m[4] * m[8] - m[5] * m[7];
    const double c01 = m[5] * m[6] - m[3] * m[8];
    const double c02 = m[3] * m[7] - m[4] * m[6];
    const double dt = m[0] * c00 + m[1] * c01 + m[2] * c02;
    if (std::fabs(dt) < 1e-300) return false;

    const double c11 = m[0] * m[8] - m[2] * m[6];
    const double c12 = m[1] * m[6] - m[0] * m[7];
    const double c22 = m[0] * m[4] - m[1] * m[3];
    const double inv = 1.0 / dt;
    // Symmetric, so the inverse is symmetric too.
    const double a00 = c00 * inv, a01 = c01 * inv, a02 = c02 * inv;
    const double a11 = c11 * inv, a12 = c12 * inv, a22 = c22 * inv;

    x[0] = a00 * q[0] + a01 * q[1] + a02 * q[2];
    x[1] = a01 * q[0] + a11 * q[1] + a12 * q[2];
    x[2] = a02 * q[0] + a12 * q[1] + a22 * q[2];
    return true;
}

// Unit eigenvector of a symmetric 3x3 for a known eigenvalue, by taking the
// largest-norm column of the adjugate of (S - lambda I).
inline void eigenvector_for(const double S[9], double lam, double v[3]) {
    double m[9];
    for (int i = 0; i < 9; ++i) m[i] = S[i];
    m[0] -= lam; m[4] -= lam; m[8] -= lam;

    // Rows of the adjugate: cross products of pairs of rows of m.
    const double cand[3][3] = {
        {m[1] * m[5] - m[2] * m[4], m[2] * m[3] - m[0] * m[5],
         m[0] * m[4] - m[1] * m[3]},
        {m[4] * m[8] - m[5] * m[7], m[5] * m[6] - m[3] * m[8],
         m[3] * m[7] - m[4] * m[6]},
        {m[7] * m[2] - m[8] * m[1], m[8] * m[0] - m[6] * m[2],
         m[6] * m[1] - m[7] * m[0]}};

    int best = 0;
    double best_n = -1.0;
    for (int k = 0; k < 3; ++k) {
        const double n2 = cand[k][0] * cand[k][0] + cand[k][1] * cand[k][1] +
                          cand[k][2] * cand[k][2];
        if (n2 > best_n) { best_n = n2; best = k; }
    }
    if (best_n <= 0.0) { v[0] = 1.0; v[1] = v[2] = 0.0; return; }
    const double inv = 1.0 / std::sqrt(best_n);
    for (int i = 0; i < 3; ++i) v[i] = cand[best][i] * inv;
}

// An orthonormal eigenbasis of a symmetric 3x3, rows of V. Built from
// eigenvector_for and then Gram-Schmidted in ascending-eigenvalue order, which
// repairs the near-degenerate case where two eigenvectors come back almost
// parallel. Deterministic: fixed order, positive-projection sign convention.
inline void orthonormal_eigenbasis(const double S[9], const double lam[3],
                                   double V[3][3]) {
    for (int i = 0; i < 3; ++i) {
        eigenvector_for(S, lam[i], V[i]);
        for (int j = 0; j < i; ++j) {
            const double d =
                V[i][0] * V[j][0] + V[i][1] * V[j][1] + V[i][2] * V[j][2];
            for (int k = 0; k < 3; ++k) V[i][k] -= d * V[j][k];
        }
        double n2 = V[i][0] * V[i][0] + V[i][1] * V[i][1] + V[i][2] * V[i][2];
        if (n2 < 1e-20) {
            // Degenerate: substitute whichever coordinate axis is most nearly
            // orthogonal to what is already fixed, then re-orthogonalise.
            for (int ax = 0; ax < 3; ++ax) {
                double c[3] = {0, 0, 0};
                c[ax] = 1.0;
                for (int j = 0; j < i; ++j) {
                    const double d =
                        c[0] * V[j][0] + c[1] * V[j][1] + c[2] * V[j][2];
                    for (int k = 0; k < 3; ++k) c[k] -= d * V[j][k];
                }
                const double cn = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
                if (cn > 1e-6) {
                    for (int k = 0; k < 3; ++k) V[i][k] = c[k];
                    n2 = cn;
                    break;
                }
            }
        }
        const double inv = 1.0 / std::sqrt(n2);
        for (int k = 0; k < 3; ++k) V[i][k] *= inv;
    }
}

// Deterministic Fibonacci-lattice directions on the sphere. Fixed count, fixed
// formula, no randomness, so two runs agree bit for bit.
inline void fibonacci_direction(int k, int n, double out[3]) {
    const double kPi = 3.14159265358979323846;
    const double golden = kPi * (3.0 - std::sqrt(5.0));
    const double z = 1.0 - 2.0 * (static_cast<double>(k) + 0.5) /
                               static_cast<double>(n);
    const double r = std::sqrt(std::fmax(0.0, 1.0 - z * z));
    const double th = golden * static_cast<double>(k);
    out[0] = r * std::cos(th);
    out[1] = r * std::sin(th);
    out[2] = z;
}

}  // namespace detail

// The exact global minimum of N over the unit sphere, plus three independent
// cross-checks.
inline NecResult nec_minimum(const FrameStress& f, int quadrature_points = 4096) {
    NecResult r;

    // Everything is done in the EIGENBASIS of S. In that basis the problem
    // separates completely -- minimise  sum_i lambda_i n_i^2 - 2 a_i n_i  with
    // a_i = q.v_i -- and the hard case becomes a visible condition (a_1 = 0)
    // rather than a numerical accident to be discovered when a bracket fails.
    double lam[3];
    symmetric3_eigenvalues(f.S, lam);

    double V[3][3];
    detail::orthonormal_eigenbasis(f.S, lam, V);

    double a[3];
    for (int i = 0; i < 3; ++i)
        a[i] = f.q[0] * V[i][0] + f.q[1] * V[i][1] + f.q[2] * V[i][2];

    double m[3] = {0, 0, 0};  // the minimiser, in eigen-coordinates
    const double scale = 1.0 + std::fabs(lam[0]) + std::fabs(lam[2]);

    // HARD CASE: q has no component along the lowest eigenvector. Then
    // (S - mu I)^{-1} q stays bounded as mu -> lambda_1, the secular equation
    // may have no root below lambda_1, and the minimiser is the bounded
    // particular solution plus whatever multiple of v_1 the unit-norm
    // constraint needs. A solver that assumed a root exists silently returns a
    // direction that is not the minimum -- which is exactly what a spherical
    // quadrature then beats.
    bool hard = std::fabs(a[0]) <= 1e-10 * scale;
    if (hard) {
        double partial2 = 0.0;
        for (int i = 1; i < 3; ++i) {
            const double gap = lam[i] - lam[0];
            m[i] = (gap > 1e-12 * scale) ? a[i] / gap : 0.0;
            partial2 += m[i] * m[i];
        }
        if (partial2 <= 1.0) {
            m[0] = std::sqrt(1.0 - partial2);   // sign chosen below
        } else {
            hard = false;  // the particular solution is already too long
        }
    }

    if (!hard) {
        // Secular equation  sum_i (a_i/(lambda_i - mu))^2 = 1  for mu < lambda_1.
        // The left side increases monotonically toward lambda_1 and tends to 0
        // as mu -> -infinity, so bisection on that interval is unconditional.
        auto norm_at = [&](double mu) {
            double s = 0.0;
            for (int i = 0; i < 3; ++i) {
                const double d = lam[i] - mu;
                if (std::fabs(d) < 1e-300) return 1e300;
                s += (a[i] / d) * (a[i] / d);
            }
            return s;
        };

        double hi = lam[0] - 1e-14 * scale;
        double lo = hi - scale;
        for (int it = 0; it < 300 && norm_at(lo) > 1.0; ++it)
            lo = lam[0] - 2.0 * (lam[0] - lo);

        for (int it = 0; it < 300; ++it) {
            const double mid = 0.5 * (lo + hi);
            if (norm_at(mid) > 1.0) hi = mid; else lo = mid;
        }
        const double mu = 0.5 * (lo + hi);
        double n2 = 0.0;
        for (int i = 0; i < 3; ++i) {
            const double d = lam[i] - mu;
            m[i] = (std::fabs(d) > 1e-300) ? a[i] / d : 0.0;
            n2 += m[i] * m[i];
        }
        const double nn = std::sqrt(n2);
        if (nn > 1e-14)
            for (int i = 0; i < 3; ++i) m[i] /= nn;
        else
            m[0] = 1.0;
    }
    r.hard_case = hard;

    // Back to Cartesian, then fix the free sign (hard case only) by value.
    auto to_cart = [&](const double mm[3], double out[3]) {
        for (int k = 0; k < 3; ++k)
            out[k] = mm[0] * V[0][k] + mm[1] * V[1][k] + mm[2] * V[2][k];
    };
    to_cart(m, r.direction);
    if (hard) {
        double m2[3] = {-m[0], m[1], m[2]};
        double alt[3];
        to_cart(m2, alt);
        if (nec_value(f, alt) < nec_value(f, r.direction))
            for (int i = 0; i < 3; ++i) r.direction[i] = alt[i];
    }

    // Stationarity residual, in the tangent space: the gradient of N restricted
    // to the sphere must vanish. Stated this way rather than as
    // |(S - mu I)n - q| so the same measure covers the hard case, where mu is
    // lambda_1 and the residual lies entirely along v_1.
    {
        double gr[3];
        for (int i = 0; i < 3; ++i) {
            double s = -2.0 * f.q[i];
            for (int j = 0; j < 3; ++j) s += 2.0 * f.S[3 * i + j] * r.direction[j];
            gr[i] = s;
        }
        const double dot = gr[0] * r.direction[0] + gr[1] * r.direction[1] +
                           gr[2] * r.direction[2];
        double tang = 0.0;
        for (int i = 0; i < 3; ++i) {
            const double t = gr[i] - dot * r.direction[i];
            tang += t * t;
        }
        r.stationarity_residual = std::sqrt(tang);
    }

    const double dn = std::sqrt(r.direction[0] * r.direction[0] +
                                r.direction[1] * r.direction[1] +
                                r.direction[2] * r.direction[2]);
    r.sphere_constraint_residual = std::fabs(dn - 1.0);
    r.n_min = nec_value(f, r.direction);
    r.direct_substitution = r.n_min;

    // Cross-check 1: deterministic spherical quadrature. Must not find anything
    // lower than the analytic minimum by more than its own resolution.
    double best_q = 1e300;
    for (int k = 0; k < quadrature_points; ++k) {
        double n[3];
        detail::fibonacci_direction(k, quadrature_points, n);
        const double v = nec_value(f, n);
        if (v < best_q) best_q = v;
    }
    r.quadrature_crosscheck = best_q - r.n_min;

    // Cross-check 2: multi-start projected descent from fixed seeds.
    double best_m = 1e300;
    for (int s = 0; s < 26; ++s) {
        double n[3];
        detail::fibonacci_direction(s, 26, n);
        double step = 0.5;
        double cur = nec_value(f, n);
        for (int it = 0; it < 400 && step > 1e-14; ++it) {
            // Gradient of N on the sphere, projected tangentially.
            double gr[3];
            for (int i = 0; i < 3; ++i) {
                double s2 = -2.0 * f.q[i];
                for (int j = 0; j < 3; ++j) s2 += 2.0 * f.S[3 * i + j] * n[j];
                gr[i] = s2;
            }
            const double dot = gr[0] * n[0] + gr[1] * n[1] + gr[2] * n[2];
            double t[3];
            for (int i = 0; i < 3; ++i) t[i] = gr[i] - dot * n[i];
            double cand[3];
            for (int i = 0; i < 3; ++i) cand[i] = n[i] - step * t[i];
            const double cn = std::sqrt(cand[0] * cand[0] + cand[1] * cand[1] +
                                        cand[2] * cand[2]);
            if (cn < 1e-14) { step *= 0.5; continue; }
            for (int i = 0; i < 3; ++i) cand[i] /= cn;
            const double vc = nec_value(f, cand);
            if (vc < cur) {
                cur = vc;
                for (int i = 0; i < 3; ++i) n[i] = cand[i];
            } else {
                step *= 0.5;
            }
        }
        if (cur < best_m) best_m = cur;
    }
    r.multistart_crosscheck = best_m - r.n_min;

    r.ok = true;
    return r;
}

// --- WEC, SEC, DEC ------------------------------------------------------------
//
// Never inferred from a single observer. For each condition the worst value over
// a deterministic sweep of the timelike velocity ball is reported, together with
// the boost that achieved it.

struct ConditionResult {
    double worst_value = 0.0;
    double worst_velocity[3] = {0, 0, 0};
    ConditionStatus status = ConditionStatus::IndeterminateWithinDiscretization;
    double error_envelope = 0.0;
    double certification_margin = 0.0;
};

namespace detail {

// A deterministic sweep of the open velocity ball: Fibonacci directions times a
// fixed ladder of speeds. Never touches |v| = 1.
template <typename Fn>
inline double worst_over_velocity_ball(Fn&& value_at, double out_v[3]) {
    static const double speeds[6] = {0.0, 0.25, 0.5, 0.75, 0.9, 0.98};
    constexpr int kDirs = 128;
    double worst = 1e300;
    for (int s = 0; s < 6; ++s) {
        const int dirs = (s == 0) ? 1 : kDirs;
        for (int k = 0; k < dirs; ++k) {
            double n[3] = {0, 0, 0};
            if (s > 0) fibonacci_direction(k, kDirs, n);
            const double v[3] = {speeds[s] * n[0], speeds[s] * n[1],
                                 speeds[s] * n[2]};
            const double val = value_at(v);
            if (val < worst) {
                worst = val;
                for (int i = 0; i < 3; ++i) out_v[i] = v[i];
            }
        }
    }
    return worst;
}

// T_ab u^a u^b for a frame velocity v, with u^a = gamma(1, v).
inline double frame_contract_uu(const FrameStress& f, const double v[3]) {
    const double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    const double g2 = 1.0 / (1.0 - v2);
    double s = f.T[0][0];
    for (int i = 0; i < 3; ++i) s += 2.0 * f.T[0][i + 1] * v[i];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) s += f.T[i + 1][j + 1] * v[i] * v[j];
    return g2 * s;
}

}  // namespace detail

// WEC: T_ab u^a u^b >= 0 for every timelike u.
inline ConditionResult wec(const FrameStress& f, double envelope) {
    ConditionResult r;
    r.error_envelope = envelope;
    r.worst_value = detail::worst_over_velocity_ball(
        [&](const double v[3]) { return detail::frame_contract_uu(f, v); },
        r.worst_velocity);
    r.status = certify(r.worst_value, envelope);
    r.certification_margin =
        std::fabs(r.worst_value) - kCertificationSafetyFactor * envelope;
    return r;
}

// SEC: (T_ab - (1/2) T eta_ab) u^a u^b >= 0, with T = eta^ab T_ab.
inline ConditionResult sec(const FrameStress& f, double envelope) {
    ConditionResult r;
    r.error_envelope = envelope;
    const double trace = -f.T[0][0] + f.S[0] + f.S[4] + f.S[8];

    FrameStress g = f;
    for (int A = 0; A < 4; ++A)
        for (int B = 0; B < 4; ++B) {
            const double eta = (A == B) ? ((A == 0) ? -1.0 : 1.0) : 0.0;
            g.T[A][B] = f.T[A][B] - 0.5 * trace * eta;
        }
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) g.S[3 * i + j] = g.T[i + 1][j + 1];

    r.worst_value = detail::worst_over_velocity_ball(
        [&](const double v[3]) { return detail::frame_contract_uu(g, v); },
        r.worst_velocity);
    r.status = certify(r.worst_value, envelope);
    r.certification_margin =
        std::fabs(r.worst_value) - kCertificationSafetyFactor * envelope;
    return r;
}

// DEC: for every future timelike u, the flux j^a = -T^a_b u^b is future causal,
// i.e. j^0 >= 0 and j.j <= 0. Both are folded into a single worst value: the
// smaller of j^0 and -j.j, so a negative worst means a violation of one or the
// other.
inline ConditionResult dec(const FrameStress& f, double envelope) {
    ConditionResult r;
    r.error_envelope = envelope;
    r.worst_value = detail::worst_over_velocity_ball(
        [&](const double v[3]) {
            const double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            const double gam = 1.0 / std::sqrt(1.0 - v2);
            const double u[4] = {gam, gam * v[0], gam * v[1], gam * v[2]};
            // j^a = -T^a{}_b u^b = -eta^{ac} T_{cb} u^b.
            //
            // u^b is already contravariant, so T_{cb} u^b needs NO metric
            // factor; the single eta raises the free index c only. Applying eta
            // to both slots -- as an earlier version did -- flips the sign of
            // the spatial flux relative to the energy flux and reports DEC as
            // violated for dust, which is how the error was caught.
            double j[4];
            for (int A = 0; A < 4; ++A) {
                double s = 0.0;
                for (int B = 0; B < 4; ++B) s += f.T[A][B] * u[B];
                const double etaA = (A == 0) ? -1.0 : 1.0;
                j[A] = -etaA * s;
            }
            const double jj = -j[0] * j[0] + j[1] * j[1] + j[2] * j[2] + j[3] * j[3];
            return std::fmin(j[0], -jj);
        },
        r.worst_velocity);
    r.status = certify(r.worst_value, envelope);
    r.certification_margin =
        std::fabs(r.worst_value) - kCertificationSafetyFactor * envelope;
    return r;
}

}  // namespace geometry
}  // namespace warp_hypad
