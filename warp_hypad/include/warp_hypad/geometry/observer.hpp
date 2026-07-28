#pragma once

#include <cmath>
#include <string>

#include "warp_hypad/core/mat4.hpp"
#include "warp_hypad/core/scalar_traits.hpp"

// ---------------------------------------------------------------------------
// Milestone A6, part 1: ADM decomposition, orthonormal tetrads, and observers.
//
// A5 produced a stress-energy tensor with coordinate indices. Nothing about it
// is an energy density. This file builds the machinery that makes an energy
// density definable: a foliation, an orthonormal frame, and explicitly typed
// observer four-velocities that are checked for admissibility before use.
//
// --- The failure this guards against -----------------------------------------
//
// The tempting move is u^mu = (1,0,0,0)/sqrt(-g_00), the coordinate-stationary
// observer, and then to call T_{mu nu} u^mu u^nu "the" energy density. For an
// Alcubierre bubble at v_s >= 1 that observer DOES NOT EXIST: g_00 reaches zero
// and then turns positive, so the normalisation is a division by zero and then
// a square root of a negative number. Taking |g_00| there produces a number --
// a plausible, plottable, entirely fictitious number.
//
// So: the coordinate-stationary observer is a separate type from the Eulerian
// observer, it is admissible only where g_00 < 0, and it is REJECTED with a
// reason elsewhere rather than regularised. The Eulerian observer, built from
// the ADM normal, remains timelike straight through that crossing, which is
// exactly why the foliation is the right structure to hang the physics on.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace geometry {

// --- ADM 3+1 decomposition ---------------------------------------------------
//
//     ds^2 = -alpha^2 dt^2 + gamma_ij (dx^i + beta^i dt)(dx^j + beta^j dt)
//
// recovered from the INVERSE metric, which is where lapse and shift appear
// algebraically:
//
//     g^{00} = -1/alpha^2        =>  alpha = 1/sqrt(-g^{00})
//     beta^i = -g^{0i}/g^{00}
//     gamma_ij = g_ij            (the spatial block of the covariant metric)
struct AdmDecomposition {
    double alpha = 0.0;      // lapse, must be > 0
    double beta[3] = {0, 0, 0};   // shift, contravariant spatial
    double gamma[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};  // 3x3 row-major, covariant

    bool valid = false;
    std::string rejection_reason;

    // Diagnostics.
    double lapse_squared_from_inverse = 0.0;
    double gamma_min_eigenvalue = 0.0;   // must be > 0 for a Riemannian slice
    double normal_norm_residual = 0.0;   // |g(n,n) + 1|
};

// Smallest eigenvalue of a symmetric 3x3, by the closed-form trigonometric
// solution. Used only as a positive-definiteness diagnostic, so operating on
// real parts is correct.
inline double symmetric3_min_eigenvalue(const double m[9]) {
    const double p1 = m[1] * m[1] + m[2] * m[2] + m[5] * m[5];
    const double tr = m[0] + m[4] + m[8];
    if (p1 <= 0.0) {
        double lo = m[0];
        if (m[4] < lo) lo = m[4];
        if (m[8] < lo) lo = m[8];
        return lo;
    }
    const double q = tr / 3.0;
    const double d0 = m[0] - q, d4 = m[4] - q, d8 = m[8] - q;
    const double p2 = d0 * d0 + d4 * d4 + d8 * d8 + 2.0 * p1;
    const double p = std::sqrt(p2 / 6.0);
    double b[9];
    for (int i = 0; i < 9; ++i) b[i] = m[i] / p;
    b[0] -= q / p; b[4] -= q / p; b[8] -= q / p;
    const double det =
        b[0] * (b[4] * b[8] - b[5] * b[7]) - b[1] * (b[3] * b[8] - b[5] * b[6]) +
        b[2] * (b[3] * b[7] - b[4] * b[6]);
    double r = det / 2.0;
    if (r < -1.0) r = -1.0;
    if (r > 1.0) r = 1.0;
    const double phi = std::acos(r) / 3.0;
    const double kPi = 3.14159265358979323846;
    // eig1 >= eig2 >= eig3; only the smallest is wanted.
    return q + 2.0 * p * std::cos(phi + (2.0 * kPi / 3.0));
}

// Decompose from the covariant metric and its inverse.
//
// Both are taken as arguments rather than inverting here, because the caller
// already has a verified inverse from the curvature pipeline and a second
// inversion could silently disagree with it.
inline AdmDecomposition adm_decompose(const Mat4<double>& g,
                                      const Mat4<double>& g_inv) {
    AdmDecomposition a;

    const double g00_up = g_inv(0, 0);
    a.lapse_squared_from_inverse = -1.0 / g00_up;

    if (!(g00_up < 0.0)) {
        a.valid = false;
        a.rejection_reason =
            "g^{00} = " + std::to_string(g00_up) +
            " is not negative, so the t = const surfaces are not spacelike and "
            "no ADM foliation exists here";
        return a;
    }

    a.alpha = 1.0 / std::sqrt(-g00_up);
    for (int i = 0; i < 3; ++i) a.beta[i] = -g_inv(0, i + 1) / g00_up;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) a.gamma[3 * i + j] = g(i + 1, j + 1);

    a.gamma_min_eigenvalue = symmetric3_min_eigenvalue(a.gamma);
    if (!(a.gamma_min_eigenvalue > 0.0)) {
        a.valid = false;
        a.rejection_reason =
            "the spatial metric is not positive definite (min eigenvalue " +
            std::to_string(a.gamma_min_eigenvalue) + ")";
        return a;
    }
    if (!(a.alpha > 0.0)) {
        a.valid = false;
        a.rejection_reason = "lapse is not positive";
        return a;
    }

    a.valid = true;
    return a;
}

// The future-directed unit normal to the t = const slices:
//
//     n^mu = ( 1/alpha, -beta^i/alpha ),   n_mu = ( -alpha, 0, 0, 0 )
//
// NOT proportional to d/dt in general. Where g_00 = 0 the coordinate vector
// d/dt is null and the coordinate-stationary observer ceases to exist, while
// this one remains a perfectly good timelike observer -- which is the whole
// reason the foliation is the right structure to build on.
inline void eulerian_normal(const AdmDecomposition& a, double n_up[4],
                            double n_down[4]) {
    n_up[0] = 1.0 / a.alpha;
    for (int i = 0; i < 3; ++i) n_up[i + 1] = -a.beta[i] / a.alpha;
    n_down[0] = -a.alpha;
    n_down[1] = n_down[2] = n_down[3] = 0.0;
}

inline double four_velocity_norm(const Mat4<double>& g, const double u[4]) {
    double s = 0.0;
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu) s += g(mu, nu) * u[mu] * u[nu];
    return s;
}

// --- Orthonormal tetrad ------------------------------------------------------
//
//     g_{mu nu} e_a^mu e_b^nu = eta_{ab}
//
// e_0 is the Eulerian normal. The spatial triad is built by Gram-Schmidt on the
// coordinate basis vectors PROJECTED orthogonal to n, in a fixed order, with a
// positive-diagonal convention on the resulting factor.
//
// Determinism matters more than elegance here. An eigen-decomposition would
// return eigenvectors whose signs and ordering flip between neighbouring grid
// points wherever eigenvalues nearly cross, and each flip becomes a spurious
// sign change in q_i and in the worst NEC direction -- a discontinuity in an
// animation that no physics produced. Fixed-order Gram-Schmidt with a positive
// diagonal is the same construction as a Cholesky factor and is continuous
// wherever the metric is.
struct Tetrad {
    // e[a][mu]: the a-th frame vector, contravariant components.
    double e[4][4] = {{0}};
    // The coframe: e_dn[a][mu] = eta_{ab} g_{mu nu} e[b][nu], so that
    // T_{ab} = e_dn[a][mu] e_dn[b][nu] T^{mu nu} and, for a covariant T,
    // T_{ab} = e[a][mu] e[b][nu] T_{mu nu}.
    double e_dn[4][4] = {{0}};

    bool valid = false;
    std::string rejection_reason;

    double orthonormality_residual = 0.0;  // max |g(e_a,e_b) - eta_ab|
    double handedness = 0.0;               // det of the frame matrix; > 0 wanted
};

inline Tetrad build_tetrad(const Mat4<double>& g, const AdmDecomposition& a) {
    Tetrad t;
    if (!a.valid) {
        t.valid = false;
        t.rejection_reason = "no ADM decomposition: " + a.rejection_reason;
        return t;
    }

    double n_up[4], n_dn[4];
    eulerian_normal(a, n_up, n_dn);
    for (int mu = 0; mu < 4; ++mu) t.e[0][mu] = n_up[mu];

    // Spatial triad: project d/dx, d/dy, d/dz orthogonal to everything already
    // fixed, then normalise. Fixed axis order, never reordered by magnitude.
    for (int i = 0; i < 3; ++i) {
        double v[4] = {0, 0, 0, 0};
        v[i + 1] = 1.0;

        // Remove the timelike component. eta_00 = -1, hence the plus sign:
        // v -> v + <v, e_0> e_0.
        double dot0 = 0.0;
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) dot0 += g(mu, nu) * v[mu] * t.e[0][nu];
        for (int mu = 0; mu < 4; ++mu) v[mu] += dot0 * t.e[0][mu];

        // Remove the previously built spatial directions.
        for (int j = 1; j <= i; ++j) {
            double dj = 0.0;
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu)
                    dj += g(mu, nu) * v[mu] * t.e[j][nu];
            for (int mu = 0; mu < 4; ++mu) v[mu] -= dj * t.e[j][mu];
        }

        const double n2 = four_velocity_norm(g, v);
        if (!(n2 > 0.0)) {
            t.valid = false;
            t.rejection_reason =
                "spatial basis vector " + std::to_string(i) +
                " degenerated during orthogonalisation (norm^2 = " +
                std::to_string(n2) + ")";
            return t;
        }
        const double inv = 1.0 / std::sqrt(n2);  // positive root: fixes the sign
        for (int mu = 0; mu < 4; ++mu) t.e[i + 1][mu] = v[mu] * inv;
    }

    // Orthonormality residual against eta.
    for (int A = 0; A < 4; ++A)
        for (int B = 0; B < 4; ++B) {
            double s = 0.0;
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu)
                    s += g(mu, nu) * t.e[A][mu] * t.e[B][nu];
            const double eta = (A == B) ? ((A == 0) ? -1.0 : 1.0) : 0.0;
            const double d = std::fabs(s - eta);
            if (d > t.orthonormality_residual) t.orthonormality_residual = d;
        }

    // Coframe, for lowering frame indices.
    for (int A = 0; A < 4; ++A) {
        const double eta = (A == 0) ? -1.0 : 1.0;
        for (int mu = 0; mu < 4; ++mu) {
            double s = 0.0;
            for (int nu = 0; nu < 4; ++nu) s += g(mu, nu) * t.e[A][nu];
            t.e_dn[A][mu] = eta * s;
        }
    }

    // Handedness: determinant of the 4x4 frame matrix. A sign flip between
    // neighbouring points would show up here first.
    Mat4<double> F;
    for (int A = 0; A < 4; ++A)
        for (int mu = 0; mu < 4; ++mu) F(A, mu) = t.e[A][mu];
    t.handedness = det(F);

    t.valid = true;
    return t;
}

// --- Observer families -------------------------------------------------------
//
// Separately typed rather than a single struct with a string tag, so a caller
// cannot pass a rejected observer where a valid one is required without saying
// so.

enum class ObserverKind {
    Eulerian,
    CoordinateStationary,
    BubbleComoving,
    TetradBoosted
};

inline const char* to_string(ObserverKind k) {
    switch (k) {
        case ObserverKind::Eulerian: return "eulerian";
        case ObserverKind::CoordinateStationary: return "coordinate_stationary";
        case ObserverKind::BubbleComoving: return "bubble_comoving";
        case ObserverKind::TetradBoosted: return "tetrad_boosted";
    }
    return "unknown";
}

struct Observer {
    ObserverKind kind = ObserverKind::Eulerian;
    double u[4] = {0, 0, 0, 0};      // contravariant four-velocity
    bool valid = false;
    std::string rejection_reason;
    double norm = 0.0;               // g(u,u), should be -1
    double normalization_residual = 0.0;  // |g(u,u) + 1|
};

inline Observer make_eulerian(const Mat4<double>& g, const AdmDecomposition& a) {
    Observer o;
    o.kind = ObserverKind::Eulerian;
    if (!a.valid) {
        o.valid = false;
        o.rejection_reason = "no ADM foliation: " + a.rejection_reason;
        return o;
    }
    double n_dn[4];
    eulerian_normal(a, o.u, n_dn);
    o.norm = four_velocity_norm(g, o.u);
    o.normalization_residual = std::fabs(o.norm + 1.0);
    o.valid = o.normalization_residual < 1e-10;
    if (!o.valid)
        o.rejection_reason = "Eulerian normal failed normalisation: g(n,n) = " +
                             std::to_string(o.norm);
    return o;
}

// u^mu = (1,0,0,0)/sqrt(-g_00), admissible ONLY where g_00 < 0.
//
// At g_00 = 0 the coordinate vector d/dt is null; for g_00 > 0 it is spacelike.
// Neither is an observer. The rejection is explicit and carries the value, and
// under no circumstance is |g_00| substituted to keep the arithmetic alive.
inline Observer make_coordinate_stationary(const Mat4<double>& g) {
    Observer o;
    o.kind = ObserverKind::CoordinateStationary;
    const double g00 = g(0, 0);

    if (g00 == 0.0) {
        o.valid = false;
        o.rejection_reason =
            "g_00 = 0 exactly: d/dt is NULL here, not an observer worldline";
        return o;
    }
    if (g00 > 0.0) {
        o.valid = false;
        o.rejection_reason =
            "g_00 = " + std::to_string(g00) +
            " > 0: d/dt is SPACELIKE here, not an observer worldline";
        return o;
    }

    const double inv = 1.0 / std::sqrt(-g00);
    o.u[0] = inv;
    o.u[1] = o.u[2] = o.u[3] = 0.0;
    o.norm = four_velocity_norm(g, o.u);
    o.normalization_residual = std::fabs(o.norm + 1.0);
    o.valid = o.normalization_residual < 1e-10;
    if (!o.valid)
        o.rejection_reason = "normalisation failed: g(u,u) = " +
                             std::to_string(o.norm);
    return o;
}

// An observer carried along with the bubble: coordinate velocity (v_s, 0, 0).
// Timelikeness is CHECKED, not assumed -- the bubble centre's coordinate
// trajectory is not automatically a physical worldline, and for a superluminal
// bubble in the exterior region it is not one.
inline Observer make_bubble_comoving(const Mat4<double>& g, double v_s) {
    Observer o;
    o.kind = ObserverKind::BubbleComoving;
    double w[4] = {1.0, v_s, 0.0, 0.0};
    const double n2 = four_velocity_norm(g, w);
    if (!(n2 < 0.0)) {
        o.valid = false;
        o.rejection_reason =
            "the bubble-centre coordinate trajectory is not timelike here: "
            "g(w,w) = " + std::to_string(n2) +
            " >= 0, so it is not a physical worldline";
        return o;
    }
    const double inv = 1.0 / std::sqrt(-n2);
    for (int mu = 0; mu < 4; ++mu) o.u[mu] = w[mu] * inv;
    o.norm = four_velocity_norm(g, o.u);
    o.normalization_residual = std::fabs(o.norm + 1.0);
    o.valid = o.normalization_residual < 1e-10;
    return o;
}

// u^a = gamma (1, v) in the ORTHONORMAL frame, |v| < 1 strictly.
inline Observer make_tetrad_boosted(const Mat4<double>& g, const Tetrad& t,
                                    const double v[3]) {
    Observer o;
    o.kind = ObserverKind::TetradBoosted;
    if (!t.valid) {
        o.valid = false;
        o.rejection_reason = "no tetrad: " + t.rejection_reason;
        return o;
    }
    const double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    if (!(v2 < 1.0)) {
        o.valid = false;
        o.rejection_reason =
            "boost speed |v| = " + std::to_string(std::sqrt(v2)) +
            " is at or beyond the unit sphere; no such observer exists";
        return o;
    }
    const double gam = 1.0 / std::sqrt(1.0 - v2);
    for (int mu = 0; mu < 4; ++mu) {
        o.u[mu] = gam * t.e[0][mu];
        for (int i = 0; i < 3; ++i) o.u[mu] += gam * v[i] * t.e[i + 1][mu];
    }
    o.norm = four_velocity_norm(g, o.u);
    o.normalization_residual = std::fabs(o.norm + 1.0);
    o.valid = o.normalization_residual < 1e-9;
    if (!o.valid)
        o.rejection_reason = "normalisation failed: g(u,u) = " +
                             std::to_string(o.norm);
    return o;
}

// --- Frame components of a stress-energy tensor -------------------------------
//
//     T_{ab} = e_a^mu e_b^nu T_{mu nu}
//
// with a covariant T. rho = T_{00}, q_i = -T_{0i}, S_ij = T_{ij}.
struct FrameStress {
    double T[4][4] = {{0}};   // T_{ab}, orthonormal frame indices
    double rho = 0.0;         // energy density measured by e_0
    double q[3] = {0, 0, 0};  // momentum density / energy flux
    double S[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};  // spatial stress, row-major
    double principal[3] = {0, 0, 0};            // eigenvalues of S, ascending
    double trace_spatial = 0.0;
    double q_norm = 0.0;
};

// All three eigenvalues of a symmetric 3x3, ascending, closed form.
inline void symmetric3_eigenvalues(const double m[9], double out[3]) {
    const double p1 = m[1] * m[1] + m[2] * m[2] + m[5] * m[5];
    const double tr = m[0] + m[4] + m[8];
    if (p1 <= 0.0) {
        out[0] = m[0]; out[1] = m[4]; out[2] = m[8];
    } else {
        const double q = tr / 3.0;
        const double d0 = m[0] - q, d4 = m[4] - q, d8 = m[8] - q;
        const double p2 = d0 * d0 + d4 * d4 + d8 * d8 + 2.0 * p1;
        const double p = std::sqrt(p2 / 6.0);
        double b[9];
        for (int i = 0; i < 9; ++i) b[i] = m[i] / p;
        b[0] -= q / p; b[4] -= q / p; b[8] -= q / p;
        const double dt =
            b[0] * (b[4] * b[8] - b[5] * b[7]) -
            b[1] * (b[3] * b[8] - b[5] * b[6]) +
            b[2] * (b[3] * b[7] - b[4] * b[6]);
        double r = dt / 2.0;
        if (r < -1.0) r = -1.0;
        if (r > 1.0) r = 1.0;
        const double phi = std::acos(r) / 3.0;
        const double kPi = 3.14159265358979323846;
        out[2] = q + 2.0 * p * std::cos(phi);
        out[0] = q + 2.0 * p * std::cos(phi + (2.0 * kPi / 3.0));
        out[1] = 3.0 * q - out[0] - out[2];
    }
    // Ascending, by fixed comparisons -- no sort with a comparator, so the
    // ordering is reproducible bit for bit.
    for (int i = 0; i < 2; ++i)
        for (int j = i + 1; j < 3; ++j)
            if (out[j] < out[i]) { const double s = out[i]; out[i] = out[j]; out[j] = s; }
}

inline FrameStress project_to_frame(const Tetrad& t, const Mat4<double>& T_cov) {
    FrameStress f;
    for (int A = 0; A < 4; ++A)
        for (int B = 0; B < 4; ++B) {
            double s = 0.0;
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu)
                    s += t.e[A][mu] * t.e[B][nu] * T_cov(mu, nu);
            f.T[A][B] = s;
        }

    f.rho = f.T[0][0];
    // Sign convention, stated once and used everywhere: q_i = -T_{0i}, so that
    // q points in the direction energy flows for ordinary matter.
    for (int i = 0; i < 3; ++i) f.q[i] = -f.T[0][i + 1];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) f.S[3 * i + j] = f.T[i + 1][j + 1];

    f.trace_spatial = f.S[0] + f.S[4] + f.S[8];
    f.q_norm = std::sqrt(f.q[0] * f.q[0] + f.q[1] * f.q[1] + f.q[2] * f.q[2]);
    symmetric3_eigenvalues(f.S, f.principal);
    return f;
}

}  // namespace geometry
}  // namespace warp_hypad
