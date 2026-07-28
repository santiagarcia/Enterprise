#pragma once

#include <string>

#include "warp_hypad/core/field.hpp"
#include "warp_hypad/core/grid.hpp"
#include "warp_hypad/core/mat4.hpp"
#include "warp_hypad/core/scalar_traits.hpp"
#include "warp_hypad/core/tensor.hpp"
#include "warp_hypad/geometry/chart.hpp"
#include "warp_hypad/geometry/finite_difference.hpp"

// ---------------------------------------------------------------------------
// Milestone A4: the covariant curvature pipeline.
//
//     g_{mu nu}  ->  Gamma^rho_{mu nu}  ->  R^rho_{sigma mu nu}  ->  R_{mu nu}
//                ->  R  ->  R_{mu nu} R^{mu nu},  R_{abcd} R^{abcd}
//
// Conventions are fixed in docs/conventions.md sections 5-7 and are NOT restated
// as choices here; this file implements a prior commitment.
//
//   Gamma^rho_{mu nu} = 1/2 g^{rho sigma}
//                       ( d_mu g_{sigma nu} + d_nu g_{sigma mu} - d_sigma g_{mu nu} )
//
//   R^rho_{sigma mu nu} = d_mu Gamma^rho_{nu sigma} - d_nu Gamma^rho_{mu sigma}
//                       + Gamma^rho_{mu lambda} Gamma^lambda_{nu sigma}
//                       - Gamma^rho_{nu lambda} Gamma^lambda_{mu sigma}
//
//   R_{sigma nu} = R^rho_{sigma rho nu},   R = g^{mu nu} R_{mu nu}
//
// MTW/Wald sign classification (+, +, +) with signature (-, +, +, +).
//
// --- On time -----------------------------------------------------------------
//
// Curvature needs d_0 g and d_0 Gamma. Time is not gridded (see
// temporal_derivative.hpp), so those come from re-evaluating the metric at
// shifted times. HOW that is done is not left implicit: every call carries a
// TemporalPolicy that states the claim being made and, where the claim is
// checkable, is checked before any curvature is produced.
//
// The frozen-time treatment is REFUSED rather than warned about. Imposing
// d_0 g = 0 on a metric that genuinely depends on t does not give the curvature
// of some other static spacetime; it gives an inconsistent input to the
// Christoffel formula. See docs/conventions.md section 11 and
// docs/time_dependence_decision.md section 3.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace geometry {

// --- Temporal derivative policy ---------------------------------------------

enum class TemporalPolicyKind {
    // Curvature must not be computed. The metric depends on t and the caller
    // has chosen to discard that dependence.
    Refused,

    // d_0 (anything built pointwise from g) is exactly zero because the metric
    // is independent of t IN THIS CHART. This is a checkable claim, not an
    // approximation, and compute_curvature() checks it: the sampled metric field
    // must be bit-identical at t and at a shifted time, or the policy is
    // downgraded to a refusal.
    ExactlyStatic,

    // d_0 obtained by re-evaluating the metric at shifted times and applying the
    // same second-order centred stencil used spatially. Makes no assumption
    // about the form of the t dependence.
    AnalyticResample
};

struct TemporalPolicy {
    TemporalPolicyKind kind = TemporalPolicyKind::Refused;

    // Step for AnalyticResample, and the staticity probe offset for
    // ExactlyStatic.
    double h = 0.0;

    // Recorded in metadata and in the refusal message. A policy without a stated
    // justification is a policy nobody can audit.
    std::string justification;

    static TemporalPolicy refused(std::string why) {
        TemporalPolicy p;
        p.kind = TemporalPolicyKind::Refused;
        p.justification = std::move(why);
        return p;
    }
    static TemporalPolicy exactly_static(std::string why, double probe_h = 1.0) {
        TemporalPolicy p;
        p.kind = TemporalPolicyKind::ExactlyStatic;
        p.h = probe_h;
        p.justification = std::move(why);
        return p;
    }
    static TemporalPolicy analytic_resample(double h, std::string why) {
        TemporalPolicy p;
        p.kind = TemporalPolicyKind::AnalyticResample;
        p.h = h;
        p.justification = std::move(why);
        return p;
    }

    const char* kind_name() const {
        switch (kind) {
            case TemporalPolicyKind::Refused: return "refused";
            case TemporalPolicyKind::ExactlyStatic: return "exactly_static";
            case TemporalPolicyKind::AnalyticResample: return "analytic_resample";
        }
        return "unknown";
    }
};

// --- Pointwise assembly ------------------------------------------------------

// Gamma^rho_{mu nu} from the inverse metric and dg(alpha, mu, nu) = d_alpha g_{mu nu}.
//
// Routed through contract_first() rather than an inlined quadruple loop so that
// the "first kind -> second kind" step is the one already verified by
// contract.contract_first_* in Milestone 1.
template <typename Scalar>
Rank3<Scalar> christoffel_at(const Mat4<Scalar>& g_inv, const Rank3<Scalar>& dg) {
    // D_{sigma mu nu} = d_mu g_{sigma nu} + d_nu g_{sigma mu} - d_sigma g_{mu nu}
    Rank3<Scalar> d;
    for (int s = 0; s < 4; ++s)
        for (int m = 0; m < 4; ++m)
            for (int n = 0; n < 4; ++n)
                d(s, m, n) = dg(m, s, n) + dg(n, s, m) - dg(s, m, n);

    Rank3<Scalar> gamma = contract_first(g_inv, d);
    const Scalar half = scalar_cast<Scalar>(0.5);
    for (int i = 0; i < 64; ++i) gamma.t[i] = gamma.t[i] * half;
    return gamma;
}

// R^rho_{sigma mu nu} from Gamma and dgamma(alpha, rho, mu, nu) = d_alpha Gamma^rho_{mu nu}.
template <typename Scalar>
Rank4<Scalar> riemann_at(const Rank3<Scalar>& gamma, const Rank4<Scalar>& dgamma) {
    Rank4<Scalar> r;
    for (int rho = 0; rho < 4; ++rho)
        for (int sig = 0; sig < 4; ++sig)
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu) {
                    Scalar acc = dgamma(mu, rho, nu, sig) - dgamma(nu, rho, mu, sig);
                    for (int lam = 0; lam < 4; ++lam) {
                        acc = acc + gamma(rho, mu, lam) * gamma(lam, nu, sig) -
                              gamma(rho, nu, lam) * gamma(lam, mu, sig);
                    }
                    r(rho, sig, mu, nu) = acc;
                }
    return r;
}

// Metric compatibility residual at a point:
//
//     nabla_lambda g_{mu nu}
//         = d_lambda g_{mu nu} - Gamma^rho_{lambda mu} g_{rho nu}
//                              - Gamma^rho_{lambda nu} g_{mu rho}
//
// Identically zero for the Levi-Civita connection. This does NOT reuse the
// Christoffel formula, so it is a genuine Level A check rather than a
// restatement of the definition -- the point made in docs/conventions.md
// section 5.
template <typename Scalar>
double metric_compatibility_residual(const Mat4<Scalar>& g, const Rank3<Scalar>& dg,
                                     const Rank3<Scalar>& gamma) {
    double worst = 0.0;
    for (int lam = 0; lam < 4; ++lam)
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) {
                Scalar acc = dg(lam, mu, nu);
                for (int rho = 0; rho < 4; ++rho) {
                    acc = acc - gamma(rho, lam, mu) * g(rho, nu) -
                          gamma(rho, lam, nu) * g(mu, rho);
                }
                const double a = real_abs(acc);
                if (a > worst) worst = a;
            }
    return worst;
}

// --- Field-level assembly ----------------------------------------------------

template <typename Scalar>
Rank3Field<Scalar> christoffel_field(const Grid& grid,
                                     const Rank2Field<Scalar>& g_inv,
                                     const Rank3Field<Scalar>& dg) {
    Rank3Field<Scalar> out(grid);
    for (std::size_t p = 0; p < grid.num_points(); ++p)
        out.set_tensor(p, christoffel_at(g_inv.matrix_at(p), dg.tensor_at(p)));
    return out;
}

template <typename Scalar>
Rank4Field<Scalar> riemann_field(const Grid& grid, const Rank3Field<Scalar>& gamma,
                                 const Rank4Field<Scalar>& dgamma) {
    Rank4Field<Scalar> out(grid);
    for (std::size_t p = 0; p < grid.num_points(); ++p)
        out.set_tensor(p, riemann_at(gamma.tensor_at(p), dgamma.tensor_at(p)));
    return out;
}

// --- Curvature result --------------------------------------------------------

template <typename Scalar>
struct CurvatureFields {
    Rank2Field<Scalar> g;
    Rank2Field<Scalar> g_inv;
    Rank3Field<Scalar> dg;       // d_alpha g_{mu nu}
    Rank3Field<Scalar> gamma;    // Gamma^rho_{mu nu}
    Rank4Field<Scalar> dgamma;   // d_alpha Gamma^rho_{mu nu}
    Rank4Field<Scalar> riemann;  // R^rho_{sigma mu nu}
    Rank2Field<Scalar> ricci;    // R_{mu nu}

    ScalarField<Scalar> ricci_scalar;  // R
    ScalarField<Scalar> ricci_square;  // R_{mu nu} R^{mu nu}
    ScalarField<Scalar> kretschmann;   // R_{abcd} R^{abcd}

    // False means nothing above may be read. The refusal is a first-class
    // outcome, not an exception: a parameter sweep that wanders into a
    // configuration curvature cannot be computed for should produce a labelled
    // gap, not a crash.
    bool ok = false;
    std::string refusal_reason;

    // Filled when the policy was ExactlyStatic: the largest difference found
    // between the metric sampled at t and at the probe time. Nonzero means the
    // staticity claim was false and the run was refused.
    double staticity_probe_residual = 0.0;

    // Worst |g g^-1 - I| over the grid, and the count of points where the
    // inverse could not be formed at all.
    double worst_inverse_residual = 0.0;
    std::size_t num_singular_points = 0;
};

// Derived invariants from riemann/ricci/g_inv. Separated so a caller that only
// wants scalars does not pay for storing them and a caller that wants both does
// not compute them twice.
template <typename Scalar>
void fill_invariants(const Grid& grid, CurvatureFields<Scalar>& c) {
    c.ricci_scalar = ScalarField<Scalar>(grid);
    c.ricci_square = ScalarField<Scalar>(grid);
    c.kretschmann = ScalarField<Scalar>(grid);

    for (std::size_t p = 0; p < grid.num_points(); ++p) {
        const Mat4<Scalar> gi = c.g_inv.matrix_at(p);
        const Mat4<Scalar> gm = c.g.matrix_at(p);
        const Mat4<Scalar> ric = c.ricci.matrix_at(p);
        const Rank4<Scalar> riem = c.riemann.tensor_at(p);

        c.ricci_scalar.at(p) = full_contract(gi, ric);

        // R^{mu nu} = g^{mu a} R_{a b} g^{b nu}; then R_{mu nu} R^{mu nu}.
        const Mat4<Scalar> ric_up = mat_mul(mat_mul(gi, ric), gi);
        c.ricci_square.at(p) = full_contract(ric_up, ric);

        // K = R_{abcd} R^{abcd}. The stored Riemann has its first slot raised,
        // so lower it before raising all four.
        const Rank4<Scalar> riem_down = lower_first(gm, riem);
        const Rank4<Scalar> riem_up = raise_all_rank4(gi, riem_down);
        c.kretschmann.at(p) = full_contract_rank4(riem_down, riem_up);
    }
}

namespace detail {

// d_0 g as a Rank3Field carrying data in the alpha = 0 slot only, formed as a
// centred difference of two metric taps separated by 2h.
template <typename Scalar>
Rank3Field<Scalar> centred_time_slot(const Grid& grid, const Rank2Field<Scalar>& minus,
                                     const Rank2Field<Scalar>& plus, double h) {
    Rank3Field<Scalar> out(grid);
    const Scalar inv = scalar_cast<Scalar>(1.0 / (2.0 * h));
    for (std::size_t p = 0; p < grid.num_points(); ++p)
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu)
                out.at(p, 0, mu, nu) = (plus.at(p, mu, nu) - minus.at(p, mu, nu)) * inv;
    return out;
}

// The same for Gamma, producing the alpha = 0 slot of dgamma.
template <typename Scalar>
Rank4Field<Scalar> centred_time_slot_rank3(const Grid& grid,
                                           const Rank3Field<Scalar>& minus,
                                           const Rank3Field<Scalar>& plus, double h) {
    Rank4Field<Scalar> out(grid);
    const Scalar inv = scalar_cast<Scalar>(1.0 / (2.0 * h));
    for (std::size_t p = 0; p < grid.num_points(); ++p)
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int c = 0; c < 4; ++c)
                    out.at(p, 0, a, b, c) =
                        (plus.at(p, a, b, c) - minus.at(p, a, b, c)) * inv;
    return out;
}

template <typename Scalar>
Rank2Field<Scalar> invert_field(const Grid& grid, const Rank2Field<Scalar>& g,
                                double& worst_residual, std::size_t& num_singular) {
    Rank2Field<Scalar> out(grid);
    worst_residual = 0.0;
    num_singular = 0;
    for (std::size_t p = 0; p < grid.num_points(); ++p) {
        const Mat4<Scalar> m = g.matrix_at(p);
        const InverseResult<Scalar> inv = inverse(m);
        if (!inv.ok) {
            ++num_singular;
            continue;
        }
        out.set_matrix(p, inv.inverse);
        const double r = identity_residual(m, inv.inverse);
        if (r > worst_residual) worst_residual = r;
    }
    return out;
}

// Largest |real| componentwise difference between two rank-2 fields.
template <typename Scalar>
double field_max_abs_diff(const Grid& grid, const Rank2Field<Scalar>& a,
                          const Rank2Field<Scalar>& b) {
    double worst = 0.0;
    for (std::size_t p = 0; p < grid.num_points(); ++p)
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) {
                const double d = real_part(a.at(p, mu, nu)) - real_part(b.at(p, mu, nu));
                const double ad = d < 0.0 ? -d : d;
                if (ad > worst) worst = ad;
            }
    return worst;
}

}  // namespace detail

// --- The driver --------------------------------------------------------------
//
// `sample` is any callable double -> Rank2Field<Scalar> returning the metric
// field over `grid` at the given coordinate time.
//
// Cost. ExactlyStatic evaluates the metric twice (once for the staticity probe)
// and differentiates twice. AnalyticResample evaluates it five times, at
// t + {-2h, -h, 0, +h, +2h}, forms d_0 g at the three inner taps, builds Gamma at
// those three taps, and centre-differences Gamma in time. That is the price of
// making no assumption about the form of the t dependence; the rigid-translation
// shortcut d_0 = -v_s d_x is available for the Alcubierre family and is used in
// the A4 suite as an INDEPENDENT check of this path rather than as the path.
template <typename Scalar, typename SampleMetric>
CurvatureFields<Scalar> compute_curvature(const Grid& grid, SampleMetric&& sample,
                                          double t, const TemporalPolicy& policy) {
    CurvatureFields<Scalar> c;

    if (policy.kind == TemporalPolicyKind::Refused) {
        c.ok = false;
        c.refusal_reason =
            "temporal policy is Refused: " + policy.justification +
            " -- curvature built from a discarded temporal derivative is not the "
            "curvature of any spacetime";
        return c;
    }

    c.g = sample(t);
    c.g_inv = detail::invert_field<Scalar>(grid, c.g, c.worst_inverse_residual,
                                           c.num_singular_points);
    if (c.num_singular_points != 0) {
        c.ok = false;
        c.refusal_reason = "metric is singular at " +
                           std::to_string(c.num_singular_points) +
                           " grid point(s); no connection exists there";
        return c;
    }

    if (policy.kind == TemporalPolicyKind::ExactlyStatic) {
        // Check the claim before relying on it. A metric that quietly depends on
        // t would otherwise be silently frozen -- exactly the failure mode the
        // refusal above exists to prevent, arriving by a different door.
        const Rank2Field<Scalar> probe = sample(t + policy.h);
        c.staticity_probe_residual = detail::field_max_abs_diff(grid, c.g, probe);
        if (c.staticity_probe_residual != 0.0) {
            c.ok = false;
            c.refusal_reason =
                "ExactlyStatic was claimed (" + policy.justification +
                ") but the metric changed by " +
                std::to_string(c.staticity_probe_residual) + " between t and t+h";
            return c;
        }

        c.dg = spacetime_gradient<Scalar>(c.g, grid, nullptr);
        c.gamma = christoffel_field(grid, c.g_inv, c.dg);
        c.dgamma = spacetime_gradient_rank3<Scalar>(c.gamma, grid, nullptr);
    } else {
        const double h = policy.h;

        const Rank2Field<Scalar> g_m2 = sample(t - 2.0 * h);
        const Rank2Field<Scalar> g_m1 = sample(t - h);
        const Rank2Field<Scalar> g_p1 = sample(t + h);
        const Rank2Field<Scalar> g_p2 = sample(t + 2.0 * h);

        // Gamma at the three inner taps. Each needs its own d_0 g, which is why
        // the outer two taps exist.
        auto gamma_at_tap = [&](const Rank2Field<Scalar>& gm,
                                const Rank2Field<Scalar>& before,
                                const Rank2Field<Scalar>& after) {
            double res = 0.0;
            std::size_t sing = 0;
            const Rank2Field<Scalar> gi =
                detail::invert_field<Scalar>(grid, gm, res, sing);
            const Rank3Field<Scalar> d0 =
                detail::centred_time_slot(grid, before, after, h);
            const Rank3Field<Scalar> dgt = spacetime_gradient<Scalar>(gm, grid, &d0);
            return christoffel_field(grid, gi, dgt);
        };

        const Rank3Field<Scalar> gamma_m1 = gamma_at_tap(g_m1, g_m2, c.g);
        const Rank3Field<Scalar> gamma_p1 = gamma_at_tap(g_p1, c.g, g_p2);

        const Rank3Field<Scalar> d0g = detail::centred_time_slot(grid, g_m1, g_p1, h);
        c.dg = spacetime_gradient<Scalar>(c.g, grid, &d0g);
        c.gamma = christoffel_field(grid, c.g_inv, c.dg);

        const Rank4Field<Scalar> d0gamma =
            detail::centred_time_slot_rank3(grid, gamma_m1, gamma_p1, h);
        c.dgamma = spacetime_gradient_rank3<Scalar>(c.gamma, grid, &d0gamma);
    }

    c.riemann = riemann_field(grid, c.gamma, c.dgamma);

    c.ricci = Rank2Field<Scalar>(grid);
    for (std::size_t p = 0; p < grid.num_points(); ++p)
        c.ricci.set_matrix(p, contract_13(c.riemann.tensor_at(p)));

    fill_invariants(grid, c);

    c.ok = true;
    return c;
}

// --- Interior/boundary reporting ---------------------------------------------
//
// Reported separately, always, for the same reason the A2 suite does it: a
// one-sided stencil has a larger error constant, and curvature uses TWO nested
// derivative passes, so the boundary shell is where a plot goes wrong first.
// A caller that averages the two is measuring its own boundary treatment.

struct FieldError {
    double interior_worst = 0.0;
    double boundary_worst = 0.0;
    std::size_t interior_points = 0;
    std::size_t boundary_points = 0;
};

// Points within `depth` of a face are counted as boundary. depth = 2 is the
// right default for curvature: the second derivative pass reaches two points in,
// so a point one step from a face is contaminated by the face's one-sided
// stencil even though its own stencil is centred.
inline bool is_shell_point(const Grid& g, std::size_t p, int depth) {
    int i, j, k;
    g.unflatten(p, i, j, k);
    return i < depth || i >= g.nx() - depth || j < depth || j >= g.ny() - depth ||
           k < depth || k >= g.nz() - depth;
}

template <typename Scalar, typename PointError>
FieldError split_error(const Grid& g, PointError&& err, int depth = 2) {
    FieldError fe;
    for (std::size_t p = 0; p < g.num_points(); ++p) {
        const double e = err(p);
        if (is_shell_point(g, p, depth)) {
            ++fe.boundary_points;
            if (e > fe.boundary_worst) fe.boundary_worst = e;
        } else {
            ++fe.interior_points;
            if (e > fe.interior_worst) fe.interior_worst = e;
        }
    }
    return fe;
}

}  // namespace geometry
}  // namespace warp_hypad
