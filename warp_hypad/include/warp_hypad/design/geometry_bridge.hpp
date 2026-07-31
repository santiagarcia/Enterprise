#pragma once

#include <cmath>
#include <string>

#include "warp_hypad/core/mat4.hpp"
#include "warp_hypad/core/tensor.hpp"
#include "warp_hypad/design/adm_design.hpp"
#include "warp_hypad/geometry/curvature.hpp"
#include "warp_hypad/geometry/einstein.hpp"

// ---------------------------------------------------------------------------
// O0C: the pointwise geometry bridge.
//
// One job: turn a design MetricJet -- g, dg, d2g, d3g, all analytic -- into
// curvature and the Einstein tensor, WITHOUT a finite difference anywhere.
//
// WHY THIS IS AN ADAPTER AND NOT A SECOND GR IMPLEMENTATION. A4 already has a
// verified Christoffel kernel, a verified Riemann kernel and verified symmetry
// residuals; A5 has a verified Einstein kernel and the source decomposition.
// What A4 does NOT have is an analytic d_alpha Gamma: its grid route forms that
// by centred differences of Gamma, which is exactly the thing O0C must not do.
// So the only genuinely new content here is the CHAIN FROM d2g TO dGamma, and
// even that is expressed through A4's own christoffel_at rather than as a fresh
// formula. Everything downstream -- riemann_at, contract_13, einstein_at,
// full_contract, lower_first, raise_all_rank4, the four symmetry residuals --
// is called, not reimplemented.
//
// THE ONE IDENTITY THAT MAKES THE REUSE WORK. A4's
//
//     christoffel_at(A, dg)^rho_{mu nu}
//         = (1/2) A^{rho sigma} [ dg(mu, sigma, nu) + dg(nu, sigma, mu)
//                                 - dg(sigma, mu, nu) ]
//
// is bilinear: linear in the matrix A and linear in the rank-3 argument. The
// Christoffel symbol is that expression at (g^{-1}, dg). Differentiating it is
// therefore a product rule over the two slots, and each resulting term is
// ANOTHER call to the same function:
//
//     d_lambda Gamma = christoffel_at( d_lambda g^{-1},  dg              )
//                    + christoffel_at( g^{-1},           d_lambda dg     )
//
//     d_lambda d_kappa Gamma
//                    = christoffel_at( d_l d_k g^{-1},   dg              )
//                    + christoffel_at( d_l g^{-1},       d_k dg          )
//                    + christoffel_at( d_k g^{-1},       d_l dg          )
//                    + christoffel_at( g^{-1},           d_l d_k dg      )
//
// No new connection algebra is written down, so no new connection algebra can
// be wrong. The inverse-metric derivatives come from differentiating
// g^{-1} g = I, which is the only place a formula is introduced:
//
//     d_l g^{rho sigma} = - g^{rho a} (d_l g_{ab}) g^{b sigma}
//     d_l d_k g^{rho sigma} =
//           g^{ra}(d_l g_{ab}) g^{bc}(d_k g_{cd}) g^{ds}
//         + g^{ra}(d_k g_{ab}) g^{bc}(d_l g_{cd}) g^{ds}
//         - g^{ra}(d_l d_k g_{ab}) g^{bs}
//
// and the gate checks both against the metric-inverse identity directly.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

// --- Inverse-metric derivatives ---------------------------------------------

// d_lambda g^{rho sigma} = - g^{rho a} (d_lambda g_{ab}) g^{b sigma}.
template <typename Scalar>
Mat4<Scalar> inverse_metric_first_derivative(const Mat4<Scalar>& g_inv,
                                             const Mat4<Scalar>& dg) {
    Mat4<Scalar> out;
    for (int r = 0; r < 4; ++r)
        for (int s = 0; s < 4; ++s) {
            Scalar acc = scalar_cast<Scalar>(0.0);
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    acc = acc + g_inv(r, a) * dg(a, b) * g_inv(b, s);
            out(r, s) = -acc;
        }
    return out;
}

// d_lambda d_kappa g^{rho sigma}, from one more derivative of the same identity.
template <typename Scalar>
Mat4<Scalar> inverse_metric_second_derivative(const Mat4<Scalar>& g_inv,
                                              const Mat4<Scalar>& dg_l,
                                              const Mat4<Scalar>& dg_k,
                                              const Mat4<Scalar>& d2g_lk) {
    Mat4<Scalar> out;
    for (int r = 0; r < 4; ++r)
        for (int s = 0; s < 4; ++s) {
            Scalar plus = scalar_cast<Scalar>(0.0);
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    for (int c = 0; c < 4; ++c)
                        for (int d = 0; d < 4; ++d)
                            plus = plus + g_inv(r, a) * dg_l(a, b) *
                                              g_inv(b, c) * dg_k(c, d) *
                                              g_inv(d, s) +
                                   g_inv(r, a) * dg_k(a, b) * g_inv(b, c) *
                                       dg_l(c, d) * g_inv(d, s);
            Scalar minus = scalar_cast<Scalar>(0.0);
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b)
                    minus = minus + g_inv(r, a) * d2g_lk(a, b) * g_inv(b, s);
            out(r, s) = plus - minus;
        }
    return out;
}

// --- Jet repacking -----------------------------------------------------------
//
// The design jet stores derivatives as arrays of CovariantMetric; A4's kernels
// want Rank3 with the derivative index leading. Repacking is not a computation
// and is kept separate so that neither side has to learn the other's layout.

template <typename Scalar>
Rank3<Scalar> pack_dg(const MetricJet<Scalar>& jet) {
    Rank3<Scalar> out;
    for (int a = 0; a < 4; ++a)
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu)
                out(a, mu, nu) = jet.dg[static_cast<std::size_t>(a)].m(mu, nu);
    return out;
}

// The lambda-slice of the second derivatives, as a Rank3 whose leading index is
// the REMAINING derivative: slice(a, mu, nu) = d_lambda d_a g_{mu nu}.
template <typename Scalar>
Rank3<Scalar> pack_d2g_slice(const MetricJet<Scalar>& jet, int lambda) {
    Rank3<Scalar> out;
    for (int a = 0; a < 4; ++a)
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu)
                out(a, mu, nu) =
                    jet.d2g[static_cast<std::size_t>(lambda)]
                           [static_cast<std::size_t>(a)].m(mu, nu);
    return out;
}

template <typename Scalar>
Rank3<Scalar> pack_d3g_slice(const MetricJet<Scalar>& jet, int lambda,
                             int kappa) {
    Rank3<Scalar> out;
    for (int a = 0; a < 4; ++a)
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu)
                out(a, mu, nu) =
                    jet.d3g[static_cast<std::size_t>(lambda)]
                           [static_cast<std::size_t>(kappa)]
                           [static_cast<std::size_t>(a)].m(mu, nu);
    return out;
}

// --- Connection jet ----------------------------------------------------------

template <typename Scalar>
struct ConnectionJet {
    Rank3<Scalar> gamma;             // Gamma^rho_{mu nu}
    Rank4<Scalar> dgamma;            // d_alpha Gamma^rho_{mu nu}
    Rank3<Scalar> d2gamma[4][4];     // d_alpha d_beta Gamma^rho_{mu nu}
    bool second_available = false;
};

// Gamma and d Gamma, analytically. `want_second` additionally fills d2gamma,
// which costs four more christoffel_at calls per index pair and is only needed
// by the differential-Bianchi diagnostic.
template <typename Scalar>
ConnectionJet<Scalar> connection_jet(const MetricJet<Scalar>& jet,
                                     bool want_second = false) {
    using geometry::christoffel_at;
    ConnectionJet<Scalar> c;

    const Mat4<Scalar>& gi = jet.g_contra.m;
    const Rank3<Scalar> dg = pack_dg(jet);
    c.gamma = christoffel_at(gi, dg);

    Mat4<Scalar> dgi[4];
    for (int l = 0; l < 4; ++l)
        dgi[l] = inverse_metric_first_derivative(
            gi, jet.dg[static_cast<std::size_t>(l)].m);

    Rank3<Scalar> d2slice[4];
    for (int l = 0; l < 4; ++l) d2slice[l] = pack_d2g_slice(jet, l);

    for (int l = 0; l < 4; ++l) {
        const Rank3<Scalar> a = christoffel_at(dgi[l], dg);
        const Rank3<Scalar> b = christoffel_at(gi, d2slice[l]);
        for (int rho = 0; rho < 4; ++rho)
            for (int mu = 0; mu < 4; ++mu)
                for (int nu = 0; nu < 4; ++nu)
                    c.dgamma(l, rho, mu, nu) =
                        a(rho, mu, nu) + b(rho, mu, nu);
    }

    if (!want_second) return c;

    for (int l = 0; l < 4; ++l)
        for (int k = 0; k < 4; ++k) {
            const Mat4<Scalar> d2gi = inverse_metric_second_derivative(
                gi, jet.dg[static_cast<std::size_t>(l)].m,
                jet.dg[static_cast<std::size_t>(k)].m,
                jet.d2g[static_cast<std::size_t>(l)]
                       [static_cast<std::size_t>(k)].m);
            const Rank3<Scalar> t1 = christoffel_at(d2gi, dg);
            const Rank3<Scalar> t2 = christoffel_at(dgi[l], d2slice[k]);
            const Rank3<Scalar> t3 = christoffel_at(dgi[k], d2slice[l]);
            const Rank3<Scalar> t4 =
                christoffel_at(gi, pack_d3g_slice(jet, l, k));
            for (int i = 0; i < 64; ++i)
                c.d2gamma[l][k].t[i] =
                    t1.t[i] + t2.t[i] + t3.t[i] + t4.t[i];
        }
    c.second_available = true;
    return c;
}

// --- Pointwise geometry ------------------------------------------------------

struct PointGeometry {
    Mat4<double> g;
    Mat4<double> g_inv;
    Rank3<double> gamma;
    Rank4<double> riemann;        // R^rho_{sigma mu nu}
    Rank4<double> riemann_down;   // R_{rho sigma mu nu}
    Mat4<double> ricci;
    Mat4<double> einstein;        // G_{mu nu}
    double ricci_scalar = 0.0;
    double ricci_square = 0.0;
    double kretschmann = 0.0;
    double einstein_trace = 0.0;

    // Diagnostics, all pointwise and all from A4's own residual functions.
    double antisymmetry_first_pair = 0.0;
    double antisymmetry_last_pair = 0.0;
    double pair_exchange = 0.0;
    double first_bianchi = 0.0;
    double metric_compatibility = 0.0;
    double inverse_residual = 0.0;
    double trace_identity = 0.0;   // |G^mu_mu + R|

    bool valid = false;
    std::string rejection;
};

// The complete pointwise chain. Every step is an existing kernel except the
// connection jet above.
inline PointGeometry point_geometry(const MetricJet<double>& jet) {
    PointGeometry p;
    if (!jet.valid) {
        p.rejection = "the design metric jet is invalid at this point";
        return p;
    }
    p.g = jet.g_cov.m;
    p.g_inv = jet.g_contra.m;

    // g g^{-1} - I, measured rather than assumed: the analytic inverse is only
    // exact if alpha and gamma are what the jet says they are.
    {
        double worst = 0.0;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) {
                double acc = 0.0;
                for (int k = 0; k < 4; ++k) acc += p.g(i, k) * p.g_inv(k, j);
                worst = std::fmax(worst, std::fabs(acc - (i == j ? 1.0 : 0.0)));
            }
        p.inverse_residual = worst;
    }

    const ConnectionJet<double> c = connection_jet(jet, /*want_second=*/false);
    p.gamma = c.gamma;
    p.riemann = geometry::riemann_at(c.gamma, c.dgamma);
    p.ricci = contract_13(p.riemann);
    p.ricci_scalar = full_contract(p.g_inv, p.ricci);
    p.einstein = geometry::einstein_at(p.g, p.ricci, p.ricci_scalar);
    p.einstein_trace = full_contract(p.g_inv, p.einstein);
    p.trace_identity = std::fabs(p.einstein_trace + p.ricci_scalar);

    const Mat4<double> ric_up = mat_mul(mat_mul(p.g_inv, p.ricci), p.g_inv);
    p.ricci_square = full_contract(ric_up, p.ricci);

    p.riemann_down = lower_first(p.g, p.riemann);
    const Rank4<double> riem_up = raise_all_rank4(p.g_inv, p.riemann_down);
    p.kretschmann = full_contract_rank4(p.riemann_down, riem_up);

    p.antisymmetry_first_pair = antisymmetry_first_pair_residual(p.riemann_down);
    p.antisymmetry_last_pair = antisymmetry_last_pair_residual(p.riemann_down);
    p.pair_exchange = pair_exchange_residual(p.riemann_down);
    p.first_bianchi = first_bianchi_residual(p.riemann_down);
    p.metric_compatibility =
        geometry::metric_compatibility_residual(p.g, pack_dg(jet), c.gamma);

    p.valid = true;
    return p;
}

// --- Differential Bianchi diagnostic ------------------------------------------
//
// nabla_mu G^{mu nu}, which the contracted Bianchi identity makes zero in the
// continuum. Here it is an ANALYTIC quantity, not a truncation statement, so a
// nonzero value means the third-derivative chain is wrong rather than that the
// grid is coarse. It is still reported as a diagnostic and never used as an
// optimizer constraint: the spec's own instruction, and the right one, because
// its convergence law has not been demonstrated.
//
// This is the only place a curvature derivative is written out longhand, since
// A4 has no analytic d Riemann to call.
struct BianchiDiagnostic {
    double divergence[4] = {0, 0, 0, 0};
    double worst = 0.0;
    double einstein_scale = 0.0;
    double relative = 0.0;
    bool available = false;
};

inline BianchiDiagnostic bianchi_diagnostic(const MetricJet<double>& jet) {
    BianchiDiagnostic d;
    if (!jet.valid) return d;

    const ConnectionJet<double> c = connection_jet(jet, /*want_second=*/true);
    const Mat4<double>& gi = jet.g_contra.m;
    const Mat4<double>& g = jet.g_cov.m;

    const Rank4<double> riem = geometry::riemann_at(c.gamma, c.dgamma);
    const Mat4<double> ric = contract_13(riem);
    const double R = full_contract(gi, ric);
    const Mat4<double> G = geometry::einstein_at(g, ric, R);

    // d_kappa R^rho_{sigma mu nu}, by the product rule on riemann_at's body.
    Rank4<double> driem[4];
    for (int k = 0; k < 4; ++k)
        for (int rho = 0; rho < 4; ++rho)
            for (int sig = 0; sig < 4; ++sig)
                for (int mu = 0; mu < 4; ++mu)
                    for (int nu = 0; nu < 4; ++nu) {
                        double acc = c.d2gamma[k][mu](rho, nu, sig) -
                                     c.d2gamma[k][nu](rho, mu, sig);
                        for (int lam = 0; lam < 4; ++lam)
                            acc += c.dgamma(k, rho, mu, lam) *
                                       c.gamma(lam, nu, sig) +
                                   c.gamma(rho, mu, lam) *
                                       c.dgamma(k, lam, nu, sig) -
                                   c.dgamma(k, rho, nu, lam) *
                                       c.gamma(lam, mu, sig) -
                                   c.gamma(rho, nu, lam) *
                                       c.dgamma(k, lam, mu, sig);
                        driem[k](rho, sig, mu, nu) = acc;
                    }

    Mat4<double> dric[4];
    for (int k = 0; k < 4; ++k) dric[k] = contract_13(driem[k]);

    Mat4<double> dgi[4];
    for (int k = 0; k < 4; ++k)
        dgi[k] = inverse_metric_first_derivative(
            gi, jet.dg[static_cast<std::size_t>(k)].m);

    double dR[4] = {0, 0, 0, 0};
    for (int k = 0; k < 4; ++k)
        dR[k] = full_contract(dgi[k], ric) + full_contract(gi, dric[k]);

    // d_kappa G_{mu nu} = d Ric - (1/2)(d g) R - (1/2) g (d R).
    Mat4<double> dG[4];
    for (int k = 0; k < 4; ++k)
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu)
                dG[k](mu, nu) =
                    dric[k](mu, nu) -
                    0.5 * jet.dg[static_cast<std::size_t>(k)].m(mu, nu) * R -
                    0.5 * g(mu, nu) * dR[k];

    // G^{mu nu} and its derivative, then the covariant divergence.
    const Mat4<double> Gup = mat_mul(mat_mul(gi, G), gi);
    Mat4<double> dGup[4];
    for (int k = 0; k < 4; ++k) {
        const Mat4<double> a = mat_mul(mat_mul(dgi[k], G), gi);
        const Mat4<double> b = mat_mul(mat_mul(gi, dG[k]), gi);
        const Mat4<double> e = mat_mul(mat_mul(gi, G), dgi[k]);
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu)
                dGup[k](mu, nu) = a(mu, nu) + b(mu, nu) + e(mu, nu);
    }

    for (int nu = 0; nu < 4; ++nu) {
        double acc = 0.0;
        for (int mu = 0; mu < 4; ++mu) {
            acc += dGup[mu](mu, nu);
            for (int lam = 0; lam < 4; ++lam)
                acc += c.gamma(mu, mu, lam) * Gup(lam, nu) +
                       c.gamma(nu, mu, lam) * Gup(mu, lam);
        }
        d.divergence[nu] = acc;
        d.worst = std::fmax(d.worst, std::fabs(acc));
    }
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu)
            d.einstein_scale = std::fmax(d.einstein_scale, std::fabs(G(mu, nu)));
    d.relative = d.einstein_scale > 0.0 ? d.worst / d.einstein_scale : d.worst;
    d.available = true;
    return d;
}

}  // namespace design
}  // namespace warp_hypad
