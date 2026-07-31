#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "warp_hypad/design/coefficient_spaces.hpp"
#include "warp_hypad/design/geometry_bridge.hpp"
#include "warp_hypad/design/shift_family.hpp"
#include "warp_hypad/design/wall_family.hpp"
#include "warp_hypad/geometry/energy_conditions.hpp"
#include "warp_hypad/geometry/observer.hpp"

// ---------------------------------------------------------------------------
// O0C: observer-frame energy, energy-condition functionals and the integrated
// objectives, plus the scalar-generic path the optimizer coordinates travel.
//
// THE SCALAR-GENERIC JET, AND WHY IT IS NOT A SECOND IMPLEMENTATION. beta^x is
// AFFINE in the coefficients:
//
//     beta = beta_base(x) + sum_k c_k B_k(x)
//
// and every spacetime derivative of beta is the same affine form with B_k
// replaced by its own derivative, because chained_derivative() is linear in the
// design-derivative table it is handed. So a Scalar-valued jet is obtained by
// calling the EXISTING double kernels once for the baseline and once per active
// basis function, then accumulating with Scalar coefficients. No derivative
// formula is rewritten; the only thing that changes is the type of the
// accumulator. That is why an OTI gradient of the energy is a gradient of the
// same arithmetic the double path runs, and not of a parallel copy of it.
//
// THE VOLUME ELEMENT. Physical spatial volume is sqrt(gamma) d^3x. In the
// axisymmetric design chart d^3x = pi R0^3 dxi ds, so
//
//     dV_gamma = pi R0^3 sqrt(gamma) dxi ds
//
// and build_quadrature() already carries the pi R0^3 factor in its weights.
// sqrt(gamma) is therefore the ONLY thing this file multiplies in, and for the
// shift-only families gamma = I makes it exactly 1 -- which the gate measures
// rather than assumes, because the identity stops holding the moment a lapse or
// spatial-metric family is added.
//
// E_neg IS A POSITIVE MAGNITUDE. Written out once so no reader has to guess:
//
//     E_net = integral rho dV          (signed)
//     E_plus = integral max(rho, 0) dV  (positive)
//     E_neg  = integral max(-rho, 0) dV (POSITIVE: the magnitude of the
//                                        negative-energy content)
//     E_abs  = integral |rho| dV = E_plus + E_neg
//     E_net  = E_plus - E_neg
//
// A5/A6 called the third of these `E_minus`, also as a positive magnitude. The
// legacy comparison record states that mapping explicitly rather than relying
// on the reader to notice.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

// --- Scalar-generic shift jet -------------------------------------------------

template <typename Scalar>
struct ShiftJetT {
    Scalar value = scalar_cast<Scalar>(0.0);
    Scalar d1[4];
    Scalar d2[4][4];
    Scalar d3[4][4][4];

    ShiftJetT() {
        const Scalar z = scalar_cast<Scalar>(0.0);
        for (int i = 0; i < 4; ++i) {
            d1[i] = z;
            for (int j = 0; j < 4; ++j) {
                d2[i][j] = z;
                for (int k = 0; k < 4; ++k) d3[i][j][k] = z;
            }
        }
    }
};

// The per-basis-function design-derivative table, isolated so the accumulation
// below can be over Scalar coefficients. Same channels correction_shift_table
// fills, same univariate kernel for the two mixed thirds.
inline DesignDerivativeTable single_mode_table(const AxisymmetricBasis& basis,
                                               const TensorBasisEvaluation& e,
                                               const DesignCoordinates& p,
                                               std::size_t m) {
    DesignDerivativeTable t;
    t.at(0, 0) = e.value[m];
    t.at(1, 0) = e.d_xi[m];
    t.at(0, 1) = e.d_s[m];
    t.at(2, 0) = e.d_xi_xi[m];
    t.at(1, 1) = e.d_xi_s[m];
    t.at(0, 2) = e.d_s_s[m];
    t.at(3, 0) = e.d3_xi[m];
    t.at(0, 3) = e.d3_s[m];
    t.at(2, 1) = basis.axial().derivative_of(e.index_xi[m], p.xi, 2) *
                 basis.radial().derivative_of(e.index_s[m], p.radial_squared, 1);
    t.at(1, 2) = basis.axial().derivative_of(e.index_xi[m], p.xi, 1) *
                 basis.radial().derivative_of(e.index_s[m], p.radial_squared, 2);
    return t;
}

template <typename Scalar>
ShiftJetT<Scalar> shift_jet_scalar(const ShiftFamilyConfig& cfg,
                                   const AxisymmetricBasis& basis,
                                   const std::vector<Scalar>& coefficients,
                                   double t, double x, double y, double z) {
    const double R0 = cfg.scale_R0;
    const double xc = cfg.bubble_center_x0 + cfg.velocity * t;
    const DesignCoordinates p = make_design_coordinates(x - xc, y, z, 0.0, R0);
    const ChainFactors F = chain_factors(y, z, cfg.velocity, R0);

    const DesignDerivativeTable base = baseline_shift_table(
        p.xi, p.radial_squared, cfg.velocity, R0, cfg.bubble_radius,
        cfg.wall_sigma);

    ShiftJetT<Scalar> j;
    j.value = scalar_cast<Scalar>(base.at(0, 0));
    for (int m = 0; m < 4; ++m) {
        const int one[1] = {m};
        j.d1[m] = scalar_cast<Scalar>(chained_derivative(base, F, one, 1));
        for (int n = 0; n < 4; ++n) {
            const int two[2] = {m, n};
            j.d2[m][n] = scalar_cast<Scalar>(chained_derivative(base, F, two, 2));
            for (int l = 0; l < 4; ++l) {
                const int three[3] = {m, n, l};
                j.d3[m][n][l] =
                    scalar_cast<Scalar>(chained_derivative(base, F, three, 3));
            }
        }
    }
    if (coefficients.empty()) return j;

    const TensorBasisEvaluation e = basis.evaluate(p, /*want_third=*/true);
    for (std::size_t k = 0; k < e.size(); ++k) {
        const Scalar& c =
            coefficients[static_cast<std::size_t>(e.flat_index[k])];
        const DesignDerivativeTable Bt = single_mode_table(basis, e, p, k);
        j.value = j.value + c * scalar_cast<Scalar>(Bt.at(0, 0));
        for (int m = 0; m < 4; ++m) {
            const int one[1] = {m};
            j.d1[m] = j.d1[m] +
                      c * scalar_cast<Scalar>(chained_derivative(Bt, F, one, 1));
            for (int n = 0; n < 4; ++n) {
                const int two[2] = {m, n};
                j.d2[m][n] =
                    j.d2[m][n] +
                    c * scalar_cast<Scalar>(chained_derivative(Bt, F, two, 2));
                for (int l = 0; l < 4; ++l) {
                    const int three[3] = {m, n, l};
                    j.d3[m][n][l] =
                        j.d3[m][n][l] +
                        c * scalar_cast<Scalar>(
                                chained_derivative(Bt, F, three, 3));
                }
            }
        }
    }
    return j;
}

// The scalar-generic F1 metric jet. Identical algebra to f1_metric_jet, which
// the gate checks by running both at the same point in double.
template <typename Scalar>
MetricJet<Scalar> f1_metric_jet_scalar(const ShiftFamilyConfig& cfg,
                                       const AxisymmetricBasis& basis,
                                       const std::vector<Scalar>& coefficients,
                                       double t, double x, double y, double z) {
    const ShiftJetT<Scalar> b =
        shift_jet_scalar<Scalar>(cfg, basis, coefficients, t, x, y, z);

    ADMDesignFields<Scalar> adm = ADMDesignFields<Scalar>::minkowski();
    adm.shift_contravariant[0] = b.value;

    MetricJet<Scalar> jet;
    jet.g_cov = assemble_covariant(adm);
    bool ok = false;
    jet.g_contra = assemble_contravariant(adm, ok);
    jet.valid = ok;
    jet.rejection =
        ok ? MetricRejection::None : MetricRejection::SingularSpatialMetric;

    const Scalar two = scalar_cast<Scalar>(2.0);
    for (int m = 0; m < 4; ++m) {
        CovariantMetric<Scalar>& d1 = jet.dg[static_cast<std::size_t>(m)];
        d1.m(0, 0) = two * b.value * b.d1[m];
        d1.m(0, 1) = b.d1[m];
        d1.m(1, 0) = b.d1[m];
        for (int n = 0; n < 4; ++n) {
            CovariantMetric<Scalar>& d2 =
                jet.d2g[static_cast<std::size_t>(m)][static_cast<std::size_t>(n)];
            d2.m(0, 0) = two * (b.d1[m] * b.d1[n] + b.value * b.d2[m][n]);
            d2.m(0, 1) = b.d2[m][n];
            d2.m(1, 0) = b.d2[m][n];
            for (int l = 0; l < 4; ++l) {
                CovariantMetric<Scalar>& d3 =
                    jet.d3g[static_cast<std::size_t>(m)]
                           [static_cast<std::size_t>(n)]
                           [static_cast<std::size_t>(l)];
                d3.m(0, 0) = two * (b.d2[m][n] * b.d1[l] + b.d1[m] * b.d2[n][l] +
                                    b.d1[n] * b.d2[m][l] +
                                    b.value * b.d3[m][n][l]);
                d3.m(0, 1) = b.d3[m][n][l];
                d3.m(1, 0) = b.d3[m][n][l];
            }
        }
    }
    return jet;
}

// --- Source semantics ---------------------------------------------------------
//
// The A5C convention, carried as data so every tensor and objective can name it.
//
//     G + Lambda g = 8 pi T^matter
//     T^matter = (G + Lambda g)/(8 pi)
//     T^Lambda = -Lambda g/(8 pi)
//     T^total  = G/(8 pi)
//
// With Lambda = 0 all three collapse, which is the canonical optimization case
// and the ONLY one an objective history is allowed to mix.
enum class SourceSemantics { MatterOnly, TotalEffective, LambdaOnly };

inline const char* to_string(SourceSemantics s) {
    switch (s) {
        case SourceSemantics::MatterOnly: return "matter_only";
        case SourceSemantics::TotalEffective: return "total_effective";
        case SourceSemantics::LambdaOnly: return "lambda_only";
    }
    return "unknown";
}

inline constexpr double kEightPi = 25.132741228718345908;

template <typename Scalar>
Mat4<Scalar> stress_energy_from_einstein(const Mat4<Scalar>& G,
                                         const Mat4<Scalar>& g, double lambda,
                                         SourceSemantics which) {
    Mat4<Scalar> out;
    const Scalar inv = scalar_cast<Scalar>(1.0 / kEightPi);
    const Scalar lam = scalar_cast<Scalar>(lambda);
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu) {
            if (which == SourceSemantics::MatterOnly)
                out(mu, nu) = (G(mu, nu) + lam * g(mu, nu)) * inv;
            else if (which == SourceSemantics::TotalEffective)
                out(mu, nu) = G(mu, nu) * inv;
            else
                out(mu, nu) = -(lam * g(mu, nu)) * inv;
        }
    return out;
}

// --- The scalar-generic Eulerian density --------------------------------------
//
// rho = T_{mu nu} n^mu n^nu with the ADM normal n^mu = (1/alpha, -beta^i/alpha).
// Written from the ADM fields rather than by inverting g, so it is differentiable
// and needs no eigenvalue test in the inner loop. The gate checks it against
// A6's make_eulerian + project_to_frame at every sample point; agreement there
// is what licenses using this form under an OTI scalar.
//
// The Einstein tensor is built from the same connection jet the double path
// uses, minus the invariants (Kretschmann in particular) that rho does not need.
template <typename Scalar>
Mat4<Scalar> einstein_from_jet(const MetricJet<Scalar>& jet) {
    const ConnectionJet<Scalar> c = connection_jet(jet, /*want_second=*/false);
    const Rank4<Scalar> riem = geometry::riemann_at(c.gamma, c.dgamma);
    const Mat4<Scalar> ric = contract_13(riem);
    const Scalar R = full_contract(jet.g_contra.m, ric);
    return geometry::einstein_at(jet.g_cov.m, ric, R);
}

template <typename Scalar>
Scalar eulerian_density(const ShiftFamilyConfig& cfg,
                        const AxisymmetricBasis& basis,
                        const std::vector<Scalar>& coefficients, double t,
                        double x, double y, double z, double lambda,
                        SourceSemantics which) {
    const MetricJet<Scalar> jet =
        f1_metric_jet_scalar<Scalar>(cfg, basis, coefficients, t, x, y, z);
    const Mat4<Scalar> G = einstein_from_jet(jet);
    const Mat4<Scalar> T =
        stress_energy_from_einstein(G, jet.g_cov.m, lambda, which);

    // The ADM normal, read off the INVERSE metric so no branch and no
    // eigenvalue test enters a differentiated path:
    //
    //     n_mu = (-alpha, 0, 0, 0),  alpha = 1/sqrt(-g^{00})
    //     n^mu = g^{mu nu} n_nu = -g^{mu 0} / sqrt(-g^{00})
    //
    // n^0 = sqrt(-g^{00}) = 1/alpha and n^i = -beta^i/alpha follow, which is
    // A6's eulerian_normal written in terms of what is already available here.
    // Valid wherever g^{00} < 0, i.e. wherever the foliation exists at all;
    // the double path checks that through adm_decompose and rejects otherwise.
    using std::sqrt;
    const Mat4<Scalar>& gi = jet.g_contra.m;
    const Scalar root = sqrt(-gi(0, 0));
    Scalar n[4];
    for (int mu = 0; mu < 4; ++mu) n[mu] = -gi(mu, 0) / root;

    Scalar rho = scalar_cast<Scalar>(0.0);
    for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu) rho = rho + T(mu, nu) * n[mu] * n[nu];
    return rho;
}

// --- Smooth surrogate ---------------------------------------------------------
//
// s_tau(u) = tau log(1 + exp(u/tau)), evaluated in the overflow-safe form
//
//     s_tau(u) = max(u, 0) + tau log(1 + exp(-|u|/tau))
//
// which is exact for both signs and never forms exp of a large positive number.
// softplus OVERESTIMATES max(u, 0) everywhere, by exactly tau log 2 at u = 0,
// so E_neg,tau is biased upward and the bias is measured, not assumed away.
template <typename Scalar>
Scalar softplus(const Scalar& u, double tau) {
    using std::exp;
    using std::log;
    const double ur = real_part(u);
    const Scalar t = scalar_cast<Scalar>(tau);
    const Scalar inv_t = scalar_cast<Scalar>(1.0 / tau);
    if (ur >= 0.0) {
        const Scalar e = exp(-u * inv_t);
        return u + t * log(scalar_cast<Scalar>(1.0) + e);
    }
    const Scalar e = exp(u * inv_t);
    return t * log(scalar_cast<Scalar>(1.0) + e);
}

inline double softplus_bias_at_zero(double tau) {
    return tau * 0.69314718055994530942;   // tau log 2
}

// --- Objectives ---------------------------------------------------------------

struct ObjectiveValues {
    double E_net = 0.0;
    double E_plus = 0.0;
    double E_neg = 0.0;    // POSITIVE magnitude
    double E_abs = 0.0;
    double J_NEC = 0.0;

    double rho_min = 0.0;
    double rho_max = 0.0;
    double N_min_global = 0.0;
    double K_max = 0.0;
    double sqrt_gamma_min = 0.0;
    double sqrt_gamma_max = 0.0;

    double net_closure_residual = 0.0;   // |E_net - (E_plus - E_neg)|
    double abs_closure_residual = 0.0;   // |E_abs - (E_plus + E_neg)|

    std::size_t points = 0;
    std::size_t invalid_points = 0;
    std::size_t nec_failures = 0;
    double volume = 0.0;                 // integral of dV_gamma, for the element check
    double lambda = 0.0;
    SourceSemantics semantics = SourceSemantics::TotalEffective;
};

// One integration set: points in Cartesian coordinates at a fixed time, with the
// dxi ds weight already carrying pi R0^3. sqrt(gamma) is applied here.
struct IntegrationSet {
    std::vector<double> x, y, z;
    std::vector<double> weight;         // pi R0^3 dxi ds
    double time = 0.0;
    std::string route;
    int quadrature_order = 8;
    int subdivision = 1;
    double xi_extent = 0.0;
    double s_extent = 0.0;
    std::size_t elements = 0;
};

// The full objective sweep. want_nec adds the A6 global null minimisation, which
// dominates the cost and is only needed where J_NEC is wanted.
inline ObjectiveValues evaluate_objectives(
    const ShiftFamilyConfig& cfg, const AxisymmetricBasis& basis,
    const std::vector<double>& coefficients, const IntegrationSet& set,
    double lambda, SourceSemantics which, bool want_nec,
    bool want_kretschmann) {
    ObjectiveValues v;
    v.lambda = lambda;
    v.semantics = which;
    v.rho_min = 1e300;
    v.rho_max = -1e300;
    v.N_min_global = 1e300;
    v.sqrt_gamma_min = 1e300;
    v.sqrt_gamma_max = -1e300;
    v.points = set.x.size();

    for (std::size_t p = 0; p < set.x.size(); ++p) {
        const double t = set.time, x = set.x[p], y = set.y[p], z = set.z[p];
        const MetricJet<double> jet =
            f1_metric_jet(cfg, basis, coefficients, t, x, y, z);
        if (!jet.valid) { ++v.invalid_points; continue; }

        const PointGeometry pg =
            want_kretschmann ? point_geometry(jet) : PointGeometry();
        Mat4<double> G;
        if (want_kretschmann) {
            G = pg.einstein;
            v.K_max = std::fmax(v.K_max, std::fabs(pg.kretschmann));
        } else {
            G = einstein_from_jet(jet);
        }
        const Mat4<double> T =
            stress_energy_from_einstein(G, jet.g_cov.m, lambda, which);

        const geometry::AdmDecomposition adm =
            geometry::adm_decompose(jet.g_cov.m, jet.g_contra.m);
        if (!adm.valid) { ++v.invalid_points; continue; }

        // sqrt(gamma): the physical volume element. Exactly 1 for gamma = I,
        // measured rather than assumed.
        const double det_gamma =
            adm.gamma[0] * (adm.gamma[4] * adm.gamma[8] - adm.gamma[5] * adm.gamma[7]) -
            adm.gamma[1] * (adm.gamma[3] * adm.gamma[8] - adm.gamma[5] * adm.gamma[6]) +
            adm.gamma[2] * (adm.gamma[3] * adm.gamma[7] - adm.gamma[4] * adm.gamma[6]);
        const double sg = std::sqrt(std::fmax(det_gamma, 0.0));
        v.sqrt_gamma_min = std::fmin(v.sqrt_gamma_min, sg);
        v.sqrt_gamma_max = std::fmax(v.sqrt_gamma_max, sg);

        const double dV = set.weight[p] * sg;
        v.volume += dV;

        const geometry::Observer obs = geometry::make_eulerian(jet.g_cov.m, adm);
        if (!obs.valid) { ++v.invalid_points; continue; }
        double rho = 0.0;
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu)
                rho += T(mu, nu) * obs.u[mu] * obs.u[nu];

        v.rho_min = std::fmin(v.rho_min, rho);
        v.rho_max = std::fmax(v.rho_max, rho);
        v.E_net += rho * dV;
        v.E_plus += std::fmax(rho, 0.0) * dV;
        v.E_neg += std::fmax(-rho, 0.0) * dV;
        v.E_abs += std::fabs(rho) * dV;

        if (!want_nec) continue;
        const geometry::Tetrad tet = geometry::build_tetrad(jet.g_cov.m, adm);
        if (!tet.valid) { ++v.nec_failures; continue; }
        const geometry::FrameStress fs = geometry::project_to_frame(tet, T);
        const geometry::NecResult nec = geometry::nec_minimum(fs, 256);
        if (!nec.ok) { ++v.nec_failures; continue; }
        v.N_min_global = std::fmin(v.N_min_global, nec.n_min);
        v.J_NEC += std::fmax(-nec.n_min, 0.0) * dV;
    }

    if (v.rho_min > 1e299) v.rho_min = 0.0;
    if (v.rho_max < -1e299) v.rho_max = 0.0;
    if (v.N_min_global > 1e299) v.N_min_global = 0.0;
    if (v.sqrt_gamma_min > 1e299) v.sqrt_gamma_min = 0.0;
    if (v.sqrt_gamma_max < -1e299) v.sqrt_gamma_max = 0.0;
    v.net_closure_residual = std::fabs(v.E_net - (v.E_plus - v.E_neg));
    v.abs_closure_residual = std::fabs(v.E_abs - (v.E_plus + v.E_neg));
    return v;
}

// The smooth surrogate, scalar-generic so it can be differentiated. This is the
// candidate optimization objective; it integrates s_tau(-rho), NOT max(-rho, 0).
template <typename Scalar>
Scalar smooth_negative_energy(const ShiftFamilyConfig& cfg,
                              const AxisymmetricBasis& basis,
                              const std::vector<Scalar>& coefficients,
                              const IntegrationSet& set, double tau,
                              double lambda, SourceSemantics which) {
    Scalar acc = scalar_cast<Scalar>(0.0);
    for (std::size_t p = 0; p < set.x.size(); ++p) {
        const Scalar rho = eulerian_density<Scalar>(
            cfg, basis, coefficients, set.time, set.x[p], set.y[p], set.z[p],
            lambda, which);
        // gamma = I for the shift-only families, so sqrt(gamma) = 1 exactly and
        // carries no coefficient dependence. The double path measures that;
        // here it would only add a constant factor of one per point.
        acc = acc + softplus(-rho, tau) * scalar_cast<Scalar>(set.weight[p]);
    }
    return acc;
}

// The signed net energy, scalar-generic. Smooth everywhere (no max), so its
// gradient is exact rather than a surrogate.
template <typename Scalar>
Scalar net_energy(const ShiftFamilyConfig& cfg, const AxisymmetricBasis& basis,
                  const std::vector<Scalar>& coefficients,
                  const IntegrationSet& set, double lambda,
                  SourceSemantics which) {
    Scalar acc = scalar_cast<Scalar>(0.0);
    for (std::size_t p = 0; p < set.x.size(); ++p) {
        const Scalar rho = eulerian_density<Scalar>(
            cfg, basis, coefficients, set.time, set.x[p], set.y[p], set.z[p],
            lambda, which);
        acc = acc + rho * scalar_cast<Scalar>(set.weight[p]);
    }
    return acc;
}

// --- Integration sets ---------------------------------------------------------
//
// Gauss-Legendre nodes at arbitrary order, by Newton iteration on P_n with the
// Chebyshev initial guess. Written here rather than tabulated because the
// convergence study needs orders 4, 6, 8, 10 and 12, and five hand-typed tables
// is five chances to lose a digit. The order-8 case is checked against O0A's
// existing gauss8_nodes() table by the gate.
inline void gauss_legendre(int n, std::vector<double>& node,
                           std::vector<double>& weight) {
    node.assign(static_cast<std::size_t>(n), 0.0);
    weight.assign(static_cast<std::size_t>(n), 0.0);
    const double kPiL = 3.14159265358979323846;
    for (int i = 0; i < n; ++i) {
        double x = std::cos(kPiL * (static_cast<double>(i) + 0.75) /
                            (static_cast<double>(n) + 0.5));
        double dp = 0.0;
        for (int it = 0; it < 100; ++it) {
            double p0 = 1.0, p1 = 0.0;
            for (int k = 0; k < n; ++k) {
                const double p2 = p1;
                p1 = p0;
                p0 = ((2.0 * k + 1.0) * x * p1 - static_cast<double>(k) * p2) /
                     (static_cast<double>(k) + 1.0);
            }
            dp = static_cast<double>(n) * (x * p0 - p1) / (x * x - 1.0);
            const double dx = -p0 / dp;
            x += dx;
            if (std::fabs(dx) < 1e-16) break;
        }
        node[static_cast<std::size_t>(i)] = x;
        weight[static_cast<std::size_t>(i)] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
}

// The elementwise analytic-jet route, at a chosen order, subdivision and extent.
// Elements are the knot spans of the basis, each optionally split `subdivision`
// ways in each direction.
inline IntegrationSet build_elementwise_set(const AxisymmetricBasis& basis,
                                            int order, int subdivision,
                                            double xi_extent, double s_extent,
                                            double time, double velocity) {
    IntegrationSet set;
    set.time = time;
    set.route = "elementwise_analytic_jet";
    set.quadrature_order = order;
    set.subdivision = subdivision;
    set.xi_extent = xi_extent;
    set.s_extent = s_extent;
    const double R0 = basis.scale_R0();
    const double kPiL = 3.14159265358979323846;
    const double vol = kPiL * R0 * R0 * R0;

    std::vector<double> gn, gw;
    gauss_legendre(order, gn, gw);

    // Element boundaries: the basis knot spans, clipped to the requested extent,
    // then subdivided. Beyond the basis range the integrand is still the
    // BASELINE, which is what the domain-enlargement study needs, so the tail is
    // tiled uniformly at the same element size.
    auto edges = [&](const QuinticBSplineBasis1D& b, double lo, double hi) {
        std::vector<double> e;
        e.push_back(lo);
        for (double t : b.distinct_interior_knots())
            if (t > lo + 1e-12 && t < hi - 1e-12) e.push_back(t);
        if (b.upper() < hi - 1e-12) {
            const double step = (hi - b.upper()) /
                                std::fmax(1.0, std::floor((hi - b.upper()) /
                                                          std::fmax(1e-9, 0.5)));
            for (double t = b.upper(); t < hi - 1e-12; t += step)
                if (t > lo + 1e-12) e.push_back(t);
        }
        e.push_back(hi);
        std::sort(e.begin(), e.end());
        e.erase(std::unique(e.begin(), e.end(),
                            [](double a, double c) {
                                return std::fabs(a - c) < 1e-12;
                            }),
                e.end());
        std::vector<double> out;
        for (std::size_t i = 0; i + 1 < e.size(); ++i)
            for (int k = 0; k < subdivision; ++k)
                out.push_back(e[i] + (e[i + 1] - e[i]) * k / subdivision);
        out.push_back(e.back());
        return out;
    };

    const std::vector<double> ex = edges(basis.axial(), -xi_extent, xi_extent);
    const std::vector<double> es = edges(basis.radial(), 0.0, s_extent);
    set.elements = (ex.size() - 1) * (es.size() - 1);

    for (std::size_t a = 0; a + 1 < ex.size(); ++a)
        for (std::size_t b = 0; b + 1 < es.size(); ++b) {
            const double xa = ex[a], xb = ex[a + 1];
            const double sa = es[b], sb = es[b + 1];
            if (!(xb > xa) || !(sb > sa)) continue;
            for (int i = 0; i < order; ++i)
                for (int j = 0; j < order; ++j) {
                    const double xi = 0.5 * (xa + xb) +
                                      0.5 * (xb - xa) * gn[static_cast<std::size_t>(i)];
                    const double s = 0.5 * (sa + sb) +
                                     0.5 * (sb - sa) * gn[static_cast<std::size_t>(j)];
                    const double w = 0.5 * (xb - xa) *
                                     gw[static_cast<std::size_t>(i)] * 0.5 *
                                     (sb - sa) * gw[static_cast<std::size_t>(j)] *
                                     vol;
                    // The design chart is comoving; the Cartesian point at time
                    // `time` follows the bubble centre.
                    set.x.push_back(xi * R0 + velocity * time);
                    set.y.push_back(R0 * std::sqrt(std::fmax(s, 0.0)));
                    set.z.push_back(0.0);
                    set.weight.push_back(w);
                }
        }
    return set;
}

// The independent reference: a UNIFORM grid in (x, rho) integrated with a
// composite rule, not with Gauss-Legendre and not on the basis elements.
//
// It evaluates the SAME analytic pointwise quantities, so it cross-checks the
// integration and nothing else -- in particular it does not reintroduce a
// finite-difference metric derivative, which is the trap the A6 grid route falls
// into and the reason that route is a legacy comparison rather than an authority.
//
// The measure: d^3x = 2 pi rho drho dx, so dV = 2 pi rho drho dx sqrt(gamma).
// Note this is a DIFFERENT parameterisation of the same measure from the
// elementwise route's pi R0^3 dxi ds; agreement between them is therefore also a
// check that the two forms of the volume element are the same measure.
inline IntegrationSet build_uniform_set(double x_half_extent,
                                        double rho_extent, int nx, int nrho,
                                        double time, double velocity) {
    IntegrationSet set;
    set.time = time;
    set.route = "uniform_axisymmetric_composite";
    set.quadrature_order = 0;
    set.subdivision = 1;
    set.xi_extent = x_half_extent;
    set.s_extent = rho_extent;
    set.elements = static_cast<std::size_t>(nx) * static_cast<std::size_t>(nrho);
    const double kPiL = 3.14159265358979323846;
    const double hx = 2.0 * x_half_extent / nx;
    const double hr = rho_extent / nrho;

    // Composite Simpson in each direction, on an even number of intervals.
    auto simpson_weight = [](int i, int n, double h) {
        if (i == 0 || i == n) return h / 3.0;
        return (i % 2 == 1 ? 4.0 : 2.0) * h / 3.0;
    };
    for (int i = 0; i <= nx; ++i)
        for (int j = 0; j <= nrho; ++j) {
            const double x = -x_half_extent + i * hx;
            const double r = j * hr;
            const double w = simpson_weight(i, nx, hx) *
                             simpson_weight(j, nrho, hr) * 2.0 * kPiL * r;
            if (w == 0.0) continue;
            set.x.push_back(x + velocity * time);
            set.y.push_back(r);
            set.z.push_back(0.0);
            set.weight.push_back(w);
        }
    return set;
}

// --- Objective classification -------------------------------------------------

struct ObjectiveClassification {
    std::string name;
    std::string family;
    std::string observer;
    std::string source_semantics;
    std::string integration_route;
    double smoothing_tau = 0.0;

    bool implemented = false;
    bool pointwise_verified = false;
    bool behavior_preserving = false;
    bool quadrature_converged = false;
    bool domain_converged = false;
    bool independent_reference_agreement = false;
    bool gradient_verified = false;
    bool resolution_certified = false;
    bool exact_gradient_available = false;
    bool smooth_surrogate_gradient_available = false;
    bool nonsmooth_or_minimizer_degenerate = false;
    bool ready_for_optimization = false;
    std::string blocking_reason;
};

}  // namespace design
}  // namespace warp_hypad
