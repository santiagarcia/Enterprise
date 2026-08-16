#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "warp_hypad/design/adm_family.hpp"

// ---------------------------------------------------------------------------
// O2: positive-energy feasibility.
//
// THE DESIGN CRITERION CHANGED, AND SO DID WHAT COUNTS AS SUCCESS.
//
// Through O1 the objective was to make the integrated negative energy SMALLER.
// From O2 the objective is to make it ZERO -- pointwise:
//
//     rho(x) >= 0   everywhere in the certified domain.
//
// E_neg is now a CONSTRAINT-VIOLATION MEASURE, not a score. A design has not
// succeeded because E_neg fell; it succeeds only when both
//
//     rho_min >= -epsilon_rho     and     E_neg <= epsilon_E
//
// with both tolerances derived from convergence and independent-reference
// studies rather than chosen once the answer is visible.
//
// WHY THE INTEGRAL ALONE CANNOT CERTIFY. A quadrature sum can be arbitrarily
// small while a narrow pocket still has rho < 0 between the nodes: shrink the
// support and E_neg falls like the volume even as the pocket DEEPENS. So every
// certification here carries rho_min, where it is, what fraction of the domain
// is negative, and a local refinement around each detected pocket. A candidate
// whose sign is resolution-dependent anywhere is not positive-energy, whatever
// the integral says.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

// The Eulerian energy density of an arbitrary ADM family, scalar-generic.
// Same chain as O0C -- metric jet, connection, Riemann, Einstein, observer --
// with the family spec deciding which fields are live.
template <typename Scalar>
Scalar family_eulerian_density(const AdmFamilySpec& spec,
                               const ShiftFamilyConfig& cfg,
                               const AxisymmetricBasis& basis,
                               const AdmFieldCoefficients<Scalar>& c,
                               double t, double x, double y, double z,
                               double lambda, SourceSemantics which) {
    const MetricJet<Scalar> jet =
        adm_metric_jet<Scalar>(spec, cfg, basis, c, t, x, y, z);
    const Mat4<Scalar> G = einstein_from_jet(jet);
    const Mat4<Scalar> T =
        stress_energy_from_einstein(G, jet.g_cov.m, lambda, which);
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

// THE FEASIBILITY RESIDUAL.
//
//     R_rho_tau = integral [ softplus_tau(-rho) ]^2 dV
//
// SQUARED, and that is the whole reason this is not just O1's objective again.
// softplus_tau(u) tends to tau log 2 rather than to zero as u -> 0 from below,
// so the linear form carries the vacuum floor O0D found and, worse, keeps
// pulling on regions that are ALREADY feasible. Squaring makes the floor
// quadratically small there while leaving deep violations dominant, so the
// gradient points at the pockets that actually matter.
//
// It is a RESIDUAL. Its reduction is not a physical success; rho_min and E_neg
// are the certification quantities.
template <typename Scalar>
Scalar feasibility_residual(const AdmFamilySpec& spec,
                            const ShiftFamilyConfig& cfg,
                            const AxisymmetricBasis& basis,
                            const AdmFieldCoefficients<Scalar>& c,
                            const IntegrationSet& set, double tau,
                            double lambda, SourceSemantics which) {
    Scalar acc = scalar_cast<Scalar>(0.0);
    for (std::size_t p = 0; p < set.x.size(); ++p) {
        const Scalar rho = family_eulerian_density<Scalar>(
            spec, cfg, basis, c, set.time, set.x[p], set.y[p], set.z[p],
            lambda, which);
        const Scalar s = softplus(-rho, tau);
        acc = acc + s * s * scalar_cast<Scalar>(set.weight[p]);
    }
    return acc;
}

// Everything a feasibility verdict needs, and nothing that could stand in for
// it. Note that rho_min carries its LOCATION: a violation you cannot point at
// is a violation you cannot refine around.
struct FeasibilityReport {
    double rho_min = 0.0;
    double rho_max = 0.0;
    double rho_min_x = 0.0, rho_min_y = 0.0, rho_min_z = 0.0;
    double E_neg = 0.0, E_plus = 0.0, E_net = 0.0, E_abs = 0.0;
    double negative_volume = 0.0;
    double total_volume = 0.0;
    double negative_volume_fraction = 0.0;
    std::size_t negative_points = 0;
    std::size_t points = 0;
    std::size_t invalid_points = 0;
    double residual = 0.0;

    // Filled by the local refinement pass.
    bool refined = false;
    double refined_rho_min = 0.0;
    bool sign_stable_under_refinement = false;

    bool feasible = false;
    std::string verdict;
};

inline FeasibilityReport evaluate_feasibility(
    const AdmFamilySpec& spec, const ShiftFamilyConfig& cfg,
    const AxisymmetricBasis& basis,
    const AdmFieldCoefficients<double>& c,
    const IntegrationSet& set, double tau, double lambda,
    SourceSemantics which) {
    FeasibilityReport r;
    r.rho_min = 1e300;
    r.rho_max = -1e300;
    r.points = set.x.size();
    for (std::size_t p = 0; p < set.x.size(); ++p) {
        const MetricJet<double> jet = adm_metric_jet<double>(
            spec, cfg, basis, c, set.time, set.x[p], set.y[p], set.z[p]);
        if (!jet.valid) { ++r.invalid_points; continue; }
        const geometry::AdmDecomposition adm =
            geometry::adm_decompose(jet.g_cov.m, jet.g_contra.m);
        if (!adm.valid) { ++r.invalid_points; continue; }
        const double rho = family_eulerian_density<double>(
            spec, cfg, basis, c, set.time, set.x[p], set.y[p], set.z[p],
            lambda, which);
        if (!std::isfinite(rho)) { ++r.invalid_points; continue; }

        const double det =
            adm.gamma[0] * (adm.gamma[4] * adm.gamma[8] - adm.gamma[5] * adm.gamma[7]) -
            adm.gamma[1] * (adm.gamma[3] * adm.gamma[8] - adm.gamma[5] * adm.gamma[6]) +
            adm.gamma[2] * (adm.gamma[3] * adm.gamma[7] - adm.gamma[4] * adm.gamma[6]);
        const double dV = set.weight[p] * std::sqrt(std::fmax(det, 0.0));
        r.total_volume += dV;
        if (rho < r.rho_min) {
            r.rho_min = rho;
            r.rho_min_x = set.x[p]; r.rho_min_y = set.y[p];
            r.rho_min_z = set.z[p];
        }
        r.rho_max = std::fmax(r.rho_max, rho);
        r.E_net += rho * dV;
        r.E_plus += std::fmax(rho, 0.0) * dV;
        r.E_neg += std::fmax(-rho, 0.0) * dV;
        r.E_abs += std::fabs(rho) * dV;
        const double s = softplus(-rho, tau);
        r.residual += s * s * set.weight[p];
        if (rho < 0.0) { ++r.negative_points; r.negative_volume += dV; }
    }
    if (r.rho_min > 1e299) r.rho_min = 0.0;
    if (r.rho_max < -1e299) r.rho_max = 0.0;
    r.negative_volume_fraction =
        r.total_volume > 0.0 ? r.negative_volume / r.total_volume : 0.0;
    return r;
}

// LOCAL REFINEMENT AROUND THE WORST POINT.
//
// A tensor-product Gauss rule can step over a narrow pocket entirely, so the
// certification re-samples a small box around the recorded minimum on a dense
// uniform lattice. If the refined minimum is materially lower than the coarse
// one, the coarse grid was not resolving the violation and the candidate's sign
// is resolution-dependent -- which disqualifies it regardless of E_neg.
inline void refine_around_minimum(FeasibilityReport& r,
                                  const AdmFamilySpec& spec,
                                  const ShiftFamilyConfig& cfg,
                                  const AxisymmetricBasis& basis,
                                  const AdmFieldCoefficients<double>& c,
                                  double time, double half_width, int n,
                                  double lambda, SourceSemantics which,
                                  double epsilon_rho) {
    double worst = r.rho_min;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j) {
            const double x = r.rho_min_x - half_width +
                             2.0 * half_width * i / n;
            const double y = std::fmax(0.0, r.rho_min_y - half_width +
                                                2.0 * half_width * j / n);
            const double rho = family_eulerian_density<double>(
                spec, cfg, basis, c, time, x, y, 0.0, lambda, which);
            if (std::isfinite(rho)) worst = std::fmin(worst, rho);
        }
    r.refined = true;
    r.refined_rho_min = worst;
    // The sign is stable when refinement does not push a nonnegative field
    // negative, and does not deepen an existing violation beyond tolerance.
    r.sign_stable_under_refinement =
        std::fabs(worst - r.rho_min) <= std::fmax(epsilon_rho,
                                                  1e-3 * std::fabs(r.rho_min));
}

// The verdict, as a conjunction. Written this way so no single number can carry
// it alone.
inline void classify_feasibility(FeasibilityReport& r, double epsilon_rho,
                                 double epsilon_E) {
    const double worst =
        r.refined ? std::fmin(r.rho_min, r.refined_rho_min) : r.rho_min;
    const bool pointwise = worst >= -epsilon_rho;
    const bool integrated = r.E_neg <= epsilon_E;
    const bool stable = !r.refined || r.sign_stable_under_refinement;
    r.feasible = pointwise && integrated && stable && r.invalid_points == 0;
    if (r.feasible) {
        r.verdict = "no negative energy detected within certified tolerance";
    } else {
        r.verdict = "negative energy present:";
        if (!pointwise)
            r.verdict += " rho_min " + std::to_string(worst) + " below -" +
                         std::to_string(epsilon_rho) + ";";
        if (!integrated)
            r.verdict += " E_neg " + std::to_string(r.E_neg) + " above " +
                         std::to_string(epsilon_E) + ";";
        if (!stable) r.verdict += " sign is resolution-dependent;";
        if (r.invalid_points) r.verdict += " invalid points present;";
    }
}

}  // namespace design
}  // namespace warp_hypad
