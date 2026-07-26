#pragma once

#include "program_b/inverse_design/passivity.hpp"

// ---------------------------------------------------------------------------
// Adaptive frequency validation.
//
// WHAT THIS IS FOR, AND WHAT IT IS NOT. B3's gate samples passivity on a fixed
// frequency grid. That establishes passivity AT THE SAMPLED POINTS, and nothing
// between them: L(w) is a rational function whose smallest eigenvalue can dip
// negative inside a narrow interval -- near a resonance, near a zero of det mu,
// or where two eigenvalue branches cross -- and a coarse grid walks straight
// past it.
//
// This utility refines the grid where a violation is most likely to hide:
//
//   * local minima of lambda_min(L)
//   * intervals of rapid eigenvalue change
//   * resonances of the model (poles of eps and mu)
//   * zeros and near-zeros of det mu, where mu^-1 blows up
//   * the immediate neighbourhood of any negative sample found
//
// It NARROWS the unverified region. It does not eliminate it, and it is NOT a
// replacement for B3's sampled gate -- it belongs to B4 validation. A genuine
// global statement would need an analytic positive-real argument, which no
// milestone here has carried out. Both results are reported, coarse and refined,
// so the difference between them is visible rather than absorbed.
// ---------------------------------------------------------------------------

namespace program_b {
namespace inverse_design {

struct AdaptiveScanResult {
    double coarse_min_eigenvalue = 0.0;
    double coarse_worst_frequency = 0.0;
    int coarse_samples = 0;

    double refined_min_eigenvalue = 0.0;
    double refined_worst_frequency = 0.0;
    int refined_samples = 0;

    // How much worse the refined scan found things. Positive means the coarse
    // grid was optimistic.
    double improvement = 0.0;
    bool coarse_says_passive = false;
    bool refined_says_passive = false;

    // Structure the refinement targeted.
    std::vector<double> resonances;
    std::vector<double> mu_near_zeros;
    int refinement_rounds = 0;
};

// Frequencies where the model has structure worth resolving: the Lorentz
// resonances themselves, and any frequency where |det mu| dips.
template <typename Scalar>
void collect_structure(const PhysicalDesignVariables<Scalar>& v, double w_min,
                       double w_max, AdaptiveScanResult& out) {
    Scalar p[kPhysicalSlots];
    physical_values(v, p);

    auto add = [&](double w) {
        if (w > w_min && w < w_max) out.resonances.push_back(w);
    };
    if (v.spec.eps_poles > 0) add(real_part(p[kSlotEpsPole0W0]));
    if (v.spec.eps_poles > 1) add(real_part(p[kSlotEpsPole1W0]));
    if (v.spec.mu_poles > 0) add(real_part(p[kSlotMuPole0W0]));
    if (v.spec.magnetoelectric) add(real_part(p[kSlotZetaW0]));

    // Zeros of mu: sweep and bracket sign changes of Re mu, which is where
    // mu^-1 diverges and the loss matrix is most delicate.
    const int n = 2000;
    double prev_w = w_min;
    double prev = real_part(evaluate_response(v, w_min).mu.re);
    for (int i = 1; i <= n; ++i) {
        const double w = w_min + (w_max - w_min) * static_cast<double>(i) /
                                     static_cast<double>(n);
        const double cur = real_part(evaluate_response(v, w).mu.re);
        if ((prev < 0.0) != (cur < 0.0)) out.mu_near_zeros.push_back(0.5 * (prev_w + w));
        prev_w = w;
        prev = cur;
    }
}

// The adaptive search itself.
template <typename Scalar>
AdaptiveScanResult adaptive_passivity_scan(const PhysicalDesignVariables<Scalar>& v,
                                          double w_min, double w_max,
                                          int coarse_samples = 21,
                                          int rounds = 6,
                                          int per_interval = 9) {
    AdaptiveScanResult out;
    out.coarse_samples = coarse_samples;

    auto lam = [&](double w) { return passivity_margin_at(v, w); };

    // --- coarse pass --------------------------------------------------------
    std::vector<double> grid;
    for (int i = 0; i < coarse_samples; ++i)
        grid.push_back(w_min + (w_max - w_min) * static_cast<double>(i) /
                                   static_cast<double>(coarse_samples - 1));
    out.coarse_min_eigenvalue = 1e300;
    std::vector<double> value;
    for (double w : grid) {
        const double l = lam(w);
        value.push_back(l);
        if (l < out.coarse_min_eigenvalue) {
            out.coarse_min_eigenvalue = l;
            out.coarse_worst_frequency = w;
        }
    }
    out.coarse_says_passive = out.coarse_min_eigenvalue >= -1e-12;

    out.refined_min_eigenvalue = out.coarse_min_eigenvalue;
    out.refined_worst_frequency = out.coarse_worst_frequency;
    out.refined_samples = coarse_samples;

    collect_structure(v, w_min, w_max, out);

    // --- seed intervals: structure first, then the coarse minimum -----------
    std::vector<std::pair<double, double>> intervals;
    auto bracket = [&](double centre, double half) {
        const double a = std::fmax(w_min, centre - half);
        const double b = std::fmin(w_max, centre + half);
        if (b > a) intervals.push_back({a, b});
    };
    const double span = (w_max - w_min) / static_cast<double>(coarse_samples - 1);
    for (double r : out.resonances) bracket(r, 2.0 * span);
    for (double z : out.mu_near_zeros) bracket(z, 2.0 * span);
    bracket(out.coarse_worst_frequency, span);

    // Intervals of rapid change, and every local minimum of the coarse samples.
    for (std::size_t i = 1; i + 1 < value.size(); ++i) {
        const bool local_min = value[i] < value[i - 1] && value[i] < value[i + 1];
        const double d1 = std::fabs(value[i] - value[i - 1]);
        const double d2 = std::fabs(value[i + 1] - value[i]);
        const double scale = std::fmax(std::fabs(value[i]), 1e-12);
        const bool rapid = std::fmax(d1, d2) > 0.5 * scale;
        if (local_min || rapid) bracket(grid[i], span);
    }

    // --- refine -------------------------------------------------------------
    for (int round = 0; round < rounds && !intervals.empty(); ++round) {
        ++out.refinement_rounds;
        std::vector<std::pair<double, double>> next;
        for (const std::pair<double, double>& iv : intervals) {
            double best = 1e300, best_w = iv.first;
            for (int i = 0; i <= per_interval; ++i) {
                const double w = iv.first + (iv.second - iv.first) *
                                                static_cast<double>(i) /
                                                static_cast<double>(per_interval);
                const double l = lam(w);
                ++out.refined_samples;
                if (l < best) { best = l; best_w = w; }
                if (l < out.refined_min_eigenvalue) {
                    out.refined_min_eigenvalue = l;
                    out.refined_worst_frequency = w;
                }
            }
            // Zoom on the local minimum for the next round.
            const double h = (iv.second - iv.first) / static_cast<double>(per_interval);
            const double a = std::fmax(w_min, best_w - h);
            const double b = std::fmin(w_max, best_w + h);
            if (b > a) next.push_back({a, b});
        }
        intervals.swap(next);
    }

    out.refined_says_passive = out.refined_min_eigenvalue >= -1e-12;
    out.improvement = out.coarse_min_eigenvalue - out.refined_min_eigenvalue;
    return out;
}

}  // namespace inverse_design
}  // namespace program_b
