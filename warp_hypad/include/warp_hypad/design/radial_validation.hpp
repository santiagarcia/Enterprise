#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "warp_hypad/design/constraint_projection.hpp"
#include "warp_hypad/design/quintic_bspline.hpp"

// ---------------------------------------------------------------------------
// F0: the one-dimensional radial validation family, and the projection of a
// known shape function into its constrained space.
//
// WHAT F0 IS FOR. Everything else in O0A is a statement about the basis. F0 is
// the one statement about whether the basis can REPRESENT something: it takes
// the Alcubierre shape function -- a function the project has verified
// independently since A3 -- and asks how well the constrained spline space
// reproduces f, f' and f''. Two derivatives, because the metric built from this
// space is differentiated twice to get curvature, and a projection that matched
// f beautifully while getting f'' wrong would produce a plausible metric with
// wrong curvature.
//
// THE RADIAL COORDINATE IS s_r = (r/R0)^2, for the same structural reason as the
// axisymmetric families: a spline in s_r is automatically even in r with zero
// slope at the origin, so f(0) = 1 with f'(0) = 0 is representable exactly
// rather than approximately.
//
// THE FLAT COLLAR IS NOT FREE. Every basis function whose support reaches beyond
// s_flat_outer is FIXED at zero, so the projection cannot reproduce f there at
// all -- it is required to be exactly zero, and f is not exactly zero. That
// truncation is a deliberate, reported error, not a fit residual: section 12 of
// the specification asks for the collar error separately for precisely this
// reason, and the wall steepness is chosen so the discarded tail is small enough
// to quote honestly rather than hidden by a favourable norm.
// ---------------------------------------------------------------------------

namespace warp_hypad {
namespace design {

// The four zones the projection error is reported in. A single global norm would
// average the wall -- where the physics is -- against the exterior, where
// everything is nearly zero and any method looks good.
enum class RadialZone { ProtectedCabin, Wall, ExteriorTransition, FlatCollar };

inline const char* to_string(RadialZone z) {
    switch (z) {
        case RadialZone::ProtectedCabin: return "protected_cabin";
        case RadialZone::Wall: return "wall";
        case RadialZone::ExteriorTransition: return "exterior_transition";
        case RadialZone::FlatCollar: return "flat_collar";
    }
    return "unknown";
}

struct RadialValidationConfig {
    int n_basis = 12;       // the mandated F0 raw space
    double s_max = 4.0;     // r up to 2 R0
    double scale_R0 = 2.0;  // metres
    double bubble_radius = 2.0;
    double wall_sigma = 3.0;
    double outer_flat_collar = 0.6;   // in s
    double cabin_radius = 0.25;       // in s
    double wall_s_lower = 0.5;        // in s
    double wall_s_upper = 1.75;       // in s

    double s_flat_outer() const { return s_max - outer_flat_collar; }
};

inline RadialZone zone_of(const RadialValidationConfig& cfg, double s) {
    if (s > cfg.s_flat_outer()) return RadialZone::FlatCollar;
    if (s < cfg.cabin_radius) return RadialZone::ProtectedCabin;
    if (s >= cfg.wall_s_lower && s <= cfg.wall_s_upper) return RadialZone::Wall;
    return RadialZone::ExteriorTransition;
}

// The reference shape function and its first two radial derivatives.
//
// f and f' come from production (metrics::shape_function and
// metrics::shape_function_radial_derivative). f'' does NOT exist in production,
// so it is written here in closed form from the same tanh identity and then
// cross-checked against a centred difference of production's f' -- the reference
// is independent of the spline being tested, which is the only property that
// matters for a projection check.
struct ShapeReference {
    double value = 0.0;
    double first = 0.0;
    double second = 0.0;
};

inline ShapeReference shape_reference(double r, double R, double sigma) {
    ShapeReference out;
    const double a = sigma * (r + R), b = sigma * (r - R);
    const double ta = std::tanh(a), tb = std::tanh(b);
    const double den = 2.0 * std::tanh(sigma * R);
    const double sa = 1.0 - ta * ta;  // sech^2(a)
    const double sb = 1.0 - tb * tb;
    out.value = (ta - tb) / den;
    out.first = sigma * (sa - sb) / den;
    // d/dr sech^2(u) = -2 sech^2(u) tanh(u) du/dr
    out.second = sigma * sigma * (-2.0 * sa * ta + 2.0 * sb * tb) / den;
    return out;
}

struct RadialFamily {
    QuinticBSplineBasis1D basis;
    std::vector<ExclusionReason> reason;
    std::vector<int> active;   // free coefficients
    std::vector<int> fixed;    // held at zero by the collar
    RadialValidationConfig config;

    RadialFamily(QuinticBSplineBasis1D b, RadialValidationConfig c)
        : basis(std::move(b)), config(c) {}
};

// A clamped basis on [0, s_max]: clamped at BOTH ends, because unlike a
// correction family F0 must represent a function that is nonzero at the origin
// (f(0) = 1) and it does not have to taper smoothly into the collar -- the
// collar coefficients are fixed at zero and the resulting truncation is reported.
inline RadialFamily make_f0_family(const RadialValidationConfig& cfg) {
    QuinticBSplineBasis1D b(0.0, cfg.s_max, cfg.n_basis);
    RadialFamily f(std::move(b), cfg);
    f.reason.assign(static_cast<std::size_t>(f.basis.size()),
                    ExclusionReason::Active);
    for (int i = 0; i < f.basis.size(); ++i) {
        if (f.basis.support_upper(i) > cfg.s_flat_outer() + 1e-12) {
            f.reason[static_cast<std::size_t>(i)] =
                ExclusionReason::OuterFlatCollar;
            f.fixed.push_back(i);
        } else {
            f.active.push_back(i);
        }
    }
    return f;
}

// The field and its radial derivatives from a coefficient vector over the FULL
// raw space (collar coefficients included and required to be zero).
//
// The chain rule from s_r to r, written once:
//     dq/dr   = q_s (2r/R0^2)
//     d2q/dr2 = q_ss (2r/R0^2)^2 + q_s (2/R0^2)
struct RadialEvaluation {
    double value = 0.0;
    double d_r = 0.0;
    double d_rr = 0.0;
    bool inside = false;
};

inline RadialEvaluation evaluate_radial(const RadialFamily& f,
                                        const std::vector<double>& c,
                                        double s) {
    RadialEvaluation out;
    const BasisEvaluation1D e = f.basis.evaluate(s);
    if (e.empty()) return out;
    out.inside = true;
    double q = 0.0, qs = 0.0, qss = 0.0;
    for (std::size_t k = 0; k < e.size(); ++k) {
        const double ck = c[static_cast<std::size_t>(e.active_indices[k])];
        q += ck * e.value[k];
        qs += ck * e.first[k];
        qss += ck * e.second[k];
    }
    const double R0 = f.config.scale_R0;
    const double r = R0 * std::sqrt(std::fmax(s, 0.0));
    const double dsdr = 2.0 * r / (R0 * R0);
    out.value = q;
    out.d_r = qs * dsdr;
    out.d_rr = qss * dsdr * dsdr + qs * 2.0 / (R0 * R0);
    return out;
}

struct ZoneError {
    double value = 0.0;
    double first = 0.0;
    double second = 0.0;
    double reference_scale = 0.0;
    std::size_t samples = 0;
};

struct RadialProjection {
    std::vector<double> coefficients;  // full raw space; collar entries zero
    int n_raw = 0;
    int n_fixed = 0;
    int n_free = 0;
    int dropped_normal_modes = 0;
    double collar_max_abs = 0.0;       // must be exactly zero
    ZoneError zone[4];
    double condition_estimate = 0.0;
};

// Least-squares projection of the reference shape function onto the family, with
// the collar coefficients held at zero.
//
// Solved as a normal-equation system over the FREE coefficients only. That is
// the elimination of section 7 applied to the simplest possible constraint set:
// the fixed coefficients never enter the unknown vector, so they cannot acquire a
// nonzero value through rounding, and the collar is exactly zero by construction
// rather than by a residual that happens to be small.
inline RadialProjection project_shape_function(const RadialFamily& f,
                                               int samples_per_span = 12) {
    RadialProjection pr;
    const RadialValidationConfig& cfg = f.config;
    pr.n_raw = f.basis.size();
    pr.n_fixed = static_cast<int>(f.fixed.size());
    pr.n_free = static_cast<int>(f.active.size());

    // Sample points: dense and uniform in s over the non-collar region, which
    // is where the fit is allowed to do anything at all.
    std::vector<double> s_samples;
    const int n = samples_per_span * 20;
    for (int k = 0; k <= n; ++k) {
        const double s = cfg.s_flat_outer() * static_cast<double>(k) /
                         static_cast<double>(n);
        s_samples.push_back(s);
    }

    // Normal equations in the free coefficients.
    std::vector<int> slot(static_cast<std::size_t>(pr.n_raw), -1);
    for (std::size_t k = 0; k < f.active.size(); ++k)
        slot[static_cast<std::size_t>(f.active[k])] = static_cast<int>(k);

    DenseMatrix A(pr.n_free, pr.n_free);
    std::vector<double> rhs(static_cast<std::size_t>(pr.n_free), 0.0);
    for (double s : s_samples) {
        const BasisEvaluation1D e = f.basis.evaluate(s);
        if (e.empty()) continue;
        const double r = cfg.scale_R0 * std::sqrt(std::fmax(s, 0.0));
        const double target =
            shape_reference(r, cfg.bubble_radius, cfg.wall_sigma).value;
        for (std::size_t a = 0; a < e.size(); ++a) {
            const int ia = slot[static_cast<std::size_t>(e.active_indices[a])];
            if (ia < 0) continue;
            rhs[static_cast<std::size_t>(ia)] += e.value[a] * target;
            for (std::size_t b = 0; b < e.size(); ++b) {
                const int ib =
                    slot[static_cast<std::size_t>(e.active_indices[b])];
                if (ib < 0) continue;
                A(ia, ib) += e.value[a] * e.value[b];
            }
        }
    }

    const EigenResult ae = jacobi_eigen(A);
    pr.condition_estimate = condition_number(ae.eigenvalues);
    const std::vector<double> z =
        spd_solve(A, rhs, 1e-13, &pr.dropped_normal_modes);

    pr.coefficients.assign(static_cast<std::size_t>(pr.n_raw), 0.0);
    for (std::size_t k = 0; k < f.active.size(); ++k)
        pr.coefficients[static_cast<std::size_t>(f.active[k])] = z[k];

    // Zone errors, on a fine independent sample set.
    const int m = 4000;
    for (int k = 0; k <= m; ++k) {
        const double s = cfg.s_max * static_cast<double>(k) /
                         static_cast<double>(m);
        const RadialEvaluation q = evaluate_radial(f, pr.coefficients, s);
        if (!q.inside) continue;
        const double r = cfg.scale_R0 * std::sqrt(std::fmax(s, 0.0));
        const ShapeReference ref =
            shape_reference(r, cfg.bubble_radius, cfg.wall_sigma);
        const int zi = static_cast<int>(zone_of(cfg, s));
        ZoneError& ze = pr.zone[zi];
        ze.value = std::fmax(ze.value, std::fabs(q.value - ref.value));
        ze.first = std::fmax(ze.first, std::fabs(q.d_r - ref.first));
        ze.second = std::fmax(ze.second, std::fabs(q.d_rr - ref.second));
        ze.reference_scale = std::fmax(ze.reference_scale, std::fabs(ref.value));
        ++ze.samples;
        if (zone_of(cfg, s) == RadialZone::FlatCollar)
            pr.collar_max_abs = std::fmax(pr.collar_max_abs,
                                          std::fabs(q.value));
    }
    return pr;
}

}  // namespace design
}  // namespace warp_hypad
