#pragma once

#include <algorithm>

#include "program_c/catalogue.hpp"

// ---------------------------------------------------------------------------
// Astrometric transformations, epoch propagation, light time and uncertainty.
//
// THE CONVERSION, written out because every sign and factor matters:
//
//     r_hat  = (cos d cos a, cos d sin a, sin d)
//     p_hat  = (-sin a, cos a, 0)              increasing right ascension
//     q_hat  = (-sin d cos a, -sin d sin a, cos d)   increasing declination
//
//     x = D r_hat
//     v = v_r r_hat  +  V(mu_alpha*, D) p_hat  +  V(mu_dec, D) q_hat
//
// with V(mu, D) = 4.740470446 * (mu[mas/yr]/1000) * D[pc]  km/s.
//
// mu_alpha* ALREADY carries cos(dec). It multiplies p_hat directly and is NOT
// divided by cos(dec) anywhere in this file. The coordinate rate of right
// ascension is d(a)/dt = mu_alpha*/cos(dec), which appears only when converting
// back to a spherical coordinate, and the gate tests that relation explicitly at
// high declination -- the one place a missing or spurious cos(dec) is visible.
// ---------------------------------------------------------------------------

namespace program_c {

// --- Observables -> barycentric Cartesian -----------------------------------

struct ConversionResult {
    BarycentricCartesian state;
    bool ok = false;
    std::string message;
};

// Requires a decided distance: the caller must have run decide_distance() and
// must pass the estimator through, so a distance can never enter the geometry
// without its provenance.
inline ConversionResult to_barycentric(const ICRSPhaseSpace& s, Parsec distance,
                                      DistanceEstimator estimator) {
    ConversionResult out;
    if (estimator == DistanceEstimator::None) {
        out.message =
            "no defensible distance estimator: refusing to place this source in "
            "Cartesian space rather than inventing a position";
        return out;
    }
    if (!(distance.value > 0.0) || !std::isfinite(distance.value)) {
        out.message = "distance is not positive and finite";
        return out;
    }

    const ICRSDirection dir{s.right_ascension, s.declination};
    const Vector3 r_hat = unit_vector(dir);
    Vector3 p_hat, q_hat;
    tangent_basis(dir, p_hat, q_hat);

    out.state.position_pc = scale(r_hat, distance.value);

    const KmPerSecond v_alpha = transverse_velocity(s.pm_ra_cosdec, distance);
    const KmPerSecond v_delta = transverse_velocity(s.pm_dec, distance);

    const Vector3 v_km_s =
        add(add(scale(r_hat, s.radial_velocity.value), scale(p_hat, v_alpha.value)),
            scale(q_hat, v_delta.value));

    out.state.velocity_pc_per_year = {
        to_parsec_per_year(KmPerSecond(v_km_s.x)).value,
        to_parsec_per_year(KmPerSecond(v_km_s.y)).value,
        to_parsec_per_year(KmPerSecond(v_km_s.z)).value};
    out.state.epoch = s.epoch;
    out.ok = true;
    out.message = std::string("converted with estimator ") + to_string(estimator);
    return out;
}

// --- Barycentric Cartesian -> observables ------------------------------------
//
// The inverse. This is where d(ra)/dt = mu_alpha*/cos(dec) would appear if a
// coordinate rate were wanted; it is not, because the catalogued quantity is the
// physical rate along p_hat.
inline ICRSPhaseSpace to_observables(const BarycentricCartesian& b) {
    ICRSPhaseSpace s;
    const double d = norm(b.position_pc);
    s.right_ascension = Radian(std::atan2(b.position_pc.y, b.position_pc.x));
    if (s.right_ascension.value < 0.0)
        s.right_ascension = Radian(s.right_ascension.value + 2.0 * kPi);
    s.declination = Radian(std::asin(b.position_pc.z / d));
    s.parallax = parallax_from_distance(Parsec(d));

    const ICRSDirection dir{s.right_ascension, s.declination};
    const Vector3 r_hat = unit_vector(dir);
    Vector3 p_hat, q_hat;
    tangent_basis(dir, p_hat, q_hat);

    const Vector3 v_km_s = {to_km_per_second(ParsecPerYear(b.velocity_pc_per_year.x)).value,
                            to_km_per_second(ParsecPerYear(b.velocity_pc_per_year.y)).value,
                            to_km_per_second(ParsecPerYear(b.velocity_pc_per_year.z)).value};

    s.radial_velocity = KmPerSecond(dot(v_km_s, r_hat));
    s.pm_ra_cosdec = proper_motion_from_transverse(KmPerSecond(dot(v_km_s, p_hat)),
                                                   Parsec(d));
    s.pm_dec = proper_motion_from_transverse(KmPerSecond(dot(v_km_s, q_hat)),
                                            Parsec(d));
    s.epoch = b.epoch;
    return s;
}

// --- Epoch propagation ------------------------------------------------------
//
// A LOCAL KINEMATIC MODEL, not a Galactic orbit integrator. Rectilinear
// barycentric motion ignores Galactic differential rotation and the vertical
// potential. The threshold below is where that starts to matter at the per-cent
// level for a typical disc star, and it is a stated modelling choice rather than
// a computed bound.
inline constexpr double kRectilinearValidityYears = 1.0e5;

struct PropagationResult {
    BarycentricCartesian state;
    double displacement_pc = 0.0;
    double radial_displacement_pc = 0.0;
    double transverse_displacement_pc = 0.0;
    double interval_years = 0.0;
    bool model_validity_warning = false;
    std::string message;
};

inline PropagationResult propagate_rectilinear(
    const BarycentricCartesian& state_at_epoch, JulianYear target_epoch) {
    PropagationResult r;
    const double dt = target_epoch.value - state_at_epoch.epoch.value;
    r.interval_years = dt;

    r.state.position_pc =
        add(state_at_epoch.position_pc, scale(state_at_epoch.velocity_pc_per_year, dt));
    r.state.velocity_pc_per_year = state_at_epoch.velocity_pc_per_year;
    r.state.epoch = target_epoch;

    const Vector3 disp = subtract(r.state.position_pc, state_at_epoch.position_pc);
    r.displacement_pc = norm(disp);

    // Split along the ORIGINAL line of sight, which is what "radial" means for a
    // catalogue epoch.
    const double d0 = norm(state_at_epoch.position_pc);
    if (d0 > 0.0) {
        const Vector3 r_hat = scale(state_at_epoch.position_pc, 1.0 / d0);
        r.radial_displacement_pc = dot(disp, r_hat);
        const Vector3 transverse =
            subtract(disp, scale(r_hat, r.radial_displacement_pc));
        r.transverse_displacement_pc = norm(transverse);
    }

    r.model_validity_warning = std::fabs(dt) > kRectilinearValidityYears;
    r.message = r.model_validity_warning
                    ? "interval exceeds the stated rectilinear validity of 1e5 yr: "
                      "Galactic differential rotation and the vertical potential "
                      "are NOT modelled, so this position is a kinematic "
                      "extrapolation and not an orbit"
                    : "rectilinear barycentric propagation within the stated "
                      "validity interval";
    return r;
}

// --- Light time -------------------------------------------------------------
//
// Three distinct states, never conflated:
//   catalogue-coordinate state at t
//   geometric state at the observer time    (where it is now)
//   apparent retarded state                 (where the arriving light left)
//
//     t_emit = t_observe - |x_star(t_emit) - x_obs(t_observe)| / c
//
// solved by fixed-point iteration. NO relativistic aberration: that is C1.

struct LightTimeResult {
    BarycentricCartesian geometric_state;   // at the observer time
    BarycentricCartesian retarded_state;    // at t_emit
    double light_time_years = 0.0;
    double emission_epoch_jyear = 0.0;
    int iterations = 0;
    double residual_years = 0.0;
    bool converged = false;
};

inline LightTimeResult solve_light_time(const BarycentricCartesian& star_at_epoch,
                                       const Vector3& observer_position_pc,
                                       JulianYear observe_epoch,
                                       double tolerance_years = 1e-12,
                                       int max_iterations = 100) {
    LightTimeResult out;
    out.geometric_state =
        propagate_rectilinear(star_at_epoch, observe_epoch).state;

    double t_emit = observe_epoch.value;
    for (out.iterations = 0; out.iterations < max_iterations; ++out.iterations) {
        const BarycentricCartesian s =
            propagate_rectilinear(star_at_epoch, JulianYear(t_emit)).state;
        const double range =
            norm(subtract(s.position_pc, observer_position_pc));
        const double next = observe_epoch.value - range / kSpeedOfLightParsecPerYear;
        out.residual_years = std::fabs(next - t_emit);
        t_emit = next;
        if (out.residual_years < tolerance_years) {
            out.converged = true;
            ++out.iterations;
            break;
        }
    }

    out.emission_epoch_jyear = t_emit;
    out.light_time_years = observe_epoch.value - t_emit;
    out.retarded_state = propagate_rectilinear(star_at_epoch, JulianYear(t_emit)).state;
    return out;
}

// --- Uncertainty propagation ------------------------------------------------
//
// Sigma_y = J Sigma_x J^T, with J the Jacobian of
//   (ra, dec, parallax, pm_ra_cosdec, pm_dec, v_r) -> (x, y, z, vx, vy, vz)
//
// METHOD: CENTRED FINITE DIFFERENCES, and the label travels with the result.
// HyPAD is deliberately not used: the available scalar carries three derivative
// directions and this Jacobian needs six independent inputs, so a batched scheme
// costs six evaluations against twelve for centred differences and gains nothing
// in accuracy on a map this smooth. That is a stated trade, not an oversight.

struct CovarianceResult {
    double sigma[36] = {};    // 6x6 over (x,y,z,vx,vy,vz), pc and pc/yr
    double position_sigma_pc[3] = {};      // principal axes of the position block
    double velocity_sigma_pc_per_year[3] = {};
    double position_rms_pc = 0.0;
    double velocity_rms_pc_per_year = 0.0;
    bool have_covariance = false;
    std::string jacobian_method = "centred finite differences";
    std::string message;

    double& operator()(int i, int j) { return sigma[6 * i + j]; }
    const double& operator()(int i, int j) const { return sigma[6 * i + j]; }
};

namespace detail {

inline void state_vector_from(const ICRSPhaseSpace& s, double out[6]) {
    out[0] = s.right_ascension.value;
    out[1] = s.declination.value;
    out[2] = s.parallax.value;
    out[3] = s.pm_ra_cosdec.value;
    out[4] = s.pm_dec.value;
    out[5] = s.radial_velocity.value;
}

inline ICRSPhaseSpace observables_from(const double v[6], JulianYear epoch) {
    ICRSPhaseSpace s;
    s.right_ascension = Radian(v[0]);
    s.declination = Radian(v[1]);
    s.parallax = Milliarcsecond(v[2]);
    s.pm_ra_cosdec = MasPerYear(v[3]);
    s.pm_dec = MasPerYear(v[4]);
    s.radial_velocity = KmPerSecond(v[5]);
    s.epoch = epoch;
    return s;
}

// The map used for the Jacobian: inverse parallax is applied UNCONDITIONALLY
// here, which is legitimate only because the caller has already passed the
// quality gate -- the covariance of a source without a defensible distance is
// not computed at all.
inline bool cartesian_from(const double v[6], JulianYear epoch, double out[6]) {
    const ICRSPhaseSpace s = observables_from(v, epoch);
    if (!(s.parallax.value > 0.0)) return false;
    const Parsec d = inverse_parallax_distance(s.parallax);
    const ConversionResult c =
        to_barycentric(s, d, DistanceEstimator::InverseParallaxHighQuality);
    if (!c.ok) return false;
    out[0] = c.state.position_pc.x;
    out[1] = c.state.position_pc.y;
    out[2] = c.state.position_pc.z;
    out[3] = c.state.velocity_pc_per_year.x;
    out[4] = c.state.velocity_pc_per_year.y;
    out[5] = c.state.velocity_pc_per_year.z;
    return true;
}

// Symmetric 3x3 eigenvalues by cyclic Jacobi, ascending. For the uncertainty
// ellipsoid only: a diagnostic, and comparisons are legitimate here.
inline void symmetric_eigenvalues3(const double a_in[3][3], double out[3]) {
    double a[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) a[i][j] = 0.5 * (a_in[i][j] + a_in[j][i]);
    for (int sweep = 0; sweep < 60; ++sweep) {
        double off = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = i + 1; j < 3; ++j) off = std::fmax(off, std::fabs(a[i][j]));
        if (off <= 1e-18) break;
        for (int p = 0; p < 2; ++p)
            for (int q = p + 1; q < 3; ++q) {
                if (std::fabs(a[p][q]) <= 1e-300) continue;
                const double th = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                const double t = (th >= 0.0 ? 1.0 : -1.0) /
                                 (std::fabs(th) + std::sqrt(th * th + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < 3; ++k) {
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k) {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
    }
    for (int i = 0; i < 3; ++i) out[i] = a[i][i];
    std::sort(out, out + 3);
}

}  // namespace detail

inline CovarianceResult propagate_covariance(const ICRSPhaseSpace& s,
                                            const Mat6& sigma_x,
                                            bool have_covariance) {
    CovarianceResult out;
    out.have_covariance = have_covariance;
    if (!have_covariance) {
        out.message =
            "no catalogue covariance for this source: uncertainty is UNKNOWN and "
            "is reported as such rather than assumed zero";
        return out;
    }

    double x0[6];
    detail::state_vector_from(s, x0);

    // Relative steps, floored, so each coordinate is differenced on its own
    // scale rather than a shared absolute step.
    double J[6][6] = {};
    for (int k = 0; k < 6; ++k) {
        const double scale_k = std::fmax(std::fabs(x0[k]), 1.0);
        const double h = 1e-6 * scale_k;
        double xp[6], xm[6], yp[6], ym[6];
        for (int i = 0; i < 6; ++i) { xp[i] = x0[i]; xm[i] = x0[i]; }
        xp[k] += h;
        xm[k] -= h;
        if (!detail::cartesian_from(xp, s.epoch, yp) ||
            !detail::cartesian_from(xm, s.epoch, ym)) {
            out.message = "the conversion failed inside the Jacobian: no "
                          "covariance is reported";
            out.have_covariance = false;
            return out;
        }
        for (int i = 0; i < 6; ++i) J[i][k] = (yp[i] - ym[i]) / (2.0 * h);
    }

    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
            double acc = 0.0;
            for (int a = 0; a < 6; ++a)
                for (int b = 0; b < 6; ++b)
                    acc += J[i][a] * sigma_x(a, b) * J[j][b];
            out(i, j) = acc;
        }

    double pos[3][3], vel[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            pos[i][j] = out(i, j);
            vel[i][j] = out(i + 3, j + 3);
        }
    double ep[3], ev[3];
    detail::symmetric_eigenvalues3(pos, ep);
    detail::symmetric_eigenvalues3(vel, ev);
    for (int i = 0; i < 3; ++i) {
        out.position_sigma_pc[i] = std::sqrt(std::fmax(ep[i], 0.0));
        out.velocity_sigma_pc_per_year[i] = std::sqrt(std::fmax(ev[i], 0.0));
    }
    out.position_rms_pc = std::sqrt(std::fmax(out(0, 0) + out(1, 1) + out(2, 2), 0.0));
    out.velocity_rms_pc_per_year =
        std::sqrt(std::fmax(out(3, 3) + out(4, 4) + out(5, 5), 0.0));
    out.message = "Sigma_y = J Sigma_x J^T, J by " + out.jacobian_method;
    return out;
}

// Uncertainty growth under rectilinear propagation. With x(t) = x0 + v dt the
// position block grows as
//     Sigma_xx(t) = Sigma_xx + dt (Sigma_xv + Sigma_vx) + dt^2 Sigma_vv
// which is exact for this model, so it is written out rather than differenced.
inline double propagated_position_rms_pc(const CovarianceResult& c, double dt_years) {
    if (!c.have_covariance) return 0.0;
    double trace = 0.0;
    for (int i = 0; i < 3; ++i)
        trace += c(i, i) + dt_years * (c(i, i + 3) + c(i + 3, i)) +
                 dt_years * dt_years * c(i + 3, i + 3);
    return std::sqrt(std::fmax(trace, 0.0));
}

}  // namespace program_c
