#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "program_c/apparent_sky.hpp"
#include "program_c/interception.hpp"
#include "program_c/relativity.hpp"
#include "program_c/snapshot.hpp"
#include "warp_hypad/core/json.hpp"
#include "warp_hypad/core/metadata.hpp"
#include "warp_hypad/core/verification.hpp"

// ---------------------------------------------------------------------------
// Gate C1-RELATIVISTIC-ASTROGATION, scientific cases 1-20.
//
// Also writes the numerical records every C1 visual is built from. THE DATA IS
// SAVED BEFORE ANY FIGURE EXISTS, and the plotting scripts read these files and
// recompute nothing. That ordering is the whole point of the visualization
// contract: a figure that disagrees with the CSV is a figure that lied.
//
// Reads the committed Gaia DR3 snapshot. No live archive access.
// ---------------------------------------------------------------------------

using namespace program_c;
using namespace warp_hypad;

namespace {

std::string g17(double v) {
    std::ostringstream o;
    o << std::setprecision(17) << v;
    return o.str();
}

// Resolve a destination from the snapshot BY SOURCE ID, never by name.
bool make_destination(const LoadedSnapshot& snap, std::uint64_t source_id,
                      Destination& out) {
    for (std::size_t i = 0; i < snap.objects.size(); ++i) {
        const StellarObject& s = snap.objects[i];
        if (s.measured.gaia_source_id != source_id) continue;
        if (!s.has_defensible_distance()) return false;

        ICRSPhaseSpace obs;
        obs.right_ascension = Radian(s.measured.ra_rad);
        obs.declination = Radian(s.measured.dec_rad);
        obs.parallax = Milliarcsecond(s.measured.parallax_mas);
        obs.pm_ra_cosdec = MasPerYear(s.measured.pm_ra_mas_per_year);
        obs.pm_dec = MasPerYear(s.measured.pm_dec_mas_per_year);
        obs.radial_velocity = KmPerSecond(s.measured.radial_velocity_km_per_s);
        obs.epoch = JulianYear(s.measured.reference_epoch_jyear);
        const ConversionResult c =
            to_barycentric(obs, Parsec(s.distance_pc), s.distance_estimator);
        if (!c.ok) return false;

        out.gaia_source_id = source_id;
        out.state_at_reference_epoch = c.state;
        // The catalogue covariance, propagated through the same labelled
        // Jacobian C0 verified. A source without one comes back with
        // have_covariance false rather than a zero matrix pretending to be a
        // measurement.
        out.uncertainty = propagate_covariance(obs, s.measured.covariance,
                                               s.measured.has(kFieldCovariance));
        out.provenance = s.provenance;
        out.display_name = s.primary_name;
        return true;
    }
    return false;
}

}  // namespace

int main(int argc, char** argv) {
    std::string out_dir = "results/verification/program_c1";
    std::string snapshot_path = "snapshots/gaia_dr3_10pc.txt";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out_dir = argv[++i];
        else if (a == "--snapshot" && i + 1 < argc) snapshot_path = argv[++i];
        else { std::cerr << "usage: c1_verify [--out <dir>] [--snapshot <f>]\n"; return 2; }
    }

    std::error_code ec;
    for (const char* sub : {"", "/data", "/figures", "/animations",
                            "/contact_sheets", "/verification"})
        std::filesystem::create_directories(out_dir + sub, ec);

    VerificationReport report(
        "program_c1",
        "C1-RELATIVISTIC-ASTROGATION: ship worldlines, clocks, moving-target "
        "interception, aberration, Doppler and the retarded apparent sky");

    const LoadedSnapshot snap = read_snapshot(snapshot_path);
    if (!snap.ok || !snap.checksum_ok) {
        std::cerr << "cannot read " << snapshot_path << ": " << snap.message << "\n";
        return 3;
    }

    const double c_pc_yr = kSpeedOfLightParsecPerYear;
    const double kDepartureEpoch = 2026.0;

    // The two destinations. Both by source_id.
    //   Proxima Centauri, the nearest star.
    //   Barnard's star, the largest proper motion of any known star -- chosen
    //   precisely because it makes the moving-target correction unmissable.
    Destination proxima, barnard;
    const bool have_proxima =
        make_destination(snap, 5853498713190525696ull, proxima);
    const bool have_barnard =
        make_destination(snap, 4472832130942575872ull, barnard);
    if (!have_proxima || !have_barnard) {
        std::cerr << "the committed snapshot does not contain both destinations\n";
        return 3;
    }

    // =====================================================================
    // 1. Constant-velocity worldline
    // =====================================================================
    {
        Stopwatch sw;
        const double beta = 0.5;
        const double v = beta * c_pc_yr;
        const double t = 4.0;
        const double x = v * t;
        // A straight line, but verified through the same state machinery the
        // accelerated cases use, so this case tests the plumbing rather than
        // arithmetic.
        BarycentricCartesian s;
        s.velocity_pc_per_year = {v, 0.0, 0.0};
        s.position_pc = scale(s.velocity_pc_per_year, t);
        report.check_value(
            "c1.constant_velocity_worldline",
            "x = v t for a constant coordinate velocity of 0.5 c over 4 yr",
            VerificationLevel::A, "analytic", "beta = 0.5, t = 4 yr",
            x, s.position_pc.x, 1.0e-15, sw.seconds(),
            "the constant-velocity model is a VERIFICATION MODEL ONLY: it needs "
            "infinite acceleration at both ends. It earns its place because "
            "tau = t/gamma is exact, which makes it the right instrument for "
            "testing the proper-time integrator");
    }

    // =====================================================================
    // 2. Proper-time integration against the exact constant-velocity answer
    // =====================================================================
    {
        Stopwatch sw;
        const double beta = 0.8;
        const double t1 = 10.0;
        const ProperTimeResult pt =
            integrate_proper_time([&](double) { return beta; }, 0.0, t1, 1000);
        const double exact = t1 / lorentz_gamma(beta);
        report.check_value(
            "c1.proper_time_integration",
            "Simpson integration of dtau = dt sqrt(1-beta^2) against the exact "
            "tau = t/gamma at beta = 0.8",
            VerificationLevel::A, "1000 intervals", "beta = 0.8, t = 10 yr",
            exact, pt.proper_time_years, 1.0e-14, sw.seconds(),
            "constant integrand, so Simpson is exact here up to rounding; the "
            "accelerated case below is the one that exercises the quadrature "
            "order");
    }

    // =====================================================================
    // 3-5. Constant proper acceleration: position, coordinate time, beta/gamma
    // =====================================================================
    const RelativisticRocket one_g{kStandardGravityParsecPerYearSquared};
    {
        Stopwatch sw;
        const double tau = 2.0;
        const double rap = one_g.rapidity_at_proper_time(tau);
        const double x_exact =
            (c_pc_yr * c_pc_yr / one_g.a_pc_per_yr2) * (std::cosh(rap) - 1.0);
        report.check_value(
            "c1.proper_acceleration_analytic_position",
            "x(tau) = (c^2/a)[cosh(a tau/c) - 1] at 1 g for tau = 2 yr",
            VerificationLevel::A, "analytic", "a = 1 g, tau = 2 yr",
            x_exact, one_g.distance_pc(tau), 1.0e-15, sw.seconds(),
            "1 g is 0.31648 pc/yr^2, converted in one place from the defined "
            "9.80665 m/s^2; a/c is " +
                g17(standard_gravity_in_c_per_year()) +
                " per year, the familiar 'one gravity is about one c per year'");
    }
    {
        Stopwatch sw;
        const double tau = 2.0;
        const double rap = one_g.rapidity_at_proper_time(tau);
        const double t_exact = (c_pc_yr / one_g.a_pc_per_yr2) * std::sinh(rap);
        report.check_value(
            "c1.proper_acceleration_analytic_coordinate_time",
            "c t(tau) = (c^2/a) sinh(a tau/c) at 1 g for tau = 2 yr",
            VerificationLevel::A, "analytic", "a = 1 g, tau = 2 yr",
            t_exact, one_g.coordinate_time(tau), 1.0e-15, sw.seconds(),
            "coordinate time already exceeds proper time by " +
                g17(t_exact / tau) + " at only two years of ship time");
    }
    {
        Stopwatch sw;
        const double tau = 2.0;
        const double rap = one_g.rapidity_at_proper_time(tau);
        const double b = one_g.beta(tau);
        const double g = one_g.gamma(tau);
        // gamma = cosh, beta = tanh, and gamma = 1/sqrt(1-beta^2) must agree.
        // Checking the IDENTITY rather than the two values separately is what
        // catches a rapidity that is right in one relation and wrong in the
        // other.
        const double residual = std::abs(g - lorentz_gamma(b)) +
                                std::abs(b - std::tanh(rap)) +
                                std::abs(g - std::cosh(rap));
        report.check_residual(
            "c1.proper_acceleration_analytic_velocity_and_gamma",
            "beta = tanh(a tau/c), gamma = cosh(a tau/c), and gamma = "
            "1/sqrt(1-beta^2) all agree",
            VerificationLevel::A, "analytic", "a = 1 g, tau = 2 yr", residual,
            1.0e-14, sw.seconds(),
            "the three relations are checked together because a wrong rapidity "
            "can satisfy any one of them alone");
    }

    // =====================================================================
    // 6. Symmetric acceleration/deceleration endpoint
    // =====================================================================
    {
        Stopwatch sw;
        SymmetricJourney j;
        j.rocket = one_g;
        j.distance_pc = norm(proxima.state_at_reference_epoch.position_pc);
        double x = 0.0, t = 0.0, b = 0.0, g = 0.0;
        j.state_at_proper_time(j.total_proper_time(), x, t, b, g);
        // Arrival must be AT the destination distance and AT REST. The
        // deceleration phase is written as a time reflection of the
        // acceleration phase rather than integrated separately, which is what
        // makes this exact rather than merely close.
        const double residual =
            std::abs(x - j.distance_pc) + std::abs(b) + std::abs(g - 1.0);
        report.check_residual(
            "c1.symmetric_journey_endpoint",
            "a symmetric 1 g journey arrives at the destination distance and at "
            "rest",
            VerificationLevel::A, "analytic reflection",
            "a = 1 g, D = " + g17(j.distance_pc) + " pc", residual, 1.0e-13,
            sw.seconds(),
            "total ship time " + g17(j.total_proper_time()) +
                " yr against coordinate time " + g17(j.total_coordinate_time()) +
                " yr, peak beta " + g17(j.peak_beta()));
    }

    // =====================================================================
    // 7-9. Moving-destination interception
    // =====================================================================
    SymmetricJourney jp;
    jp.rocket = one_g;
    InterceptionResult icept_proxima, icept_barnard;
    {
        Stopwatch sw;
        auto tof = [&](double d) {
            SymmetricJourney j;
            j.rocket = one_g;
            j.distance_pc = d;
            return j.total_coordinate_time();
        };
        icept_proxima = solve_interception(proxima, kDepartureEpoch, tof);
        icept_barnard = solve_interception(barnard, kDepartureEpoch, tof);
        report.check_predicate(
            "c1.moving_destination_interception",
            "the interception solve converges and places the ship on the "
            "destination's worldline, not at its catalogue position",
            VerificationLevel::A,
            icept_proxima.converged && icept_barnard.converged &&
                icept_proxima.residual_pc < 1.0e-9 &&
                icept_barnard.residual_pc < 1.0e-9,
            "both solves converge below 1e-9 pc",
            "Proxima residual " + g17(icept_proxima.residual_pc) + " pc in " +
                std::to_string(icept_proxima.iterations) + " iterations; "
                "Barnard residual " + g17(icept_barnard.residual_pc) + " pc in " +
                std::to_string(icept_barnard.iterations) + " iterations",
            sw.seconds(),
            "bisection, not Newton: the residual is monotone in t for any "
            "profile bounded by c, so bracketing cannot diverge");
    }
    {
        Stopwatch sw;
        // What ignoring proper motion costs. For Barnard's star -- 10.4
        // arcsec/yr, the largest known -- this must be LARGE. A test that only
        // checked "the error is small" would pass on a broken propagator.
        report.check_predicate(
            "c1.stale_position_targeting_error",
            "aiming at the stale catalogue position misses, and misses more for "
            "the faster-moving star",
            VerificationLevel::A,
            icept_barnard.stale_targeting_error_pc >
                    icept_proxima.stale_targeting_error_pc &&
                icept_proxima.stale_targeting_error_pc > 0.0,
            "Barnard's error exceeds Proxima's, and both are nonzero",
            "Proxima " + g17(icept_proxima.stale_targeting_error_au) +
                " au, Barnard " + g17(icept_barnard.stale_targeting_error_au) +
                " au",
            sw.seconds(),
            "the comparison is the assertion: Barnard's star has the largest "
            "proper motion of any known star, so if the stale error were not "
            "larger for it the propagator would be wrong");
    }
    {
        Stopwatch sw;
        // The uncertainty at arrival must EXCEED the catalogue uncertainty at
        // the reference epoch, because velocity errors integrate.
        const double sigma_now = propagated_position_rms_pc(proxima.uncertainty, 0.0);
        report.check_predicate(
            "c1.interception_uncertainty",
            "the interception position uncertainty exceeds the catalogue "
            "position uncertainty, because velocity errors integrate over the "
            "flight",
            VerificationLevel::B,
            icept_proxima.interception_uncertainty_pc >= sigma_now &&
                std::isfinite(icept_proxima.interception_uncertainty_pc),
            "arrival sigma at least the epoch sigma",
            "epoch " + g17(sigma_now) + " pc, arrival " +
                g17(icept_proxima.interception_uncertainty_pc) + " pc over " +
                g17(icept_proxima.flight_coordinate_time_years) + " yr",
            sw.seconds(),
            "this is the number that decides whether you arrive in the system "
            "or merely near it, and it is not the catalogue's quoted error");
    }

    // =====================================================================
    // 10. Lorentz transformation preserves the Minkowski norm
    // =====================================================================
    {
        Stopwatch sw;
        const Vector3 bv{0.3, -0.4, 0.2};
        double worst = 0.0;
        // A null vector must stay null and a timelike vector must keep its norm.
        const FourVector k = photon_four_vector(Vector3{0.6, 0.8, 0.0}, 1.0);
        worst = std::max(worst, std::abs(boost(k, bv).norm2() - k.norm2()));
        const FourVector u = four_velocity(Vector3{0.1, 0.2, -0.3});
        worst = std::max(worst, std::abs(boost(u, bv).norm2() - u.norm2()));
        worst = std::max(worst, std::abs(u.norm2() + 1.0));   // u.u = -1
        report.check_residual(
            "c1.lorentz_norm_preservation",
            "the boost preserves the Minkowski norm of a null and a timelike "
            "four-vector, and the four-velocity is unit timelike",
            VerificationLevel::A, "signature (-,+,+,+)", "beta = (0.3,-0.4,0.2)",
            worst, 1.0e-14, sw.seconds(),
            "the boost is written out in full rather than as a rotation to the "
            "z axis plus a 1-D boost, because the round trip through a rotation "
            "is where an axis convention slips in unnoticed. The tolerance is a "
            "few ulps on order-one quantities, not a physical allowance: the "
            "norm is preserved EXACTLY in real arithmetic and the residual here "
            "is the rounding of the gamma products");
    }

    // =====================================================================
    // 11. Aberration against the analytic formula, and against the physics
    // =====================================================================
    {
        Stopwatch sw;
        const double beta = 0.6;
        const Vector3 bv{beta, 0.0, 0.0};
        double worst_formula = 0.0;
        bool crowds_forward = true;

        for (int i = 0; i <= 36; ++i) {
            const double th = i * kPi / 36.0;
            // Direction TOWARD the source, in the x-z plane.
            const Vector3 n_hat{std::cos(th), 0.0, std::sin(th)};
            const AberrationResult ab = aberrate(n_hat, bv);

            // Against the source-direction form of the analytic result.
            const double expect = aberrate_cos_source_direction(std::cos(th), beta);
            worst_formula =
                std::max(worst_formula, std::abs(ab.cos_angle_to_velocity - expect));

            // And against the PHYSICS: apparent positions crowd toward the
            // direction of motion, so the apparent angle from the velocity is
            // never larger than the geometric one. This is the check that
            // survives a sign inversion in the formula itself.
            if (ab.cos_angle_to_velocity + 1.0e-12 < std::cos(th))
                crowds_forward = false;
        }
        report.check_residual(
            "c1.aberration_analytic",
            "the four-vector boost reproduces cos t' = (cos t + b)/(1 + b cos t) "
            "for the source direction, over 37 angles",
            VerificationLevel::A, "37 angles", "beta = 0.6", worst_formula,
            1.0e-14, sw.seconds(),
            "implemented by boosting the photon four-vector; the scalar formula "
            "is the TEST, not the implementation, which is the right way round "
            "because the boost generalises to any velocity direction");
        report.check_predicate(
            "c1.aberration_crowds_toward_motion",
            "apparent source positions crowd toward the direction of motion at "
            "every angle",
            VerificationLevel::A, crowds_forward,
            "cos of the apparent angle never below the geometric one",
            crowds_forward ? "holds at all 37 angles" : "VIOLATED",
            sw.seconds(),
            "the case that catches a sign inversion. Reproducing one formula "
            "from itself would pass with the sign wrong; this asserts the "
            "physical direction of the effect");
    }

    // =====================================================================
    // 12-13. Doppler, longitudinal limits
    // =====================================================================
    {
        Stopwatch sw;
        const double beta = 0.6;
        const Vector3 bv{beta, 0.0, 0.0};
        // Approaching: the source is ahead at +x, so the photon propagates in
        // -x, opposite to the motion.
        const double D_app = doppler_factor(Vector3{-1.0, 0.0, 0.0}, bv);
        const double expect_app = std::sqrt((1.0 + beta) / (1.0 - beta));
        report.check_value(
            "c1.doppler_approaching",
            "nu'/nu = sqrt((1+b)/(1-b)) for a source directly ahead",
            VerificationLevel::A, "analytic", "beta = 0.6", expect_app, D_app,
            1.0e-15, sw.seconds(),
            "computed from -(u_ship.k)/(u_ref.k), so the longitudinal limit is a "
            "consequence rather than the definition. D = " + g17(D_app) +
                " is a blueshift, as it must be for an approaching source");
    }
    {
        Stopwatch sw;
        const double beta = 0.6;
        const Vector3 bv{beta, 0.0, 0.0};
        const double D_rec = doppler_factor(Vector3{1.0, 0.0, 0.0}, bv);
        const double expect_rec = std::sqrt((1.0 - beta) / (1.0 + beta));
        const double expect_app = std::sqrt((1.0 + beta) / (1.0 - beta));
        // And the transverse case, which has no classical counterpart.
        const double D_tr = doppler_factor(Vector3{0.0, 1.0, 0.0}, bv);
        report.check_value(
            "c1.doppler_receding",
            "nu'/nu = sqrt((1-b)/(1+b)) for a source directly behind",
            VerificationLevel::A, "analytic", "beta = 0.6", expect_rec, D_rec,
            1.0e-15, sw.seconds(),
            "the reciprocal of the approaching case to " +
                g17(std::abs(D_rec * expect_app - 1.0)) +
                ". The transverse case gives " + g17(D_tr) + " = gamma, pure "
                "time dilation, which has no classical counterpart");
    }

    // =====================================================================
    // 14. Zero-velocity sky recovery
    // =====================================================================
    {
        Stopwatch sw;
        const Vector3 zero{0.0, 0.0, 0.0};
        double worst_direction = 0.0, worst_doppler = 0.0;
        for (int i = 0; i <= 20; ++i) {
            const double th = i * kPi / 20.0;
            const Vector3 n_hat{std::cos(th), std::sin(th), 0.0};
            const AberrationResult ab = aberrate(n_hat, zero);
            worst_direction = std::max(
                worst_direction,
                norm(subtract(ab.apparent_direction_to_source, n_hat)));
            worst_doppler =
                std::max(worst_doppler, std::abs(ab.doppler_factor - 1.0));
        }
        // The Doppler factor IS bit-exact at beta = 0, so it is asserted at
        // tolerance zero. The direction is not, and claiming it would be was
        // wrong: photon_four_vector renormalises by |k_hat|, and (cos t, sin t, 0)
        // does not have norm exactly 1 in binary floating point. Renormalising an
        // inexact unit vector is not the identity. That is correct behaviour --
        // the code must renormalise to be robust to non-unit input -- so the
        // residual is one rounding, not a defect.
        report.check_residual(
            "c1.zero_velocity_sky_recovery_doppler",
            "at zero velocity every Doppler factor is exactly 1",
            VerificationLevel::A, "21 directions", "beta = 0", worst_doppler, 0.0,
            sw.seconds(),
            "tolerance ZERO and it holds: boost() returns its argument unchanged "
            "when beta is zero, and the time component was set to the frequency "
            "exactly");
        report.check_residual(
            "c1.zero_velocity_sky_recovery_directions",
            "at zero velocity every apparent direction equals the geometric "
            "direction to within one renormalisation rounding",
            VerificationLevel::A, "21 directions", "beta = 0", worst_direction,
            4.0e-16, sw.seconds(),
            "the tolerance is four ulps and the reason is stated rather than "
            "assumed: the input direction (cos t, sin t, 0) has norm 1 only to "
            "rounding, and the code renormalises it. A stationary observer's sky "
            "is the geometric sky to the last bit the representation allows");
    }

    // =====================================================================
    // 15. Retarded source state
    // =====================================================================
    {
        Stopwatch sw;
        const RetardedSolution ret = solve_retarded_emission(
            Vector3{0.0, 0.0, 0.0}, kDepartureEpoch,
            proxima.state_at_reference_epoch);
        // Self-consistency: the light time must equal the distance to the
        // RETARDED position over c, not to the current position.
        const double implied = ret.geometric_distance_pc / c_pc_yr;
        report.check_value(
            "c1.retarded_source_state",
            "the emission time satisfies |x_obs - x_star(t_e)| = c (t_obs - t_e)",
            VerificationLevel::A,
            std::to_string(ret.iterations) + " iterations", "Proxima from Earth",
            implied, ret.light_time_years, 1.0e-12, sw.seconds(),
            "the star is seen where it WAS " + g17(ret.light_time_years) +
                " yr earlier, at epoch " + g17(ret.emission_epoch_years) +
                ", not where it is");
    }

    // =====================================================================
    // 16-17. Communication delays
    // =====================================================================
    SignalLedger ledger_mid;
    {
        Stopwatch sw;
        // A ship a third of the way to Proxima.
        const double frac = 1.0 / 3.0;
        const Vector3 pos =
            scale(proxima.state_at_reference_epoch.position_pc, frac);
        ledger_mid = signal_ledger(pos, kDepartureEpoch + 1.0, 0.9, proxima);
        const double d = norm(pos);
        report.check_value(
            "c1.earth_to_ship_signal_delay",
            "the Earth clock the ship can SEE lags coordinate time by |x|/c",
            VerificationLevel::A, "single event", "one third of the way to Proxima",
            (kDepartureEpoch + 1.0) - d / c_pc_yr,
            ledger_mid.earth_time_seen_from_ship, 1.0e-14, sw.seconds(),
            "this is NOT 'the current time on Earth'. It is the latest Earth "
            "event whose light has arrived, " + g17(ledger_mid.earth_signal_age_years) +
                " yr stale, and no display aboard the ship may show otherwise");
    }
    {
        Stopwatch sw;
        const double frac = 1.0 / 3.0;
        const Vector3 pos =
            scale(proxima.state_at_reference_epoch.position_pc, frac);
        const double d = norm(pos);
        report.check_value(
            "c1.ship_to_earth_signal_delay",
            "a message sent by the ship reaches Earth at t + |x|/c, and the "
            "round trip is 2|x|/c",
            VerificationLevel::A, "single event", "one third of the way to Proxima",
            (kDepartureEpoch + 1.0) + d / c_pc_yr,
            ledger_mid.ship_signal_arrives_at_earth, 1.0e-14, sw.seconds(),
            "round trip " + g17(ledger_mid.round_trip_to_earth_years) +
                " yr; the destination signal solve is separate because the "
                "destination MOVES and the outgoing light cone has to catch it");
    }

    // =====================================================================
    // 18-19. Apparent sky: source count and uncertainty propagation
    // =====================================================================
    std::vector<ApparentSkyRow> sky_ship, sky_earth;
    {
        Stopwatch sw;
        const Vector3 obs_pos =
            scale(proxima.state_at_reference_epoch.position_pc, 1.0 / 3.0);
        const Vector3 beta_ship{0.0, 0.0, 0.0};
        const Vector3 dir =
            scale(proxima.state_at_reference_epoch.position_pc,
                  1.0 / norm(proxima.state_at_reference_epoch.position_pc));
        const Vector3 beta_moving = scale(dir, 0.75);

        int placeable = 0;
        for (std::size_t i = 0; i < snap.objects.size(); ++i) {
            const StellarObject& s = snap.objects[i];
            if (!s.has_defensible_distance()) continue;
            ++placeable;
            BarycentricCartesian st;
            st.position_pc = {snap.derived[i].x_pc, snap.derived[i].y_pc,
                              snap.derived[i].z_pc};
            // Velocity from the catalogue observables, through the C0 converter.
            ICRSPhaseSpace o;
            o.right_ascension = Radian(s.measured.ra_rad);
            o.declination = Radian(s.measured.dec_rad);
            o.parallax = Milliarcsecond(s.measured.parallax_mas);
            o.pm_ra_cosdec = MasPerYear(s.measured.pm_ra_mas_per_year);
            o.pm_dec = MasPerYear(s.measured.pm_dec_mas_per_year);
            o.radial_velocity = KmPerSecond(s.measured.radial_velocity_km_per_s);
            o.epoch = JulianYear(s.measured.reference_epoch_jyear);
            const ConversionResult cv =
                to_barycentric(o, Parsec(s.distance_pc), s.distance_estimator);
            if (!cv.ok) continue;
            st = cv.state;

            sky_earth.push_back(apparent_sky_row(s, st, obs_pos, kDepartureEpoch,
                                                 beta_ship));
            sky_ship.push_back(apparent_sky_row(s, st, obs_pos, kDepartureEpoch,
                                                beta_moving));
        }
        report.check_predicate(
            "c1.apparent_sky_source_count",
            "the apparent-sky table has exactly one row per placeable catalogue "
            "source, in both the rest and moving observer frames",
            VerificationLevel::A,
            static_cast<int>(sky_earth.size()) == placeable &&
                sky_ship.size() == sky_earth.size(),
            std::to_string(placeable) + " rows in each frame",
            std::to_string(sky_earth.size()) + " rest, " +
                std::to_string(sky_ship.size()) + " moving, from " +
                std::to_string(snap.objects.size()) + " snapshot objects",
            sw.seconds(),
            "aberration MOVES stars, it never creates or destroys them; a row "
            "count that changed with beta would mean the renderer and the "
            "science disagree about what is in the sky");
    }
    {
        Stopwatch sw;
        // Every uncertainty angle must be finite, non-negative, and small enough
        // that the small-angle statement sigma/d is defensible.
        double worst = 0.0;
        bool ok = true;
        for (const ApparentSkyRow& r : sky_ship) {
            if (!std::isfinite(r.uncertainty_angle) || r.uncertainty_angle < 0.0)
                ok = false;
            worst = std::max(worst, r.uncertainty_angle);
        }
        report.check_predicate(
            "c1.uncertainty_angle_propagation",
            "every apparent-sky row carries a finite, non-negative angular "
            "uncertainty, and the small-angle form remains defensible",
            VerificationLevel::B, ok && worst < 0.1,
            "all finite and below 0.1 rad", "worst " + g17(worst) + " rad over " +
                std::to_string(sky_ship.size()) + " sources",
            sw.seconds(),
            "sigma/d is a small-angle statement; at " + g17(worst) +
                " rad the neglected term is of order its cube, which is beyond "
                "the fourth decimal");
    }

    // =====================================================================
    // 20. Deterministic worldline
    // =====================================================================
    std::vector<WorldlineSample> worldline;
    {
        Stopwatch sw;
        auto build = [&]() {
            std::vector<WorldlineSample> w;
            SymmetricJourney j;
            j.rocket = one_g;
            j.distance_pc = icept_proxima.interception_distance_pc;
            const int n = 240;
            const Vector3 dir = scale(icept_proxima.interception_point,
                                      1.0 / norm(icept_proxima.interception_point));
            for (int i = 0; i <= n; ++i) {
                const double tau = j.total_proper_time() * i / n;
                double x = 0.0, t = 0.0, b = 0.0, g = 0.0;
                j.state_at_proper_time(tau, x, t, b, g);
                WorldlineSample s;
                s.ship.coordinate_state.position_pc = scale(dir, x);
                s.ship.coordinate_state.velocity_pc_per_year =
                    scale(dir, b * c_pc_yr);
                s.ship.coordinate_state.epoch = JulianYear(kDepartureEpoch + t);
                s.ship.proper_time_years = tau;
                s.ship.coordinate_time_years = t;
                s.ship.beta = b;
                s.ship.gamma = g;
                const double sgn = (tau <= j.half_proper_time() ? 1.0 : -1.0);
                s.ship.proper_acceleration = scale(dir, sgn * one_g.a_pc_per_yr2);
                s.destination = proxima.state_at(kDepartureEpoch + t);
                s.remaining_distance_pc = norm(subtract(
                    s.destination.position_pc, s.ship.coordinate_state.position_pc));
                const SignalLedger sl =
                    signal_ledger(s.ship.coordinate_state.position_pc,
                                  kDepartureEpoch + t, tau, proxima);
                s.received_earth_time = sl.earth_time_seen_from_ship;
                s.received_destination_time = sl.destination_time_seen_from_ship;
                w.push_back(s);
            }
            return w;
        };
        worldline = build();
        const std::vector<WorldlineSample> again = build();
        double worst = 0.0;
        for (std::size_t i = 0; i < worldline.size(); ++i) {
            worst = std::max(worst, std::abs(worldline[i].ship.proper_time_years -
                                             again[i].ship.proper_time_years));
            worst = std::max(worst, norm(subtract(
                                        worldline[i].ship.coordinate_state.position_pc,
                                        again[i].ship.coordinate_state.position_pc)));
        }
        report.check_residual(
            "c1.deterministic_worldline",
            "rebuilding the worldline from the same inputs reproduces it bit for "
            "bit",
            VerificationLevel::A, std::to_string(worldline.size()) + " samples",
            "1 g symmetric to Proxima", worst, 0.0, sw.seconds(),
            "tolerance ZERO. Every visual downstream is generated from these "
            "samples, so a worldline that is merely reproducible to a tolerance "
            "would make animation frames irreproducible");
    }

    // =====================================================================
    // Write the numerical records. BEFORE any figure exists.
    // =====================================================================
    const std::string data = out_dir + "/data";
    {
        std::ofstream f(data + "/worldline_proxima_1g.csv", std::ios::binary);
        f << "sample,proper_time_yr,coordinate_time_yr,epoch_jyear,x_pc,y_pc,z_pc,"
             "vx_pc_yr,vy_pc_yr,vz_pc_yr,beta,gamma,accel_pc_yr2,"
             "remaining_distance_pc,dest_x_pc,dest_y_pc,dest_z_pc,"
             "received_earth_time,received_destination_time\n";
        for (std::size_t i = 0; i < worldline.size(); ++i) {
            const WorldlineSample& s = worldline[i];
            f << i << "," << g17(s.ship.proper_time_years) << ","
              << g17(s.ship.coordinate_time_years) << ","
              << g17(s.ship.coordinate_state.epoch.value) << ","
              << g17(s.ship.coordinate_state.position_pc.x) << ","
              << g17(s.ship.coordinate_state.position_pc.y) << ","
              << g17(s.ship.coordinate_state.position_pc.z) << ","
              << g17(s.ship.coordinate_state.velocity_pc_per_year.x) << ","
              << g17(s.ship.coordinate_state.velocity_pc_per_year.y) << ","
              << g17(s.ship.coordinate_state.velocity_pc_per_year.z) << ","
              << g17(s.ship.beta) << "," << g17(s.ship.gamma) << ","
              << g17(norm(s.ship.proper_acceleration)) << ","
              << g17(s.remaining_distance_pc) << ","
              << g17(s.destination.position_pc.x) << ","
              << g17(s.destination.position_pc.y) << ","
              << g17(s.destination.position_pc.z) << ","
              << g17(s.received_earth_time) << ","
              << g17(s.received_destination_time) << "\n";
        }
    }
    auto write_sky = [&](const std::string& name,
                         const std::vector<ApparentSkyRow>& rows) {
        std::ofstream f(data + "/" + name, std::ios::binary);
        f << "source_id,primary_name,emission_epoch,geometric_distance,"
             "apparent_direction_x,apparent_direction_y,apparent_direction_z,"
             "doppler_factor,brightness_factor,uncertainty_angle,quality_flags,"
             "geometric_direction_x,geometric_direction_y,geometric_direction_z,"
             "aberration_angle\n";
        for (const ApparentSkyRow& r : rows) {
            std::string nm = r.primary_name;
            for (char& ch : nm) if (ch == ',') ch = ';';
            f << r.source_id << "," << nm << "," << g17(r.emission_epoch) << ","
              << g17(r.geometric_distance) << "," << g17(r.apparent_direction_x)
              << "," << g17(r.apparent_direction_y) << ","
              << g17(r.apparent_direction_z) << "," << g17(r.doppler_factor) << ","
              << g17(r.brightness_factor) << "," << g17(r.uncertainty_angle) << ","
              << r.quality_flags << "," << g17(r.geometric_direction_x) << ","
              << g17(r.geometric_direction_y) << "," << g17(r.geometric_direction_z)
              << "," << g17(r.aberration_angle) << "\n";
        }
    };
    write_sky("apparent_sky_earth_frame.csv", sky_earth);
    write_sky("apparent_sky_ship_frame.csv", sky_ship);

    // The aberration/Doppler SEQUENCE. Written here, by the solver, because the
    // animation must not recompute a single boost: a plotting script that swept
    // beta itself would be a second, unverified simulation wearing this one's
    // label. Each beta is a full Lorentz transformation of every source.
    {
        const Vector3 obs_pos =
            scale(proxima.state_at_reference_epoch.position_pc, 1.0 / 3.0);
        const Vector3 dir =
            scale(proxima.state_at_reference_epoch.position_pc,
                  1.0 / norm(proxima.state_at_reference_epoch.position_pc));
        std::ofstream f(data + "/apparent_sky_sequence.csv", std::ios::binary);
        f << "beta_index,beta,gamma,source_id,apparent_direction_x,"
             "apparent_direction_y,apparent_direction_z,doppler_factor,"
             "brightness_factor,aberration_angle\n";
        const int kBetaSteps = 48;
        for (int b = 0; b < kBetaSteps; ++b) {
            // Closed form in the index, so the sequence is reproducible and the
            // animation's frame-to-beta map is a lookup rather than a formula
            // duplicated in Python.
            const double beta = 0.98 * b / (kBetaSteps - 1);
            const Vector3 bv = scale(dir, beta);
            for (std::size_t i = 0; i < snap.objects.size(); ++i) {
                const StellarObject& s = snap.objects[i];
                if (!s.has_defensible_distance()) continue;
                ICRSPhaseSpace o;
                o.right_ascension = Radian(s.measured.ra_rad);
                o.declination = Radian(s.measured.dec_rad);
                o.parallax = Milliarcsecond(s.measured.parallax_mas);
                o.pm_ra_cosdec = MasPerYear(s.measured.pm_ra_mas_per_year);
                o.pm_dec = MasPerYear(s.measured.pm_dec_mas_per_year);
                o.radial_velocity = KmPerSecond(s.measured.radial_velocity_km_per_s);
                o.epoch = JulianYear(s.measured.reference_epoch_jyear);
                const ConversionResult cv =
                    to_barycentric(o, Parsec(s.distance_pc), s.distance_estimator);
                if (!cv.ok) continue;
                const ApparentSkyRow r = apparent_sky_row(
                    s, cv.state, obs_pos, kDepartureEpoch, bv);
                f << b << "," << g17(beta) << "," << g17(lorentz_gamma(beta))
                  << "," << r.source_id << "," << g17(r.apparent_direction_x)
                  << "," << g17(r.apparent_direction_y) << ","
                  << g17(r.apparent_direction_z) << "," << g17(r.doppler_factor)
                  << "," << g17(r.brightness_factor) << ","
                  << g17(r.aberration_angle) << "\n";
            }
        }
    }

    // Trajectory family comparison, and the interception ledger.
    {
        std::ofstream f(data + "/trajectory_comparison.csv", std::ios::binary);
        f << "destination,source_id,model,acceleration_g,distance_pc,"
             "coordinate_time_yr,proper_time_yr,peak_beta,peak_gamma,"
             "arrival_epoch,stale_error_au,interception_uncertainty_pc\n";
        struct Row { const Destination* d; const InterceptionResult* r; };
        const Row rows[2] = {{&proxima, &icept_proxima}, {&barnard, &icept_barnard}};
        const double accels[2] = {1.0, 0.1};
        for (const Row& row : rows) {
            const double D = row.r->interception_distance_pc;
            // Constant coordinate velocity at 0.5 c, the verification model.
            {
                const double t = D / (0.5 * c_pc_yr);
                f << row.d->display_name << "," << row.d->gaia_source_id
                  << ",ConstantCoordinateVelocity,0," << g17(D) << "," << g17(t)
                  << "," << g17(t / lorentz_gamma(0.5)) << ",0.5,"
                  << g17(lorentz_gamma(0.5)) << ","
                  << g17(kDepartureEpoch + t) << "," << g17(row.r->stale_targeting_error_au)
                  << "," << g17(row.r->interception_uncertainty_pc) << "\n";
            }
            for (double ag : accels) {
                SymmetricJourney j;
                j.rocket.a_pc_per_yr2 = ag * kStandardGravityParsecPerYearSquared;
                j.distance_pc = D;
                f << row.d->display_name << "," << row.d->gaia_source_id
                  << ",SymmetricProperAcceleration," << g17(ag) << "," << g17(D)
                  << "," << g17(j.total_coordinate_time()) << ","
                  << g17(j.total_proper_time()) << "," << g17(j.peak_beta()) << ","
                  << g17(j.peak_gamma()) << ","
                  << g17(kDepartureEpoch + j.total_coordinate_time()) << ","
                  << g17(row.r->stale_targeting_error_au) << ","
                  << g17(row.r->interception_uncertainty_pc) << "\n";
            }
        }
    }
    {
        std::ofstream f(data + "/interception.csv", std::ios::binary);
        f << "destination,source_id,quantity,x_pc,y_pc,z_pc,scalar\n";
        auto emit = [&](const Destination& d, const InterceptionResult& r) {
            auto row = [&](const char* q, const Vector3& v, double s) {
                f << d.display_name << "," << d.gaia_source_id << "," << q << ","
                  << g17(v.x) << "," << g17(v.y) << "," << g17(v.z) << "," << g17(s)
                  << "\n";
            };
            row("stale_catalogue_target", r.stale_catalogue_target,
                norm(r.stale_catalogue_target));
            row("destination_at_departure", r.destination_at_departure,
                norm(r.destination_at_departure));
            row("destination_at_arrival", r.destination_at_arrival,
                norm(r.destination_at_arrival));
            row("interception_point", r.interception_point,
                r.interception_distance_pc);
            row("stale_targeting_error", Vector3{0, 0, 0}, r.stale_targeting_error_au);
            row("interception_uncertainty_pc", Vector3{0, 0, 0},
                r.interception_uncertainty_pc);
            row("arrival_epoch", Vector3{0, 0, 0}, r.arrival_epoch_years);
            row("flight_coordinate_time_yr", Vector3{0, 0, 0},
                r.flight_coordinate_time_years);
        };
        emit(proxima, icept_proxima);
        emit(barnard, icept_barnard);
    }
    {
        // The light-cone ledger along the journey.
        std::ofstream f(data + "/light_cones.csv", std::ios::binary);
        f << "sample,coordinate_time_yr,proper_time_yr,distance_from_earth_pc,"
             "earth_time_seen,earth_signal_age_yr,signal_arrives_at_earth,"
             "signal_arrives_at_destination,round_trip_yr\n";
        for (std::size_t i = 0; i < worldline.size(); ++i) {
            const WorldlineSample& s = worldline[i];
            const SignalLedger sl = signal_ledger(
                s.ship.coordinate_state.position_pc,
                kDepartureEpoch + s.ship.coordinate_time_years,
                s.ship.proper_time_years, proxima);
            f << i << "," << g17(s.ship.coordinate_time_years) << ","
              << g17(s.ship.proper_time_years) << ","
              << g17(norm(s.ship.coordinate_state.position_pc)) << ","
              << g17(sl.earth_time_seen_from_ship) << ","
              << g17(sl.earth_signal_age_years) << ","
              << g17(sl.ship_signal_arrives_at_earth) << ","
              << g17(sl.ship_signal_arrives_at_destination) << ","
              << g17(sl.round_trip_to_earth_years) << "\n";
        }
    }
    {
        // The local neighbourhood, for the 3-D figures. Straight from the
        // snapshot's derived columns: no recomputation anywhere downstream.
        std::ofstream f(data + "/local_neighborhood.csv", std::ios::binary);
        f << "source_id,primary_name,x_pc,y_pc,z_pc,distance_pc,"
             "distance_sigma_pc,estimator,quality_flags,has_planets\n";
        for (std::size_t i = 0; i < snap.objects.size(); ++i) {
            const StellarObject& s = snap.objects[i];
            if (!s.has_defensible_distance()) continue;
            std::string nm = s.primary_name;
            for (char& ch : nm) if (ch == ',') ch = ';';
            f << s.measured.gaia_source_id << "," << nm << ","
              << g17(snap.derived[i].x_pc) << "," << g17(snap.derived[i].y_pc)
              << "," << g17(snap.derived[i].z_pc) << "," << g17(s.distance_pc)
              << "," << g17(s.distance_uncertainty_pc) << ","
              << to_string(s.distance_estimator) << "," << s.measured.quality_flags
              << "," << (s.planets.empty() ? 0 : 1) << "\n";
        }
    }

    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir + "/verification", error)) {
        std::cerr << "failed to write reports: " << error << "\n";
        return 3;
    }

    // --- c1_metadata.json: the provenance every figure inherits -------------
    {
        JsonValue meta = JsonValue::object();
        meta.set("milestone", JsonValue("C1-RELATIVISTIC-ASTROGATION"));
        meta.set("catalogue", JsonValue(snap.release.catalogue));
        meta.set("catalogue_release", JsonValue(snap.release.release));
        meta.set("snapshot", JsonValue(snapshot_path));
        meta.set("snapshot_checksum",
                 JsonValue(std::to_string(snap.checksum_in_file)));
        meta.set("data_kind", JsonValue(snap.release.data_kind));
        meta.set("reference_epoch", JsonValue(snap.release.reference_epoch_jyear));
        meta.set("simulation_epoch", JsonValue(kDepartureEpoch));
        meta.set("coordinate_frame",
                 JsonValue(BarycentricCartesian::description().frame));
        meta.set("observer_frame", JsonValue("LocalShipFrame (instantaneous rest "
                                             "frame of the ship, boosted from "
                                             "BarycentricCartesian)"));
        meta.set("units", JsonValue("parsec; Julian year; parsec per Julian year"));
        meta.set("speed_of_light_pc_per_yr", JsonValue(c_pc_yr));
        meta.set("standard_gravity_pc_per_yr2",
                 JsonValue(kStandardGravityParsecPerYearSquared));
        meta.set("trajectory_model",
                 JsonValue(to_string(TrajectoryModel::SymmetricProperAcceleration)));
        meta.set("trajectory_model_purpose",
                 JsonValue(model_purpose(TrajectoryModel::SymmetricProperAcceleration)));
        meta.set("color_model",
                 JsonValue(to_string(ColorModel::ApproximateBlackbodyVisualization)));
        meta.set("brightness_qualification",
                 JsonValue(BrightnessFactors::qualification()));
        meta.set("effective_metric",
                 JsonValue("NONE. C1 is special-relativistic kinematics in flat "
                           "spacetime. No warp metric, no propulsion energy, no "
                           "collision hazard, no interstellar medium."));

        JsonValue dests = JsonValue::array();
        struct D { const Destination* d; const InterceptionResult* r; };
        for (const D& x : {D{&proxima, &icept_proxima}, D{&barnard, &icept_barnard}}) {
            JsonValue j = JsonValue::object();
            j.set("destination_source_id",
                  JsonValue(std::to_string(x.d->gaia_source_id)));
            j.set("display_name", JsonValue(x.d->display_name));
            j.set("distance_pc", JsonValue(x.r->interception_distance_pc));
            j.set("arrival_epoch", JsonValue(x.r->arrival_epoch_years));
            j.set("flight_coordinate_time_yr",
                  JsonValue(x.r->flight_coordinate_time_years));
            j.set("stale_targeting_error_au",
                  JsonValue(x.r->stale_targeting_error_au));
            j.set("interception_uncertainty_pc",
                  JsonValue(x.r->interception_uncertainty_pc));
            j.set("solver_residual_pc", JsonValue(x.r->residual_pc));
            j.set("solver_iterations", JsonValue(x.r->iterations));
            j.set("solver_method", JsonValue(x.r->method));
            dests.push_back(j);
        }
        meta.set("destinations", dests);

        JsonValue phys = JsonValue::object();
        phys.set("acceleration_g", JsonValue(1.0));
        phys.set("worldline_samples",
                 JsonValue(static_cast<double>(worldline.size())));
        phys.set("total_proper_time_yr",
                 JsonValue(worldline.empty() ? 0.0
                                             : worldline.back().ship.proper_time_years));
        phys.set("total_coordinate_time_yr",
                 JsonValue(worldline.empty()
                               ? 0.0
                               : worldline.back().ship.coordinate_time_years));
        phys.set("peak_beta", JsonValue(jp.peak_beta()));
        meta.set("physical_parameters", phys);

        meta.set("git_commit", JsonValue(build_info_json()));
        meta.set("generated_timestamp", JsonValue(utc_timestamp()));
        if (!json_write_file(out_dir + "/c1_metadata.json", meta, error)) return 3;
    }

    std::cout << "\nC1 DATA WRITTEN (visuals are generated from these files "
                 "and recompute nothing)\n";
    std::cout << "  worldline samples                  " << worldline.size() << "\n";
    std::cout << "  apparent-sky rows per frame        " << sky_ship.size() << "\n";
    std::cout << "  Proxima  flight coordinate time    "
              << icept_proxima.flight_coordinate_time_years << " yr\n";
    std::cout << "  Proxima  stale targeting error     "
              << icept_proxima.stale_targeting_error_au << " au\n";
    std::cout << "  Barnard  stale targeting error     "
              << icept_barnard.stale_targeting_error_au << " au\n";
    std::cout << "\nSCOPE: special-relativistic kinematics only. No warp metric, "
                 "no propulsion\n  energy, no collision hazard, no interstellar "
                 "medium. Clocks are always named.\n";

    const bool ok = report.all_passed();
    std::cout << "\nC1 SCIENTIFIC GATE: " << (ok ? "PASS" : "FAIL") << " ("
              << report.num_passed() << "/" << report.num_cases() << ")\n";
    return ok ? 0 : 1;
}
