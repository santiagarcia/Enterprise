#include <cmath>
#include <fstream>
#include <sstream>

#include "c0_checks.hpp"
#include "c0_fixtures.hpp"

// ---------------------------------------------------------------------------
// Gate C0-CATALOGUE-UNIVERSE.
//
// SCOPE. C0 verifies the MACHINERY of a catalogue-grounded star database:
// schema, provenance, units, coordinate transformations, epoch propagation,
// light time, uncertainty and a deterministic offline snapshot. It does NOT
// verify agreement with Gaia DR3, because no Gaia record was retrieved -- the
// fixtures are SYNTHETIC and the gate says so in its own output.
// ---------------------------------------------------------------------------

namespace program_c {
namespace checks {

using namespace program_c::fixtures;
using warp_hypad::Stopwatch;
using warp_hypad::VerificationLevel;
using warp_hypad::VerificationReport;

namespace {

std::string sci(double v) {
    std::ostringstream o;
    o.setf(std::ios::scientific);
    o.precision(4);
    o << v;
    return o.str();
}

ICRSPhaseSpace observables_of(const StellarObject& s) {
    ICRSPhaseSpace o;
    o.right_ascension = Radian(s.measured.ra_rad);
    o.declination = Radian(s.measured.dec_rad);
    o.parallax = Milliarcsecond(s.measured.parallax_mas);
    o.pm_ra_cosdec = MasPerYear(s.measured.pm_ra_mas_per_year);
    o.pm_dec = MasPerYear(s.measured.pm_dec_mas_per_year);
    o.radial_velocity = KmPerSecond(s.measured.radial_velocity_km_per_s);
    o.epoch = JulianYear(s.measured.reference_epoch_jyear);
    return o;
}

const StellarObject* by_name(const std::vector<StellarObject>& v,
                             const std::string& n) {
    for (const StellarObject& s : v)
        if (s.primary_name == n) return &s;
    return nullptr;
}

}  // namespace

void run_c0_checks(VerificationReport& report, const std::string& out_dir,
                  C0Metrics& metrics) {
    const Snapshot snap50 = build_snapshot(50.0);
    const FixtureSet fx = synthetic_catalogue();
    metrics.data_kind = snap50.release.data_kind;
    metrics.catalogue_release = snap50.release.release;
    metrics.adql_query = snap50.release.adql_query;
    metrics.rows_in = fx.report.rows_in;
    metrics.accepted = fx.report.accepted;
    metrics.rejected = static_cast<int>(fx.report.rejected.size());
    metrics.duplicates_collapsed = fx.report.duplicates_collapsed;

    // =======================================================================
    // The honesty case, first, so it cannot be missed.
    // =======================================================================
    {
        Stopwatch sw;
        bool all_tagged = true;
        for (const StellarObject& s : snap50.objects)
            if (s.provenance.is_real_catalogue()) all_tagged = false;
        report.check_predicate(
            "c0.fixtures_are_declared_synthetic",
            "every committed fixture is tagged SYNTHETIC or "
            "LITERATURE_FROM_MEMORY, and none claims to be retrieved catalogue "
            "data",
            VerificationLevel::A,
            all_tagged && snap50.release.data_kind == std::string(kSyntheticTag),
            "no fixture claims a real catalogue",
            "snapshot data_kind = " + snap50.release.data_kind + "; " +
                std::to_string(snap50.objects.size()) + " objects, all tagged",
            sw.seconds(),
            "this milestone had no network access, so no Gaia DR3 record was "
            "retrieved. Writing plausible numbers under a real catalogue's name "
            "would be fabricating records and the gate would be verifying a "
            "fiction. tools/import_gaia.py is what produces fixtures with genuine "
            "provenance, as a separate optional step");
    }

    // =======================================================================
    // 3. Units.
    // =======================================================================
    {
        Stopwatch sw;
        double worst = 0.0;
        // Round trips through every boundary conversion.
        worst = std::fmax(worst, std::fabs(to_degree(to_radian(Degree(37.5))).value - 37.5));
        worst = std::fmax(worst,
                          std::fabs(to_milliarcsecond(to_radian(Milliarcsecond(768.5))).value -
                                    768.5));
        worst = std::fmax(worst, std::fabs(to_parsec(to_metre(Parsec(3.2))).value - 3.2));
        worst = std::fmax(worst,
                          std::fabs(to_julian_year(to_second(JulianYear(11.0))).value - 11.0));
        worst = std::fmax(worst,
                          std::fabs(to_km_per_second(to_parsec_per_year(KmPerSecond(42.0))).value -
                                    42.0));
        // The two factors most often wrong, checked against their definitions.
        const double au_per_yr = kAuKm / kJulianYearSeconds;
        worst = std::fmax(worst, std::fabs(kAuPerYearInKmPerSecond - au_per_yr));
        // 1 mas/yr at 1 pc is 4.740470446 km/s.
        worst = std::fmax(worst,
                          std::fabs(transverse_velocity(MasPerYear(1000.0), Parsec(1.0)).value -
                                    au_per_yr));
        metrics.worst_unit_roundtrip = worst;

        report.check_residual(
            "c0.unit_conversions_round_trip",
            "every boundary conversion round-trips, and 1 arcsec/yr at 1 pc is "
            "4.740470446 km/s",
            VerificationLevel::A, "n/a", "6 conversions plus 2 factor checks",
            worst, 1e-12, sw.seconds(),
            "the transverse-velocity factor is the one astrometry gets wrong most "
            "often, so it is checked against its definition rather than a quoted "
            "value");

        // Speed of light in the internal system, derived not quoted: a star at
        // 1 pc must be 3.2616 light-years away.
        const double ly_per_pc = 1.0 / light_travel_time(Parsec(1.0)).value;
        report.check_value(
            "c0.speed_of_light_in_internal_units",
            "1 pc corresponds to 3.261563777 light-years, so c = 0.3066 pc/yr",
            VerificationLevel::B, "n/a", "derived from c and the Julian year",
            3.261563777, 1.0 / kSpeedOfLightParsecPerYear, 1e-6, sw.seconds());
        (void)ly_per_pc;
    }

    // =======================================================================
    // 4. Direction to unit vector.
    // =======================================================================
    {
        Stopwatch sw;
        double worst = 0.0;
        struct Case { double ra, dec; double x, y, z; };
        const Case cases[5] = {{0, 0, 1, 0, 0},
                               {90, 0, 0, 1, 0},
                               {0, 90, 0, 0, 1},
                               {180, 0, -1, 0, 0},
                               {0, -90, 0, 0, -1}};
        for (const Case& c : cases) {
            const ICRSDirection d{to_radian(Degree(c.ra)), to_radian(Degree(c.dec))};
            const Vector3 u = unit_vector(d);
            worst = std::fmax(worst, std::fabs(u.x - c.x));
            worst = std::fmax(worst, std::fabs(u.y - c.y));
            worst = std::fmax(worst, std::fabs(u.z - c.z));
            worst = std::fmax(worst, std::fabs(norm(u) - 1.0));
        }
        // The tangent basis must be orthonormal and orthogonal to r_hat.
        for (double dec : {-60.0, -10.0, 0.0, 25.0, 75.0}) {
            const ICRSDirection d{to_radian(Degree(43.0)), to_radian(Degree(dec))};
            const Vector3 r = unit_vector(d);
            Vector3 p, q;
            tangent_basis(d, p, q);
            worst = std::fmax(worst, std::fabs(dot(r, p)));
            worst = std::fmax(worst, std::fabs(dot(r, q)));
            worst = std::fmax(worst, std::fabs(dot(p, q)));
            worst = std::fmax(worst, std::fabs(norm(p) - 1.0));
            worst = std::fmax(worst, std::fabs(norm(q) - 1.0));
        }
        report.check_residual(
            "c0.ra_dec_to_unit_vector",
            "the ICRS unit vector is correct at the cardinal directions, and the "
            "tangent basis (p_hat, q_hat) is orthonormal and orthogonal to r_hat",
            VerificationLevel::A, "n/a", "5 cardinal points, 5 declinations",
            worst, 1e-15, sw.seconds());
    }

    // =======================================================================
    // 5 and 6. The distance decision.
    // =======================================================================
    {
        Stopwatch sw;
        const StellarObject* good = by_name(fx.objects, "SYN-AXIS-10PC");
        const StellarObject* poor = by_name(fx.objects, "SYN-POOR-PARALLAX");
        const StellarObject* neg = by_name(fx.objects, "SYN-NEGATIVE-PARALLAX");

        const bool good_ok =
            good && good->distance_estimator ==
                        DistanceEstimator::InverseParallaxHighQuality &&
            std::fabs(good->distance_pc - 10.0) < 1e-12;

        report.check_predicate(
            "c0.parallax_distance_high_quality",
            "a parallax of 100 mas with signal-to-noise 200 yields exactly 10 pc "
            "by inverse parallax, with the estimator recorded",
            VerificationLevel::B, good_ok,
            "estimator InverseParallaxHighQuality, distance 10 pc",
            good ? std::string(to_string(good->distance_estimator)) + ", d = " +
                       sci(good->distance_pc) + " +- " +
                       sci(good->distance_uncertainty_pc) + " pc"
                 : std::string("fixture missing"),
            sw.seconds());

        const bool poor_rejected =
            poor && poor->distance_estimator == DistanceEstimator::None &&
            poor->measured.flagged(kQualityParallaxLowSignalToNoise);
        const bool neg_rejected =
            neg && neg->distance_estimator == DistanceEstimator::None &&
            neg->measured.flagged(kQualityParallaxNonPositive);

        report.check_predicate(
            "c0.poor_parallax_is_rejected_not_inverted",
            "a source with parallax signal-to-noise 2 and one with a NEGATIVE "
            "parallax both get estimator None, with distinct reasons, and are "
            "RETAINED with a flag rather than deleted",
            VerificationLevel::A, poor_rejected && neg_rejected,
            "both rejected for inverse parallax, both retained",
            poor && neg
                ? std::string("snr=2 flagged low-signal-to-noise; negative "
                              "parallax flagged non-positive; both still in the "
                              "catalogue")
                : std::string("fixtures missing"),
            sw.seconds(),
            "d = 1/parallax is NOT a universally unbiased estimator, and the "
            "quality gate that licenses it lives in exactly one place. A source "
            "that fails it keeps its astrometry and loses only its distance");
    }

    // =======================================================================
    // 7. The proper-motion convention.
    // =======================================================================
    {
        Stopwatch sw;
        // A star at dec = 60 deg with pure mu_alpha* = 500 mas/yr. Propagating
        // and converting back must give a COORDINATE rate
        //     d(ra)/dt = mu_alpha* / cos(dec)
        // which is 2x the physical rate at dec = 60. A missing or spurious
        // cos(dec) shows up here and nowhere else.
        const StellarObject* s = by_name(fx.objects, "SYN-HIGHDEC-CONVENTION");
        double worst = 1.0;
        double measured_ratio = 0.0;
        if (s) {
            const ICRSPhaseSpace o = observables_of(*s);
            const ConversionResult c =
                to_barycentric(o, Parsec(s->distance_pc), s->distance_estimator);
            const double dt = 100.0;
            const PropagationResult p = propagate_rectilinear(
                c.state, JulianYear(o.epoch.value + dt));
            const ICRSPhaseSpace back = to_observables(p.state);

            const double dra_rad = back.right_ascension.value - o.right_ascension.value;
            const double dra_mas = to_milliarcsecond(Radian(dra_rad)).value;
            const double expected_coordinate_rate =
                o.pm_ra_cosdec.value / std::cos(o.declination.value);
            measured_ratio = (dra_mas / dt) / o.pm_ra_cosdec.value;
            worst = std::fabs(dra_mas / dt - expected_coordinate_rate) /
                    std::fabs(expected_coordinate_rate);
        }
        metrics.proper_motion_convention_error = worst;

        report.check_residual(
            "c0.proper_motion_convention_mu_alpha_star",
            "the catalogued pm_ra is mu_alpha* = mu_alpha cos(dec): at dec = 60 "
            "the coordinate rate d(ra)/dt comes out at mu_alpha*/cos(dec), twice "
            "the physical rate",
            VerificationLevel::B, "n/a",
            "dec = 60 deg, mu_alpha* = 500 mas/yr, 100 yr; measured "
            "d(ra)/dt / mu_alpha* = " + sci(measured_ratio) + " (expected 2)",
            worst, 1e-3, sw.seconds(),
            "the single test that distinguishes the two proper-motion conventions. "
            "mu_alpha* multiplies p_hat directly and is never divided by cos(dec) "
            "in the Cartesian conversion; the division appears only on the way "
            "back to a spherical coordinate");
    }

    // =======================================================================
    // 8. Radial velocity, and 9/10 the Cartesian conversion and its inverse.
    // =======================================================================
    {
        Stopwatch sw;
        const StellarObject* s = by_name(fx.objects, "SYN-RADIAL-ONLY");
        double radial_fraction = 0.0;
        bool ok = false;
        if (s) {
            const ICRSPhaseSpace o = observables_of(*s);
            const ConversionResult c =
                to_barycentric(o, Parsec(s->distance_pc), s->distance_estimator);
            const PropagationResult p =
                propagate_rectilinear(c.state, JulianYear(o.epoch.value + 1000.0));
            radial_fraction = p.transverse_displacement_pc /
                              std::fmax(p.displacement_pc, 1e-300);
            ok = c.ok && radial_fraction < 1e-12 && p.radial_displacement_pc > 0.0;
        }
        report.check_predicate(
            "c0.radial_velocity_contribution_is_purely_radial",
            "a star with zero proper motion and positive radial velocity moves "
            "only along the line of sight",
            VerificationLevel::B, ok,
            "transverse displacement fraction below 1e-12",
            "transverse/total = " + sci(radial_fraction), sw.seconds(),
            "radial velocity positive means receding, so the radial displacement "
            "must be positive: the sign is asserted, not just the magnitude");
    }

    {
        Stopwatch sw;
        double worst = 0.0;
        int tested = 0;
        for (const StellarObject& s : fx.objects) {
            if (!s.has_defensible_distance()) continue;
            const ICRSPhaseSpace o = observables_of(s);
            const ConversionResult c =
                to_barycentric(o, Parsec(s.distance_pc), s.distance_estimator);
            if (!c.ok) { worst = 1.0; continue; }
            const ICRSPhaseSpace b = to_observables(c.state);
            // Compare angles modulo 2 pi, and the rest directly.
            double dra = b.right_ascension.value - o.right_ascension.value;
            while (dra > kPi) dra -= 2.0 * kPi;
            while (dra < -kPi) dra += 2.0 * kPi;
            // Near the pole right ascension is degenerate, so it is compared
            // through the unit vector instead of the angle.
            const Vector3 u0 = unit_vector({o.right_ascension, o.declination});
            const Vector3 u1 = unit_vector({b.right_ascension, b.declination});
            worst = std::fmax(worst, norm(subtract(u0, u1)));
            worst = std::fmax(worst, std::fabs(b.parallax.value - o.parallax.value) /
                                         o.parallax.value);
            worst = std::fmax(worst, std::fabs(b.pm_ra_cosdec.value - o.pm_ra_cosdec.value) /
                                         std::fmax(std::fabs(o.pm_ra_cosdec.value), 1.0));
            worst = std::fmax(worst, std::fabs(b.pm_dec.value - o.pm_dec.value) /
                                         std::fmax(std::fabs(o.pm_dec.value), 1.0));
            worst = std::fmax(worst,
                              std::fabs(b.radial_velocity.value - o.radial_velocity.value) /
                                  std::fmax(std::fabs(o.radial_velocity.value), 1.0));
            ++tested;
        }
        metrics.worst_cartesian_roundtrip = worst;
        report.check_residual(
            "c0.astrometry_cartesian_round_trip",
            "observables -> barycentric Cartesian -> observables is the identity, "
            "for every fixture with a defensible distance",
            VerificationLevel::A, "n/a",
            std::to_string(tested) + " sources, unit vector compared near the "
            "pole where right ascension is degenerate",
            worst, 1e-12, sw.seconds());
    }

    // =======================================================================
    // 11, 12, 13. Propagation.
    // =======================================================================
    {
        Stopwatch sw;
        const StellarObject* st = by_name(fx.objects, "SYN-STATIONARY");
        const StellarObject* hp = by_name(fx.objects, "SYN-HIGH-PM");
        double stationary_disp = 1.0, hp_disp = 0.0;
        bool warn_fires = false, warn_silent_within = true;

        if (st) {
            const ICRSPhaseSpace o = observables_of(*st);
            const ConversionResult c =
                to_barycentric(o, Parsec(st->distance_pc), st->distance_estimator);
            const PropagationResult p =
                propagate_rectilinear(c.state, JulianYear(o.epoch.value + 1.0e4));
            stationary_disp = p.displacement_pc;
            warn_silent_within = !p.model_validity_warning;
        }
        if (hp) {
            const ICRSPhaseSpace o = observables_of(*hp);
            const ConversionResult c =
                to_barycentric(o, Parsec(hp->distance_pc), hp->distance_estimator);
            const PropagationResult p =
                propagate_rectilinear(c.state, JulianYear(o.epoch.value + 1000.0));
            hp_disp = p.displacement_pc;
            const PropagationResult far = propagate_rectilinear(
                c.state, JulianYear(o.epoch.value + 1.0e6));
            warn_fires = far.model_validity_warning;
        }
        metrics.stationary_displacement_pc = stationary_disp;
        metrics.high_pm_displacement_pc = hp_disp;

        report.check_residual(
            "c0.stationary_star_does_not_move",
            "a source with zero proper motion and zero radial velocity has "
            "exactly zero displacement after 10 000 yr",
            VerificationLevel::A, "n/a", "10 000 yr", stationary_disp, 0.0,
            sw.seconds(),
            "EXACT zero, not a tolerance: the propagation is x + v dt and v is "
            "identically zero, so any displacement at all would be a bug in the "
            "conversion");

        report.check_predicate(
            "c0.high_proper_motion_propagation",
            "a 5000 mas/yr source at 4 pc moves a measurable distance in 1000 yr, "
            "and the rectilinear validity warning fires past 1e5 yr while staying "
            "silent within it",
            VerificationLevel::B,
            hp_disp > 1e-3 && warn_fires && warn_silent_within,
            "measurable displacement; warning fires only outside the stated "
            "validity",
            "1000 yr displacement " + sci(hp_disp) + " pc; warning at 1e6 yr = " +
                (warn_fires ? "yes" : "no") + ", at 1e4 yr = " +
                (warn_silent_within ? "no" : "yes"),
            sw.seconds(),
            "rectilinear barycentric motion is a LOCAL KINEMATIC MODEL, not a "
            "Galactic orbit integrator, and the warning is how that limit is "
            "communicated rather than assumed known");
    }

    // =======================================================================
    // 14. Light time.
    // =======================================================================
    {
        Stopwatch sw;
        double worst_resid = 0.0;
        int worst_iter = 0;
        bool all_converged = true;
        double check_pc = 0.0, check_yr = 0.0;
        for (const StellarObject& s : fx.objects) {
            if (!s.has_defensible_distance()) continue;
            const ICRSPhaseSpace o = observables_of(s);
            const ConversionResult c =
                to_barycentric(o, Parsec(s.distance_pc), s.distance_estimator);
            const LightTimeResult lt = solve_light_time(
                c.state, {0.0, 0.0, 0.0}, JulianYear(o.epoch.value));
            if (!lt.converged) all_converged = false;
            worst_resid = std::fmax(worst_resid, lt.residual_years);
            worst_iter = std::max(worst_iter, lt.iterations);
            if (s.primary_name == "SYN-AXIS-10PC") {
                check_pc = s.distance_pc;
                check_yr = lt.light_time_years;
            }
        }
        metrics.light_time_residual_years = worst_resid;
        metrics.light_time_iterations = worst_iter;

        // A stationary star at 10 pc: light time is exactly 10 * 3.261563777 yr.
        const double expected = check_pc * (1.0 / kSpeedOfLightParsecPerYear);
        report.check_predicate(
            "c0.light_time_iteration_converges",
            "the light-time fixed point converges for every nearby source, and "
            "for a stationary star at 10 pc gives 32.6 yr",
            VerificationLevel::B,
            all_converged && worst_resid < 1e-12 &&
                std::fabs(check_yr - expected) < 1e-9,
            "converged everywhere, and the analytic value at 10 pc",
            "worst residual " + sci(worst_resid) + " yr in at most " +
                std::to_string(worst_iter) + " iterations; 10 pc gives " +
                sci(check_yr) + " yr against " + sci(expected),
            sw.seconds(),
            "three states are kept distinct: catalogue-coordinate, geometric at "
            "the observer time, and apparent retarded. NO relativistic aberration "
            "-- that is C1, and it is absent rather than approximated");
    }

    // =======================================================================
    // 15, 16. Uncertainty.
    // =======================================================================
    {
        Stopwatch sw;
        const StellarObject* s = by_name(fx.objects, "SYN-AXIS-10PC");
        bool ok = false;
        double growth_ratio = 0.0;
        std::string method;
        if (s) {
            const ICRSPhaseSpace o = observables_of(*s);
            const CovarianceResult cv = propagate_covariance(
                o, s->measured.covariance, s->measured.has(kFieldCovariance));
            method = cv.jacobian_method;
            metrics.position_sigma_pc = cv.position_rms_pc;
            metrics.velocity_sigma_pc_per_year = cv.velocity_rms_pc_per_year;

            // The dominant term for a source on the +x axis is the parallax:
            // sigma_d/d = sigma_parallax/parallax = 0.5/100, so sigma_x ~ 0.05 pc.
            const double expected_x = 10.0 * (0.5 / 100.0);
            const double at0 = propagated_position_rms_pc(cv, 0.0);
            const double at1e4 = propagated_position_rms_pc(cv, 1.0e4);
            growth_ratio = at1e4 / std::fmax(at0, 1e-300);
            ok = cv.have_covariance && std::fabs(cv(0, 0) - expected_x * expected_x) <
                                          0.05 * expected_x * expected_x &&
                 at1e4 > at0;
            metrics.uncertainty_growth_ratio = growth_ratio;
        }
        report.check_predicate(
            "c0.covariance_transformation",
            "Sigma_y = J Sigma_x J^T reproduces the expected radial uncertainty "
            "for a source on the +x axis, sigma_x = d * sigma_parallax/parallax",
            VerificationLevel::B, ok,
            "the dominant parallax term recovered within 5 per cent",
            "sigma_position = " + sci(metrics.position_sigma_pc) + " pc, "
            "sigma_velocity = " + sci(metrics.velocity_sigma_pc_per_year) +
                " pc/yr, Jacobian by " + method,
            sw.seconds(),
            "the Jacobian method is LABELLED and travels with the result. HyPAD is "
            "not used here: the scalar carries three derivative directions and "
            "this Jacobian needs six inputs, so batching costs six evaluations "
            "against twelve for centred differences with no accuracy gain");

        report.check_predicate(
            "c0.uncertainty_grows_with_epoch",
            "the propagated position uncertainty grows with the propagation "
            "interval, as Sigma_xx + dt(Sigma_xv + Sigma_vx) + dt^2 Sigma_vv",
            VerificationLevel::B, growth_ratio > 1.0,
            "strictly growing", "ratio at 1e4 yr = " + sci(growth_ratio),
            sw.seconds(),
            "written out from the model rather than differenced, because for "
            "x = x0 + v dt the expression is exact");
    }

    // =======================================================================
    // 17. Missing fields.
    // =======================================================================
    {
        Stopwatch sw;
        const StellarObject* s = by_name(fx.objects, "SYN-NO-RV");
        const StellarObject* nc = by_name(fx.objects, "SYN-POOR-PARALLAX");
        const bool ok = s && !s->measured.has(kFieldRadialVelocity) &&
                        s->measured.flagged(kQualityMissingRadialVelocity) &&
                        nc && !nc->measured.has(kFieldCovariance) &&
                        nc->measured.flagged(kQualityMissingCovariance);
        // And the covariance path must REFUSE rather than assume zero.
        bool refuses = false;
        if (nc) {
            const CovarianceResult cv = propagate_covariance(
                observables_of(*nc), nc->measured.covariance, false);
            refuses = !cv.have_covariance && !cv.message.empty();
        }
        report.check_predicate(
            "c0.missing_fields_are_flagged_not_nan",
            "a missing radial velocity and a missing covariance are expressed by "
            "CLEAR VALIDITY BITS and quality flags, and the covariance path "
            "refuses to report an uncertainty it does not have",
            VerificationLevel::A, ok && refuses,
            "validity bits clear, flags set, covariance refused",
            ok ? std::string("missing rv and missing covariance both flagged; "
                             "covariance reported as UNKNOWN rather than zero")
               : std::string("fixtures missing"),
            sw.seconds(),
            "NaN alone never expresses missing data here. An absent uncertainty "
            "reported as zero would be the most dangerous possible default");
    }

    // =======================================================================
    // 1, 2, 18, 23, 24, 25. The snapshot.
    // =======================================================================
    {
        Stopwatch sw;
        const std::string path = out_dir + "/c0_snapshot_50pc.txt";
        std::string err;
        const bool written = write_snapshot(snap50, path, err);
        const LoadedSnapshot loaded = read_snapshot(path);

        // Round trip: every field the query layer needs must survive.
        double worst = 0.0;
        bool ids_ok = loaded.objects.size() == snap50.objects.size();
        if (ids_ok)
            for (std::size_t i = 0; i < loaded.objects.size(); ++i) {
                const StellarObject& a = snap50.objects[i];
                const StellarObject& b = loaded.objects[i];
                if (a.measured.gaia_source_id != b.measured.gaia_source_id ||
                    a.primary_name != b.primary_name ||
                    a.distance_estimator != b.distance_estimator ||
                    a.measured.valid_fields != b.measured.valid_fields ||
                    a.measured.quality_flags != b.measured.quality_flags)
                    ids_ok = false;
                worst = std::fmax(worst, std::fabs(a.measured.ra_rad - b.measured.ra_rad));
                worst = std::fmax(worst, std::fabs(a.distance_pc - b.distance_pc));
            }
        metrics.snapshot_roundtrip_error = worst;
        metrics.checksum = snap50.checksum;

        report.check_predicate(
            "c0.snapshot_schema_round_trip",
            "the snapshot writes and reads back with every identifier, validity "
            "bitmask, quality bitmask, estimator and float preserved",
            VerificationLevel::A,
            written && loaded.ok && ids_ok && worst < 1e-15 &&
                loaded.schema_version == kSnapshotSchemaVersion,
            "all fields preserved, schema version " +
                std::to_string(kSnapshotSchemaVersion),
            "worst float difference " + sci(worst) + ", " +
                std::to_string(loaded.objects.size()) + " objects",
            sw.seconds(),
            "floats are written with %.17g, so the round trip is exact rather "
            "than close");

        report.check_predicate(
            "c0.snapshot_is_deterministic",
            "serialising twice gives byte-identical output and the checksum in "
            "the file matches a recomputation",
            VerificationLevel::A,
            serialize(snap50) == serialize(snap50) && loaded.checksum_ok,
            "byte-identical and checksum verified",
            "FNV-1a checksum " + std::to_string(snap50.checksum) +
                (loaded.checksum_ok ? " verified on read" : " MISMATCH"),
            sw.seconds(),
            "records are written in ascending source_id and the checksum covers "
            "the canonical serialization excluding its own line. The checksum is "
            "an integrity check against corruption and reordering, NOT a "
            "cryptographic guarantee, and is labelled so in the manifest");

        // Pin the numeric bit values of the validity and quality masks. The
        // snapshot stores these masks as INTEGERS, so any external reader --
        // tools/snapshot_to_sqlite.py, tools/import_gaia.py -- has to agree on
        // the bit order or it will decode "low parallax signal-to-noise" as
        // "non-positive parallax". That is not a hypothetical: both tools were
        // written against a guessed order and decoded the fixtures wrongly until
        // this case existed. Reordering the enum now fails here loudly.
        report.check_predicate(
            "c0.flag_bit_values_are_pinned",
            "the ValidField and QualityFlag bit values are fixed, because the "
            "snapshot stores them as integers that external tools must decode",
            VerificationLevel::A,
            kFieldRa == 1ull && kFieldDec == 2ull && kFieldParallax == 4ull &&
                kFieldPmRa == 8ull && kFieldPmDec == 16ull &&
                kFieldRadialVelocity == 32ull && kFieldCovariance == 64ull &&
                kFieldDistance == 128ull &&
                kQualityParallaxNonPositive == 1ull &&
                kQualityParallaxLowSignalToNoise == 2ull &&
                kQualityMissingRadialVelocity == 4ull &&
                kQualityMissingCovariance == 8ull &&
                kQualityDuplicatedSource == 16ull &&
                kQualityNonFiniteValue == 32ull &&
                kQualityOutsideRequestedVolume == 64ull,
            "validity bits 1..128 and quality bits 1..64 in declared order",
            "fields " + std::to_string(kFieldRa) + ".." +
                std::to_string(kFieldDistance) + ", quality " +
                std::to_string(kQualityParallaxNonPositive) + ".." +
                std::to_string(kQualityOutsideRequestedVolume),
            sw.seconds(),
            "adding a flag means appending a bit; reordering or reusing one "
            "silently reinterprets every snapshot ever written");

        // Pin the hash itself against the published FNV-1a 64-bit test vectors.
        // Without this, a mistyped offset basis still produces a stable,
        // self-consistent checksum -- the determinism case above passes happily --
        // while silently making the snapshot unreadable by any independent
        // implementation, which is exactly how tools/import_gaia.py and this
        // header drifted apart once already.
        report.check_predicate(
            "c0.checksum_matches_published_fnv1a_vectors",
            "fnv1a() reproduces the published FNV-1a 64-bit test vectors, so an "
            "independent implementation can read what this one writes",
            VerificationLevel::A,
            fnv1a("") == 0xcbf29ce484222325ull &&
                fnv1a("a") == 0xaf63dc4c8601ec8cull &&
                fnv1a("foobar") == 0x85944171f73967e8ull,
            "the three published vectors reproduce exactly",
            "empty -> " + std::to_string(fnv1a("")) + ", \"a\" -> " +
                std::to_string(fnv1a("a")) + ", \"foobar\" -> " +
                std::to_string(fnv1a("foobar")),
            sw.seconds(),
            "the offset basis is 0xcbf29ce484222325 and the prime is "
            "0x100000001b3; tools/import_gaia.py --verify-format checks the same "
            "agreement from outside the C++ build");

        bool prov_ok = !loaded.objects.empty();
        for (const StellarObject& s : loaded.objects)
            if (s.provenance.source_catalogue.empty() ||
                s.provenance.original_identifier.empty())
                prov_ok = false;
        report.check_predicate(
            "c0.provenance_is_retained_through_the_snapshot",
            "every record still carries its source catalogue, release, table, "
            "timestamp, licence and original identifier after a write and read",
            VerificationLevel::A,
            prov_ok && !loaded.release.adql_query.empty() &&
                loaded.release.release == std::string(kPinnedCatalogue),
            "provenance and the exact ADQL both survive",
            "release pinned to " + loaded.release.release + "; ADQL " +
                std::to_string(loaded.release.adql_query.size()) +
                " characters stored beside the result",
            sw.seconds());

        report.check_predicate(
            "c0.catalogue_release_is_pinned",
            "the snapshot records the pinned catalogue, table and reference "
            "epoch, so a silent migration to a later release fails here",
            VerificationLevel::A,
            loaded.release.catalogue == std::string(kPinnedCatalogue) &&
                loaded.release.source_table == std::string(kPinnedTable) &&
                std::fabs(loaded.release.reference_epoch_jyear -
                          kGaiaDR3ReferenceEpochJyear) < 1e-12,
            "Gaia DR3, gaiadr3.gaia_source, J2016.0",
            loaded.release.catalogue + " / " + loaded.release.source_table +
                " / epoch " + sci(loaded.release.reference_epoch_jyear),
            sw.seconds(),
            "migration to Gaia DR4 must be a deliberate versioned update with "
            "regression tests, exactly as the OTIlib pin is enforced elsewhere");

        // Offline queries, on the loaded snapshot only.
        const std::vector<QueryHit> near5 = within_radius(loaded, 5.0);
        const std::vector<QueryHit> near50 = within_radius(loaded, 50.0);
        const std::vector<QueryHit> hosts = exoplanet_hosts(loaded);
        const StellarObject* alias =
            find_by_identifier(loaded, "SYN-ALIAS-EXOHOST");
        metrics.count_within_5pc = static_cast<int>(near5.size());
        metrics.count_within_50pc = static_cast<int>(near50.size());
        metrics.exoplanet_hosts = static_cast<int>(hosts.size());

        report.check_predicate(
            "c0.offline_queries_after_import",
            "the 5 pc and 50 pc volume queries, the nearest-N query and the "
            "exoplanet-host query all run against the LOADED snapshot with no "
            "network access",
            VerificationLevel::A,
            !near5.empty() && near50.size() >= near5.size() && !hosts.empty() &&
                !nearest(loaded, 5).empty(),
            "all four queries return results, 5 pc a subset of 50 pc",
            std::to_string(near5.size()) + " within 5 pc, " +
                std::to_string(near50.size()) + " within 50 pc, " +
                std::to_string(hosts.size()) + " exoplanet hosts",
            sw.seconds(),
            "network access is not part of the simulation runtime: this header "
            "reads files and nothing else, and the online client lives in tools/");

        report.check_predicate(
            "c0.identifier_enrichment_and_exoplanet_association",
            "a source is findable by an alias as well as its primary name, and a "
            "planet record is associated with its host",
            VerificationLevel::A,
            alias != nullptr && alias->primary_name == "SYN-EXOHOST" &&
                !alias->planets.empty() &&
                alias->planets[0].name == "SYN-EXOHOST b",
            "alias resolves to the host, and the host carries its planet",
            alias ? alias->primary_name + " found via alias, " +
                        std::to_string(alias->planets.size()) + " planet(s)"
                  : std::string("alias lookup failed"),
            sw.seconds(),
            "SIMBAD supplies names and cross-identifications only; it is never "
            "the primary star catalogue");
    }

    // =======================================================================
    // Reconciliation: 0 silently dropped records.
    // =======================================================================
    {
        Stopwatch sw;
        // Duplicate handling: the same source id offered twice must collapse to
        // one record, and the collapse must be counted.
        FixtureSet dup = synthetic_catalogue();
        const int before = static_cast<int>(dup.objects.size());
        dup.objects.push_back(dup.objects[0]);
        dup.report.rows_in = static_cast<int>(dup.objects.size());
        Snapshot s;
        s.objects = dup.objects;
        // Collapse by id, counting.
        std::vector<StellarObject> unique;
        int collapsed = 0;
        for (const StellarObject& o : s.objects) {
            bool seen = false;
            for (const StellarObject& u : unique)
                if (u.measured.gaia_source_id == o.measured.gaia_source_id) seen = true;
            if (seen) { ++collapsed; continue; }
            unique.push_back(o);
        }
        s.objects = unique;
        s.ingestion = dup.report;
        s.ingestion.duplicates_collapsed = collapsed;
        s.ingestion.accepted = 0;
        s.ingestion.rejected.clear();
        for (const StellarObject& o : s.objects) {
            if (o.distance_estimator == DistanceEstimator::None)
                s.ingestion.rejected.push_back(
                    {o.measured.gaia_source_id, decide_distance(o.measured).reason});
            else ++s.ingestion.accepted;
        }
        finalise(s);
        // Deliberately NOT metrics.duplicates_collapsed: this is the injected
        // ledger, not the delivered snapshot's.
        metrics.injected_duplicate_rows_in = s.ingestion.rows_in;
        metrics.injected_duplicate_collapsed = collapsed;

        report.check_predicate(
            "c0.duplicate_sources_collapse_and_are_counted",
            "the same source_id offered twice yields one record, and the collapse "
            "is counted rather than silent",
            VerificationLevel::A,
            collapsed == 1 && static_cast<int>(s.objects.size()) == before,
            "one duplicate collapsed, count preserved",
            std::to_string(collapsed) + " collapsed, " +
                std::to_string(s.objects.size()) + " unique records",
            sw.seconds());

        report.check_predicate(
            "c0.no_record_is_silently_dropped",
            "rows in equals accepted plus rejected plus duplicates collapsed, and "
            "every rejection carries a reason",
            VerificationLevel::A,
            // BOTH ledgers: the delivered snapshot and the injected-duplicate
            // set. Checking only one would leave the delivered snapshot, which is
            // the artefact anyone actually consumes, unreconciled.
            snap50.ingestion.reconciles() && s.ingestion.reconciles() &&
                [&] {
                    for (const RejectionRecord& r : s.ingestion.rejected)
                        if (r.reason.empty()) return false;
                    for (const RejectionRecord& r : snap50.ingestion.rejected)
                        if (r.reason.empty()) return false;
                    return true;
                }(),
            "both ledgers balance and every rejection is explained",
            "delivered " + std::to_string(snap50.ingestion.rows_in) + " in = " +
                std::to_string(snap50.ingestion.accepted) + " + " +
                std::to_string(snap50.ingestion.rejected.size()) + " + " +
                std::to_string(snap50.ingestion.duplicates_collapsed) +
                "; with an injected duplicate " +
                std::to_string(s.ingestion.rows_in) + " in = " +
                std::to_string(s.ingestion.accepted) + " accepted + " +
                std::to_string(s.ingestion.rejected.size()) + " rejected + " +
                std::to_string(s.ingestion.duplicates_collapsed) + " collapsed",
            sw.seconds(),
            "this is the case that makes the acceptance criterion '0 silently "
            "dropped records' checkable rather than aspirational");
    }

    // =======================================================================
    // 0 unknown coordinate frames, 0 unit-ambiguous public fields.
    // =======================================================================
    {
        Stopwatch sw;
        const FrameDescription frames[5] = {
            ICRSDirection::description(), ICRSPhaseSpace::description(),
            BarycentricCartesian::description(),
            GalactocentricCartesian::description(), LocalShipFrame::description()};
        bool all_described = true;
        for (const FrameDescription& d : frames)
            if (!d.frame || !d.origin || !d.orientation || !d.epoch || !d.units ||
                !d.velocity_convention || std::string(d.frame).empty() ||
                std::string(d.units).empty())
                all_described = false;
        metrics.frames_described = 5;

        report.check_predicate(
            "c0.every_frame_declares_its_six_facts",
            "all five coordinate frames declare frame, origin, orientation, "
            "epoch, units and velocity convention",
            VerificationLevel::A, all_described,
            "five frames fully described",
            "ICRSDirection, ICRSPhaseSpace, BarycentricCartesian, "
            "GalactocentricCartesian, LocalShipFrame",
            sw.seconds(),
            "an anonymous Vec3 never crosses a frame boundary: each frame is a "
            "distinct type and every transformation is a named function between "
            "two named types, so the substitution is a compile error");

        // Unit-ambiguity: the strong types must not be constructible from a bare
        // double implicitly, and must not interconvert.
        const bool strong =
            !std::is_convertible<double, Parsec>::value &&
            !std::is_convertible<double, KmPerSecond>::value &&
            !std::is_convertible<Parsec, Metre>::value &&
            !std::is_convertible<MasPerYear, KmPerSecond>::value &&
            !std::is_convertible<Radian, Degree>::value;
        report.check_predicate(
            "c0.no_unit_ambiguous_public_field",
            "the unit types are not implicitly constructible from double and do "
            "not implicitly interconvert",
            VerificationLevel::A, strong,
            "all five conversions rejected at compile time",
            "checked with std::is_convertible on Parsec, KmPerSecond, Metre, "
            "MasPerYear and Degree",
            sw.seconds(),
            "milliarcseconds against arcseconds, parsecs against light-years, "
            "km/s against pc/yr: each has produced published errors, so each has "
            "its own type and the compiler refuses the mixture");
    }

    // =======================================================================
    // Literature sanity, at loose tolerance and clearly labelled.
    // =======================================================================
    {
        Stopwatch sw;
        double worst_rel = 0.0;
        for (const LiteratureStar& l : literature_reference()) {
            const double d = inverse_parallax_distance(Milliarcsecond(l.parallax_mas)).value;
            worst_rel = std::fmax(worst_rel,
                                  std::fabs(d - l.approximate_distance_pc) /
                                      l.approximate_distance_pc);
        }
        metrics.literature_sanity_error = worst_rel;
        report.check_residual(
            "c0.literature_reference_sanity",
            "inverse parallax reproduces the recalled approximate distances of "
            "three well-known nearby stars",
            VerificationLevel::D, "n/a",
            "3 stars, values RECALLED not retrieved, tolerance 1 per cent",
            worst_rel, 1e-2, sw.seconds(),
            "a LOOSE order-of-magnitude sanity check on the arithmetic, and "
            "nothing more. These parallaxes are recorded from memory and must be "
            "replaced by a real import before any result depending on them is "
            "quoted; no precision assertion rests on them");
    }

    // --- The manifest and results -------------------------------------------
    {
        std::ofstream os(out_dir + "/c0_manifest.txt", std::ios::binary);
        if (os) {
            os << "program_c C0 catalogue manifest\n\n";
            os << "DATA KIND: " << snap50.release.data_kind << "\n";
            os << "  This snapshot contains NO retrieved catalogue data. Every\n"
                  "  record is SYNTHETIC or LITERATURE_FROM_MEMORY. Run\n"
                  "  tools/import_gaia.py to produce a snapshot with genuine\n"
                  "  provenance.\n\n";
            os << "catalogue          " << snap50.release.catalogue << "\n";
            os << "release            " << snap50.release.release << "\n";
            os << "source table       " << snap50.release.source_table << "\n";
            os << "reference epoch    J" << snap50.release.reference_epoch_jyear << "\n";
            os << "volume radius      " << snap50.release.volume_radius_pc << " pc\n";
            os << "retrieval          " << snap50.release.retrieval_timestamp << "\n";
            os << "licence            " << snap50.release.licence << "\n";
            os << "checksum (FNV-1a)  " << snap50.checksum
               << "   integrity only, not cryptographic\n\n";
            os << "rows in            " << fx.report.rows_in << "\n";
            os << "accepted           " << fx.report.accepted << "\n";
            os << "rejected           " << fx.report.rejected.size() << "\n";
            for (const RejectionRecord& r : fx.report.rejected)
                os << "  reject " << r.source_id << ": " << r.reason << "\n";
            os << "distance estimator counts\n";
            os << "  InverseParallaxHighQuality " << fx.report.distance_inverse_parallax << "\n";
            os << "  ExternalPosterior          " << fx.report.distance_external << "\n";
            os << "  UserSupplied               " << fx.report.distance_user << "\n";
            os << "  None                       " << fx.report.distance_none << "\n";
            os << "missing radial velocity      " << fx.report.missing_radial_velocity << "\n";
            os << "missing covariance           " << fx.report.missing_covariance << "\n";
            os << "duplicates collapsed         " << fx.report.duplicates_collapsed
               << "\n";
            os << "ledger balances              "
               << (fx.report.reconciles() ? "yes" : "NO") << "   ("
               << fx.report.rows_in << " = " << fx.report.accepted << " + "
               << fx.report.rejected.size() << " + "
               << fx.report.duplicates_collapsed << ")\n";

            os << "\nuncertainty summary over the accepted records\n";
            {
                double min_sigma = 0.0, max_sigma = 0.0, sum_frac = 0.0;
                int n = 0;
                for (const StellarObject& s : snap50.objects) {
                    if (!s.has_defensible_distance()) continue;
                    const double sig = s.distance_uncertainty_pc;
                    if (n == 0 || sig < min_sigma) min_sigma = sig;
                    if (n == 0 || sig > max_sigma) max_sigma = sig;
                    if (s.distance_pc > 0.0) sum_frac += sig / s.distance_pc;
                    ++n;
                }
                os << "  records with a distance    " << n << "\n";
                os << "  distance sigma, min        " << min_sigma << " pc\n";
                os << "  distance sigma, max        " << max_sigma << " pc\n";
                os << "  mean fractional sigma      "
                   << (n > 0 ? sum_frac / n : 0.0) << "\n";
                os << "  parallax snr floor         " << kMinParallaxSignalToNoise
                   << " (below it the estimator is None, never a biased inverse)\n";
            }

            os << "\nrow counts by volume, from the offline query layer\n";
            os << "  within 5 pc                " << metrics.count_within_5pc << "\n";
            os << "  within 50 pc               " << metrics.count_within_50pc
               << "\n";
            os << "  confirmed exoplanet hosts  " << metrics.exoplanet_hosts << "\n";

            // Say plainly that this text was never sent anywhere. The query is a
            // committed SPECIFICATION of what a retrieval would ask for; calling
            // it "issued" would be a provenance claim this milestone cannot make.
            os << "\nADQL specifying the intended import. NOT ISSUED: this "
                  "milestone performed\nno retrieval. Reproduce it with "
                  "tools/import_gaia.py, whose adql_for_radius()\nemits this text "
                  "character for character.\n"
               << snap50.release.adql_query << "\n";
        }

        std::ofstream vs(out_dir + "/c0_validation_objects.csv", std::ios::binary);
        if (vs) {
            vs << "primary_name,source_id,data_kind,object_type,estimator,"
                  "distance_pc,distance_sigma_pc,x_pc,y_pc,z_pc,quality_flags,"
                  "planets\n";
            vs.setf(std::ios::scientific);
            vs.precision(10);
            for (std::size_t i = 0; i < snap50.objects.size(); ++i) {
                const StellarObject& s = snap50.objects[i];
                const Snapshot::Derived& d = snap50.derived[i];
                vs << s.primary_name << ',' << s.measured.gaia_source_id << ','
                   << s.provenance.source_catalogue << ',' << s.object_type << ','
                   << to_string(s.distance_estimator) << ',' << s.distance_pc << ','
                   << s.distance_uncertainty_pc << ',' << d.x_pc << ',' << d.y_pc
                   << ',' << d.z_pc << ',' << s.measured.quality_flags << ','
                   << s.planets.size() << '\n';
            }
        }
    }
}

}  // namespace checks
}  // namespace program_c
