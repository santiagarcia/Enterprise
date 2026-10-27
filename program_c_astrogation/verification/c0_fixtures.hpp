#pragma once

#include "program_c/snapshot.hpp"

// ---------------------------------------------------------------------------
// Recorded fixtures.
//
// READ THIS BEFORE TRUSTING ANY NUMBER IN THIS FILE.
//
// This milestone was implemented WITHOUT NETWORK ACCESS, so no Gaia DR3 record
// was retrieved. Nothing here is presented as retrieved catalogue data. There are
// exactly two kinds of record and they are tagged apart:
//
//   SYNTHETIC              constructed so the correct answer is known in closed
//                          form. These drive every numerical assertion in the
//                          gate: round trips, propagation, covariance, snapshot
//                          determinism. Their values are chosen for arithmetic
//                          convenience and correspond to no real star.
//
//   LITERATURE_FROM_MEMORY approximate, well-known values recalled rather than
//                          retrieved. Used ONLY for order-of-magnitude sanity
//                          checks at loose tolerance. They must be replaced by a
//                          real import before any result depending on them is
//                          quoted.
//
// No fixture carries a Gaia source_id, release string or retrieval timestamp it
// does not have. Writing plausible numbers under a real catalogue's name would be
// fabricating records, and the gate would then be verifying a fiction.
//
// tools/import_gaia.py is the live client. Running it is what produces fixtures
// with genuine provenance, and it is a separate optional step.
// ---------------------------------------------------------------------------

namespace program_c {
namespace fixtures {

inline constexpr const char* kSyntheticTag = "SYNTHETIC";
inline constexpr const char* kMemoryTag = "LITERATURE_FROM_MEMORY";

// The ADQL the live import issues. This is a SPECIFICATION, authored here, not
// retrieved data -- so it is legitimate to commit and to assert against. Nothing
// in C0 executes it: the text records what a retrieval WOULD ask for.
//
// THIS TEXT AND adql_for_radius() IN tools/import_gaia.py MUST BE IDENTICAL,
// character for character, including the fixed six-decimal formatting of both
// numbers. They diverged once -- the C++ text omitted the proper-motion
// correlation columns that the importer needs to assemble a covariance, so a
// snapshot built from the committed query text would have carried
// kQualityMissingCovariance on every record for no physical reason. Check them
// against each other with
//
//     python tools/import_gaia.py --print-query --radius-pc 50
//
// which contacts no network.
//
// Selection is on PARALLAX, never on a distance column: the catalogue measures
// parallax, and a distance cut would bury an estimator choice in the selection.
inline std::string adql_for_radius(double radius_pc, double min_parallax_snr) {
    const double min_parallax_mas = 1000.0 / radius_pc;
    std::ostringstream q;
    q.setf(std::ios::fixed);
    q.precision(6);
    q << "SELECT source_id, ra, ra_error, dec, dec_error, "
         "parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error, "
         "radial_velocity, radial_velocity_error, ref_epoch, "
         "ra_dec_corr, ra_parallax_corr, ra_pmra_corr, ra_pmdec_corr, "
         "dec_parallax_corr, dec_pmra_corr, dec_pmdec_corr, "
         "parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, "
         "duplicated_source, phot_g_mean_mag, ruwe, astrometric_params_solved "
      << "FROM " << kPinnedTable << " "
      << "WHERE parallax IS NOT NULL "
      << "AND parallax > " << min_parallax_mas << " "
      << "AND parallax_error > 0 "
      << "AND parallax / parallax_error >= " << min_parallax_snr << " "
      << "ORDER BY parallax DESC";
    return q.str();
}

inline Provenance synthetic_provenance(const std::string& id) {
    Provenance p;
    p.source_catalogue = kSyntheticTag;
    p.data_release = "n/a: constructed for verification, not retrieved";
    p.source_table = "n/a";
    p.query = "n/a: no query was issued";
    p.retrieval_timestamp = "n/a";
    p.licence = "n/a: not third-party data";
    p.original_identifier = id;
    return p;
}

inline Provenance memory_provenance(const std::string& id) {
    Provenance p;
    p.source_catalogue = kMemoryTag;
    p.data_release =
        "UNVERIFIED: recalled, not retrieved. Replace before quoting.";
    p.source_table = "n/a";
    p.query = "n/a: no query was issued";
    p.retrieval_timestamp = "n/a";
    p.licence = "n/a";
    p.original_identifier = id;
    return p;
}

// Build one synthetic source. Angles in degrees for readability at the call site;
// converted immediately.
inline StellarObject make_synthetic(std::uint64_t id, const std::string& name,
                                   double ra_deg, double dec_deg,
                                   double parallax_mas, double parallax_err_mas,
                                   double pmra_cosdec, double pmdec, double rv,
                                   bool have_rv = true, bool have_cov = false) {
    StellarObject s;
    CatalogueAstrometry& a = s.measured;
    a.gaia_source_id = id;
    a.ra_rad = to_radian(Degree(ra_deg)).value;
    a.dec_rad = to_radian(Degree(dec_deg)).value;
    a.parallax_mas = parallax_mas;
    a.parallax_error_mas = parallax_err_mas;
    a.pm_ra_mas_per_year = pmra_cosdec;
    a.pm_dec_mas_per_year = pmdec;
    a.radial_velocity_km_per_s = have_rv ? rv : 0.0;
    a.reference_epoch_jyear = kGaiaDR3ReferenceEpochJyear;

    a.set_valid(kFieldRa);
    a.set_valid(kFieldDec);
    a.set_valid(kFieldParallax);
    a.set_valid(kFieldPmRa);
    a.set_valid(kFieldPmDec);
    if (have_rv) a.set_valid(kFieldRadialVelocity);
    else a.flag(kQualityMissingRadialVelocity);

    if (have_cov) {
        a.set_valid(kFieldCovariance);
        // A deliberately simple diagonal covariance in the catalogue's own mixed
        // units: mas, mas, mas, mas/yr, mas/yr, km/s.
        const double sig[6] = {to_radian(Milliarcsecond(0.05)).value,
                               to_radian(Milliarcsecond(0.04)).value,
                               parallax_err_mas, 0.06, 0.05, 0.4};
        for (int i = 0; i < 6; ++i) a.covariance(i, i) = sig[i] * sig[i];
    } else {
        a.flag(kQualityMissingCovariance);
    }

    const DistanceDecision d = decide_distance(a);
    s.distance_estimator = d.estimator;
    s.distance_pc = d.distance_pc;
    s.distance_uncertainty_pc = d.uncertainty_pc;
    if (d.estimator != DistanceEstimator::None) a.set_valid(kFieldDistance);
    else if (!(a.parallax_mas > 0.0)) a.flag(kQualityParallaxNonPositive);
    else a.flag(kQualityParallaxLowSignalToNoise);

    s.primary_name = name;
    s.identifiers.push_back(name);
    s.object_type = "synthetic test source";
    s.provenance = synthetic_provenance(name);
    return s;
}

// --- The validation set -----------------------------------------------------
//
// One synthetic record per required category, each with a closed-form answer.

struct FixtureSet {
    std::vector<StellarObject> objects;
    IngestionReport report;
};

inline FixtureSet synthetic_catalogue() {
    FixtureSet f;
    auto add = [&](const StellarObject& s) { f.objects.push_back(s); };

    // 1. Exactly on the +x axis at exactly 10 pc: position (10, 0, 0).
    add(make_synthetic(1001, "SYN-AXIS-10PC", 0.0, 0.0, 100.0, 0.5,
                       0.0, 0.0, 0.0, true, true));
    // 2. On the celestial pole at exactly 5 pc: position (0, 0, 5).
    add(make_synthetic(1002, "SYN-POLE-5PC", 0.0, 90.0, 200.0, 1.0,
                       0.0, 0.0, 0.0, true, true));
    // 3. High declination with pure mu_alpha*: the cos(dec) convention test.
    add(make_synthetic(1003, "SYN-HIGHDEC-CONVENTION", 30.0, 60.0, 50.0, 0.2,
                       500.0, 0.0, 0.0, true, true));
    // 4. High proper motion, for propagation.
    add(make_synthetic(1004, "SYN-HIGH-PM", 120.0, -20.0, 250.0, 0.5,
                       4000.0, -3000.0, 0.0, true, true));
    // 5. Pure radial velocity: propagation must be purely radial.
    add(make_synthetic(1005, "SYN-RADIAL-ONLY", 200.0, 10.0, 125.0, 0.4,
                       0.0, 0.0, 50.0, true, true));
    // 6. Stationary: zero proper motion and zero radial velocity.
    add(make_synthetic(1006, "SYN-STATIONARY", 45.0, 45.0, 40.0, 0.1,
                       0.0, 0.0, 0.0, true, true));
    // 7. Missing radial velocity.
    add(make_synthetic(1007, "SYN-NO-RV", 88.0, -5.0, 66.6667, 0.3,
                       120.0, 80.0, 0.0, false, true));
    // 8. Poor parallax: MUST be rejected for inverse-parallax distance.
    add(make_synthetic(1008, "SYN-POOR-PARALLAX", 15.0, 15.0, 4.0, 2.0,
                       10.0, 10.0, 0.0, true, false));
    // 9. Negative parallax: also rejected, for a different reason.
    add(make_synthetic(1009, "SYN-NEGATIVE-PARALLAX", 300.0, -60.0, -1.5, 0.6,
                       5.0, 5.0, 0.0, true, false));
    // 10. Exoplanet host.
    {
        StellarObject s = make_synthetic(1010, "SYN-EXOHOST", 270.0, 30.0, 90.9091,
                                        0.3, 200.0, -150.0, -12.0, true, true);
        s.object_type = "synthetic exoplanet host";
        PlanetRecord p;
        p.name = "SYN-EXOHOST b";
        p.orbital_period_days = 11.186;
        p.semi_major_axis_au = 0.0485;
        p.radius_earth = 1.1;
        p.mass_earth = 1.27;
        p.provenance = synthetic_provenance("SYN-EXOHOST b");
        s.planets.push_back(p);
        s.identifiers.push_back("SYN-ALIAS-EXOHOST");
        add(s);
    }
    // 11. Binary / multiple.
    {
        StellarObject s = make_synthetic(1011, "SYN-BINARY-A", 100.0, -35.0,
                                        160.0, 0.4, 300.0, 200.0, 20.0, true, true);
        s.object_type = "synthetic binary component";
        s.identifiers.push_back("SYN-BINARY-A-ALT");
        add(s);
        StellarObject b = make_synthetic(1012, "SYN-BINARY-B", 100.001, -35.001,
                                        160.1, 0.4, 301.0, 199.0, 20.5, true, true);
        b.object_type = "synthetic binary component";
        add(b);
    }
    // 12. Outside the default 50 pc volume: retained and flagged.
    {
        StellarObject s = make_synthetic(1013, "SYN-OUTSIDE-50PC", 10.0, 10.0,
                                        16.6667, 0.1, 10.0, 10.0, 0.0, true, true);
        s.measured.flag(kQualityOutsideRequestedVolume);
        add(s);
    }

    f.report.rows_in = static_cast<int>(f.objects.size());
    for (const StellarObject& s : f.objects) {
        if (s.distance_estimator == DistanceEstimator::None) {
            const DistanceDecision d = decide_distance(s.measured);
            f.report.rejected.push_back({s.measured.gaia_source_id, d.reason});
            ++f.report.distance_none;
        } else {
            ++f.report.accepted;
            ++f.report.distance_inverse_parallax;
        }
        if (s.measured.flagged(kQualityMissingRadialVelocity))
            ++f.report.missing_radial_velocity;
        if (s.measured.flagged(kQualityMissingCovariance))
            ++f.report.missing_covariance;
    }
    return f;
}

// --- Literature reference values, for loose sanity checks only ---------------

struct LiteratureStar {
    const char* name;
    double parallax_mas;
    double total_pm_mas_per_year;
    double approximate_distance_pc;
};

// RECALLED, NOT RETRIEVED. Quoted to about three significant figures at best,
// and used only where a several-per-cent tolerance is appropriate. Every one must
// be replaced by a real import before being cited.
inline std::vector<LiteratureStar> literature_reference() {
    return {
        {"Proxima Centauri (recalled)", 768.5, 3853.0, 1.301},
        {"Barnard's Star (recalled)", 547.0, 10358.0, 1.828},
        {"Sirius (recalled)", 379.0, 1339.0, 2.638},
    };
}

// --- Snapshot assembly ------------------------------------------------------

inline Snapshot build_snapshot(double radius_pc) {
    Snapshot snap;
    const FixtureSet f = synthetic_catalogue();
    snap.objects = f.objects;
    snap.ingestion = f.report;

    snap.release.catalogue = kPinnedCatalogue;
    snap.release.release = kPinnedCatalogue;
    snap.release.source_table = kPinnedTable;
    snap.release.reference_epoch_jyear = kGaiaDR3ReferenceEpochJyear;
    snap.release.volume_radius_pc = radius_pc;
    snap.release.retrieval_timestamp =
        "n/a: no retrieval was performed; fixtures are SYNTHETIC";
    snap.release.adql_query = adql_for_radius(radius_pc, kMinParallaxSignalToNoise);
    snap.release.licence =
        "Gaia data would be ESA/Gaia/DPAC, CC BY-NC 3.0 IGO; NOT APPLICABLE here "
        "because no Gaia data is present";
    snap.release.data_kind = kSyntheticTag;

    finalise(snap);
    return snap;
}

}  // namespace fixtures
}  // namespace program_c
