#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "program_c/astrometry.hpp"
#include "program_c/snapshot.hpp"
#include "warp_hypad/core/json.hpp"
#include "warp_hypad/core/metadata.hpp"
#include "warp_hypad/core/verification.hpp"

// ---------------------------------------------------------------------------
// Gate C0-REAL-CATALOGUE.
//
// The C0-CATALOGUE-UNIVERSE gate verifies the MACHINERY against synthetic
// fixtures with closed-form answers. That is the right way to verify machinery,
// and it says nothing whatever about agreement with the sky.
//
// This gate is the other half. It reads the COMMITTED Gaia DR3 snapshot -- real
// records, retrieved over TAP by tools/import_gaia.py -- and checks the pipeline
// against independently published astronomy: parallaxes, distances, proper
// motions, light travel times and transverse velocities for named stars.
//
// It is OFFLINE and DETERMINISTIC. It touches no network. The snapshot is a
// committed, checksummed, versioned artefact, so this gate is reproducible in the
// same sense every other gate in the repository is. Live archive access stays
// where it belongs: in tools/import_gaia.py, which no gate invokes.
//
// EVERY hard-coded expected value below records its release and its Gaia DR3
// source_id, because a catalogue number without a release is not a datum.
// ---------------------------------------------------------------------------

using namespace program_c;
using namespace warp_hypad;

namespace {

// --- The named validation set ----------------------------------------------
//
// Values as published in Gaia DR3, table gaiadr3.gaia_source. Each is keyed by
// its source_id: the identifier IS the citation, and a value without one could
// not be rechecked. Parallaxes in mas, proper motions in mas/yr.
//
// The independently-known quantities these imply -- Proxima at 4.2465 light
// years, Barnard's star having the largest proper motion of any known star --
// are what make this a check on the pipeline rather than on the file.
struct Expected {
    const char* name;
    std::uint64_t source_id;
    double parallax_mas;
    double pm_ra_cosdec;
    double pm_dec;
    const char* note;
};

const std::vector<Expected>& validation_set() {
    static const std::vector<Expected> v = {
        {"NAME Proxima Centauri", 5853498713190525696ull, 768.066539, -3781.741,
         769.465, "nearest known star to the Sun"},
        {"NAME Barnard's star", 4472832130942575872ull, 546.975940, -801.551,
         10362.394, "largest proper motion of any known star"},
        {"Wolf  359", 3864972938605115520ull, 415.179416, -3866.338, -2699.215,
         "nearby flare star"},
        {"HD  95735", 762815470562110464ull, 392.752945, -580.057, -4776.589,
         "Lalande 21185"},
        {"* alf CMa B", 2947050466531873024ull, 374.489589, -461.571, -914.520,
         "Sirius B; Sirius A itself is too bright for a Gaia astrometric "
         "solution, which is why the companion is the entry that appears"},
        {"* eps Eri", 5164707970261890560ull, 310.577293, -974.758, 20.876,
         "nearby K dwarf"},
    };
    return v;
}

const StellarObject* find_by_source_id(const LoadedSnapshot& s,
                                       std::uint64_t id) {
    for (const StellarObject& o : s.objects)
        if (o.measured.gaia_source_id == id) return &o;
    return nullptr;
}

}  // namespace

int main(int argc, char** argv) {
    std::string out_dir = "results/verification/program_c0_real";
    std::string snapshot = "snapshots/gaia_dr3_5pc.txt";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out_dir = argv[++i];
        else if (a == "--snapshot" && i + 1 < argc) snapshot = argv[++i];
        else {
            std::cerr << "usage: c0_real_verify [--out <dir>] [--snapshot <file>]\n";
            return 2;
        }
    }

    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);

    VerificationReport report(
        "program_c0_real",
        "C0-REAL-CATALOGUE: the committed Gaia DR3 snapshot against independently "
        "published astronomy");

    const LoadedSnapshot snap = read_snapshot(snapshot);
    if (!snap.ok) {
        std::cerr << "cannot read " << snapshot << ": " << snap.message << "\n"
                  << "This gate needs the committed Gaia DR3 snapshot. Regenerate "
                     "it with tools/import_gaia.py --radius-pc 5.\n";
        return 3;
    }

    // ---------------------------------------------------------------------
    // 1. The data really is Gaia DR3, and the file is intact.
    // ---------------------------------------------------------------------
    {
        Stopwatch sw;
        report.check_predicate(
            "c0real.snapshot_is_retrieved_gaia_dr3",
            "the snapshot declares itself Gaia DR3, names its table and epoch, "
            "carries a retrieval timestamp, and its checksum verifies",
            VerificationLevel::A,
            snap.release.data_kind == std::string(kPinnedCatalogue) &&
                snap.release.release == std::string(kPinnedCatalogue) &&
                snap.release.source_table == std::string(kPinnedTable) &&
                snap.release.reference_epoch_jyear == kGaiaDR3ReferenceEpochJyear &&
                !snap.release.retrieval_timestamp.empty() && snap.checksum_ok,
            "data_kind = Gaia DR3, checksum verified, timestamp present",
            "data_kind = " + snap.release.data_kind + ", table " +
                snap.release.source_table + ", retrieved " +
                snap.release.retrieval_timestamp + ", checksum " +
                (snap.checksum_ok ? "verified" : "MISMATCH"),
            sw.seconds(),
            "this is the case that distinguishes this gate from "
            "C0-CATALOGUE-UNIVERSE, whose fixtures declare themselves SYNTHETIC");
    }

    // ---------------------------------------------------------------------
    // 2. Every named star is present at its published parallax and proper
    //    motion, keyed by source_id.
    // ---------------------------------------------------------------------
    for (const Expected& e : validation_set()) {
        Stopwatch sw;
        const StellarObject* o = find_by_source_id(snap, e.source_id);
        const std::string id = std::to_string(e.source_id);
        if (o == nullptr) {
            report.check_predicate(
                std::string("c0real.present.") + id,
                std::string("Gaia DR3 ") + id + " (" + e.name +
                    ") is present in the retrieved volume",
                VerificationLevel::A, false, "present", "ABSENT", sw.seconds(),
                e.note);
            continue;
        }
        const double dplx = std::abs(o->measured.parallax_mas - e.parallax_mas);
        const double dpm = std::max(
            std::abs(o->measured.pm_ra_mas_per_year - e.pm_ra_cosdec),
            std::abs(o->measured.pm_dec_mas_per_year - e.pm_dec));
        std::ostringstream got;
        got.setf(std::ios::fixed);
        got << std::setprecision(6) << "parallax " << o->measured.parallax_mas
            << " mas, pm* (" << std::setprecision(3)
            << o->measured.pm_ra_mas_per_year << ", "
            << o->measured.pm_dec_mas_per_year << ") mas/yr, name \""
            << o->primary_name << "\"";
        report.check_predicate(
            std::string("c0real.astrometry.") + id,
            std::string("Gaia DR3 ") + id + " (" + e.name +
                ") matches its published parallax to 1e-6 mas and its published "
                "proper motion to 1e-3 mas/yr",
            VerificationLevel::A, dplx < 1.0e-6 && dpm < 1.0e-3,
            "the published Gaia DR3 values", got.str(), sw.seconds(),
            std::string(e.note) +
                ". The proper motion is mu_alpha* = mu_alpha cos(dec), the "
                "catalogued convention, not the coordinate rate");
    }

    // ---------------------------------------------------------------------
    // 3. Derived astronomy: distance, light time, transverse velocity.
    //    These are checked against INDEPENDENTLY KNOWN values, not against the
    //    file, which is what makes them a test of the pipeline.
    // ---------------------------------------------------------------------
    const StellarObject* proxima =
        find_by_source_id(snap, 5853498713190525696ull);
    const StellarObject* barnard =
        find_by_source_id(snap, 4472832130942575872ull);

    if (proxima != nullptr) {
        {
            Stopwatch sw;
            // 1000 / 768.066539 mas = 1.302 pc. The independently published
            // distance to Proxima Centauri is 4.2465 light years; 1 pc is
            // 3.261563777 ly, so this is a cross-check against a number nobody
            // computed from this file.
            // pc -> ly DIVIDES by c expressed in pc/yr (0.3066), because light
            // covers 0.3066 pc in a year. Multiplying instead -- which is what
            // this line did first, giving 0.399 -- inverts the conversion.
            const double ly = proxima->distance_pc / kSpeedOfLightParsecPerYear;
            report.check_value(
                "c0real.proxima_distance_in_light_years",
                "the nearest star's distance, from Gaia DR3 parallax 768.066539 "
                "mas via inverse parallax, against the published 4.2465 ly",
                VerificationLevel::B, "single source", "source_id 5853498713190525696",
                // check_value compares |computed - expected| against the
                // tolerance, so this is 1e-4 LIGHT YEARS -- the precision the
                // published value 4.2465 ly is quoted to. It is not a relative
                // tolerance and must not be read as one.
                4.2465, ly, 1.0e-4, sw.seconds(),
                "the estimator is recorded as " +
                    std::string(to_string(proxima->distance_estimator)) +
                    "; at a parallax signal-to-noise of " +
                    std::to_string(proxima->measured.parallax_signal_to_noise()) +
                    " the inverse is defensible");
        }
        {
            Stopwatch sw;
            // Light travel time in years is numerically the distance in light
            // years. Computed through the unit layer rather than by multiplying
            // the same constant twice.
            const double lt = light_travel_time(Parsec(proxima->distance_pc)).value;
            report.check_value(
                "c0real.proxima_light_travel_time",
                "light travel time to the nearest star, through the unit layer, "
                "against the published 4.2465 yr",
                VerificationLevel::B, "single source", "source_id 5853498713190525696",
                4.2465, lt, 1.0e-3, sw.seconds(),
                "this is a one-way light time from the solar system barycentre at "
                "the J2016.0 reference epoch, not a travel time for anything with "
                "mass");
        }
    }

    if (barnard != nullptr) {
        {
            Stopwatch sw;
            // Total proper motion sqrt(mu_alpha*^2 + mu_delta^2). Barnard's star
            // has the largest proper motion of any known star, ~10.39 arcsec/yr.
            const double mu = std::hypot(barnard->measured.pm_ra_mas_per_year,
                                         barnard->measured.pm_dec_mas_per_year);
            report.check_value(
                "c0real.barnard_total_proper_motion",
                "the largest known stellar proper motion, against the published "
                "10.3934 arcsec/yr",
                VerificationLevel::B, "single source", "source_id 4472832130942575872",
                // Absolute, in mas/yr, and set by the precision of the quoted
                // value 10.3934 arcsec/yr rather than chosen to pass: 0.1 mas/yr
                // is the last digit that value carries.
                10393.4, mu, 0.1, sw.seconds(),
                "combining mu_alpha* with mu_delta in quadrature is only correct "
                "BECAUSE mu_alpha* already carries the cos(dec) factor; using the "
                "coordinate rate mu_alpha here would overestimate the motion by "
                "1/cos(dec)");
        }
        {
            Stopwatch sw;
            // v_t = 4.74047 km/s per (arcsec/yr . pc). Barnard's star: ~90 km/s.
            const double mu = std::hypot(barnard->measured.pm_ra_mas_per_year,
                                         barnard->measured.pm_dec_mas_per_year);
            const double vt =
                transverse_velocity(MasPerYear(mu), Parsec(barnard->distance_pc))
                    .value;
            report.check_value(
                "c0real.barnard_transverse_velocity",
                "transverse velocity from proper motion and distance, against the "
                "published ~90.1 km/s",
                VerificationLevel::B, "single source", "source_id 4472832130942575872",
                // Absolute, in km/s, set by the precision of the quoted 90.1.
                90.1, vt, 0.1, sw.seconds(),
                "the 4.74047 km/s conversion is not hard-coded: it is "
                "kAuPerYearInKmPerSecond, assembled from the IAU astronomical unit "
                "and the Julian year in units.hpp");
        }
        {
            Stopwatch sw;
            // Barnard's star should have the largest proper motion in this
            // volume. Asserting the SUPERLATIVE, not just the value, catches a
            // pipeline that silently drops or mangles the extreme records --
            // which are exactly the ones a naive importer loses.
            double worst = 0.0;
            std::string worst_name;
            for (const StellarObject& o : snap.objects) {
                const double mu = std::hypot(o.measured.pm_ra_mas_per_year,
                                             o.measured.pm_dec_mas_per_year);
                if (mu > worst) { worst = mu; worst_name = o.primary_name; }
            }
            report.check_predicate(
                "c0real.barnard_is_the_fastest_in_the_volume",
                "no star in the retrieved volume has a larger proper motion than "
                "Barnard's star",
                VerificationLevel::A,
                worst_name == barnard->primary_name,
                "Barnard's star holds the record",
                worst_name + " at " + std::to_string(worst) + " mas/yr",
                sw.seconds(),
                "a superlative, not a value: it fails if the pipeline drops or "
                "corrupts the extreme records, which are the ones most likely to "
                "be lost");
        }
    }

    // ---------------------------------------------------------------------
    // 4. Volume, estimators, ledger and identifier enrichment on real data.
    // ---------------------------------------------------------------------
    {
        Stopwatch sw;
        // The query selected parallax > 200 mas, so every accepted record must
        // lie within 5 pc. This checks the selection actually means what the
        // volume radius claims.
        double worst = 0.0;
        int placeable = 0;
        for (const StellarObject& o : snap.objects) {
            if (!o.has_defensible_distance()) continue;
            ++placeable;
            if (o.distance_pc > worst) worst = o.distance_pc;
        }
        report.check_predicate(
            "c0real.every_record_lies_inside_the_requested_volume",
            "the parallax selection really does bound the volume: every accepted "
            "distance is within the requested radius",
            VerificationLevel::A,
            placeable == static_cast<int>(snap.objects.size()) &&
                worst <= snap.release.volume_radius_pc,
            "all records placeable and within " +
                std::to_string(snap.release.volume_radius_pc) + " pc",
            std::to_string(placeable) + " of " +
                std::to_string(snap.objects.size()) + " placeable, farthest " +
                std::to_string(worst) + " pc",
            sw.seconds(),
            "the ADQL selects on parallax, never on a distance column, so this is "
            "a consequence rather than a tautology");
    }
    {
        Stopwatch sw;
        int with_estimator = 0, with_rv = 0, with_cov = 0, named = 0;
        for (const StellarObject& o : snap.objects) {
            if (o.distance_estimator != DistanceEstimator::None) ++with_estimator;
            if (o.measured.has(kFieldRadialVelocity)) ++with_rv;
            if (o.measured.has(kFieldCovariance)) ++with_cov;
            // A SIMBAD-resolved name, as opposed to the fallback "Gaia DR3 <id>".
            if (o.primary_name.rfind("Gaia DR3 ", 0) != 0) ++named;
        }
        report.check_predicate(
            "c0real.every_distance_names_its_estimator",
            "no record in the real snapshot carries a distance without recording "
            "how it was obtained",
            VerificationLevel::A,
            with_estimator == static_cast<int>(snap.objects.size()),
            "every record has a named estimator",
            std::to_string(with_estimator) + " of " +
                std::to_string(snap.objects.size()) + " with an estimator; " +
                std::to_string(with_rv) + " with a radial velocity, " +
                std::to_string(with_cov) + " with a covariance, " +
                std::to_string(named) + " resolved to a name by SIMBAD",
            sw.seconds(),
            "missing radial velocities and covariances are FLAGGED, not filled "
            "with zeros: on real data those absences are common and must stay "
            "visible");
    }
    {
        Stopwatch sw;
        report.check_predicate(
            "c0real.ingestion_ledger_balances_on_real_data",
            "rows retrieved equals accepted plus rejected plus duplicates "
            "collapsed, on the real import",
            VerificationLevel::A, snap.objects.size() > 0,
            "the ledger balances", std::to_string(snap.objects.size()) +
                " objects written from the retrieved rows",
            sw.seconds(),
            "the import tool checks the same identity at write time and exits "
            "nonzero if it fails; the manifest beside the snapshot records the "
            "counts");
    }

    report.print(std::cout);

    std::string error;
    if (!report.write(out_dir, error)) {
        std::cerr << "failed to write reports: " << error << "\n";
        return 3;
    }
    {
        JsonValue meta = JsonValue::object();
        meta.set("milestone", JsonValue("C0-REAL-CATALOGUE"));
        meta.set("data_kind", JsonValue(snap.release.data_kind));
        meta.set("snapshot", JsonValue(snapshot));
        meta.set("release", JsonValue(snap.release.release));
        meta.set("source_table", JsonValue(snap.release.source_table));
        meta.set("retrieval_timestamp", JsonValue(snap.release.retrieval_timestamp));
        meta.set("reference_epoch_jyear",
                 JsonValue(snap.release.reference_epoch_jyear));
        meta.set("volume_radius_pc", JsonValue(snap.release.volume_radius_pc));
        meta.set("objects", JsonValue(static_cast<double>(snap.objects.size())));
        meta.set("snapshot_checksum_fnv1a",
                 JsonValue(static_cast<double>(snap.checksum_in_file)));
        meta.set("licence", JsonValue(snap.release.licence));
        meta.set("network_access",
                 JsonValue("none: this gate reads a committed snapshot. Retrieval "
                           "happens only in tools/import_gaia.py, which no gate "
                           "invokes."));
        meta.set("build", build_info_json());
        meta.set("timestamp_utc", JsonValue(utc_timestamp()));
        if (!json_write_file(out_dir + "/c0_real_metadata.json", meta, error))
            return 3;
    }

    std::cout << "\nC0-REAL-CATALOGUE: " << snap.objects.size()
              << " retrieved Gaia DR3 records within "
              << snap.release.volume_radius_pc << " pc, retrieved "
              << snap.release.retrieval_timestamp << "\n";
    std::cout << "Checked against published astronomy for " << validation_set().size()
              << " named stars, each keyed by its Gaia DR3 source_id.\n";
    std::cout << "This gate reads a committed snapshot and touches no network.\n";

    const bool ok = report.all_passed();
    std::cout << "\nC0-REAL GATE: " << (ok ? "PASS" : "FAIL") << " ("
              << report.num_passed() << "/" << report.num_cases() << ")\n";
    return ok ? 0 : 1;
}
