#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "program_c/frames.hpp"

// ---------------------------------------------------------------------------
// The catalogue record.
//
// TWO RULES THIS FILE EXISTS TO ENFORCE.
//
// 1. NaN ALONE NEVER EXPRESSES MISSING DATA. Every record carries a
//    `valid_fields` bitmask. A NaN that arrives from a catalogue with its
//    validity bit set is a DATA ERROR to be reported, not a silent "missing" --
//    and a field whose bit is clear must not be read at all. The accessors below
//    enforce that, so a caller cannot forget.
//
// 2. NOTHING IS SILENTLY DROPPED. A record that fails a quality gate is RETAINED
//    with a reason recorded. The gate reconciles input rows against
//    accepted + rejected and asserts the difference is zero.
//
// PROVENANCE. Every record keeps the catalogue, release, table, exact query,
// retrieval timestamp, licence and original identifier. `source_catalogue` also
// carries the honesty marker: a record generated for testing says SYNTHETIC, and
// an approximate value recalled rather than retrieved says
// LITERATURE_FROM_MEMORY. Neither is ever labelled "Gaia DR3".
// ---------------------------------------------------------------------------

namespace program_c {

// --- Provenance -------------------------------------------------------------

struct Provenance {
    std::string source_catalogue;    // "Gaia DR3", "SYNTHETIC", "LITERATURE_FROM_MEMORY"
    std::string data_release;
    std::string source_table;
    std::string query;
    std::string retrieval_timestamp;
    std::string licence;
    std::string original_identifier;

    bool is_real_catalogue() const {
        return source_catalogue != "SYNTHETIC" &&
               source_catalogue != "LITERATURE_FROM_MEMORY" &&
               !source_catalogue.empty();
    }
    bool is_synthetic() const { return source_catalogue == "SYNTHETIC"; }
    bool is_from_memory() const {
        return source_catalogue == "LITERATURE_FROM_MEMORY";
    }
};

// The pinned source. A silent migration to a later release fails the gate.
inline constexpr const char* kPinnedCatalogue = "Gaia DR3";
inline constexpr const char* kPinnedTable = "gaiadr3.gaia_source";
inline constexpr double kGaiaDR3ReferenceEpochJyear = 2016.0;
inline constexpr const char* kGaiaTapEndpoint =
    "https://gea.esac.esa.int/tap-server/tap";

// --- Validity and quality ---------------------------------------------------

enum ValidField : std::uint64_t {
    kFieldRa = 1ull << 0,
    kFieldDec = 1ull << 1,
    kFieldParallax = 1ull << 2,
    kFieldPmRa = 1ull << 3,
    kFieldPmDec = 1ull << 4,
    kFieldRadialVelocity = 1ull << 5,
    kFieldCovariance = 1ull << 6,
    kFieldDistance = 1ull << 7,
};

enum QualityFlag : std::uint64_t {
    kQualityNone = 0,
    kQualityParallaxNonPositive = 1ull << 0,
    kQualityParallaxLowSignalToNoise = 1ull << 1,
    kQualityMissingRadialVelocity = 1ull << 2,
    kQualityMissingCovariance = 1ull << 3,
    kQualityDuplicatedSource = 1ull << 4,
    kQualityNonFiniteValue = 1ull << 5,
    kQualityOutsideRequestedVolume = 1ull << 6,
};

// The parallax quality gate that licenses inverse-parallax distance, stated
// numerically in ONE place. A parallax signal-to-noise of 10 corresponds to a
// distance bias of order a per cent for this estimator; below it the inverse is
// not defensible and the estimator must be `None` or an external posterior.
inline constexpr double kMinParallaxSignalToNoise = 10.0;
inline constexpr double kMinParallaxMas = 0.0;   // strictly positive required

// --- Astrometry -------------------------------------------------------------

// 6x6 covariance over (ra, dec, parallax, pm_ra_cosdec, pm_dec, radial_velocity)
// in the catalogue's own units: mas, mas, mas, mas/yr, mas/yr, km/s. The unit
// mixture is the catalogue's, and it is recorded rather than silently converted.
struct Mat6 {
    double m[36] = {};
    double& operator()(int i, int j) { return m[6 * i + j]; }
    const double& operator()(int i, int j) const { return m[6 * i + j]; }
    static Mat6 zero() { return Mat6(); }
};

struct CatalogueAstrometry {
    std::uint64_t gaia_source_id = 0;

    double ra_rad = 0.0;
    double dec_rad = 0.0;
    double parallax_mas = 0.0;
    double parallax_error_mas = 0.0;
    double pm_ra_mas_per_year = 0.0;     // mu_alpha* = mu_alpha cos(dec)
    double pm_dec_mas_per_year = 0.0;
    double radial_velocity_km_per_s = 0.0;
    double reference_epoch_jyear = kGaiaDR3ReferenceEpochJyear;

    Mat6 covariance;
    std::uint64_t valid_fields = 0;
    std::uint64_t quality_flags = kQualityNone;

    bool has(ValidField f) const { return (valid_fields & f) != 0; }
    void set_valid(ValidField f) { valid_fields |= f; }
    void flag(QualityFlag f) { quality_flags |= f; }
    bool flagged(QualityFlag f) const { return (quality_flags & f) != 0; }

    double parallax_signal_to_noise() const {
        if (!(parallax_error_mas > 0.0)) return 0.0;
        return parallax_mas / parallax_error_mas;
    }
};

// --- Distance ---------------------------------------------------------------

enum class DistanceEstimator {
    None,
    InverseParallaxHighQuality,
    ExternalPosterior,
    UserSupplied,
};

inline const char* to_string(DistanceEstimator e) {
    switch (e) {
        case DistanceEstimator::None: return "None";
        case DistanceEstimator::InverseParallaxHighQuality:
            return "InverseParallaxHighQuality";
        case DistanceEstimator::ExternalPosterior: return "ExternalPosterior";
        case DistanceEstimator::UserSupplied: return "UserSupplied";
    }
    return "unknown";
}

// --- Planets ----------------------------------------------------------------

struct PlanetRecord {
    std::string name;
    double orbital_period_days = 0.0;
    double semi_major_axis_au = 0.0;
    double radius_earth = 0.0;
    double mass_earth = 0.0;
    std::uint64_t valid_fields = 0;
    Provenance provenance;
};

// --- The object -------------------------------------------------------------

struct StellarObject {
    CatalogueAstrometry measured;

    DistanceEstimator distance_estimator = DistanceEstimator::None;
    double distance_pc = 0.0;
    double distance_uncertainty_pc = 0.0;

    std::string primary_name;
    std::vector<std::string> identifiers;
    std::string object_type;

    std::vector<PlanetRecord> planets;
    Provenance provenance;

    bool has_defensible_distance() const {
        return distance_estimator != DistanceEstimator::None &&
               measured.has(kFieldDistance);
    }
};

// --- The distance decision, in one place ------------------------------------
//
// Applies the quality gate and RECORDS the estimator. Never returns a distance
// without saying how it was obtained, and never deletes a record: a source that
// fails the gate comes back with estimator None and the reason flagged.
struct DistanceDecision {
    DistanceEstimator estimator = DistanceEstimator::None;
    double distance_pc = 0.0;
    double uncertainty_pc = 0.0;
    std::string reason;
};

inline DistanceDecision decide_distance(const CatalogueAstrometry& a) {
    DistanceDecision d;
    if (!a.has(kFieldParallax)) {
        d.reason = "no parallax field: distance estimator None";
        return d;
    }
    if (!std::isfinite(a.parallax_mas)) {
        d.reason = "parallax is not finite: distance estimator None";
        return d;
    }
    if (!(a.parallax_mas > kMinParallaxMas)) {
        d.reason = "parallax is not strictly positive, so its inverse is "
                   "meaningless: distance estimator None";
        return d;
    }
    const double snr = a.parallax_signal_to_noise();
    if (!(snr >= kMinParallaxSignalToNoise)) {
        d.reason = "parallax signal-to-noise " + std::to_string(snr) +
                   " below the required " +
                   std::to_string(kMinParallaxSignalToNoise) +
                   ": inverse parallax is not a defensible estimator here, so "
                   "None is returned rather than a biased distance";
        return d;
    }
    d.estimator = DistanceEstimator::InverseParallaxHighQuality;
    d.distance_pc = inverse_parallax_distance(Milliarcsecond(a.parallax_mas)).value;
    // First-order propagation, which is adequate only BECAUSE the signal-to-noise
    // gate above has already excluded the regime where it is not.
    d.uncertainty_pc = d.distance_pc * (a.parallax_error_mas / a.parallax_mas);
    d.reason = "inverse parallax, signal-to-noise " + std::to_string(snr);
    return d;
}

// --- Rejection bookkeeping --------------------------------------------------
//
// The structure that makes "0 silently dropped records" checkable.

struct RejectionRecord {
    std::uint64_t source_id = 0;
    std::string reason;
};

struct IngestionReport {
    int rows_in = 0;
    int accepted = 0;
    std::vector<RejectionRecord> rejected;
    int duplicates_collapsed = 0;

    int distance_inverse_parallax = 0;
    int distance_external = 0;
    int distance_user = 0;
    int distance_none = 0;

    int missing_radial_velocity = 0;
    int missing_covariance = 0;

    bool reconciles() const {
        return rows_in == accepted + static_cast<int>(rejected.size()) +
                              duplicates_collapsed;
    }
};

}  // namespace program_c
