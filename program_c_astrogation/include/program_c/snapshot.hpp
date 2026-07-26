#pragma once

#include <fstream>
#include <map>
#include <sstream>

#include "program_c/astrometry.hpp"

// ---------------------------------------------------------------------------
// The local snapshot: a portable, versioned, deterministic star database.
//
// WHY NOT SQLITE, YET. SQLite is the recommended first implementation and it
// remains the intended destination. No SQLite amalgamation is vendored here, and
// adding a ~250 000-line third-party C source file is a dependency decision for
// the project owner rather than something to slip into a milestone -- the same
// reasoning that keeps OTIlib an explicit opt-in.
//
// So the C++ core writes a dependency-free versioned text snapshot with the SAME
// LOGICAL TABLES, and `tools/snapshot_to_sqlite.py` materialises exactly those
// tables into a real .sqlite file using Python's built-in sqlite3. The SQLite
// deliverable therefore exists, and the C++ library and gate stay dependency-free
// and byte-deterministic.
//
// DETERMINISM. Records are written in ascending source_id, every float is written
// with %.17g, and the checksum is FNV-1a over the canonical serialization
// EXCLUDING the checksum line. Two runs over the same data produce byte-identical
// files, which is what the gate asserts.
//
// NETWORK ACCESS IS NOT PART OF THE RUNTIME. This header reads and writes files
// and nothing else. Online work lives in tools/, which the library never calls.
// ---------------------------------------------------------------------------

namespace program_c {

inline constexpr int kSnapshotSchemaVersion = 1;

// FNV-1a, 64-bit. Chosen for being short, deterministic and dependency-free; it
// is an integrity check against accidental corruption and reordering, not a
// cryptographic guarantee, and it is labelled as such in the manifest.
inline std::uint64_t fnv1a(const std::string& s) {
    // The standard 64-bit FNV-1a offset basis, 0xcbf29ce484222325. Spelled in
    // hexadecimal deliberately: the decimal form 14695981039346656037 is one
    // digit away from a plausible-looking wrong number, and an independent
    // reimplementation of this function must agree with it bit for bit, or a
    // snapshot written by one tool is unreadable by the other.
    std::uint64_t h = 0xcbf29ce484222325ull;
    for (unsigned char c : s) {
        h ^= static_cast<std::uint64_t>(c);
        h *= 1099511628211ull;
    }
    return h;
}

struct CatalogueRelease {
    std::string catalogue = kPinnedCatalogue;
    std::string release = kPinnedCatalogue;
    std::string source_table = kPinnedTable;
    double reference_epoch_jyear = kGaiaDR3ReferenceEpochJyear;
    double volume_radius_pc = 50.0;
    std::string retrieval_timestamp;
    std::string adql_query;
    std::string licence;
    // The honesty marker, propagated into the manifest and asserted by the gate.
    std::string data_kind = "SYNTHETIC";
};

struct Snapshot {
    int schema_version = kSnapshotSchemaVersion;
    CatalogueRelease release;
    std::vector<StellarObject> objects;
    IngestionReport ingestion;

    // Derived columns, stored rather than recomputed on read.
    struct Derived {
        double x_pc = 0.0, y_pc = 0.0, z_pc = 0.0;
        double distance_pc = 0.0;
        long long spatial_cell = 0;
        double epoch_jyear = 0.0;
    };
    std::vector<Derived> derived;

    std::uint64_t checksum = 0;
};

// A coarse spatial cell index, for a first pass at spatial queries. 5 pc cells,
// packed into one integer. Deliberately coarse and deliberately explicit: it is
// an index, not a claim about structure.
inline long long spatial_cell_index(double x, double y, double z,
                                   double cell_pc = 5.0) {
    const long long i = static_cast<long long>(std::floor(x / cell_pc));
    const long long j = static_cast<long long>(std::floor(y / cell_pc));
    const long long k = static_cast<long long>(std::floor(z / cell_pc));
    return ((i + 4096) * 8192 + (j + 4096)) * 8192 + (k + 4096);
}

namespace detail {

inline std::string g17(double v) {
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return buf;
}

// Field escaping: tab and newline are the separators, so they are escaped rather
// than assumed absent from a catalogue name.
inline std::string esc(const std::string& s) {
    std::string o;
    for (char c : s) {
        if (c == '\t') o += "\\t";
        else if (c == '\n') o += "\\n";
        else if (c == '\\') o += "\\\\";
        else o += c;
    }
    return o;
}

inline std::string unesc(const std::string& s) {
    std::string o;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '\\' && i + 1 < s.size()) {
            ++i;
            if (s[i] == 't') o += '\t';
            else if (s[i] == 'n') o += '\n';
            else o += s[i];
        } else {
            o += s[i];
        }
    }
    return o;
}

inline std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : line) {
        if (c == '\t') { out.push_back(cur); cur.clear(); }
        else cur += c;
    }
    out.push_back(cur);
    return out;
}

}  // namespace detail

// The canonical serialization. One section per logical table, so the mapping to
// SQLite is one-to-one.
inline std::string serialize(const Snapshot& snap) {
    std::ostringstream o;
    o << "# program_c snapshot\n";
    o << "schema_version\t" << snap.schema_version << "\n";

    o << "[catalogue_release]\n";
    o << "catalogue\t" << detail::esc(snap.release.catalogue) << "\n";
    o << "release\t" << detail::esc(snap.release.release) << "\n";
    o << "source_table\t" << detail::esc(snap.release.source_table) << "\n";
    o << "reference_epoch_jyear\t"
      << detail::g17(snap.release.reference_epoch_jyear) << "\n";
    o << "volume_radius_pc\t" << detail::g17(snap.release.volume_radius_pc) << "\n";
    o << "retrieval_timestamp\t" << detail::esc(snap.release.retrieval_timestamp)
      << "\n";
    o << "licence\t" << detail::esc(snap.release.licence) << "\n";
    o << "data_kind\t" << detail::esc(snap.release.data_kind) << "\n";

    o << "[import_query]\n";
    o << "adql\t" << detail::esc(snap.release.adql_query) << "\n";

    o << "[ingestion]\n";
    o << "rows_in\t" << snap.ingestion.rows_in << "\n";
    o << "accepted\t" << snap.ingestion.accepted << "\n";
    o << "duplicates_collapsed\t" << snap.ingestion.duplicates_collapsed << "\n";
    for (const RejectionRecord& r : snap.ingestion.rejected)
        o << "rejected\t" << r.source_id << "\t" << detail::esc(r.reason) << "\n";

    o << "[stellar_object]\n";
    o << "# source_id\tprimary_name\tobject_type\testimator\tdistance_pc\t"
         "distance_sigma_pc\tvalid_fields\tquality_flags\n";
    for (const StellarObject& s : snap.objects)
        o << s.measured.gaia_source_id << "\t" << detail::esc(s.primary_name)
          << "\t" << detail::esc(s.object_type) << "\t"
          << to_string(s.distance_estimator) << "\t"
          << detail::g17(s.distance_pc) << "\t"
          << detail::g17(s.distance_uncertainty_pc) << "\t"
          << s.measured.valid_fields << "\t" << s.measured.quality_flags << "\n";

    o << "[astrometry]\n";
    o << "# source_id\tra_rad\tdec_rad\tparallax_mas\tparallax_error_mas\t"
         "pmra_mas_yr\tpmdec_mas_yr\trv_km_s\tepoch_jyear\n";
    for (const StellarObject& s : snap.objects) {
        const CatalogueAstrometry& a = s.measured;
        o << a.gaia_source_id << "\t" << detail::g17(a.ra_rad) << "\t"
          << detail::g17(a.dec_rad) << "\t" << detail::g17(a.parallax_mas) << "\t"
          << detail::g17(a.parallax_error_mas) << "\t"
          << detail::g17(a.pm_ra_mas_per_year) << "\t"
          << detail::g17(a.pm_dec_mas_per_year) << "\t"
          << detail::g17(a.radial_velocity_km_per_s) << "\t"
          << detail::g17(a.reference_epoch_jyear) << "\n";
    }

    o << "[covariance]\n";
    for (const StellarObject& s : snap.objects) {
        if (!s.measured.has(kFieldCovariance)) continue;
        o << s.measured.gaia_source_id;
        for (int i = 0; i < 36; ++i) o << "\t" << detail::g17(s.measured.covariance.m[i]);
        o << "\n";
    }

    o << "[identifier]\n";
    for (const StellarObject& s : snap.objects)
        for (const std::string& id : s.identifiers)
            o << s.measured.gaia_source_id << "\t" << detail::esc(id) << "\n";

    o << "[planet]\n";
    o << "# host_source_id\tname\tperiod_days\ta_au\tradius_earth\tmass_earth\n";
    for (const StellarObject& s : snap.objects)
        for (const PlanetRecord& p : s.planets)
            o << s.measured.gaia_source_id << "\t" << detail::esc(p.name) << "\t"
              << detail::g17(p.orbital_period_days) << "\t"
              << detail::g17(p.semi_major_axis_au) << "\t"
              << detail::g17(p.radius_earth) << "\t" << detail::g17(p.mass_earth)
              << "\n";

    o << "[provenance]\n";
    o << "# source_id\tcatalogue\trelease\ttable\ttimestamp\tlicence\toriginal_id\n";
    for (const StellarObject& s : snap.objects) {
        const Provenance& p = s.provenance;
        o << s.measured.gaia_source_id << "\t" << detail::esc(p.source_catalogue)
          << "\t" << detail::esc(p.data_release) << "\t"
          << detail::esc(p.source_table) << "\t"
          << detail::esc(p.retrieval_timestamp) << "\t" << detail::esc(p.licence)
          << "\t" << detail::esc(p.original_identifier) << "\n";
    }

    o << "[derived]\n";
    o << "# source_id\tx_pc\ty_pc\tz_pc\tdistance_pc\tspatial_cell\tepoch_jyear\n";
    for (std::size_t i = 0; i < snap.objects.size() && i < snap.derived.size(); ++i) {
        const Snapshot::Derived& d = snap.derived[i];
        o << snap.objects[i].measured.gaia_source_id << "\t" << detail::g17(d.x_pc)
          << "\t" << detail::g17(d.y_pc) << "\t" << detail::g17(d.z_pc) << "\t"
          << detail::g17(d.distance_pc) << "\t" << d.spatial_cell << "\t"
          << detail::g17(d.epoch_jyear) << "\n";
    }
    return o.str();
}

// Sort into canonical order and fill the derived columns and checksum.
inline void finalise(Snapshot& snap) {
    // Stable canonical order: ascending source_id.
    std::vector<std::size_t> order(snap.objects.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return snap.objects[a].measured.gaia_source_id <
               snap.objects[b].measured.gaia_source_id;
    });
    std::vector<StellarObject> sorted;
    sorted.reserve(snap.objects.size());
    for (std::size_t i : order) sorted.push_back(snap.objects[i]);
    snap.objects.swap(sorted);

    snap.derived.clear();
    snap.derived.reserve(snap.objects.size());
    for (const StellarObject& s : snap.objects) {
        Snapshot::Derived d;
        d.epoch_jyear = s.measured.reference_epoch_jyear;
        if (s.has_defensible_distance()) {
            ICRSPhaseSpace o;
            o.right_ascension = Radian(s.measured.ra_rad);
            o.declination = Radian(s.measured.dec_rad);
            o.parallax = Milliarcsecond(s.measured.parallax_mas);
            o.pm_ra_cosdec = MasPerYear(s.measured.pm_ra_mas_per_year);
            o.pm_dec = MasPerYear(s.measured.pm_dec_mas_per_year);
            o.radial_velocity = KmPerSecond(s.measured.radial_velocity_km_per_s);
            o.epoch = JulianYear(s.measured.reference_epoch_jyear);
            const ConversionResult c =
                to_barycentric(o, Parsec(s.distance_pc), s.distance_estimator);
            if (c.ok) {
                d.x_pc = c.state.position_pc.x;
                d.y_pc = c.state.position_pc.y;
                d.z_pc = c.state.position_pc.z;
                d.distance_pc = s.distance_pc;
                d.spatial_cell = spatial_cell_index(d.x_pc, d.y_pc, d.z_pc);
            }
        }
        snap.derived.push_back(d);
    }
    snap.checksum = fnv1a(serialize(snap));
}

inline bool write_snapshot(const Snapshot& snap, const std::string& path,
                          std::string& error) {
    std::ofstream f(path, std::ios::binary);
    if (!f) { error = "cannot open " + path; return false; }
    f << serialize(snap);
    f << "checksum\t" << snap.checksum << "\n";
    return true;
}

// Read back. Only the fields the offline query layer needs are parsed; the rest
// is available in the file and is not silently discarded from it.
struct LoadedSnapshot {
    int schema_version = 0;
    CatalogueRelease release;
    std::vector<StellarObject> objects;
    std::vector<Snapshot::Derived> derived;
    std::uint64_t checksum_in_file = 0;
    std::uint64_t checksum_recomputed = 0;
    bool checksum_ok = false;
    bool ok = false;
    std::string message;
};

inline LoadedSnapshot read_snapshot(const std::string& path) {
    LoadedSnapshot out;
    std::ifstream f(path, std::ios::binary);
    if (!f) { out.message = "cannot open " + path; return out; }

    std::string all((std::istreambuf_iterator<char>(f)),
                    std::istreambuf_iterator<char>());
    // Split off the trailing checksum line for recomputation.
    const std::string marker = "\nchecksum\t";
    const std::size_t pos = all.rfind(marker);
    if (pos == std::string::npos) {
        out.message = "no checksum line";
        return out;
    }
    const std::string body = all.substr(0, pos + 1);
    out.checksum_in_file = std::strtoull(all.c_str() + pos + marker.size(), nullptr, 10);
    out.checksum_recomputed = fnv1a(body);
    out.checksum_ok = out.checksum_in_file == out.checksum_recomputed;

    std::istringstream in(body);
    std::string line, section;
    std::map<std::uint64_t, std::size_t> index;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        if (line[0] == '[') { section = line; continue; }
        const std::vector<std::string> f2 = detail::split_tabs(line);
        if (section.empty()) {
            if (f2[0] == "schema_version" && f2.size() > 1)
                out.schema_version = std::atoi(f2[1].c_str());
        } else if (section == "[catalogue_release]" && f2.size() > 1) {
            const std::string v = detail::unesc(f2[1]);
            if (f2[0] == "catalogue") out.release.catalogue = v;
            else if (f2[0] == "release") out.release.release = v;
            else if (f2[0] == "source_table") out.release.source_table = v;
            else if (f2[0] == "reference_epoch_jyear")
                out.release.reference_epoch_jyear = std::atof(v.c_str());
            else if (f2[0] == "volume_radius_pc")
                out.release.volume_radius_pc = std::atof(v.c_str());
            else if (f2[0] == "retrieval_timestamp")
                out.release.retrieval_timestamp = v;
            else if (f2[0] == "licence") out.release.licence = v;
            else if (f2[0] == "data_kind") out.release.data_kind = v;
        } else if (section == "[import_query]" && f2.size() > 1) {
            if (f2[0] == "adql") out.release.adql_query = detail::unesc(f2[1]);
        } else if (section == "[stellar_object]" && f2.size() >= 8) {
            StellarObject s;
            s.measured.gaia_source_id = std::strtoull(f2[0].c_str(), nullptr, 10);
            s.primary_name = detail::unesc(f2[1]);
            s.object_type = detail::unesc(f2[2]);
            const std::string est = f2[3];
            s.distance_estimator =
                est == "InverseParallaxHighQuality"
                    ? DistanceEstimator::InverseParallaxHighQuality
                    : est == "ExternalPosterior" ? DistanceEstimator::ExternalPosterior
                    : est == "UserSupplied" ? DistanceEstimator::UserSupplied
                                            : DistanceEstimator::None;
            s.distance_pc = std::atof(f2[4].c_str());
            s.distance_uncertainty_pc = std::atof(f2[5].c_str());
            s.measured.valid_fields = std::strtoull(f2[6].c_str(), nullptr, 10);
            s.measured.quality_flags = std::strtoull(f2[7].c_str(), nullptr, 10);
            index[s.measured.gaia_source_id] = out.objects.size();
            out.objects.push_back(s);
        } else if (section == "[astrometry]" && f2.size() >= 9) {
            const std::uint64_t id = std::strtoull(f2[0].c_str(), nullptr, 10);
            auto it = index.find(id);
            if (it == index.end()) continue;
            CatalogueAstrometry& a = out.objects[it->second].measured;
            a.ra_rad = std::atof(f2[1].c_str());
            a.dec_rad = std::atof(f2[2].c_str());
            a.parallax_mas = std::atof(f2[3].c_str());
            a.parallax_error_mas = std::atof(f2[4].c_str());
            a.pm_ra_mas_per_year = std::atof(f2[5].c_str());
            a.pm_dec_mas_per_year = std::atof(f2[6].c_str());
            a.radial_velocity_km_per_s = std::atof(f2[7].c_str());
            a.reference_epoch_jyear = std::atof(f2[8].c_str());
        } else if (section == "[covariance]" && f2.size() >= 37) {
            const std::uint64_t id = std::strtoull(f2[0].c_str(), nullptr, 10);
            auto it = index.find(id);
            if (it == index.end()) continue;
            for (int i = 0; i < 36; ++i)
                out.objects[it->second].measured.covariance.m[i] =
                    std::atof(f2[static_cast<std::size_t>(i) + 1].c_str());
        } else if (section == "[identifier]" && f2.size() >= 2) {
            const std::uint64_t id = std::strtoull(f2[0].c_str(), nullptr, 10);
            auto it = index.find(id);
            if (it != index.end())
                out.objects[it->second].identifiers.push_back(detail::unesc(f2[1]));
        } else if (section == "[planet]" && f2.size() >= 6) {
            const std::uint64_t id = std::strtoull(f2[0].c_str(), nullptr, 10);
            auto it = index.find(id);
            if (it == index.end()) continue;
            PlanetRecord p;
            p.name = detail::unesc(f2[1]);
            p.orbital_period_days = std::atof(f2[2].c_str());
            p.semi_major_axis_au = std::atof(f2[3].c_str());
            p.radius_earth = std::atof(f2[4].c_str());
            p.mass_earth = std::atof(f2[5].c_str());
            out.objects[it->second].planets.push_back(p);
        } else if (section == "[provenance]" && f2.size() >= 7) {
            const std::uint64_t id = std::strtoull(f2[0].c_str(), nullptr, 10);
            auto it = index.find(id);
            if (it == index.end()) continue;
            Provenance& p = out.objects[it->second].provenance;
            p.source_catalogue = detail::unesc(f2[1]);
            p.data_release = detail::unesc(f2[2]);
            p.source_table = detail::unesc(f2[3]);
            p.retrieval_timestamp = detail::unesc(f2[4]);
            p.licence = detail::unesc(f2[5]);
            p.original_identifier = detail::unesc(f2[6]);
            p.query = out.release.adql_query;
        } else if (section == "[derived]" && f2.size() >= 7) {
            Snapshot::Derived d;
            d.x_pc = std::atof(f2[1].c_str());
            d.y_pc = std::atof(f2[2].c_str());
            d.z_pc = std::atof(f2[3].c_str());
            d.distance_pc = std::atof(f2[4].c_str());
            d.spatial_cell = std::atoll(f2[5].c_str());
            d.epoch_jyear = std::atof(f2[6].c_str());
            out.derived.push_back(d);
        }
    }
    out.ok = true;
    out.message = "loaded " + std::to_string(out.objects.size()) + " objects";
    return out;
}

// --- Offline queries --------------------------------------------------------
//
// The whole query surface reads a loaded snapshot. No network, ever.

struct QueryHit {
    std::size_t index = 0;
    double distance_pc = 0.0;
};

inline std::vector<QueryHit> nearest(const LoadedSnapshot& s, std::size_t count) {
    std::vector<QueryHit> hits;
    for (std::size_t i = 0; i < s.objects.size() && i < s.derived.size(); ++i) {
        if (!s.objects[i].has_defensible_distance()) continue;
        hits.push_back({i, s.derived[i].distance_pc});
    }
    std::sort(hits.begin(), hits.end(), [](const QueryHit& a, const QueryHit& b) {
        if (a.distance_pc != b.distance_pc) return a.distance_pc < b.distance_pc;
        return a.index < b.index;
    });
    if (hits.size() > count) hits.resize(count);
    return hits;
}

inline std::vector<QueryHit> within_radius(const LoadedSnapshot& s,
                                          double radius_pc) {
    std::vector<QueryHit> hits;
    for (std::size_t i = 0; i < s.objects.size() && i < s.derived.size(); ++i) {
        if (!s.objects[i].has_defensible_distance()) continue;
        if (s.derived[i].distance_pc <= radius_pc)
            hits.push_back({i, s.derived[i].distance_pc});
    }
    std::sort(hits.begin(), hits.end(), [](const QueryHit& a, const QueryHit& b) {
        if (a.distance_pc != b.distance_pc) return a.distance_pc < b.distance_pc;
        return a.index < b.index;
    });
    return hits;
}

inline std::vector<QueryHit> exoplanet_hosts(const LoadedSnapshot& s) {
    std::vector<QueryHit> hits;
    for (std::size_t i = 0; i < s.objects.size() && i < s.derived.size(); ++i)
        if (!s.objects[i].planets.empty())
            hits.push_back({i, s.derived[i].distance_pc});
    std::sort(hits.begin(), hits.end(), [](const QueryHit& a, const QueryHit& b) {
        if (a.distance_pc != b.distance_pc) return a.distance_pc < b.distance_pc;
        return a.index < b.index;
    });
    return hits;
}

inline const StellarObject* find_by_identifier(const LoadedSnapshot& s,
                                             const std::string& id) {
    for (const StellarObject& o : s.objects) {
        if (o.primary_name == id) return &o;
        for (const std::string& x : o.identifiers)
            if (x == id) return &o;
    }
    return nullptr;
}

}  // namespace program_c
