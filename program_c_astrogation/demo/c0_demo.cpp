#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "c0_fixtures.hpp"

// ---------------------------------------------------------------------------
// The C0 offline demonstration.
//
// Loads the local snapshot and reports the Solar System barycentric origin, the
// nearest systems, 3D positions at the Gaia reference epoch, positions
// propagated to a chosen future epoch, nearest confirmed exoplanet hosts, and
// candidate destinations ranked by distance.
//
// THIS IS A CATALOGUE AND COORDINATE DEMONSTRATION, NOT A VOYAGE SIMULATION.
// Nothing here plans a route, models a ship, or uses Program A or Program B.
//
// It reads the snapshot from disk and makes no network access whatsoever.
// ---------------------------------------------------------------------------

using namespace program_c;

int main(int argc, char** argv) {
    std::string out_dir = "results/verification/program_c0";
    std::string snapshot_path;
    std::string report_name = "c0_demonstration.txt";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out_dir = argv[++i];
        // Point the demonstration at any snapshot, in particular one written by
        // tools/import_gaia.py from the live archive. The demonstration reads the
        // FILE and reports the file's own data_kind, so real and synthetic runs
        // are distinguishable in the output rather than by remembering which was
        // which.
        else if (a == "--snapshot" && i + 1 < argc) snapshot_path = argv[++i];
        else if (a == "--report" && i + 1 < argc) report_name = argv[++i];
        else {
            std::cerr << "usage: c0_demo [--out <dir>] [--snapshot <file>] "
                         "[--report <name>]\n";
            return 2;
        }
    }

    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);

    // Build and write a snapshot if the gate has not already, then read it back:
    // the demonstration works only from the file, never from memory.
    const std::string path =
        snapshot_path.empty() ? out_dir + "/c0_snapshot_50pc.txt" : snapshot_path;
    if (!std::filesystem::exists(path)) {
        const Snapshot s = fixtures::build_snapshot(50.0);
        std::string err;
        if (!write_snapshot(s, path, err)) {
            std::cerr << "cannot write the snapshot: " << err << "\n";
            return 3;
        }
    }

    const LoadedSnapshot snap = read_snapshot(path);
    if (!snap.ok) {
        std::cerr << "cannot read the snapshot: " << snap.message << "\n";
        return 3;
    }

    std::ostringstream o;
    o << "=== Program C, C0 offline demonstration ===\n\n";
    o << "DATA KIND: " << snap.release.data_kind << "\n";
    if (snap.release.data_kind != "Gaia DR3")
        o << "  These are NOT retrieved catalogue records. Every position below is\n"
             "  computed from SYNTHETIC fixtures, so the geometry is exercised and\n"
             "  the astronomy is not. For real records, pass\n"
             "    --snapshot snapshots/gaia_dr3_10pc.txt\n"
             "  or regenerate one with tools/import_gaia.py.\n";
    else
        o << "  Real Gaia DR3 records, retrieved "
          << snap.release.retrieval_timestamp << ".\n  " << snap.release.licence
          << "\n";
    o << "\nsnapshot            " << path << "\n";
    o << "schema version      " << snap.schema_version << "\n";
    o << "checksum verified   " << (snap.checksum_ok ? "yes" : "NO") << "\n";
    o << "pinned release      " << snap.release.release << " / "
      << snap.release.source_table << "\n";
    o << "reference epoch     J" << snap.release.reference_epoch_jyear << "\n";
    o << "volume radius       " << snap.release.volume_radius_pc << " pc\n";
    o << "objects             " << snap.objects.size() << "\n";
    // Say plainly how many can be placed in space. A record without a defensible
    // distance is KEPT in the database but cannot appear in a position listing,
    // and omitting it silently would be exactly the dropped-record failure the
    // gate exists to prevent.
    std::size_t placeable = 0, outside = 0;
    for (const StellarObject& s : snap.objects) {
        if (s.has_defensible_distance()) ++placeable;
        if (s.measured.flagged(kQualityOutsideRequestedVolume)) ++outside;
    }
    o << "placeable           " << placeable << " of " << snap.objects.size() << "  ("
      << (snap.objects.size() - placeable)
      << " retained without a defensible distance, so they cannot be positioned)\n";
    o << "outside the volume  " << outside
      << "  (retained and flagged, never deleted)\n\n";

    // --- The origin ---------------------------------------------------------
    const FrameDescription bary = BarycentricCartesian::description();
    o << "--- Solar System barycentric origin ---\n";
    o << "frame               " << bary.frame << "\n";
    o << "origin              " << bary.origin << "\n";
    o << "orientation         " << bary.orientation << "\n";
    o << "units               " << bary.units << "\n";
    o << "velocity convention " << bary.velocity_convention << "\n";
    o << "position            (0, 0, 0) pc by definition\n\n";

    // --- Nearest systems ----------------------------------------------------
    const std::vector<QueryHit> near = nearest(snap, 100);
    o << "--- nearest " << near.size()
      << " systems, 3D position at the reference epoch ---\n";
    o.setf(std::ios::fixed);
    o << std::setw(26) << std::left << "name" << std::right << std::setw(11) << "d [pc]"
      << std::setw(9) << "+-" << std::setw(11) << "x" << std::setw(11) << "y"
      << std::setw(11) << "z" << "  estimator\n";
    for (const QueryHit& h : near) {
        const StellarObject& s = snap.objects[h.index];
        const Snapshot::Derived& d = snap.derived[h.index];
        // Never print more precision than the uncertainty supports.
        const int digits = significant_digits_for(s.distance_pc,
                                                  s.distance_uncertainty_pc);
        o << std::setw(26) << std::left << s.primary_name << std::right
          << std::setw(11) << std::setprecision(std::min(digits, 6)) << s.distance_pc
          << std::setw(9) << std::setprecision(3) << s.distance_uncertainty_pc
          << std::setw(11) << std::setprecision(4) << d.x_pc << std::setw(11)
          << d.y_pc << std::setw(11) << d.z_pc << "  "
          << to_string(s.distance_estimator)
          << (s.measured.flagged(kQualityOutsideRequestedVolume)
                  ? "  [OUTSIDE REQUESTED VOLUME]" : "")
          << "\n";
    }

    // --- Propagated to a future epoch ---------------------------------------
    const double target_epoch = 12016.0;   // 10 000 yr after the reference epoch
    o << "\n--- positions propagated to J" << std::setprecision(1) << target_epoch
      << " (rectilinear, a LOCAL KINEMATIC model) ---\n";
    o << std::setw(26) << std::left << "name" << std::right << std::setw(12)
      << "moved [pc]" << std::setw(12) << "radial" << std::setw(12) << "transverse"
      << "  warning\n";
    for (const QueryHit& h : near) {
        const StellarObject& s = snap.objects[h.index];
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
        if (!c.ok) continue;
        const PropagationResult p =
            propagate_rectilinear(c.state, JulianYear(target_epoch));
        o << std::setw(26) << std::left << s.primary_name << std::right
          << std::setw(12) << std::setprecision(5) << p.displacement_pc
          << std::setw(12) << p.radial_displacement_pc << std::setw(12)
          << p.transverse_displacement_pc << "  "
          << (p.model_validity_warning ? "VALIDITY WARNING" : "within validity")
          << "\n";
    }

    // --- Exoplanet hosts ----------------------------------------------------
    const std::vector<QueryHit> hosts = exoplanet_hosts(snap);
    o << "\n--- nearest confirmed exoplanet hosts (" << hosts.size() << ") ---\n";
    for (const QueryHit& h : hosts) {
        const StellarObject& s = snap.objects[h.index];
        o << std::setw(26) << std::left << s.primary_name << std::right
          << std::setw(11) << std::setprecision(4) << s.distance_pc << " pc";
        for (const PlanetRecord& p : s.planets)
            o << "   " << p.name << " (P = " << std::setprecision(3)
              << p.orbital_period_days << " d, a = " << p.semi_major_axis_au
              << " au)";
        o << "\n";
    }

    // --- Candidate destinations ---------------------------------------------
    o << "\n--- candidate destinations ranked by distance ---\n";
    o << "Ranked by CURRENT CATALOGUE DISTANCE ONLY. This is not a route: a real\n"
         "voyage must intercept a destination's WORLDLINE, not travel to its stale\n"
         "catalogue position, and that is C1.\n";
    int rank = 0;
    for (const QueryHit& h : near) {
        if (rank++ >= 10) break;
        const StellarObject& s = snap.objects[h.index];
        const double lt = light_travel_time(Parsec(s.distance_pc)).value;
        o << std::setw(3) << std::right << rank << ". " << std::setw(24) << std::left
          << s.primary_name << std::right << std::setw(10) << std::setprecision(4)
          << s.distance_pc << " pc   light time " << std::setw(9)
          << std::setprecision(2) << lt << " yr"
          << (s.planets.empty() ? "" : "   [has planets]") << "\n";
    }

    o << "\nThis is a catalogue and coordinate demonstration. There is no ship, no\n"
         "route optimisation, no warp bubble and no rendering: those are C1 to C3.\n";

    std::cout << o.str();
    std::ofstream f(out_dir + "/" + report_name, std::ios::binary);
    if (f) f << o.str();
    return 0;
}
