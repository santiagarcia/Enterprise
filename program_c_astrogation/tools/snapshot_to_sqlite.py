#!/usr/bin/env python3
"""Materialise a Program C snapshot into a real SQLite database.

WHY THIS IS A SEPARATE TOOL AND NOT PART OF THE C++ CORE
--------------------------------------------------------
C0's authoritative store is the versioned text snapshot, for three reasons that
are recorded here rather than left implicit:

1. It adds no C++ dependency. Vendoring the SQLite amalgamation into a verified
   scientific repository is a real decision with a real licence and audit
   surface, and C0 does not need to make it to be correct.
2. It is diffable and reviewable. A reviewer can read the snapshot and see every
   number, provenance field and rejection reason.
3. It is byte-deterministic and checksummed, so "the same data" is a statement
   about bytes rather than about a binary file whose page layout may vary.

SQLite is the right *query* substrate, so this tool projects the snapshot into
it. The projection is lossless for every column the snapshot carries, and the
snapshot's checksum is stored in the database so a query result can always be
traced back to the exact input file.

This tool is NOT part of any gate: it needs only the Python standard library, but
the deterministic verification path does not depend on it.

Usage
-----
    python snapshot_to_sqlite.py c0_snapshot_50pc.txt --out stars.sqlite
    python snapshot_to_sqlite.py c0_snapshot_50pc.txt --out stars.sqlite --demo

Tables
------
    catalogue_release  one row: release, table, epoch, timestamp, licence,
                       data_kind, adql, snapshot checksum
    ingestion          one row: rows_in, accepted, duplicates_collapsed
    rejection          one row per rejected source, WITH ITS REASON
    stellar_object     one row per object, including the distance estimator
    astrometry         the measured block, with validity and quality bitmasks
    covariance         one row per object that HAS a covariance, 36 columns
    identifier         one row per (source_id, identifier)
    planet             one row per confirmed planet
    provenance         one row per object
    derived            barycentric x, y, z, distance, spatial cell, epoch
"""

import argparse
import os
import sqlite3
import sys

SCHEMA_VERSION = 1

# Must match ValidField / QualityFlag in program_c/catalogue.hpp.
FIELD_NAMES = [
    (1 << 0, "ra"), (1 << 1, "dec"), (1 << 2, "parallax"), (1 << 3, "pm_ra"),
    (1 << 4, "pm_dec"), (1 << 5, "radial_velocity"), (1 << 6, "covariance"),
    (1 << 7, "distance"),
]
# The bit ORDER here is not free: it is the order of the QualityFlag enum, and
# the C0 gate case c0.flag_bit_values_are_pinned asserts the numeric values so
# this table cannot silently drift out of step with the C++ side again.
QUALITY_NAMES = [
    (1 << 0, "parallax_nonpositive"), (1 << 1, "parallax_low_snr"),
    (1 << 2, "missing_radial_velocity"), (1 << 3, "missing_covariance"),
    (1 << 4, "duplicated_source"), (1 << 5, "nonfinite_value"),
    (1 << 6, "outside_requested_volume"),
]


def fnv1a(s: str) -> int:
    """FNV-1a 64-bit. Must match fnv1a() in snapshot.hpp; vectors are pinned by
    the C0 gate case c0.checksum_matches_published_fnv1a_vectors."""
    h = 0xCBF29CE484222325
    for b in s.encode("utf-8"):
        h ^= b
        h = (h * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return h


def unesc(s: str) -> str:
    """Inverse of detail::esc in snapshot.hpp."""
    out = []
    i = 0
    while i < len(s):
        if s[i] == "\\" and i + 1 < len(s):
            n = s[i + 1]
            out.append("\t" if n == "t" else "\n" if n == "n" else n)
            i += 2
        else:
            out.append(s[i])
            i += 1
    return "".join(out)


def decode_mask(mask, names):
    """Turn a bitmask into a readable comma-separated list.

    Stored ALONGSIDE the raw integer, never instead of it: the integer is the
    datum and the text is for humans reading query output.
    """
    return ",".join(name for bit, name in names if mask & bit)


def parse_snapshot(path):
    """Parse the snapshot and verify its checksum.

    Refuses a file whose checksum does not verify. A database silently built from
    a corrupted snapshot is worse than no database, because every downstream
    query then carries an unfalsifiable provenance claim.
    """
    with open(path, "r", encoding="utf-8", newline="") as fh:
        text = fh.read()

    marker = "\nchecksum\t"
    pos = text.rfind(marker)
    if pos < 0:
        sys.exit("%s has no checksum line: refusing to build a database from it"
                 % path)
    body = text[:pos + 1]
    stored = int(text[pos + len(marker):].strip())
    recomputed = fnv1a(body)
    if stored != recomputed:
        sys.exit("%s FAILED its checksum (stored %d, recomputed %d). Refusing to "
                 "build a database from a file whose integrity is unproven."
                 % (path, stored, recomputed))

    data = {
        "checksum": stored,
        "schema_version": 0,
        "release": {},
        "adql": "",
        "ingestion": {"rows_in": 0, "accepted": 0, "duplicates_collapsed": 0},
        "rejected": [],
        "stellar_object": [],
        "astrometry": [],
        "covariance": [],
        "identifier": [],
        "planet": [],
        "provenance": [],
        "derived": [],
    }

    section = None
    for line in body.split("\n"):
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1]
            continue
        parts = line.split("\t")
        key = parts[0]

        if section is None:
            if key == "schema_version":
                data["schema_version"] = int(parts[1])
            continue
        if section == "catalogue_release":
            data["release"][key] = unesc(parts[1]) if len(parts) > 1 else ""
        elif section == "import_query":
            if key == "adql":
                data["adql"] = unesc(parts[1]) if len(parts) > 1 else ""
        elif section == "ingestion":
            if key == "rejected":
                data["rejected"].append((parts[1], unesc(parts[2])
                                         if len(parts) > 2 else ""))
            elif key in data["ingestion"]:
                data["ingestion"][key] = int(parts[1])
        elif section in ("stellar_object", "astrometry", "covariance",
                         "identifier", "planet", "provenance", "derived"):
            data[section].append([unesc(p) for p in parts])

    if data["schema_version"] != SCHEMA_VERSION:
        sys.exit("snapshot schema version %d, this tool understands %d. A schema "
                 "change is a deliberate, versioned event: update this tool."
                 % (data["schema_version"], SCHEMA_VERSION))
    return data


DDL = """
PRAGMA foreign_keys = ON;

CREATE TABLE catalogue_release (
    schema_version         INTEGER NOT NULL,
    catalogue              TEXT    NOT NULL,
    release                TEXT    NOT NULL,
    source_table           TEXT    NOT NULL,
    reference_epoch_jyear  REAL    NOT NULL,
    volume_radius_pc       REAL    NOT NULL,
    retrieval_timestamp    TEXT,
    licence                TEXT,
    -- SYNTHETIC, LITERATURE_FROM_MEMORY or a real release string. Any query
    -- result must be read in the light of this column.
    data_kind              TEXT    NOT NULL,
    adql                   TEXT,
    snapshot_checksum      TEXT    NOT NULL,
    snapshot_path          TEXT    NOT NULL
);

CREATE TABLE ingestion (
    rows_in              INTEGER NOT NULL,
    accepted             INTEGER NOT NULL,
    rejected             INTEGER NOT NULL,
    duplicates_collapsed INTEGER NOT NULL,
    ledger_balances      INTEGER NOT NULL
);

CREATE TABLE rejection (
    source_id TEXT NOT NULL,
    reason    TEXT NOT NULL     -- never empty: a rejection without a reason is a
);                              -- silently dropped record

CREATE TABLE stellar_object (
    source_id          TEXT PRIMARY KEY,
    primary_name       TEXT,
    object_type        TEXT,
    -- Never hidden: every distance in this database says how it was obtained.
    distance_estimator TEXT NOT NULL,
    distance_pc        REAL,
    distance_sigma_pc  REAL,
    valid_fields       INTEGER NOT NULL,
    valid_fields_text  TEXT,
    quality_flags      INTEGER NOT NULL,
    quality_flags_text TEXT
);

CREATE TABLE astrometry (
    source_id          TEXT PRIMARY KEY REFERENCES stellar_object(source_id),
    ra_rad             REAL NOT NULL,
    dec_rad            REAL NOT NULL,
    parallax_mas       REAL NOT NULL,
    parallax_error_mas REAL NOT NULL,
    pm_ra_mas_per_year REAL,          -- mu_alpha* = mu_alpha cos(dec)
    pm_dec_mas_per_year REAL,
    radial_velocity_km_per_s REAL,
    epoch_jyear        REAL NOT NULL,
    parallax_snr       REAL
);

-- Present ONLY for objects whose covariance is known. Absence here means
-- unknown, which is not the same as zero.
CREATE TABLE covariance (
    source_id TEXT PRIMARY KEY REFERENCES stellar_object(source_id),
    %COVCOLS%
);

CREATE TABLE identifier (
    source_id  TEXT NOT NULL REFERENCES stellar_object(source_id),
    identifier TEXT NOT NULL,
    PRIMARY KEY (source_id, identifier)
);

CREATE TABLE planet (
    host_source_id     TEXT NOT NULL REFERENCES stellar_object(source_id),
    name               TEXT NOT NULL,
    orbital_period_days REAL,
    semi_major_axis_au REAL,
    radius_earth       REAL,
    mass_earth         REAL,
    PRIMARY KEY (host_source_id, name)
);

CREATE TABLE provenance (
    source_id           TEXT PRIMARY KEY REFERENCES stellar_object(source_id),
    source_catalogue    TEXT,
    data_release        TEXT,
    source_table        TEXT,
    retrieval_timestamp TEXT,
    licence             TEXT,
    original_identifier TEXT
);

CREATE TABLE derived (
    source_id    TEXT PRIMARY KEY REFERENCES stellar_object(source_id),
    x_pc         REAL NOT NULL,     -- barycentric ICRS-aligned Cartesian
    y_pc         REAL NOT NULL,
    z_pc         REAL NOT NULL,
    distance_pc  REAL NOT NULL,
    spatial_cell INTEGER NOT NULL,  -- 5 pc cells; an index, not a claim
    epoch_jyear  REAL NOT NULL
);

CREATE INDEX idx_derived_distance ON derived(distance_pc);
CREATE INDEX idx_derived_cell     ON derived(spatial_cell);
CREATE INDEX idx_identifier_name  ON identifier(identifier);
CREATE INDEX idx_object_name      ON stellar_object(primary_name);
"""


def covariance_columns():
    labels = ["ra", "dec", "parallax", "pmra", "pmdec", "rv"]
    return ",\n    ".join(
        "c_%s_%s REAL NOT NULL" % (labels[i], labels[j])
        for i in range(6) for j in range(6))


def build(path, out, snapshot):
    if os.path.exists(out):
        os.remove(out)
    con = sqlite3.connect(out)
    con.executescript(DDL.replace("%COVCOLS%", covariance_columns()))

    rel = snapshot["release"]
    con.execute(
        "INSERT INTO catalogue_release VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
        (snapshot["schema_version"], rel.get("catalogue", ""),
         rel.get("release", ""), rel.get("source_table", ""),
         float(rel.get("reference_epoch_jyear", 0.0)),
         float(rel.get("volume_radius_pc", 0.0)),
         rel.get("retrieval_timestamp", ""), rel.get("licence", ""),
         rel.get("data_kind", "UNKNOWN"), snapshot["adql"],
         str(snapshot["checksum"]), os.path.abspath(path)))

    ing = snapshot["ingestion"]
    n_rej = len(snapshot["rejected"])
    balances = int(ing["rows_in"] == ing["accepted"] + n_rej
                   + ing["duplicates_collapsed"])
    con.execute("INSERT INTO ingestion VALUES (?,?,?,?,?)",
                (ing["rows_in"], ing["accepted"], n_rej,
                 ing["duplicates_collapsed"], balances))
    con.executemany("INSERT INTO rejection VALUES (?,?)", snapshot["rejected"])

    for r in snapshot["stellar_object"]:
        valid, quality = int(r[6]), int(r[7])
        con.execute("INSERT INTO stellar_object VALUES (?,?,?,?,?,?,?,?,?,?)",
                    (r[0], r[1], r[2], r[3], float(r[4]), float(r[5]),
                     valid, decode_mask(valid, FIELD_NAMES),
                     quality, decode_mask(quality, QUALITY_NAMES)))

    for r in snapshot["astrometry"]:
        plx, plx_err = float(r[3]), float(r[4])
        snr = plx / plx_err if plx_err > 0.0 else None
        con.execute("INSERT INTO astrometry VALUES (?,?,?,?,?,?,?,?,?,?)",
                    (r[0], float(r[1]), float(r[2]), plx, plx_err,
                     float(r[5]), float(r[6]), float(r[7]), float(r[8]), snr))

    if snapshot["covariance"]:
        placeholders = ",".join(["?"] * 37)
        for r in snapshot["covariance"]:
            con.execute("INSERT INTO covariance VALUES (%s)" % placeholders,
                        [r[0]] + [float(v) for v in r[1:37]])

    con.executemany("INSERT INTO identifier VALUES (?,?)",
                    [(r[0], r[1]) for r in snapshot["identifier"]])
    con.executemany("INSERT INTO planet VALUES (?,?,?,?,?,?)",
                    [(r[0], r[1], float(r[2]), float(r[3]), float(r[4]),
                      float(r[5])) for r in snapshot["planet"]])
    con.executemany("INSERT INTO provenance VALUES (?,?,?,?,?,?,?)",
                    [tuple(r[:7]) for r in snapshot["provenance"]])
    con.executemany("INSERT INTO derived VALUES (?,?,?,?,?,?,?)",
                    [(r[0], float(r[1]), float(r[2]), float(r[3]),
                      float(r[4]), int(r[5]), float(r[6]))
                     for r in snapshot["derived"]])
    con.commit()

    if not balances:
        print("WARNING: the ingestion ledger in this snapshot does not balance: "
              "%d in against %d + %d + %d."
              % (ing["rows_in"], ing["accepted"], n_rej,
                 ing["duplicates_collapsed"]), file=sys.stderr)
    return con


def demo_queries(con):
    """A few queries, to show the projection is usable and honest."""
    print("\n--- catalogue_release ---")
    row = con.execute(
        "SELECT release, source_table, reference_epoch_jyear, data_kind, "
        "snapshot_checksum FROM catalogue_release").fetchone()
    print("release %s / table %s / epoch J%g" % (row[0], row[1], row[2]))
    print("data_kind %s   snapshot checksum %s" % (row[3], row[4]))
    if row[3] != row[0]:
        print("  NOTE: data_kind differs from the release string, so these rows "
              "are NOT retrieved catalogue data.")

    print("\n--- ledger ---")
    print("rows_in %d = accepted %d + rejected %d + duplicates %d (balances: %d)"
          % con.execute("SELECT rows_in, accepted, rejected, "
                        "duplicates_collapsed, ledger_balances "
                        "FROM ingestion").fetchone())

    print("\n--- nearest 10, with the estimator named ---")
    for r in con.execute(
            "SELECT o.primary_name, d.distance_pc, o.distance_sigma_pc, "
            "       d.x_pc, d.y_pc, d.z_pc, o.distance_estimator "
            "FROM derived d JOIN stellar_object o USING (source_id) "
            "WHERE d.distance_pc > 0 ORDER BY d.distance_pc LIMIT 10"):
        print("  %-24s %8.4f +- %.4f pc  (%8.4f, %8.4f, %8.4f)  %s" % r)

    print("\n--- objects with no defensible distance (KEPT, not deleted) ---")
    rows = con.execute(
        "SELECT o.source_id, o.primary_name, o.quality_flags_text, r.reason "
        "FROM stellar_object o LEFT JOIN rejection r USING (source_id) "
        "WHERE o.distance_estimator = 'None'").fetchall()
    for r in rows:
        print("  %s %-24s [%s]\n      %s" % r)
    if not rows:
        print("  none")

    print("\n--- confirmed exoplanet hosts ---")
    for r in con.execute(
            "SELECT o.primary_name, d.distance_pc, p.name, "
            "       p.orbital_period_days, p.semi_major_axis_au "
            "FROM planet p JOIN stellar_object o "
            "  ON o.source_id = p.host_source_id "
            "JOIN derived d USING (source_id) ORDER BY d.distance_pc"):
        print("  %-24s %8.4f pc   %s (P = %g d, a = %g au)" % r)

    print("\n--- objects missing a covariance (unknown, NOT zero) ---")
    n = con.execute(
        "SELECT COUNT(*) FROM stellar_object o "
        "WHERE NOT EXISTS (SELECT 1 FROM covariance c "
        "                  WHERE c.source_id = o.source_id)").fetchone()[0]
    print("  %d of %d" % (n, con.execute(
        "SELECT COUNT(*) FROM stellar_object").fetchone()[0]))


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("snapshot", help="a Program C text snapshot")
    ap.add_argument("--out", default=None,
                    help="output .sqlite path (default: alongside the snapshot)")
    ap.add_argument("--demo", action="store_true",
                    help="run a few demonstration queries afterwards")
    a = ap.parse_args(argv)

    out = a.out or (os.path.splitext(a.snapshot)[0] + ".sqlite")
    snapshot = parse_snapshot(a.snapshot)
    con = build(a.snapshot, out, snapshot)

    counts = {t: con.execute("SELECT COUNT(*) FROM %s" % t).fetchone()[0]
              for t in ("stellar_object", "astrometry", "covariance",
                        "identifier", "planet", "provenance", "derived",
                        "rejection")}
    print("%s written from %s (checksum %d verified)"
          % (out, a.snapshot, snapshot["checksum"]))
    for t in sorted(counts):
        print("  %-16s %d" % (t, counts[t]))

    if a.demo:
        demo_queries(con)
    con.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
