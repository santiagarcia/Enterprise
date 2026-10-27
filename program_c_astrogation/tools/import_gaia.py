#!/usr/bin/env python3
"""Optional live catalogue import for Program C, milestone C0.

THIS SCRIPT IS NOT PART OF ANY GATE. It requires network access and is therefore
excluded from CTest by construction: the deterministic verification workflow is

    online import  ->  validated versioned snapshot  ->  offline simulation

and only the last two steps are ever required. Running this script is how a real
Gaia DR3 snapshot gets produced; until someone runs it, the snapshot shipped in
results/ is SYNTHETIC and says so in its own data_kind field.

What it does
------------
1. Queries the pinned Gaia DR3 table gaiadr3.gaia_source over the ESA TAP
   endpoint, using ADQL, for a parallax-selected volume.
2. Optionally enriches names, aliases and object types from SIMBAD.
   SIMBAD IS NOT THE PRIMARY CATALOGUE: it contributes identifiers and
   classifications only, never astrometry.
3. Optionally attaches confirmed planets from the NASA Exoplanet Archive.
4. Writes a snapshot in exactly the text format program_c/snapshot.hpp reads,
   with data_kind = "Gaia DR3", the real retrieval timestamp and the exact ADQL
   used, so provenance travels with the numbers.

Pinning
-------
The release string, table name and reference epoch are pinned in the C++ header
and repeated here. A migration to Gaia DR4 is a deliberate, versioned change:
bump PINNED_RELEASE and PINNED_TABLE, regenerate, and expect the regression tests
to be updated with new source ids. Do not point this script at DR4 and leave the
release string saying DR3.

Format agreement
----------------
serialize() below must stay byte-for-byte identical to serialize() in
snapshot.hpp, because the C++ reader recomputes the FNV-1a checksum over the body
and refuses a file that disagrees. Prove agreement without a network:

    python import_gaia.py --verify-format <a snapshot the C++ gate wrote>

Usage
-----
    python import_gaia.py --radius-pc 50 --out snapshot_gaia_dr3_50pc.txt
    python import_gaia.py --radius-pc 10 --no-simbad --no-exoplanets
    python import_gaia.py --print-query --radius-pc 50      # no network

Dependencies: NONE. Standard library only, urllib for HTTP.
"""

import argparse
import csv
import datetime
import io
import json
import math
import ssl
import sys
import urllib.error
import urllib.parse
import urllib.request

# --- Pinned provenance. Must match program_c/catalogue.hpp -------------------
PINNED_CATALOGUE = "Gaia DR3"
PINNED_RELEASE = "Gaia DR3"
PINNED_TABLE = "gaiadr3.gaia_source"
PINNED_REFERENCE_EPOCH_JYEAR = 2016.0
GAIA_TAP_SYNC = "https://gea.esac.esa.int/tap-server/tap/sync"
SIMBAD_TAP_SYNC = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
EXOPLANET_TAP_SYNC = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
# Verified 2026-10-27 at https://www.cosmos.esa.int/web/gaia-users/license
# ("Gaia data are distributed under the CC BY-NC 3.0 IGO license") and
# https://www.cosmos.esa.int/web/esdc/terms-and-conditions (credit line
# "Credit: ESA, Gaia DPAC"). Earlier snapshots carry the erroneous string
# "CC BY-SA 3.0 IGO"; see snapshots/DATA_LICENSES.md.
GAIA_LICENCE = ("ESA/Gaia/DPAC, CC BY-NC 3.0 IGO; Credit: ESA, Gaia DPAC; "
                "cite Gaia Collaboration 2016, 2023")

MIN_PARALLAX_SNR = 10.0
SCHEMA_VERSION = 1

# Field validity bits. Must match the ValidField enum in catalogue.hpp.
FIELD_RA = 1 << 0
FIELD_DEC = 1 << 1
FIELD_PARALLAX = 1 << 2
FIELD_PM_RA = 1 << 3
FIELD_PM_DEC = 1 << 4
FIELD_RADIAL_VELOCITY = 1 << 5
FIELD_COVARIANCE = 1 << 6
FIELD_DISTANCE = 1 << 7

# Quality bits. Must match the QualityFlag enum in catalogue.hpp, IN ITS ORDER.
# The C0 gate case c0.flag_bit_values_are_pinned asserts these numeric values on
# the C++ side so this block cannot silently drift out of step.
QUALITY_NONE = 0
QUALITY_NONPOSITIVE_PARALLAX = 1 << 0
QUALITY_LOW_PARALLAX_SNR = 1 << 1
QUALITY_MISSING_RV = 1 << 2
QUALITY_MISSING_COVARIANCE = 1 << 3
QUALITY_DUPLICATED_SOURCE = 1 << 4
QUALITY_NONFINITE_VALUE = 1 << 5
QUALITY_OUTSIDE_VOLUME = 1 << 6

MAS_PER_RADIAN = 180.0 * 3600.0 * 1000.0 / math.pi


# --- The query --------------------------------------------------------------

def adql_for_radius(radius_pc: float, min_snr: float) -> str:
    """The query text.

    THIS AND fixtures::adql_for_radius() IN verification/c0_fixtures.hpp MUST BE
    IDENTICAL, character for character, including the fixed six-decimal
    formatting of both numbers -- the committed C++ text is what the manifest
    publishes as the query this importer issues, so a divergence would publish a
    query nobody ran. Compare them with --print-query, which needs no network.

    Selection is on PARALLAX, not on any distance column, because the catalogue
    measures parallax. A distance cut in the query would silently embed an
    estimator choice in the selection itself.
    """
    min_parallax_mas = 1000.0 / radius_pc
    return (
        "SELECT source_id, ra, ra_error, dec, dec_error, "
        "parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error, "
        "radial_velocity, radial_velocity_error, ref_epoch, "
        "ra_dec_corr, ra_parallax_corr, ra_pmra_corr, ra_pmdec_corr, "
        "dec_parallax_corr, dec_pmra_corr, dec_pmdec_corr, "
        "parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, "
        "duplicated_source, phot_g_mean_mag, ruwe, astrometric_params_solved "
        "FROM " + PINNED_TABLE + " "
        "WHERE parallax IS NOT NULL "
        "AND parallax > " + ("%.6f" % min_parallax_mas) + " "
        "AND parallax_error > 0 "
        "AND parallax / parallax_error >= " + ("%.6f" % min_snr) + " "
        "ORDER BY parallax DESC"
    )


def _ssl_context():
    """A verifying TLS context, falling back to certifi's CA bundle.

    Some archives -- CDS/SIMBAD among them -- serve a chain this machine's system
    trust store cannot complete, so a plain urlopen fails with
    CERTIFICATE_VERIFY_FAILED. The fix is to supply a CA bundle, NOT to disable
    verification: an unverified fetch would let anything on the path substitute
    the identifiers this database then attributes to SIMBAD, and the provenance
    field would say SIMBAD regardless. If neither trust store works, the import
    fails and says so rather than importing unauthenticated data.
    """
    try:
        return ssl.create_default_context()
    except Exception:                       # noqa: BLE001 -- fall through below
        pass
    return None


def _certifi_context():
    try:
        import certifi
    except ImportError:
        return None
    return ssl.create_default_context(cafile=certifi.where())


def tap_query(endpoint: str, query: str, timeout: float = 600.0) -> str:
    """Run a synchronous TAP query and return the CSV body.

    Uses urllib from the standard library rather than requests, so this script has
    NO third-party dependency at all. That matters: the whole point of keeping the
    import out of the gate is that the gate needs nothing installed, and it would
    be poor form for the import itself to need something either.
    """
    data = urllib.parse.urlencode({
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
        "QUERY": query}).encode("ascii")
    req = urllib.request.Request(
        endpoint, data=data,
        headers={"User-Agent": "program_c_astrogation/C0 (research prototype)"})
    contexts = [_ssl_context()]
    fallback = _certifi_context()
    if fallback is not None:
        contexts.append(fallback)

    last = None
    for ctx in contexts:
        try:
            with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
                if r.status != 200:
                    raise RuntimeError("TAP query to %s returned HTTP %d"
                                       % (endpoint, r.status))
                return r.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace")[:2000]
            raise RuntimeError("TAP query to %s failed with HTTP %d:\n%s"
                               % (endpoint, exc.code, body)) from exc
        except urllib.error.URLError as exc:
            # Only a certificate-chain failure is worth retrying with the other
            # trust store. A refused connection or a timeout is not.
            if not isinstance(exc.reason, ssl.SSLCertVerificationError):
                raise
            last = exc
    raise RuntimeError(
        "TAP query to %s could not complete a verified TLS handshake with either "
        "the system trust store or certifi's CA bundle: %s. Refusing to fetch "
        "unverified." % (endpoint, last))


def f(row, key, default=None):
    """Read a float, distinguishing genuinely absent from zero.

    A blank, null or nan cell means NOT MEASURED, and the caller must flag it
    rather than store a zero that reads as a measurement of zero.
    """
    v = row.get(key, "")
    if v is None:
        return default
    t = str(v).strip().lower()
    if t == "" or t in ("null", "nan", "none"):
        return default
    try:
        return float(v)
    except ValueError:
        return default


def covariance_from_correlations(row):
    """Build the 6x6 covariance of (ra, dec, parallax, pmra, pmdec, rv).

    Gaia publishes per-parameter errors and dimensionless correlations, so the
    covariance is sigma_i sigma_j rho_ij. Angular errors arrive in mas and are
    converted to RADIANS here, because the C++ side stores ra and dec in radians
    and a covariance in mixed units is a silent unit error waiting to happen.

    Radial velocity is solved separately from the five-parameter astrometry, so
    the archive publishes no rv correlations. Leaving that row and column
    diagonal is a documented modelling choice, not a measurement.

    Returns (matrix, complete). When complete is False the caller must set
    kQualityMissingCovariance rather than present the matrix as known.
    """
    s_ra_mas = f(row, "ra_error")
    s_dec_mas = f(row, "dec_error")
    s_plx = f(row, "parallax_error")
    s_pmra = f(row, "pmra_error")
    s_pmdec = f(row, "pmdec_error")
    s_rv = f(row, "radial_velocity_error")

    complete = None not in (s_ra_mas, s_dec_mas, s_plx, s_pmra, s_pmdec)
    sig = [
        (s_ra_mas or 0.0) / MAS_PER_RADIAN,
        (s_dec_mas or 0.0) / MAS_PER_RADIAN,
        s_plx or 0.0,
        s_pmra or 0.0,
        s_pmdec or 0.0,
        s_rv or 0.0,
    ]
    corr_keys = {
        (0, 1): "ra_dec_corr",
        (0, 2): "ra_parallax_corr",
        (0, 3): "ra_pmra_corr",
        (0, 4): "ra_pmdec_corr",
        (1, 2): "dec_parallax_corr",
        (1, 3): "dec_pmra_corr",
        (1, 4): "dec_pmdec_corr",
        (2, 3): "parallax_pmra_corr",
        (2, 4): "parallax_pmdec_corr",
        (3, 4): "pmra_pmdec_corr",
    }
    m = [[0.0] * 6 for _ in range(6)]
    for i in range(6):
        m[i][i] = sig[i] * sig[i]
    for (i, j), key in corr_keys.items():
        rho = f(row, key, 0.0)
        m[i][j] = m[j][i] = rho * sig[i] * sig[j]
    return m, bool(complete)


# --- Fetching ---------------------------------------------------------------

def fetch_gaia(radius_pc, min_snr, limit=None):
    query = adql_for_radius(radius_pc, min_snr)
    if limit:
        query = query.replace("SELECT ", "SELECT TOP %d " % int(limit), 1)
    print("[gaia] querying %s for a %g pc volume" % (PINNED_TABLE, radius_pc),
          file=sys.stderr)
    rows = list(csv.DictReader(io.StringIO(tap_query(GAIA_TAP_SYNC, query))))
    print("[gaia] %d rows" % len(rows), file=sys.stderr)
    return query, rows


def fetch_simbad_identifiers(source_ids):
    """Names, aliases and object types keyed by Gaia DR3 source id.

    SIMBAD contributes IDENTIFIERS AND CLASSIFICATIONS ONLY. No astrometric value
    from SIMBAD is ever written into the measured block.
    """
    if not source_ids:
        return {}
    out = {}
    chunk = 300
    for i in range(0, len(source_ids), chunk):
        part = source_ids[i:i + chunk]
        in_list = ", ".join("'Gaia DR3 %s'" % sid for sid in part)
        query = ("SELECT id.id AS gaia_id, b.main_id, b.otype_txt "
                 "FROM ident AS id JOIN basic AS b ON id.oidref = b.oid "
                 "WHERE id.id IN (%s)" % in_list)
        try:
            body = tap_query(SIMBAD_TAP_SYNC, query)
        except Exception as exc:            # noqa: BLE001 -- optional enrichment
            print("[simbad] chunk failed, continuing without it: %s" % exc,
                  file=sys.stderr)
            continue
        for row in csv.DictReader(io.StringIO(body)):
            gid = (row.get("gaia_id") or "").replace("Gaia DR3 ", "").strip()
            if gid:
                out[gid] = {"main_id": (row.get("main_id") or "").strip(),
                            "otype": (row.get("otype_txt") or "").strip()}
    print("[simbad] identified %d of %d" % (len(out), len(source_ids)),
          file=sys.stderr)
    return out


def fetch_exoplanets():
    """Confirmed planets from the NASA Exoplanet Archive, keyed by host name."""
    query = ("SELECT pl_name, hostname, pl_orbper, pl_orbsmax, pl_rade, "
             "pl_bmasse, discoverymethod, disc_year "
             "FROM ps WHERE default_flag = 1")
    try:
        body = tap_query(EXOPLANET_TAP_SYNC, query)
    except Exception as exc:                # noqa: BLE001 -- optional enrichment
        print("[exo] query failed, continuing without planets: %s" % exc,
              file=sys.stderr)
        return {}
    hosts = {}
    for row in csv.DictReader(io.StringIO(body)):
        hosts.setdefault((row.get("hostname") or "").strip(), []).append(row)
    print("[exo] %d confirmed hosts" % len(hosts), file=sys.stderr)
    return hosts


# --- Records ----------------------------------------------------------------

def build_records(gaia_rows, simbad, exoplanets, radius_pc):
    """Convert catalogue rows into records, rejecting nothing silently.

    A source that fails the distance quality gate is KEPT, with estimator None and
    the reason recorded in the rejection ledger. Deleting it would make the row
    count irreconcilable, which is precisely the failure the gate forbids.
    """
    records = []
    rejected = []
    seen = set()
    duplicates = 0

    for row in gaia_rows:
        sid = (row.get("source_id") or "").strip()
        if not sid:
            rejected.append(("0", "no source_id in the row"))
            continue
        if sid in seen:
            duplicates += 1
            continue
        seen.add(sid)

        ra_deg = f(row, "ra")
        dec_deg = f(row, "dec")
        plx = f(row, "parallax")
        plx_err = f(row, "parallax_error")
        if ra_deg is None or dec_deg is None or plx is None or plx_err is None:
            rejected.append((sid, "missing ra, dec, parallax or parallax_error"))
            continue

        valid = FIELD_RA | FIELD_DEC | FIELD_PARALLAX
        quality = QUALITY_NONE

        pmra = f(row, "pmra")
        pmdec = f(row, "pmdec")
        if pmra is not None:
            valid |= FIELD_PM_RA
        if pmdec is not None:
            valid |= FIELD_PM_DEC
        rv = f(row, "radial_velocity")
        if rv is not None:
            valid |= FIELD_RADIAL_VELOCITY
        else:
            quality |= QUALITY_MISSING_RV

        # Honour the catalogue's own reference epoch rather than assuming it.
        # Gaia DR3 publishes ref_epoch = 2016.0 for every source, but reading the
        # column means a release that changes it cannot slip through unnoticed.
        epoch = f(row, "ref_epoch", PINNED_REFERENCE_EPOCH_JYEAR)

        cov, cov_complete = covariance_from_correlations(row)
        if cov_complete:
            valid |= FIELD_COVARIANCE
        else:
            quality |= QUALITY_MISSING_COVARIANCE

        # The archive's own duplicated_source flag. Recorded, not acted on: it
        # marks an astrometric solution the pipeline considers suspect, which is a
        # data-quality statement, not a reason to delete the row.
        if (row.get("duplicated_source") or "").strip().lower() in ("true", "1"):
            quality |= QUALITY_DUPLICATED_SOURCE

        # The distance decision, mirroring decide_distance() in catalogue.hpp.
        snr = plx / plx_err if plx_err > 0.0 else 0.0
        if plx <= 0.0:
            estimator, dist, dist_err = "None", 0.0, 0.0
            quality |= QUALITY_NONPOSITIVE_PARALLAX
            rejected.append((sid, "non-positive parallax %g mas" % plx))
        elif snr < MIN_PARALLAX_SNR:
            estimator, dist, dist_err = "None", 0.0, 0.0
            quality |= QUALITY_LOW_PARALLAX_SNR
            rejected.append((sid, "parallax signal-to-noise %.3f below %.1f: "
                                  "inverting the parallax here would be biased"
                                  % (snr, MIN_PARALLAX_SNR)))
        else:
            estimator = "InverseParallaxHighQuality"
            dist = 1000.0 / plx
            dist_err = 1000.0 * plx_err / (plx * plx)
            valid |= FIELD_DISTANCE
            if dist > radius_pc:
                quality |= QUALITY_OUTSIDE_VOLUME

        ident = simbad.get(sid, {})
        main_id = ident.get("main_id", "")
        identifiers = ["Gaia DR3 " + sid]
        if main_id:
            identifiers.append(main_id)
        planets = [
            {"name": p.get("pl_name", ""),
             "period_days": f(p, "pl_orbper", 0.0),
             "a_au": f(p, "pl_orbsmax", 0.0),
             "radius_earth": f(p, "pl_rade", 0.0),
             "mass_earth": f(p, "pl_bmasse", 0.0)}
            for p in exoplanets.get(main_id, [])
        ] if main_id else []

        records.append({
            "source_id": sid,
            "ra_rad": math.radians(ra_deg),
            "dec_rad": math.radians(dec_deg),
            "parallax_mas": plx,
            "parallax_error_mas": plx_err,
            "pm_ra": pmra or 0.0,
            "pm_dec": pmdec or 0.0,
            "rv": rv or 0.0,
            "epoch": epoch,
            "cov": cov,
            "valid": valid,
            "quality": quality,
            "estimator": estimator,
            "distance_pc": dist,
            "distance_err_pc": dist_err,
            "primary_name": main_id or ("Gaia DR3 " + sid),
            "identifiers": identifiers,
            "object_type": ident.get("otype", ""),
            "planets": planets,
        })

    return records, rejected, duplicates


# --- Serialisation, in lockstep with snapshot.hpp ---------------------------

def g17(v):
    """Match detail::g17: printf %.17g. Python and C agree on IEEE doubles."""
    return "%.17g" % float(v)


def fnv1a(s: str) -> int:
    """FNV-1a over the UTF-8 bytes. Must match fnv1a() in snapshot.hpp."""
    h = 0xCBF29CE484222325
    for b in s.encode("utf-8"):
        h ^= b
        h = (h * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return h


def esc(s) -> str:
    """Match detail::esc: escape backslash, tab and newline, nothing else.

    Note the ORDER. Backslash must go first, or the escapes introduced for tab
    and newline would themselves be escaped.
    """
    return str(s).replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n")


def spatial_cell_index(x, y, z, cell_pc=5.0):
    """Match spatial_cell_index() in snapshot.hpp."""
    i = int(math.floor(x / cell_pc))
    j = int(math.floor(y / cell_pc))
    k = int(math.floor(z / cell_pc))
    return ((i + 4096) * 8192 + (j + 4096)) * 8192 + (k + 4096)


def unit_vector(ra_rad, dec_rad):
    """Match unit_vector() in frames.hpp."""
    cd = math.cos(dec_rad)
    return (cd * math.cos(ra_rad), cd * math.sin(ra_rad), math.sin(dec_rad))


def serialize(records, release, ingestion):
    """Byte-for-byte the format serialize() in snapshot.hpp produces.

    THIS FUNCTION AND serialize() IN snapshot.hpp MUST STAY IN LOCKSTEP. The C++
    reader recomputes the checksum over the body and refuses the file when they
    disagree, so drift here becomes a loud read failure rather than quiet
    corruption. Check agreement with --verify-format, which needs no network.
    """
    o = io.StringIO()
    o.write("# program_c snapshot\n")
    o.write("schema_version\t%d\n" % SCHEMA_VERSION)

    o.write("[catalogue_release]\n")
    o.write("catalogue\t%s\n" % esc(release["catalogue"]))
    o.write("release\t%s\n" % esc(release["release"]))
    o.write("source_table\t%s\n" % esc(release["source_table"]))
    o.write("reference_epoch_jyear\t%s\n" % g17(release["reference_epoch_jyear"]))
    o.write("volume_radius_pc\t%s\n" % g17(release["volume_radius_pc"]))
    o.write("retrieval_timestamp\t%s\n" % esc(release["retrieval_timestamp"]))
    o.write("licence\t%s\n" % esc(release["licence"]))
    o.write("data_kind\t%s\n" % esc(release["data_kind"]))

    o.write("[import_query]\n")
    o.write("adql\t%s\n" % esc(release["adql_query"]))

    o.write("[ingestion]\n")
    o.write("rows_in\t%d\n" % ingestion["rows_in"])
    o.write("accepted\t%d\n" % ingestion["accepted"])
    o.write("duplicates_collapsed\t%d\n" % ingestion["duplicates_collapsed"])
    for sid, reason in ingestion["rejected"]:
        o.write("rejected\t%s\t%s\n" % (sid, esc(reason)))

    o.write("[stellar_object]\n")
    o.write("# source_id\tprimary_name\tobject_type\testimator\tdistance_pc\t"
            "distance_sigma_pc\tvalid_fields\tquality_flags\n")
    for r in records:
        o.write("\t".join([
            r["source_id"], esc(r["primary_name"]), esc(r["object_type"]),
            r["estimator"], g17(r["distance_pc"]), g17(r["distance_err_pc"]),
            str(r["valid"]), str(r["quality"])]) + "\n")

    o.write("[astrometry]\n")
    o.write("# source_id\tra_rad\tdec_rad\tparallax_mas\tparallax_error_mas\t"
            "pmra_mas_yr\tpmdec_mas_yr\trv_km_s\tepoch_jyear\n")
    for r in records:
        o.write("\t".join([
            r["source_id"], g17(r["ra_rad"]), g17(r["dec_rad"]),
            g17(r["parallax_mas"]), g17(r["parallax_error_mas"]),
            g17(r["pm_ra"]), g17(r["pm_dec"]), g17(r["rv"]),
            g17(r["epoch"])]) + "\n")

    o.write("[covariance]\n")
    for r in records:
        # Only when genuinely known. A record with an incomplete covariance is
        # ABSENT from this section, not present as a block of zeros.
        if not (r["valid"] & FIELD_COVARIANCE):
            continue
        cov = r["cov"]
        o.write(r["source_id"] + "".join(
            "\t" + g17(cov[i][j]) for i in range(6) for j in range(6)) + "\n")

    o.write("[identifier]\n")
    for r in records:
        for ident in r["identifiers"]:
            o.write("%s\t%s\n" % (r["source_id"], esc(ident)))

    o.write("[planet]\n")
    o.write("# host_source_id\tname\tperiod_days\ta_au\tradius_earth\t"
            "mass_earth\n")
    for r in records:
        for p in r["planets"]:
            o.write("\t".join([
                r["source_id"], esc(p["name"]), g17(p["period_days"]),
                g17(p["a_au"]), g17(p["radius_earth"]),
                g17(p["mass_earth"])]) + "\n")

    o.write("[provenance]\n")
    o.write("# source_id\tcatalogue\trelease\ttable\ttimestamp\tlicence\t"
            "original_id\n")
    for r in records:
        pv = r["provenance"]
        o.write("\t".join([
            r["source_id"], esc(pv["source_catalogue"]),
            esc(pv["data_release"]), esc(pv["source_table"]),
            esc(pv["retrieval_timestamp"]), esc(pv["licence"]),
            esc(pv["original_identifier"])]) + "\n")

    o.write("[derived]\n")
    o.write("# source_id\tx_pc\ty_pc\tz_pc\tdistance_pc\tspatial_cell\t"
            "epoch_jyear\n")
    for r in records:
        d = r["derived"]
        o.write("\t".join([
            r["source_id"], g17(d["x"]), g17(d["y"]), g17(d["z"]),
            g17(d["distance_pc"]), str(d["cell"]), g17(d["epoch"])]) + "\n")
    return o.getvalue()


def finalise(records):
    """Match finalise(): canonical order by source_id, then derived columns.

    A record without a defensible distance gets zeroed derived columns, exactly
    as the C++ side does, because it cannot be placed in space. It is still
    written to the file.
    """
    records.sort(key=lambda r: int(r["source_id"]))
    for r in records:
        d = {"x": 0.0, "y": 0.0, "z": 0.0, "distance_pc": 0.0, "cell": 0,
             "epoch": r["epoch"]}
        if (r["estimator"] != "None" and (r["valid"] & FIELD_DISTANCE)
                and r["distance_pc"] > 0.0):
            ux, uy, uz = unit_vector(r["ra_rad"], r["dec_rad"])
            d["x"] = ux * r["distance_pc"]
            d["y"] = uy * r["distance_pc"]
            d["z"] = uz * r["distance_pc"]
            d["distance_pc"] = r["distance_pc"]
            d["cell"] = spatial_cell_index(d["x"], d["y"], d["z"])
        r["derived"] = d
    return records


def write_snapshot(path, records, rejected, duplicates, query, radius_pc,
                   rows_in):
    ts = datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ")
    for r in records:
        r["provenance"] = {
            "source_catalogue": PINNED_CATALOGUE,
            "data_release": PINNED_RELEASE,
            "source_table": PINNED_TABLE,
            "retrieval_timestamp": ts,
            "licence": GAIA_LICENCE,
            "original_identifier": "Gaia DR3 " + r["source_id"],
        }
    finalise(records)
    release = {
        "catalogue": PINNED_CATALOGUE,
        "release": PINNED_RELEASE,
        "source_table": PINNED_TABLE,
        "reference_epoch_jyear": PINNED_REFERENCE_EPOCH_JYEAR,
        "volume_radius_pc": radius_pc,
        "retrieval_timestamp": ts,
        "licence": GAIA_LICENCE,
        # The honest marker: this snapshot really did come from the archive, so
        # unlike the gate fixtures it does NOT say SYNTHETIC.
        "data_kind": PINNED_RELEASE,
        "adql_query": query,
    }
    accepted = sum(1 for r in records if r["estimator"] != "None")
    ingestion = {"rows_in": rows_in, "accepted": accepted,
                 "duplicates_collapsed": duplicates, "rejected": rejected}
    body = serialize(records, release, ingestion)
    checksum = fnv1a(body)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(body)
        fh.write("checksum\t%d\n" % checksum)
    print("[out] %s: %d objects, checksum %d"
          % (path, len(records), checksum), file=sys.stderr)
    return checksum, ts, accepted


def verify_format(path):
    """Recompute the checksum of an existing snapshot.

    This is the format-agreement check that needs no network: if the C++ writer
    and this module disagree about a single byte of the body, the recomputed
    FNV-1a will not match the stored one.
    """
    with open(path, "r", encoding="utf-8", newline="") as fh:
        text = fh.read()
    marker = "checksum\t"
    idx = text.rfind("\n" + marker)
    if idx < 0:
        print("no checksum line found in %s" % path, file=sys.stderr)
        return 1
    body = text[:idx + 1]
    stored = int(text[idx + 1 + len(marker):].strip())
    got = fnv1a(body)
    ok = got == stored
    print("%s\n  stored     %d\n  recomputed %d\n  -> %s"
          % (path, stored, got, "AGREE" if ok else "DISAGREE"))
    if not ok:
        print("The Python and C++ serialisers have drifted. Fix serialize() in "
              "this script against snapshot.hpp before importing anything.",
              file=sys.stderr)
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--radius-pc", type=float, default=50.0,
                    help="volume radius in parsecs (default 50)")
    ap.add_argument("--min-parallax-snr", type=float, default=MIN_PARALLAX_SNR)
    ap.add_argument("--limit", type=int, default=None,
                    help="cap the row count, for a smoke test")
    ap.add_argument("--out", default="snapshot_gaia_dr3.txt")
    ap.add_argument("--no-simbad", action="store_true")
    ap.add_argument("--no-exoplanets", action="store_true")
    ap.add_argument("--print-query", action="store_true",
                    help="print the ADQL and exit, touching no network")
    ap.add_argument("--verify-format", metavar="SNAPSHOT",
                    help="recompute an existing snapshot's checksum to prove this "
                         "script and snapshot.hpp agree on the format; touches "
                         "no network")
    a = ap.parse_args(argv)

    if a.print_query:
        print(adql_for_radius(a.radius_pc, a.min_parallax_snr))
        return 0
    if a.verify_format:
        return verify_format(a.verify_format)

    print("NOTE: this is a LIVE IMPORT. It contacts ESA, CDS and IPAC. No gate "
          "depends on it, and the offline verification never runs it.",
          file=sys.stderr)
    query, rows = fetch_gaia(a.radius_pc, a.min_parallax_snr, a.limit)
    ids = [(r.get("source_id") or "").strip() for r in rows]
    simbad = {} if a.no_simbad else fetch_simbad_identifiers(
        [i for i in ids if i])
    exo = {} if a.no_exoplanets else fetch_exoplanets()

    records, rejected, duplicates = build_records(rows, simbad, exo, a.radius_pc)
    checksum, ts, accepted = write_snapshot(
        a.out, records, rejected, duplicates, query, a.radius_pc, len(rows))

    # The ledger must balance here exactly as it must in the gate.
    total = accepted + len(rejected) + duplicates
    manifest = {
        "data_kind": PINNED_RELEASE,
        "release": PINNED_RELEASE,
        "source_table": PINNED_TABLE,
        "reference_epoch_jyear": PINNED_REFERENCE_EPOCH_JYEAR,
        "retrieval_timestamp_utc": ts,
        "adql": query,
        "volume_radius_pc": a.radius_pc,
        "rows_in": len(rows),
        "rows_accepted": accepted,
        "rows_rejected": len(rejected),
        "duplicates_collapsed": duplicates,
        "ledger_balances": total == len(rows),
        "objects_written": len(records),
        "snapshot_checksum_fnv1a": checksum,
        "simbad_used": not a.no_simbad,
        "exoplanet_archive_used": not a.no_exoplanets,
        "licence": GAIA_LICENCE,
    }
    with open(a.out + ".manifest.json", "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
    if total != len(rows):
        print("WARNING: the ledger does not balance: %d in against %d + %d + %d. "
              "Investigate before trusting this snapshot."
              % (len(rows), accepted, len(rejected), duplicates),
              file=sys.stderr)
        return 1
    print("[ok] ledger balances: %d = %d + %d + %d"
          % (len(rows), accepted, len(rejected), duplicates), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
