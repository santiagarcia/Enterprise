# Program C, milestone C0 — catalogue-grounded astrogation

    catalogues -> time-dependent celestial state -> route and ship worldline
      -> Program A metric -> Program B electromagnetic observables
      -> virtual navigation environment

C0 builds only the first link: a **local, offline, deterministic star database**
with typed coordinates, epoch propagation and honest uncertainty.

## Scope

C0 contains: the astronomical data schema; catalogue ingestion; coordinate
transformations; epoch propagation; uncertainty and data-quality handling; a
searchable local database; deterministic verification.

C0 does **not** simulate a warp bubble, optimize a voyage, render a cockpit, or
claim full Galactic dynamics.

## Separation from Programs A and B

Program C is a separate directory, a separate CMake target set and a separate
gate. It may consume the public interfaces of `warp_hypad` and
`program_b_em_inverse_design`; **nothing about catalogues, routes or UI enters
their cores.** At C0 it consumes neither — the dependency edge is created at C2,
where a trajectory adapter is written.

---

## THE DATA HONESTY PROBLEM, STATED FIRST

**This milestone was implemented without network access, so no Gaia DR3 record
was retrieved.** Every fixture committed here is one of two clearly separated
kinds, and neither is presented as retrieved catalogue data:

| kind | provenance tag | what it is | what it is used for |
|---|---|---|---|
| **synthetic** | `SYNTHETIC` | constructed so the correct answer is known in closed form | the numerical gate: round trips, propagation, covariance, snapshot determinism |
| **literature reference** | `LITERATURE_FROM_MEMORY` | approximate well-known values recorded from memory | order-of-magnitude sanity checks at loose tolerance only |

**No file in this repository claims a Gaia `source_id`, a retrieval timestamp or
a release string it does not have.** Writing plausible numbers under a real
catalogue's name would be fabricating records, and the gate would then be
verifying a fiction. Every synthetic record carries
`source_catalogue = "SYNTHETIC"` and every literature value carries
`LITERATURE_FROM_MEMORY` together with a note that it must be replaced before any
result depending on it is quoted.

`tools/import_gaia.py` is the live client. Running it against the real TAP
service is what produces fixtures with genuine provenance, and it is a **separate
optional step**, never part of CTest. Until it is run, the C0 gate verifies the
*machinery* and says so.

---

## 1. Pinned data sources

    catalogue        Gaia DR3
    table            gaiadr3.gaia_source
    reference epoch  J2016.0  (Gaia DR3 reference epoch)
    query language   ADQL over TAP
    endpoint         https://gea.esac.esa.int/tap-server/tap

**Migration to Gaia DR4 must be a deliberate, versioned catalogue update with
regression tests**, exactly as the OTIlib pin is enforced elsewhere in this
repository. The release string is stored in the snapshot header and the gate
asserts it, so a silent migration fails.

SIMBAD supplies **only** primary names, aliases, object classifications and
cross-identifications. It is never the primary star catalogue.

The NASA Exoplanet Archive supplies confirmed planet records, host-star
associations and selected planetary properties.

Every imported record retains a full `Provenance`: source catalogue, data
release, source table, the exact query text, retrieval timestamp, licence and the
original identifier.

## 2. First universe volume

A configurable sphere around the Solar System, default `R_max = 50 pc`, with 5 pc
and 10 pc verification volumes. The whole Gaia catalogue is never ingested.

The exact ADQL is stored beside the result, in `import_query`.

**`d = 1/parallax` is not treated as a universally unbiased distance estimate.**
Three estimators are supported and the one used is always recorded:

| estimator | when |
|---|---|
| `InverseParallaxHighQuality` | only where the parallax quality gate passes |
| `ExternalPosterior` | a distance supplied by an external posterior catalogue |
| `UserSupplied` | an explicit override, marked as such |
| `None` | no defensible distance; the record is retained and flagged, not deleted |

The quality gate for `InverseParallaxHighQuality` is stated numerically in
`catalogue.hpp` and applied in one place.

## 3. Missing data

**NaN alone is never used to express missing data.** Every record carries a
`valid_fields` bitmask and a `quality_flags` bitmask. A field whose validity bit
is clear must not be read, and the accessors enforce it. Records that fail
quality are **retained and flagged**, never silently dropped — the gate asserts
`0 silently dropped records` by reconciling input rows against
`accepted + rejected` with a reason for each rejection.

## 4. Typed coordinate frames

`ICRSDirection`, `ICRSPhaseSpace`, `BarycentricCartesian`,
`GalactocentricCartesian`, `LocalShipFrame` are **distinct types**. An anonymous
`Vec3` never crosses a frame boundary; every transformation is a named function
between two named frame types.

Each frame documents: frame, origin, orientation, epoch, units and velocity
convention. See `frames.hpp`, where those six facts are recorded as data on the
frame type itself and printed into the results metadata.

## 5. Units

Strong types, no implicit conversion from `double`: `Radian`, `Degree`,
`Milliarcsecond`, `MasPerYear`, `Parsec`, `Metre`, `JulianYear`, `Second`,
`KmPerSecond`, `MetrePerSecond`. A public interface never carries a bare
`double distance`.

**Internal dynamical system: parsec, Julian year, and parsec per Julian year.**
Every boundary conversion is a named function with the factor written out:

    1 AU/yr        = 4.740470446 km/s   (149597870.7 km / 31 557 600 s)
    1 pc           = 3.0856775814913673e16 m
    1 pc/yr        = 977792.2 km/s      (derived, see units.hpp)
    c              = 0.3066013938 pc/yr (derived from 1 pc = 3.261563777 ly)

## 6. Epoch propagation

    x(t) = x(t0) + v (t - t0)

**A local kinematic model, not a Galactic orbit integrator.** Rectilinear
barycentric motion ignores Galactic differential rotation and the vertical
potential, both of which matter on 10⁷–10⁸ yr timescales and at kpc distances.
`propagate_rectilinear` returns a validity warning whose threshold is stated in
`astrometry.hpp`, and the gate asserts the warning fires.

## 7. Light-time

Three distinct states, never conflated:

* **catalogue-coordinate state at t** — the model position at coordinate time t;
* **geometric state at observer time** — where the star is *now*, ignoring light travel;
* **apparent retarded state** — where the light arriving now was emitted.

    t_emit = t_observe - |x_star(t_emit) - x_obs(t_observe)| / c

solved by fixed-point iteration, with the iteration count and residual reported.

**No relativistic aberration at C0.** Aberration, Doppler shift and the apparent
sky belong to C1 and are absent rather than approximated.

## 8. Uncertainty

The Gaia 6×6 covariance is stored where available. The transformation to
Cartesian propagates it as `Σ_y = J Σ_x Jᵀ`, with **J computed by centred finite
differences and the method labelled** in the output. HyPAD is not used here: the
available scalar carries 3 derivative directions and this Jacobian needs 6
independent inputs, so a batched scheme would be six evaluations against twelve
for centred differences with no accuracy gain on a smooth analytic map. That
trade is recorded rather than left implicit.

Reported: the positional uncertainty ellipsoid (principal axes), velocity
uncertainty, uncertainty growth with epoch, and an explicit
`missing-covariance` status.

**Derived coordinates are never printed to more precision than the measurements
support.** `significant_digits_for()` derives the digit count from the reported
uncertainty and the report uses it.

## 9. Local database

**Deferred decision, stated openly.** The recommended first implementation is
SQLite. No SQLite amalgamation is vendored in this repository, and adding a
~250 000-line third-party C source file is a dependency decision for the project
owner, not one to slip into a milestone — the same reasoning applied to OTIlib.

So C0 ships:

* a **portable versioned snapshot** written and read by dependency-free C++, with
  the same logical tables (`catalogue_release`, `stellar_object`, `astrometry`,
  `identifier`, `planet`, `provenance`, `import_query`), a schema version in the
  header, and a deterministic FNV-1a checksum over the canonical serialization;
* `tools/snapshot_to_sqlite.py`, which materialises exactly those tables into a
  real `.sqlite` file using Python's built-in `sqlite3` — so the SQLite
  deliverable exists without a new C++ dependency.

Derived columns are stored: Cartesian x, y, z; distance; a coarse spatial cell
index; and the epoch.

**Network access is not part of the simulation runtime.** The workflow is
enforced by directory layout: `tools/` is online, `snapshots/` is a validated
versioned artefact, and the library and gate read only snapshots.

## 10. Known-object validation set

Resolved through catalogue identifiers, never handwritten coordinates. Categories
covered: nearest system; high-proper-motion star; bright nearby star; confirmed
exoplanet host; binary/multiple; a star lacking radial velocity; and a
deliberately poor-parallax source that **must** be rejected.

Until the live import is run these are the synthetic and literature-tagged
records described above, and the gate reports which kind each one is.

## 11. Gate

`C0-CATALOGUE-UNIVERSE`, **28 cases delivered against 25 required**, 0 failures,
0 skips, in both the OTIlib-enabled and OTIlib-disabled configurations (C0 uses
no automatic differentiation, so its case count is identical either way).
Acceptance: 0 failures, 0 silently dropped records, 0 unknown coordinate frames,
0 unit-ambiguous public fields.

Network-dependent import runs only as an optional tool, never in CTest.

The three cases beyond the required 25 exist because writing the optional Python
tools found three real defects that the required cases could not have caught, and
each defect earned a case:

| Case | The defect it now prevents |
| --- | --- |
| `c0.checksum_matches_published_fnv1a_vectors` | `fnv1a()` used a mistyped offset basis (`1469598103934665603` for `14695981039346656037`). It hashed *deterministically*, so the determinism case passed — while making every snapshot unreadable by any independent implementation. Found by `import_gaia.py --verify-format` disagreeing. |
| `c0.flag_bit_values_are_pinned` | The snapshot stores validity and quality masks as **integers**. Both Python tools were written against a guessed bit order and decoded a low-parallax-signal-to-noise source as "non-positive parallax" — a wrong data-quality statement in query output, with the correct reason string sitting right beside it. |
| `c0.no_record_is_silently_dropped` (strengthened) | It reconciled only the injected-duplicate ledger, leaving the *delivered* snapshot — the artefact anyone consumes — unreconciled. It now checks both, and the manifest prints both separately rather than mixing them into a sum that does not balance. |

## 12. Tools and the demonstration

Three artefacts sit outside the gate, in dependency order:

* [`tools/import_gaia.py`](../tools/import_gaia.py) — the optional live import.
  Queries the pinned `gaiadr3.gaia_source` over TAP, enriches identifiers from
  SIMBAD (**identifiers and classifications only, never astrometry**), attaches
  NASA Exoplanet Archive planets, and writes a snapshot whose `data_kind` says
  `Gaia DR3` because it really did come from the archive. Its `serialize()` is
  byte-for-byte the C++ `serialize()`; `--verify-format` proves that against a
  snapshot the gate wrote, and `--print-query` prints the ADQL, both without
  touching a network.
* [`tools/snapshot_to_sqlite.py`](../tools/snapshot_to_sqlite.py) — the SQLite
  materialisation. **Verifies the snapshot checksum before building anything**,
  because a database built from an unverified file gives every downstream query
  an unfalsifiable provenance claim. Stores the snapshot's checksum and
  `data_kind` in the `catalogue_release` table so a query result can always be
  traced back to its input, and decodes the bitmasks into text *alongside* the
  raw integers, never instead of them.
* `c0_demo` — the offline demonstration, run by CTest as `programc0.demo`. Reads
  only the snapshot file. Reports the barycentric origin with its six frame
  facts, the nearest systems at the reference epoch, positions propagated to
  J12016, exoplanet hosts, and candidates ranked by distance — stating in its own
  output that a distance ranking is not a route, and how many records are
  retained but not *placeable* because they have no defensible distance.

The ADQL text is committed in **one** canonical form, duplicated in
`fixtures::adql_for_radius()` and `import_gaia.py`, and the manifest publishes it
under the heading **"NOT ISSUED: this milestone performed no retrieval"**. The two
copies diverged once — the C++ text omitted the proper-motion correlation columns
the importer needs — which would have published a query that, if run, produced a
snapshot with `kQualityMissingCovariance` on every record for no physical reason.

---

## After C0 — documented, not implemented

Fully specified in [`c1_c2_c3_forward_plan.md`](c1_c2_c3_forward_plan.md),
including per-milestone gate sketches and the interfaces C0 must **not** be
altered to accommodate. In brief:

**C1 relativistic astrogation** — ship worldline, proper time, acceleration
profiles, light-time navigation, Doppler shift, stellar aberration, apparent sky,
and interception of a destination's *worldline* rather than travel to its stale
catalogue position.

**C2 effective-metric voyage** — Program A trajectory adapter, moving bubble
centre, local metric and curvature, target-star worldline intersection, parameter
sensitivities, causal and horizon diagnostics. This is where the dependency on
Programs A and B is created, and it runs **one way only**.

**C3 virtual universe and cockpit** — 3D stellar visualization,
observer-dependent sky, route playback, instrument display, exact versus
restricted Program B propagation, uncertainty visualization.
