# program_c0 verification report

**C0-CATALOGUE-UNIVERSE: astronomical data schema, catalogue ingestion, typed coordinate frames, epoch propagation, light time, uncertainty and a deterministic offline snapshot**

| | |
|---|---|
| Generated (UTC) | 2026-07-26T23:38:30Z |
| Git commit | `f7f98dae050799cd41206988dfbefb8e290febb5` |
| Compiler | MSVC 19.44.35223.0 |
| Build type | RelWithDebInfo |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 28 |
| Passed | 28 |
| Failed | 0 |
| Skipped | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | A | c0.fixtures_are_declared_synthetic | every committed fixture is tagged SYNTHETIC or LITERATURE_FROM_MEMORY, and none claims to be retrieved catalogue data | n/a | n/a | no fixture claims a real catalogue | snapshot data_kind = SYNTHETIC; 13 objects, all tagged | - | - | exact | 1.000000000000e-06 |
| PASS | A | c0.unit_conversions_round_trip | every boundary conversion round-trips, and 1 arcsec/yr at 1 pc is 4.740470446 km/s | n/a | 6 conversions plus 2 factor checks | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 2.100000000000e-06 |
| PASS | B | c0.speed_of_light_in_internal_units | 1 pc corresponds to 3.261563777 light-years, so c = 0.3066 pc/yr | n/a | derived from c and the Julian year | 3.261563777000e+00 | 3.261563777167e+00 | 1.674340666113e-10 | 5.133551819286e-11 | 1.000000000000e-06 | 1.430000000000e-05 |
| PASS | A | c0.ra_dec_to_unit_vector | the ICRS unit vector is correct at the cardinal directions, and the tangent basis (p_hat, q_hat) is orthonormal and orthogonal to r_hat | n/a | 5 cardinal points, 5 declinations | 0 | 1.224646799147e-16 | 1.224646799147e-16 | - | 1.000000000000e-15 | 2.200000000000e-06 |
| PASS | B | c0.parallax_distance_high_quality | a parallax of 100 mas with signal-to-noise 200 yields exactly 10 pc by inverse parallax, with the estimator recorded | n/a | n/a | estimator InverseParallaxHighQuality, distance 10 pc | InverseParallaxHighQuality, d = 1.0000e+01 +- 5.0000e-02 pc | - | - | exact | 7.000000000000e-07 |
| PASS | A | c0.poor_parallax_is_rejected_not_inverted | a source with parallax signal-to-noise 2 and one with a NEGATIVE parallax both get estimator None, with distinct reasons, and are RETAINED with a flag rather than deleted | n/a | n/a | both rejected for inverse parallax, both retained | snr=2 flagged low-signal-to-noise; negative parallax flagged non-positive; both still in the catalogue | - | - | exact | 8.400000000000e-06 |
| PASS | B | c0.proper_motion_convention_mu_alpha_star | the catalogued pm_ra is mu_alpha* = mu_alpha cos(dec): at dec = 60 the coordinate rate d(ra)/dt comes out at mu_alpha*/cos(dec), twice the physical rate | n/a | dec = 60 deg, mu_alpha* = 500 mas/yr, 100 yr; measured d(ra)/dt / mu_alpha* = 2.0000e+00 (expected 2) | 0 | 7.834821758479e-08 | 7.834821758479e-08 | - | 1.000000000000e-03 | 6.800000000000e-06 |
| PASS | B | c0.radial_velocity_contribution_is_purely_radial | a star with zero proper motion and positive radial velocity moves only along the line of sight | n/a | n/a | transverse displacement fraction below 1e-12 | transverse/total = 2.4799e-15 | - | - | exact | 1.100000000000e-06 |
| PASS | A | c0.astrometry_cartesian_round_trip | observables -> barycentric Cartesian -> observables is the identity, for every fixture with a defensible distance | n/a | 11 sources, unit vector compared near the pole where right ascension is degenerate | 0 | 5.620824514355e-13 | 5.620824514355e-13 | - | 1.000000000000e-12 | 7.800000000000e-06 |
| PASS | A | c0.stationary_star_does_not_move | a source with zero proper motion and zero radial velocity has exactly zero displacement after 10 000 yr | n/a | 10 000 yr | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.900000000000e-06 |
| PASS | B | c0.high_proper_motion_propagation | a 5000 mas/yr source at 4 pc moves a measurable distance in 1000 yr, and the rectilinear validity warning fires past 1e5 yr while staying silent within it | n/a | n/a | measurable displacement; warning fires only outside the stated validity | 1000 yr displacement 9.6963e-02 pc; warning at 1e6 yr = yes, at 1e4 yr = no | - | - | exact | 5.900000000000e-06 |
| PASS | B | c0.light_time_iteration_converges | the light-time fixed point converges for every nearby source, and for a stationary star at 10 pc gives 32.6 yr | n/a | n/a | converged everywhere, and the analytic value at 10 pc | worst residual 2.2737e-13 yr in at most 5 iterations; 10 pc gives 3.2616e+01 yr against 3.2616e+01 | - | - | exact | 9.200000000000e-06 |
| PASS | B | c0.covariance_transformation | Sigma_y = J Sigma_x J^T reproduces the expected radial uncertainty for a source on the +x axis, sigma_x = d * sigma_parallax/parallax | n/a | n/a | the dominant parallax term recovered within 5 per cent | sigma_position = 5.0000e-02 pc, sigma_velocity = 4.0910e-07 pc/yr, Jacobian by centred finite differences | - | - | exact | 6.700000000000e-06 |
| PASS | B | c0.uncertainty_grows_with_epoch | the propagated position uncertainty grows with the propagation interval, as Sigma_xx + dt(Sigma_xv + Sigma_vx) + dt^2 Sigma_vv | n/a | n/a | strictly growing | ratio at 1e4 yr = 1.0033e+00 | - | - | exact | 1.180000000000e-05 |
| PASS | A | c0.missing_fields_are_flagged_not_nan | a missing radial velocity and a missing covariance are expressed by CLEAR VALIDITY BITS and quality flags, and the covariance path refuses to report an uncertainty it does not have | n/a | n/a | validity bits clear, flags set, covariance refused | missing rv and missing covariance both flagged; covariance reported as UNKNOWN rather than zero | - | - | exact | 9.000000000000e-07 |
| PASS | A | c0.snapshot_schema_round_trip | the snapshot writes and reads back with every identifier, validity bitmask, quality bitmask, estimator and float preserved | n/a | n/a | all fields preserved, schema version 1 | worst float difference 0.0000e+00, 13 objects | - | - | exact | 1.568840000000e-02 |
| PASS | A | c0.snapshot_is_deterministic | serialising twice gives byte-identical output and the checksum in the file matches a recomputation | n/a | n/a | byte-identical and checksum verified | FNV-1a checksum 517333436261633138 verified on read | - | - | exact | 1.570350000000e-02 |
| PASS | A | c0.flag_bit_values_are_pinned | the ValidField and QualityFlag bit values are fixed, because the snapshot stores them as integers that external tools must decode | n/a | n/a | validity bits 1..128 and quality bits 1..64 in declared order | fields 1..128, quality 1..64 | - | - | exact | 1.691100000000e-02 |
| PASS | A | c0.checksum_matches_published_fnv1a_vectors | fnv1a() reproduces the published FNV-1a 64-bit test vectors, so an independent implementation can read what this one writes | n/a | n/a | the three published vectors reproduce exactly | empty -> 14695981039346656037, "a" -> 12638187200555641996, "foobar" -> 9625390261332436968 | - | - | exact | 1.691550000000e-02 |
| PASS | A | c0.provenance_is_retained_through_the_snapshot | every record still carries its source catalogue, release, table, timestamp, licence and original identifier after a write and read | n/a | n/a | provenance and the exact ADQL both survive | release pinned to Gaia DR3; ADQL 564 characters stored beside the result | - | - | exact | 1.693190000000e-02 |
| PASS | A | c0.catalogue_release_is_pinned | the snapshot records the pinned catalogue, table and reference epoch, so a silent migration to a later release fails here | n/a | n/a | Gaia DR3, gaiadr3.gaia_source, J2016.0 | Gaia DR3 / gaiadr3.gaia_source / epoch 2.0160e+03 | - | - | exact | 1.695030000000e-02 |
| PASS | A | c0.offline_queries_after_import | the 5 pc and 50 pc volume queries, the nearest-N query and the exoplanet-host query all run against the LOADED snapshot with no network access | n/a | n/a | all four queries return results, 5 pc a subset of 50 pc | 2 within 5 pc, 10 within 50 pc, 1 exoplanet hosts | - | - | exact | 1.696550000000e-02 |
| PASS | A | c0.identifier_enrichment_and_exoplanet_association | a source is findable by an alias as well as its primary name, and a planet record is associated with its host | n/a | n/a | alias resolves to the host, and the host carries its planet | SYN-EXOHOST found via alias, 1 planet(s) | - | - | exact | 1.698300000000e-02 |
| PASS | A | c0.duplicate_sources_collapse_and_are_counted | the same source_id offered twice yields one record, and the collapse is counted rather than silent | n/a | n/a | one duplicate collapsed, count preserved | 1 collapsed, 13 unique records | - | - | exact | 7.319000000000e-04 |
| PASS | A | c0.no_record_is_silently_dropped | rows in equals accepted plus rejected plus duplicates collapsed, and every rejection carries a reason | n/a | n/a | both ledgers balance and every rejection is explained | delivered 13 in = 11 + 2 + 0; with an injected duplicate 14 in = 11 accepted + 2 rejected + 1 collapsed | - | - | exact | 7.420000000000e-04 |
| PASS | A | c0.every_frame_declares_its_six_facts | all five coordinate frames declare frame, origin, orientation, epoch, units and velocity convention | n/a | n/a | five frames fully described | ICRSDirection, ICRSPhaseSpace, BarycentricCartesian, GalactocentricCartesian, LocalShipFrame | - | - | exact | 3.000000000000e-06 |
| PASS | A | c0.no_unit_ambiguous_public_field | the unit types are not implicitly constructible from double and do not implicitly interconvert | n/a | n/a | all five conversions rejected at compile time | checked with std::is_convertible on Parsec, KmPerSecond, Metre, MasPerYear and Degree | - | - | exact | 8.000000000000e-06 |
| PASS | D | c0.literature_reference_sanity | inverse parallax reproduces the recalled approximate distances of three well-known nearby stars | n/a | 3 stars, values RECALLED not retrieved, tolerance 1 per cent | 0 | 1.980392117639e-04 | 1.980392117639e-04 | - | 1.000000000000e-02 | 1.000000000000e-06 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
