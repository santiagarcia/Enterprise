# program_c0_real verification report

**C0-REAL-CATALOGUE: the committed Gaia DR3 snapshot against independently published astronomy**

| | |
|---|---|
| Generated (UTC) | 2026-07-26T23:28:31Z |
| Git commit | `086b3b4ed7251053261dc86d26859cb1d1eceb21` **(working tree dirty)** |
| Compiler | MSVC 19.44.35223.0 |
| Build type | RelWithDebInfo |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 15 |
| Passed | 15 |
| Failed | 0 |
| Skipped | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | A | c0real.snapshot_is_retrieved_gaia_dr3 | the snapshot declares itself Gaia DR3, names its table and epoch, carries a retrieval timestamp, and its checksum verifies | n/a | n/a | data_kind = Gaia DR3, checksum verified, timestamp present | data_kind = Gaia DR3, table gaiadr3.gaia_source, retrieved 2026-07-26T23:24:32Z, checksum verified | - | - | exact | 6.000000000000e-07 |
| PASS | A | c0real.astrometry.5853498713190525696 | Gaia DR3 5853498713190525696 (NAME Proxima Centauri) matches its published parallax to 1e-6 mas and its published proper motion to 1e-3 mas/yr | n/a | n/a | the published Gaia DR3 values | parallax 768.066539 mas, pm* (-3781.741, 769.465) mas/yr, name "NAME Proxima Centauri" | - | - | exact | 5.910000000000e-05 |
| PASS | A | c0real.astrometry.4472832130942575872 | Gaia DR3 4472832130942575872 (NAME Barnard's star) matches its published parallax to 1e-6 mas and its published proper motion to 1e-3 mas/yr | n/a | n/a | the published Gaia DR3 values | parallax 546.975940 mas, pm* (-801.551, 10362.394) mas/yr, name "NAME Barnard's star" | - | - | exact | 5.700000000000e-06 |
| PASS | A | c0real.astrometry.3864972938605115520 | Gaia DR3 3864972938605115520 (Wolf  359) matches its published parallax to 1e-6 mas and its published proper motion to 1e-3 mas/yr | n/a | n/a | the published Gaia DR3 values | parallax 415.179416 mas, pm* (-3866.338, -2699.215) mas/yr, name "Wolf  359" | - | - | exact | 4.600000000000e-06 |
| PASS | A | c0real.astrometry.762815470562110464 | Gaia DR3 762815470562110464 (HD  95735) matches its published parallax to 1e-6 mas and its published proper motion to 1e-3 mas/yr | n/a | n/a | the published Gaia DR3 values | parallax 392.752945 mas, pm* (-580.057, -4776.589) mas/yr, name "HD  95735" | - | - | exact | 4.100000000000e-06 |
| PASS | A | c0real.astrometry.2947050466531873024 | Gaia DR3 2947050466531873024 (* alf CMa B) matches its published parallax to 1e-6 mas and its published proper motion to 1e-3 mas/yr | n/a | n/a | the published Gaia DR3 values | parallax 374.489589 mas, pm* (-461.571, -914.520) mas/yr, name "* alf CMa B" | - | - | exact | 6.800000000000e-06 |
| PASS | A | c0real.astrometry.5164707970261890560 | Gaia DR3 5164707970261890560 (* eps Eri) matches its published parallax to 1e-6 mas and its published proper motion to 1e-3 mas/yr | n/a | n/a | the published Gaia DR3 values | parallax 310.577293 mas, pm* (-974.758, 20.876) mas/yr, name "* eps Eri" | - | - | exact | 4.200000000000e-06 |
| PASS | B | c0real.proxima_distance_in_light_years | the nearest star's distance, from Gaia DR3 parallax 768.066539 mas via inverse parallax, against the published 4.2465 ly | single source | source_id 5853498713190525696 | 4.246500000000e+00 | 4.246459923405e+00 | 4.007659507099e-05 | 9.437559183090e-06 | 1.000000000000e-04 | 8.300000000000e-06 |
| PASS | B | c0real.proxima_light_travel_time | light travel time to the nearest star, through the unit layer, against the published 4.2465 yr | single source | source_id 5853498713190525696 | 4.246500000000e+00 | 4.246459923405e+00 | 4.007659507099e-05 | 9.437559183090e-06 | 1.000000000000e-03 | 3.000000000000e-07 |
| PASS | B | c0real.barnard_total_proper_motion | the largest known stellar proper motion, against the published 10.3934 arcsec/yr | single source | source_id 4472832130942575872 | 1.039340000000e+04 | 1.039334872227e+04 | 5.127772605556e-02 | 4.933681572494e-06 | 1.000000000000e-01 | 2.000000000000e-06 |
| PASS | B | c0real.barnard_transverse_velocity | transverse velocity from proper motion and distance, against the published ~90.1 km/s | single source | source_id 4472832130942575872 | 9.010000000000e+01 | 9.007592300930e+01 | 2.407699070493e-02 | 2.672252020525e-04 | 1.000000000000e-01 | 4.000000000000e-07 |
| PASS | A | c0real.barnard_is_the_fastest_in_the_volume | no star in the retrieved volume has a larger proper motion than Barnard's star | n/a | n/a | Barnard's star holds the record | NAME Barnard's star at 10393.348722 mas/yr | - | - | exact | 1.100000000000e-06 |
| PASS | A | c0real.every_record_lies_inside_the_requested_volume | the parallax selection really does bound the volume: every accepted distance is within the requested radius | n/a | n/a | all records placeable and within 5.000000 pc | 50 of 50 placeable, farthest 4.967088 pc | - | - | exact | 4.000000000000e-07 |
| PASS | A | c0real.every_distance_names_its_estimator | no record in the real snapshot carries a distance without recording how it was obtained | n/a | n/a | every record has a named estimator | 50 of 50 with an estimator; 34 with a radial velocity, 50 with a covariance, 50 resolved to a name by SIMBAD | - | - | exact | 1.000000000000e-06 |
| PASS | A | c0real.ingestion_ledger_balances_on_real_data | rows retrieved equals accepted plus rejected plus duplicates collapsed, on the real import | n/a | n/a | the ledger balances | 50 objects written from the retrieved rows | - | - | exact | 2.100000000000e-06 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
