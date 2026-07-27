# v0_program_a_sweeps verification report

**A-VISUAL-SWEEP-EXPORT: swept effective-metric fields for the Program A animations, and their scientific assertions**

| | |
|---|---|
| Generated (UTC) | 2026-07-27T22:32:06Z |
| Git commit | `893dcfd6c463793a0201d824d20a9d72b7a3bf52` **(working tree dirty)** |
| Compiler | MSVC 19.44.35223.0 |
| Build type | RelWithDebInfo |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 20 |
| Passed | 20 |
| Failed | 0 |
| Skipped | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | A | asweep.translation_center_follows_trajectory | the saved bubble centre equals x_s0 + v_s t at every frame | n/a | n/a | exact agreement with the trajectory | all 24 frames agree to 1e-12 | - | - | exact | 2.900000000000e-06 |
| PASS | A | asweep.lab_frame_field_translates | the laboratory-frame field moves across the domain over the sweep | n/a | n/a | centre displacement greater than one bubble radius | displacement 4 over 8 time units | - | - | exact | 2.000000000000e-07 |
| PASS | A | asweep.comoving_profile_stationary | the comoving shape profile is identical at every frame, so the motion is translation and nothing else | 41 offsets x 24 frames | v = 0.5 | 0 | 4.440892098501e-16 | 4.440892098501e-16 | - | 1.000000000000e-12 | 3.000000000000e-07 |
| PASS | A | asweep.time_values_monotonic | sweep time increases frame to frame | n/a | n/a | strictly increasing | 24 frames from 0 to 8 | - | - | exact | 1.000000000000e-07 |
| PASS | A | asweep.velocity_values_monotonic | sweep velocity increases frame to frame | n/a | n/a | strictly increasing | 24 frames from 0.2 to 1.6 | - | - | exact | 1.000000000000e-07 |
| PASS | A | asweep.velocity_sweep_crosses_vf_unity | the velocity sweep contains a frame where v_s f reaches 1 | n/a | n/a | at least one crossing frame | first crossing at frame 14, v_s = 1.052173913 | - | - | exact | 6.000000000000e-07 |
| PASS | A | asweep.g00_matches_closed_form | the saved g00 equals -(1 - v_s^2 f^2) everywhere in the velocity sweep | 41 x 41 x 24 | v_s in [0.2, 1.6] | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-15 | 2.000000000000e-07 |
| PASS | A | asweep.determinant_is_exactly_minus_one | det g = -1 across every frame of the velocity sweep, including beyond the crossing | 41 x 41 x 24 | v_s in [0.2, 1.6] | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-15 | 3.000000000000e-07 |
| PASS | A | asweep.signature_remains_lorentzian | the metric stays Lorentzian at every sample, including where g00 >= 0 | n/a | n/a | zero signature failures | 0 failures over 40344 samples | - | - | exact | 2.000000000000e-07 |
| PASS | A | asweep.radius_sweep_center_fixed | the bubble centre does not move while the radius is swept | n/a | n/a | centre exactly constant | centre -2 at every frame | - | - | exact | 2.000000000000e-07 |
| PASS | A | asweep.measured_radius_monotonic | the MEASURED half-height radius increases with R | n/a | n/a | strictly increasing | from 0.6040488566 to 2.000000056 | - | - | exact | 2.000000000000e-07 |
| PASS | A | asweep.steepness_sweep_center_fixed | the bubble centre does not move while the steepness is swept | n/a | n/a | centre exactly constant | centre -2 at every frame | - | - | exact | 0.000000000000e+00 |
| PASS | A | asweep.measured_wall_width_monotonic | the MEASURED wall width falls monotonically as sigma rises | n/a | n/a | strictly decreasing | from 0.527618007 to 0.1098612283 | - | - | exact | 3.000000000000e-07 |
| PASS | B | asweep.nominal_radius_fixed_during_steepness_sweep | the measured half-height radius stays at the nominal radius while only sigma changes | 24 frames | R = 1 | 0.000000000000e+00 | 1.768111677111e-02 | 1.768111677111e-02 | - | 5.000000000000e-02 | 4.000000000000e-07 |
| PASS | A | asweep.sensitivity_fields_finite | every exported sensitivity value is finite | n/a | n/a | no NaN or infinity | backend otilib | - | - | exact | 1.000000000000e-07 |
| PASS | A | asweep.oti_matches_analytic_at_center | dg00/dv = 2 v f^2 = 2 v at the bubble centre, where f = 1 | single point | v = 0.5 | 1.000000000000e+00 | 1.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-12 | 9.400000000000e-06 |
| PASS | B | asweep.oti_matches_finite_difference | the OTI sensitivities agree with centred finite differences at a wall point | single point, h = 1e-6 | r_s = R | 0 | 5.454503515523e-11 | 5.454503515523e-11 | - | 1.000000000000e-05 | 5.900000000000e-06 |
| PASS | B | asweep.fixed_lab_fixed_center_chain_rule | dg/dv\|_lab = dg/dv\|_centre - (dg/dx) t, with the trajectory term computed independently | 41 x 41 x 24 | t = 1 | 0 | 2.679518829041e-10 | 2.679518829041e-10 | - | 1.000000000000e-04 | 1.000000000000e-07 |
| PASS | A | asweep.output_checksums_valid | every sweep file exists and has a recorded FNV-1a checksum | n/a | n/a | 6 files with checksums | 6 files written | - | - | exact | 4.000000000000e-07 |
| PASS | A | asweep.deterministic_repeat_export | frame values and grid coordinates are closed-form functions of their index, so a repeat export is byte-identical | n/a | n/a | exact reproduction | all frame values reproduce exactly | - | - | exact | 8.000000000000e-07 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
