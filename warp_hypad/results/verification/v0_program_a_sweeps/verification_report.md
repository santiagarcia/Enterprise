# v0_program_a_sweeps verification report

**A-VISUAL-SWEEP-EXPORT: swept effective-metric fields for the Program A animations, and their scientific assertions**

| | |
|---|---|
| Generated (UTC) | 2026-07-27T23:00:41Z |
| Git commit | `cf39f1cc37fe1590ae329224a39bc0fd61f44e5f` **(working tree dirty)** |
| Compiler | MSVC 19.44.35223.0 |
| Build type | RelWithDebInfo |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 23 |
| Passed | 23 |
| Failed | 0 |
| Skipped | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | A | asweep.translation_center_follows_trajectory | the saved bubble centre equals x_s0 + v_s t at every frame | n/a | n/a | exact agreement with the trajectory | all 24 frames agree to 1e-12 | - | - | exact | 4.000000000000e-06 |
| PASS | A | asweep.lab_frame_field_translates | the laboratory-frame field moves across the domain over the sweep | n/a | n/a | centre displacement greater than one bubble radius | displacement 4 over 8 time units | - | - | exact | 2.000000000000e-07 |
| PASS | A | asweep.comoving_profile_stationary | the comoving shape profile is identical at every frame, so the motion is translation and nothing else | 41 offsets x 24 frames | v = 0.5 | 0 | 4.440892098501e-16 | 4.440892098501e-16 | - | 1.000000000000e-12 | 5.000000000000e-07 |
| PASS | A | asweep.time_values_monotonic | sweep time increases frame to frame | n/a | n/a | strictly increasing | 24 frames from 0 to 8 | - | - | exact | 2.000000000000e-07 |
| PASS | A | asweep.velocity_values_monotonic | sweep velocity increases frame to frame | n/a | n/a | strictly increasing | 24 frames from 0.2 to 1.6 | - | - | exact | 1.000000000000e-07 |
| PASS | A | asweep.velocity_sweep_crosses_vf_unity | the velocity sweep contains a frame where v_s f reaches 1 | n/a | n/a | at least one crossing frame | first crossing at frame 12, v_s = 1 | - | - | exact | 5.000000000000e-07 |
| PASS | A | asweep.g00_matches_closed_form | the saved g00 equals -(1 - v_s^2 f^2) everywhere in the velocity sweep | 41 x 41 x 24 | v_s in [0.2, 1.6] | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-15 | 2.000000000000e-07 |
| PASS | A | asweep.determinant_is_exactly_minus_one | det g = -1 across ALL FOUR sweeps -- translation, velocity, radius and steepness -- including beyond the ergosurface crossing | 41 x 41 x 24 per sweep, 4 sweeps | translation 0, velocity 0, radius 0, steepness 0 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-15 | 9.000000000000e-07 |
| PASS | A | asweep.signature_remains_lorentzian | the metric stays Lorentzian at every sample, including where g00 >= 0 | n/a | n/a | zero signature failures | 0 failures over 40344 samples | - | - | exact | 4.400000000000e-06 |
| PASS | A | asweep.radius_sweep_center_fixed | the bubble centre does not move while the radius is swept | n/a | n/a | centre exactly constant | centre -2 at every frame | - | - | exact | 9.000000000000e-07 |
| PASS | A | asweep.measured_radius_monotonic | the MEASURED half-height radius increases with R | n/a | n/a | strictly increasing | from 0.6040488566 to 2.000000056 | - | - | exact | 3.000000000000e-07 |
| PASS | A | asweep.steepness_sweep_center_fixed | the bubble centre does not move while the steepness is swept | n/a | n/a | centre exactly constant | centre -2 at every frame | - | - | exact | 1.000000000000e-07 |
| PASS | A | asweep.measured_wall_width_monotonic | the MEASURED wall width falls monotonically as sigma rises | n/a | n/a | strictly decreasing | from 0.527618007 to 0.1098612283 | - | - | exact | 6.000000000000e-07 |
| PASS | B | asweep.nominal_radius_fixed_during_steepness_sweep | the measured half-height radius stays at the nominal radius while only sigma changes | 24 frames | R = 1 | 0.000000000000e+00 | 1.768111677111e-02 | 1.768111677111e-02 | - | 5.000000000000e-02 | 1.200000000000e-06 |
| PASS | A | asweep.sensitivity_fields_finite | every exported sensitivity value is finite | n/a | n/a | no NaN or infinity | backend otilib | - | - | exact | 1.000000000000e-07 |
| PASS | A | asweep.oti_matches_analytic_at_center | dg00/dv = 2 v f^2 = 2 v at the bubble centre, where f = 1 | single point | v = 0.5 | 1.000000000000e+00 | 1.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-12 | 1.390000000000e-05 |
| PASS | B | asweep.oti_matches_finite_difference | the OTI sensitivities agree with centred finite differences at a wall point | single point, h = 1e-6 | r_s = R | 0 | 5.454503515523e-11 | 5.454503515523e-11 | - | 1.000000000000e-05 | 9.700000000000e-06 |
| PASS | B | asweep.fixed_lab_fixed_center_chain_rule | dg/dv\|_lab = dg/dv\|_centre - (dg/dx) t, with the trajectory term computed independently | 41 x 41 x 24 | t = 1, central-difference step h = 1e-6, worst at (x = -1.6, y = 0, v_s = 0.9652173913) | 0 | 2.679518829041e-10 | 2.679518829041e-10 | - | 1.000000000000e-07 | 1.200000000000e-06 |
| PASS | A | asweep.exact_v_unity_frame_exists | the velocity sequence contains v_s = 1 exactly, and remains strictly increasing | n/a | n/a | one frame at exactly 1.0, sequence strictly increasing | frame 12, v_s = 1 | - | - | exact | 1.000000000000e-06 |
| PASS | A | asweep.g00_zero_at_center_on_exact_onset | at v_s = 1 the saved g00 is exactly zero at the bubble centre, where f = 1 | grid sample at r_s = 0 | v_s = 1, f(0) = 1 | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 3.000000000000e-07 |
| PASS | A | asweep.spatial_zero_contour_exists_above_unity | above v_s = 1 the saved grid brackets zero, so a nontrivial g00 = 0 contour is supported; at the exact onset the field only touches zero and no contour is supported | n/a | n/a | first bracketing frame strictly after the exact onset | exact onset frame 12, first bracketing frame 13 | - | - | exact | 1.100000000000e-06 |
| PASS | A | asweep.output_checksums_valid | every sweep file exists and has a recorded FNV-1a checksum | n/a | n/a | 6 files with checksums | 6 files written | - | - | exact | 4.000000000000e-07 |
| PASS | A | asweep.deterministic_repeat_export | frame values and grid coordinates are closed-form functions of their index, so a repeat export is byte-identical | n/a | n/a | exact reproduction | all frame values reproduce exactly | - | - | exact | 7.000000000000e-07 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
