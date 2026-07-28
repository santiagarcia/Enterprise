# milestone_06 verification report

**Observer frames and energy conditions: ADM decomposition, orthonormal tetrads, typed observer families, the frame stress-energy decomposition, the exact NEC minimiser on the sphere, WEC/DEC/SEC certification, and the integrated optimization-ready objectives**

| | |
|---|---|
| Generated (UTC) | 2026-07-28T20:55:41Z |
| Git commit | `f9a1921cd92656b1790341831eb61799ef575e67` **(working tree dirty)** |
| Compiler | MSVC 19.44.35228.0 |
| Build type | RelWithDebInfo |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 40 |
| Passed | 40 |
| Failed | 0 |
| Skipped | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | B | frames.adm_lapse_recovery | alpha = 1 exactly for the Alcubierre metric | 17x17x17 | v=0.5 R=2 sigma=1 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 1.333000000000e-04 |
| PASS | B | frames.adm_spatial_metric_recovery | gamma_ij = delta_ij for the Alcubierre metric | 17x17x17 | v=0.5 R=2 sigma=1 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-15 | 1.916000000000e-04 |
| PASS | A | frames.eulerian_normalization | g_{mu nu} n^mu n^nu = -1 | 17x17x17 | v=0.5 R=2 sigma=1 | 0 | 1.110223024625e-16 | 1.110223024625e-16 | - | 1.000000000000e-13 | 1.940000000000e-04 |
| PASS | A | frames.tetrad_orthonormality | g_{mu nu} e_a^mu e_b^nu = eta_{ab} | 17x17x17 | v=0.5 R=2 sigma=1 | 0 | 1.110223024625e-16 | 1.110223024625e-16 | - | 1.000000000000e-13 | 7.322000000000e-04 |
| PASS | A | frames.tetrad_handedness_uniform | the frame determinant does not change sign across the grid | n/a | n/a | one sign throughout | determinant range [1.000000, 1.000000] | - | - | exact | 7.376000000000e-04 |
| PASS | A | frames.tetrad_continuity | adjacent tetrads differ by O(h), not by O(1) | 17x17x17 | v=0.5 R=2 sigma=1 | 0 | 1.197303001732e-01 | 1.197303001732e-01 | - | 2.000000000000e-01 | 9.897000000000e-04 |
| PASS | A | observers.stationary_valid_only_below_zero | the coordinate-stationary observer exists exactly where g_00 < 0 | n/a | n/a | valid for every g_00 < 0 | valid | - | - | exact | 2.369000000000e-03 |
| PASS | A | observers.stationary_rejected_at_null | the coordinate-stationary observer is REJECTED at g_00 = 0 | n/a | n/a | rejected | rejected: d/dt is null | - | - | exact | 2.371600000000e-03 |
| PASS | A | observers.stationary_rejected_when_spacelike | the coordinate-stationary observer is REJECTED for g_00 > 0 | n/a | n/a | rejected | rejected: d/dt is spacelike | - | - | exact | 2.372500000000e-03 |
| PASS | A | observers.eulerian_survives_the_crossing | the Eulerian observer remains timelike through g_00 = 0 and beyond | n/a | n/a | valid at every v_s | valid throughout | - | - | exact | 2.376100000000e-03 |
| PASS | A | observers.rejections_counted | inadmissible observers are counted, not silently normalised | n/a | n/a | some candidates rejected | 7 rejections across the sweep | - | - | exact | 2.378400000000e-03 |
| PASS | A | observers.boost_rejected_at_and_beyond_light_speed | boosted observers require \|v\| < 1 strictly | n/a | n/a | both rejected | boost speed \|v\| = 1.000000 is at or beyond the unit sphere; no such observer exists | - | - | exact | 3.100000000000e-06 |
| PASS | B | fixtures.vacuum | known energy-condition signature of a standard source | n/a | n/a | NEC ok, WEC ok, DEC ok, SEC ok | NEC ok (N_min = 0.000000), WEC ok, DEC ok, SEC ok | - | - | exact | 9.300000000000e-05 |
| PASS | B | fixtures.dust | known energy-condition signature of a standard source | n/a | n/a | NEC ok, WEC ok, DEC ok, SEC ok | NEC ok (N_min = 1.000000), WEC ok, DEC ok, SEC ok | - | - | exact | 9.260000000000e-05 |
| PASS | B | fixtures.radiation | known energy-condition signature of a standard source | n/a | n/a | NEC ok, WEC ok, DEC ok, SEC ok | NEC ok (N_min = 1.333333), WEC ok, DEC ok, SEC ok | - | - | exact | 8.280000000000e-05 |
| PASS | B | fixtures.lambda_positive | known energy-condition signature of a standard source | n/a | n/a | NEC ok, WEC ok, DEC ok, SEC violated | NEC ok (N_min = 0.000000), WEC ok, DEC ok, SEC violated | - | - | exact | 8.430000000000e-05 |
| PASS | B | fixtures.phantom | known energy-condition signature of a standard source | n/a | n/a | NEC violated, WEC violated, DEC violated, SEC violated | NEC violated (N_min = -0.600000), WEC violated, DEC violated, SEC violated | - | - | exact | 8.450000000000e-05 |
| PASS | B | fixtures.anisotropic | known energy-condition signature of a standard source | n/a | n/a | NEC violated, WEC violated, DEC violated, SEC violated | NEC violated (N_min = -0.400000), WEC violated, DEC violated, SEC violated | - | - | exact | 1.431000000000e-04 |
| PASS | B | fixtures.momentum_flux | known energy-condition signature of a standard source | n/a | n/a | NEC violated, WEC violated, DEC violated, SEC violated | NEC violated (N_min = -0.963331), WEC violated, DEC violated, SEC violated | - | - | exact | 9.380000000000e-05 |
| PASS | A | nec.axis_only_search_is_insufficient | the six coordinate directions do not find the NEC minimum | n/a | n/a | axis-only search overstates N_min | axis-only = -0.600000, true minimum = -0.963331 at (-0.832050, -0.554700, 0.000000) | - | - | exact | 5.970000000000e-05 |
| PASS | A | nec.sphere_constraint | the reported direction is a unit vector | n/a (pointwise) | off-axis test tensor | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 2.037000000000e-04 |
| PASS | A | nec.stationarity | (S - mu I) n = q at the minimiser | n/a (pointwise) | off-axis test tensor | 0 | 2.482534153247e-16 | 2.482534153247e-16 | - | 1.000000000000e-08 | 2.082000000000e-04 |
| PASS | D | nec.quadrature_crosscheck | 16384-point deterministic spherical quadrature finds nothing lower | n/a | n/a | quadrature minimum >= analytic minimum | difference = 0.000059 | - | - | exact | 2.112000000000e-04 |
| PASS | D | nec.multistart_crosscheck | 26-seed projected gradient descent finds nothing lower | n/a | n/a | descent minimum >= analytic minimum | difference = -0.000000 | - | - | exact | 2.132000000000e-04 |
| PASS | E | nec.angular_convergence | quadrature approaches the analytic minimum from above at every resolution | n/a | n/a | no quadrature beats the analytic value | confirmed at 64, 256, 1024, 4096, 16384 points | - | - | exact | 6.320000000000e-04 |
| PASS | A | nec.direct_substitution | N(n_worst) recomputed independently equals the reported minimum | n/a (pointwise) | off-axis test tensor | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 0.000000000000e+00 |
| PASS | A | nec.minimizer_is_off_axis | the test tensor's worst direction is not a coordinate axis | n/a | n/a | all components nonzero | (-0.326600, -0.846114, 0.421217) | - | - | exact | 3.200000000000e-06 |
| PASS | B | alcubierre.eulerian_density_vs_closed_form | rho = -(v_s^2/32 pi)[(d_y f)^2 + (d_z f)^2] | 25x25x25 | v=0.5 R=2 sigma=1 | 0 | 1.457184556087e-04 | 1.457184556087e-04 | - | 4.000000000000e-04 | 1.805284800000e+00 |
| PASS | B | alcubierre.density_vanishes_on_axis | rho = 0 on the axis of travel, where y = z = 0 | 25x25x25 | v=0.5 R=2 sigma=1 | 0 | 3.530921807969e-05 | 3.530921807969e-05 | - | 4.000000000000e-04 | 1.805303000000e+00 |
| PASS | A | alcubierre.density_negative_off_axis | rho < 0 in the off-axis wall | n/a | n/a | a negative minimum exists | min rho = -0.000560 | - | - | exact | 1.805314500000e+00 |
| PASS | A | alcubierre.tetrad_residual | orthonormality across the energy field | 25x25x25 | v=0.5 R=2 sigma=1 | 0 | 1.110223024625e-16 | 1.110223024625e-16 | - | 1.000000000000e-12 | 1.805320900000e+00 |
| PASS | B | alcubierre.density_scales_as_velocity_squared | rho ~ v_s^2, so doubling v_s quadruples \|rho\| | 17x17x17 | v_s = 0.4 vs 0.8 | 4.000000000000e+00 | 4.003549280883e+00 | 3.549280882846e-03 | 8.873202207116e-04 | 6.000000000000e-02 | 8.924643000000e-01 |
| PASS | A | alcubierre.nec_violated_and_certified | the minimum NEC value is negative by more than the error envelope | n/a | n/a | certified_violated | N_min = -0.013680, envelope = 0.000400, safety factor 3.000000 | - | - | exact | 1.700597100000e+00 |
| PASS | A | integrals.e_minus_positive | E_- is strictly positive: the inferred source has negative-density regions | n/a | n/a | E_- > 0 | E_- = 0.044084 | - | - | exact | 1.149986430000e+01 |
| PASS | A | integrals.e_plus_and_e_minus_reported_separately | E_+ and E_- are both reported, not folded into E_net | n/a | n/a | E_+, E_- and E_abs separately available | E_+ = 0.000289, E_- = 0.044084, E_abs = 0.044373 | - | - | exact | 1.149987050000e+01 |
| PASS | A | integrals.nec_violation_positive | J_NEC > 0: the integrated NEC violation is nonzero | n/a | n/a | J_NEC > 0 | J_NEC = 1.424342 | - | - | exact | 1.149987420000e+01 |
| PASS | A | integrals.not_boundary_dominated | the energy integral is not dominated by the outer shell | n/a | n/a | boundary share < 35% | boundary share = 0.016324 | - | - | exact | 1.149987550000e+01 |
| PASS | A | integrals.exact_and_smooth_penalties_distinct | the softplus surrogate is not the exact objective | n/a | n/a | the two differ | exact E_- = 0.044084, softplus E_- = 2.076494 (width 200.000000) | - | - | exact | 1.000000000000e-07 |
| PASS | E | integrals.resolution_convergence | E_- changes by less than 20% between the two finest grids | n/a | n/a | relative change < 20% | E_- = 0.028093, 0.029822, 0.031099; finest relative change = 0.041061 | - | - | exact | 6.898226400000e+00 |
| PASS | E | integrals.domain_convergence | E_- is insensitive to enlarging the domain | n/a | n/a | relative change < 25% | E_- at half-width 4 -> 5: 0.030645 -> 0.030021, relative change 0.020795 | - | - | exact | 5.616658900000e+00 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
