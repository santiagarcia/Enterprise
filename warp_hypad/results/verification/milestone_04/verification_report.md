# milestone_04 verification report

**Covariant curvature: Christoffel symbols, Riemann and Ricci tensors, curvature invariants, and the temporal-derivative policy they require**

| | |
|---|---|
| Generated (UTC) | 2026-07-28T17:29:02Z |
| Git commit | `bd262fb99592277ec844c2873831a362e4f07842` **(working tree dirty)** |
| Compiler | MSVC 19.44.35228.0 |
| Build type | RelWithDebInfo |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 43 |
| Passed | 43 |
| Failed | 0 |
| Skipped | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | A | curvature.minkowski_pipeline_ok | compute_curvature accepts a genuinely static metric | n/a | n/a | ok | ok | - | - | exact | 4.005400000000e-03 |
| PASS | A | curvature.minkowski_christoffel_zero | Gamma^rho_{mu nu} = 0 for a constant metric | 9x9x9 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 4.031700000000e-03 |
| PASS | A | curvature.minkowski_riemann_zero | R^rho_{sigma mu nu} = 0 in flat Cartesian coordinates | 9x9x9 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 4.284500000000e-03 |
| PASS | A | curvature.flat_spherical_connection_nonzero | Gamma != 0 although the spacetime is flat | n/a | n/a | max\|Gamma\| > 0.5 | 5.000000 | - | - | exact | 3.113560000000e-02 |
| PASS | B | curvature.flat_spherical_christoffel_closed_form | Gamma^r_{th th} = -r, Gamma^th_{r th} = 1/r, Gamma^ph_{th ph} = cot(th), ... | 33x33x5 | r in [1,5], th in [0.5, pi-0.5] | 0 | 4.627884422951e-04 | 4.627884422951e-04 | - | 2.000000000000e-03 | 4.122000000000e-04 |
| PASS | B | curvature.flat_spherical_christoffel_closed_form_boundary | same, on the one-sided shell | 33x33x5 | shell depth 2 | 0 | 1.060977397394e-03 | 1.060977397394e-03 | - | 5.000000000000e-02 | 4.278000000000e-04 |
| PASS | B | curvature.flat_spherical_riemann_zero | R^rho_{sigma mu nu} = 0 for flat space in ANY chart | 33x33x5 | r in [1,5] | 0 | 6.464646464650e-03 | 6.464646464650e-03 | - | 1.500000000000e-02 | 1.691300000000e-03 |
| PASS | A | curvature.alcubierre_pipeline_ok | the time-dependent path produces a result | n/a | n/a | ok | ok | - | - | exact | 5.646530000000e-02 |
| PASS | A | curvature.christoffel_lower_pair_symmetry | Gamma^rho_{mu nu} = Gamma^rho_{nu mu} (torsion free) | 17x17x17 | v=0.5 R=2 sigma=1 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-15 | 1.158000000000e-04 |
| PASS | A | curvature.metric_compatibility | nabla_lambda g_{mu nu} = 0 for the Levi-Civita connection | 17x17x17 | v=0.5 R=2 sigma=1 | 0 | 3.252606517457e-17 | 3.252606517457e-17 | - | 1.000000000000e-14 | 4.013000000000e-04 |
| PASS | A | curvature.riemann_antisymmetry_last_pair | R^rho_{sigma mu nu} = -R^rho_{sigma nu mu} | 17x17x17 | v=0.5 R=2 sigma=1 | 0 | 1.040834085586e-17 | 1.040834085586e-17 | - | 1.000000000000e-15 | 2.845300000000e-03 |
| PASS | A | curvature.riemann_first_bianchi | R^rho_{[sigma mu nu]} = 0 | 17x17x17 | v=0.5 R=2 sigma=1 | 0 | 1.393199791644e-17 | 1.393199791644e-17 | - | 1.000000000000e-14 | 2.849200000000e-03 |
| PASS | A | curvature.riemann_antisymmetry_first_pair | R_{rho sigma mu nu} = -R_{sigma rho mu nu} (indices lowered) | 17x17x17 | v=0.5 R=2 sigma=1 | 0 | 7.360073493998e-03 | 7.360073493998e-03 | - | 2.000000000000e-02 | 5.662600000000e-03 |
| PASS | A | curvature.riemann_pair_exchange | R_{rho sigma mu nu} = R_{mu nu rho sigma} | 17x17x17 | v=0.5 R=2 sigma=1 | 0 | 5.072302569248e-03 | 5.072302569248e-03 | - | 2.000000000000e-02 | 5.670600000000e-03 |
| PASS | A | curvature.ricci_symmetry | R_{mu nu} = R_{nu mu} | 17x17x17 | v=0.5 R=2 sigma=1 | 0 | 4.744600297774e-03 | 4.744600297774e-03 | - | 2.000000000000e-02 | 5.677100000000e-03 |
| PASS | B | curvature.schwarzschild_ricci_flat | R_{mu nu} = 0 (vacuum Einstein equation) | 33x33x5 | M=1, r in [4,6] (horizon at r=2) | 0 | 1.670741694639e-03 | 1.670741694639e-03 | - | 4.000000000000e-03 | 2.515000000000e-04 |
| PASS | B | curvature.schwarzschild_ricci_scalar_zero | R = 0 | 33x33x5 | M=1, r in [4,6] (horizon at r=2) | 0 | 6.221584040016e-05 | 6.221584040016e-05 | - | 4.000000000000e-03 | 3.410000000000e-05 |
| PASS | B | curvature.schwarzschild_kretschmann | R_{abcd} R^{abcd} = 48 M^2 / r^6 | 33x33x5 | M=1, r in [4,6] (horizon at r=2) | 0 | 1.781993460093e-03 | 1.781993460093e-03 | - | 2.000000000000e-03 | 5.820000000000e-05 |
| PASS | A | curvature.chart_connections_differ | the two charts describe the same geometry with different Gamma^rho_{mu nu} | n/a | n/a | max\|Gamma_schw - Gamma_iso\| > 0.02 | 0.049190 | - | - | exact | 1.550000000000e-05 |
| PASS | A | curvature.kretschmann_chart_independent | K is a scalar: its value at areal radius 5 does not depend on the chart | 33x33x5 | M=1, areal radius exactly 5 | 3.067950364674e-03 | 3.067841582144e-03 | 1.087825297797e-07 | 3.545772155646e-05 | 1.000000000000e-03 | 2.060000000000e-05 |
| PASS | B | curvature.kretschmann_isotropic_closed_form | K(rho) = 48 M^2 / r(rho)^6 in isotropic coordinates | 33x33x5 | M=1, rho midpoint | 3.072000000000e-03 | 3.067841582144e-03 | 4.158417855774e-06 | 1.353651645760e-03 | 1.000000000000e-03 | 3.910000000000e-05 |
| PASS | B | curvature.desitter_ricci_scalar | R = 12/L^2 for de Sitter, everywhere | 33x33x5 | L=10, r in [1,5] (horizon at r=10) | 0 | 5.221851117329e-03 | 5.221851117329e-03 | - | 1.200000000000e-02 | 4.040000000000e-05 |
| PASS | A | curvature.desitter_ricci_scalar_sign | R > 0 for de Sitter under the MTW (+,+,+) conventions of docs/conventions.md section 6 | n/a | n/a | R > 0 | 3.000124 | - | - | exact | 5.410000000000e-05 |
| PASS | B | curvature.desitter_ricci_proportional_to_metric | R_{mu nu} = (3/L^2) g_{mu nu} | 33x33x5 | L=10, r in [1,5] (horizon at r=10) | 0 | 2.001836755930e-02 | 2.001836755930e-02 | - | 4.500000000000e-02 | 1.407000000000e-04 |
| PASS | B | curvature.desitter_kretschmann | R_{abcd} R^{abcd} = 24/L^4 | 33x33x5 | L=10, r in [1,5] (horizon at r=10) | 0 | 1.070973666375e-02 | 1.070973666375e-02 | - | 2.500000000000e-02 | 3.360000000000e-05 |
| PASS | B | curvature.desitter_ricci_square | R_{mu nu} R^{mu nu} = 36/L^4 | 33x33x5 | L=10, r in [1,5] (horizon at r=10) | 0 | 1.060657088569e-02 | 1.060657088569e-02 | - | 2.500000000000e-02 | 3.370000000000e-05 |
| PASS | A | curvature.frozen_time_refused | a Refused temporal policy yields no curvature at all | n/a | n/a | ok == false | temporal policy is Refused: TimeTreatment::FrozenTimeDiagnostic discards d_0 g -- curvature built from a discarded temporal derivative is not the curvature of any spacetime | - | - | exact | 1.800000000000e-06 |
| PASS | A | curvature.false_static_claim_refused | ExactlyStatic is verified against the sampled metric before use | n/a | n/a | ok == false | ExactlyStatic was claimed (Alcubierre (a deliberately false claim) has no t dependence in this chart; d_0 g vanishes identically rather than being discarded) but the metric changed by 0.119730 between t and t+h | - | - | exact | 7.373000000000e-04 |
| PASS | A | curvature.rigid_translation_identity | d_0 g = -v_s d_x g exactly, for x_s(t) = v_s t | 17x17x17 | v=0.5 R=2 sigma=1, h_t = h_x/v_s | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-15 | 4.247150000000e-02 |
| PASS | A | curvature.frozen_time_would_be_materially_wrong | discarding d_0 g changes the Ricci scalar by an amount comparable to the Ricci scalar itself | n/a | n/a | max\|R_frozen - R\| > 0.1 max\|R\| | max\|R_frozen - R\| = 0.151391, max\|R\| = 0.114287 | - | - | exact | 6.850860000000e-02 |
| PASS | A | curvature.alcubierre_zero_velocity_flat | v_s = 0 gives g = eta exactly, hence R^rho_{sigma mu nu} = 0 | 17x17x17 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 4.383000000000e-02 |
| PASS | A | curvature.alcubierre_reference_point_on_grid | every reference coordinate is an exact grid coordinate | n/a | n/a | found on grid | found on grid | - | - | exact | 3.590424000000e-01 |
| PASS | B | curvature.alcubierre_ricci_scalar_vs_symbolic | R against exact symbolic differentiation of the closed form | 33x33x33 | v=0.5 R=2 sigma=1, t=0 | 0 | 1.681360865602e-02 | 1.681360865602e-02 | - | 4.000000000000e-02 | 3.590449000000e-01 |
| PASS | B | curvature.alcubierre_ricci_tensor_vs_symbolic | R_{mu nu} componentwise against the symbolic reference | 33x33x33 | v=0.5 R=2 sigma=1, t=0 | 0 | 9.806013443099e-03 | 9.806013443099e-03 | - | 2.500000000000e-02 | 3.590637000000e-01 |
| PASS | B | curvature.alcubierre_kretschmann_vs_symbolic | R_{abcd} R^{abcd} against the symbolic reference | 33x33x33 | v=0.5 R=2 sigma=1, t=0 | 0 | 4.873350692668e-03 | 4.873350692668e-03 | - | 2.500000000000e-02 | 3.590658000000e-01 |
| PASS | B | curvature.alcubierre_ricci_square_vs_symbolic | R_{mu nu} R^{mu nu} against the symbolic reference | 33x33x33 | v=0.5 R=2 sigma=1, t=0 | 0 | 2.477560362196e-03 | 2.477560362196e-03 | - | 2.500000000000e-02 | 3.590671000000e-01 |
| PASS | A | curvature.alcubierre_centre_finite | curvature is finite at r_s = 0 | n/a | n/a | finite | finite | - | - | exact | 9.000000000000e-07 |
| PASS | A | curvature.alcubierre_axial_symmetry | R(x,y,z) = R(x,z,y): the bubble is axisymmetric about x | 33x33x33 | v=0.5 R=2 sigma=1, t=0 | 0 | 2.775557561563e-17 | 2.775557561563e-17 | - | 1.000000000000e-14 | 5.910000000000e-05 |
| PASS | A | curvature.alcubierre_asymptotically_flat | curvature at the far corner is small compared with the wall | n/a | n/a | \|R\| at corner < 2% of max\|R\| | corner = 0.000170, max = 0.145296 | - | - | exact | 2.182000000000e-04 |
| PASS | E | curvature.convergence_schwarzschild_kretschmann | K error falls as O(h^2) under refinement | 9,17,33 (theta refined with r) | M=1, evaluated at the pinned midpoint r=5 | 2.000000000000e+00 | 1.984385537974e+00 | 1.561446202635e-02 | 7.807231013177e-03 | 2.500000000000e-01 | 5.414000000000e-02 |
| PASS | E | curvature.convergence_flat_spherical_riemann | the spurious Riemann of flat space falls as O(h^2) | 9,17,33 | r in [1,5] | 2.000000000000e+00 | 1.997062592816e+00 | 2.937407184201e-03 | 1.468703592100e-03 | 3.000000000000e-01 | 4.354890000000e-02 |
| PASS | E | curvature.convergence_ricci_symmetry | the Ricci antisymmetry induced by the discrete product rule falls as O(h^2) | 17,25,33 | v=0.5 R=2 sigma=1 | 2.000000000000e+00 | 2.077223650842e+00 | 7.722365084178e-02 | 3.861182542089e-02 | 4.000000000000e-01 | 5.948733000000e-01 |
| PASS | E | curvature.convergence_alcubierre_ricci_scalar | R converges to the symbolic reference at O(h^2), including the resampled temporal derivative | 9,17,33 | v=0.5 R=2 sigma=1, t=0 | 2.000000000000e+00 | 1.628565302839e+00 | 3.714346971615e-01 | 1.857173485807e-01 | 4.000000000000e-01 | 5.728005000000e-01 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
