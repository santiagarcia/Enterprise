# milestone_05 verification report

**Einstein tensor and inferred stress-energy: G = R_{mu nu} - g R / 2, raw and symmetric, its trace and square, the contracted Bianchi residual, cosmological-constant treatment and the units layer**

| | |
|---|---|
| Generated (UTC) | 2026-07-28T20:32:19Z |
| Git commit | `3c81091981a06dfcdae03c2cda932a2a53cd7ef4` **(working tree dirty)** |
| Compiler | MSVC 19.44.35228.0 |
| Build type | RelWithDebInfo |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 38 |
| Passed | 38 |
| Failed | 0 |
| Skipped | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | A | einstein.minkowski_vanishes | G_{mu nu} = 0 in flat Cartesian coordinates | 9x9x9 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 7.104100000000e-03 |
| PASS | B | einstein.flat_spherical_vanishes | G_{mu nu} = 0 for flat space in ANY chart, despite Gamma != 0 | 33x33x5 | r in [1,5], theta in [1,pi-1] | 0 | 1.073681501689e-02 | 1.073681501689e-02 | - | 2.000000000000e-02 | 3.892280000000e-02 |
| PASS | B | einstein.schwarzschild_vanishes | G_{mu nu} = 0 (vacuum Einstein equation) | 33x33x5 | M=1, r in [4,6] | 0 | 1.149362154495e-03 | 1.149362154495e-03 | - | 4.000000000000e-03 | 3.257510000000e-02 |
| PASS | B | einstein.schwarzschild_still_curved | K = 48 M^2 / r^6 remains nonzero where G vanishes | 33x33x5 | M=1, r=5 | 3.072000000000e-03 | 3.067950364674e-03 | 4.049635325995e-06 | 1.318240666014e-03 | 1.000000000000e-03 | 3.259890000000e-02 |
| PASS | A | einstein.square_chart_independent | G_{mu nu} G^{mu nu} is a scalar and agrees between charts | 33x33x5 | M=1, areal radius exactly 5 | 0 | 5.603592209536e-11 | 5.603592209536e-11 | - | 1.000000000000e-05 | 8.000000000000e-06 |
| PASS | B | einstein.desitter_analytic_value | G_{mu nu} = -(3/L^2) g_{mu nu} for static de Sitter | 33x33x5 | L=2, r in [0.5,1.5] | 0 | 1.300655204787e-02 | 1.300655204787e-02 | - | 6.000000000000e-02 | 4.117760000000e-02 |
| PASS | B | einstein.desitter_vacuum_with_lambda | G_{mu nu} + Lambda g_{mu nu} = 0 with Lambda = 3/L^2 | 33x33x5 | L=2, Lambda=0.75 | 0 | 1.300655204787e-02 | 1.300655204787e-02 | - | 6.000000000000e-02 | 2.408000000000e-04 |
| PASS | A | einstein.desitter_matter_and_lambda_sources_differ | T^matter vanishes for de Sitter while T^Lambda does not | n/a | n/a | T^matter ~ 0 and T^Lambda != 0 | \|T^matter\| = 0.000118, \|T^Lambda\| = 0.039789 | - | - | exact | 1.825800000000e-03 |
| PASS | B | einstein.alcubierre_vs_symbolic | G_{mu nu} componentwise against exact symbolic differentiation | 33x33x33 | v=0.5 R=2 sigma=1, t=0 | 0 | 8.907449976072e-03 | 8.907449976072e-03 | - | 3.000000000000e-02 | 3.993248000000e-01 |
| PASS | A | einstein.g00_is_not_an_energy_density | G_00 is a coordinate component, not an observer-measured scalar | n/a | n/a | recorded as a coordinate component | G_00 = 0.074989 at the bubble centre | - | - | exact | 5.300000000000e-06 |
| PASS | A | einstein.raw_antisymmetry_measured | max \|G_{mu nu} - G_{nu mu}\| on the directly computed tensor | 17x17x17 | v=0.5 R=2 sigma=1 | 0 | 3.104141866397e-03 | 3.104141866397e-03 | - | 2.000000000000e-02 | 4.596690000000e-02 |
| PASS | A | einstein.symmetric_projection_is_symmetric | (1/2)(G + G^T) is exactly symmetric | 17x17x17 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 3.210000000000e-05 |
| PASS | A | einstein.raw_and_symmetric_both_available | the raw tensor is preserved alongside its symmetric projection | n/a | n/a | the two fields differ | max \|raw - symmetric\| = 0.001552 | - | - | exact | 4.690000000000e-05 |
| PASS | A | einstein.trace_identity | G^mu_mu = -R in four dimensions | 17x17x17 | v=0.5 R=2 sigma=1 | 0 | 2.775557561563e-17 | 2.775557561563e-17 | - | 1.000000000000e-13 | 3.130000000000e-05 |
| PASS | A | bianchi.absolute_residual_is_not_enough | the absolute divergence tracks the size of G, not the quality of the identity | n/a | n/a | the geometry with the smaller \|G\| has the smaller absolute residual | Schwarzschild: max\|G\| = 0.031160, \|div G\| = 0.000000, relative = 0.000000;  Alcubierre: max\|G\| = 0.139083, \|div G\| = 0.005993, relative = 0.043088 | - | - | exact | 3.196929000000e-01 |
| PASS | A | bianchi.contracted_residual_absolute | nabla_mu G^{mu nu} = 0, absolute | 21x21x21 | v=0.5 R=2 sigma=1 | 0 | 5.992793649573e-03 | 5.992793649573e-03 | - | 5.000000000000e-02 | 3.197122000000e-01 |
| PASS | A | bianchi.contracted_residual_relative | nabla_mu G^{mu nu} = 0, normalised by max \|G_{mu nu}\| | 21x21x21 | v=0.5 R=2 sigma=1 | 0 | 4.308790461157e-02 | 4.308790461157e-02 | - | 2.000000000000e-01 | 3.197268000000e-01 |
| PASS | E | einstein.antisymmetry_converges | the antisymmetric part of G falls as O(h^2) | 17,25,33 | v=0.5 R=2 sigma=1 | 2.000000000000e+00 | 2.077223650842e+00 | 7.722365084178e-02 | 3.861182542089e-02 | 4.000000000000e-01 | 6.169348000000e-01 |
| PASS | E | bianchi.converges | the contracted Bianchi residual falls under refinement | 17,25,33 | v=0.5 R=2 sigma=1, relative | 2.000000000000e+00 | 1.169243402642e+00 | 8.307565973581e-01 | 4.153782986791e-01 | 9.000000000000e-01 | 1.996490700000e+00 |
| PASS | A | units.geometrized_conversion | T_{mu nu} = G_{mu nu} / (8 pi) in geometrized units | n/a (pointwise) | G = c = 1 | 2.670466304858e-03 | 2.670466304858e-03 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-15 | 8.000000000000e-07 |
| PASS | A | units.si_requires_an_explicit_scale | no SI value is produced without a declared coordinate scale | n/a | n/a | NaN | NaN | - | - | exact | 6.000000000000e-07 |
| PASS | A | units.si_conversion_dimensional_consistency | T^SI = c^4/(8 pi G_N) * G_dimensionless / L0^2 | n/a (pointwise) | L0 = 100 m | 3.231946704832e+37 | 3.231946704832e+37 | 0.000000000000e+00 | 0.000000000000e+00 | 3.231946704832e+31 | 3.300000000000e-06 |
| PASS | A | units.si_scales_as_inverse_length_squared | T^SI scales as 1/L0^2 | n/a (pointwise) | L0 = 100 m vs 50 m | 4.000000000000e+00 | 4.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-09 | 1.000000000000e-07 |
| PASS | A | einstein.frozen_time_refused | no Einstein tensor is produced from a discarded temporal derivative | n/a | n/a | ok == false | curvature was refused: temporal policy is Refused: FrozenTimeDiagnostic discards d_0 g -- curvature built from a discarded temporal derivative is not the curvature of any spacetime | - | - | exact | 5.300000000000e-06 |
| PASS | A | einstein.deterministic | two identical runs produce bit-identical tensors | 9x9x9 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 6.194380000000e-02 |
| PASS | A | einstein.scope_statement_recorded | every A5 report carries its scope limitation | n/a | n/a | scope statement present | A5 computes the geometric Einstein tensor and the stress-energy tensor inferred through Einstein's equation for the selected metric. It does not establish that this source can be generated physically, that energy conditions are satisfied, or that FTL travel is realizable. | - | - | exact | 0.000000000000e+00 |
| PASS | A | lambda.matter_source_identity | 8 pi T^matter = G + Lambda g | 33x33x5 | de Sitter L=2, Lambda=0.75 | 0 | 2.775557561563e-17 | 2.775557561563e-17 | - | 1.000000000000e-14 | 1.871000000000e-04 |
| PASS | A | lambda.effective_lambda_source_identity | 8 pi T^Lambda = -Lambda g | 33x33x5 | de Sitter L=2, Lambda=0.75 | 0 | 2.220446049250e-16 | 2.220446049250e-16 | - | 1.000000000000e-14 | 1.619000000000e-04 |
| PASS | A | lambda.total_source_identity | 8 pi T^total = G | 17x17x5 | de Sitter L=2, Lambda=0.75 | 0 | 2.220446049250e-16 | 2.220446049250e-16 | - | 1.000000000000e-14 | 1.557000000000e-04 |
| PASS | A | lambda.decomposition_closure | T^total = T^matter + T^Lambda | 33x33x5 | de Sitter L=2, Lambda=0.75 | 0 | 1.387778780781e-17 | 1.387778780781e-17 | - | 1.000000000000e-15 | 1.643000000000e-04 |
| PASS | B | lambda.desitter_matter_source_vanishes | T^matter = 0 for de Sitter with Lambda = 3/L^2: it is a VACUUM | 33x33x5 | de Sitter L=2, Lambda=0.75 | 0 | 5.175142627500e-04 | 5.175142627500e-04 | - | 1.500000000000e-03 | 2.460000000000e-05 |
| PASS | A | lambda.desitter_lambda_source_nonzero | T^Lambda != 0 for de Sitter | n/a | n/a | max \|T^Lambda\| > 1e-3 | max \|T^Lambda\| = 0.068209 | - | - | exact | 4.550000000000e-05 |
| PASS | B | lambda.desitter_total_equals_lambda_source | T^total = T^Lambda for de Sitter, because T^matter vanishes | 33x33x5 | de Sitter L=2, Lambda=0.75 | 0 | 5.175142627500e-04 | 5.175142627500e-04 | - | 1.500000000000e-03 | 2.460000000000e-05 |
| PASS | A | lambda.zero_lambda_source_vanishes | T^Lambda = 0 exactly when Lambda = 0 | 13x13x13 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 3.161820000000e-02 |
| PASS | A | lambda.zero_lambda_matter_equals_total | T^matter = T^total = G/(8 pi) when Lambda = 0 | 13x13x13 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 2.450000000000e-05 |
| PASS | A | lambda.labels_match_arithmetic | the recorded convention states the equation actually implemented | n/a | n/a | convention string names the implemented decomposition | G_{mu nu} + Lambda g_{mu nu} = 8 pi T^matter_{mu nu}; T^Lambda = -Lambda g/(8 pi); T^total = T^matter + T^Lambda = G/(8 pi) | - | - | exact | 7.000000000000e-07 |
| PASS | A | lambda.no_ambiguous_policy_remains | the public API exposes no single policy-selected 'the source' field | n/a | n/a | SourceDecomposition with four separately named tensors | lambda = 0.750000, lambda_on_geometric_side = true | - | - | exact | 2.000000000000e-07 |
| PASS | A | lambda.deterministic_records | two identical runs give bit-identical decompositions | 33x33x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 2.035400000000e-03 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
