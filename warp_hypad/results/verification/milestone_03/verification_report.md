# milestone_03 verification report

**Metric generation: Minkowski and frozen-time Alcubierre, with determinant, inverse, signature and HyPAD parameter sensitivities**

| | |
|---|---|
| Generated (UTC) | 2026-07-26T14:25:39Z |
| Git commit | `611053f22b136ad570115f0535ab02e677e46fef` |
| Compiler | MSVC 19.44.35223.0 |
| Build type | RelWithDebInfo |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 31 |
| Passed | 31 |
| Failed | 0 |
| Skipped | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | A | metric.minkowski.exact_values | g_{mu nu} == diag(-1,1,1,1) at every grid point | 21x21x21 | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 9.238000000000e-04 |
| PASS | A | metric.minkowski.determinant | det g == -1 | 21x21x21 | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 9.943000000000e-04 |
| PASS | A | metric.minkowski.symmetry | g_{mu nu} == g_{nu mu} | 21x21x21 | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 9.985000000000e-04 |
| PASS | B | metric.shape.interior_value | f(0) == 1 exactly, by construction of the numerator | n/a | R=3, sigma=2 | 1.000000000000e+00 | 1.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-15 | 3.100000000000e-06 |
| PASS | B | metric.shape.exterior_decay | f(r) -> 0 as r -> infinity (asymptotic flatness of the shape function) | n/a | r = 60, R = 3, sigma = 2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 1.870000000000e-05 |
| PASS | A | metric.shape.monotone_decreasing | f is non-increasing in r on [0, 40] | n/a | 2000 samples | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-15 | 1.079000000000e-04 |
| PASS | A | metric.shape.bounded | 0 <= f <= 1 everywhere | n/a | 2000 samples | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-15 | 1.100000000000e-04 |
| PASS | A | metric.shape.radius_shifts_wall | increasing R moves the f = 1/2 wall outward | n/a | n/a | wall radius strictly increasing with R | R=2 -> 2.000335, R=4 -> 4.000000, R=6 -> 6.000000 | - | - | exact | 2.600000000000e-05 |
| PASS | A | metric.shape.sigma_steepens_wall | increasing sigma increases \|df/dr\| at the wall | n/a | n/a | wall slope strictly increasing with sigma | sigma=1 -> 0.502460, 2 -> 1.000012, 4 -> 2.000000 | - | - | exact | 6.600000000000e-05 |
| PASS | A | metric.shape.velocity_independent | f depends on R and sigma only, never on v | n/a | v = 0.1 vs 0.9 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 4.700000000000e-06 |
| PASS | A | metric.alcubierre.velocity_affects_expected_components | changing v leaves the purely spatial block g_{ij} unchanged; only g_00 and g_0x respond | n/a | v = 0.1 vs 0.9 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 4.300000000000e-06 |
| PASS | B | metric.alcubierre.determinant_exact | det g == -1 identically, for every parameter value and every point (the (t,x) block has determinant (-1 + v^2 f^2) - (v f)^2 = -1) | 41x41x41 | v=0.2, R=3, sigma=2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-15 | 3.652860000000e-02 |
| PASS | A | metric.alcubierre.symmetry | g_{mu nu} == g_{nu mu} everywhere | 41x41x41 | v=0.2, R=3, sigma=2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 3.656920000000e-02 |
| PASS | A | metric.alcubierre.inverse_identity | g_{mu alpha} g^{alpha nu} == delta_mu^nu at every point | 41x41x41 | v=0.2, R=3, sigma=2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 3.657340000000e-02 |
| PASS | A | metric.alcubierre.lorentzian_signature | every point has exactly one negative and three positive metric eigenvalues | n/a | n/a | 0 points with wrong signature | 0 of 68921 points | - | - | exact | 3.660780000000e-02 |
| PASS | A | metric.alcubierre.all_points_valid | no point is non-finite, asymmetric, singular, ill-conditioned or of wrong signature | n/a | n/a | 0 invalid points | 0 invalid | - | - | exact | 3.661400000000e-02 |
| PASS | B | metric.alcubierre.inverse_vs_analytic | the numerical adjugate inverse agrees with the closed-form inverse g^00=-1, g^01=-vf, g^11=1-v^2f^2 | 21x21x21 | v=0.2, R=3, sigma=2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 2.071600000000e-03 |
| PASS | B | metric.alcubierre.asymptotic_flatness | g -> eta far from the bubble | n/a | \|x\| = 60, R = 3, sigma = 2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 8.000000000000e-07 |
| PASS | A | metric.alcubierre.zero_velocity_is_minkowski | v = 0 gives exactly the Minkowski metric regardless of R, sigma | n/a | v = 0 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 4.900000000000e-06 |
| PASS | A | metric.signature_analyser_works | the eigenvalue signature classifier gives (4+,0-) for the identity, (3+,1-) for Minkowski, detects a zero eigenvalue, and detects a second negative eigenvalue | n/a | n/a | all four classified correctly | all four correct | - | - | exact | 4.500000000000e-06 |
| PASS | A | metric.signature_rotation_invariant | eigenvalues of R^T eta R equal those of eta | n/a | rotation by 0.7 rad in the x-y plane | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 6.700000000000e-06 |
| PASS | A | metric.generic.test_scalar_bit_exact | the metric evaluates bit-identically with double and a custom scalar type | n/a | v=0.2, R=3, sigma=2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 5.000000000000e-06 |
| PASS | D | metric.generic.long_double_agrees | long double agrees with double to double precision | n/a | v=0.2, R=3, sigma=2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 8.800000000000e-06 |
| PASS | A | metric.frozen_time_contract | in the frozen-time formulation the metric is independent of t, exactly | n/a | t = 0, 1, -5, 1000 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.100000000000e-06 |
| PASS | B | metric.sensitivity.shape_function | df/dR and df/dsigma match the hand-derived closed form, and df/dv is exactly zero, at 6 sample points | n/a | v=0.2, R=3, sigma=2 | 0 | 1.110223024625e-16 | 1.110223024625e-16 | - | 1.000000000000e-13 | 5.450000000000e-05 |
| PASS | B | metric.sensitivity.g00_vs_analytic | dg_00/dp == 2 v f^2, 2 v^2 f f_R, 2 v^2 f f_sigma | n/a | 6 points x 3 parameters | 0 | 1.110223024625e-16 | 1.110223024625e-16 | - | 1.000000000000e-13 | 3.470000000000e-04 |
| PASS | B | metric.sensitivity.g0x_vs_analytic | dg_0x/dp == -f, -v f_R, -v f_sigma, and the stored transpose agrees | n/a | 6 points x 3 parameters | 0 | 1.110223024625e-16 | 1.110223024625e-16 | - | 1.000000000000e-13 | 3.509000000000e-04 |
| PASS | A | metric.sensitivity.spatial_block_insensitive | the purely spatial block g_{ij} is the constant 1, so every sensitivity is exactly zero | n/a | 6 points x 3 parameters x 9 components | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 3.538000000000e-04 |
| PASS | A | metric.sensitivity.simultaneous_vs_separate | seeding all three parameters at once gives the same metric sensitivities as three single-seed evaluations | n/a | 6 points x 16 components | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.361000000000e-04 |
| PASS | A | metric.sensitivity.zero_seed_matches_double | with zero seeds the OTI metric real parts equal the double metric to within a few ulp | n/a | 6 points | 0 | 2.775557561563e-17 | 2.775557561563e-17 | - | 1.000000000000e-15 | 3.450000000000e-05 |
| PASS | D | metric.sensitivity.vs_finite_difference | for every point and parameter, the best centered finite-difference metric sensitivity agrees with the HyPAD value | n/a | 6 points x 3 parameters x 6 step sizes | 0 | 3.116662736531e-05 | 3.116662736531e-05 | - | 1.000000000000e-03 | 1.513200000000e-03 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
