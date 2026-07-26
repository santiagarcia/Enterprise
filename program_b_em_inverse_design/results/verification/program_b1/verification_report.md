# program_b1 verification report

**Metric to electromagnetic constitutive response: Plebanski-Tamm forward map, conformal-class inverse, and non-metric detection**

| | |
|---|---|
| Generated (UTC) | 2026-07-26T13:59:54Z |
| Git commit | `018cf5de355b1ab9785cf9f8bfce3046e4929fd6` |
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
| PASS | B | b1.minkowski_is_vacuum | eta_{mu nu} maps to eps = mu = I and zero magnetoelectric coupling | n/a | Heaviside-Lorentz, c = 1 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-15 | 2.640000000000e-05 |
| PASS | B | b1.isotropic_forward | g = diag(-1,a,a,a) gives eps = mu = sqrt(a) I | n/a | a = 0.25, 1, 4, 9 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 4.100000000000e-06 |
| PASS | B | b1.gordon_optical_metric_speed | the static isotropic Gordon metric diag(-1/n^2,1,1,1) has null coordinate speed 1/n | n/a | n = 1, 1.5, 2, 3 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-15 | 5.000000000000e-07 |
| PASS | A | b1.plain_dielectric_not_metric_representable | a non-magnetic dielectric (eps = 4I, mu = I) is rejected as metric-representable, despite having a well-defined Gordon optical metric | n/a | n/a | rejected, impedance mismatch reported | eps != mu: the medium is not impedance matched, so it is birefringent in the general case and no single pseudo-Riemannian metric reproduces full Maxwell propagation in it. An optical (Gordon) metric may still describe rays. | - | - | exact | 4.000000000000e-06 |
| PASS | B | b1.anisotropic_forward | g = diag(-1,2,3,5) gives eps^ii = sqrt(30)/g_ii | n/a | manufactured diagonal metric | 0 | 2.220446049250e-16 | 2.220446049250e-16 | - | 1.000000000000e-14 | 4.000000000000e-07 |
| PASS | A | b1.conformal_invariance | the medium is unchanged under g -> Omega^2 g, for Omega^2 = 0.25, 0.5, 2, 7 | n/a | diag(-1,2,3,5) | 0 | 4.440892098501e-16 | 4.440892098501e-16 | - | 1.000000000000e-13 | 2.900000000000e-06 |
| PASS | A | b1.roundtrip_conformal | g -> (eps,mu) -> ghat recovers the CONFORMAL CLASS of g | n/a | 4 diagonal metrics | 0 | 9.868649107779e-17 | 9.868649107779e-17 | - | 1.000000000000e-13 | 6.500000000000e-06 |
| PASS | A | b1.roundtrip_factor_is_nontrivial | at least one round trip returns a representative with Omega^2 != 1, confirming the conformal test is not a disguised equality test | n/a | n/a | some Omega^2 differs from 1 | largest \|Omega^2 - 1\| = 3.000000 | - | - | exact | 1.000000000000e-05 |
| PASS | A | b1.null_cone_equivalence | the reconstructed metric has the same null coordinate speeds along each axis as the original | n/a | diag(-1,2,3,5) | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 2.500000000000e-06 |
| PASS | A | b1.nonmetric_detected | constitutive tensors that no single pseudo-Riemannian metric can produce -- birefringent, axion/Tellegen, and skewon-like -- are each REPORTED rather than fitted | n/a | n/a | 3 of 3 detected | 3 detected. uniaxial birefringent: not_metric_representable; axion/Tellegen: not_metric_representable; asymmetric eps: not_metric_representable;  | - | - | exact | 1.000000000000e-05 |
| PASS | A | b1.valid_medium_accepted | an impedance-matched symmetric medium with no magnetoelectric coupling is accepted and reconstructed | n/a | n/a | accepted | ok | - | - | exact | 8.000000000000e-07 |
| PASS | B | b1.warp_subluminal_medium_exists | for v_s = 0.2 the Alcubierre metric has an equivalent medium at every sampled point | n/a | n/a | 65 of 65 | 65 of 65 | - | - | exact | 1.280000000000e-05 |
| PASS | B | b1.warp_is_bianisotropic | the warp medium has a nonzero magnetoelectric vector w, i.e. it is bianisotropic, not a simple dielectric | n/a | n/a | w_x != 0 at the bubble wall | w_x = -0.101011 | - | - | exact | 1.630000000000e-05 |
| PASS | A | b1.warp_inverse_reports_unsupported | reconstructing a metric from a magnetoelectric medium is reported as unsupported in B1 rather than approximated | n/a | n/a | MagnetoelectricUnsupported | magnetoelectric_unsupported_in_B1 | - | - | exact | 1.980000000000e-05 |
| PASS | B | b1.superluminal_medium_fails_reported | for v_s = 2 there are points where g_00 >= 0 and NO equivalent medium exists, and the map reports this instead of returning a value | n/a | n/a | some points rejected | 121 of 321 sample points have no equivalent medium; first at x = -3 | - | - | exact | 7.480000000000e-05 |
| PASS | A | b1.generic.test_scalar | the constitutive mapping instantiates with a custom scalar having no ordering, no equality and no implicit double conversion | n/a | diag(-1,2,3,5) | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.000000000000e-06 |
| PASS | D | b1.sensitivity.metric_parameters | d eps_yy / d(v, R, sigma) through the Plebanski map matches centred finite differences | n/a | v=0.2, R=3, sigma=2 at the wall x=3 | 0 | 9.067909012186e-08 | 9.067909012186e-08 | - | 1.000000000000e-06 | 6.100000000000e-05 |
| PASS | B | b1.sensitivity.refractive_index | for eps = mu = n^2 I the reconstruction gives h = n^4 I and dh/dn = 4n^3, both matched by HyPAD | n/a | n = 1.7 | 0 | 1.807812781804e-16 | 1.807812781804e-16 | - | 1.000000000000e-13 | 5.400000000000e-06 |
| PASS | A | b1.sensitivity.conformal_normalisation | d eps / d Omega == 0: the differential form of conformal invariance of the forward map | n/a | Omega = 1.4 on diag(-1,2,3,5) | 0 | 8.881784197001e-16 | 8.881784197001e-16 | - | 1.000000000000e-13 | 4.500000000000e-06 |
| PASS | B | b1.sensitivity.observable_null_speed | the null coordinate speed 1/n^2 of the reconstructed medium, and its derivative -2/n^3, both match HyPAD | n/a | n = 1.7 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 5.300000000000e-06 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
