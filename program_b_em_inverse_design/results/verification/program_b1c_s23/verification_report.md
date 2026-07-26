# program_b1c_s23 verification report

**B1C-S23: four-dimensional Hodge duals, Bel extraction from Schuster & Visser (A7a)-(A7d), spatial projector and block spatiality**

| | |
|---|---|
| Generated (UTC) | 2026-07-26T22:17:39Z |
| Git commit | `9bf1028688994b1fd2c9718f519998ef2ea1d7dc` |
| Compiler | MSVC 19.44.35223.0 |
| Build type | RelWithDebInfo |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 42 |
| Passed | 42 |
| Failed | 0 |
| Skipped | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | A | b1c.dual.levi_civita_normalisation | eps_{0123} = +sqrt(-det g) and eps^{0123} = -1/sqrt(-det g), both totally antisymmetric, for Minkowski AND a curved non-diagonal lab | n/a | project convention, section 11 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 1.700000000000e-05 |
| PASS | A | b1c.dual.contraction_identity | eps_{abcd} eps^{cdef} = -2(delta_a^e delta_b^f - delta_a^f delta_b^e) | n/a | Minkowski and curved lab | 0 | 2.220446049250e-16 | 2.220446049250e-16 | - | 1.000000000000e-13 | 1.370000000000e-05 |
| PASS | A | b1c.dual.two_form_double_dual_is_minus_identity | **F = -F for an arbitrary two-form, in four-dimensional Lorentzian signature | n/a | Minkowski and curved lab | 0 | 2.220446049250e-16 | 2.220446049250e-16 | - | 1.000000000000e-13 | 1.540000000000e-05 |
| PASS | A | b1c.dual.left_dual_squared_is_minus_identity | *(*Z) = -Z | n/a | random Z with Z's antisymmetries; Minkowski and curved lab | 0 | 3.330669073875e-16 | 3.330669073875e-16 | - | 1.000000000000e-12 | 4.140000000000e-05 |
| PASS | A | b1c.dual.right_dual_squared_is_minus_identity | (Z*)* = -Z | n/a | random Z with Z's antisymmetries; Minkowski and curved lab | 0 | 3.330669073875e-16 | 3.330669073875e-16 | - | 1.000000000000e-12 | 4.460000000000e-05 |
| PASS | A | b1c.dual.left_right_commute | *(Z*) = (*Z)*, since the two duals contract disjoint index pairs | n/a | Minkowski and curved lab | 0 | 6.661338147751e-16 | 6.661338147751e-16 | - | 1.000000000000e-12 | 3.650000000000e-05 |
| PASS | A | b1c.dual.double_dual_squared_is_identity | (*Z*) applied twice returns +Z: (-1) from each index pair | n/a | Minkowski and curved lab | 0 | 5.551115123126e-16 | 5.551115123126e-16 | - | 1.000000000000e-12 | 3.650000000000e-05 |
| PASS | B | b1c.dual.vacuum_analytic | the left dual of the vacuum Z equals (1/2) eps^{abcd} | n/a | Minkowski and curved lab | 0 | 5.551115123126e-17 | 5.551115123126e-17 | - | 1.000000000000e-13 | 1.460000000000e-05 |
| PASS | A | b1c.dual.transforms_tensorially | *(L.Z) computed in the new coordinates equals L.(*Z) computed in the old ones | n/a | linear map with det L = 1.61 | 0 | 1.110223024625e-15 | 1.110223024625e-15 | - | 1.000000000000e-12 | 1.618000000000e-04 |
| PASS | A | b1c.projector.annihilates_V_right | h^a_b V^b = 0 | n/a | 4 observers: lab_coordinate, boosted_beta0.6, eulerian_alcubierre_v0.2, eulerian_ergoregion_v2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 7.500000000000e-06 |
| PASS | A | b1c.projector.annihilates_V_left | V_a h^a_b = 0 | n/a | 4 observers: lab_coordinate, boosted_beta0.6, eulerian_alcubierre_v0.2, eulerian_ergoregion_v2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 1.390000000000e-05 |
| PASS | A | b1c.projector.idempotent | h^a_c h^c_b = h^a_b | n/a | 4 observers: lab_coordinate, boosted_beta0.6, eulerian_alcubierre_v0.2, eulerian_ergoregion_v2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 1.590000000000e-05 |
| PASS | A | b1c.projector.trace_is_three | h^a_a = 3 | n/a | 4 observers: lab_coordinate, boosted_beta0.6, eulerian_alcubierre_v0.2, eulerian_ergoregion_v2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 1.790000000000e-05 |
| PASS | A | b1c.projector.metric_reconstruction | g_ab = h_ab - V_a V_b | n/a | 4 observers: lab_coordinate, boosted_beta0.6, eulerian_alcubierre_v0.2, eulerian_ergoregion_v2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 1.990000000000e-05 |
| PASS | A | b1c.projector.index_forms_consistent | h^{ac} g_cb = h^a_b, so the upper, lower and mixed forms are the same object | n/a | 4 observers: lab_coordinate, boosted_beta0.6, eulerian_alcubierre_v0.2, eulerian_ergoregion_v2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 2.700000000000e-05 |
| PASS | D | b1c.bel.a7_contractions_independent | all four A7 contractions match a reference computed in the test with its own loops | n/a | random Z with Z's antisymmetries, 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 1.503000000000e-04 |
| PASS | B | b1c.bel.vacuum_epsilon | eps^{ab} = h^{ab} for vacuum | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 7.780000000000e-05 |
| PASS | B | b1c.bel.vacuum_mu_inverse | [mu^-1]^{ab} = h^{ab} for vacuum | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 8.080000000000e-05 |
| PASS | B | b1c.bel.vacuum_zeta_zero | zeta = 0 for vacuum | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 8.260000000000e-05 |
| PASS | B | b1c.bel.vacuum_zeta_adjoint_zero | zeta+ = 0 for vacuum | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 8.420000000000e-05 |
| PASS | B | b1c.bel.isotropic_rest_epsilon | eps^{ab} = eps h^{ab} from (B1a)/(B2), with eps = 2.75 | n/a | eps = 2.75, mu^-1 = 0.4, 4 observers | 0 | 8.881784197001e-16 | 8.881784197001e-16 | - | 1.000000000000e-12 | 7.660000000000e-05 |
| PASS | B | b1c.bel.isotropic_rest_mu_inverse | [mu^-1]^{ab} = mu^-1 h^{ab} from (B1b)/(B2), with mu^-1 = 0.4 | n/a | eps = 2.75, mu^-1 = 0.4, 4 observers | 0 | 3.330669073875e-16 | 3.330669073875e-16 | - | 1.000000000000e-12 | 7.960000000000e-05 |
| PASS | B | b1c.bel.isotropic_rest_zeta_zero | zeta = 0 for an isotropic medium AT REST relative to the observer | n/a | 4 observers | 0 | 4.440892098501e-16 | 4.440892098501e-16 | - | 1.000000000000e-12 | 8.420000000000e-05 |
| PASS | A | b1c.bel.left_right_dual_distinct | for a magnetoelectric medium *Z != Z*, and the two A7 duals therefore give genuinely different magnetoelectric blocks | n/a | n/a | both differences well above roundoff | max\|*Z - Z*\| = 0.200002, max\|W_zeta - W_zeta^T\| = 0.200002 | - | - | exact | 3.310000000000e-05 |
| PASS | A | b1c.spatiality.epsilon.first_index | V_a epsilon^{ab} = 0 | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 2.000000000000e-07 |
| PASS | A | b1c.spatiality.epsilon.second_index | epsilon^{ab} V_b = 0 | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 2.600000000000e-06 |
| PASS | A | b1c.projection_invariance.epsilon | h^a_c h^b_d epsilon^{cd} = epsilon^{ab} | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 4.600000000000e-06 |
| PASS | A | b1c.spatiality.mu_inverse.first_index | V_a mu_inverse^{ab} = 0 | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 1.000000000000e-07 |
| PASS | A | b1c.spatiality.mu_inverse.second_index | mu_inverse^{ab} V_b = 0 | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 2.200000000000e-06 |
| PASS | A | b1c.projection_invariance.mu_inverse | h^a_c h^b_d mu_inverse^{cd} = mu_inverse^{ab} | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 1.110000000000e-05 |
| PASS | A | b1c.spatiality.zeta.first_index | V_a zeta^{ab} = 0 | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 1.000000000000e-07 |
| PASS | A | b1c.spatiality.zeta.second_index | zeta^{ab} V_b = 0 | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 2.300000000000e-06 |
| PASS | A | b1c.projection_invariance.zeta | h^a_c h^b_d zeta^{cd} = zeta^{ab} | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 4.200000000000e-06 |
| PASS | A | b1c.spatiality.zeta_adjoint.first_index | V_a zeta_adjoint^{ab} = 0 | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 2.000000000000e-06 |
| PASS | A | b1c.spatiality.zeta_adjoint.second_index | zeta_adjoint^{ab} V_b = 0 | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 4.000000000000e-06 |
| PASS | A | b1c.projection_invariance.zeta_adjoint | h^a_c h^b_d zeta_adjoint^{cd} = zeta_adjoint^{ab} | n/a | 4 observers | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 5.700000000000e-06 |
| PASS | A | b1c.spatiality.zeta_is_nonzero | the magnetoelectric block is genuinely nonzero for the warp medium, so its spatiality and projection tests are not vacuous | n/a | n/a | max\|zeta\| well above roundoff | max\|zeta\| = 0.200002 | - | - | exact | 1.870000000000e-05 |
| PASS | D | b1c.bel.zeta_normalisation | zeta = +2 W_zeta, determined from the (B17) moving-medium benchmark | n/a | 5 boost configurations | 2.000000000000e+00 | 2.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | exact | 1.000000000000e-07 |
| PASS | B | b1c.b17.zeta_is_nonzero | the (B17) benchmark produces a genuinely nonzero magnetoelectric tensor in every configuration, so the normalisation is actually being constrained | n/a | n/a | min over configurations of max\|zeta\| above 1e-3 | min max\|zeta\| = 0.357407 | - | - | exact | 1.000000000000e-07 |
| PASS | D | b1c.b17.zeta_proportional_to_A7d | the (A7d) Bel block and the (B17) reference are proportional componentwise, with a single ratio per configuration | n/a | 5 boost configurations | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 2.000000000000e-07 |
| PASS | D | b1c.b17.normalisation_is_universal | the ratio zeta(B17) / W_zeta(A7d) is one constant, independent of boost magnitude, boost direction and the eps/mu pair | n/a | 5 configurations; ratios: 2.0000000000 2.0000000000 2.0000000000 2.0000000000 2.0000000000  | 0 | 6.661338147751e-16 | 6.661338147751e-16 | - | 1.000000000000e-11 | 6.500000000000e-06 |
| PASS | D | b1c.b17.normalisation_determined | the magnetoelectric normalisation constant c_zeta is determined from (B17), with one specific sign | n/a | n/a | a single determined constant | zeta = (2.000000000000) * W_zeta | - | - | exact | 2.100000000000e-06 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
