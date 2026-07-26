# program_b2 verification report

**B2-ABSTRACT-INVERSE: abstract constitutive inverse design. g_eff^target -> Z_target -> observer blocks -> nearest admissible constitutive design, over the A11 reconstruction verified at B1C**

| | |
|---|---|
| Generated (UTC) | 2026-07-26T22:17:42Z |
| Git commit | `9bf1028688994b1fd2c9718f519998ef2ea1d7dc` |
| Compiler | MSVC 19.44.35223.0 |
| Build type | RelWithDebInfo |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 33 |
| Passed | 33 |
| Failed | 0 |
| Skipped | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | A | b2.observer_frame_is_orthonormal_and_spatial | the spatial triad satisfies g_lab(e_i, e_j) = delta_ij and g_lab(e_i, V) = 0 | n/a | n/a | frame built | orthonormal spatial triad adapted to V. | - | - | exact | 3.900000000000e-06 |
| PASS | A | b2.observer_frame_orthonormality | g_lab(e_i, e_j) = delta_ij and g_lab(e_i, V) = 0 | n/a | g_lab = eta, V = (1,0,0,0) | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 2.170000000000e-05 |
| PASS | A | b2.parameter_vector_roundtrip | q -> coefficients -> q is the identity, and q -> blocks -> coefficients likewise, for every linear family | n/a | 4 families | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 2.440000000000e-05 |
| PASS | A | b2.symmetry_by_construction | eps and mu^-1 are symmetric and zeta_adjoint = zeta^T for every family, with no symmetrisation step applied afterwards | n/a | 5 families | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 6.690000000000e-05 |
| PASS | A | b2.spatiality_by_construction | every block is spatial with respect to V, structurally | n/a | 5 families | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 7.460000000000e-05 |
| PASS | A | b2.reconstructed_tensor_is_reciprocal | Z_design^{abcd} = Z_design^{cdab} for every family | n/a | 5 families | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 7.950000000000e-05 |
| PASS | A | b2.a11_reconstruction_roundtrip | design blocks -> Z via (A11) -> blocks via (A7) returns the design blocks | n/a | 4 perturbed magnetoelectric designs | 0 | 2.220446049250e-16 | 2.220446049250e-16 | - | 1.000000000000e-12 | 1.532000000000e-04 |
| PASS | A | b2.norms_agree_in_the_analogue_configuration | with g_lab = eta and V = (1,0,0,0) the observer Gram is the identity, so the coordinate-basis and observer-referred norms coincide exactly | n/a | magnetoelectric design | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 8.950000000000e-05 |
| PASS | A | b2.observer_referred_norm_is_invariant | \|\|Z\|\|^2_V is unchanged by a nontrivial linear coordinate change | n/a | non-orthogonal L, det = 1.2170e+00 | 0 | 1.612495088506e-16 | 1.612495088506e-16 | - | 1.000000000000e-12 | 9.360000000000e-05 |
| PASS | A | b2.coordinate_basis_norm_is_not_invariant | the coordinate-basis norm DOES change under the same coordinate change, which is why it is labelled a coordinate-basis norm | n/a | n/a | a relative change above 1e-2 | coord: 2.2032e+01 -> 5.0390e+01, observer-referred: 2.2032e+01 -> 2.2032e+01 | - | - | exact | 1.060000000000e-04 |
| PASS | A | b2.objective_zero_at_exact_target | J_Z and J_blocks vanish when the design equals the target | n/a | 3 designs | 0 | 7.893660995243e-32 | 7.893660995243e-32 | - | 1.000000000000e-28 | 1.977000000000e-04 |
| PASS | A | b2.objective_positive_away_from_target | J_Z is strictly positive for every perturbed design | n/a | n/a | smallest J_Z over 9 perturbations above 1e-6 | smallest J_Z = 2.4916e-03 | - | - | exact | 2.011000000000e-04 |
| PASS | A | b2.objective_is_exactly_quadratic | J_Z is an exact quadratic in q for the linear families: its third difference along any direction vanishes | n/a | 4 base points, h = 0.25 | 0 | 2.303712776097e-15 | 2.303712776097e-15 | - | 1.000000000000e-12 | 1.704000000000e-04 |
| PASS | B | b2.gradient_manufactured_analytic | for J = sum c_k (q_k - a_k)^2 the centred difference reproduces the hand-derived 2 c_k (q_k - a_k) | n/a | 21 parameters, 4 step sizes | 0 | 3.376410262490e-12 | 3.376410262490e-12 | - | 1.000000000000e-09 | 4.900000000000e-06 |
| PASS | B | b2.oti_gradient_manufactured_matches_analytic | the batch-assembled HyPAD gradient of J = sum c_k (q_k - a_k)^2 equals the hand-derived 2 c_k (q_k - a_k) exactly | n/a | 21 parameters in 7 batches of at most 3 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 2.630000000000e-05 |
| PASS | D | b2.gradient_finite_difference_matches_analytic | centred differences of J_Z reproduce 2 A_k . (A q - z) / s, the gradient implied by the linear design map | n/a | 21 parameters, h in {1e-3, 1e-4, 1e-5} | 0 | 1.285398870676e-12 | 1.285398870676e-12 | - | 1.000000000000e-08 | 1.255500000000e-03 |
| PASS | B | b2.oti_gradient_matches_analytic | the batch-assembled HyPAD gradient of J_Z equals the closed-form 2 A_k . (A q - z) / s implied by the linear design map | n/a | 21 parameters | 0 | 3.469446951954e-18 | 3.469446951954e-18 | - | 1.000000000000e-13 | 2.227600000000e-03 |
| PASS | A | b2.oti_gradient_is_batched | the 21-parameter gradient is assembled from batches of at most kOtiDirections = 3 active directions | n/a | n/a | 7 batches of at most 3 | 7 evaluations, active indices [0,1,2][3,4,5][6,7,8][9,10,11][12,13,14][15,16,17][18,19,20] | - | - | exact | 2.248100000000e-03 |
| PASS | D | b2.oti_gradient_matches_finite_difference | the batch-assembled HyPAD gradient equals the centred-difference gradient of the same objective | n/a | 21 parameters | 0 | 1.285398870676e-12 | 1.285398870676e-12 | - | 1.000000000000e-08 | 2.250800000000e-03 |
| PASS | A | b2.oti_gradient_batch_independence | a partial derivative is identical whether its parameter was seeded alone or alongside two others | n/a | 5 parameters, batch of 1 versus 3 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 3.376700000000e-03 |
| PASS | B | b2.recovery.isotropic | starting from a perturbed design, the optimizer recovers the known target: relative tensor residual below 1e-6 and every parameter within 1e-4 | n/a | n/a | converged onto the known answer | J: 2.8890e-02 -> 9.4689e-21 in 71 iterations, ObjectiveBelowTolerance; max\|dq\| = 2.6017e-10 | - | - | exact | 3.492500000000e-03 |
| PASS | A | b2.optimizer_descends.isotropic | the objective decreases monotonically and the line search succeeds at every accepted iteration | n/a | n/a | monotone objective history | 72 iterations recorded, 432 objective evaluations | - | - | exact | 3.498000000000e-03 |
| PASS | B | b2.recovery.anisotropic | starting from a perturbed design, the optimizer recovers the known target: relative tensor residual below 1e-6 and every parameter within 1e-4 | n/a | n/a | converged onto the known answer | J: 5.1506e-02 -> 8.4309e-21 in 200 iterations, ObjectiveBelowTolerance; max\|dq\| = 2.1562e-10 | - | - | exact | 7.727180000000e-02 |
| PASS | A | b2.optimizer_descends.anisotropic | the objective decreases monotonically and the line search succeeds at every accepted iteration | n/a | n/a | monotone objective history | 201 iterations recorded, 8844 objective evaluations | - | - | exact | 7.727710000000e-02 |
| PASS | B | b2.recovery.magnetoelectric | starting from a perturbed design, the optimizer recovers the known target: relative tensor residual below 1e-6 and every parameter within 1e-4 | n/a | n/a | converged onto the known answer | J: 7.6055e-02 -> 8.8153e-21 in 217 iterations, ObjectiveBelowTolerance; max\|dq\| = 2.8777e-10 | - | - | exact | 8.483400000000e-02 |
| PASS | A | b2.optimizer_descends.magnetoelectric | the objective decreases monotonically and the line search succeeds at every accepted iteration | n/a | n/a | monotone objective history | 218 iterations recorded, 9592 objective evaluations | - | - | exact | 8.483920000000e-02 |
| PASS | B | b2.recovery.metric_induced | a subluminal Alcubierre effective metric with g_lab = eta is recovered in the unrestricted reciprocal family, compared against Z_target itself | n/a | n/a | converged onto Z_target | target: reciprocal; pair-antisymmetric; finite; observer-decomposable; metric-induced. J: 2.7186e-02 -> 6.8098e-21, relative \|Z_design - Z_target\| = 1.4297e-10, ObjectiveBelowTolerance | - | - | exact | 1.715850000000e-02 |
| PASS | D | b2.restricted_family_is_a_projection_not_a_fit | against a magnetoelectric target, the zeta == 0 family leaves a STRICTLY NONZERO residual, and the optimizer finds the same point as the exact least-squares projection | n/a | n/a | nonzero residual, matching the normal-equation solution | optimizer J = 7.5053e-02, normal-equation J = 7.5053e-02, max\|dq\| = 7.1098e-11 | - | - | exact | 4.874600000000e-02 |
| PASS | B | b2.infeasible_bounded_reports_best_approximation | with bounds excluding the target, the result is a BEST BOUNDED APPROXIMATION: J_Z stays strictly positive, the bound is reported as active, and no exact realization is claimed | n/a | n/a | nonzero J_Z with the entry bound respected and reported | J_Z = 1.7450e-01, penalty = 2.2031e-05, max\|eps\| = 1.5001e+00 against bound 1.5, entry bound active: yes, MaxIterationsReached | - | - | exact | 1.413461500000e+00 |
| PASS | B | b2.coordinate_change_consistency | optimizing entirely in the transformed coordinates reconstructs the transformed target: the FINAL TENSOR RELATIONSHIP holds | n/a | n/a | Z'_design = Z'_target | primed target reciprocal; pair-antisymmetric; finite; observer-decomposable; not metric-induced. J -> 9.5999e-21, relative \|Z'_design - Z'_target\| = 1.0971e-10, ObjectiveBelowTolerance | - | - | exact | 2.484516000000e-01 |
| PASS | A | b2.deterministic_repeatability | two runs from the same initialization are BIT-IDENTICAL in history, termination reason and final parameters | n/a | n/a | bit-identical runs | identical over 200 iterations | - | - | exact | 1.382641000000e-01 |
| PASS | A | b2.nonfinite_is_rejected_not_absorbed | a non-finite design is rejected with reason NonFiniteEncountered, and an overflowing design is reported as not finite rather than silently scored | n/a | n/a | NonFiniteEncountered, and evaluation reports finite = false | NonFiniteEncountered; overflow evaluation finite = false | - | - | exact | 4.940000000000e-05 |
| PASS | A | b2.target_classification_is_recorded | a metric-induced target is classified as metric-induced, reciprocal, pair-antisymmetric, finite and observer-decomposable | n/a | n/a | all six classification flags set | reciprocal; pair-antisymmetric; finite; observer-decomposable; metric-induced. | - | - | exact | 3.560000000000e-05 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
