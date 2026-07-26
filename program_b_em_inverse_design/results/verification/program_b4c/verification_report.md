# program_b4c verification report

**B4C-PROPAGATION-CORRECTION: the metricity claim withdrawn and replaced, propagation-invisible constitutive pieces demonstrated, the premetric classification reported piece by piece, the layered solver made scalar-generic and differentiated through, and a root finder that sees repeated and complex roots**

| | |
|---|---|
| Generated (UTC) | 2026-07-26T22:19:10Z |
| Git commit | `9bf1028688994b1fd2c9718f519998ef2ea1d7dc` |
| Compiler | MSVC 19.44.35223.0 |
| Build type | RelWithDebInfo |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 18 |
| Passed | 18 |
| Failed | 0 |
| Skipped | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | A | b4c.metricity_claim_is_corrected | the unrestricted claim 'metricity iff no birefringence' is WITHDRAWN and replaced by the one-directional statement | n/a | n/a | the corrected wording states only the provable direction | Metric-induced constitutive tensors are non-birefringent. Generic reciprocal constitutive tensors need not be. The converse requires additional assumptions and is subject to conformal/dilaton freedom and constitutive pieces that do not affect the Fresnel polynomial. | - | - | exact | 3.400000000000e-06 |
| PASS | B | b4c.axion_changes_the_tensor_but_not_the_cone | adding alpha eps^{abcd} to a metric-induced tensor changes Z and changes the zeta block by exactly 2 alpha h, while leaving the Fresnel polynomial bit-identical | n/a | n/a | tensor and zeta differ; eps, mu^-1 and the cone do not | max\|Z_sum - Z_metric\| = 3.7000e-01; block differences eps 0.0000e+00, mu^-1 0.0000e+00, zeta 7.4000e-01 (predicted 2*alpha*h = 7.4000e-01); Fresnel coefficient difference 0.0000e+00; alpha recovered 3.7000e-01 against 0.37 | - | - | exact | 1.514000000000e-04 |
| PASS | B | b4c.conformal_freedom_is_explicit | two effective metrics differing by a conformal factor give the SAME Z and therefore the same cone: the metric is determined only up to that factor | n/a | n/a | metrics differ, Z does not | max\|g2 - g1\| = 2.3000e+00, max\|Z2 - Z1\| = 1.1102e-16 | - | - | exact | 4.200000000000e-06 |
| PASS | A | b4c.classification_reports_each_piece_separately | the metric tensor, the metric-plus-axion tensor and a birefringent anisotropic tensor are distinguished by SEPARATE flags: axion content, skewon content, single-cone, and principal-part metric-compatibility | n/a | n/a | three tensors, three distinct flag patterns | metric: reciprocal; no skewon; no axion; single Fresnel cone (principal part metric-COMPATIBLE, which does NOT establish metricity: the axion piece is invisible to the cone and the cone fixes only a conformal class) \| with axion: reciprocal; no skewon; contains axion; single Fresnel cone (principal part metric-COMPATIBLE, which does NOT establish metricity: the axion piece is invisible to the cone and the cone fixes only a conformal class) \| anisotropic: reciprocal; no skewon; no axion; birefringent: two distinct Fresnel cones, so the response is NOT metric-induced | - | - | exact | 5.416000000000e-04 |
| PASS | B | b4c.roots.distinct_repeated_and_near_split | four distinct real roots, two exact DOUBLE roots, and a pair split by 1e-7 are all recovered, with the double roots counted at multiplicity two | n/a | n/a | all roots found with the right real count | four_distinct: 4 real, 4 distinct; two_double: 4 real, 2 distinct; near_split: 4 real, 3 distinct; worst root error 9.4414e-08, worst residual 1.2212e-15 | - | - | exact | 1.580200000000e-03 |
| PASS | B | b4c.roots.complex_conjugate_pair | two real roots and one complex-conjugate pair are separated and reported, where the sign-change scan would have seen only the two real ones | n/a | n/a | 2 real roots and 1 complex pair | 2 real, 1 complex pair(s), worst residual 1.0531e-17 | - | - | exact | 1.597900000000e-03 |
| PASS | A | b4c.roots.scaling_invariance | scaling every coefficient by 1e6 leaves the roots unchanged | n/a | factor 1e6 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 1.611800000000e-03 |
| PASS | D | b4c.roots.coefficient_reconstruction | rebuilding the coefficients from the returned root set recovers the original polynomial | n/a | quartic with a complex pair | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 1.634400000000e-03 |
| PASS | B | b4c.roots.metric_double_root_is_recovered | the root finder itself returns the metric-induced Fresnel quartic's REPEATED branches, at multiplicity two, rather than relying on a separate factorization comparison | n/a | n/a | two distinct roots, multiplicity two, coincident positive branches | 2 distinct, max multiplicity 2, positive-branch separation 0.0000e+00, worst residual 6.6174e-24 (repeated root detected: attainable accuracy is about eps^(1/m), so about 1e-8 for a double root) | - | - | exact | 1.740000000000e-05 |
| PASS | B | b4c.roots.anisotropic_split_is_recovered | the anisotropic medium's two distinct propagation branches are resolved as four distinct roots | n/a | n/a | four distinct roots and a measurable split | 4 distinct, relative split 2.4309e-01 | - | - | exact | 8.110000000000e-05 |
| PASS | B | b4c.solver.matched_slab_is_analytic | the templated solver at Scalar = double gives r = 0 and phase = omega eps d for an impedance-matched slab | n/a | n/a | zero reflection and the analytic phase | R = 0.0000e+00, phase = -5.3719e-01 against -5.3719e-01 | - | - | exact | 4.990000000000e-05 |
| PASS | B | b4c.gradient.matched_slab_matches_closed_form | d(phase)/d(eps) = omega d and d(phase)/d(thickness) = omega eps for a matched slab, obtained as Im(dt/t) through the actual solver | n/a | 3 step sizes | 0 | 6.436020250078e-08 | 6.436020250078e-08 | - | 1.000000000000e-06 | 1.777900000000e-03 |
| PASS | D | b4c.gradient.finite_difference_step_independence | phase sensitivities to permeability, magnetoelectric amplitude, and a Lorentz oscillator's strength, resonance and damping agree across three step sizes | n/a | 5 parameters, h in {1e-4, 1e-5, 1e-6}: matched slab: d(phase)/d(eps) target 2.2100e+00, d(phase)/d(d) target 3.3800e+00; mu_inf d(phase) = 1.5434e+00; zeta_amp d(phase) = -3.3458e-01; F d(phase) = 1.2234e-01; w0 d(phase) = -1.5982e-01; gamma d(phase) = -5.7605e-03;  | 0 | 1.538025706793e-09 | 1.538025706793e-09 | - | 1.000000000000e-05 | 1.786200000000e-03 |
| PASS | B | b4c.hypad.through_actual_stack_assembly | a physical parameter seeded through interfaces, propagation and the Redheffer cascade gives d(phase)/d(eps) = omega d and d(phase)/d(thickness) = omega eps | n/a | n/a | both closed forms recovered | d(phase)/d(eps) = 2.2100e+00 against 2.2100e+00; d(phase)/d(d) = 3.3800e+00 against 3.3800e+00; d(phase) = Im(dt/t), no inverse trigonometry involved | - | - | exact | 5.936000000000e-04 |
| PASS | A | b4c.hypad.separate_versus_batched_seeding | the same partial derivative of the transmission amplitude comes out identically whether eps was seeded alone or alongside the thickness | n/a | batch of 1 versus 2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-15 | 1.190700000000e-03 |
| PASS | B | b4c.sensitivity.magnetoelectric_scale_matters | the propagation phase is strongly sensitive to the magnetoelectric scale in the exact profile, so the 1.6-radian mismatch of B4 is not an accident of one parameter choice | n/a | n/a | a large phase sensitivity to the magnetoelectric scale | d(phase)/d(zeta scale): exact 1.3560e+00, restricted 1.8736e+00 | - | - | exact | 3.322390000000e-02 |
| PASS | B | b4c.no_regression_in_energy_balance | after templating the solver, R + T + A = 1 still holds with A integrated independently | n/a | the same lossy slab B4 used, 3 frequencies | 0 | 6.866251789361e-09 | 6.866251789361e-09 | - | 1.000000000000e-07 | 1.660300000000e-03 |
| PASS | A | b4c.deterministic_repeatability | two identical solves of the templated stack, and two identical root solves, are bit-identical | n/a | n/a | bit-identical | identical | - | - | exact | 1.751700000000e-03 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
