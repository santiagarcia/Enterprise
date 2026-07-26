# program_b3 verification report

**B3-PHYSICAL-ADMISSIBILITY: frequency-dependent constitutive models. q_physical -> eps(w), mu^-1(w), zeta(w) -> Z(w) -> J_band, with structural causality and a derived passivity condition**

| | |
|---|---|
| Generated (UTC) | 2026-07-26T18:04:38Z |
| Git commit | `976ef4960578d7eece4aa22a861cbeb1c63c9144` |
| Compiler | MSVC 19.44.35223.0 |
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
| PASS | A | b3.complex.arithmetic | addition, multiplication and division against hand-computed values; (a/b)*b = a; a conj(a) is real and equals \|a\|^2 | n/a | a = 2-3i, b = -1.5+0.75i | 0 | 2.220446049250e-16 | 2.220446049250e-16 | - | 1.000000000000e-15 | 1.140000000000e-05 |
| PASS | A | b3.complex.transpose_and_adjoint_are_distinct | for a complex symmetric matrix the plain transpose residual vanishes while the Hermitian adjoint residual does not | n/a | n/a | \|A - A^T\| = 0 and \|A - A^dagger\| > 0 | transpose residual 0.0000e+00, adjoint residual 1.0000e+00 | - | - | exact | 1.490000000000e-05 |
| PASS | A | b3.complex.matrix_inverse_roundtrip | A A^-1 = I for a complex 3x3 | n/a | adjugate form, branch-free | 0 | 2.220446049250e-16 | 2.220446049250e-16 | - | 1.000000000000e-14 | 2.880000000000e-05 |
| PASS | B | b3.complex.real_limit_recovery | with every damping set to zero the response is exactly real and matches the undamped closed form | n/a | 21 frequencies, gamma = 0 | 0 | 4.440892098501e-16 | 4.440892098501e-16 | - | 1.000000000000e-14 | 2.850000000000e-05 |
| PASS | B | b3.convention.sign_of_loss_is_pinned | under the recorded e^{-i omega t} kernel both specified model forms have STRICTLY POSITIVE imaginary part for positive parameters, which is passive; under the opposite kernel both would be gain | n/a | n/a | Im > 0 for the Lorentz term across the band and for the Drude term | kernel exp(-i*omega*t); min Im(Lorentz) over band = 5.6073e-03; Im(Drude) at w=1 = 1.4257e-01 | - | - | exact | 1.490000000000e-05 |
| PASS | A | b3.convention.recorded | the Fourier kernel, time derivative, constitutive relation, dissipated power, Poynting vector, passive sign, causal half-plane and reality condition are all recorded | n/a | n/a | all conventions recorded | exp(-i*omega*t) \| <p> = (omega/2) Im[ E* . D + H* . B ] = (omega/2) v^dagger L v, L = (C - C^dagger)/(2i) | - | - | exact | 2.090000000000e-05 |
| PASS | B | b3.dispersion.lorentz_closed_form | the two-pole Lorentz permittivity matches a hand-written closed form | n/a | 3 frequencies | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 6.200000000000e-06 |
| PASS | A | b3.dispersion.mu_inverse_is_derived_from_mu | [mu^-1](w) mu(w) = 1 at every sampled frequency | n/a | w from 0.5 to 5 | 0 | 2.220446049250e-16 | 2.220446049250e-16 | - | 1.000000000000e-13 | 2.419000000000e-04 |
| PASS | A | b3.dispersion.reality_condition | eps(-w) = conj(eps(w)), the condition for a real time-domain response | n/a | 4 frequencies | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 2.488000000000e-04 |
| PASS | A | b3.dispersion.magnetoelectric_is_reciprocal_by_construction | zeta_adjoint = zeta^T at every frequency, exactly | n/a | 21 frequencies; max\|zeta\| = 7.5000e-02 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 3.100000000000e-05 |
| PASS | B | b3.causality.poles_in_lower_half_plane | every pole of every Lorentz term sits in the LOWER half omega-plane, including the overdamped case where both roots are purely imaginary | n/a | n/a | Im(pole) < 0 for all 8 roots | eps pole 0: 2.9986-0.0900i, -2.9986-0.0900i; eps pole 1: 4.3974-0.1500i, -4.3974-0.1500i; mu pole 0: 3.5983-0.1100i, -3.5983-0.1100i; overdamped: 0.0000-0.0134i, 0.0000-2.9866i (overdamped);  | - | - | exact | 1.840000000000e-05 |
| PASS | A | b3.causality.positive_transform_is_strictly_positive | the damping transform is strictly positive for every input, including extreme negative values, so no optimizer step can produce an acausal model | n/a | n/a | positive for all inputs | smallest value over u in {-1e6 .. 1e2} = 1.0000e-12 | - | - | exact | 2.370000000000e-05 |
| PASS | A | b3.causality.transform_roundtrip | positive_transform(positive_transform_inverse(x)) = x, so a known physical model can be expressed exactly as a design vector | n/a | 5 values | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 2.730000000000e-05 |
| PASS | D | b3.causality.kramers_kronig_diagnostic | the Kramers-Kronig residual DECREASES monotonically as the integration band widens, and is small on the widest band -- the signature of quadrature truncation rather than a causality violation | n/a | n/a | monotone decrease with band width, and below 5e-3 on the widest | [1.000e-02,1.000e+02] -> 6.089e-08  [1.000e-03,1.000e+03] -> 6.075e-11  [1.000e-04,1.000e+04] -> 6.084e-14   | - | - | exact | 2.426929000000e-01 |
| PASS | A | b3.passivity.loss_matrix_is_hermitian | L = (C - C^dagger)/(2i) is Hermitian at every band frequency | n/a | 21 frequencies | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 1.436000000000e-04 |
| PASS | B | b3.passivity.reciprocity_in_EH_representation | xi_EH = zeta_EH^T, the reciprocity relation DERIVED for this project's (E,B) sign convention | n/a | 21 frequencies | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 1.492000000000e-04 |
| PASS | B | b3.passivity.vacuum_is_lossless | the nondispersive vacuum limit has loss eigenvalues exactly zero: passive, and on the boundary rather than inside it | n/a | n/a | min eigenvalue of L is 0 | min loss eigenvalue 0.0000e+00 | - | - | exact | 9.060000000000e-05 |
| PASS | B | b3.passivity.damped_lorentz_is_passive | a damped two-pole Lorentz medium is passive across the whole band: the smallest eigenvalue of the loss matrix is non-negative | n/a | n/a | min eigenvalue of L >= 0 at every frequency | min loss eigenvalue 2.4604e-03 at w = 2.0000e+00 | - | - | exact | 1.041000000000e-04 |
| PASS | B | b3.passivity.active_oscillator_is_rejected | a Lorentz term with NEGATIVE damping is detected as active -- the loss matrix acquires a negative eigenvalue -- and its poles move into the upper half plane, so it is acausal as well | n/a | n/a | passivity violated and poles in the upper half plane | min loss eigenvalue -2.9630e+00, violation 2.9630e+00, pole Im = 9.0000e-02 | - | - | exact | 1.005000000000e-04 |
| PASS | B | b3.passivity.diagonal_imaginary_parts_are_not_sufficient | a medium with STRICTLY POSITIVE Im eps and Im mu at every frequency is nevertheless NOT passive when the magnetoelectric coupling is strong: the full loss matrix has a negative eigenvalue | n/a | n/a | positive diagonal imaginary parts, negative loss eigenvalue | min Im eps = 3.1250e-06, min Im mu = 1.1211e-06, yet min eigenvalue of L = -3.7943e+00 at w = 3.6000e+00 | - | - | exact | 1.463000000000e-04 |
| PASS | A | b3.band.a11_is_real_linear | (A11) is homogeneous linear in the blocks with real coefficients | n/a | alpha = 1.7, beta = -0.9 | 0 | 4.440892098501e-16 | 4.440892098501e-16 | - | 1.000000000000e-12 | 2.950000000000e-05 |
| PASS | A | b3.band.zero_at_exact_target_and_reported_separately | J_band vanishes at the exact target, and the centre, mean and worst frequency mismatches are reported as three distinct numbers | n/a | n/a | all mismatch measures zero at the target | centre 0.0000e+00, mean 0.0000e+00, worst 0.0000e+00 at w = 0.0000e+00 | - | - | exact | 1.111300000000e-03 |
| PASS | B | b3.gradient.analytic_one_lorentz_pole | d eps/dF = 1/D, d eps/dw0 = -2 F w0/D^2, d eps/dgamma = i F w/D^2 match centred differences | n/a | 3 step sizes, w = 1.3 | 0 | 7.099333898530e-10 | 7.099333898530e-10 | - | 1.000000000000e-06 | 7.000000000000e-07 |
| PASS | D | b3.gradient.complex_step | complex-step derivatives, from instantiating the model at Scalar = Complex<double> with a step of 1e-30, agree with centred differences | n/a | 8 parameters, h = 1e-30 against centred h = 1e-6 | 0 | 1.790479574500e-10 | 1.790479574500e-10 | - | 1.000000000000e-06 | 1.410000000000e-05 |
| PASS | A | b3.gradient.oti_is_batched | the band objective's gradient is assembled from batches of at most kOtiDirections = 3 active directions | n/a | n/a | 4 batches | 4 batches over 11 parameters, active [0,1,2][3,4,5][6,7,8][9,10] | - | - | exact | 9.186440000000e-02 |
| PASS | D | b3.gradient.oti_matches_finite_difference | the batch-assembled HyPAD gradient of J_band matches centred differences | n/a | 11 parameters | 0 | 1.452185593998e-10 | 1.452185593998e-10 | - | 1.000000000000e-06 | 9.187100000000e-02 |
| PASS | B | b3.recovery.one_pole | a manufactured single-pole target is recovered from a 20 percent perturbation, in physical parameters and not merely in objective | n/a | n/a | J_band below 1e-16 and every physical parameter within 1e-4 | J: 2.4407e-02 -> 9.9441e-19, max\|dp\| = 6.7958e-09, ObjectiveBelowTolerance | - | - | exact | 5.863817800000e+00 |
| PASS | B | b3.recovery.two_pole | a two-pole target is recovered to a band mismatch below 1e-12 | n/a | n/a | J_band below 1e-14 | J: 6.7191e-02 -> 9.9876e-19; centre 1.6712e-18, mean 9.4358e-19, worst 2.6670e-18 at w = 4.6000e+00, ObjectiveBelowTolerance | - | - | exact | 9.155838700000e+00 |
| PASS | B | b3.restricted_one_pole_approximation | a one-pole model fitted to a two-pole target leaves a STRICTLY NONZERO residual, and still improves substantially on its start | n/a | n/a | nonzero residual, and a genuine improvement | J: 1.3813e-02 -> 1.3446e-02, StepBelowTolerance | - | - | exact | 6.979785000000e-01 |
| PASS | B | b3.hard_feasible_versus_soft_penalty | the hard-feasible run satisfies the bound EXACTLY -- zero violation -- and is labelled feasible, while the soft-penalised run is reported under a different name whatever the size of its violation | n/a | n/a | hard run exactly feasible; soft run never labelled feasible | hard feasible approximation: violation 0.0000e+00, J = 4.7906e-02; soft penalized approximation: violation 2.2401e-06, J = 9.9888e-03 | - | - | exact | 2.860108600000e+01 |
| PASS | B | b3.bandwidth_tradeoff | a fixed one-pole model fits a two-pole target progressively worse as the band widens, and the degradation is monotone | n/a | n/a | monotone increase in residual with bandwidth | J at half-width 0.15 = 3.3696e-06, at 1.6 = 1.2561e-02 | - | - | exact | 2.095336420000e+01 |
| PASS | B | b3.recovery.reciprocal_magnetoelectric | a dispersive magnetoelectric target is recovered while reciprocity holds exactly at every frequency throughout | n/a | n/a | recovered, with zeta_adjoint = zeta^T exactly | J -> 9.9668e-19, reciprocity residual 0.0000e+00, min loss eigenvalue 2.4171e-03 | - | - | exact | 8.202126100000e+00 |
| PASS | A | b3.frame_orientation_matters_for_anisotropic_shapes | the same coefficients on an independently rebuilt triad give a DIFFERENT tensor, because the magnetoelectric shape K singles out particular spatial directions | n/a | n/a | a measurable difference between pushed-forward and rebuilt frames | max\|Z(pushed) - Z(rebuilt)\| = 3.9715e-03 | - | - | exact | 2.551900000000e-03 |
| PASS | A | b3.coordinate_change_consistency_per_frequency | the same physical model built in transformed coordinates gives the pushed-forward Z at EVERY band frequency, real and imaginary parts alike | n/a | 21 frequencies, non-orthogonal L | 0 | 1.332267629550e-15 | 1.332267629550e-15 | - | 1.000000000000e-11 | 2.556300000000e-03 |
| PASS | B | b3.real_limit_reduces_to_b2 | a nondispersive lossless model reproduces the B2 isotropic design's Z exactly at every frequency, with identically zero imaginary part | n/a | n/a | B3 real limit equals the B2 tensor, imaginary part exactly zero | max\|Z_B3 - Z_B2\| = 0.0000e+00, max\|Im Z\| = 0.0000e+00 | - | - | exact | 2.716000000000e-04 |
| PASS | B | b3.limit.vacuum_nondispersive | eps = mu = 1 with no poles reproduces the vacuum Z at every frequency | n/a | 21 frequencies | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 2.417000000000e-04 |
| PASS | B | b3.limit.lossless_lorentz_away_from_resonance | an undamped Lorentz medium evaluated well below its resonance is exactly lossless and sits on the passivity boundary | n/a | n/a | Im eps exactly zero and min loss eigenvalue zero | max\|Im eps\| = 0.0000e+00, min loss eigenvalue 0.0000e+00 | - | - | exact | 5.690000000000e-05 |
| PASS | A | b3.deterministic_repeatability | two band-optimizer runs from the same initialization are bit-identical | n/a | n/a | bit-identical runs | identical over 401 iterations | - | - | exact | 1.035569670000e+01 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
