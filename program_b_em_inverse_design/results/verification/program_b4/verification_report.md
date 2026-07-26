# program_b4 verification report

**B4-MAXWELL-PROPAGATION: does the B3 constitutive model produce the electromagnetic propagation the target effective metric predicts? Plane-wave operator, Fresnel polynomial, layered propagation and energy balance**

| | |
|---|---|
| Generated (UTC) | 2026-07-26T22:19:09Z |
| Git commit | `9bf1028688994b1fd2c9718f519998ef2ea1d7dc` |
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
| PASS | A | b4.convention.maxwell_signs_recorded | the frequency-domain curl equations are recorded and derived from this project's own e^{-i omega t} kernel | n/a | n/a | curl E = +i w B and curl H = -i w D + J | curl E = +i*omega*B ; curl H = -i*omega*D + J | - | - | exact | 3.400000000000e-06 |
| PASS | B | b4.vacuum.plane_wave_dispersion | \|k\| = omega in vacuum with c = 1: the analytic transverse mode is annihilated by the plane-wave operator exactly on the light cone | n/a | 3 directions, omega = 1.3 | 0 | 1.110223024625e-16 | 1.110223024625e-16 | - | 1.000000000000e-15 | 4.250000000000e-05 |
| PASS | A | b4.vacuum.off_shell_is_not_singular | at \|k\| = 1.7 omega the operator is NOT singular, so the previous case is not passing for a trivial reason | n/a | n/a | smallest singular value well away from zero off the cone | min off-shell singular value = 9.1000e-01 | - | - | exact | 1.064000000000e-04 |
| PASS | B | b4.isotropic.analytic_dispersion | \|k\|^2 = eps mu omega^2 for an isotropic medium, with H_z = sqrt(eps/mu) E_y | n/a | eps in {2,4}, mu in {1,1.5} | 0 | 8.881784197001e-16 | 8.881784197001e-16 | - | 1.000000000000e-15 | 6.500000000000e-06 |
| PASS | B | b4.isotropic.phase_velocity | v_p = 1/sqrt(eps mu) | n/a | 4 media | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-15 | 1.200000000000e-05 |
| PASS | B | b4.vacuum.transversality_and_curl | k.E = k.H = 0 and k x E = omega H for the vacuum transverse mode | n/a | n/a | both divergences and the Faraday relation vanish | max\|k.E\|,\|k.H\| = 0.0000e+00, \|k x E - w H\| = 0.0000e+00 | - | - | exact | 2.000000000000e-07 |
| PASS | B | b4.fresnel.vacuum_factorizes_into_the_light_cone | the Fresnel quartic from adj(W) = lambda k (x) k is proportional to (eta^{ab} k_a k_b)^2, with ONE constant across every wavevector | n/a | 6 wavevectors; rank-one residual 3.8722e-16, constant -1.2500e-01 | 0 | 2.627746803847e-16 | 2.627746803847e-16 | - | 1.000000000000e-12 | 8.100000000000e-06 |
| PASS | A | b4.fresnel.conformal_invariance_of_the_cone | the Fresnel scalar is unchanged under g_eff -> Omega^2 g_eff, so the light cone is a property of the conformal class | n/a | vacuum and Alcubierre v_s = 0.3, Omega^2 = 2.7 | 0 | 1.207610826990e-15 | 1.207610826990e-15 | - | 1.000000000000e-12 | 7.800000000000e-06 |
| PASS | B | b4.fresnel.metric_media_do_not_birefringe | a metric-induced medium has its two propagation branches degenerate, while a generic reciprocal anisotropic medium splits them | n/a | n/a | degenerate for the metric medium, split for the generic one | metric split = 0.0000e+00, generic split = 2.4309e-01 | - | - | exact | 4.460000000000e-05 |
| PASS | D | b4.group_velocity_matches_finite_differences | v_g = dw/dk from the dispersion curve equals c/(n + w dn/dw) | n/a | 3 frequencies, lossless Lorentz | 0 | 8.451688822033e-12 | 8.451688822033e-12 | - | 1.000000000000e-07 | 7.400000000000e-06 |
| PASS | B | b4.poynting_vacuum | <S> = (1/2) Re[E x conj(H)] gives energy flow along +x with magnitude 1/2 for a unit vacuum transverse mode | n/a | n/a | S = (0.5, 0, 0) | S = (5.0000e-01, 0.0000e+00, 0.0000e+00) | - | - | exact | 2.000000000000e-07 |
| PASS | B | b4.magnetoelectric_directional_asymmetry | in a magnetoelectric medium the forward and backward wavenumbers differ in magnitude: phase velocity is direction dependent | n/a | n/a | \|kappa_+\| != \|kappa_-\| | \|kappa_+\| = 1.3622e+00, \|kappa_-\| = 2.2022e+00, omega = 1.0000e+00 | - | - | exact | 4.900000000000e-06 |
| PASS | B | b4.layered.lossless_slab_conserves_energy | R + T = 1 for a lossless dielectric slab | n/a | eps = 4, d = 1.7, 3 frequencies | 0 | 3.330669073875e-16 | 3.330669073875e-16 | - | 1.000000000000e-13 | 1.106000000000e-04 |
| PASS | B | b4.layered.lossless_slab_absorbs_nothing | the independently integrated dissipated power is zero in a lossless slab | n/a | 3 frequencies | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 1.129000000000e-04 |
| PASS | B | b4.layered.lossy_slab_energy_balance | R + T + A = 1 with A obtained by integrating <p> = (w/2) Im[E* D + H* B] through the internal field profile | n/a | eps = 4 + 0.35i, d = 1.7, 3 frequencies, 6000 quadrature points | 0 | 6.866251789361e-09 | 6.866251789361e-09 | - | 1.000000000000e-07 | 1.584800000000e-03 |
| PASS | B | b4.layered.absorption_is_positive_for_a_passive_slab | a passive lossy slab absorbs a strictly positive amount of power | n/a | n/a | A > 0 | R = 7.6170e-02, T = 3.5058e-01, A = 5.7325e-01 | - | - | exact | 1.588000000000e-03 |
| PASS | B | b4.layered.active_slab_is_detected | a slab with Im eps < 0 shows NEGATIVE absorption and R + T > 1 | n/a | n/a | gain detected, and the balance still closes | A = -4.0014e-01, R + T = 1.4001e+00 | - | - | exact | 3.640000000000e-05 |
| PASS | A | b4.layered.thick_evanescent_stack_is_stable | a 40-unit evanescent layer returns finite R and T with T exponentially small, rather than overflowing | n/a | n/a | finite results, T below 1e-30 | R = 9.9530e-01, T = 2.5936e-49 | - | - | exact | 3.040000000000e-05 |
| PASS | B | b4.propagation.restricted_family_costs_phase | dropping the magnetoelectric term changes measurable propagation observables -- phase accumulation and group delay -- through a warp-analogue profile, and both stacks conserve energy | n/a | n/a | a measurable phase difference, with energy conserved in both | exact: R = 2.4038e-32, T = 1.0000e+00, phase = -1.9224e+00; restricted (zeta = 0): R = 5.6635e-04, T = 9.9943e-01, phase = -3.1381e-01; d(phase) = 1.6086e+00, d(group delay) = 1.6097e+00 | - | - | exact | 3.366300000000e-03 |
| PASS | B | b4.adaptive.refines_the_passivity_margin | the adaptive scan finds a strictly LOWER minimum loss eigenvalue than the coarse grid, by refining around resonances, zeros of det mu and local minima | n/a | n/a | a strictly lower minimum found by refinement | coarse (21 samples): min eig = -3.1308e+01 at w = 3.0500e+00; adaptive (981 samples, 6 rounds): min eig = -4.9965e+02 at w = 3.0558e+00; resonances targeted: 3, mu zeros: 2 | - | - | exact | 3.622200000000e-03 |
| PASS | B | b4.hypad.propagation_phase_gradient | d(phase)/d(eps) through a slab matches the analytic w d /(2 sqrt(eps)) | n/a | eps = 4, d = 1.7, omega = 1.3; finite-difference cross-check 1.7275e-11 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 1.030000000000e-05 |
| PASS | A | b4.deterministic_repeatability | two identical stack solves return bit-identical R, T and A | n/a | n/a | bit-identical | identical | - | - | exact | 7.450000000000e-05 |
| PASS | A | b4.branch_count_is_reported | the number of real propagation branches along a direction is counted and reported rather than assumed | n/a | n/a | at least two real branches for an anisotropic medium | 2 real branches along +x: 1.8708e+00 2.2361e+00  | - | - | exact | 2.920000000000e-05 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
