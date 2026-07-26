# milestone_02t verification report

**Temporal derivative operators, and the rigid-translation identity d_t g = -v_s d_x g**

| | |
|---|---|
| Generated (UTC) | 2026-07-26T18:03:11Z |
| Git commit | `976ef4960578d7eece4aa22a861cbeb1c63c9144` |
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
| PASS | E | temporal.order.first.interior | d/dt converges at order 2 with the interior stencil | n/a | n/a | >= 1.85 (theoretical 2) | 1.999963 | - | - | exact | 7.790000000000e-05 |
| PASS | E | temporal.order.second.interior | d2/dt2 converges at order 2 with the interior stencil | n/a | n/a | >= 1.85 (theoretical 2) | 1.999975 | - | - | exact | 3.500000000000e-06 |
| PASS | E | temporal.order.first.window_start | d/dt converges at order 2 with the window_start stencil | n/a | n/a | >= 1.85 (theoretical 2) | 1.943961 | - | - | exact | 1.500000000000e-06 |
| PASS | E | temporal.order.second.window_start | d2/dt2 converges at order 2 with the window_start stencil | n/a | n/a | >= 1.85 (theoretical 2) | 2.004902 | - | - | exact | 1.500000000000e-06 |
| PASS | E | temporal.order.first.window_end | d/dt converges at order 2 with the window_end stencil | n/a | n/a | >= 1.85 (theoretical 2) | 2.049729 | - | - | exact | 1.400000000000e-06 |
| PASS | E | temporal.order.second.window_end | d2/dt2 converges at order 2 with the window_end stencil | n/a | n/a | >= 1.85 (theoretical 2) | 1.993976 | - | - | exact | 1.200000000000e-06 |
| PASS | B | temporal.mixed.dt_dx.interior | d2f/dt dx against the closed form, interior points | 21x21x21 | h_t = 0.01, h_s = 0.1 | 0 | 1.477151686861e-03 | 1.477151686861e-03 | - | 5.000000000000e-03 | 1.365600000000e-03 |
| PASS | B | temporal.mixed.dt_dx.boundary | d2f/dt dx against the closed form, boundary points | 21x21x21 | h_t = 0.01, h_s = 0.1 | 0 | 2.909883438719e-03 | 2.909883438719e-03 | - | 1.000000000000e-02 | 1.436100000000e-03 |
| PASS | B | temporal.mixed.dt_dy.interior | d2f/dt dy against the closed form, interior points | 21x21x21 | h_t = 0.01, h_s = 0.1 | 0 | 6.186314770118e-05 | 6.186314770118e-05 | - | 5.000000000000e-03 | 1.460900000000e-03 |
| PASS | B | temporal.mixed.dt_dy.boundary | d2f/dt dy against the closed form, boundary points | 21x21x21 | h_t = 0.01, h_s = 0.1 | 0 | 1.473354082616e-04 | 1.473354082616e-04 | - | 1.000000000000e-02 | 1.485900000000e-03 |
| PASS | B | temporal.mixed.dt_dz.interior | d2f/dt dz against the closed form, interior points | 21x21x21 | h_t = 0.01, h_s = 0.1 | 0 | 2.393857045022e-05 | 2.393857045022e-05 | - | 5.000000000000e-03 | 1.091800000000e-03 |
| PASS | B | temporal.mixed.dt_dz.boundary | d2f/dt dz against the closed form, boundary points | 21x21x21 | h_t = 0.01, h_s = 0.1 | 0 | 6.509300063162e-05 | 6.509300063162e-05 | - | 1.000000000000e-02 | 1.119900000000e-03 |
| PASS | A | temporal.mixed_commutes | d_t(D_a f) == D_a(d_t f) at every point and on every axis | 21x21x21 | h_t = 0.01 | 0 | 5.684341886081e-13 | 5.684341886081e-13 | - | 1.000000000000e-12 | 4.852700000000e-03 |
| PASS | A | temporal.constant_field_zero | d/dt and d2/dt2 of a time-constant field are zero, for all three stencil positions | n/a | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 1.000000000000e-06 |
| PASS | A | temporal.linear_field_exact | d/dt of a field linear in t is exact and d2/dt2 vanishes | n/a | n/a | 0 | 3.637978807092e-12 | 3.637978807092e-12 | - | 1.000000000000e-11 | 6.000000000000e-07 |
| PASS | B | temporal.tensor_field | every rank-2 component is differentiated in time against its own closed form, into the alpha = 0 slot | 11x11x11 | h_t = 0.01 | 0 | 4.852773459341e-05 | 4.852773459341e-05 | - | 1.000000000000e-04 | 3.399000000000e-04 |
| PASS | A | temporal.generic.test_scalar_bit_exact | the temporal operators instantiate with a custom scalar carrying no ordering, no equality and no implicit double conversion, and give bit-identical results | n/a | n/a | bit-identical | bit-identical | - | - | exact | 6.000000000000e-07 |
| PASS | D | temporal.generic.long_double_agrees | long double agrees with double to double precision | n/a | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 2.300000000000e-06 |
| PASS | D | temporal.rigid_translation_identity | d_t g_{mu nu} == -v_s d_x g_{mu nu} for a rigidly translating bubble: the residual between analytic d_t and numerical d_x converges to zero at the order of the spatial stencil | n/a | n/a | order >= 1.70 and residual falling by more than 4x over 4x refinement | 1.824772 | - | - | exact | 6.735080000000e-02 |
| PASS | A | temporal.zero_velocity_static | v_s = 0 gives d_t g == 0 exactly | n/a | v = 0 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.480000000000e-05 |
| PASS | A | temporal.frozen_mode_reports_zero | FrozenTimeDiagnostic reports d_t g == 0 exactly | n/a | v=0.35, t=2 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 3.700000000000e-06 |
| PASS | A | temporal.frozen_differs_from_moving | at t = 2 with v_s = 0.35 the frozen slice and the time-dependent metric DISAGREE, quantifying what the diagnostic mode discards | n/a | n/a | substantial disagreement (> 0.1) | max component difference 0.199773 | - | - | exact | 6.400000000000e-06 |
| PASS | D | temporal.oti_through_time_derivative | d/dp (d_t g_00) AT FIXED BUBBLE TRAJECTORY, from OTI pushed through the temporal stencil, matches a centred parameter difference of the same quantity computed in double | n/a | v=0.35, R=3, sigma=2, t=1.1, sampled at r = R (the wall) | 0 | 4.483150132910e-08 | 4.483150132910e-08 | - | 1.000000000000e-06 | 3.680000000000e-05 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
