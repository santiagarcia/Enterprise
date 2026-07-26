# milestone_01 verification report

**Grid, contiguous storage, fixed-size tensors and 4x4 matrix algebra**

| | |
|---|---|
| Generated (UTC) | 2026-07-26T14:25:37Z |
| Git commit | `611053f22b136ad570115f0535ab02e677e46fef` |
| Compiler | MSVC 19.44.35223.0 |
| Build type | RelWithDebInfo |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 77 |
| Passed | 77 |
| Failed | 0 |
| Skipped | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | A | grid.lower_bound_exact [41x41x41] | x(0) == x_min for every axis | 41x41x41 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 6.000000000000e-07 |
| PASS | A | grid.upper_bound_exact [41x41x41] | x(n-1) == x_max for every axis | 41x41x41 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 2.000000000000e-07 |
| PASS | A | grid.spacing_uniform [41x41x41] | max_i \|(x_{i+1} - x_i) - dx\| / dx == 0 | 41x41x41 | n/a | 0 | 3.552713678801e-15 | 3.552713678801e-15 | - | 1.000000000000e-12 | 3.500000000000e-06 |
| PASS | A | grid.center_coordinate_exact [41x41x41] | x((n-1)/2) == (x_min + x_max)/2 for odd n | 41x41x41 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 2.000000000000e-06 |
| PASS | A | grid.index_roundtrip_ijk [41x41x41] | (i,j,k) -> flat -> (i,j,k) is the identity at every point | 41x41x41 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 7.756000000000e-04 |
| PASS | A | grid.index_roundtrip_flat [41x41x41] | flat -> (i,j,k) -> flat is the identity for every flat index | 41x41x41 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 8.188000000000e-04 |
| PASS | A | grid.flat_index_bijective [41x41x41] | flat() is a bijection from the index box onto [0, N) | 41x41x41 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 3.731000000000e-04 |
| PASS | A | grid.strides_match_flat [41x41x41] | stride_a() == flat(...+1...) - flat(...) for each axis | 41x41x41 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 2.000000000000e-07 |
| PASS | A | grid.dx_value [41x41x41] | dx == (x_max - x_min) / (nx - 1) | 41x41x41 | n/a | 5.000000000000e-01 | 5.000000000000e-01 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-15 | 0.000000000000e+00 |
| PASS | A | grid.cell_volume [41x41x41] | cell_volume == dx * dy * dz | 41x41x41 | n/a | 1.250000000000e-01 | 1.250000000000e-01 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-15 | 0.000000000000e+00 |
| PASS | A | grid.num_points [41x41x41] | num_points == nx * ny * nz | 41x41x41 | n/a | 6.892100000000e+04 | 6.892100000000e+04 | 0.000000000000e+00 | - | exact | 1.000000000000e-07 |
| PASS | A | grid.lower_bound_exact [17x9x5] | x(0) == x_min for every axis | 17x9x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 3.000000000000e-07 |
| PASS | A | grid.upper_bound_exact [17x9x5] | x(n-1) == x_max for every axis | 17x9x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 2.000000000000e-07 |
| PASS | A | grid.spacing_uniform [17x9x5] | max_i \|(x_{i+1} - x_i) - dx\| / dx == 0 | 17x9x5 | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 6.000000000000e-07 |
| PASS | A | grid.center_coordinate_exact [17x9x5] | x((n-1)/2) == (x_min + x_max)/2 for odd n | 17x9x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.900000000000e-06 |
| PASS | A | grid.index_roundtrip_ijk [17x9x5] | (i,j,k) -> flat -> (i,j,k) is the identity at every point | 17x9x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 8.600000000000e-06 |
| PASS | A | grid.index_roundtrip_flat [17x9x5] | flat -> (i,j,k) -> flat is the identity for every flat index | 17x9x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 9.100000000000e-06 |
| PASS | A | grid.flat_index_bijective [17x9x5] | flat() is a bijection from the index box onto [0, N) | 17x9x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 3.100000000000e-06 |
| PASS | A | grid.strides_match_flat [17x9x5] | stride_a() == flat(...+1...) - flat(...) for each axis | 17x9x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.000000000000e-07 |
| PASS | A | grid.dx_value [17x9x5] | dx == (x_max - x_min) / (nx - 1) | 17x9x5 | n/a | 6.250000000000e-01 | 6.250000000000e-01 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-15 | 0.000000000000e+00 |
| PASS | A | grid.cell_volume [17x9x5] | cell_volume == dx * dy * dz | 17x9x5 | n/a | 1.562500000000e-01 | 1.562500000000e-01 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-15 | 0.000000000000e+00 |
| PASS | A | grid.num_points [17x9x5] | num_points == nx * ny * nz | 17x9x5 | n/a | 7.650000000000e+02 | 7.650000000000e+02 | 0.000000000000e+00 | - | exact | 0.000000000000e+00 |
| PASS | A | grid.validation_rejects_malformed | validate_grid_spec rejects n<1, inverted bounds, non-finite bounds and unallocatable sizes | n/a | n/a | all four malformed specs rejected | all four rejected | - | - | exact | 9.700000000000e-06 |
| PASS | A | grid.validation_accepts_valid | validate_grid_spec accepts both reference specs | n/a | n/a | both accepted | both accepted | - | - | exact | 1.700000000000e-06 |
| PASS | A | layout.rank2_roundtrip | every rank-2 component stores and retrieves its own index tag exactly | 11x7x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 6.320000000000e-05 |
| PASS | A | layout.rank2_documented_offset | raw()[16*p + 4*mu + nu] holds component (mu,nu) at point p | 11x7x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 4.100000000000e-05 |
| PASS | A | layout.rank3_roundtrip_and_offset | rank-3 accessor and raw()[64*p + 16*a + 4*b + c] agree with the index tag everywhere | 11x7x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 2.099000000000e-04 |
| PASS | A | layout.rank4_roundtrip_and_offset | rank-4 accessor and raw()[256*p + 64*a + 16*b + 4*c + d] agree with the index tag everywhere | 11x7x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 8.930000000000e-04 |
| PASS | A | layout.scalar_field_roundtrip | ScalarField (i,j,k) and flat accessors agree and store exactly | 11x7x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 3.000000000000e-06 |
| PASS | A | layout.mat4_field_conversion | set_matrix followed by matrix_at is the identity, with no transposition | 11x7x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 4.660000000000e-05 |
| PASS | A | layout.rank3_tensor_conversion | set_tensor followed by tensor_at is the identity for rank 3 | 11x7x5 | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 8.680000000000e-05 |
| PASS | A | algebra.matmul_hand_calculated | diag(1,2,3,4)*B scales rows; B*diag(1,2,3,4) scales columns | n/a | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.060000000000e-05 |
| PASS | D | algebra.matmul_vs_independent | mat_mul agrees with an independently written naive triple loop | n/a | entries in [-2, 2] | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 1.758000000000e-04 |
| PASS | B | algebra.det_analytic | det of Minkowski (-1), identity (1), diag(2,3,5,7) (210), an upper-triangular matrix (24) and a rank-deficient matrix (0) | n/a | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 5.000000000000e-07 |
| PASS | D | algebra.det_vs_independent | the paired-2x2-minor determinant agrees with recursive Laplace expansion along row 0 | n/a | entries in [-2, 2] | 0 | 2.664535259100e-15 | 2.664535259100e-15 | - | 1.000000000000e-12 | 2.003800000000e-03 |
| PASS | A | algebra.inverse_identity_left [diagonal] | A * A^-1 == I | n/a | diagonal matrices with mixed signs | 0 | 3.330669073875e-16 | 3.330669073875e-16 | - | 1.000000000000e-09 | 2.219500000000e-03 |
| PASS | A | algebra.inverse_identity_right [diagonal] | A^-1 * A == I | n/a | diagonal matrices with mixed signs | 0 | 3.330669073875e-16 | 3.330669073875e-16 | - | 1.000000000000e-09 | 2.226100000000e-03 |
| PASS | A | algebra.inverse_identity_left [symmetric] | A * A^-1 == I | n/a | general symmetric matrices | 0 | 3.126388037344e-13 | 3.126388037344e-13 | - | 1.000000000000e-09 | 2.214700000000e-03 |
| PASS | A | algebra.inverse_identity_right [symmetric] | A^-1 * A == I | n/a | general symmetric matrices | 0 | 3.410605131648e-13 | 3.410605131648e-13 | - | 1.000000000000e-09 | 2.218100000000e-03 |
| PASS | A | algebra.inverse_identity_left [lorentzian] | A * A^-1 == I | n/a | perturbed Lorentzian metrics, signature -+++ | 0 | 6.661338147751e-16 | 6.661338147751e-16 | - | 1.000000000000e-09 | 3.304700000000e-03 |
| PASS | A | algebra.inverse_identity_right [lorentzian] | A^-1 * A == I | n/a | perturbed Lorentzian metrics, signature -+++ | 0 | 6.661338147751e-16 | 6.661338147751e-16 | - | 1.000000000000e-09 | 3.331900000000e-03 |
| PASS | A | algebra.inverse_identity_left [general] | A * A^-1 == I | n/a | general nonsingular matrices | 0 | 8.526512829121e-14 | 8.526512829121e-14 | - | 1.000000000000e-09 | 2.698300000000e-03 |
| PASS | A | algebra.inverse_identity_right [general] | A^-1 * A == I | n/a | general nonsingular matrices | 0 | 5.684341886081e-14 | 5.684341886081e-14 | - | 1.000000000000e-09 | 2.706800000000e-03 |
| PASS | D | algebra.inverse_vs_independent | the adjugate inverse agrees with Gauss-Jordan elimination | n/a | perturbed Lorentzian metrics | 0 | 6.661338147751e-16 | 6.661338147751e-16 | - | 1.000000000000e-09 | 1.791000000000e-04 |
| PASS | B | algebra.minkowski_self_inverse | eta^{mu nu} == eta_{mu nu} == diag(-1,1,1,1) | n/a | signature -+++ | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.200000000000e-06 |
| PASS | A | algebra.singular_detected | inverse() reports ok == false for exactly singular matrices (row 3 = 2*row 0 - 3*row 1) | n/a | n/a | 0 accepted | 0 of 200 exactly singular matrices were accepted | - | - | exact | 4.180000000000e-05 |
| PASS | A | algebra.near_singular_reported | a near-singular matrix is either rejected or flagged ill-conditioned; it is never quietly accepted with a bad residual | n/a | n/a | 0 quietly accepted with a bad residual | 200 rejected, 0 flagged ill-conditioned, 0 accepted quietly with a bad residual | - | - | exact | 5.080000000000e-05 |
| PASS | A | algebra.ill_conditioning_flag_fires | across perturbations from 1e-4 to 1e-14 the ill_conditioned flag is set for at least one accepted matrix, and no matrix at any scale is quietly accepted with a bad residual | n/a | n/a | at least one flagged, zero quietly accepted with a bad residual | eps=0.0001: 100 ok (4 flagged, 0 quietly bad, worst unflagged residual 2.46108e-09 at cond1 1.02718e+07), 0 rejected; eps=1e-06: 100 ok (85 flagged, 0 quietly bad, worst unflagged residual 1.86265e-09 at cond1 3.81862e+07), 0 rejected; eps=1e-08: 96 ok (96 flagged, 0 quietly bad, worst unflagged residual 0 at cond1 0), 4 rejected; eps=1e-10: 17 ok (17 flagged, 0 quietly bad, worst unflagged residual 0 at cond1 0), 83 rejected; eps=1e-12: 0 ok (0 flagged, 0 quietly bad, worst unflagged residual 0 at cond1 0), 100 rejected; eps=1e-14: 0 ok (0 flagged, 0 quietly bad, worst unflagged residual 0 at cond1 0), 100 rejected | - | - | exact | 1.742000000000e-04 |
| PASS | A | algebra.condition_number_consistent | cond1 == \|\|A\|\|_1 * \|\|A^-1\|\|_1, and cond1 >= 1 | n/a | general matrices | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 1.925000000000e-04 |
| PASS | A | contract.trace_mixed | trace of diag(1,2,3,4) == 10 | n/a | n/a | 1.000000000000e+01 | 1.000000000000e+01 | 0.000000000000e+00 | 0.000000000000e+00 | exact | 4.000000000000e-07 |
| PASS | B | contract.trace_minkowski | trace of diag(-1,1,1,1) == 2 | n/a | signature -+++ | 2.000000000000e+00 | 2.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | exact | 2.000000000000e-07 |
| PASS | A | contract.full_contract | sum_{mu,nu} (mu + nu) * 1 == 48 | n/a | n/a | 4.800000000000e+01 | 4.800000000000e+01 | 0.000000000000e+00 | 0.000000000000e+00 | exact | 2.000000000000e-07 |
| PASS | A | contract.metric_full_trace | eta^{mu nu} eta_{mu nu} == 4 (dimension of spacetime) | n/a | signature -+++ | 4.000000000000e+00 | 4.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | exact | 2.000000000000e-06 |
| PASS | A | contract.raise_lower_roundtrip | eta^{mu a} eta_{a b} T^{b}{}_{nu} == T^{mu}{}_{nu} | n/a | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 4.000000000000e-07 |
| PASS | A | contract.contract_first_diagonal | for diagonal a, (a . d)(r,m,n) == a(r,r) * d(r,m,n) at all 64 components; spot value out(1,2,3) == 369 | n/a | a = diag(2,3,5,7) | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.000000000000e-06 |
| PASS | A | contract.contract_first_spot_value | out(1,2,3) == 3 * 123 == 369 | n/a | a = diag(2,3,5,7) | 3.690000000000e+02 | 3.690000000000e+02 | 0.000000000000e+00 | 0.000000000000e+00 | exact | 8.000000000000e-07 |
| PASS | A | contract.contract_13 | out(s,n) == sum_r t(r,s,r,n); with the constructed t this is 4*(10s+n), so out(2,3) == 92 | n/a | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 9.000000000000e-07 |
| PASS | A | contract.contract_13_slot_selectivity | a rank-4 vanishing wherever slot 0 == slot 2 contracts to zero under contract_13 | n/a | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 6.000000000000e-07 |
| PASS | A | contract.lower_first_diagonal | for diagonal g, (g . t)(a,b,c,d) == g(a,a) * t(a,b,c,d) at all 256 components | n/a | g = diag(2,3,5,7) | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 4.500000000000e-06 |
| PASS | A | contract.symmetry_detector_works | lower_pair_symmetry_residual returns 0 for a b<->c symmetric tensor and 1 for a deliberately antisymmetric one | n/a | n/a | 0 for symmetric, 1 for the planted asymmetry | symmetric: 0.000000, planted: 1.000000 | - | - | exact | 5.000000000000e-07 |
| PASS | A | contract.mat4_symmetry_detector_works | symmetry_residual returns 0 for a symmetric matrix and 3 for a planted asymmetry of 3 | n/a | n/a | 0 for symmetric, 3 for the planted asymmetry | symmetric: 0.000000, planted: 3.000000 | - | - | exact | 2.000000000000e-07 |
| PASS | A | generic.compiles_with_custom_scalar | Mat4, Rank3, Rank4 and every contraction instantiate with double, long double and a custom scalar that has no ordering, no equality and no implicit conversion to double | n/a | n/a | all three instantiations compile | all three instantiations compiled and ran | - | - | exact | 2.000000000000e-07 |
| PASS | A | generic.inverse_ok_agrees | the singularity verdict does not depend on the scalar type | n/a | n/a | all three report ok == true | double=ok, long double=ok, TestScalar=ok | - | - | exact | 0.000000000000e+00 |
| PASS | A | generic.double_vs_custom_scalar_bit_exact | TestScalar reproduces the double results bit for bit, since it performs the same operations in the same order | n/a | n/a | bit-identical results | bit-identical | - | - | exact | 1.000000000000e-07 |
| PASS | D | generic.long_double_agrees | the determinant computed in long double agrees with the double result to double precision | n/a | perturbed Lorentzian test matrix | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-14 | 5.000000000000e-07 |
| PASS | A | generic.inverse_residual_custom_scalar | A * A^-1 == I holds when computed entirely in TestScalar | n/a | perturbed Lorentzian test matrix | 0 | 2.220446049250e-16 | 2.220446049250e-16 | - | 1.000000000000e-12 | 0.000000000000e+00 |
| PASS | A | generic.fields_with_custom_scalar | Rank2Field instantiates with a custom scalar and produces bit-identical results to double | n/a | n/a | bit-identical accumulations | bit-identical | - | - | exact | 1.370000000000e-05 |
| PASS | A | generic.scalar_traits_roundtrip | real_part(scalar_cast<S>(x)) == x for every supported scalar type | n/a | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.000000000000e-07 |
| PASS | A | config.reference_parses | the specified reference configuration parses, validates and produces no warnings | n/a | n/a | parsed with 0 warnings | parsed with 0 warnings | - | - | exact | 9.090000000000e-05 |
| PASS | A | config.values_correct | every field of the reference configuration parses to its specified value, exactly | n/a | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 4.000000000000e-07 |
| PASS | A | config.roundtrip_byte_exact | dump(to_json(from_json(dump(to_json(cfg))))) reproduces the identical byte stream | n/a | n/a | identical byte streams | identical | - | - | exact | 1.536000000000e-04 |
| PASS | A | config.number_roundtrip_bit_exact | 0.1+0.2, 1/3, a subnormal-scale value and -DBL_MAX all survive write-then-read bit for bit | n/a | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.141000000000e-04 |
| PASS | A | config.unknown_keys_warn | a misspelled key produces a warning rather than silently taking the default value | n/a | n/a | 2 warnings | 2 warnings: unrecognised key 'unknown_section' ignored; unrecognised key 'metric.sigmaa' ignored | - | - | exact | 1.330000000000e-05 |
| PASS | A | config.invalid_rejected | unknown metric type, wrong value types, a non-integer grid size, a non-object root and malformed JSON are all rejected | n/a | n/a | 0 of 5 accepted | 0 of 5 accepted | - | - | exact | 1.830000000000e-05 |
| PASS | A | config.validation_rejects_unsupported | negative radius, zero sigma, zero epsilon, an unimplemented signature and an unimplemented derivative order are rejected | n/a | n/a | 0 of 5 accepted | 0 of 5 accepted | - | - | exact | 5.900000000000e-06 |
| PASS | A | config.partial_inherits_defaults | a partial config overrides only what it names and inherits the reference defaults elsewhere | n/a | n/a | velocity 0.5, radius 3, nx 41 | velocity 0.500000, radius 3.000000, nx 41 | - | - | exact | 4.600000000000e-06 |
| PASS | A | config.json_parser_edge_cases | the parser accepts well-formed JSON and rejects trailing commas, trailing content, unquoted keys and unterminated strings | n/a | n/a | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.990000000000e-05 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
