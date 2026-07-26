# otilib verification report

**OTIlib onumm3n2 adapter: 3 imaginary directions, truncation order 2**

| | |
|---|---|
| Generated (UTC) | 2026-07-25T01:09:35Z |
| Git commit | `unknown` **(working tree dirty)** |
| Compiler | MSVC 19.44.35223.0 |
| Build type | RelWithDebInfo |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 23 |
| Passed | 23 |
| Failed | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | B | oti.elementary.sqrt | first and second derivatives against the closed form | n/a | x=1.3, d/dx=0.438529, d2/dx2=-0.168665 | 0 | 2.775557561563e-17 | 2.775557561563e-17 | - | 1.000000000000e-12 | 8.740000000000e-05 |
| PASS | B | oti.elementary.exp | first and second derivatives against the closed form | n/a | x=1.3, d/dx=3.6693, d2/dx2=3.6693 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 6.300000000000e-06 |
| PASS | B | oti.elementary.log | first and second derivatives against the closed form | n/a | x=1.3, d/dx=0.769231, d2/dx2=-0.591716 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 4.300000000000e-06 |
| PASS | B | oti.elementary.sin | first and second derivatives against the closed form | n/a | x=1.3, d/dx=0.267499, d2/dx2=-0.963558 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 2.700000000000e-06 |
| PASS | B | oti.elementary.cos | first and second derivatives against the closed form | n/a | x=1.3, d/dx=-0.963558, d2/dx2=-0.267499 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 2.900000000000e-06 |
| PASS | B | oti.elementary.tan | first and second derivatives against the closed form | n/a | x=1.3, d/dx=13.9751, d2/dx2=100.68 | 0 | 1.411490347209e-16 | 1.411490347209e-16 | - | 1.000000000000e-12 | 6.600000000000e-06 |
| PASS | B | oti.elementary.sinh | first and second derivatives against the closed form | n/a | x=1.3, d/dx=1.97091, d2/dx2=1.69838 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 3.000000000000e-06 |
| PASS | B | oti.elementary.cosh | first and second derivatives against the closed form | n/a | x=1.3, d/dx=1.69838, d2/dx2=1.97091 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 2.100000000000e-06 |
| PASS | B | oti.elementary.tanh | first and second derivatives against the closed form | n/a | x=1.3, d/dx=0.257433, d2/dx2=-0.443672 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 3.900000000000e-06 |
| PASS | B | oti.elementary.divide | d(a/b)/da == 1/b and d2(a/b)/da2 == 0 | n/a | a=1.3, b=0.7 constant | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 1.700000000000e-06 |
| PASS | B | oti.composite.value | f(a,b,c) = a^2 sin(b) + exp(ac), real part | n/a | a=1.3, b=0.7, c=0.45 | 2.883718877072e+00 | 2.883718877072e+00 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-14 | 7.000000000000e-07 |
| PASS | B | oti.composite.gradient | all three first derivatives of a^2 sin(b) + exp(ac) | n/a | a=1.3, b=0.7, c=0.45 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 7.000000000000e-07 |
| PASS | B | oti.composite.hessian | all six independent second derivatives of a^2 sin(b) + exp(ac) | n/a | a=1.3, b=0.7, c=0.45 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 1.900000000000e-06 |
| PASS | A | oti.hessian_symmetry | H(i,j) == H(j,i) | n/a | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 4.000000000000e-07 |
| PASS | A | oti.seed_isolation | seeding direction k leaves every other derivative channel exactly 0 | n/a | 3 runs, one seed each | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 8.000000000000e-07 |
| PASS | A | oti.simultaneous_vs_separate | 3 directions seeded at once give the same derivatives as 3 single-seed runs | n/a | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 3.000000000000e-06 |
| PASS | A | oti.zero_seed_matches_double | with all seeds zero the OTI real part equals the double result | n/a | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 4.000000000000e-07 |
| PASS | A | oti.aliasing_is_unsafe | OTIlib *_to functions produce different results when the output aliases an input, confirming the wrapper's no-alias rule is needed | n/a | n/a | aliased result differs from the clean result | differs, as expected | - | - | exact | 1.000000000000e-07 |
| PASS | A | oti.generic_matrix_inverse | Mat4/inverse/identity_residual instantiate and run with OTIScalar | n/a | n/a | 0 | 2.775557561563e-17 | 2.775557561563e-17 | - | 1.000000000000e-12 | 1.290000000000e-05 |
| PASS | D | oti.matrix_inverse_derivative_identity | d(A^-1)/dp == -A^-1 (dA/dp) A^-1 for all three parameters | n/a | a=1.3, b=0.7, c=0.45 | 0 | 1.884885476866e-16 | 1.884885476866e-16 | - | 1.000000000000e-11 | 7.800000000000e-06 |
| PASS | D | oti.vs_finite_difference.first | the best centered-FD first derivative agrees with OTI | n/a | a=1.3, b=0.7, c=0.45 | 0 | 1.615667020779e-11 | 1.615667020779e-11 | - | 1.000000000000e-09 | 3.500000000000e-06 |
| PASS | D | oti.vs_finite_difference.second | the best centered-FD second derivative agrees with OTI | n/a | a=1.3, b=0.7, c=0.45 | 0 | 3.356945284952e-09 | 3.356945284952e-09 | - | 1.000000000000e-05 | 1.600000000000e-06 |
| PASS | D | oti.fd_error_is_v_shaped | FD error falls then rises as h shrinks (truncation then roundoff), confirming FD is a genuinely independent method | n/a | n/a | error at h=1e-5 below both ends | h=1e-1: 2.726e-04, h=1e-5: 4.011e-11, h=1e-12: 1.908e-04 | - | - | exact | 3.000000000000e-06 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
