# milestone_02 verification report

**Finite-difference operators: first, second and mixed derivatives at second order, interior and boundary**

| | |
|---|---|
| Generated (UTC) | 2026-07-26T16:48:53Z |
| Git commit | `1208b5c04d3bd81a382921faef75f1eb71a3365e` |
| Compiler | MSVC 19.44.35223.0 |
| Build type | RelWithDebInfo |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 41 |
| Passed | 41 |
| Failed | 0 |
| Skipped | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | A | fd.exact.constant.first | d/dx of a constant is exactly 0 on every axis | 21x21x21 | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-13 | 2.528600000000e-03 |
| PASS | A | fd.exact.constant.second | second derivative of a constant field is exactly 0 | 21x21x21 | n/a | 0 | 5.684341886081e-14 | 5.684341886081e-14 | - | 1.000000000000e-11 | 2.643000000000e-03 |
| PASS | A | fd.exact.linear.first | d/dx of a linear field is exact on every axis | 21x21x21 | n/a | 0 | 3.197442310920e-14 | 3.197442310920e-14 | - | 1.000000000000e-13 | 2.203300000000e-03 |
| PASS | A | fd.exact.linear.second | second derivative of a linear field is exactly 0 | 21x21x21 | n/a | 0 | 9.094947017729e-13 | 9.094947017729e-13 | - | 1.000000000000e-11 | 2.210800000000e-03 |
| PASS | E | fd.order.interior.polynomial.d/dx | observed interior convergence order approaches 2 | n/a | n/a | >= 1.85 (theoretical 2) | 1.962525 | - | - | exact | 3.250000000000e-05 |
| PASS | E | fd.order.boundary.polynomial.d/dx | observed boundary convergence order approaches 2, i.e. the one-sided stencils are the same order as the interior ones | n/a | n/a | >= 1.85 (theoretical 2) | 1.972166 | - | - | exact | 4.240000000000e-05 |
| PASS | E | fd.order.interior.polynomial.d2/dx2 | observed interior convergence order approaches 2 | n/a | n/a | >= 1.85 (theoretical 2) | 2.000000 | - | - | exact | 2.190000000000e-05 |
| PASS | E | fd.order.boundary.polynomial.d2/dx2 | observed boundary convergence order approaches 2, i.e. the one-sided stencils are the same order as the interior ones | n/a | n/a | >= 1.85 (theoretical 2) | 2.000000 | - | - | exact | 3.260000000000e-05 |
| PASS | E | fd.order.interior.polynomial.d/dy | observed interior convergence order approaches 2 | n/a | n/a | >= 1.85 (theoretical 2) | 2.000000 | - | - | exact | 2.770000000000e-05 |
| PASS | E | fd.order.boundary.polynomial.d/dy | observed boundary convergence order approaches 2, i.e. the one-sided stencils are the same order as the interior ones | n/a | n/a | >= 1.85 (theoretical 2) | 2.000000 | - | - | exact | 3.620000000000e-05 |
| PASS | A | fd.exact_operator.polynomial.d2/dy2 | the stencil is exact for this field (the derivative it truncates against vanishes identically), so the error stays at roundoff on every grid | n = 11..81 | polynomial | 0 | 2.000888343900e-11 | 2.000888343900e-11 | - | 1.000000000000e-08 | 1.800000000000e-06 |
| PASS | A | fd.exact_operator.polynomial.d/dz | the stencil is exact for this field (the derivative it truncates against vanishes identically), so the error stays at roundoff on every grid | n = 11..81 | polynomial | 0 | 1.723066134218e-13 | 1.723066134218e-13 | - | 1.000000000000e-08 | 1.500000000000e-06 |
| PASS | A | fd.exact_operator.polynomial.d2/dz2 | the stencil is exact for this field (the derivative it truncates against vanishes identically), so the error stays at roundoff on every grid | n = 11..81 | polynomial | 0 | 2.728484105319e-11 | 2.728484105319e-11 | - | 1.000000000000e-08 | 2.100000000000e-06 |
| PASS | A | fd.exact_operator.polynomial.d2/dxdy | the stencil is exact for this field (the derivative it truncates against vanishes identically), so the error stays at roundoff on every grid | n = 11..81 | polynomial | 0 | 7.423839321063e-12 | 7.423839321063e-12 | - | 1.000000000000e-08 | 2.300000000000e-06 |
| PASS | A | fd.exact_operator.polynomial.d2/dxdz | the stencil is exact for this field (the derivative it truncates against vanishes identically), so the error stays at roundoff on every grid | n = 11..81 | polynomial | 0 | 1.026601026410e-11 | 1.026601026410e-11 | - | 1.000000000000e-08 | 1.400000000000e-06 |
| PASS | A | fd.exact_operator.polynomial.d2/dydz | the stencil is exact for this field (the derivative it truncates against vanishes identically), so the error stays at roundoff on every grid | n = 11..81 | polynomial | 0 | 1.015232342638e-11 | 1.015232342638e-11 | - | 1.000000000000e-08 | 2.000000000000e-06 |
| PASS | E | fd.order.interior.trigonometric.d/dx | observed interior convergence order approaches 2 | n/a | n/a | >= 1.85 (theoretical 2) | 1.999588 | - | - | exact | 2.670000000000e-05 |
| PASS | E | fd.order.boundary.trigonometric.d/dx | observed boundary convergence order approaches 2, i.e. the one-sided stencils are the same order as the interior ones | n/a | n/a | >= 1.85 (theoretical 2) | 2.025824 | - | - | exact | 3.490000000000e-05 |
| PASS | E | fd.order.interior.trigonometric.d2/dx2 | observed interior convergence order approaches 2 | n/a | n/a | >= 1.85 (theoretical 2) | 1.993405 | - | - | exact | 2.740000000000e-05 |
| PASS | E | fd.order.boundary.trigonometric.d2/dx2 | observed boundary convergence order approaches 2, i.e. the one-sided stencils are the same order as the interior ones | n/a | n/a | >= 1.85 (theoretical 2) | 1.992334 | - | - | exact | 3.640000000000e-05 |
| PASS | E | fd.order.interior.trigonometric.d/dy | observed interior convergence order approaches 2 | n/a | n/a | >= 1.85 (theoretical 2) | 1.998717 | - | - | exact | 7.230000000000e-05 |
| PASS | E | fd.order.boundary.trigonometric.d/dy | observed boundary convergence order approaches 2, i.e. the one-sided stencils are the same order as the interior ones | n/a | n/a | >= 1.85 (theoretical 2) | 1.973216 | - | - | exact | 1.017000000000e-04 |
| PASS | E | fd.order.interior.trigonometric.d2/dy2 | observed interior convergence order approaches 2 | n/a | n/a | >= 1.85 (theoretical 2) | 1.999739 | - | - | exact | 3.440000000000e-05 |
| PASS | E | fd.order.boundary.trigonometric.d2/dy2 | observed boundary convergence order approaches 2, i.e. the one-sided stencils are the same order as the interior ones | n/a | n/a | >= 1.85 (theoretical 2) | 1.912710 | - | - | exact | 4.380000000000e-05 |
| PASS | E | fd.order.interior.trigonometric.d/dz | observed interior convergence order approaches 2 | n/a | n/a | >= 1.85 (theoretical 2) | 1.971233 | - | - | exact | 2.760000000000e-05 |
| PASS | E | fd.order.boundary.trigonometric.d/dz | observed boundary convergence order approaches 2, i.e. the one-sided stencils are the same order as the interior ones | n/a | n/a | >= 1.85 (theoretical 2) | 1.978478 | - | - | exact | 3.480000000000e-05 |
| PASS | E | fd.order.interior.trigonometric.d2/dz2 | observed interior convergence order approaches 2 | n/a | n/a | >= 1.85 (theoretical 2) | 1.971204 | - | - | exact | 3.020000000000e-05 |
| PASS | E | fd.order.boundary.trigonometric.d2/dz2 | observed boundary convergence order approaches 2, i.e. the one-sided stencils are the same order as the interior ones | n/a | n/a | >= 1.85 (theoretical 2) | 1.968746 | - | - | exact | 4.190000000000e-05 |
| PASS | E | fd.order.interior.trigonometric.d2/dxdy | observed interior convergence order approaches 2 | n/a | n/a | >= 1.85 (theoretical 2) | 1.998152 | - | - | exact | 2.060000000000e-05 |
| PASS | E | fd.order.boundary.trigonometric.d2/dxdy | observed boundary convergence order approaches 2, i.e. the one-sided stencils are the same order as the interior ones | n/a | n/a | >= 1.85 (theoretical 2) | 1.990122 | - | - | exact | 3.330000000000e-05 |
| PASS | E | fd.order.interior.trigonometric.d2/dxdz | observed interior convergence order approaches 2 | n/a | n/a | >= 1.85 (theoretical 2) | 1.971066 | - | - | exact | 2.790000000000e-05 |
| PASS | E | fd.order.boundary.trigonometric.d2/dxdz | observed boundary convergence order approaches 2, i.e. the one-sided stencils are the same order as the interior ones | n/a | n/a | >= 1.85 (theoretical 2) | 1.988453 | - | - | exact | 3.450000000000e-05 |
| PASS | E | fd.order.interior.trigonometric.d2/dydz | observed interior convergence order approaches 2 | n/a | n/a | >= 1.85 (theoretical 2) | 1.970098 | - | - | exact | 2.210000000000e-05 |
| PASS | E | fd.order.boundary.trigonometric.d2/dydz | observed boundary convergence order approaches 2, i.e. the one-sided stencils are the same order as the interior ones | n/a | n/a | >= 1.85 (theoretical 2) | 1.947364 | - | - | exact | 2.800000000000e-05 |
| PASS | A | fd.mixed_commutes | D_a(D_b f) == D_b(D_a f) at every point, including corners | 21x21x21 | trigonometric field | 0 | 9.592326932761e-14 | 9.592326932761e-14 | - | 1.000000000000e-12 | 3.672500000000e-03 |
| PASS | B | fd.tensor_field.interior | every rank-2 component is differentiated against its own analytic derivative (interior points) | 41x41x41 | g_{mu nu} = (mu+1)sin x + (nu+1)y^2 + (mu nu+1)z | 0 | 1.666458345753e-03 | 1.666458345753e-03 | - | 2.500000000000e-03 | 1.848778000000e-01 |
| PASS | B | fd.tensor_field.boundary | same, at boundary points using one-sided stencils | 41x41x41 | same field | 0 | 1.904572274663e-03 | 1.904572274663e-03 | - | 5.000000000000e-03 | 1.849032000000e-01 |
| PASS | A | fd.tensor_field.time_slot_zero | with no time-derivative field supplied, d_0 g_{mu nu} is exactly 0 | 41x41x41 | n/a | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.849069000000e-01 |
| PASS | A | fd.generic.test_scalar_bit_exact | the derivative operators give bit-identical results for double and a custom scalar with no ordering and no implicit double conversion | 21x21x21 | trigonometric field | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.250600000000e-03 |
| PASS | D | fd.generic.long_double_agrees | long double agrees with double to double precision | 21x21x21 | trigonometric field | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | 1.000000000000e-12 | 1.268700000000e-03 |
| PASS | A | fd.grid_validation | a 3-point axis is accepted for first derivatives and rejected for second derivatives, whose one-sided stencil needs 4 points | n/a | n/a | first accepted, second rejected | first accepted, second rejected | - | - | exact | 1.500000000000e-05 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
