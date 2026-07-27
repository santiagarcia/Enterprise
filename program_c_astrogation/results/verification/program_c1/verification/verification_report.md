# program_c1 verification report

**C1-RELATIVISTIC-ASTROGATION: ship worldlines, clocks, moving-target interception, aberration, Doppler and the retarded apparent sky**

| | |
|---|---|
| Generated (UTC) | 2026-07-27T22:29:32Z |
| Git commit | `893dcfd6c463793a0201d824d20a9d72b7a3bf52` **(working tree dirty)** |
| Compiler | MSVC 19.44.35223.0 |
| Build type | Release |
| C++ standard | C++17 |
| Platform | Windows-AMD64 |
| Cases | 22 |
| Passed | 22 |
| Failed | 0 |
| Skipped | 0 |
| **Gate** | **PASS** |

## All cases

| Result | Level | Test | Reference | Grid | Parameters | Expected | Computed | Abs error | Rel error | Tolerance | Runtime (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PASS | A | c1.constant_velocity_worldline | x = v t for a constant coordinate velocity of 0.5 c over 4 yr | analytic | beta = 0.5, t = 4 yr | 6.132027875711e-01 | 6.132027875711e-01 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-15 | 5.000000000000e-07 |
| PASS | A | c1.proper_time_integration | Simpson integration of dtau = dt sqrt(1-beta^2) against the exact tau = t/gamma at beta = 0.8 | 1000 intervals | beta = 0.8, t = 10 yr | 6.000000000000e+00 | 6.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-14 | 2.900000000000e-06 |
| PASS | A | c1.proper_acceleration_analytic_position | x(tau) = (c^2/a)[cosh(a tau/c) - 1] at 1 g for tau = 2 yr | analytic | a = 1 g, tau = 2 yr | 8.923558363143e-01 | 8.923558363143e-01 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-15 | 5.200000000000e-06 |
| PASS | A | c1.proper_acceleration_analytic_coordinate_time | c t(tau) = (c^2/a) sinh(a tau/c) at 1 g for tau = 2 yr | analytic | a = 1 g, tau = 2 yr | 3.756289449781e+00 | 3.756289449781e+00 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-15 | 9.300000000000e-06 |
| PASS | A | c1.proper_acceleration_analytic_velocity_and_gamma | beta = tanh(a tau/c), gamma = cosh(a tau/c), and gamma = 1/sqrt(1-beta^2) all agree | analytic | a = 1 g, tau = 2 yr | 0 | 2.664535259100e-15 | 2.664535259100e-15 | - | 1.000000000000e-14 | 1.800000000000e-06 |
| PASS | A | c1.symmetric_journey_endpoint | a symmetric 1 g journey arrives at the destination distance and at rest | analytic reflection | a = 1 g, D = 1.3019705311704333 pc | 0 | 2.220446049250e-16 | 2.220446049250e-16 | - | 1.000000000000e-13 | 1.110000000000e-05 |
| PASS | A | c1.moving_destination_interception | the interception solve converges and places the ship on the destination's worldline, not at its catalogue position | n/a | n/a | both solves converge below 1e-9 pc | Proxima residual 7.4851236320228054e-13 pc in 39 iterations; Barnard residual 8.3000273320976703e-13 pc in 38 iterations | - | - | exact | 3.470000000000e-04 |
| PASS | A | c1.stale_position_targeting_error | aiming at the stale catalogue position misses, and misses more for the faster-moving star | n/a | n/a | Barnard's error exceeds Proxima's, and both are nonzero | Proxima 79.747097422293933 au, Barnard 335.42125009883409 au | - | - | exact | 3.000000000000e-07 |
| PASS | B | c1.interception_uncertainty | the interception position uncertainty exceeds the catalogue position uncertainty, because velocity errors integrate over the flight | n/a | n/a | arrival sigma at least the epoch sigma | epoch 8.4540921002773233e-05 pc, arrival 8.4613680339280109e-05 pc over 5.8713284105906496 yr | - | - | exact | 2.300000000000e-06 |
| PASS | A | c1.lorentz_norm_preservation | the boost preserves the Minkowski norm of a null and a timelike four-vector, and the four-velocity is unit timelike | signature (-,+,+,+) | beta = (0.3,-0.4,0.2) | 0 | 1.415534356397e-15 | 1.415534356397e-15 | - | 1.000000000000e-14 | 5.000000000000e-07 |
| PASS | A | c1.aberration_analytic | the four-vector boost reproduces cos t' = (cos t + b)/(1 + b cos t) for the source direction, over 37 angles | 37 angles | beta = 0.6 | 0 | 2.220446049250e-16 | 2.220446049250e-16 | - | 1.000000000000e-14 | 2.800000000000e-06 |
| PASS | A | c1.aberration_crowds_toward_motion | apparent source positions crowd toward the direction of motion at every angle | n/a | n/a | cos of the apparent angle never below the geometric one | holds at all 37 angles | - | - | exact | 5.900000000000e-06 |
| PASS | A | c1.doppler_approaching | nu'/nu = sqrt((1+b)/(1-b)) for a source directly ahead | analytic | beta = 0.6 | 2.000000000000e+00 | 2.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-15 | 2.600000000000e-06 |
| PASS | A | c1.doppler_receding | nu'/nu = sqrt((1-b)/(1+b)) for a source directly behind | analytic | beta = 0.6 | 5.000000000000e-01 | 5.000000000000e-01 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-15 | 3.100000000000e-06 |
| PASS | A | c1.zero_velocity_sky_recovery_doppler | at zero velocity every Doppler factor is exactly 1 | 21 directions | beta = 0 | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.500000000000e-06 |
| PASS | A | c1.zero_velocity_sky_recovery_directions | at zero velocity every apparent direction equals the geometric direction to within one renormalisation rounding | 21 directions | beta = 0 | 0 | 2.237726045656e-16 | 2.237726045656e-16 | - | 4.000000000000e-16 | 4.000000000000e-06 |
| PASS | A | c1.retarded_source_state | the emission time satisfies \|x_obs - x_star(t_e)\| = c (t_obs - t_e) | 5 iterations | Proxima from Earth | 4.246038798014e+00 | 4.246038798014e+00 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-12 | 6.200000000000e-06 |
| PASS | A | c1.earth_to_ship_signal_delay | the Earth clock the ship can SEE lags coordinate time by \|x\|/c | single event | one third of the way to Proxima | 2.025584513359e+03 | 2.025584513359e+03 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-14 | 2.500000000000e-06 |
| PASS | A | c1.ship_to_earth_signal_delay | a message sent by the ship reaches Earth at t + \|x\|/c, and the round trip is 2\|x\|/c | single event | one third of the way to Proxima | 2.028415486641e+03 | 2.028415486641e+03 | 0.000000000000e+00 | 0.000000000000e+00 | 1.000000000000e-14 | 1.800000000000e-06 |
| PASS | A | c1.apparent_sky_source_count | the apparent-sky table has exactly one row per placeable catalogue source, in both the rest and moving observer frames | n/a | n/a | 315 rows in each frame | 315 rest, 315 moving, from 315 snapshot objects | - | - | exact | 3.087000000000e-04 |
| PASS | B | c1.uncertainty_angle_propagation | every apparent-sky row carries a finite, non-negative angular uncertainty, and the small-angle form remains defensible | n/a | n/a | all finite and below 0.1 rad | worst 0.034074858474507121 rad over 315 sources | - | - | exact | 3.200000000000e-06 |
| PASS | A | c1.deterministic_worldline | rebuilding the worldline from the same inputs reproduces it bit for bit | 241 samples | 1 g symmetric to Proxima | 0 | 0.000000000000e+00 | 0.000000000000e+00 | - | exact | 1.855000000000e-04 |

## Verification levels

| Level | Meaning |
|---|---|
| A | exact algebraic identity |
| B | analytic solution |
| C | manufactured solution |
| D | independent numerical method |
| E | grid-convergence study |
| F | cross-code validation |
