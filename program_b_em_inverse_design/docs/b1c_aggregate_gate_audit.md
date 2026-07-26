# B1C aggregate gate: audit of the six failed assertions

The aggregate B1C gate (`b1c_verify`, CTest `programb1c.gate`) reported

```
23/29 passed
6 failed assertions
```

and was the only failing test in the repository. This document classifies every
one of the six, and records what replaced it.

**No test was deleted or weakened to make the suite pass.** Two assertions were
formally withdrawn as mathematically invalid, and both are replaced by strictly
stronger statements — including, in each case, a positive assertion of what is
actually true where the old one asserted something false. No tolerance was
relaxed. The chain-rule tolerance of 1e-5 was not relaxed; the rebuilt case runs
at 1e-14.

Corrected tally: **53/53 passed, 0 failed** with OTIlib enabled; **45/45 passed,
4 skipped** with OTIlib disabled, the four skips being HyPAD-backend cases only.

---

## Classification summary

| # | test | class |
|---|---|---|
| 1 | `b1c.observer_decomposition_epsilon_matches_plebanski` | 3 — mathematically invalid test |
| 2 | `b1c.magnetoelectric_sign_resolved` | 2 + 3 + 4 — legacy, invalid, and superseded |
| 3 | `b1c.eigenvalue_conditioning_diverges_at_ergosurface` | 3 — mathematically invalid test |
| 4 | `b1c.sensitivity_fixed_center_matches_reference` | 3 — mathematically invalid test (ill-conditioned error measure) |
| 5 | `b1c.sensitivity_fixed_lab_matches_total` | 3 — mathematically invalid test (same) |
| 6 | `b1c.total_derivative_matches_chain_rule` | 3 — valid claim, invalid numerical realisation |

Also migrated, though they were passing, because they routed through the
quarantined decomposition: `b1c.observer_decomposition_mu_matches_plebanski`,
`b1c.boosted_observer_decomposition`, `b1c.statuses_are_distinguishable`, and
the `ergosurface_conditioning.csv` / `corrected_status_map.csv` exports.

**Class 1 (implementation defect) does not appear.** That is itself a result: the
rebuilt velocity cases compare the production HyPAD sensitivities against
independently written analytic closed forms and agree to **2.2e-16**. Nothing in
the implementation was wrong.

---

## 1. `b1c.observer_decomposition_epsilon_matches_plebanski`

| | |
|---|---|
| **Original claim** | "for the coordinate-time observer in a flat lab, eps from the Bel decomposition of Z equals the closed-form Plebański eps" |
| **Computed** | 4.166666666667e-02 |
| **Tolerance** | 1e-12 |
| **Reason for failure** | The two quantities belong to **different observers** and are not equal. Derived, not surmised: for the laboratory observer V = (1,0,0,0) normalised with g_lab = η, ε<sup>ij</sup> = −g<sub>eff</sub><sup>00</sup>(γ<sup>−1</sup>)<sup>ij</sup>; the Alcubierre spatial block is δ<sub>ij</sub> and det g = −1, so g<sub>eff</sub><sup>00</sup> = det(γ)/det(g) = −1 and **ε = 𝟙 exactly**. Plebański normalises by g<sub>00</sub> of g<sub>eff</sub> instead, giving ε<sub>P</sub> = diag(1, 1/(1−v²f²), 1/(1−v²f²)). At v_s = 0.2, f = 1 the difference is 1/(1−0.04) − 1 = **1/24 = 4.16667e-2**, which is the failure to every digit. |
| **Classification** | **3 — mathematically invalid test.** It asserted equality between observer-dependent quantities belonging to two different observers. |
| **Replacement** | `b1c.warp_epsilon_closed_form` (ε<sub>lab</sub> = h exactly, all v_s), `b1c.plebanski_epsilon_closed_form` (ε<sub>P</sub> against its own closed form), `b1c.plebanski_differs_by_observer_normalisation` (records the 1/(1−v_s²f²) factor as a measured quantity), and `b1c.observer_change_reconstruction_over_profile`. |
| **Action taken** | Withdrawn as an equality. Both blocks are now asserted against their **own** closed forms, which together *prove* they cannot be equal; the discrepancy is recorded as a number rather than treated as an error. |
| **Source** | Schuster & Visser (A7a) with ε = −2W_ε; Plebański–Tamm ε<sup>ij</sup> = −√(−g) g<sup>ij</sup>/g<sub>00</sub> as implemented in `metric_medium_map.hpp`. The identity g<sup>ij</sup> − g<sup>0i</sup>g<sup>0j</sup>/g<sup>00</sup> = (γ<sup>−1</sup>)<sup>ij</sup> is what connects them. |

### Companion, passing but migrated

`b1c.observer_decomposition_mu_matches_plebanski` passed at 1e-12 and rests on
the same observer identification. It passes for a real reason:
[μ⁻¹]<sub>lab</sub> = diag(0, 1, 1−v²f², 1−v²f²) while μ<sub>P</sub> = ε<sub>P</sub>,
so (μ<sub>P</sub>)⁻¹ = diag(1, 1−v²f², 1−v²f²) — they coincide **componentwise**.
Retained as `b1c.plebanski_mu_inverse_coincides`, restated as a cross-check
between two independent routes and explicitly **not** as evidence that the two
observers agree, since ε demonstrably does not. It is special to a metric whose
spatial block is δ<sub>ij</sub>.

---

## 2. `b1c.magnetoelectric_sign_resolved`

| | |
|---|---|
| **Original claim** | "the covariant zeta agrees with one of the two Plebański magnetoelectric blocks, fixing the sign convention by derivation rather than by assumption" |
| **Computed** | predicate false; reported `matches zeta = -cross_matrix(w)` |
| **Tolerance** | 1e-12 on `min(‖ζ − ξ_P‖, ‖ζ − ζ_P‖)` |
| **Reason for failure** | Three independent defects. (a) It used the quarantined tetrad decomposition. (b) It compared against Plebański blocks belonging to a different observer, so no agreement was possible — ζ<sub>lab</sub> has magnitude v f while \|w\| = v f/(1−v²f²). (c) The `min(...)` over two candidate signs means the case **asserted no sign at all**, which is precisely the construction this project forbids; the reported string names a match while the predicate reports none. |
| **Classification** | **2 + 3 + 4.** Obsolete (legacy path), mathematically invalid (cross-observer), and superseded by B1C-S456. |
| **Replacement** | `b1c.warp_zeta_closed_form` — ζ<sup>23</sup> = −g<sub>01</sub> = v f, ζ<sup>32</sup> = +g<sub>01</sub>, all other components zero, over 4 velocities × 7 positions, residual 0.0 — and `b1c.warp_magnetoelectric_handedness` / `b1c.warp_magnetoelectric_magnitude`, which assert **one** direction: with the bubble moving along +x, B = +B₀ẑ and E = 0, the induced displacement points along **−y**, with D_y = g<sub>01</sub> exactly. In B1C-S456: `b1c.s5.b17_zeta_magnitude`, `b1c.s5.b17_zeta_sign`, `b1c.s5.handedness_direction`, `b1c.s5.handedness_matches_minkowski_relations`. |
| **Action taken** | Withdrawn. No surviving case anywhere in the suite accepts the better of two candidate signs. |
| **Source** | (A7d) with ζ = +2W_ζ, pinned against (B17); D<sup>a</sup> = ε<sup>ab</sup>E_b − ζ<sup>ab</sup>B_b as derived in `field_reconstruction.hpp`; handedness cross-checked against Minkowski's moving-medium relations. |

---

## 3. `b1c.eigenvalue_conditioning_diverges_at_ergosurface`

| | |
|---|---|
| **Original claim** | "the permittivity conditioning for the coordinate-time observer degrades sharply as v_s f → 1 from outside, while the covariant tensor stays finite" |
| **Computed** | `cond(eps) = 1` just outside the ergosurface versus `1` far away |
| **Tolerance** | `cond_near > 10 × cond_far` |
| **Reason for failure** | The asserted quantity **cannot** diverge. A tensor may scale uniformly towards infinity with its condition number pinned at 1 — and here ε does not scale at all: for the laboratory observer ε<sup>ab</sup> = h<sup>ab</sup> identically, at every point, for every v_s, so cond(ε) ≡ 1. The measured 1-versus-1 is the correct answer to the question asked, and carries no information about the physics. |
| **Classification** | **3 — mathematically invalid test.** |
| **Replacement** | Six diagnostics, all computed and exported to `ergosurface_conditioning.csv` (max and min eigenvalue magnitude, determinant, max component magnitude, magnetoelectric norm, and max\|Z<sup>abcd</sup>\|, for both ε and μ⁻¹). Divergence is asserted **only** for a quantity predicted to diverge, under the same observer and the same conventions: <br>• `b1c.ergosurface.epsilon_is_the_projector` — ε = h to 0.0, so cond(ε) ≡ 1 <br>• `b1c.ergosurface.epsilon_conditioning_does_not_diverge` — the withdrawal, kept visible in the gate output <br>• `b1c.ergosurface.mu_inverse_min_eigenvalue_closed_form` — min\|eig μ⁻¹\| = min(1, \|g<sub>00</sub>\|) to 0.0 <br>• `b1c.ergosurface.permeability_diverges` — min\|eig μ⁻¹\| = 2.46e-5, cond(μ⁻¹) = 4.07e4, \|det μ⁻¹\| = 6.04e-10 at the sampled point nearest the surface, against 1, 1, 1 at x = 5 <br>• `b1c.ergosurface.covariant_tensor_stays_bounded` — max\|Z\| = 0.500006 near versus 0.5 far, and bounded in (0.1, 2) over the whole sweep |
| **Action taken** | **Formally withdrawn.** What actually diverges at the ergosurface, for the laboratory observer, is the **permeability**: the spectrum of μ⁻¹ is {1, −g<sub>00</sub>, −g<sub>00</sub>}, so its transverse pair reaches zero exactly on v_s f = 1 and μ = (μ⁻¹)⁻¹ blows up, while ε and Z do not notice the surface. |
| **Source** | (A7a)/(A7b) with ε = −2W_ε and μ⁻¹ = +2W_μ; the closed forms are those verified in `b1c.warp_epsilon_closed_form` and `b1c.warp_mu_inverse_closed_form`. |

Two notes on scope. max\|Z\| is **bounded but not constant** over the sweep — it
rises to about 1.43 deep inside the ergoregion, where g<sub>eff</sub><sup>11</sup>
= 1 − v²f² reaches −2.98 — so boundedness is what is asserted. And the minimum
eigenvalue magnitude of μ⁻¹ is min(1, \|g<sub>00</sub>\|), **not** \|g<sub>00</sub>\|:
inside the ergoregion \|g<sub>00</sub>\| exceeds 1 and the unit longitudinal
eigenvalue becomes the smaller one. Both of these were over-general first
attempts at the replacement assertion, caught by the gate itself.

---

## 4. `b1c.sensitivity_fixed_center_matches_reference`

| | |
|---|---|
| **Original claim** | "FixedBubbleCenter reproduces dg/dv with the trajectory held fixed" |
| **Computed** | 2.908488468709e-02 |
| **Tolerance** | 1e-6 |
| **Reason for failure** | The error measure, not the physics. It computed `\|a − b\| / max(\|b\|, 1e-9)` against a centred finite difference with h = 1e-6. At the exterior point x = 6 the true dg₀₀/dv is **2.64e-11**, while a centred difference of a quantity of size 1 with h = 1e-6 carries roundoff of about 2e-10 — the reference was noise, an order of magnitude larger than the signal it was meant to represent (measured: 5.55e-11 against a true 2.64e-11). The floored denominator then converted that noise into 2.9e-2. The identical failure value in case 5 is the tell: at t = 0 the two modes coincide exactly and both cases were failing at the same point, on the same noise. |
| **Classification** | **3 — mathematically invalid test.** A relative error against a reference that may be, and at one point is, smaller than its own roundoff. |
| **Replacement** | `b1c.velocity.oti_fixed_center_matches_analytic.g00` and `.g01` — HyPAD against the independently written analytic closed form, **absolute** tolerance 1e-12, measured **2.2e-16**; plus `b1c.velocity.analytic_matches_fd_fixed_center.{g00,g01}` — analytic against centred differences of the reference, absolute, four step sizes, measured 3.4e-14 and 1.5e-14; plus `b1c.velocity.reference_reproduces_metric.{g00,g01}` at 0.0, so the reference is not merely self-consistent. |
| **Action taken** | Rebuilt from `evaluate_fixed_center_reference(v, fixed_xs, t, x, y, z)` in `verification/b1c_velocity_reference.hpp`, which calls **no** production code — no `AlcubierreMetric`, no `shape_function`, no `center_at`, no OTI. `fixed_xs` is an argument and `t` is accepted but never used in the body; that absence *is* the fixed-centre condition, so a finite-difference sweep in v cannot move the bubble. |
| **Source** | g₀₀ = −1 + v²f², g₀₁ = −v f; ∂g₀₀/∂v\|<sub>x_s</sub> = 2v f², ∂g₀₁/∂v\|<sub>x_s</sub> = −f. |

---

## 5. `b1c.sensitivity_fixed_lab_matches_total`

| | |
|---|---|
| **Original claim** | "FixedLaboratoryCoordinates reproduces the TOTAL dg/dv\|_{t,x}, including the dx_s/dv_s = t trajectory term" |
| **Computed** | 2.908488468709e-02 |
| **Tolerance** | 1e-6 |
| **Reason for failure** | Identical to case 4, at the same point and on the same finite-difference noise. |
| **Classification** | **3 — mathematically invalid test.** |
| **Replacement** | `b1c.velocity.oti_fixed_lab_matches_analytic.{g00,g01}` — absolute 1e-12, measured **2.2e-16**; `b1c.velocity.analytic_matches_fd_fixed_lab.{g00,g01}` — 3.8e-11 and 6.9e-11 against 1e-8. |
| **Action taken** | Rebuilt from `evaluate_fixed_lab_reference(v, xs0, t, x, y, z)`, which recomputes x_s(v) = x_s0 + v t **inside every call**, so a finite-difference sweep in v moves the bubble exactly as the trajectory does. |
| **Source** | dg/dv\|<sub>t,x</sub> = ∂g/∂v\|<sub>x_s</sub> + (∂g/∂f)(df/dr_s)(∂r_s/∂v). |

---

## 6. `b1c.total_derivative_matches_chain_rule`

| | |
|---|---|
| **Original claim** | "dg/dv\|_{t,x} − dg/dv\|_{r_s} = −t (dg/dx), the hand-derived chain-rule term with dx_s/dv_s = t" |
| **Computed** | 3.353573049043e-03 |
| **Tolerance** | 1e-5 |
| **Reason for failure** | **The claim is correct**; its numerical realisation was not. Same ill-conditioned relative measure — at t = 2, x = 6 the total derivative is 2.71e-8, so the floored denominator turned a ~9e-11 finite-difference residual (from h_x = 1e-6 on a quantity of size 1) into 3.35e-3. |
| **Classification** | **3 — valid claim, invalid numerical realisation.** |
| **Replacement** | Two cases, per component, replacing one. `b1c.velocity.trajectory_term_is_chain_rule.{g00,g01}` — the identity between analytic expressions, absolute tolerance **1e-14**, measured 4.5e-17 and 1.1e-16. `b1c.velocity.trajectory_equals_minus_t_dgdx.{g00,g01}` — the −t(∂g/∂x) form **measured** rather than assumed, with ∂g/∂x differenced from the *fixed-centre* reference so that x_s genuinely cannot drift while x moves; 3.2e-11 and 9.1e-11 against 1e-8. |
| **Action taken** | Rebuilt. The tolerance was **not** relaxed; it was tightened from 1e-5 to 1e-14 by removing the finite difference from the identity being tested. The trajectory contribution is assembled component by component as (∂g/∂f)(df/dr_s)(∂r_s/∂v) rather than through a global assumption. |
| **Source** | r_s = √((x−x_s)² + y² + z²) with x_s = x_s0 + v t gives ∂r_s/∂v = −((x−x_s)/r_s)·t. Then ∂g₀₀/∂f = 2v²f and ∂g₀₁/∂f = −v. At r_s = 0 the trajectory term is **exactly** zero, not regularised: it always carries a factor df/dr_s, which vanishes there because f is even, while \|(x−x_s)/r_s\| ≤ 1 stays bounded. |

### On the identity ∂g/∂x_s = −∂g/∂x

The correction task warns against assuming this globally. It does hold for this
metric, because g depends on x and x_s solely through (x − x_s) while the
explicit v-dependence is untouched — but it is now **measured** per component,
at fixed x_s and fixed v, rather than taken for granted. The two derivatives are
kept distinct throughout: `analytic_dgdv_fixed_center` is
∂g/∂v\|<sub>x_s</sub> and `analytic_dgdv_fixed_lab` is
dg/dv\|<sub>t,x, x_s = x_s0 + v t</sub>, and their difference is the trajectory
term.

### Coverage

14 points: bubble centre, interior, wall (x_s + 3.0), outer wall (x_s + 3.385),
exterior (x_s + 6.0), and off-axis (y, z ≠ 0), at t = 0, t = 2 and t = −1.5.
Four step sizes, 1e-3 to 1e-6, best taken per point and component. Errors are
reported **per metric component** (g₀₀ and g₀₁, the only two that depend on v)
with the worst point named in each case; no case reduces to a single maximum
over components and points. `b1c.velocity.constant_components_carry_no_channel`
additionally asserts that the other seven components carry **exactly** zero
velocity sensitivity in both modes, tolerance 0.0.

---

## Legacy quarantine

`legacy_unverified.hpp` is no longer `#include`d by any translation unit. It
survives as a header for historical reference and its contents are unreachable
from every gate. The four call sites in the aggregate gate were migrated, not
deleted:

| original | property it verified | replacement |
|---|---|---|
| `legacy_unverified::decompose` in the Plebański cases | ε, μ⁻¹ of the warp analogue | `decompose_A7` in `b1c.warp_*_closed_form` |
| `b1c.boosted_observer_decomposition` | blocks are observer dependent | `b1c.boosted_observer_blocks_differ`, via `decompose_A7` and `laboratory_boosted_observer` |
| `legacy_unverified::decompose` in the conditioning sweep | ε conditioning near the ergosurface | the six-diagnostic battery above, via `decompose_A7` |
| `legacy_unverified::analyse` in the status map and `b1c.statuses_are_distinguishable` | the five-way status classification | `analyse_metric_medium()` in `metric_medium_status.hpp`, which takes the laboratory observer and the effective-geometry observer as **separately typed** arguments so the two questions cannot be conflated |

---

## Consolidation

The aggregate gate no longer duplicates its sub-gates. B1C-S23 and B1C-S456 are
separate binaries and separate CTest tests, and `programb1c.gate` depends on both
through a **single** `set_tests_properties(... DEPENDS ...)` assignment —
successive calls replace the property rather than appending to it, so two of the
three dependencies were previously being dropped silently.

What remains in the aggregate binary is only what is cross-cutting: construction
of Z from a pair of metrics (vacuum, conformal invariance, tensoriality — the
sub-gates take Z as given), the warp analogue's blocks against closed forms over
the profile, the ergoregion correction and its divergence diagnostics, the two
velocity-sensitivity modes, objective invariance, the five-way status
classification, and HyPAD channels through Z. It writes
`consolidated_status.md` pointing at the sub-gate records.
