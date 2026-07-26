# B4C — propagation interpretation and differentiability correction

Three corrections to B4. Gate `B4C-PROPAGATION-CORRECTION`: **18/18 passed,
0 failed, 0 skipped** with OTIlib enabled; 16 passed and 2 skipped with it
disabled, the two being HyPAD comparisons.

---

## 1. The metricity claim — withdrawn and replaced

**Withdrawn:**

    metricity iff no birefringence

**Replacement, and the only statement the tests support:**

> Metric-induced constitutive tensors are non-birefringent. Generic reciprocal
> constitutive tensors need not be. The converse requires additional assumptions
> and is subject to conformal/dilaton freedom and constitutive pieces that do not
> affect the Fresnel polynomial.

What B4 actually established was: metric-induced tensors give a repeated
quadratic cone; **one selected** generic reciprocal anisotropic tensor
birefringes; reciprocity alone does not imply metricity. None of that is the
converse.

**The Fresnel polynomial does not uniquely reconstruct `Z^{abcd}`.** Two
independent obstructions, both now measured rather than argued.

---

## 2. Propagation-invisible constitutive pieces

### The axion term

    Z_axion^{abcd} = alpha eps^{abcd}

with `eps^{abcd}` the Levi-Civita **tensor** of `g_lab` in this project's own
convention `eps^{0123} = -1/sqrt(-det g_lab)`. It is totally antisymmetric, hence
antisymmetric in each index pair and symmetric under pair exchange — a legitimate
reciprocal constitutive tensor.

**Why it is invisible.** The Fresnel construction goes through
`W^{bd} = Z^{abcd}k_a k_c`, and for the axion piece

    W_axion^{bd} = alpha eps^{abcd} k_a k_c = 0

**identically** — a totally antisymmetric symbol contracted with the symmetric
`k_a k_c`. Not approximately, not up to a factor: the axion contributes nothing to
`W`, so nothing to its adjugate, so nothing to the cone.

Added to a metric-induced tensor with `alpha = 0.37`:

| quantity | measured |
|---|---|
| change in Fresnel coefficients | **0.0** (bit-identical) |
| change in `Z` | 3.70e-1 |
| change in the `zeta` block | 7.40e-1 = exactly `2 alpha` |
| change in `eps` and `mu^-1` blocks | **0.0** |
| `alpha` recovered by projection | 3.70e-1 |

The axion is exactly the **trace part of a symmetric zeta**, since
`zeta^{bd} = 2 V_a V_c (Z*)^{abcd} = 2 alpha h^{bd}` for the axion piece alone —
and the measured `2 alpha h` confirms it to machine precision.

### Conformal / dilaton freedom

Two effective metrics differing by `Omega² = 3.3` give `max|g2 - g1| = 2.3` and
`max|Z2 - Z1| = 1.11e-16`. `Z` is conformally invariant, so the cone cannot
distinguish them: `g_eff` is determined only up to a conformal factor. This is a
second, independent reason the cone underdetermines the constitutive tensor.

**Together:** same Fresnel cone ⇏ identical constitutive tensor.

---

## 3. The classification, piece by piece

B4's single `metric_induced` boolean is gone. It could not have distinguished the
metric tensor from the metric-plus-axion tensor at all — they share a cone and
differ in constitutive content. The replacement reports each measurable property
under its own name:

| flag | what it measures |
|---|---|
| `reciprocal` | `Z^{abcd} = Z^{cdab}` |
| `pair_antisymmetric` | `Z^{(ab)cd} = Z^{ab(cd)} = 0` |
| `contains_axion` | the projection `alpha` is nonzero |
| `contains_skewon` | the pair-antisymmetric part is nonzero |
| `fresnel_single_cone` | the quartic is a doubled quadric |
| `principal_part_metric_compatible` | **a necessary condition read off the cone, not metricity** |

Measured, three tensors, three distinct patterns:

| tensor | reciprocal | axion | skewon | single cone |
|---|---|---|---|---|
| metric-induced | yes | no | no | **yes** |
| metric + axion | yes | **yes** | no | **yes** |
| anisotropic `diag(2, 3.5, 5)` | yes | no | no | **no** |

For the reciprocal tensors this project works with, pair symmetry kills the
skewon piece and the count is `21 = 20 (principal) + 1 (axion)`.

---

## 4. The root finder — replaced

B4 bracketed **sign changes**. A root of even multiplicity does not change the
sign, so the doubled root of the metric-induced case was *never actually
recovered*; the factorization test passed only because it compared `lambda`
against a separately computed cone. Complex roots were invisible too, so
evanescent branches were silently absent rather than reported.

B4C uses **Durand–Kerner** simultaneous iteration with Newton polish, clustering
for multiplicity, and per-root residual and conditioning reporting. Every root
carries real and imaginary parts, a residual, a multiplicity estimate, and a
warning when the multiplicity makes it ill-conditioned — the attainable accuracy
for multiplicity `m` is about `eps^(1/m)`, so about 1e-8 for a double root.

| case | result |
|---|---|
| four distinct real roots | 4 real, 4 distinct |
| two real + one complex-conjugate pair | 2 real, 1 complex pair, residual 1.05e-17 |
| two exact double roots | 4 real, **2 distinct**, multiplicity 2 |
| pair split by 1e-7 | 4 real, 3 distinct (resolved) |
| worst root error over all cases | 9.44e-8 |
| **metric-induced quartic** | 2 distinct, multiplicity 2, positive-branch separation **0.0** |
| anisotropic quartic | 4 distinct, relative split 2.43e-1 |
| scaling by 1e6 | roots unchanged, **0.0** |
| coefficients rebuilt from roots | **0.0** |

The metric double root is now recovered **by the root finder itself**, not
inferred from a separate comparison.

---

## 5. The layered solver — scalar-generic

`LayerMode`, `Layer`, `SMatrix`, `StackResult`, `layer_modes`,
`interface_smatrix`, `propagation_smatrix`, `star` and `solve_stack` are all
templated on `Scalar`. The same code path runs at `double` and at the HyPAD
scalar, so a physical parameter is differentiated **through the actual assembly**
— interfaces, propagation, Redheffer cascade — and not through a surrogate.

No branch is taken on a `Scalar` value in the production arithmetic. The
comparisons live in reporting and branch-selection helpers, deliberately outside
the differentiated path, and `R`, `T`, `A` and `arg t` are reduced to `double` as
reporting quantities.

Verified: the templated solver at `double` gives `R = 0.0` and
`phase = omega eps d` for an impedance-matched slab.

---

## 6. Differentiating the phase without `atan2`

`arg t` needs `atan2`, which is branch sensitive. Its **derivative** does not:

    d(phase) = Im( dt / t )      for t != 0

    d(phase)/dq = ( t_re dt_im/dq - t_im dt_re/dq ) / |t|²

using only the real and imaginary channels of the transmission amplitude. No
inverse trigonometry appears anywhere in the differentiated path.

Against the matched-slab closed forms `d(phase)/d(eps) = omega d` and
`d(phase)/d(thickness) = omega eps`:

| route | error |
|---|---|
| HyPAD through the actual stack | `2.2100` vs `2.2100`, `3.3800` vs `3.3800` — error below 1e-12 |
| finite differences, 3 step sizes | 6.44e-8 |
| separate versus batched seeding | **0.0** |

Step-independence across `h` in {1e-4, 1e-5, 1e-6} for permeability, the
magnetoelectric amplitude, and a Lorentz oscillator's strength, resonance and
damping: **1.54e-9**.

### Exclusions, stated rather than glossed

No claim is made to differentiate through branch sorting, mode switching, poles,
zeros of transmission, or degeneracies. `phase_sensitivity()` **refuses to report
a value** when `|t|` falls below 1e-8 and says why, rather than returning a large
number. Root multiplicity carries an explicit conditioning warning.

---

## 7. Restricted-versus-exact propagation sensitivities

Sensitivities of phase, transmission, reflection and group delay to `v_s`, the
magnetoelectric scale, and the wall steepness `sigma`, for the exact profile and
the `zeta = 0` restriction, are exported to `b4c_sensitivities.csv`.

`d(phase)/d(magnetoelectric scale)` is **1.36 (exact)** and **1.87
(restricted)** — both of order unity. So B4's 1.6-radian mismatch is not an
accident of one parameter value: the phase moves by radians per unit change in
the magnetoelectric scale, and the restricted family sits on a comparably steep
part of the same landscape.

---

## 8. No regression

The lossy-slab energy balance after templating: **6.87e-9**, unchanged from B4.
Two identical solves are bit-identical.

---

## 9. What B4C does not change

The scope statement of B4 stands: a passing result demonstrates that the B3
constitutive model reproduces selected propagation observables of the target
effective metric. It is not a fabricated metamaterial, unit cell, antenna, coil
or experimental device. Causality remains global; passivity remains verified over
the tested band only.
