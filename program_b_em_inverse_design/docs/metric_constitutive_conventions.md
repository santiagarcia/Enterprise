# Metric ↔ constitutive conventions

Every convention Program B depends on, fixed in one place. Conventions in this
area differ between sources more than in most of physics, and a sign or a factor
of √(−g) silently absorbed is the most likely way for this module to be wrong
while appearing to work.

For each equation: **convention · signature · determinant factor · index
positions · units · source · domain of validity.**

---

## 0. Global choices

| | |
|---|---|
| Metric signature | **(−, +, +, +)** — same as Program A |
| Coordinate index map | 0 = t, 1 = x, 2 = y, 3 = z |
| Latin indices i, j, k | spatial, 1–3 |
| Units | **Heaviside–Lorentz with c = 1**, so vacuum is ε = μ = 1 |
| ε₀, μ₀ | both 1 by the unit choice; never appear |
| Frequency dependence | none — non-dispersive, lossless, real-valued |
| Scalar type | templated; all formulas hold for `double` and `OTIScalar` alike |

**Why Heaviside–Lorentz and c = 1.** It makes vacuum the identity, which turns
`b1.minkowski_is_vacuum` into an exact test rather than a
units-conversion test, and it matches Program A's geometrized units so a metric
can pass between the programmes with no rescaling. SI conversion is one overall
factor per block and is **not applied anywhere** in this module.

## 1. Constitutive relations

    D^i = ε^{ij} E_j + ξ^{ij}  H_j
    B^i = ζ^{ij} E_j + μ^{ij} H_j

- **Index positions:** ε, μ, ξ, ζ carry two upper spatial indices; E and H carry
  lower spatial indices. In Cartesian coordinates with a flat background 3-metric
  this distinction does not change any component, but it is stated because it
  *does* matter the moment curvilinear coordinates appear.
- **Units:** ε and μ dimensionless; ξ and ζ dimensionless (they would carry
  1/c in SI).
- **Domain:** linear, local, non-dispersive, lossless media only.
- **Source:** standard bianisotropic form; matches the layout in Hehl & Obukhov.

## 2. Plebański–Tamm forward map

    ε^{ij} = μ^{ij} = − √(−g) g^{ij} / g_00                    (2.1)
    w_i    = − g_{0i} / g_00                                    (2.2)
    D = ε E + w × H,     B = μ H − w × E                        (2.3)

equivalently ξ = [w]<sub>×</sub> and ζ = −[w]<sub>×</sub>, where
[w]<sub>×</sub>H ≡ w × H.

- **Convention:** Plebański–Tamm.
- **Signature:** (−,+,+,+). With (+,−,−,−) the signs of (2.1) and (2.2) change
  and vacuum would come out as −I; **do not mix**.
- **Determinant factor:** √(−g) with the **4-dimensional** determinant, not the
  spatial 3-determinant. Getting this wrong changes the conformal weight and
  breaks §5.
- **Index positions:** g<sup>ij</sup> in (2.1) is the **spatial block of the
  inverse 4-metric**, `g_inv(i+1, j+1)` — it is *not* the inverse of the spatial
  block g<sub>ij</sub>. These differ whenever g<sub>0i</sub> ≠ 0, i.e. for every
  warp metric. In code: `gi.inverse(i+1, j+1)`.
- **Source:** Plebański, *Phys. Rev.* **118**, 1396 (1960) — verified;
  reproduced in Leonhardt & Philbin (2006) — unverified.
- **Domain of validity:** **g<sub>00</sub> < 0** and det g < 0. Both are checked
  and reported by `metric_to_medium()`; see §7.

### 2.1 The sign of w — this project's convention

Sources differ on the sign of w and on whether (2.3) is written with `+w×H` or
`−w×H`. **This project fixes the convention above and does not claim it matches
any particular source.**

It is pinned by three checks that are insensitive to external convention:

| Check | Pins |
|---|---|
| `b1.minkowski_is_vacuum` | overall normalisation and the sign of ε |
| `b1.conformal_invariance` | the power of √(−g) |
| `b1.roundtrip_conformal` | self-consistency of forward and inverse |

None of these is sensitive to the sign of w. **Consequence: the physical
handedness of the magnetoelectric drag is not established by B1**, and no claim
depending on it may be made until a propagation test exists (B2). The magnitude
of w is well determined; its sign is a convention until then.

## 3. Cross-product matrix

    [w]_{ik} = ε_{ijk} w_j    so that   [w]_× H = w × H

explicitly

        ⎡  0   −w₃   w₂ ⎤
    [w]=⎢  w₃   0   −w₁ ⎥
        ⎣ −w₂   w₁   0  ⎦

- **Convention:** ε<sub>123</sub> = +1, right-handed Cartesian.
- **Note:** this is the *ordinary* Levi-Civita symbol, not a tensor density —
  the background here is flat Cartesian, so no √γ appears.

## 4. Inverse map (w = 0 only)

Gauge: **g<sub>00</sub> = −1** (`ConformalGauge::UnitLapse`).

    h = det(ε) · ε^{-1}                                          (4.1)
    g = diag-block(−1, h)

**Derivation.** With g<sub>00</sub> = −1 and g<sub>0i</sub> = 0, (2.1) gives
ε = √(det h)·h<sup>−1</sup>. Taking determinants,
det ε = (det h)<sup>3/2</sup>/det h = (det h)<sup>1/2</sup>, so
√(det h) = det ε, and substituting back gives (4.1).

- **Domain of validity:** w = 0, ε = μ symmetric positive-definite, det ε > 0.
  All four are checked; violations return a status, never a fitted answer.
- **Not implemented:** w ≠ 0. This *is* metric-representable (it corresponds to
  g<sub>0i</sub> ≠ 0) but the closed form above does not cover it, and the code
  returns `MagnetoelectricUnsupported` rather than approximating. Every warp
  metric falls in this case — see the limitations list.

### Sanity values

| ε | h from (4.1) | g | null speed |
|---|---|---|---|
| I | I | diag(−1,1,1,1) | 1 |
| n²I | n⁴I | diag(−1,n⁴,n⁴,n⁴) | 1/n² |

The second is *not* 1/n, and that is correct: ε = μ = n² gives
1/√(εμ) = 1/n². Compare §6.

## 5. Conformal ambiguity

Under g → Ω²g with Ω a positive function:

    √(−g) → Ω⁴ √(−g),   g^{ij} → Ω^{−2} g^{ij},   g_00 → Ω² g_00

so ε → Ω<sup>4−2−2</sup>ε = ε, and w is unchanged as a ratio of two Ω² factors.

**The medium determines the metric only up to a conformal factor.** The API
represents this explicitly rather than hiding it:

- `medium_to_metric()` returns a `ConformalGauge` naming the representative;
- `compare_up_to_conformal()` returns the recovered Ω² *and* the residual, so a
  caller sees what was divided out;
- round-trip tests assert conformal equality, and a companion case
  (`b1.roundtrip_factor_is_nontrivial`) confirms Ω² ≠ 1 so the conformal test is
  not a disguised equality test.

Conformally invariant observables — null cones, ray paths, the dispersion
relation — are the only things a round trip may be asked to preserve.

## 6. Gordon optical metric

    ḡ_{μν} = g_{μν} + (1 − 1/n²) u_μ u_ν                        (6.1)

- **Signature:** (−,+,+,+); u<sub>μ</sub>u<sup>μ</sup> = −1.
- **Static isotropic medium in flat space:** u<sub>μ</sub> = (−1,0,0,0), giving
  ḡ = diag(−1/n², 1, 1, 1) and null coordinate speed **1/n**.
- **Domain:** ray/eikonal level only. Applies to media that are *not* impedance
  matched, where no Plebański metric exists.
- **Source:** Gordon (1923) — unverified.
- **Do not conflate with §2.** Different claim, different domain, different
  speed. See [effective_metric_theory.md](effective_metric_theory.md) §4.

## 7. Validity conditions, and what happens when they fail

| Condition | Where checked | On failure |
|---|---|---|
| g<sub>00</sub> < 0 | `metric_to_medium` | `ok = false` + reason string |
| det g < 0 | `metric_to_medium` | `ok = false` + reason string |
| g invertible | `metric_to_medium` | `ok = false` + reason string |
| ε = μ | `classify` | `NotMetricRepresentable` |
| ε, μ symmetric | `classify` | `NotMetricRepresentable` |
| ξ = −ζ antisymmetric | `classify` | `NotMetricRepresentable` |
| w = 0 (inverse only) | `medium_to_metric` | `MagnetoelectricUnsupported` |
| det ε > 0 | `medium_to_metric` | `Degenerate` |

**No path returns a best-fit metric with a warning.** Either the conditions hold
and the result is exact, or a status explains why there is no result.

## 8. Null coordinate speed

For a diagonal metric, g<sub>00</sub>dt² + g<sub>ii</sub>dx<sup>i</sup>² = 0 gives

    |dx^i/dt| = √( −g_00 / g_ii )                                (8.1)

Note the ordering: **−g<sub>00</sub>/g<sub>ii</sub>**, not its reciprocal. An
inverted version of (8.1) was written in an early draft of the Gordon check and
produced n instead of 1/n; the gate caught it.

## 9. Cross-programme consistency

| | Program A | Program B |
|---|---|---|
| Signature | (−,+,+,+) | (−,+,+,+) |
| Units | geometrized, G = c = 1 | Heaviside–Lorentz, c = 1 |
| Index map | 0=t, 1=x, 2=y, 3=z | same |
| `Mat4` layout | row-major, `m(row,col)` | same type, reused |
| Failure mode at v<sub>s</sub>f > 1 | no coordinate observers | no equivalent medium |

The last row is the same condition reached independently from two directions.

---

# Part II — the covariant formulation (B1C)

Sections 0–9 above describe the 3+1 Plebański–Tamm route of B1. Sections 10–13
describe the covariant route that superseded it, following Schuster & Visser,
*Phys. Rev. D* **96**, 124019 (2017). Where the two disagree the covariant one
is authoritative; see [b1c_source_equations.md](b1c_source_equations.md).

## 10. The covariant constitutive tensor and its observer blocks

    Z^{abcd} = ½ √(det g_eff / det g_lab)
               [ g_eff^{ac} g_eff^{bd} − g_eff^{ad} g_eff^{bc} ]        (6)

    H^{ab} = Z^{abcd} F_{cd}                                            (4)

`H^{ab}` is the excitation tensor, written `G^{ab}` in the paper. Symmetries:
antisymmetric in `(ab)` and in `(cd)` (12), and pair-symmetric
`Z^{abcd} = Z^{cdab}` (13) — 21 degrees of freedom.

Blocks relative to a **laboratory** observer `V^a`, `g_lab(V,V) = −1`:

    W_ε^{bd}      = V_a V_c Z^{abcd}          ε^{bd}    = −2 W_ε^{bd}
    W_μ^{bd}      = V_a V_c (*Z*)^{abcd}      [μ⁻¹]^{bd} = +2 W_μ^{bd}
    [W_ζ^T]^{bd}  = V_a V_c (*Z)^{abcd}       [ζ†]^{bd} = +2 [W_ζ^T]^{bd}
    W_ζ^{bd}      = V_a V_c (Z*)^{abcd}       ζ^{bd}    = +2 W_ζ^{bd}

The −2/+2/+2/+2 pattern is derived and measured, not guessed. All four blocks
are spatial, `V_a X^{ab} = X^{ab} V_b = 0`, structurally.

## 11. Levi-Civita tensor and Hodge duals

    ε_{0123} = + √(−det g_lab)          ε^{0123} = − 1/√(−det g_lab)

The contravariant sign is **forced** once the covariant one and the signature
are fixed: raising four indices in (−,+,+,+) contributes
`(−1)(+1)(+1)(+1) = −1`. It is not a free choice.

The Levi-Civita **tensor** is built from `g_lab` at every point. The bare
permutation **symbol** appears only inside the two constructors; hard-coding a
Minkowski symbol into the dual operations would give wrong answers on any
non-Cartesian or curved laboratory background — the failure the
`b1c.coordinate_change.block_excitation` case is designed to catch.

    (*Z)^{abcd} = ½ ε^{ab}_{~~ef} Z^{efcd}          left dual
    ( Z*)^{abcd} = ½ Z^{abef} ε_{ef}^{~~cd}         right dual
    (*Z*)^{abcd} = left then right

In four-dimensional Lorentzian signature `**F = −F`, from
`ε_{abcd} ε^{cdef} = −2(δ_a^e δ_b^f − δ_a^f δ_b^e)`.

## 12. The 3+1 fields, and how they relate to the blocks

For a normalised laboratory observer `V^a`, applied identically to the field
strength and to the excitation:

    x^a =  X^{ab} V_b            y^a = −(*X)^{ab} V_b
    X^{ab} = V^a x^b − V^b x^a + ε^{abcd} V_c y_d

`(E,B)` from `F`, `(D,H)` from the excitation. Vacuum fixes the naming: with
`g_eff = g_lab` the excitation is `F^{ab}`, so `D = E` and `H = B`.

    D^a = ε^{ab} E_b − ζ^{ab} B_b
    H^a = [ζ†]^{ab} E_b + [μ⁻¹]^{ab} B_b

**The minus sign on ζ is this project's convention and is derived**, not
inherited: it differs from the paper's three-dimensional equation (15), which is
written in the index convention of its equation (20) rather than that of its own
covariant (A6)/(A7). One consistent transposition of ζ explains both this sign
and the fact that (A11)'s magnetoelectric slot takes ζ† rather than ζ. The
handedness is checked against a source outside the paper — Minkowski's
moving-medium relations — and asserted in one direction, not compared both ways.

Reconstruction of `Z` from the blocks is equation (A11), implemented in
`z_reconstruction.hpp`. It requires pair symmetry (13), under which
`ζ† = ζ^T`; that identity is measured at every call site, never assumed.

## 13. g_lab and g_eff are different objects

| | `g_lab` | `g_eff` |
|---|---|---|
| what it is | the physical laboratory background | the light cones the medium mimics |
| role | raises and lowers **every** index | an output of the design problem |
| normalises observers | yes | only diagnostic ones |
| builds the Levi-Civita tensor | yes | never |
| value for the analogue experiment | η | whatever is being designed |

Enforced in the type system: `LaboratoryMetric<Scalar>` and
`EffectiveMetric<Scalar>` are distinct types with explicit constructors, and
`make_laboratory_observer()` and `make_effective_geometry_observer()` return
distinct types, so a diagnostic observer of the effective geometry cannot be
fed into the constitutive decomposition. The substitution is a silent physics
error otherwise — every residual stays finite and the answer is wrong — and it
is the family of mistake behind B1's retracted conclusion.

Non-Minkowski laboratory backgrounds do appear, in the general tensor tests and
in the coordinate-change case. They are labelled as such and are **not** a claim
about a curved laboratory.
