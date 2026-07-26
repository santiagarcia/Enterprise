# Effective-metric electrodynamics: theory review for B1

**Scope.** What is known about the correspondence between a spacetime metric
g<sub>μν</sub>, the linear constitutive tensor χ<sup>μνρσ</sup>, and the 3+1
material tensors ε, μ, ξ, ζ — and, more importantly, where that correspondence
stops holding.

> **On references.** The bibliography in §11 is recorded from memory. Author,
> title, venue and year are stated because they are load-bearing; volume and
> page numbers are given only where confident. **Plebański (1960) has been
> verified against the APS record.** The rest carry a confidence marking and
> must be checked against the originals before appearing in a paper. No DOIs are
> given, because they have not been verified. Nothing here is cited for a claim
> the project has not also verified numerically.

---

## 1. Premetric electrodynamics: the general setting

Maxwell's equations split into two parts of very different character:

    dF = 0            (homogeneous; F = dA, topological)
    dH = J            (inhomogeneous; charge conservation)

Neither mentions a metric. All metric information sits in the **constitutive
relation** linking the excitation H<sub>μν</sub> to the field strength
F<sub>ρσ</sub>. For a linear, local, non-dispersive medium:

    H_{μν} = ½ χ_{μν}{}^{ρσ} F_{ρσ}

This is the premetric formulation (Post; Hehl and Obukhov). It is the right
starting point for B1 because it makes the central question precise: *the metric
is not an input to electrodynamics, it is a possible property of χ.*

**Component count.** With the antisymmetries of both index pairs, χ has 36
independent components. Hehl and Obukhov decompose it irreducibly:

| Part | Components | Metric can produce it? |
|---|---:|---|
| **principal** | 20 | yes — this is where a metric lives |
| **skewon** | 15 | no |
| **axion** | 1 | no |

A metric-induced χ is purely principal, and moreover a *specific* principal
form. This decomposition is what makes "can this medium be a metric?" a
well-posed question with a checkable answer, and it is the theoretical basis for
the `classify()` function in the code.

## 2. The metricity condition: no birefringence

The general χ produces a **quartic** dispersion relation — a Fresnel surface —
which generically factorises into two distinct light cones. That is
birefringence: two polarisations, two speeds.

A single pseudo-Riemannian metric produces exactly **one** light cone, so the
quartic must degenerate into the square of a quadric. The condition under which
this happens is a closure relation on χ; Lämmerzahl and Hehl established that
**vanishing birefringence in premetric vacuum electrodynamics forces a
Riemannian light cone**, i.e. the existence of a metric, unique up to a
conformal factor.

This is the deep reason for two things the code does:

- `ε ≠ μ` (impedance mismatch) is rejected: it is the generic route to
  birefringence;
- the recovered metric is a **conformal class**, never a single metric.

## 3. The Plebański–Tamm map: metric → medium

For a metric with g<sub>00</sub> < 0, Maxwell's equations in that spacetime are
identical in form to flat-space Maxwell in a medium with

    ε^{ij} = μ^{ij} = − √(−g) g^{ij} / g_00
    w_i    = − g_{0i} / g_00

    D = ε E + w × H
    B = μ H − w × E

Reference: Plebański (1960), building on Tamm; hence "Plebański–Tamm". The same
relations are reproduced in Leonhardt & Philbin (2006).

Three structural facts follow immediately, and all three are verified in the
code:

1. **Impedance matching is automatic.** ε = μ identically. A metric can *only*
   produce an impedance-matched medium. `b1.plain_dielectric_not_metric_representable`
2. **Bianisotropy appears exactly when g<sub>0i</sub> ≠ 0.** A metric with a
   shift vector is a *magnetoelectric* medium, not a dielectric.
   `b1.warp_is_bianisotropic`
3. **The map is conformally invariant.** Under g → Ω²g: √(−g) → Ω⁴√(−g),
   g<sup>ij</sup> → Ω<sup>−2</sup>g<sup>ij</sup>, g<sub>00</sub> → Ω²g<sub>00</sub>,
   so ε is unchanged (Ω<sup>4−2−2</sup> = Ω⁰) and w is a ratio of two Ω² factors.
   `b1.conformal_invariance`, and differentially `b1.sensitivity.conformal_normalisation`

## 4. Optical metrics are a *different* statement

The most common conflation in this area, and worth stating sharply.

- **Plebański map**: reproduces the **full Maxwell system**. Requires ε = μ.
- **Gordon optical metric**: describes **ray propagation only**, at the eikonal
  level. Exists far more broadly, including for ordinary dielectrics.

For an isotropic medium of index n with 4-velocity u in background g, Gordon's
optical metric is

    ḡ_{μν} = g_{μν} + (1 − 1/n²) u_μ u_ν

For a static medium in flat spacetime, u<sub>μ</sub> = (−1,0,0,0):

    ḡ = diag(−1/n², 1, 1, 1),   null coordinate speed 1/n

A plain dielectric (ε = n², μ = 1) **has** a Gordon metric but is **not**
metric-representable in the Plebański sense. The code holds these apart
deliberately: `gordon_metric_static_isotropic()` is separate from
`metric_to_medium()`, and `b1.plain_dielectric_not_metric_representable` asserts
the dielectric is rejected by the latter while §4's speed relation still holds.

Note the speeds differ and both are correct for their own claim:
- impedance-matched ε = μ = n²: speed 1/√(εμ) = 1/n²
- Gordon, index n: ray speed 1/n

`b1.sensitivity.observable_null_speed` verifies the first, `b1.gordon_optical_metric_speed`
the second.

## 5. Transformation optics

The engineering-facing half of the same mathematics. Maxwell's equations are
form-invariant under coordinate transformation if the change is absorbed into ε
and μ (Ward & Pendry). Pendry, Schurig & Smith and, independently, Leonhardt
turned this into a design method: choose a coordinate map, read off the material
that implements it.

The essential point for B1, and the reason §8 matters: a coordinate
transformation of *vacuum* produces a medium whose ε and μ are **the same
tensor**, and which is therefore always impedance matched and reflectionless.
This is the class of media that a metric can represent. Anything outside it is
not a coordinate transformation of vacuum, whatever else it may be.

## 6. RETRACTED: "where the correspondence fails"

> **This section previously claimed that no equivalent medium exists where
> v_s f > 1. That claim was wrong and is retracted.** It was an
> observer-dependent result stated as an invariant one. The corrected analysis
> is in §6b; the original text is kept below it only so the error is legible.

### 6b. Corrected: the v_s f = 1 surface is an ERGOSURFACE

The Alcubierre metric has **det g = −1 and signature (−,+,+,+) at every point,
for every v_s** — verified to v_s = 10, 0 of 405 sampled points non-Lorentzian.
The reason is algebraic: the (t,x) block has determinant −1 identically, so its
two eigenvalues always have opposite signs however large v_s f becomes.

A Lorentzian nondegenerate metric always admits the covariant constitutive
tensor Z^abcd. Verified directly: at v_s = 2, across **241 sampled points with
g₀₀ ≥ 0, Z is valid at all 241 and the metric is Lorentzian at all 241**, while
the coordinate-time observer of g_eff is inadmissible at all 241.

What fails at v_s f = 1 is that ∂ₜ stops being timelike. That is an
**ergosurface**, structurally identical to Kerr's, and the failure is removed by
changing observer: the Eulerian observer of g_eff and the laboratory coordinate
observer are both admissible at 201/201 ergoregion points tested.

This matches the established result that susceptibility eigenvalues diverge on
the ergosurface while the covariant tensor does not (Schuster & Visser 2017,
§11).

**Five distinct statuses** are now reported rather than one verdict:
`CovariantMetricMediumValid`, `SelectedObserverDecompositionInvalid`,
`ConstitutiveParametersSingular`, `NonMetricConstitutiveTensor`,
`MetricDegenerateOrNonLorentzian`. Only the last means "no metric medium here".

### 6c. Original text, retained as the record of the error

The following was the B1 conclusion. It is **superseded**.

The Plebański map requires g<sub>00</sub> < 0 — the coordinate time direction
must be timelike. For the Alcubierre metric,

    g_00 = −1 + v_s² f²

so the map **fails wherever v_s f > 1**.

For the reference parameters (v_s = 0.2) this never happens. For the
superluminal regime the project actually cares about, it happens in a shell
around the bubble wall: measured at v_s = 2, **121 of 321 sampled points on the
axis have no equivalent medium** (`b1.superluminal_medium_fails_reported`).

~~**This is the single most important realizability limitation found in B1.**~~
**RETRACTED — see §6b.** What is true is the narrower statement that the
*coordinate-observer* Plebański decomposition is singular there, by the same
condition that costs Program A its coordinate observers — see
[time_dependence_decision.md](../../warp_hypad/docs/time_dependence_decision.md) §7.
The two programmes hit the same wall from opposite directions, which is evidence
it is a property of the spacetime rather than of either implementation.

Whether a different frame or observer family rescues it is **open**, and is the
first thing B2 should address.

## 7. Reciprocity, passivity, causality, dispersion, loss

B1 implements the **non-dispersive, lossless, reciprocal** case only. What that
excludes:

- **Reciprocity.** A metric gives ξ = −ζ<sup>T</sup> with ξ antisymmetric. A
  Tellegen (axion) medium has ξ = ζ symmetric and is non-reciprocal; it cannot
  come from a metric. Verified as a rejection case.
- **Passivity and loss.** Everything here is real-valued and lossless. A real
  metamaterial has complex, frequency-dependent ε and μ, and passivity
  constrains the sign of the imaginary parts. **Not modelled at all.**
- **Causality and dispersion.** Kramers–Krönig relations force any medium with
  frequency-dependent response to be dispersive, and a medium with ε < 1 or
  ε < 0 over a band is necessarily dispersive. Several of the warp media
  computed in B1 have ε components far from 1. **The code makes no
  frequency-domain claim whatsoever.**
- **Magnitude realizability.** Nothing checks whether the computed ε, μ are
  physically achievable. Transformation-optics designs routinely demand ε → 0 or
  ε → ∞ somewhere, which is achievable only narrowband and with loss.

These are not oversights to be fixed by more code in B1 — they are the boundary
between a geometric-optics correspondence and an actual device, and they should
gate any claim about physical realizability.

## 8. Conditions summary

**A constitutive tensor is representable by a single pseudo-Riemannian metric
only if** (necessary conditions, all checked in `classify()`):

1. ε = μ (impedance matched) — else generically birefringent
2. ε, μ symmetric — else skewon part present
3. ξ = −ζ, both antisymmetric — else skewon and/or axion part present

and, for the *frame-dependent* map to exist at a point,

4. g<sub>00</sub> < 0.

**Sufficiency** requires in addition the full closure/no-birefringence relation
on χ, which B1 does **not** check: the conditions above are checked on the 3+1
blocks, not on the covariant χ. A response could in principle satisfy 1–3 and
still fail closure. This is a stated gap, not a claimed result.

## 9. What B1 verifies

20 cases, all passing. The load-bearing ones:

| Claim | Case | Result |
|---|---|---|
| Minkowski ↔ vacuum fixes normalisation | `b1.minkowski_is_vacuum` | exact |
| Forward map matches closed forms | `b1.isotropic_forward`, `b1.anisotropic_forward` | ≤ 2.2e-16 |
| Conformal invariance | `b1.conformal_invariance` | 4.4e-16 |
| …and differentially, dε/dΩ = 0 | `b1.sensitivity.conformal_normalisation` | 8.9e-16 |
| Round trip recovers the conformal class | `b1.roundtrip_conformal` | 9.9e-17 |
| …and the factor is genuinely ≠ 1 | `b1.roundtrip_factor_is_nontrivial` | Ω² up to 4 |
| Null cones preserved | `b1.null_cone_equivalence` | exact |
| Non-metric media rejected (3 kinds) | `b1.nonmetric_detected` | 3/3 |
| Valid media still accepted | `b1.valid_medium_accepted` | ok |
| Superluminal breakdown reported | `b1.superluminal_medium_fails_reported` | 121/321 |

## 10. What B1 does *not* establish

- no covariant χ<sup>μνρσ</sup> is constructed; only the 3+1 blocks
- no Fresnel/dispersion surface is computed, so birefringence is inferred from
  the ε ≠ μ condition rather than measured
- no closure-relation test (§8 sufficiency)
- the inverse map covers **w = 0 only**; magnetoelectric media are reported as
  unsupported rather than approximated (`b1.warp_inverse_reports_unsupported`)
- no frequency dependence, loss, passivity or causality
- the **sign convention of w is unverified** against the literature; it is fixed
  internally and pinned by Minkowski, conformal invariance and the round trip,
  none of which are sensitive to it. The physical handedness of the
  magnetoelectric drag is therefore **not** a claim B1 supports.
- no ray tracing, so "null cones agree" is checked along coordinate axes only

## 11. References

Verified:

- J. Plebański, "Electromagnetic Waves in Gravitational Fields", *Physical
  Review* **118**, 1396 (1960). *(Verified against the APS record.)*

Recorded from memory — **verify before citing**:

- E. J. Post, *Formal Structure of Electromagnetics*, North-Holland (1962).
- F. W. Hehl and Yu. N. Obukhov, *Foundations of Classical Electrodynamics:
  Charge, Flux, and Metric*, Birkhäuser (2003). — premetric formulation, the
  principal/skewon/axion decomposition.
- C. Lämmerzahl and F. W. Hehl, "Riemannian light cone from vanishing
  birefringence in premetric vacuum electrodynamics", *Physical Review D* **70**,
  105022 (2004).
- W. Gordon, "Zur Lichtfortpflanzung nach der Relativitätstheorie", *Annalen der
  Physik* **377**, 421 (1923). — the optical metric.
- A. J. Ward and J. B. Pendry, "Refraction and geometry in Maxwell's equations",
  *Journal of Modern Optics* **43**, 773 (1996).
- J. B. Pendry, D. Schurig and D. R. Smith, "Controlling Electromagnetic
  Fields", *Science* **312**, 1780 (2006).
- U. Leonhardt, "Optical Conformal Mapping", *Science* **312**, 1777 (2006).
- U. Leonhardt and T. G. Philbin, "General relativity in electrical
  engineering", *New Journal of Physics* **8**, 247 (2006).
- I. E. Tamm (1924), on electrodynamics in curved spacetime — the "Tamm" of
  Plebański–Tamm. *(Details not recalled with confidence.)*
