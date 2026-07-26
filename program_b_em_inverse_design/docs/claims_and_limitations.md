# What this implementation may and may not support

A claim-by-claim register. The purpose is that nobody — including a future
reader of a results directory — has to reconstruct from the code what is
established and what is merely computed.

**Status:** 417/417 verification cases pass across eleven gates, 0 failed,
0 skipped. 24/24 CTest with the OTIlib backend enabled, 23/23 with it disabled.
Zero compiler warnings at `/W4 /WX`.

| Gate | Cases |
|---|---:|
| Milestone 1 — grid, storage, algebra | 77 |
| Milestone A2 — spatial derivative operators | 41 |
| Milestone A2T — temporal derivatives, rigid-translation identity | 23 |
| Milestone A3 — metric generation | 31 |
| OTIlib AD backend | 23 |
| Program B1 — metric ↔ constitutive mapping | 20 |
| Program B1C-S23 — Hodge duals, Bel extraction, ζ normalisation | 42 |
| Program B1C-S456 — field reconstruction, ζ†, handedness, (A11) | 36 |
| Program B1C aggregate — ergosurface, velocity sensitivities | 53 |
| Program B2 — abstract constitutive inverse design | 33 |
| Program B3 — physical admissibility | 38 |

With the OTIlib backend disabled the HyPAD-specific cases are recorded as
skipped rather than omitted, and the totals still reconcile to the same case
counts.

---

## Claims the implementation SUPPORTS

### Numerics

- Finite-difference operators are **second order interior and boundary alike**,
  for first, second and mixed spatial derivatives. Measured 1.913–2.026 across
  24 operator/field combinations.
- Temporal derivatives are second order in all three stencil positions
  (interior, window start, window end). Measured 1.944–2.050.
- Space–time mixed derivatives **commute exactly** (5.7e-13), not merely to
  truncation order.
- 4×4 inversion, determinant and contractions agree with independently written
  references; singular matrices are rejected 200/200 and ill-conditioning is
  flagged on a derived, not chosen, threshold.

### Automatic differentiation

- `d(A⁻¹)/dp = −A⁻¹(dA/dp)A⁻¹` holds to **1.9e-16**. This is the gate for using
  OTI on inverse metrics.
- Metric parameter sensitivities match hand-derived closed forms to **1.1e-16**.
- Gradient *and* Hessian come from a single evaluation (`onumm3n2`, 3 directions,
  order 2), verified against analytic second derivatives.
- OTI channels pass through the temporal stencil correctly (4.5e-08 against an
  independent single-difference reference).

### Metric

- Minkowski is exact; Alcubierre has **det g = −1 identically**, is symmetric,
  Lorentzian at every tested point, and asymptotically flat.
- The numerical inverse matches the closed-form inverse to 1e-14.
- **∂ₜg = −v_s ∂ₓg** for a rigidly translating bubble: the residual between an
  analytic ∂ₜ and an independently computed numerical ∂ₓ converges at order 1.82
  toward zero (0.0194 → 0.0056 → 0.0016 over 4× refinement). A false identity
  would level off at order 0.

### Electromagnetic correspondence (B1)

- Minkowski ↔ vacuum exactly, in Heaviside–Lorentz units.
- The Plebański forward map matches closed forms for isotropic and diagonal
  anisotropic metrics (≤ 2.2e-16).
- **Conformal invariance** holds to 4.4e-16, and *differentially*
  dε/dΩ = 0 to 8.9e-16.
- A round trip g → medium → ĝ recovers the **conformal class** (9.9e-17), with
  the recovered factor genuinely ≠ 1, and preserves null cones exactly.
- Non-metric constitutive tensors — birefringent, axion/Tellegen, skewon-like —
  are **detected and reported**, 3/3, while a valid medium is still accepted.
- The warp medium is **bianisotropic**: w_x ≈ −0.101 at the wall. No arrangement
  of ordinary dielectrics reproduces it.

### Dispersive models (B3)

Gate B3-PHYSICAL-ADMISSIBILITY, 38/38.

- The Fourier convention is `e^{-i omega t}`, fixed once and **pinned** by the
  requirement that the specified Lorentz and Drude forms come out passive rather
  than active under it.
- **Causality is GLOBAL**, and structural: damping passes through a strictly
  positive transform, so every pole of every term lies in the lower half
  omega-plane for every parameter value. No sampling is involved and no
  optimizer step can escape it.
- **Passivity is established OVER THE VERIFIED BAND ONLY.** `min eig L >= 0` is
  checked on a finite set of sampled frequencies. **Global passivity is not
  claimed.** `L(omega)` is rational and can dip negative in a narrow interval
  between samples — especially near a resonance, a zero of `det mu`, or an
  eigenvalue crossing — and a sampled gate cannot see that. Proving it globally
  needs an analytic positive-real argument that has not been carried out here.
  The adaptive search introduced at B4 narrows the unverified region; it does not
  close it, and it is not a replacement for the B3 sampled gate.
- Kramers–Krönig is a **diagnostic**, reported as a convergence trend, and is
  explicitly not the causality proof.
- `hard feasible approximation` and `soft penalized approximation` are reported
  under different names, and a soft-penalty result is never called feasible
  merely because its violation is small.
- A recovered model reproduces a **response** over a band. That is not the same
  claim as recovering **parameters**: below resonance the model is not
  identifiable, and the gate measures both.

### Covariant correspondence (B1C, steps 1–6)

Gates B1C-S23 (42/42) and B1C-S456 (36/36), both 0 failures and 0 skips.

- Four-dimensional Hodge duals, with `**F = −F` and the double-dual identities
  verified, and the contravariant Levi-Civita sign **derived** rather than
  chosen.
- The Bel extraction of Schuster & Visser (A7a)–(A7d), with the normalisations
  ε = −2W_ε, μ⁻¹ = +2W_μ, ζ = +2W_ζ, ζ† = +2W_ζ^T — every one of them measured
  against a medium with independent ε ≠ μ⁻¹ or a nonzero magnetoelectric
  response, never inferred from vacuum.
- **Full field reconstruction**: `H^{ab} = Z^{abcd}F_{cd}` computed directly and
  rebuilt from (E,B) → (D,H) through the four blocks agree to **2.7e-15** over
  five media (vacuum, isotropic at rest with ε ≠ μ⁻¹, diagonal anisotropic,
  moving isotropic, warp analogue with g_lab = η) and twelve two-forms (six
  basis, six random), plus eight random reciprocal media.
- E, B, D and the magnetic excitation are spatial to **0.0**; F is recovered
  from (E,B) to **0.0**; the chain is linear to 1.8e-15.
- **Z is reconstructed from the observer blocks** via (A11) to **2.2e-16**,
  including random constitutive tensors carrying only the symmetries (12) and
  (13). Nothing is fitted componentwise.
- The blocks are observer dependent and the tensor is not: for a fixed warp
  Z, two admissible laboratory observers see blocks differing by 6.7e-1 and
  rebuild the same Z to 4.4e-16.
- The whole chain is covariant under a nontrivial linear coordinate change that
  takes g_lab away from Minkowski: 1.1e-14.

**Scope limit:** (A11) reconstructs exactly the **reciprocal** media, those with
`Z^{abcd} = Z^{cdab}`. That covers everything the effective-metric programme
produces, since equation (6) is manifestly pair-symmetric, but a skewon-carrying
response is outside it. The residual is measured at every call site rather than
assumed.

---

## Claims the implementation DOES NOT support

### Nothing about curvature

There is no Christoffel symbol, Riemann, Ricci, Einstein or stress-energy
tensor. **No statement about warp-drive energy requirements, negative energy, or
energy conditions is supported by anything in this repository.**

### The frozen-time mode is not physics

`TimeTreatment::FrozenTimeDiagnostic` discards ∂ₜg rather than computing it. Its
outputs are legitimate for metric components, determinant, signature, validity
masks and parameter sensitivities. They are **not** legitimate inputs to any
curvature calculation, and `valid_for_curvature()` exists so the curvature code
can refuse rather than warn. The measured discrepancy against the correct
time-dependent metric at t = 2, v_s = 0.35 is 0.20 in g₀₀ — not a small
correction.

### Sensitivity is "at fixed bubble trajectory"

`center_at()` advances the bubble using `real_part(velocity)`, so OTI computes
∂/∂v **at fixed trajectory** — the shape sensitivity, not the total lab-frame
derivative, which would add a d(x_s)/dv = t term. The two differ substantially
at t ≠ 0; an earlier version of the test conflated them and reported an 88%
"error" between two different quantities.

**This vanishes at t = 0**, and since the objective for a rigidly moving bubble
is invariant under time translation, all optimization can be done at t = 0 where
the two coincide. Revisit if sensitivities at fixed co-moving position are ever
wanted.

### The magnetoelectric sign — resolved for the covariant blocks, in B1C

This limitation applied to the **B1 3+1 route**, where the sign of w was fixed
internally and pinned only by tests insensitive to it. It no longer applies to
the covariant blocks.

For the covariant magnetoelectric tensor ζ the magnitude, the tensor sign and
the physical handedness are now all determined and gated (B1C-S23 and
B1C-S456):

* ζ = +2 W_ζ, measured against Schuster & Visser (B17) across five boost
  configurations; the opposite sign is recorded as a distance, not merely
  excluded.
* ζ† = +2 W_ζ^T, measured against the transpose of (B17).
* Handedness: for V = (1,0,0,0), a medium with ε > μ⁻¹ boosted along **+x**, and
  B = +B₀ ẑ with E = 0, the induced displacement points along **−y**. Asserted
  as one direction; a +y result fails the gate. Cross-checked at β = 1e-3
  against Minkowski's moving-medium relations, a reference outside the paper.

**What is still not established:** this is the handedness of the *constitutive
response* of the effective medium, verified algebraically. It is not a
propagation result. No claim about the direction of drag on a propagating wave
is supported until B2 exists, and the B1 3+1 `w`-sign convention in sections 2.1
and 3 of the conventions document remains unpinned on its own terms.

### B1 does not test metricity sufficiently

`classify()` checks **necessary** conditions on the 3+1 blocks (ε = μ, symmetry,
ξ = −ζ antisymmetric). It does **not** construct the covariant χ^μνρσ, compute a
Fresnel surface, or test the closure/no-birefringence relation. A response could
satisfy all the implemented conditions and still fail closure.

### The inverse map covers w = 0 only

Every warp metric has w ≠ 0, so **the inverse map does not currently apply to
the case of interest**. It returns `MagnetoelectricUnsupported` rather than
approximating.

### No frequency domain at all

Non-dispersive, lossless, real-valued, reciprocal. No Kramers–Krönig, no
passivity constraint, no bandwidth. Several computed warp media have ε far from
1, which in a real material forces dispersion and loss. **No claim about
physical realizability of any computed medium is supported.**

### No magnitude realizability check

Nothing tests whether a computed ε or μ is achievable. Transformation-optics
designs routinely demand ε → 0 or ∞ somewhere.

---

## RETRACTED CLAIM (B1 → B1C)

**B1 claimed: "no equivalent medium exists" where v_s f > 1. That is withdrawn.**

It was computed from the Plebański 3+1 map relative to the coordinate-time
observer, and stated as though it were invariant. The corrected position:

- the Alcubierre metric is **Lorentzian with det g = −1 everywhere**, for every
  v_s tested up to 10 (0 of 405 points non-Lorentzian);
- therefore the **covariant tensor Z^abcd exists everywhere** — verified at
  241/241 ergoregion points at v_s = 2;
- what fails is that **∂ₜ stops being timelike**: an ergosurface;
- **other observers work there** — Eulerian and laboratory-coordinate observers
  admissible at 201/201 ergoregion points.

The correct statement is *"the coordinate-observer Plebański decomposition is
singular on this surface"*, not *"no medium exists"*. Whether a realizable
laboratory material exists there remains **open** and is a B2 question.

### What actually diverges at the ergosurface

A second, smaller retraction, made while repairing the aggregate gate. B1C
originally asserted that the **condition number of ε** must diverge as
v_s f → 1. It cannot, and the assertion is **withdrawn**: for the laboratory
observer ε^ab = h^ab identically, at every point and every v_s, so cond(ε) ≡ 1.
More generally a tensor may scale uniformly towards infinity with its condition
number pinned at 1, so a conditioning test is the wrong instrument for detecting
a blow-up.

What does diverge, under the same observer and the same conventions, is the
**permeability**. The spectrum of [μ⁻¹] for the Alcubierre analogue is
{1, −g₀₀, −g₀₀} exactly, so its transverse pair reaches zero precisely on
v_s f = 1 and μ = (μ⁻¹)⁻¹ blows up. Measured at the sampled point nearest the
surface (|g₀₀| = 2.46e-5): min|eig μ⁻¹| = 2.46e-5, cond(μ⁻¹) = 4.07e4,
|det μ⁻¹| = 6.04e-10, against 1, 1, 1 far away — while max|Z^abcd| changes by
6.1e-6 across the same surface and stays bounded in (0.1, 2) over the whole
sweep. Divergence is asserted only for the quantity predicted to diverge.

Full record: [b1c_aggregate_gate_audit.md](b1c_aggregate_gate_audit.md).

## The observer-dependent limitation, reached from two directions

**Where v_s f > 1, g₀₀ ≥ 0**, and:

- Program A loses its coordinate observers — ∂ₜ is no longer timelike, so
  "energy density measured by the coordinate observer" is undefined.
- Program B loses the equivalent medium — the Plebański map requires g₀₀ < 0.
  Measured at v_s = 2: **121 of 321 axis points have no equivalent medium.**

For the reference parameters (v_s = 0.2) neither occurs. **For the superluminal
regime the project actually cares about, both do**, in the same region.

This is a property of the spacetime and the coordinates, not of either
implementation, and it is not fixed by switching frames — the co-moving
formulation makes it *worse*, losing coordinate observers over the entire
exterior rather than a shell near the wall.

Two programmes reaching the same *surface* from opposite directions is evidence
the surface is real — but it is an **ergosurface**, not a realizability wall.
Both failures are failures of the same observer choice, not of the geometry. It should be treated as the primary open scientific
question, ahead of any optimization work:

- Program A must choose an observer family that exists there (Milestone 7).
- Program B must establish whether any frame or generalised medium supports the
  correspondence there (B2).

---

## Open questions, in priority order

1. **Does the metric ↔ medium correspondence survive into v_s f > 1 in any
   frame?** If not, the electromagnetic-analogue route has a hard ceiling below
   the regime of interest, and that should be known before more is built on it.
2. **Which observer family** replaces the coordinate observers where they fail?
3. **Closure relation** — implement the covariant χ and a genuine birefringence
   test, so metricity is decided sufficiently rather than necessarily.
4. **w ≠ 0 inverse**, so the warp medium can be round-tripped.
5. **Co-moving cross-check** at Milestone 5: curvature invariants must agree
   between the inertial and co-moving formulations while components differ.
6. Dispersion, loss and passivity, before any realizability claim.
