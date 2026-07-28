# Sign, index, frame and unit conventions

Every convention this project will ever depend on is fixed here, in one place.
Ambiguity about index ordering or curvature sign is the single most common
source of silent error in a numerical relativity code, and it is the thing that
makes cross-code comparison against Warp Factory (Milestone 10) meaningless if
left implicit.

**Status.** Sections marked *(declared)* describe conventions that are fixed but
whose implementation does not exist yet. They are written down now so that the
implementation is checked against a prior commitment rather than the commitment
being reverse-engineered from whatever the code happens to do. Sections marked
*(implemented and verified)* are backed by passing cases in
`results/verification/milestone_01/`.

---

## 1. Signature — implemented and verified

**(−, +, +, +).**

Minkowski in Cartesian coordinates is `diag(−1, 1, 1, 1)`, so `det η = −1` and a
timelike vector has `g(u, u) < 0`.

Verified by `algebra.minkowski_self_inverse` (η is its own inverse, exactly) and
`contract.trace_minkowski` (the diagonal sums to 2, not 4 or −2).

Only this signature is implemented. `config_validate` rejects any other value of
`conventions.metric_signature` rather than ignoring it.

## 2. Coordinate index map — implemented and verified

| Index | Coordinate |
|---|---|
| 0 | t |
| 1 | x |
| 2 | y |
| 3 | z |

Greek indices μ, ν, ρ, σ, α, β, λ run 0–3 over spacetime. Latin i, j, k are
**grid** indices over x, y, z and are never tensor indices.

Recorded in `metadata.json` under `layout.coordinate_index_map` so a reader of a
results directory does not need the source.

## 3. Tensor index placement — implemented and verified

The containers (`Mat4`, `Rank3`, `Rank4`) store numbers and do not track
variance. The project-wide convention is that **the leading slot is the raised
one**, matching how the objects are written on paper:

| Object | Access | Meaning |
|---|---|---|
| metric | `g(mu, nu)` | g<sub>μν</sub> |
| inverse metric | `g_inv(mu, nu)` | g<sup>μν</sup> |
| Christoffel | `gamma(rho, mu, nu)` | Γ<sup>ρ</sup><sub>μν</sub> |
| Riemann | `riem(rho, sig, mu, nu)` | R<sup>ρ</sup><sub>σμν</sub> |
| Ricci | `ric(mu, nu)` | R<sub>μν</sub> |
| Einstein | `g_tensor(mu, nu)` | G<sub>μν</sub> |
| stress-energy | `t(mu, nu)` | T<sub>μν</sub> |

`Mat4` is row-major: `m(row, col)` is at flat offset `4*row + col`. This is
verified against a hand-calculated asymmetric product
(`algebra.matmul_hand_calculated`), which distinguishes it from the transposed
convention — a symmetric or identity test case could not.

## 4. Grid and storage layout — implemented and verified

See [memory_layout.md](memory_layout.md). Summary:

- flat grid index `p = (i·n_y + j)·n_z + k`, so a walk along +z is contiguous;
- fields are **point-major**: `raw()[16p + 4μ + ν]` for rank 2.

## 5. Christoffel symbols — implemented and verified (A4)

Levi-Civita connection (metric-compatible and torsion-free):

    Γ^ρ_{μν} = ½ g^{ρσ} ( ∂_μ g_{σν} + ∂_ν g_{σμ} − ∂_σ g_{μν} )

Implemented as `geometry::christoffel_at` / `christoffel_field`, routed through
the Milestone 1 `contract_first` so the first-kind-to-second-kind step is the one
already verified there.

Verified independently:

- symmetry in the lower pair, Γ<sup>ρ</sup><sub>μν</sub> = Γ<sup>ρ</sup><sub>νμ</sub>
  — exactly 0;
- metric compatibility, ∇<sub>λ</sub> g<sub>μν</sub> = 0 — 3.3e-17;
- against the closed form for flat space in spherical coordinates
  (Γ<sup>r</sup><sub>θθ</sub> = −r, Γ<sup>θ</sup><sub>rθ</sub> = 1/r,
  Γ<sup>φ</sup><sub>θφ</sub> = cot θ) — 3.4e-4 interior.

Metric compatibility is the important one: it is an identity the connection must
satisfy that does **not** reuse the formula above, so it is a genuine Level A
check rather than a restatement. It comes out exact because its derivation never
differentiates a product — see §7a.

## 6. Riemann tensor — implemented and verified (A4)

    R^ρ_{σμν} = ∂_μ Γ^ρ_{νσ} − ∂_ν Γ^ρ_{μσ} + Γ^ρ_{μλ} Γ^λ_{νσ} − Γ^ρ_{νλ} Γ^λ_{μσ}

This is the **MTW / Wald** convention: signature (−,+,+,+) together with this
Riemann sign and the Ricci contraction below gives positive Ricci scalar for a
sphere and positive energy density for ordinary matter.

In the Misner–Thorne–Wheeler sign classification this project is **(+, +, +)**.

**The sign is verified, not merely declared.** Static de Sitter is the only
reference geometry in the A4 suite with nonzero Ricci, and
`curvature.desitter_ricci_scalar_sign` asserts R > 0 there. A pipeline built with
the opposite Riemann sign returns −12/L² and passes every index identity in the
suite; that one case is what catches it. What it fixes is the **sign convention**
that A5 and A6 inherit — not the sign of any energy. See §7b.

## 7. Ricci tensor and scalar — implemented and verified (A4)

    R_{σν} = R^ρ_{σρν}          (contract slot 1 against slot 3)
    R      = g^{μν} R_{μν}

The contraction is `contract_13`, and `contract.contract_13_slot_selectivity`
verifies that it contracts *that* index pair and not one of the other five
possibilities.

Curvature invariants, added at A4:

    R_{μν} R^{μν}        `ricci_square`
    R_{αβγδ} R^{αβγδ}    `kretschmann`

## 7a. Which curvature identities survive discretisation — measured (A4)

A result, not a convention, and it decides how every A4 tolerance is written.

The identities whose derivation is a rearrangement of the same discrete
∂<sub>α</sub> g hold **exactly**, at roundoff:

| Identity | Residual |
|---|---|
| Γ<sup>ρ</sup><sub>μν</sub> = Γ<sup>ρ</sup><sub>νμ</sub> | 0 |
| ∇<sub>λ</sub> g<sub>μν</sub> = 0 | 3.3e-17 |
| R<sup>ρ</sup><sub>σμν</sub> = −R<sup>ρ</sup><sub>σνμ</sub> | 0 |
| first Bianchi, R<sup>ρ</sup><sub>[σμν]</sub> = 0 | ~1e-16 |
| ∂<sub>0</sub> g = −v<sub>s</sub> ∂<sub>x</sub> g (aligned step, §11) | 0 |

The identities whose textbook derivation passes through the **product rule** do
not, because a difference operator does not satisfy it: D(fg) differs from
(Df)g + f(Dg) by O(h²) f′g′. These are second-order truncation statements:

| Identity | Residual (17³, v=0.5 R=2 σ=1) | Measured order |
|---|---|---|
| R<sub>ρσμν</sub> = −R<sub>σρμν</sub> (lowered) | 7.4e-3 | — |
| R<sub>ρσμν</sub> = R<sub>μνρσ</sub> | 5.1e-3 | — |
| R<sub>μν</sub> = R<sub>νμ</sub> | 4.7e-3 | 2.08 |

Writing any of the second group with a roundoff tolerance and then widening it
until it passed would have converted a check into a calibration. They are gated
instead at twice the measured value **with** a companion convergence case, and
Milestone A5 inherits exactly this asymmetry when it forms
G<sub>μν</sub> = R<sub>μν</sub> − ½ g<sub>μν</sub> R.

## 8. Einstein tensor — implemented and verified (A5)

    G_{μν} = R_{μν} − ½ g_{μν} R

Verified: Minkowski and flat-spherical give exactly and effectively zero
(6.7e-3 interior on a large connection — a coordinate connection does *not*
become a false source); Schwarzschild gives 1.1e-3 while remaining curved;
static de Sitter gives the analytic −(3/L²)g<sub>μν</sub>; and the Alcubierre
components agree with an independent symbolic reference to 8.9e-3.

**Raw and symmetric are stored separately.** A4 established that Ricci symmetry
is O(h²), not roundoff, and G inherits that. `raw_einstein` preserves the
directly computed tensor and its antisymmetry is *measured* (3.1e-3, converging
at order 2.08); `symmetric_einstein` is the explicitly labelled projection
½(G<sub>μν</sub> + G<sub>νμ</sub>) used downstream. Symmetrising in place would
destroy the only direct measure of the discretisation error, and a large
antisymmetric part would become invisible.

| Identity | Class | Residual |
|---|---|---|
| G<sup>μ</sup><sub>μ</sub> = −R | exact | 2.8e-17 |
| ½(G + Gᵀ) is symmetric | exact | 0 |
| ∇<sub>μ</sub>G<sup>μν</sup> = 0 | truncation | 6.0e-3 absolute, 4.3e-2 relative |

The trace identity is exact because g<sup>μν</sup>g<sub>μν</sub> = 4 holds for
whatever R<sub>μν</sub> the discretisation produced — it needs no product rule.
The contracted Bianchi identity is **not** an exact finite-difference identity:
it costs a third nested derivative pass. Both an absolute and a relative
residual are reported, because a vacuum's absolute residual is tiny purely
because its G is tiny.

## 8a. Cosmological constant — the source decomposition (A5, corrected at A5C)

**A5 got this wrong and A5C fixes it.** A5 carried an enum member named
`LambdaAsEffectiveStressEnergy` whose arithmetic computed
(G<sub>μν</sub> + Λg<sub>μν</sub>)/8π. That quantity is the **matter** source.
The effective cosmological-constant source is −Λg<sub>μν</sub>/8π — a different
tensor, not a sign slip. The name asserted one physical interpretation while the
code implemented another, which is the exact failure this project's naming
discipline exists to prevent.

The fix is not a rename. A single policy-selected "the source" field cannot
express the decomposition, because there are three distinct tensors and a caller
needs all of them. `CosmologicalConstantTreatment` is **removed**, not renamed,
and replaced by `SourceDecomposition`.

Writing the field equation with Λ on the geometric side:

    G_{μν} + Λ g_{μν} = 8π T^matter_{μν}

| Field | Definition | Λ = 0 | de Sitter, Λ = 3/L² |
|---|---|---|---|
| `geometric_einstein` | G<sub>μν</sub> | G | −Λg, **nonzero** |
| `matter_stress_energy` | (G + Λg)/8π | G/8π | **0** — a vacuum |
| `lambda_effective_stress_energy` | −Λg/8π | **0** | −Λg/8π |
| `total_effective_stress_energy` | G/8π | G/8π | = T<sup>Λ</sup> |

Four identities, all verified exact at ~1e-16:

    8π T^matter = G + Λg      8π T^Λ = −Λg
    T^total = T^matter + T^Λ  8π T^total = G

The last is the closure that makes it a decomposition. `T^total` is independent
of Λ, because moving Λ between the two sides of the field equation is
bookkeeping and not physics.

de Sitter is where the distinction is sharpest and where A5's error would have
misled: it contains **no matter** (T<sup>matter</sup> → 5.2e-4) and is
nonetheless **sourced** (|T<sup>Λ</sup>| = 0.068). A5's mislabelled field
computed T<sup>matter</sup> while calling it the Λ source, so a reader would
have concluded de Sitter's Λ source was zero.

**The A5 numerical results are unaffected.** Warp-drive metrics use Λ = 0, where
the decomposition collapses to the single tensor A5 reported: T<sup>matter</sup>
= T<sup>total</sup> = G/8π exactly, and T<sup>Λ</sup> = 0 exactly.

## 9. Stress-energy tensor and units — implemented (A5)

    T_{μν}^inferred = G_{μν} / (8π)          geometrized, G = c = 1

**The name is deliberate.** The tensor is *inferred* by inserting a chosen
geometry into Einstein's equation and reading off what source would have to be
present. It is not a material model, not a physical source, not a fabricated
source and not a realizable source; no mechanism generating it has been
demonstrated. The A5 gate checks the naming.

**Coordinate components are not energy densities.** G<sub>00</sub> is positive
throughout the sampled Alcubierre points, and that is *not* a statement that
the inferred source has positive energy density. An energy density is
T<sub>μν</sub>u<sup>μ</sup>u<sup>ν</sup> for a specified normalised observer,
which A6 introduces. No A5 output may be labelled an energy density.

**SI conversion requires an explicit scale and refuses to invent one.** A grid
coordinate is a pure number, not metres:

    G^SI = G / L₀²                    1/m²
    T^SI = c⁴/(8πG_N) · G^SI          J/m³ = Pa

with `c⁴/(8πG_N) ≈ 4.82e42 Pa per 1/m²`. `si_stress_energy` returns NaN when no
`CoordinateScale` is supplied, and every converted output records the scale.
T^SI scales as 1/L₀², so the *same* dimensionless geometry means a very
different source depending on how large you declare it to be.

## 10. Observer frame and energy density — declared (not implemented)

`T_00` **is not** "energy density" in general. It is the energy density measured
by an observer at rest in the coordinate frame, which is not a physically
meaningful observer for a warp metric.

The energy density this project reports is

    ρ = T_{μν} u^μ u^ν

for an explicitly specified, normalized four-velocity `u` with
g<sub>μν</sub> u<sup>μ</sup> u<sup>ν</sup> = −1. Milestone 7 exists to make this
distinction and to verify the normalization; until then, no quantity in this
project may be labelled "energy density".

## 7b. What a curvature invariant does NOT tell you — interpretation guard

Two inferences are invalid and are called out here because both are easy to make
from an A4 figure, and both would poison every later milestone.

**A sign change in R or in the Kretschmann invariant says nothing about the sign
of an energy.** Curvature invariants characterise geometry. An energy density is
`ρ = T_{μν} u^μ u^ν` for an explicitly specified, normalised observer, and it
requires the Einstein tensor (A5), a stress-energy interpretation (A5) and a
declared observer family (A6). Neither R nor `R_{αβγδ}R^{αβγδ}` is proportional
to any energy density, and the Alcubierre invariants change sign across the wall
for reasons that are purely geometric. Nothing in A4 licenses reading those
lobes as regions of positive or negative energy.

The remark in §6 that these conventions give "positive energy density for
ordinary matter" is a statement identifying *which* sign convention this project
uses, in the standard MTW classification. It is not a rule for inferring an
energy from a curvature invariant.

**Finite curvature at g₀₀ = 0 is a statement about curvature and nothing else.**
The A4 velocity sweep establishes exactly one thing at the crossing: every
curvature invariant stays finite, and `det g = −1` identically. That is *not*:

| Not established by A4 | Why |
|---|---|
| no horizon | a horizon is a global causal-structure property of null congruences, not a pointwise curvature statement |
| no causal disconnection | requires tracing null geodesics, which does not exist in this repository |
| safe traversal | requires tidal accelerations on a worldline, i.e. geodesic deviation |
| physical realizability | requires a source, energy conditions, and a mechanism — none of which exist |

`g₀₀ = 0` in these coordinates marks where the coordinate-stationary observer
stops being timelike. Whether that surface is a horizon, and what lies causally
beyond it, are questions A4 does not ask and cannot answer.

## 10a. Coordinate charts — implemented (A4)

Up to A3 the three grid axes were always (x, y, z) and the identification was
implicit. A4 breaks that: the strongest closed-form references are written in
(r, θ, φ), and flat-space-in-curvilinear-coordinates exists precisely to separate
"the connection vanishes" from "the curvature vanishes".

Nothing about the finite-difference operators changes — a stencil differentiates
with respect to a *grid axis*, and which coordinate that axis carries is a
property of the metric. What changes is that a results directory is ambiguous
unless the chart is recorded, so `geometry::Chart` is a first-class object.

| Chart | Grid axes | Coordinate singularities to exclude |
|---|---|---|
| `CartesianTXYZ` | x, y, z | none |
| `SphericalTRThetaPhi` | r, θ, φ | r = 0, sin θ = 0, and r = 2M for Schwarzschild |
| `IsotropicTRhoThetaPhi` | ρ, θ, φ | ρ = 0, sin θ = 0, ρ = M/2 |

A grid that straddles one of those produces large finite numbers that look like
physics. Every A4 domain excludes them by construction and records the exclusion
in the case's parameters.

## 11. Time dependence — DECIDED for A3: frozen-time spatial slice

**A3 implements `TimeTreatment::FrozenTimeDiagnostic`, explicitly and by name.**

The bubble centre is held at a fixed coordinate `bubble_center_x`, so
∂<sub>t</sub> g<sub>μν</sub> = 0 **by construction rather than by assumption**.
`metric.evaluate(t, x, y, z, params)` accepts `t` and ignores it; the
`metric.frozen_time_contract` case verifies the metric is bit-identical for
t = 0, 1, −5, 1000.

**What this is NOT.** It is a spatial slice of a moving spacetime, not the
moving spacetime. The Alcubierre bubble genuinely moves, so
∂<sub>t</sub> g<sub>μν</sub> ≠ 0 in the coordinate frame, and **curvature
computed from this slice alone is not the complete physical Alcubierre
spacetime curvature.** Every A3 results directory records this in
`metric_metadata.json`, and the gate prints it on every run.

Three things keep the future path open:

- `evaluate` already takes `t`, so no call site changes when the bubble moves;
- `TimeTreatment::MovingBubble` exists as the unimplemented branch that will use
  x<sub>s</sub>(t);
- `spacetime_gradient` takes the ∂<sub>0</sub> field as a **required** argument.
  The frozen-time form is a separately named function,
  `spacetime_gradient_frozen_time`, so the approximation is spelled out at each
  call site rather than defaulting silently.

The remaining choice — co-moving frame versus an explicit time axis — is still
open and must be settled before Milestone 5 computes curvature it intends to
call physical. It no longer blocks A2 or A3.

## 12. Numerical derivative conventions — implemented and verified (A2)

- second-order centred differences in the interior;
- second-order one-sided differences on the boundary;
- interior and boundary errors reported **separately**, always.

Boundary error is reported separately because a one-sided stencil has a larger
error constant, and averaging it into a global norm hides exactly the artefact
that shows up in warp-bubble curvature plots.

The one-sided second-derivative stencil uses **four** points,
(2f₀ − 5f₁ + 4f₂ − f₃)/h². The familiar three-point form is only O(h) at an
endpoint and would silently drop the boundary to first order while the interior
still reported second.

Measured at the finest resolution pair (n = 41 → 81), across 24 operator/field
combinations: interior order 1.963–2.000, boundary 1.913–2.026.

Mixed derivatives are composed first derivatives. Because each operator touches
only its own axis index, the two orderings commute **exactly** rather than to
truncation order — verified to 7.9e-14, i.e. arithmetic roundoff.

## 13. Floating-point conventions — implemented

The build sets `/fp:precise` (MSVC) or `-ffp-contract=off` (GCC/Clang). This is
required rather than incidental: the verification suite asserts that `double`
and the custom scalar type produce **bit-identical** results, which assumes the
compiler does not reassociate or fuse floating-point expressions differently
between two template instantiations.
