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

## 5. Christoffel symbols — declared (not implemented)

Levi-Civita connection (metric-compatible and torsion-free):

    Γ^ρ_{μν} = ½ g^{ρσ} ( ∂_μ g_{σν} + ∂_ν g_{σμ} − ∂_σ g_{μν} )

Consequences that Milestone 4 must verify independently:

- symmetry in the lower pair, Γ<sup>ρ</sup><sub>μν</sub> = Γ<sup>ρ</sup><sub>νμ</sub>;
- metric compatibility, ∇<sub>λ</sub> g<sub>μν</sub> = 0.

Metric compatibility is the important one: it is an identity the connection must
satisfy that does **not** reuse the formula above, so it is a genuine Level A
check rather than a restatement.

## 6. Riemann tensor — declared (not implemented)

    R^ρ_{σμν} = ∂_μ Γ^ρ_{νσ} − ∂_ν Γ^ρ_{μσ} + Γ^ρ_{μλ} Γ^λ_{νσ} − Γ^ρ_{νλ} Γ^λ_{μσ}

This is the **MTW / Wald** convention: signature (−,+,+,+) together with this
Riemann sign and the Ricci contraction below gives positive Ricci scalar for a
sphere and positive energy density for ordinary matter.

In the Misner–Thorne–Wheeler sign classification this project is **(+, +, +)**.

## 7. Ricci tensor and scalar — declared (not implemented)

    R_{σν} = R^ρ_{σρν}          (contract slot 1 against slot 3)
    R      = g^{μν} R_{μν}

The contraction is implemented as `contract_13`, and
`contract.contract_13_slot_selectivity` already verifies that it contracts *that*
index pair and not one of the other five possibilities.

## 8. Einstein tensor — declared (not implemented)

    G_{μν} = R_{μν} − ½ g_{μν} R

Identities Milestone 6 must verify: symmetry, the trace identity
G<sup>μ</sup><sub>μ</sub> = −R in four dimensions, and the contracted Bianchi
identity ∇<sub>μ</sub> G<sup>μν</sup> = 0.

## 9. Stress-energy tensor and units — declared (not implemented)

    T_{μν} = G_{μν} / (8π)

**Geometrized units: G = c = 1.** The factor 1/(8π) is a single configurable
normalization constant so SI can be added later without touching the geometry
code; Milestone 6 verifies that changing it rescales T exactly.

`conventions.units` currently accepts only `"geometrized"`.

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
