# Time dependence: technical decision

**Question.** The Alcubierre bubble moves. Milestones A3 evaluated the metric on a
frozen slice. Before any curvature is computed, the project must decide how time
dependence enters — because a Riemann tensor built from spatial derivatives
alone is the Riemann tensor of a *different spacetime*, not an approximation to
the right one.

**Recommendation: Formulation 1 — explicit inertial-frame time dependence
through x_s(t).** Reasoning in §5. Formulation 2 is retained as an independent
cross-check, Formulation 3 is retained as a diagnostic mode that may never
support a physical claim.

---

## 0. The result that drives the decision

For a bubble translating rigidly at constant velocity, x_s(t) = v_s t, every
metric component depends on t and x **only through the combination**

    u = x − v_s t

because r_s = √((x − v_s t)² + y² + z²) and every component is a function of
r_s alone (times constants). Therefore, by the chain rule,

    ∂_t g_μν = (∂g_μν/∂u)(∂u/∂t) = −v_s ∂g_μν/∂u
    ∂_x g_μν = (∂g_μν/∂u)(∂u/∂x) = +∂g_μν/∂u

    ⟹  ∂_t g_μν = −v_s ∂_x g_μν        (EXACT, constant v_s, rigid profile)

This is not an approximation and not a small-velocity expansion. It is an
algebraic identity of rigid translation.

Verified numerically before this document was written (v = 0.35, R = 3, σ = 2,
four spacetime points, centred differences at h = 1e-6): agreement to 3e-10 …
4e-9 relative, with the single 4e-5 outlier occurring where the derivative
itself is ~4e-7 and the discrepancy is finite-difference noise rather than a
real disagreement.

**Consequence.** The temporal derivative needed for the full curvature of a
constant-velocity warp bubble costs *nothing beyond a spatial derivative the
project already has, and has already verified to second order*. The usual
argument for the co-moving frame — "time dependence is expensive" — does not
apply here.

It also hands us a free Level A verification: compute ∂_t analytically and
compare against −v_s × (the numerical ∂_x). The two routes share no code.

**Where the identity fails.** It requires constant v_s and a rigid profile. An
accelerating bubble, v_s(t), or a time-varying R(t) or σ(t), breaks it, and a
genuine temporal derivative operator is then required. That is why the temporal
operators are implemented and verified regardless (§6) — they are the general
mechanism, and the identity is the fast path plus the check on it.

---

## 1. Formulation 1 — explicit inertial-frame time dependence

### Metric

    ds² = −dt² + (dx − v_s f(r_s(t)) dt)² + dy² + dz²

    r_s(t) = √((x − x_s(t))² + y² + z²),     x_s(t) = ∫ v_s dt′

Components in the (−,+,+,+) signature:

    g_00 = −1 + v_s² f²        g_0x = g_x0 = −v_s f
    g_xx = g_yy = g_zz = 1     all others 0

### Coordinate system

Inertial Cartesian (t, x, y, z). These are the coordinates of an observer at
rest at infinity — the frame in which the bubble is seen to move.

### Required temporal derivatives

∂_t g_μν ≠ 0, and ∂_t² g_μν ≠ 0. For constant v_s both reduce to spatial
derivatives via §0:

    ∂_t g = −v_s ∂_x g              ∂_t² g = v_s² ∂_x² g
    ∂_t ∂_x g = −v_s ∂_x² g

For variable v_s(t) they must be computed with genuine temporal stencils.

### Observer interpretation

Coordinate observers (u^μ ∝ δ^μ_0) are asymptotically the inertial rest
observers. This is the physically meaningful family: "the energy density seen by
someone watching the bubble fly past" is the quantity that makes the
negative-energy requirement mean what it is usually taken to mean.

Valid only where g_00 < 0. For v_s f > 1 the coordinate ∂_t is spacelike and
coordinate observers cease to exist — see §7.

### Compatibility with the existing Cartesian grid

Complete. The spatial grid is unchanged; t enters as a scalar argument that
`evaluate(t, x, y, z, params)` already accepts. No grid rework, no time axis
required for constant v_s.

### Implications for Christoffel and curvature

Γ^ρ_μν needs the full ∂_α g_μν including α = 0. The `Rank3Field` layout already
reserves slot 0 for it, and `spacetime_gradient` already takes the ∂_0 field as
a **required** argument. Filling it is the whole change.

Curvature is then the genuine curvature of the moving spacetime, and the
existing checks (Minkowski → 0, metric compatibility, Bianchi) apply unmodified.

### Implications for OTI sensitivities

Clean. Coordinates stay `double`; parameters stay `Scalar`. The identity
∂_t g = −v_s ∂_x g is linear in v_s, and v_s is a *seeded* parameter, so the
product rule propagates automatically:

    ∂/∂v_s (∂_t g) = −∂_x g − v_s ∂_x(∂g/∂v_s)

OTI handles this without special-casing, since ∂_x is a linear operator with
constant coefficients acting on OTI-valued fields. No new seeding logic.

### Verification opportunities

- **Level A**: ∂_t g computed analytically vs −v_s × numerical ∂_x. Shares no
  code between the two sides.
- **Level A**: v_s = 0 ⟹ ∂_t g ≡ 0 exactly.
- **Level E**: temporal convergence order of the general operator, verified
  against the identity.
- **Level D**: curvature invariants must equal those of Formulation 2 (§2).
- **Level B**: asymptotically g → η, so the existing asymptotic-flatness test
  survives unchanged.

### Numerical complexity

For constant v_s: **no extra cost** — one already-computed spatial derivative
field, scaled. For variable v_s: one extra metric evaluation per temporal
stencil tap (3 for first order, 4 one-sided), i.e. ≤4× the metric evaluation
cost, which is far cheaper than the curvature contraction that follows it.

### Risk of confusing coordinate effects with physical curvature

**Low.** These are the standard coordinates in which the Alcubierre solution is
stated and in which Warp Factory works, so cross-code comparison (Milestone 10)
is a component-by-component comparison rather than a transformation exercise.

---

## 2. Formulation 2 — co-moving coordinates

### Metric

Substituting x′ = x − v_s t, t′ = t (so dx = dx′ + v_s dt′):

    ds² = −dt′² + (dx′ + v_s(1 − f(r′)) dt′)² + dy² + dz²

    r′ = √(x′² + y² + z²)

Components:

    g′_00 = −1 + v_s²(1 − f)²      g′_0x′ = v_s(1 − f)
    g′_x′x′ = g′_yy = g′_zz = 1

### Coordinate system

Cartesian, co-moving with the bubble centre.

### Required temporal derivatives

**None.** ∂_t′ g′_μν = 0 **exactly**, by construction rather than by
approximation — this is the honest version of what A3 pretended to have.

### Observer interpretation

**This is the drawback.** Coordinate observers now move *with* the bubble. At
spatial infinity f → 0, so

    g′_00 → −1 + v_s²,      g′_0x′ → v_s

The metric does **not** approach η. The spacetime is still asymptotically flat —
that is a geometric statement, invariant — but the *components* approach a
boosted Minkowski form. Consequences:

- the existing `metric.alcubierre.asymptotic_flatness` test would have to be
  rewritten to compare against boosted-Minkowski rather than η;
- "energy density measured by the coordinate observer" is now a different
  physical quantity from the usual one, and every energy-condition statement
  would need requalifying;
- for v_s ≥ 1 the asymptotic coordinate observer is not timelike at all
  (g′_00 = −1 + v_s² ≥ 0 as r → ∞), so the entire exterior loses coordinate
  observers — a far worse failure mode than Formulation 1's, where the problem is
  confined to the bubble wall.

### Compatibility with the existing Cartesian grid

Complete, and marginally simpler — the bubble stays centred in the grid for all
time, so a fixed domain suffices even for long evolutions.

### Implications for Christoffel and curvature

∂_0 g′ ≡ 0, so the Christoffel and Riemann code needs only spatial derivatives.
Cheapest possible path.

Curvature **invariants** (R, R_μν R^μν, Kretschmann) are identical to
Formulation 1 — they are scalars. Curvature **components** are not, since the
frames differ.

### Implications for OTI sensitivities

Equivalent in difficulty. One subtlety: the coordinate transformation itself
depends on v_s, so ∂/∂v_s at fixed x′ is a *different* derivative from ∂/∂v_s at
fixed x. Sensitivity maps from the two formulations are not directly comparable
without transforming, and conflating them would be a silent error.

### Verification opportunities

- **Level A**: ∂_t′ g′ = 0 exactly (trivially true, therefore weak).
- **Level D**: curvature invariants must match Formulation 1 — genuinely strong.
- Weaker overall: making time dependence vanish also removes the checks that
  would catch a temporal-derivative bug.

### Numerical complexity

Lowest of the three. No temporal derivatives at all.

### Risk of confusing coordinate effects with physical curvature

**Moderate, and specifically here.** The nonzero asymptotic g′_0x′ is a pure
coordinate artefact — flat space in a moving frame — and it appears in exactly
the component (g_0x) that carries the physical warp shift. A reader, or a plot,
that does not carry the frame declaration will misread a boost as a shift. This
is the classic way coordinate effects get mistaken for physics, and it is why
this formulation is not recommended as primary despite being the cheapest.

---

## 3. Formulation 3 — frozen-time spatial slice (diagnostic only)

### Metric

Formulation 1 with x_s held fixed and ∂_t g := 0 imposed.

### What it is

A spatial slice of a moving spacetime, with the temporal derivative **discarded
rather than computed**. It is not a solution of anything: setting ∂_t g = 0 in
a metric that genuinely depends on t does not produce the metric of a different,
static spacetime — it produces an inconsistent set of inputs to the Christoffel
formula.

### Legitimate uses

- inspecting metric components, determinant, signature and validity masks;
- exercising the parameter-sensitivity machinery, which is pointwise in t;
- regression-testing the spatial derivative operators;
- fast parameter sweeps for grid sizing and conditioning.

All of these are properties of the slice itself and are unaffected by ∂_t.

### Illegitimate uses

Any Christoffel symbol, any curvature tensor, any stress-energy tensor, any
energy condition, any statement about negative energy. Γ^ρ_μν built from a
truncated ∂_α g_μν is not the connection of any metric, so nothing downstream of
it means anything.

### Enforcement

Renamed `TimeTreatment::FrozenTimeDiagnostic`. The name appears in the type, in
`metric_metadata.json`, and on every gate run. **The curvature code, when it is
written, must refuse to run in this mode** rather than warn — a warning in a log
is not a control.

---

## 4. Comparison

| | 1 · inertial | 2 · co-moving | 3 · frozen diagnostic |
|---|---|---|---|
| ∂_t g_μν | ≠ 0, = −v_s ∂_x g | ≡ 0 exactly | discarded |
| Extra cost, constant v_s | **none** (identity) | none | none |
| Extra cost, variable v_s | ≤4 metric evals | n/a (breaks) | n/a |
| Asymptotic components → η | **yes** | no (boosted) | yes |
| Coordinate observer | inertial rest observer | co-moving | inertial rest |
| Observers exist for v_s ≥ 1 | fails only near wall | **fails everywhere outside** | n/a |
| Curvature invariants | correct | correct | **meaningless** |
| Warp Factory comparison | direct | needs transform | n/a |
| Verification strength | strong (§1) | moderate | none |
| Coordinate-artefact risk | low | moderate | n/a |

## 5. Recommendation and reasoning

**Adopt Formulation 1 as the primary formulation.**

The decision does not rest on ease of coding. Formulation 2 is strictly easier
to code — it needs no temporal derivative at all — and is nevertheless not
recommended. The reasons are physical and verificational:

1. **The cost argument for Formulation 2 evaporates.** §0 shows the temporal
   derivative for a constant-velocity bubble is an exact multiple of a spatial
   derivative the project already computes and has verified to second order.
   The main reason to prefer co-moving coordinates does not exist here.

2. **Observer interpretation is the point of the project.** The objective is a
   negative-energy penalty, and energy density is observer-dependent. In
   Formulation 1 the coordinate observers are asymptotically the inertial rest
   observers, which is the family every published Alcubierre energy-condition
   statement refers to. In Formulation 2 they co-move, and every such statement
   would need requalification.

3. **Formulation 2 degrades exactly where the project is headed.** The
   interesting regime is v_s ≥ 1. There, co-moving coordinate observers fail
   *everywhere in the exterior* (g′_00 → −1 + v_s² ≥ 0), whereas Formulation 1's
   breakdown is confined to the region where v_s f > 1, near the wall — a local,
   detectable, mappable failure rather than a global one.

4. **Formulation 1 verifies better.** It supplies an independent Level A check
   (analytic ∂_t vs −v_s × numerical ∂_x) that shares no code between its two
   sides. Formulation 2's corresponding check, ∂_t′ g′ = 0, is true by
   construction and therefore proves nothing about the derivative machinery.

5. **Cross-code validation is direct.** Warp Factory works in the inertial
   frame; matching components requires no coordinate transformation, removing a
   whole class of comparison error from Milestone 10.

**Formulation 2 is retained as an independent cross-check, not discarded.**
Implementing the co-moving metric costs little and buys the strongest available
verification of the curvature code: curvature *invariants* computed in the two
formulations must agree to truncation error, while their *components* differ.
That is a Level D check against an independently written metric, and it would
catch a whole class of frame and sign errors. It should be added at Milestone 5.

**Formulation 3 is retained as `FrozenTimeDiagnostic`**, restricted as in §3.

## 6. Consequent work, in order

1. Rename `TimeTreatment::FrozenTimeDiagnostic` → `FrozenTimeDiagnostic`; add
   `InertialTimeDependent`; keep `CoMoving` declared and unimplemented.
2. Implement temporal derivative operators with the same second-order interior
   and one-sided boundary stencils as the spatial ones, verified against
   manufactured analytic functions of (t,x,y,z). Needed for variable v_s and to
   verify the identity. **(Done — see `results/verification/milestone_02t/`.)**
3. Implement `AlcubierreMetric` with x_s(t) = v_s t and analytic ∂_t.
4. Verify ∂_t g = −v_s ∂_x g as a Level A gate.
5. Only then begin Christoffel symbols, with the curvature entry point refusing
   `FrozenTimeDiagnostic`.
6. At Milestone 5, add the co-moving metric and compare curvature invariants.

## 7. Open issue this decision does not resolve

**Coordinate observers fail where v_s f > 1.** Then g_00 = −1 + v_s²f² ≥ 0, ∂_t
is no longer timelike, and there is no observer at rest in these coordinates.
For the reference parameters (v_s = 0.2, f ≤ 1) this never happens. For the
superluminal regime the project actually cares about, it happens in a shell
around the bubble wall.

This is a property of the spacetime and the coordinates, not of any formulation
choice above, and it is **not** solved by switching frames — Formulation 2 makes
it worse. It must be handled at Milestone 7 by choosing an observer family that
exists there (for instance Eulerian observers normal to the constant-t slices,
which remain timelike where the slicing is spacelike), and by exporting a mask
marking where any given observer family is undefined.

The same condition reappears independently in Program B: the Plebański
metric-to-medium mapping also requires g_00 < 0, and fails in the same region.
See [effective_metric_theory.md](effective_metric_theory.md) §6.
