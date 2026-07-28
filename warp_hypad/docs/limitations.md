# Missing features and current limitations

Honest state of Program A as of Milestones 1, A2, A2T, A3, A4 and A5, with the
OTIlib backend integrated. Program B keeps its own register at
[../../program_b_em_inverse_design/docs/claims_and_limitations.md](../../program_b_em_inverse_design/docs/claims_and_limitations.md).

## What exists

- **Milestone 1** — Cartesian grid, contiguous storage, fixed-size tensors, 4×4
  algebra with singularity and conditioning reporting, generic-scalar seam, JSON
  config, metadata export, verification harness. **77 cases.**
- **OTIlib backend** — real hypercomplex AD (`onumm3n2`, 3 directions, order 2)
  wrapped as `OTIScalar`. **23 cases.**
- **Milestone A2** — first, second and mixed derivatives; second-order centred
  interior and one-sided boundary stencils; scalar and tensor field
  differentiation. **41 cases.**
- **Milestone A3** — Minkowski and frozen-time Alcubierre metrics, shape
  function, inverse, determinant, symmetry, Lorentzian-signature checking by
  eigenvalue, validity masks, and HyPAD sensitivities for v, R and σ. **31
  cases.**

- **Milestone A2T** — temporal derivatives at second order in three stencil
  positions, space–time mixed derivatives, and the rigid-translation identity.
  **23 cases.**
- **Program B1** — metric ↔ electromagnetic constitutive mapping. **20 cases**,
  tracked separately.

- **Milestone A5** — Einstein tensor (raw and symmetric projection kept
  separately), its trace and square, the contracted Bianchi residual, the
  cosmological-constant treatment policy, the inferred stress-energy tensor and
  the geometrized/SI units layer. **27 cases.**

**285 verification cases pass, 0 fail, 0 skipped** across eight gates (with the
OTIlib backend on; 2 skipped with it off). Full CTest: 55/55 with OTIlib, 54/54
without.

- **Milestone A4** — Christoffel symbols, Riemann and Ricci tensors, Ricci
  scalar, R<sub>μν</sub>R<sup>μν</sup> and the Kretschmann invariant; the
  coordinate-chart abstraction; the temporal-derivative policy and its two
  refusals; four closed-form reference spacetimes and an independently generated
  symbolic reference for Alcubierre. **43 cases.**

## What does not exist

**No observers, so no energy.** The Einstein tensor now exists and a
stress-energy tensor can be *inferred* from it, but there is still no observer
frame, no energy density, no energy condition, no objective and no optimizer.
Nothing in this repository may yet be called an energy: a coordinate component
of G or of T is not one, and A6 is what introduces the observer four-velocities
that would make `ρ = T_{μν}u^μu^ν` meaningful. There is therefore still nothing
to minimise.

Specifically absent:

| Milestone | Missing |
|---|---|
| ~~2~~ | *done (A2)* |
| ~~3~~ | *done (A3); time dependence now resolved, see §1* |
| ~~4~~ | *done (A4): Christoffel, Riemann, Ricci, invariants, metric compatibility, Bianchi* |
| ~~5~~ | *done (A4), which absorbed the Riemann/Ricci/invariant scope* |
| ~~6~~ | *done (A5): Einstein tensor, inferred stress-energy, contracted Bianchi, units* |
| 7 | observer frames, four-velocity normalization, energy density |
| 8 | energy conditions, smooth negative-energy objective, integration |
| 9 | HDF5/VTK export, full Python plotting package |
| 10 | *(scalar, seeding and extraction are DONE ahead of schedule — see §2; sensitivity **fields** still need the stress-energy code from Milestone 6)* |
| 11+ | optimization, scaling study, Warp Factory cross-validation. *(FD verification of the AD scalar itself is done; FD verification of full-field sensitivities is not.)* |

### A5-specific limitations

- **The contracted Bianchi identity is the weakest-converging quantity in the
  project.** It costs a *third* nested derivative pass, so at 17–33³ it carries
  more pre-asymptotic error than anything else here: measured order 1.17 against
  a design order of 2, and gated with a deliberately wide band. What the case
  establishes is that it converges at all — a non-converging Bianchi residual
  would mean the Einstein tensor is wrong.

- **`compute_einstein` with a resampling policy costs three full curvature
  computations.** The temporal slot of ∇<sub>μ</sub>G<sup>μν</sup> needs G on
  neighbouring time slices. Those are scoped so only one extra `CurvatureFields`
  is alive at a time, but the Bianchi cases still run on smaller grids than the
  rest of A5 for that reason.

- **SI values are a presentation layer, and nothing verifies them against
  physics.** Dimensional consistency and the 1/L₀² scaling are checked; whether a
  particular declared `CoordinateScale` is the *right* one for a given model is a
  modelling decision the code cannot check.

- **`inferred_stress_energy` has no mechanism behind it.** It is what Einstein's
  equation says would have to be present to source the chosen geometry. Whether
  such a source can exist is not addressed, and A5 makes no claim either way.

### A4-specific limitations

- **Curvature is `double` only.** The pipeline is templated on Scalar and would
  instantiate over `OTIScalar`, but nothing does. A `Rank4Field` is 256 scalars
  per grid point — 141 MB at 41³ in plain double — and the AD type multiplies
  that by the channel count. Curvature *sensitivities* are deliberately deferred
  until after the curvature itself is verified.

- **Memory sets the grid size, and the grid size sets the accuracy.** The
  Alcubierre reference comparison runs at 33³ with a peak working set near
  350 MB. That resolves the σ = 1 wall with about four points, which is why its
  measured convergence order is 1.63 rather than the 1.98–2.08 the smooth
  reference geometries reach. The order is real but the finest grid is only just
  entering the asymptotic regime.

- **Second-order accuracy, twice nested.** Curvature applies the derivative
  operators twice, so the boundary shell is where a figure goes wrong first.
  Every A4 residual is reported interior and boundary separately, with a shell
  depth of 2, and the boundary figure is typically an order of magnitude worse.

- **No second Bianchi identity.** ∇<sub>[λ</sub>R<sup>ρ</sup><sub>|σ|μν]</sub> = 0
  needs a covariant derivative of a rank-4 field, which needs a third derivative
  pass. Not implemented; the contracted form ∇<sub>μ</sub>G<sup>μν</sup> = 0 is
  Milestone 6's obligation and is the more useful of the two.

- **Time is still not gridded.** Temporal derivatives come from re-evaluating the
  metric at shifted times. `AnalyticResample` costs five metric-field evaluations
  and three Christoffel assemblies per curvature call. That is exact-in-principle
  and assumption-free, but it means the pipeline cannot consume a metric that
  exists only as tabulated data.

- **No co-moving cross-check.** `TimeTreatment::CoMoving` remains declared and
  unimplemented, so the planned "invariants agree while components differ" check
  between inertial and co-moving charts has not been run. The
  Schwarzschild-versus-isotropic case now provides an equivalent statement for a
  static geometry, but not for a moving one.

## Known open questions

### 1. Time dependence — RESOLVED

`TimeTreatment::InertialTimeDependent` with x_s(t) = v_s t is the chosen
formulation. See [time_dependence_decision.md](time_dependence_decision.md).

It rests on an exact identity: for a rigidly translating bubble every component
depends on t and x only through u = x − v_s t, so **∂ₜg = −v_s ∂ₓg**. The
temporal derivative therefore costs no more than a spatial one, which removes
the usual argument for co-moving coordinates. Verified at Level D: the residual
against an independently computed numerical ∂ₓ converges at order 1.82 toward
zero.

`FrozenTimeDiagnostic` is retained for metric inspection, validity masks and
parameter sensitivities. **`valid_for_curvature()` returns false for it**, and
the curvature code must refuse rather than warn.

Still open: adding the co-moving formulation as a curvature-invariant
cross-check at Milestone 5 (recommended, not required).

### 1b. Sensitivity is at fixed bubble trajectory

`center_at()` advances the bubble with `real_part(velocity)`, so OTI computes
∂/∂v **at fixed trajectory** — shape sensitivity, not the total lab-frame
derivative which would add a d(x_s)/dv = t term. The two differ substantially at
t ≠ 0.

This vanishes at t = 0, and the objective for a rigidly moving bubble is
invariant under time translation, so all optimization can be done at t = 0 where
the two coincide. Revisit only if sensitivities at fixed co-moving position are
wanted.

### 2. AD backend — RESOLVED (OTIlib integrated)

Superseded. A machine-wide search located the authentic OTI sources and led to
upstream `mauriaristi/otilib`, pinned at commit `de3f7fa2`. The real
`onumm3n2_t` type (3 imaginary directions, truncation order 2) is now wrapped as
`OTIScalar` and verified — 23/23 cases, including
`d(A⁻¹)/dp == −A⁻¹(dA/dp)A⁻¹` at 1.9e-16.

No fictional `HyPADScalar` was ever created. `TestScalar` remains what it always
was: a compile-time proof obligation with zero derivative channels.

Second-order propagation **is** supported, so Hessians are available — better
than the order-1-only Fortran variants present locally.

See [otilib_discovery_report.md](otilib_discovery_report.md).

**Licence: RESOLVED.** GPL-3.0 adopted for the research prototype. The backend
stays optional (`-DWARP_HYPAD_WITH_OTILIB=ON`) and isolated behind the adapter,
and OTIlib is a pinned external checkout, never vendored into this tree.

### 3. The project specification is truncated

The source specification cuts off mid-sentence in Milestone 10's verification
list ("Compare every major field: metric; inverse metric; …"). The
verification-gated milestone plan therefore has **no** stated requirements for
gradient verification against finite differences, gradient-based optimization,
the performance study, or Warp Factory cross-validation — all of which existed
in the earlier, non-gated revision as Milestones 8–10. Those requirements need
reconstructing before work reaches them.

## Limitations of what does exist

### Matrix inversion uses the adjugate, not LU

Chosen because pivoting requires comparing scalar magnitudes, which a
hypercomplex AD type cannot support. The cost is roughly 2× the flops of LU and
worse behaviour on ill-conditioned input.

This is mitigated, not eliminated: `InverseResult` reports `ok`, `cond1` and
`ill_conditioned`, and callers must check them. The measured relationship is
residual ≈ κ₁ × 6e-17, and the ill-conditioning threshold is derived from that
(κ₁ > 1e-8/ε ≈ 4.5e7) rather than picked as a round number.

**Callers must check `.ok`.** Nothing in the type system enforces it yet.

### Singularity detection is determinant-magnitude based

`|det| ≤ 1e-13 · ‖A‖₁⁴` — scale-invariant, but a determinant test is a weaker
singularity criterion than a rank-revealing decomposition. It rejects exactly
singular matrices reliably (verified, 200/200) and rejects near-singular ones
aggressively; a matrix that is ill-conditioned but not small-determinant could
in principle slip through, which is what the separate `cond1` diagnostic is for.

### Signature checking is diagnostic-only, in double

`analyse_signature` uses cyclic Jacobi rotation on real parts. Signature is a
count of signs, so it has no useful derivative and is deliberately not
differentiable. It is a validity diagnostic, never an input to physics.

### Alcubierre parameter set is fixed at three

`onumm3n2` carries exactly three derivative directions, matching v, R and σ. The
bubble centre is a fixed coordinate, not a differentiated parameter; seeding it
would need a fourth direction and a different OTI type.

### Grid spacing is uniform only to a few ulp

A deliberate consequence of pinning endpoints and midpoint exactly. See
[memory_layout.md](memory_layout.md).

### No unit-testing framework

The project specification suggests Catch2 or GoogleTest. Neither is used.

The reason is that neither produces what the verification policy actually
requires — `verification_summary.json`, `verification_table.csv` and
`verification_report.md`, each carrying reference, tolerance, errors, runtime,
commit and compiler per case. That machinery had to be written regardless, and
once written, a second assertion framework layered on top would add a
dependency without adding coverage. Each verification group is registered
individually with CTest, so per-area reporting is preserved.

If a framework is wanted later, the check functions are already independent and
would wrap cleanly.

### Provenance is currently incomplete

The git repository has **no commits**, so `metadata.json` and every verification
record show `git_commit: "unknown"` and `git_dirty: true`. The gate prints a
warning about this rather than hiding it. Committing the tree will resolve it;
the CMake capture is already in place.

### No performance work

Nothing is optimized, parallelized or vectorized, by design — correctness first.
There is no benchmark harness; that is Milestone 9.

### A2/A3 specific limitations

- **Mixed second derivatives of the same axis.** `mixed_derivative(f, a, a)`
  composes two centred first derivatives, giving the wider
  (f[i+2] − 2f[i] + f[i−2])/(4h²) stencil — second order, but with a four times
  larger error constant than `second_derivative`. Callers wanting ∂²/∂x² must
  use the latter. Not currently enforced by the type system.

- **Convergence order is not asserted where the stencil is exact.** For fields
  whose truncated derivative vanishes identically (most operators on the test
  polynomial), the measured error is pure roundoff and *rises* under refinement.
  Those six cases assert exactness instead. This is correct, but it means the
  polynomial field exercises fewer convergence paths than it appears to.

- **OTI real parts are not bit-identical to double when division is present.**
  OTIlib implements division as pow(−1) composed with multiply, so the shape
  function's real part rounds differently from a hardware divide — measured
  2.8e-17, about 0.12 ulp. Division-free expressions *are* bit-exact.

- **Finite differences cannot verify small metric sensitivities.** At the
  exterior sample ∂g₀₀/∂σ ≈ 1.8e-11 against g₀₀ ≈ −1, so FD loses eleven orders
  to cancellation and bottoms out near 1e-3 whatever the step. That is a limit
  of finite differences, not of HyPAD, and is the reason the project uses AD.

- **No 3D or VTK output yet.** A2/A3 export CSV slices and profiles only. HDF5
  and VTK remain Milestone 9.
