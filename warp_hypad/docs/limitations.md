# Missing features and current limitations

Honest state of Program A as of Milestones 1, A2, A2T and A3, with the OTIlib
backend integrated. Program B keeps its own register at
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

**215 verification cases pass, 0 fail, 0 skipped** across six gates (with the
OTIlib backend on; 2 skipped with it off).

## What does not exist

**No curvature.** There is no Christoffel symbol, no Riemann or Ricci tensor, no
Einstein tensor, no stress-energy tensor, no observer frame, no energy
condition, no objective and no optimizer.

Specifically absent:

| Milestone | Missing |
|---|---|
| ~~2~~ | *done (A2)* |
| ~~3~~ | *done (A3); time dependence now resolved, see §1* |
| 4 | Christoffel symbols, metric compatibility |
| 5 | Riemann, Ricci, curvature invariants, Bianchi identities |
| 6 | Einstein tensor, stress-energy, contracted Bianchi |
| 7 | observer frames, four-velocity normalization, energy density |
| 8 | energy conditions, smooth negative-energy objective, integration |
| 9 | HDF5/VTK export, full Python plotting package |
| 10 | *(scalar, seeding and extraction are DONE ahead of schedule — see §2; sensitivity **fields** still need the metric and curvature code from Milestones 3–6)* |
| 11+ | optimization, scaling study, Warp Factory cross-validation. *(FD verification of the AD scalar itself is done; FD verification of full-field sensitivities is not.)* |

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
