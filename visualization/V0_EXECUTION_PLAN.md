# Phase 0 — execution plan

Written from the committed CSV schemas, not from memory. Every "have" row below
names the file and the columns that carry the product; every "needs exporter" row
names what is missing and why no existing record substitutes.

**Repository state at planning time**

    HEAD          c1262d4  V0 foundation
    working tree  clean
    tags          7, all present
    CTest         33/33 with OTIlib

**The rule this plan enforces:** *do not add a second implementation path for a
quantity that is already saved.* Where a column exists, the renderer reads it.
Where no column exists, a C++ exporter calls the same production code the gates
call — and the Python renderer still computes nothing.

---

## Summary

| front | products | from existing records | need a new exporter | blocked |
|---|---|---|---|---|
| A static | 11 | 9 (done) | 0 | 2 (curvature → A4) |
| A animations | 9 | 0 | 9 | 0 |
| B1/B1C | 13 | 6 | 7 | 0 |
| B2 | 8 | 5 | 3 | 0 |
| B3 | 12 | 6 | 6 | 0 |
| B4/B4C | 15 | 5 | 10 | 0 |
| cross-program | 3 | 1 | 2 | 0 |

**37 of 71 products need new production exporters.** That is the real shape of
this milestone: it is majority *scientific output-writing*, not plotting. The
existing gates saved what they needed to assert, which is a narrower set than
what an animation needs — a gate checks a value, an animation needs the field.

---

## Program A

### Static — 9 complete, 2 blocked

Delivered in `c1262d4` from `milestone_03/{shape_profiles, metric_slice_xy,
sensitivity_slice_xy, hypad_vs_fd}.csv`. 28 diagnostics, 0 fired.

`ricci_scalar_slice`, `curvature_norm_slice` — **BLOCKED_BY_UNIMPLEMENTED_SCIENCE**,
pending A4. No Christoffel, Riemann, Ricci or Kretschmann exists in `warp_hypad`.

### Static — refinements required by Phase 1.1

| product | source | note |
|---|---|---|
| deviation views `Δg00 = g00+1`, `Δg11 = g11−1`, `Δdet = det+1` | `metric_slice_xy.csv` | have columns; must be titled **contrast-enhanced deviation from stated reference** and never shown as the absolute field |
| `det(g)` status panel | `metric_slice_xy.csv` | analytic −1, measured min/max, worst deviation, tolerance, pass/fail — replaces reliance on a flat heatmap |
| dynamic ergosurface title | `metric_slice_xy.csv` | already detects `g00` never crossing zero; retitle to `g00 diagnostic — no ergosurface for this parameter set` |

### Animations — all 9 need the exporter

**No swept field data exists anywhere in the repository.** `metric_slice_xy.csv`
is a single snapshot at one parameter set. Every A animation varies time,
`v_s`, `R`, `sigma` or the observer convention.

New target `a_sweep_export`, calling the production `MetricField` path:

    translation   frame,time,bubble_center_x,x,y,f,g00,g01,g11,det,signature_ok
    velocity      frame,v,x,y,f,g00,g01,det,signature_ok
    radius        frame,R,x,y,f,g00,g01,det,signature_ok
    steepness     frame,sigma,x,y,f,g00,g01,det,signature_ok
    sensitivity   frame,<param>,x,y,dg00_dv,dg00_dR,dg00_dsigma,dg0x_*

Each file records coordinate mode, fixed-lab vs fixed-centre convention, metric
parameters, grid, units, time, git commit and an output checksum.

The velocity sweep must span `v·f = 1` where the domain permits, because that is
the ergosurface condition and the current parameter set does not reach it — the
committed slice tops out at `g00 = −0.96`.

---

## Program B1 / B1C

### Have

| product | source | columns |
|---|---|---|
| `epsilon_eigenvalues` | `b1c/ergosurface_conditioning.csv` | `eps_eig_min, eps_eig_max, eps_cond, eps_det` along `x` |
| `mu_inverse_eigenvalues` | same | `mu_inv_cond, mu_inv_*` |
| `ergosurface_constitutive_behavior` | same | `g00, f, vf, covariant_valid, effective_observer_admissible` |
| `zeta_norm` | `b1c/corrected_status_map.csv` | `zeta_max` over `(v_s, x)` |
| `constitutive_vs_velocity` (animation) | same | the `v_s × x` grid is a ready-made sweep |
| `ergosurface_crossing` (animation) | same | frames indexed by `v_s`, with `g00` and the status flags |

### Need an exporter

| product | why nothing substitutes |
|---|---|
| `zeta_adjoint_norm` | only `zeta_max` is saved; ζ and ζ† are distinct objects and the project spent a whole milestone establishing that `reconstruct_Z_A11` needs **`zeta_adjoint`**, so conflating them here would undo that |
| `observer_block_comparison`, `observer_blocks_same_Z` | **no observer sweep exists.** This is the central B1C claim — `blocks(V) ≠ blocks(U)` while `Ẑ[V] = Ẑ[U] = Z` — and it needs per-observer blocks plus both reconstruction residuals in one record |
| `Z_roundtrip_residual` | the residual is asserted by the gate as a scalar; the *field* is not saved |
| `constitutive_vs_radius`, `constitutive_vs_steepness` | no `R` or `sigma` sweep |

New target `b1c_sweep_export`. The observer record must carry the laboratory
metric, the observer four-velocity, the effective metric, all four blocks per
observer, the block-difference norm and each reconstruction residual — so the
renderer can show difference and invariance in one frame without computing
either. **The laboratory observer is normalised with the laboratory metric**, per
the B1C type separation; the exporter must not normalise it with `g_eff`.

---

## Program B2

### Have

`objective_history`, `gradient_norm_history` (`objective_histories.csv`),
`parameter_recovery` (`recovered_designs.csv`), `Z_residual_matrix`
(`z_residual_by_pair.csv`), `restricted_family_comparison`
(`family_residuals.csv`).

### Need an exporter

`target_vs_recovered_blocks`, `constraint_activity`, and the
`inverse_design_optimization` animation. The saved history has
`objective, gradient_norm, step, backtracks` per iteration but **no tensor
state** — no ε eigenvalues, no μ⁻¹ eigenvalues, no |ζ|, no per-iteration
residual, no constraint status. The animation's whole purpose is to show whether
a falling objective corresponds to sensible tensor convergence, and that cannot
be drawn from a scalar history.

New target `b2_history_export`, writing per-iteration tensor state. It also
enables the required diagnostic: **objective decreasing while a physical block
residual increases catastrophically.**

---

## Program B3

### Have

`epsilon_frequency_response`, `mu_frequency_response`, `zeta_frequency_response`,
`passivity_margin`, `constitutive_frequency_sweep`,
`passivity_eigenvalue_sweep` — all from `b3_passivity_sweep.csv`
(`omega, re/im eps, re/im mu, re/im zeta, min_loss_eigenvalue`), which is a
complete frequency sweep. `target_vs_model_band` from
`b3_bandwidth_tradeoff.csv`.

### Need an exporter

`pole_map` (pole locations are proved structurally, never written out),
`kk_diagnostic` (the Kramers-Kronig residual is a gate scalar, not a curve),
`hard_vs_soft_bounds`, `resonance_damping_sweep`, `one_pole_vs_two_pole`.

The pole map matters most: causality in B3 is established **by pole placement**,
and a figure asserting poles lie in the correct half-plane must read actual pole
coordinates rather than infer them from a sampled response.

---

## Program B4 / B4C

### Have

The three `*_dispersion_surface` figures from `b4_dispersion_surfaces.csv`
(`medium, theta, kappa_branch1, kappa_branch2, branches`), `fresnel_roots`
partially from the same, `exact_vs_restricted_propagation` from
`b4c_sensitivities.csv`, and the axion assertions from `b4c_classification.csv`
(`contains_axion, fresnel_single_cone, axion_coefficient` — enough to assert the
cone is unchanged while the tensor differs).

### Need an exporter

`metric_fresnel_factorization`, `group_phase_energy_velocity`,
`layer_field_profile` (only ε, μ⁻¹, ζ against `x` are saved, not the field),
`reflection_transmission_absorption`, and five of the six animations —
`dispersion_frequency_sweep` (only θ is swept, not ω), `fresnel_root_evolution`,
`layered_wave_propagation`, `resonance_transmission`, `exact_vs_zero_zeta_wave`.

The repeated-root requirement needs care in the exporter, not the renderer: a
metric-induced quartic has a **repeated** quadratic factor, and the export must
record root multiplicity explicitly so the renderer can draw a double root as a
double root instead of two coincident branches that framewise sorting may swap.

---

## Cross-program

`catalogue_to_relativistic_journey` — **have**, entirely from committed C1
records, no recomputation.

`metric_to_constitutive` and `constitutive_to_propagation` — need the A and B1C
sweep exporters above, then join on the shared parameter.

---

## Ordering

1. `a_sweep_export` → the 9 A animations. Largest single unlock, and it proves
   the exporter pattern end to end.
2. `b1c_sweep_export` → the 7 B1 products, including the observer animation that
   carries B1C's central claim.
3. `b2_history_export`, `b3_pole_export`, `b4_field_export` → the remaining B
   products.
4. The two cross-program chains, which fall out of 1 and 2.
5. Gallery, then the `V0-CROSS-PROGRAM-VISUALS` gate with its three-state
   result: `PASS`, `FAIL`, `BLOCKED_BY_UNIMPLEMENTED_SCIENCE`.
6. Tag `v0-cross-program-visuals-complete`, **without claiming curvature
   coverage.**

Then A4, then C2A, then C2 — each with its own gate, commit and tag, and the
full suite run between phases.

## One scope observation, recorded now rather than at the end

The A4 specification asks for a de Sitter check at `R = 12H²`. That is the
correct value for **four-dimensional** de Sitter with the convention where `R > 0`
for a sphere. Since the phase requires "0 unresolved sign-convention cases", the
convention will be fixed and written down before any curvature code is compiled —
`R^a_{bcd} = ∂_c Γ^a_{bd} − ∂_d Γ^a_{bc} + …` and `R_{ab} = R^c_{acb}` — and each
reference value stated against it, so a sign disagreement surfaces as a failing
case rather than as a silently flipped figure.
