# V0 — cross-program scientific visual verification

A repository-wide visualization layer for Programs A, B and C, built on the
contract C1 established and now made global.

    production solver → saved CSV/JSON/frame records → visualization scripts
                      → PNG/SVG/GIF/MP4

**No visualization script recomputes physics.** Renderers consume saved numerical
outputs only, and production solver functions are never imported into a plotting
script. In C1 this is enforced mechanically — a gate case parses the AST of every
script and fails on a call to any physics routine — and V0 extends that check
across the whole tree.

---

## THE TWO FINDINGS THAT SHAPE THIS MILESTONE

An inventory of every committed numerical output was taken before a single
renderer was written. It produced two results that change the V0-A scope, and
both are stated here rather than discovered halfway through.

### 1. Most of what V0 needs is already saved

Program A alone already commits, with verified provenance:

| file | columns | supports |
|---|---|---|
| `milestone_03/shape_profiles.csv` | `R, sigma, r, f` | `alcubierre_shape_function` |
| `milestone_03/metric_slice_xy.csv` | `x, y, g00, g01, g11, det, inverse_residual, signature_ok` | `metric_components_slice`, `metric_determinant`, `metric_signature_map`, `ergosurface` |
| `milestone_03/sensitivity_slice_xy.csv` | `x, y, dg00_dv, dg00_dR, dg00_dsigma, dg0x_dv, dg0x_dR, dg0x_dsigma` | `sensitivity_v_slice`, `sensitivity_R_slice`, `sensitivity_sigma_slice` |
| `milestone_03/hypad_vs_fd.csv` | `parameter, step_size, hypad, finite_difference, absolute_error` | `oti_vs_fd_residual` |
| `milestone_03/centerline_x.csv` | `x, f, g00, g01` | centreline overlays |
| `milestone_02/error_map_xy.csv`, `convergence.csv` | operator errors and observed orders | operator diagnostics |

Nine of the eleven required V0-A static figures are therefore renderable from
data that already passed a numerical gate. **No new writer is needed for them,
and none will be added** — inventing a second path to a number that is already
saved is how two versions of "the truth" appear.

### 2. Program A HAS NO CURVATURE CODE — two required figures are blocked

`ricci_scalar_slice.png`, `curvature_norm_slice.png` and
`curvature_wall_evolution.gif/mp4` cannot be produced, because **nothing in
Program A computes curvature.** The word appears in the codebase only as a
*guard against* computing it. From `metrics/metric.hpp`:

> Legitimate: metric components, determinant, signature, validity masks,
> parameter sensitivities, regression-testing the spatial operators.
> **ILLEGITIMATE: any Christoffel symbol, curvature, stress-energy or energy
> condition.** The curvature code must REFUSE this mode, not warn about it.

There is no Christoffel routine, no Riemann tensor, no Ricci scalar and no
Kretschmann invariant anywhere in `warp_hypad`.

**What this milestone will not do about it.** It will not compute a Ricci scalar
inside a plotting script. That is precisely the failure the visualization
contract exists to prevent: a curvature field invented in Python, rendered as a
figure, and thereafter believed. A curvature figure whose numbers never passed a
gate is worse than no curvature figure, because it looks identical to one that
did.

**What is required instead.** Curvature is a *scientific* milestone for Program A,
not a visualization task. It needs Christoffel symbols from the verified
finite-difference operators, a Riemann and Ricci assembly, and verification
against cases with known answers — flat space giving exactly zero, Schwarzschild
giving `R = 0` with a nonzero Kretschmann invariant, and a convergence study
showing the expected order. Only then does a curvature figure mean anything.

Until that exists, the three curvature products are recorded as **BLOCKED —
pending verified curvature in Program A**, and the V0 gate lists them as blocked
rather than missing, so the distinction stays visible.

### 3. Animations need swept fields, which are not saved

Every V0-A animation varies something — time, radius, steepness, observer — and
the committed slices are single-parameter snapshots. This is the one place a
**new output writer is legitimate**, because the fields do not exist yet in any
form.

The writer is a dedicated export target that calls the same verified
`MetricField` production code the gates use. It computes nothing new: it
evaluates the existing metric over a parameter sweep and writes it. Its outputs
carry the same provenance block as every gate record.

---

## Structure

    visualization/
        common/        readers, provenance, render helpers, anomaly reporting
        themes/        scientific, scientific_print, cinematic
        program_a/     metric, shape function, sensitivity, signature
        program_b1/    covariant constitutive
        program_b2/    inverse design
        program_b3/    dispersion and admissibility
        program_b4/    propagation
        program_c/     re-exports the C1 renderers, which already comply
        regression/    the V0 gate and the anomaly-report checker
        galleries/     generated markdown indices

`common/` and `themes/` are lifted from the C1 implementation, which already
proved them: `records.py` (readers that refuse NaN at the boundary, provenance
validation, a closed colour-model vocabulary), `render.py` (deterministic frames,
GIF/MP4 with graceful encoder fallback, contact sheets, manifests) and
`themes.py`.

## Global metadata contract

Every visualization writes a sidecar containing:

    program, milestone, git_commit, git_dirty, input_files, input_checksums,
    simulation_parameters, coordinate_frame, observer, units, time_convention,
    render_theme, color_model, camera_model, frame_count, frame_rate,
    generated_timestamp

`color_model` is drawn from the closed set established in C1 —
`physical_spectral_model`, `approximate_blackbody_visualization`,
`catalogue_photometric_proxy`, `illustrative_only` — and the writer **refuses** an
undeclared or invented value. `input_checksums` extends the C1 sidecar: V0
renderers hash their inputs, so a figure can be traced to the exact bytes it was
drawn from.

## Visual anomaly report

Every renderer computes diagnostics on the data it draws and writes
`visual_anomaly_report.json`. **These are tests, not aesthetic scores**, and the
V0 gate reads them:

    unexpected_discontinuity, symmetry_violation, axis_inversion,
    center_trajectory_mismatch, unexpected_negative_determinant,
    non_lorentzian_signature, blank_field, clipped_values,
    frame_to_frame_jump, branch_reordering, unstable_exposure,
    lost_object_identifiers, changing_object_count, bound_violation

Each entry records the diagnostic, the measured value, the threshold and whether
it fired. A renderer that finds an anomaly still produces its figure — suppressing
the picture would hide the evidence — but the gate fails on it.

## Logical visual assertions

The gate does not rely on file existence. Each program contributes assertions
that tie the picture to the numbers:

**Program A** — the bubble centre in each frame matches the saved trajectory; a
radius sweep changes radius without translating the centre; a steepness sweep
changes wall width without moving the nominal radius; the determinant stays
consistent with Lorentzian signature everywhere the saved `signature_ok` flag is
set; the OTI-versus-FD residual map agrees with the numerical gate's reported
errors.

**Program B1** — `blocks(V) ≠ blocks(U)` while
`reconstruct_Z(blocks(V)) = reconstruct_Z(blocks(U))`, both shown in the same
frame, because the point is that observer-dependent blocks reconstruct one
observer-independent `Z`.

**Program B3** — poles stay in the half-plane causality requires; adaptive
passivity minima appear at the same frequencies the gate reports; a
hard-feasible run never displays a bound violation, and a soft-penalty result
labels one visibly.

**Program B4** — the axion animation shows a constitutive difference and a
ζ-block difference with an *identical* Fresnel cone; the repeated-root animation
shows the metric double root beside the split anisotropic roots.

## Gate `V0-CROSS-PROGRAM-VISUALS`

Files exist; PNG dimensions valid; SVGs parse; GIF/MP4 frame counts match their
manifests; contact sheets carry the expected frames; no frame blank; no saved
field non-finite; axes and units stated; equal spatial scaling where required;
source row counts match rendered object counts; camera paths deterministic;
timestamps monotonic; **no visual script calls production physics**; metadata
matches source checksums; re-rendering yields stable metadata; and the logical
assertions above agree with the numerical gates.

Blocked products are reported as blocked, with their reason, and are neither
counted as passes nor as silent omissions.

---

## Relationship to C2

V0 and C2 proceed together but do not depend on each other. C2 adds the
effective-metric voyage chain and its own outputs; V0 renders whatever is saved,
including C2's outputs once they exist. Neither weakens or bypasses an existing
gate: every current gate still runs, unchanged, and the full suite stands at
**33/33 with OTIlib and 32/32 without**.
