# Program A gallery index

The chain, in the order the milestones built it:

```
g_{mu nu}  ->  curvature  ->  G_{mu nu}  ->  T_{a-hat b-hat}  ->  energy conditions
   A3/V0.1A       A4             A5              A5/A6                  A6
```

Each arrow is a milestone with its own numerical gate, and each gallery below is
verified against the records that gate wrote. Nothing in a gallery is drawn by a
script that computes physics; the visual gates check that with an AST scan.

| stage | gallery | scientific gate | visual gate |
| --- | --- | --- | --- |
| metric components, determinant, signature, sensitivities | [program_a/](.) | `milestone3.gate` | `V0-PROGRAM-A-ANIMATIONS` |
| covariant curvature | [program_a4/](../program_a4/) | `milestone4.gate` | `V0-PROGRAM-A-ANIMATIONS` (A4 section) |
| Einstein tensor and the inferred source | [program_a5/](../program_a5/) | `milestone5.gate` | — |
| observer frames and energy conditions | [program_a6/](../program_a6/) | `milestone6.gate` | `V0-PROGRAM-A6-ENERGY-VISUALS` |

## What each stage may and may not be called

**A3 / V0.1A — the metric.** Components, determinant, signature masks and
parameter sensitivities. `g_00 = 0` first appears here as a coordinate
condition. It is not a metric degeneracy: `det(g)` stays negative through it.

**A4 — curvature.** Christoffel symbols, Riemann, Ricci and the invariants,
verified against flat space, Schwarzschild and de Sitter. The animations carry
the crossing through `v_s = 1` and show curvature remaining finite there, which
is why the crossing is an ergosurface-type condition and not a singularity.

**A5 — the Einstein tensor.** `G_{mu nu}` and the stress-energy tensor inferred
through Einstein's equation. Every A5 panel is a COORDINATE-BASIS component and
none of them is an energy density; the A5 frames say so explicitly, because that
is the figure most likely to be screenshotted and mislabelled.

**A6 — observers and energy conditions.** The first stage at which a quantity
may be called an energy density, because it is the first with a foliation, an
orthonormal tetrad and explicitly typed, admissibility-checked observers.

  * `rho = T_{a-hat b-hat} n^{a-hat} n^{b-hat}` for the Eulerian observer built
    from the ADM normal, which stays timelike through `g_00 = 0`.
  * The coordinate-stationary observer is a SEPARATE family, admissible only
    where `g_00 < 0`. At `g_00 = 0` the vector `d/dt` is null and beyond it is
    spacelike; the records leave its energy density empty and no figure draws
    it there.
  * Every value carries a certification status. A negative measured value is
    not a certified negative value, and the figures keep the two apart.

## The frame convention, unchanged since V0.1A

```
exact_onset_frame            = 8    (v_s = 1 exactly, g_00(0) = 0 exactly)
first_spatial_contour_frame  = 9    (the first frame whose grid brackets zero)
```

At onset the field only touches zero, so no contour is supported; the next frame
is the first that can carry one. A4, A5 and A6 all use the same velocity
sequence `v_s = k/10, k = 2..16`, so their animations are comparable frame for
frame.

## Standing scope

None of these galleries demonstrates that a warp source can be generated, that
quantum energy inequalities are satisfied, that the geometry is dynamically
stable, or that operational FTL travel is possible. They show what the selected
metric implies under the field equations, and what a named observer would
measure.

## Promotion

Galleries are promoted explicitly and never by a test:

```
cmake --build build --target promote_a_sweeps
cmake --build build --target promote_a4_curvature
cmake --build build --target promote_a5_einstein
cmake --build build --target promote_a6_energy
```

Tests read and write only build-local directories. Promotion decides which
verified build defines the committed record, and a test that made that decision
silently would let whichever suite ran last rewrite the repository.
