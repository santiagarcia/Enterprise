# Technical design

## Goal

An independent C++ implementation of the numerical physics pipeline

    p → g_μν(x; p) → ∂_α g_μν → Γ^ρ_μν → R_μν → R → G_μν → T_μν → J(p)

with hypercomplex automatic differentiation designed in from the start, so that
parameter sensitivities ∂J/∂p come out of the same code that computes J — no
separate sensitivity equations, no adjoint derivation, no divergence between the
two.

This is not a port of Warp Factory. The equations are implemented independently
from their mathematical statements and documented in
[conventions.md](conventions.md); cross-validation against Warp Factory is a
late milestone, not a development method.

## The one architectural decision

Everything follows from templating every kernel on `Scalar`:

```cpp
template <typename Scalar> ... compute_stress_energy(...);
```

so that `compute_stress_energy<double>` and `compute_stress_energy<HyPADScalar>`
are the same code. That only works if kernels never assume anything about the
scalar beyond arithmetic. Three consequences shaped the Milestone 1 code:

### 1. No branching on a scalar value

A comparison on a hypercomplex number either is ill-defined or silently discards
the derivative channels; a branch introduces a derivative discontinuity the AD
type cannot see.

So matrix inversion uses the **closed-form adjugate**, not LU with partial
pivoting — pivoting selects a row by comparing magnitudes. The adjugate is
branch-free straight-line arithmetic that differentiates cleanly.

Where an ordering decision is genuinely required — is this matrix singular? does
this residual pass? — the code goes through `real_part()` and decides on the
real component **explicitly**. Whether a matrix is invertible is a property of
the primal value, not of its sensitivities.

### 2. The seam is one small trait

`ScalarTraits<Scalar>` provides `real()`, `from_double()` and two constants.
Nothing else in the library knows anything about the scalar type. Milestone 10
specialises this trait and nothing else changes.

`HyPADScalar` was **not invented**. No HyPAD source was available, and a
fabricated interface connected to nothing would be worse than none.

### 3. Genericity is proven at compile time, not asserted

`TestScalar` deliberately provides no ordering, no equality, and no implicit
conversion to or from `double`. If any kernel compared two scalars or named
`double` internally, the `TestScalar` instantiation would fail to build.

`static_assert`s pin those absences, so the proof cannot quietly weaken if
`TestScalar` later grows a conversion operator.

The runtime check is that `double` and `TestScalar` agree **bit for bit** — a
stronger statement than a tolerance, and the one that would catch a kernel
treating the two types differently. It holds, which is why the build pins
`/fp:precise` and `-ffp-contract=off`.

## No third-party dependencies

Not minimalism for its own sake: a vendored linear-algebra library that could
not be instantiated with `HyPADScalar` would have to be removed again in
Milestone 10, and the specification explicitly conditions Eigen on custom-scalar
support. The 4×4 operations needed are small enough to own.

JSON is hand-written for one specific reason: numbers must round-trip
bit-exactly, which off-the-shelf writers at fixed precision do not guarantee.
The writer emits the shortest representation that reparses to the identical
double, so `0.2` stays `"0.2"` and a config recovered from a results directory
reproduces the run exactly.

## Verification architecture

A milestone is complete when a **process returns zero**, not when someone reads
output and judges it plausible.

`milestone1_verify` exits 0 only if every case passes. Each case records name,
mathematical reference, verification level, grid, parameters, expected, computed,
absolute error, relative error, tolerance, verdict, runtime, git commit,
compiler and build mode. Three artefacts are written per milestone:
`verification_summary.json`, `verification_table.csv`, `verification_report.md`.

Levels, in decreasing strength: **A** exact algebraic identity · **B** analytic
solution · **C** manufactured solution · **D** independent numerical method ·
**E** grid convergence · **F** cross-code validation.

**No implementation verifies itself.** `tests/reference/ref_linalg.hpp` exists to
disagree with `mat4.hpp` and shares no algorithm with it: Gauss-Jordan with
partial pivoting against the adjugate, recursive Laplace expansion against
paired 2×2 minors. Agreement to 1e-12 across 500 random matrices per family then
means something. Had it been a refactor of the code under test, it would mean
nothing.

Two design habits fall out of this:

**Test cases must be asymmetric.** The hand-calculated multiplication check uses
`diag(1,2,3,4)·B` and `B·diag(1,2,3,4)`, which differ — so it detects a
transposed row/column convention. An identity or symmetric test case could not.

**Detectors are themselves tested.** `contract.symmetry_detector_works` feeds a
deliberately antisymmetric tensor to the symmetry residual function and requires
it to return exactly the planted magnitude. A residual function that always
returned zero would otherwise make every later symmetry gate pass vacuously.

That habit already paid: the near-singular check originally reported "200
rejected, 0 flagged" — the determinant test caught every case, so the
`ill_conditioned` flag was never exercised and the check passed without proving
anything. Adding a sweep across perturbation scales exposed that the flag
threshold was wrong by four orders of magnitude. See below.

## A defect the gate caught

The ill-conditioning threshold was initially `1e12`, a round number chosen by
analogy with "fewer than four digits survive". Sweeping the perturbation from
1e-4 to 1e-14 produced matrices that were accepted, **not** flagged, and had
identity residuals up to 3e-5.

The measured relationship is residual ≈ κ₁ × 6e-17 ≈ 0.3 κ₁ε, so κ₁ε is a
conservative bound. The constant is now **derived** from an accuracy target:

```cpp
inline constexpr double kInverseAccuracyTarget = 1e-8;
inline constexpr double kConditionWarnLimit =
    kInverseAccuracyTarget / std::numeric_limits<double>::epsilon();  // ≈ 4.5e7
```

so "flagged" means exactly "the residual may exceed 1e-8". After the fix the
worst *unflagged* residual across the sweep is 2.5e-9 — just under target, which
is the intended behaviour.

## Layout

```
warp_hypad/
├── include/warp_hypad/core/   scalar_traits, mat4, tensor, grid, field,
│                              json, config, metadata, verification, test_scalar
├── src/core/                  grid, json, config, metadata, verification
├── tests/
│   ├── reference/             independent linalg — must share no algorithm
│   └── *_checks.cpp           one file per verification group
├── python/                    figure generation from saved output only
├── config/                    example configurations
└── docs/
```

Directories for `geometry/`, `metrics/`, `physics/`, `sensitivity/` and
`optimization/` are **not** created yet. Empty scaffolding for unwritten code
implies structure that has not been designed.

## Build

C++17, CMake, no dependencies. Warnings are errors (`/W4 /WX`, or
`-Wall -Wextra -Wpedantic -Werror -Wshadow -Wconversion`) — a gate that passes
while the compiler reports a signed/unsigned mismatch in an indexing expression
is not a gate.

Provenance is captured at configure time and compiled in, so a verification
record cannot drift from the revision that produced it.

## Next step

Milestone 2 (finite-difference operators) **must not begin** until the time-
dependence question in [conventions.md §11](conventions.md) is resolved, because
it determines whether a time axis is needed at all.
