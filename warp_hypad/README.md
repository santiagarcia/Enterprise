# warp_hypad

An independent C++ implementation of the numerical physics pipeline used for
warp-spacetime analysis, with hypercomplex automatic differentiation (HyPAD)
designed into the scalar type from the beginning so parameter sensitivities come
from the same code that computes the objective.

    p → g_μν(x; p) → ∂_α g_μν → Γ^ρ_μν → R_μν → R → G_μν → T_μν → J(p)

> **Status: Milestones 1, A2, A2T, A3 and A4 complete, plus a verified AD backend
> and Program B1.** **258 / 258 verification cases pass**, 43/43 CTest.
> **Curvature exists as of A4** — Christoffel, Riemann, Ricci, and the
> R<sub>μν</sub>R<sup>μν</sup> and Kretschmann invariants, checked against
> Schwarzschild, static de Sitter and flat space in spherical coordinates.
> **Still no stress-energy:** no Einstein tensor, no T<sub>μν</sub>, no observer
> frames, no energy conditions, no optimizer. Nothing here supports any claim
> about warp-drive energy requirements, and there is as yet no quantity in this
> repository that may be called an energy.
> See [docs/limitations.md](docs/limitations.md) and
> [Program B's claims register](../program_b_em_inverse_design/docs/claims_and_limitations.md).

## Build and verify

Requires a C++17 compiler and CMake ≥ 3.20. The core has no third-party
dependencies; the AD backend needs the pinned OTIlib checkout.

```sh
# core only (no GPL dependency)
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo

# with the HyPAD backend (GPL-3.0, see below)
git clone https://github.com/mauriaristi/otilib.git external/otilib
git -C external/otilib checkout de3f7fa2b1ba667dff054ad5ba1f6f41f4ba8b54
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DWARP_HYPAD_WITH_OTILIB=ON

cmake --build build
ctest --test-dir build --output-on-failure
```

Each milestone gate is a process that returns a status:

```sh
./build/milestone1_verify --out results/verification/milestone_01
./build/milestone2_verify --out results/verification/milestone_02
./build/milestone3_verify --out results/verification/milestone_03
./build/milestone4_verify --out results/verification/milestone_04
./build/otilib_verify     --out results/verification/otilib
echo $?          # 0 = the milestone passes; anything else blocks the next one
```

Figures regenerate from saved output alone, without rerunning the solver:

```sh
python python/plot_milestone1.py results/verification/milestone_01
python python/plot_milestone2.py results/verification/milestone_02
python python/plot_milestone3.py results/verification/milestone_03
```

## Current result

**258 / 258 verification cases pass**, 0 failed, 0 skipped
(MSVC 19.44, RelWithDebInfo, Windows-AMD64, `-DWARP_HYPAD_WITH_OTILIB=ON`).
43 / 43 CTest tests pass (42 / 42 with the backend off), zero compiler warnings
at `/W4 /WX`.

| Gate | Cases | Highlights |
|---|---:|---|
| 1 — grid, storage, algebra | 77 | index map bijective; singular matrices rejected 200/200 |
| A2 — spatial derivatives | 41 | observed order **1.913–2.026**, interior *and* boundary |
| A2T — temporal derivatives | 23 | `∂ₜg = −v_s ∂ₓg` residual converges at order 1.82 |
| A3 — metric generation | 31 | `det g = −1` exactly; HyPAD vs analytic **1.1e-16** |
| A4 — covariant curvature | 43 | Schwarzschild `K = 48M²/r⁶`; de Sitter `R = +12/L²` fixes the sign |
| OTIlib AD backend | 23 | `d(A⁻¹)/dp = −A⁻¹(dA/dp)A⁻¹` at **1.9e-16** |
| [Program B1](../program_b_em_inverse_design/) | 20 | conformal invariance 4.4e-16; 3/3 non-metric media rejected |

By level: 152 × A (exact identity), 57 × B (analytic), 34 × E (grid
convergence), 15 × D (independent method).

### What A4 actually establishes

Four reference geometries, chosen so that no single bug passes all of them:

- **Minkowski, Cartesian** — Γ = 0 and Riemann = 0, exactly. Passes trivially for
  a pipeline that always returns zero, which is why the next case exists.
- **Flat space, spherical** — Γ ≠ 0 while Riemann = 0. Separates "the connection
  vanishes" from "the curvature vanishes"; catches a sign error in the ΓΓ terms.
- **Schwarzschild** — `R_μν = 0` while `K = 48M²/r⁶`. Separates Ricci-flat from
  flat, and the same geometry in isotropic coordinates gives `K` agreeing to
  3e-3 across two charts.
- **Static de Sitter** — the only reference with nonzero Ricci, so the only one
  that fixes the **sign**: `R = +12/L²` under the MTW (+,+,+) conventions. A
  pipeline with the opposite Riemann sign passes every index identity here and
  fails only this.

Plus the Alcubierre metric against an independently generated symbolic reference
(exact sympy derivatives, algebra written separately in numpy), which is the only
case exercising the genuinely time-dependent path.

Measured convergence: **1.98** (Schwarzschild K), **2.00** (flat-spherical
Riemann), **2.08** (Ricci symmetry), **1.63** (Alcubierre R, where the grid only
just resolves the wall).

Two refusals are verified rather than documented: a frozen-time treatment yields
no curvature at all, and a false `ExactlyStatic` claim is caught by probing the
metric before it is used. The size of the error the first one prevents is
measured, not asserted — discarding `∂₀g` changes the Ricci scalar by more than
the Ricci scalar itself.

### Milestone 1 detail

| Group | Cases | What it establishes |
|---|---:|---|
| `grid` | 24 | endpoints and midpoint exact; spacing uniform to 3.6e-15; index map bijective and round-tripping exhaustively; malformed specs rejected |
| `layout` | 7 | every rank-2/3/4 component stores and retrieves its own index tag exactly, and at the documented raw offset |
| `algebra` | 18 | multiply/determinant/inverse against hand calculations **and** an independent implementation; singular matrices rejected 200/200; ill-conditioning flagged on a derived threshold |
| `contraction` | 12 | each contraction against a closed form, including slot-selectivity; residual detectors verified against planted violations |
| `genericity` | 7 | kernels compile and run with `double`, `long double` and a custom scalar with no ordering, no equality and no implicit `double` conversion — bit-identical results |
| `config` | 9 | reference config parses to exactly the specified values; byte-exact round trip; `0.1+0.2`, `1/3` and `-DBL_MAX` survive bit-exactly; typos warn; invalid values rejected |

By verification level: 70 × A (exact identity), 3 × B (analytic), 4 × D
(independent method). No C, E or F cases yet — manufactured solutions and grid
convergence arrive with the derivative operators in Milestone 2, cross-code
validation last.

Full record, regenerated on every run:
`results/verification/milestone_01/verification_report.md`.

## How verification works here

A milestone is complete when a process returns zero — not when the code compiles
or the figures look plausible. Every case carries a mathematical reference, a
quantitative tolerance, a verdict, and the commit and compiler that produced it.

Cases are graded **A** exact identity · **B** analytic solution · **C**
manufactured solution · **D** independent method · **E** grid convergence ·
**F** cross-code validation.

No implementation verifies itself: `tests/reference/ref_linalg.hpp` uses
Gauss-Jordan with partial pivoting and recursive Laplace expansion precisely
because the code under test uses neither.

The gate has already earned its keep — it caught an ill-conditioning threshold
that was wrong by four orders of magnitude, letting matrices through whose
inverses had lost every useful digit. See [docs/design.md](docs/design.md).

## Documentation

| | |
|---|---|
| [docs/conventions.md](docs/conventions.md) | signature, index placement, curvature signs, units, frames |
| [docs/memory_layout.md](docs/memory_layout.md) | grid indexing, point-major storage, and why |
| [docs/design.md](docs/design.md) | architecture, the generic-scalar seam, verification design |
| [docs/limitations.md](docs/limitations.md) | what is missing and what is unresolved |

## Two things to know before continuing

1. **Time dependence is DECIDED: inertial frame with x_s(t).**
   ([docs/time_dependence_decision.md](docs/time_dependence_decision.md))
   The decision rests on an exact identity — for a rigidly translating bubble
   ∂ₜg = −v_s ∂ₓg — so the correct time-dependent formulation costs no more than
   a spatial derivative. `FrozenTimeDiagnostic` remains available for metric
   inspection and sensitivities, and **the curvature code must refuse it**
   (`valid_for_curvature()`).

   Still open before Milestone 5: whether to add the co-moving formulation as a
   curvature-invariant cross-check (recommended), and which observer family to
   use where v_s f > 1.

2. **The AD backend is real, verified, and GPL-3.0.** Upstream
   [OTIlib](https://github.com/mauriaristi/otilib) (pinned at `de3f7fa2`) is
   integrated as `OTIScalar` — 3 imaginary directions, truncation order 2, so
   gradient *and* Hessian come from one evaluation. It builds natively with MSVC
   despite upstream documenting "Windows only under WSL".

   GPL-3.0 is adopted for the research prototype. The backend stays **optional
   and isolated** behind `include/warp_hypad/sensitivity/`, and OTIlib is a
   pinned external checkout — never vendored into this tree.
   ([otilib_discovery_report.md](docs/otilib_discovery_report.md))
