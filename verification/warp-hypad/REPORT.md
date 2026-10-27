# REPORT — warp_hypad feasibility round

Builder: Ada, 2026-10-27. **Builder report, not independently reviewed.** Saved by the main session from
Ada's hand-back, because the harness blocked the subagent's report write; content as delivered.
Repo: local clone of `Enterprise` @ 477141dd4a322edff4b71b6e02416dacd923a676. Nothing pushed and no
source edited. Workspace: this directory. Labels: reproduced / failed / new / not attempted.

## 1. Environment and commands
- Windows 11 10.0.26200, Intel Core Ultra 9 285K (24 threads), 127 GB RAM.
- GCC 16.1.0 (WinLibs UCRT), CMake 4.4.1, Ninja 1.13.2. RelWithDebInfo, `-ffp-contract=off`.
- OTIlib vendored with PINNED_SHA de3f7fa (pin **not verified** against upstream history).
- No matplotlib in any local Python, and PyPI was blocked.

```
FL="-Wno-error=literal-suffix -Wno-error=unused-but-set-variable -Wno-error=unused-variable -Wno-error=unused-function -Wno-error=misleading-indentation"
cmake -S Enterprise/warp_hypad -B build-oti -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo -DWARP_HYPAD_WITH_OTILIB=ON "-DCMAKE_CXX_FLAGS=$FL"   (and build-nooti with OFF)
cmake --build build-oti -j 16; ctest --test-dir build-oti -j 8 --timeout 1800 > logs/ctest-oti.log
.venv/Scripts/python.exe driver/count_cases.py build-oti/results/verification
```

**Verification driver** `driver/wh_probe.cpp`:
- Links the repo's `libwarp_hypad_core.a`, `libwarp_hypad_otilib.a` and `tests/o0d_convergence_checks.cpp`, with the same flags plus `-fopenmp`.
- Modes: `o2` (a verbatim port of the o2_verify L-BFGS with configurable working/validation sets), `fd`, and `eval` (ladder plus a 1200–1600² meridian scan with polish).
- Node values are computed in parallel and summed serially, so integrals are bit-identical to the repo's.
- Scripts: `driver/run_refine.sh`, `driver/run_it40.sh`, `driver/analyze_fd.py`.

**Why the warning flags were needed:** GCC `-Wliteral-suffix` fires in OTIlib `comm.h`, and the tests have unused variables and misleading indentation (`o0c_gradient_checks.cpp:526` checked: benign). The MSVC `/W4` build does not flag these.

## 2. Test suite
| Item | Repo claim | Here | Status |
|---|---|---|---|
| 8 gates in README group | 285/285 (reported) | 296/296 (A5 now 38, was 27) | reproduced; README count stale |
| Program B, 11 gates | 417/417 | 417/417 | reproduced |
| All 25 gate summaries, OTI | n/a | 868/868 pass | new |
| All 25 gate summaries, no OTI | n/a | 780 pass, 50 declared skips, 0 fail | new |

**CTest, OTI build: 63/94.**
- 30 failures are figure, animation or visual tests (matplotlib missing); the v0a6 visual gate fails because its SVGs are absent for the same reason.
- `o2.gate` timed out at 1800 s. It runs serially and ran under CPU contention, so it counts as **not completed**.

**CTest, no-OTI build: 63/90.** All 27 failures are matplotlib-related.

**O1:** reproduced to 12 digits (E_neg 0.046345695595 → 0.046191311534, 14 iterations, converged).

**Side effect of running the tests:** they rewrite tracked Program C1 records. The diff is archived in `results/patches/` and the files were reverted. The committed C1 record shows `git_commit 893dcfd, git_dirty true`, a Program C provenance flag.

## 3. Terminology; coordinator's facts verified
**A-series grid route (A2–A6):**
- Spatial derivatives use 2nd-order FD stencils (`geometry/finite_difference.hpp`; Γ derivatives are differenced again).
- OTI is applied only to metric sensitivities with respect to v, R, σ (A3) and to temporal checks.
- No parameter derivative runs through curvature or energy on this route.

**O-series design route (O0B–O2):**
- `adm_metric_jet` / Jet3 propagates g, ∂g, ∂²g, ∂³g analytically; `geometry_bridge.hpp` states no FD is used anywhere.
- The chain Γ → Riemann → G → T → ρ = T(n,n) is evaluated at quadrature nodes.
- p is the 108 B-spline design coefficients (3 fields × 36 active). v, R and σ are plain doubles on this route.

**What the derivatives are:** ∂J/∂p and ∂²J/∂p² are derivatives of the **discretised** objective, not the continuum one. That discretisation is a fixed-weight Gauss sum over nodes on a truncated domain, with a softplus surrogate, a finite basis, and order-2 truncated OTI arithmetic. "Exact" holds only as exact derivatives of that implemented sum, up to rounding.

**Fact (1), verified.**
- O2 minimises `feasibility_residual` = Σ softplus_τ(−ρ)²·w with τ = 1e-3; the manifest says so.
- R uses the raw weights and E_neg uses w·√det γ; for F2A γ = I, so they coincide.
- The brief's "minimising integrated negative energy" is wrong for O2 and correct for O1.

**Fact (2), verified.** Refinement on this route means quadrature order, subdivision and extent, plus the sample density for ρ_min.

**Fact (3), verified.** onumm3n2 gives 3 directions at order 2; full Hessians are not assembled.

**Stale documentation:**
- `limitations.md` still says "no observers/energy/optimizer" and "no FD verification of full-field sensitivities".
- O0C does contain a Richardson-FD gradient check of E_neg_tau (3 coordinates, gradient only).

**ε_E record:** O2's `feasibility_tolerances.json` says ε_E is 10 × the O0D uncertainty. But `o0d_state()` is never populated in the O2 process, so ε_E falls to the 1e-9 floor and becomes 1e-8. The O0D value is 1.207e-9, which would give 1.207e-8. The verdict is unaffected, but the recorded derivation is inaccurate.

## 4. Derivative verification
**Setup** (data in `results/fd/fd_all.csv` and `fd_summary.md`):
- Directions: one lapse coefficient, one axial-shift coefficient and one dense unit vector per point.
- Steps h = 1e-1…1e-9. Central 2-point FD for the gradient, 3-point for the diagonal Hessian, 4-point for mixed entries.
- Points: y = 0, the combined start, and the O2 zero-start optimum.
- Sets: working 4:1 (6,624 nodes) and validation 8:1 (246,656 nodes).

**R (the actual O2 objective), reproduced/new.**
- Gradient best relative error 4e-11 to 2e-7; Hessian blocks 6e-9 to 4e-6.
- The error sweep has the expected V-shape: the gradient is best near h ≈ 1e-5 to 1e-4 and the Hessian near h ≈ 1e-3. Round-off dominates below about 1e-7 (gradient) and 1e-5 (Hessian).

**E_neg_tau:** gradient ≤ 2e-8, Hessian ≤ 3e-6.

**E_neg (max(−ρ,0), branch on the real part).**
- At y = 0 and at the combined start: gradient ≤ 1e-9, Hessian ≤ 8e-7. Every node is negative at y = 0, so E_neg is smooth there; it is exactly quadratic in the axial coefficient.
- At the O2 optimum (negative fraction 0.99) nodes cross zero. The gradient agrees to about 1e-7, but only for h ≤ 1e-5; at h = 1e-3 it is off by 4e-3. Some Hessian entries are 5e-2 or worse.
- **Conclusion:** E_neg is piecewise smooth. OTI returns a one-sided, branch-local derivative, and the Hessian is unreliable near sign changes.

**ρ_min (min over nodes; OTI derivative of the active node).**
- Combined start: the argmin is unique, its gap to the second node is 4.5e-5, and agreement is 1e-10.
- **Baseline: nondifferentiable.**
  - The nodes at x = ±0.028 (validation set) are exactly tied by x-symmetry, so central FD error is 0.06 to 0.7 and Hessian FD diverges like 1/h.
  - The one-sided FD equals the min or max of the tied OTI gradients: forward 4.1385e-5 vs min 4.1386e-5, backward 7.3997e-5 vs max 7.3994e-5.
- **O2 optimum: still tied.** The x = ±1.348 pair matches to 2e-15 because the zero-start design stays x-symmetric. Directional derivatives along the dense direction are +4.09e-4 and −5.62e-4.
- **Conclusion:** treat ρ_min as Clarke-subdifferentiable only.

**Not attempted:** the full 108×108 Hessian, and derivatives with respect to v, R or σ (that path is not OTI-enabled on the O-route).

## 5. O2 reproduction and refinement
**Original resolution** (working 4:1, 6,624 nodes; validation 8:1, 246,656 nodes; 7.6 min at 12 threads):
- Baseline: ρ_min −9.6024e-4, E_neg 4.634569e-2.
- Zero-start best: ρ_min −1.037949e-3, E_neg 2.073573e-2.
- That is **E_neg −55.26 %, ρ_min −8.09 %, reproduced.**
- All 5 starts end at `accepted_iteration_limit` (12 iterations), so **none converged**; the gradient norm was still falling.

**Fixed-design evaluation ladder** (orders 4/6/8/12/16; subdivisions 1/2/4; extent 5/25):
- Baseline E_neg is constant to 1e-10 relative.
- The optimum's E_neg spans 2.07352e-2 to 2.07400e-2 (about 5e-5 to 2e-4 relative, non-monotone because of the max kink). No clean convergence order.
- Node ρ_min wanders between −1.019e-3 and −1.039e-3 as sample points move.
- Continuum ρ_min (scan plus polish): baseline −9.60727e-4, optimum −1.04096e-3. That is **−8.35 %: the deepening survives.**

**Working-set refinement** (zero start, 12 iterations, validated at 8:1, plus continuum):
| Working set | Nodes | E_neg | Δ | Continuum ρ_min | Δ |
|---|---|---|---|---|---|
| 2:1 (coarser) | — | 2.08135e-2 | −55.09 % | −9.8939e-4 | −3.0 % |
| 3:1 | 3726 | 2.07389e-2 | −55.25 % | −1.06247e-3 | −10.6 % |
| 4:1 (original) | 6624 | 2.07357e-2 | −55.26 % | −1.04096e-3 | −8.35 % |
| 4:2 | 26496 | 2.07366e-2 | −55.26 % | −1.04212e-3 | −8.47 % |
| 8:1 | 26496 | 2.07365e-2 | −55.26 % | −1.04212e-3 | −8.47 % |
| 8:2 (26 min) | — | 2.07365e-2 | −55.26 % | −1.04211e-3 | −8.47 % |

At the 12-iteration snapshot both numbers are converged in working resolution from 4:1 upward (ρ_min to about 1e-3 relative, E_neg to about 4e-5). The coarsest set changes the size of the deepening, not its sign.

**Iteration dependence (new; failed robustness check).**
- ρ_min is non-monotone along the trajectory. Working-set values run −9.6e-4 → −1.33e-3 (iteration 7) → −8.65e-4 (iteration 11) → −1.035e-3 (iteration 12).

| Run | E_neg | Working ρ_min | Validation ρ_min | Continuum ρ_min |
|---|---|---|---|---|
| 4:1, 60 iterations | 1.6748e-2 (−63.9 %) | −5.69e-4 ("41 % better") | −2.058e-3 | −3.867e-3 (4.0× baseline) |
| 4:1, 40 iterations | 1.6664e-2 (−64.0 %) | — | — | −2.029e-3 (−111 %) |
| 8:1, 40 iterations | 1.6376e-2 (−64.7 %) | (= validation nodes) | −6.89e-4 ("+28 % better") | −1.188e-3 (−24 %) |

- In the 60-iteration run, the gate's 25×25 refinement found −3.861e-3; the deep pocket sits on the axis between nodes (x≈0, y≈1.63).
- In the 8:1, 40-iteration run the working nodes coincide with the validation nodes, yet the optimizer still digs a pocket between them; refinement found −1.176e-3.
- **Interpretation:** the optimizer exploits gaps between quadrature nodes at every density tested, so node-sampled ρ_min is not a valid certification quantity on its own. The repo's local refinement catches this. With more iterations, "E_neg down, worst point deeper" grows stronger, but its magnitude depends on the iteration count and the sampling.

## 6. Status
| Item | Status |
|---|---|
| Build (GCC; flag-only changes) | reproduced |
| 417 cases | reproduced |
| 285 cases | stale; 296/296 now |
| Figure and visual tests | failed (environment: matplotlib) |
| o2.gate via CTest | not completed (timeout) |
| O1 | reproduced |
| O2 numbers | reproduced |
| Smooth objectives, OTI vs FD (gradient + Hessian blocks) | reproduced/new |
| E_neg and ρ_min derivatives | nonsmooth, flagged |
| 12-iteration concentration under grid refinement | survives (new) |
| Concentration as a converged-optimum property | failed (not converged; iteration-dependent; node-gap exploitation) |
| Warp Factory cross-check, novelty, v/R/σ derivatives | not attempted |

## 7. Next decisive experiment
1. Rerun O2 to a real stopping criterion with an objective or constraint on a sampling-robust worst-point measure: a continuum ρ_min from scan plus polish, or an adaptive node set that re-inserts each iterate's argmin.
2. Report E_neg and continuum ρ_min along the whole trajectory at 2 or more working densities.
3. If concentration persists at convergence, it is an objective-design result; otherwise the finding is "integral objectives on fixed quadrature are exploitable".
4. Rough cost: about 1–2 hours of machine time.

## 8. Readiness (kept separate)
- **Scientific promise: moderate.** The methods point is clean: integral objectives concentrate or hide violations, and fixed-node certification is exploitable. Novelty is unknown (see the Scout literature file).
- **Evidence readiness: low to moderate.**
  - Strong: the code and test suite reproduce, and derivatives of the actual objective are verified.
  - Weak: the headline number is a non-converged snapshot; ρ_min is nonsmooth; there is no domain review and no Warp Factory comparison.
- **Provenance, unresolved:**
  - No top-level LICENSE; the OTIlib pin is unverified.
  - AI co-authorship appears in the commits ("Co-Authored-By: Claude Opus 5").
  - Stale README and limitations docs; an inaccurate ε_E derivation record; the Program C1 record comes from a dirty commit.
- **Effort to a defensible note:** the §7 experiment, a sampling-robust certification, the novelty search, Noether/Gauss review, and doc and licence cleanup. Moderate effort; not a schedule.

## Artifacts
- `driver/`
- `logs/`
- `results/o2_w4s1_v8s1/`
- `results/refine/*/` and `results/refine/runs.txt`
- `results/eval/ladder_orig.csv` and `results/eval/ladder_refine.csv`
- `results/fd/`
- `results/patches/`
- `build-oti/` and `build-nooti/`

Limitations:
- The continuum ρ_min comes from a dense scan plus coordinate polish of the z=0 meridian, relying on axisymmetry. It is not a certified global minimum.
- The driver is Ada's own code and needs independent review.
- About 3 hours of machine time in total.

Reviewer handoff: Vera (numerics and driver fidelity), Noether (GR conventions). Reviewers should rerun `results/refine/w4s1_it60` and inspect `results/eval/ladder_refine.csv`.
