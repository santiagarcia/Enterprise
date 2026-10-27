# REPORT: warp_hypad study 2 (matched A/B optimization study)

Builder: Ada, 2026-10-27. **Builder report, not independently reviewed.** Saved by the main session from Ada's
hand-back, because the subagent's write was blocked.
- Protocol: `PROTOCOL.md`, written at 11:45 CDT before any study-2 computation.
- Repository: `Enterprise` at 477141d, unmodified, linked against the existing `build-oti/` libraries.
- Every statement is about the **Eulerian energy density ρ** (normal observer, TotalEffective, Λ = 0, t = 0). No other
  observer or energy condition was evaluated.
- Status labels: **reproduced / new / failed**.

## 1. Headline
1. **No run converged.** All 10 runs (2 formulations × 2 densities × 2–3 starts) ended at ITERATION_LIMIT (150 accepted
   iterations). The pre-registered question "does concentration persist at convergence?" is therefore
   **failed (not converged)**. Final ‖g‖/‖g0‖ ranged from 1.5e-3 to 2.0e-1; the threshold was 1e-3.
2. **Formulation A (fixed nodes, current O2) is exploitable**, by the pre-registered rule (new).
   - At every matched start and density, A's final continuum worst point is far deeper than B's: q_A = 2.7–23.5
     vs q_B = 0.69–0.98.
   - The fixed-node check under-reports the worst point by more than 10 % at 4 of 5 A final designs.
   - The deepening shrinks as the working set is refined (q_A ≈ 14–23 at 4:1, ≈ 2.7–3.7 at 6:1). This is the
     signature of pockets dug between nodes.
3. **B (adaptive argmin re-insertion, same integrand and τ) shows no concentration at its final iterates**
   (q < 1.01 in 5/5), yet it lowers E_neg by **58–64 %**, comparable to A's 56–66 % (new; not converged).
   - In this family, at the iterates reached, reducing E_neg did **not** require deepening the worst point.
   - Study 1's "E_neg down, worst point deeper" is attributable to the fixed-node formulation, not to the integral
     objective as such. This holds only for the non-converged iterates examined.

## 2. Protocol deviations
| # | Deviation | Effect |
|---|---|---|
| D1 | **Scheduling.** Two priority pairs ran at once (four processes, 5–6 OpenMP threads each) instead of one pair at a time; P5 ran next to the S3 assessments. | Wall time only. Results are deterministic per run; node sums are serial. |
| D2 | **Trace-row timing for B.** `f` and ‖g‖ in a trace row are those of the accepted iterate on the set it was accepted on. The node columns (`rhomin_w`, `nodes`) are on the set after that iterate's refinement. | None on decisions. |
| D3 | **C2/STALLED window.** Implemented over iterates since the last working-set change, for both criteria. PROTOCOL states this for C2 only. | Never binding: no run stalled. |
| D4 | **Assessment of the active region and the design-independent floor.** Polish used 16 starts each, not 64. The global D and Dw polishes used 64 as registered. | Applies only to the ρ_active and floor columns, which equal ρ*_D at every final. |
| D5 | **Concentration criterion.** The ≥ 10 % E_neg condition was moot: every run has r_E ≥ 0.56. | None. |
| D6 | **REPORT.md** could not be saved by the subagent (harness block); the text was returned to the main session. | None. |

**Budget:** total wall time 11:45–15:23 = **3.6 h**, under the 6 h cap. The pilot projection was 3.3 h, so the
budget rule was not triggered and every priority (P1–P5) ran with K_max = 150.

## 3. Checks before running
| Check | Result | Status |
|---|---|---|
| Driver base sets (4:1, 6:1, 8:1) vs `build_elementwise_set` | bit-identical | new |
| Parallel E_neg / ρ_min vs repo `evaluate_feasibility` (2 designs) | identical | new |
| R(4:1) vs repo `feasibility_residual` (2 designs) | identical to the last bit | new |
| Baseline E_neg (8:1) | 4.634569e-2 | reproduced |
| Baseline continuum ρ* | −9.60727e-4 at (0, 2.001) | reproduced (study 1: −9.60727e-4) |
| O2 zero-start iterate 3 residual | 1.48943928e-4 | reproduced |
| Axisymmetry, 2,000 random 3D points | max relative deviation ≤ 2.8e-13 on the baseline, a study-1 design and all 10 finals | new (PASS) |

- The analytic argument for axisymmetry is in PROTOCOL §5. With the numerical check passing, the meridian search
  represents the full 3D domain D up to rounding.
- The design-independent floor (baseline ρ where no active mode reaches) is **−4.9985e-4** at (−1.05, 1.34), in the
  protected interior. No design in this family can make ρ*_D shallower than this, so q ≥ 0.52.

## 4. Results (final iterates; assessment seeds 1 and 2 agree to all printed digits in every case)
b = −9.60727e-4 and E0 = 4.634569e-2. q = ρ*_D / b, and r_E = 1 − E_neg(8:1)/E0.

"Check worst" is the minimum of three fixed-node checks:
- working-node ρ_min,
- 8:1 validation-node ρ_min,
- the repository's `refine_around_minimum`.

| Run | Term. | ‖g‖/‖g0‖ | q | r_E | ρ*_D (location) | E_neg 8:1 / 16:1 | Check worst | Under-report |
|---|---|---|---|---|---|---|---|---|
| A 4:1 S1 | ITER_LIMIT | 1.3e-2 | **23.5** | 0.559 | −2.256e-2 (1.22, 1.97) | 2.046e-2 / 1.999e-2 | −1.722e-2 | yes |
| B 4:1 S1 | ITER_LIMIT | 1.7e-2 | 0.687 | 0.611 | −6.597e-4 (−1.14, 1.98) | 1.801e-2 / 1.806e-2 | −6.425e-4 | no |
| A 6:1 S1 | ITER_LIMIT | 6.0e-3 | **3.75** | 0.653 | −3.600e-3 (0.76, 1.65) | 1.610e-2 / 1.623e-2 | −1.664e-3 | yes |
| B 6:1 S1 | ITER_LIMIT | 8.5e-3 | 0.892 | 0.643 | −8.573e-4 (0.79, 2.44) | 1.654e-2 / 1.652e-2 | −8.571e-4 | no |
| A 4:1 S3 | ITER_LIMIT | 2.5e-3 | **18.6** | 0.599 | −1.788e-2 | 1.857e-2 / 1.831e-2 | −1.616e-2 | no (9.6 %) |
| B 4:1 S3 | ITER_LIMIT | 9.4e-3 | 0.827 | 0.582 | −7.947e-4 | 1.937e-2 / 1.937e-2 | −7.923e-4 | no |
| A 6:1 S3 | ITER_LIMIT | 1.5e-3 | **2.70** | 0.659 | −2.592e-3 | 1.583e-2 / 1.583e-2 | −2.184e-3 | yes |
| B 6:1 S3 | ITER_LIMIT | 4.2e-3 | 0.754 | 0.628 | −7.239e-4 | 1.723e-2 / 1.723e-2 | −7.207e-4 | no |
| A 4:1 S2 | ITER_LIMIT | 1.9e-2 | **14.0** | 0.604 | −1.342e-2 | 1.836e-2 / 1.815e-2 | −1.009e-2 | yes |
| B 4:1 S2 | ITER_LIMIT | 2.0e-1* | 0.979 | 0.581 | −9.410e-4 | 1.944e-2 / 1.945e-2 | −9.303e-4 | no |

\* B 4:1 S2 refined at iterate 150. The ‖g‖ shown is after the set change.

- For every final design, ρ*_D = ρ*_Dw = ρ*(active region): the worst point lies in the design-active part of the
  working domain.
- **Under-reporting by the working nodes alone is far larger** than the column shows. For A, the final working-node
  ρ_min is −7.3e-4 to −1.4e-3, which looks no worse, or only slightly worse, than the baseline. The 8:1 validation
  set and the local refinement catch most of the pockets.
- **Refinement in B:**
  - 9–27 element splits per run, maximum depth 3;
  - final node count 1.06–1.20× the base;
  - the node cap was never reached.

## 5. Iteration and resolution dependence (figures/)
**Iteration** (`iteration_S1/S2/S3.png`):
- A's continuum ρ* deepens steadily after about 40 iterations, most sharply at 4:1 (S1: q ≈ 1 → 5 at it 80 →
  23.5 at it 150).
- Over the same iterations, A's 4:1 validation E_neg *rises* again (r_E 0.65 → 0.56), while the working objective R
  falls monotonically. The validation integral sees the pockets that the working rule does not.

**Range over the last 20 iterates** (in-loop q; all exceed the 0.05 iteration-stability threshold):

| Formulation | Range of q | Flag |
|---|---|---|
| A | 1.9 to 24.4 | not iteration-stable |
| B | 0.60 to 1.14 | not iteration-stable |

- B shows transient excursions above 1.05 (B 4:1 S1 up to 1.055; B 4:1 S2 up to 1.138). These are pockets that
  appear at an accepted iterate and are refined away at the next check.
- B's r_E is stable to ±0.003 over the last 20 iterates; A's 4:1 r_E varies by up to 0.04.

**Resolution** (`resolution.png`):
- From 4:1 to 6:1, q_A falls from 23.5 to 3.75 (S1) and from 18.6 to 2.70 (S3).
- q_B goes from 0.69 to 0.89 (S1) and from 0.83 to 0.75 (S3); |Δq_B| = 0.21 and 0.07, both above the 0.05
  resolution criterion.
- r_E changes by less than 0.05 between densities for both formulations.

**Maps** (`meridian_<density>_<start>.png`):
- Each map shows baseline, A and B side by side, with ρ = 0 contours, the 8 deepest polished local minima, ρ*_D,
  the B split triggers, and the working domain.
- A's deepest points sit in narrow pockets at the inner edge of the positive-ρ region, near (±1.2, 1.97),
  (0, 1.62) and (0.76, 1.65).
- B's minima stay broad and baseline-like.

## 6. Decisions under PROTOCOL §7
| Rule | Outcome | Status |
|---|---|---|
| Concentration persists at convergence (objective-design result) | Not evaluable: no B run CONVERGED | **failed (not converged)** |
| Formulation A exploitable | (i) under-reporting at 4/5 A finals; (ii) q_A − q_B ≥ 1.7 at all 5 matched pairs, and B shows no concentration where A does | **new: A exploitable** |
| B concentration on non-converged trajectories | Not at final iterates (q_B < 1.01 in 5/5); transient excursions to q ≤ 1.14 inside the last 20 iterates | new (trajectory statement only) |
| Iteration stability | A and B both flagged not iteration-stable | failed |
| B resolution stability (\|Δq\| < 0.05) | 0.21 (S1) and 0.07 (S3) | failed |
| Study-1 iteration-dependent deepening under A | Recurs at 150 iterations and 3 starts; larger at the coarser working set | reproduced |

## 7. Failed checks and limitations
- **Convergence:**
  - 10/10 runs hit the iteration limit. A 4:1/6:1 S3 came closest (‖g‖/‖g0‖ = 2.5e-3 and 1.5e-3).
  - Whether the objectives have a finite minimiser in this unregularised 108-coefficient family is unknown; ‖y‖
    kept growing in study 1.
  - All final-state statements are about iterate 150, not optima.
- **Worst-point search** (sampling plus multi-start local polish):
  - It is not a certified bound. ρ*_D is an upper bound on the infimum over D.
  - Two Halton seeds and the adaptive subdivision agree, but a narrower pocket than the 0.01 lattice plus the
    3-level subdivision can resolve could still be missed.
- **Partial independence of B's check:** B's in-loop search and the final assessment use the same method family
  (lattice plus compass polish), though with different lattices, a Halton design and the adaptive stage. B was
  refined exactly where an in-loop search found pockets, so its good final assessment is partly expected. A pocket
  that both searches miss would be missed by both.
- **Search domain:** D is the O2 validation domain (|x| ≤ 7, r ≤ 7) at t = 0 only. There is no time dependence and no
  observer other than Eulerian.
- **B is one design choice:**
  - the 5 % resolution threshold, depth 6, 4× cap and in-loop lattice were fixed in advance, not tuned;
  - other adaptive or constraint formulations (for example a worst-point penalty) were not run.
- **Driver authorship:** the driver is Ada's own code (`driver/s2.cpp`), checked bit-for-bit against the
  repository's integration and residual routines but not independently reviewed.

## 8. Reproduction
```
bash study2/driver/build.sh                      # g++ 16.1, study-1 flags + -fopenmp
study2/driver/s2.exe check                       # identity checks (logs/check.log)
bash study2/driver/run_pair.sh 4 S1 baseline_zero 6     # likewise 6 S1, 4 S3/6 S3 combined_perturbation, 4 S2 lapse_perturbation
bash study2/driver/assess_run.sh A_4_S1 24       # per run: 2-seed final assessment, axi check, field, snapshots
C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I study2/driver/analyze.py
```

## 9. Artifacts (under C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad/study2/)
- `PROTOCOL.md`
- `driver/`: `s2.cpp`, `build.sh`, `run_pair.sh`, `assess_run.sh`, `analyze.py`
- `pilot/`: the timing pilot
- `runs/<F>_<d>_<S>/`, one folder per run, containing:
  - `config.txt`, `trace.csv`, `ys.txt` (every accepted y), `linesearch.csv`,
  - `refine_events.csv`, `set_changes.csv` (B only),
  - `termination.txt`, `y_final.txt`, `log.txt`;
  - plus `runs/runs.txt` with timings.
- `assess/`:
  - `assess.csv` (finals and baseline, both seeds), `snapshots.csv` (every 25 iterates);
  - `assess.csv.<run>.minima.csv`;
  - `snap/`
- `fields/`: meridian ρ grids
- `figures/`
- `results/summary.csv`
- `logs/`

## 10. Next action
- **Reviewers.**
  - Vera: rerun one matched pair, for example 4:1 S1, and the assessment, then audit `refine_B` and the stopping
    logic in `s2.cpp`.
  - Gauss: assess whether the B refinement rule and the search are adequate.
- **The convergence question needs one of:**
  - a much larger iteration budget, or
  - a bounded or regularised design space (a protocol change requiring a new pre-registration).
- **Candidate claim for the manuscript** (not yet reviewed): "with fixed-node quadrature the O2 residual is
  exploitable; with adaptive argmin re-insertion E_neg falls by a similar amount without deepening the Eulerian
  worst point at the iterates reached."
