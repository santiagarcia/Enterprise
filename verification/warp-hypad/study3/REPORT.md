# REPORT: warp_hypad study 3 (regularised, matched-progress A/B study)

Builder: Ada, 2026-10-27. **This is a builder report and has not been independently reviewed.** Saved by the main
session from Ada's hand-back, because the subagent's write was blocked.
- Protocol: `PROTOCOL.md`, written 15:59–16:00 CDT before any study-3 computation. Its SHA-256 (`05fd4082…a125bfe`)
  was recorded at 16:00:04 in `logs/protocol.sha256`, and the file was not changed afterwards.
- Repository: `Enterprise` at 477141d, unmodified, linked against the existing `build-oti/` libraries.
- Scope: every statement here is about the **Eulerian energy density ρ** (normal observer, TotalEffective source,
  Λ = 0, t = 0), family F2A, start S1, this L-BFGS port, λ1 = 3.75e-8 and λ2 = 3.75e-9, and working densities 4:1
  and 6:1. Nothing here is claimed for other observers, energy conditions, starts, or "the O2 residual" in general.
- Status labels: **reproduced / new / failed**.

## 1. Headline
1. **No run converged** (failed). All 8 runs stopped at the registered WALL_LIMIT.
   - λ1 runs: 175–245 accepted iterates. λ2 runs: 84–133.
   - C1' (‖∇F‖ ≤ 1e-5·F/‖y‖) and C2' (10-iterate relative decrease ≤ 1e-8) were never true at any iterate.
   - Final ‖∇F‖ was 45× the C1' tolerance (A 4:1 λ1) to 38,000× (B 4:1 λ2).
   - Therefore D2, D3, D4 and D6 are all **failed (not converged)**. No statement is made about minimisers.
2. **The regulariser bounded the design norm** (new).
   - ‖y‖ plateaued at 16.5–16.8 (λ1) and 37.5–39.3 (λ2) in every run. In study 2 A drifted to 95–139.
   - The non-convergence is therefore slow convergence inside a bounded region, not drift.
   - Conditioning is estimated, not measured: the logged L-BFGS γ is 3e4–4e5, so the curvature along the steps is
     about 3e-6 to 3e-5. Against λ = 3.75e-8 or 3.75e-9 this gives an apparent condition number of about 10^3–10^4.
3. **At λ1 (design norm ≈ 16.6), A and B are indistinguishable, and neither deepens the worst point after
   iterate 8** (new; not converged).
   - At the finals, q = 0.59–0.60 and r_E = 0.571–0.572 in all four runs.
   - Over iterates ≥ 40, in-loop q ≤ 0.88 for 4:1 runs and ≤ 1.04 for A 6:1; no iterate reaches 1.05.
   - Fixed-node quadrature is resolved at these designs: R on Dw at orders 4/6/8 agrees to 1e-3 relative, and E8 and
     E16 agree to 0.1 %.
4. **At λ2 (design norm ≈ 38), A deepens the worst point strongly; B does so far less, but not never** (new; not
   converged; trajectory statements only).
   - A, finals: q_A = 32.0 (4:1) and 4.73 (6:1), strongly resolution-dependent.
   - A, onset: the deepening starts after r_E has plateaued, from about iterate 40 at 4:1 and about iterate 110 at 6:1.
   - A 4:1, quadrature gap: R on the 8:1 Dw set exceeds the working-set R by 20 % (6:1: 0.08 %). This is consistent
     with between-node exploitation at 4:1 (interpretation).
   - B, finals: q_B = 1.05 (4:1; the final iterate is also a refinement iterate) and 0.87 (6:1).
   - B, trajectories: concentration at 13 (4:1) and 20 (6:1) iterates; maximum assessed q 1.67 (iterate 73) and 1.42.
5. **Matched progress (Gauss B1): by the pre-registered D5 rule, "B suppresses deepening" is not established at
   either λ.**
   - The rule is valid, but it tests the wrong phase. Every registered r_E level (0.40, 0.50, 0.60) is first crossed
     at iterates 10–15. At that point A and B are still nearly the same trajectory (assessed |q_A − q_B| ≤ 0.0011).
   - A's deepening at λ2 occurs later, at almost constant r_E.
   - Section 5 gives a descriptive late-phase comparison at overlapping r_E and ‖y‖. It is not a registered decision.
6. **Every run shows the pre-registered trajectory label "concentration along the trajectory"** (Vera B1).
   - All 8 runs share a transient at iterates 7–8: q = 1.38–1.42, r_E = 0.10–0.24, ‖y‖ ≤ 3.4. It occurs before B's
     first refinement.
   - At λ1 that transient is the only concentration in each run. Labels for each run are in section 4.

## 2. Protocol deviations and declarations
| # | Item | Effect |
|---|---|---|
| — | **Pilot passed** (Δ1 = 11.05 > Δ2 = 5.90; `logs/pilot_decision.txt`). λ1 and λ2 were unchanged. The 8:1 diagnostics cost about 7 % of an iterate, so they ran at every iterate (cost rule). | No deviation |
| — | **WALL_LIMIT reached by all 8 runs.** λ1 runs were stopped by their 2 h 15 min cap at 18:26; λ2 runs by the global deadline at 19:41:13 (pilot end + 3 h 30 min). | Registered state, not a deviation; it caused the non-convergence verdicts |
| D1 | `analyze.py report` was first run at 19:55, before the iterate assessments had finished (37 of 80 rows). It was rerun at 20:03 on the complete data, and every number here comes from the rerun. | None |
| D2 | `results/assess_list.txt` has 88 entries. 80 distinct (run, iterate) pairs were assessed; the other 8 duplicate an iterate already listed (for example, r_E 0.40 and ‖y‖ 10 are both crossed at iterate 10). All assessment priorities 1–5 were completed. | None |
| D3 | The outside-Ab floor is −1.752e-4 at (−2.10, 1.371), on the Ab boundary. It is a minimum over D \ Ab with that region's boundary included. It is far shallower than b, and design-independence outside Ab held at all 34,266 check points for every design, so the fallback whole-D search was never triggered. | None |
| D4 | PROTOCOL 3.4 asks for working-node ρ_min before and after refinement. Both were logged, and the "after" value is used for the fixed-node check. The in-loop ρ* is logged once per iterate: the design does not change at refinement, so ρ* does not change either. | None |
| D5 | Under the D5 rule, the ‖y‖ = 40 level was reached only by A 6:1 λ2 (iterate 71); B's maximum was 39.3. That level is reported as "not reached by both". | None |
| D6 | REPORT.md could not be saved by the subagent (harness block). The text was returned to the main session. | None |

No other departures. Runs used 2/2/4/4 threads (12 in total), as registered; assessments used 12 threads after the
runs had finished.

## 3. Checks
| Check | Result | Status |
|---|---|---|
| Driver base sets 4:1, 6:1, 8:1 vs `build_elementwise_set` | bit-identical | reproduced |
| E_neg and ρ_min vs `evaluate_feasibility`; R(4:1) vs `feasibility_residual` (2 designs) | identical | reproduced |
| Baseline ρ*_D (new assessment: 0.005 lattice, Halton 2^17, Nelder–Mead) | b = −9.60726526e-4 at (0, 2.00134) | reproduced (study 2: −9.60726526e-4) |
| Baseline E_neg (8:1) | E0 = 4.634569e-2 | reproduced |
| First 20 iterates of B 4:1 λ1 in the pilot and in the production run | identical F at every printed digit | new (determinism) |
| Axisymmetry, 2,000 random 3D points, baseline and 8 finals | max ratio 4.2e-13 | new (PASS) |
| Design-independence outside Ab (34,266 points per design, 88 designs) | 0 mismatches | new (PASS) |
| Metric validity and admissibility on a 0.025 lattice over Ab (21,125 points) | valid at every point in all finals | new |

## 4. Results at the final iterates (all non-converged)
b = −9.60727e-4 and E0 = 4.634569e-2. q = ρ*_D/b and r_E = 1 − E8/E0. ρ*_D is the minimum over two seeds and the
outside floor (the seeds agree in every case). "Check" is the worst of the fixed-node checks (working nodes, 8:1
nodes, `refine_around_minimum`). "Under" is (check − ρ*_D)/|ρ*_D|.

| Run | it | ‖y‖ | q | r_E | r_E on Dw | R4 / R6 / R8 on Dw (×1e-4) | E8/E16 − 1 | Under | ln N range | max\|β\| | Trajectory label (count, first–last, max q) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A 4:1 λ1 | 245 | 16.78 | 0.595 | 0.572 | 0.588 | 1.0991 / 1.0998 / 1.0999 | −0.08 % | 0.3 % | −0.11 to 4.01 | 0.92 | conc. along trajectory, not at final (2, it 7–8, 1.42) |
| B 4:1 λ1 | 175 | 16.59 | 0.601 | 0.571 | 0.586 | 1.1004 / 1.1006 / 1.1006 | −0.01 % | 0.1 % | −0.11 to 3.98 | 0.89 | same (2, 7–8, 1.42) |
| A 6:1 λ1 | 229 | 16.54 | 0.598 | 0.571 | 0.586 | 1.1006 / 1.1007 / 1.1008 | −0.05 % | 1.3 % | −0.12 to 3.97 | 0.89 | same (2, 7–8, 1.42) |
| B 6:1 λ1 | 194 | 16.50 | 0.592 | 0.571 | 0.586 | 1.1011 / 1.1010 / 1.1010 | +0.02 % | 1.0 % | −0.14 to 3.96 | 0.89 | same (2, 7–8, 1.42) |
| A 4:1 λ2 | 133 | 37.60 | **31.96** | 0.589 | 0.604 | 1.0436 / 1.1557 / **1.2545** | −0.17 % | 1.5 % | −0.02 to 6.80 | 4.15 | concentration along trajectory (105, 7–133, 31.96) |
| B 4:1 λ2 | 84 | 39.27 | 1.052 | 0.632 | 0.649 | 1.0534 / 1.0549 / 1.0551 | −0.05 % | 1.7 % | 0.00 to 6.70 | 4.33 | concentration along trajectory (13, 7–84, 1.67 at it 73) |
| A 6:1 λ2 | 126 | 37.54 | **4.73** | 0.630 | 0.647 | 1.0450 / 1.0469 / 1.0477 | −0.60 % | 2.0 % | 0.00 to 6.86 | 4.00 | concentration along trajectory (47, 7–126, 5.51 at it 117) |
| B 6:1 λ2 | 90 | 38.90 | 0.874 | 0.628 | 0.645 | 1.0515 / 1.0523 / 1.0524 | −0.01 % | 0.3 % | 0.00 to 6.67 | 4.29 | conc. along trajectory, not at final (20, 7–82, 1.42) |

Notes:
- **Objective balance at the finals.** The penalty (λ/2)‖y‖² is 5.1e-6 to 5.3e-6 at λ1 (4.7 % of F) and 2.6e-6 to
  2.9e-6 at λ2 (2.5 %).
- **Worst-point locations.** Deepest points lie at r ≈ 1.94–2.44. For A 4:1 λ2 it is a narrow pocket at
  (0.0008, 2.440), on the x ≈ 0 seam between the two positive lobes (`figures/meridian_4_L2.png`).
- **Fixed-node under-reporting** never exceeded 10 %: at most 2.0 %.
  - Unlike study 2, the 8:1 node set and `refine_around_minimum` caught A's λ2 pockets.
  - The working nodes alone did not. At A 4:1 λ2 the working-node ρ_min is about −1e-3, against ρ*_D = −3.07e-2
    (from `trace.csv`, `rhomin_w_post`).
- **Tail measures.**
  - V(ρ < 1.05·b) on the 8:1 set is nonzero only for A 4:1 λ2 (0.44).
  - L4 of ρ_− relative to the baseline: 0.48 at λ1 (all runs), 0.44–0.45 for B λ2 and A 6:1 λ2, and 4.82 for A 4:1 λ2.
- **Admissibility (Gauss I4).** All metrics are valid. At λ2 the lapse reaches N ≈ e^6.9 ≈ 950 and the shift
  |β| ≈ 4.3. These are large but finite designs. The present claims concern only their Eulerian ρ, not physical
  reasonableness.
- **Iteration stability (last 20 iterates, in-loop q).**
  - λ1: ranges are 0.002–0.036 wide, so the λ1 runs are iteration-stable by the 0.05 flag.
  - λ2: all are flagged "not iteration-stable": A 4:1 4.3–32.0, A 6:1 2.2–5.5, B 4:1 0.72–1.67, B 6:1 0.87–1.40.

## 5. Matched progress (D5) and late-phase comparison
**Pre-registered D5** (`results/d5_matched.csv`; assessed q at the first crossing):

| λ, density | r_E 0.40 (it A/B): q_A, q_B | r_E 0.50 | r_E 0.60 | ‖y‖ 10 | ‖y‖ 20 |
|---|---|---|---|---|---|
| λ1 4:1 | it 10/10: 0.9441, 0.9442 | 11/11: 0.7794, 0.7787 | not reached | 10/10: as for r_E 0.40 | not reached |
| λ1 6:1 | 10/10: 0.9438, 0.9438 | 11/11: 0.7773, 0.7773 | not reached | 10/10 | not reached |
| λ2 4:1 | 10/10: 1.0992, 1.0992 | 11/11: 0.9362, 0.9350 | 15/15: 0.5655, 0.5655 | 10/10 | 14/14: 0.7596, 0.7599 |
| λ2 6:1 | 10/10: 1.0988, 1.0988 | 11/11: 0.9332, 0.9332 | 15/15: 0.5654, 0.5654 | 10/10 | 14/14: 0.7605, 0.7605 |

**Verdict (D5):** "not established at matched E_neg progress" at both λ.
- The rule needs q_A − q_B ≥ 0.05 at ≥ 2 r_E levels at both densities. The largest difference observed was 0.0011.
- This is **failed as a discriminating test**, not evidence that A and B behave alike later on.
- The registered levels sit in the early phase, before B's refinements have changed its trajectory appreciably.

**Late-phase comparison** (descriptive, NOT pre-registered; `results/late_phase.txt`; accepted iterates ≥ 40):

| Run | r_E range | ‖y‖ range | median q | max q | fraction q ≥ 1.05 |
|---|---|---|---|---|---|
| A 4:1 λ2 | 0.588–0.631 | 35.6–40.0 | 3.13 | 31.96 | 0.99 |
| B 4:1 λ2 | 0.625–0.634 | 35.1–39.3 | 0.76 | 1.67 | 0.13 |
| A 6:1 λ2 | 0.624–0.636 | 34.2–40.9 | 1.03 | 5.51 | 0.46 |
| B 6:1 λ2 | 0.626–0.636 | 36.6–39.0 | 0.87 | 1.40 | 0.31 |
| λ1, all four | 0.569–0.581 | 16.5–18.8 | 0.60 | 0.85–1.04 | 0.00 |

Interpretation, labelled as such:
- At overlapping r_E (about 0.63) and design norm (about 35–40), A's worst point is far deeper than B's. B reaches
  equal or slightly better r_E.
- At 6:1 the difference is moderate: median q 1.03 vs 0.87, and B still exceeds 1.05 at 31 % of its late iterates.
- Two qualifications apply:
  1. A's q is iteration-unstable (it jumps between basins).
  2. The comparison pools iterates rather than matching individual pairs.
- A future pre-registration should match in the plateau phase: for example, compare at equal accepted-iterate windows
  after r_E has plateaued, or at equal R.

`figures/matched_L1.png` and `figures/matched_L2.png` plot q against r_E and against ‖y‖ for the full trajectories,
with assessed points marked.

## 6. Decisions under PROTOCOL §8
| Rule | Outcome | Status |
|---|---|---|
| D1 Converged | 0 of 8 runs (all WALL_LIMIT; C1' and C2' never true) | **failed** |
| D2 Regularisation dominance | Not evaluable (no converged B). Descriptive: B λ1 r_E = 0.571 ≥ 0.40 | failed (not converged) |
| D3 Concentration at convergence (A/B × λ1/λ2) | Not evaluable | **failed (not converged)** |
| D4 Fixed-node exploitation at convergence | Not evaluable. Descriptive, A λ2 finals: q_A(4:1) − q_A(6:1) = 27.2, but under-reporting is 1.5 % (< 10 %); E8/E16 differ by 0.17 % (< 2 %) while R8/R4 = 1.20 (> 1.02). Neither registered branch would hold as written. | failed (not converged) |
| D5 Matched E_neg progress | λ1 and λ2: "not established at matched E_neg progress" (max \|Δq\| 0.0011 at the registered levels) | new (verdict); failed as a discriminating test (section 5) |
| D6 λ dependence | Not evaluable for A or B | failed (not converged) |
| D7 Search sanity | See section 7: no deeper point than the polished value was ever found | new |
| T Trajectory labels | All 8: concentration along the trajectory. λ1 runs and B 6:1 λ2: "not at the final iterate". Details in section 4. | new |
| Iteration stability | λ1: stable (range ≤ 0.036); λ2: not stable (all four) | new / failed |

## 7. Search adequacy and independence (Gauss I1, I5; Vera I2)
**In-loop search** (I1 fix):
- On average 10.7 polishes per iterate, of which 1.3 duplicates were merged.
- 8 distinct polished minima were found at most iterates.
- B refinement events: 123 splits in total (35/43/16/29), at depths 0 (76), 1 (42) and 2 (5).
- No node cap was reached; nodes grew by at most 1.31× the base set.
- Refinement no longer resets L-BFGS memory: 0 memory resets in any run (B2 fix).

**Assessment against in-loop:**
- At 79 of the 80 assessed iterates, the assessed ρ*_D equals the in-loop ρ* within 1e-6 relative.
- At A 6:1 λ2 iterate 100, the assessment found a point 1.75 % deeper than the in-loop search.
- The assessment uses a structurally different method: a 0.005 lattice, Halton 2^17, Lipschitz cell screening, and
  Nelder–Mead.
- The agreement therefore says more than study 2's did. It is still not a certificate.

**Cell screening and D7:**
- Flagged cells appeared only in λ2 designs: A 17–29 and B 394–429 at the finals. A sub-lattice value fell below the
  pre-polish sample minimum at 11 assessments (B finals at λ2 and 9 iterates).
- In every case the polished ρ*_D was deeper than every sub-lattice value, so no D7 event changed ρ*_D.
- The heuristic cell estimate lies below ρ*_D by up to 76 % (B 4:1 λ2: −1.78e-3 vs −1.01e-3). This bounds nothing;
  it only marks where the 0.00125 resampling was done.

**Limits:**
- Pockets narrower than about 0.00125 inside unflagged cells, or narrower than 0.005 elsewhere, could be missed by
  both searches.
- ρ*_D is an attained value, and hence an upper bound on the infimum.
- The meridian reduction rests on the axisymmetry argument and check.

## 8. Limitations
- **Not converged.** All 8 runs stopped at WALL_LIMIT. Every result is a trajectory or final-iterate statement, not
  an optimum statement.
  - λ2 runs had only 84–133 iterates. Whether A's λ2 deepening continues, saturates or reverses, and whether B's
    excursions die out, is unknown.
- **The budget forced a narrow design:**
  - one start (S1);
  - no S3 runs;
  - no 1 % threshold sensitivity run for B (Gauss I2).
- **The regularisation choice frames the comparison.** The penalty is the coefficient-space Euclidean norm, which is
  basis-dependent. The results are budget-specific: no deepening at budget ≈ 17, strong deepening by A at ≈ 38.
  Nothing here extrapolates to the unregularised problem.
- **D5 was mis-targeted** (section 5). The late-phase table is descriptive and was chosen after seeing the data.
- **B is one adaptive design** (5 % rule, polished minima only). It reduces deepening but does not prevent
  intermittent concentration (λ2: q up to 1.67).
- **The driver is the builder's own code** (`driver/s3.cpp`). The shared paths are verbatim copies of study 2's
  s2.cpp and were checked bit-for-bit against the repository. The new code (search, refinement edge rule,
  regularisation, assessment) has not been independently reviewed.
- **Wall time:** 15:52 → about 20:10 CDT, about 4.3 h. The run phase took 16:11–19:41 at 12 threads.

## 9. Reproduction
```
bash study3/driver/build.sh                                  # g++ 16.1, study-2 flags
study3/driver/s3.exe check                                   # logs/check.log
OMP_NUM_THREADS=12 study3/driver/s3.exe timing               # pilot/timing_12thr.log
s3.exe run --form B --ws 4:1 --start baseline_zero --lambda 3.75e-8 --kmax 20 --out pilot/B_4_L1
bash study3/driver/lane.sh F D T DEADLINE_EPOCH 3.75e-8 3.75e-9   # lanes A4(2) B4(2) A6(4) B6(4); deadline 1793148073
bash study3/driver/assess_finals.sh 12
C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I study3/driver/analyze.py select
bash study3/driver/assess_iters.sh 12
C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I study3/driver/analyze.py report
```
- Runs are deterministic; node sums are serial. The pilot reproduced the production trajectory.
- WALL_LIMIT runs stop at a clock-dependent iterate, so a rerun reproduces their trajectories up to that iterate.
  Rerunning with `--kmax <it>` from `termination.txt` reproduces the final design.

## 10. Artifacts (under this folder)
- **Pre-registration and logs:**
  - `PROTOCOL.md`, `logs/protocol.sha256`, `logs/pilot_decision.txt`;
  - `logs/check.log`, `logs/assess_finals.log`, `logs/assess_iters.log`, `logs/analyze.log`.
- **Driver:** `driver/s3.cpp`, `build.sh`, `lane.sh`, `assess_finals.sh`, `assess_iters.sh`, `analyze.py`.
- **Pilot:** `pilot/` (timing logs, `B_4_L1` 20-iterate λ pilot, `A_6_L1` 3-iterate timing pilot).
- **Runs:** `runs/<F>_<d>_L<k>/`, one folder per run, containing:
  - `config.txt`, `trace.csv`, `ys.txt` (every accepted y), `linesearch.csv`, `inloop_minima.csv`;
  - `refine_events.csv`, plus `set_changes.csv` for B;
  - `final_diag.csv`, `termination.txt`, `y_final.txt`, `log.txt`;
  - plus `runs/runs.txt`.
- **Assessment:** `assess/assess.csv` (finals and baseline, 2 seeds, full); `assess/iters.csv` (80 iterates);
  `*.minima.csv`; `assess/snap/`.
- **Fields:** `fields/` (ρ on Ab at 0.01).
- **Results:**
  - `results/summary.csv`, `decisions.json`, `d5_matched.csv`;
  - `results/inloop_vs_assessment.csv`, `late_phase.txt`, `assess_list.txt`.
- **Figures:** `figures/trajectory_L1.png`, `trajectory_L2.png`, `matched_L1.png`, `matched_L2.png`, `resolution.png`,
  and `meridian_{4,6}_L{1,2}.png`.

## 11. Candidate manuscript wording (not reviewed)
"In family F2A (Eulerian ρ, t = 0, one start), we minimised a Tikhonov-regularised O2 residual for up to 245 L-BFGS
iterations; no run converged.
- At a design-norm budget of about 17, neither fixed-node nor adaptive-argmin quadrature deepened the continuum worst
  point beyond an early transient: q ≈ 0.6 at a 57 % E_neg reduction.
- At a budget of about 38, fixed-node minimisation deepened it 4.7–32× (strongly resolution-dependent). Adaptive
  re-insertion reached a similar E_neg reduction (63 %) with final q of 0.87–1.05, though intermittent excursions
  reached 1.67.
- A pre-registered matched-progress test was not discriminating, because both formulations crossed the registered
  E_neg levels before they diverged."

## 12. Next action
- **Independent review:**
  - Vera: rerun one λ2 pair to a fixed `--kmax` and audit the s3.cpp changes (search, edge rule, no-reset L-BFGS,
    assessment).
  - Gauss: assess the D5 mis-targeting and the late-phase comparison, and judge whether conditioning (about
    10^3–10^4) makes convergence reachable at all with L-BFGS m = 8.
- **Pre-registration for a study 4** (proposal only):
  - continue the λ2 runs from their archived `y_final` (warm start, same driver) to convergence or to a fixed
    iteration budget;
  - match in the plateau phase;
  - consider a preconditioned or Newton-type solver (the Hessian is only 108 × 108).
