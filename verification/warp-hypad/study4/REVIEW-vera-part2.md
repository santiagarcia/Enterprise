# Independent review of REPORT-v4.md Part 2 (§7–13) (Vera, 2026-10-28)

Saved by the main session from Vera's hand-back (read-only review). Text unchanged except formatting.

**Verdict: CLEARED WITH CORRECTIONS.** I found no error that changes a registered decision, a classification or a licensed claim. Every number I recomputed matches the raw outputs. There are seven corrections, all to labels, wording or declarations, plus one status update (P4).

**Role and independence**
- I am the reviewer, not a builder. I did not write Part 2, any driver or analysis script, or any run.
- Possible bias: I wrote the R3 ruling that allowed P1 to start. I therefore have a stake in this resumption looking procedurally clean, so I looked hard for procedural faults.
- Anchoring: I read Part 2 before recomputing. To limit it, every check below uses my own scripts on the raw `trace.csv`, `termination.txt`, `assess/*.csv`, `logs/assess.log` and `hessians/final_*` / `s2_*` binaries. I used the builder's code only where I was checking that code against the protocol text (`d5w.py`, `analyze_runs.py`, `v4_extra.py`).

## Findings, ranked by severity × likelihood × consequence

**F1 (medium; certain; wording only). The registered "no study-4 refinement check" wording is missing for B6L2.**
- PROTOCOL-v4 §4.3 (lines 379–380) says: "if B made no refinement check, B's D3/D1 statements read 'failed (not converged; no study-4 refinement check)'".
- B6L2 has n_check = 0 (`runs/B6L2/termination.txt`).
- REPORT-v4 §10 D1 (l. 214) and D3-B (l. 219–220) give only the terminal state and "(n_outer 0, n_check 0)". The cause is that `v4_extra.py` (l. 157–160) does not produce the phrase.
- Correction: add the registered phrase to both the D1 and D3-B statements for B6L2.

**F2 (medium; certain; no change to what may be claimed). The 4:1 D5w combined label follows `d5w.py`, not the protocol text, and this is not declared.**
- PROTOCOL-v4 §8 (l. 575): the verdict "must agree at W = 20 and W = 25; otherwise 'not established (window-sensitive)'".
- At 4:1, W = 20 gives "B lower" (4 windows) and W = 25 gives "not evaluable" (3 windows). Read literally, the combined verdict is "not established (window-sensitive)".
- `driver/d5w.py` `across_W` (l. 50–53) applies the gap-companion rule (§8 l. 628: "if either W is not evaluable, the result is not evaluable") to the verdict as well. That file was hashed before any data existed (21:35), so the choice was not tuned to the result.
- Under §8's licence table (l. 655), both labels give "No level statement". No licence changes.
- Correction: declare this as a reading of the pre-registered code (a new deviation 9), and give the literal label alongside. Part 2 l. 250–251 repeats "B lower" and "gap grows" at W = 20. The manuscript must not quote those, even descriptively.

**F3 (medium for the manuscript; certain). A required quantity, λ_min(H_R)/λ, is not reported, and it shows that minimiser status depends on the regulariser.**
- §2.3 (l. 241) requires λ_min(H_R)/λ at every Hessian design.
- From `hessians/final_*/summary.csv`, which I confirmed with my own eigen-analysis:

| Final | λ_min(H_R)/λ |
|---|---|
| A6L2m | −0.90 |
| A4L2f | −0.21 |
| A6L2f | −0.90 |
| X4L1fromL2 | −0.18 |

- So H_R is indefinite at every converged final, and H_F is PD only because of the λI term. At A 6:1, μ_min(H_F) = 3.77e-10 ≈ 0.10λ.
- Part 2 §9/§11 correctly scope the claim to "the regularised, discretised F", but without this column the dependence on λ cannot be seen.
- Correction: add the λ_min(H_R)/λ column to the §9 table.

**F4 (low-medium; certain). Hash and FAILED-line statements are stale, and one fix cannot be verified.**
- §12 deviation 7 (l. 349–350) and §13 (l. 356) say `sha256sum -c logs/driver.sha256` gives "one expected FAILED line".
- Since 23:40 there have been **two**: the run-1 hash of `fd_checks_v4.py` and the pre-fix hash of `v4_extra.py` (`be2281d9…`). I reproduced this, and `logs/p4_prerun_hashcheck.txt` (07:22) also records two.
- §7 (l. 111) is correct for the 23:38 checkpoint.
- Deviation 1 says the `v4_extra.py` fix "changed only the output encoding". This cannot be checked, because the pre-fix file was not archived (unlike `fd_checks_v4_run1.py.txt`).
- Plausibility: the failed run wrote nothing, and the run happened while P1 was running and before any run result existed.
- Correction: state 2 FAILED lines, and record that the pre-fix version was not archived.

**F5 (low-medium; certain; status). P4's status has changed since Part 2 was written.**
- `logs/execution.log` l. 42–43: P4 was started at 07:22:27 under unchanged v4 at the coordinator's request, 72 s after Part 2 was appended (07:21:15).
- The seven Part 2 result hashes are still byte-identical: `decisions.json` 1c0a047e, `v4_decisions.json` f8e4be07, `runs.csv` 11477a51, `v4_eps.csv` ef92d174, `assess.csv` 4283d5b7, `starts.csv` f6ccda6b, `iters.csv` f399a43c.
- Effect of P4 on Part 2's conclusions:
  - **D6x cannot change.** Neither arm keeps its start character, so labels 1 and 2 cannot apply whatever A4L1f does.
  - **D6 stays "not evaluable".** λ2 D3 is "failed (not converged)", and that reason alone is enough.
  - **Only λ1 D5w is pending.**
- Deviation 5 was legitimate under §5 (l. 435–436). Its session-time reasoning is now moot.
- Correction: add a dated note to l. 132, l. 266, l. 268 and §12.5 saying that P4 was started at 07:22 and will be reported in Part 3.
- Before P4's analysis, fix the CRLF writer in `analyze_runs.py select` (l. 82) and hash and log the change. Otherwise Part 3 will hit the same failure as deviation 3.

**F6 (low; certain; interpretation). The D4 results discriminate little.**
- (a) The only out-of-sample application, A6L2m, started with u_w = 72.9 %, already well above the 10 % threshold.
- (b) X4L1fromL2 (descriptive) is **shallow** (q 0.72) but has u_w = 21.7 %, so it is labelled "fixed-node exploitation". The study-3 separation (0.3–4.6 % against 87.9–98.2 %) therefore does not hold at a new design.
- A D4 label must not be read as "has a pocket", and D4 must not be described as validated out-of-sample.
- Part 2 is not wrong here, but §11 (l. 312) should state this limitation.

**F7 (low; certain or near-certain). Registered labels are not attached where the protocol says they should be.**
- §8 D1b (l. 495) says D4 statements quote the D1b label. The D4 table rows for the four converged runs (l. 234–236, 240) do not.
- §8 D7/T (l. 700–701) says "every final-state sentence carries the T label". The §11 sentences (l. 308–313) do not.
- The NV3 limitation is stated once (l. 181–182), not with each D1b use in §10 and §11.
- Lower confidence: the §8 general rule (l. 486) says *any* label derived from a non-converged run gets "(not converged)". By that rule, the T labels of A4L2m, B4L2, B6L2, X4L2fromL1 and XB4L2fromA (l. 282–290) arguably need the suffix.

**F8 (low; certain; reporting completeness).**
- The line-search-failure bullet (l. 137–138) gives the C2' value for 2 of the 3 runs. X4L2fromL1's is 1.08e-7.
- D7 events also occurred at 9 assessed trajectory iterates: B4L2 @100, 117, 125, 150 and 175; X4L1fromL2 @125; XB4L2fromA @100, 105 and 125 (`logs/assess.log`). §10 D7 lists only finals and starts.

## Verification results (all pass unless noted)

**1. Hash integrity.**
- PROTOCOL v1–v4 match their `.sha256` files (v4 is `0d77d7d5…c3ad`).
- `s4.cpp` is `ef247445…9463`. All current driver files check OK; the only FAILED lines are the two superseded ones (F4).
- The five FD `.out` hashes equal the v4 §9 table.
- Script timing, from `logs/execution.log`:
  - `v4_extra.py`: hash `be2281d9` logged 23:40:16, before its first run. The failed run and fix were at 23:40:24, and the new hash `35dc5a3e` was logged before the re-run. The file's mtime (23:40:24.76) is consistent with this.
  - `s2_newton.py`: written 06:40:50.870, `driver.sha256` updated 06:40:50.902, first `y_newton.txt` at 06:41:06. It was hashed before first use.
  - Both scripts were written before the outputs they process existed (P1 results 03:42; finals' Hessians 06:01 onward).

**2. Run classification** (my own C1'/C2' recomputation from `trace.csv`).
- My C1/C2 flags match the logged flags at every iterate of every run. The only exceptions are the two B check rows, where `nodes` is logged after the split while C1/C2 are evaluated on the set before it. That is bookkeeping, not an error.

| Run | Terminal state | Recomputed facts |
|---|---|---|
| A4L2m | LINESEARCH_FAIL (it 91) | C1' at 87–91 (5); C2' 3.93e-8 |
| B6L2 | LINESEARCH_FAIL (it 56) | C1' at 54–56 (3); C2' 2.31e-6 |
| X4L2fromL1 | LINESEARCH_FAIL (it 65) | C1' at 61–65 (5); C2' 1.08e-7 |
| B4L2 | ITERATION_LIMIT (it 200) | inner C1'∧C2' at it 70 → check → SPLIT 2 |
| XB4L2fromA | ITERATION_LIMIT (it 150) | check at it 72 → SPLIT 6 |
| A6L2m | CONVERGED | first C1'∧C2' at its final iterate (49) |
| A4L2f | CONVERGED | first C1'∧C2' at its final iterate (72) |
| A6L2f | CONVERGED | first C1'∧C2' at its final iterate (66) |
| X4L1fromL2 | CONVERGED | first C1'∧C2' at its final iterate (126) |

- None of the three LINESEARCH_FAIL runs meets C1'∧C2', so none is CONVERGED_LS_LIMITED. The classification is correct, and nothing was repaired.
- Precision-limited iterates: 0 (max ê_g/tol = 3.2e-8). STALLED: 0. Refresh triggers are all "grad". Hang guard: maximum wall/cap is 0.33. Block wall times reproduce from the logs.

**3. D1b / marginal / S2** (my own `eigh` on H_R_fwd + λI; τ_H = max(1e-9·max|μ|, 10·max(e_H, e_H^c))).
- `y.bin` equals `runs/<R>/y_final.txt` bit for bit at all four finals, and λ is correct (λ1 for X4L1fromL2).
- μ_min/τ_H = 1.399e4, 1.694e4, 1.399e4, 9.419e4. At every final there are 0 negative, 0 near-zero and 0 marginal-band eigenvalues.
- H1a mismatches are 0, the H1b gate does not fail, and κ and ε*_all match the report.
- No "(marginal)" label exists anywhere, so none was dropped.
- S2: `y.bin` equals `y_newton.txt`. The Newton step is at most 2.09e-6 relative to ‖y‖. μ_min changes by at most 6.6e-5 relative, and μ_min > τ_H at both points, so "margin positive" holds.
- NV6 not applicable: correct.

**4. Decision rules.**
- **D3 P1-A vs P2-A.** At 6:1, "same point" holds: ‖Δy‖/‖y‖ = 3.12e-6 and |Δq| = 2.7e-5. At 4:1, "P1 not converged" is correct.
- **D3 per formulation.** "failed (not converged)" for A and B is correct, apart from the F1 wording.
- **D5w.**
  - My windows reproduce Δ_w exactly: 4:1 W = 20 gives 23.39 / 60.88 / 36.54 / 23.85; 6:1 gives −11.12 / 0.47 (W = 20) and −12.81 / 0.20 (W = 25).
  - Horizons: k_c = 91 (4:1) and 56 (6:1).
  - The 6:1 padding label is correct (n_ext = 1, m = 2 at W = 25). "(no study-4 refinement check)" (6:1), "(refinement-limited: 1 of …)" (4:1) and "(not converged)" are all correct.
  - Reconstruction: Δ_0 is 1.6168 against 1.617, and 1.1457 against 1.146.
  - The P5 route is correctly not licensed.
- **D4.** u_w and the R(8:1)/R_w ratios reproduce. Every final's assessed worst point lies in Dw (maximum |x*| = 0.988 and r* = 2.41, against |x| ≤ 3.4 and r ≤ 3.22, using R0 = 2 from `s4.cpp` l. 637 and 792). Note that `analyze_runs.py` has no branch for the "outside Dw" label, so this was checked only by hand; I confirmed it.
- **D6x.** "budget-dependent (not converged)" follows from the precedence rules applied exactly.
- **T labels and assessed iterates.** Assessed q equals in-loop q to 6 significant figures at all 37 assessed iterates. Seeds 1 and 2 agree on ρ*_D at every final and start. All 10 NV4 step-0 rows match `v4_eps.csv`.

**5. The eight declared deviations.**

| # | Deviation | Ruling |
|---|---|---|
| 1 | New analysis scripts | Legitimate; see F4 |
| 2 | v3 limitation wording in `decisions.json` | Legitimate; cite the v4 wording only |
| 3 | CRLF re-run | Legitimate and verified: the list is now LF with 37 lines; there are no CRs in any CSV; there are 37 unique `iters.csv` rows; every snapshot equals its `ys.txt` line k+1; the failed run's log shows the "25\r" arithmetic errors and segfaults and no rows; no code was changed |
| 4 | X4L1fromL2 max-q iterate not assessed | Legitimate. It is a pre-registered `select()` defect (the argmax includes iterate 0). No verdict uses assessed q there, and in-loop and assessed q agreed at all 37 assessed iterates. Recommend assessing iterate 1 in Part 3 (about 14 s) |
| 5 | P4 not run | Legitimate under §5; now superseded (F5) |
| 6 | D4 descriptive for X4L1fromL2 | This is the registered scope ("on each CONVERGED A run (P1, and P2 separately)"), not a deviation |
| 7 | One FAILED line | Incorrect count (F4) |
| 8 | Other machine use | Legitimate |

None of the deviations changes a conclusion.

## What Part 2 licenses for the manuscript

**Licensed** (with full scope: Eulerian ρ, F2A, start-S1 lineage, fixed base node sets, regularised and discretised F):
1. At λ2, dense BFGS initialised from the measured Hessian reached first-order convergence (C1'∧C2') at A 4:1 (A4L2f) and A 6:1 (A6L2m, the same point as A6L2f). The second-order label is "strict local minimiser (second-order)" with μ_min/τ_H ≥ 1.4e4, quoted with:
   - the NV3 limitation (ε*_all 8.3e-6 to 1.6e-5);
   - the fact that H_R alone is indefinite there (F3).

   The same holds at λ1 for the continuation X4L1fromL2.
2. These A λ2 minimisers have deep assessed worst points: q 62.6 at 4:1 and 27.4 at 6:1. D4 labels them "fixed-node exploitation", which must be given with:
   - the registered post-hoc wording;
   - the "(start on calibration design)" label where it applies;
   - the F6 caveat.
3. D6x "budget-dependent (not converged)": the λ2-pocket start lost its pocket at λ1, and the λ1 shallow start reached the λ2 4:1 point.
4. Descriptive only, not a decision rule:
   - The A 4:1 point is reached from three different starts (‖Δy‖/‖y‖ ≤ 9.3e-7), but two of those runs are not converged under the registered rule.
   - At A4L2f's converged final, which is bit-identical to XB4L2fromA iterate 72, B's refinement rule requested 6 splits. So A's 4:1 minimiser is not a fixed point of B's refinement.

**Not licensed:**
- Any A-versus-B statement ("B lower", "B shallower", "suppresses") at either density. This includes the W = 20 "B lower / gap grows" and B4L2's final q of 4.4. B4L2 is not converged, and at 6:1 B ended at q 27.7, about the same as A.
- Any claim that B converged.
- Minimiser claims about R, about the unregularised or continuum problem, or global ones.
- "Strict" without the limitation sentence.
- Treating any of the three LINESEARCH_FAIL runs as converged.
- λ dependence (D6), λ1 D5w, or "multiple minima".
- That D4 has been validated, or that a D4 label implies a pocket.
- Convergence-rate claims (the 20-iterate slopes are about 0).
- Timing claims beyond "indicative".
- Anything from P4 until Part 3 has been reviewed.

## Reproduction (study4 directory; `C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I`)
- Hashes: `sha256sum PROTOCOL*.md; sha256sum -c logs/driver.sha256; sha256sum results/*.json results/*.csv assess/*.csv`
- Vera's scripts (session scratchpad, not saved in study4; may not persist): v1.py (C1'/C2' and terminal facts), v2.py (spectra, τ_H, S2, same-point checks), v3.py (assessed q, Dw, u_w, assessed vs in-loop q), v4.py (D5w windows). Each is run as `python -I <script> "$PWD"`.

## Limitations
- I did not rerun the optimiser, the Hessians or the assessments. Bit-reproducibility of the runs is not tested here.
- I checked D7, admissibility and the outside-Ab claims from logs only.
- The `v4_extra.py` pre-fix content cannot be verified (F4).

## Candidate lessons (proposed only)
1. Diff pre-registered analysis code against the protocol text, rule by rule, before any data exist. Here the `across_W` code used the gap rule for the verdict (F2). Scope: pre-registered decision scripts. Evidence: d5w.py l. 50–53 against §8 l. 575. Limitation: found in one study.
2. When a hashed script is fixed, archive the pre-fix copy, as was done for `fd_checks_v4` but not for `v4_extra`.
3. On Windows, Python list files read by bash must be written with `newline="\n"`.

**Next action:** the builder (Ada) applies corrections F1–F8 to Part 2 as a dated corrections section, without changing any results. The main session then sends Part 3 (P4) for separate review.

## Correction note (2026-10-28, from Vera's Part 3 review, G8)
Verification item 5, deviation 3, says "there are no CRs in any CSV". That is wrong: the s4/Python CSV and JSON outputs are CRLF throughout (for example, `assess/assess.csv` has 18 CRs in 18 lines). Effect: none. The claim that mattered, that `assess_list.txt` is LF, holds, and the readers tolerate CRLF (the mirror re-run reproduces the Part 2 records identically).
