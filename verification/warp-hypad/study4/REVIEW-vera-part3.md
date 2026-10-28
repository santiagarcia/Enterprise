# Independent review of REPORT-v4.md: "Part 2 corrections" and "Part 3" (P4 at λ1) (Vera, 2026-10-28)

Saved by the main session from Vera's hand-back (read-only review). Text unchanged except formatting.

**Verdict: CLEARED WITH CORRECTIONS.** No error changes a registered classification, a decision-rule label, or a Part 2 result. My recomputation of every Part 3 number in scope matches the raw outputs, with one descriptive exception (G5). The corrections concern wording, label completeness and traceability. One of them (G1) matters for the manuscript: the registered λ1 4:1 D5w "B lower" gloss calls the difference "inherited", which is wrong. The start gap was about zero, and the difference appeared during study 4.

**Role and independence**
- I am the reviewer, not a builder. I did not write Part 2, Part 3, any driver or analysis script, or any run.
- Possible bias:
  - I wrote F1–F8, which these corrections apply. I therefore have a stake in seeing them applied.
  - F5 recommended the `analyze_runs.py` CRLF fix whose provenance I am checking here. I recommended it but did not write it.
- Anchoring: I read the report first. Every check below uses my own scripts on the raw `trace.csv`, `termination.txt`, `outer.csv`, `assess/*.csv`, `logs/*.log` and Hessian `.bin` files. I used builder code only to compare its outputs with the protocol text.
- Affect enabled. Protective concern was moderate because these sections feed the manuscript. Its only effect was that I checked the A-vs-B licence more carefully. It did not change the evidence threshold.

## Findings, ranked by severity × likelihood × consequence

**G1 (medium for the manuscript; certain; no rule outcome changes). The licence gloss for λ1 4:1 contradicts the data.**
- P3-4 reports the registered label "B lower at matched cumulative iterates (start gap inherited) (determined by converged final(s)) (unequal prior history)". It then glosses it as "persistence of an inherited difference, with no change detected", and in the same paragraph says the start gap is near zero. Those two statements contradict each other.
- What I measured:
  - Δq at iterate 0 is −0.0063 in-loop and −0.0063 assessed. A was slightly lower at the start.
  - A moved from 0.595 to 0.723, and B from 0.601 to 0.577. So the +0.146 difference appeared during study 4.
  - Every window with Δ_w ≥ +0.05 contains extended (stationary) values. The two windows that are not extended give +0.002 and +0.028 at W = 20.
- The gap companion is undefined at λ1, so "change not detected" is true only because change was never tested. "Persistence of an inherited difference" is false.
- Correction:
  - Quote the registered label verbatim, followed by Δq(0) = −0.006 in-loop and assessed.
  - State that the difference is between two converged shallow finals: 0.723 against 0.577, both "none" (q < 1.01).
  - Drop "persistence of an inherited difference" in the report and in any manuscript text.

**G2 (low-medium; certain). Labels and a required quantity are missing for λ1 D5w.**
- (a) The 6:1 verdict needs its labels. n_ext = 5 of m = 7 at W = 20 and 5 of 6 at W = 25, which is at least m/2, so §8 N2a applies. The full label is "not established (determined by converged final(s))", plus "unequal prior history" (§5).
- (b) §8 requires Δq at iterate 0 "(in-loop and assessed)" alongside the verdict. Only the in-loop 4:1 value is given. The assessed values are −0.0063 (4:1) and +0.0055 (6:1). The 6:1 in-loop value (+0.0055) is also missing.

**G3 (low-medium; certain; risk downstream). `part3/results/decisions.json` contains λ1 strings that must not be quoted.**
- They are computed with the λ2 constants, as P3-6.2 says:
  - `L1_4.licensed` = "…; gap narrowed against both baselines …";
  - `L1_4.gap` = "gap shrinks";
  - `L1_6.gap` = "not established (baseline-sensitive)";
  - `reconstruction_deviation: true` at both densities.
- P3-6.2 correctly says these are not used. The JSON itself carries no warning, though, and the same JSON holds the registered labels.
- Correction: add one line to P3-7 naming these fields as invalid. Better still, add a short sidecar note next to the file, as a new file rather than an edit of the hashed output.

**G4 (low; certain). D3 at λ1 is a registration reading, and its labels need context.**
- v4 §8 registers D3 "per formulation, at λ2 on the P1 runs". D3 at λ1 is implied only through D6 ("if D3 verdicts exist at both λ").
- Mitigation: the λ1 D3 branch was already in `v4_extra.py`, hashed `35dc5a3e` at 23:40 on 10-27, before any P4 data existed. So the choice was not made after seeing data.
- Correction:
  - Declare it as a reading.
  - Attach the T label ("none along the trajectory") to the D3 sentence, as §8 D7/T requires.
  - State plainly that A's "resolution-dependent" comes from two "none" finals (0.723 against 0.570, |Δq| = 0.1525 ≥ 0.05). It must not be read as concentration that depends on resolution.

**G5 (low; certain; descriptive). The final-to-final distance range in P3-2 is wrong as written.**
- P3-2 says "The P4 finals otherwise differ from each other by ‖Δy‖/‖y‖ = 0.14–0.18".
- That range covers only the same-density and same-formulation pairs:

  | Pair | ‖Δy‖/‖y‖ |
  |---|---|
  | A4–B4 | 0.141 |
  | A4–A6 | 0.151 |
  | B4–B6 | 0.166 |
  | A6–B6 | 0.178 |

- The cross pairs are closer: A4L1f–B6L1f is 0.088, and **B4L1f–A6L1f is 0.036**.
- Correction: give 0.036–0.178, or name the pairs.

**G6 (low; certain). The time in C4 for the third FAILED line does not match the file times.**
- C4 says "From 07:27 there are three".
- `driver/analyze_runs.py` mtime is 07:30:48.4. The `bcdbf477` line was appended at 07:31:00 (`logs/driver.sha256` mtime 07:31:00.49; `execution.log` 07:31:00).
- Correct the time to 07:30:48 or 07:31. The order (hash logged, then first use) is not affected; see check 2.

**G7 (low; process). Part 2's byte-identity cannot be fully verified.**
- No hash of REPORT-v4.md was logged at the Part 2 append (07:21:15). The reference for the "prefix compare" in `execution.log` (10:12:01) is not archived.
- Indirect evidence that Part 2 is unchanged:
  - Every line reference in REVIEW-vera-part2.md still points to the same content: l. 111, 132, 137–138, 181–182, 214, 219–220, 234–236, 240, 250–251, 266, 268, 282–290, 308–313, 349–350 and 356.
  - The seven result hashes are unchanged.
  - Re-running the fixed `analyze_runs.py` on the mirror reproduces every Part 2 run record and both λ2 D5w records exactly. The only difference is the newly assessed X4L1fromL2 T_maxq_assessed: 29.47, previously NaN.
- I consider Part 2 unchanged on this evidence, but not proven byte-identical.
- Recommendation: log `sha256sum REPORT-v4.md` at every append.

**G8 (informational; my own error in REVIEW-vera-part2.md).**
- I wrote "there are no CRs in any CSV". That is wrong: the s4/Python CSV and JSON outputs are CRLF throughout (for example, `assess/assess.csv` has 18 CRs in 18 lines).
- Effect: none. The claim that mattered, that `assess_list.txt` is LF, holds. The readers tolerate CRLF, which the identical mirror re-run demonstrates.
- Please append a correction to my Part 2 review record.

## Verification results (all pass unless noted)

**1. Part 2 integrity and the mirror**
- The seven Part 2 hashes are unchanged: decisions.json 1c0a047e, v4_decisions.json f8e4be07, runs.csv 11477a51, v4_eps.csv ef92d174, assess.csv 4283d5b7, starts.csv f6ccda6b, iters.csv f399a43c.
- The Part 3 hashes match P3-7: 23332480, e270947a, c9382a42, 5b1e9059, 423295a1, 8575cdfe.
- Mirror:
  - 0 mismatches over all copied run files.
  - `part3/assess/assess.csv` and `iters.csv` begin with the exact Part 2 bytes, followed by 8 and 11 new lines.
  - `starts.csv` and the Hessian summaries are identical.
- Files modified after 07:22:07 outside `part3/` and the new P4 directories are only these: `driver/analyze_runs.py` (+ `.pyc`), `logs/driver.sha256`, `logs/execution.log`, `logs/post_p4.log`, `runs/runs.txt` (appended), REPORT-v4.md and REVIEW-vera-part2.md. No Part 2 Hessian, run, assessment or result file was touched.

**2. The fixed script: hashes and timing order**
- `analyze_runs_prefix.py.txt` has sha256 `a033a541…`, which is the registered hash. Its mtime (21:36:37 on 10-27) was preserved from the original.
- The diff against the current file is exactly one line (l. 82, adding `newline="\n"`). There is no CR difference.
- Order of events: the current file `bcdbf477…` was written at 07:30:48, its hash logged at 07:31:00, and it was first used (`select`) at 09:33:58 (`post_p4.log`). P4 output was first analysed only after P4 ended at 08:59:12, so the hash was logged before use.
- `p4_post.sh` is `e2b07e85…` and `p4_s2.sh` is `15cec41b…`, both matching `execution.log`. Both were logged before their first runs (08:59:23 and about 09:36:52).
- `sha256sum -c logs/driver.sha256` gives exactly 3 FAILED lines: analyze_runs, fd_checks_v4 and v4_extra, as declared.

**3. Classification**
- Method: my own C1' (‖∇F‖ ≤ 1e-5·F/max(‖y‖,1)) and C2' (10-iterate window on one node set) from `trace.csv`.
- My flags match the logged C1/C2 at every iterate of all four runs (0 mismatches).

| Run | First C1'∧C2' | C1' first true | C2' quantity at final | ‖∇F‖/tol | Check | Terminal |
|---|---|---|---|---|---|---|
| A4L1f | 79 (= final) | 5 | 6.31e-9 | 0.0188 | – | CONVERGED ✓ |
| B4L1f | 49 (= final) | 44 | 3.10e-9 | 0.0063 | it 49 PASS, n_need 0, C3 = 1 | CONVERGED ✓ |
| A6L1f | 41 (= final) | 36 | 2.68e-9 | 0.0121 | – | CONVERGED ✓ |
| B6L1f | 49 (= final) | 45 | 7.52e-9 | 0.0092 | it 49 PASS, C3 = 1 | CONVERGED ✓ |

- Both B runs converge at iterate 49. This is a genuine coincidence in the traces, not an artefact.
- B warm-start sets meet the §4.3 acceptance tests exactly: 8304 nodes / 519 leaves / ρ_min,w −5.6904060374e-4, and 16632 / 462 / −5.63314357538e-4, matching the study-3 last rows.
- Solver counts: precision-limited 0 (max ê_g/tol 1.2e-10); STALLED 0; one set per run; rates 0.185, 0.065, 0.048 and 0.060.
- Wall time: 1.61 h, 0.61 h and 0.55 h reproduce from the logs. Maximum wall/cap is 0.24.

**4. D1b (my own `eigh` on sym(H_R_fwd) + λ1·I; τ_H = max(1e-9·max|μ|, 10·max(e_H, e_H^c)))**
- `y.bin` equals `y_final.txt` bit for bit, and λ = 3.75e-8 in the metadata.

| Final | neg / nz / marg | μ_min | τ_H | μ_min/τ_H | κ | λ_min(H_R)/λ | n neg(H_R) | μ_min/λ |
|---|---|---|---|---|---|---|---|---|
| A4L1f | 0/0/0 | 3.0913e-8 | 3.282e-13 | 9.419e4 | 1.062e4 | −0.176 | 3 | 0.82 |
| B4L1f | 0/0/0 | 1.3541e-8 | 2.318e-13 | 5.840e4 | 1.712e4 | −0.639 | 3 | 0.36 |
| A6L1f | 0/0/0 | 2.2692e-8 | 2.441e-13 | 9.295e4 | 1.076e4 | −0.395 | 5 | 0.61 |
| B6L1f | 0/0/0 | 1.7251e-8 | 2.343e-13 | 7.362e4 | 1.358e4 | −0.540 | 4 | 0.46 |

- All four are strict and none is marginal.
- H_R is indefinite at all four, so positive definiteness of H_F rests on the λI term. At B4L1f, μ_min(H_F) is only 0.36λ.
- H1a mismatches and pair errors are 0, and the H1b gate does not fail.
- ε*_all values (7.7e-5, 4.5e-5, 7.2e-5, 5.7e-5) match `part3/results/v4_eps.csv`.
- S2: `y.bin` equals `y_newton.txt`. The step is ≤ 6.2e-8 relative to ‖y‖, and μ_min changes by ≤ 1.5e-6 relative. "Margin positive" holds at all four.

**5. Assessed q, D3, D4, T, D7**
- Final q: 0.7228, 0.5768, 0.5703 and 0.5708. Seeds 1 and 2 agree exactly. r_E = 0.575 at all four; every final is "none".
- D3 at λ1:
  - A: |Δq| = 0.1525, both < 1.01, so "resolution-dependent" (literal rule) ✓.
  - B: |Δq| = 0.0060, so "no concentration at convergence" ✓.
- Study-3 D2 is not carried into v4. Had it been applied, it would not trigger: r_E 0.575 ≥ 0.40.
- D4 (descriptive):

  | Run | u_w start → final |
  |---|---|
  | A4L1f | 4.56 % → 21.67 % |
  | B4L1f | 1.45 % → 0.45 % |
  | A6L1f | 1.59 % → 0.24 % |
  | B6L1f | 1.00 % → 0.77 % |

  - ‖Δy‖/‖y_start‖, Δq and R(8:1 Dw)/R_w reproduce.
  - Worst points lie in Dw (|x| ≤ 1.056, r ≤ 2.357).
- T label "none along the trajectory" ✓ for all four (maximum q 0.857, 0.699, 0.690, 0.986, all < 1.01). Assessed q equals in-loop q at all 11 new assessed iterates.
- X4L1fromL2 iterate 1: 29.47 assessed and in-loop, so Part 2 deviation 4 is closed ✓.
- D7 events: A4L1f final (seed 1) and A4L1f @75 only ✓. Outside-Ab mismatches 0 ✓.

**6. Decision-rule labels, checked against the registered text**
- D1: all four converged; B via PASS ✓.
- D5w 4:1 λ1 (my own implementation):
  - W = 20: +0.0021, +0.0277, +0.1557, +0.1460 ×4, so 5 of 7 windows ≥ 0.05, none ≤ −0.05 → "B lower".
  - W = 25: +0.0047, +0.0707, +0.1461, +0.1460 ×3, so 5 of 6 → "B lower".
  - No window lies within 0.001 of a threshold.
  - k_c = 150. n_ext 5/7 and 5/6, so "(determined by converged final(s))" ✓.
  - No refinement label, under the PASS exception ✓. "(unequal prior history)" ✓. No "(not converged)" ✓.
- D5w 6:1 λ1: every |Δ_w| ≤ 0.0025 → "not established" at both W ✓. The labels are incomplete (G2).
- λ level: none, because the densities differ ✓.
- Gap companion at λ1: "not defined" is the only coherent reading. v4 registers s0 = 84/90 and the baselines for λ2 only, and P4's unequal s0 (245/229 against 175/194) has no common pre-window. Mapping it to the "not evaluable" licence row is a reading, but every candidate row forbids "suppresses".
- D6: "not evaluable", because λ2 D3 is "failed (not converged)" ✓.
- D6x: "budget-dependent (not converged)" unchanged ✓.
  - The reference A4L1f is the same point as X4L1fromL2: ‖Δy‖/‖y‖ = 7.4e-7 and |Δq| ≈ 1e-6.
  - Neither arm keeps its start character ✓.
- "not established" and "not evaluable" are used as registered. No "(marginal)" or "(not converged)" label is missing in Part 3.

**7. Part 2 corrections C1–C8**
- C1: the registered §4.3 phrase is applied ✓.
- C2: the literal label "not established (window-sensitive)" is given, and the instruction not to quote it is included ✓.
- C3: the values match my earlier eigen-analysis ✓.
- C4: correct except the time (G6) ✓.
- C5: 72 s is confirmed (07:21:15 → 07:22:27) ✓.
- C6: ✓. P4 strengthens the caveat: A4L1f flips to "exploitation" (u_w 21.7 %) while q stays < 1.01.
- C7 and C8: ✓.
- No result was changed (check 1).

## Deviations

| # | Deviation | Ruling |
|---|---|---|
| P3-6.1 | `select()` newline fix | Legitimate and declared. Pre-fix copy archived with the registered hash; one-line diff; hash logged before use; Part 2 records reproduce identically under the fixed script |
| P3-6.2 | λ1 gap baselines not used | Legitimate and necessary. This is a gap in the registration rather than a choice. No conclusion depends on it, but see G3 |
| P3-6.3 | Mirror root | Legitimate and declared. Verified byte for byte (check 1) |
| P3-6.4 | D4 at P4 reported descriptively | Legitimate. This is the registered scope |
| P3-6.5 | S2 on B's final set | Consistent with §7. Not really a deviation |
| — | 3 FAILED lines | Verified |

None of these changes a conclusion.

## What Parts 2 and 3 together license for the manuscript

All statements below carry the full scope: Eulerian ρ, F2A, start-S1 lineage, fixed or refined node sets, and the regularised, discretised F.

**Licensed:**
1. **Part 2 items 1–4, unchanged** (REVIEW-vera-part2.md), with C1–C8 applied.
2. **λ1 second order.** All four P4 runs (A and B at 4:1 and 6:1) reach first-order convergence; B is converged via a PASS refinement check. Each final is a "strict local minimiser (second-order)" of the regularised F, with μ_min/τ_H from 5.8e4 to 9.4e4. Each must be stated with:
   - the NV3 limitation (ε*_all 4.5e-5 to 7.7e-5);
   - the fact that H_R alone is indefinite (λ_min(H_R)/λ from −0.18 to −0.64);
   - T "none along the trajectory".

   All are shallow (q 0.570–0.723; "none").
3. **D3 at λ1** (declared reading, G4). A: "resolution-dependent", from two "none" finals. B: "no concentration at convergence".
4. **D6x** "budget-dependent (not converged)". Descriptively, the A 4:1 λ1 point is reached both from the λ2 pocket and from the λ1 shallow start.
5. **D5w λ1, 4:1, quoted as a label only.** "B lower at matched cumulative iterates (start gap inherited) (determined by converged final(s)) (unequal prior history)". It must be accompanied by:
   - start Δq ≈ −0.006;
   - the statement that the difference is between two shallow converged finals (0.723 against 0.577) and appeared during study 4 (G1).

   6:1: "not established (determined by converged final(s))". There is no λ-level statement.

**Not licensed:**
- **"Suppresses", in any form, at either λ.**
  - λ2: D5w is not evaluable or window-sensitive, and the P5 route is not met.
  - λ1: the gap companion is undefined. Also, no λ1 final has a pocket, so there is no deepening to suppress.
- **Any general A-versus-B or formulation-effect claim.**
  - The only "B lower" holds at one density and one λ, and comes from the converged finals.
  - The histories are unequal, and the difference lies below the concentration threshold.
  - 6:1 shows no difference (|Δ| ≤ 0.003).
  - At λ2 nothing is licensed. This includes B4L2 q 4.4 and the W = 20 intermediates.
- "Persistence of an inherited difference" at λ1.
- Any λ1 gap-change statement, and any quotation of the λ1 `gap`/`licensed`/`reconstruction_deviation` fields in `part3/results/decisions.json`.
- **Attributing A4L1f's higher q to fixed-node exploitation.**
  - Speculative hypothesis, not licensed: A 4:1 is the outlier on its coarse base set. B4L1f and A6L1f end 0.036 apart with q 0.577 and 0.570, and A4L1f has u_w 21.7 %.
  - D4 is descriptive here and not validated out-of-sample (C6). A test that could discriminate: compare A4L1f's assessed worst point and F with a 6:1-node evaluation of the same y. It is optional and not needed for any licensed statement.
- λ dependence (D6), "multiple minima", and claims about the unregularised, continuum or global problem.
- "Strict" without the NV3 and H_R-indefinite sentences.
- Rate or timing claims beyond "indicative".

## Acceptance conditions
- The builder adds a dated "Part 3 corrections" section applying G1–G6 without changing any result file.
- The main session appends the G8 correction to REVIEW-vera-part2.md as a dated note.
- G7 is a recommendation for future appends, not a blocker.

No re-review is needed if the corrections are wording-only and the six Part 3 hashes plus the seven Part 2 hashes are unchanged.

## Reproduction (from `warp_hypad/study4`; `PY = C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I`)
```
sha256sum -c logs/driver.sha256        # 3 FAILED (analyze_runs, fd_checks_v4, v4_extra)
sha256sum results/*.json results/*.csv assess/*.csv part3/results/* part3/assess/*.csv part3/*.sh driver/analyze_runs*.py*
diff driver/analyze_runs_prefix.py.txt driver/analyze_runs.py      # one line (82)
find . -newer logs/p4_prerun_hashcheck.txt -type f | grep -v -E "^./part3|L1f"
```
Vera's scripts p3a.py (hashes, mirror fidelity, C1'/C2', in-loop q, own D5w windows), p3b.py (own eigh/τ_H/κ/λ_min(H_R), S2, assessed q, D3, D4, same-point distances) and p3c.py (decision-record equality; λ1 D5w JSON fields) are in the session scratchpad, which may not persist. They are not saved in study4.

## Limitations
- I did not rerun the optimiser, the Hessians or the assessments, so bit-reproducibility of the runs is not tested.
- I took metric validity (21,125 points) and the outside-Ab checks from the logs only.
- Part 2 byte-identity is supported only indirectly (G7).
- `council assemble` initially failed on a cp1252 console-encoding error. It succeeded with `PYTHONIOENCODING=utf-8`, and I read its output. No relevant adopted memories were retrieved.

## Candidate lessons (proposed only)
1. **A registered licence gloss can contradict the data it licenses.**
   - Scope: pre-registered decision tables whose wording embeds an assumption ("start gap inherited").
   - Evidence: the λ1 4:1 D5w, where Δq(0) = −0.006.
   - Practice: check each gloss against the measured premise before quoting it.
   - Limitation: one study.
2. **Hash the report at every append.** Otherwise a "byte-identical" claim cannot be checked later (G7). Scope: append-only reports.
3. **Reviewer self-check.** Verify blanket "none/no X anywhere" statements with a count before writing them (G8). Scope: Vera's own reviews.

**Next action:** Ada applies G1–G6 as a dated "Part 3 corrections" section. The main session appends the G8 note to REVIEW-vera-part2.md. Quill then drafts manuscript text only within the licence above, and a separate Vera pass checks the draft wording against it.
