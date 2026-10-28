# Study 4 REPORT under PROTOCOL-v4: resumption stopped at condition R3 (cross-implementation check) for review

Builder: Ada, 2026-10-27 (23:30–23:45 CDT). These are builder results and have **not been independently reviewed**.
Saved by the main session from Ada's hand-back (the subagent's write was blocked). `REPORT.md` is kept unchanged as the
v3 record. Scope: Eulerian ρ only; F2A, start-S1 lineage, λ1 = 3.75e-8 and λ2 = 3.75e-9 (PROTOCOL-v4 scope).

## 1. Outcome
- **Integrity checks pass.**
  - `PROTOCOL-v4.md` SHA-256 `0d77d7d5…c3ad` matches `logs/protocol-v4.sha256`, and v1–v3 are unchanged.
  - `driver/s4.cpp` is `ef247445…9463`, and every file in `logs/driver.sha256` checks OK.
  - The five archived `hessians/fd/*.out` hashes equal the v4 §9 table.
  - Evidence: `logs/v4_resume_hashcheck.txt`.
- **Vera's conditions R1–R3 were declared and logged before the analysis script existed**
  (`logs/v4_conditions_R1R3.txt`, SHA-256 `db1c4679…6bfc`; `logs/execution.log` 23:31:12).
- **R1 passes.** All 20 inputs pinned in REVIEW-vera-v4.md §4 match.
- **The V4 substance passes under v4 with R2.**
  - H2 (NV1) passes in all 8 columns.
  - H3 (a), (b) and (c) pass at all 5 designs.
  - At A_6_L2, (c) is "slope not resolvable by the R reference".
- **R3 is not met.** 84 of 86 compared items reproduce Vera's §1 table to the printed digits; the 2 that do not are both
  at B_4_L2 (§3).
  - R3 as declared says: "any disagreement is reported and the study STOPS for review".
  - **The study is stopped.** Gate V4 is not declared passed.
  - No P1–P5 run, Hessian refresh, assessment, D1/D1b/D3/D4/D5w/gap/D6/D6x/D7 quantity or (marginal) label was computed.

## 2. V4 re-analysis (`driver/fd_checks_v4.py`, run 2; `logs/gate_H2H3_v4.txt`; `hessians/fd/fd_checks_v4.csv`)
In every row, h₁…h₄ = 1e-2…1e-5. All five H3 designs used s = h, since |v·y| = 0.581, 0.877, 0.493, 0.196 and 0.455.

| Row | Result |
|---|---|
| H2, A_4_L1 (cols 0, 36, 72, 69) | PASS. One pair of ~100× in each column. Best e = 9.1e-9, 4.3e-9, 1.3e-8 and 6.4e-9. Two-pair reading (descriptive): col0 fails. |
| H2, A_4_L2 (cols 0, 36, 72, 51) | PASS. ~100× in each column. Best e = 3.2e-8, 1.7e-8, 1.8e-8 and 1.7e-7. Two-pair reading: col0 fails. |
| H3 (a) Hv | PASS at all 5 designs. Ratios ~100×. Best e_a = 1.27e-8 (A_4_L1), 6.37e-8 (A_4_L2), 3.41e-8 (B_4_L2), 7.12e-9 (A_6_L2), 5.42e-8 (B_6_L2). |
| H3 (b) R–∇R | PASS at all 5 designs. First ratio 100–103×. Best c/\|φ′\| = 8.2e-7, 4.0e-7, 1.4e-8, 1.2e-7, 6.6e-6, all ≤ 1e-4. |
| H3 (c) vᵀHv | PASS at all 5 designs. Min e_c = 1.13e-6, 3.35e-6, 2.12e-6, 2.48e-6, 2.81e-6. Resolvable and slope ≥ 30× (45.5, 55.3, 153.9, 71.6) at A_4_L1, A_4_L2, B_4_L2 and B_6_L2. **A_6_L2: "slope not resolvable by the R reference"**, with ρ(s₂) = 1.0e-5 > e_c(h₁)/30 = 5.6e-7 and a ratio of 6.7. |
| Descriptive | Richardson r(s₁) = 0.49996, 0.49984, 0.49981, 0.49949, 0.49991 (theory ½). In every row the minimum lies at h₂ or below, so R2 did not bind. |

**Carried over (v4 §9 step 5) without regeneration:** H1a, H1b, H1c, H4, V1–V3, V5a–V5e and V6 (REPORT.md §2).

**Strength of evidence, using Vera's correction.** There are two independent lines of evidence, not three. An OTI error
along v is excluded at ≲ 3e-8 by the gradient differences and at ≲ 1e-7 by the R values.

**Coverage limit (Vera F3).** One random direction plus 4 columns at 2 designs cannot detect an error with Ev = 0.
Structured mapping errors are covered by H1a, V5c and H2.

## 3. R3 disagreements (reason for the stop)
| Design | Item | Mine (full precision) | Vera (printed) | In R3's listed items? |
|---|---|---|---|---|
| B_4_L2 | best e_a (h₄) | 3.414538e-8, which prints as 3.41e-8 | 3.42e-8 | **yes** |
| B_4_L2 | signed projected-gradient error, best (h₄) | −1.366e-10 | 1.4e-10 | no (sign only) |

- **The e_a value is stable across inputs.** It is identical to 1e-8 relative with H_R_rev, H_R_comp or H_R_ld in place
  of H_R_fwd (3.4145379e-8 in every case), so the difference is not caused by which Hessian file is used.
- **Possible explanation, not verified.** The value agrees with Vera's if it was rounded twice
  (3.4145 → 3.415 → 3.42).
- **Effect on the gate.** Neither item affects any pass or fail: e_a = 3.4e-8 against a tolerance of 1e-4, and the
  projected error is descriptive.
- **What review is needed.** Vera, or another reviewer, rules on whether this is a printing difference. If it is, V4
  passes and P1 may start under v4 with no other change.

## 4. Protocol deviations and builder readings (declared)
1. **R1–R3 are pre-analysis conditions (Vera).** R2 was applied to the H2, H3 (a) and H3 (b) slope rules as well as to
   (c). Under H2 an empty pair set already fails, so this tightens nothing new.
2. **NV2 reading (ii) is lenient, not conservative (Vera F4).** Using ρ over k ∈ {3, 4} can only exempt pairs from the
   slope gate. It did not matter here, because ρ was set by h₄ at every design.
3. **A comparator code fix between runs 1 and 2. The analysis was not changed.**
   - Run 1 (script `a0c567be…`; archived as `logs/gate_H2H3_v4_run1.txt` and `driver/fd_checks_v4_run1.py.txt`) mapped
     Vera's elided cells "first … x" to h₄. Five of its seven mismatches came from that mapping. In Vera's table, x is the
     value at the best step; her A_4_L1 "(h₃)" annotation shows this.
   - The fix was made after the mismatches were seen. It was logged, and the new hash (`3c40f767…`) was recorded before
     run 2.
   - All gate numbers are identical in runs 1 and 2.
4. **H2 FD symmetry residual (descriptive).** It is NaN for pairs with an OTI entry H_R[i][j] = 0 exactly
   (A_4_L1: pairs with col69). It is not a gate.
5. **The script also checks** that the `.dirs` v equals `default_rng(4)` bit-for-bit, and that every step equals the
   registered h·max(1, ·).

## 5. Not run
- P1–P5, assessments, post-run Hessians and NV6: **not run**.
- Machine load was not checked, because no run started.

## 6. Reproduction (from `warp_hypad/study4`; `C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I`)
```
sha256sum PROTOCOL*.md hessians/fd/*.out; sha256sum -c logs/driver.sha256
python -I driver/fd_checks_v4.py "$PWD"     # writes logs/gate_H2H3_v4.txt, hessians/fd/fd_checks_v4.csv; exit 3 = R3 stop
```
Hashes:
- `logs/gate_H2H3_v4.txt`: `f6d89fea…caf5`
- `hessians/fd/fd_checks_v4.csv`: `b788c1aa…cd17`

---

# Part 2. Resumption after the R3 ruling: P1–P3 and P5 runs, assessments and decision rules (Ada, builder, 2026-10-27/28)

These are builder results. They have **not been independently reviewed**. Scope is as in PROTOCOL-v4:
- Eulerian ρ = T(n,n) only;
- family F2A;
- the start-S1 lineage;
- λ1 = 3.75e-8 and λ2 = 3.75e-9.

All statements concern the regularised, discretised F on the stated node sets.

## 7. Gate V4 and resumption
- **The R3 ruling is logged** in `logs/execution.log`.
  - Source: REVIEW-vera-v4.md, "R3 ruling (Vera, 2026-10-27)": "R3 SATISFIED (V4 passes; P1 may start under v4)".
  - Both B_4_L2 differences in §3 were transcription errors in Vera's §1 table. Erratum: e_a is 3.41e-8, and the
    signed projected error is −1.4e-10.
  - **Gate V4 passes under v4.** This ruling supersedes the statement in §1 that the study is stopped.
- **Hashes were re-verified before P1** (`logs/p1_prerun_hashcheck.txt`).
  - PROTOCOL-v4 (`0d77d7d5…c3ad`), `driver/s4.cpp` (`ef247445…9463`) and every file in `logs/driver.sha256` match.
  - The one FAILED line in `logs/driver.sha256` is the superseded run-1 hash of `fd_checks_v4.py` (`a0c567be…`).
    This failure is expected.
  - The five FD `.out` hashes match.
- **Machine load before P1:** 3–8 % total CPU on 24 logical CPUs. Every block used at most 12 threads.
- **Wall time.** P1 started at 23:38:48, and all work ended at 07:18. Total: 7.7 h, within the ~11 h worst case.
  - P1: 4.05 h (critical lane B4L2).
  - P2 + P3 + P5: 2.32 h (critical lane XB4L2fromA).
  - Post-run Hessians, assessments and S2 Hessians: 1.3 h.
  - No run reached its hang guard.

## 8. Runs (`runs/runs.txt`, `runs/<ID>/termination.txt`, `results/runs.csv`, `results/decisions.json`)
| Pri | Run | Terminal state | Iterates | ‖∇F‖/tol_C1' at final | ‖y‖ | Final assessed q | r_E | Checks (B) |
|---|---|---|---|---|---|---|---|---|
| P1 | A4L2m | LINESEARCH_FAIL (C1' true, C2' false) | 91 | 0.031 | 36.467 | 62.62 | 0.341 | – |
| P1 | B4L2 | ITERATION_LIMIT | 200 | 5.2e2 | 35.424 | 4.436 | 0.625 | 1 (iterate 70, SPLIT, 2 splits); n_outer 1 |
| P1 | A6L2m | **CONVERGED** | 49 | 0.064 | 35.013 | 27.38 | 0.619 | – |
| P1 | B6L2 | LINESEARCH_FAIL (C1' true, C2' false) | 56 | 0.0068 | 35.186 | 27.68 | 0.619 | none (n_check 0) |
| P2 | A4L2f | **CONVERGED** | 72 | 0.44 | 36.467 | 62.62 | 0.341 | – |
| P2 | A6L2f | **CONVERGED** | 66 | 0.12 | 35.013 | 27.38 | 0.619 | – |
| P3 | X4L1fromL2 (λ1) | **CONVERGED** | 126 | 0.13 | 16.695 | 0.723 | 0.575 | – |
| P3 | X4L2fromL1 (λ2) | LINESEARCH_FAIL (C1' true, C2' false) | 65 | 0.10 | 36.467 | 62.62 | 0.341 | – |
| P4 | A4L1f, B4L1f, A6L1f, B6L1f | **not run** (§12, deviation 5) | – | – | – | – | – | – |
| P5 | XB4L2fromA | ITERATION_LIMIT | 150 | 1.7e5 | 39.126 | 3.035 | 0.604 | 1 (iterate 72, SPLIT, 6 splits); n_outer 1 |

- **The three LINESEARCH_FAIL runs are not converged (D1).**
  - Each met C1' for its last 3–5 iterates.
  - The line search failed before the 10-iterate C2' window closed. At the last iterate, (F_{k−10} − F_k)/F_k was
    3.9e-8 (A4L2m) and 2.3e-6 (B6L2), against a threshold of 1e-8.
  - This is the registered state. No repair was made.
- **Descriptive: final designs that coincide.**
  - A4L2m, X4L2fromL1 and A4L2f end at the same design. Against A4L2f, ‖Δy‖/‖y‖ = 9.3e-7. F agrees to 12 digits, and
    assessed q is 62.62 for all three.
  - A6L2m and A6L2f end at the same design (‖Δy‖/‖y‖ = 3.1e-6).
  - B6L2 does **not** end at the A 6:1 design (‖Δy‖/‖y‖ = 6.9e-2), although its q (27.68) is close to A's (27.38).
- **B refinement.**
  - **B4L2.** Its first study-4 check, at iterate 70, split 2 leaves. On the new set F rose from 1.0678e-4 to
    1.5387e-4 (`runs/B4L2/set_changes.csv`). The second inner solve did not converge by K.
  - **XB4L2fromA.** It followed the A4L2f trajectory exactly until its first check: assessed q is identical at iterates
    25 and 50, because both start on the base set with the same W_0. It then split 6 leaves at iterate 72 and did not
    converge again by K = 150.
- **Precision and solver counts.**
  - No iterate in any run was C1' "precision-limited" (H1c). ê_g/tol_C1' ≈ 2e-11 at the converged finals.
  - STALLED iterates: 0 in every run.
  - Refreshes: 1–2 per run, all triggered by the gradient rule.
  - W resets: 0.
- **Rate (descriptive):** slope of log10‖∇F‖ per iterate over iterates 1–20.

  | Run | Slope |
  |---|---|
  | A4L2m | +0.018 |
  | B4L2 | −0.007 |
  | A6L2m | +0.020 |
  | B6L2 | +0.001 |
  | A4L2f | −0.004 |
  | A6L2f | +0.007 |
  | X4L1fromL2 | −0.011 |
  | X4L2fromL1 | +0.017 |
  | XB4L2fromA | −0.004 |

## 9. Second-order character (D1b), NV3, NV6 and S2 (`hessians/final_<run>/`, `hessians/s2_<run>/`, `results/v4_eps.csv`)
| Converged final | H_F label | n negative / near-zero / marginal | μ_min | τ_H | μ_min/τ_H | κ(H_F) | ε*_all | D1b |
|---|---|---|---|---|---|---|---|---|
| A6L2m | PD | 0 / 0 / 0 | 3.767e-10 | 2.69e-14 | 1.40e4 | 7.15e4 | 8.3e-6 | strict local minimiser (second-order) |
| A4L2f | PD | 0 / 0 / 0 | 2.946e-9 | 1.74e-13 | 1.69e4 | 5.90e4 | 1.6e-5 | strict local minimiser (second-order) |
| A6L2f | PD | 0 / 0 / 0 | 3.767e-10 | 2.69e-14 | 1.40e4 | 7.15e4 | 8.3e-6 | strict local minimiser (second-order) |
| X4L1fromL2 | PD | 0 / 0 / 0 | 3.091e-8 | 3.28e-13 | 9.42e4 | 1.06e4 | 7.7e-5 | strict local minimiser (second-order) |

- **No final is "strict (marginal)".** Every μ_min/τ_H is at least 1.4e4, far above the marginal threshold of 100.
  No "(marginal)" suffix arises anywhere in this part.
- **The H1a and H1b gates pass** at all four finals.
- **Limitation (v4, NV3), stated with every D1b label above:** "τ_H covers node-sum roundoff; a sign change by
  per-node evaluation error requires a coherent relative error ≥ ε*_all (reported)." Here ε*_all = 8.3e-6 to 7.7e-5.
- **NV6 is not applicable.** No converged final is labelled "degenerate or saddle", so no P_neg projection or
  F(y ± α·w_min) evaluation was made.
- **S2 (descriptive): H_F was re-measured at y + Newton step.**
  - ‖Newton step‖: 7.3e-5 (A6L2m), 3.4e-5 (A4L2f), 6.4e-5 (A6L2f) and 1.2e-5 (X4L1fromL2). Each is at most 2.1e-6
    relative to ‖y‖.
  - μ_min at the Newton point agrees with the final value to at least 4 digits.
  - μ_min/τ_H there is 1.40e4, 1.69e4, 1.40e4 and 9.42e4.
  - Result: **"margin positive"** at all four.
- **Step-0 labels in v4 wording (NV4), with ε* (NV3); descriptive.**
  - Every step-0 design is "indefinite (n negative), non-stationary (‖∇F‖/tol_C1' = x)":

    | Design | Negative eigenvalues | ‖∇F‖/tol_C1' |
    |---|---|---|
    | A_4_L1 | 1 | 44.7 |
    | A_4_L2 | 2 | 1.3e4 |
    | A_6_L1 | 1 | 689 |
    | A_6_L2 | 2 | 5.67e3 |
    | B_4_L1 | 1 | 4.87e3 |
    | B_4_L2 | 7 | 3.81e4 |
    | B_6_L1 | 2 | 836 |
    | B_6_L2 | 3 | 1.61e4 |
    | A_4_L2_it84 | 5 | 1.81e4 |
    | A_6_L2_it90 | 3 | 1.23e4 |

  - "Saddle" is therefore not used for any step-0 design.
  - ε*_min ranges from 1.6e-5 to 9.1e-4, and ε*_all from 2.0e-6 to 6.6e-5 (`results/v4_eps.csv`).
  - The other step-0 labels are unchanged from REPORT.md.

## 10. Decision rules
**D1 (converged).**
- Converged: A6L2m, A4L2f, A6L2f and X4L1fromL2.
- Not converged: A4L2m, B6L2 and X4L2fromL1 (LINESEARCH_FAIL); B4L2 and XB4L2fromA (ITERATION_LIMIT).
- No final is precision-limited.

**D3, per formulation at λ2 on P1:**
- A: **"failed (not converged)"**. A4L2m ended LINESEARCH_FAIL; A6L2m is CONVERGED and D1b strict.
- B: **"failed (not converged)"**. B4L2 ended ITERATION_LIMIT (n_outer 1, n_check 1); B6L2 ended LINESEARCH_FAIL
  (n_outer 0, n_check 0).

**D3, P1-A vs P2-A:**
- 4:1: **"P1 not converged"**. Descriptively, A4L2m's final is at ‖Δy‖/‖y‖ = 9.3e-7 and |Δq| = 4e-4 from A4L2f's.
- 6:1: **"same point"**. |Δq| = 2.7e-5 and ‖Δy‖/‖y‖ = 3.1e-6; both runs are CONVERGED and D1b strict, not marginal.

**D4 (fixed-node exploitation; registered on CONVERGED A runs of P1 and P2).**
- Post-hoc status, in the registered wording: "D4 statistic chosen post hoc (study-3 REPORT §6, Noether R5);
  retro-check on the motivating data is a consistency check, not a validation; P2 starts on the calibration designs
  (start u_w reported)."
- All assessed worst points lie inside Dw (|x| ≤ 0.99, r ≤ 2.41).

| Run | u_w start → final | Label | ‖Δy‖/‖y_start‖ | q_final − q_start | R(8:1 Dw)/R_w |
|---|---|---|---|---|---|
| A6L2m (P1, out-of-sample) | 72.9 % → 97.7 % | **fixed-node exploitation** | 0.50 | +25.4 | 1.027 (> 1.02 noted) |
| A4L2f (P2) | 98.2 % → 99.1 % | **fixed-node exploitation (start on calibration design)** | 0.48 | +30.7 | 4.23 (> 1.02 noted) |
| A6L2f (P2) | 87.9 % → 97.7 % | **fixed-node exploitation (start on calibration design)** | 0.43 | +22.6 | 1.027 (> 1.02 noted) |
| A4L2m (descriptive) | 74.6 % → 99.1 % | fixed-node exploitation (not converged) | 0.59 | +59.9 | 4.23 |
| B4L2 (descriptive) | 1.6 % → 84.9 % | fixed-node exploitation (not converged) (start on calibration design) | 0.58 | +3.38 | 1.003 |
| B6L2 (descriptive) | 0.3 % → 97.6 % | fixed-node exploitation (not converged) (start on calibration design) | 0.53 | +26.8 | 1.028 |
| X4L1fromL2 (descriptive, P3) | 98.2 % → 21.7 % | fixed-node exploitation (start on calibration design) | 0.86 | −31.2 | 1.001 |
| X4L2fromL1 (descriptive) | 4.6 % → 99.1 % | fixed-node exploitation (not converged) (start on calibration design) | 1.38 | +62.0 | 4.23 |
| XB4L2fromA (descriptive) | 98.2 % → 59.0 % | fixed-node exploitation (not converged) (start on calibration design) | 0.75 | −28.9 | 1.002 |

The only out-of-sample application of the registered rule is A6L2m, and it gives "fixed-node exploitation".

**D5w (matched iterate budget), λ2:**
- **4:1: "not evaluable (< 4 windows)"**, with the labels "(refinement-limited: 1 of 3 windows after the first
  check)" and "(not converged)".
  - The horizon is k_c = 91, set by A4L2m's LINESEARCH_FAIL.
  - At W = 20 there are 4 windows. Δ_w = +23.4, +60.9, +36.5 and +23.8, which gives "B lower" and "gap grows" on both
    baselines.
  - At W = 25 there are only 3 windows, so the result across W is "not evaluable".
  - **No level statement is made, and "suppresses" is not used.** k_chk1 = 70, and n_check within k_c = 1.
- **6:1: "not evaluable (< 4 windows)"**, with the labels "(determined by converged final(s))", "(no study-4
  refinement check)" and "(not converged)".
  - The padding label applies because n_ext = 1 of m = 2 at W = 25.
  - The horizon is k_c = 56, set by B6L2's LINESEARCH_FAIL. There are 2 windows at each W.
  - Δ_w = −11.12 and +0.48 at W = 20, and −12.81 and +0.20 at W = 25.
  - **No level statement is made.**
- **Gap companion:** "not evaluable (< 4 windows)" at both densities.
- **Reconstruction check.** The study-4 iterate-0 in-loop Δ is 1.6168 (4:1) and 1.1457 (6:1), against the archived
  Δ_0 of 1.617 and 1.146. Both differences are below 0.01, so there is no reconstruction deviation.
- **P5 route: not licensed.** It requires a D5w "B lower" verdict without a refinement label. Descriptively,
  XB4L2fromA's final has a pocket (q 3.03, not converged), as does A4L2f's (q 62.6). P5 therefore does not show a
  shallow B final.
- **λ1 D5w: not run**, because P4 was not run.

**D6 (λ dependence): "not evaluable".** There is no D3 verdict at λ2, and P4 was not run.

**D6x (cross warm start): "budget-dependent (not converged)".**
- X4L1fromL2 (λ1, started from the q = 32 λ2 pocket) ends shallow: q 0.72, CONVERGED, strict, ‖y‖ 16.7. It adopts λ1's
  character.
- X4L2fromL1 (λ2, started from the shallow λ1 final) ends with a pocket: q 62.6, ‖y‖ 36.5, LINESEARCH_FAIL (not
  converged). It adopts λ2's character.
- The suffix "(not converged)" is added because X4L2fromL1 determines the label and did not converge.

**T trajectory labels.** These count study-4 iterates only. The maximum in-loop q is confirmed by assessment where that
iterate was assessed.

| Run | T label | Iterates with concentration | Max in-loop q (iterate) | Assessed there | q range, last 20 (flag > 0.05) |
|---|---|---|---|---|---|
| A4L2m | concentration along the trajectory | 91 | 127.7 (30) | 127.7 | 61.5–64.0 (flagged) |
| B4L2 | concentration along the trajectory | 194 | 143.8 (117) | 143.8 | 4.44–11.8 (flagged) |
| A6L2m | concentration along the trajectory | 49 | 31.0 (24) | 31.0 | 27.3–28.2 (flagged) |
| B6L2 | concentration along the trajectory | 51 | 42.1 (20) | 42.1 | 26.7–27.9 (flagged) |
| A4L2f | concentration along the trajectory | 72 | 106.4 (41) | 106.4 | 62.5–66.8 (flagged) |
| A6L2f | concentration along the trajectory | 66 | 32.0 (29) | 32.0 | 27.1–27.6 (flagged) |
| X4L1fromL2 | concentration along the trajectory, not at the final iterate (converged) | 10 | 29.5 (1) | not assessed (deviation 4) | 0.718–0.726 |
| X4L2fromL1 | concentration along the trajectory | 63 | 110.9 (15) | 110.9 | 62.6–67.1 (flagged) |
| XB4L2fromA | concentration along the trajectory | 150 | 558.1 (105) | 558.1 | 1.34–10.3 (flagged) |

**D7 (search sanity; descriptive).**
- Sub-lattice D7 events occurred at:
  - B4L2 (both seeds);
  - XB4L2fromA (both seeds);
  - X4L1fromL2 (seed 1 only);
  - the starts A_4_L2_it84, B_4_L2 and B_6_L2.

  (A D7 event is a sub-lattice value below the pre-polish best.) ρ*_D includes these values.
- Seeds 1 and 2 agree on ρ*_D at every final and start. This is not cited as robustness.
- The outside-Ab independence check passed everywhere, so no whole-D fallback was needed.
- Metric validity holds at all 21,125 admissibility points at every assessed final and start.

**H (Question 1):** see §9 for the NV4 wording, and REPORT.md for the remaining labels. The study-3 κ statement and
Vera's I5 rewording refer to the study-3 finals at step 0, and their labels are unchanged from REPORT.md.

## 11. What the results license (builder reading; for review)
- **Converged second-order minimisers exist on the fixed base sets at λ2.**
  - A 4:1 (A4L2f) and A 6:1 (A6L2m, the same point as A6L2f) reach strict local minimisers (not marginal) of the
    regularised, discretised F.
  - Both have deep pockets: q 62.6 and 27.4.
  - The registered D4 rule labels both "fixed-node exploitation" (u_w ≈ 98–99 %).
- **No A-vs-B formulation statement is licensed.**
  - D3 failed (not converged) for both formulations.
  - D5w is not evaluable at either density.
  - B4L2's single refinement check produced a split, and the following inner solve did not converge within K. B6L2
    made no check.
- **D6x is "budget-dependent (not converged)".** The λ2-pocket start lost its pocket at λ1, and the shallow λ1 start
  gained the λ2 pocket.

## 12. Protocol deviations and builder notes (declared)
1. **New analysis scripts.**
   - `driver/v4_extra.py` computes NV3 ε*, the NV4 label, NV6, D3, D6 and D6x: the sections not covered by the
     pre-registered `analyze_runs.py`.
   - `driver/s2_newton.py` computes the S2 Newton point.
   - The hashes of both scripts were recorded in `logs/driver.sha256` before first use.
   - The first `v4_extra.py` run failed on a cp1252 encoding error before writing anything. The fix changed only the
     output encoding (to UTF-8). The script was re-hashed (`35dc5a3e…`) and re-run.
2. **D1b limitation wording.** `analyze_runs.py` was pre-registered under v3 and still writes the v3 wording into
   `decisions.json`. The registered wording is the v4 (NV3) text above, taken from `results/v4_decisions.json`.
3. **Trajectory assessment was re-run.**
   - `analyze_runs.py select` wrote `results/assess_list.txt` with CRLF line endings.
   - `assess.sh` therefore read the iterate as "25\r" and wrote empty snapshots. `s4 assess` then segfaulted on all 37,
     and no row was written.
   - The CRs were stripped from the list (the same 37 entries; no code change), and the assessment was re-run.
4. **The max-q iterate of X4L1fromL2 was not assessed.**
   - `select` takes the argmax over all iterates, including iterate 0, and then drops iterate 0.
   - `report` takes the maximum over iterates ≥ 1 with r_E ≥ 0.10, which is iterate 1 (in-loop q 29.5).
   - Iterate 1 was therefore not in the list. This affects descriptive output only.
5. **P4 was not run.**
   - At 06:40, 7.0 h into the session, I estimated the remaining cost from the measured B 4:1 lane rate in P1
     (73 s per iterate, including 2 refreshes). At that rate B4L1f (K = 150) needs about 3 h.
   - With P4's assessment and up to 4 second-order Hessians, the session would have reached about 11 h or more.
   - v4 §5 forbids truncating a started run for session time, so P4 was not started.
   - It is reported "not run". It may be run later under unchanged v4, with the date recorded.
   - Consequently λ1 D5w, D6 and the D6x "multiple minima" reference are unavailable.
6. **D4 scope.** The script prints X4L1fromL2's D4 result without "[descriptive]", because it is a converged A run.
   The registered rule covers P1 and P2 only, so the result is treated as descriptive here.
7. **`logs/driver.sha256`** keeps the superseded run-1 line for `fd_checks_v4.py`. `sha256sum -c` therefore reports one
   expected FAILED line.
8. **Other machine use.** Other user processes were running (3–8 % total CPU before P1). Results are thread-count
   independent (v4 §6), so timings are indicative only.

## 13. Reproduction (from `warp_hypad/study4`; `PY = C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I`)
```
sha256sum -c logs/driver.sha256                     # one expected FAILED line (deviation 7)
# runs: driver/lane.sh with the arguments logged in runs/<ID>/config.txt (P1, then the P2+P3+P5 block; <= 12 threads)
s4.exe hessian --y runs/<R>/y_final.txt --ws <d>:1 --lambda <L> --passes driver/hessian_passes.txt --out hessians/final_<R>
$PY driver/hess_analyze.py $PWD design hessians/final_<R>
bash driver/assess.sh finals <runs...>; bash driver/assess.sh starts; bash driver/assess.sh iters results/assess_list.txt  # list must be LF
$PY driver/analyze_runs.py $PWD select|report      # results/decisions.json, results/runs.csv
$PY driver/v4_extra.py $PWD eps <hessian dirs>; $PY driver/v4_extra.py $PWD decide   # results/v4_eps.csv, results/v4_decisions.json
$PY driver/s2_newton.py hessians/final_<R> hessians/final_<R>/y_newton.txt; s4.exe hessian --y ... --out hessians/s2_<R>
```
**Logs:** `logs/execution.log`, `logs/post.log`, `logs/assess.log` and `runs/runs.txt`.

**Result hashes (SHA-256 prefix):**

| File | Hash |
|---|---|
| `results/decisions.json` | 1c0a047e |
| `results/v4_decisions.json` | f8e4be07 |
| `results/runs.csv` | 11477a51 |
| `results/v4_eps.csv` | ef92d174 |
| `assess/assess.csv` | 4283d5b7 |
| `assess/starts.csv` | f6ccda6b |
| `assess/iters.csv` | f399a43c |

---

# Part 2 corrections (Ada, builder, 2026-10-28; applying REVIEW-vera-part2.md F1–F8)

Vera's review cleared Part 2 with corrections (`REVIEW-vera-part2.md`). The corrections below change no result,
classification or licensed claim. The Part 2 text and its result files are left unchanged; the seven Part 2 result
hashes in §13 still match (re-checked 10:09). Where this section and Part 2 differ, this section applies.

- **C1 (F1). B6L2 wording.**
  - §10 D1 now reads: B6L2 "failed (not converged; no study-4 refinement check)" (LINESEARCH_FAIL, n_check 0).
  - §10 D3-B now reads: "failed (not converged)"; B4L2 ITERATION_LIMIT; B6L2 "failed (not converged; no study-4
    refinement check)".
  - This is the registered §4.3 wording. `v4_extra.py` does not produce it.
- **C2 (F2). New deviation 9: the D5w 4:1 λ2 combined label.**
  - PROTOCOL-v4 §8 (verdict across W) gives the literal label **"not established (window-sensitive)"**: W = 20 gives
    "B lower" and W = 25 gives "not evaluable".
  - Part 2 reports "not evaluable (< 4 windows)". That is the reading in the pre-registered `driver/d5w.py`
    (`across_W`), which applies the gap-companion rule to the verdict. `d5w.py` was hashed before any data existed.
  - Both labels give "no level statement", so the licence is unchanged.
  - The W = 20 values "B lower" and "gap grows" in §10 are window-level intermediates. **Do not quote them**, even
    descriptively.
- **C3 (F3). λ_min(H_R)/λ, which §2.3 requires but §9 omitted.**

  | Converged final | λ_min(H_R)/λ |
  |---|---|
  | A6L2m | −0.90 |
  | A4L2f | −0.21 |
  | A6L2f | −0.90 |
  | X4L1fromL2 | −0.18 |

  - H_R is indefinite at every converged final, so H_F is PD only because of the λI term. At A 6:1,
    μ_min(H_F) ≈ 0.10λ.
  - Every "strict local minimiser" statement concerns the regularised F.
- **C4 (F4). Hash statements.**
  - Since 23:40, `sha256sum -c logs/driver.sha256` has given **two** expected FAILED lines: the run-1 hash of
    `fd_checks_v4.py` and the pre-fix hash of `v4_extra.py` (`be2281d9…`).
  - From 07:27 there are **three**: the registered hash of `analyze_runs.py` (`a033a541…`) is now superseded (Part 3,
    deviation P3-1).
  - The pre-fix version of `v4_extra.py` was **not archived**. Deviation 1's statement that the fix changed only the
    output encoding therefore cannot be verified. The failed run wrote nothing and ran before any run result existed.
- **C5 (F5). P4 status.**
  - P4 was started at 07:22:27 under unchanged v4, at the coordinator's request, 72 s after Part 2 was appended. It is
    reported in Part 3.
  - This supersedes "not run" in §8, §10 (λ1 D5w, D6), §12.5 and the related notes.
  - Effect on Part 2's conclusions: D6x and D6 cannot change (see Part 3).
- **C6 (F6). D4 discriminates little.**
  - The one out-of-sample application, A6L2m, started at u_w = 72.9 %, already far above 10 %.
  - X4L1fromL2 (descriptive) is shallow (q 0.72) and is still labelled "fixed-node exploitation" (u_w = 21.7 %). The
    study-3 separation (0.3–4.6 % against 87.9–98.2 %) does not hold at new designs. P4 adds A4L1f: shallow,
    u_w = 21.7 %.
  - A D4 label does not imply a pocket, and D4 has not been validated out-of-sample.
- **C7 (F7). Registered labels attached.**
  - **D4 rows for converged runs.** A6L2m, A4L2f and A6L2f are each "strict local minimiser (second-order)". The
    X4L1fromL2 row is descriptive and carries the same D1b label. Each carries the limitation: "τ_H covers node-sum
    roundoff; a sign change by per-node evaluation error requires a coherent relative error ≥ ε*_all (reported)."
  - **§11 sentences with T labels.** A4L2f, A6L2m and A6L2f: "concentration along the trajectory". X4L1fromL2:
    "concentration along the trajectory, not at the final iterate (converged)". Each D1b use in §10–11 carries the
    NV3 limitation and ε*_all (8.3e-6 to 7.7e-5).
  - **"(not converged)" suffix (general rule, §8).** I adopt this suffix for the T labels of A4L2m, B4L2, B6L2,
    X4L2fromL1 and XB4L2fromA: "concentration along the trajectory (not converged)".
- **C8 (F8). Completeness.**
  - X4L2fromL1's C2' quantity at its last iterate is 1.08e-7.
  - D7 events also occurred at 9 assessed trajectory iterates: B4L2 @100, 117, 125, 150, 175; X4L1fromL2 @125;
    XB4L2fromA @100, 105, 125 (`logs/assess.log`).

---

# Part 3. P4 (λ1 runs) under PROTOCOL-v4 (Ada, builder, 2026-10-28)

These are builder results and have **not been independently reviewed**. The scope is the same as in Part 2.

## P3-1. Start, integrity and procedure
- **Request and checks.** P4 was requested by the coordinator under unchanged v4. The same checks were made before
  the run (`logs/p4_prerun_hashcheck.txt`, 07:22):
  - PROTOCOL-v4 `0d77d7d5…c3ad`, `s4.cpp` `ef247445…9463` and all current driver files are OK;
  - the 2 FAILED lines then were the superseded `fd_checks_v4.py` and `v4_extra.py` hashes;
  - the FD `.out` hashes match;
  - machine load was 4–11 % CPU.
- **Runs.** P1 lane layout (A4 2, B4 2, A6 4, B6 4 = 12 threads), K = 150, λ1. Starts were the four archived λ1
  `y_final` files (B on its reconstructed set) with the step-0 Hessians.
- **Wall time.** P4 ran 07:22:27–08:59 (1.61 h; critical lane A4L1f, 5805 s). Post-run Hessians and assessments
  took 0.6 h, and the S2 Hessians 0.55 h. No run reached its hang guard (caps 24210–36540 s).
- **CRLF fix before any P4 analysis (coordinator; Vera F5).**
  - `driver/analyze_runs.py` `select()` line 82 now writes with `newline="\n"`. This is the only change; see
    `diff --strip-trailing-cr`.
  - The pre-fix copy is archived as `driver/analyze_runs_prefix.py.txt`; it has the registered hash `a033a541…`.
  - The new hash `bcdbf477…` was appended to `logs/driver.sha256` and logged before use.
  - The resulting list `part3/results/assess_list.txt` has 0 CRs.

## P3-2. Runs
| Run | Terminal state | Iterates | ‖∇F‖/tol_C1' | C2' quantity at final | ‖y‖ | Final assessed q | r_E | Checks (B) |
|---|---|---|---|---|---|---|---|---|
| A4L1f | **CONVERGED** | 79 | 0.019 | 6.3e-9 | 16.695 | 0.7228 | 0.575 | – |
| B4L1f | **CONVERGED** (C3: PASS) | 49 | 0.0063 | 3.1e-9 | 16.481 | 0.5768 | 0.575 | 1 (iterate 49, PASS); n_outer 0 |
| A6L1f | **CONVERGED** | 41 | 0.012 | 2.7e-9 | 16.396 | 0.5703 | 0.575 | – |
| B6L1f | **CONVERGED** (C3: PASS) | 49 | 0.0092 | 7.5e-9 | 16.375 | 0.5708 | 0.575 | 1 (iterate 49, PASS); n_outer 0 |

- **Solver counts.** Precision-limited iterates: 0. STALLED: 0. W resets: 0. BFGS skips: 4 (A4L1f), 0 for the others.
  Refreshes: A4L1f 2, the others 1, all triggered by the gradient rule.
- **Descriptive.** A4L1f ends at the same design as X4L1fromL2 (‖Δy‖/‖y‖ = 7.4e-7). The P4 finals otherwise differ
  from each other by ‖Δy‖/‖y‖ = 0.14–0.18.
- **Rate (descriptive, log10‖∇F‖ per iterate over iterates 1–20):** A4L1f +0.185, B4L1f +0.065, A6L1f +0.048,
  B6L1f +0.060.

## P3-3. Second order (D1b), NV3, S2 (`hessians/final_<run>/`, `hessians/s2_<run>/`, `part3/results/v4_eps.csv`)
| Final | H_F label | neg / near-zero / marginal | μ_min | τ_H | μ_min/τ_H | κ(H_F) | λ_min(H_R)/λ | ε*_all | D1b |
|---|---|---|---|---|---|---|---|---|---|
| A4L1f | PD | 0 / 0 / 0 | 3.091e-8 | 3.28e-13 | 9.42e4 | 1.06e4 | −0.18 | 7.7e-5 | strict local minimiser (second-order) |
| B4L1f | PD | 0 / 0 / 0 | 1.354e-8 | 2.32e-13 | 5.84e4 | 1.71e4 | −0.64 | 4.5e-5 | strict local minimiser (second-order) |
| A6L1f | PD | 0 / 0 / 0 | 2.269e-8 | 2.44e-13 | 9.30e4 | 1.08e4 | −0.39 | 7.2e-5 | strict local minimiser (second-order) |
| B6L1f | PD | 0 / 0 / 0 | 1.725e-8 | 2.34e-13 | 7.36e4 | 1.36e4 | −0.54 | 5.7e-5 | strict local minimiser (second-order) |

- **Gates and limitation.** H1a mismatches and pair errors are 0, and the H1b gate passes at all four. Each D1b label
  carries: "τ_H covers node-sum roundoff; a sign change by per-node evaluation error requires a coherent relative
  error ≥ ε*_all (reported)."
- **H_R alone is indefinite at all four** (λ_min(H_R)/λ < 0), so H_F is PD only through λI.
- **No "(marginal)" label arises** (μ_min/τ_H ≥ 5.8e4).
- **NV6 is not applicable**: no final is degenerate or a saddle.
- **S2.**
  - ‖Newton step‖ = 9.0e-7, 2.7e-7, 5.5e-7 and 1.0e-6 (≤ 6.2e-8 relative to ‖y‖).
  - μ_min/τ_H at y + step = 9.42e4, 5.84e4, 9.30e4 and 7.36e4.
  - **"Margin positive"** at all four.

## P3-4. Decision rules (`part3/results/decisions.json`, `part3/results/v4_decisions.json`)
**D1:** all four P4 runs are CONVERGED, and B is converged via PASS. No final is precision-limited.

**D3 at λ1 (P4; study-3 rule):**
- **A: "resolution-dependent".** q = 0.723 (4:1) and 0.570 (6:1) are both < 1.01, but |Δq| = 0.152 ≥ 0.05.
  D1b: strict / strict, with the NV3 limitation.
- **B: "no concentration at convergence".** q = 0.577 and 0.571 (|Δq| = 0.006). D1b: strict / strict, with the NV3
  limitation.

**D6 (λ dependence): "not evaluable".** λ2 D3 is "failed (not converged)" for both formulations (Part 2).

**D6x: "budget-dependent (not converged)", unchanged.**
- The reference A4L1f now exists. It is CONVERGED, strict and shallow (q 0.72), and it is the same point as
  X4L1fromL2.
- Neither arm keeps its start character, so labels 1 and 2 cannot apply.

**D4 (descriptive: P4 is outside the registered P1/P2 scope).** Every row carries "(start on calibration design)", the
post-hoc wording and the C6 caveat. All four finals are D1b strict, with the NV3 limitation. Every assessed worst point
lies in Dw (|x| ≤ 1.06, r ≤ 2.36).

| Run | u_w start → final | Descriptive label | ‖Δy‖/‖y_start‖ | q_final − q_start | R(8:1 Dw)/R_w |
|---|---|---|---|---|---|
| A4L1f | 4.6 % → 21.7 % | fixed-node exploitation (shallow final; see C6) | 0.65 | +0.128 | 1.0015 |
| B4L1f | 1.5 % → 0.45 % | no fixed-node exploitation | 0.63 | −0.024 | 1.0005 |
| A6L1f | 1.6 % → 0.24 % | no fixed-node exploitation | 0.58 | −0.027 | 1.0002 |
| B6L1f | 1.0 % → 0.77 % | no fixed-node exploitation | 0.59 | −0.021 | 1.0000 |

**D5w at λ1 (P4 pair; labels "unequal prior history" and "(determined by converged final(s))"):**
- **4:1: "B lower"** at both W.
  - k_c = 150. A converged at 79 and B at 49, so both are extended as stationary.
  - W = 20: Δ_w = +0.002, +0.028, +0.156, +0.146, +0.146, +0.146, +0.146 (5 of 7 ≥ +0.05; none ≤ −0.05).
  - W = 25: +0.005, +0.071, +0.146, +0.146, +0.146, +0.146 (5 of 6).
  - The padding label applies: n_ext = 5 of m = 7 (W = 20) and 5 of 6 (W = 25). The verdict is therefore a
    comparison of the converged finals (A4L1f q 0.723 against B4L1f q 0.577, both shallow, "none") with the
    trajectories. Its W = 20/25 agreement is not cited as robustness.
  - No refinement label applies, because B is CONVERGED via PASS.
- **6:1: "not established".** Every |Δ_w| ≤ 0.003, with n_ext = 5 at both W.
- **λ-level statement: none** (the two densities differ).
- **Gap companion at λ1: not defined by v4 (deviation P3-2).**
  - The registered baselines (§8 table) exist only for the λ2 P1 pair (s0 = 84/90).
  - `d5w.py` applied those λ2 constants to the λ1 pair, which produced "gap shrinks" (4:1), "not established
    (baseline-sensitive)" (6:1) and a spurious "reconstruction deviation" flag. **These outputs are not used.**
- **Licensed statement (4:1, λ1; §8 table, row "B lower, no refinement label", gap not evaluable).** "B lower at
  matched cumulative iterates (start gap inherited) (determined by converged final(s)) (unequal prior history)":
  persistence of an inherited difference, with no change detected. **"Suppresses" is not used.**
  - Context: the in-loop Δq at iterate 0 is −0.006, so the start gap is near zero.
  - Both finals are shallow (q < 1.01), and the difference is between two shallow strict minimisers (0.723 against
    0.577).
  - The assessed q at the window ends shows no D7-type sign reversal: at iterate 25, A 0.609 against B 0.580.

**T labels (all D1-converged):** A4L1f, B4L1f, A6L1f and B6L1f are all "none along the trajectory".

| Run | Max in-loop q (iterate) | Assessed there | q range, last 20 (flag > 0.05) |
|---|---|---|---|
| A4L1f | 0.857 (51) | 0.857 | 0.716–0.730 |
| B4L1f | 0.699 (18) | 0.699 | 0.577–0.577 |
| A6L1f | 0.690 (13) | 0.690 | 0.570–0.613 |
| B6L1f | 0.986 (9) | 0.986 | 0.570–0.571 |

**X4L1fromL2 iterate 1 (Part 2 deviation 4, now closed):** assessed q = 29.47, equal to the in-loop 29.47.

**D7 (descriptive):**
- Events occurred at A4L1f (final, seed 1) and A4L1f @75.
- Seeds agree on ρ*_D at all four finals.
- The outside-Ab check passed with no fallback.
- Metric validity holds at 21,125/21,125 points at each final.
- The λ1 starts were taken from Part 2 `assess/starts.csv` (unchanged).

## P3-5. What Part 3 licenses (builder reading; for review)
- **At λ1 all four formulation/density runs reach strict local minimisers** of the regularised, discretised F, with
  the NV3 limitation (ε*_all 4.5e-5 to 7.7e-5). H_R is indefinite at each. All are shallow (q 0.57–0.72; "none along
  the trajectory").
- **D3 at λ1:** A "resolution-dependent"; B "no concentration at convergence".
- **D5w at λ1, 4:1:** the statement quoted above, a comparison of two shallow converged finals. 6:1: "not
  established".
- **Not licensed:** any λ-dependence (D6) statement; "suppresses"; any λ1 gap-change statement.

## P3-6. Deviations (declared)
1. **`analyze_runs.py` select() newline fix** (P3-1). It is a code change to a pre-registered script made after the
   Part 2 data, with no effect on selection logic. The pre-fix copy is archived, and the hash was logged before use.
2. **λ1 gap companion.** `d5w.py` uses the λ2-only registered baselines for the λ1 pair. v4 registers no λ1 baselines,
   so the gap companion at λ1 is reported as not defined, and `d5w.py`'s λ1 gap labels and reconstruction flag are not
   used. The D5w verdict itself does not depend on the baselines.
3. **Separate output paths so that the Part 2 files stay byte-identical while under review.**
   - P4 assessments use the same `s4.exe assess` calls and arguments as `driver/assess.sh`, but write to
     `part3/assess/`. The mirror `part3/` also holds copies of the Part 2 run, assessment and summary files.
   - `analyze_runs.py` and `v4_extra.py` were run with root `part3/`.
   - The scripts are `part3/p4_post.sh` (sha256 logged before first run; one non-executed edit fixed and re-logged)
     and `part3/p4_s2.sh`. Both hashes are in `logs/execution.log` but not in `logs/driver.sha256`.
   - Appended, not rewritten: `runs/runs.txt`, `logs/execution.log`, `logs/post_p4.log` and `logs/driver.sha256`.
   - The 7 Part 2 result hashes are unchanged.
4. **D4 at P4** is reported descriptively. The script prints A6L1f's and A4L1f's labels without "[descriptive]";
   P4 is outside the registered P1/P2 scope.
5. **S2 at the B finals** used the run's final set (`--set runs/<R>/set_final.txt`) for the Hessian at y + Newton
   step, as for its final Hessian.

## P3-7. Reproduction (study4; `PY = C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I`)
```
sha256sum -c logs/driver.sha256        # 3 expected FAILED lines (fd_checks_v4 run 1, v4_extra pre-fix, analyze_runs pre-fix)
bash part3/p4_post.sh; bash part3/p4_s2.sh     # after P4 runs via driver/lane.sh (args in runs/<ID>/config.txt)
$PY driver/analyze_runs.py $PWD/part3 report; $PY driver/v4_extra.py $PWD/part3 decide
$PY driver/v4_extra.py $PWD/part3 eps hessians/final_A4L1f hessians/final_B4L1f hessians/final_A6L1f hessians/final_B6L1f
```
**Result hashes (SHA-256 prefix):**

| File | Hash |
|---|---|
| `part3/results/decisions.json` | 23332480 |
| `part3/results/v4_decisions.json` | e270947a |
| `part3/results/runs.csv` | c9382a42 |
| `part3/results/v4_eps.csv` | 5b1e9059 |
| `part3/assess/assess.csv` | 423295a1 |
| `part3/assess/iters.csv` | 8575cdfe |

---

# Part 3 corrections (Ada, builder, 2026-10-28; applying REVIEW-vera-part3.md G1–G6)

Vera's review cleared Part 3 and the Part 2 corrections, with corrections. This section changes no result file and no
earlier section. Where it differs from Part 3 or from "Part 2 corrections", this section applies. The REPORT-v4.md
hash was logged in `logs/execution.log` before and after this append (G7).

- **H1 (G1). λ1 4:1 D5w.**
  - The registered label, verbatim: **"B lower at matched cumulative iterates (start gap inherited) (determined by
    converged final(s)) (unequal prior history)"**.
  - Δq(0) = −0.006, both in-loop and assessed (−0.0063 and −0.0063). At the start, A was slightly lower.
  - The difference is between **two shallow converged finals**: A4L1f q 0.723 against B4L1f q 0.577. Both are "none"
    (q < 1.01).
  - The difference **appeared during study 4**: A went from 0.595 to 0.723, and B from 0.601 to 0.577.
  - Every window with Δ_w ≥ +0.05 contains extended (stationary) values. The two windows that are not extended give
    +0.002 and +0.028 (W = 20).
  - **The P3-4 gloss "persistence of an inherited difference, with no change detected" is withdrawn.** It must not be
    used in the report or in any manuscript. The gap companion is undefined at λ1, so change was not tested.
  - No "suppresses".
- **H2 (G2). λ1 6:1 D5w label and Δq(0).**
  - The full label at 6:1 is **"not established (determined by converged final(s)) (unequal prior history)"**.
    n_ext = 5 of m = 7 at W = 20 and 5 of 6 at W = 25.
  - Δq at iterate 0, in-loop / assessed:

    | Density | In-loop | Assessed |
    |---|---|---|
    | 4:1 | −0.0063 | −0.0063 |
    | 6:1 | +0.0055 | +0.0055 |
- **H3 (G3). Invalid fields in `part3/results/decisions.json`.** P3-7 gains this line: the λ1 fields `D5w.L1_4.licensed`,
  `D5w.L1_4.gap`, `D5w.L1_6.gap` and `D5w.L1_4/L1_6.reconstruction_deviation` (with the λ1 per_W G_hi, G_lo and
  single_baseline entries) were computed with λ2 constants and are **invalid; do not quote them**.
  - A new sidecar file states the same: `part3/results/decisions.README.txt`.
  - The JSON is not edited.
- **H4 (G4). D3 at λ1 is a reading.**
  - v4 §8 registers D3 "at λ2 on the P1 runs". D3 at λ1 is implied only through D6. I applied the same study-3 rule to
    the P4 pair; this is a declared reading, made in `v4_extra.py` (hash `35dc5a3e`, 23:40 on 10-27) before any P4
    data existed.
  - **With labels attached:**
    - A at λ1: "resolution-dependent". Both finals are "none" (0.723 and 0.570, |Δq| = 0.1525 ≥ 0.05). The T labels
      are "none along the trajectory". D1b is strict, with the NV3 limitation.
    - B at λ1: "no concentration at convergence". The T labels are "none along the trajectory". D1b is strict, with
      the NV3 limitation.
  - A's "resolution-dependent" **must not be read as concentration that depends on resolution**. It is a difference
    in q between two finals that both have no concentration.
- **H5 (G5). Final-to-final distances** (‖Δy‖/‖y‖, P4 finals). The P3-2 range "0.14–0.18" is replaced by
  **0.036–0.178**:

  | Pair | ‖Δy‖/‖y‖ |
  |---|---|
  | B4L1f–A6L1f | 0.036 |
  | A4L1f–B6L1f | 0.088 |
  | A4L1f–B4L1f | 0.141 |
  | A4L1f–A6L1f | 0.151 |
  | B4L1f–B6L1f | 0.166 |
  | A6L1f–B6L1f | 0.178 |

  Separately, A4L1f and X4L1fromL2 coincide at 7.4e-7.
- **H6 (G6). Time correction in "Part 2 corrections" C4.** "From 07:27 there are three" should read "from 07:31".
  - The `analyze_runs.py` mtime is 07:30:48.
  - The `bcdbf477` line was appended to `logs/driver.sha256` at 07:31:00.
  - The order was: hash logged, then first use.
- **G7 (process).** No REPORT-v4.md hash was logged at the Part 2 append, and the prefix-compare reference copy from
  10:12 is a session scratch file, not archived in study4. Part 2's byte-identity therefore rests on Vera's indirect
  evidence. From now on the report hash is logged at every append.
