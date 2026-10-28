# PROTOCOL v2: warp_hypad study 4 (measured Hessian, dense quasi-Newton, nested B, iterate-matched A/B)

Pre-registered by Ada (builder), 2026-10-27. **Version 2**, written after Noether's pre-computation review
(`study4/REVIEW-noether.md`). **No study-4 computation has been run.** No study-4 objective, gradient, Hessian or
assessment has been evaluated.
- v1 (`PROTOCOL.md`, SHA-256 `a21e0ff65283875bf98b441652e777d028999bc5bc67317921eb76e7758bf64c`,
  `logs/protocol.sha256`) is kept unchanged and is **superseded** by this file.
- The SHA-256 of this file is recorded in `study4/logs/protocol-v2.sha256` immediately after writing.
- v2 has not been re-reviewed. A reviewer (Noether or another qualified member) checks the change table (section 13)
  before any study-4 computation. Changes requested by that review produce v3 with a new hash.
- After that review, nothing here changes. Departures are declared in REPORT.md under "Protocol deviations".

**Inputs read before writing v2:** everything listed in v1, plus `REVIEW-noether.md` (R1–R11, S1–S7) and
`study3/driver/s3.cpp` lines 326–378 (`refine_B`) and 630–690 (main loop: `capped`, C3, STALLED).

**Pre-registration artifacts** (no study-4 quantity evaluated; the v2 additions read archived study-3 files only):
| File | Content | SHA-256 |
|---|---|---|
| `calibration/d5w_calibration.py` | v1 calibration of D5w on archived study-3 traces | `4c58da899875267cd5f7a3d62625bfe2dcfc81c564272156e08a18c129059d7f` |
| `calibration/d5w_calibration.txt` | Its output | `33c9164844dcdf1754804d12e91dc0b9b7f826aaf20c402470a21d78c5e204ae` |
| `calibration/d5w_calibration_v2.py` | **v2:** iterate-1 alignment and gap companion on the study-3 λ1 pairs (Noether S5, S4) | `f508d4e05ef44dde5a66f13eb4d382aa6a3885adb9c3107600fca543ebbb8b23` |
| `calibration/d5w_calibration_v2.txt` | Its output | `dfd80e2e3b8ef413edc9b8924df40c90a6acfa0bae517e98a6bd01f8f087a57e` |
| `calibration/d4_retro_check.py` | **v2:** the replacement D4 statistic on the 8 archived study-3 finals (R5) | `767719f085dda2106a886e2476f75dfdd604b4775908f73a1bc3f97e5840a51e` |
| `calibration/d4_retro_check.txt` | Its output | `bdd8a4c148284d11ac577831ad83d0326a3bcce2d8c27d7cdc9995e78351bb0a` |
| `driver/make_hessian_passes.py` | Builds the OTI Hessian pass schedule (section 2.1) | `973ff2aca51fa4594385568ffb597b487bfc76632dc64c0fdb2f68cdc17313aa` |
| `driver/hessian_passes.txt` | The schedule: 1,962 seed sets | `a50821e98765c6ef6efe55d5cba4ff9c943c97cb576bff9e002c9454da6cfc4a` |

Reproduction (from `warp_hypad/`, any Python 3.11+, `-I`):
`python -I study4/calibration/d5w_calibration_v2.py study3` and `python -I study4/calibration/d4_retro_check.py study3`.

**Scope (all statements):**
- the **Eulerian energy density ρ = T(n,n)** (ADM normal observer, TotalEffective source, Λ = 0, t = 0);
- family **F2A**, 108 coefficients, τ = 1e-3;
- the **start-S1 lineage only**: every study-4 start is an archived study-3 iterate that descends from S1
  (`baseline_zero`);
- λ1 = 3.75e-8 and λ2 = 3.75e-9 (study 3, unchanged), working densities 4:1 and 6:1, and this solver.
- **(v2, R-scope)** Every "minimiser", "stationary point", "second-order" and κ statement concerns the **regularised,
  discretised F on a stated node set**, not the continuum problem. The Euclidean ‖y‖ and the criteria C1'/C2' are
  basis-dependent (study-3 limitation, unchanged).

No other observer, energy condition, family, start or time slice is evaluated, and no statement is made about the
unregularised problem or "the O2 residual" in general.

## 0. Questions
1. **Conditioning (measured, replaces the study-3 estimate).** What is the spectrum of ∇²F at the archived study-3
   designs: definiteness, κ, and λ_min(∇²R) relative to λ?
2. **Convergence.** Does a dense quasi-Newton solver initialised from the measured Hessian reach C1' ∧ C2' (∧ C3)
   from the archived designs within a fixed iterate budget? Is a converged point a second-order local minimiser?
3. **A vs B at convergence and at matched iterate budgets (λ2).** Does B (adaptive re-insertion, refined only between
   converged inner solves) end with a shallower continuum worst point than A? Does A's 4:1 pocket persist?
4. **λ path (Gauss G-I4).** Is "deepening at λ2, none at λ1" a design-norm budget effect or a path effect?

## 1. Fixed problem (unchanged from study 3)
- Enterprise clone at 477141d, unmodified; the existing `build-oti/` libraries; build flags of `study3/driver/build.sh`.
- F(y) = R(y) + (λ/2)‖y‖², R = R_A or R_B as in study 3; integrand s(ρ)², s = softplus_τ(−ρ).
- Working domain Dw (ξ ∈ [−1.7, 1.7], s ∈ [0, 2.6]), assessment domain D (|x| ≤ 7, r ≤ 7).
- Base working sets 4:1 (6,624 nodes) and 6:1 (14,904 nodes); validation set 8:1.
- Shared code paths copied verbatim from `study3/driver/s3.cpp`: node sets, `residual`, `gradient`, `node_rho`,
  in-loop search (PROTOCOL-3 §3.2), 8:1 diagnostics, assessment.
- **(v2, R1) `refine_B` is copied with one declared change:** its split decisions (the `need` test, the depth-6 test,
  the node-cap arithmetic, edge-leaf testing, the guard of 16 passes) are verbatim, but it returns an outcome record
  instead of a split count and has **no persistent `capped` flag** (section 4.3).
- **(v2, R6b)** `gradient()` additionally keeps the per-node contributions and sums them a second time in reverse
  node order. The forward sum is unchanged, so R, F and ∇F remain bit-identical to `s3` (gate V2).
- b = −9.60726526e-4 and E0 = 4.634569e-2 (study-3 assessment, reproduced there). q = ρ*_D/b, r_E = 1 − E8/E0.

## 2. Step 0 — Hessian first (before any optimisation run)

### 2.1 Method: assembled OTI Hessian (exact up to roundoff), not finite differences
- **Fact used (from the headers):** `OTIScalar` is OTIlib `onumm3n2`: 3 imaginary directions, truncation order 2.
  Every study-3 gradient pass already propagates order-2 coefficients; `extract_hessian(acc, i, j)` reads them via
  `onumm3n2_get_deriv`. Noether confirmed from `gradient()` that softplus, square and the weighted sum are all done
  in `OTIScalar`, so the order-2 coefficients are exact second derivatives of R on the node set, and
  H_F = H_R + λI is exact up to roundoff (REVIEW R6).
- **Pass schedule.** One pass seeds three coefficients (y_i, y_j, y_k) and returns ∂²R for all six pairs among them,
  plus three gradient components. `driver/hessian_passes.txt` lists 1,962 seed sets (1,908 triples and 54 pairs) from
  a cyclic Steiner triple system STS(109) with point 108 removed. **Every unordered pair of the 108 coefficients
  occurs in exactly one pass** (independently re-checked by Noether: 5,778 pairs, each once), and every coefficient
  occurs in 54 passes.
- **Assembly.** H_R[i][j] = H_R[j][i] comes from the single pass that contains {i, j}. H_R[i][i] is the value from the
  first pass that contains i. Then H_F = H_R + λI. **Symmetry is imposed by construction and is not evidence of
  correctness** (R6).
- **Forward and reverse accumulation (v2, R6b).** In every pass, the node sum is accumulated in node order (H_fwd,
  ∇R_fwd, R_fwd; this is the value used) and also in reverse node order (H_rev, ∇R_rev, R_rev) from the same
  per-node contributions. No extra passes.
- **Why not FD:** central FD of the OTI gradient needs 216 gradients = 7,776 passes. The OTI assembly needs
  1,962 passes ≈ 54.5 gradient-equivalents, about 4× cheaper, and has no truncation error. FD is used only as an
  accuracy check.
- **Same objective as the run.** The Hessian of R_B is evaluated on B's node set at that design (section 4.3): the
  frozen-set objective of the inner solve.

### 2.2 Accuracy checks (all gates must pass before any run; results in `hessians/checks.csv`)
- **H1a (indexing and determinism; not an error estimate).** At every Hessian design, the 54 values of each H_R[i][i]
  and of each gradient component are **bit-identical**, and the gradient components are bit-identical to
  `gradient()`. Noether's argument (R6) predicts bit-identity, because the order-2 coefficient of e_i² does not
  depend on the other seeded directions. Any difference is a coding error. This check measures neither roundoff
  nor truncation, and it is not used in τ_H.
- **H1b (roundoff estimate; measured at every Hessian design, including refreshes and finals).**
  e_H = ‖H_fwd − H_rev‖₂, e_g = ‖∇F_fwd − ∇F_rev‖, e_F = |F_fwd − F_rev|.
  - **Reference and limitation:** this estimates the summation roundoff of the node sum. It does not capture
    per-node evaluation roundoff, which is identical in both orders. The per-node part is bounded only by H2/H3
    (to their 1e-4 resolution).
  - **Gate:** if 10·e_H > 1e-6·max|μ(H_F)| at any design, eigenvalue classification there is precision-limited and
    the study stops for review.
- **H1c (gradient noise floor; reported, gate-free).** e_g / tol_C1' at every Hessian design. In addition, e_g is
  logged at **every accepted iterate** (from the reverse sum in `gradient()`). If e_g > 0.1·tol_C1' at an iterate,
  C1' is labelled **"precision-limited"** there. A D1 statement on a run whose final row is precision-limited carries
  that label, and this is known before D1 is read.
- **H2 (column FD with step sweep).** At `A_4_L1/y_final` and `A_4_L2/y_final`, central differences of the OTI
  gradient on four columns: 0, 36, 72 (the first lapse, axial-shift and radial-shift coefficients) and the column
  with the largest |H_R[i][i]|. Steps h ∈ {1e-2, 1e-3, 1e-4, 1e-5}·max(1, |y_i|).
  - **Reference:** the OTI column. The full sweep e(h) = ‖FD − H col‖/‖H col‖ is reported. On the truncation side
    e(h) must fall by ≥ 30× per decade (consistent with O(h²)) for at least two consecutive steps, and the best step
    must give e ≤ 1e-4.
  - **FD symmetry residual (descriptive):** for each pair (i, j) of checked columns, |FD_ij − FD_ji| at the best step,
    relative to |H_R[i][j]|. It estimates the FD error itself, against the OTI entry as reference.
- **H3 (directional).** One seeded random unit vector v (seed 4) at `A_4_L1/y_final`, `A_4_L2/y_final`,
  `B_4_L2/y_final` (reconstructed set) and **(v2, S7) `A_6_L2/y_final` and `B_6_L2/y_final` (reconstructed set)**.
  Steps as in H2. Compare H_R v with the central FD of ∇R along v, and vᵀH_R v with the second difference of R. Same
  tolerance and slope rule.
- **H4 (sanity, not a gate).** 1/γ_last from the final study-3 trace row lies in [μ_min(H_F), μ_max(H_F)] at that
  design. γ comes from a step-averaged Hessian, so a failure is reported, not disqualifying.
- **H5 = gate V5c (assembly on a known dense Hessian, v2, R6a).** See section 9. It catches isolated off-diagonal
  mapping errors that H2/H3 cannot.
- If a gate fails, no run starts. A coding error is fixed and logged, and the checks are repeated. If the method
  itself fails (for example, step sensitivity that FD cannot resolve, or the H1b gate), the study stops for review.

### 2.3 Designs and what is reported
Hessians at 10 designs:
- the **8 archived study-3 finals** (`runs/<F>_<d>_L<k>/y_final.txt`; for B on the reconstructed final set);
- the **2 matched-start A designs** of section 5: `A_4_L2` iterate 84 and `A_6_L2` iterate 90, taken from `ys.txt`
  line k+1. Checked: line 85 and line 91 have ‖y‖ = 38.7803 and 39.8462, equal to the trace values.

The cross-λ starts reuse H_R, because R does not depend on λ: H_F = H_R + λ_new·I.

For each design (`hessians/<design>/H.bin`, `eig.csv`, `summary.csv`):
- all 108 eigenvalues μ_i of H_F and of H_R (`numpy.linalg.eigh` on the symmetric matrix);
- e_H, e_g, e_F (H1b) and e_g/tol_C1';
- counts of negative, near-zero and positive μ. "Negative" means μ < −τ_H, "near-zero" means |μ| ≤ τ_H, with
  **τ_H = max(1e-9·max|μ|, 10·e_H)** (v2, R6);
- κ(H_F) = μ_max/μ_min if H_F is PD; otherwise κ_|·| = max|μ|/min|μ| and the label "indefinite";
- λ_min(H_R)/λ and μ_max(H_F)/λ;
- the share of each eigenvector in the three field blocks (lapse, axial shift, radial shift), **including (v2, S3) the
  shares of the negative and near-zero eigenvectors**, reported separately;
- the Newton step ‖H_F⁻¹∇F‖ and the decrement ½∇Fᵀ H_F⁻¹∇F / F, if PD (Gauss O4; descriptive);
- ‖∇F‖/tol_C1' and **(v2) ‖∇F‖/(λ‖y‖)** (gradient relative to the regulariser gradient; Noether's convergence answer).

**Registered labels (descriptive, Question 1):**
- **"PD"** (with the measured κ), **"indefinite (n negative)"** or **"singular to τ_H"**.
- The study-3 statement "apparent condition number 10^3–10^4" is labelled **"consistent"** if the measured κ(H_F) of
  every PD final at that λ lies in [1e3, 1e4], **"not supported"** if at least one PD final lies outside (the measured
  range replaces it), and **(v2, R11) "not evaluable (no PD final)"** if no final at that λ is PD.
- Vera's I5 rewording ("curvature along recent steps exceeds λ by 10^2–10^4") is labelled **"consistent"** if
  μ_max(H_F)/λ at every final at that λ lies in **[1e2, 1e4]** (v2, R11), and "not supported" otherwise, with the
  measured range. (This uses all finals, PD or not, since μ_max > 0 is expected; if μ_max ≤ 0 the label is
  "not evaluable".)

## 3. Solver (all runs; replaces study-3 L-BFGS m = 8)
**Dense BFGS** on the inverse-Hessian approximation W ≈ H_F⁻¹ (108 × 108; O(n²) per iterate, negligible).
- **Initialisation from the measured Hessian.** Write H_F = QΛQᵀ and set W_0 = Q diag(1/max(|μ_i|, λ)) Qᵀ.
  - A builder change to Gauss's "clipped at λ": |μ| (saddle-free scaling) gives a negative-curvature direction a
    curvature-scaled step instead of a 1/λ step. W_0 is SPD either way. Section 2 reports how many eigenvalues were
    modified.
- **Direction and line search.** p = −W∇F, truncated to ‖p‖ ≤ trust. Unchanged from study 3:
  - Armijo c1 = 1e-4, ×0.5 backtracking, ≤ 10 trials;
  - trust 0.5 at start, ×2 up to 4, ×0.5 after a failed line search, LINESEARCH_FAIL below 1e-8;
  - the same validity probe.
  - Noether's assessment (recorded): this is a step-length safeguard, not a model trust region. It is adequate for
    descent but does not guarantee convergence on nonconvex F. A weak-Wolfe search or Powell damping (S1) is
    **not adopted** (section 13); the skip-fraction trigger below is adopted instead.
- **Update.** The standard BFGS inverse update with s = y_{k+1} − y_k and t = ∇F_{k+1} − ∇F_k, both on the current
  node set.
  - Skip the update if sᵀt ≤ 1e-14‖s‖‖t‖ (study-3 rule, registered; Vera O5). Skips are logged, and **(v2, S1) the
    skip fraction over the last 20 accepted iterates is logged at every iterate.**
  - If gᵀp ≥ 0, reset W to W_0 of the most recent measured Hessian, then log and count the reset.
- **Hessian refresh** (same rule for A and B; v2, R10 and S1). At accepted iterate k, re-measure H_F at y_k and reset W
  to its modification (as for W_0) only if all of the following hold:
  1. at least 20 accepted iterates since the last Hessian measurement;
  2. **at least 20 accepted iterates since the last node-set change** (B only; for A this is always true). So
     iterates k−20..k are all on the current set and the comparison in 4 is defined;
  3. fewer than 2 refreshes so far in the run;
  4. ‖∇F_k‖ > 0.5·‖∇F_{k−20}‖, **or** (S1) more than 50 % of the BFGS updates in iterates k−19..k were skipped.

  **A split never triggers a re-measurement by itself**, so the at-most-2 rule stays meaningful. After a split, W is
  kept (as study 3 kept its pairs). Each refresh costs one Hessian (section 10), is logged with its trigger, and is
  not an accepted iterate.
- **Logging per accepted iterate:** as study 3 (`trace.csv` columns), plus Newton/BFGS step norm, the BFGS-skip and
  W-reset flags, skip fraction, refresh flag and trigger, e_g and the precision-limited flag, the STALLED flag, and
  inner/outer index for B. Every accepted y goes to `ys.txt`.

## 4. Formulations, stopping and B's nested loop
### 4.1 Stopping criteria and terminal states
- **C1':** ‖∇F‖ ≤ 1e-5·F/max(‖y‖, 1). (tol_C1' ≈ 2.8e-11 at λ2 and ≈ 6.9e-11 at λ1 at the study-3 finals.)
- **C2':** (F_{k−10} − F_k)/F_k ≤ 1e-8 over iterates on the current node set only.
- **C3 (B only; v2, R1):** at the refinement check made at inner convergence, **no candidate test row has
  `need = 1`** (section 4.3). "No split" is not sufficient: a refused split is not convergence.
- **STALLED (v2, R9): a logged label, not a terminal state.** It is set while C2' ∧ ¬C1' has held for 30 consecutive
  iterates on one set. A runs continue to C1' ∧ C2' or K. B inner solves continue. Rationale (Noether R9): near a
  minimiser, C2' can hold while ‖∇F‖ is several times tol_C1', and at the study-3 linear rate (~0.98) closing that
  gap takes more than 30 iterates; C2' was never true in study 3, so the stop was untested. The number of STALLED
  iterates is reported per run.
- **Terminal states:**

  | State | Applies | Converged? |
  |---|---|---|
  | CONVERGED | A: C1' ∧ C2'. B: inner C1' ∧ C2' and C3 at that check | yes |
  | CONVERGED_LS_LIMITED | LINESEARCH_FAIL with C1' ∧ C2' (B also needs C3 at a check made at that point) | yes (labelled) |
  | ITERATION_LIMIT | K accepted iterates | no |
  | LINESEARCH_FAIL | trust < 1e-8 without C1' ∧ C2' | no |
  | INVALID | validity probe failure that cannot be backtracked | no |
  | EVAL_LIMIT | 3·K function evaluations (B: split re-evaluations count) | no |
  | OUTER_LIMIT | B: a split is needed at outer round 10 | no |
  | **REFINE_BLOCKED_DEPTH** (v2, R1) | B: at a check, ≥ 1 `need = 1` row, no split, all refusals at depth 6 | no |
  | **REFINE_BLOCKED_CAP** (v2, R1) | B: as above with ≥ 1 refusal at the node cap (depth refusals may also occur; both counts are logged) | no |
  | WALL_LIMIT | hang guard (section 6) | no |

### 4.2 A
One BFGS solve on the fixed base set. CONVERGED iff C1' ∧ C2'.

The in-loop search (study-3 §3.2) runs after every accepted iterate. In A it is diagnostic only: it gives in-loop q.

### 4.3 B: refinement only between converged inner solves (Gauss G-I2)
- **Inner solve.** BFGS on a **frozen** node set until inner convergence (C1' ∧ C2' on that set). The in-loop search
  runs after every accepted iterate as a diagnostic, and **no split is made inside an inner solve**.
- **Outer refinement check (v2, R1).** Only after inner convergence, `refine_B` runs on the converged design with the
  study-3 rule (in-loop candidates with v < 0, the 5 % `need` test on every leaf whose closed rectangle contains the
  candidate, depth ≤ 6, node cap 4× base, ≤ 16 passes). It returns an outcome record
  `{n_need, n_split, n_refused_depth, n_refused_cap, n_outside_Dw}` counted over all test rows of the check, and
  writes every row to `refine_events.csv` as in study 3.
  - **PASS:** n_need = 0 → C3 holds → **CONVERGED**.
  - **SPLIT:** n_split ≥ 1 → apply the splits, recompute F and ∇F on the new set, keep W, restart the C2' history,
    and begin the next inner solve. Refusals in the same check are logged. If this is outer round 10, the run ends
    **OUTER_LIMIT** instead.
  - **BLOCKED:** n_need ≥ 1 and n_split = 0 → **REFINE_BLOCKED_CAP** if n_refused_cap ≥ 1, otherwise
    **REFINE_BLOCKED_DEPTH**. The final design is reported as "inner-converged on a blocked set": stationary on its
    frozen set, but not a fixed point of (inner optimise → refine). It is never a D1 convergence.
  - Candidates outside Dw (no containing leaf) are counted in n_outside_Dw and do not affect C3, as in study 3; the
    count is reported.
- **Driver logic plan (v2, R1).** Relative to `s3.cpp`:
  - delete the persistent `bool capped` and the `if (form == 'B' && !capped)` guard. The cap is evaluated per split
    attempt, using the study-3 arithmetic, at every check;
  - `refine_B` is templated on the ρ evaluator, so gate V5d can drive it with a synthetic ρ field;
  - C3 = (n_need == 0) is computed only at a check; between checks B's logged C3 column is empty, not true;
  - `termination.txt` adds n_outer, the last outcome record and the cumulative refusal counts.
- **Outer cap.** 10 refinement rounds (Gauss's value).
- **No refinement except at inner convergence.** If an inner solve ends at K, LINESEARCH_FAIL (without C1' ∧ C2'),
  INVALID or EVAL_LIMIT, that is the run's final state, and no refinement is made at that point. A line-search
  failure with C1' ∧ C2' is an inner convergence, so the check is made (CONVERGED_LS_LIMITED only on PASS).
  - Consequence stated in advance: between refinements, B is a fixed-node optimisation on its current set and could
    exploit that set. The in-loop q trajectory reports this, labelled **"B inner exploitation"**: in-loop q ≥ 1.05
    while the working-node ρ_min is above the in-loop ρ* by more than 10 %.
- **(v2, R7) B may never refine in study 4.** B's inner solve has no iterate cap. If the first inner solve does not
  converge within K = 200 (plausible given study 3), the B run is a fixed-node BFGS run on B's archived refined set
  and ends ITERATION_LIMIT with n_outer = 0. This is accepted as an outcome, not repaired by an inner cap (an inner
  cap would make refinement happen away from inner convergence, which this design excludes). It is reported:
  - n_outer completed, and n_outer completed within the D5w horizon k_c, for every B run;
  - D5w carries the label in section 8 when n_outer within k_c is 0;
  - B's D3/D1 statements read "failed (not converged; no study-4 refinement)".
- **The B warm-start set is the archived final set.** It is reconstructed by replaying the `split` rows of
  `refine_events.csv` in file order on the base set. Acceptance requires all of:
  - the node count equals `termination.txt` (B_4_L1 8,304; B_6_L1 16,632; B_4_L2 8,688; B_6_L2 18,036);
  - the leaf count equals the last trace row;
  - the working-node ρ_min equals the last row's `rhomin_w_post` to 1e-12 relative.

  Fallback: a deterministic `s3.exe` replay with `--kmax k_B`. It is bit-identical per Vera's rerun.
- The first refinement check of a warm-started B run is at its first inner convergence, not at iterate 0.

## 5. Runs (warm starts only; iterate-matched)
Iterate counts are **accepted study-4 iterates after the warm start** (iterate 0 = the start design).

| Pri | ID | Form, density, λ | Start design | K |
|---|---|---|---|---|
| P1 | `A4L2m` | A 4:1 λ2 | `A_4_L2` **iterate 84** (`ys.txt` line 85) | 200 |
| P1 | `B4L2` | B 4:1 λ2 | `B_4_L2/y_final` (iterate 84), reconstructed set | 200 |
| P1 | `A6L2m` | A 6:1 λ2 | `A_6_L2` **iterate 90** (`ys.txt` line 91) | 200 |
| P1 | `B6L2` | B 6:1 λ2 | `B_6_L2/y_final` (iterate 90), reconstructed set | 200 |
| P2 | `A4L2f` | A 4:1 λ2 | `A_4_L2/y_final` (iterate 133) | 200 |
| P2 | `A6L2f` | A 6:1 λ2 | `A_6_L2/y_final` (iterate 126) | 200 |
| P3 | `X4L1fromL2` | A 4:1 **λ1** | `A_4_L2/y_final` (q = 32 pocket, ‖y‖ 37.6) | 150 |
| P3 | `X4L2fromL1` | A 4:1 **λ2** | `A_4_L1/y_final` (q 0.60, ‖y‖ 16.8) | 150 |
| P4 | `A4L1f`, `B4L1f`, `A6L1f`, `B6L1f` | A/B × 4:1/6:1 λ1 | the four λ1 `y_final` (A 245/229, B 175/194) | 150 |
| P5 (optional) | `XB4L2fromA` | B 4:1 λ2 | `A_4_L2/y_final` on the **base** set | 150 |

**Matched budgets.**
- **P1 is the iterate-matched A/B comparison.** Each A run starts from A's own study-3 iterate at the count where B
  stopped (84 at 4:1, 90 at 6:1). A and B therefore enter study 4 with equal cumulative accepted iterates (Gauss
  G-B1, Vera B1) and receive the same K.
  - This is a declared departure from "warm start from y_final" for A. The y_final A continuations are kept as P2.
  - **(v2, R8) Declared confounds** (Noether): (1) the starts differ, inherited from study 3: in-loop q at the start
    is A 4:1 2.67 vs B 1.05 and A 6:1 2.02 vs B 0.87; (2) the cumulative count mixes study-3 L-BFGS (m = 8) iterates
    with study-4 dense-BFGS iterates; (3) B's study-3 history used a different refinement schedule (every iterate).
    P1 is a fair equal-budget comparison of **continuations**; it is not a test of the formulation from a common
    point. D5w asks whether the inherited difference **persists or changes**; the gap companion (section 8) and P5
    address change and common start.
  - The 6:1 A start (iterate 90) precedes A's study-3 onset (~105–110), so the windows cover the onset.
- **P4 starts from unequal study-3 histories** (A 245/229 vs B 175/194). D5w at λ1 is reported, labelled "unequal prior
  history".
- **P5** tests suppression from a common point: does B's refinement remove A's pocket? It is optional and descriptive,
  and is not part of D5w; it may support the word "suppresses" at 4:1 (section 8).

**Priority and budget rule (v2, R11).**
- Runs are started in priority order, except that **P5 may start in the 2 threads left free by P2 + P3** (section 6).
- **No started run is truncated because the session's time budget ends.** The per-run hang guard (section 6) takes
  precedence over this rule: a run that reaches its guard ends WALL_LIMIT and is a declared deviation.
- A priority that is not started in this session is reported as "not run". It may be run later under this unchanged
  protocol, with the date recorded.
- Within P1, the four runs start concurrently; likewise within P2 + P3 (+ P5).

## 6. Budget, threads and hang guard
- **The budget is in accepted iterates (K above).** A and B in a matched pair always get the same K.
- **Threads:** ≤ 12 in total. P1 uses lanes A4 2, B4 2, A6 4, B6 4; P2 + P3 use A4 2, A6 4, X 2 + 2, and P5 (if run)
  2; P4 uses the P1 layout. Step-0 Hessians and post-run Hessians use 12 threads, one design at a time; in-lane
  refreshes use the lane's threads. Results are thread-count independent (study-3 rerun: bit-identical at 2 vs
  6 threads).
- **Hang guard (v2):** each run's wall cap is 3 × (K × its upper s/iterate in section 10 + 2 × its in-lane refresh
  time). Reaching it is WALL_LIMIT and a declared deviation; D5w treats it as a non-converged terminal (section 8).

## 7. Assessment (unchanged method from study 3, as executed)
- The study-3 assessment code is copied verbatim:
  - region split with the outside-Ab independence check and the floor;
  - h = 0.005 lattice on Ab, 2^17 Halton points (seeded Cranley–Patterson shift);
  - Lipschitz cell screening and D7 sub-lattice (0.00125, ≤ 4096 cells);
  - ≤ 128 starts;
  - Nelder–Mead polish.
- **Declared simplex value.** Study-3 PROTOCOL §7 registered an initial simplex of 0.005. The study-3 code ran
  **0.5h = 0.0025** with two restarts at 0.25h = 0.00125 (s3.cpp lines 866–867; Vera I4). Study 4 registers the
  **as-executed values: 0.0025, restarts 0.00125, stop at diameter 1e-11, ≤ 4000 iterations**. No other change.
- **Assessed designs (priority):**
  1. every run's final (`--full 1`, two seeds; includes R on Dw at orders 4/6/8, E16, axisymmetry, admissibility);
  2. every start design (two seeds; matched-start A designs and the cross starts);
  3. for each run, the in-loop max-q iterate and the last iterate of every D5w window (W = 25);
  4. every 25th iterate.
- **Use of q.** D5w and T use **in-loop q** at every iterate. In-loop q is an attained value, so it is a lower bound
  on q; in study 3 it agreed with the assessment at 79 of 80 iterates. Assessed q at the window ends is reported
  next to it.
  - If assessed q at a window end would reverse that window's sign of Δ by more than 0.05, the event is reported as
    a D7-type event. The registered verdict still uses in-loop q.
  - **D3, D4 and D6x use final assessed q** (v2, R3).
- **Second-order check at convergence.** For every CONVERGED run, compute the Hessian (section 2.1, on the run's final
  set), its spectrum and H1b.
- **(v2, S2; descriptive, lowest priority after P1–P3 assessment)** At each converged final labelled D1b "strict",
  measure H_F once more at y + Newton step and report μ_min there and ‖Newton step‖. "Margin positive" if
  μ_min > τ_H at both points. Not run → "not done".

## 8. Decision rules
Concentration is defined as in study 3: r_E ≥ 0.10 and q ≥ 1.05; marginal: 1.01 ≤ q < 1.05; none: q < 1.01. Below,
"pocket" means final assessed q ≥ 1.05 and "shallow" means final assessed q < 1.05.

**General (v2, R3):** any label derived from a non-converged run carries the suffix **"(not converged)"**. "Same point"
means |Δq| < 0.01 and ‖Δy‖/‖y‖ < 1e-3 between two finals.

- **H (Question 1; descriptive):** the section-2.3 labels at the 10 designs. They are reported first, regardless of
  how the runs go.
- **D1 Converged:** terminal state CONVERGED or CONVERGED_LS_LIMITED (section 4.1; for B this includes C3 as defined
  in v2). REFINE_BLOCKED_*, OUTER_LIMIT and every other state are not converged. A precision-limited final (H1c)
  carries that label.
- **D1b Second-order:** for a converged run, "strict local minimiser (second-order)" if every μ(H_F) > τ_H at the
  final; otherwise "degenerate or saddle (n negative, m near-zero)". D3, D4 and D6x statements quote this label.
- **D3 Concentration at convergence** (per formulation, at λ2 on the P1 runs): as in study 3. It needs both densities
  converged; otherwise "failed (not converged)", with the terminal states and, for B, n_outer.
  - **P1-A vs P2-A (v2, R4)**, per density:
    - both CONVERGED and the same point → "same point";
    - both CONVERGED, not the same point, both D1b strict → **"distinct local minimisers"**;
    - both CONVERGED, not the same point, not both strict → **"distinct stationary points"**;
    - otherwise → "P1 not converged", "P2 not converged" or both. No multiple-stationary-point wording is used.
- **D4 Fixed-node exploitation at convergence (v2, R5; replaces the study-3 rule).** Changed because study-3 REPORT
  section 6 showed that the study-3 rule could not flag the A 4:1 λ2 pocket (q 32 vs 4.7): its under-reporting branch
  took the worst of working, 8:1 and polish nodes (≈ 1.5 %), and its E8/E16 branch needed > 2 %.
  - **Statistic:** working-node under-reporting u_w = (ρ_min,w − ρ*_D)/|ρ*_D|, with ρ_min,w the minimum of ρ over
    the run's working node set at the final and ρ*_D the final assessed worst value (the same quantity as section 4.3
    "B inner exploitation", with the assessed instead of the in-loop worst value).
  - **Rule:** on each CONVERGED A run (P1, and P2 separately), "fixed-node exploitation" if u_w > 10 %, otherwise
    "no fixed-node exploitation". If the assessed worst point lies outside Dw, the label is "under-reporting outside
    Dw (domain, not node exploitation)". Computed descriptively for B finals and non-converged finals, with
    "(not converged)".
  - **Companion (descriptive):** R on the 8:1 Dw set / working R; values > 1.02 are noted.
  - **Retrospective check, made before registration on archived data (`calibration/d4_retro_check.txt`):** on the
    8 study-3 finals, u_w is 98.2 % (A 4:1 λ2) and 87.9 % (A 6:1 λ2), and 0.3–4.6 % for the other six. The rule would
    have flagged both A λ2 pockets and nothing else. The threshold (10 %) is the study-3 value and was not tuned;
    the check shows behaviour on non-converged study-3 finals only.
- **D5w Matched iterate budget, plateau phase (replaces D5; Gauss G-I3 (b); v2 R2, R7, R8, S4, S6).** For each λ and
  density, on the A/B pair of P1 (P4 at λ1):
  - **Terminal handling and horizon (v2, R2).**
    - A run that ends CONVERGED (or CONVERGED_LS_LIMITED) at k_f < K is **extended as stationary**: its in-loop
      q(k) = q(k_f) for k_f < k ≤ K.
    - Every other terminal state (ITERATION_LIMIT, LINESEARCH_FAIL, INVALID, EVAL_LIMIT, OUTER_LIMIT,
      REFINE_BLOCKED_*, WALL_LIMIT) **truncates** at its last accepted iterate. STALLED is not terminal (section 4.1).
    - **k_c = min(K, the terminal iterate of each non-converged run in the pair).** If both converged, k_c = K.
    - The number of windows that contain extended (stationary) values is reported.
  - **Windows.** Consecutive **complete** windows [1, W], [W+1, 2W], … with end ≤ k_c; a partial final window is
    dropped. Δ_w = median in-loop q_A − median in-loop q_B over the window.
  - **Verdict for one W.** "B lower" if Δ_w ≥ +0.05 in ≥ 2/3 of the windows and no window has Δ_w ≤ −0.05. "A lower"
    is the mirror. Otherwise "not established". With fewer than 4 windows: "not evaluable (< 4 windows)".
  - **Window-alignment sensitivity (v2, S6).** The verdict must agree at W = 20 and W = 25; otherwise "not established
    (window-sensitive)". The windows overlap heavily, so agreement is reported as **alignment sensitivity, not as
    replication**.
  - **Verdict wording (v2, R8).** Per density: **"B lower at matched cumulative iterates (start gap inherited)"**, or
    the mirror "A lower at matched cumulative iterates (start gap inherited)", "not established" or "not evaluable".
    Δq at iterate 0 (in-loop and assessed) is reported alongside.
  - **No-refinement label (v2, R7).** If B completed 0 outer rounds within k_c, the verdict carries **"(no study-4
    refinement: fixed refined set vs base set)"** and is not cited as a formulation effect.
  - **Gap companion (v2, S4).** G_w = Δ_w − Δ_pre, where Δ_pre = median in-loop q_A − median in-loop q_B over the
    archived study-3 iterates [s0 − W + 1, s0] of the two start trajectories (s0 = 84 at 4:1, 90 at 6:1; both A and
    B traces contain those iterates). Same rule as above: "gap grows", "gap shrinks", "not established",
    "not evaluable", with W = 20/25 agreement.
    - A pre-window baseline is registered rather than the single iterate-0 difference because the v2 calibration
      found that a single-iterate baseline gives 2 false "gap shrinks" verdicts in 14 λ1 alignments, the pre-window
      baseline 0 (`calibration/d5w_calibration_v2.txt`).
  - **"Suppresses" (v2, R8).** The word "B suppresses deepening" at a density is used only if that density's D5w
    verdict is "B lower" without the no-refinement label **and** either the gap companion gives "gap grows", or (4:1
    only) P5 agrees: `XB4L2fromA` final is shallow while `A4L2f` final has a pocket (suffix "(not converged)" if
    either is not converged). Otherwise the R8 wording stands.
  - **λ level.** A λ-level statement is made only if both densities give the same verdict; the wording and labels of
    the weaker density apply.
  - **Equal ‖y‖ levels are dropped** (both formulations sit at the same norm plateau). First-crossing r_E levels are
    dropped (study-3 D5 failure).
  - **Calibration, done before registration (v1, `calibration/d5w_calibration.txt`):** on archived study-3 in-loop q
    with windows from A's r_E plateau onset: at λ1 no false positive; at λ2 (2–3 windows) the verdict depended on
    alignment; hence ≥ 4 windows and W = 20/25 agreement. Applied to study 3: λ1 "not established", λ2 "not
    evaluable".
  - **Calibration v2 (`calibration/d5w_calibration_v2.txt`; Noether S5):** the registered **iterate-1 alignment**
    (windows start at the iterate after the start s0) on the study-3 λ1 pairs, at s0 ∈ {0, 20, 40, 60, 80, 84, 100}
    (4:1) and {0, 20, 40, 60, 80, 90, 100} (6:1), horizon min(k_A, k_B) = 175 / 194:
    - D5w gave 0 false positives in 14 alignments (all "not established" or "not evaluable"); at the registered
      s0 (84 / 90) the late |Δ_w| ≤ 0.015.
    - Gap companion: 2 false "gap shrinks" (4:1, s0 = 40 and 60, where the single-iterate Δ_0 caught a transient of
      +0.14 / +0.08) with a single-iterate baseline, 0 with the registered pre-window baseline.
    - Limitation: this is a null-like case (A ≈ B at λ1 late). It shows the false-positive behaviour of the
      registered alignment, not its power at λ2 on study-4 trajectories.
  - Descriptive companions, reported per window: max q and frac(q ≥ 1.05).
- **D6 λ dependence:** as in study 3, if D3 verdicts exist at both λ.
- **D6x Cross warm start (Gauss section 2(4); P3; v2, R3).** Both arms use **final assessed q** only; in-loop and
  trajectory excursions (study 3: 1.38–1.42 transients on every S1-lineage run) are descriptive only.
  - Arm character: `X4L1fromL2` (λ1, from a λ2 pocket) **keeps** its start character if its final has a pocket and
    **adopts** λ1's character if shallow. `X4L2fromL1` (λ2, from a λ1 shallow final) **keeps** if shallow and
    **adopts** if it has a pocket.
  - Labels, in precedence order (first that applies; mutually exclusive by construction):
    1. **"multiple minima at λ_k"** — an arm keeps its start character, is CONVERGED and D1b strict, and the reference
       continuation at the same λ and density (λ1: `A4L1f` from P4; λ2: `A4L2f` from P2) is CONVERGED, D1b strict,
       of the opposite character, and not the same point. Both minima are named.
    2. **"path-dependent (trajectory)"** — at least one arm keeps its start character, and 1 does not hold.
    3. **"budget-dependent"** — both arms adopt the character of their λ.
    4. **"inconclusive"** — an arm has no assessed final (not run, or INVALID before iterate 1).
  - Labels 2 and 3 carry "(not converged)" if an arm that determines them did not converge. If the reference run
    for label 1 was not run (P4 is low priority), the label is "path-dependent (trajectory)" with "reference not run".
  - The λ1 runs use the same design-norm sense as study 3 (budget ≈ 17 vs ≈ 38). ‖y‖ at the finals is reported.
- **D7 Search sanity, T trajectory labels and iteration stability:** as in study 3, applied to the study-4
  trajectories. The T label counts study-4 iterates only, and the study-3 history is cited separately. Every
  final-state sentence carries the T label.
- **Rate (Gauss O1; descriptive):** the log‖∇F‖ slope over the first 20 iterates after each warm start.

## 9. Pre-run verification gates (after review; before any run)
| Gate | Check |
|---|---|
| V1 | Enterprise HEAD 477141d with a clean tree. Build `s4.exe` with the study-3 flags (`driver/build.sh`). |
| V2 | At `A_4_L2/y_final` (base set) and `B_4_L2/y_final` (reconstructed set), `s4` reproduces `s3` R, F and ∇F **bit-identically** (forward sum unchanged). |
| V3 | B set reconstruction acceptance (section 4.3) for all four B finals. |
| V4 | H1a, H1b gate, H2, H3 (section 2.2). |
| V5a | BFGS unit checks on synthetic problems (no warp objective): a 108-dimensional SPD quadratic with W_0 = H⁻¹ converges in one step; the inverse update satisfies the secant equation to 1e-12; the skip, reset, skip-fraction and refresh-eligibility rules (including the 20-iterates-since-set-change condition) fire on constructed cases. |
| V5b | **Softplus smoothness (R6).** Inspect the Enterprise `softplus` implementation used by `OTIScalar`. If it has a branch (cutoff), record the threshold and evaluate analytically the jump in its first and second derivatives there; pass if each is ≤ 1e-12 relative. Otherwise stop for review. |
| V5c | **Assembly on a known dense Hessian (R6a).** f(y) = ½ yᵀMy + Σ_{m=1..5} softplus_1(w_mᵀy)², with M symmetric with seeded, pairwise-distinct entries and seeded dense w_m, at a seeded random y, pushed through the same 1,962-pass schedule, the same forward/reverse accumulation and the same assembly code. Every one of the 5,886 independent entries must match the analytic Hessian to 1e-12·max|H|, and the gradient likewise. |
| V5d | **Refinement outcome logic (R1).** `refine_B` driven by a synthetic ρ field on the 4:1 base set: a resolved candidate gives PASS; a candidate needing a split at depth 6 gives BLOCKED-DEPTH with C3 false; a set at the node cap gives BLOCKED-CAP with C3 false; a needed split below both limits gives SPLIT; a check after a cap refusal is still executed (no persistent flag). |

A gate failure caused by code is fixed and logged, and the gate is repeated. A failure that needs a protocol change
stops the study pending review.

## 10. Cost estimate (ESTIMATE, not a measurement; v2 corrects the wall times)
**Basis:**
- Gradient timing from `study3/pilot/timing_*.log`. 4:1: 28.5 s at 2 threads, 5.09 s at 12 threads. 6:1: 64.4 s at 2
  threads, 10.96 s at 12 threads. Thread scaling is near-linear (5.6× and 5.9× for 6× the threads).
- Per-iterate wall time from the study-3 traces (including in-loop search and 8:1 diagnostics):

  | Run | s per iterate | Threads |
  |---|---|---|
  | A 4:1 | 33 | 2 |
  | B 4:1 | 43–54 (43 ≈ 33 × 1.31 if B rarely refines) | 2 |
  | A 6:1 | 35–36 | 4 |
  | B 6:1 | 42–50 | 4 |

- **One OTI Hessian** = 1,962 passes ≈ 54.5 gradient-equivalents. At 12 threads: about 4.6 min at 4:1 and 10 min at
  6:1, × the B set ratio (1.12–1.31). **In-lane refresh:** 4:1 at 2 threads 26 min (A) to 34 min (B); 6:1 at
  4 threads 29–31 min (A) to 35–37 min (B).
- **v2 correction (Noether):** v1 divided thread-hours by 12, which assumes perfect packing. Lanes are fixed, so a
  block's wall time is its critical-path lane.

| Item | Critical path | Wall |
|---|---|---|
| Step 0: 10 Hessians (≈ 79 min at 12 threads) + H2/H3 (5 designs) + V5 | — | ≈ 1.4–1.6 h |
| **P1** (4 × 200, no early convergence) | B4L2: 200 × 43–54 s at 2 threads | **≈ 2.4–3.0 h** |
| P1 in-lane refreshes (≤ 2 on the critical lane) | B4L2: 2 × 34 min | + ≤ 1.1 h |
| P2 + P3 (+ P5 in the spare 2 threads) | A6L2f: 200 × 36 s at 4 threads (P5: 150 × 33–54 s ≈ 1.4–2.3 h) | ≈ 2.0 h (≈ 2.3 h with P5) |
| P2 + P3 in-lane refreshes | A6L2f: 2 × 30 min | + ≤ 1.0 h |
| Assessment P1–P3 (finals full; starts; window ends) | 12 threads | ≈ 0.5 h |
| Second-order Hessians at converged finals P1–P3 (≤ 8) | 12 threads | 0–1.0 h |
| S2 Hessians at y + Newton step (≤ 8, lowest priority) | 12 threads | 0–1.0 h |
| **Core P0–P3** | | **≈ 6.3–7.0 h without refreshes or post-run Hessians; ≈ 10 h with all refreshes and second-order Hessians; ≈ 11 h with S2** |
| P4 (4 × 150, P1 layout) + its assessment | B4L1f: 150 × 43–54 s | ≈ 2.0–2.6 h |

- If the runs converge in tens of iterates (hoped for, unverified), P1–P3 run time falls to roughly a quarter and
  the core is ≈ 3.5–5 h including second-order Hessians.
- **Not included:** writing and verifying `s4.cpp` (V1–V5, now including V5b–V5d), about 2–2.5 h of builder time
  before step 0.
- **Uncertainty:** about ×1.5. B's per-iterate cost after the change is unknown, because it now refines rarely.
- Noether's figure (core upper ≈ 6.5–8.5 h) used v1's P2 + P3 value of 1.5 h. With the same packing correction
  applied to P2 + P3 (critical lane A6L2f ≈ 2.0 h), the core is as tabulated above.

## 11. Outputs and status labels
- `study4/driver/` (s4.cpp, build and lane scripts, analysis);
- `hessians/<design>/`, `runs/<ID>/` (study-3 file set plus `bfgs.csv` and, for B, `outer.csv` with every outcome
  record), `assess/`, `results/`, `figures/`, `logs/`, `REPORT.md`.
- Every analysis number comes from archived code (Vera I1).
- Labels: **reproduced / new / failed**. Builder results are not independently reviewed until a reviewer does so.

## 12. Builder changes to Gauss's draft (study3/REVIEW-gauss.md §3), declared (v1, unchanged)
1. **The Hessian is assembled from OTI order-2 passes (1,962 STS passes), not FD.** FD is used only for checks
   H2/H3. The gradient already uses the order-2 type, so a pass costs one gradient pass.
2. **Hessians at all 8 finals plus 2 matched-start designs**, not 2 designs.
3. **W_0 uses max(|μ|, λ)** rather than clipping at λ.
4. **B has no inner iterate cap.** Refinement happens strictly only at inner convergence. Outer cap 10. (v2 states the
   consequence, R7.)
5. **K = 200 (λ2 A/B), 150 (cross and λ1)** instead of 300.
6. **Matched-start A runs (P1)** from A's iterate at B's final count. The A y_final continuations are kept (P2).
7. **D5w replaces Gauss's (b)/(c):** ≥ 4 windows, a no-opposite-window condition and W-alignment sensitivity. Equal-R_ref
   matching (c) is not registered, because the study-3 traces do not log R_ref.
8. **λ1 continuations (P4)** and the optional formulation cross start (P5) are added at lower priority.

## 13. Changes from v1 (Noether review, REVIEW-noether.md)
### Required changes
| ID | Finding | Resolution in v2 | Where |
|---|---|---|---|
| R1 | C3 = "no split" is true when a split is refused at depth 6 or the node cap; `capped` disables the check permanently → false CONVERGED for B | C3 = no `need = 1` row at the check. New non-converged terminals REFINE_BLOCKED_DEPTH / REFINE_BLOCKED_CAP ("inner-converged on a blocked set", never D1). Persistent `capped` flag removed; cap evaluated per attempt at every check. Driver plan: outcome record from `refine_B`; gate V5d tests all four outcomes | §1, §4.1, §4.3, §8 D1, §9 V5d |
| R2 | D5w horizon k_c = min(k_A, k_B) contradicts the converged-run clause | CONVERGED runs extended as stationary to K; k_c = min(K, terminal iterates of non-converged runs); complete windows only; every terminal state's effect stated (non-converged ones truncate); count of extended windows reported | §8 D5w |
| R3 | D6x labels overlap, use different statistics, over-claim "multiple minima" | Final assessed q for both arms; trajectory excursions descriptive; explicit precedence (multiple minima > path-dependent (trajectory) > budget-dependent > inconclusive); "multiple minima" only with two CONVERGED, D1b-strict, distinct finals at the same λ; "(not converged)" suffix | §8 general, D6x |
| R4 | D3 "start-dependent (multiple stationary points)" fires on non-convergence | Labels only when P1-A and P2-A both CONVERGED: "same point", "distinct local minimisers" (both strict), "distinct stationary points"; otherwise "P1/P2 not converged" | §8 D3 |
| R5 | D4 inherited although study 3 showed it misses A's pocket | Replaced by working-node under-reporting u_w > 10 %, change justified by study-3 REPORT §6; outside-Dw case labelled; 8:1/working R ratio kept as descriptive companion; archived-data check shows it flags exactly the two A λ2 pockets | §8 D4, `calibration/d4_retro_check.*` |
| R6 | H1 measures determinism, not error; τ_H rests on ≈ 0; mapping errors can escape | H1 split: H1a (indexing/determinism, bit-identity expected, not in τ_H); H1b forward/reverse node-order accumulation in the same passes with stated reference and limitation, τ_H = max(1e-9 max|μ|, 10 e_H), precision gate; H1c gradient noise floor vs tol_C1' at every iterate with "precision-limited" label; H2 step-sweep slope rule and FD symmetry residual against the OTI reference; V5c synthetic dense-Hessian assembly test (every entry); V5b softplus smoothness | §1, §2.1, §2.2, §2.3, §9 |
| R7 | B may never refine within K = 200 | Stated as an accepted outcome (ITERATION_LIMIT, n_outer = 0; no inner cap added, with reason); n_outer and n_outer within k_c reported; D5w label "(no study-4 refinement: fixed refined set vs base set)", not cited as a formulation effect; D3 wording | §4.3, §8 D5w |
| R8 | D5w wording asserts causation | Verdict renamed "B lower at matched cumulative iterates (start gap inherited)"; "suppresses" only with gap growth (S4) or, at 4:1, P5 agreement; confounds (start gap, L-BFGS/BFGS mix, B's study-3 schedule) declared | §5, §8 D5w |
| R9 | STALLED can stop a still-converging run and truncate D5w | STALLED is a logged label, not terminal, for A and B; count reported | §4.1 |
| R10 | Refresh rule undefined across set changes | Eligible only ≥ 20 iterates after the last set change (and the last measurement); a split never triggers a refresh | §3 |
| R11 | Label completeness | κ label "not evaluable (no PD final)"; I5 interval [1e2, 1e4]; hang guard takes precedence over "never truncated" (which now refers to the session budget) | §2.3, §5, §6 |
| Cost | P1 wall underestimated | P1 ≈ 2.4–3.0 h (critical lane B4L2); packing correction also applied to P2 + P3 (≈ 2.0 h); core ≈ 6.3–7.0 h, ≈ 10–11 h worst case | §10 |
| Scope | — | Statements concern the regularised, discretised F on a stated node set; ‖∇F‖/(λ‖y‖) reported | Scope, §2.3 |

### Suggestions
| ID | Adopted? | Note |
|---|---|---|
| S1 | **Partly.** Skip-fraction logging and the > 50 %-in-20 refresh trigger adopted. Weak-Wolfe / Powell damping not adopted | Keeps the line search identical to study 3, so solver changes are limited to the measured-Hessian initialisation and the dense update. Noether judged the current method well posed as a descent method |
| S2 | **Yes, lowest priority** | Descriptive; "not done" if not reached; costed |
| S3 | **Yes** | Shares of negative and near-zero eigenvectors; no extra cost |
| S4 | **Yes, with a pre-window baseline** | Single-iterate baseline produced false positives in calibration v2 |
| S5 | **Yes, done before registration** | `calibration/d5w_calibration_v2.*`; 0 D5w false positives in 14 λ1 alignments |
| S6 | **Yes** | W = 20/25 agreement reported as alignment sensitivity, not replication |
| S7 | **Yes** | H3 adds `A_6_L2/y_final` and `B_6_L2/y_final` |

### Not changed, and why
- The node-cap and depth limits themselves, the outer cap 10 and K are unchanged; R1 changes only how a refusal is
  classified.
- No inner iterate cap for B (R7 alternative not taken; see section 4.3).
- One new scheduling change not requested by the review: P5 may use the 2 spare threads during P2 + P3, so that the
  S4/R8 "suppresses" route through P5 is not starved by priority order. P5 remains descriptive.
