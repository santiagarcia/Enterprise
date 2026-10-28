# PROTOCOL: warp_hypad study 4 (measured Hessian, dense quasi-Newton, nested B, iterate-matched A/B)

Pre-registered by Ada (builder), 2026-10-27, written 20:20–20:51 CDT (hash time in logs/protocol.sha256). **No study-4 computation has been run.**
No objective, gradient, Hessian or assessment has been evaluated for study 4.
- The SHA-256 of this file is recorded in `study4/logs/protocol.sha256` immediately after writing.
- An independent reviewer checks this protocol before any study-4 computation starts. Changes requested by that
  review produce a new protocol version with a new hash; the superseded file and hash are kept.
- After review, nothing here changes. Departures are declared in REPORT.md under "Protocol deviations".

**Inputs read before writing:**
- study3/PROTOCOL.md, REPORT.md, REVIEW-gauss.md (his section 3 is the draft for this study), REVIEW-vera.md;
- study3/driver/s3.cpp (gradient, L-BFGS loop, refine_B, assessment) and the pilot timing logs;
- the archived study-3 traces (`runs/*/trace.csv`, `set_changes.csv`, `termination.txt`, `ys.txt`);
- the OTI headers `Enterprise/warp_hypad/include/warp_hypad/sensitivity/otilib_{scalar,extract}.hpp`.

**Pre-registration artifacts made before this file** (no objective evaluated):
| File | Content | SHA-256 |
|---|---|---|
| `calibration/d5w_calibration.py` | Applies the candidate window rule (section 8, D5w) to the archived study-3 traces | `4c58da899875267cd5f7a3d62625bfe2dcfc81c564272156e08a18c129059d7f` |
| `calibration/d5w_calibration.txt` | Its output | `33c9164844dcdf1754804d12e91dc0b9b7f826aaf20c402470a21d78c5e204ae` |
| `driver/make_hessian_passes.py` | Builds the OTI Hessian pass schedule (section 2.1) | `973ff2aca51fa4594385568ffb597b487bfc76632dc64c0fdb2f68cdc17313aa` |
| `driver/hessian_passes.txt` | The schedule: 1,962 seed sets | `a50821e98765c6ef6efe55d5cba4ff9c943c97cb576bff9e002c9454da6cfc4a` |

**Scope (all statements):**
- the **Eulerian energy density ρ = T(n,n)** (ADM normal observer, TotalEffective source, Λ = 0, t = 0);
- family **F2A**, 108 coefficients, τ = 1e-3;
- the **start-S1 lineage only**: every study-4 start is an archived study-3 iterate that descends from S1
  (`baseline_zero`);
- λ1 = 3.75e-8 and λ2 = 3.75e-9 (study 3, unchanged), working densities 4:1 and 6:1, and this solver.

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
  in-loop search (PROTOCOL-3 §3.2), `refine_B` (study-3 rule with edge-leaf testing), 8:1 diagnostics, assessment.
- b = −9.60726526e-4 and E0 = 4.634569e-2 (study-3 assessment, reproduced there). q = ρ*_D/b, r_E = 1 − E8/E0.

## 2. Step 0 — Hessian first (before any optimisation run)

### 2.1 Method: assembled OTI Hessian (exact up to roundoff), not finite differences
- **Fact used (from the headers):** `OTIScalar` is OTIlib `onumm3n2`: 3 imaginary directions, truncation order 2.
  Every study-3 gradient pass already propagates order-2 coefficients; `extract_hessian(acc, i, j)` reads them via
  `onumm3n2_get_deriv`. The header records a check of this extraction against an analytic Hessian (prior evidence
  only; re-checked in 2.2).
- **Pass schedule.** One pass seeds three coefficients (y_i, y_j, y_k) and returns ∂²R for all six pairs among them,
  plus three gradient components. `driver/hessian_passes.txt` lists 1,962 seed sets (1,908 triples and 54 pairs).
  They come from a cyclic Steiner triple system STS(109) with point 108 removed. **Every unordered pair of the 108
  coefficients occurs in exactly one pass**, and every coefficient occurs in 54 passes. The generator asserts both
  properties.
- **Assembly.** H_R[i][j] = H_R[j][i] comes from the single pass that contains {i, j}. H_R[i][i] is the value from the
  first pass that contains i; the other 53 values are kept for check H1. Then H_F = H_R + λI.
- **Why not FD:** central FD of the OTI gradient needs 216 gradients = 7,776 passes. The OTI assembly needs
  1,962 passes ≈ 54.5 gradient-equivalents, about 4× cheaper, and has no truncation error. FD is used only as the
  accuracy check.
- **Same objective as the run.** The Hessian of R_B is evaluated on B's node set at that design (section 4.3).
  Accumulation over nodes is serial in node order, as for `gradient()`, so it is deterministic.

### 2.2 Accuracy checks (all must pass before any run; results in `hessians/checks.csv`)
- **H1 (repeatability).** At every Hessian design, the 54 values of each H_R[i][i] and of each gradient component
  agree with each other. The gradient components also agree with `gradient()`. Tolerance: 1e-12 relative to
  max|H_R| (or to ‖∇R‖ for the gradient).
- **H2 (column FD).** At `A_4_L1/y_final` and `A_4_L2/y_final`, central differences of the OTI gradient on four
  columns: 0, 36, 72 (the first lapse, axial-shift and radial-shift coefficients) and the column with the largest
  |H_R[i][i]|. Steps h ∈ {1e-2, 1e-3, 1e-4, 1e-5}·max(1, |y_i|). Pass if the best step gives
  ‖FD − H col‖/‖H col‖ ≤ 1e-4, and the error falls with h on the truncation side of the sweep (at least two
  consecutive steps).
- **H3 (directional).** At the same two designs and at `B_4_L2/y_final` on its reconstructed set: one seeded random
  unit vector v (seed 4). Steps as in H2. Compare H_R v with the central FD of ∇R along v, and vᵀH_R v with the
  second difference of R. Same tolerance.
- **H4 (sanity, not a gate).** Check that 1/γ_last from the final study-3 trace row lies in [μ_min(H_F), μ_max(H_F)]
  at that design. γ comes from a step-averaged Hessian, so a failure is reported, not disqualifying.
- If H1–H3 fail, no run starts. A coding error is fixed and logged, and the checks are repeated. If the method
  itself fails (for example, step sensitivity that FD cannot resolve), the study stops for review.

### 2.3 Designs and what is reported
Hessians at 10 designs:
- the **8 archived study-3 finals** (`runs/<F>_<d>_L<k>/y_final.txt`; for B on the reconstructed final set);
- the **2 matched-start A designs** of section 5: `A_4_L2` iterate 84 and `A_6_L2` iterate 90, taken from `ys.txt`
  line k+1. Checked: line 85 and line 91 have ‖y‖ = 38.7803 and 39.8462, equal to the trace values.

The cross-λ starts reuse H_R, because R does not depend on λ: H_F = H_R + λ_new·I.

For each design (`hessians/<design>/H.bin`, `eig.csv`, `summary.csv`):
- all 108 eigenvalues μ_i of H_F and of H_R (`numpy.linalg.eigh` on the symmetric matrix);
- counts of negative, near-zero and positive μ. "Negative" means μ < −τ_H, "near-zero" means |μ| ≤ τ_H,
  with τ_H = max(1e-9·max|μ|, 10 × the H1 discrepancy);
- κ(H_F) = μ_max/μ_min if H_F is PD; otherwise κ_|·| = max|μ|/min|μ| and the label "indefinite";
- λ_min(H_R)/λ;
- the share of each eigenvector in the three field blocks (lapse, axial shift, radial shift);
- the Newton step ‖H_F⁻¹∇F‖ and the decrement ½∇Fᵀ H_F⁻¹∇F / F, if PD (Gauss O4; descriptive distance to
  stationarity);
- ‖∇F‖/tol_C1'.

**Registered labels (descriptive, Question 1):**
- **"PD"** (with the measured κ), **"indefinite (n negative)"** or **"singular to τ_H"**.
- The study-3 statement "apparent condition number 10^3–10^4" is labelled **"consistent"** if the measured κ(H_F) of
  every PD final at that λ lies in [1e3, 1e4]. Otherwise it is **"not supported"**, and the measured range replaces
  it.
- Vera's I5 rewording ("curvature along recent steps exceeds λ by 10^2–10^4") is checked against μ_max/λ in the same
  way.

## 3. Solver (all runs; replaces study-3 L-BFGS m = 8)
**Dense BFGS** on the inverse-Hessian approximation W ≈ H_F⁻¹ (108 × 108; O(n²) per iterate, negligible).
- **Initialisation from the measured Hessian.** Write H_F = QΛQᵀ and set W_0 = Q diag(1/max(|μ_i|, λ)) Qᵀ.
  - This is a builder change to Gauss's "clipped at λ". It uses |μ| (saddle-free scaling): a negative-curvature
    direction gets a curvature-scaled step instead of a 1/λ step.
  - Either choice gives a descent direction. Section 2 reports how many eigenvalues were modified.
- **Direction and line search.** p = −W∇F, truncated to ‖p‖ ≤ trust. The rest is unchanged from study 3 (Gauss:
  "the same Armijo/trust logic"):
  - Armijo c1 = 1e-4, ×0.5 backtracking, ≤ 10 trials;
  - trust 0.5 at start, ×2 up to 4, ×0.5 after a failed line search, LINESEARCH_FAIL below 1e-8;
  - the same validity probe.
- **Update.** The standard BFGS inverse update with s = y_{k+1} − y_k and t = ∇F_{k+1} − ∇F_k, both on the current
  node set.
  - Skip the update if sᵀt ≤ 1e-14‖s‖‖t‖. This is the study-3 rule, now registered (Vera O5). Skips are logged.
  - If gᵀp ≥ 0, reset W to W_0 of the most recent measured Hessian, then log and count the reset.
- **Hessian refresh** (same rule for A and B). At accepted iterate k, re-measure H_F at y_k and reset W to its
  modification (as for W_0) only if all of the following hold:
  - ‖∇F_k‖ > 0.5·‖∇F_{k−20}‖ on the current node set;
  - at least 20 accepted iterates since the last measurement;
  - fewer than 2 refreshes so far in the run.

  Each refresh costs one Hessian (section 10). Refreshes are logged and reported as Hessian evaluations; they are
  not accepted iterates.
- **Logging per accepted iterate:** as study 3 (`trace.csv` columns), plus Newton/BFGS step norm, the BFGS-skip and
  W-reset flags, refresh flag, and inner/outer index for B. Every accepted y goes to `ys.txt`.

## 4. Formulations, stopping and B's nested loop
### 4.1 Stopping criteria (unchanged from study 3)
- **C1':** ‖∇F‖ ≤ 1e-5·F/max(‖y‖, 1).
- **C2':** (F_{k−10} − F_k)/F_k ≤ 1e-8 over iterates on the current node set only.
- **C3 (B only):** the refinement check made no split.
- **STALLED:** C2' ∧ ¬C1' for 30 consecutive iterates on one set.
- Other states: LINESEARCH_FAIL (CONVERGED_LS_LIMITED only if C1' ∧ C2'), INVALID, EVAL_LIMIT (3·K function
  evaluations), and **ITERATION_LIMIT (K, section 6)**.
- WALL_LIMIT exists only as a hang guard (section 6).

### 4.2 A
One BFGS solve on the fixed base set. CONVERGED iff C1' ∧ C2'.

The in-loop search (study-3 §3.2) runs after every accepted iterate. In A it is diagnostic only: it gives in-loop q.

### 4.3 B: refinement only between converged inner solves (Gauss G-I2)
- **Inner solve.** BFGS on a **frozen** node set until inner convergence (C1' ∧ C2' on that set). The in-loop search
  runs after every accepted iterate as a diagnostic, and **no split is made inside an inner solve**.
- **Outer step.** Only after inner convergence, run the study-3 refinement check on the converged design (in-loop
  candidates, 5 % rule, edge-leaf testing, depth ≤ 6, node cap 4× base).
  - No split: C3 holds and the run is **CONVERGED**.
  - Otherwise apply the splits and recompute F and ∇F on the new set. Keep W (as study 3 kept its pairs), restart
    the C2' history, and begin the next inner solve.
- **Outer cap.** 10 refinement rounds. A split at round 10 ends the run as **OUTER_LIMIT** (non-converged).
- **No refinement except at inner convergence.** If an inner solve ends at K, STALLED or LINESEARCH_FAIL, that is
  the run's final state, and no refinement is made at that point.
  - Consequence stated in advance: between refinements, B is a fixed-node optimisation on its current set and could
    exploit that set. The in-loop q trajectory reports this, labelled **"B inner exploitation"**: in-loop q ≥ 1.05
    while the working-node ρ_min is above the in-loop ρ* by more than 10 %.
  - A refinement repairs it only at inner convergence. "B at convergence" means a fixed point of
    (inner optimise → refine).
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
  - The starts already differ, which is inherited from study 3: in-loop q at the start is A 4:1 2.67 vs B 1.05 and
    A 6:1 2.02 vs B 0.87; assessed q at the starts is reported. The P1 window test (D5w) therefore asks whether that difference
    **persists or changes** under continued optimisation at equal budgets. It is not a test from a common point.
- **P4 starts from unequal study-3 histories** (A 245/229 vs B 175/194). At λ1 A and B were indistinguishable after
  iterate ~70 (calibration). D5w at λ1 is reported, labelled "unequal prior history".
- **P5** tests suppression from a common point: does B's refinement remove A's pocket? It is optional and run only
  if budget allows. It is descriptive and not part of D5w.

**Priority and budget rule.**
- Runs are started in priority order. **A started run is never truncated by wall time.**
- A priority that is not started in this session is reported as "not run". It may be run later under this unchanged
  protocol, with the date recorded.
- Within P1, the four runs start concurrently; likewise within P2 + P3.

## 6. Budget, threads and hang guard
- **The budget is in accepted iterates (K above).** A and B in a matched pair always get the same K.
- **Threads:** ≤ 12 in total. P1 uses lanes A4 2, B4 2, A6 4, B6 4; P2 + P3 use A4 2, A6 4, X 2 + 2; P4 uses the P1
  layout. Hessians use 12 threads, one design at a time. Results are thread-count independent (study-3 rerun:
  bit-identical at 2 vs 6 threads).
- **Hang guard:** each run gets a wall cap of 3× its section-10 estimate. Hitting it is WALL_LIMIT and a declared
  deviation. D5w then uses the common iterate window min(k_A, k_B).

## 7. Assessment (unchanged method from study 3, as executed)
- The study-3 assessment code is copied verbatim:
  - region split with the outside-Ab independence check and the floor;
  - h = 0.005 lattice on Ab, 2^17 Halton points (seeded Cranley–Patterson shift);
  - Lipschitz cell screening and D7 sub-lattice (0.00125, ≤ 4096 cells);
  - ≤ 128 starts;
  - Nelder–Mead polish.
- **Declared simplex value.** Study-3 PROTOCOL §7 registered an initial simplex of 0.005. The study-3 code ran
  **0.5h = 0.0025** with two restarts at 0.25h = 0.00125 (s3.cpp lines 866–867; Vera I4). Study 4 registers the
  **as-executed values: 0.0025, restarts 0.00125, stop at diameter 1e-11, ≤ 4000 iterations**. This keeps every
  study-4 q comparable with the study-3 assessed numbers and with b. No other change.
- **Assessed designs (priority):**
  1. every run's final (`--full 1`, two seeds; includes R on Dw at orders 4/6/8, E16, axisymmetry, admissibility);
  2. every start design (two seeds; matched-start A designs and the cross starts);
  3. for each run, the in-loop max-q iterate and the last iterate of every D5w window (W = 25);
  4. every 25th iterate.
- **Use of q.** D5w and T use **in-loop q** at every iterate. In-loop q is an attained value, so it is a lower bound
  on q; in study 3 it agreed with the assessment at 79 of 80 iterates. Assessed q at the window ends is reported
  next to it.
  - If assessed q at a window end would reverse that window's sign of Δmedian by more than 0.05, the event is
    reported as a D7-type event. The registered verdict still uses in-loop q.
- **Second-order check at convergence.** For every CONVERGED run, compute the Hessian (section 2.1, on the run's
  final set) and its spectrum.

## 8. Decision rules
Concentration is defined as in study 3: r_E ≥ 0.10 and q ≥ 1.05; marginal: 1.01 ≤ q < 1.05; none: q < 1.01.

- **H (Question 1; descriptive):** the section-2.3 labels at the 10 designs. They are reported first, regardless of
  how the runs go.
- **D1 Converged:** C1' ∧ C2' (∧ C3 at the outer level for B).
- **D1b Second-order:** for a converged run, "strict local minimiser (second-order)" if every μ(H_F) > τ_H at the
  final; otherwise "degenerate or saddle (n negative, m near-zero)". D3 and D4 statements quote this label.
- **D3 Concentration at convergence** (per formulation, at λ2 on the P1 runs): as in study 3. It needs both densities
  converged; otherwise "failed (not converged)".
  - P2 (A from y_final) is evaluated the same way. If P1 and P2 A converge at one density to finals with
    |Δq| < 0.01 and ‖Δy‖/‖y‖ < 1e-3, they are the same point.
  - Otherwise the label **"start-dependent (multiple stationary points)"** is attached to A's D3 verdict.
- **D4 Fixed-node exploitation at convergence:** as in study 3, on converged A runs (P1, and P2 separately).
- **D5w Matched iterate budget, plateau phase (replaces D5; Gauss G-I3 (b)).** For each λ and density, on the A/B
  pair of P1 (P4 at λ1):
  - **Windows.** Consecutive windows of W accepted iterates, [1, W], [W+1, 2W], …, up to the common horizon
    k_c = min(k_A, k_B). A run that CONVERGED before k_c contributes its converged q to the later windows (it is
    stationary). Δ_w = median in-loop q_A − median in-loop q_B over the window.
  - **Verdict for one W.** "B lower" if Δ_w ≥ +0.05 in ≥ 2/3 of the windows and no window has Δ_w ≤ −0.05.
    "A lower" is the mirror. Otherwise "not established". With fewer than 4 windows: "not evaluable".
  - **Window-alignment robustness.** The verdict must agree at W = 20 and W = 25; otherwise "not established
    (window-sensitive)".
  - **λ level.** "B suppresses deepening at matched iterate budgets" at λ holds only if both densities give
    "B lower". Δq at iterate 0 (inherited) is reported alongside.
  - **Equal ‖y‖ levels are dropped** (Gauss: both formulations sit at the same norm plateau). First-crossing r_E levels
    are dropped (study-3 D5 failure).
  - **Calibration, done before registration** (`calibration/d5w_calibration.txt`, archived study-3 in-loop q, windows
    from A's r_E plateau onset, defined as the first iterate after which r_E stays in a 0.01 band for 10 iterates):
    1. At λ1 every variant gives "not established": late-phase |Δ_w| ≤ 0.011 after iterate ~70, so there was no
       false positive. One early window at 4:1 reached +0.054 (W = 25) and +0.051 (W = 20); it is isolated, so it
       does not decide.
    2. At λ2, with only 2–3 windows available, the verdict **depended on window alignment**. At 6:1, W = 25 from
       the plateau onset gave "B lower", but W = 20 gave "not established", because B > A in the window
       ending at iterate 85 (Δ = −0.090). The reviewers' common-window reading is "indistinguishable at 6:1".
    3. This is why the registered form needs ≥ 4 windows and W = 20/25 agreement. Applied to study 3, it gives
       λ1 "not established" and λ2 "not evaluable (< 4 windows)". With K = 200, study 4 has 10 windows (W = 20) or
       8 windows (W = 25).
    4. Limitation: the calibration shows the rule's behaviour on study-3 trajectories, not its power on study-4
       trajectories.
  - Descriptive companions, reported per window: max q and frac(q ≥ 1.05).
- **D6 λ dependence:** as in study 3, if D3 verdicts exist at both λ.
- **D6x Cross warm start (Gauss section 2(4); P3).** The q values are final assessed q.
  - "**budget-dependent**" if X4L1fromL2 ends with q < 1.05 and X4L2fromL1 reaches q ≥ 1.05 (at its final, or at an
    assessed trajectory iterate);
  - "**path-dependent (multiple minima)**" if X4L1fromL2 is CONVERGED with q ≥ 1.05, or X4L2fromL1 is CONVERGED with
    q < 1.05;
  - otherwise, including when the run that keeps its starting character did not converge, "**inconclusive**".
  - The λ1 runs use the same design-norm sense as study 3 (budget ≈ 17 vs ≈ 38). ‖y‖ at the finals is reported.
- **D7 Search sanity, T trajectory labels and iteration stability:** as in study 3, applied to the study-4
  trajectories. The T label counts study-4 iterates only, and the study-3 history is cited separately. Every
  final-state sentence carries the T label.
- **Rate (Gauss O1; descriptive):** the log‖∇F‖ slope over the first 20 iterates after each warm start.

## 9. Pre-run verification gates (after review; before any run)
| Gate | Check |
|---|---|
| V1 | Enterprise HEAD 477141d with a clean tree. Build `s4.exe` with the study-3 flags (`driver/build.sh`). |
| V2 | At `A_4_L2/y_final` (base set) and `B_4_L2/y_final` (reconstructed set), `s4` reproduces `s3` R, F and ∇F **bit-identically** (the shared code is verbatim). |
| V3 | B set reconstruction acceptance (section 4.3) for all four B finals. |
| V4 | H1–H3 (section 2.2). |
| V5 | BFGS unit checks on synthetic problems (no warp objective): a 108-dimensional SPD quadratic with W_0 = H⁻¹ converges in one step; the inverse update satisfies the secant equation to 1e-12; the skip and reset rules fire on constructed cases. |

A gate failure caused by code is fixed and logged, and the gate is repeated. A failure that needs a protocol change
stops the study pending review.

## 10. Cost estimate (ESTIMATE, not a measurement)
**Basis:**
- Gradient timing from `study3/pilot/timing_*.log`. 4:1: 28.5 s at 2 threads, 5.09 s at 12 threads. 6:1: 64.4 s at 2
  threads, 10.96 s at 12 threads. Thread scaling is therefore measured and near-linear (5.6× and 5.9× for 6× the
  threads).
- Per-iterate wall time from the study-3 traces, which include the in-loop search and 8:1 diagnostics:

  | Run | s per iterate | Threads |
  |---|---|---|
  | A 4:1 | 33 | 2 |
  | B 4:1 | 46–54 | 2 |
  | A 6:1 | 35–36 | 4 |
  | B 6:1 | 42–50 | 4 |

  BFGS bookkeeping is negligible.
- **One OTI Hessian** = 1,962 passes ≈ 54.5 gradient-equivalents. At 12 threads: about 4.6 min at 4:1 and 10 min at
  6:1. B sets with 1.12–1.31× the nodes cost proportionally more.

| Item | Thread-hours | Wall at 12 threads |
|---|---|---|
| Step 0: 10 Hessians + H2/H3 checks | ≈ 17 | ≈ 1.4–1.5 h |
| P1 (4 runs × 200 iterates, upper bound: no early convergence) | ≈ 29 | ≈ 2.4 h |
| P2 + P3 (2 × 200 + 2 × 150) | ≈ 17 | ≈ 1.5 h |
| Assessment P1–P3 (finals full; starts; window ends) | ≈ 6 | ≈ 0.5 h |
| Refreshes (≤ 2 per run) and second-order Hessians at converged finals, P1–P3 | 0–20 | 0–1.7 h |
| **Core P0–P3** | | **≈ 5.8 h upper without refreshes or final Hessians; ≈ 7.5 h worst case** |
| P4 (4 × 150) + its assessment | ≈ 23 | ≈ 1.9 h |
| P5 (optional) | ≈ 4.5 | ≈ 0.4 h |

- If the runs converge superlinearly in tens of iterates (the hoped-for outcome, unverified), P1–P3 run time falls to
  roughly a quarter. Second-order Hessians are then added for each converged final. Core ≈ 3.5–4.5 h.
- **Not included:** writing and verifying `s4.cpp` (V1–V5), about 1.5–2 h of builder time before step 0.
- Gauss's 5–7 h estimate assumed FD Hessians at 2 designs and 6 runs. This plan measures 10 Hessians more cheaply
  (OTI) and adds the matched-start A runs.
- **Uncertainty:** about ×1.5. B's per-iterate cost after the change is unknown, because it now refines rarely.

## 11. Outputs and status labels
- `study4/driver/` (s4.cpp, build and lane scripts, analysis);
- `hessians/<design>/`, `runs/<ID>/` (study-3 file set plus `bfgs.csv` and, for B, `outer.csv`), `assess/`,
  `results/`, `figures/`, `logs/`, `REPORT.md`.
- Every analysis number comes from archived code (Vera I1).
- Labels: **reproduced / new / failed**. Builder results are not independently reviewed.

## 12. Builder changes to Gauss's draft (study3/REVIEW-gauss.md §3), declared
1. **The Hessian is assembled from OTI order-2 passes (1,962 STS passes), not FD.** FD is used only for checks
   H2/H3. Gauss judged OTI "probably comparable or dearer", assuming a costlier order-2 type. The gradient already
   uses the order-2 type, so a pass costs one gradient pass.
2. **Hessians at all 8 finals plus 2 matched-start designs**, not 2 designs.
3. **W_0 uses max(|μ|, λ)** rather than clipping at λ.
4. **B has no inner iterate cap.** Refinement happens strictly only at inner convergence, as the task requires.
   Outer cap 10 (Gauss's value).
5. **K = 200 (λ2 A/B), 150 (cross and λ1)** instead of 300. Rationale: with a measured initial Hessian, a run that
   has not converged by 200 iterates is not near a minimiser. D5w then still has ≥ 8 windows. The refresh rule
   covers model degradation.
6. **Matched-start A runs (P1)** from A's iterate at B's final count, so cumulative iterate budgets match. The A
   y_final continuations are kept (P2).
7. **D5w replaces Gauss's (b)/(c):** ≥ 4 windows, a no-opposite-window condition and W-alignment robustness, added
   after the calibration. Equal-R_ref matching (c) is not registered, because the study-3 traces do not log R_ref,
   so it could not be calibrated before registration.
8. **λ1 continuations (P4)** and the optional formulation cross start (P5) are added at lower priority.
