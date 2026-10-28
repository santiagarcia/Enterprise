# PROTOCOL: warp_hypad study 3 (regularised, matched-progress A/B study)

Pre-registered by Ada (builder), 2026-10-27, written 15:59–16:25 CDT, before any study-3 computation.
Starting point: Gauss's draft (study2/REVIEW-gauss.md section 6), adopted with the builder changes listed in section 10.
Inputs read before writing: study2 PROTOCOL/REPORT/driver, REVIEW-gauss.md and REVIEW-vera.md, and archived study-2
traces (scale of f, ‖g‖ and per-iterate cost only). No study-3 number existed when this was written.
After the pilot (section 6) nothing here changes; departures are declared in REPORT.md under "Protocol deviations".
The SHA-256 of this file is written to `study3/logs/protocol.sha256` before the pilot.

**Scope (all statements):** the Eulerian energy density ρ = T(n,n) (ADM normal observer, TotalEffective source,
Λ = 0, t = 0), family F2A, start S1 (`baseline_zero`) only, this L-BFGS port, these two λ values and two working
densities. No other observer, energy condition, family or time slice is evaluated, and no statement about
"the O2 residual" in general is made.

## 0. Question
For a coercive (regularised) objective, does the continuum worst point of ρ deepen **at convergence** under the
fixed-node formulation A and under adaptive argmin re-insertion B, and does B suppress deepening **at matched
progress** (equal validation E_neg reduction, and separately equal design norm), not just at equal iteration count?

## 1. Fixed problem (unchanged from study 2 unless stated)
- Enterprise clone at 477141d, unmodified; existing `build-oti/` libraries; family F2A, 108 coefficients, τ = 1e-3,
  integrand s(ρ)², s = softplus_τ(−ρ); working domain Dw (ξ ∈ [−1.7, 1.7], s ∈ [0, 2.6]); assessment domain
  D (|x| ≤ 7, r ≤ 7).
- Base working sets 4:1 (6,624 nodes) and 6:1 (14,904 nodes); validation set 8:1 on the O2 extent.
- L-BFGS port of study 2 (memory 8, Armijo 1e-4, ×0.5 backtracking, ≤ 10 trials, trust 0.5 → ×2 up to 4, ×0.5 after a
  failed line search, same validity probe at (0.7, 1.3, 0.8)), except the changes in section 3.

## 2. Regularised objective (the bounded design space)
F(y) = R(y) + (λ/2)‖y‖², with R = R_A or R_B exactly as in study 2.
- λ_k = 0.1^k · f0 / y_ref², f0 = 1.499e-4 (baseline R at 4:1), y_ref = 20: **λ1 = 3.75e-8, λ2 = 3.75e-9.** Fixed;
  not tuned. Only the pilot rule (section 6) may change them, once, by ×10.
- **Why Tikhonov rather than a box:** R ≥ 0 is continuous and finite for every y (the lapse is exp of a field and
  γ = I, so the metric is always valid), hence F is continuous and coercive and has a global minimiser; F stays
  smooth, so the L-BFGS port and its line search are unchanged; a hard box would need a projected method and would
  make the active set part of the result.
- **Effect on interpretation (stated in advance):** if y* is a local minimiser of F, then y* is a local minimiser of R
  over the ball ‖y‖ ≤ ‖y*‖ (for nearby y in that ball, R(y) ≥ R(y*) + (λ/2)(‖y*‖² − ‖y‖²) ≥ R(y*)). Converged
  results are therefore statements about **R-optimal designs within a design-norm budget ‖y*‖ in B-spline coefficient
  space** (a basis-dependent norm), not about the unregularised infimum. The two λ values give two budgets.
  A claim about the unregularised problem is made only if the λ1 and λ2 verdicts agree (D6), and even then it is
  labelled an extrapolation.

## 3. Optimiser and formulation changes (fixes for Gauss B2, I1, I2; Vera)
1. **No wasteful L-BFGS resets (B2).** B keeps its curvature pairs across a working-set change; after a change f and
   ∇F are recomputed on the new set and the next pair uses the new-set gradient. Memory is cleared only if the
   L-BFGS direction is not a descent direction; the first step after such a clear is p = −γ_last·∇F (γ_last = s'y/y'y
   of the last accepted pair) instead of p = −∇F. Iterate 1 keeps the study-2 behaviour for both A and B (matched).
   Because refinement no longer resets, every accepted iterate counts toward K_max; resets are logged and counted.
2. **In-loop search (I1)**, run after every accepted iterate in both A (diagnostic) and B (refinement driver):
   241 × 116 lattice on Dw; all discrete 8-neighbour lattice local minima with ρ < 0 are ranked; they are polished
   (compass to 1e-9, as before) in rank order, in batches of 8, and polished results within 1e-3 of each other are
   merged, until **8 distinct polished minima** exist or 64 polishes have been done. The 8 deepest distinct minima are
   the candidates. Duplicates are logged.
3. **B refinement (I2):** the study-2 rule (split leaf L when ρ* < ρ_L − 0.05|ρ*|, re-test the child, depth ≤ 6,
   node cap 4× base, never coarsen) applied to every candidate, with one change: **every** leaf whose closed
   rectangle contains the point is tested (a minimum on a leaf edge is tested on both sides). Correct description:
   B resolves **polished local minima**, not "every pocket of ρ < 0".
4. **Logging (I3, O1):** per accepted iterate F, R, (λ/2)‖y‖², ‖∇F‖, ‖∇R‖, ‖y‖, step, trust, γ_k, reset flag,
   nfev/ngev, nodes, splits, in-loop ρ* and location **before** refinement, and for B the in-loop ρ* after
   refinement is not recomputed (the design is unchanged; only the nodes change), working-node ρ_min **before and
   after** refinement, 8:1 validation E_neg (whole D and restricted to Dw), negative-volume fraction, the tail
   volume V(ρ < 1.05·b_s2) and the L4 norm of ρ_− on the 8:1 set (b_s2 = −9.60726526e-4, the study-2 baseline),
   8:1 node ρ_min, `refine_around_minimum`, and every accepted y.

## 4. Stopping (Gauss section 6 item 3)
- **C1':** ‖∇F‖ ≤ 1e-5 · F / max(‖y‖, 1) (scale-aware, start-independent).
- **C2':** (F_{k−10} − F_k)/F_k ≤ 1e-8, using only iterates since the last working-set change (stated explicitly
  this time; study-2 deviation D3).
- **C3 (B only):** the refinement check at iterate k made no split.
- **CONVERGED** iff C1' ∧ C2' (∧ C3). This is the only pass state.
- Failure states (archived, assessed, labelled non-converged): ITERATION_LIMIT (K_max = 600), EVAL_LIMIT (3000),
  STALLED (C2' true with C1' false at 30 consecutive iterates on one set), LINESEARCH_FAIL (trust < 1e-8; relabelled
  CONVERGED_LS_LIMITED only if C1' ∧ C2' hold), INVALID (non-finite F or ∇F), and **WALL_LIMIT** (section 5 deadline;
  new, required by the 6 h budget).

## 5. Runs, budget, priority
Runs (start S1 only), priority order: P1 λ1 4:1 A,B; P2 λ1 6:1 A,B; P3 λ2 4:1 A,B; P4 λ2 6:1 A,B.
The optional S3 runs of Gauss's draft are **dropped in advance** (budget; section 10).
- Thread budget ≤ 12 in total (other agents may run Abaqus). Four lanes run concurrently: A-4:1 and B-4:1 with
  2 threads each, A-6:1 and B-6:1 with 4 threads each. Each lane runs its λ1 run, then its λ2 run.
- Per-run wall cap: a λ1 run is stopped (WALL_LIMIT) at lane start + 2 h 15 min; a λ2 run at the global run deadline.
- **Global run deadline: pilot end + 3 h 30 min, and in no case later than 20:50 CDT.** Hard stop for all
  computation (runs, assessment, figures): 21:40 CDT (≈ 6 h after the 15:52 task start). Whatever is unfinished is
  reported as not done.
- **Validation diagnostics cost rule:** if the pilot shows the per-iterate 8:1 diagnostics exceed 25 % of the
  per-iterate wall time, they are computed at every 5th accepted iterate (and at the final) in the run, and the
  intermediate iterates needed by D5/T are evaluated afterwards from `ys.txt` with the same code. Otherwise every
  iterate. Recorded in the REPORT either way.

## 6. Pilot (before the main runs)
20 accepted iterates of B 4:1 S1 at λ1 with 12 threads, plus 3 iterates of A 6:1 for timing. Let Δ1 = ‖y_10‖ − ‖y_0‖
and Δ2 = ‖y_20‖ − ‖y_10‖. Pass if Δ2 < Δ1. If not, λ1 and λ2 are both multiplied by 10 once (a declared deviation)
and the pilot is not repeated. The pilot also records the per-iterate time split (gradient, in-loop search, 8:1
diagnostics) used by the section-5 cost rule. Pilot outputs are archived in `study3/pilot/` and not used in decisions.

## 7. Independent assessment (Gauss I4–I6, D7; Vera I2): structurally different from the in-loop search
Shared with the in-loop search: only the ρ evaluator. Different: grid (0.005 vs 0.028), random design, cell-wise
Lipschitz screening, start selection, and local method (Nelder–Mead instead of compass search).
- **Region split.** The design changes ρ only on its active region (study 2: |x| ≤ 1.97, r ≤ 2.99). The search box
  is Ab = [−2.1, 2.1] × [0, 3.1]. For every assessed design, design-independence outside Ab is verified on a 0.05
  lattice over D \ Ab (ρ(y) = ρ(0) exactly). If it fails anywhere, that design is searched on the whole of D with the
  same method and the event is reported. The design-independent floor outside Ab is computed once, on the baseline,
  with the same method on D \ Ab. ρ*_D = min(ρ*_Ab, floor).
- **Search on Ab.** (1) uniform lattice h = 0.005 (841 × 621); (2) 2^17 Halton points with a seeded Cranley–Patterson
  shift; (3) **cell screening (I5a):** for every lattice cell, est = min corner − 2·L_cell·h/√2 with L_cell the
  largest central-difference |∇ρ| at its corners (a heuristic, not certified); **D7:** every cell whose est is more
  than 10 % below the current best value (up to the 4096 lowest-est cells) is resampled on a 0.00125 sub-lattice;
  (4) starts: all discrete 8-neighbour local minima of the lattice, sub-lattice minima and the 64 best Halton
  points, deduplicated at 0.01, best 128; (5) Nelder–Mead polish (initial simplex 0.005, two restarts from the best
  vertex, stop at simplex diameter 1e-11), confined to Ab.
- **Reported:** ρ*_D and location, the 8 deepest distinct minima, sample minimum, minimum cell estimate, number of
  flagged cells, D7 events (sub-lattice found a value below the pre-polish best), agreement with the in-loop ρ*.
  E_neg on 8:1 (D and Dw), negative-volume fraction, tail volume V(ρ < 1.05·b), L4(D) of ρ_− relative to baseline
  (all on the 8:1 set). Finals and baseline additionally: E_neg 16:1, R on Dw at orders 4, 6, 8 (D4), the
  axisymmetry check (2,000 random 3D points, tolerance 1e-9·max|ρ| + 1e-14), and **admissibility (I4):** min/max ln N,
  max |β| and metric validity on a 0.025 lattice over Ab.
- Finals and baseline use two Halton seeds. Seed agreement is reported but is **not** cited as robustness (Vera I2).
- **Assessment priority** (stop at the hard stop): (1) baseline and all finals; (2) the D5 matched-progress iterates;
  (3) the trajectory max-q iterate of each run (T); (4) snapshots every 50 iterates; (5) every 25.
- Limitations stated in advance: sampling plus local search; not a certified bound; ρ*_D is an attained value, hence
  an upper bound on the infimum; the meridian reduction rests on axisymmetry.

## 8. Decision rules (q, b, r_E as in study 2; b and E0 recomputed by the new assessment)
q = ρ*_D/b, r_E = 1 − E_neg(8:1)/E0. A design shows **concentration** if r_E ≥ 0.10 and q ≥ 1.05; marginal if
1.01 ≤ q < 1.05; none if q < 1.01. In-loop q (an attained value of ρ, so a lower bound on the true q) is used along
trajectories; assessed q at assessed iterates.
- **D1 Converged:** C1' ∧ C2' (∧ C3). Only converged runs enter D2, D3, D4, D6.
- **D2 Regularisation dominance:** report R, ‖y‖ and the penalty at the optimum. If B's converged r_E at λ1 is
  below 0.40, the regulariser dominates and D3/D6 rest on λ2 only.
- **D3 Concentration at convergence** (per formulation and λ): q ≥ 1.05 and r_E ≥ 0.10 at both densities with
  |q(4:1) − q(6:1)| < 0.05 → "concentration at convergence"; q < 1.01 at both densities with |Δq| < 0.05 → "no
  concentration at convergence"; otherwise "resolution-dependent". Needs both densities converged; else
  "failed (not converged)".
- **D4 Fixed-node exploitation at convergence:** A is exploitable if q_A(4:1) − q_A(6:1) ≥ 0.05 and the fixed-node
  check (worst of working-node, 8:1-node and `refine_around_minimum`) under-reports ρ*_D by more than 10 % at A's 4:1
  optimum; or if E8(A) and E16(A) differ by more than 2 % while R_A on the 8:1 Dw set exceeds R_A on the working set
  by more than 2 %.
- **D5 Matched progress (Gauss B1; applies to all runs, converged or not; builder change):** for each λ and density,
  find the first accepted iterates where r_E ≥ 0.40, 0.50 and 0.60, and where ‖y‖ ≥ 10, 20 and 40, in A and in B.
  "B suppresses deepening at matched E_neg progress" holds for a λ only if q_A − q_B ≥ 0.05 (assessed q) at ≥ 2 of
  the r_E levels reached by both, at both densities. "A deepens more at matched norm" is reported the same way for
  the ‖y‖ levels (descriptive, not decisive). The final-iterate comparison is reported but not decisive. If fewer
  than 2 r_E levels are reached by both at a density: "not evaluable at matched progress".
- **D6 λ dependence:** if D3 verdicts differ between λ1 and λ2 → "regularisation-dependent"; no claim about the
  unregularised problem.
- **D7 Search sanity:** as in section 7; any deeper point found replaces ρ*_D and is reported.
- **T Trajectory labels (Vera B1; every run):** over all accepted iterates, the count, first and last iterate of
  concentration (in-loop q ≥ 1.05 with r_E ≥ 0.10), the maximum q and its iterate (confirmed by assessment), and one
  label: "concentration along the trajectory" (≥ 1 iterate), "marginal along the trajectory" (max q in [1.01, 1.05)),
  or "none along the trajectory". For a converged run whose optimum shows no concentration but whose trajectory
  does: "concentration along the trajectory, not at convergence". These labels are reported in the decision table
  next to every final-state verdict, and no final-state sentence omits them.
- **Iteration stability** (reported, as in study 2): range of in-loop q over the last 20 accepted iterates
  (flag > 0.05).

## 9. Outputs and status labels
`study3/driver/` (s3.cpp, scripts, analysis), `pilot/`, `runs/<F>_<d>_L<k>/` (config, trace, ys, linesearch,
refine_events, set_changes, inloop_minima, termination, y_final, log), `assess/`, `figures/` (made with Anaconda
Python), `results/`, `logs/`, `REPORT.md`. Status labels: **reproduced** (a study-1/2 number recovered), **new**,
**failed** (a check or criterion not met). Builder results are not independently reviewed.

## 10. Builder changes to Gauss's draft (declared in advance)
1. K_max 600 and eval cap 3000 kept, but a WALL_LIMIT state and per-run/global deadlines are added: at 12 threads
   Gauss's ≈ 70,000 process-seconds plan (4 processes at 5–6 threads) does not fit 6 h. Non-converged runs are
   reported as such.
2. Optional S3 runs dropped.
3. D5 applies to trajectories of all runs (it is a trajectory statement); D1's "only converged runs" applies to
   D2, D3, D4, D6. D5 decides on assessed q at the matched iterates; in-loop q is reported alongside.
4. Reset iterates count toward K_max because refinement no longer resets the memory (section 3.1).
5. Edge-leaf testing in B (section 3.3).
6. The assessment is restricted to the verified active box plus a once-computed floor (section 7), with a finer
   lattice, Lipschitz screening and Nelder–Mead, to make it both affordable and structurally different.
7. STALLED is redefined for the new C2' (30 consecutive iterates).
8. Trajectory labels T (Vera) and the E_neg-on-Dw, tail-volume and L4 measures (Gauss I6) are added.
