# REVIEW (Gauss, numerical methods): warp_hypad study 2

Reviewer: Gauss, 2026-10-27. Independent of the builder (Ada): I did not write or modify the driver, runs or report.
Same-model separate context; this reduces anchoring but is not proof of independence. Vera is rerunning one pair
separately; I did **not** rerun any optimisation or assessment. Everything below comes from reading
`driver/s2.cpp` or post-processing the archived `runs/*/trace.csv`, `ys.txt`, `refine_events.csv`,
`set_changes.csv` and `assess/snapshots.csv` (each check under 1 min). No dispatch packet was supplied and Council
`assemble` was not run for this review; validated memory retrieval was not used.

Scope: Eulerian rho only, as in the protocol. Classification: **blocking** (must be resolved before the manuscript
claim in REPORT section 10 can be used), **important**, **optional**.

## 0. Summary verdict
- Execution matches the protocol and the code does what PROTOCOL section 3 says, with the caveats below. I found no
  arithmetic bug in the objective, refinement or stopping logic.
- **"A is exploitable" is supported**, but the strongest evidence is internal to A (working f falls while the 8:1
  validation E_neg rises and continuum rho* deepens, with working-node rho_min flat), not the A-vs-B final-iterate gap.
- **The A-vs-B comparison is not matched in optimisation progress** (B1). At matched E_neg reduction, A's q is
  similar to B's in 3 of 5 pairs. The headline "q_A = 2.7-23.5 vs q_B = 0.69-0.98" compares A deep in an
  over-fitting phase with a B that never reached that phase. The candidate manuscript claim must be weakened (section 3).
- Non-convergence is most plausibly structural (drift without a demonstrated finite minimiser; I3), so a larger
  budget alone is unlikely to help. A regularised re-run is proposed in section 6 with pre-registered rules.

## 1. Findings, ranked by severity, likelihood and consequence

### B1 (blocking). A and B are compared at equal iteration count, not at equal progress
Evidence (method in section 7):

| Pair | B final: r_E, q_B, norm(y_B) | A at first iterate with r_E >= r_E(B final): it, q_A, norm(y) | A at first iterate with norm(y) >= norm(y_B): it, q_A |
|---|---|---|---|
| 4:1 S1 | 0.611, 0.687, 27.8 | it 16, **0.645**, 32.5 | it 15, 0.687 |
| 6:1 S1 | 0.643, 0.892, 49.8 | it 37, **1.681**, 56.1 | it 22, 0.971 |
| 4:1 S3 | 0.582, 0.827, 19.0 | it 20, **0.748**, 27.8 | it 18, 1.394 |
| 6:1 S3 | 0.628, 0.753, 36.5 | it 28, **0.759**, 58.8 | it 23, 0.748 |
| 4:1 S2 | 0.580, 0.979, 19.3 | it 14, **1.134**, 22.1 | it 13, 1.095 |

q_A here is from the in-loop Dw search. It agrees with the independent assessment rho*_D to all printed digits at
59 of 60 snapshot iterates, so it is a valid proxy.
- A reaches B's final E_neg reduction within 14-37 iterations, with q = 0.65-1.7. A's deepening happens **after**
  its validation r_E has plateaued (figures/iteration_S1.png): the working objective keeps falling but the true
  integral does not.
- Final norm(y): A 95-139, B 19-50. Median accepted step: A 0.64-1.27, B 0.08-0.16. B travelled 3-5x less far in
  design space.
- Consequence: the data support "A over-fits its nodes late in the run". They do **not** show that B reaches the same
  or a better E_neg without deepening *where A deepens*, because B never reached A's best r_E (4:1: 0.61 vs A max
  0.645) or A's distance in design space. Whether B would deepen if driven further is open. REPORT 1.3 ("attributable
  to the fixed-node formulation, not to the integral objective as such") goes beyond the evidence. Rule 7(ii) was
  satisfied as written, so this is a weakness of the pre-registered rule, not a protocol violation.

### B2 (blocking for the comparison, important otherwise). L-BFGS restarts after every set change waste B's budget
- `mode_run` clears S and Y after each refinement. The next direction is -g with alpha = min(1, trust/norm(g)) = 1,
  so the step equals norm(g), about 2e-6 to 1e-5 (median post-change step 1.7e-6 to 9.8e-6). B has 9-17 such
  near-null iterations per run, and 44-74 of its 150 iterations lie within 5 iterates of a reset. A had none after
  iteration 1.
- Even 6 or more iterates after a reset, B's median step (0.09-0.40) is below A's (0.64-0.83), so the landscape also
  slows B. The two causes are confounded.
- Fix: scale the first post-reset step to the trust radius (or reuse the last gamma = s'y/y'y), or keep curvature
  pairs across a refinement (each set change moves f by at most about 1e-5 relative; see `set_changes.csv`), and do
  not count reset iterations toward K_max.

### I1 (important). The in-loop search tests only 1-5 distinct minima per iterate, and can miss the global one
- `pick_distinct` enforces separation >= 0.1 on the **lattice starts**, not on the polished results. Inside one broad
  basin the 8 starts polish to the same point (e.g. B_4_S1 iterate 0: 8 identical events at (0, 2.0013)). The
  refinement is therefore driven by very few minima.
- At B_6_S1 iterate 75 the in-loop search returned rho* = -7.27e-4 at (0.75, 2.41), while the independent assessment
  found -7.63e-4 at (+-1.194, 1.979): a 4.7 % deeper minimum that the refinement rule never tested at that iterate.
  This is the only miss among 30 B snapshots, but it shows that minima outside the polished set are not controlled.
- Fix: deduplicate after polishing (merge within 1e-3) and keep polishing further lattice minima until 8
  **distinct** polished minima exist, or polish every discrete 8-neighbour lattice local minimum below a threshold.
  This costs little compared with one gradient.

### I2 (important). The B rule controls the depth of polished minima, not the integral error
- The split test (rho* below rho_L minus 5 % of |rho*|) compares one point value with one leaf's node minimum. It does
  not test neighbouring leaves (a minimum on a leaf edge is refined on one side only: `refine_B` takes the first leaf
  containing the point), nor shoulders and secondary pockets that are not among the polished minima.
- Such gaps can be exploited only down to the current global-minimum depth: anything deeper becomes the global
  minimum and is refined at the next accepted iterate. This explains B's transient in-loop excursions (q up to 1.42
  at B_4_S1 iterate 8 and 1.70 at B_4_S2 iterate 20; REPORT quotes only the last-20 window). The rule therefore makes
  deepening **detected and reverted**, not impossible: a step can still be accepted *because of* an unresolved pocket.
  E_neg aliasing in B is small (E8 vs E16 within 0.3 % at all B finals), so the leak is minor at these iterates.
- Depth 6 and the 4x cap never bound: maximum depth 3, at most 1.20x base nodes. Tested depths over all B events:
  4324 at depth 0, 1356 at 1, 405 at 2, 42 at 3. The 5 % threshold is the only active parameter and was not varied;
  one sensitivity run at 1 % would show whether conclusions depend on it.
- Wording: PROTOCOL says every detected pocket of rho < 0 is resolved. With a negative-volume fraction of about
  0.98, rho < 0 almost everywhere; the rule resolves **local minima**, not pockets. Correct the description.

### I3 (important). Non-convergence is structural, not a budget issue
- norm(y) grows throughout every run (A: 95-139 at iterate 150, still +10 to +30 per 50 iterates; B: 19-50). For A
  the displacement over consecutive 20-iterate windows has cosine 0.91-0.96: a steady drift in one direction.
- f still drops by 1.2e-3 to 7.8e-3 (relative) per 20 iterates; norm(g)/norm(g0) is not decreasing (A_4_S1 minimum
  7.2e-3, final 1.3e-2). C2 (1e-6 over 5 iterates) was never met in any run; C1 was met at only 3 and 2 iterates in
  A_4_S3 and A_6_S3.
- Mechanism (hypothesis consistent with the data, not proven): s = softplus_tau(-rho) with tau = 1e-3 comparable to
  |rho| (about 1e-3), so s^2 is positive even for rho >= 0 and vanishes only as rho grows without bound. Wherever the
  design can make rho positive, the integrand keeps rewarding larger coefficients, and the infimum need not be
  attained.
- C1 is relative to norm(g0), which depends on the start (S1 starts at a nearly stationary baseline), so C1 is not
  comparable across starts. C2 at 1e-6 per 5 iterates is about 100x tighter than the observed decrease.
- Conditioning was not measured (no Hessian or gamma logging). The trust radius stayed at its 4.0 cap and A took
  unit L-BFGS steps, so the line search was not the limiter; the drift is.

### I4 (important). Physical admissibility of the final designs is unchecked
- Final coefficient magnitudes reach 24-38 (A) and 6-10 (B); A's log-lapse block norms are 52-76. If the quintic
  B-spline coefficients map roughly to values of ln N, these imply lapses of order exp(30) or exp(-30). The only
  validity test is one `adm_metric_jet` probe at (0.7, 1.3, 0.8). Report max |ln N|, max |shift| and metric validity
  on the assessment lattice for each final design; A's worst points may sit in a degenerate regime.

### I5 (important). Worst-point search: credible as an attained value, weak as an independent check
- rho*_D is an attained value, so it is a valid upper bound on the infimum over D by construction. The open question
  is tightness.
- The assessment adds little beyond the in-loop search. The two agree to all digits at 59/60 snapshots. The two
  Halton seeds agree trivially because the 64 polish starts are dominated by the 0.01 lattice (Halton spacing, about
  sqrt(98/65536) = 0.039, is coarser than the lattice). Subdivision only refines cells whose corner minimum is within
  10 % of the current global minimum, so a pocket narrower than about 0.01, with shallow corner values and away from
  the incumbent, is never sampled. The compass search uses only 8 axis/diagonal directions.
- Cheap strengthening: (a) a per-cell Lipschitz estimate, min corner minus 2*L_cell*h/sqrt(2), with L_cell from
  |grad rho| sampled at the corners, reported as a heuristic (not certified) lower estimate; (b) a 0.0025 lattice on
  A's pocket neighbourhoods and on every B refined leaf.

### I6 (important). Metrics
- q (single worst point over baseline worst point) is the right quantity for "the worst point deepens", but it jumps
  when the argmin switches basin, which is why both formulations fail the iteration-stability rule. Also report the
  location jumps and a continuous tail measure, e.g. the volume where rho < 1.05*b and the L4(D) norm of rho_-
  relative to the baseline. These distinguish a measure-zero spike from a growing region.
- The design-independent floor bounds q from below by 0.52; state this when calling q < 1 an improvement.
- r_E uses E0 over D (|x|, r up to 7) while the design acts on Dw, so the design-independent part caps r_E. Also
  report r_E restricted to Dw. E8 vs E16 agree within 2.3 % (A) and 0.3 % (B), so E8 is adequate.
- Neither formulation minimises E_neg (the integral of rho_-): both minimise the integral of softplus^2, which also
  rewards raising rho above 0. This is O2's design; say so when interpreting r_E.

### O1 (optional)
- REPORT D2: the trace f and norm(g) precede refinement, while node columns follow it, so B's rhomin_w plots mix sets.
- Log gamma_k and norm(y_k) in the trace, and log duplicate polished minima.
- B's iteration instability is driven by refinement-induced transients; log q before and after refinement separately.

## 2. Answers to the assigned questions
1. **B rule.** Depth and cap are irrelevant at these iterates. The 5 % test is reasonable for point depth, but it is
   applied to only 1-5 distinct minima, after acceptance, one leaf at a time. It can still be exploited between
   checks and outside the polished set (one observed miss; transient excursions to q 1.4-1.7). Adequate as a
   detector, not as a guarantee (I1, I2).
2. **Worst-point search.** rho*_D is a valid upper bound on the infimum (it is attained). It is probably tight for
   features at least 0.03 wide; pockets under 0.01 away from the incumbent minimum are untested. The in-loop and
   "independent" searches behave as one method (I5).
3. **Convergence.** Structural drift of norm(y) under an objective without a demonstrated minimiser; C1 is
   start-relative and C2 too tight for the decay rate (I3). Proposal in section 6.
4. **q and r_E.** Appropriate as headline numbers, with the additions in I6.
5. **A-vs-B fairness.** Not matched in progress (B1), and B is handicapped by restarts (B2).

## 3. What the current data do support
- A's working objective is exploitable. In A 4:1 S1, the working f falls monotonically, the 8:1 validation r_E falls
  from 0.65 to 0.56 after about iterate 80, and the working-node rho_min stays near baseline while rho*_D reaches
  23x baseline. The deepening shrinks with working density (4:1 much more than 6:1). This conclusion does not need B.
- Over 150 iterations, B's assessed snapshot q ranges from 0.58 to 1.12 (above 1.05 at B_4_S3@100, B_4_S2@75 and
  B_4_S2@100); its in-loop q reached 1.42-1.70 transiently. Final B q is below 1.01 in 5/5.
- Suggested replacement for the candidate manuscript claim: "With fixed-node quadrature the O2 residual is
  exploitable: continued minimisation lowers the working residual while the validation integral rises and the
  continuum worst point deepens. Adaptive argmin re-insertion suppressed this within 150 iterations, but it also
  progressed less far; whether it prevents deepening at comparable progress or at convergence is not established."

## 4. Cost basis
From runs/runs.txt (4 concurrent processes, 5-6 threads each), wall time per accepted iterate including diagnostics:
A 4:1 about 9-18 s, B 4:1 10-20 s, A 6:1 30-33 s, B 6:1 35-37 s.

## 5. Why not simply a larger budget
The evidence in I3 points to drift without a minimiser. 600 unregularised iterations would most likely produce more
drift and further A deepening, not CONVERGED states. A larger budget is not recommended as the next study.

## 6. PROPOSED study 3 (pre-registration draft by Gauss; not adopted; needs builder review and Santiago's approval)
**Question.** For a well-posed (coercive) objective, does the worst point deepen at convergence under A and under B?

**Changes (only these).**
1. Objective F(y) = R(y) + (lambda/2) norm(y)^2, with R as in study 2 (A or B). Since R is non-negative, F is
   coercive and a minimiser exists. lambda_k = 0.1^k * f0 / y_ref^2 with f0 = 1.499e-4 (baseline R at 4:1),
   y_ref = 20 (about B's final scale), k in {1, 2}: lambda1 = 3.75e-8, lambda2 = 3.75e-9. Fixed before running;
   no tuning.
2. B fixes: keep L-BFGS curvature pairs across refinement; scale the first step after any reset to the trust radius;
   reset iterations do not count toward K_max. The in-loop search keeps 8 **distinct polished** minima (merge
   radius 1e-3). The same search runs as a diagnostic in A.
3. Stopping: C1' norm(grad F) <= 1e-5 * F / max(norm(y), 1) (scale-aware and start-independent); C2' relative
   decrease over 10 iterates <= 1e-8; C3 as before. K_max = 600, evaluation cap 3000. Log gamma_k, norm(y),
   norm(grad F) and norm(grad R).
4. Runs, in priority order: S1 at 4:1 and 6:1, A and B, lambda1 then lambda2 (8 runs). Optional S3 at lambda1
   (4 runs) only if the budget allows.
5. Assessment as in study 2, plus I5(a) cell estimate, I6 tail volume and L4 measure, and I4 field admissibility.

**Cost estimate.** Assuming at most 400 iterates to C1' (unverified): 8 runs x 400 iterates x about 22 s mean is
about 70,000 process-seconds; with 4 concurrent processes, about 5 h wall, plus about 1 h of assessment. Hard stop at
8 h wall. Pilot first: 20 iterates of B 4:1 S1 at lambda1 (about 7 min) to confirm that the growth of norm(y) per 10
iterates is decreasing. If not, multiply lambda by 10 once and record it as a deviation.

**Decision rules (proposed).** q, b and r_E as in study 2.
- D1 Converged: C1' and C2' (and C3 for B). Only converged runs enter D3-D6.
- D2 Regularisation dominance: report R and norm(y) at the optimum. If B's converged r_E at lambda1 is below 0.40,
  the regulariser dominates and conclusions rest on lambda2 only.
- D3 Concentration at convergence (per formulation and lambda): q >= 1.05 and r_E >= 0.10 at both densities with
  |q(4:1) - q(6:1)| < 0.05 gives "concentration at convergence". q < 1.01 at both densities with |delta q| < 0.05
  gives "no concentration at convergence". Anything else is "resolution-dependent".
- D4 Fixed-node exploitation at convergence: A is exploitable if q_A(4:1) - q_A(6:1) >= 0.05 and the fixed-node
  check under-reports by more than 10 % at A's 4:1 optimum; or if E8(A) and E16(A) differ by more than 2 % while
  R_A on the 8:1 set exceeds R_A on the working set by more than 2 %.
- D5 Matched-progress comparison (fixes B1): along trajectories, compare q_A and q_B at the first iterates where
  r_E >= 0.40, 0.50 and 0.60, and at norm(y) = 10, 20 and 40. "B suppresses deepening" holds only if
  q_A - q_B >= 0.05 at 2 or more of the r_E levels reached by both, at both densities. The final-iterate comparison
  is reported but not decisive.
- D6 lambda dependence: if D3 verdicts differ between lambda1 and lambda2, report "regularisation-dependent" and make
  no claim about the unregularised problem.
- D7 Search sanity: if the I5(a) estimate in any cell is more than 10 % below rho*_D, run a 0.0025 local lattice
  there. If it finds a deeper point, it replaces rho*_D and the event is reported.

## 7. Reproduction of my numbers
Python (`C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I`) on `study2/runs/<run>/trace.csv` and `ys.txt`:
- norm(y) = norm of row k of `ys.txt` (row k is iterate k);
- q_in = rhostar_in / (-9.60727e-4); r_E = 1 - Eneg_v8 / 4.634569e-2;
- "A at matched r_E" = first A row with r_E >= the last B r_E; likewise for norm(y);
- snapshot cross-check: column 4 (rho*_D) of `assess/snapshots.csv` (no header) vs `rhostar_in` at the same iterate;
- reset statistics: rows with splits > 0 and the `step` of the following row; window = rows within 5 after a split;
- depth counts: column 6 of `runs/B_*/refine_events.csv`; duplicates: repeated (x, r) within one iterate.

## 8. Limitations of this review
- No reruns. Conclusions rest on archived traces and code reading. I3's mechanism and I4's admissibility concern are
  hypotheses with supporting evidence, not demonstrations.
- The matched-progress comparison uses in-loop q, validated against the assessment only at 25-iterate snapshots.
- The study-3 cost estimate assumes convergence within 400 iterates; this is unverified.
- Same-model reviewer: separate context, not a guarantee of independent judgement.
