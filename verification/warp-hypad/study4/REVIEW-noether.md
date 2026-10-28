# REVIEW (Noether): warp_hypad study 4 PROTOCOL, pre-computation

Reviewer: Noether (mathematical foundations). 2026-10-27. **This is a review of a protocol. It contains no study-4 result.**

## Independence and what was examined
- I had no part in drafting PROTOCOL.md (Ada built it from Gauss's study-3 section 3). Separate context, not
  independent verification. Nothing was repaired.
- PROTOCOL.md SHA-256 `a21e0ff6...f758bf64c` matches `logs/protocol.sha256`. The four pre-registration artifact hashes
  in the header table also match.
- Read: study3/REPORT.md, REVIEW-gauss.md (sections 1-3), REVIEW-vera.md (outline), study3/PROTOCOL.md section 8,
  `study3/driver/s3.cpp` (`residual`, `gradient`, `refine_B`, main loop), and `calibration/d5w_calibration.txt`.
- Ran one check, with no objective evaluation: `driver/hessian_passes.txt` has 1,962 sets (1,908 triples, 54 pairs).
  All 5,778 unordered pairs of 108 indices occur exactly once. The schedule is correct.
- Not examined: the OTIlib headers (not at the expected path), `make_hessian_passes.py` source, `d5w_calibration.py`
  source, and the Enterprise `softplus` implementation.

## Overall
Most of the design is sound. The OTI assembly is a good idea: the Hessian is that of the discretised objective on a
frozen node set. The nested B loop is the right fix for G-I2. The matched-start A runs address G-B1.

However, several decision rules are incomplete or overlap, and one inherited code path can make a B run report
CONVERGED when it is not. These must be fixed before running.

## Required changes before running (ranked by severity x likelihood x consequence)

### R1 (high). C3 can hold when refinement was *blocked*, giving a false CONVERGED for B
In `s3.cpp` `refine_B`, a candidate that needs a split but sits at depth 6 is logged `max_depth` and skipped. One
that would exceed the cap is logged `node_cap`, sets `capped = true`, and from then on the loop skips `refine_B`
entirely (`if (form == 'B' && !capped)`). In both cases `splits == 0`, so `C3 = (splits == 0)` is true. PROTOCOL
section 4.3 says "No split: C3 holds and the run is CONVERGED".
- Study 3 never hit the cap (at most 1.31x), but study 4 allows up to 10 outer rounds on sets already 1.12-1.31x.
- Consequence: the "fixed point of (inner optimise -> refine)" can be declared on a set that the refinement rule
  itself rejects.
- **Fix:** define C3 as "no candidate with `need = 1`", not "no split". Add the non-converged terminal states
  DEPTH_LIMIT and CAP_LIMIT, or "converged on a capped set" as a separately labelled, non-D1 state. A `capped` flag
  must never silently turn the check off.

### R2 (high). D5w horizon contradicts its own converged-run clause
k_c = min(k_A, k_B) is the common horizon, but "a run that CONVERGED before k_c contributes its converged q to the
later windows" can never apply. The converged run's k defines k_c, so nothing lies after it.
- If A converges at iterate 30 and B runs to 200, k_c = 30 gives fewer than 4 windows: "not evaluable". That is the
  case where the rule is most informative.
- **Fix:**
  - Extend a CONVERGED run as stationary to K.
  - Let k_c be the minimum over non-CONVERGED terminal iterates, or K.
  - Use complete windows only. A partial final window is dropped.
  - State what STALLED, LINESEARCH_FAIL, INVALID, EVAL_LIMIT and OUTER_LIMIT do to the horizon. With the fix above,
    they truncate it.

### R3 (high). D6x labels are not mutually exclusive, are asymmetric, and over-claim "multiple minima"
- "Budget-dependent" accepts X4L2fromL1 reaching q >= 1.05 "at an assessed trajectory iterate". "Path-dependent"
  fires if X4L2fromL1 is CONVERGED with q < 1.05. A run with a transient excursion that then converges shallow
  satisfies both.
  - Every S1-lineage run in study 3 showed a 1.38-1.42 transient, so this is likely, not hypothetical.
- The two arms use different statistics: the lambda1 arm uses its final, while the lambda2 arm can use any iterate.
- "Path-dependent (multiple minima)" needs two converged minima at the same lambda. One run converging with a pocket
  at lambda1 shows one minimum, not two, unless the shallow lambda1 point (A4L1f, P4) is also CONVERGED. P4 is lower
  priority and may be "not run".
- "Budget-dependent" does not require convergence, yet reads as a statement about the objective.
- **Fix:**
  - Use the final assessed q for both arms; trajectory excursions are descriptive only.
  - Give an explicit precedence order.
  - Write "multiple minima" only when both minima at that lambda are CONVERGED with D1b "strict"; otherwise write
    "path-dependent (trajectory)".
  - Append "(not converged)" to any label from a non-converged arm.

### R4 (medium-high). D3 "start-dependent (multiple stationary points)" fires on non-convergence
"Otherwise" covers P2 not converging, which would attach a multiple-stationary-points label with no second
stationary point.
- **Fix:** apply the label only when P1-A and P2-A are both CONVERGED and differ. If both are D1b-strict, write
  "distinct local minimisers"; otherwise write "distinct stationary points". If P2 did not converge, write
  "P2 not converged".

### R5 (medium-high). D4 is inherited unchanged although study 3 showed it misses the case it exists for
Study-3 REPORT section 6 found that at A 4:1 lambda2 (q 32 vs 4.7, R8/R4 = 1.20) "neither registered branch would
hold as written".
- Branch 1 checks under-reporting by the worst of working, 8:1 and `refine_around_minimum`. The 8:1 and polish
  checks catch the pocket, so that value is about 1.5 %.
- Branch 2 needs E8/E16 to differ by more than 2 %.
- So the rule cannot flag the exploitation already observed. Re-registering a rule known to be blind to that case is
  honest only if declared. It is better replaced.
- **Fix:** use working-node under-reporting, (rho_min,w - rho*_D)/|rho*_D| > 10 %. This is the quantity section 4.3
  already uses for "B inner exploitation". Alternatively use "R on the 8:1 Dw set exceeds working R by > 2 %" alone.
  Either way, state that the change was made because of study-3 REPORT section 6.

### R6 (medium). The Hessian has no accuracy estimate; tau_H and the D1b label rest on a quantity that is probably zero
- **Exactness: yes, modulo roundoff.** In `gradient()`, softplus, square and the weighted sum are all done in
  `OTIScalar`, so the order-2 coefficients of `acc` are exact second derivatives of R on that node set.
  The Hessian of (lambda/2)|y|^2 is lambda*I, so H_F = H_R + lambda*I is exact.
  - For A this is the Hessian of the run's objective.
  - For B it is the Hessian of the frozen-set objective, which is the inner solve's objective. The piecewise
    structure enters only through set changes. Assembling on the current set is therefore correct.
  - Softplus is smooth, provided the Enterprise implementation has no hard cutoff with an O(1) derivative jump. I did
    not check this; one line in V5 should confirm it.
- **H1 does not measure error.** In truncated OTI arithmetic, the e_i^2 coefficient of a product or composition is
  formed from the same operands in the same order whatever the other two seeded directions are. The 54 diagonal values
  will therefore almost certainly be bit-identical. H1 checks indexing and determinism, not roundoff, so
  "10 x H1 discrepancy" in tau_H will be 0.
- **Off-diagonals are each measured once.** An assembly-mapping error (pass slot -> (i, j)) would escape H1 entirely.
  It would be caught by H2 only if it lies in the four checked columns, and by H3 only if it is above 1e-4 of |Hv|.
- **Symmetry is imposed by construction, so it is not evidence.**
- **Fix (cheap):**
  - (a) Add a V5 synthetic test of the assembly code itself: a smooth 108-variable function with a known dense Hessian
    (e.g. a quartic with all pairwise couplings), pushed through the same pass schedule and assembly. It must
    reproduce every entry to roundoff.
  - (b) In the same passes, accumulate the node sum also in reverse node order. This costs nothing extra in passes.
    Report |H_fwd - H_rev|_2, |gradF_fwd - gradF_rev| and |F_fwd - F_rev| at every design.
  - Set tau_H = max(1e-9 max|mu|, 10 |H_fwd - H_rev|_2).
  - Report the gradient noise floor relative to tol_C1'. If it exceeds 0.1 tol, C1' is precision-limited, and that
    must be known before D1 is read.

### R7 (medium). B in P1 may never refine; the D5w verdict would then not test adaptive refinement
B's inner solve has no cap. If the first inner solve does not converge within K = 200, which study 3 makes plausible,
"B" in study 4 is a fixed-node BFGS run on B's archived refined set. D5w would compare A on the base set with a
fixed-node run on a finer, B-shaped set. A "B lower" verdict would then not be evidence about adaptive re-insertion.
- **Fix:** register the reporting of the number of completed outer rounds within k_c. If it is 0, the D5w verdict
  carries the label "no study-4 refinement: fixed refined set vs base set", and it is not cited as a formulation
  effect.

### R8 (medium). D5w verdict wording asserts causation that P1 cannot support
The starts differ: in-loop q is 2.67 vs 1.05 at 4:1 and 2.02 vs 0.87 at 6:1. A "B lower" verdict could therefore be
nothing but persistence of the inherited gap. Section 5 says this, but section 8's verdict text "B suppresses
deepening at matched iterate budgets" does not.
- **Fix:** rename the verdict to "B lower at matched cumulative iterates (start gap inherited)". The word
  "suppresses" is used only if P5 (common start) agrees or the gap grows (see S4).

### R9 (medium). STALLED can terminate a run that is still converging, and it truncates D5w
- Near a minimiser, F - F* ~ |g|^2/(2 mu). Take mu ~ lambda2, F ~ 1.05e-4 and |y| ~ 38, so tol_C1' ~ 2.8e-11.
  Then C2' (a 10-iterate decrease of at most ~1e-12) already holds when |g| is several times tol_C1'.
- Under linear convergence at the study-3 rate (0.98), going from about 5x tol to tol takes more than 30 iterates.
  The run then ends STALLED before C1' and before K.
- C2' was never true in study 3, so this path is untested.
- **Fix:** make STALLED a logged label that does not stop A runs (they go to K or convergence). For B, record it and
  continue the inner solve. Alternatively, justify the stop with a computed estimate. Either way, register it.

### R10 (low-medium). Refresh rule is undefined across B set changes
"|gradF_k| > 0.5 |gradF_{k-20}| on the current node set" is undefined when a split occurred within the last 20
iterates.
- **Fix:** count only iterates since the last set change, so a refresh first becomes eligible 20 iterates after a
  change. State whether a split itself may trigger a re-measurement. I suggest not, so that the at-most-2 rule stays
  meaningful.

### R11 (low). Label completeness
- The "consistent" label for 10^3-10^4 is vacuously true when no final at that lambda is PD. Add "not evaluable
  (no PD final)".
- The mu_max/lambda check of Vera's I5 rewording says "in the same way" but gives no interval. State it:
  [1e2, 1e4].
- "Never truncated by wall time" (section 5) and WALL_LIMIT (section 6) conflict. State that the 3x hang guard takes
  precedence.

## Answers to the specific questions
- **Hessian assembly:** exact for R on the frozen set, and the Tikhonov term is handled correctly (see R6). The FD
  checks H2/H3 are adequate for scale and convention errors, such as a factor of 2 on the diagonal: H3's random
  direction weights every entry. They are not adequate for isolated off-diagonal mapping errors (R6a). Symmetry is
  not a check.
- **Fixed node set:** yes. For B the Hessian is taken on the reconstructed/current set (sections 2.1 and 3 refresh),
  which is the inner solve's smooth objective. F_B is piecewise only across outer steps, and the Hessian is never
  used across one: W is kept as an SPD approximation, which is legitimate.
- **Dense BFGS with indefinite H:** well posed as a descent method. W_0 = Q |Lambda|_lambda^-1 Q^T is SPD, and
  skipping when s^T t <= 0 preserves SPD, so p is always descent and the g^T p >= 0 reset only guards against
  roundoff. The trust cap plus Armijo backtracking is a step-length safeguard, not a model trust region. That is
  adequate for descent, but it does not guarantee convergence on nonconvex F (Armijo without Wolfe; skipping discards
  negative-curvature information). See S1. A full trust-region method is not required for well-posedness.
- **Convergence criteria and lambda:** C1' and C2' are invariant to scaling F. C1' relative to the regulariser
  gradient, tol/(lambda |y|), is about 1.9e-4 at lambda2 and 1.1e-4 at lambda1. This is comparable by coincidence of
  F/penalty, not by design; report |gradF|/(lambda |y|) as well. Neither criterion is invariant to reparametrising
  y, and the Euclidean |y| is basis-dependent (study-3 limitation).
  - Reachability: plausible for lambda1 A.
  - Uncertain for lambda2 A: 12,986x tol at the final, the iterate is still entering a pocket, and the tau = 1e-3
    softplus gives a small Newton region.
  - For B it depends on R1 and R7. Also see R9.
- **B termination and cycling:** splits only add leaves and never coarsen. The node set therefore increases strictly
  through a finite lattice (depth 6, cap 4x), and the outer cap is 10. No cycling is possible. The defect is R1, a
  false C3, not cycling.
- **Matched budgets:** the cumulative accepted iterates are equal, but:
  1. the starts differ (R8);
  2. the cumulative count mixes study-3 L-BFGS m = 8 iterates with study-4 dense-BFGS iterates;
  3. B's study-3 history used a different refinement schedule.

  This is a fair equal-budget comparison of continuations, but it is confounded as a test of the formulation. The
  6:1 A start (iterate 90) precedes A's onset (about 105-110), which is a good property: the windows cover the onset.
- **D5w:**
  - It is fixed before study-4 data, so it cannot be tuned to them.
  - The calibration explored four variants on study-3 data from the same lineage, while knowing that A was deeper at
    lambda2. The changes went in the conservative direction, which limits the harm.
  - The calibration aligned windows to A's plateau onset, while the registered rule aligns to study-4 iterate 1. That
    alignment was never exercised (S5).
  - It does target the divergence phase, because every study-4 iterate is post-divergence.
  - Defects: R2, R7, R8.
- **Decision rules:** not complete or mutually exclusive as written (R1-R5, R11). On non-convergence they are honest
  (D3/D4 "failed (not converged)"), except where noted in R3 and R4.
- **Scope:** consistent: Eulerian rho, F2A, S1 lineage, two lambda values, two densities. Add that every
  "minimiser", "second-order" and kappa statement concerns the regularised, discretised F on a stated node set, not
  the continuum problem.
- **Cost:**
  - The per-Hessian arithmetic checks out: 1,962/36 = 54.5; 54.5 x 5.09 s = 4.6 min; 54.5 x 10.96 s = 10 min. Ten
    Hessians plus H2/H3 take about 1.4 h.
  - **P1 wall time is underestimated.** 2.4 h is thread-hours/12, which assumes perfect packing. Lanes are fixed, so
    the critical path is B4L2 at 2 threads: 200 x 46-54 s, about 2.6-3.0 h. If B rarely refines, use 33 s x 1.31 (set
    size), about 2.4 h.
  - In-lane refreshes add about 26-34 min each at 2 threads for 4:1 and about 30-36 min at 4 threads for 6:1, up to
    two per run.
  - Core upper bound: about 6.5-8.5 h, not 5.8-7.5 h. This does not affect validity.

## Suggestions (not blocking)
- **S1.** Use a weak-Wolfe line search, or Powell-damped BFGS, so that s^T t > 0 is guaranteed and the skip rule
  rarely fires. Log the skip fraction and register a trigger: more than 50 % skips in 20 iterates triggers a refresh,
  if one remains.
- **S2.** At each converged final, report the Newton step against a ball where mu_min(H) stays positive (H measured
  at y + Newton step), as a crude second-order sufficiency margin for D1b.
- **S3.** In the eigenvector field-block shares, also report the share of the negative and near-zero eigenvectors.
  That is where the "weak-influence coefficient" explanation of kappa is tested.
- **S4.** Add a change-based companion to D5w: Delta_w - Delta_0, with the inherited gap at iterate 0 subtracted. A
  growing gap is evidence for "suppression"; a constant gap is persistence.
- **S5.** Re-run `d5w_calibration.py` with iterate-1 alignment on the study-3 lambda1 tails (from iterate 84/90
  onward), to show the registered alignment's false-positive behaviour. This can be done before registration,
  without study-4 data.
- **S6.** Because the windows overlap heavily, W = 20/25 agreement is a weak robustness check. Report it as alignment
  sensitivity, not as replication.
- **S7.** H3: also test one direction at a B 6:1 set and one at a 6:1 A design, because all three H3 designs are
  currently 4:1.

## Reproduction of my checks
```
cd pub-feasibility/warp_hypad
sha256sum study4/PROTOCOL.md study4/driver/* study4/calibration/*; cat study4/logs/protocol.sha256
python -I -c "from itertools import combinations; from collections import Counter; S=[tuple(map(int,l.split(','))) for l in open('study4/driver/hessian_passes.txt') if l.strip() and not l.startswith('#')]; c=Counter(p for s in S for p in combinations(sorted(s),2)); print(len(S),len(c),set(c.values()))"
#   -> 1962 5778 {1}
grep -n "max_depth\|node_cap\|!capped\|C3 =" study3/driver/s3.cpp      # R1
```

## Candidate lessons (for the main session; proposed only)
- When a refinement routine can refuse a split (depth or cap), a "no split" convergence test must count refusals as
  non-convergence (scope: adaptive-quadrature solvers; evidence: s3.cpp `refine_B` and PROTOCOL-4 section 4.3).
- Repeatability across OTI passes does not estimate roundoff, because the order-2 coefficients do not depend on the
  other seeded directions. Use reversed-order accumulation instead (scope: OTI Hessian assembly; evidence: the
  structure of truncated multiplication; argued, not run).

---

## v2 re-check (Noether, 2026-10-27; protocol review only, no study-4 result)

### Independence, integrity and what was run
- I did not draft v2 (Ada did). I did write R1-R11 and S1-S7, so this re-check checks my own requested fixes. I did not
  repair anything.
- `PROTOCOL-v2.md` SHA-256 is `b2c0af84...176885da` and matches `logs/protocol-v2.sha256`. `PROTOCOL.md` is unchanged:
  `a21e0ff6...f758bf64c` matches `logs/protocol.sha256`. All eight artifact hashes in the v2 header table match.
- Both calibrations reproduce byte-for-byte, ignoring CR line endings. They read only `study3/runs/*/trace.csv` and
  `study3/assess/assess.csv`.
- I computed one extra quantity from archived study-3 traces: the registered λ2 pre-window gap baseline (finding N1).
  No objective, gradient or Hessian was evaluated.

### Verdicts on R1-R11
| ID | Verdict | Reason |
|---|---|---|
| R1 | **resolved** | C3 = (n_need == 0) at a check. BLOCKED is never D1. The persistent `capped` flag is gone. V5d tests all four outcomes and that a check still runs after a cap refusal. The case analysis covers every check: PASS when n_need = 0; SPLIT when n_split ≥ 1 (OUTER_LIMIT at round 10); BLOCKED when n_need ≥ 1 and n_split = 0. Two low clarifications are listed under N5. |
| R2 | **resolved** (fairness: see below; one labelling addition, N2) | The horizon and terminal handling are now consistent. Complete windows only. |
| R3 | **resolved** | Final assessed q is used for both arms. The four labels are exhaustive and exclusive under the stated precedence. "Multiple minima" needs two CONVERGED, strict, distinct finals at the same λ. Missing reference → "path-dependent (trajectory)". |
| R4 | **resolved** | — |
| R5 | **resolved as a rule**; reporting change required (N4) | The statistic is right for the target. Its retro-check carries a selection effect (see below). |
| R6 | **partially resolved** (N3) | H1a/H1b split, V5b and V5c are right. Forward-vs-reverse measures only one part of the error (see below). |
| R7 | **partially resolved** (N2b) | The label condition counts splits, when it should count refinement checks and where they fall in the windows. |
| R8 | **resolved in wording; defeated at 6:1 by the gap-companion baseline** (N1, high) | — |
| R9 | **resolved** | — |
| R10 | **resolved** | — |
| R11 | **resolved** | — |

### The specific questions
**REFINE_BLOCKED and the driver gate plan.** The logic is sound.
- A refusal that coexists with a split is logged and re-tested at the next check. Persistent depth refusals therefore
  end as BLOCKED_DEPTH, not as a false CONVERGED.
- Cap refusals are evaluated per attempt. Within one check the result depends on the order of the attempts, but it is
  deterministic. That is acceptable, and it is the study-3 order.
- Gap: after a CONVERGED_LS_LIMITED-type check that returns SPLIT, the trust radius is still below 1e-8. In `s3.cpp`
  (lines 617-635) trust is carried, so the next inner solve would end LINESEARCH_FAIL at once. That is an artefact,
  not a property of B (N5).

**R2: is holding converged runs at their final value fair?** Yes, on the matched-budget question, which asks what
each method has produced after k cumulative iterates.
- A CONVERGED run has reached a fixed point of its own iteration. More iterates would not move it (for B, this
  includes C3). Its q is therefore its q at every later budget.
- Truncating at non-converged terminals is conservative: it can only produce "not evaluable".
- **The new problem is padding.** Suppose both runs converge early, say at 20 and 30. Then k_c = K = 200, and 8-10
  windows are constants. The verdict reduces to one comparison of finals, but it passes the "≥ 4 windows" and
  W = 20/25 robustness tests mechanically. Those tests exist to guard against trajectory noise; stationary windows
  manufacture agreement.

**R6: forward vs reverse summation.** It measures the right thing only for summation roundoff. It is not a bound.
1. **It does not measure per-node evaluation error.** That error is identical in both orders. For this objective it
   is plausibly the dominant term: ρ is built from metric second derivatives and may cancel. Above 1e-4 relative it
   is bounded by H2/H3; below that it is unmeasured.
   - τ_H and the H1b gate are therefore floors, not error bars.
   - A D1b "strict" with μ_min just above τ_H is not established to roundoff.
2. **e_H is the difference of two estimates that carry comparable error.** It is not the error of the value used.
   - With heavy cancellation both orders can err in a similar direction, so the difference can underestimate the
     error.
   - The direct, equally free measure is the forward sum minus a compensated (Neumaier, or `long double`) sum of the
     same stored contributions.
   - The standard a-priori bound (n−1)·u·Σ|contrib| is also free.
3. By Weyl, |Δμ_i| ≤ ‖ΔH‖₂, so τ_H = 10·e_H is the correct *form* once e_H estimates ‖ΔH‖₂.

**R7: the label.** It fires on "0 outer rounds within k_c". Two cases go wrong:
- B passes its first check (CONVERGED, n_outer = 0). Its set then *has* been tested by the refinement rule, and that
  is a legitimate B fixed point. The label wrongly says "not a formulation effect".
- B splits once at iterate 190 of 200. The label is dropped, but almost every window compares a fixed set.

The condition should depend on **refinement checks** (not splits) and on **how many verdict windows end after the
first check**.

**D4 retro-calibration: a tuning risk?** The threshold choice is low risk, but the statistic choice carries a
selection effect. Study 4 should report it as follows.
- **Threshold.** The separation on the archive is 0.3-4.6 % vs 87.9-98.2 %. Any threshold in about (5 %, 85 %) gives
  the same classification, so 10 % is not tuned in a way that matters.
- **Statistic.** I proposed u_w (R5) after seeing that study 3's A λ2 finals hid the pocket off the working nodes. The
  retro-check confirms that the rule flags the cases it was chosen for. It is a consistency check, not validation.
- **Not independent of the archive.** P2 `A4L2f`/`A6L2f` start *at* the archived designs with u_w = 98 % / 88 %. A
  "fixed-node exploitation" label on P2 is therefore close to inherited, unless the run moves far. Only P1 finals and
  any design that moves substantially are out-of-sample applications.
- **Required wording (N4):** "D4 statistic chosen post hoc (study-3 REPORT §6, Noether R5); retro-check on the
  motivating data is a consistency check, not a validation; P2 starts on the calibration designs (start u_w
  reported)". Report u_w at the start next to u_w at the final for every run.
- Also correct the v2 text "the threshold (10 %) is the study-3 value". The 10 % was study 4 v1's "B inner
  exploitation" threshold. If the study-3 rule had a 10 % branch, cite it. Either way, state where the number comes
  from.

### New problems introduced by v2 (required changes)
**N1 (high). The pre-window gap baseline voids the R8 safeguard at 6:1.**
From the archived traces (in-loop q, s0 = 90):
- A rises from 0.71 at iterate 87 to 2.02 at iterate 90. The start iterate is A's onset spike.
- Δ_pre = median_A − median_B over 66-90 is **−0.093 (W = 20) / −0.084 (W = 25)**. The point gap at the start is
  **+1.15**.

Consequence:
- G_w = Δ_w + 0.09, so every window that supports "B lower" (Δ_w ≥ 0.05) has G_w ≥ 0.14.
- "No window with Δ_w ≤ −0.05" implies no G_w ≤ −0.05.
- Therefore at 6:1, **"B lower" implies "gap grows"**, and the word "suppresses" is licensed by the D5w verdict alone.
  It would apply even if the inherited +1.15 gap shrank to +0.2.

At 4:1 the bias runs the other way: Δ_pre ≈ 2.28 vs point 1.62. That is conservative for "grows" and liberal for
"shrinks".

The λ1 calibration could not show this. It is a null case with stationary trajectories, and it tested only "shrinks"
false positives.

**Fix:**
- Require "gap grows" against **both** baselines: G_w computed with max(Δ_pre, Δ_0).
- Require "gap shrinks" against both: computed with min(Δ_pre, Δ_0).
- Otherwise "not established (baseline-sensitive)".
- Report both baselines and the numbers above in the protocol.
- State that the pre-window lies in study-3 L-BFGS iterates under B's per-iterate refinement.

**N2 (medium). D5w padding and the no-refinement condition.**
- (a) Report the number of windows in which either run is extended. If at least half of the verdict windows (at
  either W) contain extended values, append "(determined by converged final(s))". The ≥ 4-window and W-agreement
  conditions are then not cited as robustness.
- (b) Replace the R7 condition. Use "(no study-4 refinement check)" if B made 0 refinement checks within k_c. Use
  "(refinement-limited: n of m windows after the first check)" if fewer than 4 complete windows end after B's first
  check. A B run that is CONVERGED via PASS is not labelled.

**N3 (medium). R6 error measure.**
- Add the compensated-sum reference: e_H^c = ‖H_fwd − H_comp‖₂, likewise for e_g and e_F. Also report the
  (n−1)u Σ|contrib| bound. Use e_H^c, or max(e_H, e_H^c), in τ_H and in the H1b gate.
- Add to D1b: "τ_H covers node-sum roundoff only; per-node evaluation error is unmeasured below the H2/H3 resolution".
- Label "strict" finals with μ_min < 100·τ_H as **"strict (marginal)"**.
- V2 bit-identity is unaffected, because the forward sum is unchanged.

**N4 (low-medium). D4 reporting:** as above.

**N5 (low). Clarifications.**
- (a) Trust is reset to 0.5 at every SPLIT.
- (b) V5d adds a mixed case (split plus a depth refusal in one check → SPLIT, and the refusal is re-tested at the next
  check) and the OUTER_LIMIT case.
- (c) In D6x, check label 4 (inconclusive) before 2/3 when an arm is missing. Under the present ordering, one missing
  arm plus one "keeps" arm gives "path-dependent", which rests on a single arm. If that is intended, say so.

### Reproduction
```
cd pub-feasibility/warp_hypad
sha256sum study4/PROTOCOL.md study4/PROTOCOL-v2.md study4/calibration/* study4/driver/*; cat study4/logs/*.sha256
PY=<python 3.11+>
$PY -I study4/calibration/d4_retro_check.py study3   | diff --strip-trailing-cr - study4/calibration/d4_retro_check.txt
$PY -I study4/calibration/d5w_calibration_v2.py study3 | diff --strip-trailing-cr - study4/calibration/d5w_calibration_v2.txt
```
N1: the λ2 point gap and the pre-window baseline come from archived traces only. For d in {4, 6} with s0 = 84 / 90
and W in {20, 25}, compute q = rhostar_in / b over `study3/runs/{A,B}_d_L2/trace.csv`. Then compute
q_A(s0) − q_B(s0), and median(q_A) − median(q_B) over iterates s0−W+1..s0. Result (point, pre-window):
4:1 W=20 (1.617, 2.287); 4:1 W=25 (1.617, 2.273); 6:1 W=20 (1.146, −0.093); 6:1 W=25 (1.146, −0.084).

### Candidate lessons (proposed only)
- A change-based companion whose baseline is a window median can be dominated by the level rule when the start iterate
  sits on a transient. Test the companion against the point baseline on the *actual* start trajectories, not only on
  a null pair. (Scope: pre-registered window comparisons; evidence: N1 numbers above.)
- Padding converged runs as stationary is fair for matched-budget questions, but it manufactures robustness windows.
  Label verdicts that the padding decides. (Scope: matched-iterate designs; argued.)

### Verdict: **NOT CLEARED**
Remaining required changes, to be made in v3 with a new hash:
- **N1** (gap companion uses both baselines);
- **N2a/N2b** (padding label; R7 condition based on checks);
- **N3** (compensated-sum reference and D1b limitation wording);
- **N4** (D4 post-hoc reporting);
- **N5a** (trust reset at SPLIT).

N5b and N5c are recommended. None of these needs new calibration data. A v3 re-check can be limited to these items.
Step-0 code (V1-V5) may be written meanwhile, but it should include N3's compensated sum before V4 is run.

---

## v3 re-check (Noether, 2026-10-27; protocol review only, no study-4 result)

### Independence, integrity and what was run
- I did not draft v3 (Ada did). N1-N5 are my own findings, so this re-check tests fixes I requested. Nothing was
  repaired.
- `PROTOCOL-v3.md` SHA-256 is `768b6036...e70a3` and matches `logs/protocol-v3.sha256`. v1 (`a21e0ff6...f758bf64c`)
  and v2 (`b2c0af84...176885da`) are unchanged against `logs/protocol.sha256` and `logs/protocol-v2.sha256`.
- I read v3 sections 0-9 and 14, `diff PROTOCOL-v2.md PROTOCOL-v3.md` (every removed v2 line maps to a row in
  section 14), and `study3/PROTOCOL.md` lines 130-140 (the 10 % origin in N4 is correct).
- I checked the section-8 arithmetic by hand. Nothing was evaluated, and no archived trace was re-read. The N1 table
  matches my v2 numbers digit for digit.

### Verdicts on N1-N5
| ID | Verdict | Reason |
|---|---|---|
| N1 | **resolved** | "Grows" is tested on G^hi = Δ_w - max(Δ_pre, Δ_0) and "shrinks" on G^lo = Δ_w - min(Δ_pre, Δ_0). The thresholds are correct: 6:1 grows ≥ 1.196, shrinks ≤ -0.143 / -0.134; 4:1 grows ≥ 2.337 / 2.323, shrinks ≤ 1.567. Exclusivity holds (G^hi ≤ G^lo), even with one window. The licence table removes the 6:1 "B lower ⇒ grows" route, and "suppresses" now needs both baselines or P5. The 6:1 onset-spike correction in section 5 is accurate and honest. |
| N2a | **resolved** | n_ext(W) ≥ m(W)/2 at either W gives "(determined by converged final(s))", and the robustness tests are then not cited. |
| N2b | **resolved** | The labels are based on checks and k_chk1, with the PASS exception; n_check, k_chk1 and `outer.csv` are logged. Reading: "CONVERGED (via PASS)" includes CONVERGED_LS_LIMITED via PASS, as D1 does. |
| N3 | **resolved** | Neumaier reference, ê = max(e, e^c) in τ_H and the H1b/H1c gates, the a-priori bound, the limitation wording, the "strict (marginal)" band at 100·τ_H, and the V5e test. The forward sum is untouched, so V2 still holds. |
| N4 | **resolved** | The registered post-hoc wording is used. Start u_w is reported with final u_w. "(start on calibration design)" is applied correctly, including to the P1 B runs. Only `A4L2m` and `A6L2m` are out-of-sample. The 10 % origin is corrected (verified at study-3 PROTOCOL line 135). |
| N5a | **resolved** | Trust is reset to 0.5 at every SPLIT, the pre-split radius is logged, and V5d tests it. It affects B only, and it is declared. |
| N5b | **resolved** | V5d includes the mixed case and OUTER_LIMIT. |
| N5c | **resolved** | "Inconclusive" now comes first. Labels 0-3 are exhaustive and exclusive. |

### Builder choices
1. **"(marginal)" propagates without demotion: acceptable.**
   - D1b already says that a marginal "strict" is one where unmeasured per-node error could change the
     classification. Keeping the label with a suffix keeps more information than demoting it, and does not claim more.
   - Condition, already met by the text: the suffix must never be dropped. That includes λ-level statements (the
     weaker density's labels apply), summaries and the abstract. μ_min/τ_H must be reported next to it.
   - A "multiple minima (marginal)" must not be paraphrased as "multiple minima" anywhere.
2. **Neumaier in `double` as the reference: sound for what N3 asked.**
   - N3 did **not** ask the reference to capture per-node roundoff. No summation scheme can. Every sum (forward,
     reverse, Neumaier, `long double`) uses the same stored per-node contributions, so per-node evaluation error is
     identical in all of them and cancels in every difference.
   - N3 asked for two things: (a) a reference whose own error is negligible next to the forward-sum error, and (b)
     wording that per-node error is unmeasured below the H2/H3 resolution.
   - (a) holds. The Neumaier error is ≤ 2u·Σ|contrib| + O(n u²)Σ|contrib|. That is about 2/(n-1) ≤ 3e-4 of the
     a-priori forward bound for n ≥ 6,624, so e^c measures the error of the value used. The stated bound is
     conservative (the sharper form is u|S| + O(n u²)Σ|x|).
   - (b) is stated at H1b and with every D1b label.
   - Measuring per-node error would need per-node evaluation in higher precision (for example a `long double` OTI
     build at one design). That is a possible future descriptive check, not a requirement.
   - Proviso, already implied by V2/V5e: no `-ffast-math` or reassociation. The `build.sh` flags (`-O2
     -ffp-contract=off`) preserve compensated summation.
3. **N1 baselines are mine, not recomputed by the builder; V6 recomputes them and stops on mismatch: sound.**
   - The driver code is a separate implementation. V6 catches transcription errors, arithmetic errors and
     indexing-convention differences, and the 87/90 spike check pins the iterate indexing.
   - It fails safe: a mismatch stops the study before any threshold is used.
   - The tolerances (±0.001 on 3-decimal values, ±0.005 on 2-decimal values) match the quoted rounding.
   - **Reading (no change needed):** the registered table values are the ones used in G^hi and G^lo. V6 only
     verifies them. A window within 0.001 of a threshold is reported as such. If the builder intends the recomputed
     full-precision values instead, that must be declared before any D5w computation.

### Does v3 reopen R1-R11?
**No.**
- R1 (the C3 and outcome logic) is unchanged. The N5a reset acts only after SPLIT and does not touch C3.
- The R2 horizon and truncation are unchanged. N2a only adds a label.
- In R3, precedence is reordered by N5c, and the labels remain exclusive and exhaustive.
- R4 wording only gains "(marginal)". R5 now carries the N4 wording.
- R6 is strengthened (ê replaces e). R7 is replaced by N2b, which is strictly finer.
- R8 is tightened by the licence table.
- R9 and R10 are unchanged. The refresh rule is untouched, and the trust reset is not a refresh.
- R11 is unchanged. The unrequested P5 scheduling change keeps the total at ≤ 12 threads (10 + 2). Results are
  thread-independent, and P5 remains descriptive.
- The gradient now also accumulates reverse, Neumaier and `long double` sums. This does not alter the forward sum,
  and V2 tests that.

### Reproduction
```
cd pub-feasibility/warp_hypad/study4
sha256sum PROTOCOL.md PROTOCOL-v2.md PROTOCOL-v3.md; cat logs/*.sha256
diff PROTOCOL-v2.md PROTOCOL-v3.md | grep '^<'      # every removed v2 line is covered by section 14
sed -n 130,140p ../study3/PROTOCOL.md               # D4 10 % origin (line 135)
```

### Candidate lessons (proposed only)
- A compensated sum is the right reference for the error of a sum, but it is blind to errors in the summands. When
  a review asks for a reference, also state which error component it cannot see, so that the request is not read as
  covering that component. (Scope: roundoff estimation for assembled sums; evidence: the framing of the v3 builder
  question.)
- When registered thresholds depend on reviewer-computed numbers, a fail-safe recomputation gate plus an explicit
  rule for which values enter the rule (table or recomputed) avoids a post-hoc choice. (Scope: pre-registration.)

### Verdict: **CLEARED TO RUN**
- N1-N5a are resolved, and so are the recommended N5b and N5c. No v4 is required.
- The two readings above are interpretations of existing text, not changes:
  - the table values are used in G^hi and G^lo;
  - the PASS exception includes CONVERGED_LS_LIMITED.

  If the builder disagrees with either reading, the disagreement must be declared before any D5w computation.
- The gates V1-V6 still apply, and they must all pass before any run.
