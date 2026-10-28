# REVIEW (Vera): warp_hypad study 3

Reviewer: Vera. Date: 2026-10-27, 20:09–20:45 CDT. Builder: Ada. I did not write or modify any study-3 code,
data, or report text. I worked in a separate context from Ada. That reduces anchoring but is not a certificate.
The only file I wrote in `study3/` is this one. Reviewer scratch outputs are in the session scratchpad.

Inputs: PROTOCOL.md, REPORT.md (as saved by the main session), driver/ (s3.cpp, analyze.py, scripts), runs/,
assess/, results/, figures/ (trajectory_L2.png inspected), logs/; study2/driver/s2.cpp for the verbatim claim;
Enterprise at 477141d (read-only; `git status` clean).

## Verdict
- The **negative and procedural results stand**:
  - no run converged;
  - D2, D3, D4 and D6 are "failed (not converged)";
  - D5 is "not established" (non-discriminating);
  - the lambda1 results;
  - the T labels;
  - the reproduction of b and E0.
- I found no defect in the new driver code (section 3).
- The rerun reproduced the archived lambda2 4:1 pair bit-for-bit.
- The **A-versus-B contrast at lambda2 holds at 4:1 only**. At 6:1 it is an artefact of unequal iteration counts
  (blocking B1). The candidate manuscript wording in REPORT section 11 cannot be used as drafted (B2).

## 1. Pre-registration integrity
| Check | Result |
|---|---|
| sha256 of PROTOCOL.md vs logs/protocol.sha256 | equal (05fd4082...a125bfe), recorded 16:00:04 |
| PROTOCOL.md mtime | 16:00:00.88. Never modified afterwards |
| First study-3 computation | logs/check.log 16:04:50; s3.cpp written 16:03:39; pilot 16:07–16:11; lanes 16:11:25 |
| Pilot decision | Delta1 = 11.05 > Delta2 = 5.90, PASS; lambda unchanged |

Pre-registration preceded all computation. One inconsistency (optional O1): PROTOCOL's header says "written
15:59–16:25 CDT", which cannot be right for a file hashed at 16:00:04. REPORT says 15:59–16:00, which matches the
file system.

## 2. Reproduction (requested rerun)
Commands, run concurrently at 6 + 6 = 12 threads, with a safety deadline (now + 55 min) that was never reached:

    OMP_NUM_THREADS=6 study3/driver/s3.exe run --form A --ws 4:1 --start baseline_zero --lambda 3.75e-9 --kmax 133 --maxfev 3000 --deadline DL --diag-every 1 --out SCRATCH/A_4_L2
    OMP_NUM_THREADS=6 study3/driver/s3.exe run --form B --ws 4:1 --start baseline_zero --lambda 3.75e-9 --kmax 84 --maxfev 3000 --deadline DL --diag-every 1 --out SCRATCH/B_4_L2
    OMP_NUM_THREADS=12 study3/driver/s3.exe assess --y SCRATCH/RUN/y_final.txt --seed 1|2 --tag RUN --full 1 --out SCRATCH/assess.csv

Wall time: 20:09:32 to 20:35:18 (A) and 20:35:43 (B), so 26 min. The assessment of the 4 designs took another ~1 min.

| Item | A 4:1 lambda2 | B 4:1 lambda2 |
|---|---|---|
| y_final.txt, ys.txt (all iterates) | byte-identical | byte-identical |
| trace.csv columns 1–36, including C1, C2, C3 (only the timing columns differ) | identical | identical |
| inloop_minima.csv, refine_events.csv, set_changes.csv | identical | identical |
| termination counters (it, nfev, ngev, nodes, resets, gF, F, R, norm y) | identical; state ITERATION_LIMIT vs archived WALL_LIMIT, as expected | identical; same |
| assessment rows (both seeds; every field except path and seconds) and minima lists | identical | identical |
| q, r_E | 31.9598, 0.5885 | 1.0521, 0.6321 |

- **Thread-count independence.** The archive used 2 threads and the rerun 6. In addition, the pilot (12 threads)
  and production B 4:1 lambda1 (2 threads) have byte-identical first 21 ys.txt rows. The determinism claim holds.
- **Limitation.** I used the archived s3.exe and did not rebuild it, to stay within the thread budget. Its mtime
  (16:04:29) is after s3.cpp's (16:03:39), and check.log passes. That the binary corresponds to the source is
  inferred, not shown.

## 3. Code audit of s3.cpp (new code)
Verbatim claim. Diffing Elem, NodeSet, edges_of, rebuild_nodes, base_set, node_rho, residual, gradient and polish
against study2/driver/s2.cpp shows comment differences only.

- **In-loop search (lines 277–324).** Matches PROTOCOL 3.2:
  - 241x116 lattice on Dw (h about 0.028);
  - 8-neighbour minima with rho < 0, ranked;
  - compass polish to 1e-9 in batches of 8;
  - 1e-3 merge done serially in rank order (deterministic);
  - stops at 8 distinct or 64 polishes.

  Ties count as minima (non-strict test). Harmless.
- **Refinement edge rule (lines 327–375).**
  - Every leaf whose closed (xi, s) rectangle contains the point is tested, with the 5 % rule, depth <= 6 and the
    4x node cap.
  - Children are appended, so the untouched leaf indices stay valid. Re-testing happens through the guard loop.
  - A resolved neighbour is re-logged on each guard pass (log noise only). Correct.
- **No-reset L-BFGS (lines 551–633).**
  - The two-loop recursion is standard. Pairs are kept across set changes.
  - After a change, f and grad F are recomputed on the new set, so the next pair is consistent within that set.
  - Memory is cleared only on a non-descent direction, followed by p = -gamma_last g. Iterate 1 uses p = -g.
  - Pairs with s'y <= 1e-14 |s||y| are skipped (inherited behaviour, not mentioned in PROTOCOL).
  - The reset flag is 0 at every iterate of all 8 traces. Correct.
- **Tikhonov term.**
  - F = R + (lambda/2)|y|^2 and grad F = grad R + lambda y are correct.
  - The line search and the Armijo test use F on the current set.
  - The pen column is evaluated at the accepted y.
- **Stopping.**
  - C1' and C2' are as registered. The C2' history restarts after a set change with the new-set F.
  - STALLED is 30 consecutive iterates and resets on a change. C3 = no split.
  - In 8/8 traces C1' and C2' are never true. The final gradient/tolerance ratios (45x to 38,000x) are confirmed.
- **WALL_LIMIT.**
  - The deadline is checked only at the top of the loop, so the iterate in progress (search plus refinement)
    completes. The overrun is at most one iterate: B 4:1 lambda2 finished 47 s after the 19:41:13 deadline.
  - The final design is fully reproducible with --kmax (section 2).
- **Assessment (lines 723–974).**
  - Halton bases 2/3 with a seeded Cranley–Patterson shift, 2^17 points.
  - The Lipschitz estimate min-corner - sqrt(2) L h is 2x more conservative than a nearest-corner bound, and it is
    a heuristic, as stated.
  - Flagging (est < best - 10 %, at most 4096 cells) and the 0.00125 sub-lattice match PROTOCOL 7.
  - Starts are deduplicated at 0.01, best 128. Nelder–Mead is confined by +inf outside Ab.
  - Outside-Ab independence is tested by exact equality at 34,266 points.
  - The lattice endpoints (2.1, 3.1) evaluate inside Ab in floating point (checked).
- **Undeclared deviation.** The Nelder–Mead initial simplex is 0.5h = 0.0025, with restarts at 0.25h. PROTOCOL 7
  says 0.005. See I4.

## 4. Findings

### Blocking (for the claims and the manuscript wording, not for the negative verdicts)
**B1. The lambda2 A-versus-B contrast at 6:1 is confounded by unequal iteration counts.**
- B's iterates cost more (refined nodes, extra line-search evaluations). At equal wall time, B reached iterate 90 at
  6:1 and A reached 126; at 4:1 the counts were 84 and 133.
- results/late_phase.txt pools A over iterates 40–126 and B over 40–90.
- The table below gives a reviewer-computed comparison over a common iterate window, from trace.csv and using
  in-loop q. It is also post hoc.

| iterate window | run | median q | max q | fraction q >= 1.05 | q at window end |
|---|---|---|---|---|---|
| 40–84 | A 4:1 | 2.69 | 3.36 | 0.98 | 2.67 |
| 40–84 | B 4:1 | 0.76 | 1.67 | 0.13 | 1.05 |
| 40–90 | A 6:1 | 0.92 | 2.02 | 0.20 | 2.02 |
| 40–90 | B 6:1 | 0.87 | 1.40 | 0.31 | 0.87 |
| 60–84 | A 6:1 vs B 6:1 | 0.99 vs 1.07 | – | 0.24 vs 0.64 | – |

- At 4:1 the contrast survives matching.
- At 6:1, within the common window, B has more iterates with q >= 1.05 than A, and the medians are close.
- A's 6:1 deepening (q 2–5.5) comes almost entirely from iterates 89–126, which B never reached.

These REPORT statements are therefore not supported at 6:1:
- sec. 1.4 ("B does so far less");
- the 6:1 sentence in sec. 5;
- sec. 11's "final q of 0.87–1.05" for B against "4.7–32x" for A.

Fix: compare at equal accepted iterates. Restrict any A > B statement to 4:1 lambda2, with the qualifiers one
start, non-converged and post hoc.

**B2. The candidate manuscript wording (REPORT sec. 11) is not usable as drafted.**
1. The 6:1 contrast (B1).
2. The A-versus-B contrast is not pre-registered. D5 was the registered test and did not discriminate. The sentence
   that states the contrast must itself say "exploratory, not pre-registered". Placing the D5 caveat in a later
   bullet is not enough.
3. "Deepened it 4.7–32x" quotes the last iterate of a quantity that was still rising: A 4:1 in-loop q went 3.1 to
   32 over its last 25 iterates, and REPORT itself flags lambda2 as iteration-unstable. Say "at the last iterate
   (q still increasing; range over the last 20 iterates 4.3–32)".
4. "Final q of 0.87–1.05" hides that B 4:1's final iterate (q = 1.052, r_E = 0.632) meets the registered
   per-design concentration criterion.
5. "Reached a similar E_neg reduction (63 %)": A 4:1 ended at 59 %, with r_E falling while q rose (an
   exploitation signature). State both values.

Suggested replacement for the lambda2 bullets: "At a budget of about 38, after 84–133 non-converged iterations,
fixed-node minimisation at the coarser working set produced worst points 2.7–32x deeper than the baseline over
iterations 40–133. Adaptive re-insertion stayed near the baseline depth (median 0.76x) over the same iteration
window, with excursions to 1.67x. At the finer working set no consistent difference was seen over the iterations
both runs reached. These comparisons were not pre-registered: the registered matched-progress test was not
discriminating, because both formulations crossed its E_neg levels before they diverged."

### Important
- **I1. late_phase.txt provenance and label.**
  - No archived script produces it (analyze.py has no late-phase code), and it was written at 19:55, during the
    incomplete first analysis.
  - The file has no "post hoc / not pre-registered" label. The window (iterates >= 40) was chosen after seeing the
    data.
  - Its q is in-loop q (a lower bound), but the REPORT sec. 5 column says only "median q".
  - REPORT labels the comparison post hoc in sec. 1.5, 5 and 8. Headline sec. 1.4 and sec. 11, which rest on the
    same comparison, do not.
  - Fix: archive the generating code, add a header line to the file, and label every A-versus-B sentence.
- **I2. The onset statement (sec. 1.4) is inaccurate.**
  - A 6:1: in-loop q >= 1.05 at iterates 32–34, 42–43, 64–65 and 70–73, then sustained over 89–97 and 104–126. The
    rise starts at about 88, not about 110.
  - A 4:1: q >= 1.05 at 27–34, before 40.
  - All four lambda2 runs, B included, share concentration at iterates 10 and 21 (q 1.37–1.49 at r_E about 0.61,
    norm y about 32). This is not mentioned anywhere. Sec. 1.6 is correct only because it limits itself to lambda1.
- **I3. "Slow convergence inside a bounded region" (sec. 1.2) is evidenced at lambda1 only.**
  - Median gradient/C1' tolerance over iterates 30–60, then over the last 30 iterates:

    | lambda | run | iterates 30–60 | last 30 iterates |
    |---|---|---|---|
    | lambda1 | A 4:1 | 16,152 | 111 |
    | lambda1 | other three runs | falls 7–17x | |
    | lambda2 | A 4:1 | 15.8k | 9.5k |
    | lambda2 | B 4:1 | 16.8k | 18.4k |
    | lambda2 | A 6:1 | 15.3k | 8.7k |
    | lambda2 | B 6:1 | 15.5k | 16.7k |

  - So at lambda2 the ratio is flat or rising.
  - A bounded norm y does not imply slow convergence. Wandering between basins (A 4:1's q jumps) fits the data
    equally well.
  - Say: "not converged; no clear progress toward C1' at lambda2".
- **I4. Undeclared Nelder–Mead simplex deviation** (section 3). Low consequence: the in-loop search and the
  assessment agree at 79 of 80 iterates. It should still be declared in REPORT sec. 2.
- **I5. The conditioning sentence (sec. 1.2) overreaches.**
  - gamma spans 1.07e3 to 2.45e6 over all iterates; 3e4–4e5 is the late-phase range only.
  - 1/(gamma lambda) is about 10^2–10^3 at lambda1 and 10^3–10^4 at lambda2, not 10^3–10^4 for both.
  - lambda is not a lower bound on the eigenvalues of the Hessian of F when R is non-convex, so this is not a
    condition-number estimate.
  - Say instead: "the curvature along recent steps exceeds lambda by 10^2–10^4".

### Optional
- **O1.** The PROTOCOL header says "written 15:59–16:25" (section 1).
- **O2.** REPORT sec. 4 note: the working-node rho_min at the A 4:1 lambda2 final is -5.6e-4 (rhomin_w_post,
  final_diag.csv), not "about -1e-3". This strengthens the point: the working nodes see q about 0.59 where the
  assessment finds 32.
- **O3.** REPORT sec. 7, "Flagged cells appeared only in lambda2 designs", holds at the finals only. The lambda1
  iterates A 4:1 @50 and @100 had 20 and 18 flagged cells.
- **O4.** At lambda1, A 6:1 has one in-loop iterate with q = 1.041 (iterate 55), which is "marginal" by the
  registered labels. REPORT gives "<= 1.04", which is accurate; the manuscript's "beyond an early transient" could
  mention it.
- **O5.** List the L-BFGS skipped-pair rule in config.txt and in any future protocol.

## 5. Scope check
REPORT sec. 1–8 keep to:
- Eulerian rho;
- start S1;
- these lambda values and densities;
- non-converged iterates.

The final-state sentences carry the T labels. The only scope slips are:
- the 6:1 A/B generalisation (B1);
- the unlabelled post-hoc contrast in sec. 1.4 and 11 (I1, B2);
- "slow convergence" at lambda2 (I3).

Sec. 11 opens with the correct scope phrase, and "no run converged" is right. The admissibility caveat in sec. 4 is
appropriate.

## 6. Verified as stated
- 0/8 converged; C1' and C2' were never true.
- No resets.
- Splits 35/43/16/29.
- The D5 table, recomputed from the traces. A and B are identical up to the first B split: iterate 9 at 4:1, 12 at
  6:1.
- T labels and counts.
- The late-phase numbers as tabulated (recomputed from the traces).
- In-loop versus assessment: 79/80 agree within 1e-6; A 6:1 lambda2 @100 is 1.75 % deeper.
- 11 D7 events.
- Under-reporting <= 2.0 %.
- R4/R6/R8 and E8/E16 at the finals.
- Seed agreement.
- Zero outside-Ab mismatches.
- The study-2 norm figure ("A 95–139"), which comes from study2/REVIEW-gauss.md.

## 7. Candidate lessons (for the main session; proposed only)
1. When formulations differ in per-iterate cost, a fixed wall budget yields unequal iteration counts. Trajectory
   comparisons should be made over a common iterate window, or registered as such.
   - Evidence: B1.
   - Scope: optimisation studies with adaptive quadrature.
2. Matched-progress levels should be placed using a pilot of where the trajectories diverge.
   - Evidence: D5, where |dq| <= 0.0011 at every registered level.
