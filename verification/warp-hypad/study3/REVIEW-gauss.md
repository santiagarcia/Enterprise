# REVIEW (Gauss): warp_hypad study 3, numerical assessment

Reviewer: Gauss (numerical methods), 2026-10-27. Builder: Ada. Affect layer: enabled, low influence (concern about
my own share in the D5 design prompted the extra self-check in section 2).

**Independence.** I did not write `s3.cpp`, run anything in study 3 or edit any study-3 artifact apart from this file.
I am **not** independent of the design: PROTOCOL sections 2, 4 and 8 (including D5's r_E levels 0.40/0.50/0.60 and
the 400-iterate cost estimate) come from my study-2 review (`../study2/REVIEW-gauss.md` section 6). Where study 3
failed because of that design, the error is mine and I say so. I did not rerun the study and did not audit the code
changes (Vera's lane). All my numbers come from the archived `runs/*/trace.csv` and `linesearch.csv` (section 5).

## 0. Summary
- **Why nothing converged.** The regulariser worked: the design norm is bounded and the λ1 A runs are near-stationary
  (A 4:1 λ1 is 45× above the C1' tolerance and falling by a factor 0.973–0.979 per iterate). The failure is a budget
  and solver problem, not drift. L-BFGS m = 8 could plausibly finish the λ1 A runs (about 180–450 more iterates). It is
  not a credible route for the λ2 runs (gradient 5,700–13,000× above tolerance while q is still moving) or for any B
  run (B splits every few iterates, so C2'/C3 keep resetting). Next step: measure the Hessian. Then use a dense quasi-Newton
  solver initialised from that Hessian, and give B a nested outer loop.
- **Condition estimate.** "10^3–10^4" is a plausible order of magnitude, not an estimate the γ data support. γ
  only brackets the curvature. The rigorous lower bound from the logged γ is κ ≳ 10–50. The 10^3–10^4 figure assumes
  λ_min(∇²F) ≈ λ, which nobody has checked.
- **D5 mis-targeting** was foreseeable from my own study-2 evidence (section 2). The REPORT handles it honestly.
- **The late-phase comparison is not iterate-matched (G-B1).** At λ2 the A runs had 126–133 accepted iterates and the
  B runs 84–90. At 6:1, A and B are indistinguishable inside their common iterate window 40–90, because A's 6:1
  deepening starts after iterate ~105, beyond B's horizon. The 4:1 separation does survive iterate matching and
  equal-working-R matching.
- **λ dependence.** "No deepening at budget ≈17" is robust for these near-stationary trajectories. "Strong A
  deepening at ≈38" holds only as a 4:1 trajectory statement (the size is a wall-clock snapshot). Calling ‖y‖ the
  cause is an interpretation. A cheap cross-warm-start test can discriminate (section 2(4)).

## 1. Findings ranked (severity, likelihood, consequence)

### G-B1 (blocking for the candidate manuscript wording; the REPORT's decision table is unaffected)
**The λ2 A-vs-B contrast at 6:1 is an iterate-horizon artefact.** It holds at 4:1 only.

In-loop q, recomputed from `trace.csv`:

| λ2, density | window | A: median q, max q, frac q ≥ 1.05 | B: median q, max q, frac q ≥ 1.05 |
|---|---|---|---|
| 4:1 | it 40–84 (common) | 2.69, 3.36, 0.98 | 0.76, 1.67, 0.13 |
| 4:1 | it 55–84 | 2.97, 3.36, 1.00 | 0.79, 1.67, 0.20 |
| 6:1 | it 40–90 (common) | 0.92, 2.02, 0.20 | 0.87, 1.40, 0.31 |
| 6:1 | it 61–90 | 0.98, 2.02, 0.27 | 1.06, 1.40, 0.53 |

**Equal working R** (first iterate with R ≤ level). These use different quadrature functionals for A and B; see O3.
- 4:1, at B's final R = 1.0552e-4: A reaches it at it 64 with q 2.63; B at it 84 with q 1.05.
- 6:1, R ≤ 1.055e-4: A it 74, q 1.03; B it 84, q 0.98.

**Consequences:**
- REPORT §1.4, §5 and §11 present "A deepens strongly, B far less" at λ2 as a two-density result. At 6:1 it compares
  A's iterates 105–126 with a B that stopped at iterate 90.
- At 6:1 the honest statement is: "indistinguishable within the common window; A deepened later (from iterate ≈105),
  beyond B's horizon".
- The same applies to "final q_B 0.87–1.05" against "q_A 4.7–32" in §11, and to D4's descriptive
  q_A(4:1) − q_A(6:1) = 27.2. Both compare non-stationary values taken at different iterates.

**Cause:** the run budget was set by wall clock, and B costs 45–54 s per iterate against A's 33–35 s.

**Minimum fix:** report the common-window table above (or recompute it from the assessed iterates) and restrict the
λ2 contrast to 4:1. Study 4 should budget by accepted iterates, not wall time.

### G-I1 (important). The condition estimate is stated more firmly than the data allow
**What γ bounds.** For a pair (s, y), γ = s'y/y'y satisfies 1/γ ∈ [λ_min(H̄), λ_max(H̄)], where H̄ is the averaged
Hessian along the step. γ is weighted towards the large eigenvalues, because y'y = s'H̄²s.

**What the logged γ supports.** Over the last 50 iterates of each run, γ spans 1.9e4–2.9e5 at λ1 and 2.5e4–1.5e6 at
λ2. This gives:
- λ_max ≳ 3–5e-5;
- λ_min ≲ 1–4e-6;
- hence κ ≳ γ_max/γ_min ≈ 10–50. That is the only rigorous statement.

**What the REPORT assumes.** "10^3–10^4" divides 1/γ by λ, i.e. it assumes some direction has R-curvature ≈ 0
(λ_min ≈ λ).
- This is plausible: 108 B-spline coefficients, many with weak influence on ρ < 0.
- It is unverified: λ_min could be much larger than λ, or ∇²R could be indefinite where the iterate is not near a
  minimiser (likely for A 4:1 λ2, where q is still growing).

**A better empirical proxy** is the contraction of ‖∇F‖ on runs without set changes. Over the last 50–100 iterates:
- A 4:1 λ1: r = 0.973–0.979;
- A 6:1 λ1: 0.986;
- A λ2: 0.977–0.992.

Reading r through the CG bound r ≈ 1 − 2/√κ gives κ_eff ≈ 6,000–9,000 for A 4:1 λ1; the steepest-descent bound
r ≈ 1 − 2/κ gives about 90. L-BFGS sits between these, so the data are consistent with 10^2–10^4.

**Fix:** reword as "apparent; bounds only", and measure the Hessian (section 3, step 0).

### G-I2 (important). B cannot meet C1' ∧ C2' ∧ C3 under the current loop
- **Splits keep coming.** Splits occur up to the last iterate at λ2: B 4:1 had 28 split iterates (last at 84 of 84);
  B 6:1 had 19 (last at 80 of 90). At λ1 they are sparse but persistent: B 4:1 last split at 166, B 6:1 at 146.
- **Each split resets C2'.** C2' counts only iterates since the last set change, so it needs ≥ 10 split-free iterates
  and then a 1e-8 decrease.
- **The gradient stops falling after splits.** The ‖∇F‖ contraction over the last 50 iterates is ≥ 1.0 for B 4:1 λ1
  (1.021) and B λ2 (1.003). B 4:1 λ1 ended 4,873× above tolerance, against 45× for A 4:1 λ1.
- **This is structural.** B minimises a piecewise-defined objective that changes whenever the in-loop search finds a
  new polished minimum. A larger budget alone does not fix it.
- **Remedy:** a nested formulation:
  1. converge the inner problem on a frozen node set;
  2. run the refinement check;
  3. repeat until no split (C3), with an outer cap.

  The fixed point of this loop is what "B at convergence" should mean.

### G-I3 (important). D5 was mis-targeted, and this was foreseeable. The error is in my study-2 proposal
Section 2(2).

### G-I4 (important). The λ dependence is confounded with the path
Section 2(4).

### Optional
- **O1. My study-2 cost estimate was wrong.** "≤ 400 iterates to C1'" was unverified. The pilot criterion (Δ2 < Δ1
  for ‖y‖) tested only boundedness, not the convergence rate. Study 4 should pilot the rate:
  log‖∇F‖ slope over 20 iterates after the warm start.
- **O2. Onset wording.** REPORT §1.4 says A 4:1 λ2 deepening starts "from about iterate 40". In-loop q is 1.45 at 40,
  1.15 at 55 and ≥ 2.7 from 60 onwards. "Sustained from ≈ 57–60" is more accurate.
- **O3. Equal-R matching compares different functionals.** A's R and B's R use different node sets (B's set is
  refined). Match instead on a formulation-neutral functional evaluated on a common set, e.g. R_ref = R on the 8:1 Dw
  validation set, which costs about 7 % of an iterate already (cost rule).
- **O4. C1' is strict relative to the decision quantity.** At A 4:1 λ1:
  - in-loop q was stable to 0.036 over the last 20 iterates;
  - ‖∇F‖ = 3.1e-9 against ‖∇R‖ = 6.3e-7 ≈ λ‖y‖, i.e. the penalty and the data gradient balance to 0.5 %;
  - with κ unknown, the distance to y* cannot be bounded.

  Once a Hessian is available, also report the Newton step ‖H⁻¹∇F‖ as a pre-registered distance-to-stationarity
  measure. Do not relax C1' post hoc.

## 2. Answers to the questions

### (1) Why no run converged; can L-BFGS m = 8 reach C1'/C2'?
**Not drift.** ‖y‖ is bounded (16.5–18.8 at λ1, 34–41 at λ2), and F decreases monotonically.

**Status at WALL_LIMIT:**

| Run | ‖∇F‖/tol_C1' | ‖∇F‖ contraction/it (last 50) | Naive extra iterates to C1' | ≈ hours at the study-3 rate |
|---|---|---|---|---|
| A 4:1 λ1 | 45 | 0.973 | ~180 | 1.6 |
| A 6:1 λ1 | 689 | 0.986 | ~450 | 4.4 |
| A 4:1 λ2 | 12,986 | 0.982 | ~530 | 4.9 |
| A 6:1 λ2 | 5,671 | 0.977 | ~380 | 3.7 |
| B (all) | 836–38,112 | ≥ 1.0 (3 of 4); B 6:1 λ1 0.983 (but 11 split iterates, last at 146) | not defined (splits) | — |

**C2'.** A 4:1 λ1 has a 10-iterate relative decrease of 1.4e-7 against 1e-8 required. F − F* scales like ‖∇F‖², so
this needs about 50 iterates and is not the binding criterion.

**Caveats on the projections:**
- They assume the asymptotic linear rate already holds.
- That is reasonable for λ1 A, where q and ‖y‖ are stationary.
- It is not credible for λ2 A: A 4:1 q tripled in the last 15 iterates while R kept falling 6e-5 per iterate, so the
  iterate is still in transit into a pocket.

**Verdict:**
- L-BFGS m = 8 is feasible for A 4:1 λ1, marginal for A 6:1 λ1, and not a sound plan for λ2 or for B.
- m = 8 in 108 dimensions is an unnecessary handicap. A dense BFGS inverse (108×108) costs O(n²) ≈ 10^4 flops per
  iterate, negligible next to a 28–37 s gradient. It keeps the full curvature history.

**Newton or preconditioning:**
1. **Exact Hessian by FD of the OTI gradient** (central, with a step sweep, section 3): 216 gradient evaluations.
   That is about 100 min at 2 threads for 4:1 (t_grad ≈ 28 s). Thread scaling is unmeasured; if it is near-linear,
   expect about 20–35 min at 6 threads. This is too expensive to repeat at every Newton step, but cheap as an
   initialiser or as a refresh every ~25 iterates.
2. **OTI order-2 Hessian blocks.** The gradient uses 36 passes of 3 directions. All 108×108 second derivatives need
   36 diagonal blocks plus 630 cross-block passes, each with 6 order-2 directions. Per-pass cost is unknown (several×
   an order-1 pass), so this is probably comparable to or dearer than FD. Its advantage is no truncation error.
   It is worth it only if the FD sweep shows step-size sensitivity.
3. **Cheap PSD part.** R = Σ w s(ρ)² has the exact part H₁ = 2 Σ w (s'² + s s'') ∇ρ∇ρᵀ, which is PSD because s² is
   convex in ρ. It comes from the per-node ρ derivatives that the existing gradient passes already compute and then
   discard (the `rho[p]` OTI values in `gradient()`).
   - It omits Σ w s s' ∇²ρ. That term is not small: the residuals do not vanish, and ρ is nonlinear in y through the
     lapse exponential.
   - Use H₁ + λI as a **preconditioner / BFGS initialiser**, not as a Newton Hessian.
   - Note: s s'' is comparable to s'² near ρ ≈ 0 for τ = 1e-3, so plain Gauss–Newton (s'² only) would be a poor
     approximation.

**Recommendation.** First measure the Hessian once at two archived finals (step 0). Then run dense BFGS initialised
with the measured H⁻¹ (made positive definite by eigenvalue clipping at λ), with an FD-Hessian refresh if the BFGS
model degrades. Expect superlinear convergence in tens of iterates if the Hessian is PD there. If the Hessian at
A 4:1 λ2 is indefinite, that iterate is not near a minimiser, which is itself a reportable result.

### (2) Was the D5 mis-targeting foreseeable? A better design
**Yes, and I should have caught it.**
- My study-2 B1 table showed A reaching B's final r_E within 14–37 iterates, with deepening only after r_E had
  plateaued.
- A and B share the trajectory exactly until B's first split (iterate 9 at 4:1, 11 at 6:1 in study 3).
- First-crossing levels at 0.40/0.50/0.60 were therefore bound to fall in the shared fast-descent phase (iterates
  10–15). In study 3, |Δq| ≤ 0.0011 there.
- First crossing is also the wrong statistic for a non-monotone r_E: A 4:1 λ2 r_E falls from 0.631 to 0.589 while q
  deepens.
- Ada applied the rule as registered; this is a design error, not an execution error.

**Matched-progress design for study 4.** All rules apply only after divergence, i.e. after the warm start.
- **(a) Primary, if convergence is reached:** D3/D4 at convergence. Matching then becomes secondary.
- **(b) Equal-iterate windows** (registered trajectory statistic):
  - budget by accepted iterates, identical K for A and B;
  - compare median q, max q and frac(q ≥ 1.05) over fixed windows [k, k+W] (e.g. W = 25) of the warm-started runs;
  - decision: "B suppresses" if median q_A − median q_B ≥ 0.05 in ≥ 2 of 3 windows at both densities.
- **(c) Equal reference-integral progress:**
  - match on R_ref (R on the common 8:1 Dw set) at levels defined relative to the warm-start value, e.g. R_ref ≤
    (1 − {0.5, 1, 2} %)·R_ref(y_final);
  - compare q at the first iterate reaching each level, but only for levels reached by both;
  - this is formulation-neutral. A's exploitation makes its R_ref non-monotone, so also report best-so-far R_ref.
- **Equal ‖y‖:** drop it. Both formulations sit at the same norm plateau (λ-bounded), so it does not discriminate.

### (3) Is the late-phase comparison sound and fairly qualified?
- **Fairly labelled:** post hoc, descriptive, pooled, with A's iteration instability noted. The REPORT is candid that
  the table was chosen after seeing the data.
- **Missing qualifier:** unequal iterate horizons (G-B1). This changes the 6:1 conclusion, so as worded the comparison
  is not sound at 6:1. The 4:1 conclusion survives both iterate matching and equal-working-R matching.
- **Smaller points:**
  - The "λ1, all four" row pools very different run lengths (175–245); harmless, because nothing deepens.
  - The REPORT says "B reaches equal or slightly better r_E". That is true at the finals (0.632 vs 0.589 at 4:1), but
    those finals are again at different iterates.

### (4) λ dependence: robust or artefact?
**λ1 "no deepening" is robust within scope.**
- All four runs reach near-stationarity with q = 0.59–0.60, iteration-stable to ≤ 0.036.
- They agree across densities and formulations (both lattices) and agree with the assessment.
- R4/R6/R8 agree to 1e-3. A 4:1 λ1 is within 45× of C1' with ‖∇F‖ still falling.
- A local minimiser with q ≈ 0.6 near these iterates is likely, but not proven.

**λ2 "strong A deepening" is a 4:1 trajectory statement.**
- It is real in the sense that F decreases while q deepens and the 8:1 quadrature gap grows to 20 %. That pattern is
  exploitation, not a line-search accident.
- Its size (q = 32) is a snapshot at a wall-clock cutoff: q went 3.4 → 32 in the last ~15 iterates.
- At 6:1 the onset is later (≈105). That fits exploitation needing narrower inter-node pockets (interpretation), but
  study 3 cannot say whether 6:1 converges deep or shallow.

**The causal reading "budget ≈38 permits deepening" is confounded with the path.**
- The λ1 and λ2 trajectories separate by iterate ~12 (‖y‖ 15 vs 26 at iterate 15).
- So the two runs reach different regions, not one region under two budgets.
- Discriminating test (cheap, uses the archives): warm-start A 4:1 at λ1 from `A_4_L2/y_final`, and A 4:1 at λ2 from
  `A_4_L1/y_final`.
  - If the first loses the pocket (‖y‖ → ~17, q → ~0.6) and the second develops one, the budget reading is supported.
  - If the first keeps its pocket, or the second stays shallow, there are multiple minimisers and the dependence is
    partly path-dependent.

**Wording I would accept:** "for this start, A trajectories at λ2 = 3.75e-9 developed node-exploiting pockets at 4:1
(from iterate ≈ 57) and later at 6:1 (≈ 105); no λ1 trajectory did before near-stationarity".

## 3. Proposed study 4 (minimal; to be pre-registered by the builder)
**Question:** at convergence (or at an iterate-matched budget) at λ2, does B suppress the pocket A finds, and is the λ
effect a budget effect?

**Solver and protocol changes (declared):**
- dense BFGS (no m limit), initialised from the step-0 Hessian (clipped PD at λ);
- the same Armijo/trust logic;
- B nested: inner C1' ∧ C2' on a frozen set, then a refinement check, outer cap 10;
- budget in accepted iterates (K = 300 per run), with a wall cap as a safety net only;
- C1', C2' and C3 unchanged.

**Steps:**

| Step | Content | Cost (assumptions: t_grad 28 s at 2 thr for 4:1 and 37 s at 4 thr for 6:1, from study 3; thread scaling unmeasured) |
|---|---|---|
| 0 | FD Hessian (central differences; step sweep h ∈ {1e-3, 1e-4, 1e-5}·max(1, \|y_i\|) on 4 columns first) at `A_4_L1/y_final` and `A_4_L2/y_final`. Report eigenvalues, κ, λ_min/λ, indefiniteness. Diagnostic, not decisional | ≈ 2 × 216 gradients; ≈ 1–1.5 h wall at 6 thr each |
| 1 | Warm-started λ2 runs A/B × 4:1/6:1 from archived `y_final` (4 runs), each with its own step-0-style initial Hessian (6:1 runs: 108 forward gradients) | Hessians ≈ 0.5–1.5 h; ≤ 300 iterates ≈ 2.5–4 h at 3 thr each, concurrent. Probably fewer if convergence is superlinear |
| 2 | Cross warm starts A 4:1: λ1 from `A_4_L2`, λ2 from `A_4_L1` (λ-path test, section 2(4)) | 2 runs × ≤ 150 iterates ≈ 1.5 h, concurrent with step 1 if threads allow |
| 3 | Assessment (study-3 method) of finals plus window/level iterates | ≈ 1 h at 12 thr |

**Total:** about 5–7 h wall at ≤ 12 threads. Optional: λ1 B to convergence (the λ1 A runs are near-done and could be
finished with plain L-BFGS from `y_final` in 1.6–4.4 h).

**Decision rules:**
- D1/D3/D4 as in study 3 for converged runs.
- D5' as in section 2(2)(b, c).
- D6' cross-warm-start outcome:
  - "budget-dependent" if λ1-from-pocket returns to q < 1.05 and λ2-from-shallow reaches q ≥ 1.05;
  - "path-dependent / multiple minima" if either keeps its starting character;
  - otherwise "inconclusive".
- Trajectory labels T kept.

## 4. What I did not check
- The `s3.cpp` changes (search, edge rule, no-reset L-BFGS, assessment). The only code I read is `gradient()` and the
  L-BFGS loop, to confirm the 3-direction OTI chunking and the γ definition (γ = s'y/y'y of the newest pair).
- No Hessian was computed: the compute budget is shared with Vera's rerun and possibly Abaqus.
- Thread scaling of t_grad is unmeasured, so all cost figures carry about a ×2 uncertainty.

## 5. Reproduction
From `study3/runs/` with Anaconda Python (`-I`):
- per run, read `trace.csv`;
- q = `rhostar_in` / (−9.60726526e-4), r_E = 1 − `Eneg_v8` / 4.634569e-2;
- tol_C1' = 1e-5 · F / max(`ynorm`, 1);
- the contraction rate is exp(slope) of a least-squares fit of ln `gF` against `it` over the last 50/100 rows;
- split iterates are rows with `splits` > 0;
- common-window statistics restrict both runs to it ∈ [40, min(n_A, n_B) − 1];
- equal-R uses the first row with R ≤ level.

The scripts used are in the reviewer's session scratchpad and are not archived. The procedure above is complete
enough to re-derive every number in this review.
