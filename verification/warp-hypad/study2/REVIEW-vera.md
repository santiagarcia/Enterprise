# REVIEW (Vera): warp_hypad study 2

Reviewer: Vera, 2026-10-27. Builder: Ada.
- **Independence:**
  - I worked in a separate context. I received the task, PROTOCOL.md, REPORT.md and the artifacts, but none of Ada’s
    reasoning beyond REPORT.md.
  - I did not read REVIEW-gauss.md.
  - I made no edits to the builder’s artifacts. My scratch builds and outputs are in the session scratchpad and are
    not part of the study.
  - For the rerun I used the builder’s driver source and assessment code, rebuilt by me. A separate context reduces
    anchoring, but it does not by itself make the verification independent.
- **Scope:** everything below concerns the Eulerian ρ = T(n,n) at t = 0, family F2A, and non-converged iterates.

## Verdict
- **The numbers reproduce bit-for-bit, and the code does what PROTOCOL §3 describes.**
- **"A exploitable" is supported**, as a statement about iterate 150.
- **The candidate manuscript claim (REPORT §10) is NOT supported as worded.** B’s own trajectories meet the
  pre-registered concentration criterion at non-final iterates, and §5/§6 under-report this (finding B1).
- The claim can be supported after it is narrowed to the final (150th) iterates and the B trajectory result is
  disclosed.

## 1. Rerun of the 4:1 S1 matched pair (about 21 min of runs plus 3 min of assessment)
- **Build:** I rebuilt `s2.exe` from `driver/s2.cpp` with the flags in `build.sh`, against `Enterprise` at 477141d
  (`git status` clean) and the existing `build-oti/` libraries, into the scratchpad.
- **Runs:** A and B ran concurrently with OMP_NUM_THREADS = 11 each. The originals used 6 threads, next to 3 other
  processes.
- **Bit-identical to the originals:** `y_final.txt`, `ys.txt`, `linesearch.csv`, `trace.csv` (all columns except
  wall time), `refine_events.csv` and `set_changes.csv` (B), and `termination.txt` (except wall_s).
  - A: ITERATION_LIMIT, it = 150, nfev = 176.
  - B: ITERATION_LIMIT, it = 150, nfev = 171, nodes = 7872, last_set_change = 147.
  - This also confirms D1: results do not depend on the thread count or on scheduling.
- **Assessment** (seeds 1 and 2, plus the baseline):

| Quantity | A 4:1 S1 | B 4:1 S1 |
|---|---|---|
| ρ*_D | −2.25594096e-2 | −6.59701955e-4 |
| q | 23.4816 | 0.6867 |
| r_E | 0.5586 | 0.6114 |
| check worst | −1.72208e-2 | −6.42484e-4 |
| under-report fraction | 0.237 | 0.026 |
| under-report flagged | yes | no |
| axi check | PASS (ratio 7.0e-14) | PASS (ratio 2.8e-13) |

  - Baseline: b = −9.60726526e-4 at (0, 2.00134), E0 = 4.634569e-2.
  - Every value matches `results/summary.csv` and `logs/assess.log`.
- **Extra falsification search** (my own probe, not part of the study design):
  - Method: a lattice on Dw at h = 0.004 (1701 x 807, 2.5x finer than the assessment lattice), then compass polish
    from 512 distinct starts (separation 0.01) down to a step of 1e-11.
  - It found the same global minimum for **all 5 B finals** and for A 4:1 S1, to 1e-10 relative. In each case the
    lattice minimum was shallower than the polished value.
  - The search found no deeper pocket that the study missed.

Reproduction (scratch paths). Build as in `build.sh` with outputs redirected, then:

    OMP_NUM_THREADS=11 s2.exe run --form A|B --ws 4:1 --start baseline_zero --kmax 150 --maxfev 1500 --out DIR
    s2.exe assess --y DIR/y_final.txt --seed 1|2 --tag T --out assess.csv ; s2.exe axi --y DIR/y_final.txt
    cmp DIR/<file> study2/runs/<run>/<file>   # y_final.txt ys.txt linesearch.csv refine_events.csv set_changes.csv

## 2. Audit of s2.cpp (refine_B, matching, stopping)
**Matched conditions: confirmed.**
- Both formulations use the same `residual`/`gradient` (s = softplus_tau(−ρ), τ from `o0d_settings`).
- `check` shows R(4:1) bit-identical to `feasibility_residual`.
- Starts match `o2_verify.cpp` lines 387–405.
- The L-BFGS and line search are a line-for-line port of `o2_verify.cpp` lines 450–548: trust 0.5, x2 up to 4,
  Armijo 1e-4, 10 halvings, probe at (0.7, 1.3, 0.8), and the 1e-14 curvature guard.
- K_max and the evaluation cap come from the same command line. B’s extra evaluations at set changes count toward
  nfev; the cap never binds (nfev ≤ 176).
- A and B 4:1 S1 trajectories are identical until B’s first split (in-loop q is identical at iterations 2–8).

**refine_B matches PROTOCOL §3.**
- It finds the leaf by containment in (ξ, s = (r/R0)^2).
- ρ_L is recomputed for the current design.
- A leaf is split when ρ* < ρ_L − 0.05|ρ*|, and the containing child is re-tested.
- It guards depth 6 and the 4x node cap.
- `rebuild_nodes` applies the same n x n Gauss rule to the leaf partition, so the integration domain is conserved.
- After a split it recomputes f and g and clears the L-BFGS memory and the f history.
- Event counts agree with REPORT: 9–27 splits, maximum depth 3, cap never reached.
- Only the 8 best polished in-loop minima are candidates. This is as registered, but it means pockets ranked 9th or
  lower are never refined.

**Stopping logic.**
- CONVERGED requires C1 ∧ C2 ∧ C3. C1 and C2 are evaluated on the set on which the iterate was accepted (before
  refinement), and C3 = no split.
- ITERATION_LIMIT and EVAL_LIMIT are tested before each iteration.
- LINESEARCH_FAIL is relabelled only if C1 ∧ C2 hold.
- No defects found.
- Note: C2 was never true at any iterate of any run. C1 was briefly true for A 4:1 S3 (3 iterates) and A 6:1 S3
  (2 iterates). The binding failure is therefore the objective-change criterion, not only the gradient. "No run
  converged" is robust.

**Independence and leakage: no information leak from the final assessment into B.**
- `run` never reads assessment outputs.
- `assess` reads only the y_final and snapshot files.
- File times show every assessment ran after its run finished.

The practical independence of the evaluation from B’s in-loop detector is weak, however (finding I2).

## 3. Deviations D1–D6
None of them changes a conclusion.
- **D1:** confirmed harmless by my bit-identical rerun.
- **D2:** affects only B’s node columns, which do not enter any decision rule.
- **D3:** I applied the protocol-literal STALLED window (all iterates, ignoring set changes) to every trace. It is
  never satisfied, so D3 had no effect.
- **D4:** has no effect on the finals, because ρ_active = ρ*_D. The floor bound q ≥ 0.52 is valid regardless:
  ρ*_D ≤ (true inactive minimum) ≤ the reported floor.
- **D5, D6:** no effect.

Two deviations are **undeclared**; neither changes a conclusion:
- **D7:** the iteration-stability ranges (§5 and the summary column qrange_last20) use the in-loop Dw ρ*, not the
  assessed ρ*_D that PROTOCOL §7 defines q with. At every final the two coincide, so the flags stand.
- **D8:** `mode_assess` truncates the first adaptive level at 4096 cells. The baseline had 4822 qualifying cells, so
  726 were dropped for b. b still reproduces study 1, and every final has ≤ 83 cells.

## 4. Axisymmetry check: meaningful, with limits
- The deviations are small but non-zero: for example, 2.5e-15 relative on the baseline. This shows that ρ really is
  computed from the 3D (y, z) arguments rather than being reduced to r internally. The check therefore exercises the
  real code path.
- I also evaluated ρ at the reported worst points of A 4:1 S1, A 6:1 S1, B 4:1 S1 and the baseline over 720 angles
  φ, and at (x, −r, 0) and (x, 0, r). The maximum deviation was ≤ 1.6e-14.
- Limits:
  - The check confirms that the implementation respects a symmetry the construction guarantees.
  - It cannot detect an asymmetric physics error, and there is no positive control (no non-axisymmetric family is
    available).
  - The analytic argument carries the weight.
- The REPORT states "≤ 2.8e-13". The logs show up to 3.8e-13 (B 4:1 S3); see O1. The check still passes against the
  1e-9 tolerance.

## 5. Findings
### Blocking
**B1. B shows concentration along its trajectories. The pre-registered partial-outcome label is missing, and the
§10 claim contradicts the data.**
- The PROTOCOL §7 concentration criterion is r_E ≥ 0.10 and q ≥ 1.05. It is met at many accepted, non-final B
  iterates. The q values below are in-loop values (from runs/*/trace.csv); these are real point values of ρ, so the
  true q is at least as large.

| Run | Iterate | q | r_E |
|---|---|---|---|
| B 4:1 S1 | 7–8 | 1.42 | 0.235 |
| B 4:1 S1 | 39–41 | 1.152 | 0.559 |
| B 4:1 S2 | 20 | 1.70 | 0.452 |
| B 4:1 S3 | 12 | 1.72 | 0.162 |
| B 6:1 S3 | 15 | 1.83 | 0.412 |
| B 6:1 S1 | 97 | 1.372 | 0.627 |

- Independent confirmation of the late case (B 6:1 S1, iterate 97):
  - I ran the study’s own assessment on the iterate-97 design (row 98 of ys.txt) and got ρ*_D = −1.31830e-3 at
    (−0.796, 2.444), so q = 1.372 and r_E = 0.627.
  - My h = 0.004 dense search gives the same value.
- The study’s own snapshot assessments (assess/snapshots.csv) also show it: B 4:1 S2 at 75/100 gives
  q = 1.086/1.117, and B 4:1 S3 at 100 gives q = 1.079.
- PROTOCOL §7 pre-registers the statement "concentration along B trajectories, not established at convergence" for
  exactly this case. REPORT §6 instead reports only "transient excursions to q ≤ 1.14 inside the last 20 iterates",
  and §5 calls them pockets that are "refined away at the next check".
  - This is not accurate for every excursion. The B 6:1 S1 iterate-97 excursion was not refined: the run has no split
    between iterates 94 and 101. It disappeared through further optimization.
- Required changes:
  - add the pre-registered label to §6;
  - report the full-trajectory maxima;
  - restrict the §1.3 and §10 wording to the final iterates.

### Important
**I1. Scope of the manuscript claim (§10).**
- "The O2 residual is exploitable" generalizes beyond the evidence. The evidence is:
  - family F2A;
  - the Eulerian ρ at t = 0;
  - this L-BFGS port;
  - working densities 4:1 and 6:1;
  - iterate 150 of non-converged runs.
- In addition:
  - q_A falls 6–7x from 4:1 to 6:1.
  - B fails the registered resolution criterion (|Δq_B| = 0.21 and 0.07).
  - B is one adaptive design, and the convergence question failed.
- Suggested wording: "In family F2A (Eulerian ρ, t = 0), after 150 non-converged L-BFGS iterations, minimising the O2
  residual on fixed working nodes deepened the continuum worst point 2.7–23.5x relative to the baseline, and the
  fixed-node checks under-reported it in 4 of 5 designs. With adaptive argmin re-insertion (same integrand and τ),
  E_neg fell by a similar 58–64 %, and the worst point at iterate 150 was not deeper than the baseline
  (q = 0.69–0.98), although intermediate iterates did deepen it by up to 1.8x."
- The final clause is required by B1.

**I2. Independence of the evaluation from B’s detector is weak in practice.**
- In 10/10 finals the assessed ρ*_D equals the final in-loop ρ* to all printed digits. The added lattice, Halton and
  adaptive stages found nothing beyond the in-loop result, and B was refined at exactly those minima.
- The two assessment seeds differ only in the shift of 65,536 Halton points. The 1401 x 701 lattice, the adaptive
  stage and the polish are seed-independent. "Seeds agree to all printed digits" is therefore weak evidence and
  should not be cited as robustness.
- My finer independent search (h = 0.004, 512 starts) confirms the B finals and A 4:1 S1. Pockets narrower than
  about 0.004, or pockets missed by compass polish, remain possible for both formulations. §7 already states this.

**I3. Declare D7 and D8** (section 3) in REPORT §2.

**I4. §1 item 2 and the §6 "A exploitable" row lack a non-converged qualifier at the point of the statement.**
- The Eulerian scope is stated once in the preamble. The non-converged qualifier is stated only in §1.3 and §7.
- Add "at iterate 150 (not converged)" to §1.2 and the §6 row.
- "This is the signature of pockets dug between nodes" is interpretive and should be labelled as an interpretation.
  The resolution trend and the gap between A working nodes and 8:1 nodes support it, but they do not prove it.

### Optional
- **O1.** Axisymmetry maximum: 3.8e-13, not 2.8e-13 (B 4:1 S3, logs/assess.log).
- **O2.** Pre-registration timing.
  - The PROTOCOL.md mtime is 11:43:57. That precedes `s2.cpp` (11:47:12), check.log (11:47:57), the pilot
    (11:48–11:49) and the first run (11:49:56).
  - The protocol therefore precedes all computation, consistent with REPORT.
  - Limits:
    - There is no VCS record, and mtimes are not tamper-evident.
    - The 38 kB driver was saved 3 minutes after the protocol, so it was probably drafted concurrently. The
      pre-registration precedes the runs, not necessarily the code design.
  - Recommend committing PROTOCOL.md (or recording its hash) before future runs.
- **O3.** §1.1 lists the final ‖g‖/‖g0‖ range. Add that C2 was never met at any iterate, and that C1 was briefly met
  by A 4:1 S3 and A 6:1 S3.
- **O4.** §6 "q_A − q_B ≥ 1.7" is true but loose. The minimum is 1.94 (6:1 S3).

## 6. Statements checked in REPORT §1 and §6
- Correct and reproduced:
  - 10/10 ITERATION_LIMIT;
  - q_A 2.70–23.5 and q_B 0.687–0.979 at the finals;
  - r_E for A 0.559–0.659 and for B 0.581–0.643;
  - under-reporting at 4/5 A finals (A 4:1 S3 at 9.6 %);
  - the B resolution and iteration-instability flags;
  - A working-node ρ_min of −7.3e-4 to −1.4e-3.
- The Eulerian/t = 0 restriction holds throughout; no other observer is evaluated anywhere in s2.cpp.
- Statements that need the B1/I4 qualifiers:
  - "B shows no concentration" is true only at the final iterates.
  - "At the iterates reached, reducing E_neg did not require deepening the worst point" is true only as an existence
    statement about the iterate-150 designs.

## 7. Candidate lessons (for the main session; not adopted)
- When a decision rule is evaluated "at final iterates" of non-converged runs, also tabulate the criterion over the
  whole trajectory. Pre-registered partial-outcome labels are easy to omit when they cut against the headline.
  - Evidence: B1.
  - Scope: optimization studies with snapshot or trajectory data.
- An assessment that shares the optimizer-side detector’s polish method, and coincides with it at every final, adds
  little independent evidence. Use a structurally different or finer probe.
  - Evidence: I2.

## 8. Next action
- The builder or the main session revises REPORT §5, §6 and §10 for B1, I1, I3 and I4.
- Gauss’s adequacy review of the refinement and search remains useful, especially for pockets narrower than 0.004.
