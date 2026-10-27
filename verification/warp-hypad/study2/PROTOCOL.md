# PROTOCOL: warp_hypad study 2 (optimization to a stopping criterion, matched A/B)

Pre-registered by Ada (builder), 2026-10-27, 11:45 CDT. Written before any study-2 computation.
Nothing in this file is changed after runs start; departures go in REPORT.md under "Protocol deviations".
Repository: the local `Enterprise` clone at 477141d, unmodified. Libraries are the existing `build-oti/` ones.
All statements concern the **Eulerian energy density ρ = T(n,n)** (ADM normal observer, TotalEffective
source, Λ = 0, t = 0). No other observer or energy condition is evaluated.

## 0. Question
Study 1 found that minimising the O2 residual R on fixed quadrature nodes lowers E_neg while the
continuum worst point of ρ deepens, but the runs were not converged and the deepening was iteration- and
sampling-dependent. Study 2 asks:
1. Does the deepening persist **at convergence** when the optimizer **cannot exploit node gaps** (an
   objective-design property), or
2. does it disappear when gaps are closed (the fixed-node formulation is exploitable)?

## 1. Fixed problem (identical for every run)
- Family F2A (lapse, axial shift, radial shift), 108 coefficients (36 active per field), the F1W
  support-excluded basis, config v = 0.6, R = 2, σ = 1, R0 = 2.
- τ = 1e-3. Integrand s(ρ)² with s = softplus_τ(−ρ), as in `feasibility_residual`.
- Optimizer: the study-1 port of the O2 L-BFGS (memory 8, Armijo c1 = 1e-4, backtracking ×0.5,
  ≤ 10 trials per trust level, initial trust radius 0.5, ×2 expansion after a first-trial acceptance up to
  4, ×0.5 shrink after a failed line search, same validity probe). Only the stopping rules (§2) change.
- Starts:
  - S1 `baseline_zero`
  - S3 `combined_perturbation`
  - S2 `lapse_perturbation`
  
  All three are defined exactly as in O2.
- Working (optimizer) domain Dw: ξ ∈ [−1.7, 1.7], s ∈ [0, 2.6], i.e. x ∈ [−3.4, 3.4], r ∈ [0, 3.225].
- Assessment domain D: the O2 validation domain ξ ∈ [−3.5, 3.5], s ∈ [0, 12.25], i.e.
  x ∈ [−7, 7], r ∈ [0, 7], with r the cylindrical radius.

## 2. Convergence criteria (same for A and B)
For accepted iterate k, let g_k be the gradient on the current working set and g_0 the gradient at the start
on the initial working set.
- **C1 (gradient):** ‖g_k‖₂ ≤ 1e-3 · ‖g_0‖₂.
- **C2 (objective change):** (f_{k−5} − f_k) / f_k ≤ 1e-6, using only iterates since the last
  working-set change.
- **C3 (B only, set stability):** the in-loop refinement check at iterate k made no change.
- **CONVERGED** iff C1 ∧ C2 (∧ C3 for B). This is the only pass state.
- **Failure states** (none count as converged):
  - **ITERATION_LIMIT:** K_max = 150 accepted iterations reached.
  - **EVAL_LIMIT:** 1500 objective evaluations reached.
  - **STALLED:** (f_{k−20} − f_k) / f_k ≤ 1e-6 with C1 false. The run stops.
  - **LINESEARCH_FAIL:** trust radius < 1e-8 after repeated failed line searches. If C1 ∧ C2 hold at
    that point, the state is "converged (line-search limited)" and reported as such.
  - **INVALID:** a non-finite objective or gradient at an accepted point.
- Failed runs are archived and reported. Their final iterates are assessed but labelled non-converged.

## 3. Formulations (the matched comparison)
**A, fixed-node (current O2).**
- Objective: R_A(y) = Σ_p w_p s(ρ(x_p; y))² on the fixed elementwise Gauss set at working density.
- Fixed-node check (the O2 certification path):
  - node ρ_min on the 8:1 validation set,
  - `refine_around_minimum` (half-width 0.35, 24×24),
  - node ρ_min on the working set.

**B, adaptive argmin re-insertion.**
- Same integrand, same τ, same base set.
- The quadrature is h-refined so that every detected pocket of ρ < 0 in Dw is resolved by nodes; the
  elements are the leaf rectangles in (ξ, s), each carrying the same n×n Gauss rule.
- Before the first iteration and after every accepted iterate, run the **in-loop search** on Dw:
  1. a uniform 241 × 116 lattice in (x, r);
  2. compass-search polish of the 8 best distinct lattice minima (separation ≥ 0.1) down to a step of 1e-9.
- For every polished minimum x* with ρ* < 0:
  - locate the leaf element L containing it, and let ρ_L be the minimum of ρ over L's nodes;
  - if ρ* < ρ_L − 0.05·|ρ*|, split L 2×2 and re-test the child containing x*;
  - the maximum depth is 6 below the base element.
- Splitting conserves the integration domain, so R_B remains a quadrature of the same integral ∫_Dw s² dV.
- After any change, recompute f and g on the new set and clear the L-BFGS memory. Set changes are not
  counted as iterations.
- Node cap: 4× the base node count. If the cap is reached, refinement stops and the run is flagged.
- B never coarsens.

**Matched conditions.** A and B use the same starts, basis, τ, line search, K_max, evaluation cap,
stopping rules, and base working set. For both, the in-loop search is also run and logged as a diagnostic
after every accepted iterate; in A it does not affect the optimizer.

## 4. Resolution study
Two working densities (Gauss order : subdivision):
- **4:1** (6,624 base nodes, the O2 setting)
- **6:1** (14,904 base nodes)

The extent is unchanged.

**Run priority** (pre-registered, executed in this order; A and B of each pair run concurrently):

| Priority | Density | Start |
|---|---|---|
| P1 | 4:1 | S1 |
| P2 | 6:1 | S1 |
| P3 | 4:1 | S3 |
| P4 | 6:1 | S3 |
| P5 | 4:1 | S2 |

**Budget rule:**
- If a 3-iteration timing pilot of A at each density projects P1–P4 above 5 h of wall time, drop P5 first,
  then lower K_max to 100 for P3–P4. Any such change is recorded as a deviation.
- Stop all work and report if the cumulative wall time passes 6 h.

## 5. Independent worst-point assessment (applied identically to every final design, the baseline, and trajectory snapshots)
It is independent of every optimizer node set and of the B in-loop lattice.

**Axisymmetry.**
- Analytic argument: F2A fields are functions of (ξ, s) with s = (y² + z²)/R0², and the radial shift is
  (y, z)/R0 · Q. With γ = I, the metric is invariant under rotations about the x-axis, so the scalar ρ is
  axisymmetric.
- Numerical check: ρ is evaluated at 2,000 random 3D points (x, r cos φ, r sin φ) and compared with
  ρ(x, r, 0). Agreement must be within 1e-9·max|ρ| + 1e-14. If the check fails, the search is repeated in
  3D and the deviation is reported.
- Given axisymmetry, the search is run in the half-plane (x, r).

**Search on D.**
1. **Space-filling:**
   - a 2^16-point Halton design (bases 2, 3), with a fixed random shift;
   - a uniform 1401 × 701 lattice (Δ = 0.01).
2. **Adaptive subdivision:**
   - every lattice cell whose corner minimum lies within 10 % of the current global minimum is resampled
     on an 8 × 8 sub-lattice;
   - the 32 best sub-cells are subdivided again, 3 levels.
3. **Multi-start local polish:**
   - compass search from the 64 best distinct candidates (separation ≥ 0.05) down to a step of 1e-10;
   - the Halton shift is repeated with a second seed, and the agreement of the two global minima is
     reported.

**Reported quantities:**
- ρ*_D and its location;
- ρ*_Dw;
- ρ* on the design-active region and the design-independent floor ρ_ext, defined as the minimum of
  baseline ρ over the part of D where the design cannot change ρ, measured as |ρ(y_rand) − ρ(0)| = 0
  exactly;
- E_neg (8:1 validation and 16:1 reference);
- negative-volume fraction.

**Limitations (stated in advance):**
- This is a sampling and local-search procedure.
- A denser grid or more starts is not proof of the global minimum.
- No certified (Lipschitz or interval) lower bound is implemented, so ρ*_D is an upper bound on the true
  infimum over D.
- The meridian reduction rests on the axisymmetry argument and check above.

## 6. Outputs
Per run, under `study2/runs/<form>_<density>_<start>/`:
- the configuration;
- per-iterate trace:
  - f, ‖g‖, step, trust radius, evaluation counts;
  - node count, refinement events;
  - E_neg on the 8:1 validation set;
  - node ρ_min on the working and validation sets, plus the `refine_around_minimum` value;
  - in-loop ρ* and its location;
- every accepted y;
- the termination state;
- the timing;
- the log.

Failed runs stay in place.

**Figures:**
- ρ on the meridian (x, r) for the baseline, A and B (matched start and density), with worst points and
  pockets marked, using symmetric-log colour;
- iteration dependence: E_neg, in-loop ρ*, node ρ_min against iteration, for A vs B;
- resolution dependence: final q and r_E against density, for A vs B.

## 7. Decision rules
Let b = ρ*_D of the baseline, E0 = baseline E_neg (8:1), q = ρ*_D / b (q > 1 means deeper than the
baseline), and r_E = 1 − E_neg / E0.
- **Concentration** at a design: r_E ≥ 0.10 and q ≥ 1.05.
  - **Marginal:** 1.01 ≤ q < 1.05.
  - **None:** q < 1.01.
- **Fixed-node under-reporting** (per design): the worst of the fixed-node check values is shallower than
  ρ*_D by more than 10 % of |ρ*_D|.
- **Formulation A exploitable** if both hold:
  - (i) fixed-node under-reporting occurs at an A final design, or along an A trajectory at the iterate
    where A's working node ρ_min is reported; and
  - (ii) at a matched start and density, q_A − q_B ≥ 0.05 at the final iterates, or B shows no
    concentration where A does.
- **Concentration persists at convergence (objective-design result)** if all hold:
  - at least one start has B CONVERGED at both densities;
  - every converged B run shows concentration;
  - |q(4:1) − q(6:1)| < 0.05 for that start.
- **Partial outcomes:**
  - B shows concentration only on non-converged trajectories: "concentration along B trajectories, not
    established at convergence".
  - No B run converges: the convergence question is "failed (not converged)", and only trajectory
    statements are made.
- Both outcomes may hold at once: A can be exploitable while B still concentrates. Both would be reported.
- Iteration dependence is reported as the range of q and r_E over the last 20 accepted iterates. A range
  in q > 0.05 is flagged as "not iteration-stable".

## 8. Status labels
These are used in REPORT.md:
- **reproduced:** a study-1 number recovered;
- **new:** a study-2 result;
- **failed:** a check or criterion not met.

Builder results are not independently reviewed.
