# PROTOCOL v3: warp_hypad study 4 (measured Hessian, dense quasi-Newton, nested B, iterate-matched A/B)

Pre-registered by Ada (builder), 2026-10-27. **Version 3**, written after Noether's v2 re-check
(`study4/REVIEW-noether.md`, section "v2 re-check", N1–N5). **No study-4 computation has been run.** No study-4
objective, gradient, Hessian or assessment has been evaluated, and the study-4 driver has not been built or run.
- v1 (`PROTOCOL.md`, SHA-256 `a21e0ff65283875bf98b441652e777d028999bc5bc67317921eb76e7758bf64c`,
  `logs/protocol.sha256`) and v2 (`PROTOCOL-v2.md`, SHA-256
  `b2c0af84b8c9449d97864cf518ad911b5bf51f5c66d4f068dcbf16cb176885da`, `logs/protocol-v2.sha256`) are kept unchanged
  and are **superseded** by this file. Both hashes were re-checked when v3 was written.
- The SHA-256 of this file is recorded in `study4/logs/protocol-v3.sha256` immediately after writing.
- v3 has not been re-reviewed. Noether's re-check says a v3 re-check "can be limited to" N1–N5. The v3 change table
  is section 14. Changes requested by that review produce v4 with a new hash.
- After that review, nothing here changes. Departures are declared in REPORT.md under "Protocol deviations".
- **Marking convention.** Text new in v3 is marked "(v3, N…)". Unmarked text is v2 text, except for the corrections
  listed in section 14.

**Inputs read before writing v2:** everything listed in v1, plus `REVIEW-noether.md` (R1–R11, S1–S7) and
`study3/driver/s3.cpp` lines 326–378 (`refine_B`) and 630–690 (main loop: `capped`, C3, STALLED).
**Inputs read before writing v3:**
- `REVIEW-noether.md`, section "v2 re-check";
- `study3/PROTOCOL.md` section 8 D4 (lines 134–137), for the origin of the 10 % threshold;
- `study3/driver/s3.cpp` lines 610–640 (trust handling);
- `study3/driver/build.sh` (compiler: MinGW `g++` with `-O2 -ffp-contract=off`);
- `calibration/d4_retro_check.txt`.

v3 adds no calibration and reads no new data. Noether computed the N1 baseline numbers in section 8 from archived
study-3 traces. The builder has not recomputed them. The driver recomputes them before any run (gate V6).

**Pre-registration artifacts** (no study-4 quantity evaluated; the v2 additions read archived study-3 files only;
all eight hashes re-checked unchanged when v3 was written):
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
  node order. **(v3, N3)** It also sums them a third time with Neumaier compensated summation, coefficient by
  coefficient. The forward sum is unchanged, so R, F and ∇F remain bit-identical to `s3` (gate V2).
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
- **Compensated reference (v3, N3).** In the same passes, the same stored per-node contributions are also summed with
  **Neumaier compensated summation**, separately for each OTI coefficient (real part, 3 first-order and 6 second-order
  coefficients). This gives H_comp, ∇R_comp and R_comp, all in `double`. No extra passes.
  - Neumaier in `double` is the registered reference because it is portable. Its error is at most about
    2u·Σ|contrib| + O(n u²)·Σ|contrib|, independent of n to first order.
  - **Secondary (descriptive).** If `std::numeric_limits<long double>::digits ≥ 64` in the build (expected for MinGW
    x86-64 `g++`; the value is logged), the forward sum is also accumulated in `long double`. ‖H_fwd − H_ld‖₂ is
    reported next to the Neumaier figure and is not used in any gate.
  - **A-priori bound (reported, free).** For each coefficient, b = (n − 1)·u·Σ_nodes |contrib|, with u = 2⁻⁵³ and n the
    node count. B_H is the matrix of these entrywise bounds on the assembled H_R, and the bound reported is
    ‖B_H‖_F ≥ ‖B_H‖₂. b_g and b_F are reported likewise.
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
- **H1b (node-sum roundoff estimate; measured at every Hessian design, including refreshes and finals).**
  - v2 measures: e_H = ‖H_fwd − H_rev‖₂, e_g = ‖∇F_fwd − ∇F_rev‖, e_F = |F_fwd − F_rev|.
  - **(v3, N3) Compensated measures:** e_H^c = ‖H_fwd − H_comp‖₂, e_g^c = ‖∇F_fwd − ∇F_comp‖ and
    e_F^c = |F_fwd − F_comp|. These measure the error of the value actually used, against a reference whose own
    summation error is second order. e_H is the difference of two estimates with comparable errors and can
    underestimate under cancellation (Noether, v2 re-check).
  - **Registered error measure:** ê_H = max(e_H, e_H^c); likewise ê_g = max(e_g, e_g^c) and ê_F = max(e_F, e_F^c).
    All of e, e^c, ê and the a-priori bound b (section 2.1) are reported at every Hessian design, together with the
    ratio e_H^c/e_H. If e_H^c > 10·e_H at any design, this is noted as "reverse-order estimate understates node-sum
    roundoff".
  - **Comparison of the measures (v3, N3; descriptive).** A table in `hessians/checks.csv` lists e_H, e_H^c,
    ‖H_fwd − H_ld‖₂ (if available) and ‖B_H‖_F per design. Two consistency checks, logged, not gates:
    - e_H^c ≤ ‖B_H‖_F (the a-priori bound must hold);
    - |e_H^c − ‖H_fwd − H_ld‖₂| ≤ 0.1·e_H^c + 1e-2·‖B_H‖_F where `long double` is available. The second term covers
      the two references' own errors: for n ≥ 6,624 nodes these are below about 2/n + 2⁻¹¹ ≈ 8e-4 of the a-priori
      bound, entrywise.

    A violation of either indicates a coding error in the accumulation. It is fixed, and the check is repeated
    (section 2.2, last bullet).
  - **Reference and limitation (v2, strengthened in v3):** every measure here estimates the **summation roundoff of
    the node sum only**. Per-node evaluation roundoff is identical in all three sums and is not measured by any of
    them. Above about 1e-4 relative it is bounded by H2/H3; **below that resolution it is unmeasured**. τ_H and this
    gate are therefore floors, not error bars.
  - **Gate:** if 10·ê_H > 1e-6·max|μ(H_F)| at any design, eigenvalue classification there is precision-limited and the
    study stops for review.
- **H1c (gradient noise floor; reported, gate-free).** ê_g / tol_C1' at every Hessian design. In addition,
  **(v3, N3) ê_g = max(e_g, e_g^c)** is logged at **every accepted iterate** (from the reverse and compensated sums in
  `gradient()`). If ê_g > 0.1·tol_C1' at an iterate, C1' is labelled **"precision-limited"** there. A D1 statement on a
  run whose final row is precision-limited carries that label, and this is known before D1 is read.
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
- e_H, e_g, e_F, **(v3, N3) e_H^c, e_g^c, e_F^c, ê_H, ê_g, ê_F and the a-priori bounds** (H1b), and ê_g/tol_C1';
- counts of negative, near-zero and positive μ. "Negative" means μ < −τ_H, "near-zero" means |μ| ≤ τ_H, with
  **τ_H = max(1e-9·max|μ|, 10·ê_H)** (v2, R6; **v3, N3: ê_H = max(e_H, e_H^c)** replaces e_H). By Weyl's
  inequality, |Δμ_i| ≤ ‖ΔH‖₂, so this is the correct form for the node-sum part of the error. It does not cover
  per-node evaluation error (H1b limitation);
- **(v3, N3)** the number of eigenvalues with τ_H < |μ| < 100·τ_H ("marginal band"), with their signs;
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
  kept (as study 3 kept its pairs). **(v3, N5a) After every split, the trust radius is reset to its start value 0.5**
  (section 4.3). Each refresh costs one Hessian (section 10), is logged with its trigger, and is
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
  - **SPLIT:** n_split ≥ 1 → apply the splits, recompute F and ∇F on the new set, keep W, **(v3, N5a) reset the
    trust radius to 0.5**, restart the C2' history, and begin the next inner solve. Refusals in the same check are
    logged. If this is outer round 10, the run ends **OUTER_LIMIT** instead.
    - **Why the reset (v3, N5a; Noether v2 re-check).** An inner convergence can be reached by a line-search failure
      with C1' ∧ C2' (section 4.1). The trust radius is then below 1e-8. `s3.cpp` (lines 610–640) carries `trust`
      across iterates. Without a reset, the next inner solve would end LINESEARCH_FAIL at once. That would be an
      artefact of the carried radius, not a property of B. The reset applies after every SPLIT, whatever the radius
      was; 0.5 is the registered start value (section 3). The pre-split radius is logged in `outer.csv`.
    - Nothing else is reset at a SPLIT (W kept; refresh counters unchanged; section 3).
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
  - **(v3, N5a)** on SPLIT, `trust = 0.5` is set before the next inner solve;
  - `termination.txt` adds n_outer, the last outcome record and the cumulative refusal counts, and **(v3, N2b)** the
    number of refinement checks n_check and the accepted iterate k_chk1 of the first check (empty if none).
  - **(v3, N2b)** `outer.csv` has one row per check: the accepted iterate at which it was made, its outcome (PASS,
    SPLIT, BLOCKED_DEPTH, BLOCKED_CAP), the outcome record, and the trust radius before any reset.
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
  - **(v3, N2b)** the number of refinement checks n_check (PASS, SPLIT or BLOCKED), n_check within k_c, the iterate
    k_chk1 of the first check, and the number of complete D5w windows that end after k_chk1, at W = 20 and W = 25;
  - D5w carries the labels in section 8. **(v3, N2b)** These are based on refinement **checks** and their position
    relative to the windows, not on splits: a check that returns PASS has tested B's set, and is a legitimate B fixed
    point;
  - **(v3, N2b)** if B made no refinement check, B's D3/D1 statements read "failed (not converged; no study-4
    refinement check)".
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
  - **(v3, N1; corrects v1/v2) A 6:1 starts on its onset spike.** v1 and v2 stated that the 6:1 A start (iterate 90)
    precedes A's study-3 onset (~105–110, the value in Noether's v1 review), "so the windows cover the onset". That
    statement is incomplete. In the archived in-loop q trace, A 6:1 rises from 0.71 at iterate 87 to 2.02 at
    iterate 90 (Noether v2 re-check, N1). The start iterate therefore sits on A's in-loop q onset spike.
    - The ~105–110 onset was not re-derived for v3, and its defining quantity is not restated here. The study-4
      windows still begin after iterate 90 and may cover later onset behaviour.
    - Iterate 0 of `A6L2m` is not a pre-onset baseline. The 6:1 point gap (+1.15) includes the spike, while the
      pre-window gap (−0.09) lies almost entirely before it. Section 8 (gap companion) uses both baselines for this
      reason.
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
  - **(v3, N3) Marginal strict.** If every μ(H_F) > τ_H but μ_min(H_F) < 100·τ_H, the label is **"strict local
    minimiser (second-order) (marginal)"**, short form **"strict (marginal)"**. μ_min/τ_H is reported for every
    converged final.
  - **(v3, N3) Limitation, stated with every D1b label:** "τ_H covers node-sum roundoff only; per-node evaluation
    error is unmeasured below the H2/H3 resolution (about 1e-4 relative)". A "strict" label is established to
    node-sum roundoff, not to all roundoff. A "strict (marginal)" label is one where unmeasured per-node error could
    plausibly change the classification.
  - **(v3, N3) Downstream use.** Where D3, D6x or S2 require "D1b strict", "strict (marginal)" satisfies the
    condition. The resulting label then carries the suffix **"(marginal)"**: for example, "distinct local minimisers
    (marginal)" or "multiple minima at λ_k (marginal)". It carries this suffix if any final that determines it is
    marginal. This is a builder choice, declared for review (section 14): the alternative was to demote marginal
    finals to "distinct stationary points" or "path-dependent (trajectory)", which would hide the measured spectrum.
- **D3 Concentration at convergence** (per formulation, at λ2 on the P1 runs): as in study 3. It needs both densities
  converged; otherwise "failed (not converged)", with the terminal states and, for B, n_outer and (v3, N2b) n_check.
  - **P1-A vs P2-A (v2, R4)**, per density:
    - both CONVERGED and the same point → "same point";
    - both CONVERGED, not the same point, both D1b strict → **"distinct local minimisers"** (v3, N3: with
      "(marginal)" if either is "strict (marginal)");
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
    have flagged both A λ2 pockets and nothing else. The check covers non-converged study-3 finals only.
  - **(v3, N4) Post-hoc status: registered wording, reported with every D4 result.** "D4 statistic chosen post hoc
    (study-3 REPORT §6, Noether R5); retro-check on the motivating data is a consistency check, not a validation;
    P2 starts on the calibration designs (start u_w reported)."
    - **The statistic was chosen after study-3 data were seen.** Noether proposed u_w (R5) after study 3 showed that
      the A λ2 finals hid the pocket off the working nodes. The retro-check shows only that the rule flags the cases
      it was chosen to flag. It is not evidence of the rule's behaviour on new designs.
    - **P2 starts on the calibration designs.** `A4L2f` and `A6L2f` start at `A_4_L2/y_final` and `A_6_L2/y_final`,
      the two designs used in the retro-check (u_w 98.2 % and 87.9 %). A "fixed-node exploitation" label on a P2 final
      is therefore close to inherited unless the run moves far from its start.
    - **Which runs start on calibration designs.** Every run except `A4L2m` and `A6L2m` starts at one of the 8
      retro-check designs: P2, P3, P4, P5, and the P1 B runs (`B_4_L2/y_final`, `B_6_L2/y_final`). Every D4 result
      (registered or descriptive) from such a run carries **"(start on calibration design)"**. It is reported with
      ‖y_final − y_start‖/‖y_start‖ and q_final − q_start. No threshold for "moved far" is registered, so these results
      are not called out-of-sample.
    - **Out-of-sample applications** of the registered rule are therefore the P1 A finals (`A4L2m`, `A6L2m`) only.
    - **u_w at the start is reported next to u_w at the final for every run** (starts are assessed; section 7, item 2),
      including B and non-converged runs (descriptive).
  - **(v3, N4; corrects v2) Origin of the 10 % threshold.** v2 said "the threshold (10 %) is the study-3 value". More
    precisely, the number 10 % comes from the study-3 D4 under-reporting branch (`study3/PROTOCOL.md` line 135: the
    worst of the working-node, 8:1-node and `refine_around_minimum` checks "under-reports ρ*_D by more than 10 %"). The
    same number is used in study-3 PROTOCOL line 170 and in study-4 v1 §4.3 for "B inner exploitation". In study 3 the
    10 % was applied to a **different statistic** (the worst of three node sets, not the working set alone). It was not
    chosen for u_w. **Threshold sensitivity:** on the archive, the separation is 0.3–4.6 % against 87.9–98.2 %, so any
    threshold in about (5 %, 85 %) gives the same classification there (Noether, v2 re-check). The threshold was not
    tuned in a way that matters for the archive. That says nothing about study-4 designs.
- **D5w Matched iterate budget, plateau phase (replaces D5; Gauss G-I3 (b); v2 R2, R7, R8, S4, S6; v3 N1, N2a, N2b).** For each λ and
  density, on the A/B pair of P1 (P4 at λ1):
  - **Terminal handling and horizon (v2, R2).**
    - A run that ends CONVERGED (or CONVERGED_LS_LIMITED) at k_f < K is **extended as stationary**: its in-loop
      q(k) = q(k_f) for k_f < k ≤ K.
    - Every other terminal state (ITERATION_LIMIT, LINESEARCH_FAIL, INVALID, EVAL_LIMIT, OUTER_LIMIT,
      REFINE_BLOCKED_*, WALL_LIMIT) **truncates** at its last accepted iterate. STALLED is not terminal (section 4.1).
    - **k_c = min(K, the terminal iterate of each non-converged run in the pair).** If both converged, k_c = K.
    - The number of windows that contain extended (stationary) values is reported.
    - **(v3, N2a) Padding label.** For each W, n_ext(W) is the number of verdict windows in which either run has at
      least one extended value, and m(W) is the number of verdict windows. If n_ext(W) ≥ m(W)/2 at W = 20 **or**
      W = 25, the D5w verdict (and the gap companion, which uses the same Δ_w) carries **"(determined by converged
      final(s))"**. It is then reported as a comparison of the converged final(s) with the other run's trajectory.
      The ≥ 4-window condition and the W = 20/25 agreement are **not cited as robustness**, because stationary windows
      agree mechanically. n_ext and m are reported at both W.
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
  - **Refinement labels (v3, N2b; replace the v2 R7 label "0 outer rounds within k_c").** They are based on B's
    refinement **checks** (PASS, SPLIT or BLOCKED; section 4.3), not on splits. k_chk1 is the accepted iterate at which
    B's first study-4 check is made.
    - **"(no study-4 refinement check)"**: B made no refinement check within k_c (k_chk1 is undefined or > k_c).
      Every window then compares A with a fixed-node run on B's archived refined set.
    - **"(refinement-limited: n of m windows after the first check)"**: B made at least one check within k_c, but
      fewer than 4 complete verdict windows end after k_chk1 (window end > k_chk1), at either W. n and m are given
      for the W that sets the label (the smaller n).
    - **Exception:** a B run that is CONVERGED (via PASS) receives neither label. Its set has been tested by the
      refinement rule, and its final is a legitimate B fixed point. The padding label (N2a) may still apply.
    - A verdict with either refinement label is **not cited as a formulation effect**, and "suppresses" is not used
      for it.
  - **Gap companion (v2, S4; rule replaced in v3, N1).** It asks whether the inherited A−B gap **changed** during
    study 4.
    - **Two baselines, both registered (v3, N1).** For each density and W:
      - **pre-window gap** Δ_pre(W) = median in-loop q_A − median in-loop q_B over the archived study-3 iterates
        [s0 − W + 1, s0] of the two start trajectories (s0 = 84 at 4:1, 90 at 6:1). **These iterates are study-3
        L-BFGS (m = 8) iterates, and B's were made under B's study-3 per-iterate refinement**, not under study-4 B;
      - **start (point) gap** Δ_0 = in-loop q_A(s0) − in-loop q_B(s0), from the same archived traces.
    - **Registered values** (Noether v2 re-check, from `study3/runs/{A,B}_{4,6}_L2/trace.csv`, q = rhostar_in/b;
      recomputed by gate V6 before any run):

      | Density | s0 | Δ_0 | Δ_pre (W = 20) | Δ_pre (W = 25) | max(Δ_pre, Δ_0) | min(Δ_pre, Δ_0) |
      |---|---|---|---|---|---|---|
      | 4:1 | 84 | +1.617 | +2.287 | +2.273 | Δ_pre | Δ_0 |
      | 6:1 | 90 | +1.146 | −0.093 | −0.084 | Δ_0 | Δ_pre |

      At 6:1 the two baselines differ in sign. The start iterate is A's onset spike (section 5): A's in-loop q rises
      from 0.71 at iterate 87 to 2.02 at iterate 90. Δ_pre measures the gap before the spike, and Δ_0 measures it on
      the spike.
    - **Why both (v3, N1).** v2 used Δ_pre alone. At 6:1, Δ_pre ≈ −0.09 makes G_w = Δ_w + 0.09, so every window that
      supports "B lower" (Δ_w ≥ +0.05) has G_w ≥ +0.14, and "B lower" implied "gap grows". The word "suppresses" would
      then have been licensed by the D5w verdict alone, even if the inherited +1.15 gap had shrunk to +0.2. At 4:1 the
      bias runs the other way (Δ_pre > Δ_0): conservative for "grows", liberal for "shrinks". The v2 λ1 calibration
      could not show this, because it was a null case that tested only "shrinks" false positives.
    - **Series.** G_w^hi = Δ_w − max(Δ_pre, Δ_0) and G_w^lo = Δ_w − min(Δ_pre, Δ_0). Since G_w^hi ≤ G_w^lo, "grows" is
      tested on the conservative series for growth, and "shrinks" on the conservative series for shrinking.
    - **Rule for one W** (same windows and thresholds as the D5w verdict):
      - **"gap grows"**: G_w^hi ≥ +0.05 in ≥ 2/3 of the windows and no window has G_w^hi ≤ −0.05;
      - **"gap shrinks"**: G_w^lo ≤ −0.05 in ≥ 2/3 of the windows and no window has G_w^lo ≥ +0.05;
      - **"not established (baseline-sensitive)"**: neither of the above holds, but the v2 rule applied with one
        baseline alone (G_w = Δ_w − Δ_pre, or G_w = Δ_w − Δ_0) gives "grows" or "shrinks";
      - **"not established"**: otherwise;
      - **"not evaluable (< 4 windows)"**.

      "Grows" and "shrinks" are mutually exclusive: with ≥ 4 windows, a window with G_w^hi ≥ +0.05 also has
      G_w^lo ≥ +0.05.
    - **Across W.** If either W is "not evaluable", the result is "not evaluable". If the two W give the same label,
      that label holds. Otherwise the result is "not established (window-sensitive)".
    - **Reported in all cases:** the single-baseline results (Δ_pre alone and Δ_0 alone) at both W, all G_w^hi and
      G_w^lo, and the study-4 iterate-0 in-loop Δ next to the archived Δ_0. A difference > 0.01 between those two is
      reported as a reconstruction deviation (for B, this tests the reconstructed set).
    - **Consequences of the registered numbers** (arithmetic on the table above, not a prediction):
      - 6:1 "grows" needs Δ_w ≥ +1.196 in ≥ 2/3 of the windows. 6:1 "shrinks" needs Δ_w ≤ −0.143 (W = 20) or
        −0.134 (W = 25): A lower than B, which is a reversal.
      - 4:1 "grows" needs Δ_w ≥ +2.337 (W = 20) or +2.323 (W = 25). 4:1 "shrinks" needs Δ_w ≤ +1.567.
      - The two-baseline rule is therefore conservative for "shrinks" at 6:1. A narrowing that stops short of
        reversal reads "not established (baseline-sensitive)" or "not established" there.
    - The v2 calibration result (single-iterate baseline: 2 false "gap shrinks" in 14 λ1 alignments; pre-window
      baseline: 0) is unchanged. The "shrinks" test now uses the smaller baseline, which is the pre-window one at 6:1
      and the point one at 4:1. A single-iterate baseline produced the calibration's false positives (λ1, 4:1,
      s0 = 40 and 60). Requiring "shrinks" against both baselines means that a point baseline alone can no longer
      produce it.
  - **What each outcome licenses (v3, N1).** Per density. Any refinement label (N2b), padding label (N2a),
    "(not converged)" or "(marginal)" applies as stated above, and every claim keeps the confounds of section 5.

    | D5w verdict | Gap companion | Licensed statement |
    |---|---|---|
    | B lower, no refinement label | gap grows | **"B suppresses deepening at matched cumulative iterates"**. The A−B in-loop q gap grew beyond both the start-point gap and the pre-start-window gap. (Continuation from inherited, unequal starts; not a common-start test.) |
    | B lower, no refinement label | not established (baseline-sensitive) | "B lower at matched cumulative iterates (start gap inherited)"; "change in gap not established (baseline-sensitive)". **No "suppresses".** At 6:1, growth against Δ_pre only means the gap stayed above its pre-onset level, which is what persistence of the onset-spike gap would also give. At 4:1, growth against Δ_0 only means the gap exceeds the start point but not the pre-start window. |
    | B lower, no refinement label | gap shrinks | "B lower at matched cumulative iterates (start gap inherited); gap narrowed against both baselines". No "suppresses". The inherited difference persists but decreased. |
    | B lower, no refinement label | not established / window-sensitive / not evaluable | "B lower at matched cumulative iterates (start gap inherited)": persistence of the inherited difference; change not detected. No "suppresses". |
    | B lower, with a refinement label | any | "B lower at matched cumulative iterates (start gap inherited)" plus the label. Not a formulation effect. No "suppresses". The gap companion is reported descriptively. |
    | A lower | any | "A lower at matched cumulative iterates (start gap inherited)". The gap companion is reported. No B-suppression statement. |
    | not established / not evaluable | any | No level statement. The gap companion is reported descriptively only. |

  - **"Suppresses", P5 route (v2, R8; unchanged; 4:1 only).** "B suppresses deepening" is also licensed at 4:1 if the
    D5w verdict is "B lower" without a refinement label **and** P5 agrees: the `XB4L2fromA` final is shallow while the
    `A4L2f` final has a pocket. The suffix "(not converged)" is added if either run is not converged. This is a
    common-start statement and is worded as such: "(common start, P5)". Otherwise the table above applies.
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
  - Labels, in precedence order (first that applies; mutually exclusive by construction). **(v3, N5c; recommended by
    Noether, adopted)** "inconclusive" is now checked first. In v2 it came last, so one missing arm plus one "keeps"
    arm gave "path-dependent", which rested on a single arm. That was not intended.
    0. **"inconclusive"**: an arm has no assessed final (not run, or INVALID before iterate 1). The available arm's
       character is reported descriptively.
    1. **"multiple minima at λ_k"**: an arm keeps its start character, is CONVERGED and D1b strict, and the reference
       continuation at the same λ and density (λ1: `A4L1f` from P4; λ2: `A4L2f` from P2) is CONVERGED, D1b strict,
       of the opposite character, and not the same point. Both minima are named. (v3, N3: "(marginal)" if either is
       "strict (marginal)".)
    2. **"path-dependent (trajectory)"**: at least one arm keeps its start character, and 1 does not hold.
    3. **"budget-dependent"**: both arms adopt the character of their λ.
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
| V5d | **Refinement outcome logic (R1).** `refine_B` driven by a synthetic ρ field on the 4:1 base set: a resolved candidate gives PASS; a candidate needing a split at depth 6 gives BLOCKED-DEPTH with C3 false; a set at the node cap gives BLOCKED-CAP with C3 false; a needed split below both limits gives SPLIT; a check after a cap refusal is still executed (no persistent flag). **(v3, N5b; recommended, adopted)** Also: a mixed check (one split plus one depth refusal) gives SPLIT, and the refused row is re-tested at the next check; a SPLIT at outer round 10 gives OUTER_LIMIT. **(v3, N5a)** After a SPLIT made with a trust radius below 1e-8, the next inner solve starts with trust = 0.5. |
| V5e | **(v3, N3) Compensated summation.** (i) Unit cases with a known exact sum: [1, 1e100, 1, −1e100] gives 2; a seeded sequence of 10⁵ terms with heavy cancellation, summed exactly with Python `fractions.Fraction`, matches the Neumaier result to ≤ 2u·Σ\|x\|. (ii) In V5c, H_comp, ∇_comp and F_comp match the analytic values to the V5c tolerance, and e_H^c ≤ ‖B_H‖_F. (iii) The build's `std::numeric_limits<long double>::digits` is logged. |
| V6 | **(v3, N1) Gap baselines from archived traces.** Recompute Δ_0 and Δ_pre(W) for 4:1 and 6:1 at W = 20 and 25 from `study3/runs/{A,B}_{4,6}_L2/trace.csv` (q = rhostar_in/b, b as in section 1). They must equal the section-8 table to ±0.001 (the rounding of the quoted values). Also confirm A 6:1 in-loop q at iterates 87 and 90 (0.71 and 2.02, ±0.005). A mismatch stops the study for review, because the registered thresholds depend on these numbers. This reads archived study-3 files only. |

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
  before step 0. (v3) V5e and V6 add an estimated 0.3 h of builder time.
- **(v3, N3)** The compensated and `long double` sums reuse the stored per-node contributions and need no extra
  passes. Their cost relative to a pass is expected to be small: one summation of n_nodes × 10 coefficients per pass,
  against n_nodes OTI evaluations of ρ. This is unmeasured and is checked at V4.
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

## 13. Changes from v1 to v2 (Noether review, REVIEW-noether.md; kept as written in v2)
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

## 14. Changes from v2 to v3 (Noether v2 re-check, REVIEW-noether.md, "v2 re-check")
No new data was read, no calibration was added and no study-4 quantity was evaluated. The v2 re-check verdicts R1,
R3, R4 and R8–R11 ("resolved") are carried over without change, except where a row below names a section.

### Required changes
| ID | Finding (Noether) | Resolution in v3 | Where |
|---|---|---|---|
| N1 (high) | The pre-window gap baseline voids the R8 safeguard at 6:1. Δ_pre ≈ −0.09 there, while the point gap is +1.15, so "B lower" implied "gap grows". The A 6:1 start iterate is its onset spike (0.71 at iterate 87 → 2.02 at iterate 90). At 4:1 the bias is reversed. | Gap companion tests "grows" against max(Δ_pre, Δ_0) and "shrinks" against min(Δ_pre, Δ_0). A result against one baseline only is "not established (baseline-sensitive)". Both baselines registered with values at both W. Single-baseline results reported. The pre-window is stated to lie in study-3 L-BFGS iterates under B's per-iterate refinement. A table states exactly which statement each D5w × gap outcome licenses; "suppresses" needs "gap grows" against both baselines (or the 4:1 P5 route). The v1/v2 sentence "the 6:1 A start precedes the onset" is corrected: iterate 0 at 6:1 is on A's onset spike. Arithmetic consequences of the registered numbers stated (6:1 "shrinks" requires reversal). Gate V6 recomputes the baselines from archived traces before any run. | §5 Matched budgets; §8 D5w gap companion and licence table; §9 V6 |
| N2a (medium) | Padding converged runs as stationary manufactures robustness windows. | n_ext(W) and m(W) reported. If n_ext ≥ m/2 at either W, the D5w verdict and gap companion carry "(determined by converged final(s))". The ≥ 4-window and W-agreement conditions are then not cited as robustness. | §8 D5w terminal handling |
| N2b (medium) | The R7 label counted splits, not refinement checks, and ignored where checks fall in the windows. | Labels based on refinement checks (PASS, SPLIT or BLOCKED): "(no study-4 refinement check)" if there is no check within k_c; "(refinement-limited: n of m windows after the first check)" if fewer than 4 complete windows end after k_chk1 at either W. A B run CONVERGED via PASS is not labelled. Either label blocks formulation-effect citation and "suppresses". n_check, k_chk1 and the per-check `outer.csv` rows logged. D3/D1 wording for B updated. | §4.3; §8 D3, D5w |
| N3 (medium) | Forward-vs-reverse is the difference of two comparably erroneous estimates, and can underestimate under cancellation. It does not measure per-node evaluation error. | Neumaier compensated sum of the same stored contributions in the same passes gives e_H^c, e_g^c, e_F^c. ê = max(e, e^c) is used in τ_H, the H1b gate and H1c. The a-priori bound (n−1)u·Σ\|contrib\| is reported. A `long double` sum is reported descriptively where the build has a ≥ 64-bit mantissa. Two logged consistency checks compare the measures. Limitation wording is added to H1b and to every D1b label. "Strict" finals with μ_min < 100·τ_H are labelled "strict (marginal)", and the suffix propagates to D3/D6x labels. Gate V5e tests the summation. The forward sum is unchanged, so V2 bit-identity is unaffected. | §1; §2.1; §2.2 H1b, H1c; §2.3; §8 D1b, D3, D6x; §9 V5e; §10 |
| N4 (low-medium) | The D4 statistic was chosen after seeing study-3 data. P2 starts on the calibration designs. v2's origin statement for the 10 % was imprecise. | Registered wording reported with every D4 result: "D4 statistic chosen post hoc (study-3 REPORT §6, Noether R5); retro-check on the motivating data is a consistency check, not a validation; P2 starts on the calibration designs (start u_w reported)". u_w at the start is reported next to u_w at the final for every run. Every run starting on a retro-check design (all except `A4L2m`, `A6L2m`) carries "(start on calibration design)". Only the P1 A finals are called out-of-sample. Origin of the 10 % corrected: study-3 D4 under-reporting branch (`study3/PROTOCOL.md` line 135), applied there to a different statistic. Threshold insensitivity on the archive is stated, with the limit that it says nothing about study-4 designs. | §8 D4 |
| N5a (low) | After a split that followed a CONVERGED_LS_LIMITED-type check, the carried trust radius (< 1e-8) would end the next inner solve at once (`s3.cpp` 610–640 carries `trust`). | Trust reset to 0.5 after every SPLIT. The pre-split radius is logged. The reset is tested in V5d. | §3; §4.3; §9 V5d |

### Recommended changes (adopted)
| ID | Resolution | Where |
|---|---|---|
| N5b | V5d adds the mixed case (split plus depth refusal → SPLIT; the refusal is re-tested at the next check) and the OUTER_LIMIT case. | §9 V5d |
| N5c | D6x checks "inconclusive" first, so that no label rests on a single arm. | §8 D6x |

### Builder choices declared for the v3 re-check
- **"(marginal)" propagates; it does not demote.** A "strict (marginal)" final satisfies "D1b strict" in D3/D6x/S2,
  and the derived label carries "(marginal)". The alternative (demotion to "distinct stationary points" or
  "path-dependent (trajectory)") was not chosen, because it would discard the measured spectrum. The reviewer may
  reverse this.
- **Neumaier in `double` is the registered reference**, not `long double`, for portability. The MinGW `long double`
  sum is descriptive only.
- **"Moved far" is not defined.** Instead, every run starting on a calibration design is labelled as such, and only
  the P1 A finals are treated as out-of-sample D4 applications.
- **The N1 baseline values are Noether's, not recomputed by the builder.** Gate V6 recomputes them before any run, and
  a mismatch stops the study.
- **D5w windows and thresholds are unchanged** (0.05, 2/3, ≥ 4 windows, W = 20/25). Only the gap baselines and the
  labels changed.

### Not changed
- D5w level rule, K, the outer cap, node-cap and depth limits, the line search and the BFGS update. Trust reset (N5a)
  is the only solver change.
- No new calibration was run. The v2 calibrations are unchanged, and their hashes were re-checked.
