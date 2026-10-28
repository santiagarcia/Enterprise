# Study 4 REPORT: execution under PROTOCOL-v3, stopped at gate V4 (H3) for review

Builder: Ada, 2026-10-27 (21:15–23:10 CDT, about 1.9 h wall time). These are builder results and have **not been
independently reviewed**. Saved by the main session from Ada's hand-back (the subagent's write was blocked).
Scope: Eulerian ρ only; F2A, start-S1 lineage, λ1 = 3.75e-8 and λ2 = 3.75e-9 (as in PROTOCOL-v3).

## 1. Outcome
- `PROTOCOL-v3.md` SHA-256 `768b6036…e70a3` matched `logs/protocol-v3.sha256` before work began. All eight
  pre-registration artifact hashes also matched.
- Gates V1, V2, V3, V5a–V5e and V6 **pass**. Step-0 Hessians were measured at all 10 designs. Within V4, H1a, H1b, H2 and
  the H_R·v part of H3 **pass**.
- **H3 fails the registered slope rule** in one row: vᵀH_R v against the second difference of R at `A_6_L2/y_final`.
  PROTOCOL-v3 §2.2: "If the method itself fails (for example, step sensitivity that FD cannot resolve …), the study stops
  for review." **The study is stopped.**
  - No optimisation run (P1–P5) was started, and no Hessian refresh was made.
  - No assessment was run, and **no D5w, gap-companion, D1, D1b, D3, D4 or D6x quantity was computed.**
- The descriptive Question-1 labels (H) are reported in §4, because the protocol reports them "first, regardless of
  how the runs go". They are not a run result.

## 2. Gates (`logs/`)
| Gate | Result | Evidence |
|---|---|---|
| V1 | PASS. Enterprise HEAD 477141d with a clean tree. `s4.exe` built with the study-3 flags (MinGW g++ 16.1.0, `-O2 -ffp-contract=off`); no `s4.cpp` warnings. | `logs/build.log`, `driver/build.sh` |
| V2 | PASS. R, F, ∇R and ∇F are bit-identical (hex comparison) between `s4` and the archived `s3.cpp` (compiled unmodified in `driver/v2_s3harness.cpp`) at `A_4_L2/y_final` (base set) and `B_4_L2/y_final` (reconstructed set). F = 1.070119697e-4 matches the archived trace. | `logs/gate_V2.txt`, `results/v2/` |
| V3 | PASS for all 4 B finals. Node counts (8304, 16632, 8688, 18036) and leaf counts (519, 462, 543, 501) are exact. Working ρ_min agrees with `rhomin_w_post` to ≤ 7e-14 relative. The archive prints 12 significant digits, so the 1e-12 test resolves to ~5e-13. | `logs/gate_V3_recon.txt`, `results/sets/` |
| V4 H1a | PASS. All 10 designs: 0 diagonal and 0 gradient mismatches across the 54 passes; pass gradients bit-identical to `gradient()`; pair_errors 0. | `hessians/checks.csv` |
| V4 H1b | PASS. ê_H = 1.6e-20 to 1.2e-19 against max\|μ\| = 1.7e-5 to 7.5e-5. e_H^c ≤ ‖B_H‖_F everywhere; the `long double` consistency check holds; no "reverse-order estimate understates" event (e_H^c/e_H = 0.24–1.27). | `hessians/checks.csv` |
| V4 H1c | Reported. ê_g/tol_C1' = 1.7e-11 to 5.3e-11. | same |
| V4 H2 | PASS under the registered primary reading. All 8 columns at A_4_L1 and A_4_L2: a fall of ≥ 30× per decade on the truncation side (~100×); best e = 4e-9 to 1.7e-7. Under the stricter reading (two consecutive decades), col0 fails at both designs. | `hessians/fd/fd_checks.csv`, `logs/gate_H2H3.txt` |
| V4 H3, H_R v | PASS at all 5 designs (best e 7e-9 to 6e-8; slope ~100×). | same |
| **V4 H3, vᵀH v** | **FAIL at A_6_L2.** e(h) = 1.67e-5, 2.48e-6, 9.2e-5, 1.0e-1 for h = 1e-2…1e-5. The first decade falls only 6.7×; roundoff dominates from h = 1e-3. The other four designs pass with a single qualifying decade (45–154×). Best e ≤ 3.4e-6 everywhere, within the 1e-4 tolerance. | same |
| V4 H4 | 1/γ_last lies within [μ_min, μ_max] at all 8 finals. | `results/step0_summary.json` |
| V5a | PASS (13 checks: one-step convergence on an SPD quadratic with the Jacobi-built W_0, secant error 2.6e-14, and the skip, reset, skip-fraction and refresh-eligibility rules). | `logs/gate_V5_selftest.txt` |
| V5b | PASS. The analytic branch jump in f′ and f″ at u = 0 is exactly 0; the numeric OTI check also gives 0. | `logs/execution.log` |
| V5c | PASS. All 5,886 entries are within 4.4e-16 of the analytic values (max\|H\| = 1.56); so is the gradient. | selftest |
| V5d | PASS (all registered cases, including the mixed check, OUTER_LIMIT and the trust reset). | selftest |
| V5e | PASS. [1, 1e100, 1, −1e100] → 2. For 10⁵ terms with heavy cancellation, the Neumaier error is 4.8e-14 against an exact `Fraction` sum, with bound 2u·Σ\|x\| = 0.30 (forward-sum error 0.14). e_H^c ≤ ‖B_H‖_F. `long double` digits: 64. | `logs/gate_V5_selftest.txt` |
| V6 | PASS. Δ_0, Δ_pre(20) and Δ_pre(25) at 4:1 and 6:1 match the table to ±0.001. A 6:1 in-loop q is 0.7145 at iterate 87 and 2.0196 at iterate 90. | `logs/gate_V6.txt` |

### Diagnosis of the H3 failure (descriptive; no new evaluation; corrected 2026-10-27, see §7)
- The observed vᵀHv error scales as 1/h² for h ≤ 1e-3 at every H3 design, which is the signature of roundoff. It sits
  15–150× above the ideal floor 4u·R/(h²\|vᵀHv\|), consistent with node-sum roundoff in R of ~10u relative.
- **(Corrected.)** The step scale is not the cause. The H3 step is s = h·max(1, \|v·y\|), and \|v·y\| < 1 at all five H3
  designs (0.581, 0.877, 0.493, 0.196, 0.455), so s = h everywhere and A_6_L2 used the same steps as the four
  designs that pass.
- **(Corrected.)** The cause is a small truncation constant meeting the s⁻² roundoff of a second difference of R
  (Noether, REVIEW-noether-V4stop.md §1.2). At A_6_L2 the truncation part at s = 1e-2 is 1.67e-5 relative, 3–20×
  smaller than at the other four designs (5.1e-5 to 3.3e-4). The roundoff part is about 1e-11/s² relative. The two
  cross at s* ≈ 3e-3, leaving about half a decade of truncation-dominated range, while a 30× fall needs about
  0.74 decade. The error therefore falls only 6.7× between 1e-2 and 1e-3. The same structure appears at every design
  (the strict two-decade reading fails at all five); A_6_L2 is the extreme case.
- The OTI quantities agree with FD: Hv to 7e-9, vᵀHv to 2.5e-6, inside the 1e-4 tolerance. Noether's analysis of the
  archived data bounds an OTI error along v at about 1e-8 relative (REVIEW-noether-V4stop.md §1.3–1.5).
- No coding error was found. The failing element is the FD reference's ability to show the slope within the registered
  steps. The protocol classifies this as a stop, and it was not overridden.
- **Options for the reviewer** (none chosen, none run; any is a protocol change needing review):
  1. a v4 that applies the slope rule only where the step range is above the reference's roundoff floor, or exempts
     the second-difference comparison from it;
  2. larger registered steps for the vᵀHv comparison;
  3. a different step scale for directions. (Corrected 2026-10-27: this would address nothing, because s = h at all
     five designs; see §7.)

## 3. Protocol deviations and builder readings (declared)
1. **Noether's v3 readings: agreed, and logged before any D5w computation.** (a) The registered table values enter
   G^hi/G^lo, as implemented in `driver/d5w.py`. (b) The PASS exception includes CONVERGED_LS_LIMITED. No D5w was
   computed.
2. **H2/H3 slope-rule reading** (logged before data):
   - Primary reading: at least one consecutive step pair on the truncation side falls ≥ 30×.
   - The stricter two-decade reading is also reported. It fails more rows: H2 col0 at both designs and vᵀHv at all
     5 designs.
   - The study stops under either reading.
3. **H3 step scale:** h·max(1, \|v·y\|) for a unit direction ("steps as in H2"), logged before data. It contributes to
   the failure (option 3). (Corrected 2026-10-27: it does not contribute; max(1, \|v·y\|) = 1 at all five designs.
   See §7.)
4. **W_0 and refreshes use a C++ cyclic Jacobi eigensolver.** Reported spectra use `numpy.linalg.eigh`. This was
   exercised only in V5a.
5. **V2 harness:** `s3.cpp` was compiled unmodified, with `main` renamed by the preprocessor.
6. **V5c synthetic accumulation:** the "nodes" are the 108 row terms of ½yᵀMy plus the 5 softplus terms.
7. **OUTER_LIMIT counting** (unexercised): the outer round is the index of the refinement check, and the 10th check
   that needs a split ends the run with the split not applied.
8. **Assessment of the start designs** was not run, because the stop came before it.

## 4. Question 1 (descriptive H labels; step 0 only)
Every label is to node-sum roundoff only. τ_H covers node-sum roundoff. (Corrected 2026-10-27, see §7: per-node
evaluation error is not measured directly; the achieved H2/H3 agreement is 4e-9 to 2e-7, not the 1e-4 tolerance, and a
sign change of μ_min would need a coherent relative per-node error ≥ ε* = 1.6e-5 to 9.1e-4, per Noether.) No design
satisfies C1', so every label below is "indefinite, non-stationary", not "saddle".

| Design | Nodes | Label | κ_\|·\| | μ_min | μ_max/λ | λ_min(H_R)/λ | ‖∇F‖/tol_C1' | ‖∇F‖/(λ‖y‖) | Eigenvalues modified in W_0 |
|---|---|---|---|---|---|---|---|---|---|
| A_4_L1 | 6624 | indefinite (1 neg) | 3.8e4 | −1.8e-9 | 1783 | −1.05 | 45 | 0.0049 | 7 |
| A_6_L1 | 14904 | indefinite (1 neg) | 7.2e4 | −1.8e-8 | 1724 | −1.49 | 689 | 0.077 | 7 |
| B_4_L1 | 8304 | indefinite (1 neg) | 1.9e4 | −7.2e-9 | 1680 | −1.19 | 4873 | 0.54 | 7 |
| B_6_L1 | 16632 | indefinite (2 neg) | 1.3e4 | −2.4e-8 | 1745 | −1.63 | 836 | 0.094 | 8 |
| A_4_L2 | 6624 | indefinite (2 neg) | 1.4e5 | −1.6e-9 | 20104 | −1.43 | 12990 | 2.6 | 9 |
| A_6_L2 | 14904 | indefinite (2 neg) | 1.4e4 | −1.7e-9 | 4426 | −1.45 | 5671 | 1.15 | 7 |
| B_4_L2 | 8688 | indefinite (7 neg) | 1.2e5 | −2.1e-8 | 5748 | −6.48 | 38110 | 7.1 | 15 |
| B_6_L2 | 18036 | indefinite (3 neg) | 5.7e4 | −6.9e-9 | 5555 | −2.83 | 16080 | 3.1 | 14 |
| A_4_L2 it 84 | 6624 | indefinite (5 neg) | 4.9e5 | −5.8e-9 | 5430 | −2.55 | 18090 | 3.5 | 14 |
| A_6_L2 it 90 | 14904 | indefinite (3 neg) | 1.1e5 | −9.2e-9 | 4767 | −3.46 | 12330 | 2.2 | 12 |

- No eigenvalue is near-zero (\|μ\| ≤ τ_H), and none lies in the marginal band. All negative eigenvalues are far outside
  τ_H: \|μ_min\|/τ_H ≥ 2.7e4.
- **Study-3 "apparent condition number 10^3–10^4": "not evaluable (no PD final)" at λ1 and λ2.** Every archived final is
  indefinite on its own node set.
- **Vera's I5 rewording ("curvature along recent steps exceeds λ by 10^2–10^4"):** "consistent" at λ1
  (μ_max/λ = 1680–1783); "not supported" at λ2 (4426–20104; A_4_L2 lies outside [1e2, 1e4]).
- Eigenvector block shares for the negative eigenvectors: `hessians/<design>/eig.csv` and `summary.csv`.
- Measured Hessian cost at 12 threads, against the estimates: 4:1 267–271 s (4.6 min estimated); 6:1 595–600 s
  (10 min estimated); B sets 364–765 s. Load from other users was ~43 % of 24 logical CPUs.
- Figure: `figures/spectra.png`.

What these labels mean, as fact rather than run results:
- None of the study-3 finals is a second-order local minimiser of the regularised, discretised F on its node set. Each
  has 1–7 directions of negative curvature.
- At λ1 the A 4:1 final has a small gradient relative to the regulariser (‖∇F‖/(λ‖y‖) = 0.005) but is a saddle to
  node-sum roundoff. **(Corrected 2026-10-27, see §7: read "has an indefinite H_F at a non-stationary point
  (‖∇F‖/tol_C1' = 45)".)**

## 5. Not run
- P1–P5, all assessments, post-run and S2 Hessians: **not run** (stopped at V4).
- D1, D1b, D3, D4, D5w and the gap companion, D6, D6x, D7, the T labels and the rate: **not computed.** They may be run
  later under an unchanged or revised protocol, with the date recorded.
- The (marginal) rule does not arise, because no D1b label exists.

## 6. Reproduction (from `warp_hypad/study4`; Python = Anaconda `-I`)
```
sha256sum PROTOCOL-v3.md; cat logs/protocol-v3.sha256
driver/build.sh                                         # V1 (assembles s4.cpp from s3.cpp + driver/src_pieces)
python -I driver/v6_baselines.py ..                     # V6
driver/s4.exe selftest --passes driver/hessian_passes.txt --out results/v5; python -I driver/v5_check.py .   # V5
driver/s4.exe recon --run ../study3/runs/B_4_L2 --ws 4:1 --out results/sets/B_4_L2.elems.txt               # V3
driver/s4.exe v2dump ... ; driver/v2_s3harness.exe ...  # V2 (commands in logs/execution.log)
driver/step0_hessians.sh; python -I driver/hess_analyze.py $PWD step0                                     # V4 H1, H4, Q1
python -I driver/fd_checks.py $PWD prepare; driver/fd_run.sh; python -I driver/fd_checks.py $PWD analyse  # V4 H2/H3
```
Driver hashes: `logs/driver.sha256` (`s4.cpp` `ef247445…9463`). Prepared but unused run tooling (`driver/lane.sh`,
`assess.sh`, `analyze_runs.py`, `d5w.py`) is untested on real runs.

## 7. Corrections (2026-10-27, Ada, after Noether's ruling REVIEW-noether-V4stop.md)
These corrections were made after the stop and after Noether's review. Corrected text in §2–§4 is marked
"(Corrected 2026-10-27)". No result was recomputed by the builder; the replacement figures are Noether's, from the
archived data. The original wording is kept here:
1. **§2 Diagnosis, step scale (error in my mechanism).** Original: "The registered step scale for a unit direction is
   h·max(1, |v·y|). At A_6_L2, |v·y| = 0.196, so the largest step is 1e-2. The O(h²) region is visible only between
   1e-2 and 1e-3, where the error falls 6.7× instead of ≥ 30×." The implied claim that the step scale compressed the
   sweep at A_6_L2 is wrong: max(1, |v·y|) = 1 at all five H3 designs (|v·y| = 0.581, 0.877, 0.493, 0.196, 0.455), so
   every design used s = h. I stated the mechanism without checking it numerically against the other designs.
2. **§2 Diagnosis, actual mechanism (added).** A small truncation constant at A_6_L2 (1.67e-5 relative at s = 1e-2,
   3–20× below the other designs) meets the s⁻² roundoff of the R second difference (≈ 1e-11/s² relative), crossing
   at s* ≈ 3e-3. This is a generic limitation of R second differences over decade steps, extreme at A_6_L2.
3. **§2 option 3 and §3 item 3.** Original: "It contributes to the failure (option 3)." Withdrawn; it does not
   contribute, and option 3 addresses nothing. Noether rejected options 2 and 3 and adopted a modified option 1
   (PROTOCOL-v4, NV2).
4. **§4 limitation sentence.** Original: "per-node evaluation error is unmeasured below the H2/H3 resolution (~1e-4
   relative), which was itself only partly established (H3 failure)." 1e-4 is the tolerance, not the achieved
   agreement (4e-9 to 2e-7); ε* = 1.6e-5 (A_4_L1) to 9.1e-4 (Noether §3).
5. **§4 "saddle" (NV4 erratum).** Original: "is a saddle to node-sum roundoff". No step-0 design satisfies C1'
   (‖∇F‖/tol_C1' = 45 to 38,110). Read "indefinite H_F at a non-stationary point". "None is a second-order local
   minimiser" remains true.

Status: PROTOCOL-v4.md (SHA-256 `0d77d7d52abf819d19c9eb8cf66cf0dbdbfb6d7261e0070f1341486499d7c3ad`,
`logs/protocol-v4.sha256`) applies NV1–NV5 and NV6, and was written after the V4/H3 data were seen. **Nothing has been
re-run.** The study remains stopped until an independent reviewer (Vera) clears v4.
