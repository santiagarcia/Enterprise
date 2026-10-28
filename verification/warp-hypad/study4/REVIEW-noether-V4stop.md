# Review of the study-4 stop at gate V4/H3 (Noether, 2026-10-27)

Reviewer: Noether (mathematical foundations). Builder: Ada. I reviewed and cleared PROTOCOL-v1 to v3. I did not
write the driver, compute any Hessian or run any FD evaluation. **Independence limit:** the H3 tolerance and slope
rule I am now asked to revise are rules I approved in v3. The change I propose was written **after** I saw the
failing row and after I recomputed the replacement quantities below, so it is not a blind pre-registration. It is
justified on general grounds (section 2), and it is applied uniformly to every H3 row. A second reviewer (Vera or
Gauss) should check that justification before the study resumes.

**What I ran:** only archived-data analysis. I ran no solver and no new FD evaluation (the 15-minute allowance was not
used). The inputs were `hessians/fd/*.out` (the archived R(y±sv) and ∇R(y±sv)), `hessians/<design>/H_R_fwd.bin`, `B_H.bin`,
`gR_fwd.bin`, `y.bin` and `meta.txt`, together with `logs/gate_H2H3.txt`, `hessians/fd/fd_checks.csv`, `results/step0_summary.json`,
`driver/fd_checks.py` and `REPORT.md`. The scripts are in the reproduction section.

## 1. Diagnosis of the H3 vᵀHv row at A_6_L2

**Verdict:** this is a limitation of the FD reference, and the evidence for that is strong. An OTI error along v is
excluded far below the tolerance. However, part of Ada's stated mechanism is incorrect.

1. **The step-scale explanation is wrong.** The H3 step is s = h·max(1, |v·y|). At all five H3 designs |v·y| < 1
   (0.581, 0.877, 0.493, 0.196, 0.455), so s = h everywhere. A_6_L2's steps are identical to those of the designs
   that pass, so the step scale did not compress anything there. Option 3 (a different step scale) therefore
   addresses nothing.
2. **The actual mechanism is a small truncation constant against a larger roundoff floor.** For φ(t) = R(y+tv), the
   central second difference has error s²φ''''/12 + O(ε_R/s²).
   - At A_6_L2 the truncation part at s = 1e-2 is 1.67e-5 relative. That is 3–20× smaller than at the other four
     designs (5.1e-5 to 3.3e-4).
   - The roundoff part, inferred from the h = 1e-5 row, is about 1e-11/s² relative, so ≈ 1e-5 at s = 1e-3. In the
     second-difference numerator this is 50–500 ulp(R). That is well inside the a-priori node-sum bound R_bound =
     1.7e-16 ≈ 1.3e4 ulp(R).
   - Truncation and roundoff therefore cross at s* ≈ (1e-11/0.167)^(1/4) ≈ 3e-3. The registered sweep has only
     about half a decade of truncation-dominated range. A 30× fall needs about 0.74 decade, so it cannot be shown with
     decade steps.
   - The same structural problem appears at every design: the strict (two-decade) reading of the vᵀHv rows fails at
     all five. A_6_L2 is the extreme case of a generic limitation of R second differences, not an outlier.
3. **Signed reconstruction from the archived FD outputs (A_6_L2, relative to vᵀH_R v = 7.0274e-7):**

   | s | Second difference of R | Projected gradient difference vᵀ(∇R₊−∇R₋)/(2s) | R–∇R residual c(s) |
   |---|---|---|---|
   | 1e-2 | +1.667e-5 | +3.337e-5 | 7.87e-12 |
   | 1e-3 | −2.477e-6 | +3.337e-7 | 7.66e-14 |
   | 1e-4 | +9.206e-5 | +2.408e-9 | 6.29e-16 |
   | 1e-5 | +1.000e-1 | −1.30e-8 | 1.7e-13 |

   - The projected gradient difference shows a clean O(s²) fall (100×, 139×) down to 2.4e-9. **So the OTI vᵀH_R v
     equals the directional second derivative of the computed ∇R to ≈ 1e-9 relative.**
   - The OTI value is vᵀ(H_R v), and the H_R v row passed (7.1e-9). Since ‖H_R v‖/|vᵀH_R v| = 3.8, an OTI error
     along v is bounded by ≈ 3e-8 relative. That is 3,000× below the 1e-4 tolerance.
4. **R and ∇R are the same function along v.** The residual c(s) = |(R₊−R₋)/(2s) − ½vᵀ(∇R₊+∇R₋)| has truncation
   (s²/3)|φ'''| and roundoff of order s⁻¹. It falls 103× and 122× per decade at A_6_L2, and c_min/|φ'| = 1.2e-7. All
   five designs pass the primary slope reading (ratios 100/61, 100/77, 100/433, 103/122, 100/19.5), and c_min/|φ'| ≤
   6.6e-6 everywhere.
5. **Richardson identity (independent tie between the R values and the Hessian).** Where both quantities are
   truncation-dominated, the second-difference error must be exactly half the projected gradient-difference error,
   because both are governed by φ'''' (1/12 against 1/6). At s = 1e-2 the measured ratio is 0.4996 (A_6_L2) and
   0.4998–0.5000 at the other four designs.
   - An error in vᵀH_R v of size δ would shift both errors by the same δ. That would move the ratio away from ½, by
     about 0.5·δ/1.7e-5 at A_6_L2. A ratio of 0.4996 bounds δ at roughly ≤ 1.4e-8 relative.

**Conclusion (1):** an OTI Hessian error along v is excluded at about the 1e-8 level, by three independent
derivations from the archived data. The failed row shows that the R second-difference reference has too little
truncation-dominated range in the registered sweep. Under v3 §2.2 the stop was correct. The H3 data do not need to be
regenerated.

## 2. Ruling on the options, and the replacement rule

- **Option 2 (larger steps) is rejected.**
  - Choosing the step from this row's failure would be tuning.
  - At s ≥ 1e-1 the O(s⁴) terms and the softplus branch at u = 0 come into play: V5b establishes continuity only of
    f′ and f″, and the second difference needs smooth φ''''.
  - The R second difference would still be the weakest reference.
- **Option 3 (a different step scale) is rejected.** As shown in 1.1, it does not address the cause.
- **Option 1 is adopted in a modified form.**
  - The second-difference row keeps its accuracy gate. Its slope rule is applied only where the reference can
    resolve it, according to a resolvability test fixed in advance and computed from each row's own sweep.
  - The slope evidence that the row was meant to supply (that H_R is the Hessian of the *value* R, not only the
    Jacobian of ∇R) is moved to a quantity that is not divided by s²: the R–∇R consistency residual. This becomes a
    new gate.
  - The general justification is that the roundoff of a second difference of an n-node sum scales as s⁻², while
    first-difference quantities scale as s⁻¹. With fixed decade steps, the truncation-dominated window of an R second
    difference is generically about one decade or less (section 1.2). A slope rule applied there tests the step
    grid, not the Hessian.

Applied to the archived data, the resolvability test gives:

| Row | Resolvable? | Result |
|---|---|---|
| A_4_L1 | yes (ρ = 7.7e-7 ≤ 1.7e-6) | pass, 45.5× |
| A_4_L2 | yes (1.5e-6 ≤ 6.2e-6) | pass, 55.3× |
| B_4_L2 | yes (4.1e-6 ≤ 1.1e-5) | pass, 154× |
| B_6_L2 | yes (5.1e-6 ≤ 6.7e-6, close) | pass, 71.6× |
| A_6_L2 | **no** (1.0e-5 > 5.6e-7) | "slope not resolvable by the R reference" |

- The new R–∇R gate (b) passes at all five designs.
- I state plainly that I computed these outcomes before writing the rule. The rule's validity rests on the error
  model, not on blindness.
- What would refute this ruling: a resolvable row failing its slope test, or failing (b), under an unchanged OTI
  code.

## 3. Step-0 finding: all 10 designs indefinite

**Robust to roundoff, yes.**

- Let S = B_H/((n−1)u) be the entrywise sum of the |per-node contributions| (u = 2⁻⁵³, which is conservative).
- If every per-node Hessian contribution carries a relative error of at most ε, then ‖ΔH‖₂ ≤ ε‖S‖₂ (Weyl's
  inequality). The sign of μ_min can flip only if ε ≥ ε* = |μ_min|/‖S‖₂.
- **ε* = 1.6e-5 (A_4_L1, smallest), 2.0e-5 (A_4_L2), and 6.8e-5 to 9.1e-4 elsewhere.** The first-order estimate
  along the eigenvector, |μ_min|/(|w|ᵀS|w|), is 1.4e-4 to 5.2e-2.
- Floating-point per-node evaluation roundoff is about (operation count)·u. Even 10⁴u ≈ 1e-12 is ≥ 10⁷ below ε*.
- The negative counts are also robust: min|μ|/‖S‖₂ ≥ 2.0e-6 (A_4_L2_it84).

Two qualifications:

1. **The "~1e-4 resolution" in v3 overstates the uncertainty. I wrote that wording (N3), and I correct it here.**
   - 1e-4 is the H2/H3 tolerance. The achieved FD agreement is 4e-9 to 2e-7 on the sampled columns and directions.
   - Errors that FD cannot see are errors that are consistent with the computed R, which means a modelling
     difference in R itself. By the R-scope rule that is outside these labels: they describe the discretised F
     implemented in s3, to which V2 is bit-identical.
   - A coherent error of 1.6e-5 to 1e-4 could flip A_4_L1 or A_4_L2, so only a code defect could do it, not roundoff.
2. **"Saddle" is the wrong word.** No design satisfies C1' (‖∇F‖/tol_C1' = 45 to 38,110). These are non-stationary
   points with an indefinite H_F, which is unremarkable for a nonconvex F. REPORT §4's "a saddle to node-sum
   roundoff" (A_4_L1) should read "indefinite H_F at a non-stationary point". The statement "none is a second-order
   local minimiser" remains true.

**Implications for the design.**

- The negative curvature is weak, |μ_min|/μ_max = 2e-5 to 1e-3, and |μ_min|/λ is 0.05 to 5.5.
- W_0 = Q diag(1/max(|μ|, λ)) Qᵀ is SPD. In a local quadratic model it makes the gradient component along a negative
  eigenvector grow by a factor of about 1 + |μ|/max(|μ|, λ) per iterate: 1.05 to 2.0.
- *Speculative (local quadratic model; to be checked by the runs):* at A_4_L1, the start of P4 `A4L1f`, the gradient
  share in the negative subspace is 1.2e-3. That is ≈ 3.6e-12, already about 19× below tol_C1' (≈ 6.9e-11 at λ1),
  with a growth factor of 1.047. In the model, C1' ∧ C2' could therefore be met in the indefinite region before an
  escape. This is plausible for A4L1f; at the other starts the share is 1.2e-2 to 2.4e-1.
- Under v3 such a run would be CONVERGED, and D1b would label it "saddle". **That is the correct, informative
  outcome.** The convergence criteria should stay first-order: the study asks *whether* converged points are
  second-order (Question 2), and the second-order check at convergence answers that question.
- Adding negative-curvature steps would change the solver and the comparison with study 3, and would need new gates.
  It is **not required**.
- What is required: the protocol must state that CONVERGED is first-order only, and that a saddle final is an
  anticipated outcome, not a failure.

## 4. PROTOCOL-v4 change specification (dated 2026-10-27; Noether; post-data, see the independence note)

Each text below replaces the cited v3 text verbatim. Mark each change "(v4, NV…)". The H2/H3 analysis is re-run on
the archived `hessians/fd/*.out` files; no new evaluation is made. Record a new SHA-256 for v4.

**NV1 (required). §2.2 H2, the sentence beginning "On the truncation side".** Replace with:
> On the truncation side (the registered steps up to and including the step with the smallest e), at least one pair of
> consecutive registered steps must satisfy e(h_k)/e(h_{k+1}) ≥ 30 (one decade, consistent with O(h²)), and the best
> step must give e ≤ 1e-4. (v4, NV1) The reading that requires two consecutive qualifying pairs is reported
> descriptively and is not a gate.

**NV2 (required). §2.2, the H3 bullet.** Replace it entirely with:
> - **H3 (directional).** One seeded random unit vector v (seed 4) at `A_4_L1/y_final`, `A_4_L2/y_final`,
>   `B_4_L2/y_final` (reconstructed set), `A_6_L2/y_final` and `B_6_L2/y_final` (reconstructed set). Steps
>   s = h·max(1, |v·y|), h ∈ {1e-2, 1e-3, 1e-4, 1e-5}. Write R± = R(y ± sv), g± = ∇R(y ± sv). The same evaluations
>   give three comparisons:
>   - **(a) H_R v (gate):** e_a(h) = ‖(g₊ − g₋)/(2s) − H_R v‖/‖H_R v‖. Tolerance and slope rule as in H2.
>   - **(b) (v4, NV2) R–∇R consistency (gate):** c(h) = |(R₊ − R₋)/(2s) − ½vᵀ(g₊ + g₋)|. Its truncation term is
>     (s²/3)|d³R/dt³| and its roundoff scales as s⁻¹. The slope rule is as in H2, and min_h c(h) ≤ 1e-4·|½vᵀ(g₊ + g₋)|
>     at the minimising step.
>   - **(c) vᵀH_R v (tolerance gate; slope rule where resolvable):** e_c(h) = |(R₊ − 2R(y) + R₋)/s² − vᵀH_R v| /
>     |vᵀH_R v|. Gate: min_h e_c ≤ 1e-4.
>     - (v4, NV2) Because the roundoff of a second difference of an n-node sum scales as s⁻², it is estimated from
>       the two smallest steps: ρ(s) = max_{k ∈ {3,4}} e_c(h_k)·(s_k/s)².
>     - A consecutive pair (h_k, h_{k+1}) on the truncation side is **resolvable** if ρ(s_{k+1}) ≤ e_c(h_k)/30.
>     - If at least one pair is resolvable, at least one resolvable pair must satisfy e_c(h_k)/e_c(h_{k+1}) ≥ 30
>       (gate).
>     - If no pair is resolvable, the row is reported as "slope not resolvable by the R reference", with ρ and e_c.
>       This is not a method failure, and (a) and (b) carry the slope evidence.
>   - Reported descriptively: the projected gradient difference vᵀ(g₊ − g₋)/(2s) against vᵀH_R v, and, at the
>     largest step, the ratio of the (c) error to that projected error (theory ½ when both are truncation-dominated).

**NV3 (required). §2.2 H1b "Reference and limitation", and the §8 D1b "Limitation" wording.**
- In H1b, replace "Above about 1e-4 relative it is bounded by H2/H3; **below that resolution it is unmeasured**."
  with:
  > (v4, NV3) It is not measured directly. FD differentiates the computed R and so cannot see errors that are
  > consistent with R. The achieved H2/H3 agreement (best e) is reported, not the 1e-4 tolerance. At each Hessian
  > design ε*_min = |μ_min(H_F)|/‖S‖₂ and ε*_all = min_i|μ_i|/‖S‖₂ are reported, with S = B_H/((n−1)u): the smallest
  > coherent relative per-node error that could change the sign of μ_min, or of any μ (Weyl).
- In D1b, the limitation sentence becomes:
  > "τ_H covers node-sum roundoff; a sign change by per-node evaluation error requires a coherent relative error
  > ≥ ε*_all (reported)."
- Add ε*_min and ε*_all to the §2.3 reported list. These are descriptive only; the marginal-band rule is unchanged.

**NV4 (required). §2.3 Registered labels, after the "singular to τ_H" label.** Insert:
> (v4, NV4) These labels describe H_F only. The word "saddle" is used only at a design that satisfies C1'. Otherwise
> the label is written "indefinite (n negative), non-stationary (‖∇F‖/tol_C1' = x)".

The correction to REPORT §4 is required as well (a report erratum, not a protocol change).

**NV5 (required). §4.1, after the terminal-state table.** Insert:
> (v4, NV5) C1' ∧ C2' are first-order criteria, and CONVERGED asserts first-order stationarity only.
> Second-order character is established only by D1b. All ten step-0 designs have indefinite H_F, so a CONVERGED
> final that D1b labels "degenerate or saddle" is an anticipated outcome and is reported as a result, not as a solver
> failure. No negative-curvature step is added (it would change the solver and its comparison with study 3).

**NV6 (recommended, descriptive).** §7, "Second-order check at convergence". At each converged final labelled
"degenerate or saddle", report the following:
- ‖P_neg∇F‖/tol_C1';
- the function-value confirmation F(y ± α·w_min) − F(y), with α = (2·10³·b_F/|μ_min|)^(1/2). This costs two
  evaluations, and the predicted decrease is ½μ_minα².

## 5. Reproduction (from `warp_hypad/study4`; Anaconda Python with `-I`; archived data only)
The scripts were run from the session scratchpad (not saved in the study folder). Their essential lines:

    v = default_rng(4).standard_normal(108); v /= norm(v); H = fromfile(H_R_fwd.bin).reshape(108,108)
    for each "v4 s R+ R- g+ g-" line of hessians/fd/<d>.out (R0 = first line):
      sd = (R+ - 2R0 + R-)/s^2 ;  pg = v.(g+ - g-)/(2s) ;  c = |(R+ - R-)/(2s) - v.(g+ + g-)/2|
      report (sd - vHv)/vHv, (pg - vHv)/vHv, c, and the ratio (sd - vHv)/(pg - vHv) at s = 1e-2
    S = B_H/((n-1)*2**-53); mu, W = eigh(sym(H_R) + lam*I); eps* = |mu_min|/||S||_2; |w|.S.|w|
    g = gR_fwd + lam*y; ||P_neg g||/||g||

## 6. Candidate lessons (proposed only; for the main session)
- A slope rule on an R second difference over fixed decade steps is generically under-resolved: roundoff scales as
  s⁻² for a node sum. Use first-difference consistency residuals for slope evidence, and keep second differences as
  tolerance checks. (Scope: FD verification of OTI/AD Hessians of sums; evidence: this stop.)
- Verify the stated mechanism numerically. Here the claimed step-scale compression was max(1, 0.196) = 1.
- "Indefinite" and "saddle" are different labels; saddle requires stationarity.

## v4 changes required
NV1, NV2, NV3, NV4 (with the REPORT §4 erratum) and NV5. NV6 is recommended. After v4 is applied, re-run the H2/H3
analysis from the archived `hessians/fd/*.out` files and confirm that V4 passes under NV1/NV2. My re-check is limited
to NV1–NV5. Because this rule was written post-data, I recommend a second reviewer check section 2's justification
before P1 starts.

**CLEARED TO RESUME UNDER v4 once NV1–NV5 are applied** (and the H2/H3 re-analysis under v4 passes from the archived
data, with no change to the OTI code).
