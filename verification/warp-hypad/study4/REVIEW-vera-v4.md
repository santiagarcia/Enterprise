# Second independent review of PROTOCOL-v4 (Vera, 2026-10-27)

**Role:** reviewer only. I did not design study 4, write the V4 ruling, write v4, or run any solver, FD or OTI
evaluation. I edited nothing in this folder except this file.
**Independence limits:**
- I read Noether's ruling before I recomputed anything, so anchoring is possible. To reduce it, I recomputed from the
  raw `hessians/fd/*.out` lines with my own script. I did not reuse `fd_checks.py` or Noether's scratch code.
- PROTOCOL-v4 cites earlier Vera items from study 3 (I4, I5, O5). None of them concerns H2/H3. I have no memory of
  them in this session.

**Inputs:**
- `PROTOCOL-v4.md`, `PROTOCOL-v3.md` (I diffed v3 against v4 in full), `REPORT.md` (including §7),
  `REVIEW-noether-V4stop.md`, `logs/gate_H2H3.txt` and `hessians/fd/fd_checks.csv`;
- the archived `hessians/fd/*.out` and `*.dirs`, and `hessians/<d>/H_R_fwd.bin` and `y.bin`;
- `driver/fd_checks.py` and `logs/{driver.sha256,execution.log,build.log}`.

## 0. Integrity checks (all pass)
- `sha256sum` of `PROTOCOL.md`, `PROTOCOL-v2.md`, `PROTOCOL-v3.md` and `PROTOCOL-v4.md` equals the hash in
  `logs/protocol.sha256`, `protocol-v2.sha256`, `protocol-v3.sha256` and `protocol-v4.sha256` respectively. v4 is
  `0d77d7d5...c3ad`, and v1-v3 are unchanged.
- The full v3->v4 diff touches only the places that §15 declares: the header, H1b, H2, H3, the §2.3 eps*/NV4 lines,
  the §4.1 NV5 paragraph, the §7 NV6 item, the §8 D1b limitation sentence, the §9 V4 row and resumption procedure,
  and §15. No D-rule, threshold, solver, budget or run definition changed.
- The `driver/s4.cpp` SHA-256 is `ef247445...9463`. This equals `logs/driver.sha256`, `logs/build.log` and v4 §9.
  `fd_checks.py` (`797b2cbe...`) is unchanged.
- The five `hessians/fd/*.out` SHA-256 values equal the v4 §9 table exactly.
- The direction v regenerated from `default_rng(4)` equals the v in the `.dirs` files bit-for-bit. All five H3
  designs used s = h (|v.y| < 1).

## 1. Is Noether's diagnosis supported? Recomputation from the archived FD outputs
Quantities are relative to vᵀH_R v, except c, which is relative to |½vᵀ(g₊+g₋)|. Steps are h = 1e-2...1e-5.

| Design | (a) e_a | (b) c/\|φ'\| | (c) e_c | ρ(s₂) vs e_c(h₁)/30 | Signed projected-gradient error | Richardson ratio at s₁ |
|---|---|---|---|---|---|---|
| A_4_L1 | 1.28e-4, 1.28e-6, 1.27e-8, 2.05e-8 | 5.0e-3, 5.0e-5, 8.2e-7, 1.1e-6 | 5.15e-5, 1.13e-6, 7.66e-5, 7.66e-3 | 7.7e-7 ≤ 1.7e-6: resolvable, 45.5x | 1.03e-4 ... 9.99e-9 (h₃) | 0.49996 |
| A_4_L2 | 1.23e-3 ... 6.37e-8 | 3.1e-3 ... 4.0e-7 | 1.85e-4, 3.35e-6, 1.47e-4, 1.08e-2 | 1.5e-6 ≤ 6.2e-6: resolvable, 55.3x | ... 5.1e-9 | 0.49984 |
| B_4_L2 | 1.04e-3 ... 3.42e-8 | 5.9e-4 ... 1.4e-8 | 3.26e-4, 2.12e-6, 1.94e-4, 4.12e-2 | 4.1e-6 ≤ 1.1e-5: resolvable, 153.9x | ... 1.4e-10 | 0.49981 |
| **A_6_L2** | 4.93e-5, 4.93e-7, **7.12e-9**, 9.61e-8 | 1.45e-3, 1.42e-5, **1.16e-7**, 3.2e-5 | 1.67e-5, 2.48e-6, 9.21e-5, 1.00e-1 | **1.0e-5 > 5.6e-7: not resolvable** | 3.34e-5, 3.34e-7, **2.41e-9**, -1.30e-8 | **0.49949** |
| B_6_L2 | 5.38e-4 ... 5.42e-8 | 1.3e-2 ... 6.6e-6 | 2.01e-4, 2.81e-6, 2.35e-4, 5.08e-2 | 5.1e-6 ≤ 6.7e-6: resolvable, 71.6x | ... 9.9e-10 | 0.49991 |

**Slope ratios on the truncation side:**
- (a): about 100x at every design.
- (b): 100/61, 100/77, 100/433, 103/122 and 100/19.5/1.0.
- **Result:** under v4 as written, all five rows pass (a), (b) and (c). Every number in Noether's §1 and §2 tables
  reproduces, apart from the two points below.

**Verdict on the diagnosis: supported.** The failing v3 row is a limitation of the R-second-difference reference. It is
not evidence of an OTI error. Three corrections to how strongly the evidence is stated (none changes the conclusion):

1. **"Three independent derivations" is overstated. There are two independent lines of evidence.**
   - The projected-gradient error is vᵀ times the (a) error vector. It is the same data as (a), not an independent
     check.
   - The bound "‖H_R v‖/|vᵀH_R v| = 3.8 => ≈ 3e-8" is also derived from (a).
   - The R-∇R residual (b) does not involve H_R at all. It ties ∇R to R, which is what makes the chain
     H_R -> ∇R -> R valid, but it bounds no Hessian error by itself.
   - The two independent lines are therefore (i) the gradient finite differences, with best (a) 7.1e-9 and a projected
     error of 2.4e-9 at A_6_L2, and (ii) the Richardson identity, which is the only one that uses the R values.
2. **The Richardson bound is a factor of 2 too tight, and it is limited by roundoff.**
   - With sd_err = T/2 - δ and pg_err = T - δ, the ratio deviates from ½ by -δ/(2(T-δ)). So |δ| ≈ 2|r-½||pg_err(s₁)|.
     Noether's "0.5·δ/1.7e-5" uses half the correct denominator.
   - At A_6_L2 I get r = 0.49949, not 0.4996, which gives δ ≈ 3.4e-8.
   - The roundoff of the second difference at s₁ is itself ρ(s₁) ≈ 1.0e-7, which exceeds the measured residual
     |sd_err - pg_err/2| = 1.7e-8. The Richardson check is therefore consistent with δ = 0 at the ~1e-7 level. It does
     not establish 1.4e-8.
   - The ~1e-8 statement holds through line (i) alone.
3. **Corrected statement:** an OTI error along v is excluded at ≲ 3e-8 (gradient-difference evidence) and at ≲ 1e-7
   (independent R-value evidence) at all five designs. That is ≥ 1000x below the 1e-4 tolerance.
   - REPORT §2 repeats "≈ 1e-8 ... three independent derivations". It should cite this wording.
   - Neither correction blocks resumption.

The mechanism holds:
- the truncation constant at A_6_L2 is 3-20x smaller than at the other designs;
- the roundoff is ∝ s⁻² (1e-1 at s = 1e-5, 9.2e-5 at s = 1e-4, so ρ is set by h₄);
- the strict two-decade reading fails at all five designs, so this is a generic limitation.

## 2. Are NV1-NV6 general, or tailored? Counterexample test
**Generality of NV2.**
- The resolvability test follows from a textbook error model (second differences ∝ ε/s², first differences ∝ ε/s),
  applied uniformly. It is not tuned to A_6_L2:
  - A_6_L2 is non-resolvable by a factor of 18;
  - it would become resolvable only if the factor 30 were reduced to below 1.67;
  - raising the factor above 39.6 makes B_6_L2 non-resolvable, which only exempts it.
- So **every factor above 1.67 gives the same archived outcome.** The threshold is not knife-edge. The other side of
  this result is that (c)'s slope gate has little power, so the real protection comes from (a). See below.
- **NV1** registers the reading that was logged before the data. `fd_checks.py`, whose docstring states both readings,
  has mtime 21:34. The FD outputs were written from 22:57 to 23:05. mtime is weak provenance, but it is consistent.
  The v3 phrase "for at least two consecutive steps" also reads most naturally as one ratio. NV1 changes no archived
  outcome under that primary reading.

**Counterexample test** (scratch script; archived FD data held fixed, and the OTI Hessian replaced by H_R + E). There
were three error types, each at ±δ with δ in {1e-3, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 1e-7}:
- `vv`: E = ±δ|vᵀHv|vvᵀ, so vᵀHv is wrong by δ relative;
- `sym_vw`: E = ±δ‖Hv‖(vwᵀ+wvᵀ) with w ⊥ v, so Hv is wrong and vᵀHv is exact;
- `ww_perp`: E = ±δ‖H‖₂wwᵀ, which is invisible along v.

Results, as the v3 / v4 / v4-with-R2 pass outcome for each wrong Hessian (P = the wrong Hessian passes):

| Case | v3 | v4 | v4 + R2 |
|---|---|---|---|
| Any design, `vv` or `sym_vw`, δ ≥ 1e-3 | F | F | F |
| A_4_L1, A_6_L2, B_6_L2: `vv` δ = ±1e-4 | F | F | F |
| **A_4_L2 `vv` δ = +1e-4** | F | **P** | F |
| A_6_L2 `vv` δ = +3e-6 | F | P | P |
| `sym_vw` δ ≤ 3e-5 at A_4_L2 and B_4_L2, δ ≤ 1e-5 at B_6_L2 | P | P | P |
| Any `ww_perp` (even 1e-3‖H‖) | P | P | P |
| True Hessian at A_6_L2 | **F** | P | P |

**Findings, ranked by severity x likelihood x consequence:**

- **F1 (medium severity; low likelihood; no consequence for the archived data). A loophole in NV2(c).**
  - If e_c is smallest at the largest step h₁, the truncation side has no pairs. "No pair is resolvable" then holds
    vacuously, and the row passes on tolerance alone.
  - A minimum at h₁ is the signature of a bias, not of roundoff. Under v3, and under H2/(a)/(b), the same pattern is a
    FAIL (the slope rule over an empty set is false).
  - Concretely, a vᵀHv error of +1e-4 at A_4_L2 passes v4. It fails v3, and it fails v4 with R2 below. (a) misses it
    there because ‖Hv‖/|vᵀHv| = 5.2 dilutes it to 1.9e-5, and e_a(h₁) = 1.2e-3 still falls 38x against that floor.
  - In the archived data every row has its minimum at h₂, so this cannot affect the V4 decision. H3 is applied only to
    these archived step-0 data.
- **F2 (low).** v4's (c) slope gate rarely binds. The protection against vᵀHv errors comes from (a), through
  ‖ΔHv‖ ≥ |Δ(vᵀHv)| (vᵀHv is computed from the same matrix).
  - At A_6_L2, v4 still catches every `vv` error ≥ 1e-5 and every `sym_vw` error ≥ 3e-6.
  - Compared with v3, v4 gives up only detection of vᵀHv errors of about 3e-6 at A_6_L2. That is 30x below the
    tolerance, and the same v3 rule rejects the true Hessian there. **This is an acceptable, error-model-based
    relaxation, not a rescue.**
- **F3 (low; pre-existing; unchanged by v4). Coverage.** One random direction plus four columns at two designs cannot
  see errors with Ev = 0 (`ww_perp` passes everything under v3 and v4).
  - The defence against structured mapping errors is H1a, V5c (a known dense Hessian) and H2.
  - v4 §2.2 should keep this limitation explicit. It does say that H2/H3 "cannot" catch isolated off-diagonal mapping
    errors (H5).
- **F4 (low). Builder reading (ii) is framed as conservative, but it is lenient.** Using k = 3 can only enlarge ρ. That
  makes pairs *less* resolvable, which exempts them from the slope gate. The archived outcome is unaffected (ρ is set by
  h₄ at every design). The REPORT should state that the reading is lenient.

## 3. Do NV3-NV6 change any registered decision in a way that favours a conclusion?
No.
- **NV4:**
  - D1b is defined only for converged runs.
  - CONVERGED and CONVERGED_LS_LIMITED both require C1' (§4.1).
  - So "saddle" stays available wherever D1b applies. NV4 changes only the descriptive step-0 wording, and it does so
    in the cautious direction.
  - D1, D1b, D3, D4, D5w, D6x, the marginal rule, τ_H and every threshold are textually unchanged.
- **NV5** is interpretive. No v3 rule counted a saddle final as a solver failure, and "distinct stationary points" (D3)
  already covered it. NV5 asserts less: CONVERGED means first-order only.
- **NV3** replaces a loose limitation sentence with a quantity, ε*, that is a correct sufficient bound. ‖ΔH‖₂ ≤ ‖|ΔH|‖₂
  ≤ ε‖S‖₂ for nonnegative S, and Weyl's inequality then applies. NV3 is descriptive only.
- **NV6** is descriptive and "not acted on". One technical point for it:
  - at α = (2·10³b_F/|μ_min|)^½, the first-order term |∇F·w|α ≤ tol_C1'·α can exceed the predicted -10³b_F;
  - for example, at A_6_L2-like values, ≈ 4e-13 against 1.7e-13;
  - each one-sided value can therefore show the wrong sign at a genuine saddle;
  - **recommended:** also report the symmetric value ½[F(y+αw)+F(y-αw)] - F(y), which uses the same two evaluations
    and cancels the odd-order terms.

## 4. Are the §9 resumption steps sound?
Mostly sound:
- the protocol hash and the `s4.cpp` hash are checked;
- the five `.out` hashes are pinned, and they match;
- the analysis goes into a new hashed script, and the v3 script is kept;
- the gates that v4 does not touch are carried over, which is justified because their inputs and code are unchanged.

**Gap (F5, medium):** the v4 script reads `H_R_fwd.bin`, `y.bin` and `meta.txt`. The OTI Hessian is the object under
test, but these files are not hash-pinned in §9. Their mtimes (22:21) predate the FD outputs, so there is no evidence of
change. I record their current hashes here, before any v4 analysis:

| File | SHA-256 |
|---|---|
| A_4_L1/H_R_fwd.bin | 7a697703b416f007dc90499c3c980ea3c5a154cd6616b334fe4a7d6781c6c1b1 |
| A_4_L1/y.bin | dfb00f93a87dce39cbefbb03c9dfb2b775577886eb8d0895759ee4c73150bd5a |
| A_4_L1/meta.txt | 41e1b8a93a180630a66c97176d0f4c2c78088346e215af5b4e189ee44276ba82 |
| A_4_L2/H_R_fwd.bin | 11a70bae2aa4cd7e48cefcd5e06ba17bfc65f52242b6bdd5f15d3c224285b27a |
| A_4_L2/y.bin | 526ce11b6b5d77108f1f252407e3d257cc956acb65cab4d259f698baec8a8162 |
| A_4_L2/meta.txt | e706826a6682c0955cd40f2808d195b1b6afd4287ee4cb6d466354567259040d |
| A_6_L2/H_R_fwd.bin | 8624a733e390a005b46a831ec7ad7434648f68beacddb6521063782455dad64e |
| A_6_L2/y.bin | a0afb16db415b3a370634375cca278eec3b3b879fc03cc82cf2caa6ffb0d6d64 |
| A_6_L2/meta.txt | 463b2091f413a5a293237ea3ef32e52dcc52fd1997b8d7c07b6abda5745c7190 |
| B_4_L2/H_R_fwd.bin | 9b048b6fc38ec921e9f70a98185cec91763547af7f598f8abe8b3caba6f6b8c5 |
| B_4_L2/y.bin | 68256a4a8e47a9600b1ca21a0dd1106f7a03afa8e035bb118f356598e7332765 |
| B_4_L2/meta.txt | 041740e4f61dd8c6884cc480328cfd63aea4ef869972096c41577174170b0c3a |
| B_6_L2/H_R_fwd.bin | 400b9570a8e8c3921cae43a8bb15e8c32bb531511c3913effcf3e05b55f9aaa8 |
| B_6_L2/y.bin | 89a91b9debbd114b6bc2b80ecf06475804234763c5d02997d09b3777f05359ec |
| B_6_L2/meta.txt | 8d23b0d1c9147e3806c4908d330bf8f24d9358df672306afcf0801640aadfbd7 |
| fd/A_4_L1.dirs | 2448f82d0d61dfe1bf1a9c885cac982dc1c4c578c409e349b88f965ebd5890dd |
| fd/A_4_L2.dirs | 027f49938231decdcdafaa6bb457cca64e9fb5f2354331a231ffcd91d0287790 |
| fd/{A_6_L2,B_4_L2,B_6_L2}.dirs | cb10441fe6b18c691e57c1f85ceaa8235a1672052485b635eac39d5cc3287a4d |

## 5. Conditions for resumption
These are pre-analysis declarations under v4's "departures are declared in REPORT.md". They do not edit PROTOCOL-v4.
Each condition only tightens a gate, and none changes the archived outcome.
- **R1 (F5).** Before `fd_checks_v4.py` runs, it (or the operator) verifies the §4 hashes above, and the result is
  logged in `logs/gate_H2H3_v4.txt`. A mismatch stops the study.
- **R2 (F1).** Declare this reading of NV2(c): *if e_c is smallest at h₁ (no truncation-side pair), the row fails its
  slope gate, as in H2; "not resolvable" requires at least one truncation-side pair, all failing the resolvability
  test.* On the archived data the minimum is at h₂ in all five rows, so the outcome is unchanged.
- **R3 (cross-implementation).** `logs/gate_H2H3_v4.txt` must reproduce the §1 table above to the printed digits:
  e_a, c, e_c, ρ, the resolvability status and the pass/fail. A disagreement is reported and the study stops for review.
  This is the independent check for a script that the builder will write after seeing the data.

**Recommended, not blocking:**
- correct the REPORT §2 wording in line with §1.3;
- relabel builder reading (ii) as lenient (F4);
- add the symmetric NV6 value (§3);
- state the single-direction coverage limit (F3) next to the H3 result.

## 6. Reproduction (from `warp_hypad/study4`; `C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I`)
`sha256sum PROTOCOL*.md hessians/fd/*.out driver/s4.cpp` was compared against `logs/*.sha256` and v4 §9.

My scripts are in the session scratchpad and are not saved here: `vera_v4.py` (recompute and perturbations) and
`vera_v4b.py` (v3 vs v4 vs v4 + R2). Their essential lines:

    v = default_rng(4).standard_normal(108); v /= norm(v)        # == .dirs v, bit-exact
    for each "v4 s R+ R- g+(108) g-(108)" line of hessians/fd/<d>.out (R0 = first line):
      e_a = ||(g+-g-)/(2s) - Hv||/||Hv||;  c = |(R+-R-)/(2s) - v.(g++g-)/2| / |v.(g++g-)/2|
      sd = (R+-2R0+R-)/s^2; pg = v.(g+-g-)/(2s); e_c = |sd-vHv|/|vHv|; ratio = (sd-vHv)/(pg-vHv) at s1
    truncation side = steps up to argmin; slope ok = any consecutive ratio >= 30; best <= 1e-4
    rho(s) = max(e_c[3]*(s3/s)^2, e_c[4]*(s4/s)^2); pair k resolvable iff rho(s_{k+1}) <= e_c[k]/30
    perturb: H + E for E in {±d|vHv| vv^T, ±d||Hv||(vw^T+wv^T), ±d||H||_2 ww^T}, w ⊥ v (rng 12345)

## 7. Candidate lessons (proposed only; for the main session)
- A rule that exempts a check when "no pair is resolvable" must not be satisfied vacuously by an empty pair set. Test
  every exemption rule on a deliberately wrong input, not only on the data that motivated it. (Scope: post-data
  verification-rule changes. Evidence: F1. Limitation: found on one constructed perturbation.)
- When a ruling counts "independent derivations", check whether they share data. Here the projected gradient is a
  component of the (a) check.

**CLEARED TO RESUME UNDER v4**, conditional on R1-R3 being declared in REPORT.md and logged before
`driver/fd_checks_v4.py` runs. If any of R1-R3 fails, the study stops for review again.

## R3 ruling (Vera, 2026-10-27)
I reran my scratch scripts `vera_v4.py` and `vera_v4b.py` on the archived inputs (both exit 0, B_4_L2 still PASSES in both).
I also printed B_4_L2 at full precision (`repr`, same code path as `vera_v4.py`):
- e_a = 1.0362395688548397e-3, 1.0359193768504914e-5, 1.0399859770426354e-7, **3.4145378758701976e-8** (argmin h₄).
- signed pg error = 6.521180800267384e-4, 6.517414727109013e-6, 6.623609867042947e-8, **−1.366498121696414e-10**.

Both values agree with Ada's `logs/gate_H2H3_v4.txt` (3.414538e-08 and −1.366498e-10) to all the digits printed there.
**Neither difference is a computational discrepancy. Both are errors in my §1 table.**
- **(a) e_a best (h₄): double rounding in my transcription.** My script printed `3.415e-08`, and I rounded that again to
  3.42e-8 by hand. The value is 3.4145e-8, which rounds to **3.41e-8**. Ada's value is correct.
- **(b) signed projected-gradient error (h₄): a dropped sign in my transcription.** My script printed `-1.366e-10`, but
  the §1 column, which is labelled "signed", shows 1.4e-10. It should read **−1.4e-10**. Ada's value is correct.
  This is descriptive only and is not a gate.

Erratum to §1, B_4_L2 row: e_a "… 3.41e-8"; signed projected-gradient error "… −1.4e-10". No other cell changes.

**Comparator fix between runs 1 and 2.** `diff --strip-trailing-cr driver/fd_checks_v4_run1.py.txt driver/fd_checks_v4.py`
shows that only two things changed: the `VERA` reference dictionary (cells tagged `best:` and moved to the best-step slot)
and the R3 comparison branch (lines 248–253). Excluding the header and the `R3 ` lines, `logs/gate_H2H3_v4_run1.txt` and
`logs/gate_H2H3_v4.txt` are identical apart from the R3 summary count. **The analysis did not change.**
The run-1 archive hash is a0c567be…, which matches `logs/execution.log`. The mapping is the correct reading of my
table: the elided "first … x" cell is the value at the best step.
Remaining limitation: the archived run-1 copy has LF line endings and the live script has CRLF.

**Reproduction:** `python -I vera_v4.py STUDY4` and `python -I vera_v4b.py STUDY4` (session scratchpad), plus a `repr`
print of `o['ea']` and `o['pg_s']` for B_4_L2.

**R3 SATISFIED (V4 passes; P1 may start under v4).** All 86 items agree with my own recomputation. The 2 mismatches came
from my mis-transcribed reference cells, corrected by the erratum above, not from Ada's implementation.
Candidate lesson (proposed): generate reference tables in code from full-precision values, and never re-round printed
output by hand. Make the comparator honour signs.
