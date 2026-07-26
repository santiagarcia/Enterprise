# B2 — abstract constitutive inverse design

Milestone B2 builds the first inverse-design layer:

    g_eff^target  ->  Z_target^{abcd}  ->  observer blocks  ->  nearest admissible design

**B2 is not a hardware optimizer.** It works entirely in an abstract
constitutive-design space and answers one question:

> Given a target effective metric or target constitutive tensor, what reciprocal,
> finite, bounded electromagnetic constitutive response most closely realizes it
> under stated constraints?

No Maxwell solver, no coils, no antennas, no metamaterial unit cells, no
dispersion fitting, no experimental control model. Those are later milestones.

Everything here rests on B1C, which is complete and tagged
`program-b-b1c-complete`. B2 adds no new physics conventions and reuses the
verified reconstruction path unchanged.

---

## 1. Conventions inherited from B1C

Signature (−,+,+,+); Heaviside–Lorentz with c = 1; Z carries four upper indices;
`H^{ab} = Z^{abcd} F_{cd}`. Levi-Civita:
ε<sub>0123</sub> = +√(−det g_lab), ε<sup>0123</sup> = −1/√(−det g_lab).

Observer blocks, from Schuster & Visser (A7a)–(A7d), with the normalisations
**fixed** at B1C and not reopened here:

    ε^{ab} = −2 W_ε^{ab}          [μ⁻¹]^{ab} = +2 W_μ^{ab}
    ζ^{ab} = +2 W_ζ^{ab}          ζ†^{ab}    = +2 [W_ζ^T]^{ab}

    D^a = ε^{ab} E_b − ζ^{ab} B_b
    H^a = ζ†^{ab} E_b + [μ⁻¹]^{ab} B_b

`g_lab` and `g_eff` remain distinct types. Every observer in B2 is a
**laboratory** observer, normalised with `g_lab`.

---

## 2. Design space

### 2.1 The observer frame

A design is specified relative to **one fixed laboratory observer**. The blocks
must be spatial with respect to it, so the parameterization is written in an
orthonormal spatial triad e<sub>(i)</sub><sup>a</sup> adapted to V:

    g_lab(e_(i), e_(j)) = δ_ij,    g_lab(e_(i), V) = 0

    ε^{ab} = Σ_ij ε_ij e_(i)^a e_(j)^b        (and likewise μ⁻¹, ζ)

Spatiality and symmetry are then **structural**: no projection or symmetrisation
is applied after the fact, so no residual can appear.

The triad is a **choice of basis**, built by real-part Gram–Schmidt from `g_lab`
and V. That is legitimate here in a way it was not in the quarantined legacy
decomposition: there a basis choice was used to *extract physics from a tensor*;
here it only supplies *coordinates on the design space*. The set of tensors
reachable is independent of the triad, and the triad depends only on the fixed
inputs (g_lab, V), never on the design variables, so no derivative channel flows
through it.

For the analogue configuration g_lab = η, V = (1,0,0,0), the triad is the
coordinate triad and ε<sub>ij</sub> = ε<sup>(i+1)(j+1)</sup> exactly.

### 2.2 The unrestricted reciprocal family — 21 degrees of freedom

| block | entries | count |
|---|---|---|
| ε<sub>ij</sub> | symmetric 3×3 | 6 |
| [μ⁻¹]<sub>ij</sub> | symmetric 3×3 | 6 |
| ζ<sub>ij</sub> | general 3×3 | 9 |
| | **total** | **21** |

with ζ† = ζ<sup>T</sup> imposed structurally. 21 is exactly the degree-of-freedom
count of a Z carrying Schuster & Visser's symmetries (12) and (13), so this
family is the **whole** reciprocal design space — nothing is reduced prematurely
to diagonal or isotropic form.

Parameter ordering, fixed and machine-checked by the round-trip case:

    q[0..5]   eps   (11, 12, 13, 22, 23, 33)
    q[6..11]  muinv (11, 12, 13, 22, 23, 33)
    q[12..20] zeta  (11, 12, 13, 21, 22, 23, 31, 32, 33)   row-major

### 2.3 Restricted families, for verification and for stated-constraint design

| family | parameters | what it fixes |
|---|---|---|
| `UnrestrictedReciprocal` | 21 | the full reciprocal space |
| `PositiveDefinite` | 21 | ε = LLᵀ + λ_min 𝟙, μ⁻¹ likewise; ζ free |
| `ZeroMagnetoelectric` | 12 | ζ ≡ 0 |
| `DiagonalAnisotropic` | 6 | ε, μ⁻¹ diagonal; ζ ≡ 0 |
| `Isotropic` | 2 | ε = ε₀𝟙, μ⁻¹ = ν₀𝟙, ζ ≡ 0 |

**The positive-definite family is a restricted physical-design family, not the
general search space.** ε = LLᵀ + λ_min𝟙 with L lower triangular guarantees
ε ≻ 0, which **excludes**:

* indefinite media;
* negative-index responses;
* ergoregion-crossing designs.

The last exclusion is concrete rather than hypothetical. B1C established that for
the Alcubierre analogue with the laboratory observer,
spec([μ⁻¹]) = {1, −g₀₀, −g₀₀}. Past the ergosurface −g₀₀ = 1 − v_s²f² turns
negative, so **[μ⁻¹] is genuinely indefinite there** and no positive-definite
parameterization can represent it. The unrestricted family therefore remains
available and is the default.

A `MetricInduced` design family is **deliberately not implemented** in B2. A
family parameterized by g_eff has 9 conformal degrees of freedom against the
design space's 21, so fitting inside it is the consistency-condition /
no-birefringence question, and that belongs to a later milestone rather than
being smuggled in as an optimizer restriction. Metric-induced targets *are*
supported — see §3 — which is what verification problem 4 needs.

---

## 3. Targets

Two entry points, both recording a full classification.

**Effective-metric target.** From `EffectiveMetric<Scalar> g_target`,
`LaboratoryMetric<Scalar> g_lab` and a laboratory observer, using the verified
covariant map of Schuster & Visser eq (6).

**Covariant constitutive target.** Direct input of Z<sub>target</sub><sup>abcd</sup>,
which permits manufactured tests and non-metric targets.

Every target records whether it is: metric-induced; reciprocal
(Z<sup>abcd</sup> = Z<sup>cdab</sup>); pair-antisymmetric
(Z<sup>(ab)cd</sup> = Z<sup>ab(cd)</sup> = 0); finite; observer-decomposable; and
inside the supported real, lossless, non-dispersive model. A target failing
reciprocity is **reported, not silently fitted** — the design space is reciprocal
by construction, so the best achievable residual is then bounded below by the
target's own skewon part.

---

## 4. Reconstruction path

For each design vector q:

1. build the observer blocks from q and the triad;
2. reconstruct Z<sub>design</sub><sup>abcd</sup> with the **verified A11
   implementation**, `reconstruct_Z_A11`;
3. evaluate symmetries and spatiality as diagnostics;
4. compare against the target.

There is **no second production reconstruction formula.** The optimizer drives
the same code path B1C verified to 2.2e-16 on random reciprocal media.

---

## 5. Objectives

Reported separately, never silently summed into one number.

### 5.1 Two norms, and which is which

For a rank-4 difference A<sup>abcd</sup>:

**Coordinate-basis norm** — plainly labelled as such, because it is *not*
invariant:

    ‖A‖²_coord = Σ_{abcd} (A^{abcd})²

**Observer-referred norm** — contracted, and a genuine scalar:

    ‖A‖²_V = A^{abcd} A^{efgh} G_ae G_bf G_cg G_dh,    G_ab = g_lab_{ab} + 2 V_a V_b

G is used rather than g_lab itself for a mathematical reason: contracting with a
Lorentzian metric is **indefinite**, so ‖A‖²_g can vanish or go negative for
A ≠ 0 and is unusable as an objective. G = h + V⊗V is positive definite — it is
g_lab with the observer's timelike direction flipped — so ‖A‖²_V is a positive
definite quadratic form and a scalar under coordinate change, given that V
transforms with the coordinates.

The two coincide exactly when g_lab = η and V = (1,0,0,0), because then
G = diag(1,1,1,1). That is the analogue configuration, which is why the
coordinate-basis norm is usable there and why calling it invariant would still be
wrong. The gate asserts both facts: agreement in that configuration, and that
the observer-referred norm is invariant under a nontrivial linear coordinate
change while the coordinate-basis norm is not.

### 5.2 Terms

    J_Z          = ‖Z_design − Z_target‖² / max(‖Z_target‖², δ_Z)

    J_blocks     = w_ε J_ε + w_μ J_μ + w_ζ J_ζ,
                   J_X = ‖X_design − X_target‖² / max(‖X_target‖², δ_X)

    J_reciprocity = ‖Z^{abcd} − Z^{cdab}‖²

`J_Z` and `J_blocks` are kept apart deliberately: covariant tensor agreement and
observer-dependent block agreement answer different questions, and B1C's
withdrawn Plebański case is what happens when they are conflated.

`J_reciprocity` is a **diagnostic**, expected to be zero to roundoff for this
structurally reciprocal parameterization. It is not a useful optimizer term and
is not included in the total by default.

    J_total = J_primary + J_penalty

with `J_primary` selected as `J_Z` or `J_blocks`.

### 5.3 Constraints, and why they are smooth

Configurable bounds on entry magnitude, Frobenius norm, determinant floor, and
magnetoelectric norm, per block.

**No branch is taken on a differentiable scalar inside the objective.** A
one-sided penalty uses the smooth hinge

    hinge(x) = ½ ( x + √(x² + τ²) ) − τ/2

which is C<sup>∞</sup>, needs only `sqrt`, and requires no comparison on a Scalar
value. Penalties are `hinge(...)²`.

Eigenvalues and condition numbers are **reported as diagnostics on real parts**,
not used as penalties, because computing them needs an iterative algorithm with
comparisons. Where a bound on them is wanted, the AD-safe surrogates are used
and labelled as such: the Frobenius norm bounds the spectral radius from above,
and the closed-form 3×3 determinant bounds the product of eigenvalues from below.
An upper bound on the condition number is available as
‖A‖_F ‖A⁻¹‖_F, computed from the adjugate, and is reported as **a bound on**
cond, never as cond.

---

## 6. Gradients

The design variables are seeded through the complete chain

    q  ->  (ε, μ⁻¹, ζ)  ->  Z_design  ->  J

There is **no separate OTI implementation of the objective**; `evaluate_design`
is scalar-generic and the same instantiation serves `double`, the test scalar,
and the HyPAD scalar.

### 6.1 Batching is not optional

The available OTI scalar is `onumm3n2`: **3 derivative directions**, order 2.
It cannot seed 21 variables at once, and pretending otherwise would be a lie
about the backend. The interface is therefore explicitly batched:

```cpp
GradientBatchResult evaluate_gradient_batch(variables, target, constraints,
                                            active_parameter_indices);
```

with at most `kOtiDirections` active indices per evaluation, and the full
gradient assembled across ⌈n/3⌉ = **7 batches** for the 21-parameter family.
Reported: number of evaluations, active indices per batch, the assembled
gradient, worst analytic error, worst finite-difference error.

### 6.2 Verification

1. **Hand-derived**, on a manufactured quadratic J = Σ c_k(q_k − a_k)² whose
   gradient 2c_k(q_k − a_k) is known in closed form. This tests the seeding and
   batch-assembly machinery in isolation from the objective.
2. **Centred finite differences** over multiple step sizes.
3. **Simultaneous versus separate seeding** — the same partial derivative must
   come out identically whether a parameter was seeded alone or alongside two
   others.

One useful exactness property is exploited and asserted: A11 is **homogeneous
linear** in the blocks and the blocks are linear in q, so for the unrestricted
family `J_Z` is an exact **quadratic** in q. Centred differences of a quadratic
have zero truncation error, so the finite-difference gradient is exact to
roundoff and the comparison can be tight rather than step-size limited.

---

## 7. Optimizer

Gradient descent with backtracking Armijo line search — deterministic, local,
and dependency-free. No external optimization library.

Recorded per run: reproducible initialization, objective history, gradient-norm
history, line-search status per iteration, termination reason, a
maximum-iteration guard, rejection of non-finite objectives or gradients, and the
final design parameters together with the reconstructed tensor.

Termination reasons are an enum, not a string comparison:
`GradientNormBelowTolerance`, `ObjectiveBelowTolerance`,
`StepBelowTolerance`, `LineSearchFailed`, `MaxIterationsReached`,
`NonFiniteEncountered`.

---

## 8. Test plan

| # | problem | what it proves |
|---|---|---|
| 1 | isotropic exact recovery | the chain closes on a known answer |
| 2 | anisotropic reciprocal exact recovery | distinct diagonal entries plus an off-diagonal symmetric entry |
| 3 | magnetoelectric exact recovery | nonsymmetric ζ with ζ† = ζᵀ; a zero-ζ solution cannot pass |
| 4 | metric-induced target | subluminal Alcubierre with g_lab = η, recovered in the **unrestricted** family, compared against Z_target and **not** against another observer's Plebański blocks |
| 5 | restricted-family nearest approximation | ζ ≡ 0 family against a magnetoelectric target: residual must stay nonzero and must match the independently computed least-squares projection |
| 6 | infeasible bounded target | bounds exclude the target; the result must be reported as a best bounded approximation, not as exact realization |
| 7 | coordinate-change consistency | the **final tensor relationship** under a nontrivial linear change, not identical iteration histories — a coordinate-basis norm need not preserve the trajectory |

Gate `B2-ABSTRACT-INVERSE`, required categories: parameter-vector round trip;
symmetry by construction; A11 reconstruction; objective zero at the exact target;
objective positive away from it; gradient against analytic; gradient against
finite differences; batched OTI gradient; optimizer descent; the five recoveries
above; restricted-family projection; infeasible bounded result; coordinate-change
consistency; deterministic repeatability; non-finite handling.

Acceptance: **0 failed, 0 unresolved mathematical skips.** With OTIlib disabled
only the OTI gradient cases may skip; every optimization and double-precision
inverse-design case must still run.

---

## 8a. Measured results

Gate `B2-ABSTRACT-INVERSE`: **33/33 passed, 0 failed, 0 skipped** with OTIlib
enabled; **28 passed, 5 skipped** with it disabled — the five being HyPAD
gradient cases only, and the totals reconciling to the same 33 so nothing
vanishes silently.

| quantity | measured |
|---|---|
| parameter-vector round trip | 0.0 |
| symmetry residual (structural) | 0.0 |
| spatiality residual (structural) | 0.0 |
| objective at the exact target | 7.9e-32 |
| (A11) round trip, design blocks → Z → blocks | 2.2e-16 |
| HyPAD gradient vs hand-derived (manufactured) | **0.0** |
| HyPAD gradient vs closed form (J_Z) | 3.5e-18 |
| HyPAD vs centred differences | 1.3e-12 |
| batch independence (1 seeded vs 3) | 0.0 |
| isotropic recovery, relative ‖Z−Z_t‖ | 9.5e-11 |
| anisotropic recovery | 7.0e-11 |
| magnetoelectric recovery | 9.3e-11 |
| metric-induced recovery (Alcubierre v_s = 0.4) | 1.4e-10 |
| coordinate-change consistency | 1.1e-10 |
| restricted family vs normal-equation projection | 7.1e-11 |

The 21-parameter gradient assembles from exactly **7 HyPAD batches** of 3 active
directions: `[0,1,2][3,4,5][6,7,8][9,10,11][12,13,14][15,16,17][18,19,20]`.

### The cost of restricting the family

Least-squares `J_Z` against the same magnetoelectric target, from the normal
equations rather than the optimizer:

| family | parameters | J_Z |
|---|---|---|
| `UnrestrictedReciprocal` | 21 | 6.5e-34 |
| `ZeroMagnetoelectric` | 12 | 7.5e-2 |
| `DiagonalAnisotropic` | 6 | 8.7e-2 |
| `Isotropic` | 2 | 1.4e-1 |

The unrestricted family reaches zero, as it must: 21 parameters against a
21-dimensional reciprocal design space. Every restriction leaves a residual that
is a property of the family, not an optimizer failure — which is the point of
problem 5.

### Two honest caveats

**The bound penalty is soft.** A smooth hinge buys differentiability at the price
of admitting a small violation, so in problem 6 `max|ε|` settles at 1.5001
against a bound of 1.5. The gate allows that explicitly rather than pretending
the constraint is hard. A genuinely hard bound needs a projection or a
structural parameterization, neither of which is in B2.

**Problem 6 terminates on `MaxIterationsReached`.** That is the expected outcome
for a penalized infeasible problem — there is no stationary point where the
penalty gradient and the objective gradient balance to below tolerance at this
weight — and the case asserts what matters (a strictly positive `J_Z`, the bound
reported active) rather than asserting convergence that did not happen.

---

## 9. What B2 does not claim

* Nothing about physical realizability of any computed design. No dispersion,
  loss, passivity or bandwidth model exists.
* Nothing about hardware. There is no geometry, no unit cell, no fabrication
  constraint.
* Nothing about measurable gravitational curvature. The laboratory target is an
  effective electromagnetic metric.
* A recovered design is a point in an abstract constitutive space that
  reproduces a target tensor. Whether any material realizes it is not addressed.
