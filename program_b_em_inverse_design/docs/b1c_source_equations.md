# B1C: primary-source equations, retrieved

Retrieved 2026-07-25 from the full text of

> **S. Schuster and M. Visser, "Effective metrics and a fully covariant
> description of constitutive tensors in electrodynamics", *Physical Review D*
> **96**, 124019 (2017). DOI 10.1103/PhysRevD.96.124019. arXiv:1706.06280v2
> [gr-qc], 20 Dec 2017. Victoria University of Wellington.**

Full extracted text is committed at
[`sources/schuster_visser_prd96_124019_extracted.txt`](sources/schuster_visser_prd96_124019_extracted.txt)
so every equation below can be checked against the source without re-fetching.

This supersedes the "recorded from memory" status for this reference.

---

## Conventions, stated in the paper's NOTATION section

> "This paper follows the sign conventions of [3] and [4]. Specifically, our
> metrics will have signature (−+++). Symmetrisation and antisymmetrisation on
> indices is indicated by enclosing these indices in round or square brackets,
> respectively. **Raising and lowering of indices shall always be done employing
> the physical background metric g.**"

| | Paper | This project | Match |
|---|---|---|---|
| Signature | (−,+,+,+) | (−,+,+,+) | ✅ |
| Index raising/lowering | physical background metric g | g_lab | ✅ |
| Spatial projector | h^ab = g^ab + V^aV^b | same | ✅ |
| Temporal projector | t^ab = −V^aV^b | — | (not yet used) |

Projector definition, from the text near eq (41):

> "we can use the temporal and spatial projection operators, respectively
> t^{ab} = −V^aV_b and h^{ab} = g^{ab} + V^aV^b, to write any vector as
> X^a = δ^a_b X^b = t^a_b X^b + h^a_b X^b. (41)"

This matches the projector the correction task specifies.

---

## The Bel extraction formulas — equations (A7a)–(A7d)

**This is the piece B1C was missing.**

    W^{bd}_ε      = V_a V_c Z^{abcd}           (A7a)
    W^{bd}_μ      = V_a V_c (*Z*)^{abcd}       (A7b)
    [W^T_ζ]^{bd}  = V_a V_c (*Z)^{abcd}        (A7c)
    W^{bd}_ζ      = V_a V_c (Z*)^{abcd}        (A7d)

where `*` is the Hodge dual taken on the **first** index pair (`*Z`), the
**second** index pair (`Z*`), or **both** (`*Z*`).

### Two consequences for the current implementation

1. **The ε formula in `covariant_constitutive.hpp` is CORRECT**, and the
   factor of −2 is right. Verified against the paper's own vacuum case:
   equation (B1a) gives ε^{ab} = ε·h^{ab}, so vacuum ε = 1 must yield h^{ab};
   substituting the vacuum Z into (A7a) gives W_ε = −½h^{bd}, hence
   ε^{bd} = −2 W_ε^{bd} = −2 V_aV_c Z^{abcd}, which in a V-adapted frame is
   exactly the implemented `ε^{ij} = −2 Z^{0i0j}`.

   Note W_ε is the **Bel piece**, not ε itself. The paper is explicit that the
   Bel names are misleading here:

   > "the names given to the three independent matrices W_ε, W_ζ, and W_μ
   > encountered in this decomposition in the GR community are very misleading
   > in the present context … what goes under the name of 'electric tensors' in
   > [39] corresponds to **both** the permittivity and the permeability tensors,
   > while the 'magnetic tensors' here are the magneto-electric tensor and its
   > transpose."

2. **The μ⁻¹ and ζ formulas in the implementation are NOT the paper's.** The
   paper obtains them by Hodge-dualising Z before contracting with V; the
   implementation instead contracts Levi-Civita symbols with frame components of
   Z. For μ⁻¹ the two happen to agree (verified numerically to 2.2e-16 against
   the inverse of the Plebański μ), but **the ζ normalisation was guessed and
   is not the paper's** — which is exactly why the magnetoelectric sign case
   still fails. It must be rebuilt from (A7c)/(A7d).

## Z in terms of the blocks — equations (A8) and (A11)

Equation (A8), transcribed from the extracted text:

    Z^{abcd} = ½ [ V^aV^d ε^{bc} + V^bV^c ε^{ad} − V^bV^d ε^{ac} − V^aV^c ε^{bd} ]
             + ⅛ ε^{abef} ε^{cdgh} [ V_f (μ⁻¹)_{eg} V_h + V_e (μ⁻¹)_{fh} V_g
                                     − V_e (μ⁻¹)_{fg} V_h − V_f (μ⁻¹)_{eh} V_g ]
             + ¼ ε^{abef} [ ζ_{fc} V^d V_e + ζ_{ed} V^c V_f
                            − ζ_{ec} V^d V_f − ζ_{fd} V^c V_e ]
             + ¼ ε^{cdgh} [ (ζ†)^{bg} V^a V^h + (ζ†)^{ah} V^b V^g
                            − (ζ†)^{ag} V^b V^h − (ζ†)^{bh} V^a V^g ]

Equation (A11) is an equivalent form written with the projector h:

    Z^{abcd} = ½ [ V^dV^a ε^{bc} − V^cV^a ε^{bd} + V^cV^b ε^{ad} − V^dV^b ε^{ac}
                 + h^{ad}(μ⁻¹)^{cb} − h^{ac}(μ⁻¹)^{db} + h^{bc}(μ⁻¹)^{ad}
                 − h^{bd}(μ⁻¹)^{ac} + (h^{bd}h^{ac} − h^{bc}h^{ad})(μ⁻¹)^e_e
                 + ε^{fabe}(V^d ζ_e{}^c − V^c ζ_e{}^d) V_f
                 + ε^{fcde}(V^b ζ_e{}^a − V^a ζ_e{}^b) V_f ]

**(A11) is the equation gate case 20 needs** — reconstruction of Z from the
decomposed blocks. It also fixes the ζ sign and normalisation by construction,
which is how case 9 should be closed: derive ζ from (A7c)/(A7d), substitute into
(A11), and require the original Z back.

## Consistency conditions — equation (7)

    ε^{ab} = μ^{ab},    ζ^{ab} = 0                                    (7)

i.e. the metricity conditions **as seen by a particular observer**. This is the
covariant statement of what `classify()` checks on the 3+1 blocks.

## Effective metric reconstruction — equations (8a), (8b), (54a)–(54c)

    (g_eff)_ab = − sqrt( −det(g) / sqrt(det ε) ) V_a V_b
                 + sqrt( sqrt(det ε) / −det(g) ) [ε]^#_ab            (8a)

with `[ε]^#` the **Moore–Penrose pseudo-inverse** and `det ε` the associated
**pseudo-determinant** — necessary because ε is spatial and therefore
singular as a 4×4 object. Equation (8b) is the same with μ, and (54c) the same
with μ⁻¹ directly.

The current `medium_to_metric()` uses a 3×3 inverse in an ad-hoc g₀₀ = −1 gauge.
It gives the right answer for the diagonal w = 0 case it is tested on, but
**(8a) is the general form** and is what the w ≠ 0 inverse should be built from.

## Moving isotropic media — equations (B1)–(B3), (B16), (B17)

A ready-made verification benchmark with a **known nonzero magnetoelectric
coupling**, which is precisely what correction task §4 asks for:

    ε^{ab}      = ε (g^{ab} + V^aV^b) = ε h^{ab}                      (B1a)
    (μ⁻¹)^{ab}  = μ⁻¹ (g^{ab} + V^aV^b) = μ⁻¹ h^{ab}                  (B1b)

    Z^{abcd} = −(ε/2)( V^aV^c h^{bd} + V^bV^d h^{ac}
                       − V^aV^d h^{bc} − V^bV^c h^{ad} )
               + (μ⁻¹/2)( h^{ac}h^{bd} − h^{ad}h^{bc} )               (B2)

    Z^{abcd} = (μ⁻¹/2)[ (h^{ac} − εμ V^aV^c)(h^{bd} − εμ V^bV^d)
                        − (h^{ad} − εμ V^aV^d)(h^{bc} − εμ V^bV^c) ]  (B3)

For a medium moving with four-velocity W relative to the observer:

    [μ⁻¹_W]^{bd} = h_W^{bd} μ + (μ⁻¹ − ε)[ (V·W)² h_W^{bd}
                                            − h_W^{be} h_{ef} h_W^{fd} ]   (B16)

    ζ_W^{ac} = (ε − μ⁻¹)(V·W) ε^{acef} W_e V_f                        (B17)

**(B17) is the magnetoelectric benchmark.** It is nonzero exactly when the
medium moves relative to the observer (V·W ≠ 0) **and** the medium is not
impedance matched (ε ≠ μ⁻¹) — and it gives the sign explicitly. This closes
correction task §4 without any `min(error_to_xi, error_to_zeta)` construction.

---

## What this changes in the plan

| B1C case | Status | Action |
|---|---|---|
| 4, 5 — projector and spatiality | **done**, B1C-S23 42/42 | h^a_b built; V_a ε^{ab} = 0 exactly |
| 6 — full F → H reconstruction | **done**, B1C-S456 | direct and block paths agree to 2.7e-15 over 5 media × 12 two-forms |
| 7 — same-observer Plebański | **withdrawn**, aggregate B1C 53/53 | the equality was between two different observers' blocks; ε_lab = 𝟙 while ε_P = diag(1, 1/(1−v²f²), …), differing by exactly the observer normalisation |
| 9 — magnetoelectric sign | **resolved** | ζ = +2W_ζ against (B17); handedness asserted in one direction |
| 12 — divergence diagnostic | **withdrawn**, aggregate B1C 53/53 | cond(ε) ≡ 1 because ε = h identically; what diverges is the *permeability*, min\|eig μ⁻¹\| = min(1, \|g₀₀\|) → 0 |
| 20 — reconstruct Z from blocks | **done**, B1C-S456 | (A11), round trip 2.2e-16 including random reciprocal media |

Cases 7 and 12 are now closed as well. Both were tests written wrongly rather
than physics, and both were withdrawn and replaced with assertions of what is
actually true — see
[b1c_aggregate_gate_audit.md](b1c_aggregate_gate_audit.md) for the derivation of
each. The whole repository now passes: B1 20/20, B1C-S23 42/42, B1C-S456 36/36,
aggregate B1C 53/53.

## Still outstanding for the reference list

Primary sources verified so far: **Plebański (1960)** (APS record) and
**Schuster & Visser (2017)** (full text, this document). The rest of the
bibliography in [effective_metric_theory.md](effective_metric_theory.md) §11
remains "recorded from memory" and must be retrieved before any preprint —
in particular a primary source for the Plebański 3+1 equations that states
**which observer or foliation** its three-dimensional tensors belong to, since
case 7 depends on that and it must not be assumed.

---

# Derivation of the Bel → physical-block normalisation (B1C step 2)

Not inferred from the vacuum case. Derived by contracting (A7a) into (A8), and
independently confirmed against the paper's own worked example (B2), which
carries **arbitrary** ε and μ⁻¹ and therefore pins a factor that vacuum cannot.

## ε from A7a

Apply `W_ε^{bd} = V_a V_c Z^{abcd}` to (A8), using `V·V = −1`, the spatiality of
every block (`V_a ε^{ab} = 0`), and the total antisymmetry of ε^{abcd}:

**ε-term of (A8):** `½[V^aV^d ε^{bc} + V^bV^c ε^{ad} − V^bV^d ε^{ac} − V^aV^c ε^{bd}]`

| sub-term | contraction | value |
|---|---|---|
| `V_aV_c V^aV^d ε^{bc}` | `V_c ε^{bc} = 0` (spatial) | 0 |
| `V_aV_c V^bV^c ε^{ad}` | `V_a ε^{ad} = 0` (spatial) | 0 |
| `−V_aV_c V^bV^d ε^{ac}` | `V_aV_c ε^{ac} = 0` (spatial) | 0 |
| `−V_aV_c V^aV^c ε^{bd}` | `(V·V)² = +1` | `−ε^{bd}` |

giving `−½ ε^{bd}`.

**μ⁻¹-term:** every sub-term carries `ε^{abef}V_aV_f` or `ε^{abef}V_aV_e` — a
totally antisymmetric symbol contracted twice with the *same* vector — hence 0.

**ζ and ζ† terms:** likewise each contains `ε^{abef}V_aV_e`, `ε^{abef}V_aV_f`,
`V_a ζ^{ab} = 0`, or `ε^{cdgh}V_cV^h` (antisymmetric on a symmetric pair) — all 0.

Therefore

    W_ε^{bd} = −½ ε^{bd}        ⟹        ε^{bd} = −2 V_a V_c Z^{abcd}      (D1)

### Independent confirmation from (B2), with arbitrary ε

(B2) gives Z for an isotropic medium at rest with **independent** ε and μ⁻¹:

    Z^{abcd} = −(ε/2)(V^aV^c h^{bd} + V^bV^d h^{ac} − V^aV^d h^{bc} − V^bV^c h^{ad})
               + (μ⁻¹/2)(h^{ac}h^{bd} − h^{ad}h^{bc})

Contracting with `V_aV_c`: only the first sub-term survives (`V_a h^{ac} = 0`
kills the rest), giving `W_ε = −(ε/2) h^{bd}`. With (B1a) `ε^{bd} = ε h^{bd}`,
this is `W_ε = −½ ε^{bd}` — **(D1) confirmed for arbitrary ε**, not merely at
ε = 1.

## μ⁻¹ from A7b

For vacuum, `(*Z*) = −Z` exactly (shown by `ε^{abef}ε_{ef}{}^{cd} =
−2(g^{ac}g^{bd} − g^{ad}g^{bc})`, the identity verified as
`b1c.dual.contraction_identity`). Hence `W_μ = −W_ε`, and since vacuum has
ε = μ⁻¹ = 1:

    (μ⁻¹)^{bd} = +2 V_a V_c (*Z*)^{abcd}                                   (D2)

**Note the sign differs from (D1).** This is checked against (B2) with
**ε ≠ μ⁻¹**, which is what distinguishes the factor and the sign — a test with
ε = μ⁻¹ could not.

## ζ, ζ† from A7c/A7d — NORMALISATION RESOLVED

    [ζ†]^{bd} = c_ζ · V_a V_c (*Z)^{abcd}                                  (D3)
    ζ^{bd}    = c_ζ · V_a V_c (Z*)^{abcd}                                  (D4)

    c_ζ = +2

The paper states this directly in equation (A6b), `2W_ζ = ζ`, and it was
independently **measured** against the (B17) moving-medium benchmark before that
line was located: across five boost configurations (three axes plus a diagonal,
four ε/μ⁻¹ pairs) the (A7d) block and the (B17) closed form are proportional
componentwise with residual exactly 0, and the single ratio is +2 to within
6.7e-16. The provisional value was **−2**, taken by symmetry with (D1), and it
was wrong: the pattern is +2 for both μ⁻¹ (D2) and ζ (D4), with ε (D1) the sole
−2.

`c_ζ` for **ζ†** is likewise +2, and that is no longer carried by symmetry with
(A6b) either. It is measured at B1C-S456 against the transpose of (B17) —
(B17) is antisymmetric, so the prediction ζ† = −ζ has a definite sign — as case
`b1c.s5.zeta_adjoint_matches_analytic`, residual 1.1e-16.

Nothing here is a `min(error_to_ξ, error_to_ζ)` construction: each case asserts
one signed value and additionally records how far the opposite sign would be.

## ζ† ≠ ζ^T in general

`ζ` and `ζ†` are stored as **separate blocks**. They are related by transposition
exactly when the pair symmetry (13) `Z^{abcd} = Z^{cdab}` holds, which is a
property of the medium — a skewon-carrying response breaks it. Under (13) the
identity is provable rather than conventional:

    [W^T_ζ]^{bd} = V_a V_c (*Z)^{abcd}
                 = V_a V_c (Z*)^{cdab}        [pair symmetry]
                 = W_ζ^{db}                   [relabel a ↔ c]

The generic code does not assume it; `reciprocity_residual()` measures it and
every (A11) call site in the gate records it.

---

# Step 4: the 3+1 field relations, DERIVED

The paper's equation (15), `D = εE + ζB`, `H = ζ†E + μ⁻¹B`, belongs to its
three-dimensional matrices of equation (20). Those are **not** in the same index
convention as its own covariant blocks (A6)/(A7). Rather than guess which, the
relations were derived from scratch, from (A1) plus `H^{ab} = Z^{abcd}F_{cd}`
plus the (A7) definitions — see the header comment of
[`field_reconstruction.hpp`](../effective_metric/include/program_b/effective_metric/field_reconstruction.hpp).
The result is

    D^a = ε^{ab} E_b − ζ^{ab} B_b                                        (S4-D)
    H^a = [ζ†]^{ab} E_b + [μ⁻¹]^{ab} B_b                                 (S4-H)

**Note the minus in (S4-D)** relative to (15). It is forced, not chosen: with
ε^{0123} = −1/√(−det g) fixed by raising four indices in (−+++), and ζ = +2W_ζ
measured against (B17), no freedom remains. D is read off the excitation
directly while H is read off its dual, and that single extra minus sign is the
whole difference; the same asymmetry is `D = +∂L/∂E`, `H = −∂L/∂B` for
`L = ½[E·εE − 2E·ζB − B·μ⁻¹B]`, whose quadratic form is symmetric, i.e.
`Z^{abcd} = Z^{cdab}`, as required.

Checked against a source **outside** the paper: Minkowski's moving-medium
relations `D + v×H = ε(E + v×B)` give, to first order, `D = εE + (ε−μ⁻¹)(v×B)`.
At β = 1e-3 the computed D_y agrees to 1.0e-6 relative
(`b1c.s5.handedness_matches_minkowski_relations`). That is what rules out a
self-consistently mirrored convention chain.

# Step 6: which magnetoelectric block goes into (A11)

(A11) writes its magnetoelectric terms with `ζ_e{}^c`. In this project's
conventions that slot is filled by **ζ†**, not ζ — the same single transposition
that explains the (S4-D) sign. Determined, not chosen, by a test no sign
convention could pass; measured `max|Ẑ − Z|`:

| medium | ζ(+) | ζ(−) | ζ†(+) | ζ†(−) |
|---|---|---|---|---|
| generic electric block | 0 | 0 | 0 | 0 |
| generic magnetic block | 1.1e-16 | 1.1e-16 | 1.1e-16 | 1.1e-16 |
| ζ diagonal | 0 | 1.0 | 0 | 1.0 |
| ζ purely off-diagonal | 5.0e-1 | 5.0e-1 | **0** | 1.0 |
| moving isotropic, boost x | 1.1 | 2.8e-17 | **2.8e-17** | 1.1 |
| random reciprocal 6×6 | 1.8 | 2.4 | **2.2e-16** | 3.0 |

Only the ζ†(+) column is exact everywhere. The rows that discriminate are the
ones whose ζ is neither symmetric nor antisymmetric: for a diagonal ζ the
transpose is the identity, and for the moving-medium ζ it is a sign flip, so
either alone would license the wrong reading — a sign constant of −1 reproduces
the moving medium exactly and is still wrong. Recorded as
`b1c.a11.magnetoelectric_slot`.
