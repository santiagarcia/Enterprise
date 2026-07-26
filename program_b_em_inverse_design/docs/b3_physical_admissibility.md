# B3 — physical admissibility

Milestone B3 turns the abstract constitutive design of B2 into a
frequency-dependent physical model:

    q_physical  ->  eps(w), mu^-1(w), zeta(w)  ->  Z(w)  ->  J_band

## Scope

> **A passing result is a model that is CAUSAL GLOBALLY and PASSIVE OVER THE
> VERIFIED FREQUENCY BAND, and that approximates the target response over that
> band. It is not yet a fabricated medium or an experimental device.**

The two adjectives have different scopes, and conflating them would overstate the
result:

| property | scope | how it is established |
|---|---|---|
| **causal** | **global**, all w | structural: every pole is in the lower half-plane for every parameter value, because damping passes through a strictly positive transform. No sampling involved. |
| **passive** | **the tested band only** | `min eig L(w) >= 0` evaluated on a finite set of sampled frequencies. A violation between samples, or outside the band, is NOT excluded by this gate. |

**Global passivity is claimed only where it is proved analytically**, and for the
model families here it is not. `L(w)` is a rational function of `w` whose
eigenvalues can dip negative in a narrow interval between grid points --
particularly near a resonance, near a zero of `det mu`, or where two eigenvalue
branches cross. The sampled gate cannot see that.

What partially closes the gap is the **adaptive passivity search** of
`adaptive_scan.hpp`, which refines intervals around local minima of
`lambda_min(L)`, rapid eigenvalue changes, resonances, zeros of `det mu`, and
poles. It is part of **B4 validation** and does not replace the B3 sampled gate;
it narrows the unverified region rather than eliminating it. A genuine global
proof would need an analytic argument -- a positive-real / Herglotz
representation of the whole constitutive matrix -- which no milestone here has
carried out.

No unit cells, no antennas, no coils, no fabrication geometry, no full-wave
Maxwell simulation. B3 is still not a geometry or hardware optimizer.

Gate `B3-PHYSICAL-ADMISSIBILITY`: **38/38 passed, 0 failed, 0 skipped** with
OTIlib enabled; **36 passed, 2 skipped** with it disabled — the two being HyPAD
gradient cases, totals reconciling to the same 38.

---

## 1. The Fourier convention, and everything that follows from it

    f(t) = Re[ f-hat(w) e^{-i w t} ]

**Fixed once, in `dispersion.hpp`, and nothing may test passivity before it.**

| quantity | value in this convention |
|---|---|
| time derivative | `d/dt -> -i w` |
| constitutive relation | `D = eps E - zeta B`, `H = zeta^T E + [mu^-1] B` (from B1C) |
| reciprocity | `zeta_adjoint = zeta^T`, **unconjugated** transpose |
| dissipated power | `<p> = (w/2) Im[E* . D + H* . B] = (w/2) v† L v` |
| Poynting vector | `<S> = (1/2) Re[E x H*]` |
| sign of Im eps | **passive ⇒ Im eps ≥ 0** for w > 0 |
| causal half-plane | analytic in the **upper** half w-plane; poles in the **lower** |
| reality | `eps(-w) = conj(eps(w))` for real w |

**How the convention is pinned rather than assumed.** The model forms this
milestone specifies — `F/(w0² − w² − i γ w)` and `−wp²/(w² + i γ w)` — both have
strictly positive imaginary part for positive parameters, which is passive. Under
`e^{+iwt}` both would be **gain**. That is the check, and it is a gate case
(`b3.convention.sign_of_loss_is_pinned`).

---

## 2. Complex scalar support

`Complex<Scalar>` is a paired-real type, **not** `std::complex`. Two reasons:
`std::complex<T>` is undefined for `T` outside the built-in floating types, and
its division and `abs` branch on magnitudes — which a hypercomplex AD scalar
cannot tolerate. Written from the definitions, with no comparison on a `Scalar`.

**No assumption is made that OTIlib supplies a complex AD scalar.** Pairing two
real AD scalars gives complex arithmetic whose derivative channels flow through
real and imaginary parts independently, which is what a real-valued objective
differentiated with respect to real physical parameters needs.

`transpose` and `hermitian_adjoint` are kept as separate functions, because
**reciprocity is the unconjugated transpose and passivity is the conjugate
one**. Using the adjoint for reciprocity is a physics error that only appears
once the medium is lossy — which is exactly where B3 begins. Verified:
`b3.complex.transpose_and_adjoint_are_distinct` (0.0 vs 1.0).

`Complex<double>` is also registered as a `Scalar` in its own right, which is
what makes **complex-step differentiation** possible: instantiating the model at
`Scalar = Complex<double>` gives `Complex<Complex<double>>`, the inner unit
physical and the outer one an infinitesimal perturbation.

---

## 3. Dispersive models

    eps(w)     = eps_inf + sum_n F_n/(w0_n² − w² − i γ_n w)  [+ Drude]
    mu(w)      = mu_inf  + sum_m G_m/(W0_m² − w² − i Γ_m w)
    [mu^-1](w) = 1 / mu(w)
    zeta(w)    = z(w) K,   K real antisymmetric, K_23 = +1

**`mu` is parameterized and `mu^-1` is derived**, not the other way round. A
Lorentz resonance is a resonance of the permeability; parameterizing `mu^-1` as a
sum of Lorentz terms would place poles of `mu^-1` — that is, **zeros of `mu`** —
at the resonance frequencies, which is a different and generally non-passive
medium. Invisible for one weak pole, decisive for a strong one. Verified:
`b3.dispersion.mu_inverse_is_derived_from_mu`.

The magnetoelectric response uses a **fixed real shape** with a scalar dispersive
amplitude, so `zeta_adjoint = zeta^T` is structural (residual exactly **0.0**).
Independent arbitrary complex entries are not used, because they would break
reciprocity by construction.

---

## 4. Passivity — derived, not assumed

Converting the project's (E,B) blocks to the (E,H) representation:

    eps_EH  = eps + zeta mu zeta^T      xi_EH   = − zeta mu
    zeta_EH = − mu zeta^T               mu_EH   = mu

    C = [[eps_EH, xi_EH], [zeta_EH, mu_EH]]
    L = (C − C†) / (2i)              Hermitian by construction

**Passivity is `L` positive semidefinite for every w > 0**, obtained from
`<p> = (w/2) v† L v`. What the gate verifies is that condition **on the sampled
band**; see the scope table above for what that does and does not establish. Eigenvalues come from the real 12×12 embedding of the
Hermitian 6×6, so no complex Jacobi rotation had to be written.

**Reciprocity in this convention is `xi_EH = zeta_EH^T`** — derived, not
imported. Many references write `xi = −zeta^T`; that sign traces to the opposite
convention for `zeta` in `D = eps E − zeta B`, and this project's sign was
derived at B1C and cross-checked against Minkowski's moving-medium relations.

### "All imaginary parts positive" is not passivity

Gate case `b3.passivity.diagonal_imaginary_parts_are_not_sufficient` exhibits a
medium with

    min Im eps = 1.29e-6 > 0,   min Im mu = 1.70e-7 > 0

that is nevertheless **not passive**: `min eig(L) = −1.19e+2`. With near-lossless
diagonal blocks the off-diagonal `L_12 = (xi_EH − zeta_EH†)/(2i)`, carrying the
magnetoelectric coupling, sets the spectrum on its own. This is precisely the
error the milestone warned against, demonstrated rather than described.

The suite also **fails** an intentionally active oscillator (negative damping):
`min eig(L) = −1.11e-2`, with its poles in the upper half-plane. A passivity test
that cannot fail a medium says nothing about the ones it passes.

---

## 5. Causality — structural

The poles of `F/(w0² − w² − i γ w)` are `± sqrt(w0² − γ²/4) − i γ/2`, so
`Im w = −γ/2 < 0` for **both** roots and every parameter value with positive
damping, including the overdamped branch. Damping is expressed through a
**strictly positive transform** of an unconstrained variable, so no optimizer
step can produce an acausal model. This is hard feasibility, not a penalty that
can be traded away.

Kramers–Kronig appears **only as a diagnostic**, and is reported as a convergence
trend rather than a fixed tolerance — because a fixed threshold there measures
the quadrature, not the physics:

| integration band | residual |
|---|---|
| [1e-2, 1e2] | 6.09e-8 |
| [1e-3, 1e3] | 6.08e-11 |
| [1e-4, 1e4] | 6.08e-14 |

Three orders per decade of band width. A first implementation simply skipped the
samples nearest the pole on a log grid; that plateaued at **1.5e-2** regardless of
band width, because principal-value cancellation needs a grid symmetric about the
pole and a log grid is not. Singularity subtraction — splitting off
`g(w)/(w'²−w²)` with its closed-form principal value — fixed it. The plateau is
how it was caught: a truncation error must fall when the band widens.

---

## 6. Band objective

    J_band = (1/N) sum_k ||Z_model(w_k) − Z_target(w_k)||_G² / max(||Z_target||_G², δ)

`Z(w) = A11(Re blocks) + i A11(Im blocks)`, exactly, because (A11) is homogeneous
linear in the blocks with real coefficients — the geometry carries no frequency
and no imaginary part. The **verified** real reconstruction is reused rather than
a second complex implementation being written; the linearity is asserted
(`b3.band.a11_is_real_linear`), not assumed.

Reported separately, never collapsed: centre mismatch, band mean, worst frequency
and where it occurred, passivity penalty, parameter-bound penalty, robustness
term.

---

## 7. Hard feasibility

B2 used smooth penalties and admitted a small bound violation. B3 adds genuine
feasibility by two mechanisms used together: **bounded parameter transforms**
(the feasible set is the whole domain of the parameterization) and **projection**
after each accepted step.

| mode | bound violation |
|---|---|
| `hard feasible approximation` | **0.0 exactly** |
| `soft penalized approximation` | 2.24e-6 |

**A soft-penalty result is never labelled feasible**, whatever the size of its
violation. The two are reported under different names.

---

## 8. Optimizer

L-BFGS (two-loop recursion, 8 pairs) with Armijo backtracking, implemented
locally — no external dependency. Steepest descent was tried first and is
retained for B2, but it stalls here: unlike B2's exactly quadratic objective,
`J_band` is **rational** in the physical parameters, the resonance parameters are
far more sensitive than the backgrounds, and the feasibility transform curves the
landscape further. Plain descent reached only ~1e-7 in 400 iterations; L-BFGS
reaches ~1e-19.

---

## 9. Measured results

| quantity | value |
|---|---|
| min loss eigenvalue over the sampled band (passive case) | 2.46e-3 |
| max passivity violation | **0.0** |
| centre / mean / worst mismatch | 1.67e-18 / 9.44e-19 / 2.67e-18 |
| worst mismatch frequency | 4.60 |
| hard-mode bound violation | **0.0** |
| soft-mode bound violation | 2.24e-6 |
| gradient vs analytic (one Lorentz pole) | 7.10e-10 |
| gradient vs complex step | 1.79e-10 |
| batched HyPAD gradient (4 batches) | 1.45e-10 |
| one-pole / two-pole recovery `J_band` | 9.94e-19 / 9.99e-19 |
| worst parameter recovery error | 6.80e-9 |
| restricted one-pole residual | **1.34e-2** (nonzero by construction) |
| coordinate-change residual | 1.33e-15 |
| real-limit residual against B2 | **0.0** |

---

## 10. Two findings worth recording

### The model is not identifiable from a band far from its resonances

The band was originally 0.6–1.6 against resonances at 3.0 and 3.6. The response
was recovered to `J_band ~ 1e-10` while the physical parameters were still
**9.6e-2** away. That is not an optimizer failure: far below resonance a Lorentz
term is `≈ F/(w0² − w²)`, so `F` and `w0` trade against each other at nearly
constant response. Moving the band to **span** the resonances (2.0–5.0) dropped
the parameter error to **6.8e-9**.

Recovering a *response* and recovering *parameters* are different claims, and
which one holds depends on where you probe.

### The magnetoelectric shape is orientation dependent

Covariance under a coordinate change holds frequency by frequency (residual
1.33e-15) **when the observer frame is pushed forward**. Rebuilding the triad
independently by Gram–Schmidt in the new coordinates gives a measurably different
tensor, because the isotropic parts are triad-independent —
`sum_i e_i^a e_i^b = h^{ab}` for any triad — but the anisotropic shape `K` singles
out particular spatial directions. Transforming the setup means transforming the
frame. Both facts are gate cases.

---

## 11. What B3 does not claim

* Nothing about fabrication. There is no geometry, unit cell, or manufacturing
  constraint anywhere.
* Nothing about a specific material. A recovered model is a point in a parametric
  family that is passive and causal by construction and reproduces a target
  response over a band. Whether any physical medium realizes it is untouched.
* Nothing about propagation. There is no Maxwell solver, no dispersion of a
  wavepacket, no scattering.
* Nothing about measurable gravitational curvature. The laboratory target remains
  an effective electromagnetic metric.
