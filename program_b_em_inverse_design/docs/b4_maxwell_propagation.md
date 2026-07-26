# B4 — Maxwell propagation

Does the B3 constitutive model actually produce the electromagnetic propagation
the target effective metric predicts?

    q_physical -> C(w, x) -> Maxwell operator -> (E, H, k, w) -> observables

## Scope

> **A passing result demonstrates that the B3 constitutive model reproduces
> selected electromagnetic propagation observables of the target effective
> metric. It still does not establish a fabricated metamaterial, unit cell,
> antenna, coil, or experimental device.**

No metamaterial geometry, no antennas, no coils, no fabrication constraints.
B4 is still not a fabrication or unit-cell optimizer.

Gate `B4-MAXWELL-PROPAGATION`: **23/23 passed, 0 failed, 0 skipped** with OTIlib
enabled; **22 passed, 1 skipped** with it disabled, the one being the HyPAD
sensitivity case.

---

## 1. Maxwell signs, derived not copied

From this project's own kernel `f(t) = Re[f-hat e^{-iwt}]`, so `d/dt -> -iw`:

| equation | form |
|---|---|
| Faraday | `curl E = +i w B` |
| Ampère | `curl H = -i w D + J` |
| Gauss (magnetic) | `div B = 0` |
| Gauss (electric) | `div D = rho` |
| plane wave | `grad -> i k`, so `k x E = w B` and `k x H = -w D` |
| Poynting | `<S> = (1/2) Re[E x conj(H)]` |

**Signs from a source using `e^{+iwt}` are not copied.** Under that convention
both curl equations flip, and mixing them yields a medium that looks passive but
propagates as if it had gain — the equations stay self-consistent, so no residual
would catch it. The vacuum plane-wave case is what pins them independently:
`k x E = omega H` holds to **0.0** for the transverse mode, and would fail by a
factor of −1 under the wrong kernel.

---

## 2. The plane-wave operator

    [ K - w zeta_EH      -w mu_EH     ] [E]
    [ w eps_EH           K + w xi_EH  ] [H] = 0

with `K v = k x v`. Propagation requires this 6×6 to be singular.

**Residuals are measured as `||M v||` for the analytically known mode, not as a
smallest singular value.** The singular-value route goes through `M† M`, whose
smallest eigenvalue is accurate to about `eps ||M||² ~ 1e-16`; the singular value
is its **square root**, so `sigma_min` carries about **1e-8** of error even when
the mode is exact. That is a property of the method, and it is why the first
attempt at these cases failed at 1.6e-9 against a 1e-13 tolerance. The direct
residual is exact to machine precision:

| case | residual |
|---|---|
| vacuum, `\|k\| = omega`, 3 directions | **1.11e-16** |
| isotropic, `\|k\|² = eps mu omega²` | 1.11e-16 |
| transversality `k·E = k·H = 0` | **0.0** |
| off-shell at `\|k\| = 1.7 omega` (must NOT be singular) | 9.1e-1 |

---

## 3. The Fresnel polynomial

From `F_ab = k_a A_b - k_b A_a` and `k_a H^{ab} = 0`:

    W^{bd} = Z^{abcd} k_a k_c,      W A = 0

`W` annihilates `k` on both sides — by antisymmetry in `(cd)` and in `(ab)` — so
`det W` is **identically zero** and carries no information. The physics is in the
adjugate, which must then be rank one:

    adj(W)_{de} = lambda k_d k_e

with `lambda` a homogeneous **quartic** in `k`: the Tamm–Rubilar scalar.
`lambda(k) = 0` is the propagation condition.

**The root finder was replaced at B4C.** B4 located branches by bracketing sign
changes, which is blind to roots of **even multiplicity** — precisely the
metric-induced case. The doubled root was never actually recovered; the
factorization test passed only because it compared `lambda` against a separately
computed cone. B4C uses Durand–Kerner with Newton polish, clustering for
multiplicity, and a per-root residual and conditioning warning, and it recovers
the metric double root directly at separation **0.0**.

| test | result |
|---|---|
| vacuum factorizes as `(eta^{ab}k_a k_b)²`, one constant for all `k` | **2.63e-16** |
| conformal invariance of the cone, `Omega² = 2.7` | 1.21e-15 |
| metric medium branch split (must vanish) | **0.0** |
| axion term: change in Fresnel coefficients | **0.0** (B4C) |
| axion term: change in Z, and in the zeta block | 3.7e-1, 7.4e-1 (B4C) |
| generic anisotropic branch split (must not) | **2.43e-1** |

**Not every reciprocal medium is metric-induced.** `eps = diag(2, 3.5, 5)` is
reciprocal, passive and perfectly respectable, and its quartic does **not**
factorize into a doubled quadric: its two branches split by 24%.

**CORRECTED at B4C.** An earlier draft of this section claimed *metricity ⟺ no
birefringence*. That is withdrawn. What the tests establish is one direction
only:

> Metric-induced constitutive tensors are non-birefringent. Generic reciprocal
> constitutive tensors need not be. The converse requires additional assumptions
> and is subject to conformal/dilaton freedom and constitutive pieces that do not
> affect the Fresnel polynomial.

**The Fresnel polynomial does not uniquely reconstruct `Z^{abcd}`.** Two
independent reasons, both demonstrated at B4C: an axion term `alpha eps^{abcd}`
changes `Z` and changes the `zeta` block by exactly `2 alpha h` while leaving the
Fresnel coefficients **bit-identical**, because `W^{bd} = Z^{abcd}k_a k_c`
annihilates it outright; and the cone is conformally invariant, so `g_eff` is
fixed only up to a factor. See
[b4c_propagation_correction.md](b4c_propagation_correction.md).

---

## 4. Energy, group velocity, and directions that differ

`v_g = dw/dk` from the dispersion curve agrees with the independent
`c/(n + w dn/dw)` to **8.45e-12**.

`<S> = (1/2)Re[E x conj(H)]` gives `(0.5, 0, 0)` for a unit vacuum transverse
mode.

In a magnetoelectric medium the forward and backward wavenumbers differ in
**magnitude** — `|kappa_+| = 1.362` against `|kappa_-| = 2.202`. That is Fresnel
drag, not a sign error: `zeta = z K` is the moving-medium form. It is also why
wavevector, phase-velocity and energy-flow directions must not be assumed to
coincide.

---

## 5. The layered solver

**Scattering matrices with the Redheffer star product, not transfer matrices.**
A transfer matrix carries `e^{+i k_- d}` for the backward wave, growing as
`e^{|Im k| d}`, and overflows long before the physics goes wrong. The scattering
form uses only `e^{+i k_+ d}` and `e^{-i k_- d}`, both of modulus ≤ 1. A
**40-unit evanescent layer** returns `R = 0.995` and `T = 2.6e-49`, finite and
stable.

The layer dispersion is derived from this project's own (E,B) form:

    mu^-1 kappa² + 2 z kappa - eps = 0,    kappa_± = (-z ± sqrt(z² + eps mu^-1))/mu^-1

reducing to `kappa = ± sqrt(eps mu)` at `z = 0`.

**Absorption is computed twice, on purpose.** Once as `1 - R - T`, which is a
definition and proves nothing, and once by integrating the local dissipated power
`<p> = (w/2) Im[E* D + H* B]` through the internal field profile. Their agreement
is the genuine balance:

| case | result |
|---|---|
| lossless slab, `R + T = 1` | 3.33e-16 |
| lossless slab, integrated absorption | **0.0** |
| lossy slab, `R + T + A = 1` with `A` integrated | **6.87e-9** |
| lossy slab R / T / A | 7.62e-2 / 3.51e-1 / 5.73e-1 |
| active slab (`Im eps < 0`) | `A = -0.400`, `R + T = 1.400`, balance closes |

The active case matters: the balance closing on a **negative** absorption shows
the sign is physical rather than clamped.

---

## 6. Target-metric propagation benchmark

A 24-layer stack sampled from a subluminal Alcubierre slice (`v_s = 0.3`),
propagated exactly and with the magnetoelectric term dropped:

| | exact | restricted (`zeta = 0`) |
|---|---|---|
| R | 2.40e-32 | 5.66e-4 |
| T | 1.000 | 0.99943 |
| phase | −1.922 | −0.314 |

**Phase mismatch 1.609 rad, group-delay mismatch 1.610.** Both stacks conserve
energy. This is the point of the milestone: the tensor-norm residual of dropping
`zeta` looked modest at B2, but the **propagation** cost is a 1.6 radian phase
error — a tensor norm does not measure that.

---

## 7. Adaptive frequency validation

A coarse 21-point grid reported a minimum loss eigenvalue of **−31.3**. The
adaptive scan — refining around resonances, zeros of `det mu`, local minima and
intervals of rapid change — found **−499.6**, at a frequency 0.006 away from the
coarse worst point, using 981 samples over 6 rounds.

**The coarse grid was optimistic by a factor of 16.** Both numbers are reported,
always, because the difference between them is the measure of what a fixed grid
misses.

This narrows the unverified region. It does **not** prove global passivity, which
would need an analytic positive-real argument this project has not carried out,
and it is **not** a replacement for B3's sampled gate.

---

## 8. HyPAD sensitivities

`d(phase)/d(eps)` through a slab matches the analytic `w d / (2 sqrt(eps))` to
**0.0**, evaluated away from any eigenvalue degeneracy or branch crossing, where
the observable is smooth. Differentiating through branch sorting or discrete mode
selection is **not** done, and would have to be documented as non-smooth if it
were. The branch count along a direction is reported rather than assumed: 2 real
branches at `kappa = 1.871, 2.236` for the anisotropic medium.

---

## 9. What B4 does not claim

* Nothing about fabrication, geometry, unit cells, antennas or coils.
* Nothing about a full-wave solution in more than one dimension. The solver is 1D
  layered at normal incidence, which is exact for the model family B3 produces
  but is not a general Maxwell solver.
* Nothing about global passivity — see §7.
* Nothing about measurable gravitational curvature. The laboratory target remains
  an effective electromagnetic metric.
