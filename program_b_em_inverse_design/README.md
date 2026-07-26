# Program B — electromagnetic inverse design

Establishing the correspondence

    target effective metric  <->  electromagnetic constitutive response

with correct conventions, stated realizability limits, and verification.

> **Status: B4 Maxwell propagation complete. B1 20/20, B1C-S23 42/42,
> B1C-S456 36/36, aggregate B1C 53/53, B2 33/33, B3 38/38, B4 23/23 -- 0 failures
> across the whole repository, with the OTIlib backend both enabled and disabled
> (26/26 and 25/25 CTest, HyPAD-only skips).**
>
> B4 shows the B3 constitutive model reproduces selected **propagation**
> observables of the target effective metric. It is NOT a fabricated
> metamaterial, unit cell, antenna, coil or experimental device.
>
> Causality is **global** and structural. Passivity is verified **over the tested
> band only** -- the adaptive scan of B4 narrows the unverified region but does
> not close it, and global passivity is not claimed.
>
> Forward (Plebanski-Tamm) map implemented and verified; inverse implemented for
> the non-magnetoelectric case with the conformal ambiguity explicit in the API;
> non-metric media detected and reported rather than fitted. The covariant
> route -- Z^{abcd}, its Bel decomposition with respect to a laboratory
> observer, full field reconstruction, and reconstruction of Z from the blocks
> -- supersedes the 3+1 route wherever they disagree.
>
> **Not supported:** any claim about physical realizability. No dispersion, no
> loss, no passivity, no frequency dependence. Nothing here generates measurable
> gravitational curvature; the laboratory target is an effective electromagnetic
> metric. See [docs/claims_and_limitations.md](docs/claims_and_limitations.md).

## Layout

```
constitutive/       Mat3, ConstitutiveResponse, metricity classification
inverse_design/     maxwell           Maxwell signs, plane-wave op, Fresnel
                    layered           1D scattering-matrix layered solver
                    adaptive_scan     adaptive passivity / resonance search
                    complex_scalar    paired-real Complex<Scalar>, AD-safe
                    dispersion        Fourier convention, Lorentz/Drude models
                    passivity         C(w), Hermitian loss matrix, spectrum
                    band_objective    Z(w), J_band, hard-feasible L-BFGS
                    design_space      21-parameter reciprocal design space
                    design_target     target definition and classification
                    design_objective  norms, smooth constraints, evaluate_design
                    design_gradient   finite-difference and manufactured cases
                    design_gradient_oti  batched HyPAD gradient (3 directions)
                    optimizer         gradient descent, backtracking Armijo
effective_metric/   metric_types      LaboratoryMetric vs EffectiveMetric
                    metric_medium_map Plebanski forward map, conformal-class
                                      inverse, Gordon metric
                    covariant_constitutive  Z^{abcd} from (6), observers
                    hodge_dual        4D Levi-Civita tensor and duals
                    bel_extraction    the (A7) observer blocks and projector
                    field_reconstruction  E,B,D,H and H^{ab} = Z^{abcd}F_{cd}
                    z_reconstruction  (A11), Z rebuilt from the blocks
                    legacy_unverified QUARANTINED; not on any gate path
verification/       b1_verify        B1 gate           (20 cases)
                    b1c_s23_verify   B1C-S23 gate      (42 cases)
                    b1c_s456_verify  B1C-S456 gate     (36 cases)
                    b1c_verify       full B1C gate     (53 cases)
                    b2_verify        B2-ABSTRACT-INVERSE (33 cases)
                    b3_verify        B3-PHYSICAL-ADMISSIBILITY (38 cases)
                    b4_verify        B4-MAXWELL-PROPAGATION (23 cases)
python/             figures, regenerated from saved CSV alone
docs/               theory review, conventions, claims register
examples/           (empty; the gates exercise everything)
```

## Build and verify

Built as part of the warp_hypad CMake project, which finds this module
automatically:

```sh
cd ../warp_hypad
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DWARP_HYPAD_WITH_OTILIB=ON
cmake --build build
ctest --test-dir build --output-on-failure -R programb1
```

Or a single gate, writing its record into this repository:

```sh
R=../program_b_em_inverse_design/results/verification
./build/b1_verify        --out $R/program_b1
./build/b1c_s23_verify   --out $R/program_b1c_s23
./build/b1c_s456_verify  --out $R/program_b1c_s456
python ../program_b_em_inverse_design/python/plot_b1.py $R/program_b1
```

The module depends on `warp_hypad_core` for `Mat4` and the verification harness,
and on nothing else. It does **not** depend on Program A's curvature pipeline —
that is deliberate, so B1 proceeds in parallel rather than waiting.

## The headline result, and the claim that was retracted

**Retracted.** B1 reported that the Plebanski map requires g_00 < 0, i.e.
v_s f < 1, and concluded that above it **no equivalent medium exists**. That
conflated two different statements. What B1 actually established is that the
3+1 decomposition *relative to the coordinate-time observer* is singular there.

The corrected result, established in B1C:

* the Alcubierre metric stays Lorentzian and nondegenerate everywhere, for every
  v_s, because its (t,x) block has determinant −1 and therefore always exactly
  one negative eigenvalue;
* the covariant constitutive tensor Z^{abcd} therefore remains defined and
  finite across v_s f = 1;
* v_s f = 1 is an **ergosurface**, structurally the same as Kerr's — the
  coordinate-time observer stops being timelike there, and other observers
  (the Eulerian one, for instance) do not;
* what fails is a *selected observer decomposition*, not the covariant
  metric-medium correspondence.

Susceptibility eigenvalues diverging on an ergosurface while the covariant
tensor does not is a known result. Whether the correspondence is *useful* into
v_s f > 1 — which observer family to use, and whether the resulting medium is
anywhere near realizable — remains open and is not claimed here.

## Documentation

| | |
|---|---|
| [docs/effective_metric_theory.md](docs/effective_metric_theory.md) | premetric electrodynamics, metricity conditions, optical vs Plebanski metrics, references |
| [docs/metric_constitutive_conventions.md](docs/metric_constitutive_conventions.md) | every equation with convention, signature, determinant factor, index positions, units, source, domain |
| [docs/claims_and_limitations.md](docs/claims_and_limitations.md) | what may and may not be claimed, and the retraction |
| [docs/b1c_source_equations.md](docs/b1c_source_equations.md) | Schuster & Visser equations as retrieved, with every derived normalisation and how it was measured |
| [docs/b1c_aggregate_gate_audit.md](docs/b1c_aggregate_gate_audit.md) | the six failed aggregate-gate assertions, classified, with withdrawals and replacements |
| [docs/b2_abstract_inverse_design.md](docs/b2_abstract_inverse_design.md) | B2 conventions, design space, objectives, constraints, test plan and measured results |
| [docs/b3_physical_admissibility.md](docs/b3_physical_admissibility.md) | B3 Fourier convention, dispersive models, passivity, causality, hard feasibility and measured results |
| [docs/b4_maxwell_propagation.md](docs/b4_maxwell_propagation.md) | B4 Maxwell signs, Fresnel polynomial, layered solver, energy balance, adaptive scan |
