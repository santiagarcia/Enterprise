# Program C, milestone C1 — relativistic astrogation

C1 computes a ship trajectory between the Solar System and a **real catalogue
destination**, keeping rigorously distinct every quantity that casual language
fuses together: coordinate time, ship proper time, Earth-frame distance, the
destination's worldline, a star's geometric position, its retarded apparent
position, and the interception point.

**C1 uses no effective warp metric.** It is the special-relativistic navigation
baseline that C2 will later compare an effective-metric voyage against. Nothing
here models propulsion energy, reaction mass, collision hazards, radiation, or
the interstellar medium. A 1 g journey to Proxima Centauri appears below as
kinematics, and the energy required to perform it is not addressed.

---

## 0. Conventions, fixed once

### Signature and units

Signature `(−,+,+,+)`, matching Programs A and B. Internal dynamical units are
those C0 established: **parsec, Julian year, parsec per Julian year**.

    c = 0.30660139...  pc/yr        (kSpeedOfLightParsecPerYear)
    1 g = 0.31648...   pc/yr^2      (9.80665 m/s^2, converted in one place)
    a/c at 1 g = 1.0323 per year    (the familiar "1 g is about 1 c per year")

Working in pc/yr rather than natural units is deliberate: it keeps every number
in the CSV outputs directly readable against the catalogue, and it means a unit
error shows up as a factor of 3.26 rather than hiding inside `c = 1`.

### Four-vectors

A four-vector is written `(t, x, y, z)` with the time component first. For a
photon, `k^a = ω (1, k̂)` where **`k̂` is the direction the photon PROPAGATES**,
not the direction toward its source. Those differ by a sign and confusing them is
the single most common error in aberration code, so the codebase names them
apart: `photon_direction` versus `direction_to_source`.

For an observer with coordinate velocity `β⃗`, the four-velocity is
`u^a = γ (1, β⃗)`, `γ = 1/√(1−β²)`, normalised to `u·u = −1`.

### Aberration sign convention — stated explicitly

Let `θ` be the angle between the **photon propagation direction `k̂`** and the
observer's velocity. Then

    cos θ' = (cos θ − β) / (1 − β cos θ)

If instead `θ` denotes the angle to the **direction toward the source**
`n̂ = −k̂`, the same physics reads

    cos θ' = (cos θ + β) / (1 + β cos θ)

Both are implemented and both are gate-tested against each other, because the
only defence against a sign error here is to write down both forms and assert
they agree when related by `n̂ = −k̂`. The physical content is the same in either:
**apparent source positions crowd toward the direction of motion.** A test that
merely reproduces one formula from itself would pass with the sign wrong; the
test that catches it is the one asserting that stars pile up ahead.

### Doppler

    ν' / ν = − (u_ship · k) / (u_reference · k)

With `u_reference = (1, 0, 0, 0)` and `k^a = ω(1, k̂)`, this reduces to

    ν' / ν = γ (1 − β⃗ · k̂)

For an **approaching** source the photon propagates opposite to the ship's
motion, `k̂ = −v̂`, so `β⃗·k̂ = −β` and

    ν'/ν = γ(1 + β) = √((1+β)/(1−β))       blueshift

and the reciprocal for a receding source. Both longitudinal limits are gate
cases. The transverse case `β⃗·k̂ = 0` gives `ν'/ν = γ`, the transverse Doppler
effect, which is pure time dilation and has no classical counterpart.

---

## 1. Destination representation

    struct Destination {
        std::uint64_t gaia_source_id;
        BarycentricCartesian state_at_reference_epoch;
        CovarianceResult uncertainty;
        Provenance provenance;
    };

**A destination is never identified by display name.** Names are resolved to a
`source_id` at selection time and the identifier is used internally thereafter.
"Proxima Centauri" is a label SIMBAD attached to Gaia DR3 5853498713190525696;
the number is the thing that can be rechecked.

## 2. Ship state and worldline

    struct ShipState {
        BarycentricCartesian coordinate_state;   // typed frame, not a bare Vec3
        double proper_time_years;                // tau
        double coordinate_time_years;            // t
        double beta, gamma;
        Vector3 proper_acceleration;             // pc/yr^2, in the ship frame
    };

    struct WorldlineSample {
        ShipState ship;
        BarycentricCartesian destination;   // where the target IS at this t
        double remaining_distance_pc;
        double received_earth_time;         // the Earth clock the ship can SEE
        double received_destination_time;
    };

`received_earth_time` is not "the current time on Earth". It is the latest Earth
event whose light has reached the ship. The distinction is the subject of §10.

## 3. Trajectory families

### Constant coordinate velocity — a verification model only

Unphysical as a mission profile (infinite acceleration at both ends), but it has
exact closed-form proper time `τ = t/γ`, which makes it the right instrument for
testing the integrator.

### Symmetric constant proper acceleration

Accelerate at constant proper acceleration `a` for half the coordinate distance,
then decelerate. For one-dimensional motion from rest the analytic relations are

    c t(τ) = (c²/a) sinh(aτ/c)
    x(τ)   = (c²/a) [cosh(aτ/c) − 1]
    β(τ)   = tanh(aτ/c)
    γ(τ)   = cosh(aτ/c)

The half-journey matching condition is applied to the **coordinate distance**:
each phase covers `D/2`, and the second phase is the time-reverse of the first,
so the ship arrives at rest. Total coordinate time and total proper time are
twice the single-phase values at `x = D/2`.

These are exact and are asserted as such, not as approximations.

### Bounded acceleration profile

A smooth sigmoid-ramped profile, so acceleration is continuous and jerk is
bounded. Provided for visual and engineering comparison.

**This is not claimed to be more physically realizable than constant proper
acceleration.** Without a propulsion model — mass ratio, exhaust velocity, power
density — "realizable" is not a statement anything here can support. It is a
different curve, and that is all that is claimed.

## 4. Moving-destination interception

The destination is not where the catalogue says by the time the ship arrives.
Solve

    x_ship(t_arrival) = x_star(t_arrival)

Reported side by side:

| quantity | meaning |
|---|---|
| stale catalogue target | the position at the J2016.0 reference epoch |
| destination at departure | propagated to the departure epoch |
| destination at arrival | propagated to the solved arrival epoch |
| interception point | where the two worldlines actually meet |
| targeting error | the miss distance from aiming at the stale position |
| interception uncertainty | the catalogue covariance propagated to arrival |

The solver is a bracketed bisection plus secant refinement on the scalar residual
`|x_ship(t)| − |x_star(t)|` along the interception direction, with the residual
and iteration count reported. Bisection is chosen over Newton because the
residual is monotone in `t` for these profiles and bracketing cannot diverge.

## 5. Clocks

    dτ = dt √(1 − v²/c²)

Four clocks are reported separately and **never** as an unqualified "travel
time":

    coordinate time         Earth / solar-system barycentric frame
    ship proper time        what the crew's clock reads
    destination coordinate time
    light-signal age        how old the light the ship currently sees is

## 6–9. Aberration, Doppler, brightness, the retarded sky

At each ship event, for each displayed source:

1. solve the emission time `t_e` from `|x_ship(t) − x_star(t_e)| = c (t − t_e)`;
2. propagate the source to `t_e`;
3. form the incoming photon direction;
4. Lorentz-transform it into the ship frame;
5. compute the Doppler factor;
6. write the record.

Brightness is decomposed rather than fused, because the exponents differ and a
single fudge factor would hide which one is wrong:

| factor | scaling |
|---|---|
| frequency shift | `D` |
| photon arrival rate | `D` |
| solid-angle aberration | `D²` |
| **bolometric specific intensity** | **`D⁴`** |

The `I_ν/ν³` invariant fixes the exponent; it is derived in the header rather
than transcribed. Rendered brightness is labelled a **physical approximation**
(bolometric, grey), never a photorealistic prediction.

**Colour.** No frequency shift is converted to RGB without a documented spectral
model. Every figure declares a `color_model` from the fixed set:

    physical_spectral_model
    approximate_blackbody_visualization
    catalogue_photometric_proxy
    illustrative_only

C1 uses `approximate_blackbody_visualization` for Doppler-shifted sky views —
a Planck spectrum at the catalogue-implied temperature, shifted by `D`, mapped
through CIE colour-matching functions. That is an approximation with a name, not
a prediction of what an eye would see.

## 10. Communication light cones

For each ship event, the latest Earth event whose signal could have arrived, and
the arrival times of ship-originated signals at Earth and at the destination.

**Earth's "current time" is never displayed aboard the ship as though information
arrived instantaneously.** The visuals show the received Earth clock alongside
the coordinate time, and the gap between them is the point.

---

## The visualization contract

A permanent requirement from this milestone onward:

    production solver → saved CSV/JSON/frame records → visualization scripts
                      → PNG/SVG/GIF/MP4

**Visualization scripts never recompute physics.** No trajectory, aberration,
Doppler factor, catalogue propagation or physical observable is calculated in a
plotting script. They read columns and draw them. This is enforced by a gate case
that greps the visualization sources for the physics symbols they must not
contain, and by the fact that every plotted number appears in a committed CSV.

Every script takes an explicit input path and output directory and **never
depends on the current working directory**.

### Themes

A **scientific** theme prioritising measurable geometry, uncertainty and labels,
and a **cinematic** theme which may add trails, glow and easing but **may not
alter positions, directions, timing or physical quantities**. The cinematic theme
draws the same numbers differently; it does not draw different numbers.

### Provenance on every figure

Each figure and animation carries an adjacent metadata file recording
`git_commit`, `git_dirty`, `catalogue`, `catalogue_release`, `snapshot_checksum`,
`data_kind`, `reference_epoch`, `simulation_epoch`, `coordinate_frame`,
`observer_frame`, `units`, `destination_source_id`, `trajectory_model`,
`physical_parameters`, `render_parameters`, `color_model` and
`generated_timestamp`.

---

## Gate `C1-RELATIVISTIC-ASTROGATION`

Twenty scientific cases and ten visualization cases. Acceptance: 0 scientific
failures, 0 visualization failures, 0 unknown observer frames, 0 ambiguous
clocks, 0 unlabelled illustrative colour models.

The gate uses the **committed Gaia DR3 snapshot**. No live archive access is
required at C1 runtime, and none occurs.

## After C1 — documented, not begun

**C2 effective-metric voyage** — connect the verified Program A metric to this
destination and worldline framework. This is where the dependency on Programs A
and B is created, and it runs one way only.

**C3 virtual cockpit** — interactive route playback, apparent sky, instruments.

**C4 mission environment** — interstellar medium, radiation, dust, collision risk
and uncertainty-aware route planning. This is where "is this survivable" first
becomes a question the code can address; C1 does not address it.

---

## Delivered state

**`C1-RELATIVISTIC-ASTROGATION`: 22 scientific cases + 10 visualization cases,
0 failures, 0 skips.** Twenty scientific cases were required; twenty-two exist
because two required cases split when a claim turned out to be false (below).

### The first journey

1 g symmetric proper acceleration to **Proxima Centauri** (Gaia DR3
5853498713190525696) and **Barnard's star** (4472832130942575872) — the latter
chosen because its 10.39 arcsec/yr proper motion, the largest of any known star,
makes the moving-target correction unmissable.

| | Proxima | Barnard |
|---|---|---|
| interception distance | 1.3020 pc | 1.8280 pc |
| coordinate time, 1 g | 5.8713 yr | 6.8035 yr |
| **ship proper time, 1 g** | **3.5417 yr** | 3.9200 yr |
| peak β | 0.9481 | 0.9632 |
| stale-position targeting error | 79.7 au | 335.4 au |
| interception uncertainty | 17.4 au | 21.8 au |

The stale-position error is the headline: aiming a 1 g ship at Proxima's
catalogue position misses by roughly twice Pluto's orbital radius, and at
Barnard's star by four times that. Neither is a rounding effect — it is the star
moving during a flight that takes years.

### Three tolerances that were wrong, and what fixing them showed

1. **Simpson at 1e-14** was unattainable: naive summation over 1001 weighted
   terms loses ~114 ulps, measured. Fixed by adding **Neumaier compensated
   summation** to `integrate_proper_time` — the integrator got better rather than
   the test getting weaker.
2. **`c1.zero_velocity_sky_recovery` at tolerance zero** was a false claim. At
   β = 0 the Doppler factor *is* bit-exact, but the apparent direction is not:
   `photon_four_vector` renormalises by `|k̂|`, and `(cos θ, sin θ, 0)` does not
   have norm exactly 1 in binary floating point. Renormalising an inexact unit
   vector is not the identity, and that is correct behaviour. The case split in
   two: Doppler at tolerance **0**, direction at **four ulps** with the reason
   stated.
3. **Lorentz norm preservation at 1e-15** was a few ulps too tight on order-one
   quantities. Loosened to 1e-14 with the ulp count named, since the norm is
   preserved *exactly* in real arithmetic and the residual is the rounding of the
   γ products.

### The aberration case that would have caught a sign error

`c1.aberration_analytic` checks the four-vector boost against the closed-form
formula — but a formula reproduced from itself passes with the sign inverted. The
case that actually protects the convention is
`c1.aberration_crowds_toward_motion`, which asserts the *physics*: apparent
source positions never move away from the direction of travel, at any of 37
angles. That is a statement no sign error survives.

### The visualization contract, mechanically enforced

`c1vis.visuals_derived_from_saved_output_only` parses the **AST** of every script
in the visualization tree and fails on a call to any hyperbolic function, Lorentz
boost, aberration, Doppler, catalogue-propagation or light-time routine. A regex
was tried first and rejected: it could not distinguish *recomputing* a Doppler
factor from *reading a saved column named* `doppler_factor`, which is precisely
what these scripts are supposed to do.

The check found one real hit — `themes.doppler_rgb` used `tanh` as a colour
compression curve. Rather than carve an exemption, the sigmoid was replaced with
the algebraic `u/(1+|u|)`. It does the same job, and "no hyperbolic function
anywhere in the visualization tree" stays a bright line needing no judgment call.

The β sweep behind the aberration animation is likewise computed **by the
solver** — 48 blocks × 315 sources, each a full Lorentz transformation, written
to `apparent_sky_sequence.csv`. Sweeping β in Python would have been a second,
unverified simulation wearing this one's label.
