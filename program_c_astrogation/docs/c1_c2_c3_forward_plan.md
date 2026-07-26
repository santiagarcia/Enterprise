# Program C: forward plan for C1, C2 and C3

**Status: DOCUMENTATION ONLY. None of this is implemented.** C0 is complete and
gated; this file records what the later milestones must do, so that C0's
interfaces are shaped by a known destination rather than guessed at. Nothing here
is a claim about existing capability.

Read this alongside [`c0_architecture.md`](c0_architecture.md), which states the
data-honesty position governing the whole program.

---

## What C0 established, and what it deliberately did not

C0 delivers a catalogue-grounded, unit-safe, provenance-carrying local star
database with deterministic offline verification. Concretely: strong unit types,
five named coordinate frames each declaring its six facts, an explicit distance
estimator on every record, validity and quality bitmasks instead of bare NaN, a
checksummed versioned snapshot, and a 28-case gate.

C0 does **not** contain, and must not be read as containing:

| Absent from C0 | Belongs to |
| --- | --- |
| Relativistic aberration and light deflection | C1 |
| A ship worldline, proper time, or any observer in motion | C1 |
| Galactic potential, orbital integration, differential rotation | C1 (as an explicitly separate model) |
| Voyage planning, target interception, delta-v or energy budgets | C2 |
| Any coupling to the Program A/B effective-metric machinery | C2 |
| Rendering, a cockpit, a sky view, any graphics | C3 |

The rectilinear propagation in C0 is a **local kinematic model** with a stated
validity horizon of 10^5 yr and a warning flag that fires beyond it. It is not a
Galactic dynamics claim and the code says so at the point of use.

---

## C1 — relativistic observer kinematics

**Goal.** Given a ship worldline through the C0 catalogue volume, compute what an
observer on that worldline measures, as opposed to what the barycentric catalogue
records.

**Scope.**

1. **Ship worldline as a first-class typed object.** A `ShipWorldline` carrying
   proper time, coordinate time, position and four-velocity, in a frame that
   declares its six facts exactly as the C0 frames do. `LocalShipFrame` already
   exists in C0 as a *static* frame; C1 makes it a worldline.
2. **Aberration.** The relativistic transformation of apparent direction between
   the barycentric frame and the instantaneous ship rest frame. C0 has none: its
   directions are catalogue directions. The C1 test that matters is the
   low-velocity limit reducing to the classical `v/c` formula and the
   high-velocity limit reproducing the standard forward-beaming angle, both
   against closed-form references.
3. **Doppler and beaming.** Frequency shift and specific-intensity
   transformation, with the exponent stated and derived rather than transcribed
   (the `I_nu / nu^3` invariant fixes it).
4. **Light-time to a moving observer.** C0 solves light time to a *static*
   barycentric observer by fixed-point iteration. C1 must solve the retarded
   intersection between a star's worldline and the ship's past light cone, and
   must report convergence exactly as C0 does.
5. **Proper time along the worldline**, and the resulting distinction between
   catalogue epoch, coordinate arrival time, and elapsed ship time. This is where
   the "the star is not where the catalogue says by the time you arrive" problem
   becomes explicit rather than a caveat in a demo footer.
6. **A separate Galactic dynamics model**, if wanted at all, as an *alternative*
   propagator behind the same interface as C0's rectilinear one — never as a
   silent upgrade to it. A snapshot propagated by a different model must record
   which model.

**What C1 must not do.** It must not modify C0's catalogue schema to carry
observer-dependent quantities. Aberrated directions are *derived for an
observer*, not properties of a star, and storing them on `StellarObject` would
repeat exactly the laboratory/effective metric conflation that Program B had to
be repaired to avoid.

**Gate sketch: `C1-RELATIVISTIC-OBSERVER`.** Closed-form aberration limits;
Doppler invariant; retarded-time convergence and its residual; proper-time
integral against an analytic constant-acceleration case; the identity that a
zero-velocity ship reproduces C0's answers *exactly*, not approximately.

---

## C2 — voyage formulation

**Goal.** Pose and solve the problem "reach a destination worldline", with the
budgets and uncertainties stated.

**Scope.**

1. **Interception, not arrival at a coordinate.** The target is a worldline. The
   solved quantity is a rendezvous event. C0's demonstration ranks candidates by
   present catalogue distance and says in its own output that this is not a
   route; C2 is where that becomes a solved problem.
2. **Uncertainty propagated to the arrival event.** C0 already propagates a
   6x6 covariance with a labelled Jacobian and reports the growth ratio with
   epoch. C2 must carry that through to a position uncertainty *at the
   rendezvous*, because a 10^4 yr propagation of a 0.4 % parallax error is the
   dominant term in any realistic voyage, not a footnote.
3. **Budgets.** Whatever cost model is used — proper time, energy, an effective
   metric requirement — must be reported as separate named terms, never a single
   fused scalar. This is the same discipline B3's band objective and B2's
   objective already follow.
4. **The Program A/B coupling, and its honest limits.** This is the only
   milestone where Program C touches the effective-metric libraries, and the
   coupling must run **one way**: C2 may *query* A/B public interfaces for what
   an effective geometry would require, and must not push catalogue or route
   concerns into their cores. Programs A and B are verified libraries; a route
   planner is not permitted to become a reason to weaken one of their gates.
5. **The claim boundary, restated.** Program B's laboratory target is an
   effective electromagnetic or optical metric. It is not measurable
   gravitational curvature and it is not a propulsion result. A C2 voyage
   formulation is therefore a *formal* optimisation over a stated model, and any
   output must say so in the same breath as it reports a number.

**Gate sketch: `C2-VOYAGE-FORMULATION`.** Rendezvous solved against an analytic
straight-line case; a zero-velocity target reducing to the static problem;
covariance at arrival against a closed-form linear-growth reference; each budget
term reported separately; a hard-feasibility flag distinct from any soft penalty,
following B3.

---

## C3 — presentation

**Goal.** Render the C0/C1/C2 state without becoming the place where physics
quietly happens.

**Scope.** A sky view from the ship frame using C1's aberrated directions and
C0's photometry; a 3D volume view of the catalogue; a route view of a C2
solution. Every rendered quantity must be traceable to a computed one, and the
presentation layer must own **no** model of its own.

**Two rules carried forward from the rest of the project.**

1. **Precision honesty.** C0 already refuses to print more significant digits
   than the measurement uncertainty supports (`significant_digits_for`). A
   renderer that draws a star at a pixel implies a precision; where that exceeds
   the data, the uncertainty must be visible, not smoothed away.
2. **Provenance visible at the surface.** The C0 snapshot's `data_kind` field
   travels into the manifest, the metadata JSON, the SQLite projection and the
   demonstration header. It must reach the screen too. A view built from
   `SYNTHETIC` fixtures that looks identical to one built from Gaia DR3 is a
   presentation bug of the most serious kind.

**Gate sketch: `C3-PRESENTATION`.** Projection round trips; a rendered direction
agreeing with C1's computed direction to within a stated pixel tolerance; the
`data_kind` marker present in every rendered frame's metadata.

---

## Ordering and preconditions

C1 before C2: interception needs the observer kinematics. C2 before C3: there is
nothing to render a route from otherwise. Each milestone gets its own gate, and
none begins before the previous one is green — the same rule that governed B1C
through B4C.

The immediate precondition for C1 is a real snapshot. C0's machinery is verified
against synthetic fixtures with closed-form answers, which is the right way to
verify machinery; but C1's aberration and light-time work becomes scientifically
meaningful only against `data_kind = "Gaia DR3"`. Running
[`tools/import_gaia.py`](../tools/import_gaia.py) is therefore step zero of C1,
not an optional extra.
