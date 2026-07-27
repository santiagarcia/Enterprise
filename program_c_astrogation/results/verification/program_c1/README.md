# C1-RELATIVISTIC-ASTROGATION records

Two halves, both required: **22 scientific cases** (`c1_verify`) and **10
visualization cases** (`c1_visual_verify.py`).

Real Gaia DR3 records throughout, from the committed `snapshots/gaia_dr3_10pc.txt`.
No live archive access at any point.

## The chain

    c1_verify  ->  data/*.csv + c1_metadata.json
               ->  visualization/make_c1_visuals.py
               ->  figures/ animations/ contact_sheets/

CTest enforces the order with `DEPENDS`, so a visual can never be checked against
a stale record. **No visualization script recomputes physics** — gate case
`c1vis.visuals_derived_from_saved_output_only` walks the AST of every script in
the visualization tree and fails on a call to any hyperbolic function, Lorentz
boost, aberration, Doppler, propagation or light-time routine.

## Layout

| Path | Contents |
| --- | --- |
| `data/` | the solver's records: worldline, both apparent skies, the 48-block β sweep, interception ledger, light cones, neighbourhood |
| `figures/` | PNG + SVG, each with a `.meta.json` sidecar and a `.plotted.csv` of the exact numbers drawn |
| `animations/` | GIF + MP4 per family, plus per-family manifests |
| `contact_sheets/` | 9 representative frames per animation, so a regression is visible without opening a video |
| `verification/` | the case-by-case reports for both halves |
| `c1_metadata.json` | the provenance every figure inherits |
| `animation_manifest.json` | consolidated: codec, resolution, frame rate, frame count, duration, source checksum, first- and final-frame checksums |

Raw animation frames are **not** committed: they are rendered into a temporary
directory and discarded, leaving the GIF, the MP4 and the contact sheet.

## The journey

1 g symmetric proper acceleration to Proxima Centauri (Gaia DR3
5853498713190525696), departing J2026.0:

| quantity | value |
| --- | --- |
| coordinate (Earth-frame) time | 5.8713 yr |
| **ship proper time** | **3.5417 yr** |
| peak β | 0.9481 |
| information lag at arrival | 4.2453 yr |
| stale-position targeting error | 79.7 au |
| Barnard's star, same error | 335.4 au |

The last two rows are the point of the interception solve: aiming at the
catalogue position misses, and misses further for a faster-moving star.

## Regenerating

```
c1_verify --snapshot ../../snapshots/gaia_dr3_10pc.txt --out <this dir>
python ../../../visualization/make_c1_visuals.py --results <this dir>
python ../../../visualization/tests/c1_visual_verify.py \
    --results <this dir> --visualization ../../../visualization
```

MP4 encoding needs ffmpeg (or the self-contained `imageio-ffmpeg` pip package).
It degrades gracefully: without an encoder the manifest records
`"written": false` with a reason, and **frame generation and PNG output are
unaffected and stay deterministic**.

## What these figures do not show

Special-relativistic kinematics only. No warp metric, no propulsion energy, no
reaction mass, no collision hazard, no radiation, no interstellar medium. Colour
in the sky views is `approximate_blackbody_visualization` — a perceptual cue for
the direction and rough size of the Doppler shift, **not** a prediction of
observed colour. Brightness is a bolometric grey approximation.
