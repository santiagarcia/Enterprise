#!/usr/bin/env python3
"""The apparent sky, in the Earth frame and the ship frame.

    python visualization/static/render_apparent_sky.py \
        results/verification/program_c1/data/apparent_sky_ship_frame.csv \
        --metadata results/verification/program_c1/c1_metadata.json \
        --out results/verification/program_c1/figures --name sky_ship_frame

READS COLUMNS AND DRAWS THEM. Every apparent direction, Doppler factor and
brightness factor was computed by the solver's Lorentz boost and saved. This
script projects saved unit vectors onto plot axes and maps saved scalars to
colour and marker area. It performs no aberration and no Doppler calculation.

COLOUR MODEL: approximate_blackbody_visualization. The Doppler-to-RGB map is a
perceptual cue for the direction and rough size of the shift. It is NOT a
prediction of observed colour, which would need the source spectrum, the shifted
Planck function, CIE colour matching, and an adaptation model.
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import cli, records
from common import render as rnd
from themes import themes


def main(argv=None) -> int:
    p = cli.parser(__doc__.splitlines()[0])
    p.add_argument("--name", default="sky_frame", help="output basename")
    p.add_argument("--title", default=None)
    p.add_argument("--show-shift", action="store_true",
                   help="draw a segment from the geometric to the apparent "
                        "direction for the most-displaced stars")
    a = p.parse_args(argv)

    rows = records.read_csv(a.input)
    meta = records.RunMetadata.load(a.metadata)
    out_dir = records.ensure_dir(a.out)
    theme = themes.get(a.theme)
    themes.apply(theme)

    ids = records.text_column(rows, "source_id")
    ax_ = records.column(rows, "apparent_direction_x")
    ay = records.column(rows, "apparent_direction_y")
    az = records.column(rows, "apparent_direction_z")
    gx = records.column(rows, "geometric_direction_x")
    gy = records.column(rows, "geometric_direction_y")
    gz = records.column(rows, "geometric_direction_z")
    dopp = records.column(rows, "doppler_factor")
    bright = records.column(rows, "brightness_factor")
    aberr = records.column(rows, "aberration_angle")
    names = records.text_column(rows, "primary_name")

    # Deterministic order, keyed by source_id: star identities are stable and a
    # rerun draws them in the same z-order.
    order = sorted(range(len(ids)), key=lambda i: ids[i])

    def to_lonlat(x, y, z):
        """Equirectangular projection of a SAVED unit direction. Presentation."""
        return math.atan2(y, x), math.asin(max(-1.0, min(1.0, z)))

    fig, ax = plt.subplots(figsize=(13.2, 5.8),
                           subplot_kw={"projection": "mollweide"})
    # Leave room at the left for the provenance stamp: on a Mollweide the data
    # fills an ellipse and an in-axes box lands on the graticule labels.
    fig.subplots_adjust(left=0.28, right=0.98, top=0.88, bottom=0.10)
    ax.set_facecolor(theme.background)

    lons, lats, cols, areas = [], [], [], []
    for i in order:
        lo, la = to_lonlat(ax_[i], ay[i], az[i])
        lons.append(lo); lats.append(la)
        cols.append(themes.doppler_rgb(dopp[i]))
    # Marker area from the SAVED bolometric brightness factor, log-compressed:
    # D^4 spans orders of magnitude at high beta and a linear map would leave one
    # visible dot.
    areas = themes.marker_areas([math.log10(max(b, 1e-9)) for b in bright], theme)

    ax.scatter(lons, lats, s=areas, c=cols, edgecolors="none", alpha=0.95,
               zorder=4)

    if a.show_shift:
        # The aberration displacement itself, for the most-shifted stars.
        top = sorted(order, key=lambda i: -aberr[i])[:40]
        for i in top:
            lo0, la0 = to_lonlat(gx[i], gy[i], gz[i])
            lo1, la1 = to_lonlat(ax_[i], ay[i], az[i])
            # Skip segments that wrap the projection seam: drawing them would
            # paint a spurious line across the whole sky.
            if abs(lo1 - lo0) > math.pi:
                continue
            ax.plot([lo0, lo1], [la0, la1], color=theme.grid, lw=0.6,
                    alpha=0.75, zorder=2)

    ax.grid(True, alpha=theme.grid_alpha)
    ax.set_xticklabels([])
    frame_label = "ship frame" if "ship" in os.path.basename(a.input) else \
                  "solar-system barycentric frame"
    ax.set_title(a.title or ("Apparent sky, %s — %d Gaia DR3 sources"
                             % (frame_label, len(order))), fontsize=10, pad=16)

    max_ab = max(aberr) if aberr else 0.0
    rnd.figure_stamp(fig, theme, lines=[
        "frame        %s" % frame_label,
        "sources      %d" % len(order),
        "Doppler      %.4f to %.4f" % (min(dopp), max(dopp)),
        "max aberr    %.4f rad (%.2f deg)" % (max_ab, math.degrees(max_ab)),
        "colour       approximate_blackbody_visualization",
        "             NOT a predicted observed colour",
        "marker       area ~ log10(D^4 brightness factor)",
    ])
    rnd.scope_note(fig, theme,
                   "Directions are retarded (each star is shown where it WAS "
                   "when the light left) and Lorentz-transformed into the stated "
                   "frame. Brightness is a bolometric grey approximation.")

    prov = meta.figure_provenance(
        color_model="approximate_blackbody_visualization",
        render_parameters=dict(theme.describe(), projection="mollweide",
                               marker_source="log10 of the saved D^4 brightness "
                                             "factor",
                               shift_segments=bool(a.show_shift)),
        observer_frame=("LocalShipFrame (instantaneous rest frame)"
                        if "ship" in os.path.basename(a.input)
                        else "BarycentricCartesian"))

    for ext in ("png", "svg"):
        pth = os.path.join(out_dir, "%s.%s" % (a.name, ext))
        fig.savefig(pth, facecolor=fig.get_facecolor())
        records.write_sidecar(pth, prov)
        print("wrote %s" % pth)
    plt.close(fig)

    records.write_plotted_csv(
        os.path.join(out_dir, "%s.plotted.csv" % a.name),
        ["source_id", "primary_name", "lon_rad", "lat_rad", "doppler_factor",
         "brightness_factor", "aberration_angle_rad", "marker_area_pt2"],
        [[ids[i], names[i], lons[k], lats[k], dopp[i], bright[i], aberr[i],
          areas[k]] for k, i in enumerate(order)])
    return 0


if __name__ == "__main__":
    sys.exit(main())
