#!/usr/bin/env python3
"""3-D local stellar neighbourhood, from saved catalogue records.

    python visualization/static/plot_local_neighborhood.py \
        results/verification/program_c1/data/local_neighborhood.csv \
        --metadata results/verification/program_c1/c1_metadata.json \
        --out results/verification/program_c1/figures

READS COLUMNS AND DRAWS THEM. Positions come from the snapshot's derived columns,
computed and verified by C0; nothing here converts a parallax, propagates an
epoch or recomputes a coordinate.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import cli, records
from themes import themes


def main(argv=None) -> int:
    p = cli.parser(__doc__.splitlines()[0])
    p.add_argument("--max-stars", type=int, default=400)
    a = p.parse_args(argv)

    rows = records.read_csv(a.input)
    meta = records.RunMetadata.load(a.metadata)
    out_dir = records.ensure_dir(a.out)
    theme = themes.get(a.theme)
    themes.apply(theme)

    x = records.column(rows, "x_pc")
    y = records.column(rows, "y_pc")
    z = records.column(rows, "z_pc")
    d = records.column(rows, "distance_pc")
    sig = records.column(rows, "distance_sigma_pc")
    names = records.text_column(rows, "primary_name")
    ids = records.text_column(rows, "source_id")
    planets = records.column(rows, "has_planets")

    # Stable, deterministic ordering: nearest first, ties broken by source_id.
    # Sorting by a float alone would let two equal distances swap between runs.
    order = sorted(range(len(x)), key=lambda i: (d[i], ids[i]))[:a.max_stars]

    fig = plt.figure(figsize=(9.5, 8.0))
    ax = fig.add_subplot(111, projection="3d")
    themes.style_3d(ax, theme)

    # Marker area from 1/distance -- a saved quantity. Stars are never all the
    # same size, and the mapping is recorded in render_parameters so a reader
    # knows a big dot means near, not large.
    areas = themes.marker_areas([1.0 / d[i] for i in order], theme)
    ax.scatter([x[i] for i in order], [y[i] for i in order],
               [z[i] for i in order], s=areas, c=theme.star_color,
               edgecolors="none", alpha=0.9, depthshade=True)

    # The Sun, at the barycentric origin by definition.
    ax.scatter([0], [0], [0], s=140, marker="*", c=theme.accent_two,
               edgecolors="none", zorder=10)
    ax.text(0, 0, 0, "  Sun", color=theme.accent_two, fontsize=8)

    # Label the nearest few and any exoplanet host, at positions DERIVED FROM THE
    # DATA -- never a random offset, which would move between renders.
    # Few labels, chosen deterministically. Crowding the origin with fifteen
    # overlapping names makes the figure less readable, not more informative.
    labelled = [i for i in order[:5]] + [i for i in order if planets[i] > 0][:5]
    for i in dict.fromkeys(labelled):
        ax.text(x[i], y[i], z[i], "  " + names[i], color=theme.foreground,
                fontsize=6.5)
        if planets[i] > 0:
            ax.scatter([x[i]], [y[i]], [z[i]], s=70, facecolors="none",
                       edgecolors=theme.accent, linewidths=0.9)

    lim = max(max(abs(v) for v in (x[i], y[i], z[i])) for i in order) * 1.05
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_xlabel("x [pc]"); ax.set_ylabel("y [pc]"); ax.set_zlabel("z [pc]")
    ax.set_title("Local stellar neighbourhood, %s\n%d sources, barycentric "
                 "ICRS-aligned, epoch J%s"
                 % (meta.get("catalogue_release", "?"), len(order),
                    meta.get("reference_epoch", "?")),
                 fontsize=10, color=theme.foreground)
    ax.grid(True, alpha=theme.grid_alpha)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=22, azim=-58)      # fixed camera: reproducible

    from common import render as rnd
    rnd.stamp(ax, theme, lines=[
        "catalogue   %s" % meta.get("catalogue_release", "?"),
        "data_kind   %s" % meta.get("data_kind", "?"),
        "checksum    %s" % str(meta.get("snapshot_checksum", "?"))[:18],
        "frame       barycentric ICRS-aligned",
        "marker      area ~ 1/distance (presentation)",
    ])
    rnd.scope_note(fig, theme,
                   "Catalogue geometry only. No trajectory, no relativity, no "
                   "warp metric. Positions are the C0 derived columns.")

    prov = meta.figure_provenance(
        color_model="catalogue_photometric_proxy",
        render_parameters=dict(theme.describe(), view_elev=22, view_azim=-58,
                               max_stars=a.max_stars,
                               marker_source="inverse distance from the CSV"),
        observer_frame="BarycentricCartesian (no observer motion in this figure)")

    written = []
    for ext in ("png", "svg"):
        p_out = os.path.join(out_dir, "local_neighborhood_3d." + ext)
        fig.savefig(p_out, facecolor=fig.get_facecolor())
        records.write_sidecar(p_out, prov)
        written.append(p_out)
    plt.close(fig)

    # The exact numbers plotted, so the picture can be checked against them.
    records.write_plotted_csv(
        os.path.join(out_dir, "local_neighborhood_3d.plotted.csv"),
        ["source_id", "primary_name", "x_pc", "y_pc", "z_pc", "distance_pc",
         "distance_sigma_pc", "marker_area_pt2"],
        [[ids[i], names[i], x[i], y[i], z[i], d[i], sig[i], areas[k]]
         for k, i in enumerate(order)])

    for p_out in written:
        print("wrote %s" % p_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
