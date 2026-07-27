#!/usr/bin/env python3
"""Rotating 3-D view of the local neighbourhood, from saved positions.

    python visualization/animation/animate_star_neighborhood.py \
        results/verification/program_c1/data/local_neighborhood.csv \
        --metadata results/verification/program_c1/c1_metadata.json \
        --out results/verification/program_c1/animations \
        --contact-sheets results/verification/program_c1/contact_sheets

Produces local_neighborhood_rotation (GIF + MP4 + contact sheet).

READS COLUMNS AND DRAWS THEM. Positions are the C0 derived columns. Only the
CAMERA moves; no coordinate is ever recomputed, and the azimuth is a closed-form
function of the frame index so the loop closes exactly.
"""

from __future__ import annotations

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import cli, records
from common import render as rnd
from themes import themes


def main(argv=None) -> int:
    p = cli.parser(__doc__.splitlines()[0])
    p.add_argument("--contact-sheets", required=True)
    p.add_argument("--frames", type=int, default=90)
    p.add_argument("--fps", type=int, default=20)
    p.add_argument("--max-stars", type=int, default=320)
    a = p.parse_args(argv)

    rows = records.read_csv(a.input)
    meta = records.RunMetadata.load(a.metadata)
    out_dir = records.ensure_dir(a.out)
    sheet_dir = records.ensure_dir(a.contact_sheets)
    theme = themes.get(a.theme)
    themes.apply(theme)

    x = records.column(rows, "x_pc")
    y = records.column(rows, "y_pc")
    z = records.column(rows, "z_pc")
    d = records.column(rows, "distance_pc")
    ids = records.text_column(rows, "source_id")
    names = records.text_column(rows, "primary_name")
    planets = records.column(rows, "has_planets")

    # Stable identities: sort by distance with source_id breaking ties, once,
    # outside the frame loop. Every frame then draws the same stars in the same
    # order.
    order = sorted(range(len(x)), key=lambda i: (d[i], ids[i]))[:a.max_stars]
    areas = themes.marker_areas([1.0 / d[i] for i in order], theme)
    lim = max(max(abs(v) for v in (x[i], y[i], z[i])) for i in order) * 1.05
    hosts = [i for i in order if planets[i] > 0]

    def frame(i):
        fig = plt.figure(figsize=(8.0, 7.0))
        ax = fig.add_subplot(111, projection="3d")
        themes.style_3d(ax, theme)
        ax.scatter([x[j] for j in order], [y[j] for j in order],
                   [z[j] for j in order], s=areas, c=theme.star_color,
                   edgecolors="none", alpha=0.9, depthshade=True)
        ax.scatter([0], [0], [0], s=150, marker="*", c=theme.accent_two,
                   edgecolors="none")
        for j in hosts[:6]:
            ax.scatter([x[j]], [y[j]], [z[j]], s=70, facecolors="none",
                       edgecolors=theme.accent, linewidths=0.9)
        # Limits frozen: no frame autoscales, so apparent motion is the camera.
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel("x [pc]"); ax.set_ylabel("y [pc]"); ax.set_zlabel("z [pc]")
        # Closed form in i, and the full 360 is divided by the frame count so the
        # loop closes exactly rather than jumping at the wrap.
        ax.view_init(elev=18.0, azim=-180.0 + 360.0 * i / a.frames)
        ax.set_title("Local neighbourhood, %s — %d sources"
                     % (meta.get("catalogue_release", "?"), len(order)),
                     fontsize=10)
        rnd.stamp(ax, theme, lines=[
            "frame       %3d / %d" % (i, a.frames),
            "azimuth     %7.2f deg  (camera only)" % (-180.0 + 360.0 * i / a.frames),
            "epoch       J%s" % meta.get("reference_epoch", "?"),
            "data_kind   %s" % meta.get("data_kind", "?"),
            "exoplanet hosts circled: %d" % len(hosts),
        ], loc="upper left")
        rnd.scope_note(fig, theme,
                       "Catalogue geometry at a single epoch. Only the camera "
                       "moves; no position is recomputed.")
        return fig

    with tempfile.TemporaryDirectory(prefix="c1_neigh_") as tmp:
        frames = rnd.save_frames(frame, a.frames, tmp, prefix="neigh")
        gif = rnd.write_gif(
            frames, os.path.join(out_dir, "local_neighborhood_rotation.gif"),
            fps=a.fps)
        mp4 = rnd.write_mp4(
            frames, os.path.join(out_dir, "local_neighborhood_rotation.mp4"),
            fps=a.fps, frame_dir=tmp)
        rnd.contact_sheet(
            frames, os.path.join(sheet_dir, "local_neighborhood_rotation.png"),
            tiles=9, label_prefix="f")
        from PIL import Image
        res = Image.open(frames[0]).size
        entries = [
            rnd.animation_entry("local_neighborhood_rotation", mp4,
                                codec="libx264", resolution=res, fps=a.fps,
                                frames=frames, duration_s=len(frames) / a.fps,
                                source_csv=records.resolve(a.input)),
            rnd.animation_entry("local_neighborhood_rotation_gif", gif,
                                codec="gif", resolution=res, fps=a.fps,
                                frames=frames, duration_s=len(frames) / a.fps,
                                source_csv=records.resolve(a.input)),
        ]

    prov = meta.figure_provenance(
        color_model="catalogue_photometric_proxy",
        render_parameters=dict(theme.describe(), frames=a.frames, fps=a.fps,
                               camera="closed form azimuth, elev fixed at 18 deg",
                               axis_limits="frozen", exposure="fixed"),
        observer_frame="BarycentricCartesian (camera orbit only; no observer "
                       "motion is modelled)")
    manifest = os.path.join(out_dir, "neighborhood_manifest.json")
    rnd.write_manifest(manifest, entries, prov)
    print("wrote %s" % manifest)
    for e in entries:
        print("  %-34s %s" % (e["filename"] or "(not written)",
                              "ok" if e["written"] else e["note"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
