#!/usr/bin/env python3
"""Aberration and Doppler sky animations, from the solver's saved sweep.

    python visualization/animation/animate_relativistic_sky.py \
        results/verification/program_c1/data/apparent_sky_sequence.csv \
        --metadata results/verification/program_c1/c1_metadata.json \
        --out results/verification/program_c1/animations \
        --contact-sheets results/verification/program_c1/contact_sheets

Produces relativistic_aberration and doppler_sky (GIF + MP4 + contact sheets).

THIS SCRIPT PERFORMS NO LORENTZ TRANSFORMATION. Every apparent direction and
Doppler factor at every beta was computed by c1_verify's four-vector boost and
written to the CSV. Sweeping beta in Python would be a second, unverified
simulation wearing this one's label -- which is exactly what the visualization
contract forbids. The frame index selects a saved beta block; that is all.

STAR IDENTITIES ARE STABLE across frames because rows are grouped by source_id
and drawn in one fixed sorted order for every frame.
"""

from __future__ import annotations

import math
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

SCOPE = ("Directions are Lorentz-transformed by the solver and read from file. "
         "Brightness is a bolometric grey approximation; colour is "
         "approximate_blackbody_visualization, not a predicted observed colour.")


def main(argv=None) -> int:
    p = cli.parser(__doc__.splitlines()[0])
    p.add_argument("--contact-sheets", required=True)
    p.add_argument("--fps", type=int, default=12)
    a = p.parse_args(argv)

    rows = records.read_csv(a.input)
    meta = records.RunMetadata.load(a.metadata)
    out_dir = records.ensure_dir(a.out)
    sheet_dir = records.ensure_dir(a.contact_sheets)
    theme = themes.get(a.theme)
    themes.apply(theme)

    # Group the saved sweep by beta index. Order within each block is by
    # source_id, so the same star occupies the same draw order in every frame.
    blocks = {}
    for r in rows:
        blocks.setdefault(int(r["beta_index"]), []).append(r)
    indices = sorted(blocks)
    for k in indices:
        blocks[k].sort(key=lambda r: r["source_id"])

    n_frames = len(indices)
    n_stars = len(blocks[indices[0]])
    # Every frame must contain the same stars: aberration moves sources, it never
    # creates or destroys them.
    for k in indices:
        if len(blocks[k]) != n_stars:
            raise records.RecordError(
                "beta block %d has %d stars against %d in the first block. "
                "Aberration moves sources; it does not create or destroy them, "
                "so this is a data defect." % (k, len(blocks[k]), n_stars))

    dest = meta.destination_name(0)

    def lonlat(r):
        x = float(r["apparent_direction_x"])
        y = float(r["apparent_direction_y"])
        z = float(r["apparent_direction_z"])
        return math.atan2(y, x), math.asin(max(-1.0, min(1.0, z)))

    # Exposure is FIXED across frames: the marker-area mapping is computed once
    # from the whole sweep, not per frame. A per-frame rescale would make the sky
    # appear to brighten and dim when only the mapping changed.
    all_bright = [math.log10(max(float(r["brightness_factor"]), 1e-12))
                  for k in indices for r in blocks[k]]
    b_lo, b_hi = min(all_bright), max(all_bright)
    lo_a, hi_a = theme.marker_area_range

    def area(bf):
        v = math.log10(max(bf, 1e-12))
        t = (v - b_lo) / (b_hi - b_lo) if b_hi > b_lo else 0.5
        return lo_a + (hi_a - lo_a) * (t ** 0.5)

    entries = []

    def sky_frame(i, color_by_doppler):
        blk = blocks[indices[i]]
        beta = float(blk[0]["beta"])
        gamma = float(blk[0]["gamma"])
        fig, ax = plt.subplots(figsize=(12.6, 5.6),
                               subplot_kw={"projection": "mollweide"})
        fig.subplots_adjust(left=0.26, right=0.98, top=0.88, bottom=0.10)
        ax.set_facecolor(theme.background)
        lo, la, cols, sz = [], [], [], []
        for r in blk:
            u, v = lonlat(r)
            lo.append(u); la.append(v)
            d = float(r["doppler_factor"])
            cols.append(themes.doppler_rgb(d) if color_by_doppler
                        else theme.star_color)
            sz.append(area(float(r["brightness_factor"])))
        ax.scatter(lo, la, s=sz, c=cols, edgecolors="none", alpha=0.95, zorder=4)
        ax.grid(True, alpha=theme.grid_alpha)
        ax.set_xticklabels([])
        ax.set_title("Apparent sky in the ship frame, toward %s" % dest,
                     fontsize=10, pad=14)
        dopp = [float(r["doppler_factor"]) for r in blk]
        aberr = [float(r["aberration_angle"]) for r in blk]
        rnd.figure_stamp(fig, theme, lines=[
            "beta          %8.5f" % beta,
            "gamma         %8.5f" % gamma,
            "sources       %8d  (fixed)" % len(blk),
            "Doppler       %.4f to %.4f" % (min(dopp), max(dopp)),
            "max aberr     %.3f deg" % math.degrees(max(aberr)),
            "colour        %s" % ("approximate_blackbody" if color_by_doppler
                                  else "none (geometry only)"),
            "              NOT a predicted observed colour",
            "exposure      fixed across all frames",
            "catalogue     %s" % meta.get("catalogue_release", "?"),
        ])
        rnd.scope_note(fig, theme, SCOPE)
        return fig

    for name, by_doppler in (("relativistic_aberration", False),
                             ("doppler_sky", True)):
        with tempfile.TemporaryDirectory(prefix="c1_sky_") as tmp:
            frames = rnd.save_frames(lambda i: sky_frame(i, by_doppler),
                                     n_frames, tmp, prefix=name)
            gif = rnd.write_gif(frames, os.path.join(out_dir, name + ".gif"),
                                fps=a.fps)
            mp4 = rnd.write_mp4(frames, os.path.join(out_dir, name + ".mp4"),
                                fps=a.fps, frame_dir=tmp)
            rnd.contact_sheet(frames, os.path.join(sheet_dir, name + ".png"),
                              tiles=9, label_prefix="beta#")
            from PIL import Image
            res = Image.open(frames[0]).size
            entries.append(rnd.animation_entry(
                name, mp4, codec="libx264", resolution=res, fps=a.fps,
                frames=frames, duration_s=len(frames) / a.fps,
                source_csv=records.resolve(a.input)))
            entries.append(rnd.animation_entry(
                name + "_gif", gif, codec="gif", resolution=res, fps=a.fps,
                frames=frames, duration_s=len(frames) / a.fps,
                source_csv=records.resolve(a.input)))

    prov = meta.figure_provenance(
        color_model="approximate_blackbody_visualization",
        render_parameters=dict(theme.describe(), frames=n_frames, fps=a.fps,
                               projection="mollweide",
                               exposure="fixed across frames; the marker mapping "
                                        "is computed once from the whole sweep",
                               beta_source="saved solver sweep, 48 blocks",
                               recomputed_in_python="nothing"),
        observer_frame="LocalShipFrame (instantaneous rest frame)")
    manifest = os.path.join(out_dir, "sky_manifest.json")
    rnd.write_manifest(manifest, entries, prov)
    print("wrote %s" % manifest)
    for e in entries:
        print("  %-30s %s" % (e["filename"] or "(not written)",
                              "ok" if e["written"] else e["note"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
