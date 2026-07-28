#!/usr/bin/env python3
"""Animate Program A4 curvature as the bubble velocity crosses v_s = 1.

Reads only the CSVs written by tools/a4_curvature_sweep.cpp. This module does
not evaluate a metric, invert anything, differentiate anything or know what a
Christoffel symbol is -- a Python reimplementation of the curvature pipeline
would drift from the verified C++ silently, and the whole point of the exporter
is that it cannot.

The animation exists to answer one question that static figures cannot:

    does anything singular happen to the curvature when g_00 reaches zero?

The answer the records give is no. det g = -1 identically for every v_s, and the
invariants grow smoothly straight through the crossing. What DOES happen at
v_s = 1 is that g_00 touches zero at the bubble centre, and one frame later the
saved grid brackets zero and a g_00 = 0 contour can be drawn. Those are two
different events and the overlay names both.

Usage:
    python animate_curvature.py --sweep <sweep dir> --out <gallery dir>
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

from common import render  # noqa: E402
from themes import themes  # noqa: E402


# --- Reading ---------------------------------------------------------------

def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def load_frames(sweep_dir):
    """Every frame_NNN.csv, in numeric order, as (xs, ys, {field: grid})."""
    names = sorted(n for n in os.listdir(sweep_dir)
                   if n.startswith("frame_") and n.endswith(".csv"))
    frames = []
    for name in names:
        rows = read_rows(os.path.join(sweep_dir, name))
        xs = sorted({float(r["x"]) for r in rows})
        ys = sorted({float(r["y"]) for r in rows})
        xi = {v: i for i, v in enumerate(xs)}
        yi = {v: i for i, v in enumerate(ys)}
        fields = {}
        for key in ("g00", "ricci_scalar", "kretschmann", "riemann_max"):
            g = [[0.0] * len(xs) for _ in ys]
            for r in rows:
                g[yi[float(r["y"])]][xi[float(r["x"])]] = float(r[key])
            fields[key] = g
        frames.append((xs, ys, fields))
    return frames


def global_range(frames, key):
    """One colour scale for the whole animation.

    A per-frame scale is the classic way to make an animation lie: every frame
    renormalises, so a field that grows by a factor of sixty across the sweep
    looks constant and only its noise appears to move. The A4 sweep grows
    max|R| from 0.022 to 1.33, which is exactly the signal a per-frame scale
    would destroy.
    """
    lo = min(min(min(row) for row in f[2][key]) for f in frames)
    hi = max(max(max(row) for row in f[2][key]) for f in frames)
    return lo, hi


# --- Rendering -------------------------------------------------------------

PANELS = [
    ("g00", r"$g_{00}$"),
    ("ricci_scalar", r"$R$"),
    ("kretschmann", r"$R_{abcd}R^{abcd}$"),
]


def make_factory(frames, summary, meta, theme, ranges):
    onset = int(meta.get("exact_onset_frame", -1))
    contour = int(meta.get("first_spatial_contour_frame", -1))
    radius = float(meta.get("radius", 2.0))

    def factory(i):
        xs, ys, fields = frames[i]
        row = summary[i]
        v_s = float(row["v_s"])
        extent = [min(xs), max(xs), min(ys), max(ys)]

        fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.6),
                                 facecolor=theme.background)
        for ax, (key, label) in zip(axes, PANELS):
            lo, hi = ranges[key]
            m = max(abs(lo), abs(hi))
            im = ax.imshow(fields[key], origin="lower", extent=extent,
                           aspect="equal", cmap="RdBu_r", vmin=-m, vmax=m,
                           interpolation="nearest")
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
            cb.ax.tick_params(labelsize=7, colors=theme.foreground)
            cb.outline.set_visible(False)

            # The nominal wall, from the manifest -- not fitted to the data.
            ax.add_patch(plt.Circle((0, 0), radius, fill=False,
                                    color=theme.foreground, linewidth=0.7,
                                    linestyle="--", alpha=0.5))

            # The g_00 = 0 contour, drawn only where the SAVED GRID brackets
            # zero. Below onset there is nothing to draw, and matplotlib
            # silently draws nothing -- which is the correct behaviour and the
            # visual form of the onset/contour distinction.
            if key == "g00":
                flat = [v for r in fields[key] for v in r]
                if min(flat) < 0.0 < max(flat):
                    ax.contour(xs, ys, fields[key], levels=[0.0],
                               colors=[theme.accent_two], linewidths=1.6)

            ax.set_title(label, fontsize=10, color=theme.foreground)
            ax.set_xlabel("x", fontsize=8, color=theme.foreground)
            ax.set_ylabel("y", fontsize=8, color=theme.foreground)
            ax.tick_params(labelsize=7, colors=theme.foreground)
            ax.set_facecolor(theme.background)
            for s in ax.spines.values():
                s.set_color(theme.grid)

        if i < onset:
            status = "below onset — $g_{00} < 0$ everywhere"
        elif i == onset:
            status = "EXACT ONSET — $v_s = 1$, $g_{00}(0) = 0$ exactly"
        elif i == contour:
            status = "first frame whose grid brackets zero — contour drawable"
        else:
            status = "past onset — $g_{00} > 0$ inside the bubble"

        fig.suptitle(
            "Alcubierre curvature through the $g_{00}=0$ crossing     "
            f"frame {i}/{len(frames) - 1}     $v_s$ = {v_s:.1f}",
            fontsize=12, color=theme.foreground, x=0.012, ha="left", y=0.985)
        fig.text(0.012, 0.932, status, fontsize=9, color=theme.accent_two,
                 ha="left")
        fig.text(0.012, 0.018,
                 "Fixed colour scale across all frames. Dashed circle is the "
                 "nominal radius R = 2. det g = −1 identically for every "
                 "$v_s$, and every invariant stays finite through the crossing "
                 "— which is a statement about CURVATURE ONLY. It does not "
                 "establish the absence of a horizon, the absence of causal "
                 "disconnection, safe traversal or physical realizability. "
                 "Sign changes in R are geometric and imply nothing about the "
                 "sign of any energy: no stress-energy, no observer, no energy "
                 "condition exists yet.",
                 fontsize=6.5, color=theme.foreground, alpha=0.75, ha="left",
                 va="bottom")
        fig.tight_layout(rect=(0, 0.05, 1, 0.90))
        return fig

    return factory


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fps", type=int, default=4)
    args = ap.parse_args(argv)

    frames = load_frames(args.sweep)
    if not frames:
        print("no frame CSVs found in " + args.sweep, file=sys.stderr)
        return 1
    summary = read_rows(os.path.join(args.sweep, "curvature_velocity_sweep.csv"))
    with open(os.path.join(args.sweep, "a4_curvature_manifest.json"),
              encoding="utf-8") as fh:
        meta = json.load(fh)

    if len(summary) != len(frames):
        print("summary rows (%d) and frame files (%d) disagree"
              % (len(summary), len(frames)), file=sys.stderr)
        return 1

    ranges = {key: global_range(frames, key) for key, _ in PANELS}

    os.makedirs(args.out, exist_ok=True)
    entries = []
    for theme_name in ("scientific", "cinematic"):
        theme = themes.get(theme_name)
        frame_dir = os.path.join(args.out, "_frames_" + theme_name)
        paths = render.save_frames(
            make_factory(frames, summary, meta, theme, ranges),
            len(frames), frame_dir, prefix="curvature", dpi=110)

        blank = [p for p in paths if render.frame_is_blank(p)]
        if blank:
            print("blank frames rendered: %s" % blank[:3], file=sys.stderr)
            return 1

        stem = "curvature_vs_velocity_" + theme_name
        gif = render.write_gif(paths, os.path.join(args.out, stem + ".gif"),
                               fps=args.fps)
        mp4 = render.write_mp4(paths, os.path.join(args.out, stem + ".mp4"),
                               fps=args.fps, frame_dir=frame_dir)
        sheet = render.contact_sheet(
            paths, os.path.join(args.out, "contact_sheets", stem + ".png"),
            tiles=9, label_prefix="v_s frame")

        src = os.path.join(args.sweep, "curvature_velocity_sweep.csv")
        for path, codec in ((gif, "gif"), (mp4, "h264")):
            entries.append(render.animation_entry(
                stem, path, codec=codec, resolution=[1430, 506], fps=args.fps,
                frames=paths, duration_s=len(frames) / float(args.fps),
                source_csv=src))
        print("wrote %s  gif=%s mp4=%s sheet=%s"
              % (stem, bool(gif), bool(mp4), bool(sheet)))

        # Raw frames are scratch and are never committed: what survives is the
        # GIF, the MP4 and the contact sheet.
        for p in paths:
            os.remove(p)
        os.rmdir(frame_dir)

    render.write_manifest(
        os.path.join(args.out, "animation_manifest.json"), entries,
        {
            "source_sweep": os.path.basename(os.path.normpath(args.sweep)),
            "git_commit": meta.get("git_commit"),
            "exact_onset_frame": meta.get("exact_onset_frame"),
            "first_spatial_contour_frame":
                meta.get("first_spatial_contour_frame"),
            "curvature_finite_through_crossing":
                meta.get("curvature_finite_through_crossing"),
            "colour_scale": "fixed globally across frames, per field",
            "renderer_computes_physics": False,
            "scope": meta.get("scope"),
        })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
