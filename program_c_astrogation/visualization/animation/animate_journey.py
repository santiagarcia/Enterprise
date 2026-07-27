#!/usr/bin/env python3
"""Journey animations, rendered from saved records only.

    python visualization/animation/animate_journey.py \
        results/verification/program_c1/data/worldline_proxima_1g.csv \
        --metadata results/verification/program_c1/c1_metadata.json \
        --out results/verification/program_c1/animations \
        --contact-sheets results/verification/program_c1/contact_sheets

Produces journey_worldline and communication_delay (GIF + MP4 + contact sheet).

READS COLUMNS AND DRAWS THEM. Every position, clock, beta, gamma and signal time
comes from the CSV. This script integrates nothing and transforms nothing.

DETERMINISM: the frame index maps to a row index by a closed-form formula; axis
limits are computed once from the whole dataset and frozen; the overlay sits at a
fixed corner; nothing autoscales; there is no random state anywhere.
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

SCOPE = ("Special-relativistic kinematics. No warp metric, no propulsion energy, "
         "no collision hazard, no interstellar medium.")


def main(argv=None) -> int:
    p = cli.parser(__doc__.splitlines()[0])
    p.add_argument("--contact-sheets", required=True)
    p.add_argument("--frames", type=int, default=120)
    p.add_argument("--fps", type=int, default=24)
    p.add_argument("--light-cones", default=None)
    a = p.parse_args(argv)

    rows = records.read_csv(a.input)
    meta = records.RunMetadata.load(a.metadata)
    out_dir = records.ensure_dir(a.out)
    sheet_dir = records.ensure_dir(a.contact_sheets)
    theme = themes.get(a.theme)
    themes.apply(theme)

    tau = records.column(rows, "proper_time_yr")
    t = records.column(rows, "coordinate_time_yr")
    beta = records.column(rows, "beta")
    gamma = records.column(rows, "gamma")
    xs = records.column(rows, "x_pc")
    ys = records.column(rows, "y_pc")
    zs = records.column(rows, "z_pc")
    remaining = records.column(rows, "remaining_distance_pc")
    dist = [(xs[i] ** 2 + ys[i] ** 2 + zs[i] ** 2) ** 0.5 for i in range(len(xs))]

    lc_path = a.light_cones or os.path.join(
        os.path.dirname(records.resolve(a.input)), "light_cones.csv")
    lc = records.read_csv(lc_path)
    seen = records.column(lc, "earth_time_seen")
    age = records.column(lc, "earth_signal_age_yr")
    epoch0 = float(meta.get("simulation_epoch") or 0.0)

    dest = meta.destination_name(0)
    sid = meta.destination_source_id(0)
    n_rows = len(t)
    n_frames = min(a.frames, n_rows)

    # Closed-form frame -> row map. Never an incremental counter: an accumulator
    # would drift and make the last frame depend on the frame count.
    def row_of(frame):
        return round(frame * (n_rows - 1) / max(n_frames - 1, 1))

    # Limits computed ONCE from the whole dataset, then frozen. No frame
    # autoscales, so nothing appears to move because the axes moved.
    XLIM = (-0.02 * t[-1], 1.02 * t[-1])
    YLIM = (-0.02 * dist[-1], 1.06 * dist[-1])
    entries = []

    def stamp_lines(k):
        return [
            "coordinate time  %8.4f yr" % t[k],
            "ship proper time %8.4f yr" % tau[k],
            "beta             %8.6f" % beta[k],
            "gamma            %8.6f" % gamma[k],
            "remaining        %8.5f pc" % remaining[k],
            "Earth clock seen %8.4f yr  (lag %.4f yr)" % (seen[k] - epoch0, age[k]),
            "destination      %s" % dest,
            "source_id        %s" % sid,
            "catalogue        %s" % meta.get("catalogue_release", "?"),
        ]

    # --- 1. journey_worldline ----------------------------------------------
    def journey_frame(i):
        k = row_of(i)
        fig, ax = plt.subplots(figsize=(8.6, 5.2))
        ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
        ax.set_xlabel("coordinate time since departure [yr]")
        ax.set_ylabel("distance from the Sun [pc]")
        ax.set_title("Journey to %s — 1 g symmetric" % dest, fontsize=10)
        ax.grid(True, alpha=theme.grid_alpha)
        # The full path, faint: context that does not change between frames.
        ax.plot(t, dist, color=theme.grid, lw=0.9, alpha=0.55)
        # The trail is a SLICE of the same saved arrays, never a smoothed or
        # resampled curve, so the cinematic theme cannot move the ship.
        themes.glow(ax, t[:k + 1], dist[:k + 1], theme, color=theme.accent,
                    base_lw=2.2)
        ax.scatter([t[k]], [dist[k]], s=70, c=theme.accent_two, zorder=6,
                   edgecolors="none")
        ax.axhline(dist[-1], color=theme.grid, lw=0.7, ls=":")
        rnd.stamp(ax, theme, lines=stamp_lines(k), loc="upper left")
        rnd.scope_note(fig, theme, SCOPE)
        return fig

    with tempfile.TemporaryDirectory(prefix="c1_journey_") as tmp:
        # Frames go to a temporary directory. Thousands of raw frames are never
        # committed; what survives is the GIF, the MP4 and the contact sheet.
        frames = rnd.save_frames(journey_frame, n_frames, tmp, prefix="journey")
        gif = rnd.write_gif(frames, os.path.join(out_dir, "journey_worldline.gif"),
                            fps=min(a.fps, 20))
        mp4 = rnd.write_mp4(frames, os.path.join(out_dir, "journey_worldline.mp4"),
                            fps=a.fps, frame_dir=tmp)
        rnd.contact_sheet(frames,
                          os.path.join(sheet_dir, "journey_worldline.png"),
                          tiles=9, label_prefix="f")
        from PIL import Image
        res = Image.open(frames[0]).size
        entries.append(rnd.animation_entry(
            "journey_worldline", mp4, codec="libx264", resolution=res,
            fps=a.fps, frames=frames, duration_s=len(frames) / a.fps,
            source_csv=records.resolve(a.input)))
        entries.append(rnd.animation_entry(
            "journey_worldline_gif", gif, codec="gif", resolution=res,
            fps=min(a.fps, 20), frames=frames,
            duration_s=len(frames) / min(a.fps, 20),
            source_csv=records.resolve(a.input)))

    # --- 2. communication_delay ---------------------------------------------
    SLIM = (min(min(seen) - epoch0, 0.0) * 1.05, max(t) * 1.05)

    def comms_frame(i):
        k = row_of(i)
        fig, ax = plt.subplots(figsize=(8.6, 5.2))
        ax.set_xlim(*XLIM); ax.set_ylim(*SLIM)
        ax.set_xlabel("ship coordinate time since departure [yr]")
        ax.set_ylabel("Earth time since departure [yr]")
        ax.set_title("What the ship can know about Earth — %s" % dest,
                     fontsize=10)
        ax.grid(True, alpha=theme.grid_alpha)
        ax.plot(t, t, color=theme.grid, lw=0.9, ls="--")
        ax.plot(t[:k + 1], [seen[j] - epoch0 for j in range(k + 1)],
                color=theme.accent, lw=2.0)
        ax.fill_between(t[:k + 1], [seen[j] - epoch0 for j in range(k + 1)],
                        t[:k + 1], color=theme.uncertainty_color, alpha=0.12)
        ax.scatter([t[k]], [seen[k] - epoch0], s=64, c=theme.accent_two,
                   zorder=6, edgecolors="none")
        ax.annotate("information lag %.4f yr" % age[k],
                    xy=(t[k], 0.5 * (t[k] + seen[k] - epoch0)),
                    xytext=(-8, 0), textcoords="offset points", ha="right",
                    fontsize=7.4, color=theme.uncertainty_color)
        rnd.stamp(ax, theme, lines=stamp_lines(k), loc="lower right")
        rnd.scope_note(fig, theme,
                       "Earth's current time is never shown aboard the ship: "
                       "only the latest Earth event whose light has arrived.")
        return fig

    with tempfile.TemporaryDirectory(prefix="c1_comms_") as tmp:
        frames = rnd.save_frames(comms_frame, n_frames, tmp, prefix="comms")
        gif = rnd.write_gif(frames,
                            os.path.join(out_dir, "communication_delay.gif"),
                            fps=min(a.fps, 20))
        mp4 = rnd.write_mp4(frames,
                            os.path.join(out_dir, "communication_delay.mp4"),
                            fps=a.fps, frame_dir=tmp)
        rnd.contact_sheet(frames,
                          os.path.join(sheet_dir, "communication_delay.png"),
                          tiles=9, label_prefix="f")
        from PIL import Image
        res = Image.open(frames[0]).size
        entries.append(rnd.animation_entry(
            "communication_delay", mp4, codec="libx264", resolution=res,
            fps=a.fps, frames=frames, duration_s=len(frames) / a.fps,
            source_csv=records.resolve(lc_path)))
        entries.append(rnd.animation_entry(
            "communication_delay_gif", gif, codec="gif", resolution=res,
            fps=min(a.fps, 20), frames=frames,
            duration_s=len(frames) / min(a.fps, 20),
            source_csv=records.resolve(lc_path)))

    prov = meta.figure_provenance(
        color_model="illustrative_only",
        render_parameters=dict(theme.describe(), frames=n_frames, fps=a.fps,
                               frame_to_row="closed form, round(i*(N-1)/(F-1))",
                               axis_limits="computed once, frozen for all frames",
                               exposure="fixed"),
        observer_frame="BarycentricCartesian; ship proper time labelled "
                       "separately in the overlay")
    manifest = os.path.join(out_dir, "journey_manifest.json")
    rnd.write_manifest(manifest, entries, prov)
    print("wrote %s" % manifest)
    for e in entries:
        print("  %-28s %s" % (e["filename"] or "(not written)",
                              "ok" if e["written"] else e["note"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
