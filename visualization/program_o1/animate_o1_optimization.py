#!/usr/bin/env python3
"""The O1 optimization animations, from saved C++ records only.

    python animate_o1_optimization.py --records <o1 dir> --out <gallery dir>

THREE THINGS VARY ACROSS FRAMES HERE AND NONE OF THEM IS PHYSICAL TIME:

    accepted iteration      the optimizer took a step
    design amplitude        a baseline-to-final interpolation at frozen time
    camera angle            nothing changes but the viewpoint

Every entry names its own axis.  Python reads records; it computes no GR.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

from program_o0c import o0c_common as O                          # noqa: E402
import plot_o1_optimization as P                                 # noqa: E402
from common import records, render                               # noqa: E402
from themes import themes                                        # noqa: E402

ANIMATIONS = (
    ("objective_progression", "accepted iteration"),
    ("design_trajectory_progression", "accepted iteration"),
    ("tradeoff_progression", "accepted iteration"),
)


def _render(fig):
    fig.canvas.draw()
    import numpy as np
    arr = np.asarray(fig.canvas.buffer_rgba()).copy()
    plt.close(fig)
    return arr


def _write(imgs, out, stem, fps=3):
    import numpy as np
    from PIL import Image
    d = records.ensure_dir(os.path.join(out, "_frames_" + stem))
    paths = []
    for k, a in enumerate(imgs):
        p = os.path.join(d, "o1_%05d.png" % k)
        Image.fromarray(np.asarray(a)).save(p)
        paths.append(p)
    gif = render.write_gif(paths, os.path.join(out, stem + ".gif"), fps=fps,
                           max_frames=len(paths), max_width=1000)
    mp4 = render.write_mp4(paths, os.path.join(out, stem + ".mp4"), fps=fps,
                           frame_dir=d)
    sheet = render.contact_sheet(
        paths, os.path.join(out, "contact_sheets", stem + ".png"),
        tiles=min(9, len(paths)), label_prefix="accepted iteration")
    res = "%dx%d" % (imgs[0].shape[1], imgs[0].shape[0])
    for p in paths:
        os.remove(p)
    os.rmdir(d)
    return {"gif": gif, "mp4": mp4, "sheet": sheet, "frames": len(paths),
            "resolution": res}


def _axis_note(fig, theme, text):
    fig.text(0.5, 0.915, text, ha="center", fontsize=8,
             color=theme.uncertainty_color)


def anim_objective(rec, out, theme, acc, cert):
    it = [float(r["iteration"]) for r in acc]
    ex = [float(r["exact_E_neg"]) for r in acc]
    ren = [float(r["renormalized_validation_E_neg_tau"]) for r in acc]
    imgs = []
    for k in range(len(acc)):
        fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.0),
                                 facecolor=theme.background)
        axes[0].axhline(ex[0], color=theme.grid, ls="--", lw=1.1)
        axes[0].plot(it[:k + 1], ex[:k + 1], "o-", color=P.RED, lw=2.2, ms=6)
        axes[0].set_xlim(-0.4, max(it) + 0.4)
        pad = (max(ex) - min(ex)) * 0.3 or 1e-6
        axes[0].set_ylim(min(ex) - pad, max(ex) + pad)
        O.style_axes(axes[0], theme, "accepted iteration",
                     r"exact $E_{\rm neg}$")
        axes[1].plot(it[:k + 1], ren[:k + 1], "o-", color=P.VIOLET, lw=2.2,
                     ms=6)
        axes[1].set_xlim(-0.4, max(it) + 0.4)
        pad2 = (max(ren) - min(ren)) * 0.3 or 1e-6
        axes[1].set_ylim(min(ren) - pad2, max(ren) + pad2)
        O.style_axes(axes[1], theme, "accepted iteration",
                     r"renormalized $E_{\rm neg,\tau}$")
        fig.text(0.5, 0.962,
                 "Exact energy and the optimized surrogate, iterate by iterate",
                 ha="center", fontsize=10.5, color=theme.foreground)
        _axis_note(fig, theme,
                   "FRAME AXIS: ACCEPTED ITERATION -- NOT PHYSICAL TIME")
        fig.text(0.5, 0.035,
                 "iterate %d of %d:  exact %.9f   (%.4f%% below baseline)"
                 % (k, len(acc) - 1, ex[k],
                    100.0 * (ex[0] - ex[k]) / ex[0] if ex[0] else 0.0),
                 ha="center", fontsize=8, color=theme.foreground)
        fig.tight_layout(rect=[0, 0.06, 1, 0.90])
        imgs.append(_render(fig))
    return _write(imgs, out, "objective_progression")


def anim_trajectory(rec, out, theme, acc, cert):
    """The accepted-iterate cloud filling in, in its own principal plane."""
    fy = O.read_csv(os.path.join(rec, "final_design_y.csv"))
    yfin = [float(r["y"]) for r in fy]
    norms = [float(r["y_norm"]) for r in acc]
    gnorm = [max(float(r["grad_norm"]), 1e-20) for r in acc]
    ex = [float(r["exact_E_neg"]) for r in acc]
    # Without a per-iterate y record the projection collapses to the norm
    # history, which is stated rather than dressed up as a 2-D cloud.
    imgs = []
    for k in range(len(acc)):
        fig, ax = plt.subplots(figsize=(11.0, 5.6), facecolor=theme.background)
        sc = ax.scatter(norms[:k + 1], gnorm[:k + 1], c=ex[:k + 1],
                        cmap="magma_r", s=110, edgecolor=theme.grid,
                        linewidth=0.6, zorder=3)
        ax.plot(norms[:k + 1], gnorm[:k + 1], "-", color=theme.grid, lw=1.0,
                zorder=1)
        ax.set_yscale("log")
        ax.set_xlim(-0.02 * max(norms) - 1e-6, max(norms) * 1.1 + 1e-6)
        ax.set_ylim(min(gnorm) * 0.4, max(gnorm) * 2.5)
        O.style_axes(ax, theme,
                     r"$\|y\|_2$   (the Sobolev design norm)",
                     r"$\|\nabla_y F_\tau\|_2$")
        cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
        cb.set_label(r"exact $E_{\rm neg}$", color=theme.foreground,
                     fontsize=8)
        cb.ax.tick_params(labelsize=6.5, colors=theme.foreground)
        fig.text(0.5, 0.962,
                 "Descent in DESIGN space: how far the design moved against "
                 "how much gradient is left", ha="center", fontsize=10.5,
                 color=theme.foreground)
        _axis_note(fig, theme,
                   "FRAME AXIS: ACCEPTED ITERATION -- NOT PHYSICAL TIME. "
                   "Axes are design-space norms, not physical space.")
        fig.tight_layout(rect=[0, 0.04, 1, 0.90])
        imgs.append(_render(fig))
    _ = yfin
    return _write(imgs, out, "design_trajectory_progression")


def anim_tradeoff(rec, out, theme, acc, cert):
    keys = [("exact_E_neg", "E_neg", True), ("E_abs", "E_abs", True),
            ("J_NEC_working", "J_NEC", True), ("K_max", "K_max", True),
            ("N_min_working", "|N_min|", True)]
    series = {}
    for k, lab, _ in keys:
        v = [abs(float(r[k])) for r in acc]
        base = v[0] if v[0] != 0 else 1.0
        series[lab] = [q / base for q in v]
    it = [float(r["iteration"]) for r in acc]
    imgs = []
    cols = [P.RED, P.AMBER, P.VIOLET, P.TEAL, P.BLUE]
    for k in range(len(acc)):
        fig, ax = plt.subplots(figsize=(11.6, 5.4), facecolor=theme.background)
        for (lab, c) in zip(series, cols):
            ax.plot(it[:k + 1], series[lab][:k + 1], "o-", color=c, lw=1.8,
                    ms=5, label=lab)
        ax.axhline(1.0, color=theme.grid, ls="--", lw=1.0)
        ax.set_xlim(-0.4, max(it) + 0.4)
        allv = [q for s in series.values() for q in s]
        ax.set_ylim(min(allv) * 0.98, max(allv) * 1.02)
        O.style_axes(ax, theme, "accepted iteration", "value / baseline")
        ax.legend(fontsize=7.4, facecolor=theme.background,
                  edgecolor=theme.grid, labelcolor=theme.foreground, ncol=5)
        fig.text(0.5, 0.962,
                 "Only E_neg was optimized: what happened to everything else",
                 ha="center", fontsize=10.5, color=theme.foreground)
        _axis_note(fig, theme,
                   "FRAME AXIS: ACCEPTED ITERATION -- NOT PHYSICAL TIME. "
                   "Above 1.0 means worse than baseline for these metrics.")
        fig.tight_layout(rect=[0, 0.04, 1, 0.90])
        imgs.append(_render(fig))
    return _write(imgs, out, "tradeoff_progression")


BUILDERS = (anim_objective, anim_trajectory, anim_tradeoff)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rec, out = records.resolve(a.records), records.ensure_dir(a.out)
    cert = O.read_json(os.path.join(rec, "final_certification.json"))
    acc = O.read_csv(os.path.join(rec, "accepted_iterates.csv"))
    theme = themes.get("scientific")
    themes.apply(theme)

    entries = []
    for fn, (stem, axis) in zip(BUILDERS, ANIMATIONS):
        p = fn(rec, out, theme, acc, cert)
        e = {"name": stem, "frame_axis": axis,
             "frame_axis_is_physical_time": False,
             "gif_written": bool(p["gif"]), "mp4_written": bool(p["mp4"]),
             "frames": p["frames"], "resolution": p["resolution"],
             "contact_sheet": (os.path.basename(p["sheet"]) if p["sheet"]
                               else None),
             "renderer_computes_physics": False, "fixed_scales": True,
             "axes_are_physical_space": False}
        if not p["mp4"]:
            e["mp4_unwritten_reason"] = (
                "no ffmpeg encoder available" if render.find_ffmpeg() is None
                else "encoder call failed")
        entries.append(e)
        print("wrote", stem, "gif" if e["gif_written"] else "",
              "mp4" if e["mp4_written"] else "(no mp4)")

    with open(os.path.join(out, "animation_manifest.json"), "w",
              encoding="utf-8") as f:
        json.dump({"animations": entries, "count": len(entries),
                   "scope": P.SCOPE,
                   "frame_axis_note":
                       "Optimization iteration is NOT physical time. A "
                       "baseline-to-final interpolation is NOT physical time. "
                       "A camera orbit is NOT physical time.",
                   "renderer_computes_physics": False}, f, indent=2)
    print("%d animations" % len(entries))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
