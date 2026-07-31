#!/usr/bin/env python3
"""The six O0C animations, from saved C++ records only.

    python animate_o0c_energy.py --records <o0c_energy dir> --out <gallery dir>

EVERY ANIMATION DECLARES WHAT ITS FRAME AXIS IS.  Three different things can
vary from frame to frame here and only one of them is time:

    physical time        the bubble translates because it moves
    correction amplitude a design coordinate is scaled, at frozen time
    quadrature refinement the integration is refined, at frozen everything
    camera angle         nothing physical changes at all

A camera orbit is not physical time, and the manifest says so per entry rather
than leaving it to a caption.
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
from matplotlib.collections import PolyCollection                # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

import o0c_common as O                                           # noqa: E402
import plot_o0c_energy as P                                      # noqa: E402
from common import records, render                               # noqa: E402
from themes import themes                                        # noqa: E402

RED, AMBER, GREEN, VIOLET, BLUE = (P.RED, P.AMBER, P.GREEN, P.VIOLET, P.BLUE)


def _write(fig_frames, out, stem, theme, fps=8):
    """Write GIF, MP4 and a contact sheet, through the shared render helpers.

    Frames go to disk as PNGs because that is what ffmpeg consumes and what the
    contact sheet is built from; the frame directory is removed afterwards, so
    no raw frame tree is ever committed.
    """
    import numpy as np
    from PIL import Image

    frame_dir = records.ensure_dir(os.path.join(out, "_frames_" + stem))
    paths = []
    for k, arr in enumerate(fig_frames):
        fp = os.path.join(frame_dir, "o0c_%05d.png" % k)
        Image.fromarray(np.asarray(arr)).save(fp)
        paths.append(fp)

    n = len(paths)
    gif = render.write_gif(paths, os.path.join(out, stem + ".gif"),
                           fps=fps, max_frames=n, max_width=900)
    mp4 = render.write_mp4(paths, os.path.join(out, stem + ".mp4"),
                           fps=fps, frame_dir=frame_dir)
    sheet = render.contact_sheet(
        paths, os.path.join(out, "contact_sheets", stem + ".png"),
        tiles=min(9, n), label_prefix="frame")
    blanks = [q for q in paths if render.frame_is_blank(q)]

    result = {"frames": n, "gif": gif, "mp4": mp4, "sheet": sheet,
              "blanks": len(blanks),
              "resolution": "%dx%d" % (fig_frames[0].shape[1],
                                       fig_frames[0].shape[0])}
    for q in paths:
        os.remove(q)
    os.rmdir(frame_dir)
    return result


def _render(fig):
    fig.canvas.draw()
    buf = fig.canvas.buffer_rgba()
    import numpy as np
    arr = np.asarray(buf).copy()
    plt.close(fig)
    return arr


def _shell_frame(mer, scales, field, key, colour, theme, azim, elev, title,
                 subtitle, with_cabin=True, extra_text=None):
    faces, nseg = P._shell_faces(mer, field, scales["isovalues"][key])
    fig, ax = plt.subplots(figsize=(9.6, 6.0), facecolor=theme.background)
    P._setup_3d_axes(ax, theme)
    P._draw_shell(ax, faces, theme, colour, azim, elev, alpha=0.32)
    if with_cabin:
        P._draw_shell(ax, P._cabin_faces(
            scales["protected_cabin_radius_comoving"]), theme, GREEN, azim,
            elev, alpha=0.55, lw=0.15)
    fig.text(0.5, 0.955, title, ha="center", fontsize=10,
             color=theme.foreground)
    fig.text(0.5, 0.908, subtitle, ha="center", fontsize=7.6,
             color=theme.uncertainty_color)
    fig.text(0.5, 0.872, O.RECONSTRUCTION_NOTE, ha="center", fontsize=6.8,
             color=theme.grid)
    if extra_text:
        fig.text(0.5, 0.04, extra_text, ha="center", fontsize=7.0,
                 color=theme.foreground)
    fig.text(0.985, 0.962, "NOT READY FOR OPTIMIZATION", ha="right",
             fontsize=6.6, color=theme.uncertainty_color)
    return _render(fig)


def anim_time_translation(rec, out, man, theme, frames):
    scales = O.read_json(os.path.join(rec, "field_scales.json"))
    mman = O.read_json(os.path.join(rec, "meridional_manifest.json"))
    imgs = []
    for k in range(frames):
        mer = O.meridional(rec, k)
        t = -3.0 + 6.0 * k / (frames - 1)
        imgs.append(_shell_frame(
            mer, scales, "energy_density", "negative_energy", RED, theme,
            0.62, 0.34,
            "Negative-energy shell translating with the bubble",
            "FRAME AXIS: PHYSICAL TIME",
            extra_text="t = %+.2f, bubble centre x_c = %+.2f m; the structure "
                       "moves because the bubble does" % (t, mman["velocity"] * t)))
    return _write(imgs, out, "physical_time_3d_translation", theme)


def anim_scientific_shell(rec, out, man, theme, frames):
    scales = O.read_json(os.path.join(rec, "field_scales.json"))
    imgs = []
    for k in range(frames):
        mer = O.meridional(rec, k)
        imgs.append(_shell_frame(
            mer, scales, "n_min", "nec_violation", VIOLET, theme, 0.45, 0.28,
            "NEC-violation region, fixed camera and fixed isovalue",
            "FRAME AXIS: PHYSICAL TIME",
            extra_text="scientific style: the camera, the isovalue and the "
                       "colour scale are all frozen, so every change is the "
                       "field itself"))
    return _write(imgs, out, "scientific_3d_energy_shell", theme)


def anim_camera_orbit(rec, out, man, theme, frames):
    scales = O.read_json(os.path.join(rec, "field_scales.json"))
    mer = O.meridional(rec, 6)
    imgs = []
    n = max(frames, 16)
    for k in range(n):
        az = 2.0 * math.pi * k / n
        imgs.append(_shell_frame(
            mer, scales, "energy_density", "negative_energy", RED, theme,
            az, 0.30 + 0.10 * math.sin(2.0 * math.pi * k / n),
            "Cinematic orbit of the negative-energy shell",
            "FRAME AXIS: CAMERA ANGLE -- NOT PHYSICAL TIME",
            extra_text="the spacetime is FROZEN at t = 0.0; only the camera "
                       "moves. azimuth %.2f rad" % az))
    return _write(imgs, out, "cinematic_3d_camera_orbit", theme, fps=12)


def anim_correction_amplitude(rec, out, man, theme, frames):
    """The F1W protected regions, as the correction amplitude is scaled.

    The amplitude axis is read from the adversarial record: production already
    evaluated Delta beta^x in the protected regions for coefficient vectors
    spanning six orders of magnitude, and every one of them is exactly zero.
    """
    rows = [r for r in O.read_csv(os.path.join(rec, "family_mode_accounting.csv"))
            if r["family"] == "F1W"]
    adv = O.read_csv(os.path.join(rec, "adversarial_candidate_checks.csv"))
    scales = O.read_json(os.path.join(rec, "field_scales.json"))
    mer = O.meridional(rec, 6)
    imgs = []
    for k in range(len(adv)):
        cand = adv[k]
        faces, _ = P._shell_faces(mer, "energy_density",
                                  scales["isovalues"]["negative_energy"])
        fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.6),
                                 facecolor=theme.background)
        P._setup_3d_axes(axes[0], theme)
        P._draw_shell(axes[0], faces, theme, RED, 0.62, 0.34, alpha=0.28)
        P._draw_shell(axes[0], P._cabin_faces(
            scales["protected_cabin_radius_comoving"]), theme, GREEN, 0.62,
            0.34, alpha=0.6, lw=0.15)

        names = [r["candidate"] for r in adv]
        vals = [abs(float(r["interior_delta_beta"])) + 1e-20 for r in adv]
        cols = [GREEN if i == k else theme.grid for i in range(len(adv))]
        axes[1].barh(range(len(adv)), vals, color=cols,
                     edgecolor=theme.grid, linewidth=0.5)
        axes[1].set_xscale("log")
        axes[1].set_yticks(range(len(adv)))
        axes[1].set_yticklabels(names, fontsize=6.4)
        O.style_axes(axes[1], theme,
                     r"$|\Delta\beta^x|$ inside the protected ball", "")
        axes[1].set_xlim(1e-21, 1e-15)

        fig.text(0.5, 0.955,
                 "Protected-region lock under adversarial correction amplitudes",
                 ha="center", fontsize=10, color=theme.foreground)
        fig.text(0.5, 0.908, "FRAME AXIS: CORRECTION AMPLITUDE, at frozen time",
                 ha="center", fontsize=7.6, color=theme.uncertainty_color)
        fig.text(0.5, 0.045,
                 "candidate '%s': max |coefficient| = %s, interior "
                 "|Delta beta^x| = %s -- exactly zero, structurally"
                 % (cand["candidate"], cand["max_coefficient"],
                    cand["interior_delta_beta"]),
                 ha="center", fontsize=7.2, color=theme.foreground)
        fig.text(0.985, 0.962, "NOT READY FOR OPTIMIZATION", ha="right",
                 fontsize=6.6, color=theme.uncertainty_color)
        imgs.append(_render(fig))
    _ = rows
    return _write(imgs, out, "correction_amplitude_3d", theme, fps=2)


def anim_quadrature(rec, out, man, theme, frames):
    order = O.read_csv(os.path.join(rec, "quadrature_order_convergence.csv"))
    sub = O.read_csv(os.path.join(rec, "subdivision_convergence.csv"))
    dom = O.read_csv(os.path.join(rec, "domain_convergence.csv"))
    stages = []
    for rows, key, label in ((order, "order", "quadrature order"),
                             (sub, "subdivision", "element subdivision"),
                             (dom, "r_max_comoving", "domain radius [m]")):
        sel = [r for r in rows if r["objective"] == "E_neg"]
        for i in range(1, len(sel) + 1):
            stages.append((label, key, sel[:i], sel))
    imgs = []
    for label, key, shown, full in stages:
        fig, ax = plt.subplots(figsize=(9.6, 5.4), facecolor=theme.background)
        xs = [float(r[key]) for r in full]
        ax.set_xlim(min(xs) * 0.92, max(xs) * 1.08)
        ys_all = [float(r["value"]) for r in full]
        pad = (max(ys_all) - min(ys_all)) or 1e-6
        ax.set_ylim(min(ys_all) - 0.4 * pad, max(ys_all) + 0.4 * pad)
        ax.plot([float(r[key]) for r in shown],
                [float(r["value"]) for r in shown], "o-", color=RED, lw=2.0,
                ms=6)
        O.style_axes(ax, theme, label, r"$E_{\rm neg}$")
        fig.text(0.5, 0.955, "Convergence of the negative-energy objective",
                 ha="center", fontsize=10, color=theme.foreground)
        fig.text(0.5, 0.905, "FRAME AXIS: QUADRATURE REFINEMENT, at frozen "
                 "time and frozen design", ha="center", fontsize=7.6,
                 color=theme.uncertainty_color)
        last = shown[-1]
        fig.text(0.5, 0.03,
                 "%s = %s, E_neg = %.9f, change from previous %.3e"
                 % (label, last[key], float(last["value"]),
                    float(last["absolute_change"])),
                 ha="center", fontsize=7.2, color=theme.foreground)
        fig.text(0.985, 0.962, "NOT READY FOR OPTIMIZATION", ha="right",
                 fontsize=6.6, color=theme.uncertainty_color)
        imgs.append(_render(fig))
    return _write(imgs, out, "quadrature_convergence", theme, fps=2)


def anim_decomposition(rec, out, man, theme, frames):
    scales = O.read_json(os.path.join(rec, "field_scales.json"))
    imgs = []
    for k in range(frames):
        mer = O.meridional(rec, k)
        fig, axes = plt.subplots(1, 3, figsize=(15.4, 4.6),
                                 facecolor=theme.background)
        xs = [mer["x"][i][0] for i in range(mer["nx"])]
        rs = [mer["rho"][0][j] for j in range(mer["nr"])]
        for ax, field, lo, hi, cm, lab in (
                (axes[0], "beta_x", scales["beta_x"]["min"], 0.0, "Blues_r",
                 r"$\beta^x$"),
                (axes[1], "energy_density", scales["energy_density"]["min"],
                 0.0, "magma", r"$\rho$ (Eulerian)"),
                (axes[2], "n_min", scales["n_min"]["min"], 0.0, "viridis",
                 r"$N_{\min}$")):
            grid = [[mer[field][i][j] for i in range(mer["nx"])]
                    for j in range(mer["nr"])]
            im = ax.imshow(grid, origin="lower",
                           extent=[min(xs), max(xs), min(rs), max(rs)],
                           aspect="auto", cmap=cm, vmin=lo, vmax=hi,
                           interpolation="nearest")
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
            cb.ax.tick_params(labelsize=6.2, colors=theme.foreground)
            O.style_axes(ax, theme, "x  [m]", r"$\rho$  [m]")
            ax.set_title(lab, fontsize=9, color=theme.foreground)
        fig.text(0.5, 0.962,
                 "Shift, Eulerian energy density and the null minimum, "
                 "on one fixed scale", ha="center", fontsize=10,
                 color=theme.foreground)
        fig.text(0.5, 0.918, "FRAME AXIS: PHYSICAL TIME", ha="center",
                 fontsize=7.6, color=theme.uncertainty_color)
        fig.text(0.985, 0.968, "NOT READY FOR OPTIMIZATION", ha="right",
                 fontsize=6.6, color=theme.uncertainty_color)
        fig.tight_layout(rect=[0, 0.02, 1, 0.90])
        imgs.append(_render(fig))
    return _write(imgs, out, "objective_decomposition", theme)


BUILDERS = (
    ("physical_time_3d_translation", anim_time_translation, "physical time",
     "scientific"),
    ("scientific_3d_energy_shell", anim_scientific_shell, "physical time",
     "scientific"),
    ("cinematic_3d_camera_orbit", anim_camera_orbit, "camera angle",
     "cinematic"),
    ("correction_amplitude_3d", anim_correction_amplitude,
     "correction amplitude", "scientific"),
    ("quadrature_convergence", anim_quadrature, "quadrature refinement",
     "scientific"),
    ("objective_decomposition", anim_decomposition, "physical time",
     "scientific"),
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rec = records.resolve(a.records)
    out = records.ensure_dir(a.out)
    man = O.read_json(os.path.join(rec, "geometry_bridge_manifest.json"))
    mman = O.read_json(os.path.join(rec, "meridional_manifest.json"))
    frames = int(mman["frames"])

    entries = []
    for stem, fn, axis, style in BUILDERS:
        theme = themes.get("cinematic" if style == "cinematic"
                           else "scientific")
        themes.apply(theme)
        paths = fn(rec, out, man, theme, frames)
        entry = {
            "name": stem,
            "frames": paths["frames"],
            "resolution": paths["resolution"],
            "blank_frames": paths["blanks"],
            "contact_sheet": (os.path.basename(paths["sheet"])
                              if paths["sheet"] else None),
            "frame_axis": axis,
            "frame_axis_is_physical_time": axis == "physical time",
            "style": style,
            "gif_written": bool(paths["gif"]),
            "mp4_written": bool(paths["mp4"]),
            "axisymmetric_reconstruction": "3d" in stem or "shell" in stem,
            "reconstruction_statement": (
                O.RECONSTRUCTION_NOTE if ("3d" in stem or "shell" in stem)
                else "two-dimensional production data, drawn directly"),
            "renderer_computes_physics": False,
            "fixed_scales": True,
        }
        if not paths["mp4"]:
            entry["mp4_unwritten_reason"] = (
                "no ffmpeg encoder is available in this environment"
                if render.find_ffmpeg() is None else "encoder call failed")
        entries.append(entry)
        print("wrote", stem, "gif" if entry["gif_written"] else "",
              "mp4" if entry["mp4_written"] else "(no mp4)")

    with open(os.path.join(out, "animation_manifest.json"), "w",
              encoding="utf-8") as f:
        json.dump({"animations": entries, "count": len(entries),
                   "scope": O.SCOPE,
                   "index_ordering": O.INDEX_ORDER,
                   "frame_axis_note":
                       "A camera orbit is not physical time. Each entry names "
                       "its own frame axis.",
                   "not_an_energy_claim": O.NOT_AN_ENERGY_CLAIM}, f, indent=2)
    print("%d animations" % len(entries))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
