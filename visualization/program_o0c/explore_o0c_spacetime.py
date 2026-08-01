#!/usr/bin/env python3
"""The O0C spacetime explorer: a worldtube in (x, rho, t) with synchronized
spatial hypersurfaces.

    python explore_o0c_spacetime.py --records <o0c_energy dir> --out <gallery>

WHAT THIS IS, SAID PRECISELY. The main view is a SPACETIME WORLDTUBE: two
spatial axes (x and the transverse radius rho) and one TIME axis, drawn as a
three-dimensional block. It is not a view of four spatial dimensions, and
nothing here projects one. The secondary view is a spatial hypersurface at one
instant, revolved about the x axis from the same production data.

The third axis of the main view is COORDINATE TIME. That is the single most
important label on these products, because a time axis rendered into a screen
direction looks exactly like a spatial one. Every frame, legend and sidecar says
which axis is which.

FOUR THINGS CAN VARY, AND ONLY ONE OF THEM IS TIME:

    physical time         the slicing plane sweeps a FIXED worldtube
    camera angle          nothing physical changes; the camera moves
    correction amplitude  frozen time, different coordinate vectors
    quadrature refinement frozen everything, the integration is refined

The first two are produced here; the second two already exist in the O0C
animation set. All four are distinguished in the manifest.

THE CAUSAL OVERLAY IS LOCAL ALGEBRA AND NOTHING MORE. The null cones are read
from null_cone_samples.csv, which production computed as
gamma_ij (v^i + beta^i)(v^j + beta^j) = alpha^2 at each event. Cone tipping is a
property of the metric coefficients at a single point. No geodesic is
integrated, and no statement is made about causal reachability, travel time or
faster-than-light motion.
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
from matplotlib.collections import PolyCollection, LineCollection  # noqa: E402
from matplotlib.lines import Line2D                              # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

import o0c_common as O                                           # noqa: E402
import plot_o0c_energy as P                                      # noqa: E402
from common import records, render                               # noqa: E402
from themes import themes                                        # noqa: E402

# The encoding, fixed once and used by both views so the eye can carry a shape
# from the worldtube to the hypersurface without relearning the palette.
WALL_COLOUR = "#9fb3c8"      # neutral translucent
ENERGY_COLOUR = "#d1495b"    # strongly contrasting
NEC_COLOUR = "#b07cff"       # visually distinct from negative rho
CABIN_COLOUR = "#7ed957"
CENTRELINE_COLOUR = "#ffb454"
CONE_COLOUR = "#ffd166"
CURVATURE_COLOUR = "#5ec8b0"

WORLDTUBE_NOTE = (
    "spacetime worldtube: two spatial axes and one TIME axis, not a view of "
    "four spatial dimensions"
)

CONE_NOTE = (
    "null cones are local metric algebra at each event: no geodesic is "
    "integrated and no claim is made about causal reachability, travel time or "
    "faster-than-light motion"
)


# --- worldtube geometry -------------------------------------------------------

def _stack(rec, frames):
    """The (x, rho, t) production block, plus its axes."""
    planes = [O.meridional(rec, k) for k in range(frames)]
    nx, nr = planes[0]["nx"], planes[0]["nr"]
    xs = [planes[0]["x"][i][0] for i in range(nx)]
    rs = [planes[0]["rho"][0][j] for j in range(nr)]
    return planes, xs, rs


def _tube_quads(planes, times, field, level, rs, v_s):
    """Loft the per-time isocontours into quadrilaterals in (x, rho, t).

    THE X GRID IS COMOVING, WHICH IS EASY TO GET WRONG. Production writes each
    meridional plane on a grid whose x values already include the bubble offset
    v_s t, so frame k's column i is at a DIFFERENT x from frame 0's column i.
    Using frame 0's axis for every frame draws a straight worldtube next to a
    leaning centreline -- visibly inconsistent, and wrong: the structure
    translates, so its worldtube must lean by exactly v_s.

    So each frame's contour is extracted on ITS OWN x axis, and the
    frame-to-frame correspondence is made in the comoving coordinate x - v_s t,
    where the curve is nearly stationary and nearest-midpoint matching is
    unambiguous.
    """
    contours = []
    for k, pl in enumerate(planes):
        xs_k = [pl["x"][i][0] for i in range(pl["nx"])]
        contours.append(O.isocontour(pl[field], xs_k, rs, level))

    quads = []
    for k in range(len(planes) - 1):
        a, b = contours[k], contours[k + 1]
        if not a or not b:
            continue
        t0, t1 = times[k], times[k + 1]
        off0, off1 = v_s * t0, v_s * t1
        for (p0, p1) in a:
            mid = (0.5 * (p0[0] + p1[0]) - off0, 0.5 * (p0[1] + p1[1]))
            best, bd = None, 1e30
            for (q0, q1) in b:
                m2 = (0.5 * (q0[0] + q1[0]) - off1, 0.5 * (q0[1] + q1[1]))
                d = (m2[0] - mid[0]) ** 2 + (m2[1] - mid[1]) ** 2
                if d < bd:
                    bd, best = d, (q0, q1)
            if best is None:
                continue
            q0, q1 = best
            quads.append(((p0[0], p0[1], t0), (p1[0], p1[1], t0),
                          (q1[0], q1[1], t1), (q0[0], q0[1], t1)))
    return quads


def _cabin_tube(radius, times, v_s, n=26):
    """The protected cabin's worldtube.

    In comoving coordinates the cabin is w = xi^2 + s <= w_cabin, which in
    metres is (x - x_c)^2 + rho^2 <= r_cabin^2: a DISC in the meridional
    half-plane, centred on the moving bubble centre. Swept through time it is a
    leaning cylinder, and it must lean by the same v_s as everything else.
    """
    def rim(t):
        xc = v_s * t
        return [(xc + radius * math.cos(math.pi * i / n),
                 radius * math.sin(math.pi * i / n), t) for i in range(n + 1)]
    quads = []
    for k in range(len(times) - 1):
        a, b = rim(times[k]), rim(times[k + 1])
        for i in range(n):
            quads.append((a[i], a[i + 1], b[i + 1], b[i]))
    return quads


def _project_st(pts, azim, elev, t_scale, distance=30.0):
    """Project (x, rho, t) with the TIME axis scaled independently.

    Time is given its own scale because x and rho are metres and t is seconds
    (geometrised, so also metres, but the ranges differ by a factor of a few and
    an unscaled block reads as a thin sliver).  The scale is printed on every
    figure so nobody mistakes the aspect for physics.
    """
    return O.project([(p[0], p[1], p[2] * t_scale) for p in pts], azim, elev,
                     distance)


def _draw_quads(ax, quads, colour, azim, elev, t_scale, alpha, lw=0.2,
                edge=(1, 1, 1, 0.05)):
    if not quads:
        return
    keyed = []
    for q in quads:
        pr = _project_st(q, azim, elev, t_scale)
        keyed.append((sum(p[2] for p in pr) / len(pr),
                      [(p[0], p[1]) for p in pr]))
    keyed.sort(key=lambda t: -t[0])
    depths = [d for d, _ in keyed]
    lo, hi = min(depths), max(depths)
    import matplotlib.colors as mc
    base = mc.to_rgb(colour)
    shades = [0.45 + 0.55 * (0.0 if hi == lo else (hi - d) / (hi - lo))
              for d in depths]
    pc = PolyCollection([p for _, p in keyed], closed=True, linewidths=lw,
                        edgecolors=edge)
    pc.set_facecolors([(base[0] * s, base[1] * s, base[2] * s, alpha)
                       for s in shades])
    ax.add_collection(pc)


def _draw_polyline(ax, pts, colour, azim, elev, t_scale, lw=1.8, ls="-"):
    pr = _project_st(pts, azim, elev, t_scale)
    ax.plot([p[0] for p in pr], [p[1] for p in pr], color=colour, lw=lw, ls=ls,
            zorder=6)


def _slice_plane(t, xr, rr, n=1):
    """The moving time-slice plane: a constant-t rectangle in the block."""
    _ = n
    return [((xr[0], rr[0], t), (xr[1], rr[0], t),
             (xr[1], rr[1], t), (xr[0], rr[1], t))]


def _cone_wireframe(centre_vx, radius, apex, dt, n=12):
    """A local null cone as a sparse wireframe, from recorded numbers only.

    The cross-section at coordinate time apex_t + dt is the circle of radius
    `radius * dt` centred at x = apex_x + centre_vx * dt.  Both numbers come
    straight from null_cone_samples.csv; this function does no metric algebra.
    """
    ax0, ar0, at0 = apex
    cx = ax0 + centre_vx * dt
    rr = radius * dt
    rim = []
    for k in range(n + 1):
        a = 2.0 * math.pi * k / n
        rim.append((cx + rr * math.cos(a), ar0 + rr * math.sin(a), at0 + dt))
    lines = [rim]
    for k in range(0, n, 3):
        lines.append([(ax0, ar0, at0), rim[k]])
    return lines


# --- the explorer figure -------------------------------------------------------

def _worldtube_axes(ax, theme, xr, rr, tr, azim, elev, t_scale, label_axes=True):
    ax.set_facecolor(theme.background)
    ax.axis("off")
    ax.set_aspect("equal")
    if not label_axes:
        return
    origin = (xr[0], rr[0], tr[0])

    # Arrows reach 72% of each range so their labels land inside the canvas
    # instead of colliding with the header text.
    def _to(frac, lo, hi):
        return lo + frac * (hi - lo)

    arrows = (((_to(0.72, xr[0], xr[1]), rr[0], tr[0]), "x  [m]", "#8fa7bf"),
              ((xr[0], _to(0.72, rr[0], rr[1]), tr[0]), r"$\rho$  [m]",
               "#8fa7bf"),
              ((xr[0], rr[0], _to(0.72, tr[0], tr[1])),
               "t  [m, geometrised]  <- TIME AXIS", CENTRELINE_COLOUR))
    for tip, lab, col in arrows:
        pr = _project_st([origin, tip], azim, elev, t_scale)
        ax.annotate("", xy=(pr[1][0], pr[1][1]), xytext=(pr[0][0], pr[0][1]),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.3))
        ax.text(pr[1][0], pr[1][1], "  " + lab, fontsize=7.0, color=col,
                va="center")


def _legend(fig, theme, x=0.012, y=0.30, cones=True):
    handles = [
        Line2D([], [], lw=6, color=WALL_COLOUR, alpha=0.5,
               label="bubble wall (shift isosurface)"),
        Line2D([], [], lw=6, color=ENERGY_COLOUR, alpha=0.6,
               label=r"negative Eulerian energy  $\rho<0$"),
        Line2D([], [], lw=6, color=NEC_COLOUR, alpha=0.6,
               label=r"NEC violation  $N_{\min}<0$"),
        Line2D([], [], lw=6, color=CABIN_COLOUR, alpha=0.8,
               label="protected cabin (F1W, support-excluded)"),
        Line2D([], [], lw=2, color=CENTRELINE_COLOUR,
               label="bubble-centre trajectory"),
        Line2D([], [], lw=1.4, color=CURVATURE_COLOUR, ls="--",
               label="Kretschmann contour"),
    ]
    if cones:
        handles.append(Line2D([], [], lw=1.2, color=CONE_COLOUR,
                              label="local null cone (metric algebra only)"))
    leg = fig.legend(handles=handles, loc="center left",
                     bbox_to_anchor=(x, y), fontsize=6.8,
                     facecolor=theme.background, edgecolor=theme.grid,
                     labelcolor=theme.foreground, framealpha=0.9)
    leg.set_zorder(20)


def build_explorer(rec, theme, frames, *, slice_frame, azim=0.86, elev=0.30,
                   show_cones=True, show_curvature=True, title=None,
                   axis_label=None):
    planes, xs, rs = _stack(rec, frames)
    scales = O.read_json(os.path.join(rec, "field_scales.json"))
    mman = O.read_json(os.path.join(rec, "meridional_manifest.json"))
    cones = O.read_csv(os.path.join(rec, "null_cone_samples.csv"))
    v_s = float(mman["velocity"])
    times = [-3.0 + 6.0 * k / (frames - 1) for k in range(frames)]
    iso = scales["isovalues"]

    all_x = [pl["x"][i][0] for pl in planes for i in range(pl["nx"])]
    xr = (min(all_x), max(all_x))
    rr = (0.0, max(rs))
    tr = (times[0], times[-1])
    # Time gets its own scale so the block is legible; printed on the figure.
    t_scale = (xr[1] - xr[0]) / max(1e-9, (tr[1] - tr[0])) * 0.85

    fig = plt.figure(figsize=(17.4, 8.0), facecolor=theme.background)
    axw = fig.add_axes([0.155, 0.075, 0.47, 0.80])
    axs_ = fig.add_axes([0.655, 0.135, 0.33, 0.68])

    _worldtube_axes(axw, theme, xr, rr, tr, azim, elev, t_scale)

    # Back-to-front: wall, then energy, then NEC, then cabin, then overlays.
    _draw_quads(axw, _tube_quads(planes, times, "beta_x", iso["wall_shift"],
                                 rs, v_s),
                WALL_COLOUR, azim, elev, t_scale, alpha=0.17)
    _draw_quads(axw, _tube_quads(planes, times, "energy_density",
                                 iso["negative_energy"], rs, v_s),
                ENERGY_COLOUR, azim, elev, t_scale, alpha=0.30)
    _draw_quads(axw, _tube_quads(planes, times, "n_min", iso["nec_violation"],
                                 rs, v_s),
                NEC_COLOUR, azim, elev, t_scale, alpha=0.20)
    cabin_r = scales["protected_cabin_radius_comoving"]
    _draw_quads(axw, _cabin_tube(cabin_r, times, v_s), CABIN_COLOUR, azim, elev,
                t_scale, alpha=0.50, lw=0.3)

    if show_curvature:
        for k in range(0, frames, 3):
            xs_k = [planes[k]["x"][i][0] for i in range(planes[k]["nx"])]
            segs = O.isocontour(planes[k]["kretschmann"], xs_k, rs,
                                iso["curvature_shell"])
            pts = []
            for (a, b) in segs:
                pr = _project_st([(a[0], a[1], times[k]),
                                  (b[0], b[1], times[k])], azim, elev, t_scale)
                pts.append([(pr[0][0], pr[0][1]), (pr[1][0], pr[1][1])])
            if pts:
                axw.add_collection(LineCollection(
                    pts, colors=CURVATURE_COLOUR, linewidths=0.7, alpha=0.55))

    # The bubble-centre trajectory: a straight line of slope v_s in (x, t).
    _draw_polyline(axw, [(v_s * t, 0.0, t) for t in times], CENTRELINE_COLOUR,
                   azim, elev, t_scale, lw=2.2)

    # The moving slice plane.
    t_sel = times[slice_frame]
    _draw_quads(axw, _slice_plane(t_sel, xr, rr), "#ffffff", azim, elev,
                t_scale, alpha=0.10, lw=1.0, edge=(1, 1, 1, 0.55))

    # Sparse null cones along the centreline, at the selected time.
    if show_cones:
        sel = [c for c in cones if int(c["frame"]) == slice_frame]
        dt = 0.9
        for c in sel:
            if c["site"].startswith("centerline") or "shoulder" in c["site"]:
                apex = (float(c["x"]), float(c["rho"]), t_sel)
                for line in _cone_wireframe(float(c["cone_center_vx"]),
                                            float(c["cone_radius"]), apex, dt):
                    pr = _project_st(line, azim, elev, t_scale)
                    axw.plot([p[0] for p in pr], [p[1] for p in pr],
                             color=CONE_COLOUR, lw=0.7, alpha=0.75, zorder=7)

    axw.autoscale_view()
    axw.set_xlim(axw.get_xlim())
    axw.set_ylim(axw.get_ylim())

    # --- the synchronized spatial hypersurface -----------------------------
    mer = planes[slice_frame]
    P._setup_3d_axes(axs_, theme, extent=6.2)
    for field, key, colour, alpha in (
            ("beta_x", "wall_shift", WALL_COLOUR, 0.14),
            ("energy_density", "negative_energy", ENERGY_COLOUR, 0.30),
            ("n_min", "nec_violation", NEC_COLOUR, 0.18)):
        faces, _ = P._shell_faces(mer, field, iso[key])
        P._draw_shell(axs_, faces, theme, colour, 0.7, 0.32, alpha=alpha)
    P._draw_shell(axs_, P._cabin_faces(cabin_r), theme, CABIN_COLOUR, 0.7,
                  0.32, alpha=0.6, lw=0.15)

    # Sparse ADM shift arrows on this hypersurface.
    sv = O.read_csv(os.path.join(rec, "shift_vectors_sparse.csv"))
    tails, tips = [], []
    for r in sv:
        b = float(r["beta_x"])
        if abs(b) < 0.06:
            continue
        p = (float(r["x"]), float(r["y"]), float(r["z"]))
        tails.append(p)
        tips.append((p[0] + 2.2 * b, p[1], p[2]))
    for a, b in zip(O.project(tails, 0.7, 0.32), O.project(tips, 0.7, 0.32)):
        axs_.annotate("", xy=(b[0], b[1]), xytext=(a[0], a[1]),
                      arrowprops=dict(arrowstyle="-|>", color=CENTRELINE_COLOUR,
                                      lw=0.7, alpha=0.7))

    # --- labels --------------------------------------------------------------
    fig.text(0.5, 0.968, title or
             "O0C spacetime explorer: worldtube in (x, rho, t) with a "
             "synchronized spatial hypersurface",
             ha="center", fontsize=11, color=theme.foreground)
    fig.text(0.39, 0.928, WORLDTUBE_NOTE, ha="center", fontsize=7.4,
             color=theme.uncertainty_color)
    fig.text(0.39, 0.898,
             "time axis scaled by %.2f for legibility; x and rho are to scale "
             "with each other" % t_scale, ha="center", fontsize=6.6,
             color=theme.grid)
    fig.text(0.82, 0.928, O.RECONSTRUCTION_NOTE, ha="center", fontsize=7.2,
             color=theme.uncertainty_color)
    fig.text(0.82, 0.898, "spatial hypersurface at t = %+.2f m" % t_sel,
             ha="center", fontsize=7.4, color=theme.foreground)
    fig.text(0.5, 0.045, axis_label or
             "FRAME AXIS: PHYSICAL TIME -- the slicing plane sweeps a fixed "
             "worldtube; the worldtube itself never moves",
             ha="center", fontsize=7.6, color=theme.uncertainty_color)
    fig.text(0.5, 0.016, CONE_NOTE + ". " + O.NOT_AN_ENERGY_CLAIM,
             ha="center", fontsize=6.3, color=theme.grid, wrap=True)
    fig.text(0.988, 0.964, "NOT READY FOR OPTIMIZATION", ha="right",
             fontsize=6.8, color=theme.uncertainty_color)
    _legend(fig, theme, cones=show_cones)
    return fig


def _render(fig):
    fig.canvas.draw()
    import numpy as np
    arr = np.asarray(fig.canvas.buffer_rgba()).copy()
    plt.close(fig)
    return arr


def _write(imgs, out, stem, fps):
    import numpy as np
    from PIL import Image
    frame_dir = records.ensure_dir(os.path.join(out, "_frames_" + stem))
    paths = []
    for k, arr in enumerate(imgs):
        fp = os.path.join(frame_dir, "o0c_%05d.png" % k)
        Image.fromarray(np.asarray(arr)).save(fp)
        paths.append(fp)
    gif = render.write_gif(paths, os.path.join(out, stem + ".gif"), fps=fps,
                           max_frames=len(paths), max_width=1000)
    mp4 = render.write_mp4(paths, os.path.join(out, stem + ".mp4"), fps=fps,
                           frame_dir=frame_dir)
    sheet = render.contact_sheet(
        paths, os.path.join(out, "contact_sheets", stem + ".png"),
        tiles=min(9, len(paths)), label_prefix="frame")
    res = "%dx%d" % (imgs[0].shape[1], imgs[0].shape[0])
    for q in paths:
        os.remove(q)
    os.rmdir(frame_dir)
    return {"gif": gif, "mp4": mp4, "sheet": sheet, "frames": len(paths),
            "resolution": res}


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
    theme = themes.get("scientific")
    themes.apply(theme)

    inputs = [os.path.join(rec, "meridional_field_%03d.csv" % k)
              for k in range(frames)]
    inputs += [os.path.join(rec, "field_scales.json"),
               os.path.join(rec, "null_cone_samples.csv"),
               os.path.join(rec, "shift_vectors_sparse.csv")]

    # --- the still ----------------------------------------------------------
    fig = build_explorer(rec, theme, frames, slice_frame=frames // 2)
    O.save_figure(fig, out, "spacetime_explorer_4d", manifest=man,
                  inputs=inputs, theme=theme, reconstruction=True,
                  camera="perspective; main view is a (x, rho, t) spacetime "
                         "worldtube, secondary view is a revolved spatial "
                         "hypersurface",
                  extra={"view_kind": "spacetime_worldtube",
                         "third_axis": "coordinate time",
                         "is_four_spatial_dimensions": False,
                         "worldtube_note": WORLDTUBE_NOTE,
                         "null_cone_note": CONE_NOTE,
                         "geodesics_integrated": False,
                         "causal_reachability_claimed": False},
                  frame_axis="physical time")
    print("wrote spacetime_explorer_4d.png")

    # A cone-focused still, so the causal overlay can be read on its own.
    cones = O.read_csv(os.path.join(rec, "null_cone_samples.csv"))
    # THE CENTRELINE ONLY, SORTED BY x. The record also carries off-axis sites,
    # several of which share an x with a centreline site; plotting all of them
    # against x produces a zigzag that reads as structure and is nothing but
    # unsorted rows.
    sel = sorted((c for c in cones
                  if int(c["frame"]) == frames // 2
                  and c["site"].startswith("centerline")),
                 key=lambda c: float(c["x"]))
    fig2, ax = plt.subplots(figsize=(12.6, 5.4), facecolor=theme.background)
    xs = [float(c["x"]) for c in sel]
    fwd = [float(c["slope_forward"]) for c in sel]
    back = [float(c["slope_backward"]) for c in sel]
    ax.fill_between(xs, back, fwd, color=CONE_COLOUR, alpha=0.22,
                    label="local null cone opening  dx/dt")
    ax.plot(xs, fwd, color=CONE_COLOUR, lw=1.6)
    ax.plot(xs, back, color=CONE_COLOUR, lw=1.6)
    ax.plot(xs, [float(c["cone_center_vx"]) for c in sel], color=CENTRELINE_COLOUR,
            lw=1.8, label=r"cone centre  $-\beta^x = v_s f$")
    ax.axhline(0.0, color=theme.foreground, lw=0.9, ls=":")
    ax.text(min(xs), 0.06, "  v = 0: the coordinate-stationary observer",
            fontsize=6.8, color=theme.foreground)
    O.style_axes(ax, theme, "x  [m]", r"$dx/dt$")
    ax.legend(fontsize=7.4, facecolor=theme.background, edgecolor=theme.grid,
              labelcolor=theme.foreground)
    fig2.text(0.5, 0.962,
              "Local null-cone tipping along the axis: the bubble TRANSLATES "
              "the cone, it does not open it", ha="center", fontsize=10,
              color=theme.foreground)
    fig2.text(0.5, 0.915,
              r"$\gamma_{ij}(v^i+\beta^i)(v^j+\beta^j)=\alpha^2$ -- the cone "
              r"half-width stays $\alpha$ while its centre moves to $-\beta^x$",
              ha="center", fontsize=7.6, color=theme.uncertainty_color)
    O.footer(fig2, theme, CONE_NOTE + ". " + O.NOT_AN_ENERGY_CLAIM)
    O.readiness_stamp(fig2, theme)
    fig2.tight_layout(rect=[0, 0.05, 1, 0.89])
    O.save_figure(fig2, out, "null_cone_tipping", manifest=man,
                  inputs=[os.path.join(rec, "null_cone_samples.csv")],
                  theme=theme,
                  extra={"null_cone_note": CONE_NOTE,
                         "geodesics_integrated": False,
                         "causal_reachability_claimed": False,
                         "travel_time_claimed": False,
                         "ftl_motion_claimed": False},
                  plotted_header=["x", "slope_forward", "slope_backward",
                                  "cone_center"],
                  plotted_rows=[[c["x"], c["slope_forward"],
                                 c["slope_backward"], c["cone_center_vx"]]
                                for c in sel])
    print("wrote null_cone_tipping.png")

    entries = []

    # --- physical time: the slice sweeps a FIXED worldtube -------------------
    imgs = [_render(build_explorer(rec, theme, frames, slice_frame=k))
            for k in range(frames)]
    p = _write(imgs, out, "spacetime_slice_sweep", fps=4)
    entries.append({
        "name": "spacetime_slice_sweep", "frame_axis": "physical time",
        "frame_axis_is_physical_time": True, "style": "scientific",
        "what_moves": "the slicing plane and the synchronized hypersurface; "
                      "the worldtube itself is FIXED",
        "camera_moves": False,
    })
    print("wrote spacetime_slice_sweep")

    # --- camera angle: nothing physical changes ------------------------------
    cin = themes.get("cinematic")
    themes.apply(cin)
    n = 18
    imgs = []
    for k in range(n):
        az = 0.45 + 2.0 * math.pi * k / n
        imgs.append(_render(build_explorer(
            rec, cin, frames, slice_frame=frames // 2, azim=az,
            elev=0.24 + 0.12 * math.sin(2.0 * math.pi * k / n),
            title="Cinematic orbit of the spacetime worldtube",
            axis_label="FRAME AXIS: CAMERA ANGLE -- NOT PHYSICAL TIME. The "
                       "spacetime is frozen; only the camera moves. "
                       "azimuth %.2f rad" % az)))
    p2 = _write(imgs, out, "spacetime_worldtube_orbit", fps=12)
    entries.append({
        "name": "spacetime_worldtube_orbit", "frame_axis": "camera angle",
        "frame_axis_is_physical_time": False, "style": "cinematic",
        "what_moves": "the camera only; the worldtube and the selected "
                      "hypersurface are frozen",
        "camera_moves": True,
    })
    print("wrote spacetime_worldtube_orbit")

    for e, paths in zip(entries, (p, p2)):
        e.update({"gif_written": bool(paths["gif"]),
                  "mp4_written": bool(paths["mp4"]),
                  "frames": paths["frames"],
                  "resolution": paths["resolution"],
                  "contact_sheet": (os.path.basename(paths["sheet"])
                                    if paths["sheet"] else None),
                  "axisymmetric_reconstruction": True,
                  "reconstruction_statement": O.RECONSTRUCTION_NOTE,
                  "view_kind": "spacetime_worldtube",
                  "third_axis": "coordinate time",
                  "is_four_spatial_dimensions": False,
                  "renderer_computes_physics": False,
                  "fixed_scales": True})
        if not paths["mp4"]:
            e["mp4_unwritten_reason"] = (
                "no ffmpeg encoder is available in this environment"
                if render.find_ffmpeg() is None else "encoder call failed")

    with open(os.path.join(out, "explorer_manifest.json"), "w",
              encoding="utf-8") as f:
        json.dump({
            "figures": ["spacetime_explorer_4d", "null_cone_tipping"],
            "animations": entries,
            "view_kind": "spacetime_worldtube",
            "worldtube_note": WORLDTUBE_NOTE,
            "null_cone_note": CONE_NOTE,
            "is_four_spatial_dimensions": False,
            "distinguished_frame_axes": [
                "physical time", "camera angle", "correction amplitude",
                "quadrature refinement"],
            "geodesics_integrated": False,
            "causal_reachability_claimed": False,
            "travel_time_claimed": False,
            "ftl_motion_claimed": False,
            "renderer_computes_physics": False,
            "scope": O.SCOPE,
            "not_an_energy_claim": O.NOT_AN_ENERGY_CLAIM}, f, indent=2)
    print("2 figures, %d animations" % len(entries))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
