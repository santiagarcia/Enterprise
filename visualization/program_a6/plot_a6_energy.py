#!/usr/bin/env python3
"""The twelve A6 observer-energy figures, drawn from saved records only.

    python plot_a6_energy.py --sweep <a6 sweep dir> --out <gallery dir>

Every panel here is an OBSERVER-MEASURED quantity, which is what distinguishes
A6 from A5: A5's figures showed coordinate-basis Einstein-tensor components and
said on every frame that they were not energy densities. These are energy
densities -- for a NAMED, NORMALISED observer, stated on the figure, with the
certification status of every sample carried alongside the value.

The two distinctions the figures exist to keep visible:

  * a negative MEASURED value is not a certified negative value. A sample whose
    magnitude sits inside three numerical error envelopes is reported as
    indeterminate and drawn as indeterminate, never as a definite sign.

  * the coordinate-stationary observer is not the Eulerian observer. It ceases
    to exist at g_00 = 0 and the records leave its energy density EMPTY there.
    No figure draws it past that point.
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
from matplotlib.colors import ListedColormap, BoundaryNorm       # noqa: E402
from matplotlib.lines import Line2D                              # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

import a6_common as A6                                           # noqa: E402
from common import anomaly, records                              # noqa: E402
from themes import themes                                        # noqa: E402


# --- Small drawing helpers ---------------------------------------------------

def _style_axes(ax, theme, xlabel="x", ylabel="y"):
    ax.set_xlabel(xlabel, fontsize=8, color=theme.foreground)
    ax.set_ylabel(ylabel, fontsize=8, color=theme.foreground)
    ax.tick_params(labelsize=7, colors=theme.foreground)
    ax.set_facecolor(theme.background)
    for s in ax.spines.values():
        s.set_color(theme.grid)


def _field(ax, fig, grid2d, extent, limit, theme, *, label, xs=None, ys=None,
           zero_contour=False, radius=None):
    im = ax.imshow(grid2d, origin="lower", extent=extent, aspect="equal",
                   cmap="RdBu_r", vmin=-limit, vmax=limit,
                   interpolation="nearest")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=7, colors=theme.foreground)
    cb.outline.set_visible(False)
    if radius:
        ax.add_patch(plt.Circle((0, 0), radius, fill=False,
                                color=theme.foreground, linewidth=0.7,
                                linestyle="--", alpha=0.55))
    if zero_contour and xs is not None:
        flat = [v for row in grid2d for v in row]
        # Drawn only where the SAVED grid brackets zero. Below the onset there
        # is correctly nothing to draw, and inventing a contour there would be
        # drawing a level the data does not reach.
        if min(flat) < 0.0 < max(flat):
            ax.contour(xs, ys, grid2d, levels=[0.0],
                       colors=[theme.accent_two], linewidths=1.4)
    ax.set_title(label, fontsize=9.5, color=theme.foreground)
    _style_axes(ax, theme)
    return im


def _status_panel(ax, fig, grid2d, extent, theme, *, label, radius=None):
    cmap = ListedColormap([A6.STATUS_COLOR[s] for s in A6.STATUS_ORDER])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    ax.imshow(grid2d, origin="lower", extent=extent, aspect="equal",
              cmap=cmap, norm=norm, interpolation="nearest")
    if radius:
        ax.add_patch(plt.Circle((0, 0), radius, fill=False,
                                color=theme.foreground, linewidth=0.7,
                                linestyle="--", alpha=0.55))
    ax.set_title(label, fontsize=9.5, color=theme.foreground)
    _style_axes(ax, theme)


def _status_legend(fig, theme, labels, y=0.012):
    handles = [Line2D([], [], marker="s", linestyle="none", markersize=7,
                      markerfacecolor=A6.STATUS_COLOR[s],
                      markeredgecolor=theme.grid, label=labels[s])
               for s in A6.STATUS_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7.2,
               frameon=False, bbox_to_anchor=(0.5, y),
               labelcolor=theme.foreground)


def _footer(fig, theme, text, y=0.002):
    fig.text(0.5, y, text, ha="center", va="bottom", fontsize=6.2,
             color=theme.foreground, alpha=0.72, wrap=True)


# --- The figures -------------------------------------------------------------

def fig_energy_density(sweep, out, meta, theme, frames, ref, sweep_rows,
                       limits):
    rows = frames[ref]
    last = len(frames) - 1
    xs, ys = A6.axes_of(rows)
    extent = [min(xs), max(xs), min(ys), max(ys)]
    rho = A6.to_grid(rows, "rho", xs, ys)
    sign = A6.status_grid(rows, "rho_sign_status", xs, ys)
    sign_last = A6.status_grid(frames[last], "rho_sign_status", xs, ys)
    s = sweep_rows[ref]
    sl = sweep_rows[last]

    # The reference frame's OWN symmetric limit. A static figure has one frame,
    # so there is nothing to hold fixed across; using the sweep-wide limit here
    # would render this frame nearly flat because |rho| grows by an order of
    # magnitude by v_s = 1.6. The animations, where a drifting scale really
    # would mislead, use the global limit instead and say so.
    local = A6.global_symmetric_limit([rows], "rho")

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.0),
                             facecolor=theme.background)
    _field(axes[0], fig, rho, extent, local, theme,
           label=r"Eulerian $\rho = T_{\hat a\hat b}\,n^{\hat a}n^{\hat b}$"
                 "\n$v_s$ = %.1f" % float(s["v_s"]),
           xs=xs, ys=ys, zero_contour=True, radius=float(meta["radius"]))
    mx, my = float(s["slice_rho_min_x"]), float(s["slice_rho_min_y"])
    axes[0].plot([mx], [my], marker="x", markersize=9, markeredgewidth=1.8,
                 color=theme.accent_two, linestyle="none")
    axes[0].annotate("min %.3e" % float(s["slice_rho_min"]), (mx, my),
                     textcoords="offset points", xytext=(8, 6), fontsize=7,
                     color=theme.accent_two)

    _status_panel(axes[1], fig, sign, extent, theme,
                  label="sign certification, SAME samples, $v_s$ = %.1f"
                        % float(s["v_s"]),
                  radius=float(meta["radius"]))
    _status_panel(axes[2], fig, sign_last, extent, theme,
                  label="sign certification at $v_s$ = %.1f"
                        % float(sl["v_s"]),
                  radius=float(meta["radius"]))

    n_cert = _count(rows, "rho_sign_status", A6.CERTIFIED_VIOLATED)
    n_ind = _count(rows, "rho_sign_status", A6.INDETERMINATE)
    n_cert_last = _count(frames[last], "rho_sign_status",
                         A6.CERTIFIED_VIOLATED)

    fig.suptitle("A6  Eulerian energy density, Alcubierre bubble     "
                 "z = 0 slice", fontsize=12, color=theme.foreground, x=0.010,
                 ha="left", y=0.985)
    fig.text(0.010, 0.930,
             "A NEGATIVE MEASURED VALUE IS NOT A CERTIFIED NEGATIVE VALUE. At "
             "$v_s$ = %.1f every one of the %d slice samples is measured but "
             "NOT certified: %d certified negative, %d indeterminate. At "
             "$v_s$ = %.1f, %d samples are certified negative."
             % (float(s["v_s"]), len(rows), n_cert, n_ind, float(sl["v_s"]),
                n_cert_last),
             fontsize=8, color=theme.uncertainty_color, ha="left")
    fig.text(0.010, 0.895,
             "numerical error estimate %.3e   safety factor %s   "
             "slice minimum %.4e at (%.3f, %.3f, 0)   frame %s status: %s"
             % (float(s["numerical_error_estimate"]),
                meta["certification_safety_factor"],
                float(s["slice_rho_min"]), mx, my, ref, s["rho_sign_status"]),
             fontsize=7.6, color=theme.foreground, ha="left", alpha=0.9)
    _status_legend(fig, theme, A6.SIGN_LABEL, y=0.028)
    _footer(fig, theme, "Dashed circle: nominal bubble wall R = %s. Orange "
                        "contour: rho = 0 where the saved grid brackets it. %s"
            % (meta["radius"], A6.SCOPE))
    fig.tight_layout(rect=(0, 0.085, 1, 0.875))

    plotted = [[r["x"], r["y"], r["rho"], r["rho_sign_status"],
                r["rho_certification_margin"], r["numerical_error_estimate"]]
               for r in rows]
    return A6.save_figure(
        fig, out, "eulerian_energy_density", meta=meta,
        inputs=[os.path.join(sweep, "energy_frame_%03d.csv" % ref),
                os.path.join(sweep, "energy_frame_%03d.csv" % last),
                os.path.join(sweep, "energy_velocity_sweep.csv")],
        theme=theme, observer="Eulerian (ADM normal n^mu)",
        extra={"reference_frame": ref, "v_s": float(s["v_s"]),
               "comparison_frame": last,
               "colour_scale": "symmetric about zero, computed from the drawn "
                               "frame; the animations use the sweep-wide "
                               "limit instead",
               "colour_limit": local,
               "certified_negative_samples": n_cert,
               "indeterminate_samples": n_ind,
               "certified_negative_samples_at_final_frame": n_cert_last},
        plotted_header=["x", "y", "rho", "rho_sign_status",
                        "rho_certification_margin",
                        "numerical_error_estimate"],
        plotted_rows=plotted)


def fig_momentum_density(sweep, out, meta, theme, frames, ref, sweep_rows,
                         limits):
    rows = frames[ref]
    xs, ys = A6.axes_of(rows)
    extent = [min(xs), max(xs), min(ys), max(ys)]
    qx = A6.to_grid(rows, "q_x", xs, ys)
    qy = A6.to_grid(rows, "q_y", xs, ys)
    qz = A6.to_grid(rows, "q_z", xs, ys)
    s = sweep_rows[ref]
    # The drawn frame's own limit; see fig_energy_density for why a static
    # figure does not inherit the animation's sweep-wide scale.
    local = max(A6.global_symmetric_limit([rows], k)
                for k in ("q_x", "q_y", "q_z"))
    arrow = 0.35 / local if local > 0 else 1.0

    fig, axes = plt.subplots(1, 3, figsize=(15.4, 4.9),
                             facecolor=theme.background)
    for ax, g, lab in ((axes[0], qx, r"$q_x$"), (axes[1], qy, r"$q_y$"),
                       (axes[2], qz, r"$q_z$")):
        _field(ax, fig, g, extent, local, theme, label=lab,
               radius=float(meta["radius"]))
    # The saved components, drawn as a field of arrows. Positions and vectors
    # both come from the CSV; only the arrow scale is a rendering choice.
    step = 2
    ax = axes[0]
    for j in range(0, len(ys), step):
        for i in range(0, len(xs), step):
            ax.arrow(xs[i], ys[j], qx[j][i] * arrow,
                     qy[j][i] * arrow, width=0.008,
                     head_width=0.10, length_includes_head=True,
                     color=theme.foreground, alpha=0.65)

    fig.suptitle("A6  Eulerian momentum density $q_i = -T_{\\hat 0\\hat i}$ "
                 "at $v_s$ = %.1f     z = 0 slice" % float(s["v_s"]),
                 fontsize=12, color=theme.foreground, x=0.010, ha="left",
                 y=0.985)
    fig.text(0.010, 0.925,
             "Measured by the EULERIAN observer. Arrows on the left panel are "
             "the saved $(q_x, q_y)$ at every second sample, scaled by a fixed "
             "factor %.1f; the arrow length is a rendering choice, its "
             "direction is not." % arrow,
             fontsize=7.8, color=theme.foreground, ha="left", alpha=0.9)
    _footer(fig, theme, "Fixed colour limits across the whole velocity sweep. "
                        + A6.SCOPE)
    fig.tight_layout(rect=(0, 0.055, 1, 0.905))

    plotted = [[r["x"], r["y"], r["q_x"], r["q_y"], r["q_z"]] for r in rows]
    return A6.save_figure(
        fig, out, "momentum_density", meta=meta,
        inputs=[os.path.join(sweep, "energy_frame_%03d.csv" % ref)],
        theme=theme, observer="Eulerian (ADM normal n^mu)",
        extra={"reference_frame": ref, "v_s": float(s["v_s"]),
               "arrow_scale": arrow, "colour_limit": local,
               "colour_scale": "symmetric about zero, computed from the drawn "
                               "frame; the animations use the sweep-wide "
                               "limit instead",
               "sign_convention": "q_i = -T_{0i}"},
        plotted_header=["x", "y", "q_x", "q_y", "q_z"], plotted_rows=plotted)


def fig_principal_stresses(sweep, out, meta, theme, frames, ref, sweep_rows,
                           limits):
    rows = frames[ref]
    s = sweep_rows[ref]
    cut_y = float(s["slice_rho_min_y"])
    line = sorted((r for r in rows if float(r["y"]) == cut_y),
                  key=lambda r: float(r["x"]))
    if not line:
        raise records.RecordError("no slice row at y = %r" % cut_y)

    x = [float(r["x"]) for r in line]
    p1 = [float(r["principal_pressure_1"]) for r in line]
    p2 = [float(r["principal_pressure_2"]) for r in line]
    p3 = [float(r["principal_pressure_3"]) for r in line]
    rho = [float(r["rho"]) for r in line]

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8),
                             facecolor=theme.background)
    ax = axes[0]
    for vals, lab, col in ((p1, r"$p_1$ (smallest)", theme.accent),
                           (p2, r"$p_2$", theme.accent_two),
                           (p3, r"$p_3$ (largest)", "#7ed957")):
        ax.plot(x, vals, lw=1.5, color=col, label=lab)
    ax.plot(x, rho, lw=1.4, color=theme.uncertainty_color, linestyle="--",
            label=r"$\rho$")
    ax.axhline(0.0, color=theme.grid, lw=0.8)
    ax.legend(fontsize=7.5, frameon=False, labelcolor=theme.foreground)
    ax.set_title("principal pressures along y = %.3f, z = 0" % cut_y,
                 fontsize=9.5, color=theme.foreground)
    _style_axes(ax, theme, "x", "eigenvalue of $S_{ij}$")

    # rho + p_i: the combination every pointwise condition is governed by for a
    # diagonal source. Plotted from the saved columns, one addition per point.
    ax = axes[1]
    for vals, lab, col in ((p1, r"$\rho + p_1$", theme.accent),
                           (p2, r"$\rho + p_2$", theme.accent_two),
                           (p3, r"$\rho + p_3$", "#7ed957")):
        ax.plot(x, [a + b for a, b in zip(rho, vals)], lw=1.5, color=col,
                label=lab)
    ax.plot(x, [float(r["N_min"]) for r in line], lw=1.6,
            color=theme.uncertainty_color,
            label=r"$N_{\min}$ (saved, minimised over the sphere)")
    ax.axhline(0.0, color=theme.grid, lw=0.8)
    ax.legend(fontsize=7.5, frameon=False, labelcolor=theme.foreground)
    ax.set_title(r"$\rho + p_i$ against the saved $N_{\min}$",
                 fontsize=9.5, color=theme.foreground)
    _style_axes(ax, theme, "x", "value")

    fig.suptitle("A6  principal stresses in the Eulerian orthonormal frame, "
                 "$v_s$ = %.1f" % float(s["v_s"]),
                 fontsize=12, color=theme.foreground, x=0.010, ha="left",
                 y=0.985)
    fig.text(0.010, 0.925,
             "The right panel is NOT a derivation of $N_{\\min}$: it shows the "
             "axis-aligned combinations beside the value production actually "
             "minimised over the whole unit sphere. Where they differ, the "
             "worst null direction is off the coordinate axes.",
             fontsize=7.6, color=theme.foreground, ha="left", alpha=0.9)
    _footer(fig, theme, A6.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.905))

    plotted = [[r["x"], r["y"], r["principal_pressure_1"],
                r["principal_pressure_2"], r["principal_pressure_3"],
                r["rho"], r["N_min"]] for r in line]
    return A6.save_figure(
        fig, out, "principal_stresses", meta=meta,
        inputs=[os.path.join(sweep, "energy_frame_%03d.csv" % ref)],
        theme=theme, observer="Eulerian (ADM normal n^mu)",
        extra={"reference_frame": ref, "cut_y": cut_y},
        plotted_header=["x", "y", "principal_pressure_1",
                        "principal_pressure_2", "principal_pressure_3", "rho",
                        "N_min"], plotted_rows=plotted)


def fig_minimum_nec(sweep, out, meta, theme, frames, ref, sweep_rows, limits):
    rows = frames[ref]
    xs, ys = A6.axes_of(rows)
    extent = [min(xs), max(xs), min(ys), max(ys)]
    n = A6.to_grid(rows, "N_min", xs, ys)
    st = A6.status_grid(rows, "NEC_status", xs, ys)
    s = sweep_rows[ref]

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0),
                             facecolor=theme.background)
    local = A6.global_symmetric_limit([rows], "N_min")
    _field(axes[0], fig, n, extent, local, theme,
           label=r"$N_{\min} = \min_{|\mathbf{n}|=1}"
                 r"\,T_{\hat a\hat b}k^{\hat a}k^{\hat b}$",
           xs=xs, ys=ys, zero_contour=True, radius=float(meta["radius"]))
    mx, my = float(s["slice_N_min_x"]), float(s["slice_N_min_y"])
    axes[0].plot([mx], [my], marker="x", markersize=9, markeredgewidth=1.8,
                 color=theme.accent_two, linestyle="none")

    _status_panel(axes[1], fig, st, extent, theme,
                  label="NEC certification of the SAME samples",
                  radius=float(meta["radius"]))

    # The independent reproduction, read from the record rather than redone:
    # the exporter saved N(n_min) recomputed from scratch by production, so the
    # agreement can be quoted without this script evaluating N at all.
    worst = max(abs(float(r["N_min"]) - float(r["nec_direct_substitution"]))
                for r in rows)

    fig.suptitle("A6  minimum null-energy functional at $v_s$ = %.1f     "
                 "z = 0 slice" % float(s["v_s"]),
                 fontsize=12, color=theme.foreground, x=0.010, ha="left",
                 y=0.985)
    fig.text(0.010, 0.930,
             "slice minimum %.5e at (%.3f, %.3f, 0), worst direction "
             "(%.3f, %.3f, %.3f)     status %s"
             % (float(s["slice_N_min"]), mx, my,
                float(s["worst_direction_x"]), float(s["worst_direction_y"]),
                float(s["worst_direction_z"]), s["NEC_status"]),
             fontsize=8, color=theme.accent_two, ha="left")
    fig.text(0.010, 0.895,
             "INDEPENDENT REPRODUCTION, read from the record: production "
             "recomputed $N(\\mathbf{n}_{\\min})$ from scratch at every sample "
             "and the worst disagreement with the reported minimum is %.2e. "
             "The quadrature cross-check is %.2e. numerical error estimate %.3e."
             % (worst,
                max(abs(float(r["nec_crosscheck_residual"])) for r in rows),
                float(s["numerical_error_estimate"])),
             fontsize=7.4, color=theme.foreground, ha="left", alpha=0.9)
    _status_legend(fig, theme, A6.STATUS_LABEL, y=0.028)
    _footer(fig, theme, A6.SCOPE)
    fig.tight_layout(rect=(0, 0.085, 1, 0.875))

    plotted = [[r["x"], r["y"], r["N_min"], r["nec_direct_substitution"],
                r["nec_crosscheck_residual"], r["nec_stationarity_residual"],
                r["NEC_status"], r["nec_certification_margin"]] for r in rows]
    return A6.save_figure(
        fig, out, "minimum_nec", meta=meta,
        inputs=[os.path.join(sweep, "energy_frame_%03d.csv" % ref),
                os.path.join(sweep, "energy_velocity_sweep.csv")],
        theme=theme, observer="Eulerian orthonormal frame (null directions)",
        extra={"reference_frame": ref, "colour_limit": local,
               "colour_scale": "symmetric about zero, computed from the drawn "
                               "frame; the animations use the sweep-wide "
                               "limit instead",
               "worst_direct_substitution_disagreement": worst},
        plotted_header=["x", "y", "N_min", "nec_direct_substitution",
                        "nec_crosscheck_residual", "nec_stationarity_residual",
                        "NEC_status", "nec_certification_margin"],
        plotted_rows=plotted)


def fig_worst_null_direction_field(sweep, out, meta, theme, ref, sweep_rows):
    all_rows = records.read_csv(os.path.join(sweep,
                                             "worst_null_directions.csv"))
    rows = [r for r in all_rows if int(r["frame"]) == ref]
    s = sweep_rows[ref]

    fig, ax = plt.subplots(figsize=(7.6, 6.6), facecolor=theme.background)
    lim = max(abs(float(r["N_min"])) for r in all_rows) or 1.0
    xs = [float(r["x"]) for r in rows]
    ys = [float(r["y"]) for r in rows]
    nx = [float(r["worst_direction_x"]) for r in rows]
    ny = [float(r["worst_direction_y"]) for r in rows]
    nz = [float(r["worst_direction_z"]) for r in rows]
    nv = [float(r["N_min"]) for r in rows]

    q = ax.quiver(xs, ys, nx, ny, nv, cmap="RdBu_r", clim=(-lim, lim),
                  angles="xy", scale=26, width=0.005, pivot="mid")
    cb = fig.colorbar(q, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(r"$N_{\min}$ at that sample", fontsize=8,
                 color=theme.foreground)
    cb.ax.tick_params(labelsize=7, colors=theme.foreground)
    cb.outline.set_visible(False)

    # Out-of-plane content, marked rather than silently projected away: a 2-D
    # arrow of a 3-D unit vector is shorter exactly where |n_z| is large.
    outer = [i for i, v in enumerate(nz) if abs(v) > 0.2]
    if outer:
        ax.plot([xs[i] for i in outer], [ys[i] for i in outer],
                marker="o", markersize=3.0, linestyle="none",
                markerfacecolor="none", markeredgecolor=theme.accent_two,
                markeredgewidth=0.8)
    ax.add_patch(plt.Circle((0, 0), float(meta["radius"]), fill=False,
                            color=theme.foreground, linewidth=0.8,
                            linestyle="--", alpha=0.55))
    ax.set_aspect("equal")
    _style_axes(ax, theme)
    ax.set_title("A6  worst null direction field, $v_s$ = %.1f, z = 0"
                 % float(s["v_s"]), fontsize=11, color=theme.foreground)
    fig.text(0.012, 0.045,
             "Arrows are the SAVED minimising directions projected onto the "
             "slice. Circled samples have $|n_z| > 0.2$, so their in-plane "
             "arrow is foreshortened. Max $|n_z|$ here is %.3f."
             % max((abs(v) for v in nz), default=0.0),
             fontsize=7.2, color=theme.foreground, alpha=0.9)
    _footer(fig, theme, A6.SCOPE)
    fig.tight_layout(rect=(0, 0.10, 1, 1))

    plotted = [[r["x"], r["y"], r["N_min"], r["worst_direction_x"],
                r["worst_direction_y"], r["worst_direction_z"],
                r["NEC_status"]] for r in rows]
    return A6.save_figure(
        fig, out, "worst_null_direction_field", meta=meta,
        inputs=[os.path.join(sweep, "worst_null_directions.csv")],
        theme=theme, observer="Eulerian orthonormal frame (null directions)",
        extra={"reference_frame": ref, "arrow_scale": 26},
        plotted_header=["x", "y", "N_min", "worst_direction_x",
                        "worst_direction_y", "worst_direction_z",
                        "NEC_status"], plotted_rows=plotted)


def fig_energy_condition_status(sweep, out, meta, theme, frames, ref,
                                sweep_rows):
    rows = frames[ref]
    xs, ys = A6.axes_of(rows)
    extent = [min(xs), max(xs), min(ys), max(ys)]
    s = sweep_rows[ref]

    fig, axes = plt.subplots(1, 4, figsize=(15.6, 4.6),
                             facecolor=theme.background)
    for ax, key, lab in ((axes[0], "NEC_status", "NEC"),
                         (axes[1], "WEC_status", "WEC"),
                         (axes[2], "DEC_status", "DEC"),
                         (axes[3], "SEC_status", "SEC")):
        counts = {st: _count(rows, key, st) for st in A6.STATUS_ORDER}
        _status_panel(ax, fig, A6.status_grid(rows, key, xs, ys), extent,
                      theme, radius=float(meta["radius"]),
                      label="%s\nviolated %d   indeterminate %d   satisfied %d"
                            % (lab, counts[A6.CERTIFIED_VIOLATED],
                               counts[A6.INDETERMINATE],
                               counts[A6.CERTIFIED_SATISFIED]))

    fig.suptitle("A6  pointwise energy-condition certification at $v_s$ = %.1f"
                 % float(s["v_s"]), fontsize=12, color=theme.foreground,
                 x=0.010, ha="left", y=0.985)
    fig.text(0.010, 0.925,
             "There is no \"pass\". A value within %s error envelopes of zero "
             "is INDETERMINATE, so finite-difference noise in a near-vacuum "
             "cannot become a false violation. Envelope %.3e."
             % (meta["certification_safety_factor"],
                float(s["numerical_error_estimate"])),
             fontsize=7.8, color=theme.uncertainty_color, ha="left")
    _status_legend(fig, theme, A6.STATUS_LABEL, y=0.048)
    _footer(fig, theme, A6.SCOPE)
    fig.tight_layout(rect=(0, 0.115, 1, 0.905))

    plotted = [[r["x"], r["y"], r["NEC_status"], r["WEC_status"],
                r["DEC_status"], r["SEC_status"],
                r["nec_certification_margin"], r["wec_certification_margin"],
                r["dec_certification_margin"], r["sec_certification_margin"]]
               for r in rows]
    return A6.save_figure(
        fig, out, "energy_condition_status", meta=meta,
        inputs=[os.path.join(sweep, "energy_frame_%03d.csv" % ref)],
        theme=theme,
        observer="worst case over a deterministic sweep of the timelike "
                 "velocity ball, per condition",
        extra={"reference_frame": ref},
        plotted_header=["x", "y", "NEC_status", "WEC_status", "DEC_status",
                        "SEC_status", "nec_certification_margin",
                        "wec_certification_margin", "dec_certification_margin",
                        "sec_certification_margin"], plotted_rows=plotted)


def fig_observer_comparison(sweep, out, meta, theme, obs_rows, sweep_rows):
    sites = sorted({r["site"] for r in obs_rows})
    fig, axes = plt.subplots(2, len(sites), figsize=(13.2, 8.0),
                             facecolor=theme.background, sharex=True)
    onset = int(meta["exact_onset_frame"])
    onset_v = float(sweep_rows[onset]["v_s"])

    for c, site in enumerate(sites):
        ax = axes[0][c]
        for fam in A6.OBSERVER_FAMILIES:
            col, ls, mk = A6.OBSERVER_STYLE[fam]
            pts = [(r["v_s"], r["measured"]) for r in obs_rows
                   if r["site"] == site and r["observer_type"] == fam
                   and r["measured"] is not None]
            bad = [r["v_s"] for r in obs_rows
                   if r["site"] == site and r["observer_type"] == fam
                   and r["measured"] is None]
            if pts:
                ax.plot([p[0] for p in pts], [p[1] for p in pts], ls,
                        color=col, marker=mk, markersize=3.4, lw=1.4,
                        label=fam)
            for v in bad:
                # An inadmissible observer is drawn as an absence with a
                # reason, never as a point. There is no value to plot: the
                # record deliberately contains none.
                ax.axvline(v, color=col, lw=6, alpha=0.16)
            if bad:
                ax.plot([], [], color=col, lw=6, alpha=0.30,
                        label="%s INADMISSIBLE" % fam)
        ax.axvline(onset_v, color=theme.uncertainty_color, lw=1.0, ls=":")
        ax.axhline(0.0, color=theme.grid, lw=0.8)
        ax.set_title("%s   measured energy density" % site, fontsize=9.5,
                     color=theme.foreground)
        ax.legend(fontsize=6.6, frameon=False, labelcolor=theme.foreground,
                  ncol=2)
        _style_axes(ax, theme, "", r"$T_{\mu\nu}u^\mu u^\nu$")

        ax = axes[1][c]
        g = sorted({(r["v_s"], r["g00"]) for r in obs_rows
                    if r["site"] == site})
        ax.plot([p[0] for p in g], [p[1] for p in g], color=theme.foreground,
                lw=1.5)
        ax.axhline(0.0, color=theme.uncertainty_color, lw=1.0, ls="--")
        ax.axvline(onset_v, color=theme.uncertainty_color, lw=1.0, ls=":")
        ax.set_title(r"$g_{00}$ at the same site", fontsize=9.5,
                     color=theme.foreground)
        _style_axes(ax, theme, r"$v_s$", r"$g_{00}$")

    fig.suptitle("A6  four observer families, and where two of them cease to "
                 "exist", fontsize=12, color=theme.foreground, x=0.010,
                 ha="left", y=0.988)
    fig.text(0.010, 0.947,
             "THE COORDINATE-STATIONARY OBSERVER IS ADMISSIBLE ONLY WHERE "
             "$g_{00} < 0$. At $g_{00} = 0$ the vector $\\partial_t$ is NULL; "
             "beyond it is SPACELIKE. Neither is a worldline, the record "
             "leaves its energy density empty, and $|g_{00}|$ is never "
             "substituted. The Eulerian observer is built from the ADM normal "
             "and stays timelike straight through.",
             fontsize=7.8, color=theme.uncertainty_color, ha="left")
    fig.text(0.010, 0.905,
             "Dotted vertical: exact onset, frame %d, $v_s$ = %.1f. Shaded "
             "bands: a family that is inadmissible at that $v_s$ and site."
             % (onset, onset_v),
             fontsize=7.4, color=theme.foreground, ha="left", alpha=0.9)
    _footer(fig, theme, A6.CROSSING_NOTE + "  " + A6.SCOPE)
    fig.tight_layout(rect=(0, 0.04, 1, 0.895))

    plotted = [[r["frame"], r["v_s"], r["site"], r["observer_type"],
                int(r["valid"]),
                "" if r["measured"] is None else r["measured"],
                r["g00"], r["norm_residual"], r["reason"]] for r in obs_rows]
    return A6.save_figure(
        fig, out, "observer_comparison", meta=meta,
        inputs=[os.path.join(sweep, "observer_comparison.csv")],
        theme=theme,
        observer="all four families, kept separate and separately labelled",
        extra={"observer_families": list(A6.OBSERVER_FAMILIES),
               "inadmissible_rows": sum(1 for r in obs_rows if not r["valid"])},
        plotted_header=["frame", "v_s", "site", "observer_type",
                        "observer_valid", "rho_measured", "g00",
                        "observer_normalization_residual",
                        "observer_rejection_reason"],
        plotted_rows=plotted)


def fig_integrated_energy(sweep, out, meta, theme, integ, sweep_rows):
    v = [float(r["v_s"]) for r in integ]
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8),
                             facecolor=theme.background)

    ax = axes[0]
    for key, lab, col in (("E_plus", r"$E_+$", theme.accent),
                          ("E_minus", r"$E_-$", theme.uncertainty_color),
                          ("E_abs", r"$E_{\rm abs}$", theme.accent_two)):
        ax.semilogy(v, [abs(float(r[key])) for r in integ], lw=1.6, color=col,
                    marker="o", markersize=3.2, label=lab)
    ax.semilogy(v, [abs(float(r["E_net"])) for r in integ], lw=1.4,
                color="#7ed957", ls="--", marker="s", markersize=3.0,
                label=r"$|E_{\rm net}|$")
    ax.legend(fontsize=7.5, frameon=False, labelcolor=theme.foreground)
    ax.set_title("integrated energy objectives", fontsize=9.5,
                 color=theme.foreground)
    _style_axes(ax, theme, r"$v_s$", "magnitude")

    ax = axes[1]
    ax.semilogy(v, [float(r["E_minus"]) for r in integ], lw=1.7,
                color=theme.uncertainty_color, marker="o", markersize=3.2,
                label=r"EXACT $E_- = \int\sqrt{\gamma}\max(0,-\rho)$")
    ax.semilogy(v, [float(r["smooth_surrogate_value"]) for r in integ], lw=1.5,
                color=theme.accent, ls=":", marker="D", markersize=3.0,
                label="softplus surrogate (width %s)" % integ[0]["softplus_width"])
    ax.legend(fontsize=7.0, frameon=False, labelcolor=theme.foreground)
    ax.set_title("exact objective vs smooth surrogate -- NEVER interchanged",
                 fontsize=9.0, color=theme.foreground)
    _style_axes(ax, theme, r"$v_s$", "value")

    ax = axes[2]
    ax.plot(v, [float(r["e_minus_interior"]) for r in integ], lw=1.6,
            color=theme.accent, marker="o", markersize=3.2,
            label=r"$E_-$ interior")
    ax.plot(v, [float(r["e_minus_shell"]) for r in integ], lw=1.6,
            color=theme.accent_two, marker="s", markersize=3.2,
            label=r"$E_-$ boundary shell")
    ax.plot(v, [float(r["tail_estimate"]) for r in integ], lw=1.4,
            color="#7ed957", ls="--", marker="^", markersize=3.2,
            label="estimated exterior tail")
    ax2 = ax.twinx()
    ax2.plot(v, [float(r["boundary_fraction"]) for r in integ], lw=1.2,
             color=theme.uncertainty_color, ls=":")
    ax2.axhline(0.35, color=theme.uncertainty_color, lw=0.9, ls="--",
                alpha=0.6)
    ax2.set_ylabel("boundary fraction (dotted)", fontsize=7.5,
                   color=theme.uncertainty_color)
    ax2.tick_params(labelsize=7, colors=theme.uncertainty_color)
    ax.legend(fontsize=7.0, frameon=False, labelcolor=theme.foreground)
    ax.set_title("interior, boundary shell and exterior tail", fontsize=9.5,
                 color=theme.foreground)
    _style_axes(ax, theme, r"$v_s$", "contribution")

    closure = max(max(float(r["net_closure_residual"]),
                      float(r["abs_closure_residual"])) for r in integ)
    fig.suptitle("A6  integrated energy objectives against $v_s$",
                 fontsize=12, color=theme.foreground, x=0.010, ha="left",
                 y=0.985)
    fig.text(0.010, 0.925,
             "Closure verified on the saved values before plotting: "
             "$|E_{\\rm net}-(E_+-E_-)|$ and $|E_{\\rm abs}-(E_++E_-)|$ are at "
             "most %.2e. $E_+$ and $E_-$ are reported separately so "
             "$E_{\\rm net}$ cannot hide the cancellation. %d of %d velocities "
             "are marked optimization-ready."
             % (closure, sum(1 for r in integ if r["optimization_ready"] == "1"),
                len(integ)),
             fontsize=7.6, color=theme.foreground, ha="left", alpha=0.92)
    _footer(fig, theme, "An objective dominated by the numerical boundary is "
                        "measuring the domain truncation, not the bubble. " +
            A6.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.905))

    hdr = ["v_s", "E_plus", "E_minus", "E_net", "E_abs",
           "smooth_surrogate_value", "softplus_width", "net_closure_residual",
           "abs_closure_residual", "e_minus_interior", "e_minus_shell",
           "tail_estimate", "boundary_fraction", "optimization_ready"]
    return A6.save_figure(
        fig, out, "integrated_energy_vs_velocity", meta=meta,
        inputs=[os.path.join(sweep, "integrated_objectives.csv")],
        theme=theme, observer="Eulerian (ADM normal n^mu)",
        extra={"worst_closure_residual": closure,
               "exact_and_smooth_objectives_distinct": True},
        plotted_header=hdr,
        plotted_rows=[[r[k] for k in hdr] for r in integ])


def fig_integrated_nec(sweep, out, meta, theme, integ, necsweep):
    v = [float(r["v_s"]) for r in integ]
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.8),
                             facecolor=theme.background)

    ax = axes[0]
    ax.plot(v, [float(r["J_NEC"]) for r in integ], lw=1.7,
            color=theme.uncertainty_color, marker="o", markersize=3.4,
            label=r"EXACT $J_{\rm NEC}=\int\sqrt{\gamma}\max(0,-N_{\min})$")
    ax.plot(v, [float(r["J_NEC_smooth"]) for r in integ], lw=1.5,
            color=theme.accent, ls=":", marker="D", markersize=3.0,
            label="softplus surrogate (width %s)" % integ[0]["softplus_width"])
    ax.legend(fontsize=7.0, frameon=False, labelcolor=theme.foreground)
    ax.set_title("integrated NEC violation -- exact and surrogate",
                 fontsize=9.5, color=theme.foreground)
    _style_axes(ax, theme, r"$v_s$", "value")

    ax = axes[1]
    ax.plot(v, [float(r["j_nec_interior"]) for r in integ], lw=1.6,
            color=theme.accent, marker="o", markersize=3.2, label="interior")
    ax.plot(v, [float(r["j_nec_shell"]) for r in integ], lw=1.6,
            color=theme.accent_two, marker="s", markersize=3.2,
            label="boundary shell")
    ax2 = ax.twinx()
    ax2.plot(v, [abs(float(r["N_min"])) for r in necsweep], lw=1.3,
             color="#7ed957", ls="--")
    ax2.set_ylabel(r"$|N_{\min}|$ on the slice (dashed)", fontsize=7.5,
                   color="#7ed957")
    ax2.tick_params(labelsize=7, colors="#7ed957")
    ax.legend(fontsize=7.2, frameon=False, labelcolor=theme.foreground)
    ax.set_title("where the violation integral comes from", fontsize=9.5,
                 color=theme.foreground)
    _style_axes(ax, theme, r"$v_s$", "contribution")

    fig.suptitle("A6  integrated NEC violation against $v_s$", fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    fig.text(0.010, 0.925,
             "The surrogate is provided for a future gradient-based optimizer "
             "and is NEVER used for pass/fail: softplus systematically "
             "overestimates $\\max(0,x)$ near zero, which is exactly where "
             "certification happens. Both curves are plotted so the gap is "
             "visible rather than assumed small.",
             fontsize=7.5, color=theme.foreground, ha="left", alpha=0.92)
    _footer(fig, theme, A6.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.905))

    hdr = ["v_s", "J_NEC", "J_NEC_smooth", "softplus_width", "j_nec_interior",
           "j_nec_shell"]
    return A6.save_figure(
        fig, out, "integrated_nec_vs_velocity", meta=meta,
        inputs=[os.path.join(sweep, "integrated_objectives.csv"),
                os.path.join(sweep, "nec_velocity_sweep.csv")],
        theme=theme, observer="Eulerian orthonormal frame (null directions)",
        extra={"exact_and_smooth_objectives_distinct": True},
        plotted_header=hdr, plotted_rows=[[r[k] for k in hdr] for r in integ])


def fig_resolution_certification(sweep, out, meta, theme):
    rows = records.read_csv(os.path.join(sweep,
                                         "resolution_certification.csv"))
    h = [float(r["h"]) for r in rows]
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.6),
                             facecolor=theme.background)

    ax = axes[0]
    ax.plot(h, [float(r["E_minus"]) for r in rows], lw=1.7,
            color=theme.uncertainty_color, marker="o", markersize=4,
            label=r"$E_-$")
    ax.plot(h, [float(r["J_NEC"]) / 40.0 for r in rows], lw=1.4,
            color=theme.accent, ls="--", marker="s", markersize=3.4,
            label=r"$J_{\rm NEC}/40$")
    ax.legend(fontsize=7.4, frameon=False, labelcolor=theme.foreground)
    ax.invert_xaxis()
    ax.set_title("integrated objectives against grid spacing", fontsize=9.5,
                 color=theme.foreground)
    _style_axes(ax, theme, "h (finer to the right)", "value")

    ax = axes[1]
    ax.plot(h, [abs(float(r["min_rho"])) for r in rows], lw=1.7,
            color=theme.accent_two, marker="o", markersize=4,
            label=r"$|\rho_{\min}|$")
    ax.plot(h, [float(r["numerical_error_estimate"]) for r in rows], lw=1.4,
            color=theme.grid, ls="--", label="numerical error estimate")
    ax.plot(h, [float(r["numerical_error_estimate"]) *
                float(meta["certification_safety_factor"]) for r in rows],
            lw=1.4, color=theme.uncertainty_color, ls=":",
            label="certification threshold")
    ax.legend(fontsize=7.0, frameon=False, labelcolor=theme.foreground)
    ax.invert_xaxis()
    ax.set_title(r"$|\rho_{\min}|$ against what would certify it",
                 fontsize=9.5, color=theme.foreground)
    _style_axes(ax, theme, "h (finer to the right)", "magnitude")

    ax = axes[2]
    for k, (key, lab) in enumerate((("rho_sign_status", r"$\rho$ sign"),
                                    ("NEC_status", "NEC"))):
        for i, r in enumerate(rows):
            ax.plot([i], [k], marker="s", markersize=15, linestyle="none",
                    markerfacecolor=A6.STATUS_COLOR[r[key]],
                    markeredgecolor=theme.grid)
    ax.set_yticks([0, 1])
    ax.set_yticklabels([r"$\rho$ sign", "NEC"], fontsize=8,
                       color=theme.foreground)
    ax.set_xticks(list(range(len(rows))))
    ax.set_xticklabels(["n=%s\nwall %s" % (r["n"], r["wall_cells"])
                        for r in rows], fontsize=7, color=theme.foreground)
    ax.set_ylim(-0.6, 1.6)
    ax.set_title("what each grid CERTIFIES", fontsize=9.5,
                 color=theme.foreground)
    _style_axes(ax, theme, "", "")

    fig.suptitle("A6  resolution certification at $v_s$ = %s"
                 % meta["reference_velocity"], fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    fig.text(0.010, 0.925,
             "The same physical claim can be certified on one grid and "
             "indeterminate on a coarser one, which is why the status is drawn "
             "per resolution rather than quoted once. Cells across the bubble "
             "wall are measured from the saved ADM shift, not assumed.",
             fontsize=7.5, color=theme.foreground, ha="left", alpha=0.92)
    _status_legend(fig, theme, A6.STATUS_LABEL, y=0.045)
    _footer(fig, theme, A6.SCOPE, y=0.004)
    fig.tight_layout(rect=(0, 0.115, 1, 0.905))

    hdr = ["n", "h", "wall_cells", "E_minus", "J_NEC", "min_rho", "N_min",
           "numerical_error_estimate", "certification_margin",
           "rho_sign_status", "NEC_status"]
    return A6.save_figure(
        fig, out, "resolution_certification", meta=meta,
        inputs=[os.path.join(sweep, "resolution_certification.csv")],
        theme=theme, observer="Eulerian (ADM normal n^mu)",
        plotted_header=hdr, plotted_rows=[[r[k] for k in hdr] for r in rows])


def fig_domain_convergence(sweep, out, meta, theme):
    rows = records.read_csv(os.path.join(sweep, "domain_convergence.csv"))
    hw = [float(r["half_width"]) for r in rows]
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.6),
                             facecolor=theme.background)

    ax = axes[0]
    ax.plot(hw, [float(r["E_minus"]) for r in rows], lw=1.7,
            color=theme.uncertainty_color, marker="o", markersize=4,
            label=r"$E_-$")
    ax.plot(hw, [float(r["E_abs"]) for r in rows], lw=1.5,
            color=theme.accent_two, marker="s", markersize=3.6,
            label=r"$E_{\rm abs}$")
    ax.plot(hw, [float(r["E_plus"]) for r in rows], lw=1.4,
            color=theme.accent, ls="--", marker="^", markersize=3.6,
            label=r"$E_+$")
    ax.legend(fontsize=7.4, frameon=False, labelcolor=theme.foreground)
    ax.set_title("integrated energy against domain half-width", fontsize=9.5,
                 color=theme.foreground)
    _style_axes(ax, theme, "half-width", "value")

    ax = axes[1]
    ax.plot(hw, [float(r["boundary_fraction"]) for r in rows], lw=1.7,
            color=theme.accent, marker="o", markersize=4,
            label="boundary fraction")
    ax.axhline(0.35, color=theme.uncertainty_color, lw=1.0, ls="--",
               label="boundary-dominated threshold")
    ax.legend(fontsize=7.2, frameon=False, labelcolor=theme.foreground)
    ax.set_title("share of the integral in the outer shell", fontsize=9.5,
                 color=theme.foreground)
    _style_axes(ax, theme, "half-width", "fraction")

    ax = axes[2]
    ax.plot(hw, [float(r["tail_estimate"]) for r in rows], lw=1.7,
            color="#7ed957", marker="o", markersize=4,
            label="estimated exterior tail")
    ax.plot(hw, [float(r["E_minus"]) * 0.25 for r in rows], lw=1.2,
            color=theme.grid, ls="--", label=r"25% of $E_-$")
    ax.legend(fontsize=7.2, frameon=False, labelcolor=theme.foreground)
    ax.set_title("what enlarging the domain adds", fontsize=9.5,
                 color=theme.foreground)
    _style_axes(ax, theme, "half-width", "value")

    fig.suptitle("A6  domain convergence at $v_s$ = %s"
                 % meta["reference_velocity"], fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    fig.text(0.010, 0.925,
             "The point count is held fixed while the extent grows, so h grows "
             "with it -- the same convention the A6 gate's domain study uses, "
             "stated rather than hidden. If the tail grew with the domain the "
             "integral would not exist and no optimizer could use it.",
             fontsize=7.5, color=theme.foreground, ha="left", alpha=0.92)
    _footer(fig, theme, A6.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.905))

    hdr = ["half_width", "h", "E_plus", "E_minus", "E_abs", "J_NEC",
           "boundary_fraction", "tail_estimate", "min_rho", "N_min",
           "wall_cells"]
    return A6.save_figure(
        fig, out, "domain_convergence", meta=meta,
        inputs=[os.path.join(sweep, "domain_convergence.csv")],
        theme=theme, observer="Eulerian (ADM normal n^mu)",
        plotted_header=hdr, plotted_rows=[[r[k] for k in hdr] for r in rows])


def fig_null_direction_sphere(sweep, out, meta, theme):
    samples = records.read_csv(os.path.join(sweep,
                                            "null_direction_sphere.csv"))
    points = records.read_csv(os.path.join(
        sweep, "null_direction_sphere_points.csv"))

    labels = [p["point_label"] for p in points]
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 8.2),
                             facecolor=theme.background,
                             subplot_kw={"projection": "mollweide"})
    flat_axes = [a for row in axes for a in row]

    for ax, p in zip(flat_axes, points):
        lab = p["point_label"]
        rows = [r for r in samples if r["point_label"] == lab]
        # Mollweide wants (longitude, latitude). phi and theta are both saved
        # by the exporter; latitude is pi/2 - theta, one subtraction.
        lon = [float(r["phi"]) for r in rows]
        lat = [math.pi / 2.0 - float(r["theta"]) for r in rows]
        val = [float(r["N_of_n"]) for r in rows]
        m = max(abs(v) for v in val) or 1.0
        sc = ax.scatter(lon, lat, c=val, cmap="RdBu_r", vmin=-m, vmax=m, s=9,
                        linewidths=0)
        cb = fig.colorbar(sc, ax=ax, fraction=0.030, pad=0.04)
        cb.ax.tick_params(labelsize=6.5, colors=theme.foreground)
        cb.outline.set_visible(False)

        ax.plot([float(p["production_phi"])],
                [math.pi / 2.0 - float(p["production_theta"])], marker="*",
                markersize=15, color=theme.accent_two, linestyle="none",
                label="production minimiser")
        ax.plot([float(p["quadrature_phi"])],
                [math.pi / 2.0 - float(p["quadrature_theta"])], marker="o",
                markersize=8, markerfacecolor="none",
                markeredgecolor=theme.uncertainty_color, markeredgewidth=1.4,
                linestyle="none", label="quadrature minimum")
        # The coordinate axes, so "off-axis" is visible rather than asserted.
        ax.plot([0.0, math.pi, math.pi / 2.0, -math.pi / 2.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, math.pi / 2.0, -math.pi / 2.0],
                marker="+", markersize=6, linestyle="none",
                color=theme.grid, label="coordinate axes")
        ax.set_title("%s   $N_{\\min}$ = %.4e   flat fraction %.1f%%"
                     % (lab, float(p["N_min"]),
                        100.0 * float(p["near_flat_fraction"])),
                     fontsize=8.6, color=theme.foreground, pad=12)
        ax.grid(True, color=theme.grid, alpha=0.35, linewidth=0.5)
        ax.tick_params(labelsize=5.5, colors=theme.foreground)
        ax.set_facecolor(theme.background)

    handles, lgl = flat_axes[0].get_legend_handles_labels()
    fig.legend(handles, lgl, loc="lower center", ncol=3, fontsize=7.4,
               frameon=False, bbox_to_anchor=(0.5, 0.028),
               labelcolor=theme.foreground)

    gap = max(abs(float(p["N_min"]) - float(p["quadrature_minimum"]))
              for p in points)
    fig.suptitle("A6  $N(\\mathbf{n})$ over the WHOLE unit sphere at four "
                 "representative points, $v_s$ = %s" % meta["reference_velocity"],
                 fontsize=12, color=theme.foreground, x=0.010, ha="left",
                 y=0.988)
    fig.text(0.010, 0.948,
             "%s deterministic Fibonacci directions per point, evaluated by "
             "production. The + markers are the six coordinate axes: testing "
             "only those samples six points of a continuum, and wherever the "
             "star does not land on one, an axis-only search would have "
             "reported a less negative worst case. Worst "
             "production-vs-quadrature gap %.2e."
             % (meta["sphere_sample_count"], gap),
             fontsize=7.6, color=theme.foreground, ha="left", alpha=0.92)
    _footer(fig, theme, "\"flat fraction\" is the share of sampled directions "
                        "within 1% of the minimum: a large value means the "
                        "minimiser is weakly determined. " + A6.SCOPE,
            y=0.004)
    fig.tight_layout(rect=(0, 0.075, 1, 0.925))

    hdr = ["point_label", "x", "y", "z", "rho", "N_min",
           "production_direction_x", "production_direction_y",
           "production_direction_z", "quadrature_minimum",
           "near_flat_fraction", "nec_hard_case", "sample_count"]
    return A6.save_figure(
        fig, out, "null_direction_sphere", meta=meta,
        inputs=[os.path.join(sweep, "null_direction_sphere.csv"),
                os.path.join(sweep, "null_direction_sphere_points.csv")],
        theme=theme, observer="Eulerian orthonormal frame (null directions)",
        extra={"probe_points": labels,
               "worst_production_quadrature_gap": gap},
        plotted_header=hdr, plotted_rows=[[p[k] for k in hdr] for p in points])


# --- Diagnostics -------------------------------------------------------------

def _count(rows, key, value):
    return sum(1 for r in rows if r[key] == value)


def run_diagnostics(rep, sweep, meta, frames, sweep_rows, necsweep, integ,
                    obs_rows, before, after):
    """The A6 anomaly checks this renderer is in a position to make."""
    flat = [r for rows in frames for r in rows]

    # 1. Nonfinite observer-frame values.
    vals = []
    for key in ("rho", "q_x", "q_y", "q_z", "N_min", "principal_pressure_1",
                "principal_pressure_3", "g00", "sqrt_gamma"):
        vals.extend(float(r[key]) for r in flat)
    rep.check_finite("observer frame fields", vals)

    # 2. Eulerian normalisation.
    resid = [float(r["observer_normalization_residual"]) for r in flat]
    resid += [r["norm_residual"] for r in obs_rows
              if r["observer_type"] == "eulerian"]
    rep.record("eulerian_normalization_failure",
               fired=max(resid) > 1e-10, measured=max(resid), threshold=1e-10,
               detail="worst |g(n,n)+1| over every drawn sample")

    # 3. Tetrad orthonormality.
    ortho = [float(r["tetrad_orthonormality_residual"]) for r in flat]
    rep.record("tetrad_orthonormality_failure",
               fired=max(ortho) > 1e-12, measured=max(ortho), threshold=1e-12,
               detail="worst max|g(e_a,e_b) - eta_ab| over every drawn sample")

    # 4. Handedness. A frame determinant that changed sign between neighbouring
    # samples would appear as a spurious sign flip in q_i and in the worst null
    # direction, which no physics produced.
    hand = [float(r["tetrad_handedness"]) for r in flat]
    rep.record("tetrad_handedness_change",
               fired=not (min(hand) > 0.0) == (max(hand) > 0.0),
               measured=[min(hand), max(hand)], threshold="one sign throughout",
               detail="frame determinant range over every drawn sample")

    # 5. Arbitrary sign flip. Fixed-order Gram-Schmidt with a positive root
    # makes every diagonal frame component positive BY CONSTRUCTION.
    diag = [float(r["tetrad_min_diagonal"]) for r in flat]
    rep.record("tetrad_sign_flip", fired=min(diag) <= 0.0, measured=min(diag),
               threshold="> 0",
               detail="min over samples of min_a e_a^a; the positive-diagonal "
                      "convention makes a non-positive value a sign flip")

    # 6. An inadmissible coordinate-stationary observer presented as valid.
    bad = [(r["frame"], r["site"], r["g00"]) for r in obs_rows
           if r["observer_type"] == "coordinate_stationary" and r["valid"]
           and r["g00"] >= 0.0]
    rep.record("invalid_stationary_observer_rendered_valid",
               fired=bool(bad), measured=len(bad), threshold=0,
               detail="coordinate-stationary rows marked valid where "
                      "g_00 >= 0: %r" % bad[:3])

    # 8. Off-axis negative energy. rho vanishes on the axis of travel because
    # it depends on the TRANSVERSE derivatives only, so a slice that showed
    # negative energy only on the axis would mean the figure had cropped the
    # torus away.
    misses = []
    for idx, rows in enumerate(frames):
        off = sum(1 for r in rows
                  if float(r["rho"]) < 0.0 and float(r["y"]) != 0.0)
        drawn_min = min(float(r["rho"]) for r in rows)
        recorded = float(sweep_rows[idx]["slice_rho_min"])
        if off == 0 or abs(drawn_min - recorded) > 1e-12:
            misses.append((idx, off, drawn_min, recorded))
    rep.record("missed_off_axis_negative_energy", fired=bool(misses),
               measured=len(misses), threshold=0,
               detail="frames whose drawn slice has no off-axis negative "
                      "sample or whose minimum disagrees with the sweep "
                      "record: %r" % misses[:3])

    # 9. The reported direction must reproduce N_min. Compared between two
    # SAVED columns: production recomputed N at its own reported direction.
    worst = max(abs(float(r["N_min"]) - float(r["nec_direct_substitution"]))
                for r in flat)
    rep.record("nec_direction_does_not_reproduce_n_min",
               fired=worst > 1e-12, measured=worst, threshold=1e-12,
               detail="worst |N_min - N(n_min) recomputed by production|")

    # 11. Status against margin, per condition. margin = |value| - k*envelope,
    # so a certified status requires a positive margin and an indeterminate one
    # requires a non-positive margin. Two production outputs, cross-checked.
    incons = []
    for r in flat:
        for st, mg in (("NEC_status", "nec_certification_margin"),
                       ("WEC_status", "wec_certification_margin"),
                       ("DEC_status", "dec_certification_margin"),
                       ("SEC_status", "sec_certification_margin"),
                       ("rho_sign_status", "rho_certification_margin")):
            definite = r[st] in (A6.CERTIFIED_SATISFIED, A6.CERTIFIED_VIOLATED)
            if definite != (float(r[mg]) > 0.0):
                incons.append((r["frame"], r["x"], r["y"], st, r[st], r[mg]))
    rep.record("status_inconsistent_with_certification_margin",
               fired=bool(incons), measured=len(incons), threshold=0,
               detail="samples whose status and margin disagree: %r"
                      % incons[:3])

    # 12. Near-zero uncertified values drawn as definite. The figures colour
    # every sample by its SAVED status, so this asks whether any sample placed
    # in a definite category has a non-positive margin.
    definite_bad = sum(
        1 for r in flat
        if r["rho_sign_status"] in (A6.CERTIFIED_SATISFIED,
                                    A6.CERTIFIED_VIOLATED)
        and float(r["rho_certification_margin"]) <= 0.0)
    rep.record("uncertified_value_rendered_definite",
               fired=definite_bad > 0, measured=definite_bad, threshold=0,
               detail="samples drawn with a definite sign whose certification "
                      "margin is not positive")

    # 13. Closure of the integrated objectives, on the saved values.
    net = max(float(r["net_closure_residual"]) for r in integ)
    ab = max(float(r["abs_closure_residual"]) for r in integ)
    rep.record("integrated_objective_closure_failure",
               fired=max(net, ab) > 1e-9, measured=[net, ab], threshold=1e-9,
               detail="worst |E_net-(E_+-E_-)| and |E_abs-(E_++E_-)|")

    # 14. A boundary-dominated integral is measuring the domain truncation.
    bf = max(float(r["boundary_fraction"]) for r in integ)
    rep.record("boundary_dominated_integral", fired=bf > 0.35, measured=bf,
               threshold=0.35,
               detail="worst share of E_abs sitting in the outer shell")

    # 15. Cells across the bubble wall.
    wc = min(int(r["wall_cells"]) for r in integ)
    rep.record("insufficient_cells_across_wall", fired=wc < 4, measured=wc,
               threshold=4,
               detail="fewest samples resolving the wall, measured from the "
                      "saved ADM shift")

    # 18. Source rows.
    defects = []
    if len(frames) != len(sweep_rows) or len(frames) != len(necsweep) \
            or len(frames) != len(integ):
        defects.append(("row counts", len(frames), len(sweep_rows),
                        len(necsweep), len(integ)))
    n = int(meta["grid_resolution"])
    for idx, rows in enumerate(frames):
        if len(rows) != n * n:
            defects.append(("frame %d has %d rows" % (idx, len(rows))))
        if int(rows[0]["frame"]) != idx:
            defects.append(("frame %d mislabelled %s" % (idx,
                                                         rows[0]["frame"])))
        if float(rows[0]["v_s"]) != float(sweep_rows[idx]["v_s"]):
            defects.append(("frame %d v_s disagrees" % idx))
    rep.record("missing_or_mismatched_source_rows", fired=bool(defects),
               measured=len(defects), threshold=0, detail=repr(defects[:3]))

    # 19. Source stability across the render.
    changed = [k for k in before if before[k] != after.get(k)]
    changed += [k for k in after if k not in before]
    rep.record("source_checksum_mismatch", fired=bool(changed),
               measured=len(changed), threshold=0,
               detail="records whose SHA-256 changed while the figures were "
                      "being drawn: %r" % changed[:3])

    # 20. Recomputed physics in this package.
    offenders = A6.scan_python_physics(A6.package_python_files())
    rep.record("python_physics_recomputation", fired=bool(offenders),
               measured=len(offenders), threshold=0,
               detail="AST scan of the A6 renderer package: %r"
                      % offenders[:3])


# --- Driver ------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--theme", default="scientific")
    args = ap.parse_args(argv)

    sweep = records.resolve(args.sweep)
    out = records.ensure_dir(args.out)
    theme = themes.get(args.theme)
    themes.apply(theme)

    meta = A6.load_manifest(sweep)
    before = A6.source_checksums(sweep)

    frames = A6.load_frames(sweep)
    sweep_rows = records.read_csv(os.path.join(sweep,
                                               "energy_velocity_sweep.csv"))
    necsweep = records.read_csv(os.path.join(sweep, "nec_velocity_sweep.csv"))
    integ = records.read_csv(os.path.join(sweep, "integrated_objectives.csv"))
    obs_rows = A6.observer_rows(sweep)
    ref = int(meta["reference_frame"])

    # Colour limits computed ONCE from the whole sweep and then frozen, so the
    # static figures and the animations share a scale and a field that grows by
    # an order of magnitude across the sweep looks like it grew.
    limits = {
        "rho": A6.global_symmetric_limit(frames, "rho"),
        "q": max(A6.global_symmetric_limit(frames, "q_x"),
                 A6.global_symmetric_limit(frames, "q_y"),
                 A6.global_symmetric_limit(frames, "q_z")),
        "N_min": A6.global_symmetric_limit(frames, "N_min"),
    }
    limits["q_arrow"] = 0.35 / limits["q"] if limits["q"] > 0 else 1.0
    with open(os.path.join(out, "global_scales.json"), "w",
              encoding="utf-8") as fh:
        json.dump({"limits": limits,
                   "computed_from": "every frame of the velocity sweep",
                   "fixed_across_frames": True}, fh, indent=2, sort_keys=True)
    with open(os.path.join(out, "source_checksums.json"), "w",
              encoding="utf-8") as fh:
        json.dump(before, fh, indent=2, sort_keys=True)

    made = []
    made.append(fig_energy_density(sweep, out, meta, theme, frames, ref,
                                   sweep_rows, limits))
    made.append(fig_momentum_density(sweep, out, meta, theme, frames, ref,
                                     sweep_rows, limits))
    made.append(fig_principal_stresses(sweep, out, meta, theme, frames, ref,
                                       sweep_rows, limits))
    made.append(fig_minimum_nec(sweep, out, meta, theme, frames, ref,
                                sweep_rows, limits))
    made.append(fig_worst_null_direction_field(sweep, out, meta, theme, ref,
                                               sweep_rows))
    made.append(fig_energy_condition_status(sweep, out, meta, theme, frames,
                                            ref, sweep_rows))
    made.append(fig_observer_comparison(sweep, out, meta, theme, obs_rows,
                                        sweep_rows))
    made.append(fig_integrated_energy(sweep, out, meta, theme, integ,
                                      sweep_rows))
    made.append(fig_integrated_nec(sweep, out, meta, theme, integ, necsweep))
    made.append(fig_resolution_certification(sweep, out, meta, theme))
    made.append(fig_domain_convergence(sweep, out, meta, theme))
    made.append(fig_null_direction_sphere(sweep, out, meta, theme))

    after = A6.source_checksums(sweep)
    rep = anomaly.AnomalyReport("A", "A6", "a6_static")
    run_diagnostics(rep, sweep, meta, frames, sweep_rows, necsweep, integ,
                    obs_rows, before, after)
    path = rep.write(out)

    for p in made:
        print("wrote %s" % os.path.basename(p))
    print("%d figures, %d fired anomalies -> %s"
          % (len(made), len(rep.fired), os.path.basename(path)))
    # A fired anomaly does NOT suppress the figure: the evidence survives and
    # the visual gate is what fails.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
