#!/usr/bin/env python3
"""The eight A6 observer-energy animations, drawn from saved records only.

    python animate_a6_energy.py --sweep <a6 sweep dir> --out <gallery dir>

Every frame is one velocity of the established sequence v_s = k/10, k = 2..16.
Frame 8 is v_s = 1 EXACTLY, where g_00 reaches zero at the bubble centre, and
frame 9 is the first whose saved grid brackets zero and can therefore carry a
contour. The animation is what makes that crossing legible, and it is also where
the crossing is most likely to be misread, so every frame states:

  * which observer measured the quantity on screen, and whether that observer
    still exists at this v_s;
  * that g_00 = 0 is a statement about the coordinate vector d/dt and NOT a
    metric degeneracy, signature failure, curvature singularity or horizon;
  * whether the value shown is certified or merely measured.

Colour limits are computed once from every frame and frozen. The cinematic theme
changes typography, easing and glow. It does not change the data, the positions,
the timing, the contours, the certification, the normalisation or the colour
limits, and the anomaly report checks the last of those directly.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402
from matplotlib.colors import ListedColormap, BoundaryNorm       # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

import a6_common as A6                                           # noqa: E402
from common import anomaly, records, render                      # noqa: E402
from themes import themes                                        # noqa: E402


FPS = 4

# Every frame of every animation carries this block. The list is fixed here so
# an animation cannot quietly drop a field.
def overlay_lines(meta, sweep_rows, necsweep, integ, i, observer,
                  observer_state):
    s, n, g = sweep_rows[i], necsweep[i], integ[i]
    return [
        "v_s              %.1f      frame %d/%d" % (float(s["v_s"]), i,
                                                    len(sweep_rows) - 1),
        "observer         %s" % observer,
        "validity         %s" % observer_state,
        "rho minimum      %+.4e  at (%.2f, %.2f, %.2f)"
        % (float(s["rho_min"]), float(s["rho_min_x"]), float(s["rho_min_y"]),
           float(s["rho_min_z"])),
        "N_min            %+.4e" % float(s["N_min"]),
        "E_minus          %.5e   E_abs %.5e"
        % (float(s["E_minus"]), float(s["E_abs"])),
        "J_NEC            %.5e" % float(s["J_NEC"]),
        "NEC %s" % s["NEC_status"],
        "WEC %s   DEC %s" % (s["WEC_status"], s["DEC_status"]),
        "SEC %s   rho sign %s" % (s["SEC_status"], s["rho_sign_status"]),
        "grid             %s^3 over [-%s, %s]^3, %s cells across the wall"
        % (s["grid_resolution"], s["domain_extent"], s["domain_extent"],
           s["wall_cells"]),
        "error estimate   %.3e   certification margin %+.3e"
        % (float(s["numerical_error_estimate"]),
           float(s["certification_margin"])),
        "Bianchi relative %.4f  (diagnostic only)"
        % float(s["bianchi_relative_residual"]),
        "NEC crosscheck   %.2e   stationarity %.2e"
        % (float(n["nec_crosscheck_residual"]),
           float(n["nec_stationarity_residual"])),
        "boundary share   %.4f   tail %.3e"
        % (float(g["boundary_fraction"]), float(g["tail_estimate"])),
    ]


def crossing_status(meta, i):
    onset = int(meta["exact_onset_frame"])
    contour = int(meta["first_spatial_contour_frame"])
    if i < onset:
        return ("below onset: $g_{00} < 0$ everywhere, every observer family "
                "admissible at both sites")
    if i == onset:
        return ("EXACT ONSET: $v_s = 1$, $g_{00}(0) = 0$ exactly. The "
                "coordinate-stationary observer is NULL here and ceases to "
                "exist; the Eulerian observer does not")
    if i == contour:
        return ("first frame whose saved grid brackets zero, so a contour is "
                "drawable; the coordinate-stationary observer is SPACELIKE "
                "inside")
    return ("past onset: $g_{00} > 0$ inside the bubble. The "
            "coordinate-stationary observer is spacelike there and is not "
            "rendered")


def stamp(fig, theme, lines, x=0.010, y=0.40):
    fig.text(x, y, "\n".join(lines), ha="left", va="center", fontsize=6.1,
             family="monospace", color=theme.foreground, linespacing=1.5,
             bbox=dict(facecolor=theme.background, edgecolor=theme.grid,
                       alpha=0.85 if theme.trail else 0.72,
                       boxstyle="round,pad=0.5", linewidth=0.6))


def banner(fig, theme, title, status, warning):
    # The cinematic theme differs here only in weight and glow; the STRINGS are
    # identical, because a caveat that fades with the styling is not a caveat.
    fig.suptitle(title, fontsize=13 if theme.trail else 12,
                 color=theme.foreground, x=0.010, ha="left", y=0.978,
                 fontweight="bold" if theme.trail else "normal")
    fig.text(0.010, 0.925, status, fontsize=8.4, color=theme.accent_two,
             ha="left")
    fig.text(0.010, 0.878, warning, fontsize=7.8,
             color=theme.uncertainty_color, ha="left")


def footer(fig, theme):
    fig.text(0.010, 0.008,
             A6.CROSSING_NOTE + "  " + A6.SCOPE, fontsize=6.0,
             color=theme.foreground, alpha=0.72, ha="left", va="bottom")


def png_resolution(path):
    """Read the real pixel size out of the PNG header.

    Measured rather than asserted: a resolution field copied from a figsize
    guess is a metadata field that quietly stops being true the first time a
    layout changes.
    """
    with open(path, "rb") as fh:
        head = fh.read(24)
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        return [0, 0]
    return [int.from_bytes(head[16:20], "big"),
            int.from_bytes(head[20:24], "big")]


def _style(ax, theme, xlabel="x", ylabel="y"):
    ax.set_xlabel(xlabel, fontsize=8, color=theme.foreground)
    ax.set_ylabel(ylabel, fontsize=8, color=theme.foreground)
    ax.tick_params(labelsize=7, colors=theme.foreground)
    ax.set_facecolor(theme.background)
    for s in ax.spines.values():
        s.set_color(theme.grid)


def _field(ax, fig, g2, extent, limit, theme, label, radius, xs=None, ys=None,
           zero=False, glow=0):
    im = ax.imshow(g2, origin="lower", extent=extent, aspect="equal",
                   cmap="RdBu_r", vmin=-limit, vmax=limit,
                   interpolation="nearest")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=7, colors=theme.foreground)
    cb.outline.set_visible(False)
    ax.add_patch(plt.Circle((0, 0), radius, fill=False,
                            color=theme.foreground, linewidth=0.7,
                            linestyle="--", alpha=0.55))
    if zero and xs is not None:
        flat = [v for row in g2 for v in row]
        if min(flat) < 0.0 < max(flat):
            # Glow is extra strokes of the SAME contour at the SAME level.
            for k in range(glow, 0, -1):
                ax.contour(xs, ys, g2, levels=[0.0],
                           colors=[theme.accent_two],
                           linewidths=1.4 + 1.8 * k, alpha=0.10)
            ax.contour(xs, ys, g2, levels=[0.0], colors=[theme.accent_two],
                       linewidths=1.4)
    ax.set_title(label, fontsize=9.5, color=theme.foreground)
    _style(ax, theme)
    return (-limit, limit)


def _status(ax, g2, extent, theme, label, radius):
    cmap = ListedColormap([A6.STATUS_COLOR[s] for s in A6.STATUS_ORDER])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    ax.imshow(g2, origin="lower", extent=extent, aspect="equal", cmap=cmap,
              norm=norm, interpolation="nearest")
    ax.add_patch(plt.Circle((0, 0), radius, fill=False,
                            color=theme.foreground, linewidth=0.7,
                            linestyle="--", alpha=0.55))
    ax.set_title(label, fontsize=9.5, color=theme.foreground)
    _style(ax, theme)


def _running(ax, theme, xs, ys, i, label, color, logy=False):
    """A curve revealed up to the current frame, on FIXED axes."""
    plot = ax.semilogy if logy else ax.plot
    plot(xs, ys, lw=1.0, color=theme.grid, alpha=0.55)
    plot(xs[:i + 1], ys[:i + 1], lw=1.8, color=color)
    ax.plot([xs[i]], [ys[i]], marker="o", markersize=5, color=color,
            linestyle="none")
    ax.set_title(label, fontsize=9.5, color=theme.foreground)
    _style(ax, theme, r"$v_s$", "")


# --- The animations ----------------------------------------------------------

def make_energy_factory(ctx, theme, scales):
    meta, frames = ctx["meta"], ctx["frames"]
    xs, ys = ctx["xs"], ctx["ys"]
    extent = ctx["extent"]
    R = float(meta["radius"])
    v = [float(r["v_s"]) for r in ctx["sweep"]]
    emin = [float(r["E_minus"]) for r in ctx["sweep"]]

    def factory(i):
        rows = frames[i]
        fig, axes = plt.subplots(1, 3, figsize=(16.4, 6.1),
                                 facecolor=theme.background)
        lim = _field(axes[0], fig, A6.to_grid(rows, "rho", xs, ys), extent,
                     ctx["limits"]["rho"], theme,
                     r"Eulerian $\rho$  (fixed scale across all frames)", R,
                     xs=xs, ys=ys, zero=True, glow=theme.glow_passes)
        scales.setdefault("rho", []).append(lim)
        _status(axes[1], A6.status_grid(rows, "rho_sign_status", xs, ys),
                extent, theme, "sign certification of the SAME samples", R)
        _running(axes[2], theme, v, emin, i,
                 r"$E_-=\int\sqrt{\gamma}\max(0,-\rho)\,d^3x$ (exact)",
                 theme.uncertainty_color, logy=True)

        banner(fig, theme,
               "A6  Eulerian energy density through the $g_{00}=0$ crossing",
               crossing_status(meta, i),
               "THIS IS AN OBSERVER-MEASURED ENERGY DENSITY FOR THE EULERIAN "
               "OBSERVER. A NEGATIVE MEASURED VALUE IS NOT A CERTIFIED "
               "NEGATIVE VALUE.")
        stamp(fig, theme, overlay_lines(meta, ctx["sweep"], ctx["nec"],
                                        ctx["integ"], i, "Eulerian (ADM normal)",
                                        "timelike at every v_s in this sweep"))
        footer(fig, theme)
        fig.tight_layout(rect=(0.205, 0.065, 1, 0.855))
        return fig
    return factory


def make_nec_factory(ctx, theme, scales):
    meta, frames = ctx["meta"], ctx["frames"]
    xs, ys = ctx["xs"], ctx["ys"]
    extent = ctx["extent"]
    R = float(meta["radius"])
    v = [float(r["v_s"]) for r in ctx["sweep"]]
    j = [float(r["J_NEC"]) for r in ctx["integ"]]

    def factory(i):
        rows = frames[i]
        fig, axes = plt.subplots(1, 3, figsize=(16.4, 6.1),
                                 facecolor=theme.background)
        lim = _field(axes[0], fig, A6.to_grid(rows, "N_min", xs, ys), extent,
                     ctx["limits"]["N_min"], theme,
                     r"$N_{\min}$ over the whole unit sphere (fixed scale)", R,
                     xs=xs, ys=ys, zero=True, glow=theme.glow_passes)
        scales.setdefault("N_min", []).append(lim)
        _status(axes[1], A6.status_grid(rows, "NEC_status", xs, ys), extent,
                theme, "NEC certification of the SAME samples", R)
        _running(axes[2], theme, v, j, i,
                 r"$J_{\rm NEC}=\int\sqrt{\gamma}\max(0,-N_{\min})\,d^3x$",
                 theme.uncertainty_color)

        banner(fig, theme,
               "A6  null-energy-condition violation through the "
               "$g_{00}=0$ crossing",
               crossing_status(meta, i),
               "MINIMISED OVER THE WHOLE UNIT SPHERE, NOT SAMPLED AT SIX AXES. "
               "AN AXIS-ONLY SEARCH CAN CERTIFY A VIOLATED NEC AS SATISFIED.")
        stamp(fig, theme, overlay_lines(meta, ctx["sweep"], ctx["nec"],
                                        ctx["integ"], i,
                                        "Eulerian orthonormal frame",
                                        "null directions; frame valid at every "
                                        "v_s"))
        footer(fig, theme)
        fig.tight_layout(rect=(0.205, 0.065, 1, 0.855))
        return fig
    return factory


def make_direction_factory(ctx, theme, scales):
    meta = ctx["meta"]
    R = float(meta["radius"])
    wnd = ctx["wnd"]
    lim = max(abs(float(r["N_min"])) for r in wnd) or 1.0
    v = [float(r["v_s"]) for r in ctx["sweep"]]
    nmins = [float(r["N_min"]) for r in ctx["nec"]]

    def factory(i):
        rows = [r for r in wnd if int(r["frame"]) == i]
        fig, axes = plt.subplots(1, 2, figsize=(14.4, 6.1),
                                 facecolor=theme.background)
        ax = axes[0]
        q = ax.quiver([float(r["x"]) for r in rows],
                      [float(r["y"]) for r in rows],
                      [float(r["worst_direction_x"]) for r in rows],
                      [float(r["worst_direction_y"]) for r in rows],
                      [float(r["N_min"]) for r in rows],
                      cmap="RdBu_r", clim=(-lim, lim), angles="xy", scale=26,
                      width=0.005, pivot="mid")
        scales.setdefault("direction", []).append((-lim, lim))
        cb = fig.colorbar(q, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=7, colors=theme.foreground)
        cb.outline.set_visible(False)
        ax.add_patch(plt.Circle((0, 0), R, fill=False, color=theme.foreground,
                                linewidth=0.8, linestyle="--", alpha=0.55))
        ax.set_xlim(-float(meta["domain_extent"]), float(meta["domain_extent"]))
        ax.set_ylim(-float(meta["domain_extent"]), float(meta["domain_extent"]))
        ax.set_aspect("equal")
        ax.set_title("worst null direction, coloured by $N_{\\min}$ "
                     "(fixed scale)", fontsize=9.5, color=theme.foreground)
        _style(ax, theme)

        _running(axes[1], theme, v, [abs(x) for x in nmins], i,
                 r"$|N_{\min}|$ on the slice", theme.accent_two)

        banner(fig, theme,
               "A6  the direction that minimises the null energy",
               crossing_status(meta, i),
               "THE MINIMISER IS OFF THE COORDINATE AXES WHEREVER THE MOMENTUM "
               "DENSITY IS NONZERO.")
        stamp(fig, theme, overlay_lines(meta, ctx["sweep"], ctx["nec"],
                                        ctx["integ"], i,
                                        "Eulerian orthonormal frame",
                                        "null directions; frame valid at every "
                                        "v_s"))
        footer(fig, theme)
        fig.tight_layout(rect=(0.235, 0.065, 1, 0.855))
        return fig
    return factory


def make_observer_factory(ctx, theme, scales):
    meta = ctx["meta"]
    obs = ctx["obs"]
    sites = sorted({r["site"] for r in obs})
    vals = [r["measured"] for r in obs if r["measured"] is not None]
    lo, hi = min(vals), max(vals)
    pad = 0.08 * (hi - lo or 1.0)
    v = [float(r["v_s"]) for r in ctx["sweep"]]
    g00c = [float(r["g00_center"]) for r in ctx["sweep"]]

    def factory(i):
        fig, axes = plt.subplots(1, 3, figsize=(15.8, 6.1),
                                 facecolor=theme.background)
        for c, site in enumerate(sites):
            ax = axes[c]
            for k, fam in enumerate(A6.OBSERVER_FAMILIES):
                col = A6.OBSERVER_STYLE[fam][0]
                row = next((r for r in obs if r["site"] == site
                            and r["observer_type"] == fam and r["frame"] == i),
                           None)
                if row is None:
                    continue
                if row["measured"] is None:
                    # No bar, because there is no value. The reason is drawn
                    # instead, in the uncertainty colour.
                    ax.text(k, 0.0, "INADMISSIBLE", rotation=90, fontsize=7,
                            ha="center", va="center",
                            color=theme.uncertainty_color)
                    ax.text(k, lo, row["reason"].split(":")[0], rotation=90,
                            fontsize=5.4, ha="center", va="bottom",
                            color=theme.uncertainty_color, alpha=0.9)
                else:
                    ax.bar([k], [row["measured"]], width=0.62, color=col,
                           alpha=0.55 if theme.trail else 0.85,
                           edgecolor=col, linewidth=1.4)
            ax.axhline(0.0, color=theme.grid, lw=0.9)
            ax.set_ylim(lo - pad, hi + pad)
            ax.set_xticks(range(len(A6.OBSERVER_FAMILIES)))
            ax.set_xticklabels([f.replace("_", "\n")
                                for f in A6.OBSERVER_FAMILIES], fontsize=6.6,
                               color=theme.foreground)
            ax.set_title("%s   $T_{\\mu\\nu}u^\\mu u^\\nu$" % site,
                         fontsize=9.5, color=theme.foreground)
            _style(ax, theme, "", "measured energy density")

        _running(axes[2], theme, v, g00c, i,
                 r"$g_{00}$ at the bubble centre", theme.foreground)
        axes[2].axhline(0.0, color=theme.uncertainty_color, lw=1.0, ls="--")

        banner(fig, theme,
               "A6  four observers, one geometry, four different energy "
               "densities",
               crossing_status(meta, i),
               "AN INADMISSIBLE OBSERVER MEASURES NOTHING. THE RECORD CONTAINS "
               "NO VALUE FOR IT AND NONE IS DRAWN. $|g_{00}|$ IS NEVER "
               "SUBSTITUTED.")
        state = "eulerian valid; coordinate-stationary valid at %s of 2 sites" \
                % ctx["sweep"][i]["stationary_valid_sites"]
        stamp(fig, theme, overlay_lines(meta, ctx["sweep"], ctx["nec"],
                                        ctx["integ"], i,
                                        "all four families, separately",
                                        state))
        footer(fig, theme)
        fig.tight_layout(rect=(0.240, 0.065, 1, 0.855))
        scales.setdefault("observer", []).append((lo - pad, hi + pad))
        return fig
    return factory


def make_integrated_factory(ctx, theme, scales):
    meta = ctx["meta"]
    integ = ctx["integ"]
    v = [float(r["v_s"]) for r in integ]
    series = {
        "E_plus": [float(r["E_plus"]) for r in integ],
        "E_minus": [float(r["E_minus"]) for r in integ],
        "E_abs": [float(r["E_abs"]) for r in integ],
        "J_NEC": [float(r["J_NEC"]) for r in integ],
    }
    smooth = [float(r["smooth_surrogate_value"]) for r in integ]
    ylo = min(min(s) for s in series.values() if min(s) > 0)
    yhi = max(max(s) for s in series.values())

    def factory(i):
        fig, axes = plt.subplots(1, 3, figsize=(15.8, 6.1),
                                 facecolor=theme.background)
        ax = axes[0]
        for key, col in (("E_plus", theme.accent),
                         ("E_minus", theme.uncertainty_color),
                         ("E_abs", theme.accent_two),
                         ("J_NEC", "#7ed957")):
            ax.semilogy(v, series[key], lw=1.0, color=theme.grid, alpha=0.5)
            ax.semilogy(v[:i + 1], series[key][:i + 1], lw=1.8, color=col,
                        label=key)
            ax.plot([v[i]], [series[key][i]], marker="o", markersize=5,
                    color=col, linestyle="none")
        ax.set_ylim(ylo * 0.5, yhi * 2.0)
        ax.legend(fontsize=7.2, frameon=False, labelcolor=theme.foreground)
        ax.set_title("exact integrated objectives (fixed axes)", fontsize=9.5,
                     color=theme.foreground)
        _style(ax, theme, r"$v_s$", "value")

        ax = axes[1]
        ax.semilogy(v, series["E_minus"], lw=1.0, color=theme.grid, alpha=0.5)
        ax.semilogy(v, smooth, lw=1.0, color=theme.grid, alpha=0.5)
        ax.semilogy(v[:i + 1], series["E_minus"][:i + 1], lw=1.8,
                    color=theme.uncertainty_color, label=r"EXACT $E_-$")
        ax.semilogy(v[:i + 1], smooth[:i + 1], lw=1.6, ls=":",
                    color=theme.accent, label="softplus surrogate")
        ax.set_ylim(min(ylo, min(smooth)) * 0.5, max(yhi, max(smooth)) * 2.0)
        ax.legend(fontsize=7.2, frameon=False, labelcolor=theme.foreground)
        ax.set_title("exact vs smooth -- NEVER interchanged", fontsize=9.5,
                     color=theme.foreground)
        _style(ax, theme, r"$v_s$", "value")

        ax = axes[2]
        r = integ[i]
        parts = [float(r["e_minus_interior"]), float(r["e_minus_shell"]),
                 float(r["tail_estimate"])]
        ax.bar([0, 1, 2], parts, width=0.6,
               color=[theme.accent, theme.accent_two, "#7ed957"],
               alpha=0.55 if theme.trail else 0.85)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["interior", "boundary\nshell", "exterior\ntail"],
                           fontsize=7, color=theme.foreground)
        ax.set_ylim(0.0, max(max(float(q["e_minus_interior"]) for q in integ),
                             1e-12) * 1.15)
        ax.set_title("where $E_-$ comes from  (boundary share %.4f)"
                     % float(r["boundary_fraction"]), fontsize=9.5,
                     color=theme.foreground)
        _style(ax, theme, "", "contribution")
        scales.setdefault("integrated", []).append((ylo * 0.5, yhi * 2.0))

        banner(fig, theme, "A6  integrated objectives an optimizer would "
                           "consume", crossing_status(meta, i),
               "$E_+$ AND $E_-$ ARE REPORTED SEPARATELY SO $E_{\\rm net}$ "
               "CANNOT HIDE THE CANCELLATION. CLOSURE RESIDUAL %.1e."
               % max(float(r["net_closure_residual"]),
                     float(r["abs_closure_residual"])))
        stamp(fig, theme, overlay_lines(meta, ctx["sweep"], ctx["nec"], integ,
                                        i, "Eulerian (ADM normal)",
                                        "optimization-ready: %s"
                                        % ("yes" if r["optimization_ready"] ==
                                           "1" else "NO")))
        footer(fig, theme)
        fig.tight_layout(rect=(0.215, 0.065, 1, 0.855))
        return fig
    return factory


def make_certification_factory(ctx, theme, scales):
    meta, frames = ctx["meta"], ctx["frames"]
    xs, ys = ctx["xs"], ctx["ys"]
    extent = ctx["extent"]
    R = float(meta["radius"])

    def factory(i):
        rows = frames[i]
        fig, axes = plt.subplots(1, 4, figsize=(16.8, 6.1),
                                 facecolor=theme.background)
        for ax, key, lab in ((axes[0], "NEC_status", "NEC"),
                             (axes[1], "WEC_status", "WEC"),
                             (axes[2], "DEC_status", "DEC"),
                             (axes[3], "SEC_status", "SEC")):
            c = {st: sum(1 for r in rows if r[key] == st)
                 for st in A6.STATUS_ORDER}
            _status(ax, A6.status_grid(rows, key, xs, ys), extent, theme,
                    "%s\nviolated %d  indeterminate %d  satisfied %d"
                    % (lab, c[A6.CERTIFIED_VIOLATED], c[A6.INDETERMINATE],
                       c[A6.CERTIFIED_SATISFIED]), R)
        scales.setdefault("certification", []).append((-0.5, 3.5))

        banner(fig, theme,
               "A6  pointwise energy-condition certification against $v_s$",
               crossing_status(meta, i),
               "THERE IS NO \"PASS\". A VALUE WITHIN %s ERROR ENVELOPES OF "
               "ZERO IS REPORTED AS INDETERMINATE."
               % meta["certification_safety_factor"])
        stamp(fig, theme, overlay_lines(meta, ctx["sweep"], ctx["nec"],
                                        ctx["integ"], i,
                                        "worst case over the timelike "
                                        "velocity ball",
                                        "per condition, deterministic sweep"))
        footer(fig, theme)
        fig.tight_layout(rect=(0.200, 0.065, 1, 0.845))
        return fig
    return factory


# --- Driver ------------------------------------------------------------------

ANIMATIONS = [
    ("energy_density_vs_velocity", make_energy_factory, ("scientific",
                                                         "cinematic"),
     "energy_velocity_sweep.csv"),
    ("nec_violation_vs_velocity", make_nec_factory, ("scientific",
                                                     "cinematic"),
     "nec_velocity_sweep.csv"),
    ("worst_null_direction_vs_velocity", make_direction_factory,
     ("scientific",), "worst_null_directions.csv"),
    ("observer_energy_comparison", make_observer_factory, ("scientific",),
     "observer_comparison.csv"),
    ("integrated_energy_objectives", make_integrated_factory, ("scientific",),
     "integrated_objectives.csv"),
    ("energy_condition_certification", make_certification_factory,
     ("scientific",), "energy_velocity_sweep.csv"),
]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fps", type=int, default=FPS)
    args = ap.parse_args(argv)

    sweep = records.resolve(args.sweep)
    out = records.ensure_dir(args.out)

    meta = A6.load_manifest(sweep)
    before = A6.source_checksums(sweep)

    frames = A6.load_frames(sweep)
    sweep_rows = records.read_csv(os.path.join(sweep,
                                               "energy_velocity_sweep.csv"))
    necsweep = records.read_csv(os.path.join(sweep, "nec_velocity_sweep.csv"))
    integ = records.read_csv(os.path.join(sweep, "integrated_objectives.csv"))
    wnd = records.read_csv(os.path.join(sweep, "worst_null_directions.csv"))
    obs = A6.observer_rows(sweep)

    if not (len(frames) == len(sweep_rows) == len(necsweep) == len(integ)):
        print("frame files and summary rows disagree: %d/%d/%d/%d"
              % (len(frames), len(sweep_rows), len(necsweep), len(integ)),
              file=sys.stderr)
        return 1

    xs, ys = A6.axes_of(frames[0])
    ctx = {
        "meta": meta, "frames": frames, "sweep": sweep_rows, "nec": necsweep,
        "integ": integ, "wnd": wnd, "obs": obs, "xs": xs, "ys": ys,
        "extent": [min(xs), max(xs), min(ys), max(ys)],
        # Computed ONCE from every frame and frozen. Both themes receive the
        # same object, so the cinematic theme cannot change a colour limit even
        # by accident.
        "limits": {
            "rho": A6.global_symmetric_limit(frames, "rho"),
            "N_min": A6.global_symmetric_limit(frames, "N_min"),
        },
    }

    entries = []
    scales = {}
    blanks = []
    encoder = render.find_ffmpeg()
    for stem, factory_of, theme_names, source in ANIMATIONS:
        for tname in theme_names:
            theme = themes.get(tname)
            themes.apply(theme)
            name = stem if len(theme_names) == 1 else "%s_%s" % (stem, tname)
            frame_dir = os.path.join(out, "_frames_" + name)
            paths = render.save_frames(factory_of(ctx, theme, scales),
                                       len(frames), frame_dir,
                                       prefix="a6", dpi=110)
            blanks += [p for p in paths if render.frame_is_blank(p)]

            gif = render.write_gif(paths, os.path.join(out, name + ".gif"),
                                   fps=args.fps, max_frames=len(frames),
                                   max_width=900)
            mp4 = render.write_mp4(paths, os.path.join(out, name + ".mp4"),
                                   fps=args.fps, frame_dir=frame_dir)
            sheet = render.contact_sheet(
                paths, os.path.join(out, "contact_sheets", name + ".png"),
                tiles=9, label_prefix="v_s frame")

            src = os.path.join(sweep, source)
            res_src = paths[0]
            for path, codec in ((gif, "gif"), (mp4, "h264")):
                e = render.animation_entry(
                    name, path, codec=codec, resolution=png_resolution(res_src),
                    fps=args.fps, frames=paths,
                    duration_s=len(frames) / float(args.fps), source_csv=src)
                e["theme"] = tname
                e["colour_scale"] = ("fixed globally across all frames, "
                                     "computed once from the whole sweep")
                e["contact_sheet"] = os.path.basename(sheet) if sheet else None
                if path is None:
                    # Required wording: an absent encoder is a recorded fact,
                    # not a silent gap, and it must not fail the scientific
                    # gate. The GIF, the frames and the contact sheet are the
                    # deterministic products.
                    e["reason"] = ("encoder unavailable" if codec == "h264"
                                   and not encoder else
                                   "writer unavailable for this format")
                entries.append(e)

            # No raw frame directory may survive the render.
            for p in paths:
                os.remove(p)
            os.rmdir(frame_dir)
            print("wrote %-40s gif=%s mp4=%s sheet=%s"
                  % (name, bool(gif), bool(mp4), bool(sheet)))

    after = A6.source_checksums(sweep)
    prov = {
        "program": "A",
        "milestone": "A6",
        "source_sweep": os.path.basename(os.path.normpath(sweep)),
        "git_commit": meta.get("git_commit"),
        "git_dirty": meta.get("git_dirty"),
        "frame_count": len(frames),
        "frame_rate": args.fps,
        "velocity_sequence": meta.get("velocity_sequence"),
        "exact_onset_frame": meta.get("exact_onset_frame"),
        "first_spatial_contour_frame": meta.get("first_spatial_contour_frame"),
        "colour_scale": ("fixed globally across all frames, computed once from "
                         "the whole velocity sweep, shared by both themes"),
        "colour_limits": ctx["limits"],
        "cinematic_theme_may_change": ["typography", "easing", "glow"],
        "cinematic_theme_may_not_change": [
            "source data", "positions", "timing", "contours", "certification",
            "normalisation", "colour limits"],
        "observer_families": list(A6.OBSERVER_FAMILIES),
        "observer_admissibility": meta.get("observer_admissibility"),
        "g00_zero_meaning": meta.get("g00_zero_meaning"),
        "exact_vs_smooth": meta.get("exact_vs_smooth"),
        "renderer_computes_physics": False,
        "encoder_available": bool(encoder),
        "numerical_error_envelope_floor":
            meta.get("numerical_error_envelope_floor"),
        "certification_safety_factor":
            meta.get("certification_safety_factor"),
        "scope": A6.SCOPE,
    }
    render.write_manifest(os.path.join(out, "animation_manifest.json"),
                          entries, prov)

    # --- Diagnostics ---------------------------------------------------------
    rep = anomaly.AnomalyReport("A", "A6", "a6_animation")

    # 7. Observer identity must not change between frames.
    per_frame = [sorted({r["observer_type"] for r in obs
                         if r["frame"] == i}) for i in range(len(frames))]
    rep.check_identifiers_stable("observer families per frame", per_frame)
    changed_identity = sum(1 for f in per_frame if f != per_frame[0])
    rep.record("observer_identity_changed_between_frames",
               fired=changed_identity > 0, measured=changed_identity,
               threshold=0,
               detail="frames whose observer family set differs from frame 0")

    # 10. A worst-direction jump between frames is only acceptable where the
    # stress eigenvalues are nearly degenerate: with a clear gap the minimiser
    # is unique and cannot move discontinuously.
    jumps = []
    for i in range(1, len(necsweep)):
        a = [float(necsweep[i - 1]["worst_direction_%s" % k])
             for k in "xyz"]
        b = [float(necsweep[i]["worst_direction_%s" % k]) for k in "xyz"]
        # |a - b| for two saved unit vectors; a difference, not a computation.
        d = max(abs(p - q) for p, q in zip(a, b))
        gap = min(float(necsweep[i - 1]["nec_eigen_gap"]),
                  float(necsweep[i]["nec_eigen_gap"]))
        scale = max(abs(float(necsweep[i]["N_min"])), 1e-30)
        if d > 0.5 and gap > 0.05 * scale:
            jumps.append((i, d, gap))
    rep.record("nec_direction_jump_without_degeneracy", fired=bool(jumps),
               measured=len(jumps), threshold=0,
               detail="frame-to-frame direction jumps larger than 0.5 with a "
                      "non-degenerate eigenvalue gap: %r" % jumps[:3])

    # 16. The domain must not change between frames.
    doms = sorted({(r["domain_extent"], r["grid_resolution"])
                   for r in sweep_rows})
    rep.record("domain_changed_between_frames", fired=len(doms) > 1,
               measured=doms[:3], threshold=1,
               detail="distinct (domain_extent, grid_resolution) pairs across "
                      "the sweep")

    # 17. Per-frame colour rescaling. The limits are recorded by every factory
    # on every frame; any spread at all is a defect.
    worst_spread = 0.0
    for key, vals in scales.items():
        spread = max(max(abs(a - vals[0][0]), abs(b - vals[0][1]))
                     for a, b in vals)
        worst_spread = max(worst_spread, spread)
    rep.record("per_frame_color_rescaling", fired=worst_spread > 0.0,
               measured=worst_spread, threshold=0.0,
               detail="largest deviation of any frame's colour limits from "
                      "frame 0, over %d recorded limit sets" % len(scales))

    # 21. Blank frames.
    rep.record("blank_frame", fired=bool(blanks), measured=len(blanks),
               threshold=0,
               detail="frames with fewer than 12 distinct colours: %r"
                      % [os.path.basename(p) for p in blanks[:3]])

    # 19. Source integrity, against the checksums the figure pass recorded.
    recorded = {}
    cpath = os.path.join(out, "source_checksums.json")
    if os.path.exists(cpath):
        with open(cpath, encoding="utf-8") as fh:
            recorded = json.load(fh)
    mismatched = [k for k, h in recorded.items() if before.get(k) != h]
    mismatched += [k for k in before if before[k] != after.get(k)]
    rep.record("source_checksum_mismatch", fired=bool(mismatched),
               measured=len(mismatched), threshold=0,
               detail="records whose SHA-256 differs from the value the figure "
                      "pass recorded, or which changed during this render: %r"
                      % sorted(set(mismatched))[:3])

    # 20. Recomputed physics in this package.
    offenders = A6.scan_python_physics(A6.package_python_files())
    rep.record("python_physics_recomputation", fired=bool(offenders),
               measured=len(offenders), threshold=0,
               detail="AST scan of the A6 renderer package: %r"
                      % offenders[:3])

    path = rep.write(out)
    left = [n for n in os.listdir(out) if n.startswith("_frames_")]
    print("%d animation entries, %d fired anomalies, %d raw frame dirs left "
          "-> %s" % (len(entries), len(rep.fired), len(left),
                     os.path.basename(path)))
    if left:
        print("raw frame directories survived: %r" % left, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
