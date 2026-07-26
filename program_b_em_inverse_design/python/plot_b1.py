#!/usr/bin/env python3
"""Program B1 figures: metric <-> constitutive mapping.

Reads only files written by b1_verify; regenerates every figure without
rerunning the solver.

Usage:
    python plot_b1.py <results_dir> [--out <figures_dir>]

Inputs:  warp_medium_profile.csv, realizability_map.csv
         medium_sensitivity.csv (optional; OTIlib backend only)
Outputs: b1_mapping_diagram.png
         b1_warp_medium.png
         b1_realizability.png
         b1_medium_sensitivity.png (when available)
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "warp_hypad" / "python"))
import wh_style as S

UNITS = ("Heaviside-Lorentz, c = 1 (vacuum: eps = mu = I)  |  signature -+++  |  "
         "Alcubierre v = 0.2, R = 3, sigma = 2 unless noted")


def read_csv(p):
    with p.open(newline="") as fh:
        return list(csv.DictReader(fh))


def plot_mapping_diagram(out):
    """The B1 mapping, including which arrows are implemented and which are not."""
    fig, ax = plt.subplots(figsize=(10.0, 4.8), facecolor=S.SURFACE)
    ax.set_xlim(0, 10); ax.set_ylim(0, 4.8); ax.axis("off")
    S.title_block(fig, "B1: the metric / constitutive correspondence",
                  "Solid arrows are implemented and verified. Dashed arrows are "
                  "mathematically valid but out of B1 scope. Crossed arrows do "
                  "not exist in general.", y=0.965)

    def box(x, y, w, h, text, colour, fill):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                                    boxstyle="round,pad=0.10,rounding_size=0.10",
                                    linewidth=1.6, edgecolor=colour,
                                    facecolor=fill, zorder=2))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=9.5, color=S.INK, zorder=3, linespacing=1.5)

    box(0.3, 2.6, 2.6, 1.15, r"$g_{\mu\nu}$" + "\nspacetime metric",
        S.SERIES[0], "#eaf2fd")
    box(3.9, 2.6, 2.6, 1.15, r"$\varepsilon,\ \mu,\ \xi,\ \zeta$" + "\nconstitutive response",
        S.SERIES[1], "#fdeee7")
    box(7.4, 2.6, 2.3, 1.15, "conformal class\n" + r"$[\,\hat g_{\mu\nu}\,]$",
        S.SERIES[2], "#e7f7f1")
    box(3.9, 0.5, 2.6, 1.15, "non-metric media\nskewon / axion /\nbirefringent",
        S.STATUS_CRIT, "#fdecec")

    def arrow(x1, y1, x2, y2, style, colour):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
                                     arrowstyle="-|>", mutation_scale=16,
                                     linewidth=1.8, color=colour,
                                     linestyle=style, zorder=4,
                                     shrinkA=2, shrinkB=2))

    arrow(2.95, 3.35, 3.85, 3.35, "-", S.SERIES[0])
    ax.text(3.4, 3.55, "Plebanski-Tamm", ha="center", fontsize=7.5, color=S.INK2)
    ax.text(3.4, 3.02, "verified", ha="center", fontsize=7, color=S.STATUS_GOOD)

    arrow(6.55, 3.05, 7.35, 3.05, "-", S.SERIES[2])
    ax.text(6.95, 3.28, r"$w=0$", ha="center", fontsize=7.5, color=S.INK2)
    ax.text(6.95, 2.72, "verified", ha="center", fontsize=7, color=S.STATUS_GOOD)

    arrow(5.4, 2.55, 5.4, 1.75, (0, (4, 3)), S.STATUS_WARN)
    ax.text(5.6, 2.15, "detected and reported,\nnever fitted", ha="left",
            fontsize=7.5, color=S.INK2)

    ax.text(0.3, 2.15, "conformal freedom\n" + r"$g \to \Omega^2 g$" +
            " leaves the\nmedium unchanged", fontsize=8, color=S.INK2,
            va="top", linespacing=1.5)

    ax.text(7.4, 2.15, "recovered only up to\n" + r"$\Omega^2$" +
            "; compare null cones,\nnot components", fontsize=8, color=S.INK2,
            va="top", linespacing=1.5)

    ax.text(0.3, 0.95, "NOT IMPLEMENTED in B1:\n"
                       r"  $w \neq 0$ inverse (every warp metric)" "\n"
                       "  covariant chi, Fresnel surface, closure relation\n"
                       "  dispersion, loss, passivity, causality",
            fontsize=7.5, color=S.STATUS_CRIT, va="top", linespacing=1.6)

    S.footer(fig, "Program B1  |  " + UNITS)
    S.save(fig, out, rect=(0, 0.06, 1, 0.89))


def plot_warp_medium(rows, out):
    xs = [float(r["x"]) for r in rows if r["medium_exists"] == "1"]
    exx = [float(r["eps_xx"]) for r in rows if r["medium_exists"] == "1"]
    eyy = [float(r["eps_yy"]) for r in rows if r["medium_exists"] == "1"]
    wx = [float(r["w_x"]) for r in rows if r["medium_exists"] == "1"]

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0), facecolor=S.SURFACE)
    S.title_block(fig, "The warp metric as an electromagnetic medium",
                  "Centreline through the bubble. The medium is bianisotropic: "
                  "no arrangement of ordinary dielectrics reproduces it.",
                  y=0.955)

    ax = axes[0]
    S.style_axes(ax)
    ax.plot(xs, exx, color=S.SERIES[0], linewidth=2.0,
            label=r"$\varepsilon^{xx} = \mu^{xx}$", zorder=3)
    ax.plot(xs, eyy, color=S.SERIES[1], linewidth=2.0,
            label=r"$\varepsilon^{yy} = \mu^{yy}$", zorder=3)
    ax.axhline(1.0, color=S.INK3, linestyle=(0, (4, 3)), linewidth=1.2, zorder=2)
    ax.annotate("vacuum", xy=(xs[0], 1.0), xytext=(3, 4),
                textcoords="offset points", fontsize=7, color=S.INK2)
    ax.set_xlabel("x (grid length unit)", fontsize=8, color=S.INK2)
    ax.set_ylabel("permittivity = permeability", fontsize=8, color=S.INK2)
    ax.set_title("impedance matched by construction", fontsize=9.5,
                 color=S.INK, loc="left")
    ax.legend(fontsize=8, frameon=False, labelcolor=S.INK2)

    ax = axes[1]
    S.style_axes(ax)
    ax.plot(xs, wx, color=S.SERIES[2], linewidth=2.0, zorder=3)
    ax.axhline(0.0, color=S.INK3, linewidth=0.8, zorder=1)
    ax.set_xlabel("x (grid length unit)", fontsize=8, color=S.INK2)
    ax.set_ylabel(r"$w_x$  (magnetoelectric)", fontsize=8, color=S.INK2)
    ax.set_title("bianisotropy, nonzero only at the wall", fontsize=9.5,
                 color=S.INK, loc="left")

    S.footer(fig, "Program B1  |  " + UNITS +
             "\nSign of w is this project's internal convention and is not "
             "verified against the literature; magnitude is well determined.")
    S.save(fig, out, rect=(0, 0.13, 1, 0.87))


def plot_realizability(rows, out):
    vs = sorted({float(r["v_s"]) for r in rows})
    xs = sorted({float(r["x"]) for r in rows})
    vi = {v: i for i, v in enumerate(vs)}
    xi = {v: i for i, v in enumerate(xs)}
    grid = [[0.0] * len(xs) for _ in vs]
    for r in rows:
        grid[vi[float(r["v_s"])]][xi[float(r["x"])]] = float(r["medium_exists"])

    fig, ax = plt.subplots(figsize=(7.4, 5.0), facecolor=S.SURFACE)
    S.title_block(fig, "Where an equivalent medium exists",
                  r"The Plebanski map needs $g_{00} < 0$, i.e. $v_s f < 1$. "
                  "Dark = a medium exists; light = none, in this frame.",
                  y=0.96)
    S.style_axes(ax)
    S.imshow_field(ax, xs, vs, grid, signed=False, label="medium exists (1/0)")
    ax.axhline(1.0, color=S.STATUS_CRIT, linestyle=(0, (5, 3)), linewidth=1.8,
               zorder=5)
    ax.annotate(r"$v_s = 1$", xy=(xs[3], 1.0), xytext=(4, 5),
                textcoords="offset points", fontsize=9, color=S.STATUS_CRIT)
    ax.set_xlabel("x (grid length unit)", fontsize=8, color=S.INK2)
    ax.set_ylabel(r"bubble velocity $v_s$", fontsize=8, color=S.INK2)

    S.footer(fig, "Program B1  |  " + UNITS +
             "\nTHE central realizability limitation: the superluminal regime "
             "of interest is exactly where the correspondence fails. This is "
             "the same v_s f > 1 condition that costs Program A its coordinate "
             "observers.")
    S.save(fig, out, rect=(0, 0.13, 1, 0.88))


def plot_sensitivity(rows, out):
    xs = [float(r["x"]) for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0), facecolor=S.SURFACE)
    S.title_block(fig, "HyPAD sensitivities of the effective medium",
                  "Derivatives of the constitutive response with respect to the "
                  "metric parameters, from one evaluation with OTIScalar.",
                  y=0.955)

    for ax, keys, title, ylab in (
        (axes[0], ["deps_yy_dv", "deps_yy_dR", "deps_yy_dsigma"],
         r"$\partial \varepsilon^{yy} / \partial p$", "sensitivity"),
        (axes[1], ["dw_dv", "dw_dR", "dw_dsigma"],
         r"$\partial w_x / \partial p$", "sensitivity"),
    ):
        S.style_axes(ax)
        for i, (k, lbl) in enumerate(zip(keys, ["v", "R", r"$\sigma$"])):
            ax.plot(xs, [float(r[k]) for r in rows], color=S.SERIES[i],
                    linewidth=1.8, label=lbl, zorder=3)
        ax.axhline(0.0, color=S.INK3, linewidth=0.8, zorder=1)
        ax.set_xlabel("x (grid length unit)", fontsize=8, color=S.INK2)
        ax.set_ylabel(ylab, fontsize=8, color=S.INK2)
        ax.set_title(title, fontsize=9.5, color=S.INK, loc="left")
        ax.legend(fontsize=8, frameon=False, labelcolor=S.INK2, title="w.r.t.",
                  title_fontsize=7)

    S.footer(fig, "Program B1  |  " + UNITS +
             "\nVerified against centred finite differences to 9.1e-08 and, for "
             "the conformal parameter, against the exact statement d eps/d Omega = 0.")
    S.save(fig, out, rect=(0, 0.13, 1, 0.87))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("results_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()

    rd: Path = a.results_dir
    if not rd.is_dir():
        print(f"error: not a directory: {rd}", file=sys.stderr)
        return 2
    out = a.out or (rd / "figures")
    out.mkdir(parents=True, exist_ok=True)

    plot_mapping_diagram(out / "b1_mapping_diagram.png")

    for name, fn, png in (
        ("warp_medium_profile.csv", plot_warp_medium, "b1_warp_medium.png"),
        ("realizability_map.csv", plot_realizability, "b1_realizability.png"),
    ):
        p = rd / name
        if not p.is_file():
            print(f"error: missing {p}; run b1_verify first", file=sys.stderr)
            return 2
        fn(read_csv(p), out / png)

    sp = rd / "medium_sensitivity.csv"
    if sp.is_file():
        plot_sensitivity(read_csv(sp), out / "b1_medium_sensitivity.png")
    else:
        print("skipped sensitivity figure (OTIlib backend disabled)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
