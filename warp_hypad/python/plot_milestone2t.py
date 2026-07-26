#!/usr/bin/env python3
"""Milestone A2T figures: temporal derivative operators.

Usage:
    python plot_milestone2t.py <results_dir> [--out <figures_dir>]

Inputs:  temporal_convergence.csv, identity_profile.csv
Outputs: milestone02t_convergence.png, milestone02t_identity.png
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh_style as S


def read_csv(p):
    with p.open(newline="") as fh:
        return list(csv.DictReader(fh))


def plot_convergence(rows, out):
    groups = defaultdict(list)
    for r in rows:
        groups[(r["quantity"], r["stencil_position"])].append(r)

    quantities = ["d/dt", "d2/dt2"]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), facecolor=S.SURFACE)
    S.title_block(fig, "Temporal derivative convergence",
                  "Three stencil positions per operator. 'Interior' samples both "
                  "sides; the window forms are the one-sided stencils a caller "
                  "restricted to a finite time window must use.", y=0.955)

    positions = ["interior", "window_start", "window_end"]
    for ax, q in zip(axes, quantities):
        S.style_axes(ax, logx=True, logy=True)
        anchor = 0.0
        for i, pos in enumerate(positions):
            v = sorted(groups[(q, pos)], key=lambda r: float(r["h"]))
            if not v:
                continue
            hs = [float(r["h"]) for r in v]
            es = [float(r["max_error"]) for r in v]
            anchor = max(anchor, max(es))
            ax.plot(hs, es, marker="o", markersize=3.6, linewidth=1.6,
                    color=S.SERIES[i], label=pos.replace("_", " "), zorder=3)
        h0, h1 = 0.008, 0.1
        top = anchor * 1.6
        ax.plot([h0, h1], [top * (h0 / h1) ** 2, top], color=S.INK2,
                linestyle=(0, (5, 3)), linewidth=1.6, zorder=2)
        ax.annotate("slope 2", xy=(h1, top), xytext=(-6, -12),
                    textcoords="offset points", ha="right", fontsize=8,
                    color=S.INK2)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(q, fontsize=10, color=S.INK, loc="left")
        ax.set_xlabel("time step h (dimensionless)", fontsize=8, color=S.INK2)
        ax.set_ylabel("max |error|", fontsize=8, color=S.INK2)
        ax.legend(fontsize=7.5, frameon=False, labelcolor=S.INK2)

    S.footer(fig, "Milestone A2T  |  manufactured field "
                  "f = sin(wt+p)cos(ax+p)exp(by+cz), dimensionless. "
                  "Time is not gridded: the operator re-evaluates the closed "
                  "form at shifted times.")
    S.save(fig, out, rect=(0, 0.11, 1, 0.87))


def plot_identity(rows, out):
    xs = [float(r["x"]) for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0), facecolor=S.SURFACE)
    S.title_block(fig, "Time dependence of the moving warp metric",
                  r"Rigid translation gives $\partial_t g = -v_s\,\partial_x g$ "
                  "exactly, which is why the inertial-frame formulation costs no "
                  "more than a spatial derivative.", y=0.955)

    ax = axes[0]
    S.style_axes(ax)
    ax.plot(xs, [float(r["dt_g00_analytic"]) for r in rows], color=S.SERIES[0],
            linewidth=2.0, label=r"$\partial_t g_{00}$", zorder=3)
    ax.plot(xs, [float(r["dt_g0x_analytic"]) for r in rows], color=S.SERIES[1],
            linewidth=2.0, label=r"$\partial_t g_{0x}$", zorder=3)
    ax.axhline(0.0, color=S.INK3, linewidth=0.8, zorder=1)
    ax.set_xlabel("x (grid length unit)", fontsize=8, color=S.INK2)
    ax.set_ylabel("temporal derivative", fontsize=8, color=S.INK2)
    ax.set_title("analytic time derivative at t = 2", fontsize=9.5,
                 color=S.INK, loc="left")
    ax.legend(fontsize=8, frameon=False, labelcolor=S.INK2)

    ax = axes[1]
    S.style_axes(ax)
    ax.plot(xs, [float(r["g00_moving"]) for r in rows], color=S.SERIES[0],
            linewidth=2.0, label="time-dependent (correct)", zorder=3)
    ax.plot(xs, [float(r["g00_frozen"]) for r in rows], color=S.STATUS_CRIT,
            linewidth=2.0, linestyle=(0, (5, 3)),
            label="frozen diagnostic (wrong at t>0)", zorder=3)
    ax.set_xlabel("x (grid length unit)", fontsize=8, color=S.INK2)
    ax.set_ylabel(r"$g_{00}$", fontsize=8, color=S.INK2)
    ax.set_title("what the diagnostic mode discards, at t = 2", fontsize=9.5,
                 color=S.INK, loc="left")
    ax.legend(fontsize=7.5, frameon=False, labelcolor=S.INK2)

    S.footer(fig, "Milestone A2T  |  Alcubierre v = 0.35, R = 3, sigma = 2, "
                  "geometrized units. The bubble has moved 0.7 length units by "
                  "t = 2; the frozen mode does not know that, which is why "
                  "curvature must refuse it.")
    S.save(fig, out, rect=(0, 0.12, 1, 0.87))


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

    for n in ("temporal_convergence.csv", "identity_profile.csv"):
        if not (rd / n).is_file():
            print(f"error: missing {rd/n}; run milestone2t_verify first",
                  file=sys.stderr)
            return 2

    plot_convergence(read_csv(rd / "temporal_convergence.csv"),
                     out / "milestone02t_convergence.png")
    plot_identity(read_csv(rd / "identity_profile.csv"),
                  out / "milestone02t_identity.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
