#!/usr/bin/env python3
"""Milestone A2 figures: finite-difference operators.

Reads only files written by milestone2_verify; regenerates every figure without
rerunning the solver.

Usage:
    python plot_milestone2.py <results_dir> [--out <figures_dir>]

Inputs:  convergence.csv, profile_x.csv, error_map_xy.csv
Outputs: milestone02_profile.png
         milestone02_convergence.png
         milestone02_error_map.png
         milestone02_order_table.png
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh_style as S


def read_csv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def fnum(s):
    return float(s) if s not in ("", None) else float("nan")


def plot_profile(rows, out):
    xs = [float(r["x"]) for r in rows]
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 6.4), facecolor=S.SURFACE)
    S.title_block(
        fig,
        "Analytic versus numerical derivatives along the central x line",
        "Trigonometric field, 41x41x41 grid, y = z = 0. Boundary points (open "
        "markers) use one-sided stencils of the same order.",
    )

    for col, (num_k, ana_k, err_k, name) in enumerate([
        ("d1_numeric", "d1_analytic", "d1_error", r"$\partial f/\partial x$"),
        ("d2_numeric", "d2_analytic", "d2_error", r"$\partial^2 f/\partial x^2$"),
    ]):
        ax = axes[0][col]
        S.style_axes(ax)
        ax.plot(xs, [float(r[ana_k]) for r in rows], color=S.SERIES[0],
                linewidth=2.0, label="analytic", zorder=2)
        ax.plot(xs, [float(r[num_k]) for r in rows], color=S.SERIES[1],
                linewidth=0, marker="o", markersize=3.0, label="numerical",
                zorder=3)
        bx = [float(r["x"]) for r in rows if r["is_boundary"] == "1"]
        by = [float(r[num_k]) for r in rows if r["is_boundary"] == "1"]
        ax.plot(bx, by, linewidth=0, marker="o", markersize=8,
                markerfacecolor="none", markeredgecolor=S.STATUS_CRIT,
                markeredgewidth=1.6, label="boundary (one-sided)", zorder=4)
        ax.set_title(name, fontsize=10, color=S.INK, loc="left")
        ax.set_xlabel("x", fontsize=8, color=S.INK2)
        ax.set_ylabel("derivative (dimensionless)", fontsize=8, color=S.INK2)
        ax.legend(fontsize=7, frameon=False, labelcolor=S.INK2)

        ax = axes[1][col]
        S.style_axes(ax)
        ax.axhline(0.0, color=S.INK3, linewidth=0.8, zorder=1)
        ax.plot(xs, [float(r[err_k]) for r in rows], color=S.SERIES[0],
                linewidth=1.6, zorder=2)
        ax.set_title("pointwise error", fontsize=9, color=S.INK, loc="left")
        ax.set_xlabel("x", fontsize=8, color=S.INK2)
        ax.set_ylabel("numerical - analytic", fontsize=8, color=S.INK2)

    S.footer(fig, "Milestone A2 - dimensionless test field, no physical units. "
                  "Grid 41x41x41 on [-1,1]^3, slice y = z = 0.")
    S.save(fig, out, rect=(0, 0.045, 1, 0.90))


def plot_convergence(rows, out):
    groups = defaultdict(list)
    for r in rows:
        groups[(r["field"], r["operator"])].append(r)

    # Only the operators with a real truncation error converge; the rest are
    # exact for that field and sit at roundoff. Plotting them together on a
    # log-log axis would show a rising line and read as a bug, so they are
    # separated and labelled.
    conv, exact = OrderedDict(), OrderedDict()
    for k, v in sorted(groups.items()):
        v.sort(key=lambda r: int(r["n"]))
        (exact if float(v[0]["interior_max"]) < 1e-10 else conv)[k] = v

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.4), facecolor=S.SURFACE)
    S.title_block(
        fig, "Grid convergence of the finite-difference operators",
        f"All {len(conv)} converging operator/field pairs, four resolutions each. "
        "Interior and boundary are shown separately and never averaged.",
        y=0.965,
    )

    # One hue for every curve, deliberately. There are 12 series and identity is
    # not the job of this figure -- the message is that they are all PARALLEL to
    # slope 2. Cycling four categorical hues across twelve series would imply an
    # identity encoding that does not survive a colour-vision check and that no
    # reader could decode anyway. Per-operator detail lives in the order table.
    for ax, key, label in ((axes[0], "interior_max", "interior (centred stencils)"),
                           (axes[1], "boundary_max", "boundary (one-sided stencils)")):
        S.style_axes(ax, logx=True, logy=True)

        anchor = 0.0
        for k, v in conv.items():
            hs = [float(r["h"]) for r in v]
            es = [float(r[key]) for r in v]
            anchor = max(anchor, es[0])
            ax.plot(hs, es, marker="o", markersize=3.2, linewidth=1.3,
                    color=S.SERIES[0], alpha=0.55, zorder=3)

        h0, h1 = 0.02, 0.24
        top = anchor * 2.2
        ax.plot([h0, h1], [top * (h0 / h1) ** 2, top], color=S.INK2,
                linestyle=(0, (5, 3)), linewidth=1.8, zorder=4)
        ax.annotate("slope 2", xy=(h1, top), xytext=(-6, -12),
                    textcoords="offset points", ha="right", fontsize=8,
                    color=S.INK2)

        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(label, fontsize=10, color=S.INK, loc="left")
        ax.set_xlabel("grid spacing h (dimensionless)", fontsize=8, color=S.INK2)
        ax.set_ylabel("max |error|", fontsize=8, color=S.INK2)

    S.footer(fig,
             f"Milestone A2 - {len(exact)} further operator/field pairs are not "
             "plotted: the stencil is exact for those fields, so their error is "
             "roundoff and rises under refinement. See the order table.")
    S.save(fig, out, rect=(0, 0.075, 1, 0.88))


def plot_order_table(rows, out):
    groups = defaultdict(list)
    for r in rows:
        groups[(r["field"], r["operator"])].append(r)

    entries = []
    for k, v in sorted(groups.items()):
        v.sort(key=lambda r: int(r["n"]))
        if float(v[0]["interior_max"]) < 1e-10:
            entries.append((k[0], k[1], None, None, float(v[-1]["interior_max"])))
        else:
            entries.append((k[0], k[1], fnum(v[-1]["observed_order_interior"]),
                            fnum(v[-1]["observed_order_boundary"]), None))

    fig, ax = plt.subplots(figsize=(9.5, 0.34 * len(entries) + 1.7),
                           facecolor=S.SURFACE)
    ax.axis("off")
    S.title_block(fig, "Observed convergence order (finest resolution pair)",
                  "Theoretical order is 2 for every operator, interior and "
                  "boundary alike.", y=0.97)

    hdr = ["field", "operator", "interior order", "boundary order", "status"]
    colx = [0.02, 0.20, 0.40, 0.58, 0.76]
    ytop = 0.86
    for c, h in enumerate(hdr):
        ax.text(colx[c], ytop, h, transform=ax.transAxes, fontsize=8.5,
                color=S.INK, fontweight="bold")
    dy = 0.80 / max(1, len(entries))
    for i, (fld, op, oi, ob, exact_err) in enumerate(entries):
        y = ytop - (i + 1) * dy
        vals = [fld, op,
                "-" if oi is None else f"{oi:.3f}",
                "-" if ob is None else f"{ob:.3f}",
                "exact (roundoff)" if oi is None else "converged"]
        ok = oi is None or (oi >= 1.85 and ob >= 1.85)
        for c, v in enumerate(vals):
            ax.text(colx[c], y, v, transform=ax.transAxes, fontsize=7.5,
                    color=(S.INK2 if ok else S.STATUS_CRIT),
                    family="monospace" if c >= 2 else None)

    S.footer(fig, "Milestone A2 - 'exact' rows are field/operator pairs whose "
                  "truncated derivative vanishes identically; no order is "
                  "defined and the residual error is roundoff.")
    fig.savefig(out, dpi=200, facecolor=S.SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def plot_error_map(rows, out):
    xs = sorted({float(r["x"]) for r in rows})
    ys = sorted({float(r["y"]) for r in rows})
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}
    grid = [[0.0] * len(xs) for _ in ys]
    for r in rows:
        grid[yi[float(r["y"])]][xi[float(r["x"])]] = float(r["error"])

    fig, ax = plt.subplots(figsize=(6.4, 5.0), facecolor=S.SURFACE)
    S.title_block(fig, r"Error in $\partial f/\partial x$, central z slice",
                  "Signed field on a zero-centred diverging scale. The boundary "
                  "rows at x = +-1 carry the one-sided stencil error.", y=0.965)
    S.style_axes(ax)
    S.imshow_field(ax, xs, ys, grid, signed=True,
                   label="numerical - analytic")
    ax.set_xlabel("x", fontsize=8, color=S.INK2)
    ax.set_ylabel("y", fontsize=8, color=S.INK2)
    S.footer(fig, "Milestone A2 - trigonometric test field, 41x41x41 on "
                  "[-1,1]^3, slice at central z. Dimensionless.")
    S.save(fig, out, rect=(0, 0.06, 1, 0.87))


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

    need = ["convergence.csv", "profile_x.csv", "error_map_xy.csv"]
    for n in need:
        if not (rd / n).is_file():
            print(f"error: missing {rd/n}; run milestone2_verify first",
                  file=sys.stderr)
            return 2

    conv = read_csv(rd / "convergence.csv")
    plot_profile(read_csv(rd / "profile_x.csv"),
                 out / "milestone02_profile.png")
    plot_convergence(conv, out / "milestone02_convergence.png")
    plot_order_table(conv, out / "milestone02_order_table.png")
    plot_error_map(read_csv(rd / "error_map_xy.csv"),
                   out / "milestone02_error_map.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
