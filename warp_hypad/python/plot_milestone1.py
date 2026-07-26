#!/usr/bin/env python3
"""Milestone 1 verification figures.

Reads only files written by the C++ solver and regenerates every Milestone 1
figure without rerunning it. That is a hard requirement of the project: a figure
that cannot be rebuilt from a results directory is not reproducible.

Usage:
    python plot_milestone1.py <results_dir> [--out <figures_dir>]

Inputs (all produced by milestone1_verify):
    <results_dir>/inverse_residuals.csv
    <results_dir>/verification_summary.json

Outputs:
    <figures_dir>/milestone01_inverse_residual_histogram.png
    <figures_dir>/milestone01_condition_vs_residual.png
    <figures_dir>/milestone01_memory_layout.png

Design notes
------------
Both data figures are faceted small multiples rather than one panel with four
coloured series. The families differ by two orders of magnitude in conditioning,
so overplotting them would hide exactly the structure the figure exists to show;
and faceting means a single validated hue carries the whole figure, so no
categorical pair has to survive a colour-vision check.

Milestone 1 quantities are dimensionless algebraic diagnostics. There is no
metric, no grid and no slice, so those annotations are stated as not applicable
rather than silently omitted.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import OrderedDict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# --- Palette (validated reference instance, light mode) ---------------------
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8983"
SERIES_1 = "#2a78d6"
STATUS_WARNING = "#fab219"
STATUS_CRITICAL = "#d03b3b"
GRID_COLOR = "#e4e3df"

# Residuals of exactly zero cannot be drawn on a log axis. They are real results,
# so they are clamped to this floor and the clamping is stated on the figure
# rather than dropped silently.
LOG_FLOOR = 1e-18

FAMILY_TITLES = OrderedDict(
    [
        ("diagonal", "Diagonal, mixed signs"),
        ("symmetric", "General symmetric"),
        ("lorentzian", "Perturbed Lorentzian ($-{+}{+}{+}$)"),
        ("general", "General nonsingular"),
    ]
)

FOOTER = (
    "Milestone 1 — dimensionless algebraic diagnostics. "
    "Units: dimensionless (no physical quantity; geometrized units not yet in play). "
    "Metric: n/a · Grid: n/a (pointwise 4×4) · Slice: n/a"
)


def style_axes(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, which="major", color=GRID_COLOR, linewidth=0.8, zorder=0)
    ax.grid(True, which="minor", color=GRID_COLOR, linewidth=0.4, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(INK_MUTED)
        ax.spines[spine].set_linewidth(0.8)
    ax.tick_params(colors=INK_SECONDARY, labelsize=8, width=0.8)


def add_footer(fig, extra=""):
    text = FOOTER if not extra else FOOTER + "\n" + extra
    fig.text(
        0.01,
        0.01,
        text,
        ha="left",
        va="bottom",
        fontsize=6.5,
        color=INK_MUTED,
        linespacing=1.5,
    )


def read_samples(path: Path):
    families = OrderedDict((k, []) for k in FAMILY_TITLES)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            fam = row["family"]
            if fam not in families:
                families[fam] = []
            families[fam].append(
                {
                    "cond1": float(row["cond1"]),
                    "residual_left": float(row["residual_left"]),
                    "residual_right": float(row["residual_right"]),
                    "det_rel_error": float(row["det_rel_error"]),
                }
            )
    return families


def clamp(values):
    """Clamp to the log floor, and report how many were clamped."""
    clamped = [max(v, LOG_FLOOR) for v in values]
    n_zero = sum(1 for v in values if v <= 0.0)
    return clamped, n_zero


def plot_residual_histogram(families, out_path: Path, tolerance: float):
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.4), facecolor=SURFACE)
    fig.suptitle(
        "Matrix-inverse identity residual, $\\max|AA^{-1}-I|$",
        fontsize=13,
        color=INK_PRIMARY,
        x=0.01,
        ha="left",
        y=0.975,
    )
    fig.text(
        0.01,
        0.925,
        "500 random matrices per family, adjugate inverse, double precision. "
        "Lower is better; the acceptance tolerance is marked.",
        fontsize=8.5,
        color=INK_SECONDARY,
        ha="left",
    )

    total_zero = 0
    for ax, (family, title) in zip(axes.flat, FAMILY_TITLES.items()):
        style_axes(ax)
        samples = families.get(family, [])
        values = [s["residual_left"] for s in samples]
        if not values:
            ax.text(0.5, 0.5, "no samples", transform=ax.transAxes,
                    ha="center", color=INK_MUTED, fontsize=9)
            ax.set_title(title, fontsize=10, color=INK_PRIMARY, loc="left")
            continue

        drawn, n_zero = clamp(values)
        total_zero += n_zero

        lo = math.floor(math.log10(min(drawn)))
        hi = math.ceil(math.log10(max(drawn)))
        if hi <= lo:
            hi = lo + 1
        bins = [10.0 ** (lo + i * (hi - lo) / 40.0) for i in range(41)]

        ax.hist(drawn, bins=bins, color=SERIES_1, edgecolor=SURFACE, linewidth=0.5, zorder=2)
        ax.set_xscale("log")

        ax.axvline(tolerance, color=STATUS_CRITICAL, linewidth=2.0, linestyle=(0, (4, 2)), zorder=3)
        ax.annotate(
            f"tolerance {tolerance:g}",
            xy=(tolerance, ax.get_ylim()[1]),
            xytext=(-4, -6),
            textcoords="offset points",
            ha="right",
            va="top",
            fontsize=7.5,
            color=STATUS_CRITICAL,
        )

        worst = max(values)
        ax.set_title(title, fontsize=10, color=INK_PRIMARY, loc="left")
        ax.text(
            0.98,
            0.94,
            f"n = {len(values)}\nworst = {worst:.2e}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7.5,
            color=INK_SECONDARY,
        )
        ax.set_xlabel("identity residual $\\max_{\\mu\\nu}|(AA^{-1}-I)_{\\mu\\nu}|$ (dimensionless)",
                      fontsize=8, color=INK_SECONDARY)
        ax.set_ylabel("count of matrices", fontsize=8, color=INK_SECONDARY)

    extra = ""
    if total_zero:
        extra = (f"{total_zero} residuals were exactly 0 and are drawn at the "
                 f"log-axis floor {LOG_FLOOR:g}.")
    add_footer(fig, extra)
    fig.tight_layout(rect=(0, 0.055, 1, 0.90))
    fig.savefig(out_path, dpi=200, facecolor=SURFACE)
    plt.close(fig)


def plot_condition_vs_residual(families, out_path: Path, warn_limit: float):
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.4), facecolor=SURFACE)
    fig.suptitle(
        "Inverse residual versus 1-norm condition number",
        fontsize=13,
        color=INK_PRIMARY,
        x=0.01,
        ha="left",
        y=0.975,
    )
    fig.text(
        0.01,
        0.925,
        "The guide line is the textbook bound residual $\\approx \\kappa_1 \\varepsilon$. "
        "Points tracking it confirm the loss of accuracy is conditioning, not a bug.",
        fontsize=8.5,
        color=INK_SECONDARY,
        ha="left",
    )

    eps = 2.220446049250313e-16
    total_zero = 0

    for ax, (family, title) in zip(axes.flat, FAMILY_TITLES.items()):
        style_axes(ax)
        samples = families.get(family, [])
        if not samples:
            ax.text(0.5, 0.5, "no samples", transform=ax.transAxes,
                    ha="center", color=INK_MUTED, fontsize=9)
            ax.set_title(title, fontsize=10, color=INK_PRIMARY, loc="left")
            continue

        conds = [max(s["cond1"], 1.0) for s in samples]
        raw_res = [s["residual_left"] for s in samples]
        res, n_zero = clamp(raw_res)
        total_zero += n_zero

        ax.scatter(conds, res, s=14, color=SERIES_1, alpha=0.55,
                   edgecolors=SURFACE, linewidths=0.4, zorder=2)

        c_lo, c_hi = min(conds), max(conds)
        if c_hi <= c_lo:
            c_hi = c_lo * 10.0
        guide_x = [c_lo, c_hi]
        guide_y = [c * eps for c in guide_x]
        ax.plot(guide_x, guide_y, color=INK_SECONDARY, linewidth=1.5,
                linestyle=(0, (5, 3)), zorder=3)
        ax.annotate(
            "$\\kappa_1\\,\\varepsilon$",
            xy=(guide_x[-1], guide_y[-1]),
            xytext=(-6, 6),
            textcoords="offset points",
            ha="right",
            fontsize=8,
            color=INK_SECONDARY,
        )

        if c_lo <= warn_limit <= c_hi:
            ax.axvline(warn_limit, color=STATUS_WARNING, linewidth=2.0,
                       linestyle=(0, (4, 2)), zorder=3)
            ax.annotate(
                "ill-conditioned\nflag fires",
                xy=(warn_limit, min(res)),
                xytext=(5, 2),
                textcoords="offset points",
                ha="left",
                va="bottom",
                fontsize=7.5,
                color="#8a6410",
            )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(title, fontsize=10, color=INK_PRIMARY, loc="left")
        ax.set_xlabel("condition number $\\kappa_1(A)=\\|A\\|_1\\|A^{-1}\\|_1$ (dimensionless)",
                      fontsize=8, color=INK_SECONDARY)
        ax.set_ylabel("identity residual (dimensionless)", fontsize=8, color=INK_SECONDARY)
        ax.text(
            0.02,
            0.94,
            f"n = {len(samples)}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.5,
            color=INK_SECONDARY,
        )

    extra = (f"Ill-conditioning flag fires at κ₁ > {warn_limit:.3g} "
             f"(= 1e-8 / ε), chosen so the flag means the residual may exceed 1e-8.")
    if total_zero:
        extra += (f" {total_zero} residuals were exactly 0 and are drawn at the "
                  f"log-axis floor {LOG_FLOOR:g}.")
    add_footer(fig, extra)
    fig.tight_layout(rect=(0, 0.065, 1, 0.90))
    fig.savefig(out_path, dpi=200, facecolor=SURFACE)
    plt.close(fig)


def plot_memory_layout(out_path: Path):
    """The memory-layout diagram required by the Milestone 1 gate.

    Shows why point-major was chosen: the 16 components of one grid point are
    contiguous, so the pointwise 4x4 algebra that dominates the pipeline reads
    sequentially, while a finite-difference stencil pays a regular stride.
    """
    fig, ax = plt.subplots(figsize=(9.6, 3.5), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.axis("off")

    fig.suptitle(
        "Rank-2 tensor field storage — point-major layout",
        fontsize=13,
        color=INK_PRIMARY,
        x=0.01,
        ha="left",
        y=0.98,
    )
    fig.text(
        0.01,
        0.885,
        "raw()[16p + 4μ + ν]   where   p = (i·n$_y$ + j)·n$_z$ + k",
        fontsize=10,
        color=INK_SECONDARY,
        ha="left",
        family="monospace",
    )

    cell_w, cell_h = 0.052, 0.26
    x0, y0 = 0.045, 0.50

    # Two grid points, 16 components each.
    for block in range(2):
        for c in range(16):
            idx = block * 16 + c
            x = x0 + idx * cell_w + block * 0.03
            shade = SERIES_1 if block == 0 else INK_MUTED
            alpha = 0.18 + 0.05 * (c // 4)
            ax.add_patch(
                Rectangle((x, y0), cell_w * 0.94, cell_h,
                          facecolor=shade, alpha=alpha,
                          edgecolor=shade, linewidth=1.0,
                          transform=ax.transAxes, clip_on=False)
            )
            ax.text(x + cell_w * 0.47, y0 + cell_h * 0.58,
                    f"{c // 4}{c % 4}",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=6.8, color=INK_PRIMARY, family="monospace")
            ax.text(x + cell_w * 0.47, y0 + cell_h * 0.18,
                    f"{idx}",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=5.6, color=INK_MUTED, family="monospace")

    for block, label in enumerate(["grid point p = 0", "grid point p = 1"]):
        span_x = x0 + block * (16 * cell_w + 0.03)
        ax.text(span_x + 8 * cell_w, y0 + cell_h + 0.055, label,
                transform=ax.transAxes, ha="center", fontsize=9,
                color=INK_PRIMARY)
        ax.annotate(
            "", xy=(span_x, y0 + cell_h + 0.03),
            xytext=(span_x + 16 * cell_w * 0.94, y0 + cell_h + 0.03),
            xycoords="axes fraction", textcoords="axes fraction",
            arrowprops=dict(arrowstyle="<->", color=INK_MUTED, linewidth=0.9),
        )

    ax.text(x0, y0 - 0.09,
            "cell label = (μν) component index      "
            "small number = offset into raw()",
            transform=ax.transAxes, fontsize=7.5, color=INK_MUTED)

    ax.text(
        x0, y0 - 0.42,
        "Chosen because pointwise 4×4 algebra (inverse, contraction) dominates the pipeline and reads one\n"
        "point contiguously. A finite-difference stencil reads the same component at neighbouring points, a\n"
        "regular stride of 16·(grid stride). Revisited with measurements in Milestone 9, not before.",
        transform=ax.transAxes, fontsize=8, color=INK_SECONDARY, linespacing=1.6,
    )

    add_footer(fig)
    fig.savefig(out_path, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("--out", type=Path, default=None,
                        help="figures directory (default: <results_dir>/figures)")
    args = parser.parse_args()

    results_dir: Path = args.results_dir
    if not results_dir.is_dir():
        print(f"error: not a directory: {results_dir}", file=sys.stderr)
        return 2

    residual_csv = results_dir / "inverse_residuals.csv"
    if not residual_csv.is_file():
        print(f"error: missing {residual_csv}; run milestone1_verify first",
              file=sys.stderr)
        return 2

    out_dir: Path = args.out or (results_dir / "figures")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Read the tolerance actually used by the gate rather than hard-coding it,
    # so the figure cannot disagree with the verification record.
    tolerance = 1e-9
    warn_limit = 1e-8 / 2.220446049250313e-16
    summary_path = results_dir / "verification_summary.json"
    if summary_path.is_file():
        with summary_path.open() as handle:
            summary = json.load(handle)
        for case in summary.get("cases", []):
            if case.get("name", "").startswith("algebra.inverse_identity_left"):
                tolerance = case.get("tolerance", tolerance)
                break

    families = read_samples(residual_csv)
    n = sum(len(v) for v in families.values())
    if n == 0:
        print(f"error: no samples in {residual_csv}", file=sys.stderr)
        return 2

    hist_path = out_dir / "milestone01_inverse_residual_histogram.png"
    scatter_path = out_dir / "milestone01_condition_vs_residual.png"
    layout_path = out_dir / "milestone01_memory_layout.png"

    plot_residual_histogram(families, hist_path, tolerance)
    plot_condition_vs_residual(families, scatter_path, warn_limit)
    plot_memory_layout(layout_path)

    print(f"read {n} samples from {residual_csv}")
    for path in (hist_path, scatter_path, layout_path):
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
