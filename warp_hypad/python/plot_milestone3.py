#!/usr/bin/env python3
"""Milestone A3 figures: metric generation and HyPAD sensitivities.

Reads only files written by milestone3_verify; regenerates every figure without
rerunning the solver.

Usage:
    python plot_milestone3.py <results_dir> [--out <figures_dir>]

Inputs (sensitivity inputs are optional -- present only when the OTIlib backend
was enabled):
    shape_profiles.csv, metric_slice_xy.csv, centerline_x.csv
    sensitivity_slice_xy.csv, hypad_vs_analytic.csv, hypad_vs_fd.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh_style as S

FROZEN = ("FROZEN-TIME SPATIAL SLICE: the bubble centre is fixed and "
          "d_t g = 0 by construction. NOT the complete Alcubierre spacetime.")
PARAMS = "Metric: Alcubierre  |  v = 0.2, R = 3, sigma = 2  |  geometrized units (G = c = 1)"


def read_csv(p):
    with p.open(newline="") as fh:
        return list(csv.DictReader(fh))


def to_grid(rows, xk, yk, vk):
    xs = sorted({float(r[xk]) for r in rows})
    ys = sorted({float(r[yk]) for r in rows})
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}
    g = [[0.0] * len(xs) for _ in ys]
    for r in rows:
        g[yi[float(r[yk])]][xi[float(r[xk])]] = float(r[vk])
    return xs, ys, g


def plot_shape(rows, out):
    by_R = defaultdict(lambda: defaultdict(list))
    for r in rows:
        by_R[float(r["R"])][float(r["sigma"])].append((float(r["r"]), float(r["f"])))

    Rs = sorted(by_R)
    fig, axes = plt.subplots(1, len(Rs), figsize=(9.6, 3.6), facecolor=S.SURFACE,
                             sharey=True)
    S.title_block(fig, "Alcubierre shape function f(r)",
                  "f(0) = 1 exactly and f -> 0 as r grows, so the interior is a "
                  "flat boosted region and the exterior is asymptotically flat. "
                  "Larger sigma steepens the wall; larger R moves it outward.",
                  y=0.965)

    for ax, R in zip(axes, Rs):
        S.style_axes(ax)
        for i, sg in enumerate(sorted(by_R[R])):
            pts = sorted(by_R[R][sg])
            ax.plot([p[0] for p in pts], [p[1] for p in pts],
                    color=S.SERIES[i], linewidth=2.0,
                    label=f"$\\sigma$ = {sg:g}", zorder=3)
        ax.axvline(R, color=S.INK3, linestyle=(0, (4, 3)), linewidth=1.2, zorder=2)
        ax.annotate(f"r = R = {R:g}", xy=(R, 1.02), fontsize=7, color=S.INK2,
                    ha="left", va="bottom", rotation=90)
        ax.set_title(f"R = {R:g}", fontsize=10, color=S.INK, loc="left")
        ax.set_xlabel("r (grid length unit)", fontsize=8, color=S.INK2)
        ax.set_ylim(-0.05, 1.15)
        ax.legend(fontsize=7, frameon=False, labelcolor=S.INK2)
    axes[0].set_ylabel("f (dimensionless)", fontsize=8, color=S.INK2)

    S.footer(fig, "Milestone A3  |  " + PARAMS + "\n" + FROZEN)
    S.save(fig, out, rect=(0, 0.11, 1, 0.86))


def plot_metric_slices(rows, out):
    panels = [
        ("g00", r"$g_{00} = -1 + v^2 f^2$", True),
        ("g01", r"$g_{0x} = -v f$", True),
        ("g11", r"$g_{xx} = 1$", False),
        ("det", r"$\det g$  (exactly $-1$ everywhere)", False),
        ("inverse_residual", r"$\max|g_{\mu\alpha}g^{\alpha\nu} - \delta|$", False),
        ("signature_ok", "Lorentzian signature mask (1 = valid)", False),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.6), facecolor=S.SURFACE)
    S.title_block(fig, "Metric components and validity, central z slice",
                  "Signed components use a zero-centred diverging scale. The "
                  "determinant and signature panels are the pointwise validity "
                  "checks, not derived quantities.", y=0.97)

    for ax, (key, label, signed) in zip(axes.flat, panels):
        S.style_axes(ax)
        xs, ys, g = to_grid(rows, "x", "y", key)
        flat = [v for row in g for v in row]
        if max(flat) - min(flat) < 1e-14:
            # A constant panel would give a degenerate colour bar; say the value.
            ax.text(0.5, 0.5, f"constant\n{flat[0]:.3g}", transform=ax.transAxes,
                    ha="center", va="center", fontsize=11, color=S.INK2)
            ax.set_xticks([]); ax.set_yticks([])
        else:
            S.imshow_field(ax, xs, ys, g, signed=signed, label="")
        ax.set_title(label, fontsize=9.5, color=S.INK, loc="left")
        ax.set_xlabel("x", fontsize=7.5, color=S.INK2)
        ax.set_ylabel("y", fontsize=7.5, color=S.INK2)

    S.footer(fig, "Milestone A3  |  " + PARAMS + "  |  81x81 slice at central z\n" + FROZEN)
    S.save(fig, out, rect=(0, 0.07, 1, 0.90))


def plot_sensitivity(rows, out):
    panels = [
        ("dg00_dv", r"$\partial g_{00}/\partial v$"),
        ("dg00_dR", r"$\partial g_{00}/\partial R$"),
        ("dg00_dsigma", r"$\partial g_{00}/\partial \sigma$"),
        ("dg0x_dv", r"$\partial g_{0x}/\partial v$"),
        ("dg0x_dR", r"$\partial g_{0x}/\partial R$"),
        ("dg0x_dsigma", r"$\partial g_{0x}/\partial \sigma$"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.6), facecolor=S.SURFACE)
    S.title_block(fig, "HyPAD parameter-sensitivity maps, central z slice",
                  "All six fields come from ONE metric evaluation with OTIScalar "
                  "seeded in three directions. Panels that change sign use a "
                  "zero-centred diverging scale; single-signed panels use a "
                  "sequential ramp with zero at the light end.", y=0.97)

    for ax, (key, label) in zip(axes.flat, panels):
        S.style_axes(ax)
        xs, ys, g = to_grid(rows, "x", "y", key)
        S.imshow_field(ax, xs, ys, g, signed=True, label="")
        ax.set_title(label, fontsize=10, color=S.INK, loc="left")
        ax.set_xlabel("x", fontsize=7.5, color=S.INK2)
        ax.set_ylabel("y", fontsize=7.5, color=S.INK2)

    S.footer(fig, "Milestone A3  |  " + PARAMS +
             "  |  derivative channels: velocity, radius, sigma\n" + FROZEN)
    S.save(fig, out, rect=(0, 0.07, 1, 0.90))


def plot_hypad_vs_analytic(rows, out):
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), facecolor=S.SURFACE)
    S.title_block(fig, "HyPAD versus analytic metric sensitivities",
                  "Closed forms derived by hand from the definition of f, not by "
                  "differentiating the implementation -- so agreement means "
                  "something.", y=0.965)

    ana = [abs(float(r["analytic"])) for r in rows]
    err = [abs(float(r["hypad"]) - float(r["analytic"])) for r in rows]

    ax = axes[0]
    S.style_axes(ax, logx=True, logy=True)
    floor = 1e-20
    ax.scatter([max(a, floor) for a in ana], [max(e, floor) for e in err],
               s=26, color=S.SERIES[0], alpha=0.7, edgecolors=S.SURFACE,
               linewidths=0.5, zorder=3)
    lo, hi = 1e-18, 1e2
    ax.plot([lo, hi], [lo * 2.2e-16, hi * 2.2e-16], color=S.INK2,
            linestyle=(0, (5, 3)), linewidth=1.5, zorder=2)
    ax.annotate(r"$\varepsilon\,|{\rm analytic}|$", xy=(hi, hi * 2.2e-16),
                xytext=(-6, 6), textcoords="offset points", ha="right",
                fontsize=8, color=S.INK2)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1e-18, 1e2)
    ax.set_title("absolute error versus derivative magnitude", fontsize=9.5,
                 color=S.INK, loc="left")
    ax.set_xlabel("|analytic derivative|", fontsize=8, color=S.INK2)
    ax.set_ylabel("|HyPAD - analytic|", fontsize=8, color=S.INK2)

    ax = axes[1]
    S.style_axes(ax)
    order = sorted(range(len(rows)), key=lambda i: -err[i])[:14]
    labels = [f"{rows[i]['component']} d/d{rows[i]['parameter'][:3]} "
              f"@{rows[i]['location'][:12]}" for i in order]
    ax.barh(range(len(order)), [err[i] for i in order], color=S.SERIES[0],
            height=0.68, zorder=3)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(labels, fontsize=6.5, color=S.INK2)
    ax.invert_yaxis()
    ax.set_xlabel("absolute error", fontsize=8, color=S.INK2)
    ax.set_title("largest 14 discrepancies", fontsize=9.5, color=S.INK, loc="left")

    S.footer(fig, "Milestone A3  |  " + PARAMS +
             f"  |  {len(rows)} comparisons across 6 sample points\n" + FROZEN)
    S.save(fig, out, rect=(0, 0.10, 1, 0.88))


def plot_hypad_vs_fd(rows, out):
    by = defaultdict(list)
    for r in rows:
        by[r["parameter"]].append(r)

    params = sorted(by)
    fig, axes = plt.subplots(1, len(params), figsize=(10.5, 3.8),
                            facecolor=S.SURFACE, sharey=True)
    S.title_block(fig, "HyPAD versus centered finite differences",
                  "Error falls as O(h^2) then rises as O(eps/h): the V-shape is "
                  "the evidence that FD is a genuinely independent method. FD "
                  "cannot do better; HyPAD has no step size at all.", y=0.965)

    for ax, pname in zip(axes, params):
        S.style_axes(ax, logx=True, logy=True)
        locs = OrderedDict()
        for r in by[pname]:
            locs.setdefault(r["location"], []).append(r)
        for i, (loc, rs) in enumerate(list(locs.items())[:3]):
            rs.sort(key=lambda r: float(r["step_size"]))
            ax.plot([float(r["step_size"]) for r in rs],
                    [max(float(r["relative_error"]), 1e-18) for r in rs],
                    marker="o", markersize=3.4, linewidth=1.5,
                    color=S.SERIES[i], label=loc, zorder=3)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(f"d/d{pname}", fontsize=10, color=S.INK, loc="left")
        ax.set_xlabel("FD step size h", fontsize=8, color=S.INK2)
        ax.legend(fontsize=6.5, frameon=False, labelcolor=S.INK2)
    axes[0].set_ylabel("scaled error vs HyPAD", fontsize=8, color=S.INK2)

    S.footer(fig, "Milestone A3  |  " + PARAMS +
             "  |  error scaled by max(|HyPAD|, 1e-9) so it stays meaningful "
             "where the derivative is negligible\n" + FROZEN)
    S.save(fig, out, rect=(0, 0.12, 1, 0.87))


def plot_centerline(rows, out):
    xs = [float(r["x"]) for r in rows]
    fig, ax = plt.subplots(figsize=(7.6, 4.0), facecolor=S.SURFACE)
    S.title_block(fig, "Centreline profile through the bubble (y = z = 0)",
                  "Shape function and the two nontrivial metric components.",
                  y=0.955)
    S.style_axes(ax)
    for i, (k, lbl) in enumerate([("f", "f(r)"), ("g00", r"$g_{00}$"),
                                  ("g01", r"$g_{0x}$")]):
        ax.plot(xs, [float(r[k]) for r in rows], color=S.SERIES[i],
                linewidth=2.0, label=lbl, zorder=3)
    ax.axhline(0.0, color=S.INK3, linewidth=0.8, zorder=1)
    ax.set_xlabel("x (grid length unit)", fontsize=8, color=S.INK2)
    ax.set_ylabel("value (dimensionless)", fontsize=8, color=S.INK2)
    ax.legend(fontsize=8, frameon=False, labelcolor=S.INK2)
    S.footer(fig, "Milestone A3  |  " + PARAMS + "\n" + FROZEN)
    S.save(fig, out, rect=(0, 0.12, 1, 0.88))


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

    for name in ("shape_profiles.csv", "metric_slice_xy.csv", "centerline_x.csv"):
        if not (rd / name).is_file():
            print(f"error: missing {rd/name}; run milestone3_verify first",
                  file=sys.stderr)
            return 2

    plot_shape(read_csv(rd / "shape_profiles.csv"),
               out / "milestone03_shape_profiles.png")
    plot_metric_slices(read_csv(rd / "metric_slice_xy.csv"),
                       out / "milestone03_metric_slices.png")
    plot_centerline(read_csv(rd / "centerline_x.csv"),
                    out / "milestone03_centerline.png")

    # Sensitivity outputs exist only when the OTIlib backend was enabled.
    opt = {
        "sensitivity_slice_xy.csv": (plot_sensitivity,
                                     "milestone03_sensitivity_maps.png"),
        "hypad_vs_analytic.csv": (plot_hypad_vs_analytic,
                                  "milestone03_hypad_vs_analytic.png"),
        "hypad_vs_fd.csv": (plot_hypad_vs_fd, "milestone03_hypad_vs_fd.png"),
    }
    missing = []
    for name, (fn, png) in opt.items():
        p = rd / name
        if p.is_file():
            fn(read_csv(p), out / png)
        else:
            missing.append(name)
    if missing:
        print("skipped sensitivity figures (OTIlib backend disabled): "
              + ", ".join(missing))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
