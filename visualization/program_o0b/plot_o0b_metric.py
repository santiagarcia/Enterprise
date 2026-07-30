#!/usr/bin/env python3
"""The twelve O0B ADM-metric figures, drawn from saved C++ records only.

    python plot_o0b_metric.py --records <o0b_metric dir> --out <gallery dir>
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402
from matplotlib.lines import Line2D                              # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

import o0b_common as O                                           # noqa: E402
from common import anomaly, records                              # noqa: E402
from themes import themes                                        # noqa: E402

BLUE, AMBER, RED, GREEN = "#4ea1ff", "#ffb454", "#d1495b", "#7ed957"


def fig_adm_components(rec, out, man, theme):
    fr = O.read_csv(os.path.join(rec, "frames_amplitude_sweep.csv"))
    last = [r for r in fr if r["frame"] == "12"]
    xs, ys, g00 = O.grid_from(last, "g00")
    _, _, g01 = O.grid_from(last, "g01")
    _, _, det = O.grid_from(last, "det_tx")
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.9),
                             facecolor=theme.background)
    for ax, gg, lab, cm in ((axes[0], g00, r"$g_{00}=-\alpha^2+\gamma_{ij}"
                                          r"\beta^i\beta^j$", "RdBu_r"),
                            (axes[1], g01, r"$g_{01}=\gamma_{1j}\beta^j$",
                             "RdBu_r"),
                            (axes[2], det, r"$\det$ of the $(t,x)$ block",
                             "RdBu_r")):
        m = max(abs(v) for row in gg for v in row) or 1.0
        im = ax.imshow(gg, origin="lower",
                       extent=[min(xs), max(xs), min(ys), max(ys)],
                       aspect="auto", cmap=cm, vmin=-m, vmax=m,
                       interpolation="nearest")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=6.5, colors=theme.foreground)
        cb.outline.set_visible(False)
        if gg is g00:
            flat = [v for row in gg for v in row]
            if min(flat) < 0.0 < max(flat):
                ax.contour(xs, ys, gg, levels=[0.0], colors=[AMBER],
                           linewidths=1.6)
        ax.set_title(lab, fontsize=9.5, color=theme.foreground)
        O.style_axes(ax, theme, "x (m)", "transverse radius (m)")
    fig.suptitle("O0B  ADM metric components of the F1 shift-only family",
                 fontsize=12, color=theme.foreground, x=0.010, ha="left",
                 y=0.985)
    O.banner(fig, theme,
             "alpha = 1 and gamma = I here, so g_00 = -1 + (beta^x)^2 and "
             "g_01 = beta^x. The amber contour is g_00 = 0. " + O.CROSSING_NOTE)
    O.footer(fig, theme, "The (t,x) determinant is exactly -1 everywhere, which "
                         "is what makes the previous sentence checkable. " +
             O.SCOPE)
    fig.tight_layout(rect=(0, 0.055, 1, 0.885))
    return O.save_figure(
        fig, out, "adm_metric_components", manifest=man,
        inputs=[os.path.join(rec, "frames_amplitude_sweep.csv")], theme=theme,
        source_rows="frame 12 of frames_amplitude_sweep.csv")


def fig_baseline_vs_zero(rec, out, man, theme):
    rows = O.read_csv(os.path.join(rec, "alcubierre_baseline_comparison.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.0),
                             facecolor=theme.background)
    ax = axes[0]
    ax.plot(O.col(rows, "g00_alcubierre"), O.col(rows, "g00_o0b"), marker="o",
            markersize=4, linestyle="none", color=BLUE, label=r"$g_{00}$")
    ax.plot(O.col(rows, "g01_alcubierre"), O.col(rows, "g01_o0b"), marker="s",
            markersize=4, linestyle="none", color=AMBER, label=r"$g_{01}$")
    lim = max(max(abs(v) for v in O.col(rows, "g00_alcubierre")),
              max(abs(v) for v in O.col(rows, "g01_alcubierre")), 1.0)
    ax.plot([-lim, lim], [-lim, lim], color=theme.grid, lw=0.9, ls=":",
            label="exact agreement")
    ax.legend(fontsize=7.2, frameon=False, labelcolor=theme.foreground)
    ax.set_title("O0B ADM path against production Alcubierre", fontsize=9.5,
                 color=theme.foreground)
    O.style_axes(ax, theme, "production value", "O0B value")

    ax = axes[1]
    d = O.col(rows, "max_component_difference")
    ax.semilogy(range(len(d)), [max(v, 1e-20) for v in d], marker="o",
                markersize=3.4, color=RED, linestyle="none")
    ax.axhline(1e-15, color=AMBER, lw=1.2, ls="--", label="1e-15")
    ax.legend(fontsize=7.2, frameon=False, labelcolor=theme.foreground)
    ax.set_title("worst component difference, per sample point", fontsize=9.5,
                 color=theme.foreground)
    O.style_axes(ax, theme, "sample", "|difference|")

    fig.suptitle("O0B  zero correction reproduces the production metric",
                 fontsize=12, color=theme.foreground, x=0.010, ha="left",
                 y=0.985)
    O.banner(fig, theme,
             "Worst difference over %d point-times is %.2e. The sign "
             "beta^x = -v_s f is FIXED by production, not chosen: "
             "AlcubierreMetric writes g_01 = -v f, and the ADM form "
             "g_0i = gamma_ij beta^j with gamma = I forces it."
             % (len(rows), max(d)))
    O.footer(fig, theme, O.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.885))
    hdr = ["t", "x", "y", "z", "g00_o0b", "g00_alcubierre", "g01_o0b",
           "g01_alcubierre", "max_component_difference"]
    return O.save_figure(
        fig, out, "baseline_vs_zero_correction", manifest=man,
        inputs=[os.path.join(rec, "alcubierre_baseline_comparison.csv")],
        theme=theme, extra={"worst_difference": max(d)},
        plotted_header=hdr, plotted_rows=[[r[k] for k in hdr] for r in rows])


def fig_single_mode(rec, out, man, theme):
    rows = O.read_csv(os.path.join(rec, "frames_single_shift_mode.csv"))
    by = O.group_by(rows, "mode", int)
    k = sorted(by)[len(by) // 2]
    sub = by[k]
    xs, ys, dg = O.grid_from(sub, "beta_delta")
    _, _, bx = O.grid_from(sub, "beta_x")
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.0),
                             facecolor=theme.background)
    m = max(abs(v) for row in dg for v in row) or 1.0
    for ax, gg, lab, lim in ((axes[0], bx, r"$\beta^x$ with the mode active",
                              None),
                             (axes[1], dg,
                              r"$\Delta\beta^x$: the mode's own contribution",
                              m)):
        mm = lim or (max(abs(v) for row in gg for v in row) or 1.0)
        im = ax.imshow(gg, origin="lower",
                       extent=[min(xs), max(xs), min(ys), max(ys)],
                       aspect="auto", cmap="RdBu_r", vmin=-mm, vmax=mm,
                       interpolation="nearest")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=6.5, colors=theme.foreground)
        cb.outline.set_visible(False)
        ax.set_title(lab, fontsize=9.5, color=theme.foreground)
        O.style_axes(ax, theme, "x (m)", "transverse radius (m)")
    sm = O.read_csv(os.path.join(rec, "single_mode_metric_checks.csv"))[0]
    fig.suptitle("O0B  a single local correction mode in the metric",
                 fontsize=12, color=theme.foreground, x=0.010, ha="left",
                 y=0.985)
    O.banner(fig, theme,
             "Mode %d. The change is confined to the mode's own support, which "
             "is why the collar stays exactly at the baseline. The analytic jet "
             "of this metric agrees with centred differences of the COMPLETE "
             "metric to %s." % (k, sm["max_first_derivative_residual"]))
    O.footer(fig, theme, O.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.885))
    return O.save_figure(
        fig, out, "single_mode_metric_perturbation", manifest=man,
        inputs=[os.path.join(rec, "frames_single_shift_mode.csv"),
                os.path.join(rec, "single_mode_metric_checks.csv")],
        theme=theme, extra={"mode_shown": k},
        source_rows="rows of frames_single_shift_mode.csv with mode = %d" % k)


def fig_inverse_residual(rec, out, man, theme):
    rows = O.read_csv(os.path.join(rec, "inverse_metric_samples.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 4.9),
                             facecolor=theme.background)
    ax = axes[0]
    r = [max(float(x["inverse_residual"]), 1e-20) for x in rows]
    ax.semilogy(range(len(r)), r, marker="o", markersize=5, color=RED,
                linestyle="none")
    ax.axhline(1e-14, color=AMBER, lw=1.2, ls="--", label="gate tolerance 1e-14")
    ax.legend(fontsize=7.2, frameon=False, labelcolor=theme.foreground)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([x["label"] for x in rows], fontsize=6.2, rotation=30,
                       ha="right", color=theme.foreground)
    ax.set_title(r"$\max|g_{\mu\lambda}g^{\lambda\nu}-\delta|$", fontsize=9.5,
                 color=theme.foreground)
    O.style_axes(ax, theme, "", "residual")

    ax = axes[1]
    ax.plot(O.col(rows, "gi00"), marker="o", markersize=4, color=BLUE,
            label=r"$g^{00}=-1/\alpha^2$")
    ax.plot(O.col(rows, "gi01"), marker="s", markersize=4, color=AMBER,
            label=r"$g^{01}=\beta^x/\alpha^2$")
    ax.legend(fontsize=7.2, frameon=False, labelcolor=theme.foreground)
    ax.set_title("analytic inverse components", fontsize=9.5,
                 color=theme.foreground)
    O.style_axes(ax, theme, "sample", "value")

    fig.suptitle("O0B  the analytic inverse and its residual", fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    O.banner(fig, theme,
             "The inverse is assembled in CLOSED FORM from the ADM variables, "
             "not by inverting the 4x4. That is what makes this residual an "
             "independent check: a numerical inverse of the same matrix could "
             "not disagree with it.")
    O.footer(fig, theme, O.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.885))
    hdr = ["label", "gi00", "gi01", "gi11", "gi22", "gi33", "inverse_residual"]
    return O.save_figure(
        fig, out, "metric_inverse_residual", manifest=man,
        inputs=[os.path.join(rec, "inverse_metric_samples.csv")], theme=theme,
        plotted_header=hdr, plotted_rows=[[r2[k] for k in hdr] for r2 in rows])


def fig_determinant(rec, out, man, theme):
    fr = O.read_csv(os.path.join(rec, "frames_amplitude_sweep.csv"))
    det = O.read_csv(os.path.join(rec, "determinant_checks.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 4.9),
                             facecolor=theme.background)
    ax = axes[0]
    dv = [abs(float(r["det_tx"]) + 1.0) for r in fr]
    ax.semilogy(range(len(dv)), [max(v, 1e-20) for v in dv], lw=0.8, color=RED)
    ax.axhline(1e-15, color=AMBER, lw=1.2, ls="--")
    ax.set_title(r"$|\det g_{(t,x)} + 1|$ over every sampled point and "
                 "amplitude", fontsize=9.2, color=theme.foreground)
    O.style_axes(ax, theme, "sample", "residual")

    ax = axes[1]
    labels = [r["case"] for r in det]
    vals = [max(float(r["residual"]), 1e-20) for r in det]
    ax.barh(range(len(vals)), vals, color=[BLUE, AMBER])
    ax.set_xscale("log")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7.4, color=theme.foreground)
    ax.set_title(r"$\det g = -\alpha^2\det\gamma$", fontsize=9.5,
                 color=theme.foreground)
    O.style_axes(ax, theme, "residual", "")

    fig.suptitle("O0B  the determinant identity", fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    O.banner(fig, theme,
             "det g = -alpha^2 det gamma is an IDENTITY of the decomposition, "
             "so it holds for every design vector rather than to a tolerance. "
             "For F1 it gives exactly -1, including where g_00 = 0. " +
             O.CROSSING_NOTE)
    O.footer(fig, theme, O.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.875))
    hdr = ["case", "det_expected", "det_computed", "residual"]
    return O.save_figure(
        fig, out, "metric_determinant_identity", manifest=man,
        inputs=[os.path.join(rec, "frames_amplitude_sweep.csv"),
                os.path.join(rec, "determinant_checks.csv")], theme=theme,
        plotted_header=hdr, plotted_rows=[[r[k] for k in hdr] for r in det])


def fig_spatial_eigenvalues(rec, out, man, theme):
    rows = O.read_csv(os.path.join(rec, "metric_samples.csv"))
    val = O.read_csv(os.path.join(rec, "metric_validity_checks.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 4.9),
                             facecolor=theme.background)
    ax = axes[0]
    for key, lab, c in (("g11", r"$\gamma_{xx}$", BLUE),
                        ("g22", r"$\gamma_{yy}$", AMBER),
                        ("g33", r"$\gamma_{zz}$", GREEN)):
        ax.plot(O.col(rows, key), marker="o", markersize=5, color=c, label=lab)
    ax.axhline(1.0, color=theme.grid, lw=0.9, ls=":")
    ax.set_ylim(0.0, 1.6)
    ax.legend(fontsize=7.2, frameon=False, labelcolor=theme.foreground)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([r["label"] for r in rows], fontsize=6.2, rotation=30,
                       ha="right", color=theme.foreground)
    ax.set_title(r"F1 spatial metric: $\gamma = I$, eigenvalues all exactly 1",
                 fontsize=9.2, color=theme.foreground)
    O.style_axes(ax, theme, "", "eigenvalue")

    ax = axes[1]
    pdrow = [r for r in val if r["check"] ==
             "spatial_metric_positive_definite"]
    ok = pdrow[0]["result"] if pdrow else "unknown"
    ax.text(0.5, 0.62, r"$\gamma = LL^T$", ha="center", va="center",
            fontsize=22, color=theme.foreground, transform=ax.transAxes)
    ax.text(0.5, 0.40,
            "positive definite for every finite L\nby construction: "
            r"$\det\gamma=e^{2(\ell_1+\ell_2+\ell_3)}>0$",
            ha="center", va="center", fontsize=9.5, color=theme.foreground,
            transform=ax.transAxes)
    ax.text(0.5, 0.16, "dormant F2/F3 path, exercised at nine synthetic\n"
                       "factors in the gate: %s" % ok,
            ha="center", va="center", fontsize=9, color=AMBER,
            transform=ax.transAxes)
    ax.axis("off")
    fig.suptitle("O0B  spatial-metric eigenvalues and the structural guarantee",
                 fontsize=12, color=theme.foreground, x=0.010, ha="left",
                 y=0.985)
    O.banner(fig, theme,
             "F1 fixes gamma = I, so this figure has nothing to plot but ones "
             "-- which is the honest picture. The parameterisation that will "
             "carry a nontrivial spatial metric is verified in the gate BEFORE "
             "it is switched on, so a dormant path is not silently wrong.")
    O.footer(fig, theme, O.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.885))
    return O.save_figure(
        fig, out, "spatial_metric_eigenvalues", manifest=man,
        inputs=[os.path.join(rec, "metric_samples.csv"),
                os.path.join(rec, "metric_validity_checks.csv")], theme=theme,
        source_rows="all rows of metric_samples.csv and "
                    "metric_validity_checks.csv")


def fig_derivative_channels(rec, out, man, theme):
    d1 = O.read_csv(os.path.join(rec, "metric_first_derivatives.csv"))
    d2 = O.read_csv(os.path.join(rec, "metric_second_derivatives.csv"))
    sym = O.read_csv(os.path.join(rec, "derivative_symmetry_checks.csv"))
    lab = sorted({r["label"] for r in d1})[0]
    s1 = [r for r in d1 if r["label"] == lab]
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.9),
                             facecolor=theme.background)
    ax = axes[0]
    g = [[0.0] * 4 for _ in range(4)]
    for r in s1:
        if r["A"] == "0" and r["B"] == "1":
            g[int(r["mu"])][0] = float(r["value"])
    ax.bar(range(4), [g[m][0] for m in range(4)],
           color=[BLUE, AMBER, GREEN, RED])
    ax.set_xticks(range(4))
    ax.set_xticklabels(["t", "x", "y", "z"], fontsize=9,
                       color=theme.foreground)
    ax.set_title(r"$\partial_\mu g_{01}$ at %s" % lab, fontsize=9.2,
                 color=theme.foreground)
    O.style_axes(ax, theme, "index (%s)" % O.INDEX_ORDER, "value")

    ax = axes[1]
    m2 = [[0.0] * 4 for _ in range(4)]
    for r in d2:
        if r["label"] == lab and r["A"] == "0" and r["B"] == "1":
            m2[int(r["mu"])][int(r["nu"])] = float(r["value"])
    mm = max(abs(v) for row in m2 for v in row) or 1.0
    im = ax.imshow(m2, cmap="RdBu_r", vmin=-mm, vmax=mm,
                   interpolation="nearest")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=6.5, colors=theme.foreground)
    cb.outline.set_visible(False)
    ax.set_xticks(range(4)); ax.set_yticks(range(4))
    ax.set_xticklabels(["t", "x", "y", "z"], fontsize=8,
                       color=theme.foreground)
    ax.set_yticklabels(["t", "x", "y", "z"], fontsize=8,
                       color=theme.foreground)
    ax.set_title(r"$\partial_\mu\partial_\nu g_{01}$", fontsize=9.2,
                 color=theme.foreground)
    O.style_axes(ax, theme, r"$\nu$", r"$\mu$")

    ax = axes[2]
    vals = [max(float(r["worst_residual"]), 1e-20) for r in sym]
    ax.barh([0, 1], vals, color=[BLUE, AMBER])
    ax.set_xscale("log")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["2nd order", "3rd order"], fontsize=8.5,
                       color=theme.foreground)
    ax.set_title("permutation-symmetry residual", fontsize=9.2,
                 color=theme.foreground)
    O.style_axes(ax, theme, "residual", "")

    fig.suptitle("O0B  analytic metric derivative channels", fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    O.banner(fig, theme,
             "Index ordering %s, recorded in every derivative record. Mixed "
             "partials commute, so the symmetry residual is a real check on the "
             "chain rule rather than a restatement of it: a channel assembled "
             "in the wrong order would not be symmetric." % O.INDEX_ORDER)
    O.footer(fig, theme, "Production derivatives are analytic throughout; "
                         "finite differences appear only as the independent "
                         "verification path. " + O.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.885))
    return O.save_figure(
        fig, out, "metric_derivative_channels", manifest=man,
        inputs=[os.path.join(rec, "metric_first_derivatives.csv"),
                os.path.join(rec, "metric_second_derivatives.csv"),
                os.path.join(rec, "derivative_symmetry_checks.csv")],
        theme=theme, extra={"probe_label": lab},
        source_rows="rows with label = %s" % lab)


def fig_temporal_channels(rec, out, man, theme):
    rows = O.read_csv(os.path.join(rec, "temporal_derivative_checks.csv"))
    fig, ax = plt.subplots(figsize=(11.6, 5.0), facecolor=theme.background)
    names = [r["channel"] for r in rows]
    res = [max(float(r["residual"]), 1e-20) for r in rows]
    tol = [float(r["tolerance"]) for r in rows]
    ax.barh(range(len(names)), res, color=BLUE, label="measured residual")
    ax.plot(tol, range(len(names)), marker="|", markersize=18, linestyle="none",
            color=RED, label="gate tolerance")
    ax.set_xscale("log")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9, family="monospace",
                       color=theme.foreground)
    ax.legend(fontsize=7.4, frameon=False, labelcolor=theme.foreground)
    O.style_axes(ax, theme, "residual against centred differences", "channel")
    nz = rows[0]["nonzero_magnitude"]
    fig.suptitle("O0B  temporal and mixed derivative channels", fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    O.banner(fig, theme,
             "EVERY ONE OF THESE CHANNELS VANISHES IDENTICALLY UNDER A "
             "FROZEN-TIME TREATMENT. They are nonzero here -- largest "
             "|d_t g| = %s -- because d/dt = -(v_s/R0) d/dxi is applied "
             "analytically. The frozen-time diagnostic path is not used "
             "anywhere in O0B." % nz)
    O.footer(fig, theme, O.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.880))
    hdr = ["channel", "residual", "tolerance", "nonzero_magnitude",
           "frozen_time_used"]
    return O.save_figure(
        fig, out, "temporal_derivative_channels", manifest=man,
        inputs=[os.path.join(rec, "temporal_derivative_checks.csv")],
        theme=theme, extra={"frozen_time_used": False},
        plotted_header=hdr, plotted_rows=[[r[k] for k in hdr] for r in rows])


def fig_gradient_checks(rec, out, man, theme):
    rows = O.read_csv(os.path.join(rec, "coefficient_gradient_checks.csv"))
    by = O.group_by(rows, "method")
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.0),
                             facecolor=theme.background)
    ax = axes[0]
    for i, (meth, sub) in enumerate(sorted(by.items())):
        d = [max(abs(float(r["absolute_difference"])), 1e-20) for r in sub]
        ax.semilogy([i] * len(d), d, marker="o", markersize=4,
                    linestyle="none", alpha=0.7,
                    color=[BLUE, AMBER, GREEN][i % 3], label=meth)
    ax.axhline(1e-14, color=RED, lw=1.2, ls="--", label="1e-14")
    ax.set_xticks(range(len(by)))
    ax.set_xticklabels(sorted(by), fontsize=7.6, color=theme.foreground)
    ax.legend(fontsize=6.8, frameon=False, labelcolor=theme.foreground)
    ax.set_title("|computed - exact reference| per coefficient", fontsize=9.2,
                 color=theme.foreground)
    O.style_axes(ax, theme, "method", "residual")

    ax = axes[1]
    targets = sorted({r["target"] for r in rows})
    for i, tg in enumerate(targets):
        sub = [r for r in rows if r["target"] == tg]
        ax.plot([float(r["reference"]) for r in sub],
                [float(r["computed"]) for r in sub], marker="o", markersize=4,
                linestyle="none", label=tg)
    lim = max(abs(float(r["reference"])) for r in rows) or 1.0
    ax.plot([-lim, lim], [-lim, lim], color=theme.grid, lw=0.9, ls=":")
    ax.legend(fontsize=6.6, frameon=False, labelcolor=theme.foreground)
    ax.set_title("computed against the exact reference", fontsize=9.2,
                 color=theme.foreground)
    O.style_axes(ax, theme, "exact reference", "computed")

    fig.suptitle("O0B  coefficient gradients of the design", fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    O.banner(fig, theme,
             "beta^x and g_01 are LINEAR in the coefficients, so d/dc_k = B_k "
             "is an identity. g_00 is QUADRATIC, so the reference is "
             "2 beta^x B_k -- a stronger test, because the baseline enters the "
             "derivative. " + O.HESSIAN_NOTE)
    O.footer(fig, theme, "No full design Hessian is implemented in O0B. " +
             O.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.880))
    hdr = ["method", "target", "coefficient_or_coordinate", "batch",
           "direction", "reference", "computed", "absolute_difference"]
    return O.save_figure(
        fig, out, "coefficient_gradient_checks", manifest=man,
        inputs=[os.path.join(rec, "coefficient_gradient_checks.csv")],
        theme=theme, extra={"full_hessian_implemented": False},
        plotted_header=hdr, plotted_rows=[[r[k] for k in hdr] for r in rows])


def fig_grid_cost(rec, out, man, theme):
    c = O.read_json(os.path.join(rec, "uniform_grid_cost_estimate.json"))
    q = O.read_json(os.path.join(rec, "elementwise_quadrature_manifest.json"))
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.0),
                             facecolor=theme.background)
    ax = axes[0]
    ax.bar([0, 1], [c["current_spacing_m"], c["required_spacing_m"]],
           color=[BLUE, RED], width=0.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["current A6 grid\n%.3f m" % c["current_spacing_m"],
                        "six-cell requirement\n%.4f m"
                        % c["required_spacing_m"]],
                       fontsize=8, color=theme.foreground)
    ax.set_title("grid spacing", fontsize=9.5, color=theme.foreground)
    O.style_axes(ax, theme, "", "h (m)")

    ax = axes[1]
    ax.bar([0, 1], [c["total_cells"], q["total_nodes"]], color=[RED, GREEN],
           width=0.5)
    ax.set_yscale("log")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Route A cells\n%.3g (%.1f GB)"
                        % (c["total_cells"], c["total_gigabytes"]),
                        "Route B nodes\n%d in %d elements"
                        % (q["total_nodes"], q["elements"])],
                       fontsize=8, color=theme.foreground)
    ax.set_title("evaluation points required", fontsize=9.5,
                 color=theme.foreground)
    O.style_axes(ax, theme, "", "points")

    fig.suptitle("O0B  cost of the two future evaluation routes", fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    O.banner(fig, theme,
             "Route A needs h <= %.4f m, %.3g cells and %.1f GB -- %.0fx the "
             "current cost. THE GRID WAS NOT RUN. The six-cell rule is a "
             "finite-difference criterion and applies to Route A because Route "
             "A differences; it does NOT transfer to Route B."
             % (c["required_spacing_m"], c["total_cells"],
                c["total_gigabytes"], c["cost_factor_versus_current"]))
    O.footer(fig, theme, q["caveat"] + " " + O.SCOPE)
    fig.tight_layout(rect=(0, 0.055, 1, 0.880))
    return O.save_figure(
        fig, out, "uniform_grid_cost", manifest=man,
        inputs=[os.path.join(rec, "uniform_grid_cost_estimate.json"),
                os.path.join(rec, "elementwise_quadrature_manifest.json")],
        theme=theme, extra={"grid_was_run": False, "energy_evaluated": False},
        source_rows="the two cost manifests in full")


def fig_quadrature_elements(rec, out, man, theme):
    rows = O.read_csv(os.path.join(rec, "elementwise_quadrature_elements.csv"))
    q = O.read_json(os.path.join(rec, "elementwise_quadrature_manifest.json"))
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.4),
                             facecolor=theme.background)
    ax = axes[0]
    for r in rows:
        xl, xh = float(r["xi_lower"]), float(r["xi_upper"])
        sl, sh = float(r["s_lower"]), float(r["s_upper"])
        collar = r["in_outer_collar"] == "1"
        ax.add_patch(plt.Rectangle((xl, sl), xh - xl, sh - sl, fill=True,
                                   facecolor=(RED if collar else BLUE),
                                   alpha=0.20, edgecolor=theme.grid,
                                   linewidth=0.5))
    ax.set_xlim(min(O.col(rows, "xi_lower")), max(O.col(rows, "xi_upper")))
    ax.set_ylim(min(O.col(rows, "s_lower")), max(O.col(rows, "s_upper")))
    ax.set_title("%d elements, one per unique tensor-product knot span"
                 % len(rows), fontsize=9.2, color=theme.foreground)
    O.style_axes(ax, theme, r"$\xi$ (dimensionless)", r"$s$ (dimensionless)")
    ax.legend(handles=[
        Line2D([], [], color=BLUE, lw=6, alpha=0.4, label="interior element"),
        Line2D([], [], color=RED, lw=6, alpha=0.4,
               label="element meeting the outer collar")],
        fontsize=6.8, frameon=False, labelcolor=theme.foreground)

    ax = axes[1]
    ax.plot(O.col(rows, "volume"), O.col(rows, "active_functions"), marker="o",
            markersize=4, linestyle="none", color=GREEN)
    ax.set_title("active basis functions per element", fontsize=9.2,
                 color=theme.foreground)
    O.style_axes(ax, theme, "element volume", "active functions")

    fig.suptitle("O0B  elementwise quadrature scaffold (Route B)", fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    O.banner(fig, theme,
             "Elements tile the design region exactly: their volumes sum to "
             "%.8f against the exact %.8f, residual %.1e. %d-node "
             "Gauss-Legendre per direction, exact to degree %d on each element. "
             "NO ENERGY IS EVALUATED HERE."
             % (q["total_volume"], q["exact_volume"], q["volume_residual"],
                q["nodes_per_direction"], q["polynomial_exactness_degree"]))
    O.footer(fig, theme, q["convergence_contract"] + " " + O.SCOPE)
    fig.tight_layout(rect=(0, 0.055, 1, 0.880))
    hdr = ["index", "xi_lower", "xi_upper", "s_lower", "s_upper", "volume",
           "active_functions", "in_outer_collar", "in_protected_cabin"]
    return O.save_figure(
        fig, out, "quadrature_elements", manifest=man,
        inputs=[os.path.join(rec, "elementwise_quadrature_elements.csv"),
                os.path.join(rec, "elementwise_quadrature_manifest.json")],
        theme=theme, plotted_header=hdr,
        plotted_rows=[[r[k] for k in hdr] for r in rows])


def fig_readiness(rec, out, man, theme):
    rd = O.read_json(os.path.join(rec, "readiness.json"))
    cons = O.read_json(os.path.join(rec, "o0a_resolution_consistency.json"))
    fam = rd["families"][0]
    flags = [("implemented", fam["implemented"]),
             ("metric construction verified",
              fam["metric_construction_verified"]),
             ("analytic jets verified", fam["analytic_jets_verified"]),
             ("gradient verified", fam["gradient_verified"]),
             ("resolution certified", fam["resolution_certified"]),
             ("ready for energy optimization",
              fam["ready_for_energy_optimization"])]
    fig, ax = plt.subplots(figsize=(12.4, 5.6), facecolor=theme.background)
    for i, (name, v) in enumerate(flags):
        ax.barh([i], [1.0], color=(GREEN if v else RED), alpha=0.75)
        ax.text(0.02, i, name, va="center", fontsize=10,
                color=theme.background if v else "#ffffff")
        ax.text(0.97, i, "TRUE" if v else "FALSE", va="center", ha="right",
                fontsize=10, family="monospace",
                color=theme.background if v else "#ffffff")
    ax.set_xlim(0, 1); ax.set_ylim(-0.6, len(flags) - 0.4)
    ax.set_yticks([]); ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_facecolor(theme.background)
    fig.suptitle("O0B  readiness of the F1 shift-only family", fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    O.readiness_stamp(fig, theme)
    O.banner(fig, theme,
             "Blocking reason: %s\nF1 minimum %.3f cells and F0 minimum %.3f "
             "cells against a six-cell contract (canonical O0A record)."
             % (fam["blocking_reason"],
                cons["f1_minimum_cells_across_narrowest_feature"],
                cons["f0_minimum_cells_across_narrowest_feature"]),
             y=0.855)
    O.footer(fig, theme, "Energy evaluated: %s. Optimizer run: %s. %s"
             % (rd["energy_evaluated"], rd["optimizer_run"], O.SCOPE))
    fig.tight_layout(rect=(0, 0.055, 1, 0.760))
    return O.save_figure(
        fig, out, "o0b_readiness", manifest=man,
        inputs=[os.path.join(rec, "readiness.json"),
                os.path.join(rec, "o0a_resolution_consistency.json")],
        theme=theme,
        extra={"required_statement": rd["required_statement"],
               "blocking_reason": fam["blocking_reason"]},
        source_rows="readiness.json and o0a_resolution_consistency.json in full")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--records", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--theme", default="scientific")
    args = ap.parse_args(argv)

    rec = records.resolve(args.records)
    out = records.ensure_dir(args.out)
    theme = themes.get(args.theme)
    themes.apply(theme)

    missing = [n for n in O.REQUIRED_RECORDS
               if not os.path.exists(os.path.join(rec, n))]
    if missing:
        print("missing required records: %r" % missing, file=sys.stderr)
        return 1

    man = O.read_json(os.path.join(rec, "metric_design_manifest.json"))
    before = O.source_checksums(rec)
    with open(os.path.join(out, "source_checksums.json"), "w",
              encoding="utf-8") as fh:
        json.dump(before, fh, indent=2, sort_keys=True)

    made = [
        fig_adm_components(rec, out, man, theme),
        fig_baseline_vs_zero(rec, out, man, theme),
        fig_single_mode(rec, out, man, theme),
        fig_inverse_residual(rec, out, man, theme),
        fig_determinant(rec, out, man, theme),
        fig_spatial_eigenvalues(rec, out, man, theme),
        fig_derivative_channels(rec, out, man, theme),
        fig_temporal_channels(rec, out, man, theme),
        fig_gradient_checks(rec, out, man, theme),
        fig_grid_cost(rec, out, man, theme),
        fig_quadrature_elements(rec, out, man, theme),
        fig_readiness(rec, out, man, theme),
    ]

    after = O.source_checksums(rec)
    rep = anomaly.AnomalyReport("O", "O0B", "o0b_static")
    changed = [k for k in before if before[k] != after.get(k)]
    rep.record("source_checksum_mismatch", fired=bool(changed),
               measured=len(changed), threshold=0, detail=repr(changed[:3]))
    blanks = [os.path.basename(p) for p in made if _blank(p)]
    rep.record("blank_figure", fired=bool(blanks), measured=len(blanks),
               threshold=0, detail=repr(blanks[:3]))
    rd = O.read_json(os.path.join(rec, "readiness.json"))
    ready = [f["family"] for f in rd["families"]
             if f["ready_for_energy_optimization"]]
    rep.record("false_optimization_readiness_claim", fired=bool(ready),
               measured=len(ready), threshold=0, detail=repr(ready))
    rep.record("full_hessian_overclaim",
               fired=bool(man.get("full_hessian_implemented")),
               measured=man.get("full_hessian_implemented"), threshold=False,
               detail=O.HESSIAN_NOTE)
    rep.record("energy_evaluated_in_o0b",
               fired=bool(man.get("energy_evaluated")) or
                     bool(man.get("optimizer_run")),
               measured=[man.get("energy_evaluated"), man.get("optimizer_run")],
               threshold=False,
               detail="O0B must evaluate no energy and run no optimizer")
    path = rep.write(out)
    for p in made:
        print("wrote %s" % os.path.basename(p))
    print("%d figures, %d fired anomalies -> %s"
          % (len(made), len(rep.fired), os.path.basename(path)))
    return 0


def _blank(path, min_colors=64):
    try:
        from PIL import Image
    except ImportError:
        return False
    im = Image.open(path).convert("RGB")
    c = im.getcolors(maxcolors=1 << 20)
    return c is not None and len(c) < min_colors


if __name__ == "__main__":
    raise SystemExit(main())
