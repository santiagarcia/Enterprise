#!/usr/bin/env python3
"""The thirteen O0C energy figures, drawn from saved C++ records only.

    python plot_o0c_energy.py --records <o0c_energy dir> --out <gallery dir>

Nine two-dimensional scientific figures and four three-dimensional views.  The
3-D views are the start of the real volumetric architecture: production writes
sampled fields, this renderer extracts isosurfaces, revolves meridional data,
projects a camera and draws.  It never evaluates physics.
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
from matplotlib.lines import Line2D                              # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

import o0c_common as O                                           # noqa: E402
from common import records                                       # noqa: E402
from themes import themes                                        # noqa: E402

BLUE, AMBER, RED, GREEN = "#4ea1ff", "#ffb454", "#d1495b", "#7ed957"
VIOLET = "#b07cff"


# --- 2-D scientific figures ---------------------------------------------------

def fig_density_reference(rec, out, man, theme):
    rows = O.read_csv(os.path.join(rec, "pointwise_objective_samples.csv"))
    r = O.col(rows, "r")
    rho = O.col(rows, "rho")
    ref = O.col(rows, "rho_analytic")
    sam = O.read_csv(os.path.join(rec, "alcubierre_density_reference.csv"))

    fig, axes = plt.subplots(1, 2, figsize=(13.4, 4.9),
                             facecolor=theme.background)
    axes[0].plot(r, rho, color=BLUE, lw=2.0,
                 label=r"$\rho$ from $G_{\mu\nu}n^\mu n^\nu/8\pi$")
    axes[0].plot(r, ref, color=AMBER, lw=1.1, ls="--",
                 label=r"$-\frac{v_s^2}{32\pi}\left[(\partial_yf)^2"
                       r"+(\partial_zf)^2\right]$")
    O.style_axes(axes[0], theme, r"transverse radius $\rho$  [m]",
                 r"Eulerian energy density  $\rho$")
    axes[0].legend(fontsize=7.6, facecolor=theme.background,
                   edgecolor=theme.grid, labelcolor=theme.foreground)

    diff = [abs(a - b) for a, b in zip(rho, ref)]
    axes[1].semilogy(r, [max(d, 1e-22) for d in diff], color=RED, lw=1.6)
    axes[1].axhline(2.2e-16, color=theme.grid, lw=0.8, ls=":")
    axes[1].text(0.02, 2.6e-16, "double epsilon", fontsize=6.4,
                 color=theme.grid)
    O.style_axes(axes[1], theme, r"transverse radius $\rho$  [m]",
                 "absolute difference")

    worst = max((abs(float(s["absolute_difference"])) for s in sam), default=0.0)
    O.banner(fig, theme,
             "The analytic Alcubierre density, reproduced by the full "
             "curvature chain to %.1e at every sampled site" % worst)
    O.footer(fig, theme,
             "Left: the density computed through metric jet -> connection -> "
             "Riemann -> Einstein -> Eulerian observer, against the closed form. "
             "Right: their difference, at roundoff. " + O.NOT_AN_ENERGY_CLAIM)
    O.readiness_stamp(fig, theme)
    fig.tight_layout(rect=[0, 0.03, 1, 0.92])
    return O.save_figure(
        fig, out, "pointwise_alcubierre_density_reference", manifest=man,
        inputs=[os.path.join(rec, "pointwise_objective_samples.csv"),
                os.path.join(rec, "alcubierre_density_reference.csv")],
        theme=theme, plotted_header=["r", "rho", "rho_analytic", "difference"],
        plotted_rows=[[a, b, c, d] for a, b, c, d in zip(r, rho, ref, diff)])


def fig_objective_decomposition(rec, out, man, theme):
    rows = O.read_csv(os.path.join(rec, "objective_values.csv"))
    v = rows[0]
    e_plus, e_neg = float(v["E_plus"]), float(v["E_neg"])
    e_net, e_abs = float(v["E_net"]), float(v["E_abs"])
    samples = O.read_csv(os.path.join(rec, "pointwise_objective_samples.csv"))

    fig, axes = plt.subplots(1, 2, figsize=(13.4, 4.9),
                             facecolor=theme.background)
    names = [r"$E_+$", r"$E_{\rm neg}$", r"$E_{\rm abs}$", r"$E_{\rm net}$"]
    vals = [e_plus, e_neg, e_abs, e_net]
    cols = [GREEN, RED, AMBER, BLUE]
    axes[0].bar(names, vals, color=cols, edgecolor=theme.grid, linewidth=0.6)
    axes[0].axhline(0.0, color=theme.foreground, lw=0.8)
    O.style_axes(axes[0], theme, "", "integrated value")
    axes[0].text(0.5, 0.92,
                 r"$E_{\rm net}=E_+-E_{\rm neg}$   "
                 r"$E_{\rm abs}=E_++E_{\rm neg}$",
                 transform=axes[0].transAxes, ha="center", fontsize=8.4,
                 color=theme.foreground)
    axes[0].text(0.5, 0.83,
                 r"$E_{\rm neg}$ is the POSITIVE magnitude",
                 transform=axes[0].transAxes, ha="center", fontsize=7.4,
                 color=theme.uncertainty_color)

    r = O.col(samples, "r")
    axes[1].plot(r, O.col(samples, "integrand_E_neg"), color=RED, lw=1.8,
                 label=r"$\max(-\rho,0)$")
    axes[1].plot(r, O.col(samples, "integrand_J_NEC"), color=VIOLET, lw=1.4,
                 label=r"$\max(-N_{\min},0)$")
    O.style_axes(axes[1], theme, r"transverse radius $\rho$  [m]",
                 "integrand")
    axes[1].legend(fontsize=7.6, facecolor=theme.background,
                   edgecolor=theme.grid, labelcolor=theme.foreground)

    O.banner(fig, theme,
             r"$E_+ = 0$ exactly for the shift-only baseline: $\rho \leq 0$ "
             r"everywhere, so $E_{\rm net} = -E_{\rm neg}$")
    O.footer(fig, theme,
             "Closure residuals are exactly zero: these are accumulation "
             "identities, not convergent statements. Source semantics: %s, "
             "Lambda = %s. %s"
             % (v.get("route", "elementwise"), man.get("lambda", 0.0),
                O.NOT_AN_ENERGY_CLAIM))
    O.readiness_stamp(fig, theme)
    fig.tight_layout(rect=[0, 0.03, 1, 0.92])
    return O.save_figure(
        fig, out, "energy_objective_decomposition", manifest=man,
        inputs=[os.path.join(rec, "objective_values.csv"),
                os.path.join(rec, "pointwise_objective_samples.csv")],
        theme=theme, plotted_header=["quantity", "value"],
        plotted_rows=[["E_plus", e_plus], ["E_neg", e_neg],
                      ["E_abs", e_abs], ["E_net", e_net]])


def _convergence_figure(rec, out, man, theme, csv_name, stem, xkey, xlabel,
                        title, note):
    rows = O.read_csv(os.path.join(rec, csv_name))
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 4.7),
                             facecolor=theme.background)
    plotted = []
    for ax, obj, c in ((axes[0], "E_net", BLUE), (axes[1], "E_neg", RED)):
        sel = [r for r in rows if r["objective"] == obj]
        xs = [float(r[xkey]) for r in sel]
        ys = [float(r["value"]) for r in sel]
        ch = [float(r["absolute_change"]) for r in sel]
        ax.plot(xs, ys, "o-", color=c, lw=1.8, ms=5)
        O.style_axes(ax, theme, xlabel, obj)
        tw = ax.twinx()
        tw.semilogy(xs[1:], [max(v, 1e-20) for v in ch[1:]], "s--",
                    color=AMBER, lw=1.1, ms=3.6)
        tw.set_ylabel("|change| from previous", color=AMBER, fontsize=8)
        tw.tick_params(colors=AMBER, labelsize=7)
        plotted += [[obj, a, b, d] for a, b, d in zip(xs, ys, ch)]
    O.banner(fig, theme, title)
    O.footer(fig, theme, note + " " + O.NOT_AN_ENERGY_CLAIM)
    O.readiness_stamp(fig, theme)
    fig.tight_layout(rect=[0, 0.03, 1, 0.92])
    return O.save_figure(
        fig, out, stem, manifest=man,
        inputs=[os.path.join(rec, csv_name)], theme=theme,
        plotted_header=["objective", xkey, "value", "absolute_change"],
        plotted_rows=plotted)


def fig_order(rec, out, man, theme):
    return _convergence_figure(
        rec, out, man, theme, "quadrature_order_convergence.csv",
        "quadrature_order_convergence", "order", "Gauss-Legendre order",
        "Convergence in quadrature order, 4 through 12",
        "The integrand is not polynomial -- exp and tanh of splines divided by "
        "a determinant -- so the rule's formal degree-15 exactness says nothing "
        "here and this measurement is the only evidence.")


def fig_subdivision(rec, out, man, theme):
    return _convergence_figure(
        rec, out, man, theme, "subdivision_convergence.csv",
        "subdivision_convergence", "subdivision", "element subdivision",
        "Convergence under element subdivision, 1, 2 and 4 ways per direction",
        "Order refinement and h refinement fail differently: a rule can be "
        "order-converged on each element and still be resolving the wrong "
        "elements.")


def fig_domain(rec, out, man, theme):
    return _convergence_figure(
        rec, out, man, theme, "domain_convergence.csv",
        "domain_convergence", "r_max_comoving", "domain radius  [m]",
        "Convergence under domain enlargement, three deterministic extents",
        "The verdict is taken on the FINAL refinement. Requiring every "
        "consecutive pair to satisfy the criterion would demand that the "
        "coarsest step already be converged.")


def fig_surrogate_bias(rec, out, man, theme):
    rows = [r for r in O.read_csv(os.path.join(rec, "smooth_surrogate_bias.csv"))
            if r.get("E_neg_exact")]
    tau = [float(r["tau"]) for r in rows]
    bias = [float(r["absolute_bias"]) for r in rows]
    bound = [float(r["predicted_upper_bound"]) for r in rows]
    exact = float(rows[0]["E_neg_exact"]) if rows else 0.0

    fig, axes = plt.subplots(1, 2, figsize=(13.4, 4.7),
                             facecolor=theme.background)
    axes[0].loglog(tau, [max(b, 1e-20) for b in bias], "o-", color=RED, lw=1.9,
                   label="measured bias")
    axes[0].loglog(tau, bound, "s--", color=AMBER, lw=1.2,
                   label=r"$\tau\log 2 \times$ volume")
    O.style_axes(axes[0], theme, r"smoothing width $\tau$",
                 r"$E_{\rm neg,\tau} - E_{\rm neg}$")
    axes[0].legend(fontsize=7.6, facecolor=theme.background,
                   edgecolor=theme.grid, labelcolor=theme.foreground)

    axes[1].semilogx(tau, [float(r["E_neg_smooth"]) for r in rows], "o-",
                     color=VIOLET, lw=1.9, label=r"$E_{\rm neg,\tau}$")
    axes[1].axhline(exact, color=BLUE, lw=1.4, ls="--",
                    label=r"exact $E_{\rm neg}$")
    O.style_axes(axes[1], theme, r"smoothing width $\tau$", "value")
    axes[1].legend(fontsize=7.6, facecolor=theme.background,
                   edgecolor=theme.grid, labelcolor=theme.foreground)

    O.banner(fig, theme,
             "The softplus surrogate approaches the exact objective FROM "
             "ABOVE, and the bias is measured rather than assumed negligible")
    O.footer(fig, theme,
             "softplus exceeds max(u,0) by exactly tau log 2 at u = 0, which is "
             "where certification happens, so it is differentiated but never "
             "substituted for the exact value. " + O.NOT_AN_ENERGY_CLAIM)
    O.readiness_stamp(fig, theme)
    fig.tight_layout(rect=[0, 0.03, 1, 0.92])
    return O.save_figure(
        fig, out, "smooth_surrogate_bias", manifest=man,
        inputs=[os.path.join(rec, "smooth_surrogate_bias.csv")], theme=theme,
        plotted_header=["tau", "bias", "bound"],
        plotted_rows=[[a, b, c] for a, b, c in zip(tau, bias, bound)])


def fig_gradient_agreement(rec, out, man, theme):
    rows = O.read_csv(os.path.join(rec, "objective_gradient_checks.csv"))
    batch = O.read_csv(os.path.join(rec, "gradient_batch_schedule.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 4.9),
                             facecolor=theme.background)

    refs = [float(r["reference"]) for r in rows]
    gots = [float(r["computed"]) for r in rows]
    if refs:
        lim = max(max(abs(v) for v in refs), max(abs(v) for v in gots)) * 1.15
        axes[0].plot([-lim, lim], [-lim, lim], color=theme.grid, lw=0.9, ls=":")
        for r, g, row in zip(refs, gots, rows):
            c = (BLUE if "richardson" in row["method"] else
                 (GREEN if "separate" in row["method"] else AMBER))
            axes[0].plot(r, g, "o", color=c, ms=6, alpha=0.85)
        axes[0].set_xlim(-lim, lim)
        axes[0].set_ylim(-lim, lim)
    O.style_axes(axes[0], theme, "reference derivative", "OTI derivative")
    axes[0].legend(handles=[
        Line2D([], [], marker="o", ls="", color=BLUE,
               label="Richardson difference"),
        Line2D([], [], marker="o", ls="", color=GREEN,
               label="separate vs batched seeding"),
        Line2D([], [], marker="o", ls="", color=AMBER,
               label="directional derivative")],
        fontsize=7.2, facecolor=theme.background, edgecolor=theme.grid,
        labelcolor=theme.foreground, loc="upper left")

    sizes = [int(r["coordinates_in_batch"]) for r in batch]
    axes[1].bar(range(len(sizes)), sizes, color=VIOLET,
                edgecolor=theme.grid, linewidth=0.6)
    O.style_axes(axes[1], theme, "OTI batch index", "coordinates seeded")
    axes[1].set_ylim(0, 4)
    axes[1].text(0.5, 0.88,
                 "%d whitened coordinates, %d batches, each covered once"
                 % (int(batch[0]["whitened_dimension"]), len(sizes)),
                 transform=axes[1].transAxes, ha="center", fontsize=8,
                 color=theme.foreground)

    O.banner(fig, theme,
             r"Gradients in the WHITENED coordinates $y$, through "
             r"$c = c_0 + PWy$ -- not $z$-gradients through $P$")
    O.footer(fig, theme, O.HESSIAN_NOTE + " " + O.NOT_AN_ENERGY_CLAIM)
    O.readiness_stamp(fig, theme)
    fig.tight_layout(rect=[0, 0.03, 1, 0.92])
    return O.save_figure(
        fig, out, "objective_gradient_agreement", manifest=man,
        inputs=[os.path.join(rec, "objective_gradient_checks.csv"),
                os.path.join(rec, "gradient_batch_schedule.csv")],
        theme=theme, plotted_header=["method", "reference", "computed"],
        plotted_rows=[[r["method"], r["reference"], r["computed"]]
                      for r in rows])


def fig_behavior_lock(rec, out, man, theme):
    rows = O.read_csv(os.path.join(rec, "family_mode_accounting.csv"))
    adv = O.read_csv(os.path.join(rec, "adversarial_candidate_checks.csv"))
    contract = O.read_json(os.path.join(rec, "behavior_contract.json"))
    f1w = O.read_json(os.path.join(rec, "f1w_basis_manifest.json"))

    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.1),
                             facecolor=theme.background)

    f1 = [r for r in rows if r["family"] == "F1"]
    f1w_rows = [r for r in rows if r["family"] == "F1W"]
    labels = ["F1 interior", "F1W protected"]
    f1_worst = max((abs(float(r["delta_beta"])) for r in f1
                    if r["inside_protected_claim"] == "1"), default=0.0)
    f1w_worst = max((abs(float(r["delta_beta"])) for r in f1w_rows), default=0.0)
    axes[0].bar(labels, [max(f1_worst, 1e-20), max(f1w_worst, 1e-20)],
                color=[RED, GREEN], edgecolor=theme.grid, linewidth=0.6)
    axes[0].set_yscale("log")
    O.style_axes(axes[0], theme, "",
                 r"worst $|\Delta\beta^x|$ in the protected region")
    axes[0].text(0.5, 0.9,
                 "F1: %.2e (volume mean over a CYLINDER)\n"
                 "F1W: %.0f, EXACTLY zero (support exclusion)"
                 % (f1_worst, f1w_worst),
                 transform=axes[0].transAxes, ha="center", fontsize=7.8,
                 color=theme.foreground)

    verdicts = {}
    for r in adv:
        verdicts[r["verdict"]] = verdicts.get(r["verdict"], 0) + 1
    ks = list(verdicts.keys())
    axes[1].bar(ks, [verdicts[k] for k in ks],
                color=[GREEN if k == "comparable" else RED for k in ks],
                edgecolor=theme.grid, linewidth=0.6)
    O.style_axes(axes[1], theme, "", "adversarial candidates")
    axes[1].tick_params(axis="x", labelrotation=12)
    axes[1].text(0.5, 0.9,
                 "every candidate leaves the protected regions untouched;\n"
                 "a lower value from a non-comparable design is not an "
                 "improvement",
                 transform=axes[1].transAxes, ha="center", fontsize=7.4,
                 color=theme.foreground)

    O.banner(fig, theme,
             "F1's cabin is a CYLINDER constrained in the mean; F1W's is a "
             "BALL enforced by support exclusion")
    O.footer(fig, theme,
             "F1W removed %d modes for cabin overlap and %d for collar overlap, "
             "leaving %d active. %s %s"
             % (f1w["removed_cabin_overlap"], f1w["removed_collar_overlap"],
                f1w["remaining_active_modes"], contract["verdict"],
                O.NOT_AN_ENERGY_CLAIM))
    O.readiness_stamp(fig, theme)
    fig.tight_layout(rect=[0, 0.04, 1, 0.92])
    return O.save_figure(
        fig, out, "behavior_lock_verification", manifest=man,
        inputs=[os.path.join(rec, "family_mode_accounting.csv"),
                os.path.join(rec, "adversarial_candidate_checks.csv"),
                os.path.join(rec, "behavior_contract.json"),
                os.path.join(rec, "f1w_basis_manifest.json")],
        theme=theme, plotted_header=["family", "worst_delta_beta"],
        plotted_rows=[["F1", f1_worst], ["F1W", f1w_worst]])


def fig_readiness(rec, out, man, theme):
    rd = O.read_json(os.path.join(rec, "o0c_objective_readiness.json"))
    flags = ["implemented", "pointwise_verified", "behavior_preserving",
             "quadrature_converged", "domain_converged",
             "independent_reference_agreement", "gradient_verified",
             "resolution_certified", "ready_for_optimization"]
    objs = rd["objectives"]

    fig, ax = plt.subplots(figsize=(13.0, 4.6), facecolor=theme.background)
    for j, o in enumerate(objs):
        for i, fl in enumerate(flags):
            v = bool(o.get(fl))
            ax.add_patch(plt.Rectangle((i, j), 0.92, 0.86,
                                       facecolor=GREEN if v else RED,
                                       edgecolor=theme.grid, linewidth=0.6,
                                       alpha=0.85))
            ax.text(i + 0.46, j + 0.43, "T" if v else "F", ha="center",
                    va="center", fontsize=8, color="#0b0d12")
    ax.set_xlim(0, len(flags))
    ax.set_ylim(0, len(objs))
    ax.set_xticks([i + 0.46 for i in range(len(flags))])
    ax.set_xticklabels([f.replace("_", "\n") for f in flags], fontsize=6.6)
    ax.set_yticks([j + 0.43 for j in range(len(objs))])
    ax.set_yticklabels([o["objective"] for o in objs], fontsize=8)
    O.style_axes(ax, theme)
    ax.grid(False)

    O.banner(fig, theme,
             "O1 may not begin: no objective has every required flag true")
    O.footer(fig, theme, rd["decision"] + " " + O.NOT_AN_ENERGY_CLAIM)
    O.readiness_stamp(fig, theme)
    fig.tight_layout(rect=[0, 0.05, 1, 0.91])
    return O.save_figure(
        fig, out, "objective_readiness", manifest=man,
        inputs=[os.path.join(rec, "o0c_objective_readiness.json")], theme=theme,
        extra={"o1_may_begin": rd["o1_may_begin"]},
        plotted_header=["objective"] + flags,
        plotted_rows=[[o["objective"]] + [bool(o.get(f)) for f in flags]
                      for o in objs])


# --- 3-D views ---------------------------------------------------------------

def _shell_faces(mer, field, level):
    xs = [mer["x"][i][0] for i in range(mer["nx"])]
    rs = [mer["rho"][0][j] for j in range(mer["nr"])]
    segs = O.isocontour(mer[field], xs, rs, level)
    return O.revolve(segs, 40), len(segs)


def _draw_shell(ax, faces, theme, color, azim, elev, alpha=0.30, lw=0.25):
    ordered = O.painter_order(faces, azim, elev)
    polys = [pts for _, pts in ordered]
    if not polys:
        return
    depths = [d for d, _ in ordered]
    lo, hi = min(depths), max(depths)
    shades = []
    for d in depths:
        t = 0.0 if hi == lo else (hi - d) / (hi - lo)
        shades.append((0.45 + 0.55 * t))
    pc = PolyCollection(polys, closed=True, linewidths=lw,
                        edgecolors=(1, 1, 1, 0.06))
    import matplotlib.colors as mc
    base = mc.to_rgb(color)
    pc.set_facecolors([(base[0] * s, base[1] * s, base[2] * s, alpha)
                       for s in shades])
    ax.add_collection(pc)


def _setup_3d_axes(ax, theme, extent=7.0):
    ax.set_facecolor(theme.background)
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent * 0.68, extent * 0.68)
    ax.set_aspect("equal")
    ax.axis("off")


def _cabin_faces(radius):
    faces = []
    n = 22
    for i in range(n):
        for j in range(n // 2):
            t0, t1 = 2 * math.pi * i / n, 2 * math.pi * (i + 1) / n
            p0 = math.pi * j / (n // 2)
            p1 = math.pi * (j + 1) / (n // 2)
            def pt(t, p):
                return (radius * math.cos(p), radius * math.sin(p) * math.cos(t),
                        radius * math.sin(p) * math.sin(t))
            faces.append((pt(t0, p0), pt(t1, p0), pt(t1, p1), pt(t0, p1)))
    return faces


def _fig_3d(rec, out, man, theme, stem, field, colour, title, note,
            with_cabin=True, with_arrows=False, azim=0.62, elev=0.34):
    mer = O.meridional(rec, 6)
    scales = O.read_json(os.path.join(rec, "field_scales.json"))
    key = {"energy_density": "negative_energy", "n_min": "nec_violation",
           "kretschmann": "curvature_shell", "beta_x": "wall_shift"}[field]
    level = scales["isovalues"][key]
    faces, nseg = _shell_faces(mer, field, level)

    fig, ax = plt.subplots(figsize=(10.6, 6.6), facecolor=theme.background)
    _setup_3d_axes(ax, theme)
    _draw_shell(ax, faces, theme, colour, azim, elev, alpha=0.34)
    if with_cabin:
        _draw_shell(ax, _cabin_faces(scales["protected_cabin_radius_comoving"]),
                    theme, GREEN, azim, elev, alpha=0.55, lw=0.15)
    if with_arrows:
        rows = O.read_csv(os.path.join(rec, "shift_vectors_sparse.csv"))
        pts, tips = [], []
        for r in rows:
            b = float(r["beta_x"])
            if abs(b) < 0.04:
                continue
            p = (float(r["x"]), float(r["y"]), float(r["z"]))
            pts.append(p)
            tips.append((p[0] + 2.4 * b, p[1], p[2]))
        for a, b in zip(O.project(pts, azim, elev),
                        O.project(tips, azim, elev)):
            ax.annotate("", xy=(b[0], b[1]), xytext=(a[0], a[1]),
                        arrowprops=dict(arrowstyle="-|>", color=AMBER,
                                        lw=0.8, alpha=0.75))

    O.banner(fig, theme, title, y=0.955)
    fig.text(0.5, 0.905, O.RECONSTRUCTION_NOTE, ha="center", fontsize=7.4,
             color=theme.uncertainty_color)
    fig.text(0.5, 0.868, "isovalue %.4g, %d meridional segments revolved"
             % (level, nseg), ha="center", fontsize=6.8, color=theme.grid)
    O.footer(fig, theme, note + " " + O.NOT_AN_ENERGY_CLAIM)
    O.readiness_stamp(fig, theme)
    return O.save_figure(
        fig, out, stem, manifest=man,
        inputs=[os.path.join(rec, "meridional_field_006.csv"),
                os.path.join(rec, "field_scales.json")],
        theme=theme, reconstruction=True,
        camera="perspective, azimuth %.2f rad, elevation %.2f rad"
               % (azim, elev),
        extra={"isovalue": level, "meridional_segments": nseg},
        plotted_header=["isovalue", "segments"], plotted_rows=[[level, nseg]])


def fig_wall_energy_3d(rec, out, man, theme):
    return _fig_3d(rec, out, man, theme, "wall_energy_3d", "energy_density",
                   RED,
                   "Negative Eulerian energy region, with the protected cabin",
                   "The shell is an isosurface of the Eulerian energy density "
                   "computed by production; the green ball is the F1W protected "
                   "interior.", with_arrows=True)


def fig_nec_3d(rec, out, man, theme):
    return _fig_3d(rec, out, man, theme, "nec_violation_3d", "n_min", VIOLET,
                   "Null-energy-condition violation region",
                   "An isosurface of the exact global null minimum N_min over "
                   "the sphere of directions, minimised by production's "
                   "trust-region solver at every sample point.")


def fig_curvature_3d(rec, out, man, theme):
    return _fig_3d(rec, out, man, theme, "curvature_shell_3d", "kretschmann",
                   AMBER, "Kretschmann curvature shell",
                   "K = R_abcd R^abcd, an invariant, so this surface is "
                   "coordinate-independent in a way the metric components are "
                   "not.")


def fig_combined(rec, out, man, theme):
    mer = O.meridional(rec, 6)
    scales = O.read_json(os.path.join(rec, "field_scales.json"))
    fig = plt.figure(figsize=(15.0, 6.4), facecolor=theme.background)
    ax3 = fig.add_axes([0.02, 0.06, 0.52, 0.80])
    _setup_3d_axes(ax3, theme)
    azim, elev = 0.7, 0.32
    for field, key, colour, alpha in (
            ("energy_density", "negative_energy", RED, 0.30),
            ("kretschmann", "curvature_shell", AMBER, 0.20)):
        faces, _ = _shell_faces(mer, field, scales["isovalues"][key])
        _draw_shell(ax3, faces, theme, colour, azim, elev, alpha=alpha)
    _draw_shell(ax3, _cabin_faces(scales["protected_cabin_radius_comoving"]),
                theme, GREEN, azim, elev, alpha=0.55, lw=0.15)

    ax2 = fig.add_axes([0.60, 0.13, 0.37, 0.72])
    grid = [[mer["energy_density"][i][j] for i in range(mer["nx"])]
            for j in range(mer["nr"])]
    xs = [mer["x"][i][0] for i in range(mer["nx"])]
    rs = [mer["rho"][0][j] for j in range(mer["nr"])]
    lo = scales["energy_density"]["min"]
    im = ax2.imshow(grid, origin="lower",
                    extent=[min(xs), max(xs), min(rs), max(rs)],
                    aspect="auto", cmap="magma", vmin=lo, vmax=0.0,
                    interpolation="nearest")
    cb = fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=6.5, colors=theme.foreground)
    O.style_axes(ax2, theme, "x  [m]", r"$\rho$  [m]")
    ax2.set_title("meridional slice, same fixed scale", fontsize=8.4,
                  color=theme.foreground)

    O.banner(fig, theme,
             "Three-dimensional shells beside the meridional data they were "
             "revolved from", y=0.955)
    fig.text(0.27, 0.905, O.RECONSTRUCTION_NOTE, ha="center", fontsize=7.4,
             color=theme.uncertainty_color)
    O.footer(fig, theme,
             "Left: energy (red) and curvature (amber) shells with the "
             "protected cabin. Right: the production half-plane itself. " +
             O.NOT_AN_ENERGY_CLAIM)
    O.readiness_stamp(fig, theme)
    return O.save_figure(
        fig, out, "combined_3d_and_slice", manifest=man,
        inputs=[os.path.join(rec, "meridional_field_006.csv"),
                os.path.join(rec, "field_scales.json")],
        theme=theme, reconstruction=True,
        camera="perspective, azimuth %.2f rad, elevation %.2f rad" % (azim, elev),
        plotted_header=["field", "isovalue"],
        plotted_rows=[["energy_density", scales["isovalues"]["negative_energy"]],
                      ["kretschmann", scales["isovalues"]["curvature_shell"]]])


FIGURES = (
    fig_density_reference, fig_objective_decomposition, fig_order,
    fig_subdivision, fig_domain, fig_surrogate_bias, fig_gradient_agreement,
    fig_behavior_lock, fig_readiness,
    fig_wall_energy_3d, fig_nec_3d, fig_curvature_3d, fig_combined,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--theme", default="scientific")
    a = ap.parse_args()

    rec = records.resolve(a.records)
    out = records.ensure_dir(a.out)
    man = O.read_json(os.path.join(rec, "geometry_bridge_manifest.json"))
    theme = themes.get(a.theme)
    themes.apply(theme)

    made = []
    for fn in FIGURES:
        made.append(fn(rec, out, man, theme))
        print("wrote", os.path.basename(made[-1]))

    with open(os.path.join(out, "figure_manifest.json"), "w",
              encoding="utf-8") as f:
        json.dump({"figures": [os.path.basename(p) for p in made],
                   "count": len(made),
                   "two_dimensional": list(O.STATICS),
                   "three_dimensional": list(O.STATICS_3D),
                   "scope": O.SCOPE,
                   "renderer_computes_physics": False,
                   "not_an_energy_claim": O.NOT_AN_ENERGY_CLAIM}, f, indent=2)
    print("%d figures" % len(made))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
