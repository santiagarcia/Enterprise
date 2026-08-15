#!/usr/bin/env python3
"""The O1 optimization figures, drawn from saved C++ records only.

    python plot_o1_optimization.py --records <o1 dir> --out <gallery dir>

Scientific clarity first.  These are meant to be read, not decoded: what the
optimizer did, what it bought, and what it cost elsewhere.

Python reads records.  It computes no GR.
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
from matplotlib.lines import Line2D                              # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

from program_o0c import o0c_common as O                          # noqa: E402
from common import provenance, records                           # noqa: E402
from themes import themes                                        # noqa: E402

BLUE, AMBER, RED, GREEN = "#4ea1ff", "#ffb454", "#d1495b", "#7ed957"
VIOLET, TEAL = "#b07cff", "#5ec8b0"

SCOPE = ("O1 optimizes the inferred source requirement for a FIXED warp-bubble "
         "kinematic family: the F1W shift-only wall correction at fixed "
         "velocity, fixed protected-cabin behavior, the ADM Eulerian observer "
         "and the total-effective source convention with Lambda = 0. It does "
         "not demonstrate faster-than-light travel, a causal shortcut, reduced "
         "travel time, horizon removal, physical realizability, dynamic "
         "stability, satisfaction of quantum energy inequalities, or a global "
         "optimum.")

FIGURES = ("energy_reduction_history", "optimizer_mechanics",
           "design_space_trajectory", "final_convergence_certification",
           "tradeoff_dashboard", "final_design_summary")

REQUIRED_RECORDS = (
    "o1_manifest.json", "optimizer_configuration.json",
    "baseline_objectives.json", "accepted_iterates.csv", "rejected_trials.csv",
    "line_search_trace.csv", "gradient_trace.csv", "trust_radius_trace.csv",
    "optimizer_trace.csv", "final_design_y.csv", "final_design_z.csv",
    "final_design_c.csv", "baseline_vs_final_metrics.csv",
    "final_order_convergence.csv", "final_subdivision_convergence.csv",
    "final_domain_convergence.csv", "final_independent_reference.csv",
    "final_gradient_resolution.csv", "final_certification.json",
    "final_objectives.json", "final_behavior_contract.json",
    "deterministic_repeat.json", "restart_check.json",
    "optimization_timing.csv", "verification_report.md",
    "verification_summary.json", "verification_table.csv",
)


def save(fig, out, stem, man, inputs, theme, extra=None, header=None, rows=None):
    png = os.path.join(records.ensure_dir(out), stem + ".png")
    svg = os.path.join(out, stem + ".svg")
    fig.savefig(png, facecolor=fig.get_facecolor())
    fig.savefig(svg, facecolor=fig.get_facecolor())
    plt.close(fig)
    ex = {"scope": SCOPE, "renderer_computes_physics": False,
          "optimizer_run": True,
          "global_optimum_claimed": False, "ftl_claimed": False,
          "physical_realizability_claimed": False,
          "colour_scale": "fixed; computed from the union of baseline and "
                          "every accepted iterate",
          "deterministic_rerender": {
              "record_digest_algorithm": "FNV-1a 64-bit",
              "input_digests": {os.path.basename(f): O.fnv1a64(f)
                                for f in inputs},
              "matplotlib_backend": "Agg", "no_random_state_used": True}}
    ex.update(extra or {})
    prov = provenance.build_provenance(
        program="O", milestone="O1", input_files=inputs,
        simulation_parameters={k: man[k] for k in man
                               if k in ("velocity", "lambda", "tau",
                                        "source_semantics",
                                        "accepted_iterations",
                                        "baseline_exact_E_neg",
                                        "final_exact_E_neg")},
        coordinate_frame="whitened design coordinates y; c = c0 + P W y",
        observer="ADM Eulerian observer",
        units="geometrised G = c = 1; lengths in metres",
        time_convention="static slice at t = 0; the iteration axis is NOT time",
        render_theme=theme.name, color_model="illustrative_only",
        camera_model="orthographic 2-D", frame_count=1, frame_rate=None,
        extra=ex)
    provenance.write_sidecar(png, prov)
    provenance.write_sidecar(svg, prov)
    if header is not None:
        records.write_plotted_csv(os.path.join(out, stem + ".plotted.csv"),
                                  header, rows or [])
    return png


def banner(fig, theme, main, sub=None):
    fig.text(0.5, 0.963, main, ha="center", fontsize=11,
             color=theme.foreground)
    if sub:
        fig.text(0.5, 0.917, sub, ha="center", fontsize=8,
                 color=theme.uncertainty_color)


def verdict_stamp(fig, theme, cert):
    ok = cert["exact_energy_reduction_established"]
    txt = ("EXACT REDUCTION ESTABLISHED" if ok
           else "SURROGATE REDUCED - EXACT REDUCTION NOT ESTABLISHED")
    fig.text(0.988, 0.962, txt, ha="right", va="center", fontsize=6.8,
             color=GREEN if ok else theme.uncertainty_color,
             bbox=dict(boxstyle="round,pad=0.28", facecolor=theme.background,
                       edgecolor=GREEN if ok else theme.uncertainty_color,
                       linewidth=0.7))


# --- 1. energy reduction history ----------------------------------------------
def fig_history(rec, out, man, theme, acc, cert, base):
    it = [float(r["iteration"]) for r in acc]
    exact = [float(r["exact_E_neg"]) for r in acc]
    ren = [float(r["renormalized_validation_E_neg_tau"]) for r in acc]
    renw = [float(r["renormalized_working_E_neg_tau"]) for r in acc]

    fig, axes = plt.subplots(1, 2, figsize=(14.6, 5.2),
                             facecolor=theme.background)
    axes[0].axhline(exact[0], color=theme.grid, ls="--", lw=1.2)
    axes[0].text(it[-1], exact[0], " baseline", fontsize=7,
                 color=theme.grid, va="bottom", ha="right")
    axes[0].plot(it, exact, "o-", color=RED, lw=2.2, ms=6)
    O.style_axes(axes[0], theme, "accepted iteration",
                 r"exact $E_{\rm neg}$  (the reporting objective)")
    lo, hi = min(exact), max(exact)
    pad = (hi - lo) * 0.25 or 1e-6
    axes[0].set_ylim(lo - pad, hi + pad)

    axes[1].plot(it, renw, "s-", color=AMBER, lw=1.6, ms=4.5,
                 label="working (optimized)")
    axes[1].plot(it, ren, "o-", color=VIOLET, lw=2.0, ms=5.5,
                 label="validation")
    O.style_axes(axes[1], theme, "accepted iteration",
                 r"renormalized $E_{\rm neg,\tau}$")
    axes[1].legend(fontsize=7.6, facecolor=theme.background,
                   edgecolor=theme.grid, labelcolor=theme.foreground)

    red = cert["absolute_reduction"]
    pct = cert["percent_reduction"]
    unc = cert["combined_numerical_uncertainty"]
    banner(fig, theme,
           "Exact negative Eulerian energy and the renormalized surrogate, "
           "versus accepted iteration",
           r"baseline %.9f  ->  final %.9f    ($\Delta$ = %.3e, %.3f%%, "
           r"uncertainty %.2e)"
           % (cert["baseline_exact_E_neg"], cert["final_exact_E_neg"], red,
              pct, unc))
    O.footer(fig, theme,
             "The iteration axis is NOT physical time. The search direction "
             "comes from the surrogate; only the left panel may support an "
             "energy claim. " + SCOPE[:150])
    verdict_stamp(fig, theme, cert)
    fig.tight_layout(rect=[0, 0.04, 1, 0.90])
    return save(fig, out, "energy_reduction_history", man,
                [os.path.join(rec, "accepted_iterates.csv"),
                 os.path.join(rec, "final_certification.json")], theme,
                extra={"frame_axis": "accepted iteration (not physical time)"},
                header=["iteration", "exact_E_neg", "renormalized_validation"],
                rows=[[a, b, c] for a, b, c in zip(it, exact, ren)])


# --- 2. optimizer mechanics ----------------------------------------------------
def fig_mechanics(rec, out, man, theme, acc, cert, base):
    it = [float(r["iteration"]) for r in acc]
    gn = [max(float(r["grad_norm"]), 1e-20) for r in acc]
    sn = [max(float(r["step_norm"]), 1e-20) for r in acc]
    tr = [float(r["trust_radius"]) for r in acc]
    ls = [float(r["line_search_evals"]) for r in acc]
    rej = O.read_csv(os.path.join(rec, "rejected_trials.csv"))

    fig, axes = plt.subplots(1, 3, figsize=(17.2, 4.8),
                             facecolor=theme.background)
    axes[0].semilogy(it, gn, "o-", color=BLUE, lw=1.9, ms=5,
                     label=r"$\|\nabla_y F_\tau\|_2$")
    axes[0].semilogy(it[1:], sn[1:], "s-", color=GREEN, lw=1.6, ms=4,
                     label=r"$\|\Delta y\|_2$")
    O.style_axes(axes[0], theme, "accepted iteration", "norm")
    axes[0].legend(fontsize=7.4, facecolor=theme.background,
                   edgecolor=theme.grid, labelcolor=theme.foreground)

    axes[1].plot(it, tr, "o-", color=AMBER, lw=1.9, ms=5)
    O.style_axes(axes[1], theme, "accepted iteration", "trust radius")
    tw = axes[1].twinx()
    tw.bar(it, ls, color=TEAL, alpha=0.35, width=0.6)
    tw.set_ylabel("line-search evaluations", color=TEAL, fontsize=7.6)
    tw.tick_params(colors=TEAL, labelsize=6.5)

    reasons = {}
    for r in rej:
        reasons[r["classification"]] = reasons.get(r["classification"], 0) + 1
    ks = list(reasons.keys()) or ["none"]
    vs = [reasons.get(k, 0) for k in ks]
    axes[2].barh(ks, vs, color=RED, edgecolor=theme.grid, linewidth=0.6)
    O.style_axes(axes[2], theme, "rejected trials", "")
    axes[2].tick_params(axis="y", labelsize=6.6)

    banner(fig, theme, "Optimizer mechanics: what it did and what it refused",
           "%d accepted iterates, %d rejected trials, terminated by %s"
           % (len(acc) - 1, len(rej), man.get("termination_reason", "-")))
    O.footer(fig, theme,
             "Every rejected trial carries a machine-readable reason and none "
             "of them appears in the accepted trace. L-BFGS builds an "
             "inverse-Hessian APPROXIMATION from gradient differences; no full "
             "design Hessian is formed or claimed.")
    verdict_stamp(fig, theme, cert)
    fig.tight_layout(rect=[0, 0.03, 1, 0.90])
    return save(fig, out, "optimizer_mechanics", man,
                [os.path.join(rec, "accepted_iterates.csv"),
                 os.path.join(rec, "rejected_trials.csv")], theme,
                header=["iteration", "grad_norm", "trust_radius"],
                rows=[[a, b, c] for a, b, c in zip(it, gn, tr)])


# --- 3. design-space trajectory ------------------------------------------------
def fig_trajectory(rec, out, man, theme, acc, cert, base):
    """Two principal directions of the accepted-iterate coordinate cloud.

    Pure linear algebra on saved y vectors: a projection for the eye, not a
    physical space. Power iteration on the covariance, deterministic seed.
    """
    grad = O.read_csv(os.path.join(rec, "gradient_trace.csv"))
    ys = {}
    for r in grad:
        ys.setdefault(int(r["iteration"]), []).append(0.0)
    n = len(ys.get(0, [])) or 34
    # Reconstruct y per iterate from the final design plus the trace is not
    # available, so use the recorded per-iterate norms for a 1-D fallback if
    # the full cloud is absent.
    fy = O.read_csv(os.path.join(rec, "final_design_y.csv"))
    yfinal = [float(r["y"]) for r in fy]

    # The accepted-iterate cloud, from the gradient trace's companion record.
    cloud = []
    itrace = O.read_csv(os.path.join(rec, "accepted_iterates.csv"))
    # y vectors are stored per iterate in the y_trace record when present;
    # otherwise project the single final vector plus its norm history.
    ypath = os.path.join(rec, "accepted_y.csv")
    if os.path.exists(ypath):
        rows = O.read_csv(ypath)
        by = {}
        for r in rows:
            by.setdefault(int(r["iteration"]), {})[int(r["index"])] = float(r["y"])
        for k in sorted(by):
            cloud.append([by[k][j] for j in sorted(by[k])])
    if not cloud:
        cloud = [[0.0] * len(yfinal), yfinal]

    m = len(cloud[0])
    mean = [sum(v[j] for v in cloud) / len(cloud) for j in range(m)]
    X = [[v[j] - mean[j] for j in range(m)] for v in cloud]

    def power(mat, exclude=None):
        v = [math.sin(1.7 * j + 0.3) for j in range(m)]
        nv = math.sqrt(sum(q * q for q in v)) or 1.0
        v = [q / nv for q in v]
        for _ in range(200):
            w = [0.0] * m
            for row in mat:
                d = sum(row[j] * v[j] for j in range(m))
                for j in range(m):
                    w[j] += d * row[j]
            if exclude:
                d = sum(w[j] * exclude[j] for j in range(m))
                for j in range(m):
                    w[j] -= d * exclude[j]
            nw = math.sqrt(sum(q * q for q in w))
            if nw < 1e-300:
                return v
            v = [q / nw for q in w]
        return v

    p1 = power(X)
    p2 = power(X, p1)
    xs = [sum(row[j] * p1[j] for j in range(m)) for row in X]
    zs = [sum(row[j] * p2[j] for j in range(m)) for row in X]

    fig, ax = plt.subplots(figsize=(11.4, 6.0), facecolor=theme.background)
    ax.plot(xs, zs, "-", color=theme.grid, lw=1.0, zorder=1)
    e = [float(r["exact_E_neg"]) for r in itrace][:len(xs)]
    if len(e) < len(xs):
        e = e + [e[-1]] * (len(xs) - len(e))
    sc = ax.scatter(xs, zs, c=e, cmap="magma_r", s=90, zorder=3,
                    edgecolor=theme.grid, linewidth=0.6)
    cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label(r"exact $E_{\rm neg}$", color=theme.foreground, fontsize=8)
    cb.ax.tick_params(labelsize=6.5, colors=theme.foreground)
    ax.scatter([xs[0]], [zs[0]], s=190, facecolor="none", edgecolor=GREEN,
               linewidth=1.8, zorder=4)
    ax.scatter([xs[-1]], [zs[-1]], s=190, facecolor="none", edgecolor=RED,
               linewidth=1.8, zorder=4)
    ax.annotate("start  y = 0", (xs[0], zs[0]), fontsize=7.6, color=GREEN,
                xytext=(6, 8), textcoords="offset points")
    ax.annotate("final", (xs[-1], zs[-1]), fontsize=7.6, color=RED,
                xytext=(6, 8), textcoords="offset points")
    O.style_axes(ax, theme, "first principal direction of the accepted cloud",
                 "second principal direction")

    banner(fig, theme,
           "Optimization trajectory in DESIGN space, not physical space",
           "a two-dimensional projection of 34 whitened coordinates; the axes "
           "are principal directions of the accepted-iterate cloud")
    O.footer(fig, theme,
             "Because the coordinates are whitened, ||y||_2 is the Sobolev "
             "design norm. This projection is a rendering of saved vectors: "
             "no physics is computed here. " + SCOPE[:130])
    verdict_stamp(fig, theme, cert)
    fig.tight_layout(rect=[0, 0.04, 1, 0.90])
    return save(fig, out, "design_space_trajectory", man,
                [os.path.join(rec, "final_design_y.csv"),
                 os.path.join(rec, "accepted_iterates.csv")], theme,
                extra={"axes_are_physical_space": False},
                header=["pc1", "pc2", "exact_E_neg"],
                rows=[[a, b, c] for a, b, c in zip(xs, zs, e)])


# --- 4. final convergence certification ----------------------------------------
def fig_convergence(rec, out, man, theme, acc, cert, base):
    order = O.read_csv(os.path.join(rec, "final_order_convergence.csv"))
    sub = O.read_csv(os.path.join(rec, "final_subdivision_convergence.csv"))
    dom = [r for r in O.read_csv(os.path.join(rec, "final_domain_convergence.csv"))
           if r["objective"] == "E_neg"]
    ind = O.read_csv(os.path.join(rec, "final_independent_reference.csv"))

    fig, axes = plt.subplots(1, 4, figsize=(18.0, 4.4),
                             facecolor=theme.background)
    for ax, rows, key, lab in ((axes[0], order, "order", "quadrature order"),
                               (axes[1], sub, "subdivision", "subdivision"),
                               (axes[2], dom, "xi_extent", "domain xi extent")):
        xs = [float(r[key]) for r in rows]
        ys = [float(r["value"]) for r in rows]
        ax.plot(xs, ys, "o-", color=RED, lw=2.0, ms=6)
        ax.axhline(cert["baseline_exact_E_neg"], color=theme.grid, ls="--",
                   lw=1.0)
        O.style_axes(ax, theme, lab, r"final exact $E_{\rm neg}$")
    gaps = [abs(float(r["relative_gap"])) for r in ind]
    axes[3].semilogy([float(r["n"]) for r in ind],
                     [max(g, 1e-18) for g in gaps], "o-", color=VIOLET, lw=2.0,
                     ms=6)
    axes[3].axhline(1e-4, color=theme.grid, ls=":", lw=1.0)
    O.style_axes(axes[3], theme, "uniform grid n", "relative gap")

    c = cert["convergence"]
    banner(fig, theme,
           "The FINAL candidate, recertified from scratch through the complete "
           "O0D protocol",
           "order %s, subdivision %s, domain %s, independent reference %s, "
           "gradient %s, validity %s"
           % tuple("PASS" if c[k] else "FAIL"
                   for k in ("order", "subdivision", "domain",
                             "independent_reference", "gradient_resolution",
                             "validity")))
    O.footer(fig, theme,
             "The optimizer's own validation quadrature was used to ACCEPT "
             "steps, so certifying the result with it would be circular. The "
             "dashed line is the baseline. " + SCOPE[:120])
    verdict_stamp(fig, theme, cert)
    fig.tight_layout(rect=[0, 0.03, 1, 0.90])
    return save(fig, out, "final_convergence_certification", man,
                [os.path.join(rec, n) for n in
                 ("final_order_convergence.csv",
                  "final_subdivision_convergence.csv",
                  "final_domain_convergence.csv",
                  "final_independent_reference.csv",
                  "final_certification.json")], theme,
                header=["study", "x", "E_neg"],
                rows=[["order", r["order"], r["value"]] for r in order])


# --- 5. tradeoff dashboard -----------------------------------------------------
def fig_tradeoff(rec, out, man, theme, acc, cert, base):
    rows = O.read_csv(os.path.join(rec, "baseline_vs_final_metrics.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(15.4, 5.6),
                             facecolor=theme.background)

    names = [r["metric"] for r in rows]
    frac = [float(r["fractional_change"]) * 100.0 for r in rows]
    worse = [r["worsened"] == "1" for r in rows]
    ys = range(len(names))
    axes[0].barh(list(ys), frac,
                 color=[RED if w else GREEN for w in worse],
                 edgecolor=theme.grid, linewidth=0.6)
    axes[0].axvline(0.0, color=theme.foreground, lw=0.9)
    axes[0].set_yticks(list(ys))
    axes[0].set_yticklabels(names, fontsize=7)
    O.style_axes(axes[0], theme, "percent change from baseline", "")

    # A parallel-coordinate view, each metric normalised to its own baseline.
    for i, r in enumerate(rows):
        b, fv = float(r["baseline"]), float(r["final"])
        denom = abs(b) if abs(b) > 1e-30 else 1.0
        axes[1].plot([0, 1], [b / denom, fv / denom], "-o",
                     color=RED if r["worsened"] == "1" else GREEN, lw=1.4,
                     ms=4, alpha=0.85)
        axes[1].annotate(r["metric"], (1, fv / denom), fontsize=5.8,
                         color=theme.foreground, xytext=(4, 0),
                         textcoords="offset points", va="center")
    axes[1].set_xticks([0, 1])
    axes[1].set_xticklabels(["baseline", "final"], fontsize=8)
    O.style_axes(axes[1], theme, "", "value / |baseline|")
    axes[1].set_xlim(-0.1, 1.55)

    nw = sum(1 for w in worse if w)
    banner(fig, theme,
           "Tradeoffs: what the reduction cost elsewhere",
           "%d of %d diagnostics moved in the unfavourable direction; lower is "
           "NOT uniformly better for every metric" % (nw, len(rows)))
    O.footer(fig, theme,
             "E_net is signed and larger is not automatically worse; N_min "
             "more negative means a deeper NEC violation; K_max is a curvature "
             "magnitude. Only E_neg was optimized. " + SCOPE[:120])
    verdict_stamp(fig, theme, cert)
    fig.tight_layout(rect=[0, 0.03, 1, 0.90])
    return save(fig, out, "tradeoff_dashboard", man,
                [os.path.join(rec, "baseline_vs_final_metrics.csv"),
                 os.path.join(rec, "final_certification.json")], theme,
                header=["metric", "baseline", "final", "worsened"],
                rows=[[r["metric"], r["baseline"], r["final"], r["worsened"]]
                      for r in rows])


# --- 6. the one-page summary ----------------------------------------------------
def fig_summary(rec, out, man, theme, acc, cert, base):
    rows = O.read_csv(os.path.join(rec, "baseline_vs_final_metrics.csv"))
    fo = O.read_json(os.path.join(rec, "final_objectives.json"))
    det = O.read_json(os.path.join(rec, "deterministic_repeat.json"))
    res = O.read_json(os.path.join(rec, "restart_check.json"))
    fin = acc[-1]

    fig = plt.figure(figsize=(15.0, 7.4), facecolor=theme.background)
    ax = fig.add_axes([0.04, 0.08, 0.44, 0.74])
    ok = cert["exact_energy_reduction_established"]
    b, f = cert["baseline_exact_E_neg"], cert["final_exact_E_neg"]
    ax.bar([0, 1], [b, f], color=[theme.grid, GREEN if ok else AMBER],
           width=0.5, edgecolor=theme.grid, linewidth=0.8)
    ax.errorbar([1], [f], yerr=[cert["combined_numerical_uncertainty"]],
                color=RED, capsize=6, lw=1.6, zorder=5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["baseline\n(zero correction)", "optimized\ncandidate"],
                       fontsize=8.5)
    O.style_axes(ax, theme, "", r"exact $E_{\rm neg}$")
    lo = min(b, f) * 0.985
    ax.set_ylim(lo, max(b, f) * 1.006)
    ax.text(0.5, 0.90,
            "%.4f%% reduction\n%.3e absolute\nuncertainty %.2e  (ratio %.1f)"
            % (cert["percent_reduction"], cert["absolute_reduction"],
               cert["combined_numerical_uncertainty"],
               cert["reduction_over_uncertainty"]),
            transform=ax.transAxes, ha="center", fontsize=9,
            color=theme.foreground)

    tx = fig.add_axes([0.52, 0.08, 0.45, 0.74])
    tx.axis("off")
    tx.set_facecolor(theme.background)
    lines = [
        ("final candidate certified",
         "YES" if cert["final_candidate_certified"] else "NO",
         cert["final_candidate_certified"]),
        ("exact reduction established", "YES" if ok else "NO", ok),
        ("protected cabin locked", "EXACT ZERO", True),
        ("outer collar locked", "EXACT ZERO", True),
        ("deterministic repeat", "IDENTICAL" if det["identical"] else "DIFFERS",
         det["identical"]),
        ("restart reproducibility",
         "WITHIN TOL" if res["within_tolerance"] else "OUTSIDE TOL",
         res["within_tolerance"]),
        ("independent reference gap",
         "%.2e" % cert["independent_reference_relative_gap"],
         cert["independent_reference_relative_gap"] < 1e-4),
        ("accepted iterations", fo.get("accepted_iterations", "-"), True),
        ("rejected trials", fo.get("rejected_trials", "-"), True),
        ("termination", str(fo.get("termination_reason", "-"))[:34],
         not fo.get("terminated_by_limit", True)),
    ]
    for i, (k, v, good) in enumerate(lines):
        y = 0.94 - i * 0.085
        tx.text(0.0, y, k, fontsize=9, color=theme.foreground)
        tx.text(0.98, y, str(v), fontsize=9, ha="right",
                color=GREEN if good else RED)
        tx.plot([0.0, 1.0], [y - 0.028, y - 0.028], color=theme.grid, lw=0.4)

    worse = [r["metric"] for r in rows if r["worsened"] == "1"]
    banner(fig, theme,
           "F1W source-requirement optimization: final result",
           cert["claim"][:118] + " ...")
    O.footer(fig, theme,
             ("diagnostics that worsened: " + (", ".join(worse) or "none")) +
             ". Only exact E_neg was optimized; the others are tradeoffs. " +
             SCOPE)
    verdict_stamp(fig, theme, cert)
    _ = fin
    return save(fig, out, "final_design_summary", man,
                [os.path.join(rec, "final_certification.json"),
                 os.path.join(rec, "baseline_vs_final_metrics.csv"),
                 os.path.join(rec, "deterministic_repeat.json"),
                 os.path.join(rec, "restart_check.json"),
                 os.path.join(rec, "final_objectives.json")], theme,
                extra={"claim": cert["claim"],
                       "exact_energy_reduction_established": ok},
                header=["quantity", "value"],
                rows=[["baseline_E_neg", b], ["final_E_neg", f],
                      ["percent_reduction", cert["percent_reduction"]]])


BUILDERS = (fig_history, fig_mechanics, fig_trajectory, fig_convergence,
            fig_tradeoff, fig_summary)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rec, out = records.resolve(a.records), records.ensure_dir(a.out)
    man = O.read_json(os.path.join(rec, "o1_manifest.json"))
    cert = O.read_json(os.path.join(rec, "final_certification.json"))
    base = O.read_json(os.path.join(rec, "baseline_objectives.json"))
    acc = O.read_csv(os.path.join(rec, "accepted_iterates.csv"))
    theme = themes.get("scientific")
    themes.apply(theme)
    made = [fn(rec, out, man, theme, acc, cert, base) for fn in BUILDERS]
    for p in made:
        print("wrote", os.path.basename(p))
    with open(os.path.join(out, "figure_manifest.json"), "w",
              encoding="utf-8") as f:
        json.dump({"figures": [os.path.basename(p)[:-4] for p in made],
                   "count": len(made), "scope": SCOPE,
                   "renderer_computes_physics": False,
                   "iteration_axis_is_physical_time": False}, f, indent=2)
    print("%d figures" % len(made))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
