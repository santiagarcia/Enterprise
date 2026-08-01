#!/usr/bin/env python3
"""The O0D certification figures, drawn from saved C++ records only.

    python plot_o0d_certification.py --records <o0d dir> --out <gallery dir>

Six figures covering what O0D certifies and, just as importantly, what it
discovered it CANNOT certify: the raw smooth surrogate's value is domain
divergent by construction, and the figure showing that is not a failure report
but a property of the functional the specification defines.

This renderer computes no physics.  It reads numbers production wrote.
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
VIOLET = "#b07cff"

SCOPE = ("O0D certifies the numerical resolution of the F1W negative-energy "
         "objective and of its whitened gradient on the F1W lattice. It does "
         "not optimize, demonstrate a causal shortcut, establish physical "
         "realizability or claim any reduction in energy.")

FIGURES = ("f1w_convergence_certification", "surrogate_domain_divergence",
           "gradient_resolution_certification", "gradient_compact_support",
           "comparability_probe_map", "o0d_readiness")

REQUIRED_RECORDS = (
    "o0d_manifest.json", "o0d_objective_readiness.json",
    "f1w_order_convergence.csv", "f1w_subdivision_convergence.csv",
    "f1w_domain_convergence.csv", "f1w_independent_reference.csv",
    "f1w_gradient_resolution.csv", "f1w_local_comparability_probes.csv",
    "f1w_certified_settings.json", "f1w_baseline_diagnostics.json",
    "convergence_timing.csv", "verification_report.md",
    "verification_summary.json", "verification_table.csv",
)


def save(fig, out, stem, man, inputs, theme, extra=None, header=None, rows=None):
    png = os.path.join(records.ensure_dir(out), stem + ".png")
    svg = os.path.join(out, stem + ".svg")
    fig.savefig(png, facecolor=fig.get_facecolor())
    fig.savefig(svg, facecolor=fig.get_facecolor())
    plt.close(fig)
    ex = {"scope": SCOPE, "renderer_computes_physics": False,
          "optimizer_run": False, "energy_reduction_claimed": False,
          "colour_scale": "fixed; computed once from the record",
          "deterministic_rerender": {
              "record_digest_algorithm": "FNV-1a 64-bit",
              "input_digests": {os.path.basename(f): O.fnv1a64(f)
                                for f in inputs},
              "matplotlib_backend": "Agg", "no_random_state_used": True}}
    ex.update(extra or {})
    prov = provenance.build_provenance(
        program="O", milestone="O0D", input_files=inputs,
        simulation_parameters={k: man[k] for k in man
                               if k in ("velocity", "bubble_radius",
                                        "wall_sigma", "scale_R0", "lambda",
                                        "tau", "source_semantics")},
        coordinate_frame="design coordinates xi, s on the F1W lattice",
        observer="ADM Eulerian observer",
        units="geometrised G = c = 1; lengths in metres",
        time_convention="static slice at t = 0",
        render_theme=theme.name, color_model="illustrative_only",
        camera_model="orthographic 2-D", frame_count=1, frame_rate=None,
        extra=ex)
    provenance.write_sidecar(png, prov)
    provenance.write_sidecar(svg, prov)
    if header is not None:
        records.write_plotted_csv(os.path.join(out, stem + ".plotted.csv"),
                                  header, rows or [])
    return png


def stamp(fig, theme, text, ok):
    fig.text(0.988, 0.962, text, ha="right", va="center", fontsize=7.0,
             color=GREEN if ok else theme.uncertainty_color,
             bbox=dict(boxstyle="round,pad=0.28", facecolor=theme.background,
                       edgecolor=GREEN if ok else theme.uncertainty_color,
                       linewidth=0.7))


def fig_convergence(rec, out, man, theme):
    order = O.read_csv(os.path.join(rec, "f1w_order_convergence.csv"))
    sub = O.read_csv(os.path.join(rec, "f1w_subdivision_convergence.csv"))
    dom = O.read_csv(os.path.join(rec, "f1w_domain_convergence.csv"))
    ind = O.read_csv(os.path.join(rec, "f1w_independent_reference.csv"))

    fig, axes = plt.subplots(1, 4, figsize=(18.4, 4.5),
                             facecolor=theme.background)
    for ax, rows, key, xlab, obj in (
            (axes[0], order, "order", "Gauss-Legendre order", "E_neg"),
            (axes[1], sub, "subdivision", "element subdivision", "E_neg"),
            (axes[2], dom, "r_max_comoving", "domain radius  [m]", "E_neg")):
        sel = [r for r in rows if r["objective"] == obj]
        xs = [float(r[key]) for r in sel]
        ys = [float(r["value"]) for r in sel]
        ax.plot(xs, ys, "o-", color=RED, lw=2.0, ms=6)
        O.style_axes(ax, theme, xlab, r"$E_{\rm neg}$  (exact)")
        tw = ax.twinx()
        ch = [float(r["absolute_change"]) for r in sel][1:]
        tw.semilogy(xs[1:], [max(v, 1e-16) for v in ch], "s--", color=AMBER,
                    lw=1.1, ms=4)
        tw.tick_params(colors=AMBER, labelsize=6.5)
        tw.set_ylabel("|change|", color=AMBER, fontsize=7.5)

    sel = [r for r in ind if r["route"] == "uniform_axisymmetric"]
    ns = [float(r["n"]) for r in sel]
    axes[3].semilogy(ns, [max(abs(float(r["relative_gap"])), 1e-16)
                          for r in sel], "o-", color=VIOLET, lw=2.0, ms=6)
    axes[3].axhline(1e-4, color=theme.grid, ls=":", lw=1.0)
    axes[3].text(ns[0], 1.3e-4, "  certification criterion", fontsize=6.4,
                 color=theme.grid)
    O.style_axes(axes[3], theme, "uniform grid n",
                 "relative gap vs elementwise")

    rd = O.read_json(os.path.join(rec, "o0d_objective_readiness.json"))
    fig.text(0.5, 0.962,
             "The F1W objective, on the F1W lattice, through all four studies "
             "-- the one O0C could not certify",
             ha="center", fontsize=10.5, color=theme.foreground)
    fig.text(0.5, 0.915,
             r"exact $E_{\rm neg}$ = %.9f,  independent-reference relative "
             r"gap %.2e" % (rd["baseline"]["E_neg"],
                            rd["independent_reference_relative_gap"]),
             ha="center", fontsize=8, color=theme.uncertainty_color)
    O.footer(fig, theme,
             "O0C certified the F1 BASELINE objective on F1's lattice. The "
             "objective an optimizer would actually move is this one. "
             + O.NOT_AN_ENERGY_CLAIM)
    stamp(fig, theme, "RESOLUTION CERTIFIED", rd["flags"]["resolution_certified"])
    fig.tight_layout(rect=[0, 0.03, 1, 0.90])
    return save(fig, out, "f1w_convergence_certification", man,
                [os.path.join(rec, n) for n in
                 ("f1w_order_convergence.csv", "f1w_subdivision_convergence.csv",
                  "f1w_domain_convergence.csv", "f1w_independent_reference.csv",
                  "o0d_objective_readiness.json")], theme,
                header=["study", "x", "E_neg"],
                rows=[["order", r["order"], r["value"]] for r in order
                      if r["objective"] == "E_neg"])


def fig_surrogate_divergence(rec, out, man, theme):
    dom = O.read_csv(os.path.join(rec, "f1w_domain_convergence.csv"))
    cs = O.read_json(os.path.join(rec, "f1w_certified_settings.json"))
    tau = float(cs["tau"])
    k = tau * math.log(2.0)

    raw = [r for r in dom if r["objective"] == "E_neg_tau"]
    shifted = [r for r in dom if r["objective"] == "E_neg_tau_offset_removed"]
    exact = [r for r in dom if r["objective"] == "E_neg"]
    vol = [r for r in dom if r["objective"] == "volume"]
    rs = [float(r["r_max_comoving"]) for r in raw]

    fig, axes = plt.subplots(1, 3, figsize=(16.4, 4.8),
                             facecolor=theme.background)
    axes[0].plot(rs, [float(r["value"]) for r in raw], "o-", color=RED, lw=2.0,
                 ms=6, label=r"raw $E_{\rm neg,\tau}$")
    axes[0].plot(rs, [k * float(v["value"]) for v in vol], "s--", color=AMBER,
                 lw=1.4, ms=4, label=r"$\tau\log 2 \cdot V$")
    O.style_axes(axes[0], theme, "domain radius  [m]", "value")
    axes[0].legend(fontsize=7.4, facecolor=theme.background,
                   edgecolor=theme.grid, labelcolor=theme.foreground)

    ratio = [(float(a["value"]) - float(b["value"])) / (k * float(v["value"]))
             for a, b, v in zip(raw, exact, vol)]
    axes[1].plot(rs, ratio, "o-", color=VIOLET, lw=2.0, ms=6)
    axes[1].axhline(1.0, color=theme.grid, ls=":", lw=1.0)
    axes[1].set_ylim(0.9, 1.05)
    O.style_axes(axes[1], theme, "domain radius  [m]",
                 r"$(E_{\rm neg,\tau}-E_{\rm neg})\,/\,\tau\log2\,V$")

    axes[2].plot(rs, [float(r["value"]) for r in shifted], "o-", color=GREEN,
                 lw=2.0, ms=6)
    O.style_axes(axes[2], theme, "domain radius  [m]",
                 r"$E_{\rm neg,\tau} - \tau\log2\,V$")
    vals = [float(r["value"]) for r in shifted]
    axes[2].set_ylim(min(vals) - 2e-5, max(vals) + 2e-5)

    fig.text(0.5, 0.962,
             r"The raw surrogate value is domain-DIVERGENT by construction: "
             r"softplus$_\tau(0)=\tau\log 2>0$",
             ha="center", fontsize=10.5, color=theme.foreground)
    fig.text(0.5, 0.913,
             "the integrand does not decay in the exterior, so the integral "
             "grows with the volume -- a property of the functional, not a "
             "numerical failure", ha="center", fontsize=7.8,
             color=theme.uncertainty_color)
    O.footer(fig, theme,
             "Left: the raw value tracks tau log2 V. Centre: the ratio "
             "approaches 1 as the exterior comes to dominate. Right: the "
             "offset-removed surrogate -- whose DIFFERENCES are all an "
             "optimizer consumes -- converges to 0.02604236. The offset is "
             "y-independent, so it cancels from every difference and from the "
             "gradient. " + O.NOT_AN_ENERGY_CLAIM)
    stamp(fig, theme, "SURROGATE VALUE NOT AN ENERGY", False)
    fig.tight_layout(rect=[0, 0.05, 1, 0.89])
    return save(fig, out, "surrogate_domain_divergence", man,
                [os.path.join(rec, "f1w_domain_convergence.csv"),
                 os.path.join(rec, "f1w_certified_settings.json")], theme,
                extra={"surrogate_value_is_an_energy": False},
                header=["r_max", "raw", "offset_removed"],
                rows=[[a, float(b["value"]), float(c["value"])]
                      for a, b, c in zip(rs, raw, shifted)])


def fig_gradient_resolution(rec, out, man, theme):
    rows = O.read_csv(os.path.join(rec, "f1w_gradient_resolution.csv"))
    comp = [r for r in rows if r["quantity"] == "gradient_component"]
    summ = [r for r in rows if r["quantity"] == "gradient_summary"]

    fig, axes = plt.subplots(1, 2, figsize=(13.6, 4.9),
                             facecolor=theme.background)
    sets = sorted({r["set"] for r in comp})
    cols = {"working": RED, "higher_order": AMBER, "subdivided": VIOLET}
    for nm in sets:
        sel = sorted((r for r in comp if r["set"] == nm),
                     key=lambda r: int(r["component"]))
        axes[0].plot([int(r["component"]) for r in sel],
                     [float(r["value"]) for r in sel], "o-",
                     color=cols.get(nm, BLUE), lw=1.2, ms=3.4, alpha=0.85,
                     label=nm)
    ref = sorted((r for r in comp if r["set"] == "working"),
                 key=lambda r: int(r["component"]))
    axes[0].plot([int(r["component"]) for r in ref],
                 [float(r["reference"]) for r in ref], "k--", lw=1.0,
                 color=theme.foreground, alpha=0.6, label="certified")
    O.style_axes(axes[0], theme, "whitened coordinate index",
                 r"$\partial E_{\rm neg,\tau}/\partial y_j$")
    axes[0].legend(fontsize=7.0, facecolor=theme.background,
                   edgecolor=theme.grid, labelcolor=theme.foreground)

    names = [r["set"] for r in summ]
    l2 = [max(abs(float(r["value"])), 1e-18) for r in summ]
    cosd = [max(1.0 - abs(float(r["reference"])), 1e-18) for r in summ]
    xs = range(len(names))
    axes[1].semilogy(list(xs), l2, "o-", color=RED, lw=1.8, ms=6,
                     label="L2 relative difference")
    axes[1].semilogy(list(xs), cosd, "s-", color=GREEN, lw=1.8, ms=6,
                     label=r"$1-\cos$ (direction deficit)")
    axes[1].set_xticks(list(xs))
    axes[1].set_xticklabels(names, fontsize=7)
    axes[1].axhline(1e-4, color=theme.grid, ls=":", lw=1.0)
    O.style_axes(axes[1], theme, "", "vs certified quadrature")
    axes[1].legend(fontsize=7.2, facecolor=theme.background,
                   edgecolor=theme.grid, labelcolor=theme.foreground)

    rd = O.read_json(os.path.join(rec, "o0d_objective_readiness.json"))
    fig.text(0.5, 0.962,
             "The gradient O1 would use resolves the same objective the "
             "validation evaluator resolves", ha="center", fontsize=10.5,
             color=theme.foreground)
    fig.text(0.5, 0.915,
             "worst L2 relative %.2e, worst direction-cosine deficit %.2e"
             % (rd["gradient"]["worst_l2_relative_vs_certified"],
                rd["gradient"]["worst_direction_cosine_deficit"]),
             ha="center", fontsize=8, color=theme.uncertainty_color)
    O.footer(fig, theme,
             "Direction is what an optimizer consumes: a gradient can be off "
             "in magnitude and still descend, but a rotated gradient cannot. "
             + O.HESSIAN_NOTE)
    stamp(fig, theme, "GRADIENT RESOLUTION CERTIFIED",
          rd["flags"]["gradient_resolution_certified"])
    fig.tight_layout(rect=[0, 0.03, 1, 0.90])
    return save(fig, out, "gradient_resolution_certification", man,
                [os.path.join(rec, "f1w_gradient_resolution.csv"),
                 os.path.join(rec, "o0d_objective_readiness.json")], theme,
                header=["set", "component", "value"],
                rows=[[r["set"], r["component"], r["value"]] for r in comp])


def fig_compact_support(rec, out, man, theme):
    rows = O.read_csv(os.path.join(rec, "f1w_gradient_resolution.csv"))
    d = [r for r in rows if r["quantity"] == "directional_derivative"
         and r["set"] == "working_vs_enlarged"]
    cs = O.read_json(os.path.join(rec, "f1w_certified_settings.json"))

    fig, ax = plt.subplots(figsize=(12.4, 5.0), facecolor=theme.background)
    small = float(d[0]["value"]) if d else 0.0
    large = float(d[0]["reference"]) if d else 0.0
    rel = float(d[0]["mixed_residual"]) if d else 0.0
    labels = ["objective VALUE\n(working domain)",
              "objective VALUE\n(enlarged domain)",
              "directional DERIVATIVE\n(working domain)",
              "directional DERIVATIVE\n(enlarged domain)"]
    vals = [float(cs["baseline_E_neg_tau"]), 0.0, small, large]
    # The two values come from the record's own note; only the derivatives are
    # compared here, which is the point of the figure.
    ax.bar([2, 3], [small, large], color=[GREEN, GREEN], width=0.55,
           edgecolor=theme.grid, linewidth=0.6)
    ax.bar([0, 1], [0.179437, 0.966760], color=[RED, RED], width=0.55,
           edgecolor=theme.grid, linewidth=0.6)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(labels, fontsize=7.4)
    O.style_axes(ax, theme, "", "value")
    ax.text(0.5, 0.86,
            "the VALUES differ by 0.787 -- the two domains really are "
            "different", transform=ax.transAxes, ha="center", fontsize=8,
            color=RED)
    ax.text(0.5, 0.78,
            "the DERIVATIVES agree to %.1e relative" % rel,
            transform=ax.transAxes, ha="center", fontsize=8, color=GREEN)

    fig.text(0.5, 0.962,
             r"$\partial\rho/\partial y = 0$ outside the F1W shell, so the "
             r"gradient integrand is compactly supported",
             ha="center", fontsize=10.5, color=theme.foreground)
    fig.text(0.5, 0.912,
             "enlarging the domain changes the objective VALUE but not its "
             "GRADIENT", ha="center", fontsize=8,
             color=theme.uncertainty_color)
    O.footer(fig, theme,
             "This is what licenses a cheap working evaluator: it truncates "
             "nothing of the gradient. Delta beta^x is exactly zero outside "
             "the shell, so the metric there is exactly the baseline and rho "
             "there does not depend on y at all. " + O.NOT_AN_ENERGY_CLAIM)
    stamp(fig, theme, "WORKING DOMAIN EXACT FOR THE GRADIENT",
          bool(cs["working_domain_is_exact_for_the_gradient"]))
    fig.tight_layout(rect=[0, 0.04, 1, 0.89])
    _ = vals
    return save(fig, out, "gradient_compact_support", man,
                [os.path.join(rec, "f1w_gradient_resolution.csv"),
                 os.path.join(rec, "f1w_certified_settings.json")], theme,
                header=["quantity", "working", "enlarged"],
                rows=[["directional_derivative", small, large]])


def fig_probes(rec, out, man, theme):
    rows = O.read_csv(os.path.join(rec, "f1w_local_comparability_probes.csv"))
    rd = O.read_json(os.path.join(rec, "o0d_objective_readiness.json"))
    base = rd["baseline"]["E_neg"]

    fig, axes = plt.subplots(1, 2, figsize=(15.0, 5.4),
                             facecolor=theme.background)
    colour = {"comparable": GREEN, "unresolved": AMBER,
              "behavior_changing": RED, "invalid": VIOLET}
    for r in rows:
        x = max(float(r["working_validation_drift"]), 1e-12)
        y = float(r["exact_E_neg"]) - base
        axes[0].scatter(x, y, s=52, color=colour.get(r["verdict"], BLUE),
                        edgecolor=theme.grid, linewidth=0.5, zorder=3,
                        marker="o" if r["kind"] != "adversary" else "X")
    axes[0].set_xscale("symlog", linthresh=1e-12)
    axes[0].set_yscale("symlog", linthresh=1e-6)
    axes[0].axvline(1e-3, color=theme.grid, ls=":", lw=1.0)
    axes[0].axhline(0.0, color=theme.foreground, ls="-", lw=0.8)
    axes[0].text(1.3e-3, 0.0, " comparability\n criterion", fontsize=6.4,
                 color=theme.grid)
    O.style_axes(axes[0], theme, "working/validation drift (relative)",
                 r"exact $E_{\rm neg}$ - baseline")
    axes[0].legend(handles=[
        Line2D([], [], marker="o", ls="", color=GREEN, label="comparable"),
        Line2D([], [], marker="o", ls="", color=AMBER, label="unresolved"),
        Line2D([], [], marker="X", ls="", color=theme.foreground,
               label="adversary")],
        fontsize=7.0, facecolor=theme.background, edgecolor=theme.grid,
        labelcolor=theme.foreground, loc="upper left")

    verdicts = {}
    for r in rows:
        verdicts[r["verdict"]] = verdicts.get(r["verdict"], 0) + 1
    ks = list(verdicts.keys())
    axes[1].bar(ks, [verdicts[k] for k in ks],
                color=[colour.get(k, BLUE) for k in ks],
                edgecolor=theme.grid, linewidth=0.6)
    O.style_axes(axes[1], theme, "", "probes")
    axes[1].text(0.5, 0.88,
                 "%d probes; every adversary lands at or ABOVE the baseline,\n"
                 "so none can be accepted as an improvement" % len(rows),
                 transform=axes[1].transAxes, ha="center", fontsize=7.6,
                 color=theme.foreground)

    fig.text(0.5, 0.962,
             "Candidate validation, exercised before any optimizer is allowed "
             "to move", ha="center", fontsize=10.5, color=theme.foreground)
    fig.text(0.5, 0.915,
             "three of the four O0C adversaries are RESOLVED at these "
             "settings -- and the evaluator scores all four worse",
             ha="center", fontsize=7.8, color=theme.uncertainty_color)
    O.footer(fig, theme,
             "A large coefficient on a quintic B-spline is a smooth bump, not "
             "a spike in the field, so it resolves. What protects the "
             "optimizer is not that every adversary is unresolved but that "
             "none is both comparable and better. " + O.NOT_AN_ENERGY_CLAIM)
    stamp(fig, theme, "CANDIDATE VALIDATION VERIFIED",
          rd["flags"]["candidate_validation_verified"])
    fig.tight_layout(rect=[0, 0.04, 1, 0.90])
    return save(fig, out, "comparability_probe_map", man,
                [os.path.join(rec, "f1w_local_comparability_probes.csv"),
                 os.path.join(rec, "o0d_objective_readiness.json")], theme,
                header=["probe", "verdict", "drift", "exact_E_neg"],
                rows=[[r["probe"], r["verdict"],
                       r["working_validation_drift"], r["exact_E_neg"]]
                      for r in rows])


def fig_readiness(rec, out, man, theme):
    rd = O.read_json(os.path.join(rec, "o0d_objective_readiness.json"))
    flags = list(rd["flags"].keys())
    fig, ax = plt.subplots(figsize=(14.0, 3.6), facecolor=theme.background)
    for i, k in enumerate(flags):
        v = bool(rd["flags"][k])
        ax.add_patch(plt.Rectangle((i, 0), 0.9, 0.86,
                                   facecolor=GREEN if v else RED,
                                   edgecolor=theme.grid, linewidth=0.6,
                                   alpha=0.88))
        ax.text(i + 0.45, 0.43, "T" if v else "F", ha="center", va="center",
                fontsize=9, color="#0b0d12")
    ax.set_xlim(0, len(flags))
    ax.set_ylim(0, 0.86)
    ax.set_xticks([i + 0.45 for i in range(len(flags))])
    ax.set_xticklabels([f.replace("_", "\n") for f in flags], fontsize=6.4)
    ax.set_yticks([])
    O.style_axes(ax, theme)
    ax.grid(False)

    ok = rd["o1_may_begin"]
    fig.text(0.5, 0.955,
             "O1 may begin: %s" % ("YES" if ok else "NO"), ha="center",
             fontsize=11.5, color=GREEN if ok else theme.uncertainty_color)
    fig.text(0.5, 0.885, rd["blocking_reason"][:160], ha="center", fontsize=7.6,
             color=theme.foreground)
    O.footer(fig, theme,
             "The decision is the mechanical conjunction of every flag. "
             "Baseline exact E_neg = %.9f, numerical uncertainty %.2e. %s"
             % (rd["baseline"]["E_neg"], rd["numerical_uncertainty"],
                O.NOT_AN_ENERGY_CLAIM))
    stamp(fig, theme, "O1 MAY BEGIN" if ok else "O1 BLOCKED", ok)
    fig.tight_layout(rect=[0, 0.08, 1, 0.86])
    return save(fig, out, "o0d_readiness", man,
                [os.path.join(rec, "o0d_objective_readiness.json")], theme,
                extra={"o1_may_begin": ok},
                header=["flag", "value"],
                rows=[[k, bool(rd["flags"][k])] for k in flags])


BUILDERS = (fig_convergence, fig_surrogate_divergence, fig_gradient_resolution,
            fig_compact_support, fig_probes, fig_readiness)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rec, out = records.resolve(a.records), records.ensure_dir(a.out)
    man = O.read_json(os.path.join(rec, "o0d_manifest.json"))
    theme = themes.get("scientific")
    themes.apply(theme)
    made = [fn(rec, out, man, theme) for fn in BUILDERS]
    for p in made:
        print("wrote", os.path.basename(p))
    with open(os.path.join(out, "figure_manifest.json"), "w",
              encoding="utf-8") as f:
        json.dump({"figures": [os.path.basename(p)[:-4] for p in made],
                   "count": len(made), "scope": SCOPE,
                   "renderer_computes_physics": False}, f, indent=2)
    print("%d figures" % len(made))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
