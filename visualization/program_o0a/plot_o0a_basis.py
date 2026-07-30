#!/usr/bin/env python3
"""The twelve O0A spline-basis figures, drawn from saved C++ records only.

    python plot_o0a_basis.py --records <o0a_basis dir> --out <gallery dir>

Nothing here evaluates a spline, a constraint, a Gram matrix or a projection. See
o0a_common for the contract and for why the tensor-product atlas is read from a
record rather than formed from two univariate columns.
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

import o0a_common as O                                           # noqa: E402
from common import anomaly, records                              # noqa: E402
from themes import themes                                        # noqa: E402

ACTIVE_COLOR = "#4ea1ff"
REMOVED_COLOR = "#8c8f99"
PAD_COLOR = "#4a4f5c"


# --- helpers ---------------------------------------------------------------

def _basis_panels(rec):
    """The three bases, each with its own knot semantics. Kept as a list so a
    figure cannot accidentally imply they share a boundary convention."""
    prof = O.read_csv(os.path.join(rec, "axisymmetric_basis_profiles.csv"))
    f0v = O.read_csv(os.path.join(rec, "univariate_basis_profiles.csv"))
    f0d = O.read_csv(os.path.join(rec, "univariate_basis_derivatives.csv"))
    f0c = {int(r["index"]): r for r in
           O.read_csv(os.path.join(rec, "f0_coefficients.csv"))}
    axial = [r for r in prof if r["basis"] == "axial_xi"]
    radial = [r for r in prof if r["basis"] == "radial_s"]

    # F0's profiles carry no status column, so the status comes from
    # f0_coefficients.csv, which records the exclusion reason per mode.
    f0rows = []
    dmap = {}
    for r in f0d:
        dmap[(r["s"], r["index"])] = r
    for r in f0v:
        k = int(r["index"])
        d = dmap.get((r["s"], r["index"]), {})
        f0rows.append({
            "coordinate": r["s"], "index": r["index"], "value": r["value"],
            "first": d.get("first", "0"), "second": d.get("second", "0"),
            "third": d.get("third", "0"),
            "exclusion_reason": f0c.get(k, {}).get("status", "active"),
            "support_lower": f0c.get(k, {}).get("support_lower", "0"),
            "support_upper": f0c.get(k, {}).get("support_upper", "0"),
            "semantics": "clamped",
        })
    return [
        ("F0 radial validation", "clamped", r"$s_r=(r/R_0)^2$", f0rows),
        ("F1 axial correction", "uniform", r"$\xi=(x-x_c)/R_0$", axial),
        ("F1 radial correction", "lower_clamped_open_upper",
         r"$s=(y^2+z^2)/R_0^2$", radial),
    ]


def _plot_modes(ax, rows, key, theme, *, show_support=True):
    by = O.group_by(rows, "index", int)
    for idx in sorted(by):
        pts = sorted(by[idx], key=lambda r: float(r["coordinate"]))
        reason = pts[0]["exclusion_reason"]
        active = reason == "active"
        col = ACTIVE_COLOR if active else (
            PAD_COLOR if reason == "outside_design_region" else REMOVED_COLOR)
        ax.plot(O.col(pts, "coordinate"), O.col(pts, key), lw=1.5 if active else 0.9,
                color=col, alpha=1.0 if active else 0.55,
                ls="-" if active else "--")
        if show_support and active:
            lo, hi = float(pts[0]["support_lower"]), float(pts[0]["support_upper"])
            ax.plot([lo, hi], [0, 0], lw=3.0, color=col, alpha=0.35,
                    solid_capstyle="butt")
    style_legend(ax, theme)


def style_legend(ax, theme):
    handles = [
        Line2D([], [], color=ACTIVE_COLOR, lw=1.6, label="active design mode"),
        Line2D([], [], color=REMOVED_COLOR, lw=1.0, ls="--",
               label="removed (collar / cabin)"),
        Line2D([], [], color=ACTIVE_COLOR, lw=3.0, alpha=0.35,
               label="local support of an active mode"),
    ]
    ax.legend(handles=handles, fontsize=6.4, frameon=False,
              labelcolor=theme.foreground, loc="upper right")


def _derivative_figure(rec, out, manifest, theme, key, stem, title, note):
    panels = _basis_panels(rec)
    fig, axes = plt.subplots(1, 3, figsize=(16.4, 4.9),
                             facecolor=theme.background)
    for ax, (name, sem, xlabel, rows) in zip(axes, panels):
        _plot_modes(ax, rows, key, theme, show_support=False)
        knots = set()
        for r in rows:
            knots.add(float(r["support_lower"]))
            knots.add(float(r["support_upper"]))
        O.draw_knots(ax, knots, theme)
        ax.axhline(0.0, color=theme.grid, lw=0.8)
        ax.set_title("%s\n%s" % (name, O.SEMANTICS_LABEL[sem]), fontsize=9,
                     color=theme.foreground)
        O.style_axes(ax, theme, xlabel, key)
    fig.suptitle(title, fontsize=12, color=theme.foreground, x=0.010,
                 ha="left", y=0.985)
    O.banner(fig, theme, note, y=0.930)
    O.footer(fig, theme, "Dotted verticals are knot positions. The three bases "
                         "use THREE DIFFERENT boundary knot conventions, shown "
                         "per panel. " + O.SCOPE)
    fig.tight_layout(rect=(0, 0.055, 1, 0.885))
    inputs = [os.path.join(rec, "axisymmetric_basis_profiles.csv"),
              os.path.join(rec, "univariate_basis_profiles.csv"),
              os.path.join(rec, "univariate_basis_derivatives.csv"),
              os.path.join(rec, "f0_coefficients.csv")]
    return O.save_figure(fig, out, stem, manifest=manifest, inputs=inputs,
                         theme=theme,
                         extra={"derivative_channel": key,
                                "production_derivative": "analytic recurrence",
                                "finite_difference_role": "verification only",
                                "continuity_label":
                                    "C4 at every ordinary interior knot"},
                         source_rows="every row of axisymmetric_basis_profiles."
                                     "csv with is_raw_support = 1, and every row "
                                     "of univariate_basis_profiles.csv / "
                                     "univariate_basis_derivatives.csv")


# --- figures ---------------------------------------------------------------

def fig_basis_functions(rec, out, manifest, theme):
    panels = _basis_panels(rec)
    fig, axes = plt.subplots(1, 3, figsize=(16.4, 5.2),
                             facecolor=theme.background)
    for ax, (name, sem, xlabel, rows) in zip(axes, panels):
        _plot_modes(ax, rows, "value", theme)
        knots = set()
        for r in rows:
            knots.add(float(r["support_lower"]))
            knots.add(float(r["support_upper"]))
        O.draw_knots(ax, knots, theme)
        n_active = len({int(r["index"]) for r in rows
                        if r["exclusion_reason"] == "active"})
        n_total = len({int(r["index"]) for r in rows})
        ax.set_title("%s\n%s\n%d of %d modes active" %
                     (name, O.SEMANTICS_LABEL[sem], n_active, n_total),
                     fontsize=9, color=theme.foreground)
        O.style_axes(ax, theme, xlabel, r"$N_i$")
    fig.suptitle("O0A  quintic B-spline basis functions, degree 5, C4 at every "
                 "ordinary interior knot", fontsize=12, color=theme.foreground,
                 x=0.010, ha="left", y=0.985)
    O.banner(fig, theme,
             "Endpoint multiplicities: clamped = 6 at both ends; uniform = 1 at "
             "both ends; lower-clamped/open-upper = 6 then 1. A simple interior "
             "knot leaves C4, verified by the scaling of the fourth-derivative "
             "jump, not by inspection.", y=0.930)
    O.footer(fig, theme,
             "A removed mode is one whose SUPPORT meets the outer flat collar or "
             "the protected cabin; it has no coefficient at all, which is what "
             "makes those regions exactly zero. " + O.SCOPE)
    fig.tight_layout(rect=(0, 0.055, 1, 0.875))
    inputs = [os.path.join(rec, "axisymmetric_basis_profiles.csv"),
              os.path.join(rec, "univariate_basis_profiles.csv"),
              os.path.join(rec, "f0_coefficients.csv"),
              os.path.join(rec, "axial_knot_vector.csv"),
              os.path.join(rec, "radial_squared_knot_vector.csv")]
    return O.save_figure(fig, out, "quintic_basis_functions", manifest=manifest,
                         inputs=inputs, theme=theme,
                         extra={"knot_semantics_shown":
                                    list(O.SEMANTICS_LABEL.values()),
                                "continuity_label":
                                    "C4 at every ordinary (multiplicity 1) "
                                    "interior knot; reduced to C3 at a "
                                    "multiplicity 2 knot",
                                "endpoint_multiplicities":
                                    {"clamped": 6, "uniform": 1,
                                     "lower_clamped_open_upper": [6, 1]}},
                         source_rows="all rows of axisymmetric_basis_profiles."
                                     "csv with is_raw_support = 1 plus all of "
                                     "univariate_basis_profiles.csv")


def fig_atlas(rec, out, manifest, theme):
    rows = O.read_csv(os.path.join(rec, "tensor_product_atlas.csv"))
    by = O.group_by(rows, "flat_index", int)
    keys = sorted(by)[:6]
    fig, axes = plt.subplots(2, 3, figsize=(15.6, 8.0),
                             facecolor=theme.background)
    flat = [a for r in axes for a in r]
    lim = max(abs(float(r["value"])) for r in rows) or 1.0
    for ax, k in zip(flat, keys):
        pts = by[k]
        xs = sorted({float(r["xi"]) for r in pts})
        ss = sorted({float(r["s"]) for r in pts})
        xi_i = {v: i for i, v in enumerate(xs)}
        s_i = {v: i for i, v in enumerate(ss)}
        grid = [[0.0] * len(xs) for _ in ss]
        for r in pts:
            grid[s_i[float(r["s"])]][xi_i[float(r["xi"])]] = float(r["value"])
        im = ax.imshow(grid, origin="lower",
                       extent=[min(xs), max(xs), min(ss), max(ss)],
                       aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim,
                       interpolation="nearest")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=6.5, colors=theme.foreground)
        cb.outline.set_visible(False)
        ax.set_title(r"$B_{ij}$  flat %d   $(i,j)=(%s,%s)$" %
                     (k, pts[0]["index_xi"], pts[0]["index_s"]),
                     fontsize=9, color=theme.foreground)
        O.style_axes(ax, theme, r"$\xi=(x-x_c)/R_0$  (dimensionless)",
                     r"$s=(y^2+z^2)/R_0^2$  (dimensionless)")
    for ax in flat[len(keys):]:
        ax.axis("off")
    fig.suptitle(r"O0A  tensor-product atlas  $B_{ij}(\xi,s)=N_i(\xi)M_j(s)$",
                 fontsize=12, color=theme.foreground, x=0.010, ha="left",
                 y=0.988)
    O.banner(fig, theme,
             "THE VERTICAL AXIS IS s, NOT A RADIUS IN METRES. "
             r"$s=(y^2+z^2)/R_0^2$ is dimensionless; a physical radius is "
             r"$\rho=R_0\sqrt{s}$, so equal steps in s are UNEQUAL steps in "
             "$\\rho$ and narrow as $\\rho$ grows. Cartesian interpretation: a "
             "band in s is an annulus about the axis of travel, and a band in "
             r"$\xi$ is a slab perpendicular to it.", y=0.958)
    O.footer(fig, theme, "Values read from tensor_product_atlas.csv, which the "
                         "C++ gate wrote; the renderer never forms the product "
                         "itself. " + O.SCOPE)
    fig.tight_layout(rect=(0, 0.035, 1, 0.900))
    return O.save_figure(fig, out, "tensor_product_basis_atlas",
                         manifest=manifest,
                         inputs=[os.path.join(rec, "tensor_product_atlas.csv")],
                         theme=theme,
                         extra={"radial_axis_quantity": "s, dimensionless",
                                "modes_shown": keys},
                         source_rows="all rows of tensor_product_atlas.csv")


def fig_support_map(rec, out, manifest, theme):
    rows = O.read_csv(os.path.join(rec, "tensor_product_supports.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 5.6),
                             facecolor=theme.background)
    ax = axes[0]
    for r in rows:
        xl, xh = float(r["xi_lower"]), float(r["xi_upper"])
        sl, sh = float(r["s_lower"]), float(r["s_upper"])
        ax.add_patch(plt.Rectangle((xl, sl), xh - xl, sh - sl, fill=False,
                                   edgecolor=ACTIVE_COLOR, linewidth=0.8,
                                   alpha=0.6))
    reg = manifest["design_region"]
    ax.add_patch(plt.Rectangle(
        (reg["xi_min"], 0.0), reg["xi_max"] - reg["xi_min"], reg["s_max"],
        fill=False, edgecolor=theme.foreground, lw=1.2, ls="--"))
    ax.axvline(reg["xi_min"] + reg["outer_flat_collar"], color="#d1495b", lw=1.2)
    ax.axvline(reg["xi_max"] - reg["outer_flat_collar"], color="#d1495b", lw=1.2)
    ax.axhline(reg["s_max"] - reg["outer_flat_collar"], color="#d1495b", lw=1.2)
    ax.axhline(reg["protected_cabin_radius"], color="#ffb454", lw=1.2, ls="-.")
    ax.set_xlim(reg["xi_min"] - 0.2, reg["xi_max"] + 0.2)
    ax.set_ylim(-0.15, reg["s_max"] + 0.15)
    ax.set_title("supports of the %d active tensor products" % len(rows),
                 fontsize=9.5, color=theme.foreground)
    O.style_axes(ax, theme, r"$\xi$ (dimensionless)", r"$s$ (dimensionless)")
    ax.legend(handles=[
        Line2D([], [], color=ACTIVE_COLOR, lw=0.9, label="active support box"),
        Line2D([], [], color=theme.foreground, lw=1.2, ls="--",
               label="design region"),
        Line2D([], [], color="#d1495b", lw=1.2, label="outer flat collar edge"),
        Line2D([], [], color="#ffb454", lw=1.2, ls="-.",
               label="protected cabin radius")],
        fontsize=6.6, frameon=False, labelcolor=theme.foreground,
        loc="upper right")

    ax = axes[1]
    ax.scatter(O.col(rows, "rho_lower"), O.col(rows, "rho_upper"), s=22,
               color=ACTIVE_COLOR, alpha=0.8)
    ax.plot([0, max(O.col(rows, "rho_upper"))],
            [0, max(O.col(rows, "rho_upper"))], color=theme.grid, lw=0.8,
            ls=":")
    ax.set_title("the same supports in PHYSICAL radius, metres", fontsize=9.5,
                 color=theme.foreground)
    O.style_axes(ax, theme, r"$\rho$ lower (m)", r"$\rho$ upper (m)")

    fig.suptitle("O0A  basis support map", fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    O.banner(fig, theme,
             "Every active support lies strictly inside the collar edges, which "
             "is why the collar is exactly zero: a mode whose support met it was "
             "removed rather than damped.", y=0.935)
    O.footer(fig, theme, O.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.885))
    hdr = ["flat_index", "index_xi", "index_s", "xi_lower", "xi_upper",
           "s_lower", "s_upper", "rho_lower", "rho_upper"]
    return O.save_figure(
        fig, out, "basis_support_map", manifest=manifest,
        inputs=[os.path.join(rec, "tensor_product_supports.csv")], theme=theme,
        plotted_header=hdr, plotted_rows=[[r[k] for k in hdr] for r in rows])


def fig_axis_modes(rec, out, manifest, theme):
    rows = O.read_csv(os.path.join(rec, "axis_regular_modes.csv"))
    fams = ["scalar", "radial_vector", "axial_radial_tensor",
            "transverse_anisotropy"]
    formulas = {
        "scalar": r"$q=Q(\xi,s)$",
        "radial_vector": r"$q_y=\frac{y}{R_0}Q$,   $q_z=\frac{z}{R_0}Q$",
        "axial_radial_tensor":
            r"$q_{xy}=\frac{y}{R_0}Q$,   $q_{xz}=\frac{z}{R_0}Q$",
        "transverse_anisotropy":
            r"$q_{yy}-q_{zz}=\frac{y^2-z^2}{R_0^2}Q$,   "
            r"$q_{yz}=\frac{yz}{R_0^2}Q$"}
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 8.4),
                             facecolor=theme.background)
    flat = [a for r in axes for a in r]
    for ax, fam in zip(flat, fams):
        sub = [r for r in rows if r["family"] == fam]
        ax.plot(O.col(sub, "component_expected"), O.col(sub, "component_actual"),
                marker="o", markersize=4, linestyle="none",
                color=ACTIVE_COLOR, label="sample")
        lim = max([abs(v) for v in O.col(sub, "component_expected")] + [1e-30])
        ax.plot([-lim, lim], [-lim, lim], color=theme.grid, lw=0.9, ls=":",
                label="exact agreement")
        worst = max(abs(v) for v in O.col(sub, "difference"))
        drho = max(abs(v) for v in O.col(sub, "d_rho_on_axis"))
        order = sub[0]["axis_vanishing_order"]
        ax.set_title("%s   vanishing order %s\n%s" %
                     (fam.replace("_", " "), order, formulas[fam]),
                     fontsize=9, color=theme.foreground)
        ax.set_xlabel("closed-form component   (worst |diff| = %.1e, "
                      r"max $|\partial B/\partial\rho|$ on axis = %.1e)"
                      % (worst, drho), fontsize=7.4, color=theme.foreground)
        O.style_axes(ax, theme, ax.get_xlabel(), "component from the basis")
        ax.legend(fontsize=6.6, frameon=False, labelcolor=theme.foreground)
    fig.suptitle("O0A  structurally regular field families: Cartesian "
                 "components against their closed forms", fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.988)
    O.banner(fig, theme,
             r"NO DIVISION BY $\rho$ OCCURS ANYWHERE. Projecting the "
             r"$\rho$-component onto $y$ gives "
             r"$q_\rho\,(y/\rho)=(\rho/R_0)Q\,(y/\rho)=(y/R_0)Q$: the $\rho$ "
             "cancels ALGEBRAICALLY, so each Cartesian component is a polynomial "
             r"in $y,z$ times a spline. There is no axis limit to take, no "
             r"small-$\rho$ branch, and no regularity penalty.", y=0.955)
    O.footer(fig, theme,
             "Samples include the axis itself and both signs of y and z. " +
             O.SCOPE)
    fig.tight_layout(rect=(0, 0.035, 1, 0.895))
    hdr = ["family", "axis_vanishing_order", "x", "y", "z", "rho", "Q",
           "component_expected", "component_actual", "difference",
           "d_rho_on_axis"]
    return O.save_figure(
        fig, out, "axis_regular_modes", manifest=manifest,
        inputs=[os.path.join(rec, "axis_regular_modes.csv")], theme=theme,
        extra={"divides_by_rho": False},
        plotted_header=hdr, plotted_rows=[[r[k] for k in hdr] for r in rows])


def fig_constraint_projection(rec, out, manifest, theme):
    proj = O.read_csv(os.path.join(rec, "constraint_projection.csv"))
    sv = O.read_csv(os.path.join(rec, "constraint_singular_values.csv"))
    summary = O.read_json(os.path.join(rec, "verification_summary.json"))
    con = manifest["constraints"]

    def matrix(name):
        sub = [r for r in proj if r["matrix"] == name]
        if not sub:
            return None, 0, 0
        nr = max(int(r["row"]) for r in sub) + 1
        nc = max(int(r["col"]) for r in sub) + 1
        m = [[0.0] * nc for _ in range(nr)]
        for r in sub:
            m[int(r["row"])][int(r["col"])] = float(r["value"])
        return m, nr, nc

    P, pr_, pc = matrix("P")
    PW, wr_, wc = matrix("PW")

    fig, axes = plt.subplots(1, 3, figsize=(15.8, 5.4),
                             facecolor=theme.background)
    lim = max(max(abs(v) for v in row) for row in P) or 1.0
    for ax, (m, name, cols, xlabel) in zip(
            axes[:2],
            [(P, r"$P$: raw $c$ $\leftarrow$ constrained $z$", pc,
              r"constrained coordinate $z_j$  (%d free)" % pc),
             (PW, r"$PW$: raw $c$ $\leftarrow$ whitened $y$", wc,
              r"whitened coordinate $y_j$  (%d retained)" % wc)]):
        im = ax.imshow(m, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim,
                       interpolation="nearest")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=6.5, colors=theme.foreground)
        cb.outline.set_visible(False)
        ax.set_title(name, fontsize=9.5, color=theme.foreground)
        O.style_axes(ax, theme, xlabel, r"raw coefficient $c_i$  (%d)" % pr_)

    ax = axes[2]
    idx = [int(r["index"]) for r in sv]
    vals = [max(float(r["singular_value"]), 1e-20) for r in sv]
    keep = [r["retained"] == "1" for r in sv]
    ax.semilogy([i for i, k in zip(idx, keep) if k],
                [v for v, k in zip(vals, keep) if k], marker="o", markersize=6,
                linestyle="none", color="#d1495b", label="retained (rank)")
    ax.semilogy([i for i, k in zip(idx, keep) if not k],
                [v for v, k in zip(vals, keep) if not k], marker="x",
                markersize=5, linestyle="none", color=REMOVED_COLOR,
                label="below the rank tolerance")
    ax.axhline(float(sv[0]["tolerance"]), color=theme.accent_two, lw=1.1,
               ls="--", label="rank tolerance")
    ax.legend(fontsize=6.8, frameon=False, labelcolor=theme.foreground)
    ax.set_title("constraint singular values: the rank gap", fontsize=9.5,
                 color=theme.foreground)
    O.style_axes(ax, theme, "index", "singular value")

    aff = null = orth = "n/a"
    for c in summary["cases"]:
        if c["name"] == "o0a.affine_constraint_residual":
            aff = "%.2e" % float(c["computed"])
        if c["name"] == "o0a.null_space_residual":
            null = "%.2e" % float(c["computed"])
    fig.suptitle("O0A  constraint elimination:  "
                 r"$c=c_0+Pz$,   $z=Wy$", fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    O.banner(fig, theme,
             "raw coefficients c: %d   |   constraint rows: %d, RANK %d, %d "
             "redundant   |   constrained coordinates z: %d free   |   whitened "
             "coordinates y: %d retained\n"
             "affine residual |C c0 - d| = %s      null-space residual |C P| = "
             "%s      the three coordinate systems are DISTINCT and are never "
             "interchanged"
             % (manifest["f1_tensor"]["active"], con["rows"], con["rank"],
                con["redundant"], con["free"], wc, aff, null), y=0.945)
    O.footer(fig, theme,
             "Constraints are ELIMINATED, not penalised: every iterate c0 + P z "
             "satisfies every constraint to machine precision by construction. " +
             O.WHITENING_NOTE + " " + O.SCOPE)
    fig.tight_layout(rect=(0, 0.055, 1, 0.870))
    return O.save_figure(
        fig, out, "constraint_projection", manifest=manifest,
        inputs=[os.path.join(rec, "constraint_projection.csv"),
                os.path.join(rec, "constraint_singular_values.csv"),
                os.path.join(rec, "verification_summary.json")],
        theme=theme,
        extra={"raw_dimension": manifest["f1_tensor"]["active"],
               "constraint_rank": con["rank"], "redundant_rows": con["redundant"],
               "free_coordinates": con["free"], "whitened_coordinates": wc,
               "whitening_note": O.WHITENING_NOTE},
        source_rows="all rows of constraint_projection.csv and "
                    "constraint_singular_values.csv")


def fig_gram_spectrum(rec, out, manifest, theme):
    raw = O.read_csv(os.path.join(rec, "gram_spectrum_raw.csv"))
    con = O.read_csv(os.path.join(rec, "gram_spectrum_constrained.csv"))
    whi = O.read_csv(os.path.join(rec, "gram_spectrum_whitened.csv"))
    fix = O.read_csv(os.path.join(rec, "linalg_symmetric_fixtures.csv"))

    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.4),
                             facecolor=theme.background)
    ax = axes[0]
    for rows, name, colr, mk in ((raw, "raw  M", "#4ea1ff", "o"),
                                 (con, r"constrained  $G=P^TMP$", "#ffb454", "s"),
                                 (whi, r"whitened  $W^TGW$", "#7ed957", "^")):
        vals = [max(float(r["eigenvalue"]), 1e-20) for r in rows]
        ax.semilogy(range(len(vals)), vals, marker=mk, markersize=4,
                    color=colr, lw=1.3,
                    label="%s   cond %.6g" % (name,
                                              float(rows[0]["condition_number"])))
    ax.axhline(float(raw[0]["condition_number"]) * 0.0 + 1.0,
               color=theme.grid, lw=0.8, ls=":")
    ax.legend(fontsize=7.0, frameon=False, labelcolor=theme.foreground)
    ax.set_title("Sobolev Gram spectra: the whitened spectrum is flat at 1",
                 fontsize=9.5, color=theme.foreground)
    O.style_axes(ax, theme, "mode index", "eigenvalue")

    ax = axes[1]
    names = [r["fixture"] for r in fix]
    disc = []
    for r in fix:
        lo = float(r["eigenvalue_min"])
        hi = float(r["eigenvalue_max"])
        disc.append(1 if lo <= 1e-12 * hi else 0)
    xs = range(len(names))
    ax.bar(xs, [float(r["largest_group"]) for r in fix], width=0.55,
           color=REMOVED_COLOR, alpha=0.75, label="largest degenerate group")
    ax.bar(xs, disc, width=0.28, color="#d1495b",
           label="modes discarded as numerically null")
    ax.set_xticks(list(xs))
    ax.set_xticklabels([n.replace("_", "\n") for n in names], fontsize=6.4,
                       color=theme.foreground)
    ax.legend(fontsize=6.8, frameon=False, labelcolor=theme.foreground)
    ax.set_title("INDEPENDENT FIXTURES, not production modes", fontsize=9.5,
                 color=theme.foreground)
    O.style_axes(ax, theme, "", "count")

    n_disc = sum(1 for r in whi if r["retained"] != "1")
    fig.suptitle("O0A  Gram conditioning and whitening", fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    O.banner(fig, theme,
             "PRODUCTION CASE: %d modes retained, %d discarded. The rank cutoff "
             "sits far below the smallest production eigenvalue, so nothing is "
             "dropped here. The right panel is the exact-null and nearly-singular "
             "FIXTURES, which each lose exactly one mode by design -- they are "
             "test matrices, not basis modes."
             % (len(whi) - n_disc, n_disc), y=0.940)
    O.footer(fig, theme,
             "A numerically null mode is DISCARDED and counted, never divided "
             "by: dividing by sqrt(lambda) there manufactures an enormous "
             "optimization direction out of rounding error. " + O.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.880))
    return O.save_figure(
        fig, out, "gram_spectrum", manifest=manifest,
        inputs=[os.path.join(rec, "gram_spectrum_raw.csv"),
                os.path.join(rec, "gram_spectrum_constrained.csv"),
                os.path.join(rec, "gram_spectrum_whitened.csv"),
                os.path.join(rec, "linalg_symmetric_fixtures.csv")],
        theme=theme,
        extra={"production_discarded_modes": n_disc,
               "fixture_panel_is_not_production": True},
        source_rows="all rows of the three gram_spectrum_*.csv records and "
                    "linalg_symmetric_fixtures.csv")


def fig_f0_projection(rec, out, manifest, theme):
    rows = O.read_csv(os.path.join(rec, "f0_baseline_projection.csv"))
    errs = {r["zone"]: r for r in
            O.read_csv(os.path.join(rec, "f0_projection_errors.csv"))}
    r_ = O.col(rows, "radius")
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.0),
                             facecolor=theme.background)
    trip = (("f_reference", "f_spline", r"$f$", "value"),
            ("f_first_reference", "f_first_spline", r"$f'$", "first"),
            ("f_second_reference", "f_second_spline", r"$f''$", "second"))
    for ax, (a, b, lab, _k) in zip(axes, trip):
        ax.plot(r_, O.col(rows, a), lw=2.0, color=theme.foreground,
                label="Alcubierre reference")
        ax.plot(r_, O.col(rows, b), lw=1.5, color="#d1495b", ls="--",
                label="constrained spline projection")
        ax.axvspan(float(manifest["scale_R0"]) *
                   math.sqrt(manifest["design_region"]["s_max"] -
                             manifest["design_region"]["outer_flat_collar"]),
                   max(r_), color=REMOVED_COLOR, alpha=0.18)
        ax.legend(fontsize=6.8, frameon=False, labelcolor=theme.foreground)
        ax.set_title(lab, fontsize=10, color=theme.foreground)
        O.style_axes(ax, theme, "radius r (m)", lab)
    fig.suptitle("O0A  F0 baseline projection of the Alcubierre shape function "
                 r"($\sigma$ = %s, %s free of %s coefficients)"
                 % (manifest["f0"]["wall_sigma"], manifest["f0"]["free"],
                    manifest["f0"]["raw"]),
                 fontsize=12, color=theme.foreground, x=0.010, ha="left",
                 y=0.985)
    O.f0_stamp(fig, theme)
    O.banner(fig, theme,
             "Shaded band: the outer flat collar, where 7 of 12 coefficients are "
             "FIXED AT ZERO, so the projection is required to be exactly zero and "
             "cannot represent f there. That truncation is deliberate and is "
             "reported as the collar error (%s), not hidden as a fit residual."
             % errs["flat_collar"]["max_value_error"], y=0.845)
    O.footer(fig, theme, O.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.800))
    hdr = ["s_r", "radius", "f_reference", "f_spline", "f_first_reference",
           "f_first_spline", "f_second_reference", "f_second_spline", "zone"]
    return O.save_figure(
        fig, out, "f0_baseline_projection", manifest=manifest,
        inputs=[os.path.join(rec, "f0_baseline_projection.csv"),
                os.path.join(rec, "f0_projection_errors.csv")], theme=theme,
        extra={"f0_resolution_certified": False,
               "f0_ready_for_optimization": False,
               "status_banner": O.F0_BANNER},
        plotted_header=hdr, plotted_rows=[[r[k] for k in hdr] for r in rows])


def fig_f0_error(rec, out, manifest, theme):
    rows = O.read_csv(os.path.join(rec, "f0_baseline_projection.csv"))
    errs = O.read_csv(os.path.join(rec, "f0_projection_errors.csv"))
    r_ = O.col(rows, "radius")
    # Three SEPARATE panels on their own linear scales. Combining them into one
    # normalised curve, or choosing a scale that flattens 1.4265, is exactly the
    # misrepresentation the specification forbids.
    trip = (("f_reference", "f_spline", r"$|\Delta f|$", "max_value_error"),
            ("f_first_reference", "f_first_spline", r"$|\Delta f'|$",
             "max_first_derivative_error"),
            ("f_second_reference", "f_second_spline", r"$|\Delta f''|$",
             "max_second_derivative_error"))
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.0),
                             facecolor=theme.background)
    for ax, (a, b, lab, key) in zip(axes, trip):
        diff = [abs(x - y) for x, y in zip(O.col(rows, a), O.col(rows, b))]
        ax.plot(r_, diff, lw=1.6, color="#d1495b")
        top = max(diff)
        ax.set_ylim(0.0, top * 1.15)
        wall = [e for e in errs if e["zone"] == "wall"][0]
        ax.axhline(float(wall[key]), color=theme.accent_two, lw=1.1, ls="--")
        ax.annotate("wall maximum %.4g" % float(wall[key]),
                    (r_[len(r_) // 2], float(wall[key])),
                    textcoords="offset points", xytext=(4, 5), fontsize=7.6,
                    color=theme.accent_two)
        ax.set_title("%s   (own linear scale, axis top %.4g)" % (lab, top * 1.15),
                     fontsize=9.2, color=theme.foreground)
        O.style_axes(ax, theme, "radius r (m)", lab)
    fig.suptitle("O0A  F0 projection error, value / first / second derivative, "
                 "SEPARATELY", fontsize=12, color=theme.foreground, x=0.010,
                 ha="left", y=0.985)
    O.f0_stamp(fig, theme)
    O.banner(fig, theme,
             "The second-derivative error reaches %s in the wall. It is plotted "
             "on its own linear axis and NOT normalised against the other two: "
             "curvature is built from second derivatives, so an error of order "
             "one here would be an error of order one in the geometry. This is "
             "the measurement that blocks optimization readiness."
             % [e for e in errs if e["zone"] == "wall"][0][
                 "max_second_derivative_error"], y=0.845)
    O.footer(fig, theme, "Five free coefficients cannot resolve a wall narrower "
                         "than one knot span. O0A records that rather than "
                         "refining the basis. " + O.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.800))
    hdr = ["zone", "samples", "max_value_error", "max_first_derivative_error",
           "max_second_derivative_error", "reference_scale"]
    return O.save_figure(
        fig, out, "f0_projection_error", manifest=manifest,
        inputs=[os.path.join(rec, "f0_baseline_projection.csv"),
                os.path.join(rec, "f0_projection_errors.csv")], theme=theme,
        extra={"separate_axes_per_derivative_order": True,
               "f0_resolution_certified": False,
               "f0_ready_for_optimization": False,
               "status_banner": O.F0_BANNER},
        plotted_header=hdr, plotted_rows=[[r[k] for k in hdr] for r in errs])


def fig_resolution(rec, out, manifest, theme):
    rows = O.read_csv(os.path.join(rec, "resolution_coverage.csv"))
    f1 = [r for r in rows if r["family"] == "f1_scalar"]
    f0 = [r for r in rows if r["family"] == "f0_radial"]
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.4),
                             facecolor=theme.background)
    ax = axes[0]
    for sub, name, colr, mk in ((f1, "F1 scalar", "#4ea1ff", "o"),
                                (f0, "F0 radial", "#ffb454", "s")):
        ax.plot(range(len(sub)), O.col(sub, "cells_across_narrowest_feature"),
                marker=mk, markersize=5, linestyle="none", color=colr,
                label=name)
    ax.axhline(6.0, color="#d1495b", lw=1.6,
               label="six-cell acceptance threshold")
    ax.axhspan(0, 6.0, color="#d1495b", alpha=0.10)
    ax.legend(fontsize=7.0, frameon=False, labelcolor=theme.foreground)
    ax.set_title("cells across the narrowest derivative feature, per active mode",
                 fontsize=9.2, color=theme.foreground)
    O.style_axes(ax, theme, "active mode (ordinal)", "grid cells")

    ax = axes[1]
    f1min = min(O.col(f1, "cells_across_narrowest_feature"))
    f0min = min(O.col(f0, "cells_across_narrowest_feature"))
    bars = ax.bar([0, 1], [f1min, f0min], width=0.5,
                  color=["#4ea1ff", "#ffb454"])
    ax.axhline(6.0, color="#d1495b", lw=1.8)
    for b, v, lab in zip(bars, [f1min, f0min], ["F1", "F0"]):
        ax.annotate("%s\n%.3f cells\nFAILS" % (lab, v),
                    (b.get_x() + b.get_width() / 2, v),
                    textcoords="offset points", xytext=(0, 6), ha="center",
                    fontsize=8.4, color=theme.uncertainty_color)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["F1 axisymmetric kernel", "F0 radial validation"],
                       fontsize=7.6, color=theme.foreground)
    ax.set_ylim(0, 7.2)
    ax.set_title("minimum over modes, against the contract", fontsize=9.2,
                 color=theme.foreground)
    O.style_axes(ax, theme, "", "grid cells")

    fig.suptitle("O0A  resolution coverage on the accepted grid "
                 "(21^3 over [-4,4]^3, spacing 0.4 m)", fontsize=12,
                 color=theme.foreground, x=0.010, ha="left", y=0.985)
    O.banner(fig, theme,
             "BOTH FAMILIES FAIL THE RESOLUTION CONTRACT. F1 %.3f cells, F0 "
             "%.3f cells, against a six-cell threshold. Neither is "
             "resolution-certified and neither is ready for optimization. Equal "
             "steps in s are unequal steps in rho, so the tightest radial "
             "feature sits at the OUTER edge of the active region."
             % (f1min, f0min), y=0.940)
    O.footer(fig, theme, "Evaluated per active mode, not on an average support "
                         "width. " + O.SCOPE)
    fig.tight_layout(rect=(0, 0.045, 1, 0.880))
    hdr = ["family", "flat_index", "narrowest_feature_width", "grid_spacing",
           "cells_across_support", "cells_across_narrowest_feature",
           "meets_six_cell_contract"]
    return O.save_figure(
        fig, out, "resolution_coverage", manifest=manifest,
        inputs=[os.path.join(rec, "resolution_coverage.csv")], theme=theme,
        extra={"six_cell_threshold": 6.0, "f1_min_cells": f1min,
               "f0_min_cells": f0min, "f0_resolution_certified": False,
               "f1_resolution_certified": False,
               "f0_ready_for_optimization": False,
               "f1_ready_for_optimization": False},
        plotted_header=hdr, plotted_rows=[[r[k] for k in hdr] for r in rows])


# --- driver ----------------------------------------------------------------

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

    manifest = O.read_json(os.path.join(rec, "basis_manifest.json"))
    before = O.source_checksums(rec)
    with open(os.path.join(out, "source_checksums.json"), "w",
              encoding="utf-8") as fh:
        json.dump(before, fh, indent=2, sort_keys=True)

    made = []
    made.append(fig_basis_functions(rec, out, manifest, theme))
    made.append(_derivative_figure(
        rec, out, manifest, theme, "first", "basis_first_derivatives",
        "O0A  first derivatives of the quintic basis",
        "Production derivative: ANALYTIC RECURRENCE. Finite differences appear "
        "in O0A only as an independent verification path and are never plotted "
        "here as production data."))
    made.append(_derivative_figure(
        rec, out, manifest, theme, "second", "basis_second_derivatives",
        "O0A  second derivatives of the quintic basis",
        "Production derivative: ANALYTIC RECURRENCE. These are the derivatives "
        "curvature will be built from, which is why degree 5 is required: simple "
        "interior knots leave C4, so the second derivative is C2."))
    made.append(_derivative_figure(
        rec, out, manifest, theme, "third", "basis_third_derivatives",
        "O0A  third derivatives of the quintic basis",
        "Production derivative: ANALYTIC RECURRENCE.   Finite-difference role: "
        "VERIFICATION ONLY.   No third finite-difference derivative is plotted "
        "here. The gate's independent check needed a halving step ladder, "
        "because a centred third difference divides by 2h^3 and its rounding "
        "floor grows like eps/h^3."))
    made.append(fig_atlas(rec, out, manifest, theme))
    made.append(fig_support_map(rec, out, manifest, theme))
    made.append(fig_axis_modes(rec, out, manifest, theme))
    made.append(fig_constraint_projection(rec, out, manifest, theme))
    made.append(fig_gram_spectrum(rec, out, manifest, theme))
    made.append(fig_f0_projection(rec, out, manifest, theme))
    made.append(fig_f0_error(rec, out, manifest, theme))
    made.append(fig_resolution(rec, out, manifest, theme))

    after = O.source_checksums(rec)
    rep = anomaly.AnomalyReport("O", "O0A", "o0a_static")
    changed = [k for k in before if before[k] != after.get(k)]
    rep.record("source_checksum_mismatch", fired=bool(changed),
               measured=len(changed), threshold=0,
               detail="records whose digest changed while the figures were "
                      "drawn: %r" % changed[:3])
    blanks = [os.path.basename(p) for p in made
              if _png_is_blank(p)]
    rep.record("blank_figure", fired=bool(blanks), measured=len(blanks),
               threshold=0, detail="figures with too few distinct colours: %r"
                                   % blanks[:3])
    rep.record("full_hessian_overclaim",
               fired=bool(manifest.get("full_hessian_implemented")),
               measured=manifest.get("full_hessian_implemented"),
               threshold=False,
               detail="the manifest must not claim a full Hessian; " +
                      O.HESSIAN_NOTE)
    path = rep.write(out)

    for p in made:
        print("wrote %s" % os.path.basename(p))
    print("%d figures, %d fired anomalies -> %s"
          % (len(made), len(rep.fired), os.path.basename(path)))
    return 0


def _png_is_blank(path, min_colors=64):
    try:
        from PIL import Image
    except ImportError:
        return False
    im = Image.open(path).convert("RGB")
    cols = im.getcolors(maxcolors=1 << 20)
    return cols is not None and len(cols) < min_colors


if __name__ == "__main__":
    raise SystemExit(main())
