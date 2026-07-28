#!/usr/bin/env python3
"""Milestone A5 figures: the Einstein tensor and the inferred stress-energy.

Redrawn entirely from the CSVs the A5 gate wrote. Nothing here computes a
curvature, an Einstein tensor or a stress-energy tensor.

SCOPE, carried on every figure: A5 computes the geometric Einstein tensor and
the stress-energy inferred through Einstein's equation for the selected metric.
It does not establish that this source can be generated physically, that energy
conditions are satisfied, or that FTL travel is realizable. No quantity plotted
here is an energy density -- that requires an observer, which A6 introduces.

Usage:
    python python/plot_milestone5.py results/verification/milestone_05
"""

from __future__ import annotations

import csv
import math
import os
import sys
import textwrap

import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wh_style as S  # noqa: E402

NOT_ENERGY = ("These are coordinate-basis components of the geometric Einstein "
              "tensor. They are NOT observer-measured energy densities.")


def wrap(text, width=150):
    return "\n".join(textwrap.wrap(" ".join(text.split()), width))


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def col(rows, name):
    return [float(r[name]) for r in rows]


def save_both(fig, out, stem, rect=(0, 0.05, 1, 0.90)):
    """PNG for review, SVG for print. Same figure, written twice."""
    fig.tight_layout(rect=rect)
    for ext in ("png", "svg"):
        path = os.path.join(out, "%s.%s" % (stem, ext))
        fig.savefig(path, dpi=200, facecolor=S.SURFACE)
        print("wrote " + path)
    plt.close(fig)


def order_of(hs, es):
    if len(hs) < 2 or es[0] <= 0 or es[1] <= 0:
        return float("nan")
    return math.log(es[0] / es[1]) / math.log(hs[0] / hs[1])


# ---------------------------------------------------------------------------
# 1. The geometries whose Einstein tensor must vanish
# ---------------------------------------------------------------------------

def fig_reference(results, out):
    mink = read_csv(os.path.join(results, "einstein_minkowski.csv"))
    flat = read_csv(os.path.join(results, "einstein_flat_spherical.csv"))
    schw = read_csv(os.path.join(results, "einstein_schwarzschild.csv"))

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0), facecolor=S.SURFACE)

    ax = axes[0]
    ax.plot(col(mink, "x"), [abs(v) for v in col(mink, "G00")],
            color=S.SERIES[0], linewidth=1.4, label="$|G_{00}|$")
    ax.plot(col(mink, "x"), [abs(v) for v in col(mink, "trace")],
            color=S.SERIES[1], linewidth=1.4, linestyle="--", label="$|G^\\mu{}_\\mu|$")
    ax.set_xlabel("x"); ax.set_ylabel("magnitude")
    ax.set_title("Minkowski: exactly zero", fontsize=9, color=S.INK2)
    S.style_axes(ax)
    ax.legend(fontsize=7.5, frameon=False)

    ax = axes[1]
    ax.plot(col(flat, "r"), col(flat, "max_abs_G"), color=S.SERIES[2],
            linewidth=1.4)
    ax.set_yscale("log")
    ax.set_xlabel("r"); ax.set_ylabel(r"max$_{\mu\nu}$ $|G_{\mu\nu}|$")
    ax.set_title("flat spherical: no false source", fontsize=9, color=S.INK2)
    S.style_axes(ax, logy=True)

    ax = axes[2]
    ax.plot(col(schw, "r"), col(schw, "max_abs_G"), color=S.SERIES[2],
            linewidth=1.4, label=r"max $|G_{\mu\nu}|$")
    ax.plot(col(schw, "r"), col(schw, "kretschmann"), color=S.SERIES[1],
            linewidth=1.4, label="$R_{abcd}R^{abcd}$")
    ax.set_yscale("log")
    ax.set_xlabel("r"); ax.set_ylabel("magnitude")
    ax.set_title("Schwarzschild: vacuum, not flat", fontsize=9, color=S.INK2)
    S.style_axes(ax, logy=True)
    ax.legend(fontsize=7.5, frameon=False)

    S.title_block(
        fig, "Three geometries whose Einstein tensor must vanish",
        "The middle panel is the one that matters: a large connection must not "
        "become a false source.")
    S.footer(fig, wrap(
        "Minkowski passes trivially for any pipeline that returns zero, which "
        "is why the other two exist. Flat space in spherical coordinates has "
        "Gamma of order 1 everywhere and must still give G = 0. Schwarzschild "
        "is a vacuum whose Kretschmann invariant is emphatically nonzero, so a "
        "vanishing G there is not a vanishing curvature. " + S.FOOTER_UNITS, 210))
    save_both(fig, out, "einstein_reference_geometries")


# ---------------------------------------------------------------------------
# 2. Alcubierre components
# ---------------------------------------------------------------------------

def fig_alcubierre_components(results, out):
    rows = read_csv(os.path.join(results, "alcubierre_einstein_centerline.csv"))
    x = col(rows, "x")

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2), facecolor=S.SURFACE)

    ax = axes[0]
    ax.plot(x, col(rows, "g00"), color=S.SERIES[0], linewidth=1.4, label="$g_{00}$")
    ax.plot(x, col(rows, "g01"), color=S.SERIES[1], linewidth=1.4, label="$g_{01}$")
    ax.set_xlabel("x"); ax.set_ylabel("metric component")
    ax.set_title("the geometry", fontsize=9, color=S.INK2)
    S.style_axes(ax)
    ax.legend(fontsize=7.5, frameon=False)

    ax = axes[1]
    for key, lab, c in (("G00_symmetric", "$G_{00}$", 0),
                        ("G01_symmetric", "$G_{01}$", 1),
                        ("G11_symmetric", "$G_{11}$", 2),
                        ("G22_symmetric", "$G_{22}$", 3)):
        ax.plot(x, col(rows, key), color=S.SERIES[c], linewidth=1.4, label=lab)
    ax.axhline(0.0, color=S.INK3, linewidth=0.8)
    ax.set_xlabel("x"); ax.set_ylabel("coordinate component")
    ax.set_title("$G_{\\mu\\nu}$ — NOT an energy density", fontsize=9,
                 color=S.INK2)
    S.style_axes(ax)
    ax.legend(fontsize=7.5, frameon=False, ncol=2)

    ax = axes[2]
    ax.plot(x, col(rows, "einstein_square"), color=S.SERIES[2], linewidth=1.4,
            label=r"$G_{\mu\nu}G^{\mu\nu}$")
    ax.plot(x, col(rows, "einstein_trace"), color=S.SERIES[1], linewidth=1.4,
            label=r"$G^\mu{}_\mu$")
    ax.axhline(0.0, color=S.INK3, linewidth=0.8)
    ax.set_xlabel("x"); ax.set_ylabel("scalar")
    ax.set_title("scalar contractions", fontsize=9, color=S.INK2)
    S.style_axes(ax)
    ax.legend(fontsize=7.5, frameon=False)

    for ax in axes:
        for xv in (-2.0, 2.0):
            ax.axvline(xv, color=S.INK3, linewidth=0.8, linestyle="--")

    S.title_block(fig, "The Einstein tensor of an Alcubierre bubble, along the axis",
                  NOT_ENERGY)
    S.footer(fig, wrap(
        "y = z = 0, t = 0, v_s = 0.5, R = 2, sigma = 1; dashed lines mark the "
        "nominal wall. G_00 is POSITIVE throughout, which says nothing about "
        "the sign of any energy density: an energy density is T_{mu nu} u^mu "
        "u^nu for a specified normalised observer, and A5 has no observers. "
        "The familiar statement that an Alcubierre bubble needs negative energy "
        "is about a contraction A6 performs, not about this component. "
        + S.FOOTER_UNITS, 210))
    save_both(fig, out, "alcubierre_einstein_components")


# ---------------------------------------------------------------------------
# 3. Trace identity
# ---------------------------------------------------------------------------

def fig_trace_identity(results, out):
    rows = read_csv(os.path.join(results, "alcubierre_einstein_centerline.csv"))
    ds = read_csv(os.path.join(results, "einstein_desitter.csv"))
    x = col(rows, "x")
    tr = col(rows, "einstein_trace")
    minus_R = [-v for v in col(rows, "ricci_scalar")]
    resid = [abs(a - b) for a, b in zip(tr, minus_R)]

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0), facecolor=S.SURFACE)

    ax = axes[0]
    ax.plot(x, minus_R, color=S.SERIES[1], linewidth=3.0, alpha=0.45, label="$-R$")
    ax.plot(x, tr, color=S.SERIES[0], linewidth=1.4, label=r"$G^\mu{}_\mu$")
    ax.set_xlabel("x"); ax.set_ylabel("trace")
    ax.set_title("Alcubierre", fontsize=9, color=S.INK2)
    S.style_axes(ax)
    ax.legend(fontsize=7.5, frameon=False)

    ax = axes[1]
    ax.semilogy(x, [max(v, 1e-20) for v in resid], color=S.SERIES[2], linewidth=1.4)
    ax.axhline(1e-13, color=S.STATUS_CRIT, linewidth=1.0, linestyle="--")
    ax.text(x[1], 1.5e-13, "tolerance 1e-13", fontsize=7, color=S.STATUS_CRIT)
    ax.set_xlabel("x"); ax.set_ylabel(r"$|G^\mu{}_\mu + R|$")
    ax.set_title("residual: roundoff, not truncation", fontsize=9, color=S.INK2)
    S.style_axes(ax, logy=True)

    ax = axes[2]
    ax.plot(col(ds, "r"), col(ds, "trace"), color=S.SERIES[0],
            linewidth=1.4, label=r"$G^\mu{}_\mu$")
    ax.plot(col(ds, "r"), [-v for v in col(ds, "ricci_scalar")],
            color=S.SERIES[1], linewidth=3.0, alpha=0.45, label="$-R = -12/L^2$")
    ax.set_xlabel("r"); ax.set_ylabel("trace")
    ax.set_title("static de Sitter", fontsize=9, color=S.INK2)
    S.style_axes(ax)
    ax.legend(fontsize=7.5, frameon=False)

    S.title_block(fig, r"$G^\mu{}_\mu = -R$ holds to roundoff",
                  "An identity that survives discretisation exactly, unlike most "
                  "of its neighbours.")
    S.footer(fig, wrap(
        "g^{mu nu}(R_{mu nu} - g_{mu nu} R/2) = R - 4R/2 = -R, using g^{mu nu} "
        "g_{mu nu} = 4. The derivation never differentiates a product, so it "
        "holds for whatever R_{mu nu} the discretisation produced -- which is "
        "exactly the criterion A4 established for telling exact identities from "
        "O(h^2) ones. Measured worst residual 2.8e-17. " + S.FOOTER_UNITS, 210))
    save_both(fig, out, "einstein_trace_identity")


# ---------------------------------------------------------------------------
# 4 and 5. Convergence
# ---------------------------------------------------------------------------

def _convergence_axes(ax, hs, es, label, colour, ref=True):
    ax.plot(hs, es, "o-", color=colour, linewidth=1.4, markersize=4, label=label)
    if ref and hs and es:
        ax.plot([min(hs), max(hs)],
                [es[0] * (h / hs[0]) ** 2 for h in (min(hs), max(hs))],
                "--", color=S.INK3, linewidth=1.0, label="slope 2")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("grid spacing h")
    S.style_axes(ax, logx=True, logy=True)
    ax.legend(fontsize=7.5, frameon=False)


def fig_antisymmetry_convergence(results, out):
    rows = read_csv(os.path.join(results, "einstein_convergence.csv"))
    pts = sorted((float(r["h"]), float(r["antisymmetry_at_probe"])) for r in rows)
    hs = [p[0] for p in pts]
    es = [p[1] for p in pts]

    fig, ax = plt.subplots(figsize=(7.0, 4.8), facecolor=S.SURFACE)
    _convergence_axes(ax, hs, es,
                      "max $|G_{\\mu\\nu} - G_{\\nu\\mu}|$  (order %.2f)"
                      % order_of(hs, es), S.SERIES[0])
    ax.set_ylabel("antisymmetry at a fixed physical point")

    S.title_block(fig, "The antisymmetry of the raw Einstein tensor converges",
                  "This is what licenses the symmetric projection.")
    S.footer(fig, wrap(
        "G inherits the O(h^2) Ricci asymmetry that A4 traced to the discrete "
        "product rule: a difference operator does not satisfy D(fg) = (Df)g + "
        "f(Dg). The raw tensor is kept and this residual is MEASURED rather "
        "than projected away -- symmetrising in place would have destroyed the "
        "only direct measure of how far the discretisation is from the "
        "continuum. Probed at fixed integer coordinates, not at the bubble "
        "centre, where the antisymmetry vanishes by symmetry and the measured "
        "order would be 0/0. " + S.FOOTER_UNITS, 120))
    save_both(fig, out, "einstein_antisymmetry_convergence",
              rect=(0, 0.07, 1, 0.89))


def fig_bianchi_convergence(results, out):
    rows = read_csv(os.path.join(results, "bianchi_convergence.csv"))
    pts = sorted((float(r["h"]), float(r["divergence_absolute"]),
                  float(r["divergence_relative"]), float(r["einstein_scale"]))
                 for r in rows)
    hs = [p[0] for p in pts]
    absr = [p[1] for p in pts]
    relr = [p[2] for p in pts]

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), facecolor=S.SURFACE)
    _convergence_axes(axes[0], hs, absr,
                      r"$|\nabla_\mu G^{\mu\nu}|$  (order %.2f)"
                      % order_of(hs, absr), S.SERIES[0])
    axes[0].set_ylabel("absolute residual")
    axes[0].set_title("absolute", fontsize=9, color=S.INK2)

    _convergence_axes(axes[1], hs, relr,
                      r"$|\nabla_\mu G^{\mu\nu}| / \max|G|$  (order %.2f)"
                      % order_of(hs, relr), S.SERIES[2])
    axes[1].set_ylabel("relative residual")
    axes[1].set_title("relative — the comparable one", fontsize=9, color=S.INK2)

    S.title_block(fig, "The contracted Bianchi residual converges",
                  r"$\nabla_\mu G^{\mu\nu} = 0$ is not an exact "
                  r"finite-difference identity.")
    S.footer(fig, wrap(
        "It costs a THIRD nested derivative pass, so at these grids it carries "
        "more pre-asymptotic error than anything else in the project and is "
        "gated with a deliberately wide band. What matters is that it converges "
        "at all: a non-converging Bianchi residual would mean the Einstein "
        "tensor is wrong. Both measures are shown because a vacuum's absolute "
        "residual is tiny purely because its G is tiny -- only the relative "
        "figure is comparable across geometries. " + S.FOOTER_UNITS, 185))
    save_both(fig, out, "contracted_bianchi_convergence", rect=(0, 0.07, 1, 0.89))


# ---------------------------------------------------------------------------
# 6. Raw versus symmetric
# ---------------------------------------------------------------------------

def fig_raw_vs_symmetric(results, out):
    rows = read_csv(os.path.join(results, "alcubierre_einstein_centerline.csv"))
    x = col(rows, "x")

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0), facecolor=S.SURFACE)

    ax = axes[0]
    ax.plot(x, col(rows, "G01_raw"), color=S.SERIES[0], linewidth=2.4,
            alpha=0.5, label="raw $G_{01}$")
    ax.plot(x, col(rows, "G01_symmetric"), color=S.SERIES[1], linewidth=1.2,
            label="symmetric $G_{01}$")
    ax.set_xlabel("x"); ax.set_ylabel("component")
    ax.set_title("indistinguishable at this scale", fontsize=9, color=S.INK2)
    S.style_axes(ax)
    ax.legend(fontsize=7.5, frameon=False)

    ax = axes[1]
    diff = [abs(a - b) for a, b in zip(col(rows, "G01_raw"),
                                       col(rows, "G01_symmetric"))]
    ax.plot(x, diff, color=S.SERIES[2], linewidth=1.4)
    ax.set_xlabel("x"); ax.set_ylabel("|raw − symmetric|")
    ax.set_title("what the projection removes", fontsize=9, color=S.INK2)
    S.style_axes(ax)

    ax = axes[2]
    ax.plot(x, col(rows, "antisymmetry_residual"), color=S.SERIES[1],
            linewidth=1.4)
    ax.set_xlabel("x"); ax.set_ylabel(r"max $|G_{\mu\nu} - G_{\nu\mu}|$")
    ax.set_title("the discarded antisymmetric part", fontsize=9, color=S.INK2)
    S.style_axes(ax)

    for ax in axes:
        for xv in (-2.0, 2.0):
            ax.axvline(xv, color=S.INK3, linewidth=0.8, linestyle="--")

    S.title_block(fig, "Raw and symmetric Einstein tensors are kept separately",
                  "The middle and right panels are what an in-place "
                  "symmetrisation would have made invisible.")
    S.footer(fig, wrap(
        "The raw tensor is the directly computed one; the symmetric projection "
        "(1/2)(G + G^T) is what goes downstream, exactly symmetric by "
        "construction. Both remain reachable, and a verification case fails if "
        "they ever become the same object. The antisymmetric part is largest on "
        "the wall, where the third derivatives of the shape function are "
        "largest -- consistent with an O(h^2) truncation origin rather than a "
        "defect. " + S.FOOTER_UNITS, 210))
    save_both(fig, out, "raw_vs_symmetric_einstein")


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: plot_milestone5.py <results/verification/milestone_05>",
              file=sys.stderr)
        return 2
    results = sys.argv[1]
    out = os.path.join(results, "figures")
    os.makedirs(out, exist_ok=True)

    fig_reference(results, out)
    fig_alcubierre_components(results, out)
    fig_trace_identity(results, out)
    fig_antisymmetry_convergence(results, out)
    fig_bianchi_convergence(results, out)
    fig_raw_vs_symmetric(results, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
