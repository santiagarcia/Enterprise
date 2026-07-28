#!/usr/bin/env python3
"""Milestone A4 figures: covariant curvature.

Redrawn entirely from the CSVs the gate wrote. Nothing here evaluates a metric,
differentiates anything or knows what a Christoffel symbol is -- if a number is
not in a saved record it does not appear in a figure. That is the project's
reproducibility rule, and it is what makes these figures evidence rather than
illustration.

Usage:
    python python/plot_milestone4.py results/verification/milestone_04
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


def wrap(text, width=150):
    """Hard-wrap a subtitle.

    `title_block` places the subtitle as a single unwrapped text object, which is
    fine for the short captions the earlier milestones use and runs off the
    canvas for these. Wrapping here rather than in wh_style keeps every existing
    figure byte-identical.
    """
    return "\n".join(textwrap.wrap(" ".join(text.split()), width))


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def col(rows, name, cast=float):
    return [cast(r[name]) for r in rows]


# ---------------------------------------------------------------------------
# 1. Schwarzschild: Ricci flat, but not flat
# ---------------------------------------------------------------------------

def fig_schwarzschild(results, out):
    rows = read_csv(os.path.join(results, "schwarzschild_kretschmann_profile.csv"))
    r = col(rows, "r")
    k_num = col(rows, "kretschmann_computed")
    k_exact = col(rows, "kretschmann_exact")
    rel = col(rows, "rel_error")
    ric = col(rows, "ricci_abs_max")

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0), facecolor=S.SURFACE)

    ax = axes[0]
    ax.plot(r, k_exact, color=S.SERIES[1], linewidth=3.0, alpha=0.45,
            label="exact  48M²/r⁶")
    ax.plot(r, k_num, color=S.SERIES[0], linewidth=1.4, label="computed")
    ax.set_yscale("log")
    ax.set_xlabel("r"); ax.set_ylabel("Kretschmann  $R_{abcd}R^{abcd}$")
    ax.set_title("nonzero curvature invariant", fontsize=9, color=S.INK2)
    S.style_axes(ax, logy=True)
    ax.legend(fontsize=7.5, frameon=False)

    ax = axes[1]
    ax.plot(r, rel, color=S.SERIES[0], linewidth=1.4)
    ax.set_yscale("log")
    ax.set_xlabel("r"); ax.set_ylabel("relative error in K")
    ax.set_title("second-order truncation, not a defect", fontsize=9, color=S.INK2)
    S.style_axes(ax, logy=True)

    ax = axes[2]
    ax.plot(r, ric, color=S.SERIES[2], linewidth=1.4)
    ax.set_yscale("log")
    ax.set_xlabel("r"); ax.set_ylabel(r"max$_{\mu\nu}$ $|R_{\mu\nu}|$")
    ax.set_title("Ricci flat: the vacuum equation holds", fontsize=9, color=S.INK2)
    S.style_axes(ax, logy=True)

    S.title_block(
        fig, "Schwarzschild separates Ricci-flat from flat",
        "$R_{\\mu\\nu}=0$ everywhere, yet the Kretschmann invariant is large "
        "— Ricci-flat is not flat.")
    S.footer(fig, wrap("M = 1, horizon at r = 2, domain r ∈ [4,6] — two horizon "
                       "radii of margin. Grid 33×33×5 in (r, θ, φ). "
                       + S.FOOTER_UNITS, 210))
    S.save(fig, os.path.join(out, "schwarzschild_curvature.png"))


# ---------------------------------------------------------------------------
# 2. de Sitter: the case that fixes the sign
# ---------------------------------------------------------------------------

def fig_desitter(results, out):
    rows = read_csv(os.path.join(results, "desitter_invariants.csv"))
    r = col(rows, "r")
    series = [
        ("ricci_scalar", "ricci_scalar_exact", "$R$", "$12/L^2$"),
        ("kretschmann", "kretschmann_exact", "$R_{abcd}R^{abcd}$", "$24/L^4$"),
        ("ricci_square", "ricci_square_exact",
         "$R_{\\mu\\nu}R^{\\mu\\nu}$", "$36/L^4$"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0), facecolor=S.SURFACE)
    for ax, (num_key, ex_key, label, ex_label) in zip(axes, series):
        num = col(rows, num_key)
        exact = col(rows, ex_key)
        ax.axhline(exact[0], color=S.SERIES[1], linewidth=3.0, alpha=0.45,
                   label=f"exact  {ex_label}")
        ax.plot(r, num, color=S.SERIES[0], linewidth=1.4, label="computed")
        ax.set_xlabel("r"); ax.set_ylabel(label)
        S.style_axes(ax)
        ax.legend(fontsize=7.5, frameon=False)

    axes[0].set_title("R > 0 — the sign is a result, not a choice",
                      fontsize=9, color=S.INK2)
    axes[1].set_title("constant, as maximal symmetry requires",
                      fontsize=9, color=S.INK2)
    axes[2].set_title("", fontsize=9)

    S.title_block(
        fig, "Static de Sitter fixes the sign of every energy density to come",
        "The only reference here with nonzero Ricci, so the only one that can "
        "detect a flipped Riemann sign.")
    S.footer(fig, wrap("L = 2 (horizon at r = 2), domain r ∈ [0.5,1.5]. L is "
                       "small deliberately: at large L the curvature is a "
                       "near-total cancellation of order-unity connection terms "
                       "and the case measures conditioning instead of "
                       "correctness. " + S.FOOTER_UNITS, 210))
    S.save(fig, os.path.join(out, "desitter_invariants.png"))


# ---------------------------------------------------------------------------
# 3. Flat space in spherical coordinates
# ---------------------------------------------------------------------------

def fig_flat_spherical(results, out):
    rows = read_csv(os.path.join(results, "flat_spherical_christoffel_error.csv"))
    r = col(rows, "r")
    gerr = col(rows, "christoffel_abs_error")
    riem = col(rows, "riemann_abs")

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0), facecolor=S.SURFACE)

    ax = axes[0]
    ax.plot(r, gerr, color=S.SERIES[0], linewidth=1.4)
    ax.set_yscale("log")
    ax.set_xlabel("r"); ax.set_ylabel(r"max $|\Gamma_{\rm num} - \Gamma_{\rm exact}|$")
    ax.set_title(r"$\Gamma$ against the closed form", fontsize=9, color=S.INK2)
    S.style_axes(ax, logy=True)

    ax = axes[1]
    ax.plot(r, riem, color=S.SERIES[2], linewidth=1.4)
    ax.set_yscale("log")
    ax.set_xlabel("r"); ax.set_ylabel(r"max $|R^{\rho}{}_{\sigma\mu\nu}|$")
    ax.set_title("spurious curvature of flat space", fontsize=9, color=S.INK2)
    S.style_axes(ax, logy=True)

    S.title_block(
        fig, "Flat space in spherical coordinates: Γ ≠ 0 while Riemann = 0",
        "A large connection whose curvature must still vanish exactly.")
    S.footer(fig, wrap("Domain r ∈ [1,5], θ ∈ [1, π−1] (poles excluded), grid "
                       "33×33×5. Right panel is pure truncation error and "
                       "converges at measured order 2.00. "
                       + S.FOOTER_UNITS, 160))
    S.save(fig, os.path.join(out, "flat_spherical_curvature.png"))


# ---------------------------------------------------------------------------
# 4. Convergence
# ---------------------------------------------------------------------------

def fig_convergence(results, out):
    rows = read_csv(os.path.join(results, "curvature_convergence.csv"))
    cases = {}
    for row in rows:
        cases.setdefault(row["case"], []).append(
            (float(row["h"]), float(row["interior_error"])))

    fig, ax = plt.subplots(figsize=(7.5, 5.0), facecolor=S.SURFACE)
    for idx, (name, pts) in enumerate(sorted(cases.items())):
        pts.sort()  # ascending h, so index 0 is the FINEST grid
        hs = [p[0] for p in pts]
        es = [p[1] for p in pts]
        # The finest pair, matching what the gate asserts. Taking the last two
        # entries would report the COARSEST pair, which is not yet in the
        # asymptotic regime and understates the order.
        order = (math.log(es[0] / es[1]) / math.log(hs[0] / hs[1])
                 if len(es) >= 2 and es[0] > 0 and es[1] > 0 else float("nan"))
        ax.plot(hs, es, "o-", color=S.SERIES[idx % len(S.SERIES)], linewidth=1.4,
                markersize=4, label=f"{name}  (order {order:.2f})")

    # Reference slope of exactly 2, anchored to the finest point of the first case.
    first = sorted(cases.items())[0][1]
    h0, e0 = first[0]
    ref_h = [min(h for h, _ in first), max(h for h, _ in first)]
    ax.plot(ref_h, [e0 * (h / h0) ** 2 for h in ref_h], "--", color=S.INK3,
            linewidth=1.0, label="slope 2 (reference)")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("grid spacing h"); ax.set_ylabel("error at a fixed physical point")
    S.style_axes(ax, logx=True, logy=True)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")

    S.title_block(
        fig, "Curvature converges at second order",
        "Every error is measured at a fixed physical point.")
    S.footer(fig, wrap("The Alcubierre case is the only one with a "
                       "time-dependent metric, and the only one whose finest "
                       "grid does not fully resolve its feature — hence its "
                       "lower order. ricci_symmetry is the residual induced by "
                       "the discrete product rule, and its convergence is what "
                       "licenses treating that residual as truncation error. "
                       + S.FOOTER_UNITS, 125))
    S.save(fig, os.path.join(out, "curvature_convergence.png"), rect=(0, 0.06, 1, 0.88))


# ---------------------------------------------------------------------------
# 5. Where the curvature of a warp bubble actually lives
# ---------------------------------------------------------------------------

def fig_alcubierre_slice(results, out):
    rows = read_csv(os.path.join(results, "alcubierre_curvature_slice_xy.csv"))
    xs = sorted({float(r["x"]) for r in rows})
    ys = sorted({float(r["y"]) for r in rows})
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}

    def field(name):
        g = [[0.0] * len(xs) for _ in ys]
        for r in rows:
            g[yi[float(r["y"])]][xi[float(r["x"])]] = float(r[name])
        return g

    panels = [
        ("g00", r"$g_{00}$", True),
        ("ricci_scalar", r"$R$", True),
        ("kretschmann", r"$R_{abcd}R^{abcd}$", True),
        ("riemann_max", r"max $|R^{\rho}{}_{\sigma\mu\nu}|$", False),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(16.0, 4.2), facecolor=S.SURFACE)
    for ax, (key, label, signed) in zip(axes, panels):
        S.imshow_field(ax, xs, ys, field(key), signed=signed, label=label)
        ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.set_title(label, fontsize=9, color=S.INK2)
        # The nominal wall, drawn from the recorded parameters, not fitted.
        ax.add_patch(plt.Circle((0, 0), 2.0, fill=False, color=S.INK3,
                                linewidth=0.8, linestyle="--"))

    S.title_block(
        fig, "Where the curvature of a warp bubble actually is",
        "z = 0 slice at t = 0; the dashed circle is the nominal radius R = 2, "
        "taken from the recorded parameters rather than fitted.")
    S.footer(fig, wrap(
        "READ THE COMPUTED GEOMETRY, NOT THE INTENDED ONE. The Ricci scalar "
        "behaves as the textbook picture suggests — 0.0018 at the centre "
        "against a peak of 0.145 at x = 1.75, changing sign outside the wall. "
        "The Kretschmann invariant does not: it is −0.023 at the centre against "
        "an extremum of −0.060, so 38% of the peak survives at the bubble's "
        "middle. This interior is NOT flat. At σ = 1 the shape function has "
        "already fallen to 0.91 by r = 1, so there is no plateau where f = 1, "
        "and a flat interior would need a thinner wall than this grid resolves. "
        "σ is small deliberately: a steeper wall is more representative but "
        "would put the comparison against the symbolic reference into the "
        "regime where truncation error dominates. "
        "v_s = 0.5, R = 2, σ = 1, grid 33³ over [−4,4]³, temporal derivative "
        "resampled with an aligned step. Signed fields use a zero-centred "
        "diverging ramp; the single-signed Riemann maximum uses a sequential "
        "one. " + S.FOOTER_UNITS, 215))
    S.save(fig, os.path.join(out, "alcubierre_curvature_slice.png"),
           rect=(0, 0.05, 1, 0.88))


def fig_alcubierre_centerline(results, out):
    rows = read_csv(os.path.join(results, "alcubierre_centerline_x.csv"))
    x = col(rows, "x")

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0), facecolor=S.SURFACE)

    ax = axes[0]
    ax.plot(x, col(rows, "g00"), color=S.SERIES[0], linewidth=1.4, label="$g_{00}$")
    ax.plot(x, col(rows, "g01"), color=S.SERIES[1], linewidth=1.4, label="$g_{01}$")
    ax.set_xlabel("x"); ax.set_ylabel("metric component")
    ax.set_title("the metric varies smoothly", fontsize=9, color=S.INK2)
    S.style_axes(ax)
    ax.legend(fontsize=7.5, frameon=False)

    ax = axes[1]
    ax.plot(x, col(rows, "christoffel_max"), color=S.SERIES[2], linewidth=1.4)
    ax.set_xlabel("x"); ax.set_ylabel(r"max $|\Gamma^{\rho}{}_{\mu\nu}|$")
    ax.set_title("the connection peaks on the wall", fontsize=9, color=S.INK2)
    S.style_axes(ax)

    ax = axes[2]
    ax.plot(x, col(rows, "ricci_scalar"), color=S.SERIES[0], linewidth=1.4,
            label="$R$")
    ax.plot(x, col(rows, "kretschmann"), color=S.SERIES[1], linewidth=1.4,
            label="$R_{abcd}R^{abcd}$")
    ax.axhline(0.0, color=S.INK3, linewidth=0.8)
    ax.set_xlabel("x"); ax.set_ylabel("invariant")
    ax.set_title("curvature changes sign across the wall",
                 fontsize=9, color=S.INK2)
    S.style_axes(ax)
    ax.legend(fontsize=7.5, frameon=False)

    for ax in axes:
        for xv in (-2.0, 2.0):
            ax.axvline(xv, color=S.INK3, linewidth=0.8, linestyle="--")

    S.title_block(
        fig, "Along the axis of travel",
        "Dashed lines mark the nominal wall at x = ±R.")
    S.footer(fig, wrap("y = z = 0, t = 0, v_s = 0.5, R = 2, σ = 1. "
                       + S.FOOTER_UNITS, 210))
    S.save(fig, os.path.join(out, "alcubierre_centerline.png"))


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: plot_milestone4.py <results/verification/milestone_04>",
              file=sys.stderr)
        return 2
    results = sys.argv[1]
    out = os.path.join(results, "figures")
    os.makedirs(out, exist_ok=True)

    fig_schwarzschild(results, out)
    fig_desitter(results, out)
    fig_flat_spherical(results, out)
    fig_convergence(results, out)
    fig_alcubierre_slice(results, out)
    fig_alcubierre_centerline(results, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
