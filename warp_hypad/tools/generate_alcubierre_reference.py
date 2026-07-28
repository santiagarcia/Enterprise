#!/usr/bin/env python3
"""Generate the Milestone A4 Alcubierre curvature reference table.

WHY THIS EXISTS
---------------
Minkowski, flat-spherical, Schwarzschild and de Sitter pin down the curvature
pipeline against closed forms, but none of them is time dependent. The Alcubierre
metric is the one case in the A4 suite where ``d_0 g`` is genuinely nonzero, so
it is the only case that exercises the ``AnalyticResample`` temporal path -- and
there is no textbook closed form for its Riemann tensor to check that path
against.

This script supplies one. It is deliberately NOT a port of the C++ code:

  * the metric derivatives are EXACT, obtained by symbolic differentiation of the
    closed form with sympy, where the C++ pipeline obtains them from second-order
    finite differences on a grid;
  * the tensor algebra is written independently here, in numpy, from the same
    published formulae in docs/conventions.md sections 5-7;
  * the temporal derivative is exact rather than resampled, so agreement is
    evidence about the resampling path rather than a restatement of it.

Agreement between the two is therefore a Level B check with an O(h^2) expected
gap, not a tautology.

WHAT IT DOES NOT PROVE
----------------------
It shares the conventions. If docs/conventions.md fixed the Riemann sign the
other way round, both this script and the C++ would be consistently wrong. The
de Sitter case is what guards that, because its Ricci scalar has a known
published sign.

Usage:
    python tools/generate_alcubierre_reference.py \
        --out tests/reference/alcubierre_curvature_reference.hpp
"""

from __future__ import annotations

import argparse
import datetime
import pathlib

import numpy as np
import sympy as sp

# --- Parameters of the reference configuration -------------------------------
#
# Chosen so the sample points below land exactly on a 33-point grid spanning
# [-4, 4]: that grid's coordinates are -4 + i/4, so every multiple of 0.25 in
# [-4, 4] is an exact grid coordinate with no interpolation anywhere.

V_S = 0.5
RADIUS = 2.0

# sigma = 1 gives a wall of width ~1, which a grid of spacing 0.25 resolves with
# about four points. A steeper wall is more interesting physically and is what
# Program A's sweeps use, but it would put the reference comparison in the regime
# where the second-order truncation error is the dominant term and the tolerance
# would have to be loosened until it stopped testing anything. The convergence
# study is what demonstrates the order; this table is what pins the value.
SIGMA = 1.0
T_EVAL = 0.0

# The bubble centre r_s = 0 is DELIBERATELY ABSENT. The shape function is even in
# r_s, so d_x f = f'(r_s) * x / r_s is a removable 0/0 there: the limit exists and
# is finite, but the symbolic expression evaluates to nan. Rather than paper over
# that with a hand-substituted limit -- which would make the reference agree with
# the C++ by construction at the one point where the C++ has its own special case
# -- the centre is excluded here and checked separately in the A4 suite for
# finiteness and for the expected axial symmetry.
#
# Every coordinate is a multiple of 0.25, which is an exact grid coordinate of the
# 33-point grid spanning [-4, 4] that the A4 suite uses. The two integer-valued
# points are additionally exact on the 9- and 17-point grids, so the convergence
# study can be evaluated at the same physical location on all three.
SAMPLE_POINTS = [
    (0.75, 0.00, 0.00),   # inside the bubble, on the x axis
    (1.00, 0.50, 0.00),   # inside the wall, off axis
    (1.00, 1.00, 0.00),   # convergence sample, off axis
    (1.50, 1.00, 0.50),   # on the wall, fully three-dimensional
    (2.00, 0.00, 0.00),   # wall crest on the x axis, convergence sample
    (2.50, 0.75, -0.50),  # outside the wall
]


def build_symbolic_derivatives():
    """Exact d_a g_{mu nu} and d_b d_a g_{mu nu} as callables of (t, x, y, z)."""
    t, x, y, z = sp.symbols("t x y z", real=True)
    v, R, sig = sp.Rational(str(V_S)), sp.Rational(str(RADIUS)), sp.Rational(str(SIGMA))

    r_s = sp.sqrt((x - v * t) ** 2 + y**2 + z**2)
    f = (sp.tanh(sig * (r_s + R)) - sp.tanh(sig * (r_s - R))) / (2 * sp.tanh(sig * R))
    vf = v * f

    g = sp.zeros(4, 4)
    g[0, 0] = -1 + vf**2
    g[0, 1] = -vf
    g[1, 0] = -vf
    g[1, 1] = 1
    g[2, 2] = 1
    g[3, 3] = 1

    coords = (t, x, y, z)

    dg = [[[sp.diff(g[m, n], coords[a]) for n in range(4)] for m in range(4)]
          for a in range(4)]
    ddg = [[[[sp.diff(dg[a][m][n], coords[b]) for n in range(4)] for m in range(4)]
            for a in range(4)] for b in range(4)]

    g_fn = sp.lambdify(coords, g, "numpy")
    dg_fn = sp.lambdify(coords, sp.Array(dg), "numpy")
    ddg_fn = sp.lambdify(coords, sp.Array(ddg), "numpy")
    return g_fn, dg_fn, ddg_fn


def curvature_at(g_fn, dg_fn, ddg_fn, t, x, y, z):
    """Riemann, Ricci, R and the invariants, assembled independently in numpy.

    Index order matches docs/conventions.md: the leading slot is the raised one.
    """
    g = np.array(g_fn(t, x, y, z), dtype=float).reshape(4, 4)
    dg = np.array(dg_fn(t, x, y, z), dtype=float).reshape(4, 4, 4)      # dg[a][m][n]
    ddg = np.array(ddg_fn(t, x, y, z), dtype=float).reshape(4, 4, 4, 4)  # ddg[b][a][m][n]

    gi = np.linalg.inv(g)

    # d_b g^{rho sigma} = -g^{rho alpha} (d_b g_{alpha beta}) g^{beta sigma}
    dgi = np.array([-gi @ dg[b] @ gi for b in range(4)], dtype=float)

    # Gamma^rho_{mu nu} = 1/2 g^{rho sigma}
    #                     (d_mu g_{sigma nu} + d_nu g_{sigma mu} - d_sigma g_{mu nu})
    D = np.empty((4, 4, 4))  # D[s][m][n]
    for s in range(4):
        for m in range(4):
            for n in range(4):
                D[s, m, n] = dg[m][s][n] + dg[n][s][m] - dg[s][m][n]
    Gamma = 0.5 * np.einsum("rs,smn->rmn", gi, D, optimize=True)

    # d_l D_{sigma mu nu}
    dD = np.empty((4, 4, 4, 4))  # dD[l][s][m][n]
    for l in range(4):
        for s in range(4):
            for m in range(4):
                for n in range(4):
                    dD[l, s, m, n] = (ddg[l][m][s][n] + ddg[l][n][s][m]
                                      - ddg[l][s][m][n])

    dGamma = 0.5 * (np.einsum("lrs,smn->lrmn", dgi, D, optimize=True)
                    + np.einsum("rs,lsmn->lrmn", gi, dD, optimize=True))

    # R^rho_{sigma mu nu} = d_mu Gamma^rho_{nu sigma} - d_nu Gamma^rho_{mu sigma}
    #                     + Gamma^rho_{mu lam} Gamma^lam_{nu sigma}
    #                     - Gamma^rho_{nu lam} Gamma^lam_{mu sigma}
    Riem = np.empty((4, 4, 4, 4))
    for rho in range(4):
        for sig in range(4):
            for mu in range(4):
                for nu in range(4):
                    acc = dGamma[mu][rho][nu][sig] - dGamma[nu][rho][mu][sig]
                    for lam in range(4):
                        acc += (Gamma[rho][mu][lam] * Gamma[lam][nu][sig]
                                - Gamma[rho][nu][lam] * Gamma[lam][mu][sig])
                    Riem[rho, sig, mu, nu] = acc

    Ric = np.einsum("rsrn->sn", Riem, optimize=True)
    R = float(np.einsum("mn,mn->", gi, Ric, optimize=True))

    # G_{mu nu} = R_{mu nu} - (1/2) g_{mu nu} R, assembled here from THIS
    # script's own Ricci rather than read back from the C++ -- the whole point
    # of the reference is that the two routes never touch.
    Ein = Ric - 0.5 * g * R

    Ric_up = gi @ Ric @ gi
    ricci_square = float(np.einsum("mn,mn->", Ric_up, Ric, optimize=True))

    Riem_down = np.einsum("ae,ebcd->abcd", g, Riem, optimize=True)
    Riem_up = np.einsum("ae,bf,cg,dh,efgh->abcd", gi, gi, gi, gi, Riem_down,
                        optimize=True)
    kretschmann = float(np.einsum("abcd,abcd->", Riem_down, Riem_up, optimize=True))

    return Ric, R, ricci_square, kretschmann, Ein


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    g_fn, dg_fn, ddg_fn = build_symbolic_derivatives()

    rows = []
    for (x, y, z) in SAMPLE_POINTS:
        Ric, R, rsq, K, Ein = curvature_at(g_fn, dg_fn, ddg_fn, T_EVAL, x, y, z)
        rows.append((x, y, z, Ric, R, rsq, K, Ein))

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    lines = []
    lines.append("#pragma once")
    lines.append("")
    lines.append("// GENERATED FILE -- do not edit by hand.")
    lines.append("//")
    lines.append("// Produced by tools/generate_alcubierre_reference.py, which")
    lines.append("// differentiates the Alcubierre metric symbolically with sympy and")
    lines.append("// assembles the curvature independently in numpy. See that script's")
    lines.append("// docstring for what this does and does not prove.")
    lines.append("//")
    lines.append(f"// generated: {datetime.date.today().isoformat()}")
    lines.append(f"// sympy:     {sp.__version__}")
    lines.append(f"// numpy:     {np.__version__}")
    lines.append("")
    lines.append("namespace warp_hypad {")
    lines.append("namespace reference {")
    lines.append("")
    lines.append("struct AlcubierreCurvaturePoint {")
    lines.append("    double x, y, z;")
    lines.append("    double ricci_scalar;")
    lines.append("    double ricci_square;")
    lines.append("    double kretschmann;")
    lines.append("    double ricci[16];    // R_{mu nu}, row-major")
    lines.append("    double einstein[16]; // G_{mu nu} = R_{mu nu} - g R / 2")
    lines.append("};")
    lines.append("")
    lines.append(f"inline constexpr double kAlcubierreRefVelocity = {V_S!r};")
    lines.append(f"inline constexpr double kAlcubierreRefRadius = {RADIUS!r};")
    lines.append(f"inline constexpr double kAlcubierreRefSigma = {SIGMA!r};")
    lines.append(f"inline constexpr double kAlcubierreRefTime = {T_EVAL!r};")
    lines.append("")
    lines.append("inline const AlcubierreCurvaturePoint kAlcubierreCurvature[] = {")
    for (x, y, z, Ric, R, rsq, K, Ein) in rows:
        ric_flat = ", ".join(f"{Ric[m, n]:.17e}" for m in range(4) for n in range(4))
        ein_flat = ", ".join(f"{Ein[m, n]:.17e}" for m in range(4) for n in range(4))
        lines.append("    {")
        lines.append(f"        {x!r}, {y!r}, {z!r},")
        lines.append(f"        {R:.17e},")
        lines.append(f"        {rsq:.17e},")
        lines.append(f"        {K:.17e},")
        lines.append(f"        {{{ric_flat}}},")
        lines.append(f"        {{{ein_flat}}},")
        lines.append("    },")
    lines.append("};")
    lines.append("")
    lines.append("inline constexpr int kAlcubierreCurvatureCount = "
                 f"{len(rows)};")
    lines.append("")
    lines.append("}  // namespace reference")
    lines.append("}  // namespace warp_hypad")
    lines.append("")

    out.write_text("\n".join(lines), encoding="utf-8")

    print(f"wrote {out} with {len(rows)} reference points")
    for (x, y, z, _Ric, R, rsq, K, Ein) in rows:
        print(f"  ({x:5.2f},{y:5.2f},{z:5.2f})  R={R: .8e}  "
              f"K={K: .8e}  G00={Ein[0, 0]: .8e}  G11={Ein[1, 1]: .8e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
