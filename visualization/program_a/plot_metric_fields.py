#!/usr/bin/env python3
"""Program A static figures, from committed milestone_03 records.

    python visualization/program_a/plot_metric_fields.py \
        --results warp_hypad/results/verification/milestone_03 \
        --out visualization/galleries/program_a

Produces:
    alcubierre_shape_function      metric_components_slice
    metric_determinant             metric_signature_map
    ergosurface                    sensitivity_v_slice
    sensitivity_R_slice            sensitivity_sigma_slice
    oti_vs_fd_residual

READS COLUMNS AND DRAWS THEM. Every value was computed by warp_hypad and passed
the milestone_03 gate. Nothing here evaluates a metric, differentiates anything,
or imports a production solver function.

BLOCKED, and deliberately absent: ricci_scalar_slice and curvature_norm_slice.
Program A contains no curvature code -- see visualization/V0_ARCHITECTURE.md.
Computing a Ricci scalar in this script would be exactly the failure the
visualization contract exists to prevent.
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import anomaly, provenance, records
from common import render as rnd
from themes import themes

PROGRAM = "A (warp_hypad)"
MILESTONE = "milestone_03"
SCOPE = ("Effective metric components and parameter sensitivities. No curvature, "
         "no stress-energy, no energy conditions: Program A does not compute "
         "them. Not a claim of measurable gravitational curvature.")


def grid_from(xs, ys, vals):
    """Reshape a saved (x, y, value) table into a 2-D array for imshow.

    PRESENTATION ONLY: this is a reindexing of saved samples, and it asserts the
    table really is a complete rectangular grid rather than silently filling
    gaps. A missing sample becomes an error, not an invisible hole.
    """
    ux = sorted(set(xs))
    uy = sorted(set(ys))
    ix = {v: i for i, v in enumerate(ux)}
    iy = {v: i for i, v in enumerate(uy)}
    grid = [[None] * len(ux) for _ in range(len(uy))]
    for x, y, v in zip(xs, ys, vals):
        grid[iy[y]][ix[x]] = v
    holes = sum(1 for row in grid for v in row if v is None)
    if holes:
        raise records.RecordError(
            "the saved slice is not a complete rectangular grid: %d of %d cells "
            "have no sample. Interpolating across them would invent data."
            % (holes, len(ux) * len(uy)))
    return ux, uy, grid


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results", required=True,
                    help="the milestone_03 results directory")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--theme", default="scientific",
                    choices=["scientific", "scientific_print", "cinematic"])
    a = ap.parse_args(argv)

    res = records.resolve(a.results)
    out_dir = records.ensure_dir(a.out)
    theme = themes.get(a.theme)
    themes.apply(theme)
    rep = anomaly.AnomalyReport(PROGRAM, MILESTONE, "plot_metric_fields")

    f_slice = os.path.join(res, "metric_slice_xy.csv")
    f_shape = os.path.join(res, "shape_profiles.csv")
    f_sens = os.path.join(res, "sensitivity_slice_xy.csv")
    f_ofd = os.path.join(res, "hypad_vs_fd.csv")

    slice_rows = records.read_csv(f_slice)
    x = records.column(slice_rows, "x")
    y = records.column(slice_rows, "y")
    g00 = records.column(slice_rows, "g00")
    g01 = records.column(slice_rows, "g01")
    g11 = records.column(slice_rows, "g11")
    det = records.column(slice_rows, "det")
    sig_ok = records.column(slice_rows, "signature_ok")

    # --- Anomalies on the metric slice, before anything is drawn ------------
    rep.check_finite("metric_slice g00", g00)
    rep.check_finite("metric_slice det", det)
    rep.check_not_blank("metric_slice g00", g00)
    # A Lorentzian metric in four dimensions has det(g) < 0 EVERYWHERE. A
    # positive determinant is not an ugly plot, it is a broken spacetime.
    rep.check_sign("metric determinant", det, expect_negative=True)
    rep.check_flag_column("saved signature_ok flag", sig_ok)

    def base_prov(color_model, camera, inputs, extra=None):
        return provenance.build_provenance(
            program=PROGRAM, milestone=MILESTONE, input_files=inputs,
            simulation_parameters={
                "note": "as recorded in the milestone_03 gate metadata",
                "slice": "xy plane through the bubble",
            },
            coordinate_frame="inertial Cartesian, comoving bubble coordinates",
            observer="coordinate (no observer boost applied)",
            units="geometrised: c = 1, lengths in bubble-radius units",
            time_convention="static slice at the gate's evaluation time",
            render_theme=theme.name, color_model=color_model,
            camera_model=camera, frame_count=1, extra=extra)

    def save(fig, name, prov, plotted=None, header=None):
        for ext in ("png", "svg"):
            p = os.path.join(out_dir, "%s.%s" % (name, ext))
            fig.savefig(p, facecolor=fig.get_facecolor())
            provenance.write_sidecar(p, prov)
        plt.close(fig)
        if plotted is not None:
            records.write_plotted_csv(
                os.path.join(out_dir, "%s.plotted.csv" % name), header, plotted)
        print("wrote %s.png / .svg" % name)

    ux, uy, G00 = grid_from(x, y, g00)
    extent = [min(ux), max(ux), min(uy), max(uy)]

    def field_fig(grid, title, cbar_label, cmap, name, prov, symmetric=False):
        fig, ax = plt.subplots(figsize=(7.6, 6.2))
        flat = [v for row in grid for v in row]
        if symmetric:
            m = max(abs(min(flat)), abs(max(flat))) or 1.0
            vmin, vmax = -m, m
        else:
            vmin, vmax = min(flat), max(flat)
        im = ax.imshow(grid, origin="lower", extent=extent, cmap=cmap,
                       vmin=vmin, vmax=vmax, aspect="equal")
        # EQUAL SPATIAL SCALING: x and y are the same physical quantity, so a
        # stretched aspect would misrepresent the bubble's shape.
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(cbar_label, color=theme.foreground)
        cb.ax.tick_params(colors=theme.foreground)
        ax.set_xlabel("x [bubble radii]")
        ax.set_ylabel("y [bubble radii]")
        ax.set_title(title, fontsize=10)
        rep.check_clipped(name, flat, vmin, vmax)
        rnd.stamp(ax, theme, lines=[
            "program     %s" % PROGRAM,
            "milestone   %s" % MILESTONE,
            "grid        %d x %d" % (len(ux), len(uy)),
            "range       %.4g to %.4g" % (min(flat), max(flat)),
            "aspect      equal (x and y are the same quantity)",
        ], loc="upper left")
        rnd.scope_note(fig, theme, SCOPE)
        save(fig, name, prov)
        return fig

    # --- 1. Shape function ---------------------------------------------------
    shape_rows = records.read_csv(f_shape)
    sr = records.column(shape_rows, "r")
    sf = records.column(shape_rows, "f")
    sR = records.column(shape_rows, "R")
    ssig = records.column(shape_rows, "sigma")
    rep.check_finite("shape f", sf)
    # f must fall monotonically from ~1 at the centre to ~0 outside; a step
    # larger than this means the sampling missed the wall entirely.
    families = {}
    for i in range(len(sr)):
        families.setdefault((sR[i], ssig[i]), []).append((sr[i], sf[i]))
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    for (R, s), pts in sorted(families.items()):
        pts.sort()
        rr = [p[0] for p in pts]
        ff = [p[1] for p in pts]
        rep.check_no_discontinuity("shape f (R=%g, sigma=%g)" % (R, s), ff, 0.35)
        ax.plot(rr, ff, lw=1.8, label="R = %g, sigma = %g" % (R, s))
    ax.set_xlabel("r [bubble radii]")
    ax.set_ylabel("f(r)")
    ax.set_title("Alcubierre shape function — %d parameter families"
                 % len(families), fontsize=10)
    ax.grid(True, alpha=theme.grid_alpha)
    ax.legend(fontsize=7, facecolor=theme.background, edgecolor=theme.grid,
              labelcolor=theme.foreground)
    rnd.stamp(ax, theme, lines=[
        "program     %s" % PROGRAM,
        "families    %d" % len(families),
        "f -> 1 inside the bubble, 0 outside",
        "larger sigma is a thinner wall, not a different radius",
    ], loc="upper right")
    rnd.scope_note(fig, theme, SCOPE)
    save(fig, "alcubierre_shape_function",
         base_prov("illustrative_only", "none (1-D profile)", [f_shape]),
         [[sR[i], ssig[i], sr[i], sf[i]] for i in range(len(sr))],
         ["R", "sigma", "r", "f"])

    # --- 2. Metric components ------------------------------------------------
    _, _, G01 = grid_from(x, y, g01)
    _, _, G11 = grid_from(x, y, g11)
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0))
    for ax, grid, lab in ((axes[0], G00, "$g_{00}$"), (axes[1], G01, "$g_{01}$"),
                          (axes[2], G11, "$g_{11}$")):
        flat = [v for row in grid for v in row]
        m = max(abs(min(flat)), abs(max(flat))) or 1.0
        im = ax.imshow(grid, origin="lower", extent=extent, cmap="RdBu_r",
                       vmin=-m, vmax=m, aspect="equal")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title(lab, fontsize=11)
        ax.set_xlabel("x [bubble radii]")
        axes[0].set_ylabel("y [bubble radii]")
    fig.suptitle("Effective metric components on the xy slice — %s" % MILESTONE,
                 fontsize=11)
    rnd.scope_note(fig, theme, SCOPE)
    save(fig, "metric_components_slice",
         base_prov("illustrative_only", "orthographic, fixed", [f_slice]),
         [[x[i], y[i], g00[i], g01[i], g11[i]] for i in range(len(x))],
         ["x", "y", "g00", "g01", "g11"])

    # --- 3. Determinant and 4. signature ------------------------------------
    _, _, DET = grid_from(x, y, det)
    # det(g) = -1 EXACTLY for the Alcubierre metric, everywhere, analytically.
    # A flat image here is the correct answer rather than a dull one, so the
    # figure states the analytic value and the deviation from it instead of
    # leaving a reader to wonder whether the field failed to populate.
    worst_det = max(abs(v + 1.0) for v in det)
    field_fig(DET, "Metric determinant: analytically hBc1$ everywhere",
              "det(g)", "magma", "metric_determinant",
              base_prov("illustrative_only", "orthographic, fixed", [f_slice],
                        extra={"analytic_determinant": -1.0,
                               "worst_deviation": worst_det}))
    rep.record("determinant_matches_analytic", fired=worst_det > 1.0e-12,
               measured=worst_det, threshold=1.0e-12,
               detail="|det(g) + 1| over the slice; the Alcubierre metric has "
                      "unit determinant by construction")

    _, _, SIG = grid_from(x, y, sig_ok)
    field_fig(SIG, "Lorentzian signature check — 1 where the saved flag holds",
              "signature_ok", "Greens", "metric_signature_map",
              base_prov("illustrative_only", "orthographic, fixed", [f_slice]))

    # --- 5. Ergosurface ------------------------------------------------------
    # The g00 = 0 contour. Drawn with matplotlib's contour of the SAVED g00
    # field: locating a level set of saved samples is presentation, and the
    # surface is not recomputed from any metric formula.
    fig, ax = plt.subplots(figsize=(7.6, 6.2))
    m = max(abs(min(g00)), abs(max(g00))) or 1.0
    im = ax.imshow(G00, origin="lower", extent=extent, cmap="RdBu_r",
                   vmin=-m, vmax=m, aspect="equal")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label("$g_{00}$")
    cs = ax.contour(ux, uy, G00, levels=[0.0], colors=[theme.accent_two],
                    linewidths=2.0)
    ax.clabel(cs, fmt={0.0: "$g_{00}=0$"}, fontsize=8)
    crosses = min(g00) < 0.0 < max(g00)
    ax.set_xlabel("x [bubble radii]")
    ax.set_ylabel("y [bubble radii]")
    # The title states what the samples SHOW, not what the figure was named
    # after. Calling a figure "ergosurface" when g00 never reaches zero would
    # assert a surface that is not there.
    ax.set_title("$g_{00}=0$ ergosurface diagnostic" if crosses
                 else "$g_{00}$ diagnostic — no ergosurface for this "
                      "parameter set", fontsize=10)
    rnd.stamp(ax, theme, lines=[
        "g00 range   %.4g to %.4g" % (min(g00), max(g00)),
        "crosses 0   %s" % ("yes, an ergosurface exists" if crosses
                            else "NO -- no ergosurface in this slice"),
        "this is a coordinate statement about g00,",
        "not a horizon and not a curvature claim",
    ], loc="upper left")
    rnd.scope_note(fig, theme, SCOPE)
    save(fig, "ergosurface",
         base_prov("illustrative_only", "orthographic, fixed", [f_slice],
                   extra={"ergosurface_present": crosses}))

    # --- 6-8. Sensitivities --------------------------------------------------
    sens = records.read_csv(f_sens)
    sx = records.column(sens, "x")
    sy = records.column(sens, "y")
    for param, col in (("v", "dg00_dv"), ("R", "dg00_dR"),
                       ("sigma", "dg00_dsigma")):
        vals = records.column(sens, col)
        rep.check_finite("sensitivity %s" % param, vals)
        rep.check_not_blank("sensitivity %s" % param, vals)
        _, _, S = grid_from(sx, sy, vals)
        field_fig(S, r"Sensitivity $\partial g_{00}/\partial %s$"
                  % ("\\sigma" if param == "sigma" else param),
                  r"$\partial g_{00}/\partial %s$"
                  % ("\\sigma" if param == "sigma" else param),
                  "RdBu_r", "sensitivity_%s_slice" % param,
                  base_prov("illustrative_only", "orthographic, fixed", [f_sens]),
                  symmetric=True)

    # --- 9. OTI versus finite differences ------------------------------------
    ofd = records.read_csv(f_ofd)
    steps = records.column(ofd, "step_size")
    abserr = records.column(ofd, "absolute_error")
    params = records.text_column(ofd, "parameter")
    rep.check_finite("oti_vs_fd absolute_error", abserr)
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    for p in sorted(set(params)):
        pts = sorted((steps[i], abserr[i]) for i in range(len(steps))
                     if params[i] == p)
        ax.loglog([q[0] for q in pts], [max(q[1], 1e-18) for q in pts],
                  marker="o", ms=3.5, lw=1.5, label=p)
    ax.set_xlabel("finite-difference step size")
    ax.set_ylabel("|OTI - FD|")
    ax.set_title("Automatic differentiation against finite differences",
                 fontsize=10)
    ax.grid(True, alpha=theme.grid_alpha, which="both")
    ax.legend(fontsize=7.5, facecolor=theme.background, edgecolor=theme.grid,
              labelcolor=theme.foreground)
    rnd.stamp(ax, theme, lines=[
        "worst |OTI - FD|  %.3e" % max(abserr),
        "the V shape is expected: truncation error falls",
        "with h, rounding error rises, and OTI has neither",
    ], loc="lower left")
    rnd.scope_note(fig, theme, SCOPE)
    save(fig, "oti_vs_fd_residual",
         base_prov("illustrative_only", "none (convergence plot)", [f_ofd]),
         [[params[i], steps[i], abserr[i]] for i in range(len(steps))],
         ["parameter", "step_size", "absolute_error"])

    # --- 10. Deviation views -------------------------------------------------
    #
    # g00 spans -1 to -0.96 and g11 is exactly 1: on an absolute scale both are
    # flat, and the structure that matters is invisible. A residual against a
    # STATED reference makes it visible -- but a contrast-enhanced residual must
    # never be presented as though it were the absolute field, so every panel
    # says so in its title and carries the reference value.
    deviations = [
        ("g00", g00, -1.0, r"$\Delta g_{00} = g_{00} + 1$", "reference $-1$ "
         "(Minkowski)"),
        ("g11", g11, 1.0, r"$\Delta g_{11} = g_{11} - 1$", "reference $+1$ "
         "(Minkowski)"),
        ("det", det, -1.0, r"$\Delta \det g = \det g + 1$",
         "reference $-1$ (analytic, exact for any $v$ and $f$)"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.0))
    dev_rows = []
    for ax, (nm, vals, ref, lab, refnote) in zip(axes, deviations):
        dev = [v - ref for v in vals]
        _, _, D = grid_from(x, y, dev)
        m = max(abs(min(dev)), abs(max(dev)))
        if m == 0.0:
            m = 1.0e-16      # an exactly-zero residual still needs a finite scale
        im = ax.imshow(D, origin="lower", extent=extent, cmap="PuOr_r",
                       vmin=-m, vmax=m, aspect="equal")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.formatter.set_powerlimits((-2, 2))
        ax.set_title("%s\n%s" % (lab, refnote), fontsize=9.5)
        ax.set_xlabel("x [bubble radii]")
        dev_rows.append([nm, ref, min(dev), max(dev),
                         max(abs(v) for v in dev)])
    axes[0].set_ylabel("y [bubble radii]")
    fig.suptitle("CONTRAST-ENHANCED DEVIATION FROM STATED REFERENCE — "
                 "these are residuals, not the absolute fields", fontsize=11)
    rnd.scope_note(fig, theme, SCOPE)
    save(fig, "metric_deviation_views",
         base_prov("illustrative_only", "orthographic, fixed", [f_slice],
                   extra={"view_kind": "contrast_enhanced_deviation",
                          "references": {"g00": -1.0, "g11": 1.0, "det": -1.0},
                          "is_absolute_field": False}),
         dev_rows,
         ["field", "reference", "min_deviation", "max_deviation",
          "worst_absolute_deviation"])

    # --- 11. Constant-field status panel -------------------------------------
    #
    # det(g) = -1 exactly, everywhere, for ANY v and f. A flat heatmap conveys
    # nothing; the useful statement is the comparison against the analytic value
    # and whether it holds.
    det_min, det_max = min(det), max(det)
    det_worst = max(abs(v + 1.0) for v in det)
    det_tol = 1.0e-12
    det_pass = det_worst <= det_tol
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    ax.axis("off")
    rows = [
        ("analytic value", "-1 exactly, for any $v$ and $f$"),
        ("measured minimum", "%.17g" % det_min),
        ("measured maximum", "%.17g" % det_max),
        ("worst |det g + 1|", "%.3e" % det_worst),
        ("tolerance", "%.1e" % det_tol),
        ("result", "PASS" if det_pass else "FAIL"),
    ]
    for k, (lab, val) in enumerate(rows):
        ax.text(0.03, 0.88 - 0.15 * k, lab, fontsize=10, family="monospace",
                color=theme.foreground, transform=ax.transAxes)
        ax.text(0.46, 0.88 - 0.15 * k, val, fontsize=10, family="monospace",
                color=(theme.accent if lab != "result"
                       else (theme.accent_two if det_pass
                             else theme.uncertainty_color)),
                transform=ax.transAxes)
    ax.set_title("det(g) status — an analytically constant field, reported "
                 "rather than shaded", fontsize=10.5, color=theme.foreground)
    rnd.scope_note(fig, theme,
                   "The (t,x) block determinant is (-1 + v^2 f^2) - (v f)^2 = -1 "
                   "identically, so this is an exact identity, not a tolerance "
                   "result. " + SCOPE)
    save(fig, "metric_determinant_status",
         base_prov("illustrative_only", "none (status panel)", [f_slice],
                   extra={"analytic_value": -1.0, "worst_deviation": det_worst,
                          "tolerance": det_tol, "passed": det_pass}),
         [[det_min, det_max, det_worst, det_tol, int(det_pass)]],
         ["measured_min", "measured_max", "worst_absolute_deviation",
          "tolerance", "passed"])

    path = rep.write(out_dir)
    fired = rep.fired
    print("\nanomaly report: %s" % path)
    print("  %d diagnostics, %d fired" % (len(rep.entries), len(fired)))
    for e in fired:
        print("  FIRED %s: measured %s against %s (%s)"
              % (e["diagnostic"], e["measured"], e["threshold"], e["detail"]))

    # BLOCKED products, recorded so they are visibly blocked rather than
    # silently absent.
    with open(os.path.join(out_dir, "BLOCKED.md"), "w", encoding="utf-8") as fh:
        fh.write("# Blocked V0-A products\n\n"
                 "`ricci_scalar_slice`, `curvature_norm_slice` and\n"
                 "`curvature_wall_evolution` are NOT generated.\n\n"
                 "Program A contains no curvature code. `warp_hypad` computes\n"
                 "metric components, determinants, signature masks and parameter\n"
                 "sensitivities; the word 'curvature' appears in it only as a\n"
                 "guard forbidding the frozen-time mode from being fed to a\n"
                 "Christoffel routine that does not yet exist.\n\n"
                 "Computing a Ricci scalar inside this plotting script would be\n"
                 "exactly the failure the visualization contract prevents: a\n"
                 "field invented in Python, rendered as a figure, and thereafter\n"
                 "believed. It would look identical to a verified one.\n\n"
                 "Unblocking requires a Program A scientific milestone:\n"
                 "Christoffel symbols from the verified finite-difference\n"
                 "operators, Riemann and Ricci assembly, and verification\n"
                 "against known answers -- flat space exactly zero, Schwarzschild\n"
                 "R = 0 with a nonzero Kretschmann invariant, and a convergence\n"
                 "study. See visualization/V0_ARCHITECTURE.md.\n")
    return 1 if fired else 0


if __name__ == "__main__":
    sys.exit(main())
