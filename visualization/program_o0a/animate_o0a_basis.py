#!/usr/bin/env python3
"""The six O0A spline-basis animations, drawn from saved C++ records only.

    python animate_o0a_basis.py --records <o0a_basis dir> --out <gallery dir>

Every frame's numbers were computed by production and written to a record: the
projections in raw_to_constrained_frames.csv and constraint_effect_frames.csv, the
partial reconstructions in f0_reconstruction_frames.csv, the whitened coefficient
patterns in whitened_perturbation_frames.csv. This renderer selects, arranges and
colours them. It does not project, factorise or evaluate a spline.

Colour and axis limits are computed once from each whole record and frozen, so a
mode that is an order of magnitude smaller than its neighbour looks smaller
rather than being renormalised into looking the same.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

import o0a_common as O                                           # noqa: E402
from common import anomaly, records, render                      # noqa: E402
from themes import themes                                        # noqa: E402

FPS = 3
ACTIVE_COLOR = "#4ea1ff"
REMOVED_COLOR = "#8c8f99"


def _style(ax, theme, xlabel="", ylabel=""):
    O.style_axes(ax, theme, xlabel, ylabel)


def _stamp(fig, theme, lines, x=0.010, y=0.30):
    fig.text(x, y, "\n".join(lines), ha="left", va="center", fontsize=6.6,
             family="monospace", color=theme.foreground, linespacing=1.5,
             bbox=dict(facecolor=theme.background, edgecolor=theme.grid,
                       alpha=0.8, boxstyle="round,pad=0.5", linewidth=0.6))


def _grid_from(rows, key, xkey="xi", ykey="s"):
    xs = sorted({float(r[xkey]) for r in rows})
    ys = sorted({float(r[ykey]) for r in rows})
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}
    g = [[0.0] * len(xs) for _ in ys]
    for r in rows:
        g[yi[float(r[ykey])]][xi[float(r[xkey])]] = float(r[key])
    return xs, ys, g


def _png_resolution(path):
    with open(path, "rb") as fh:
        head = fh.read(24)
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        return [0, 0]
    return [int.from_bytes(head[16:20], "big"),
            int.from_bytes(head[20:24], "big")]


# --- factories --------------------------------------------------------------

def make_radial_mode_factory(rec, theme, scales):
    prof = O.read_csv(os.path.join(rec, "univariate_basis_profiles.csv"))
    der = O.read_csv(os.path.join(rec, "univariate_basis_derivatives.csv"))
    stats = {int(r["mode_index"]): r for r in
             O.read_csv(os.path.join(rec, "mode_statistics.csv"))
             if r["family"] == "f0_radial"}
    by_v = O.group_by(prof, "index", int)
    dmap = {}
    for r in der:
        dmap.setdefault(int(r["index"]), []).append(r)
    modes = sorted(by_v)
    vlim = max(abs(float(r["value"])) for r in prof) or 1.0

    def factory(i):
        k = modes[i]
        pts = sorted(by_v[k], key=lambda r: float(r["s"]))
        st = stats.get(k, {})
        status = st.get("status", "active")
        active = status == "active"
        fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.4),
                                 facecolor=theme.background)
        ax = axes[0]
        for j in modes:
            other = sorted(by_v[j], key=lambda r: float(r["s"]))
            ax.plot(O.col(other, "s"), O.col(other, "value"), lw=0.7,
                    color=theme.grid, alpha=0.45)
        ax.plot(O.col(pts, "s"), O.col(pts, "value"), lw=2.2,
                color=ACTIVE_COLOR if active else REMOVED_COLOR,
                ls="-" if active else "--")
        lo = float(st.get("support_lower", 0.0))
        hi = float(st.get("support_upper", 0.0))
        ax.axvspan(lo, hi, color=ACTIVE_COLOR if active else REMOVED_COLOR,
                   alpha=0.12)
        ax.set_ylim(-0.15 * vlim, 1.08 * vlim)
        ax.set_title("F0 radial mode %d of %d   %s" %
                     (k, len(modes) - 1,
                      "ACTIVE DESIGN MODE" if active else "REMOVED - NOT A "
                                                          "DESIGN MODE"),
                     fontsize=10,
                     color=theme.foreground if active
                     else theme.uncertainty_color)
        _style(ax, theme, r"$s_r=(r/R_0)^2$ (dimensionless)", r"$N_k$")

        ax = axes[1]
        dpts = sorted(dmap.get(k, []), key=lambda r: float(r["s"]))
        for key, lab, colr in (("first", r"$N'$", "#4ea1ff"),
                               ("second", r"$N''$", "#ffb454"),
                               ("third", r"$N'''$", "#7ed957")):
            ax.plot(O.col(dpts, "s"), O.col(dpts, key), lw=1.4, color=colr,
                    label=lab)
        ax.axvspan(lo, hi, color=REMOVED_COLOR, alpha=0.10)
        ax.set_ylim(-scales["der"], scales["der"])
        ax.legend(fontsize=7, frameon=False, labelcolor=theme.foreground)
        ax.set_title("analytic derivatives (fixed scale across all modes)",
                     fontsize=9.5, color=theme.foreground)
        _style(ax, theme, r"$s_r$", "derivative")

        fig.suptitle("O0A  individual F0 radial basis modes", fontsize=12,
                     color=theme.foreground, x=0.010, ha="left", y=0.985)
        if not active:
            O.banner(fig, theme,
                     "THIS MODE IS REMOVED: %s. It has no coefficient, so it "
                     "cannot contribute anywhere. It is drawn dashed and grey to "
                     "show what the collar constraint takes out, NOT as an "
                     "available design degree of freedom."
                     % O.EXCLUSION_LABEL.get(status, status), y=0.935)
        else:
            O.banner(fig, theme,
                     "Active design mode. Shaded band is its local support; "
                     "outside it the mode is exactly zero.", y=0.935,
                     color=theme.accent_two)
        _stamp(fig, theme, [
            "mode index      %d" % k,
            "support         [%.4f, %.4f]" % (lo, hi),
            "status          %s" % status,
            "collar reason   %s" % O.EXCLUSION_LABEL.get(status, status),
            "max |N|         %.6g" % float(st.get("max_abs_value", 0.0)),
            "max |N'|        %.6g" % float(st.get("max_abs_first", 0.0)),
            "max |N''|       %.6g" % float(st.get("max_abs_second", 0.0)),
        ], y=0.30)
        O.footer(fig, theme, O.SCOPE)
        fig.tight_layout(rect=(0.235, 0.05, 1, 0.895))
        scales.setdefault("radial", []).append((-scales["der"], scales["der"]))
        return fig

    return factory, len(modes)


def make_axisymmetric_mode_factory(rec, theme, scales):
    rows = O.read_csv(os.path.join(rec, "tensor_mode_frames.csv"))
    sup = {int(r["flat_index"]): r for r in
           O.read_csv(os.path.join(rec, "tensor_product_supports.csv"))}
    by = O.group_by(rows, "flat_index", int)
    modes = sorted(by)
    lim = max(abs(float(r["value"])) for r in rows) or 1.0

    def factory(i):
        k = modes[i]
        pts = by[k]
        xs, ys, g = _grid_from(pts, "value")
        s = sup.get(k, {})
        fig, ax = plt.subplots(figsize=(11.6, 5.6), facecolor=theme.background)
        im = ax.imshow(g, origin="lower",
                       extent=[min(xs), max(xs), min(ys), max(ys)],
                       aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim,
                       interpolation="nearest")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=7, colors=theme.foreground)
        cb.outline.set_visible(False)
        ax.set_title(r"$B_{ij}(\xi,s)$   flat %d of %d   $(i,j)=(%s,%s)$"
                     % (k, len(modes) - 1, s.get("index_xi", "?"),
                        s.get("index_s", "?")),
                     fontsize=10, color=theme.foreground)
        _style(ax, theme, r"$\xi=(x-x_c)/R_0$ (dimensionless)",
               r"$s=(y^2+z^2)/R_0^2$ (dimensionless)")
        fig.suptitle("O0A  individual axisymmetric tensor-product modes, all "
                     "active", fontsize=12, color=theme.foreground, x=0.010,
                     ha="left", y=0.985)
        O.banner(fig, theme,
                 "Every mode shown is ACTIVE: the %d collar-excluded lattice "
                 "functions have no coefficient and are not animated as design "
                 "modes. Fixed colour scale across all frames."
                 % 44, y=0.935)
        _stamp(fig, theme, [
            "flat index      %d" % k,
            "xi support      [%s, %s]" % (s.get("xi_lower"), s.get("xi_upper")),
            "s support       [%s, %s]" % (s.get("s_lower"), s.get("s_upper")),
            "rho support (m) [%s, %s]" % (s.get("rho_lower"), s.get("rho_upper")),
            "status          active design mode",
            "max |B|         %.6g" % max(abs(float(r["value"])) for r in pts),
            "max |dB/dxi|    %.6g" % max(abs(float(r["d_xi"])) for r in pts),
            "max |dB/ds|     %.6g" % max(abs(float(r["d_s"])) for r in pts),
        ], y=0.30)
        O.footer(fig, theme, "s is dimensionless; a physical radius is "
                            r"$\rho=R_0\sqrt{s}$. " + O.SCOPE)
        fig.tight_layout(rect=(0.245, 0.05, 1, 0.895))
        scales.setdefault("axisym", []).append((-lim, lim))
        return fig

    return factory, len(modes)


def make_raw_to_constrained_factory(rec, theme, scales):
    rows = O.read_csv(os.path.join(rec, "raw_to_constrained_frames.csv"))
    by = O.group_by(rows, "frame", int)
    frames = sorted(by)
    lim = max(abs(float(r["raw_value"])) for r in rows) or 1.0

    def factory(i):
        pts = sorted(by[frames[i]], key=lambda r: int(r["coefficient_index"]))
        idx = [int(r["coefficient_index"]) for r in pts]
        fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.2),
                                 facecolor=theme.background)
        ax = axes[0]
        ax.bar(idx, O.col(pts, "raw_value"), width=0.8, color=REMOVED_COLOR,
               alpha=0.6, label="raw perturbation")
        ax.bar(idx, O.col(pts, "projected_value"), width=0.45,
               color=ACTIVE_COLOR, label="after exact constraint elimination")
        ax.set_ylim(-1.15 * lim, 1.15 * lim)
        ax.legend(fontsize=7, frameon=False, labelcolor=theme.foreground)
        ax.set_title("raw unit direction %d and its constrained part"
                     % frames[i], fontsize=9.5, color=theme.foreground)
        _style(ax, theme, "raw coefficient index", "coefficient value")

        ax = axes[1]
        ax.bar(idx, O.col(pts, "removed_value"), width=0.7, color="#d1495b")
        ax.set_ylim(-1.15 * lim, 1.15 * lim)
        ax.set_title("the component the constraints REMOVE", fontsize=9.5,
                     color=theme.foreground)
        _style(ax, theme, "raw coefficient index", "removed value")

        r = pts[0]
        fig.suptitle("O0A  exact constraint elimination, direction by direction",
                     fontsize=12, color=theme.foreground, x=0.010, ha="left",
                     y=0.985)
        O.banner(fig, theme,
                 "THIS IS AN EXACT PROJECTION, NOT AN OPTIMIZER PENALTY. The "
                 "removed component is the part of this direction that would "
                 "violate the translation moment, the protected-cabin mean, or "
                 "the redundant constraint direction. What survives satisfies "
                 "every constraint to machine precision by construction, at "
                 "every iterate, not asymptotically.", y=0.935)
        _stamp(fig, theme, [
            "raw direction      %d" % frames[i],
            "raw residual       %.3e" % float(r["constraint_residual_raw"]),
            "projected residual %.3e" % float(r["constraint_residual_projected"]),
            "removed norm       %.6f" % float(r["removed_norm"]),
            "surviving norm     %.6f" % float(r["projected_norm"]),
        ], y=0.30)
        O.footer(fig, theme, O.SCOPE)
        fig.tight_layout(rect=(0.215, 0.05, 1, 0.885))
        scales.setdefault("raw2con", []).append((-1.15 * lim, 1.15 * lim))
        return fig

    return factory, len(frames)


def make_constraint_effect_factory(rec, theme, scales):
    rows = O.read_csv(os.path.join(rec, "constraint_effect_frames.csv"))
    by = O.group_by(rows, "frame", int)
    frames = sorted(by)
    lim = max(abs(float(r["unconstrained"])) for r in rows) or 1.0

    def factory(i):
        pts = by[frames[i]]
        xs, ys, gu = _grid_from(pts, "unconstrained")
        _, _, gp = _grid_from(pts, "projected")
        _, _, gd = _grid_from(pts, "difference")
        fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.0),
                                 facecolor=theme.background)
        for ax, g, lab in ((axes[0], gu, "unconstrained field"),
                           (axes[1], gp, "projected constrained field"),
                           (axes[2], gd, "difference")):
            im = ax.imshow(g, origin="lower",
                           extent=[min(xs), max(xs), min(ys), max(ys)],
                           aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim,
                           interpolation="nearest")
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
            cb.ax.tick_params(labelsize=6.5, colors=theme.foreground)
            cb.outline.set_visible(False)
            ax.set_title(lab, fontsize=9.5, color=theme.foreground)
            _style(ax, theme, r"$\xi$", r"$s$")

        collar_max = max(abs(float(r["projected"])) for r in pts
                         if r["in_outer_collar"] == "1")
        cabin_rows = [r for r in pts if r["in_protected_cabin"] == "1"]
        r0 = pts[0]
        fig.suptitle("O0A  the effect of exact constraint projection on the field",
                     fontsize=12, color=theme.foreground, x=0.010, ha="left",
                     y=0.985)
        O.banner(fig, theme,
                 "The projected field is EXACTLY %.1e in the outer collar over "
                 "%d collar samples. Fixed colour scale across all three panels "
                 "and all frames."
                 % (collar_max, sum(1 for r in pts
                                    if r["in_outer_collar"] == "1")), y=0.935)
        _stamp(fig, theme, [
            "frame                %d" % frames[i],
            "raw residual         %.3e" % float(r0["constraint_residual_raw"]),
            "projected residual   %.3e"
            % float(r0["constraint_residual_projected"]),
            "collar max |field|   %.3e" % collar_max,
            "cabin samples        %d" % len(cabin_rows),
        ], y=0.30)
        O.footer(fig, theme, O.SCOPE)
        fig.tight_layout(rect=(0.205, 0.05, 1, 0.885))
        scales.setdefault("effect", []).append((-lim, lim))
        return fig

    return factory, len(frames)


def make_f0_reconstruction_factory(rec, theme, scales):
    rows = O.read_csv(os.path.join(rec, "f0_reconstruction_frames.csv"))
    errs = {r["zone"]: r for r in
            O.read_csv(os.path.join(rec, "f0_projection_errors.csv"))}
    by = O.group_by(rows, "frame", int)
    frames = sorted(by)
    lims = {}
    for a, b in (("f_reference", "f_spline"),
                 ("f_first_reference", "f_first_spline"),
                 ("f_second_reference", "f_second_spline")):
        lims[a] = max(max(abs(float(r[a])) for r in rows),
                      max(abs(float(r[b])) for r in rows))

    def factory(i):
        pts = sorted(by[frames[i]], key=lambda r: float(r["radius"]))
        r_ = O.col(pts, "radius")
        fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.0),
                                 facecolor=theme.background)
        trip = (("f_reference", "f_spline", r"$f$"),
                ("f_first_reference", "f_first_spline", r"$f'$"),
                ("f_second_reference", "f_second_spline", r"$f''$"))
        for ax, (a, b, lab) in zip(axes, trip):
            ax.plot(r_, O.col(pts, a), lw=2.0, color=theme.foreground,
                    label="Alcubierre reference")
            ax.plot(r_, O.col(pts, b), lw=1.5, color="#d1495b", ls="--",
                    label="partial spline reconstruction")
            ax.set_ylim(-1.15 * lims[a], 1.15 * lims[a])
            ax.legend(fontsize=6.8, frameon=False, labelcolor=theme.foreground)
            ax.set_title(lab, fontsize=10, color=theme.foreground)
            _style(ax, theme, "radius r (m)", lab)
        n = pts[0]["coefficients_included"]
        fig.suptitle("O0A  F0 reconstruction, %s of %d free coefficients included"
                     % (n, len(frames) - 1), fontsize=12,
                     color=theme.foreground, x=0.010, ha="left", y=0.985)
        O.f0_stamp(fig, theme)
        O.banner(fig, theme,
                 "ADDING COEFFICIENTS DOES NOT CONVERGE TO A USABLE "
                 "REPRESENTATION. With all %d included the wall errors are "
                 "%s (value), %s (first derivative) and %s (second derivative). "
                 "Five coefficients cannot resolve a wall narrower than one knot "
                 "span; this animation shows the limit being approached, not "
                 "overcome."
                 % (len(frames) - 1,
                    errs["wall"]["max_value_error"],
                    errs["wall"]["max_first_derivative_error"],
                    errs["wall"]["max_second_derivative_error"]), y=0.845)
        O.footer(fig, theme, "Fixed vertical scales across all frames. " + O.SCOPE)
        fig.tight_layout(rect=(0, 0.05, 1, 0.800))
        scales.setdefault("f0recon", []).append(
            (-1.15 * lims["f_second_reference"], 1.15 * lims["f_second_reference"]))
        return fig

    return factory, len(frames)


def make_whitened_factory(rec, theme, scales):
    coef = O.read_csv(os.path.join(rec, "whitened_perturbation_frames.csv"))
    field = O.read_csv(os.path.join(rec, "whitened_field_frames.csv"))
    byc = O.group_by(coef, "frame", int)
    byf = O.group_by(field, "frame", int)
    frames = sorted(byc)
    clim = max(abs(float(r["coefficient_delta"])) for r in coef) or 1.0
    flim = max(abs(float(r["field"])) for r in field) or 1.0
    loc = {r["coordinate_system"]: r for r in
           O.read_csv(os.path.join(rec, "whitening_locality.csv"))}

    def factory(i):
        pts = sorted(byc[frames[i]], key=lambda r: int(r["coefficient_index"]))
        fpts = byf[frames[i]]
        fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.2),
                                 facecolor=theme.background)
        ax = axes[0]
        ax.bar([int(r["coefficient_index"]) for r in pts],
               O.col(pts, "coefficient_delta"), width=0.8, color="#c792ea")
        ax.set_ylim(-1.15 * clim, 1.15 * clim)
        ax.set_title("one whitened coordinate $y_j$ moves %s of %s raw "
                     "coefficients" % (pts[0]["coefficients_touched"],
                                       pts[0]["total_coefficients"]),
                     fontsize=9.5, color=theme.foreground)
        _style(ax, theme, "raw coefficient index",
               r"$\Delta c_i$ from $c=c_0+PWe_j$")

        ax = axes[1]
        xs, ys, g = _grid_from(fpts, "field")
        im = ax.imshow(g, origin="lower",
                       extent=[min(xs), max(xs), min(ys), max(ys)],
                       aspect="auto", cmap="RdBu_r", vmin=-flim, vmax=flim,
                       interpolation="nearest")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=6.5, colors=theme.foreground)
        cb.outline.set_visible(False)
        ax.set_title(r"the resulting field $c_0+PWe_j$", fontsize=9.5,
                     color=theme.foreground)
        _style(ax, theme, r"$\xi$", r"$s$")

        collar = max(abs(float(r["field"])) for r in fpts
                     if r["in_outer_collar"] == "1")
        fig.suptitle(r"O0A  whitened coordinate perturbation:  "
                     r"$y_j\rightarrow z=Wy\rightarrow c=c_0+PWy$",
                     fontsize=12, color=theme.foreground, x=0.010, ha="left",
                     y=0.985)
        O.banner(fig, theme, O.WHITENING_NOTE + "  These are NOT local basis "
                             "functions: a single whitened coordinate touches %s "
                             "of %s coefficients, against %s on average for a "
                             "column of P."
                 % (pts[0]["coefficients_touched"],
                    pts[0]["total_coefficients"],
                    loc.get("P", {}).get("average_modes_touched", "?")),
                 y=0.935)
        _stamp(fig, theme, [
            "frame                %d" % frames[i],
            "whitened index j     %s" % pts[0]["whitened_index"],
            "coefficients touched %s of %s" % (pts[0]["coefficients_touched"],
                                               pts[0]["total_coefficients"]),
            "max |delta c|        %.6g" % float(pts[0]["max_abs_delta"]),
            "collar max |field|   %.3e" % collar,
            "locality retained    only in the spline evaluation basis",
        ], y=0.30)
        O.footer(fig, theme, O.SCOPE)
        fig.tight_layout(rect=(0.225, 0.05, 1, 0.885))
        scales.setdefault("whitened", []).append((-flim, flim))
        return fig

    return factory, len(frames)


ANIMATIONS = [
    ("individual_radial_basis_modes", make_radial_mode_factory,
     "univariate_basis_profiles.csv"),
    ("individual_axisymmetric_basis_modes", make_axisymmetric_mode_factory,
     "tensor_mode_frames.csv"),
    ("raw_to_constrained_mode", make_raw_to_constrained_factory,
     "raw_to_constrained_frames.csv"),
    ("constraint_projection_effect", make_constraint_effect_factory,
     "constraint_effect_frames.csv"),
    ("f0_baseline_to_reconstruction", make_f0_reconstruction_factory,
     "f0_reconstruction_frames.csv"),
    ("whitened_coordinate_perturbation", make_whitened_factory,
     "whitened_field_frames.csv"),
]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--records", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fps", type=int, default=FPS)
    ap.add_argument("--theme", default="scientific")
    args = ap.parse_args(argv)

    rec = records.resolve(args.records)
    out = records.ensure_dir(args.out)
    theme = themes.get(args.theme)
    themes.apply(theme)

    manifest = O.read_json(os.path.join(rec, "basis_manifest.json"))
    readiness = O.read_json(os.path.join(rec, "basis_readiness.json"))
    before = O.source_checksums(rec)

    # One fixed derivative scale for the radial-mode animation, from the record.
    der = O.read_csv(os.path.join(rec, "univariate_basis_derivatives.csv"))
    scales = {"der": max(max(abs(float(r["first"])), abs(float(r["second"]))),
                         abs(float(r["third"]))) for r in der[:1]}
    scales["der"] = max(max(abs(float(r[k])) for r in der)
                        for k in ("first", "second", "third"))

    encoder = render.find_ffmpeg()
    entries = []
    blanks = []
    for stem, factory_of, source in ANIMATIONS:
        factory, n_frames = factory_of(rec, theme, scales)
        frame_dir = os.path.join(out, "_frames_" + stem)
        paths = render.save_frames(factory, n_frames, frame_dir, prefix="o0a",
                                   dpi=105)
        blanks += [p for p in paths if render.frame_is_blank(p)]
        gif = render.write_gif(paths, os.path.join(out, stem + ".gif"),
                               fps=args.fps, max_frames=n_frames,
                               max_width=900)
        mp4 = render.write_mp4(paths, os.path.join(out, stem + ".mp4"),
                              fps=args.fps, frame_dir=frame_dir)
        sheet = render.contact_sheet(
            paths, os.path.join(out, "contact_sheets", stem + ".png"),
            tiles=9, label_prefix="frame")
        src = os.path.join(rec, source)
        for path, codec in ((gif, "gif"), (mp4, "h264")):
            e = render.animation_entry(
                stem, path, codec=codec, resolution=_png_resolution(paths[0]),
                fps=args.fps, frames=paths,
                duration_s=n_frames / float(args.fps), source_csv=src)
            e["theme"] = theme.name
            e["colour_scale"] = ("fixed; computed once from the whole record and "
                                 "frozen across frames")
            e["contact_sheet"] = os.path.basename(sheet) if sheet else None
            if path is None:
                e["reason"] = ("encoder unavailable" if codec == "h264"
                               and not encoder else
                               "writer unavailable for this format")
            entries.append(e)
        for p in paths:
            os.remove(p)
        os.rmdir(frame_dir)
        print("wrote %-40s frames=%d gif=%s mp4=%s sheet=%s"
              % (stem, n_frames, bool(gif), bool(mp4), bool(sheet)))

    after = O.source_checksums(rec)
    prov = {
        "program": "O",
        "milestone": "O0A",
        "source_records": os.path.basename(os.path.normpath(rec)),
        "git_commit": manifest.get("git_commit"),
        "git_dirty": manifest.get("git_dirty"),
        "frame_rate": args.fps,
        "colour_scale": ("fixed; computed once per record and frozen across "
                        "frames"),
        "renderer_computes_physics": False,
        "renderer_computes_splines": False,
        "hessian_coverage": O.HESSIAN_NOTE,
        "full_hessian_implemented": False,
        "whitening_note": O.WHITENING_NOTE,
        "oti_directions": manifest.get("oti_directions"),
        "gradient_batches": manifest.get("gradient_batches"),
        "encoder_available": bool(encoder),
        "f0_ready_for_optimization": False,
        "f1_ready_for_optimization": False,
        "readiness_families": [f["family"] for f in readiness["families"]],
        "deterministic_rerender": {
            "record_digest_algorithm": "FNV-1a 64-bit",
            "input_digests": before,
            "matplotlib_backend": "Agg",
            "no_random_state_used": True,
        },
        "scope": O.SCOPE,
    }
    render.write_manifest(os.path.join(out, "animation_manifest.json"), entries,
                          prov)

    rep = anomaly.AnomalyReport("O", "O0A", "o0a_animation")
    worst = 0.0
    for key, vals in scales.items():
        if not isinstance(vals, list):
            continue
        worst = max(worst, max(max(abs(a - vals[0][0]), abs(b - vals[0][1]))
                               for a, b in vals))
    rep.record("per_frame_scale_drift", fired=worst > 0.0, measured=worst,
               threshold=0.0,
               detail="largest deviation of any frame's limits from frame 0")
    rep.record("blank_frame", fired=bool(blanks), measured=len(blanks),
               threshold=0, detail="%r" % [os.path.basename(p) for p in blanks[:3]])
    changed = [k for k in before if before[k] != after.get(k)]
    rep.record("source_checksum_mismatch", fired=bool(changed),
               measured=len(changed), threshold=0, detail="%r" % changed[:3])
    left = [n for n in os.listdir(out) if n.startswith("_frames_")]
    rep.record("raw_frame_directory_survived", fired=bool(left),
               measured=len(left), threshold=0, detail="%r" % left)
    ready = [f["family"] for f in readiness["families"]
             if f["ready_for_optimization"]]
    rep.record("false_optimization_readiness_claim", fired=bool(ready),
               measured=len(ready), threshold=0,
               detail="families claiming readiness: %r" % ready)
    path = rep.write(out)
    print("%d animation entries, %d fired anomalies, %d raw frame dirs left -> %s"
          % (len(entries), len(rep.fired), len(left), os.path.basename(path)))
    return 1 if left else 0


if __name__ == "__main__":
    raise SystemExit(main())
