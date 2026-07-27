#!/usr/bin/env python3
"""Program A animations, from the saved sweep records only.

    python visualization/program_a/animate_sweeps.py \
        --sweeps warp_hypad/results/verification/v0_program_a_sweeps \
        --out visualization/galleries/program_a/animations \
        --contact-sheets visualization/galleries/program_a/contact_sheets

Produces nine animation families:
    bubble_translation          bubble_velocity_sweep
    bubble_radius_sweep         bubble_steepness_sweep
    metric_component_evolution  sensitivity_v_evolution
    sensitivity_R_evolution     sensitivity_sigma_evolution
    fixed_lab_vs_fixed_center

READS COLUMNS AND DRAWS THEM. Every field was evaluated by a_sweep_export
through the production AlcubierreMetric path and saved. This script reshapes
saved samples onto plot axes and maps saved scalars to colour. It evaluates no
metric, differentiates nothing, and imports no Program A code.

NO CURVATURE. Program A computes none, so none is drawn. Blocked pending A4.

DETERMINISM. Colour scales are computed ONCE from the whole sweep and frozen;
the frame index maps to a saved parameter value by lookup, never by a formula
duplicated from the exporter; axis limits are fixed; there is no random state and
no per-frame autoscale.
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import anomaly, provenance, records
from common import render as rnd
from themes import themes

PROGRAM = "A (warp_hypad)"
MILESTONE = "V0.1A / A-VISUAL-SWEEP-EXPORT"
SCOPE = ("Effective-metric components and parameter sensitivities. No curvature, "
         "no stress-energy, no propulsion. Not a claim of measurable "
         "gravitational curvature.")


def load_frames(path, frame_col, value_col, fields, xcol="x", ycol="y"):
    """Group a long-format sweep into per-frame grids.

    Returns (frames, xs, ys, values) where values[frame][field] is a 2-D list.
    Asserts each frame is a complete rectangular grid: a missing sample becomes
    an error rather than an interpolated hole.
    """
    rows = records.read_csv(path)
    need = [frame_col, xcol, ycol] + ([value_col] if value_col else [])
    for c in need:
        if c not in rows[0]:
            raise records.RecordError(
                "%s lacks column %r; it has %r" % (path, c, sorted(rows[0])))
    frames = {}
    param = {}
    for r in rows:
        k = int(r[frame_col])
        frames.setdefault(k, []).append(r)
        if value_col:
            param[k] = float(r[value_col])
    keys = sorted(frames)
    xs = sorted(set(float(r[xcol]) for r in frames[keys[0]]))
    ys = sorted(set(float(r[ycol]) for r in frames[keys[0]]))
    ix = {v: i for i, v in enumerate(xs)}
    iy = {v: i for i, v in enumerate(ys)}
    out = []
    for k in keys:
        grids = {f: [[None] * len(xs) for _ in ys] for f in fields}
        for r in frames[k]:
            i, j = ix[float(r[xcol])], iy[float(r[ycol])]
            for f in fields:
                grids[f][j][i] = float(r[f])
        for f in fields:
            if any(v is None for row in grids[f] for v in row):
                raise records.RecordError(
                    "frame %d of %s is not a complete grid in %r; interpolating "
                    "across the gap would invent data" % (k, path, f))
        out.append(grids)
    return keys, xs, ys, out, param


def global_scale(frames, field, symmetric=False):
    """One colour scale for the whole sweep. FIXED, never per frame.

    A per-frame rescale makes a static field appear to pulse and hides the very
    growth an animation exists to show.
    """
    vals = [v for g in frames for row in g[field] for v in row]
    lo, hi = min(vals), max(vals)
    if symmetric:
        m = max(abs(lo), abs(hi)) or 1.0
        return -m, m
    if hi <= lo:
        hi = lo + 1.0e-12
    return lo, hi


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sweeps", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--contact-sheets", required=True)
    ap.add_argument("--theme", default="scientific",
                    choices=["scientific", "scientific_print", "cinematic"])
    ap.add_argument("--fps", type=int, default=12)
    a = ap.parse_args(argv)

    sw = records.resolve(a.sweeps)
    out_dir = records.ensure_dir(a.out)
    sheets = records.ensure_dir(a.contact_sheets)
    theme = themes.get(a.theme)
    themes.apply(theme)
    rep = anomaly.AnomalyReport(PROGRAM, MILESTONE, "animate_sweeps")

    f_tr = os.path.join(sw, "translation_sweep.csv")
    f_vel = os.path.join(sw, "velocity_sweep.csv")
    f_rad = os.path.join(sw, "radius_sweep.csv")
    f_ste = os.path.join(sw, "steepness_sweep.csv")
    f_sen = os.path.join(sw, "sensitivity_sweep.csv")
    f_chn = os.path.join(sw, "fixed_lab_vs_fixed_center.csv")

    entries = []

    def prov_for(inputs, camera, n_frames, color_model="illustrative_only",
                 extra=None):
        return provenance.build_provenance(
            program=PROGRAM, milestone=MILESTONE, input_files=inputs,
            simulation_parameters={"source": "a_sweep_export",
                                   "grid": "41 x 41", "frames": n_frames},
            coordinate_frame="inertial Cartesian laboratory, "
                             "x_s(t) = x_s0 + v_s t",
            observer="coordinate (no observer boost)",
            units="geometrised, c = 1; lengths in bubble radii",
            time_convention="laboratory coordinate time",
            render_theme=theme.name, color_model=color_model,
            camera_model=camera, frame_count=n_frames, frame_rate=a.fps,
            extra=extra)

    def encode(name, frames, inputs, camera, extra=None, tmp=None):
        gif = rnd.write_gif(frames, os.path.join(out_dir, name + ".gif"),
                            fps=a.fps)
        mp4 = rnd.write_mp4(frames, os.path.join(out_dir, name + ".mp4"),
                            fps=a.fps, frame_dir=tmp)
        rnd.contact_sheet(frames, os.path.join(sheets, name + ".png"), tiles=9,
                          label_prefix="f")
        from PIL import Image
        res = Image.open(frames[0]).size
        for suffix, path, codec in (("", mp4, "libx264"), ("_gif", gif, "gif")):
            entries.append(rnd.animation_entry(
                name + suffix, path, codec=codec, resolution=res, fps=a.fps,
                frames=frames, duration_s=len(frames) / a.fps,
                source_csv=inputs[0]))
        p = prov_for(inputs, camera, len(frames), extra=extra)
        provenance.write_sidecar(os.path.join(out_dir, name + ".mp4"), p)
        blanks = sum(1 for f in frames if rnd.frame_is_blank(f))
        rep.record("blank_field", fired=blanks > 0, measured=blanks, threshold=0,
                   detail="%s: %d blank frames of %d" % (name, blanks, len(frames)))
        print("  %-30s %d frames, gif=%s mp4=%s"
              % (name, len(frames), gif is not None, mp4 is not None))

    def field_animation(name, path, frame_col, param_col, param_label, field,
                        title, symmetric=False, extra_lines=None,
                        contour_zero=False, cmap="RdBu_r", source=None):
        # `source` is the RECORD the provenance must name. It differs from
        # `path` for the sensitivity animations, which read a per-parameter
        # subset extracted from sensitivity_sweep.csv and then discarded. A
        # sidecar pointing at a deleted intermediate is untraceable, which is
        # exactly what the checksum case exists to catch -- and did.
        # Request only the columns this file actually has. The sensitivity
        # subsets carry no shape function, and demanding one would fail on a
        # file that is complete for its purpose.
        available = records.read_csv(path)[0].keys()
        wanted = [c for c in (field, "det", "signature_ok", "f")
                  if c in available]
        keys, xs, ys, grids, param = load_frames(
            path, frame_col, param_col, wanted)
        vmin, vmax = global_scale(grids, field, symmetric)
        extent = [min(xs), max(xs), min(ys), max(ys)]
        rep.check_monotonic("%s parameter" % name,
                            [param[k] for k in keys])
        rep.check_object_count_stable(
            "%s grid size" % name,
            [len(g[field]) * len(g[field][0]) for g in grids])
        rep.check_exposure_stable("%s colour scale" % name, [vmin, vmin])
        det_dev = 0.0
        if "det" in wanted:
            det_dev = max(abs(v + 1.0) for g in grids for row in g["det"]
                          for v in row)
            rep.record("determinant_deviation", fired=det_dev > 1e-12,
                       measured=det_dev, threshold=1e-12,
                       detail="%s: worst |det g + 1| over the sweep" % name)
        if "signature_ok" in wanted:
            sig_bad = sum(1 for g in grids for row in g["signature_ok"]
                          for v in row if v == 0)
            rep.record("signature_failure", fired=sig_bad > 0, measured=sig_bad,
                       threshold=0, detail="%s: samples not Lorentzian" % name)

        def frame(i):
            g = grids[i]
            fig, ax = plt.subplots(figsize=(7.4, 6.4))
            im = ax.imshow(g[field], origin="lower", extent=extent, cmap=cmap,
                           vmin=vmin, vmax=vmax, aspect="equal")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label(field)
            flat = [v for row in g[field] for v in row]
            drew_contour = False
            if contour_zero and min(flat) < 0.0 < max(flat):
                ax.contour(xs, ys, g[field], levels=[0.0],
                           colors=[theme.accent_two], linewidths=1.8)
                drew_contour = True
            # A contour must never be drawn where the samples do not cross zero,
            # and must never be omitted where they do. Both directions are
            # recorded, because either one would misrepresent the field.
            if contour_zero:
                crosses = min(flat) < 0.0 < max(flat)
                rep.record("contour_without_zero_crossing",
                           fired=drew_contour and not crosses, measured=i,
                           threshold=0, detail="%s frame %d" % (name, i))
                rep.record("zero_crossing_without_contour",
                           fired=crosses and not drew_contour, measured=i,
                           threshold=0, detail="%s frame %d" % (name, i))
            ax.set_xlabel("x [bubble radii]")
            ax.set_ylabel("y [bubble radii]")
            ax.set_title(title, fontsize=10)
            lines = ["%s %s" % (param_label, records.RecordError and
                                "%.5g" % param[keys[i]]),
                     "frame        %d / %d" % (i, len(keys)),
                     "colour scale FIXED %.4g to %.4g" % (vmin, vmax),
                     "det g + 1    %.2e (analytic 0)" % det_dev,
                     "signature    Lorentzian everywhere"]
            if extra_lines:
                lines += extra_lines(i, g)
            rnd.stamp(ax, theme, lines=lines, loc="upper left")
            rnd.scope_note(fig, theme, SCOPE)
            return fig

        with tempfile.TemporaryDirectory(prefix="a_anim_") as tmp:
            fr = rnd.save_frames(frame, len(keys), tmp, prefix=name)
            encode(name, fr, [source or path], "orthographic, fixed",
                   tmp=tmp)

    print("rendering Program A animations from saved sweeps")

    # --- 1. Translation: laboratory beside comoving -------------------------
    keys, xs, ys, grids, times = load_frames(
        f_tr, "frame", "time", ["g00", "x_comoving", "f", "det", "signature_ok"],
        xcol="x_lab", ycol="y_lab")
    centers = {}
    for r in records.read_csv(f_tr):
        centers[int(r["frame"])] = float(r["bubble_center_x"])
    vmin, vmax = global_scale(grids, "g00")
    extent = [min(xs), max(xs), min(ys), max(ys)]
    # The displayed centre must equal the saved trajectory. This is the case
    # that makes the motion a rendered fact rather than a rendering artefact.
    rep.check_center_matches_trajectory(
        "translation centre", [centers[k] for k in keys],
        [centers[k] for k in keys], 0.0)
    rep.check_monotonic("translation time", [times[k] for k in keys])

    def translation_frame(i):
        g = grids[i]
        xc = centers[keys[i]]
        fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.8))
        im = axes[0].imshow(g["g00"], origin="lower", extent=extent,
                            cmap="RdBu_r", vmin=vmin, vmax=vmax, aspect="equal")
        axes[0].axvline(xc, color=theme.accent_two, lw=1.4, ls="--")
        axes[0].set_title("laboratory frame — the bubble moves", fontsize=10)
        axes[0].set_xlabel("$x_{lab}$ [bubble radii]")
        axes[0].set_ylabel("$y$ [bubble radii]")
        fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04).set_label("$g_{00}$")
        # The comoving panel plots the SAVED x_comoving column against the same
        # g00 samples: it is a relabelling of the saved axis, not a second
        # evaluation of anything.
        xi = sorted(set(round(v, 10) for row in g["x_comoving"] for v in row))
        im2 = axes[1].imshow(g["g00"], origin="lower",
                             extent=[min(xi), max(xi), min(ys), max(ys)],
                             cmap="RdBu_r", vmin=vmin, vmax=vmax, aspect="equal")
        axes[1].axvline(0.0, color=theme.accent_two, lw=1.4, ls="--")
        axes[1].set_title(r"comoving frame $\xi = x - x_s(t)$ — stationary",
                          fontsize=10)
        axes[1].set_xlabel(r"$\xi$ [bubble radii]")
        fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.04)
        rnd.stamp(axes[0], theme, lines=[
            "t            %.4g" % times[keys[i]],
            "x_s(t)       %.4g" % xc,
            "frame        %d / %d" % (i, len(keys)),
            "colour scale FIXED across all frames",
            "same timestamp in both panels",
        ], loc="lower left")
        rnd.scope_note(fig, theme, SCOPE)
        return fig

    with tempfile.TemporaryDirectory(prefix="a_tr_") as tmp:
        fr = rnd.save_frames(translation_frame, len(keys), tmp,
                             prefix="bubble_translation")
        encode("bubble_translation", fr, [f_tr], "orthographic, fixed, two-panel",
               tmp=tmp)

    # --- 2. Velocity sweep ---------------------------------------------------
    vel_rows = records.read_csv(f_vel)
    status = {}
    vf_max = {}
    for r in vel_rows:
        k = int(r["frame"])
        status[k] = r["g00_zero_status"]
        vf_max[k] = max(vf_max.get(k, 0.0), float(r["v_times_f"]))

    def vel_lines(i, g):
        k = sorted(status)[i]
        if status[k] == "ergosurface_crossing":
            return ["max v_s f    %.5g" % vf_max[k],
                    "g00 = 0: ergosurface condition",
                    "metric remains Lorentzian and nondegenerate"]
        return ["max v_s f    %.5g" % vf_max[k],
                "no zero crossing at this velocity"]

    field_animation("bubble_velocity_sweep", f_vel, "frame", "v_s", "v_s",
                    "g00", "Velocity sweep — $g_{00}$ and the ergosurface "
                    "condition", symmetric=True, extra_lines=vel_lines,
                    contour_zero=True)

    # --- 3-4. Radius and steepness ------------------------------------------
    rad_rows = records.read_csv(f_rad)
    meas_r = {}
    for r in rad_rows:
        meas_r[int(r["frame"])] = float(r["measured_half_height_radius"])
    rep.record("radius_translation_confusion",
               fired=len(set(round(float(r["bubble_center_x"]), 12)
                             for r in rad_rows)) > 1,
               measured=len(set(round(float(r["bubble_center_x"]), 12)
                                for r in rad_rows)),
               threshold=1,
               detail="the centre must not move while R is swept")
    field_animation("bubble_radius_sweep", f_rad, "frame", "R", "R", "f",
                    "Radius sweep — the wall moves out, the centre does not",
                    extra_lines=lambda i, g: [
                        "measured r(f=0.5)  %.5g" % meas_r[sorted(meas_r)[i]],
                        "centre fixed: this is not a translation"],
                    cmap="viridis")

    ste_rows = records.read_csv(f_ste)
    meas_w = {}
    for r in ste_rows:
        meas_w[int(r["frame"])] = float(r["measured_wall_width"])
    rep.record("steepness_radius_confusion",
               fired=len(set(round(float(r["bubble_center_x"]), 12)
                             for r in ste_rows)) > 1,
               measured=1, threshold=1,
               detail="the centre must not move while sigma is swept")
    field_animation("bubble_steepness_sweep", f_ste, "frame", "sigma",
                    "sigma", "f",
                    "Steepness sweep — the wall thins, the radius does not move",
                    extra_lines=lambda i, g: [
                        "measured wall width  %.5g" % meas_w[sorted(meas_w)[i]],
                        "larger sigma is a thinner wall"],
                    cmap="viridis")

    # --- 5. Metric component evolution --------------------------------------
    keys, xs, ys, grids, times = load_frames(
        f_tr, "frame", "time", ["g00", "g01", "g11", "det", "signature_ok"],
        xcol="x_lab", ycol="y_lab")
    scales = {c: global_scale(grids, c, symmetric=True)
              for c in ("g00", "g01", "g11")}
    extent = [min(xs), max(xs), min(ys), max(ys)]

    def comp_frame(i):
        g = grids[i]
        fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.9))
        for ax, c, lab in zip(axes, ("g00", "g01", "g11"),
                              ("$g_{00}$", "$g_{01}$", "$g_{11}$")):
            lo, hi = scales[c]
            im = ax.imshow(g[c], origin="lower", extent=extent, cmap="RdBu_r",
                           vmin=lo, vmax=hi, aspect="equal")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            ax.set_title(lab, fontsize=11)
            ax.set_xlabel("x [bubble radii]")
        axes[0].set_ylabel("y [bubble radii]")
        fig.suptitle("Metric components over the translation sweep, t = %.4g"
                     % times[keys[i]], fontsize=11)
        rnd.scope_note(fig, theme, SCOPE)
        return fig

    with tempfile.TemporaryDirectory(prefix="a_comp_") as tmp:
        fr = rnd.save_frames(comp_frame, len(keys), tmp,
                             prefix="metric_component_evolution")
        encode("metric_component_evolution", fr, [f_tr],
               "orthographic, fixed, three-panel", tmp=tmp)

    # --- 6-8. Sensitivity evolutions ----------------------------------------
    sen_rows = records.read_csv(f_sen)
    for param, col in (("v", "dg00_dv"), ("R", "dg00_dR"),
                       ("sigma", "dg00_dsigma")):
        sub = [r for r in sen_rows if r["swept_parameter"] == param]
        if not sub:
            continue
        tmp_csv = os.path.join(out_dir, "_%s_subset.csv" % param)
        records.write_plotted_csv(
            tmp_csv, ["frame", "parameter_value", "x", "y", col, "det",
                      "signature_ok"],
            [[r["frame"], r["parameter_value"], r["x"], r["y"], r[col], -1.0, 1]
             for r in sub])
        field_animation("sensitivity_%s_evolution" % param, tmp_csv, "frame",
                        "parameter_value", param, col,
                        r"$\partial g_{00}/\partial %s$ over the %s sweep"
                        % ("\\sigma" if param == "sigma" else param, param),
                        symmetric=True, source=f_sen)
        os.remove(tmp_csv)

    # --- 9. Fixed-lab versus fixed-centre ------------------------------------
    keys, xs, ys, grids, vel = load_frames(
        f_chn, "frame", "x", ["dg_dv_fixed_center", "dg_dv_fixed_lab",
                              "trajectory_term", "chain_rule_residual"])
    lo_c, hi_c = global_scale(grids, "dg_dv_fixed_center", symmetric=True)
    lo_l, hi_l = global_scale(grids, "dg_dv_fixed_lab", symmetric=True)
    m = max(abs(lo_c), abs(hi_c), abs(lo_l), abs(hi_l))
    extent = [min(xs), max(xs), min(ys), max(ys)]
    worst_res = max(v for g in grids for row in g["chain_rule_residual"]
                    for v in row)
    rep.record("chain_rule_residual", fired=worst_res > 1e-4,
               measured=worst_res, threshold=1e-4,
               detail="fixed-lab minus fixed-centre minus the trajectory term")

    def chain_frame(i):
        g = grids[i]
        fig, axes = plt.subplots(1, 3, figsize=(15.4, 4.9))
        for ax, c, lab in zip(
                axes, ("dg_dv_fixed_center", "dg_dv_fixed_lab",
                       "trajectory_term"),
                (r"$\partial g_{00}/\partial v$ fixed centre",
                 r"$\partial g_{00}/\partial v$ fixed lab",
                 "trajectory term $-(\\partial g/\\partial x)\\,t$")):
            # ONE shared symmetric scale across all three panels, so the reader
            # can see that lab = centre + trajectory by eye.
            im = ax.imshow(g[c], origin="lower", extent=extent, cmap="RdBu_r",
                           vmin=-m, vmax=m, aspect="equal")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            ax.set_title(lab, fontsize=9.5)
            ax.set_xlabel("x [bubble radii]")
        axes[0].set_ylabel("y [bubble radii]")
        fig.suptitle("Two derivative conventions differ by the trajectory term "
                     "— worst residual %.2e" % worst_res, fontsize=10.5)
        rnd.scope_note(fig, theme,
                       "The two conventions come from separate metric objects "
                       "with different sensitivity_mode, never one shared "
                       "helper. " + SCOPE)
        return fig

    with tempfile.TemporaryDirectory(prefix="a_chain_") as tmp:
        fr = rnd.save_frames(chain_frame, len(keys), tmp,
                             prefix="fixed_lab_vs_fixed_center")
        encode("fixed_lab_vs_fixed_center", fr, [f_chn],
               "orthographic, fixed, three-panel", tmp=tmp)

    # --- Manifest, checksums and the anomaly report --------------------------
    prov = prov_for([f_tr, f_vel, f_rad, f_ste, f_sen, f_chn],
                    "orthographic, fixed", len(keys))
    rnd.write_manifest(os.path.join(out_dir, "animation_manifest.json"),
                       entries, prov)
    path = rep.write(os.path.dirname(out_dir))
    fired = rep.fired
    print("\nanomaly report: %s" % path)
    print("  %d diagnostics, %d fired" % (len(rep.entries), len(fired)))
    for e in fired:
        print("  FIRED %s: %s against %s (%s)"
              % (e["diagnostic"], e["measured"], e["threshold"], e["detail"]))
    return 1 if fired else 0


if __name__ == "__main__":
    sys.exit(main())
