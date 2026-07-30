#!/usr/bin/env python3
"""The six O0B ADM-metric animations, drawn from saved C++ records only.

    python animate_o0b_metric.py --records <o0b_metric dir> --out <gallery dir>

Every frame's numbers were computed by production and written to a record. This
renderer selects, arranges and colours them; it assembles no metric and applies
no coefficient transform. Colour limits are computed once per record and frozen.
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

import o0b_common as O                                           # noqa: E402
from common import anomaly, records, render                      # noqa: E402
from themes import themes                                        # noqa: E402

FPS = 3
BLUE, AMBER, RED, GREEN = "#4ea1ff", "#ffb454", "#d1495b", "#7ed957"


def _stamp(fig, theme, lines, x=0.010, y=0.30):
    fig.text(x, y, "\n".join(lines), ha="left", va="center", fontsize=6.4,
             family="monospace", color=theme.foreground, linespacing=1.5,
             bbox=dict(facecolor=theme.background, edgecolor=theme.grid,
                       alpha=0.8, boxstyle="round,pad=0.5", linewidth=0.6))


def _png_res(path):
    with open(path, "rb") as fh:
        head = fh.read(24)
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        return [0, 0]
    return [int.from_bytes(head[16:20], "big"),
            int.from_bytes(head[20:24], "big")]


def _field(ax, fig, xs, ys, g, theme, lab, lim):
    im = ax.imshow(g, origin="lower",
                   extent=[min(xs), max(xs), min(ys), max(ys)], aspect="auto",
                   cmap="RdBu_r", vmin=-lim, vmax=lim, interpolation="nearest")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=6.5, colors=theme.foreground)
    cb.outline.set_visible(False)
    ax.set_title(lab, fontsize=9.5, color=theme.foreground)
    O.style_axes(ax, theme, "x (m)", "transverse radius (m)")


def make_single_mode(rec, theme, scales):
    rows = O.read_csv(os.path.join(rec, "frames_single_shift_mode.csv"))
    by = O.group_by(rows, "mode", int)
    modes = sorted(by)
    lim = max(abs(float(r["beta_delta"])) for r in rows) or 1.0

    def factory(i):
        sub = by[modes[i]]
        xs, ys, g = O.grid_from(sub, "beta_delta")
        fig, ax = plt.subplots(figsize=(11.6, 5.4), facecolor=theme.background)
        _field(ax, fig, xs, ys, g, theme,
               r"$\Delta\beta^x$ from a single active mode", lim)
        scales.setdefault("single", []).append((-lim, lim))
        fig.suptitle("O0B  one local correction mode at a time, mode %d of %d"
                     % (modes[i], len(modes) - 1), fontsize=12,
                     color=theme.foreground, x=0.010, ha="left", y=0.985)
        O.banner(fig, theme,
                 "The correction has LOCAL SUPPORT, so each mode changes the "
                 "shift only where its own basis function lives. Outside that "
                 "patch the metric is exactly the Alcubierre baseline.")
        _stamp(fig, theme, [
            "mode            %d" % modes[i],
            "amplitude       0.08",
            "max |d beta|    %.6g" % max(abs(float(r["beta_delta"]))
                                         for r in sub),
            "colour scale    fixed across all frames",
        ])
        O.footer(fig, theme, O.SCOPE)
        fig.tight_layout(rect=(0.215, 0.055, 1, 0.885))
        return fig

    return factory, len(modes)


def make_amplitude(rec, theme, scales, key, stem_lab):
    rows = O.read_csv(os.path.join(rec, "frames_amplitude_sweep.csv"))
    by = O.group_by(rows, "frame", int)
    frames = sorted(by)
    lim = max(abs(float(r[key])) for r in rows) or 1.0

    def factory(i):
        sub = by[frames[i]]
        xs, ys, g = O.grid_from(sub, key)
        fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.2),
                                 facecolor=theme.background)
        _field(axes[0], fig, xs, ys, g, theme, stem_lab, lim)
        ax = axes[1]
        amps = [float(by[f][0]["amplitude"]) for f in frames]
        worst = [max(abs(float(r["det_tx"]) + 1.0) for r in by[f])
                 for f in frames]
        ax.semilogy(amps, [max(v, 1e-20) for v in worst], lw=1.0,
                    color=theme.grid)
        ax.semilogy(amps[:i + 1], [max(v, 1e-20) for v in worst[:i + 1]],
                    lw=1.8, color=RED)
        ax.plot([amps[i]], [max(worst[i], 1e-20)], marker="o", markersize=6,
                color=RED)
        ax.set_ylim(1e-20, 1e-10)
        ax.set_title(r"$|\det g_{(t,x)}+1|$ stays exactly zero", fontsize=9.5,
                     color=theme.foreground)
        O.style_axes(ax, theme, "correction amplitude", "residual")
        scales.setdefault(key, []).append((-lim, lim))
        fig.suptitle("O0B  baseline to corrected shift, amplitude %.2f"
                     % float(sub[0]["amplitude"]), fontsize=12,
                     color=theme.foreground, x=0.010, ha="left", y=0.985)
        O.banner(fig, theme,
                 "The correction is added to the ANALYTIC Alcubierre baseline, "
                 "not to the F0 projection: F0 is a representation test that "
                 "O0A found is not resolution-certified. At amplitude 0 the "
                 "metric is exactly the production one.")
        _stamp(fig, theme, [
            "amplitude       %.4f" % float(sub[0]["amplitude"]),
            "max |g01|       %.6g" % max(abs(float(r["g01"])) for r in sub),
            "det residual    %.2e" % worst[i],
            "inverse resid   %.2e" % max(float(r["inverse_residual"])
                                         for r in sub),
        ])
        O.footer(fig, theme, O.CROSSING_NOTE + " " + O.SCOPE)
        fig.tight_layout(rect=(0.205, 0.055, 1, 0.885))
        return fig

    return factory, len(frames)


def make_time_translation(rec, theme, scales):
    rows = O.read_csv(os.path.join(rec, "frames_time_translation.csv"))
    by = O.group_by(rows, "frame", int)
    frames = sorted(by)
    lim = max(abs(float(r["beta_x"])) for r in rows) or 1.0
    dlim = max(abs(float(r["d_t_g01"])) for r in rows) or 1.0

    def factory(i):
        sub = by[frames[i]]
        xs, ys, g = O.grid_from(sub, "beta_x")
        _, _, dg = O.grid_from(sub, "d_t_g01")
        fig, axes = plt.subplots(1, 2, figsize=(14.4, 5.2),
                                 facecolor=theme.background)
        _field(axes[0], fig, xs, ys, g, theme, r"$\beta^x$", lim)
        _field(axes[1], fig, xs, ys, dg, theme,
               r"$\partial_t g_{01}$ (analytic)", dlim)
        xc = float(sub[0]["bubble_center_x"])
        for ax in axes:
            ax.axvline(xc, color=AMBER, lw=1.4, ls="--")
        scales.setdefault("time", []).append((-lim, lim))
        fig.suptitle("O0B  analytic time translation, t = %+.3f"
                     % float(sub[0]["t"]), fontsize=12,
                     color=theme.foreground, x=0.010, ha="left", y=0.985)
        O.banner(fig, theme,
                 "The bubble centre moves as x_c = x_c0 + v_s t, and the time "
                 "derivative is ANALYTIC: d/dt = -(v_s/R0) d/dxi. The right "
                 "panel is nonzero, which is precisely what a frozen-time "
                 "treatment would have set identically to zero.")
        _stamp(fig, theme, [
            "t               %+.4f" % float(sub[0]["t"]),
            "bubble centre   %+.4f m" % xc,
            "max |d_t g01|   %.6g" % max(abs(float(r["d_t_g01"]))
                                         for r in sub),
            "max |d_tt g01|  %.6g" % max(abs(float(r["d_tt_g01"]))
                                         for r in sub),
            "frozen time     NOT used",
        ])
        O.footer(fig, theme, O.SCOPE)
        fig.tight_layout(rect=(0.200, 0.055, 1, 0.885))
        return fig

    return factory, len(frames)


def make_coordinate_systems(rec, theme, scales):
    rows = O.read_csv(os.path.join(rec, "frames_coordinate_systems.csv"))
    by = O.group_by(rows, "frame", int)
    frames = sorted(by)
    lim = max(abs(float(r["beta_delta"])) for r in rows) or 1.0

    def factory(i):
        sub = by[frames[i]]
        fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.0),
                                 facecolor=theme.background)
        for ax, sysname, lab in (
                (axes[0], "raw_c", r"raw $c_k$"),
                (axes[1], "constrained_z", r"constrained $z_j$"),
                (axes[2], "whitened_y", r"whitened $y_j$")):
            part = [r for r in sub if r["system"] == sysname]
            if not part:
                ax.axis("off")
                continue
            xs, ys, g = O.grid_from(part, "beta_delta")
            _field(ax, fig, xs, ys, g, theme,
                   "%s, index %d\ntouches %s of %s coefficients"
                   % (lab, frames[i], part[0]["coefficients_touched"],
                      part[0]["total_coefficients"]), lim)
        scales.setdefault("coords", []).append((-lim, lim))
        fig.suptitle("O0B  the same index in three coordinate systems: "
                     r"$c=c_0+Pz$, $z=Wy$", fontsize=12,
                     color=theme.foreground, x=0.010, ha="left", y=0.985)
        O.banner(fig, theme,
                 "These are DIFFERENT FIELDS. A raw index moves one local "
                 "coefficient; a whitened index moves nearly all of them. " +
                 O.WHITENED_NOTE)
        O.footer(fig, theme,
                 "Whitening is a dense transform: locality survives only in the "
                 "underlying spline evaluation basis, never in the optimization "
                 "coordinates. " + O.SCOPE)
        fig.tight_layout(rect=(0, 0.055, 1, 0.875))
        return fig

    return factory, len(frames)


def make_sampling(rec, theme, scales):
    rows = O.read_csv(os.path.join(rec, "frames_sampling_comparison.csv"))
    by = O.group_by(rows, "frame", int)
    frames = sorted(by)
    cost = O.read_json(os.path.join(rec, "uniform_grid_cost_estimate.json"))
    q = O.read_json(os.path.join(rec, "elementwise_quadrature_manifest.json"))

    def factory(i):
        sub = by[frames[i]]
        el = [r for r in sub if r["route"] == "B"]
        fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.2),
                                 facecolor=theme.background)
        ax = axes[0]
        if el:
            ax.plot(O.col(el, "xi"), O.col(el, "s"), marker="s", markersize=7,
                    linestyle="none", color=GREEN, alpha=0.8)
        ax.set_xlim(-3.2, 3.2); ax.set_ylim(-0.2, 4.2)
        ax.set_title("Route B: %d of %d elements revealed"
                     % (len(el), q["elements"]), fontsize=9.5,
                     color=theme.foreground)
        O.style_axes(ax, theme, r"$\xi$", r"$s$")

        ax = axes[1]
        ax.bar([0, 1], [cost["total_cells"], q["total_nodes"]],
               color=[RED, GREEN], width=0.5)
        ax.set_yscale("log")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Route A\n%.3g cells" % cost["total_cells"],
                            "Route B\n%d nodes" % q["total_nodes"]],
                           fontsize=8.5, color=theme.foreground)
        ax.set_title("evaluation points required", fontsize=9.5,
                     color=theme.foreground)
        O.style_axes(ax, theme, "", "points")
        scales.setdefault("sampling", []).append((0.0, 1.0))

        fig.suptitle("O0B  refined uniform grid against elementwise quadrature",
                     fontsize=12, color=theme.foreground, x=0.010, ha="left",
                     y=0.985)
        O.banner(fig, theme,
                 "Route A must resolve the narrowest feature with six CELLS "
                 "because it takes finite differences. Route B takes none -- "
                 "the derivatives are analytic -- so the six-cell rule does not "
                 "transfer to it, and what must converge instead is the "
                 "integral. NEITHER ROUTE IS RUN IN O0B.")
        O.footer(fig, theme, q["convergence_contract"] + " " + O.SCOPE)
        fig.tight_layout(rect=(0, 0.055, 1, 0.875))
        return fig

    return factory, len(frames)


ANIMATIONS = [
    ("single_shift_mode", lambda r, t, s: make_single_mode(r, t, s),
     "frames_single_shift_mode.csv"),
    ("baseline_to_corrected_shift",
     lambda r, t, s: make_amplitude(r, t, s, "g01", r"$g_{01}=\beta^x$"),
     "frames_amplitude_sweep.csv"),
    ("metric_components_vs_mode_amplitude",
     lambda r, t, s: make_amplitude(r, t, s, "g00", r"$g_{00}$"),
     "frames_amplitude_sweep.csv"),
    ("analytic_time_translation", lambda r, t, s: make_time_translation(r, t, s),
     "frames_time_translation.csv"),
    ("raw_constrained_whitened_metric_effect",
     lambda r, t, s: make_coordinate_systems(r, t, s),
     "frames_coordinate_systems.csv"),
    ("uniform_grid_vs_elementwise_sampling",
     lambda r, t, s: make_sampling(r, t, s),
     "frames_sampling_comparison.csv"),
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

    man = O.read_json(os.path.join(rec, "metric_design_manifest.json"))
    rd = O.read_json(os.path.join(rec, "readiness.json"))
    before = O.source_checksums(rec)
    encoder = render.find_ffmpeg()
    entries, blanks, scales = [], [], {}

    for stem, factory_of, source in ANIMATIONS:
        factory, n = factory_of(rec, theme, scales)
        frame_dir = os.path.join(out, "_frames_" + stem)
        paths = render.save_frames(factory, n, frame_dir, prefix="o0b", dpi=105)
        blanks += [p for p in paths if render.frame_is_blank(p)]
        gif = render.write_gif(paths, os.path.join(out, stem + ".gif"),
                               fps=args.fps, max_frames=n, max_width=900)
        mp4 = render.write_mp4(paths, os.path.join(out, stem + ".mp4"),
                               fps=args.fps, frame_dir=frame_dir)
        sheet = render.contact_sheet(
            paths, os.path.join(out, "contact_sheets", stem + ".png"), tiles=9,
            label_prefix="frame")
        for path, codec in ((gif, "gif"), (mp4, "h264")):
            e = render.animation_entry(
                stem, path, codec=codec, resolution=_png_res(paths[0]),
                fps=args.fps, frames=paths, duration_s=n / float(args.fps),
                source_csv=os.path.join(rec, source))
            e["theme"] = theme.name
            e["colour_scale"] = ("fixed; computed once from the whole record "
                                 "and frozen across frames")
            e["contact_sheet"] = os.path.basename(sheet) if sheet else None
            if path is None:
                e["reason"] = ("encoder unavailable" if codec == "h264"
                               and not encoder else
                               "writer unavailable for this format")
            entries.append(e)
        for p in paths:
            os.remove(p)
        os.rmdir(frame_dir)
        print("wrote %-42s frames=%d gif=%s mp4=%s sheet=%s"
              % (stem, n, bool(gif), bool(mp4), bool(sheet)))

    after = O.source_checksums(rec)
    render.write_manifest(
        os.path.join(out, "animation_manifest.json"), entries,
        {"program": "O", "milestone": "O0B",
         "source_records": os.path.basename(os.path.normpath(rec)),
         "git_commit": man.get("git_commit"), "frame_rate": args.fps,
         "colour_scale": "fixed; computed once per record and frozen",
         "renderer_computes_physics": False,
         "renderer_assembles_metrics": False,
         "index_ordering": O.INDEX_ORDER,
         "hessian_coverage": O.HESSIAN_NOTE,
         "full_hessian_implemented": False,
         "energy_evaluated": False, "optimizer_run": False,
         "required_statement": rd["required_statement"],
         "f1_resolution_certified": False,
         "ready_for_energy_optimization": False,
         "whitened_index_identity": O.WHITENED_NOTE,
         "encoder_available": bool(encoder),
         "deterministic_rerender": {"input_digests": before,
                                    "matplotlib_backend": "Agg",
                                    "no_random_state_used": True},
         "scope": O.SCOPE})

    rep = anomaly.AnomalyReport("O", "O0B", "o0b_animation")
    worst = 0.0
    for vals in scales.values():
        worst = max(worst, max(max(abs(a - vals[0][0]), abs(b - vals[0][1]))
                               for a, b in vals))
    rep.record("per_frame_scale_drift", fired=worst > 0.0, measured=worst,
               threshold=0.0, detail="deviation of any frame's limits from "
                                     "frame 0")
    rep.record("blank_frame", fired=bool(blanks), measured=len(blanks),
               threshold=0, detail=repr([os.path.basename(p)
                                         for p in blanks[:3]]))
    changed = [k for k in before if before[k] != after.get(k)]
    rep.record("source_checksum_mismatch", fired=bool(changed),
               measured=len(changed), threshold=0, detail=repr(changed[:3]))
    left = [n for n in os.listdir(out) if n.startswith("_frames_")]
    rep.record("raw_frame_directory_survived", fired=bool(left),
               measured=len(left), threshold=0, detail=repr(left))
    ready = [f["family"] for f in rd["families"]
             if f["ready_for_energy_optimization"]]
    rep.record("false_optimization_readiness_claim", fired=bool(ready),
               measured=len(ready), threshold=0, detail=repr(ready))
    path = rep.write(out)
    print("%d entries, %d fired anomalies, %d raw frame dirs -> %s"
          % (len(entries), len(rep.fired), len(left), os.path.basename(path)))
    return 1 if left else 0


if __name__ == "__main__":
    raise SystemExit(main())
