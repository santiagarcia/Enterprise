"""Shared reading, styling and provenance for the O0B ADM-metric visuals.

    C++ gate -> saved records -> these scripts -> PNG/SVG/GIF/MP4

Nothing here assembles a metric, inverts one, applies a chain rule or forms a
coefficient transform. Every number plotted is a column `o0b_verify` wrote,
including every animation frame: a Python renderer could multiply P and W and
sum basis values to get a metric, and that would be a second, unverified
implementation of the design layer. The visual gate scans this package for
production names to keep that true.

WHAT EVERY READINESS PRODUCT MUST SAY. O0B verified metric CONSTRUCTION. It did
not certify the spline family for energy evaluation, and the underlying O0A F1
family is still short of the six-cell resolution contract. Those are different
claims and the banner keeps them apart.
"""

from __future__ import annotations

import csv
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

from common import provenance, records          # noqa: E402


SCOPE = (
    "O0B maps verified spline coordinates into structurally valid ADM metrics "
    "and analytic metric derivatives. It does not yet minimize energy, certify "
    "the spline family for energy evaluation, demonstrate a causal shortcut or "
    "establish physical realizability.")

HESSIAN_NOTE = (
    "The complete design gradient is available through twelve batches. Mixed "
    "Hessian entries between distinct batches are not available.")

READINESS_BANNER = ("metric construction verified\n"
                    "current F1 family not yet resolution-certified\n"
                    "not ready for energy optimization")

CROSSING_NOTE = (
    "g_00 = 0 is where the coordinate-stationary vector becomes null. It is NOT "
    "a metric degeneracy, NOT a signature change, NOT a curvature singularity "
    "and NOT automatically a horizon: det g stays exactly -1 there.")

WHITENED_NOTE = (
    "The whitened Gram is the identity, so its spectrum is maximally degenerate "
    "and every orthonormal rotation of the whitened axes is equally valid. A "
    "whitened coordinate index has NO unique physical mode identity.")

INDEX_ORDER = "0=t, 1=x, 2=y, 3=z"

REQUIRED_RECORDS = (
    "metric_design_manifest.json", "coefficient_space_manifest.json",
    "readiness.json", "o0a_resolution_consistency.json",
    "coefficient_transforms.csv", "adm_field_samples.csv",
    "metric_samples.csv", "inverse_metric_samples.csv",
    "metric_first_derivatives.csv", "metric_second_derivatives.csv",
    "metric_third_derivatives.csv", "temporal_derivative_checks.csv",
    "derivative_symmetry_checks.csv", "alcubierre_baseline_comparison.csv",
    "single_mode_metric_checks.csv", "metric_validity_checks.csv",
    "determinant_checks.csv", "coefficient_gradient_checks.csv",
    "uniform_grid_cost_estimate.json", "elementwise_quadrature_elements.csv",
    "elementwise_quadrature_manifest.json", "verification_summary.json",
    "frames_single_shift_mode.csv", "frames_amplitude_sweep.csv",
    "frames_time_translation.csv", "frames_coordinate_systems.csv",
    "frames_sampling_comparison.csv",
)


def read_csv(path):
    with open(records.resolve(path), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def read_json(path):
    with open(records.resolve(path), encoding="utf-8") as fh:
        return json.load(fh)


def fnv1a64(path):
    h = 0xCBF29CE484222325
    with open(records.resolve(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            for b in chunk:
                h ^= b
                h = (h * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return "%016x" % h


def source_checksums(rec):
    out = {}
    for n in sorted(os.listdir(records.resolve(rec))):
        p = os.path.join(rec, n)
        if os.path.isfile(records.resolve(p)):
            out[n] = fnv1a64(p)
    return out


def col(rows, key, cast=float):
    return [cast(r[key]) for r in rows]


def group_by(rows, key, cast=str):
    out = {}
    for r in rows:
        out.setdefault(cast(r[key]), []).append(r)
    return out


def grid_from(rows, key, xkey="x", ykey="y"):
    xs = sorted({float(r[xkey]) for r in rows})
    ys = sorted({float(r[ykey]) for r in rows})
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}
    g = [[0.0] * len(xs) for _ in ys]
    for r in rows:
        g[yi[float(r[ykey])]][xi[float(r[xkey])]] = float(r[key])
    return xs, ys, g


def style_axes(ax, theme, xlabel="", ylabel=""):
    ax.set_xlabel(xlabel, fontsize=8, color=theme.foreground)
    ax.set_ylabel(ylabel, fontsize=8, color=theme.foreground)
    ax.tick_params(labelsize=7, colors=theme.foreground)
    ax.set_facecolor(theme.background)
    for s in ax.spines.values():
        s.set_color(theme.grid)


def banner(fig, theme, text, y=0.930, color=None):
    fig.text(0.010, y, text, ha="left", va="top", fontsize=8.0,
             color=color or theme.uncertainty_color)


def footer(fig, theme, text, y=0.004):
    fig.text(0.5, y, text, ha="center", va="bottom", fontsize=6.2,
             color=theme.foreground, alpha=0.75, wrap=True)


def readiness_stamp(fig, theme, x=0.985, y=0.955):
    """The three-line O0B status, boxed, on every readiness product."""
    fig.text(x, y, READINESS_BANNER, ha="right", va="top", fontsize=8.2,
             color=theme.uncertainty_color, family="monospace",
             linespacing=1.5,
             bbox=dict(facecolor=theme.background,
                       edgecolor=theme.uncertainty_color, alpha=0.9,
                       boxstyle="round,pad=0.5", linewidth=1.1))


def simulation_parameters(manifest):
    return {k: manifest.get(k) for k in (
        "family", "baseline", "shift_sign_convention",
        "lapse_parameterisation", "spatial_metric_parameterisation",
        "velocity", "bubble_radius", "wall_sigma", "scale_R0",
        "raw_dimension", "constraint_rank", "constrained_dimension",
        "whitened_dimension", "oti_directions", "gradient_batches",
        "index_ordering", "time_treatment")}


def save_figure(fig, out_dir, stem, *, manifest, inputs, theme, extra=None,
                plotted_header=None, plotted_rows=None, source_rows=None):
    import matplotlib.pyplot as plt

    out = records.ensure_dir(out_dir)
    png = os.path.join(out, stem + ".png")
    svg = os.path.join(out, stem + ".svg")
    fig.savefig(png, facecolor=fig.get_facecolor())
    fig.savefig(svg, facecolor=fig.get_facecolor())
    plt.close(fig)

    ex = {
        "scope": SCOPE,
        "renderer_computes_physics": False,
        "renderer_assembles_metrics": False,
        "hessian_coverage": HESSIAN_NOTE,
        "full_hessian_implemented": False,
        "energy_evaluated": False,
        "optimizer_run": False,
        "f1_resolution_certified": False,
        "ready_for_energy_optimization": False,
        "index_ordering": INDEX_ORDER,
        "g00_zero_meaning": CROSSING_NOTE,
        "colour_scale": "fixed; computed once from the whole record",
        "deterministic_rerender": {
            "record_digest_algorithm": "FNV-1a 64-bit",
            "input_digests": {os.path.basename(records.resolve(f)): fnv1a64(f)
                              for f in inputs},
            "matplotlib_backend": "Agg",
            "no_random_state_used": True,
        },
    }
    if source_rows is not None:
        ex["source_rows"] = source_rows
    ex.update(extra or {})

    prov = provenance.build_provenance(
        program="O", milestone="O0B", input_files=inputs,
        simulation_parameters=simulation_parameters(manifest),
        coordinate_frame="Alcubierre coordinates (t, x, y, z); design "
                         "coordinates xi = (x-x_c(t))/R0, s = (y^2+z^2)/R0^2",
        observer="none: O0B produces a metric, not an observer measurement",
        units="geometrised G = c = 1; R0 and lengths in metres",
        time_convention="analytic: d/dt = -(v_s/R0) d/dxi, not frozen time",
        render_theme=theme.name, color_model="illustrative_only",
        camera_model="orthographic 2-D slice",
        frame_count=1, frame_rate=None, extra=ex)
    provenance.write_sidecar(png, prov)
    provenance.write_sidecar(svg, prov)

    if plotted_header is not None:
        records.write_plotted_csv(os.path.join(out, stem + ".plotted.csv"),
                                  plotted_header, plotted_rows or [])
    return png
