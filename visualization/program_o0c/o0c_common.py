"""Shared loading, styling and provenance for the O0C energy visuals.

THE DIVISION OF LABOUR IS THE POINT OF THIS MODULE.

Production writes sampled scalar and vector fields.  This module and the two
renderers that use it read those fields, extract isosurfaces, revolve
meridional half-planes, project a camera and draw.  Nothing here evaluates a
metric, a connection, a curvature, an Einstein tensor, a stress-energy, an
observer quantity, a null-energy minimum or an energy objective.  The visual
gate scans this directory's abstract syntax tree for exactly those operations,
so the separation is enforced rather than promised.

The one arithmetic this module DOES perform is geometry of the drawing:
rotating a half-plane about the axis, interpolating an isocontour, projecting
3-D points onto a 2-D canvas.  That is rendering, not physics, and every
product built by revolution carries the statement that says so.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))

from common import provenance, records          # noqa: E402


SCOPE = (
    "O0C evaluates and verifies curvature, inferred stress-energy and "
    "converged energy-condition objectives for behavior-preserving metric "
    "families. It does not yet optimize those families, demonstrate a causal "
    "shortcut or establish physical realizability."
)

RECONSTRUCTION_NOTE = (
    "axisymmetric 3D reconstruction from production meridional data"
)

SAMPLED_NOTE = (
    "independently sampled (x,y,z) production volume; no revolution"
)

HESSIAN_NOTE = (
    "The complete design gradient is recovered through deterministic OTI "
    "batches. Mixed Hessian entries between coordinates assigned to different "
    "batches are not available."
)

NOT_AN_ENERGY_CLAIM = (
    "These are scientific views of a computed spacetime field. They are not "
    "evidence of physical realizability, of a causal shortcut, or of a "
    "reduction in warp-drive energy."
)

INDEX_ORDER = "0=t, 1=x, 2=y, 3=z"

STATICS = (
    "pointwise_alcubierre_density_reference",
    "energy_objective_decomposition",
    "quadrature_order_convergence",
    "subdivision_convergence",
    "domain_convergence",
    "smooth_surrogate_bias",
    "objective_gradient_agreement",
    "behavior_lock_verification",
    "objective_readiness",
)

STATICS_3D = (
    "wall_energy_3d",
    "nec_violation_3d",
    "curvature_shell_3d",
    "combined_3d_and_slice",
)

ANIMATIONS = (
    ("physical_time_3d_translation", "physical time"),
    ("scientific_3d_energy_shell", "physical time"),
    ("cinematic_3d_camera_orbit", "camera angle"),
    ("correction_amplitude_3d", "correction amplitude"),
    ("quadrature_convergence", "quadrature refinement"),
    ("objective_decomposition", "physical time"),
)

REQUIRED_RECORDS = (
    "geometry_bridge_manifest.json",
    "reused_api_map.json",
    "behavior_contract.json",
    "family_mode_accounting.csv",
    "f1w_basis_manifest.json",
    "pointwise_geometry_checks.csv",
    "source_semantics.json",
    "observer_frame_samples.csv",
    "alcubierre_density_reference.csv",
    "velocity_scaling.csv",
    "energy_condition_samples.csv",
    "null_minimizer_checks.csv",
    "objective_definitions.json",
    "pointwise_objective_samples.csv",
    "objective_values.csv",
    "quadrature_order_convergence.csv",
    "subdivision_convergence.csv",
    "domain_convergence.csv",
    "independent_integration_comparison.csv",
    "legacy_a6_comparison.csv",
    "smooth_surrogate_bias.csv",
    "objective_gradient_checks.csv",
    "gradient_batch_schedule.csv",
    "adversarial_candidate_checks.csv",
    "o0c_objective_readiness.json",
    "volume_field_3d.csv",
    "volume_field_3d.meta.json",
    "meridional_manifest.json",
    "field_scales.json",
    "shift_vectors_sparse.csv",
    "verification_report.md",
    "verification_summary.json",
    "verification_table.csv",
)


# --- loading -----------------------------------------------------------------

def read_csv(path):
    with open(records.resolve(path), newline="", encoding="utf-8") as f:
        text = [ln for ln in f if not ln.startswith("#")]
    return list(csv.DictReader(text))


def read_json(path):
    with open(records.resolve(path), encoding="utf-8") as f:
        return json.load(f)


def fnv1a64(path):
    h = 0xCBF29CE484222325
    with open(records.resolve(path), "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            for b in chunk:
                h ^= b
                h = (h * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return "%016x" % h


def col(rows, key, cast=float):
    return [cast(r[key]) for r in rows if r.get(key) not in (None, "")]


def meridional(rec, frame):
    """One meridional half-plane, as (nx, nr) nested lists of each field."""
    rows = read_csv(os.path.join(rec, "meridional_field_%03d.csv" % frame))
    nx = max(int(r["i"]) for r in rows) + 1
    nr = max(int(r["j"]) for r in rows) + 1
    out = {k: [[0.0] * nr for _ in range(nx)]
           for k in ("x", "rho", "beta_x", "energy_density", "n_min",
                     "kretschmann", "g00")}
    for r in rows:
        i, j = int(r["i"]), int(r["j"])
        for k in out:
            out[k][i][j] = float(r[k])
    out["nx"], out["nr"] = nx, nr
    return out


def volume(rec):
    """The sampled (x, y, z) volume, as a dict of flat arrays plus its shape."""
    rows = read_csv(os.path.join(rec, "volume_field_3d.csv"))
    n = max(int(r["i"]) for r in rows) + 1
    fields = {k: [0.0] * (n * n * n)
              for k in ("x", "y", "z", "beta_x", "rho", "kretschmann", "g00")}
    for r in rows:
        idx = (int(r["i"]) * n + int(r["j"])) * n + int(r["k"])
        for k in fields:
            fields[k][idx] = float(r[k])
    fields["n"] = n
    return fields


# --- rendering geometry (drawing, not physics) --------------------------------

def isocontour(grid, xs, ys, level):
    """Marching-squares line segments of `grid` at `level`.

    Pure interpolation on a saved array.  It knows nothing about what the
    numbers mean, which is exactly the property the visual gate checks for.
    """
    segs = []
    nx, ny = len(xs), len(ys)
    for i in range(nx - 1):
        for j in range(ny - 1):
            c = (grid[i][j], grid[i + 1][j], grid[i + 1][j + 1], grid[i][j + 1])
            px = (xs[i], xs[i + 1], xs[i + 1], xs[i])
            py = (ys[j], ys[j], ys[j + 1], ys[j + 1])
            pts = []
            for e in range(4):
                a, b = c[e], c[(e + 1) % 4]
                if (a - level) * (b - level) >= 0 or a == b:
                    continue
                t = (level - a) / (b - a)
                pts.append((px[e] + t * (px[(e + 1) % 4] - px[e]),
                            py[e] + t * (py[(e + 1) % 4] - py[e])))
            if len(pts) == 2:
                segs.append(pts)
    return segs


def revolve(segs, n_theta=48):
    """Revolve meridional (x, rho) segments into a 3-D surface strip set.

    Every product built with this MUST carry RECONSTRUCTION_NOTE.
    """
    tris = []
    for (x0, r0), (x1, r1) in segs:
        for k in range(n_theta):
            a = 2.0 * math.pi * k / n_theta
            b = 2.0 * math.pi * (k + 1) / n_theta
            ca, sa, cb, sb = math.cos(a), math.sin(a), math.cos(b), math.sin(b)
            tris.append(((x0, r0 * ca, r0 * sa), (x1, r1 * ca, r1 * sa),
                         (x1, r1 * cb, r1 * sb), (x0, r0 * cb, r0 * sb)))
    return tris


def project(points, azimuth, elevation, distance=26.0):
    """Perspective projection onto the canvas. Camera geometry only."""
    ca, sa = math.cos(azimuth), math.sin(azimuth)
    ce, se = math.cos(elevation), math.sin(elevation)
    out = []
    for (x, y, z) in points:
        xr = x * ca + y * sa
        yr = -x * sa + y * ca
        zr = z
        yv = yr * ce - zr * se
        zv = yr * se + zr * ce
        depth = distance - yv
        if depth < 1e-3:
            depth = 1e-3
        f = distance / depth
        out.append((xr * f, zv * f, depth))
    return out


def painter_order(faces, azimuth, elevation, distance=26.0):
    """Back-to-front ordering, so translucent shells composite correctly."""
    keyed = []
    for f in faces:
        pr = project(f, azimuth, elevation, distance)
        d = sum(p[2] for p in pr) / len(pr)
        keyed.append((d, [(p[0], p[1]) for p in pr]))
    keyed.sort(key=lambda t: -t[0])
    return keyed


# --- styling -----------------------------------------------------------------

def style_axes(ax, theme, xlabel="", ylabel=""):
    ax.set_facecolor(theme.background)
    for sp in ax.spines.values():
        sp.set_color(theme.grid)
    ax.tick_params(colors=theme.foreground, labelsize=8)
    ax.xaxis.label.set_color(theme.foreground)
    ax.yaxis.label.set_color(theme.foreground)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=9)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9)
    ax.grid(True, color=theme.grid, alpha=theme.grid_alpha, linewidth=0.5)


def banner(fig, theme, text, y=0.945, color=None):
    fig.text(0.5, y, text, ha="center", va="center", fontsize=9,
             color=color or theme.foreground)


def footer(fig, theme, text, y=0.004):
    fig.text(0.5, y, text, ha="center", va="bottom", fontsize=6.4,
             color=theme.grid, wrap=True)


def readiness_stamp(fig, theme, x=0.986, y=0.962):
    fig.text(x, y, "NOT READY FOR OPTIMIZATION", ha="right", va="center",
             fontsize=7.2, color=theme.uncertainty_color,
             bbox=dict(boxstyle="round,pad=0.28", facecolor=theme.background,
                       edgecolor=theme.uncertainty_color, linewidth=0.7))


def simulation_parameters(manifest):
    return {k: manifest[k] for k in manifest
            if k in ("velocity", "bubble_radius", "wall_sigma", "scale_R0",
                     "lambda", "source_semantics", "f1w_active_modes",
                     "f1w_whitened_dimension", "f1w_oti_batches",
                     "converged_E_neg", "converged_E_net",
                     "production_finite_differences_in_geometry")}


def save_figure(fig, out_dir, stem, *, manifest, inputs, theme, extra=None,
                plotted_header=None, plotted_rows=None, reconstruction=False,
                camera=None, frames=1, frame_axis=None):
    import matplotlib.pyplot as plt

    out = records.ensure_dir(out_dir)
    png = os.path.join(out, stem + ".png")
    fig.savefig(png, facecolor=fig.get_facecolor())
    svg = None
    if stem in STATICS:
        svg = os.path.join(out, stem + ".svg")
        fig.savefig(svg, facecolor=fig.get_facecolor())
    plt.close(fig)

    ex = {
        "scope": SCOPE,
        "renderer_computes_physics": False,
        "renderer_computes_curvature": False,
        "renderer_computes_stress_energy": False,
        "renderer_computes_energy_objectives": False,
        "hessian_coverage": HESSIAN_NOTE,
        "energy_evaluated_by_production": True,
        "optimizer_run": False,
        "ready_for_optimization": False,
        "not_an_energy_claim": NOT_AN_ENERGY_CLAIM,
        "index_ordering": INDEX_ORDER,
        "colour_scale": "fixed; computed once by production over every frame",
        "source_semantics": manifest.get("source_semantics", "total_effective"),
        "lambda": manifest.get("lambda", 0.0),
        "axisymmetric_reconstruction": bool(reconstruction),
        "reconstruction_statement":
            RECONSTRUCTION_NOTE if reconstruction else SAMPLED_NOTE,
        "deterministic_rerender": {
            "record_digest_algorithm": "FNV-1a 64-bit",
            "input_digests": {os.path.basename(records.resolve(f)): fnv1a64(f)
                              for f in inputs},
            "matplotlib_backend": "Agg",
            "no_random_state_used": True,
        },
    }
    if frame_axis:
        ex["frame_axis"] = frame_axis
    ex.update(extra or {})

    prov = provenance.build_provenance(
        program="O", milestone="O0C", input_files=inputs,
        simulation_parameters=simulation_parameters(manifest),
        coordinate_frame="Alcubierre coordinates (t, x, y, z); design "
                         "coordinates xi = (x-x_c(t))/R0, s = (y^2+z^2)/R0^2",
        observer="ADM Eulerian observer n^mu = (1/alpha, -beta^i/alpha)",
        units="geometrised G = c = 1; lengths in metres",
        time_convention="analytic: d/dt = -(v_s/R0) d/dxi, not frozen time",
        render_theme=theme.name, color_model="illustrative_only",
        camera_model=camera or "orthographic 2-D slice",
        frame_count=frames, frame_rate=None, extra=ex)
    provenance.write_sidecar(png, prov)
    if svg:
        provenance.write_sidecar(svg, prov)

    if plotted_header is not None:
        records.write_plotted_csv(os.path.join(out, stem + ".plotted.csv"),
                                  plotted_header, plotted_rows or [])
    return png
