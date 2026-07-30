"""Shared reading, styling and provenance for the O0A spline-basis visuals.

THE CONTRACT FOR THIS MILESTONE.

    C++ gate -> saved records -> these scripts -> PNG/SVG/GIF/MP4

Nothing here evaluates a B-spline, forms a tensor product, assembles a Gram
matrix, factorises anything or projects a constraint. Every number plotted is a
column that `o0a_verify` wrote. Where a figure needs B_ij(xi, s), it reads
`tensor_product_atlas.csv` rather than multiplying two saved univariate columns:
the tensor product IS the basis definition, case 11 is what verifies it, and a
Python multiplication would be a second, unverified implementation of it.

The only arithmetic performed here is presentation -- picking colour limits,
taking the maximum of a saved column for an axis range, and converting a saved
`s` to a radius for a label. The visual gate scans this package's AST for
production names to keep that true as the package grows.

WHAT THESE FIGURES MUST NOT SAY. O0A concludes that neither the F0 validation
family nor the F1 kernel is resolution-certified, and therefore that neither is
ready for optimization. That is a measured result, not a defect, and the banner
helpers below put it on every F0 and resolution product so a reader cannot take a
clean-looking projection for a usable one.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

from common import provenance, records          # noqa: E402


SCOPE = (
    "O0A constructs and verifies a constrained spline coordinate system for "
    "future metric design. It does not yet optimize a spacetime metric, reduce "
    "energy requirements, demonstrate a causal shortcut or establish physical "
    "realizability.")

# Required verbatim on every product that touches the OTI gradient.
HESSIAN_NOTE = (
    "Disjoint three-direction OTI batches recover the complete coefficient "
    "gradient. They do not recover mixed Hessian entries between coefficients "
    "assigned to different batches.")

WHITENING_NOTE = (
    "Whitening makes optimization coordinates globally supported in coefficient "
    "space. Locality is retained only in the underlying spline evaluation basis.")

F0_BANNER = ("verified representation test\n"
             "not resolution-certified\n"
             "not ready for optimization")

# The three knot semantics in play. Named so a figure cannot imply that every
# family shares one boundary convention.
SEMANTICS_LABEL = {
    "clamped": "clamped (multiplicity 6 both ends)",
    "uniform": "uniform (multiplicity 1 both ends)",
    "lower_clamped_open_upper": "lower-clamped, open upper (6 then 1)",
}

EXCLUSION_LABEL = {
    "active": "active design mode",
    "outer_flat_collar": "removed: support meets the outer flat collar",
    "protected_cabin": "removed: support meets the protected cabin",
    "outside_design_region": "lattice padding: outside the design region",
    "symmetry_constraint": "removed: symmetry constraint",
    "moment_constraint": "removed: moment constraint",
    "linearly_dependent": "removed: linearly dependent",
}


# --- Records ----------------------------------------------------------------

REQUIRED_RECORDS = (
    "basis_manifest.json",
    "basis_readiness.json",
    "axial_knot_vector.csv",
    "radial_squared_knot_vector.csv",
    "univariate_basis_profiles.csv",
    "univariate_basis_derivatives.csv",
    "axisymmetric_basis_profiles.csv",
    "tensor_product_atlas.csv",
    "tensor_product_supports.csv",
    "axis_regular_modes.csv",
    "sparse_stencil_statistics.csv",
    "constraint_matrix.csv",
    "constraint_projection.csv",
    "constraint_singular_values.csv",
    "gram_spectrum_raw.csv",
    "gram_spectrum_constrained.csv",
    "gram_spectrum_whitened.csv",
    "basis_normalization.csv",
    "whitening_locality.csv",
    "linalg_fixtures.csv",
    "linalg_symmetric_fixtures.csv",
    "f0_baseline_projection.csv",
    "f0_projection_errors.csv",
    "f0_coefficients.csv",
    "resolution_coverage.csv",
    "coefficient_derivative_checks.csv",
    "tensor_derivative_probes.csv",
    "verification_summary.json",
)


def read_csv(path):
    with open(records.resolve(path), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def read_json(path):
    with open(records.resolve(path), encoding="utf-8") as fh:
        return json.load(fh)


def fnv1a64(path):
    """The same content digest the C++ manifest writes, so the two can be
    compared. FNV-1a is not a cryptographic checksum and is not presented as
    one; it detects a changed record, which is what is needed here."""
    h = 0xCBF29CE484222325
    with open(records.resolve(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            for b in chunk:
                h ^= b
                h = (h * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return "%016x" % h


def sha256(path):
    h = hashlib.sha256()
    with open(records.resolve(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def source_checksums(record_dir):
    out = {}
    for name in sorted(os.listdir(records.resolve(record_dir))):
        p = os.path.join(record_dir, name)
        if os.path.isfile(records.resolve(p)):
            out[name] = fnv1a64(p)
    return out


def group_by(rows, key, cast=str):
    out = {}
    for r in rows:
        out.setdefault(cast(r[key]), []).append(r)
    return out


def col(rows, key, cast=float):
    return [cast(r[key]) for r in rows]


# --- Figure output ----------------------------------------------------------

def simulation_parameters(manifest):
    return {
        "milestone": "O0A",
        "degree": manifest.get("degree"),
        "minimum_basis_count_per_direction":
            manifest.get("minimum_basis_count_per_direction"),
        "design_region": manifest.get("design_region"),
        "scale_R0": manifest.get("scale_R0"),
        "radial_basis_coordinate": manifest.get("radial_basis_coordinate"),
        "f1_axial": manifest.get("f1_axial"),
        "f1_radial": manifest.get("f1_radial"),
        "f1_tensor": manifest.get("f1_tensor"),
        "constraints": manifest.get("constraints"),
        "f0": manifest.get("f0"),
        "oti_directions": manifest.get("oti_directions"),
        "oti_order": manifest.get("oti_order"),
        "gradient_batches": manifest.get("gradient_batches"),
    }


def save_figure(fig, out_dir, stem, *, manifest, inputs, theme, extra=None,
                plotted_header=None, plotted_rows=None, source_rows=None):
    """PNG + SVG + a V0 sidecar for each, plus the numbers that were drawn.

    `plotted_rows` writes the exact values; `source_rows` is the alternative the
    specification allows -- an explicit declaration of which record rows a figure
    consumed, for figures whose plotted set is the whole record.
    """
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
        "renderer_computes_splines": False,
        "hessian_coverage": HESSIAN_NOTE,
        "full_hessian_implemented": False,
        "f0_resolution_certified": False,
        "f0_ready_for_optimization": False,
        "f1_resolution_certified": False,
        "f1_ready_for_optimization": False,
        "colour_scale": "fixed; computed once from the whole record",
        "deterministic_rerender": {
            "record_digest_algorithm": "FNV-1a 64-bit, as written by the C++ "
                                       "manifest",
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
        program="O", milestone="O0A", input_files=inputs,
        simulation_parameters=simulation_parameters(manifest),
        coordinate_frame="design coordinates xi = (x-x_c)/R0 and "
                         "s = (y^2+z^2)/R0^2, both dimensionless",
        observer="none: O0A is a coordinate system, not a spacetime measurement",
        units="R0 in metres; xi and s dimensionless",
        time_convention="static: O0A has no time dependence",
        render_theme=theme.name, color_model="illustrative_only",
        camera_model="orthographic 2-D",
        frame_count=1, frame_rate=None, extra=ex)
    provenance.write_sidecar(png, prov)
    provenance.write_sidecar(svg, prov)

    if plotted_header is not None:
        records.write_plotted_csv(os.path.join(out, stem + ".plotted.csv"),
                                  plotted_header, plotted_rows or [])
    return png


def footer(fig, theme, text, y=0.004):
    fig.text(0.5, y, text, ha="center", va="bottom", fontsize=6.2,
             color=theme.foreground, alpha=0.75, wrap=True)


def banner(fig, theme, text, y=0.945, color=None):
    fig.text(0.010, y, text, ha="left", va="top", fontsize=8.2,
             color=color or theme.uncertainty_color)


def f0_stamp(fig, theme, x=0.985, y=0.955):
    """The three-line F0 status, on every F0 product, in the uncertainty colour
    and inside a box so it cannot be cropped out of a screenshot silently."""
    fig.text(x, y, F0_BANNER, ha="right", va="top", fontsize=8.4,
             color=theme.uncertainty_color, family="monospace",
             linespacing=1.5,
             bbox=dict(facecolor=theme.background,
                       edgecolor=theme.uncertainty_color, alpha=0.9,
                       boxstyle="round,pad=0.5", linewidth=1.1))


def style_axes(ax, theme, xlabel="", ylabel=""):
    ax.set_xlabel(xlabel, fontsize=8, color=theme.foreground)
    ax.set_ylabel(ylabel, fontsize=8, color=theme.foreground)
    ax.tick_params(labelsize=7, colors=theme.foreground)
    ax.set_facecolor(theme.background)
    for s in ax.spines.values():
        s.set_color(theme.grid)


def draw_knots(ax, knots, theme, label=None, y=None):
    for k, t in enumerate(sorted(set(knots))):
        ax.axvline(t, color=theme.grid, lw=0.7, ls=":", alpha=0.8,
                   label=label if k == 0 else None)
    if y is not None:
        for t in sorted(set(knots)):
            ax.plot([t], [y], marker="|", markersize=7, color=theme.grid,
                    linestyle="none")
