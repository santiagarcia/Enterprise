"""Shared reading, scaling and diagnostics for the A6 observer-energy visuals.

THE CONTRACT, RESTATED FOR A6
-----------------------------
    production solver -> saved records -> visualization scripts -> PNG/SVG/GIF

Nothing in this package computes an energy density, a null-energy functional, a
tetrad or a certification. Every number is read from a CSV that
`warp_hypad/tools/a6_energy_sweep.cpp` wrote by calling
`geometry::project_to_frame`, `geometry::nec_minimum`, `geometry::wec/dec/sec`
and `geometry::certify`. The only arithmetic here is PRESENTATION: choosing
colour limits, arranging a grid of saved samples into a 2-D array for imshow,
and subtracting a saved polar angle from pi/2 to place a marker on a Mollweide
axis. The test is whether a plotted quantity could differ from the CSV. None can.

OBSERVER SEMANTICS ARE PART OF THE DATA, NOT THE STYLING
--------------------------------------------------------
Four observer families appear in these records and they are never merged. The
coordinate-stationary observer is admissible only where g_00 < 0; at g_00 = 0 it
is null and beyond it is spacelike, and in both cases the exporter left its
energy-density cell EMPTY. `observer_rows` refuses to return a measured value for
an inadmissible observer, so a renderer cannot plot one by accident -- it would
have to invent the number itself.

CERTIFICATION IS ALSO PART OF THE DATA
--------------------------------------
A measured negative value and a CERTIFIED negative value are different claims.
Every field figure therefore carries the saved per-sample status column, and the
`STATUS_*` tables below are the only vocabulary a renderer may use to say what a
colour means.
"""

from __future__ import annotations

import ast
import hashlib
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

from common import provenance, records          # noqa: E402


# --- The closed status vocabulary -------------------------------------------
#
# Exactly the four values geometry::to_string(ConditionStatus) can produce.
# There is no "pass": a value inside the numerical error envelope is
# indeterminate, and a renderer that drew it as a definite sign would be
# asserting something the science declined to assert.

CERTIFIED_SATISFIED = "certified_satisfied"
CERTIFIED_VIOLATED = "certified_violated"
INDETERMINATE = "indeterminate_within_discretization"
OBSERVER_INVALID = "observer_or_tetrad_invalid"

STATUS_ORDER = (CERTIFIED_SATISFIED, INDETERMINATE, CERTIFIED_VIOLATED,
                OBSERVER_INVALID)

STATUS_INDEX = {s: i for i, s in enumerate(STATUS_ORDER)}

STATUS_COLOR = {
    CERTIFIED_SATISFIED: "#3f7fd0",
    INDETERMINATE: "#8c8f99",
    CERTIFIED_VIOLATED: "#d1495b",
    OBSERVER_INVALID: "#2a2f3a",
}

STATUS_LABEL = {
    CERTIFIED_SATISFIED: "certified satisfied",
    INDETERMINATE: "indeterminate (within envelope)",
    CERTIFIED_VIOLATED: "certified violated",
    OBSERVER_INVALID: "observer or tetrad invalid",
}

# For a SIGN rather than a condition, the same four states read differently:
# certify(rho, envelope) returns CertifiedSatisfied when rho is certifiably
# positive and CertifiedViolated when it is certifiably negative.
SIGN_LABEL = {
    CERTIFIED_SATISFIED: "certified positive",
    INDETERMINATE: "measured, NOT certified",
    CERTIFIED_VIOLATED: "certified negative",
    OBSERVER_INVALID: "observer or tetrad invalid",
}

OBSERVER_FAMILIES = ("eulerian", "coordinate_stationary", "bubble_comoving",
                     "tetrad_boosted")

OBSERVER_STYLE = {
    "eulerian": ("#4ea1ff", "-", "o"),
    "coordinate_stationary": ("#ffb454", "--", "s"),
    "bubble_comoving": ("#7ed957", "-.", "^"),
    "tetrad_boosted": ("#c792ea", ":", "D"),
}

# The standing caveat. On every A6 product, because this is the milestone at
# which a quantity may finally be called an energy density and therefore the one
# most likely to be screenshotted without its context.
SCOPE = (
    "A6 evaluates observer-dependent energy density and classical pointwise "
    "energy conditions for the stress-energy tensor INFERRED from the selected "
    "metric. It does not demonstrate that the source can be generated, that "
    "quantum energy inequalities are satisfied, that the geometry is "
    "dynamically stable, or that operational FTL travel is possible.")

# What g_00 = 0 is, and the four things it is not. Repeated verbatim on the
# frames, because every one of those four is a stronger claim than the data
# supports.
CROSSING_NOTE = (
    "g_00 = 0 is where the coordinate vector d/dt becomes null. It is NOT a "
    "metric degeneracy, NOT a signature failure, NOT a curvature singularity "
    "and NOT automatically a horizon. The Eulerian observer, built from the ADM "
    "normal, remains timelike through it.")


# --- Source files ------------------------------------------------------------

SWEEP_FILES = (
    "a6_energy_manifest.json",
    "energy_velocity_sweep.csv",
    "nec_velocity_sweep.csv",
    "worst_null_directions.csv",
    "observer_comparison.csv",
    "integrated_objectives.csv",
    "resolution_certification.csv",
    "domain_convergence.csv",
    "null_direction_sphere.csv",
    "null_direction_sphere_points.csv",
)


def sha256(path):
    h = hashlib.sha256()
    with open(records.resolve(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def source_checksums(sweep_dir):
    """SHA-256 of every source record, frames included, in a fixed order."""
    out = {}
    for name in SWEEP_FILES:
        p = os.path.join(sweep_dir, name)
        if os.path.exists(p):
            out[name] = sha256(p)
    for name in sorted(os.listdir(records.resolve(sweep_dir))):
        if name.startswith("energy_frame_") and name.endswith(".csv"):
            out[name] = sha256(os.path.join(sweep_dir, name))
    return out


def load_manifest(sweep_dir):
    p = os.path.join(records.resolve(sweep_dir), "a6_energy_manifest.json")
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def frame_files(sweep_dir):
    d = records.resolve(sweep_dir)
    return [os.path.join(d, n) for n in sorted(os.listdir(d))
            if n.startswith("energy_frame_") and n.endswith(".csv")]


def load_frames(sweep_dir):
    """Every z = 0 slice record, in frame order, as lists of row dicts."""
    return [records.read_csv(p) for p in frame_files(sweep_dir)]


# --- Presentation geometry ---------------------------------------------------

def axes_of(rows):
    xs = sorted({float(r["x"]) for r in rows})
    ys = sorted({float(r["y"]) for r in rows})
    return xs, ys


def to_grid(rows, key, xs, ys, cast=float, fill=0.0):
    """Arrange saved samples into a 2-D array for imshow.

    Rearrangement, not computation: every cell is one CSV value placed at the
    coordinates the same row carries.
    """
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}
    g = [[fill] * len(xs) for _ in ys]
    for r in rows:
        g[yi[float(r["y"])]][xi[float(r["x"])]] = cast(r[key])
    return g


def status_grid(rows, key, xs, ys):
    return to_grid(rows, key, xs, ys, cast=lambda s: STATUS_INDEX.get(s, 3),
                   fill=3)


def global_symmetric_limit(frames, key):
    """One colour limit for the WHOLE sweep, symmetric about zero.

    Computed once from every frame and then frozen. A per-frame limit
    renormalises each frame, so a field that grows by an order of magnitude
    across the sweep looks constant while only its noise appears to move.
    """
    m = 0.0
    for rows in frames:
        for r in rows:
            v = abs(float(r[key]))
            if v > m:
                m = v
    return m if m > 0.0 else 1.0


def global_range(frames, key):
    lo, hi = float("inf"), float("-inf")
    for rows in frames:
        for r in rows:
            v = float(r[key])
            lo = min(lo, v)
            hi = max(hi, v)
    return lo, hi


# --- Observer rows -----------------------------------------------------------

def observer_rows(sweep_dir):
    """observer_comparison.csv, with inadmissible observers kept but blank.

    Returns dicts carrying `measured` as a float ONLY when the exporter wrote
    one. An inadmissible observer measures nothing, so `measured` is None and a
    renderer that wanted to plot it would have to invent the value itself --
    which is the whole point of leaving the cell empty in the record.
    """
    rows = records.read_csv(os.path.join(sweep_dir, "observer_comparison.csv"))
    out = []
    for r in rows:
        valid = r["observer_valid"] == "1"
        raw = r["rho_measured"].strip()
        measured = float(raw) if (valid and raw) else None
        if valid and measured is None:
            raise records.RecordError(
                "observer %r at frame %s is marked valid but carries no "
                "measured energy density" % (r["observer_type"], r["frame"]))
        if not valid and raw:
            raise records.RecordError(
                "observer %r at frame %s is marked INVALID yet carries a "
                "measured energy density %r. An inadmissible observer measures "
                "nothing." % (r["observer_type"], r["frame"], raw))
        out.append({
            "frame": int(r["frame"]),
            "v_s": float(r["v_s"]),
            "site": r["site"],
            "g00": float(r["g00"]),
            "observer_type": r["observer_type"],
            "valid": valid,
            "reason": r["observer_rejection_reason"],
            "norm_residual": float(r["observer_normalization_residual"]),
            "tetrad_residual": float(r["tetrad_orthonormality_residual"]),
            "measured": measured,
            "crosscheck": (float(r["eulerian_rho_crosscheck_residual"])
                           if r["eulerian_rho_crosscheck_residual"].strip()
                           else None),
        })
    return out


# --- Self-scan for recomputed physics ---------------------------------------
#
# The V0 gate runs this check independently over the whole visualization tree.
# It ALSO runs here, so the anomaly report carries a diagnostic that can
# actually fire rather than a self-declaration nobody tested.

_FORBIDDEN_CALLS = frozenset({
    "evaluate", "shape_function", "center_at", "center_at_scalar",
    "analytic_inverse", "analyse_signature", "seed_all", "seed_parameter",
    "extract_gradient", "cosh", "sinh", "tanh", "atanh",
    "nec_value", "nec_minimum", "project_to_frame", "adm_decompose",
    "build_tetrad", "certify", "wec", "dec", "sec",
})


def scan_python_physics(paths):
    """Names of files that call production physics or import warp_hypad."""
    offenders = []
    for p in paths:
        with open(records.resolve(p), encoding="utf-8") as fh:
            src = fh.read()
        try:
            tree = ast.parse(src, filename=os.path.basename(p))
        except SyntaxError as exc:
            offenders.append("%s does not parse (%s)"
                             % (os.path.basename(p), exc.msg))
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                mod = getattr(node, "module", "") or ""
                names = [al.name for al in node.names]
                if "warp_hypad" in mod or any("warp_hypad" in n for n in names):
                    offenders.append("%s:%d imports warp_hypad"
                                     % (os.path.basename(p), node.lineno))
            if isinstance(node, ast.Call):
                f = node.func
                nm = (f.attr if isinstance(f, ast.Attribute)
                      else f.id if isinstance(f, ast.Name) else None)
                if nm and nm in _FORBIDDEN_CALLS:
                    offenders.append("%s:%d calls %s"
                                     % (os.path.basename(p), node.lineno, nm))
    return offenders


def package_python_files():
    return [os.path.join(_HERE, n) for n in sorted(os.listdir(_HERE))
            if n.endswith(".py")]


# --- Figure output -----------------------------------------------------------

def simulation_parameters(meta):
    return {
        "metric": "Alcubierre",
        "radius": meta.get("radius"),
        "sigma": meta.get("sigma"),
        "grid_resolution": meta.get("grid_resolution"),
        "domain_extent": meta.get("domain_extent"),
        "velocity_sequence": meta.get("velocity_sequence"),
        "exact_onset_frame": meta.get("exact_onset_frame"),
        "first_spatial_contour_frame":
            meta.get("first_spatial_contour_frame"),
        "numerical_error_envelope_floor":
            meta.get("numerical_error_envelope_floor"),
        "certification_safety_factor":
            meta.get("certification_safety_factor"),
        "softplus_width": meta.get("softplus_width"),
        "lambda": 0.0,
    }


def save_figure(fig, out_dir, stem, *, meta, inputs, theme, observer,
                extra=None, plotted_header=None, plotted_rows=None,
                frame_count=1, frame_rate=None):
    """PNG + SVG + the V0 sidecar for both, and the numbers that were drawn."""
    import matplotlib.pyplot as plt

    out = records.ensure_dir(out_dir)
    png = os.path.join(out, stem + ".png")
    svg = os.path.join(out, stem + ".svg")
    fig.savefig(png, facecolor=fig.get_facecolor())
    fig.savefig(svg, facecolor=fig.get_facecolor())
    plt.close(fig)

    ex = {
        "scope": SCOPE,
        "g00_zero_meaning": CROSSING_NOTE,
        "colour_scale": "fixed globally across the whole velocity sweep",
        "certification_status_displayed": True,
        "exact_and_smooth_objectives_distinct": True,
        "renderer_computes_physics": False,
        "source_sweep": os.path.basename(os.path.normpath(
            os.path.dirname(records.resolve(inputs[0])))),
    }
    ex.update(extra or {})

    prov = provenance.build_provenance(
        program="A", milestone="A6", input_files=inputs,
        simulation_parameters=simulation_parameters(meta),
        coordinate_frame="Alcubierre coordinates (t, x, y, z), z = 0 slice",
        observer=observer,
        units="geometrised G = c = 1",
        time_convention="t = 0 slice; the sweep parameter is v_s, not time",
        render_theme=theme.name, color_model="illustrative_only",
        camera_model="orthographic 2-D slice",
        frame_count=frame_count, frame_rate=frame_rate, extra=ex)
    provenance.write_sidecar(png, prov)
    provenance.write_sidecar(svg, prov)

    if plotted_header is not None:
        records.write_plotted_csv(os.path.join(out, stem + ".plotted.csv"),
                                  plotted_header, plotted_rows or [])
    return png
