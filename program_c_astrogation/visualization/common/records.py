"""Reading saved Program C records. THE ONLY WAY VISUALS GET NUMBERS.

THE VISUALIZATION CONTRACT
--------------------------
    production solver -> saved CSV/JSON/frame records -> visualization scripts
                      -> PNG/SVG/GIF/MP4

Plotting and animation scripts read output files. They MUST NOT independently
recompute trajectories, aberration, Doppler shifts, catalogue propagation or any
physical observable. A visualization that recomputes physics is a second,
unverified simulation wearing the first one's label, and when the two disagree
the picture is what people believe.

This module is deliberately the only import path a plotting script needs for
data. It offers readers and no physics. There is no Lorentz factor here, no
aberration, no propagation -- and the gate greps for exactly those symbols in the
visualization tree to keep it that way.

The one arithmetic this module does allow is PRESENTATION geometry: projecting a
3-D unit direction that the solver already computed onto 2-D plot axes, scaling a
marker, mapping a number to a colour. That is drawing, not physics. The test is
whether the quantity could differ from the CSV: a projection of a saved direction
cannot, because the direction came from the file.
"""

from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


# --- Errors -----------------------------------------------------------------

class RecordError(RuntimeError):
    """A record file is missing, malformed, or not what the caller asked for."""


# --- Paths ------------------------------------------------------------------
#
# NO SCRIPT MAY DEPEND ON THE CURRENT WORKING DIRECTORY. Every path a script
# touches arrives as an explicit argument and is resolved to an absolute path
# here, once, so a relative path in an argument still means what the caller meant
# rather than what the shell happened to be in.

def resolve(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path))


def ensure_dir(path: str) -> str:
    p = resolve(path)
    os.makedirs(p, exist_ok=True)
    return p


# --- CSV --------------------------------------------------------------------

def read_csv(path: str) -> List[Dict[str, str]]:
    """Read a solver-written CSV into a list of dicts, strings unconverted."""
    p = resolve(path)
    if not os.path.exists(p):
        raise RecordError(
            "%s does not exist. Visuals are generated from saved solver output, "
            "so run c1_verify before any plotting script." % p)
    with open(p, "r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise RecordError("%s has no rows" % p)
    return rows


def column(rows: List[Dict[str, str]], name: str) -> List[float]:
    """Extract a float column, refusing NaN and infinity at the boundary.

    A nonfinite value must never reach a renderer: matplotlib will silently drop
    the point, the figure will look plausible, and the missing data will be
    invisible. Failing here instead makes it loud.
    """
    if name not in rows[0]:
        raise RecordError("column %r not in %r" % (name, sorted(rows[0])))
    out = []
    for i, r in enumerate(rows):
        v = float(r[name])
        if v != v or v in (float("inf"), float("-inf")):
            raise RecordError(
                "column %r row %d is %r. A nonfinite value must not reach a "
                "renderer: it would be silently dropped and the gap would be "
                "invisible in the figure." % (name, i, v))
        out.append(v)
    return out


def text_column(rows: List[Dict[str, str]], name: str) -> List[str]:
    if name not in rows[0]:
        raise RecordError("column %r not in %r" % (name, sorted(rows[0])))
    return [r[name] for r in rows]


# --- Metadata ---------------------------------------------------------------

# The fields every figure and animation must carry. Fixed here so a script
# cannot quietly omit one.
REQUIRED_PROVENANCE_FIELDS = (
    "git_commit",
    "git_dirty",
    "catalogue",
    "catalogue_release",
    "snapshot_checksum",
    "data_kind",
    "reference_epoch",
    "simulation_epoch",
    "coordinate_frame",
    "observer_frame",
    "units",
    "destination_source_id",
    "trajectory_model",
    "physical_parameters",
    "render_parameters",
    "color_model",
    "generated_timestamp",
)

# The fixed vocabulary. An RGB mapping is NEVER presented as a physically
# predicted observed colour, and the only way to say what a figure's colour means
# is to pick one of these.
VALID_COLOR_MODELS = (
    "physical_spectral_model",
    "approximate_blackbody_visualization",
    "catalogue_photometric_proxy",
    "illustrative_only",
)


@dataclass
class RunMetadata:
    """The solver's c1_metadata.json, plus helpers to derive figure provenance."""

    path: str
    raw: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def load(cls, path: str) -> "RunMetadata":
        p = resolve(path)
        if not os.path.exists(p):
            raise RecordError(
                "%s does not exist. Every figure inherits its provenance from the "
                "solver metadata; there is no default." % p)
        with open(p, "r", encoding="utf-8") as fh:
            return cls(path=p, raw=json.load(fh))

    def get(self, key: str, default: Any = None) -> Any:
        return self.raw.get(key, default)

    def git_commit(self) -> str:
        # The solver nests build info; accept either shape rather than guessing.
        v = self.raw.get("git_commit")
        if isinstance(v, dict):
            return str(v.get("git_commit", "unknown"))
        return str(v) if v is not None else "unknown"

    def git_dirty(self) -> bool:
        v = self.raw.get("git_commit")
        if isinstance(v, dict):
            return bool(v.get("git_dirty", True))
        return bool(self.raw.get("git_dirty", True))

    def destination_source_id(self, index: int = 0) -> str:
        d = self.raw.get("destinations") or []
        if index < len(d):
            return str(d[index].get("destination_source_id", "unknown"))
        return "unknown"

    def destination_name(self, index: int = 0) -> str:
        d = self.raw.get("destinations") or []
        if index < len(d):
            return str(d[index].get("display_name", "unknown"))
        return "unknown"

    def figure_provenance(self, *, color_model: str,
                          render_parameters: Dict[str, Any],
                          destination_source_id: Optional[str] = None,
                          observer_frame: Optional[str] = None,
                          trajectory_model: Optional[str] = None
                          ) -> Dict[str, Any]:
        """Build the sidecar dict for one figure or animation.

        Refuses an unlabelled or invented colour model. "Never present an
        illustrative RGB mapping as a physically predicted observed colour" is
        only enforceable if the vocabulary is closed, so it is.
        """
        if color_model not in VALID_COLOR_MODELS:
            raise RecordError(
                "color_model %r is not one of %r. An RGB mapping must declare "
                "what it means." % (color_model, list(VALID_COLOR_MODELS)))
        return {
            "git_commit": self.git_commit(),
            "git_dirty": self.git_dirty(),
            "catalogue": self.get("catalogue", "unknown"),
            "catalogue_release": self.get("catalogue_release", "unknown"),
            "snapshot_checksum": self.get("snapshot_checksum", "unknown"),
            "data_kind": self.get("data_kind", "unknown"),
            "reference_epoch": self.get("reference_epoch"),
            "simulation_epoch": self.get("simulation_epoch"),
            "coordinate_frame": self.get("coordinate_frame", "unknown"),
            "observer_frame": observer_frame or self.get("observer_frame",
                                                         "unknown"),
            "units": self.get("units", "unknown"),
            "destination_source_id": (destination_source_id
                                      or self.destination_source_id(0)),
            "trajectory_model": (trajectory_model
                                 or self.get("trajectory_model", "unknown")),
            "physical_parameters": self.get("physical_parameters", {}),
            "render_parameters": render_parameters,
            "color_model": color_model,
            "generated_timestamp": _utc_now(),
            "generated_from": os.path.basename(self.path),
            "recomputes_physics": False,
        }


def _utc_now() -> str:
    import datetime
    return datetime.datetime.now(datetime.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ")


def write_sidecar(figure_path: str, provenance: Dict[str, Any]) -> str:
    """Write <figure>.meta.json beside a figure. Returns the sidecar path."""
    missing = [k for k in REQUIRED_PROVENANCE_FIELDS if k not in provenance]
    if missing:
        raise RecordError("figure provenance is missing %r" % missing)
    p = resolve(figure_path) + ".meta.json"
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(provenance, fh, indent=2, sort_keys=True)
    return p


def write_plotted_csv(path: str, header: List[str],
                      rows: List[List[Any]]) -> str:
    """Write the EXACT numbers a figure plotted.

    Required by the contract: a reader must be able to check the picture against
    the numbers without rerunning anything. When a figure subsets or reorders the
    solver output, this is what records which subset was drawn.
    """
    p = resolve(path)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    return p
