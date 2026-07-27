"""The V0 global metadata contract.

C1's sidecar proved the idea; V0 makes it repository-wide and adds the two fields
that matter across programs:

  * `program` and `milestone`, so a figure states which codebase produced it;
  * `input_files` and `input_checksums`, so a figure can be traced to the exact
    BYTES it was drawn from, not merely to a filename that may since have been
    regenerated.

The checksum is the important addition. C1 could rely on a single solver run
writing a single metadata file; across five programs a figure may be drawn from
records produced months apart, and a filename is not identity.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import os
import subprocess
from typing import Any, Dict, Iterable, List, Optional

from . import records

# The full field list the V0 contract requires. Kept here so a renderer cannot
# quietly omit one: build_provenance() fills every entry and write_sidecar()
# refuses a dict that is missing any.
REQUIRED_FIELDS = (
    "program",
    "milestone",
    "git_commit",
    "git_dirty",
    "input_files",
    "input_checksums",
    "simulation_parameters",
    "coordinate_frame",
    "observer",
    "units",
    "time_convention",
    "render_theme",
    "color_model",
    "camera_model",
    "frame_count",
    "frame_rate",
    "generated_timestamp",
)

# The closed vocabulary, unchanged from C1. An RGB mapping is never presented as
# a physically predicted observed colour, and the only way to say what a figure's
# colour means is to pick one of these.
VALID_COLOR_MODELS = records.VALID_COLOR_MODELS


def file_sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(records.resolve(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def git_state(repo_dir: Optional[str] = None) -> Dict[str, Any]:
    """The commit and dirty flag, read from git itself.

    Read here rather than inherited from a solver's metadata, because a V0
    renderer may draw records committed long ago while the working tree has moved
    on. The figure should say what the tree looked like when the FIGURE was made.
    """
    cwd = repo_dir or os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
    def run(*args):
        try:
            r = subprocess.run(["git", *args], cwd=cwd, capture_output=True,
                               text=True, timeout=30)
            return r.stdout.strip() if r.returncode == 0 else ""
        except Exception:                  # noqa: BLE001 -- git may be absent
            return ""
    commit = run("rev-parse", "HEAD") or "unknown"
    status = run("status", "--porcelain")
    return {"git_commit": commit, "git_dirty": bool(status.strip())}


def build_provenance(*, program: str, milestone: str, input_files: Iterable[str],
                     simulation_parameters: Dict[str, Any],
                     coordinate_frame: str, observer: str, units: str,
                     time_convention: str, render_theme: str, color_model: str,
                     camera_model: str, frame_count: int = 1,
                     frame_rate: Optional[float] = None,
                     extra: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Assemble a complete V0 sidecar.

    Refuses an undeclared colour model outright. "Never present an illustrative
    RGB mapping as a physically predicted observed colour" is only enforceable if
    the vocabulary is closed, so it is closed here rather than by convention.
    """
    if color_model not in VALID_COLOR_MODELS:
        raise records.RecordError(
            "color_model %r is not one of %r. An RGB mapping must declare what "
            "it means." % (color_model, list(VALID_COLOR_MODELS)))

    files: List[str] = [records.resolve(f) for f in input_files]
    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        raise records.RecordError(
            "cannot checksum inputs that do not exist: %r" % missing)

    prov: Dict[str, Any] = {
        "program": program,
        "milestone": milestone,
        "input_files": [os.path.basename(f) for f in files],
        "input_checksums": {os.path.basename(f): file_sha256(f) for f in files},
        "simulation_parameters": simulation_parameters,
        "coordinate_frame": coordinate_frame,
        "observer": observer,
        "units": units,
        "time_convention": time_convention,
        "render_theme": render_theme,
        "color_model": color_model,
        "camera_model": camera_model,
        "frame_count": int(frame_count),
        "frame_rate": frame_rate,
        "generated_timestamp": datetime.datetime.now(
            datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        # Stated as data, so the V0 gate can assert it rather than trust a
        # comment. The AST check is what actually enforces it.
        "recomputes_physics": False,
        "contract": ("production solver -> saved records -> visualization "
                     "scripts -> PNG/SVG/GIF/MP4"),
    }
    prov.update(git_state())
    if extra:
        prov.update(extra)
    return prov


def write_sidecar(figure_path: str, provenance: Dict[str, Any]) -> str:
    """Write <figure>.meta.json, refusing an incomplete contract."""
    missing = [k for k in REQUIRED_FIELDS if k not in provenance]
    if missing:
        raise records.RecordError(
            "V0 figure provenance is missing %r. Every field is required: a "
            "figure whose frame or units are unstated cannot be checked against "
            "anything." % missing)
    path = records.resolve(figure_path) + ".meta.json"
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(provenance, fh, indent=2, sort_keys=True)
    return path


def verify_inputs_unchanged(sidecar_path: str,
                            search_dirs: Iterable[str]) -> Dict[str, Any]:
    """Re-hash a sidecar's inputs and report agreement.

    Used by the V0 gate: a figure whose input bytes have changed since it was
    drawn is stale, and stale is indistinguishable from correct by eye.
    """
    with open(records.resolve(sidecar_path), encoding="utf-8") as fh:
        prov = json.load(fh)
    stored = prov.get("input_checksums", {})
    result = {"checked": 0, "matched": 0, "missing": [], "changed": []}
    dirs = [records.resolve(d) for d in search_dirs]
    for name, digest in stored.items():
        found = None
        for d in dirs:
            cand = os.path.join(d, name)
            if os.path.exists(cand):
                found = cand
                break
        result["checked"] += 1
        if found is None:
            result["missing"].append(name)
            continue
        if file_sha256(found) == digest:
            result["matched"] += 1
        else:
            result["changed"].append(name)
    result["ok"] = not result["missing"] and not result["changed"]
    return result
