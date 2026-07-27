"""Visual anomaly diagnostics. THESE ARE TESTS, NOT AESTHETIC SCORES.

Every renderer computes diagnostics on the data it draws and writes
`visual_anomaly_report.json`. The V0 gate reads those reports and fails on a
fired anomaly.

A renderer that finds an anomaly STILL PRODUCES ITS FIGURE. Suppressing the
picture would hide the evidence: the point of drawing a field is to see what is
in it, and a figure withheld because it looked wrong is the one you most needed
to look at. The gate fails; the image survives.

Each diagnostic records its measured value, its threshold and whether it fired,
so a reader can see how close a passing check came to failing. A bare boolean
would hide a margin that had quietly shrunk to nothing.
"""

from __future__ import annotations

import json
import math
import os
from typing import Any, Dict, List, Optional, Sequence

from . import records


class AnomalyReport:
    """Accumulates diagnostics for one renderer, then writes them."""

    def __init__(self, program: str, milestone: str, renderer: str):
        self.program = program
        self.milestone = milestone
        self.renderer = renderer
        self.entries: List[Dict[str, Any]] = []

    def record(self, name: str, *, fired: bool, measured: Any, threshold: Any,
               detail: str = "") -> bool:
        self.entries.append({
            "diagnostic": name,
            "fired": bool(fired),
            "measured": measured,
            "threshold": threshold,
            "detail": detail,
        })
        return bool(fired)

    # --- Field diagnostics --------------------------------------------------

    def check_finite(self, name: str, values: Sequence[float]) -> bool:
        """No NaN, no infinity. Checked at the renderer as well as the reader.

        records.column already refuses non-finite values on read; this catches a
        value produced by a transformation INSIDE the renderer, which the reader
        never sees.
        """
        bad = sum(1 for v in values
                  if v != v or v in (float("inf"), float("-inf")))
        return self.record("blank_field" if not values else "nonfinite_value",
                           fired=bad > 0, measured=bad, threshold=0,
                           detail="%s: %d non-finite of %d"
                                  % (name, bad, len(values)))

    def check_not_blank(self, name: str, values: Sequence[float],
                        min_distinct: int = 4) -> bool:
        """A field that is constant carries no information.

        The classic silent failure: the pipeline runs, the file appears, the
        figure renders, and it is one flat colour because the field was never
        populated.
        """
        distinct = len(set(round(v, 12) for v in values))
        return self.record("blank_field", fired=distinct < min_distinct,
                           measured=distinct, threshold=min_distinct,
                           detail="%s has %d distinct values" % (name, distinct))

    def check_clipped(self, name: str, values: Sequence[float],
                      lo: float, hi: float, max_fraction: float = 0.02) -> bool:
        """Too much of the field sitting exactly at a colour-scale limit.

        Clipping is sometimes intended, so this reports a FRACTION against a
        threshold rather than any occurrence: a colour bar that saturates on 40%
        of the domain is hiding the structure it claims to show.
        """
        if not values:
            return self.record("clipped_values", fired=False, measured=0.0,
                               threshold=max_fraction, detail=name)
        # A CONSTANT field has no clipping, because it has no range. With
        # hi == lo every sample trivially satisfies both "v <= lo" and "v >= hi",
        # so the naive count reports 100% clipped on a perfectly good field --
        # which it did, on the Alcubierre determinant, whose analytic value is
        # exactly -1 everywhere. Constancy is a real property and is checked by
        # check_not_blank; it is not clipping, and conflating the two turns a
        # correct result into a spurious failure.
        if not (hi > lo):
            return self.record("clipped_values", fired=False, measured=0.0,
                               threshold=max_fraction,
                               detail="%s: constant field (value %g), so "
                                      "clipping is undefined; constancy is "
                                      "checked separately" % (name, lo))
        n = sum(1 for v in values if v <= lo or v >= hi)
        frac = n / len(values)
        return self.record("clipped_values", fired=frac > max_fraction,
                           measured=round(frac, 6), threshold=max_fraction,
                           detail="%s: %.2f%% at or beyond the scale limits"
                                  % (name, 100 * frac))

    def check_symmetry(self, name: str, values_a: Sequence[float],
                       values_b: Sequence[float], tolerance: float) -> bool:
        """Two sequences that a symmetry says must agree.

        The caller supplies the pairing, because only the caller knows which
        symmetry the field is supposed to have.
        """
        if len(values_a) != len(values_b):
            return self.record("symmetry_violation", fired=True,
                               measured="length mismatch", threshold=tolerance,
                               detail=name)
        worst = max((abs(a - b) for a, b in zip(values_a, values_b)),
                    default=0.0)
        return self.record("symmetry_violation", fired=worst > tolerance,
                           measured=worst, threshold=tolerance, detail=name)

    def check_no_discontinuity(self, name: str, values: Sequence[float],
                               max_jump: float) -> bool:
        """No step between adjacent samples larger than max_jump.

        On a smooth field a large adjacent-sample jump means either the grid is
        too coarse to resolve what is there, or the ordering is wrong. Both are
        worth failing on.
        """
        worst, where = 0.0, -1
        for i in range(1, len(values)):
            d = abs(values[i] - values[i - 1])
            if d > worst:
                worst, where = d, i
        return self.record("unexpected_discontinuity", fired=worst > max_jump,
                           measured=worst, threshold=max_jump,
                           detail="%s: largest step at index %d" % (name, where))

    def check_sign(self, name: str, values: Sequence[float], *,
                   expect_negative: bool) -> bool:
        """A field whose sign is fixed by geometry.

        For a Lorentzian metric in four dimensions det(g) < 0 everywhere. A
        positive determinant is not an ugly plot; it is a broken spacetime.
        """
        if expect_negative:
            bad = sum(1 for v in values if v >= 0.0)
            label = "unexpected_negative_determinant"
        else:
            bad = sum(1 for v in values if v < 0.0)
            label = "unexpected_sign"
        return self.record(label, fired=bad > 0, measured=bad, threshold=0,
                           detail="%s: %d of %d samples on the wrong side"
                                  % (name, bad, len(values)))

    def check_flag_column(self, name: str, flags: Sequence[float]) -> bool:
        """A saved validity flag that must be set everywhere it is drawn."""
        bad = sum(1 for f in flags if f == 0)
        return self.record("non_lorentzian_signature", fired=bad > 0,
                           measured=bad, threshold=0,
                           detail="%s: %d of %d samples flagged invalid"
                                  % (name, bad, len(flags)))

    # --- Sequence and animation diagnostics ---------------------------------

    def check_monotonic(self, name: str, values: Sequence[float]) -> bool:
        """Timestamps and frame indices must increase."""
        bad = sum(1 for i in range(1, len(values)) if values[i] <= values[i - 1])
        return self.record("non_monotonic_sequence", fired=bad > 0,
                           measured=bad, threshold=0,
                           detail="%s: %d non-increasing steps" % (name, bad))

    def check_center_matches_trajectory(self, name: str,
                                        rendered: Sequence[float],
                                        trajectory: Sequence[float],
                                        tolerance: float) -> bool:
        """The drawn bubble centre against the saved trajectory.

        This is the assertion that a translation animation is showing motion that
        was computed rather than motion introduced by the renderer.
        """
        if len(rendered) != len(trajectory):
            return self.record("center_trajectory_mismatch", fired=True,
                               measured="length mismatch", threshold=tolerance,
                               detail=name)
        worst = max((abs(a - b) for a, b in zip(rendered, trajectory)),
                    default=0.0)
        return self.record("center_trajectory_mismatch", fired=worst > tolerance,
                           measured=worst, threshold=tolerance, detail=name)

    def check_object_count_stable(self, name: str,
                                  counts: Sequence[int]) -> bool:
        """Object count constant across frames.

        Aberration moves stars; a sweep moves a wall. Neither creates or destroys
        anything, so a changing count means the renderer and the data disagree
        about what is in the scene.
        """
        distinct = sorted(set(counts))
        return self.record("changing_object_count", fired=len(distinct) > 1,
                           measured=distinct[:5], threshold=1, detail=name)

    def check_identifiers_stable(self, name: str,
                                 per_frame_ids: Sequence[Sequence[str]]) -> bool:
        """The same identifiers present in every frame, so identity is stable."""
        if not per_frame_ids:
            return self.record("lost_object_identifiers", fired=False,
                               measured=0, threshold=0, detail=name)
        first = set(per_frame_ids[0])
        lost = max((len(first - set(f)) for f in per_frame_ids), default=0)
        return self.record("lost_object_identifiers", fired=lost > 0,
                           measured=lost, threshold=0,
                           detail="%s: up to %d identifiers missing from a frame"
                                  % (name, lost))

    def check_frame_jump(self, name: str, per_frame_scalar: Sequence[float],
                         max_jump: float) -> bool:
        """No sudden frame-to-frame jump in a quantity that should evolve
        smoothly. Catches a dropped, duplicated or misordered frame."""
        return self.check_no_discontinuity(name, per_frame_scalar, max_jump)

    def check_exposure_stable(self, name: str,
                              per_frame_scale: Sequence[float],
                              tolerance: float = 0.0) -> bool:
        """The colour or marker mapping must not drift between frames.

        A per-frame autoscale makes a static field appear to pulse. Tolerance
        defaults to ZERO because the mapping is supposed to be computed once from
        the whole dataset and frozen; any variation at all is a defect.
        """
        if not per_frame_scale:
            return self.record("unstable_exposure", fired=False, measured=0.0,
                               threshold=tolerance, detail=name)
        spread = max(per_frame_scale) - min(per_frame_scale)
        return self.record("unstable_exposure", fired=spread > tolerance,
                           measured=spread, threshold=tolerance, detail=name)

    def check_no_branch_reordering(self, name: str,
                                   per_frame_sorted: Sequence[Sequence[float]],
                                   tolerance: float) -> bool:
        """Branches (roots, eigenvalues) keep their order between frames.

        A crossing that is real must be labelled as such; a crossing that is an
        artefact of re-sorting each frame independently is a rendering bug that
        makes two curves appear to swap identity.
        """
        worst = 0.0
        for i in range(1, len(per_frame_sorted)):
            a, b = per_frame_sorted[i - 1], per_frame_sorted[i]
            if len(a) != len(b):
                return self.record("branch_reordering", fired=True,
                                   measured="branch count changed",
                                   threshold=tolerance, detail=name)
            worst = max(worst, max((abs(x - y) for x, y in zip(a, b)),
                                   default=0.0))
        return self.record("branch_reordering", fired=worst > tolerance,
                           measured=worst, threshold=tolerance, detail=name)

    def check_bound(self, name: str, values: Sequence[float], *,
                    lower: Optional[float] = None,
                    upper: Optional[float] = None) -> bool:
        """A physical bound that must hold, e.g. a passivity margin.

        A HARD-feasible result must never violate this. A soft-penalty result
        may, and when it does the figure must label the violation rather than
        smooth it away -- so this fires, the gate records it, and the renderer
        draws the violation in the uncertainty colour.
        """
        bad = 0
        for v in values:
            if lower is not None and v < lower:
                bad += 1
            if upper is not None and v > upper:
                bad += 1
        return self.record("bound_violation", fired=bad > 0, measured=bad,
                           threshold=0,
                           detail="%s: %d of %d outside [%s, %s]"
                                  % (name, bad, len(values), lower, upper))

    # --- Output -------------------------------------------------------------

    @property
    def fired(self) -> List[Dict[str, Any]]:
        return [e for e in self.entries if e["fired"]]

    def write(self, out_dir: str,
              filename: str = "visual_anomaly_report.json") -> str:
        """Write or MERGE into the report.

        Merges, because several renderers write into one results directory and a
        plain overwrite would leave only the last one's diagnostics -- silently
        discarding the others, which is exactly the class of failure this module
        exists to catch.
        """
        path = os.path.join(records.ensure_dir(out_dir), filename)
        existing: Dict[str, Any] = {"renderers": {}}
        if os.path.exists(path):
            try:
                with open(path, encoding="utf-8") as fh:
                    existing = json.load(fh)
            except (OSError, ValueError):
                existing = {"renderers": {}}
        existing.setdefault("renderers", {})[self.renderer] = {
            "program": self.program,
            "milestone": self.milestone,
            "diagnostics": self.entries,
            "num_fired": len(self.fired),
        }
        total = sum(r["num_fired"] for r in existing["renderers"].values())
        existing["total_fired"] = total
        existing["all_clear"] = total == 0
        existing["note"] = ("These are TESTS, not aesthetic scores. A fired "
                            "diagnostic fails the V0 gate. Figures are still "
                            "written when one fires, because suppressing the "
                            "picture would hide the evidence.")
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(existing, fh, indent=2, sort_keys=True)
        return path


def summarise(path: str) -> Dict[str, Any]:
    """Read a report back, for the gate."""
    with open(records.resolve(path), encoding="utf-8") as fh:
        return json.load(fh)
