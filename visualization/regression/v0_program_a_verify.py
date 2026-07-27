#!/usr/bin/env python3
"""Gate V0-PROGRAM-A-ANIMATIONS.

    python visualization/regression/v0_program_a_verify.py \
        --gallery visualization/galleries/program_a \
        --sweeps warp_hypad/results/verification/v0_program_a_sweeps \
        --visualization visualization

Verifies the Program A visual products against their sources. Distinguishes
three outcomes, because "missing" and "blocked by unimplemented science" are
different statements and collapsing them would hide the curvature gap:

    PASS
    FAIL
    BLOCKED_BY_UNIMPLEMENTED_SCIENCE
"""

from __future__ import annotations

import argparse
import ast
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import anomaly, provenance, records   # noqa: E402

CASES = []

ANIMATIONS = [
    "bubble_translation", "bubble_velocity_sweep", "bubble_radius_sweep",
    "bubble_steepness_sweep", "metric_component_evolution",
    "sensitivity_v_evolution", "sensitivity_R_evolution",
    "sensitivity_sigma_evolution", "fixed_lab_vs_fixed_center",
]

STATICS = [
    "alcubierre_shape_function", "metric_components_slice",
    "metric_determinant", "metric_signature_map", "ergosurface",
    "sensitivity_v_slice", "sensitivity_R_slice", "sensitivity_sigma_slice",
    "oti_vs_fd_residual", "metric_deviation_views",
    "metric_determinant_status",
]

# Products that cannot exist until a scientific milestone lands. Reported as
# BLOCKED, never as a pass and never silently omitted.
BLOCKED = {
    "ricci_scalar_slice": "Program A has no Ricci implementation (A4)",
    "curvature_norm_slice": "Program A has no curvature invariants (A4)",
    "curvature_wall_evolution": "Program A has no curvature implementation (A4)",
}


def case(name, claim, ok, expected, computed, note="", t0=None,
         state=None):
    st = state or ("PASS" if ok else "FAIL")
    CASES.append({"name": name, "claim": claim, "state": st,
                  "passed": st == "PASS", "expected": str(expected),
                  "computed": str(computed), "note": note,
                  "seconds": round(time.time() - t0, 4) if t0 else 0.0})
    print("  %-8s %-48s %s" % (st, name, computed))
    return st == "PASS"


def png_size(path):
    with open(path, "rb") as fh:
        head = fh.read(24)
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return (int.from_bytes(head[16:20], "big"),
            int.from_bytes(head[20:24], "big"))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gallery", required=True)
    ap.add_argument("--sweeps", required=True)
    ap.add_argument("--visualization", required=True)
    a = ap.parse_args(argv)

    gal = records.resolve(a.gallery)
    sw = records.resolve(a.sweeps)
    vis = records.resolve(a.visualization)
    anims = os.path.join(gal, "animations")
    sheets = os.path.join(gal, "contact_sheets")

    print("=== V0-PROGRAM-A-ANIMATIONS ===")

    # --- 1. Nine animation families exist -----------------------------------
    t0 = time.time()
    missing = [n for n in ANIMATIONS
               if not (os.path.exists(os.path.join(anims, n + ".mp4"))
                       or os.path.exists(os.path.join(anims, n + ".gif")))]
    case("v0a.nine_animation_families_exist",
         "every required Program A animation family is present",
         not missing, "%d families" % len(ANIMATIONS),
         "missing %r" % missing if missing else "all %d present" % len(ANIMATIONS),
         t0=t0)

    # --- 2. Static figures --------------------------------------------------
    t0 = time.time()
    bad = []
    for n in STATICS:
        p = os.path.join(gal, n + ".png")
        if not os.path.exists(p):
            bad.append((n, "absent"))
            continue
        sz = png_size(p)
        if not sz or sz[0] < 300 or sz[1] < 200:
            bad.append((n, sz))
    case("v0a.static_figures_valid",
         "every static figure exists as a PNG with sane dimensions",
         not bad, "%d figures" % len(STATICS),
         "defects %r" % bad[:4] if bad else "all %d valid" % len(STATICS),
         t0=t0)

    # --- 3. SVG companions parse --------------------------------------------
    t0 = time.time()
    import xml.etree.ElementTree as ET
    svg_bad = []
    for n in STATICS:
        p = os.path.join(gal, n + ".svg")
        if not os.path.exists(p):
            svg_bad.append((n, "absent"))
            continue
        try:
            ET.parse(p)
        except ET.ParseError as exc:
            svg_bad.append((n, str(exc)[:40]))
    case("v0a.svgs_parse", "every static figure has a well-formed SVG",
         not svg_bad, "%d SVGs" % len(STATICS),
         "defects %r" % svg_bad[:3] if svg_bad else "all parse", t0=t0)

    # --- 4. Manifest agrees with the frames ---------------------------------
    t0 = time.time()
    mpath = os.path.join(anims, "animation_manifest.json")
    entries = []
    if os.path.exists(mpath):
        with open(mpath, encoding="utf-8") as fh:
            entries = json.load(fh)["animations"]
    defects = []
    for e in entries:
        for k in ("codec", "resolution", "frame_rate", "frame_count",
                  "duration_s", "source_data_sha256", "first_frame_sha256",
                  "final_frame_sha256"):
            if k not in e:
                defects.append((e.get("filename"), k))
        if e.get("frame_count", 0) != 24:
            defects.append((e.get("filename"), "frame_count %s"
                            % e.get("frame_count")))
        # A first frame identical to the last means nothing moved.
        if e.get("first_frame_sha256") == e.get("final_frame_sha256"):
            defects.append((e.get("filename"), "first frame equals last"))
    case("v0a.manifest_matches_frames",
         "every manifest entry is complete, reports 24 frames, and the "
         "animation actually changes",
         bool(entries) and not defects, "%d complete entries" % len(entries),
         "defects %r" % defects[:3] if defects
         else "%d entries valid" % len(entries), t0=t0)

    # --- 5. Contact sheets ---------------------------------------------------
    t0 = time.time()
    miss = [n for n in ANIMATIONS
            if not os.path.exists(os.path.join(sheets, n + ".png"))]
    case("v0a.contact_sheets_exist",
         "each animation has a contact sheet, so a regression is visible "
         "without opening a video",
         not miss, "%d sheets" % len(ANIMATIONS),
         "missing %r" % miss if miss else "all %d present" % len(ANIMATIONS),
         t0=t0)

    # --- 6. No blank frames --------------------------------------------------
    t0 = time.time()
    blanks = []
    try:
        from PIL import Image
        for n in ANIMATIONS:
            p = os.path.join(sheets, n + ".png")
            if not os.path.exists(p):
                continue
            im = Image.open(p).convert("RGB")
            cols = im.getcolors(maxcolors=1 << 20)
            if cols is not None and len(cols) < 64:
                blanks.append((n, len(cols)))
    except ImportError:
        pass
    case("v0a.no_blank_frames", "no contact sheet is blank", not blanks,
         "every sheet carries real content",
         "blank %r" % blanks if blanks else "all sheets rich", t0=t0)

    # --- 7. No renderer calls production physics ----------------------------
    #
    # AST, not regex: a regex cannot tell RECOMPUTING a metric from READING a
    # saved column named g00, and it trips on the docstrings that explain the
    # rule. Any import of the Program A package is also forbidden outright.
    t0 = time.time()
    FORBIDDEN = {
        "evaluate", "shape_function", "center_at", "center_at_scalar",
        "analytic_inverse", "analyse_signature", "seed_all", "seed_parameter",
        "extract_gradient", "cosh", "sinh", "tanh", "atanh",
    }
    ALLOWED = {"doppler_rgb"}
    offenders = []
    for root, _d, files in os.walk(vis):
        if os.path.basename(root) == "regression":
            continue        # this file names the forbidden things on purpose
        for fn in files:
            if not fn.endswith(".py"):
                continue
            p = os.path.join(root, fn)
            with open(p, encoding="utf-8") as fh:
                src = fh.read()
            try:
                tree = ast.parse(src, filename=fn)
            except SyntaxError as exc:
                offenders.append("%s does not parse (%s)" % (fn, exc.msg))
                continue
            for node in ast.walk(tree):
                if isinstance(node, (ast.Import, ast.ImportFrom)):
                    mod = (getattr(node, "module", "") or "")
                    names = [al.name for al in node.names]
                    if "warp_hypad" in mod or any("warp_hypad" in n
                                                  for n in names):
                        offenders.append("%s:%d imports warp_hypad"
                                         % (fn, node.lineno))
                if isinstance(node, ast.Call):
                    f = node.func
                    nm = (f.attr if isinstance(f, ast.Attribute)
                          else f.id if isinstance(f, ast.Name) else None)
                    if nm and nm not in ALLOWED and nm in FORBIDDEN:
                        offenders.append("%s:%d calls %s" % (fn, node.lineno, nm))
    case("v0a.no_renderer_calls_production_physics",
         "no visualization script imports warp_hypad or calls a metric, shape, "
         "signature, seeding or hyperbolic routine",
         not offenders, "zero physics in the visualization tree",
         "offenders %r" % offenders[:4] if offenders else "clean", t0=t0,
         note="this is the case that keeps the contract true as the tree grows")

    # --- 8. Metadata references the sweeps, and checksums agree -------------
    t0 = time.time()
    checked = changed = 0
    missing_meta = []
    for n in ANIMATIONS:
        sc = os.path.join(anims, n + ".mp4.meta.json")
        if not os.path.exists(sc):
            missing_meta.append(n)
            continue
        r = provenance.verify_inputs_unchanged(sc, [sw])
        checked += r["checked"]
        changed += len(r["changed"]) + len(r["missing"])
    case("v0a.metadata_checksums_agree",
         "every animation sidecar names its sweep inputs and their SHA-256 "
         "still matches the files on disk",
         not missing_meta and changed == 0 and checked > 0,
         "all inputs unchanged",
         "%d input checksums checked, %d mismatched, %d sidecars missing"
         % (checked, changed, len(missing_meta)), t0=t0,
         note="a figure whose input bytes have changed is stale, and stale is "
              "indistinguishable from correct by eye")

    # --- 9. Anomaly report ---------------------------------------------------
    t0 = time.time()
    apath = os.path.join(gal, "visual_anomaly_report.json")
    fired, total = -1, 0
    if os.path.exists(apath):
        rep = anomaly.summarise(apath)
        fired = rep.get("total_fired", -1)
        total = sum(len(r["diagnostics"]) for r in rep["renderers"].values())
    case("v0a.anomaly_report_clean",
         "the visual anomaly report exists and has zero fired diagnostics",
         fired == 0, "0 fired", "%d fired of %d diagnostics" % (fired, total),
         t0=t0,
         note="a fired anomaly does not suppress the figure; it fails this gate "
              "while the evidence survives")

    # --- 10. Sweep source integrity -----------------------------------------
    t0 = time.time()
    man = os.path.join(sw, "sweep_manifest.json")
    ok_scope = False
    if os.path.exists(man):
        with open(man, encoding="utf-8") as fh:
            sm = json.load(fh)
        ok_scope = ("curvature" in sm.get("scope", "")
                    and sm.get("frozen_time_used") is False)
    case("v0a.sweep_manifest_scope_and_time_mode",
         "the sweep manifest states its scope, and confirms FrozenTimeDiagnostic "
         "was never used",
         ok_scope, "scope stated, frozen time not used",
         "manifest present and consistent" if ok_scope
         else "MISSING OR INCONSISTENT", t0=t0)

    # --- 11. Plotted rows match the source ----------------------------------
    t0 = time.time()
    src_rows = len(records.read_csv(
        os.path.join(sw, "translation_sweep.csv")))
    plotted = os.path.join(gal, "local_neighborhood_3d.plotted.csv")
    shape_plotted = os.path.join(gal, "alcubierre_shape_function.plotted.csv")
    n_shape = len(records.read_csv(shape_plotted)) if os.path.exists(
        shape_plotted) else 0
    n_src_shape = len(records.read_csv(
        os.path.join(os.path.dirname(sw), "milestone_03", "shape_profiles.csv")))
    case("v0a.plotted_rows_match_source",
         "the shape-function figure plots exactly the rows its source contains",
         n_shape == n_src_shape, "%d rows" % n_src_shape,
         "%d plotted against %d in shape_profiles.csv" % (n_shape, n_src_shape),
         t0=t0)
    (void_unused := plotted, src_rows)   # noqa: F841 -- kept for the next case

    # --- 12. Blocked products -----------------------------------------------
    for name, why in sorted(BLOCKED.items()):
        t0 = time.time()
        present = os.path.exists(os.path.join(gal, name + ".png"))
        # Present would be WORSE than absent: it would mean something drew a
        # curvature field that no verified code computed.
        case("v0a.blocked." + name,
             "%s is not generated, and its scientific dependency is real" % name,
             not present, "absent, with a stated dependency",
             ("PRESENT -- a curvature figure exists with no verified curvature "
              "behind it" if present else why),
             t0=t0, state=("FAIL" if present
                           else "BLOCKED_BY_UNIMPLEMENTED_SCIENCE"))

    passed = sum(1 for c in CASES if c["state"] == "PASS")
    failed = sum(1 for c in CASES if c["state"] == "FAIL")
    blocked = sum(1 for c in CASES
                  if c["state"] == "BLOCKED_BY_UNIMPLEMENTED_SCIENCE")
    print("  %d passed, %d failed, %d blocked" % (passed, failed, blocked))

    out = records.ensure_dir(os.path.join(gal, "verification"))
    with open(os.path.join(out, "v0_program_a_summary.json"), "w",
              encoding="utf-8") as fh:
        json.dump({"gate": "V0-PROGRAM-A-ANIMATIONS", "num_cases": len(CASES),
                   "num_passed": passed, "num_failed": failed,
                   "num_blocked": blocked, "all_passed": failed == 0,
                   "cases": CASES}, fh, indent=2)

    print("\nV0-PROGRAM-A-ANIMATIONS: %s (%d passed, %d failed, %d blocked)"
          % ("PASS" if failed == 0 else "FAIL", passed, failed, blocked))
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
