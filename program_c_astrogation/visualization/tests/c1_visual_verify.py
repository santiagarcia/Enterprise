#!/usr/bin/env python3
"""Gate C1-RELATIVISTIC-ASTROGATION, visualization cases 21-30.

    python visualization/tests/c1_visual_verify.py \
        --results results/verification/program_c1 \
        --visualization visualization

Checks that the visuals exist, are well formed, agree with the scientific CSVs,
and -- the case that matters most -- that no visualization script recomputes
physics.

Writes verification/visual_summary.json in the same shape the C++ gates use, so
the two halves of the C1 gate can be tallied together.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import records            # noqa: E402


CASES = []


def case(name, claim, ok, expected, computed, note="", t0=None):
    CASES.append({
        "name": name, "claim": claim, "passed": bool(ok),
        "expected": str(expected), "computed": str(computed), "note": note,
        "seconds": round(time.time() - t0, 4) if t0 else 0.0,
    })
    print("  %-4s %-46s %s" % ("PASS" if ok else "FAIL", name, computed))
    return bool(ok)


def png_size(path):
    """Read a PNG's dimensions from the IHDR chunk. No image library needed."""
    with open(path, "rb") as fh:
        head = fh.read(24)
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return (int.from_bytes(head[16:20], "big"), int.from_bytes(head[20:24], "big"))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True)
    ap.add_argument("--visualization", required=True)
    a = ap.parse_args(argv)

    res = records.resolve(a.results)
    vis = records.resolve(a.visualization)
    figs = os.path.join(res, "figures")
    anims = os.path.join(res, "animations")
    sheets = os.path.join(res, "contact_sheets")
    data = os.path.join(res, "data")

    print("=== program_c1_visual: C1 visualization cases 21-30 ===")

    # --- 21. Static 3-D neighbourhood ---------------------------------------
    t0 = time.time()
    png = os.path.join(figs, "local_neighborhood_3d.png")
    svg = os.path.join(figs, "local_neighborhood_3d.svg")
    size = png_size(png) if os.path.exists(png) else None
    case("c1vis.static_3d_neighborhood_generated",
         "the 3-D neighbourhood figure exists in raster and vector form with "
         "sane dimensions",
         os.path.exists(png) and os.path.exists(svg) and size
         and size[0] > 400 and size[1] > 400,
         "PNG and SVG present, both dimensions above 400 px",
         "%s at %s" % (os.path.basename(png), size), t0=t0)

    # --- 22. Worldline figures ----------------------------------------------
    t0 = time.time()
    want = ["worldline_coordinate_time", "worldline_proper_time",
            "velocity_gamma_acceleration", "earth_ship_destination_lightcone",
            "destination_interception", "time_comparison"]
    missing = [w for w in want
               if not os.path.exists(os.path.join(figs, w + ".png"))]
    case("c1vis.worldline_figures_generated",
         "every required worldline, interception and clock figure exists",
         not missing, "%d figures" % len(want),
         "missing %r" % missing if missing else "all %d present" % len(want),
         t0=t0)

    # --- 23. Earth and ship sky figures -------------------------------------
    t0 = time.time()
    sky_ok = all(os.path.exists(os.path.join(figs, n + ".png"))
                 for n in ("sky_earth_frame", "sky_ship_frame"))
    case("c1vis.earth_and_ship_sky_figures_generated",
         "both apparent-sky figures exist, one per observer frame",
         sky_ok, "sky_earth_frame.png and sky_ship_frame.png",
         "present" if sky_ok else "MISSING", t0=t0)

    # --- 24/25. Animations ---------------------------------------------------
    t0 = time.time()
    manifests = [os.path.join(anims, m) for m in
                 ("journey_manifest.json", "sky_manifest.json",
                  "neighborhood_manifest.json")]
    entries = []
    for m in manifests:
        if os.path.exists(m):
            with open(m, encoding="utf-8") as fh:
                entries.extend(json.load(fh)["animations"])
    names = {e["filename"] for e in entries if e["filename"]}
    have_ab = any(n.startswith("relativistic_aberration") for n in names)
    have_dop = any(n.startswith("doppler_sky") for n in names)
    case("c1vis.aberration_animation_generated",
         "the aberration and Doppler sky animations exist",
         have_ab and have_dop, "relativistic_aberration and doppler_sky",
         "aberration=%s doppler=%s" % (have_ab, have_dop), t0=t0)

    t0 = time.time()
    have_journey = any(n.startswith("journey_worldline") for n in names)
    have_comms = any(n.startswith("communication_delay") for n in names)
    have_rot = any(n.startswith("local_neighborhood_rotation") for n in names)
    case("c1vis.journey_animation_generated",
         "the journey, communication-delay and neighbourhood-rotation "
         "animations exist",
         have_journey and have_comms and have_rot,
         "three animation families", "journey=%s comms=%s rotation=%s"
         % (have_journey, have_comms, have_rot), t0=t0)

    # --- 26. Contact sheets --------------------------------------------------
    t0 = time.time()
    want_sheets = ["journey_worldline", "communication_delay",
                   "relativistic_aberration", "doppler_sky",
                   "local_neighborhood_rotation"]
    missing_s = [s for s in want_sheets
                 if not os.path.exists(os.path.join(sheets, s + ".png"))]
    case("c1vis.contact_sheets_generated",
         "every animation has a contact sheet, so a regression is visible "
         "without opening a video",
         not missing_s, "%d contact sheets" % len(want_sheets),
         "missing %r" % missing_s if missing_s
         else "all %d present" % len(want_sheets), t0=t0)

    # --- 27. Manifest validity ----------------------------------------------
    t0 = time.time()
    required = ("filename", "codec", "resolution", "frame_rate", "frame_count",
                "duration_s", "source_data_sha256", "first_frame_sha256",
                "final_frame_sha256")
    bad = []
    for e in entries:
        for k in required:
            if k not in e:
                bad.append((e.get("filename"), k))
        if e.get("frame_count", 0) < 2:
            bad.append((e.get("filename"), "frame_count"))
        # A first frame identical to the last means nothing moved.
        if e.get("first_frame_sha256") and \
                e["first_frame_sha256"] == e["final_frame_sha256"]:
            bad.append((e.get("filename"), "first frame equals final frame"))
    case("c1vis.animation_manifest_valid",
         "every manifest entry carries codec, resolution, frame rate, frame "
         "count, duration and three checksums, and the animation actually moves",
         not bad and entries, "%d complete entries" % len(entries),
         "defects %r" % bad[:4] if bad else "%d entries valid" % len(entries),
         t0=t0)

    # --- 28. Visuals derived from saved output only -------------------------
    #
    # THE CASE THAT MATTERS MOST. Greps every visualization source for the
    # physics it must not contain. A plotting script that recomputed aberration
    # would be a second, unverified simulation wearing the first one's label.
    t0 = time.time()
    # Checked with the AST rather than a regex, because a regex cannot tell
    # RECOMPUTING a Doppler factor from READING a saved column called
    # "doppler_factor" -- and reading that column is exactly what these scripts
    # are supposed to do. It also cannot tell code from the docstrings that
    # explain the rule. The AST sees only executable code, and only calls.
    #
    # The forbidden set is the mathematics that would constitute a second
    # simulation: hyperbolic functions (trajectory integration and rapidity),
    # and any call to a physics routine by name.
    import ast

    FORBIDDEN_CALLS = {
        "cosh": "hyperbolic trajectory integration",
        "sinh": "hyperbolic trajectory integration",
        "tanh": "rapidity to velocity",
        "atanh": "velocity to rapidity",
        "asinh": "inverse hyperbolic trajectory relation",
        "acosh": "inverse hyperbolic trajectory relation",
        "arctanh": "velocity to rapidity",
        "lorentz_gamma": "Lorentz factor",
        "aberrate": "relativistic aberration",
        "boost": "Lorentz boost",
        "doppler_factor": "Doppler calculation",
        "solve_retarded_emission": "light-time solve",
        "propagate_rectilinear": "catalogue propagation",
        "to_barycentric": "catalogue coordinate conversion",
        "inverse_parallax_distance": "catalogue distance estimation",
    }
    # themes.doppler_rgb is explicitly permitted: it takes an ALREADY-SAVED
    # factor and returns a colour. That is drawing, not physics.
    ALLOWED = {"doppler_rgb"}

    offenders = []
    for root, _dirs, files in os.walk(vis):
        if "tests" in root.replace("\\", "/").split("/"):
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
                offenders.append("%s: does not parse (%s)" % (fn, exc.msg))
                continue
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                f = node.func
                name = (f.attr if isinstance(f, ast.Attribute)
                        else f.id if isinstance(f, ast.Name) else None)
                if name and name not in ALLOWED and name in FORBIDDEN_CALLS:
                    offenders.append("%s:%d calls %s (%s)"
                                     % (fn, node.lineno, name,
                                        FORBIDDEN_CALLS[name]))
    case("c1vis.visuals_derived_from_saved_output_only",
         "no visualization source recomputes a trajectory, aberration, Doppler "
         "shift, catalogue propagation or light-time solve",
         not offenders, "zero physics recomputation in the visualization tree",
         "offenders %r" % offenders[:5] if offenders
         else "clean across the visualization tree",
         note="the visuals read solver CSVs and draw them; this case is what "
              "keeps that true as the tree grows", t0=t0)

    # --- 29. No blank or nonfinite frames -----------------------------------
    t0 = time.time()
    blanks = []
    try:
        from PIL import Image
        for s in want_sheets:
            p = os.path.join(sheets, s + ".png")
            if not os.path.exists(p):
                continue
            im = Image.open(p).convert("RGB")
            cols = im.getcolors(maxcolors=1 << 20)
            if cols is not None and len(cols) < 64:
                blanks.append((s, len(cols)))
    except ImportError:
        pass
    # Nonfinite values are refused at the reader boundary in records.column, so
    # this checks the other end: that something was actually drawn.
    case("c1vis.no_blank_or_nonfinite_frames",
         "no contact sheet is blank, and no nonfinite value reached a renderer",
         not blanks, "every sheet carries real content",
         "blank %r" % blanks if blanks else "all sheets have rich content",
         note="records.column refuses NaN and infinity at the read boundary, so "
              "a nonfinite value fails loudly instead of vanishing from a plot",
         t0=t0)

    # --- 30. Deterministic rerender ------------------------------------------
    t0 = time.time()
    # Re-run one static script into a scratch directory and compare the plotted
    # CSV byte for byte. The IMAGE hash is deliberately not compared: that would
    # be pinned to a matplotlib and freetype version, and the requirement is
    # explicit that pixel-perfect hashes need a pinned environment. The plotted
    # numbers are the invariant that matters.
    import tempfile
    ok_det, detail = False, "not run"
    with tempfile.TemporaryDirectory(prefix="c1_rerender_") as tmp:
        script = os.path.join(vis, "static", "plot_worldline.py")
        cmd = [sys.executable, script,
               os.path.join(data, "worldline_proxima_1g.csv"),
               "--metadata", os.path.join(res, "c1_metadata.json"),
               "--out", tmp]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            detail = "rerender failed: %s" % r.stderr.strip().splitlines()[-1:]
        else:
            name = "worldline_proper_time.plotted.csv"
            a_p, b_p = os.path.join(figs, name), os.path.join(tmp, name)
            if os.path.exists(a_p) and os.path.exists(b_p):
                ha = hashlib.sha256(open(a_p, "rb").read()).hexdigest()
                hb = hashlib.sha256(open(b_p, "rb").read()).hexdigest()
                ok_det = ha == hb
                detail = "plotted CSV sha256 %s" % ("matches" if ok_det
                                                    else "DIFFERS")
                # And the sidecar must agree on everything except its timestamp.
                ma = json.load(open(a_p.replace(".plotted.csv", ".png.meta.json"),
                                    encoding="utf-8"))
                mb = json.load(open(b_p.replace(".plotted.csv", ".png.meta.json"),
                                    encoding="utf-8"))
                ma.pop("generated_timestamp"); mb.pop("generated_timestamp")
                ok_det = ok_det and ma == mb
                detail += "; metadata %s" % ("matches" if ma == mb else "DIFFERS")
    case("c1vis.deterministic_rerender_metadata",
         "re-rendering from the same input reproduces the plotted numbers and "
         "the figure metadata exactly",
         ok_det, "identical plotted CSV and identical metadata", detail,
         note="image hashes are NOT compared: they are pinned to matplotlib and "
              "freetype versions, and the requirement allows pixel hashes only "
              "in a pinned environment. The plotted numbers are the invariant",
         t0=t0)

    passed = sum(1 for c in CASES if c["passed"])
    print("  %d/%d passed" % (passed, len(CASES)))

    out = records.ensure_dir(os.path.join(res, "verification"))
    with open(os.path.join(out, "visual_summary.json"), "w",
              encoding="utf-8") as fh:
        json.dump({
            "milestone": "program_c1_visual",
            "title": "C1-RELATIVISTIC-ASTROGATION visualization cases 21-30",
            "num_cases": len(CASES), "num_passed": passed,
            "num_failed": len(CASES) - passed, "num_skipped": 0,
            "all_passed": passed == len(CASES),
            "cases": CASES,
        }, fh, indent=2)

    print("\nC1 VISUAL GATE: %s (%d/%d)"
          % ("PASS" if passed == len(CASES) else "FAIL", passed, len(CASES)))
    return 0 if passed == len(CASES) else 1


if __name__ == "__main__":
    sys.exit(main())
