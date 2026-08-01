#!/usr/bin/env python3
"""V0-O0C-CONVERGED-ENERGY-VISUALS: the O0C visual gate.

    python v0_o0c_verify.py --gallery <dir> --records <dir> --visualization <dir>

What this gate is actually for.  O0C is the first milestone in the optimization
program that computes an ENERGY, and the first that draws in three dimensions.
Both are places where a picture can say more than the data supports, so the gate
checks the two claims that matter:

    the renderers compute no physics -- enforced by an AST scan, not promised;
    every product says what its frame axis is and whether it was revolved.

Plus the ordinary discipline: the records exist, the scientific gate passed, the
figures are not blank, the scales are fixed, the digests have not drifted, and
nothing anywhere claims a reduction in warp-drive energy.
"""

from __future__ import annotations

import argparse
import ast
import json
import os
import struct
import sys
import time
import xml.etree.ElementTree as ET

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

from program_o0c import o0c_common as O                          # noqa: E402


GATE = "V0-O0C-CONVERGED-ENERGY-VISUALS"

# Names a renderer must never call or define. Every one of them is an operation
# that would make Python a second, unverified physics implementation.
FORBIDDEN_CALLS = (
    "christoffel", "riemann", "ricci", "einstein_tensor", "kretschmann_of",
    "stress_energy", "energy_density_of", "nec_minimum", "compute_curvature",
    "assemble_metric", "metric_at", "eulerian_density", "connection",
    "einstein_at", "adm_decompose", "build_tetrad", "project_to_frame",
    "softplus", "evaluate_objectives", "smooth_negative_energy",
)

# Names that are legitimate rendering geometry and must NOT be flagged.
ALLOWED_GEOMETRY = ("isocontour", "revolve", "project", "painter_order")

_RESULTS = []


def case(name, reference, ok, expected, computed, t0=None, note=""):
    _RESULTS.append({
        "name": name, "reference": reference, "passed": bool(ok),
        "expected": expected, "computed": computed, "note": note,
        "runtime_seconds": (time.time() - t0) if t0 else 0.0,
    })
    print("  %-6s %-48s %s" % ("PASS" if ok else "FAIL", name, computed))
    return bool(ok)


def png_size(path):
    with open(path, "rb") as f:
        head = f.read(33)
    if len(head) < 24 or head[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    w, h = struct.unpack(">II", head[16:24])
    return int(w), int(h)


def png_is_rich(path, min_bytes=24000):
    return os.path.getsize(path) >= min_bytes


def read_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def sidecar(gal, stem):
    """One figure's provenance record.

    build_provenance FLATTENS the `extra` dict into the top level, so a reader
    that looks for a nested "extra" key finds nothing and silently passes. Every
    sidecar read goes through this one function so that cannot happen in some
    cases and not others.
    """
    p = os.path.join(gal, stem + ".png.meta.json")
    return read_json(p) if os.path.exists(p) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gallery", required=True)
    ap.add_argument("--records", required=True)
    ap.add_argument("--visualization", required=True)
    a = ap.parse_args()

    gal, rec, vis = a.gallery, a.records, a.visualization
    print("=== %s ===" % GATE)

    # --- records and the scientific gate ------------------------------------
    t0 = time.time()
    missing = [n for n in O.REQUIRED_RECORDS
               if not os.path.exists(os.path.join(rec, n))]
    case("v0o0c.records_exist", "every required O0C record is present",
         not missing, "all %d present" % len(O.REQUIRED_RECORDS),
         "all %d present" % len(O.REQUIRED_RECORDS) if not missing
         else "missing: " + ", ".join(missing[:6]), t0)

    t0 = time.time()
    summary = read_json(os.path.join(rec, "verification_summary.json"))
    man = read_json(os.path.join(rec, "geometry_bridge_manifest.json"))
    case("v0o0c.scientific_gate_passed",
         "the O0C scientific gate passed with zero failures",
         summary["num_failed"] == 0, "0 failed",
         "%d/%d passed, %d failed, %d skipped"
         % (summary["num_passed"], summary["num_cases"],
            summary["num_failed"], summary["num_skipped"]), t0)

    # --- the claims the milestone must not make ------------------------------
    t0 = time.time()
    case("v0o0c.no_optimizer_run",
         "no optimizer ran anywhere in O0C",
         man["optimizer_run"] is False, "optimizer_run false",
         "optimizer_run = %s" % man["optimizer_run"], t0)

    t0 = time.time()
    case("v0o0c.no_production_finite_difference_geometry",
         "the production geometry route contains no finite difference",
         man["production_finite_differences_in_geometry"] == 0,
         "zero", "%d, with the note: %s"
         % (man["production_finite_differences_in_geometry"],
            man["finite_difference_note"][:90]), t0)

    t0 = time.time()
    rd = read_json(os.path.join(rec, "o0c_objective_readiness.json"))
    case("v0o0c.readiness_not_overclaimed",
         "no objective is marked ready while a required flag is false",
         rd["any_ready_for_optimization"] is False and
         rd["o1_may_begin"] is False and
         all(not o["ready_for_optimization"] for o in rd["objectives"]),
         "nothing marked ready, O1 may not begin",
         "%d objectives classified, none ready; blocking reason: %s"
         % (len(rd["objectives"]), rd["objectives"][0]["blocking_reason"][:80]),
         t0)

    t0 = time.time()
    adv = O.read_csv(os.path.join(rec, "adversarial_candidate_checks.csv"))
    bad = [r for r in adv if r["verdict"] == "comparable"
           and r["behavior_preserved"] != "1"]
    case("v0o0c.no_behavior_changing_improvement",
         "no behavior-changing or unresolved candidate is treated as an "
         "improvement",
         not bad, "zero such candidates",
         "%d adversarial candidates, %d comparable, %d flagged"
         % (len(adv), sum(1 for r in adv if r["verdict"] == "comparable"),
            sum(1 for r in adv if r["verdict"] != "comparable")), t0)

    t0 = time.time()
    interior = max((abs(float(r["interior_delta_beta"])) for r in adv),
                   default=1.0)
    collar = max((abs(float(r["collar_delta_beta"])) for r in adv), default=1.0)
    case("v0o0c.protected_regions_locked",
         "every adversarial candidate leaves Delta beta^x exactly zero in the "
         "protected interior and the outer collar",
         interior == 0.0 and collar == 0.0, "exactly zero",
         "worst interior %g, worst collar %g over %d candidates"
         % (interior, collar, len(adv)), t0)

    t0 = time.time()
    hay = json.dumps(man) + json.dumps(rd)
    for n in O.STATICS + O.STATICS_3D:
        p = os.path.join(gal, n + ".png.meta.json")
        if os.path.exists(p):
            hay += open(p, encoding="utf-8").read()
    lowered = hay.lower()
    forbidden_phrases = ("energy reduced", "reduced energy", "lower energy "
                         "requirement", "energy requirement reduced",
                         "warp drive is feasible", "faster than light travel "
                         "is possible")
    hits = [p for p in forbidden_phrases if p in lowered]
    case("v0o0c.no_reduced_energy_claim",
         "nothing in the manifests or figure sidecars claims a reduction in "
         "warp-drive energy",
         not hits, "no such phrase",
         "scanned %d characters, 0 hits" % len(hay) if not hits
         else "found: " + ", ".join(hits), t0)

    # --- the renderers compute no physics -------------------------------------
    t0 = time.time()
    offenders = []
    scanned = 0
    for mod in ("program_o0c/o0c_common.py", "program_o0c/plot_o0c_energy.py",
                "program_o0c/animate_o0c_energy.py",
                "program_o0c/explore_o0c_spacetime.py"):
        path = os.path.join(vis, mod)
        if not os.path.exists(path):
            continue
        scanned += 1
        tree = ast.parse(open(path, encoding="utf-8").read())
        for node in ast.walk(tree):
            name = None
            if isinstance(node, ast.Call):
                f = node.func
                name = (f.attr if isinstance(f, ast.Attribute)
                        else getattr(f, "id", None))
            elif isinstance(node, ast.FunctionDef):
                name = node.name
            if not name:
                continue
            low = name.lower()
            if any(g in low for g in ALLOWED_GEOMETRY):
                continue
            for bad_name in FORBIDDEN_CALLS:
                if bad_name in low:
                    offenders.append("%s:%s" % (mod, name))
    case("v0o0c.no_python_gr_computation",
         "no renderer defines or calls a metric, curvature, stress-energy, "
         "observer, NEC or objective evaluation",
         not offenders and scanned == 4, "clean across all 4 modules",
         "clean across %d modules" % scanned if not offenders
         else "offenders: " + ", ".join(sorted(set(offenders))[:6]), t0)

    t0 = time.time()
    sidecars = [os.path.join(gal, n + ".png.meta.json")
                for n in O.STATICS + O.STATICS_3D]
    declared = 0
    for p in sidecars:
        if not os.path.exists(p):
            continue
        e = read_json(p)
        if e.get("renderer_computes_physics") is False:
            declared += 1
    case("v0o0c.renderers_declare_no_physics",
         "every figure sidecar declares that the renderer computed no physics",
         declared == len(sidecars), "all %d declare it" % len(sidecars),
         "%d of %d" % (declared, len(sidecars)), t0)

    # --- products -------------------------------------------------------------
    t0 = time.time()
    absent, small = [], []
    for n in O.STATICS + O.STATICS_3D:
        p = os.path.join(gal, n + ".png")
        if not os.path.exists(p):
            absent.append(n)
            continue
        sz = png_size(p)
        if not sz or sz[0] < 600 or sz[1] < 300:
            small.append(n)
    case("v0o0c.static_figures_exist",
         "all thirteen figures are present at a usable size",
         not absent and not small,
         "%d present" % len(O.STATICS + O.STATICS_3D),
         "%d present, %d absent, %d undersized"
         % (len(O.STATICS + O.STATICS_3D) - len(absent), len(absent),
            len(small)), t0)

    t0 = time.time()
    blanks = [n for n in O.STATICS + O.STATICS_3D
              if os.path.exists(os.path.join(gal, n + ".png"))
              and not png_is_rich(os.path.join(gal, n + ".png"))]
    case("v0o0c.pngs_nonblank", "no figure is blank", not blanks,
         "all rich", "all %d rich" % len(O.STATICS + O.STATICS_3D)
         if not blanks else "thin: " + ", ".join(blanks), t0)

    t0 = time.time()
    bad_svg = []
    for n in O.STATICS:
        p = os.path.join(gal, n + ".svg")
        if not os.path.exists(p):
            bad_svg.append(n + " (absent)")
            continue
        try:
            ET.parse(p)
        except Exception as exc:                                 # noqa: BLE001
            bad_svg.append("%s (%s)" % (n, exc))
    case("v0o0c.svgs_parse",
         "every two-dimensional figure has a well-formed SVG",
         not bad_svg, "all %d parse" % len(O.STATICS),
         "all %d parse" % len(O.STATICS) if not bad_svg
         else ", ".join(bad_svg[:4]), t0)

    # --- animations -----------------------------------------------------------
    t0 = time.time()
    am = read_json(os.path.join(gal, "animation_manifest.json"))
    entries = am["animations"]
    case("v0o0c.animation_manifest_parses",
         "the animation manifest lists every animation",
         len(entries) == len(O.ANIMATIONS),
         "%d entries" % len(O.ANIMATIONS), "%d entries" % len(entries), t0)

    t0 = time.time()
    unwritten = [e["name"] for e in entries
                 if not e["gif_written"] or not e["mp4_written"]]
    reasons = [e.get("mp4_unwritten_reason") for e in entries
               if not e["mp4_written"]]
    case("v0o0c.animation_status_honest",
         "every animation is written, or carries a reason it is not",
         all(e["gif_written"] for e in entries) and
         all(e["mp4_written"] or e.get("mp4_unwritten_reason")
             for e in entries),
         "written, or explained",
         "%d written, %d unwritten%s"
         % (len(entries) - len(unwritten), len(unwritten),
            (" (" + str(reasons[0])[:50] + ")") if reasons else ""), t0)

    t0 = time.time()
    expected_axes = dict(O.ANIMATIONS)
    wrong = [e["name"] for e in entries
             if e.get("frame_axis") != expected_axes.get(e["name"])]
    orbit = [e for e in entries if e["name"] == "cinematic_3d_camera_orbit"]
    orbit_ok = bool(orbit) and orbit[0]["frame_axis"] == "camera angle" and \
        orbit[0]["frame_axis_is_physical_time"] is False
    case("v0o0c.animation_axis_labels",
         "every animation names its frame axis, and the camera orbit is NOT "
         "labelled physical time",
         not wrong and orbit_ok, "all axes correct, orbit is camera angle",
         "%d of %d correct; camera orbit frame_axis = '%s', physical time = %s"
         % (len(entries) - len(wrong), len(entries),
            orbit[0]["frame_axis"] if orbit else "?",
            orbit[0]["frame_axis_is_physical_time"] if orbit else "?"), t0)

    t0 = time.time()
    revolved = [e for e in entries if e["axisymmetric_reconstruction"]]
    stated = [e for e in revolved
              if e["reconstruction_statement"] == O.RECONSTRUCTION_NOTE]
    case("v0o0c.reconstruction_statement_present",
         "every revolved product carries the axisymmetric-reconstruction "
         "statement",
         len(stated) == len(revolved) and revolved,
         "all revolved products state it",
         "%d of %d revolved animations state it" % (len(stated), len(revolved)),
         t0)

    t0 = time.time()
    vmeta = read_json(os.path.join(rec, "volume_field_3d.meta.json"))
    sampled_sidecars = 0
    for n in O.STATICS_3D:
        p = os.path.join(gal, n + ".png.meta.json")
        if not os.path.exists(p):
            continue
        e = read_json(p)
        if e.get("reconstruction_statement") in (O.RECONSTRUCTION_NOTE,
                                                 O.SAMPLED_NOTE):
            sampled_sidecars += 1
    case("v0o0c.sampled_volume_distinguished",
         "the independently sampled production volume is recorded and is not "
         "confused with a revolved reconstruction",
         vmeta["axisymmetric_reconstruction"] is False and
         sampled_sidecars == len(O.STATICS_3D),
         "the volume record declares itself sampled; every 3-D figure declares "
         "which kind it is",
         "volume kind '%s', %d of %d 3-D figures declare their kind"
         % (vmeta["kind"], sampled_sidecars, len(O.STATICS_3D)), t0)

    t0 = time.time()
    sheets = [e for e in entries if e.get("contact_sheet")]
    present = [e for e in sheets
               if os.path.exists(os.path.join(gal, "contact_sheets",
                                              e["contact_sheet"]))]
    case("v0o0c.contact_sheets_exist",
         "each animation has a contact sheet a reviewer can read in a diff",
         len(present) == len(entries), "all %d present" % len(entries),
         "%d of %d present" % (len(present), len(entries)), t0)

    t0 = time.time()
    leftovers = [d for d in os.listdir(gal)
                 if d.startswith("_frames_") and
                 os.path.isdir(os.path.join(gal, d))]
    case("v0o0c.no_raw_frame_directory",
         "no raw frame directory is left in the gallery",
         not leftovers, "none", "none" if not leftovers
         else ", ".join(leftovers[:4]), t0)

    # --- fixed scales and reproducibility -------------------------------------
    t0 = time.time()
    scales = read_json(os.path.join(rec, "field_scales.json"))
    checked = [sidecar(gal, n) for n in O.STATICS + O.STATICS_3D]
    checked = [c for c in checked if c is not None]
    fixed_declared = (len(checked) == len(O.STATICS + O.STATICS_3D) and
                      all(c.get("colour_scale", "").startswith("fixed")
                          for c in checked))
    case("v0o0c.fixed_scales_declared",
         "colour scales and isovalues are computed once by production and "
         "frozen, and every sidecar says so",
         scales["policy"].startswith("fixed") and fixed_declared and
         all(e.get("fixed_scales") for e in entries),
         "fixed everywhere",
         "production policy: '%s'; all sidecars and animation entries agree: %s"
         % (scales["policy"][:44], fixed_declared), t0)

    t0 = time.time()
    drift = []
    for n in O.STATICS + O.STATICS_3D:
        p = os.path.join(gal, n + ".png.meta.json")
        if not os.path.exists(p):
            continue
        e = read_json(p)
        digests = e.get("deterministic_rerender", {}).get("input_digests", {})
        if not digests:
            drift.append("%s:no digests recorded" % n)
        for base, digest in digests.items():
            src = os.path.join(rec, base)
            if os.path.exists(src) and O.fnv1a64(src) != digest:
                drift.append("%s:%s" % (n, base))
    case("v0o0c.checksums_agree",
         "every figure's recorded input digest still matches the record it was "
         "drawn from",
         not drift, "no drift",
         "no drift" if not drift else "drifted: " + ", ".join(drift[:5]), t0)

    t0 = time.time()
    case("v0o0c.index_ordering_recorded",
         "the spacetime index ordering is stated with the data",
         man.get("hessian_limitation") and
         scales.get("index_ordering") == O.INDEX_ORDER,
         "0=t, 1=x, 2=y, 3=z", "scales: '%s'; animations: '%s'"
         % (scales.get("index_ordering"), am.get("index_ordering")), t0)

    t0 = time.time()
    case("v0o0c.hessian_limitation_carried",
         "the Hessian limitation travels with the manifest and the figures",
         man["hessian_limitation"] == O.HESSIAN_NOTE,
         "the exact sentence", man["hessian_limitation"][:80] + "...", t0)

    t0 = time.time()
    case("v0o0c.scope_statement_displayed",
         "the O0C scope statement appears in the manifest and in every "
         "figure sidecar",
         man["scope"] == O.SCOPE and
         len(checked) == len(O.STATICS + O.STATICS_3D) and
         all(c.get("scope") == O.SCOPE for c in checked),
         "present in the manifest and in all %d sidecars"
         % len(O.STATICS + O.STATICS_3D),
         "manifest and %d of %d sidecars carry it verbatim"
         % (sum(1 for c in checked if c.get("scope") == O.SCOPE),
            len(O.STATICS + O.STATICS_3D)), t0)

    # --- the spacetime explorer -----------------------------------------------
    t0 = time.time()
    ex = read_json(os.path.join(gal, "explorer_manifest.json"))
    figs = [n for n in ex["figures"]
            if os.path.exists(os.path.join(gal, n + ".png"))]
    anims = [e for e in ex["animations"]
             if e["gif_written"] and
             os.path.exists(os.path.join(gal, e["name"] + ".gif"))]
    case("v0o0c.explorer_products_exist",
         "the spacetime explorer wrote its figures and animations",
         len(figs) == len(ex["figures"]) and len(anims) == len(ex["animations"]),
         "%d figures, %d animations" % (len(ex["figures"]),
                                        len(ex["animations"])),
         "%d of %d figures, %d of %d animations"
         % (len(figs), len(ex["figures"]), len(anims),
            len(ex["animations"])), t0)

    t0 = time.time()
    # Only the WORLDTUBE view carries worldtube properties. null_cone_tipping
    # is an ordinary two-dimensional plot and correctly declares neither, so
    # requiring it to would be testing the wrong thing.
    worldtube_figs = ["spacetime_explorer_4d"]
    side = [sidecar(gal, n) for n in worldtube_figs]
    side = [c for c in side if c is not None]
    declared = [c for c in side
                if c.get("is_four_spatial_dimensions") is False]
    case("v0o0c.worldtube_is_not_four_spatial_dimensions",
         "the explorer is described as a spacetime worldtube with a TIME axis, "
         "never as a view of four spatial dimensions",
         ex["is_four_spatial_dimensions"] is False and
         ex["view_kind"] == "spacetime_worldtube" and
         "TIME axis" in ex["worldtube_note"] and
         all(e.get("is_four_spatial_dimensions") is False
             for e in ex["animations"]) and
         len(declared) == len(worldtube_figs),
         "declared a worldtube everywhere, four spatial dimensions nowhere",
         "manifest view_kind '%s'; %d of %d worldtube figure sidecars and "
         "%d of %d animation entries declare it"
         % (ex["view_kind"], len(declared), len(worldtube_figs),
            sum(1 for e in ex["animations"]
                if e.get("is_four_spatial_dimensions") is False),
            len(ex["animations"])), t0)

    t0 = time.time()
    third = [e for e in ex["animations"] if e.get("third_axis") ==
             "coordinate time"]
    case("v0o0c.time_axis_named",
         "the third axis of the worldtube is named as coordinate time",
         len(third) == len(ex["animations"]) and
         all(c.get("third_axis") == "coordinate time" for c in side),
         "coordinate time, everywhere",
         "%d of %d animation entries and %d of %d worldtube sidecars name it"
         % (len(third), len(ex["animations"]),
            sum(1 for c in side if c.get("third_axis") == "coordinate time"),
            len(side)), t0)

    t0 = time.time()
    cone_man = read_json(os.path.join(rec, "null_cone_manifest.json"))
    disclaimed = (cone_man["geodesics_integrated"] is False and
                  cone_man["causal_reachability_claimed"] is False and
                  cone_man["travel_time_claimed"] is False and
                  cone_man["ftl_motion_claimed"] is False and
                  ex["geodesics_integrated"] is False and
                  ex["causal_reachability_claimed"] is False and
                  ex["travel_time_claimed"] is False and
                  ex["ftl_motion_claimed"] is False)
    cone_side = sidecar(gal, "null_cone_tipping")
    case("v0o0c.null_cones_make_no_causal_claim",
         "the causal overlay is local metric algebra, with geodesic "
         "integration and every causal claim explicitly disclaimed",
         disclaimed and cone_side is not None and
         cone_side.get("geodesics_integrated") is False and
         cone_side.get("ftl_motion_claimed") is False,
         "no geodesics, no reachability, no travel time, no FTL",
         "production and renderer manifests and the figure sidecar all "
         "disclaim it; production derives the cone as %s"
         % cone_man["derivation"][:64], t0)

    t0 = time.time()
    all_axes = set()
    for e in entries + ex["animations"]:
        all_axes.add(e.get("frame_axis"))
    required = {"physical time", "camera angle", "correction amplitude",
                "quadrature refinement"}
    orbit2 = [e for e in ex["animations"]
              if e["name"] == "spacetime_worldtube_orbit"]
    case("v0o0c.four_frame_axes_distinguished",
         "physical time, camera motion, correction amplitude and quadrature "
         "refinement are all present and all distinguished",
         required.issubset(all_axes) and
         bool(orbit2) and orbit2[0]["frame_axis_is_physical_time"] is False and
         orbit2[0]["camera_moves"] is True,
         "all four axes present and separately labelled",
         "axes present: %s; worldtube orbit is camera-angle with "
         "physical_time=%s" % (", ".join(sorted(all_axes)),
                               orbit2[0]["frame_axis_is_physical_time"]
                               if orbit2 else "?"), t0)

    # --- report ---------------------------------------------------------------
    passed = sum(1 for r in _RESULTS if r["passed"])
    failed = len(_RESULTS) - passed
    print("  %d passed, %d failed" % (passed, failed))
    print("  acceptance: 0 Python GR computations, 0 optimizer runs, "
          "0 readiness overclaims, 0 unlabelled frame axes, "
          "0 energy-reduction claims")

    out = os.path.join(gal, "verification")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "visual_summary.json"), "w",
              encoding="utf-8") as f:
        json.dump({"gate": GATE, "cases": _RESULTS, "num_cases": len(_RESULTS),
                   "num_passed": passed, "num_failed": failed,
                   "scope": O.SCOPE,
                   "not_an_energy_claim": O.NOT_AN_ENERGY_CLAIM}, f, indent=2)

    print("\n%s: %s (%d passed, %d failed)"
          % (GATE, "PASS" if failed == 0 else "FAIL", passed, failed))
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
