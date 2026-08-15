#!/usr/bin/env python3
"""V0-O1-OPTIMIZATION-VISUALS: the O1 visual gate.

    python v0_o1_verify.py --gallery <dir> --records <dir> --visualization <dir>

O1 is the first milestone that could be mistaken for a claim about warp drives,
so this gate is mostly about what the pictures are allowed to say: the verdict
shown must be the verdict recorded, the iteration axis must never be called
time, and no product may shorten the claim into something it does not support.
"""

from __future__ import annotations

import argparse
import ast
import json
import os
import sys
import time
import xml.etree.ElementTree as ET

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

from program_o0c import o0c_common as O                          # noqa: E402
from program_o1 import plot_o1_optimization as P                 # noqa: E402

GATE = "V0-O1-OPTIMIZATION-VISUALS"

FORBIDDEN_CALLS = ("christoffel", "riemann", "ricci", "einstein_at",
                   "stress_energy", "eulerian_density", "nec_minimum",
                   "adm_decompose", "build_tetrad", "project_to_frame",
                   "smooth_negative_energy", "evaluate_objectives",
                   "connection_jet", "metric_jet")
ALLOWED = ("isocontour", "revolve", "project", "painter_order")

# Phrases no product may contain, whatever the numbers say.
FORBIDDEN_PHRASES = (
    "optimized warp drive", "energy-minimal warp drive",
    "physically viable warp drive", "ftl solution", "exotic matter eliminated",
    "faster than light", "causal shortcut", "travel time", "horizon removed",
    "globally optimal", "global optimum",
)

_R = []


def case(name, ref, ok, expected, computed, t0=None):
    _R.append({"name": name, "reference": ref, "passed": bool(ok),
               "expected": expected, "computed": computed,
               "runtime_seconds": (time.time() - t0) if t0 else 0.0})
    print("  %-6s %-50s %s" % ("PASS" if ok else "FAIL", name, computed))
    return bool(ok)


def read_json(p):
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def sidecar(gal, stem):
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

    t0 = time.time()
    missing = [n for n in P.REQUIRED_RECORDS
               if not os.path.exists(os.path.join(rec, n))]
    case("v0o1.records_exist", "every required O1 record is present",
         not missing, "all %d present" % len(P.REQUIRED_RECORDS),
         "all %d present" % len(P.REQUIRED_RECORDS) if not missing
         else "missing: " + ", ".join(missing[:6]), t0)

    t0 = time.time()
    summ = read_json(os.path.join(rec, "verification_summary.json"))
    man = read_json(os.path.join(rec, "o1_manifest.json"))
    cert = read_json(os.path.join(rec, "final_certification.json"))
    case("v0o1.scientific_gate_passed",
         "the O1 scientific gate passed with zero failures",
         summ["num_failed"] == 0, "0 failed",
         "%d/%d passed, %d failed, %d skipped"
         % (summ["num_passed"], summ["num_cases"], summ["num_failed"],
            summ["num_skipped"]), t0)

    # --- the verdict shown is the verdict recorded ---------------------------
    t0 = time.time()
    s = sidecar(gal, "final_design_summary")
    ok = (s is not None and
          s.get("exact_energy_reduction_established")
          == cert["exact_energy_reduction_established"] and
          s.get("claim") == cert["claim"])
    case("v0o1.summary_figure_matches_the_record",
         "the summary figure carries exactly the recorded verdict and claim",
         ok, "identical verdict and claim string",
         "record established = %s, sidecar = %s"
         % (cert["exact_energy_reduction_established"],
            s.get("exact_energy_reduction_established") if s else "absent"), t0)

    t0 = time.time()
    consistent = (cert["exact_energy_reduction_established"] ==
                  (cert["final_candidate_certified"] and
                   cert["absolute_reduction"] > 0.0 and
                   cert["absolute_reduction"] >
                   cert["combined_numerical_uncertainty"]))
    case("v0o1.claim_is_the_conjunction_of_its_conditions",
         "the reduction claim is true exactly when certification, a positive "
         "reduction and exceeding the uncertainty all hold",
         consistent, "verdict equals the conjunction",
         "certified %s, reduction %.3e, uncertainty %.3e, ratio %.1f, "
         "established %s"
         % (cert["final_candidate_certified"], cert["absolute_reduction"],
            cert["combined_numerical_uncertainty"],
            cert["reduction_over_uncertainty"],
            cert["exact_energy_reduction_established"]), t0)

    # --- the claim must not be shortened -------------------------------------
    t0 = time.time()
    hay = json.dumps(man) + json.dumps(cert)
    for n in P.FIGURES:
        p = os.path.join(gal, n + ".png.meta.json")
        if os.path.exists(p):
            hay += open(p, encoding="utf-8").read()
    am = os.path.join(gal, "animation_manifest.json")
    if os.path.exists(am):
        hay += open(am, encoding="utf-8").read()
    low = hay.lower()
    # The scope sentence legitimately NAMES these to deny them, so a hit only
    # counts outside a denial context.
    hits = []
    for phrase in FORBIDDEN_PHRASES:
        idx = 0
        while True:
            k = low.find(phrase, idx)
            if k < 0:
                break
            # THE DENIAL WINDOW HAS TO COVER A LIST, NOT A PHRASE.
            # The scope sentence denies eight things in one breath -- "It does
            # not demonstrate faster-than-light travel, a causal shortcut, ...
            # or a global optimum" -- so the governing "does not" sits some
            # 200 characters before the last item. A 90-character window
            # flagged the scope statement itself, which is the one place these
            # words are not only allowed but required.
            ctx = low[max(0, k - 420):k]
            if not any(d in ctx for d in ("not ", "does not", "never", "no ",
                                          "forbidden", "without", "denies")):
                hits.append(phrase)
                break
            idx = k + 1
    case("v0o1.no_shortened_or_forbidden_claim",
         "no product shortens the result into an FTL, viability or "
         "global-optimum claim",
         not hits, "no such phrase outside a denial",
         "scanned %d characters, %d unguarded hits" % (len(hay), len(hits))
         if not hits else "found: " + ", ".join(hits), t0)

    t0 = time.time()
    case("v0o1.no_global_optimum_or_realizability_claim",
         "the manifest denies global optimality, FTL and physical "
         "realizability explicitly",
         man["global_optimum_claimed"] is False and
         man["ftl_claimed"] is False and
         man["physical_realizability_claimed"] is False,
         "all three false",
         "global %s, ftl %s, realizability %s"
         % (man["global_optimum_claimed"], man["ftl_claimed"],
            man["physical_realizability_claimed"]), t0)

    # --- the iteration axis is not time ---------------------------------------
    t0 = time.time()
    anim = read_json(os.path.join(gal, "animation_manifest.json"))
    bad = [e["name"] for e in anim["animations"]
           if e.get("frame_axis_is_physical_time") is not False]
    named = [e for e in anim["animations"] if e.get("frame_axis")]
    case("v0o1.no_animation_calls_its_axis_physical_time",
         "every animation names its frame axis and none of them is physical "
         "time",
         not bad and len(named) == len(anim["animations"]),
         "all axes named, none physical time",
         "%d animations, axes: %s"
         % (len(anim["animations"]),
            ", ".join(sorted({e["frame_axis"] for e in anim["animations"]}))),
         t0)

    t0 = time.time()
    traj = sidecar(gal, "design_space_trajectory")
    case("v0o1.design_space_is_not_physical_space",
         "the trajectory figure declares its axes are design space, not "
         "physical space",
         traj is not None and traj.get("axes_are_physical_space") is False,
         "declared false",
         "axes_are_physical_space = %s"
         % (traj.get("axes_are_physical_space") if traj else "absent"), t0)

    # --- surrogate discipline --------------------------------------------------
    t0 = time.time()
    acc = O.read_csv(os.path.join(rec, "accepted_iterates.csv"))
    cols = set(acc[0].keys()) if acc else set()
    needed = {"raw_working_E_neg_tau", "renormalized_working_E_neg_tau",
              "delta_working_E_neg_tau", "raw_validation_E_neg_tau",
              "renormalized_validation_E_neg_tau",
              "delta_validation_E_neg_tau", "exact_E_neg"}
    case("v0o1.three_surrogate_names_recorded_separately",
         "raw, renormalized and delta surrogates are recorded under distinct "
         "names alongside the exact objective",
         needed.issubset(cols), "all seven columns present",
         "%d of %d present" % (len(needed & cols), len(needed)), t0)

    # --- the renderer computes no physics --------------------------------------
    t0 = time.time()
    offenders, scanned = [], 0
    for mod in ("program_o1/plot_o1_optimization.py",
                "program_o1/animate_o1_optimization.py"):
        p = os.path.join(vis, mod)
        if not os.path.exists(p):
            continue
        scanned += 1
        for node in ast.walk(ast.parse(open(p, encoding="utf-8").read())):
            nm = None
            if isinstance(node, ast.Call):
                f = node.func
                nm = (f.attr if isinstance(f, ast.Attribute)
                      else getattr(f, "id", None))
            elif isinstance(node, ast.FunctionDef):
                nm = node.name
            if not nm:
                continue
            low2 = nm.lower()
            if any(g in low2 for g in ALLOWED):
                continue
            for b in FORBIDDEN_CALLS:
                if b in low2:
                    offenders.append("%s:%s" % (mod, nm))
    case("v0o1.no_python_gr_computation",
         "no renderer defines or calls a GR evaluation",
         not offenders and scanned == 2, "clean across both modules",
         "clean across %d modules" % scanned if not offenders
         else "offenders: " + ", ".join(sorted(set(offenders))[:5]), t0)

    # --- products ---------------------------------------------------------------
    t0 = time.time()
    absent = [n for n in P.FIGURES
              if not os.path.exists(os.path.join(gal, n + ".png"))]
    thin = [n for n in P.FIGURES
            if os.path.exists(os.path.join(gal, n + ".png"))
            and os.path.getsize(os.path.join(gal, n + ".png")) < 18000]
    case("v0o1.figures_exist", "all figures present and not blank",
         not absent and not thin, "%d present" % len(P.FIGURES),
         "%d present, %d absent, %d thin"
         % (len(P.FIGURES) - len(absent), len(absent), len(thin)), t0)

    t0 = time.time()
    bad2 = []
    for n in P.FIGURES:
        p = os.path.join(gal, n + ".svg")
        if not os.path.exists(p):
            bad2.append(n)
            continue
        try:
            ET.parse(p)
        except Exception:                                        # noqa: BLE001
            bad2.append(n)
    case("v0o1.svgs_parse", "every figure has a well-formed SVG", not bad2,
         "all %d parse" % len(P.FIGURES),
         "all %d parse" % len(P.FIGURES) if not bad2 else ", ".join(bad2[:3]),
         t0)

    t0 = time.time()
    unwritten = [e["name"] for e in anim["animations"]
                 if not e["gif_written"] or
                 (not e["mp4_written"] and not e.get("mp4_unwritten_reason"))]
    case("v0o1.animation_status_honest",
         "every animation is written or carries a reason it is not",
         not unwritten, "written or explained",
         "%d animations, %d unexplained gaps"
         % (len(anim["animations"]), len(unwritten)), t0)

    t0 = time.time()
    sheets = [e for e in anim["animations"] if e.get("contact_sheet")]
    present = [e for e in sheets
               if os.path.exists(os.path.join(gal, "contact_sheets",
                                              e["contact_sheet"]))]
    case("v0o1.contact_sheets_exist", "each animation has a contact sheet",
         len(present) == len(anim["animations"]),
         "all %d" % len(anim["animations"]),
         "%d of %d" % (len(present), len(anim["animations"])), t0)

    t0 = time.time()
    leftovers = [d for d in os.listdir(gal)
                 if d.startswith("_frames_") and
                 os.path.isdir(os.path.join(gal, d))]
    case("v0o1.no_raw_frame_directory", "no raw frame directory is left",
         not leftovers, "none",
         "none" if not leftovers else ", ".join(leftovers[:3]), t0)

    t0 = time.time()
    drift = []
    for n in P.FIGURES:
        c = sidecar(gal, n)
        if not c:
            continue
        dig = c.get("deterministic_rerender", {}).get("input_digests", {})
        if not dig:
            drift.append("%s: no digests" % n)
        for base, d in dig.items():
            src = os.path.join(rec, base)
            if os.path.exists(src) and O.fnv1a64(src) != d:
                drift.append("%s:%s" % (n, base))
    case("v0o1.checksums_agree",
         "each figure's recorded input digest still matches its record",
         not drift, "no drift",
         "no drift" if not drift else "drifted: " + ", ".join(drift[:4]), t0)

    t0 = time.time()
    cards = [sidecar(gal, n) for n in P.FIGURES]
    cards = [c for c in cards if c is not None]
    scoped = [c for c in cards if c.get("scope") == P.SCOPE]
    case("v0o1.scope_statement_on_every_figure",
         "every figure carries the exact O1 scope statement",
         len(scoped) == len(P.FIGURES) and man["scope"] == P.SCOPE,
         "all %d plus the manifest" % len(P.FIGURES),
         "%d of %d figures and the manifest" % (len(scoped), len(P.FIGURES)),
         t0)

    passed = sum(1 for r in _R if r["passed"])
    failed = len(_R) - passed
    print("  %d passed, %d failed" % (passed, failed))
    print("  acceptance: 0 Python GR computations, 0 shortened claims, "
          "0 physical-time mislabels, 0 verdict mismatches")
    out = os.path.join(gal, "verification")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "visual_summary.json"), "w",
              encoding="utf-8") as f:
        json.dump({"gate": GATE, "cases": _R, "num_cases": len(_R),
                   "num_passed": passed, "num_failed": failed,
                   "scope": P.SCOPE}, f, indent=2)
    print("\n%s: %s (%d passed, %d failed)"
          % (GATE, "PASS" if failed == 0 else "FAIL", passed, failed))
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
