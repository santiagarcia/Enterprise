#!/usr/bin/env python3
"""V0-O0D-F1W-CERTIFICATION-VISUALS: the O0D visual gate.

    python v0_o0d_verify.py --gallery <dir> --records <dir> --visualization <dir>

O0D is the milestone that decides whether an optimizer may run at all, so the
gate's job is narrow and specific: check that the readiness verdict shown in the
figures is the one the record actually contains, that the surrogate's value is
never presented as an energy, and that the renderer computed no physics.
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
from program_o0d import plot_o0d_certification as P              # noqa: E402

GATE = "V0-O0D-F1W-CERTIFICATION-VISUALS"

FORBIDDEN = ("christoffel", "riemann", "ricci", "einstein_at", "stress_energy",
             "eulerian_density", "nec_minimum", "adm_decompose",
             "build_tetrad", "project_to_frame", "smooth_negative_energy",
             "evaluate_objectives", "connection_jet", "metric_jet")
ALLOWED = ("isocontour", "revolve", "project", "painter_order")

_R = []


def case(name, ref, ok, expected, computed, t0=None):
    _R.append({"name": name, "reference": ref, "passed": bool(ok),
               "expected": expected, "computed": computed,
               "runtime_seconds": (time.time() - t0) if t0 else 0.0})
    print("  %-6s %-52s %s" % ("PASS" if ok else "FAIL", name, computed))
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
    case("v0o0d.records_exist", "every required O0D record is present",
         not missing, "all %d present" % len(P.REQUIRED_RECORDS),
         "all %d present" % len(P.REQUIRED_RECORDS) if not missing
         else "missing: " + ", ".join(missing[:6]), t0)

    t0 = time.time()
    summ = read_json(os.path.join(rec, "verification_summary.json"))
    man = read_json(os.path.join(rec, "o0d_manifest.json"))
    rd = read_json(os.path.join(rec, "o0d_objective_readiness.json"))
    case("v0o0d.scientific_gate_passed",
         "the O0D scientific gate passed with zero failures",
         summ["num_failed"] == 0, "0 failed",
         "%d/%d passed, %d failed, %d skipped"
         % (summ["num_passed"], summ["num_cases"], summ["num_failed"],
            summ["num_skipped"]), t0)

    t0 = time.time()
    case("v0o0d.no_optimizer_run", "no optimizer ran in O0D",
         man["optimizer_run"] is False, "false",
         "optimizer_run = %s" % man["optimizer_run"], t0)

    t0 = time.time()
    case("v0o0d.no_production_finite_difference_geometry",
         "the production geometry route contains no finite difference",
         man["production_finite_differences_in_geometry"] == 0, "zero",
         "%d" % man["production_finite_differences_in_geometry"], t0)

    t0 = time.time()
    case("v0o0d.no_energy_reduction_claim",
         "O0D claims no reduction in energy",
         man["energy_reduction_claimed"] is False and
         rd["energy_reduction_claimed"] is False, "false in both records",
         "manifest %s, readiness %s" % (man["energy_reduction_claimed"],
                                        rd["energy_reduction_claimed"]), t0)

    # --- the verdict shown must be the verdict recorded ----------------------
    t0 = time.time()
    ready = sidecar(gal, "o0d_readiness")
    consistent = (ready is not None and
                  ready.get("o1_may_begin") == rd["o1_may_begin"])
    case("v0o0d.readiness_figure_matches_the_record",
         "the readiness figure's verdict is the one the record contains",
         consistent, "identical",
         "record o1_may_begin = %s, figure sidecar = %s"
         % (rd["o1_may_begin"],
            ready.get("o1_may_begin") if ready else "absent"), t0)

    t0 = time.time()
    all_true = all(bool(v) for v in rd["flags"].values())
    case("v0o0d.readiness_is_the_conjunction_of_its_flags",
         "o1_may_begin is true only when every required flag is true",
         rd["o1_may_begin"] == all_true, "verdict equals the conjunction",
         "%d flags, all true = %s, o1_may_begin = %s; %s"
         % (len(rd["flags"]), all_true, rd["o1_may_begin"],
            rd["blocking_reason"][:70]), t0)

    # --- the surrogate must never be presented as an energy -------------------
    t0 = time.time()
    sd = sidecar(gal, "surrogate_domain_divergence")
    case("v0o0d.surrogate_value_not_presented_as_an_energy",
         "the figure showing the surrogate's domain divergence declares that "
         "its value is not an energy",
         sd is not None and sd.get("surrogate_value_is_an_energy") is False,
         "declared false",
         "surrogate_value_is_an_energy = %s"
         % (sd.get("surrogate_value_is_an_energy") if sd else "absent"), t0)

    t0 = time.time()
    dom = O.read_csv(os.path.join(rec, "f1w_domain_convergence.csv"))
    raw = [float(r["value"]) for r in dom if r["objective"] == "E_neg_tau"]
    shifted = [float(r["value"]) for r in dom
               if r["objective"] == "E_neg_tau_offset_removed"]
    grows = len(raw) >= 2 and raw[-1] > raw[0] * 1.5
    settles = (len(shifted) >= 2 and
               abs(shifted[-1] - shifted[-2]) < 1e-6 * max(abs(shifted[-1]), 1e-9)
               or abs(shifted[-1] - shifted[-2]) < 1e-9)
    case("v0o0d.surrogate_divergence_is_recorded_not_hidden",
         "the raw surrogate's growth and the offset-removed surrogate's "
         "convergence are both in the record",
         grows and settles, "raw grows, offset-removed settles",
         "raw %.4f -> %.4f (grows), offset-removed %.8f -> %.8f (settles)"
         % (raw[0], raw[-1], shifted[0], shifted[-1]), t0)

    # --- the renderer computes no physics --------------------------------------
    t0 = time.time()
    offenders, scanned = [], 0
    for mod in ("program_o0d/plot_o0d_certification.py",):
        p = os.path.join(vis, mod)
        if not os.path.exists(p):
            continue
        scanned += 1
        tree = ast.parse(open(p, encoding="utf-8").read())
        for node in ast.walk(tree):
            nm = None
            if isinstance(node, ast.Call):
                f = node.func
                nm = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", None)
            elif isinstance(node, ast.FunctionDef):
                nm = node.name
            if not nm:
                continue
            low = nm.lower()
            if any(g in low for g in ALLOWED):
                continue
            for bad in FORBIDDEN:
                if bad in low:
                    offenders.append("%s:%s" % (mod, nm))
    case("v0o0d.no_python_gr_computation",
         "the renderer defines and calls no GR evaluation",
         not offenders and scanned == 1, "clean",
         "clean across %d module" % scanned if not offenders
         else "offenders: " + ", ".join(sorted(set(offenders))[:5]), t0)

    # --- products ---------------------------------------------------------------
    t0 = time.time()
    absent = [n for n in P.FIGURES
              if not os.path.exists(os.path.join(gal, n + ".png"))]
    thin = [n for n in P.FIGURES
            if os.path.exists(os.path.join(gal, n + ".png"))
            and os.path.getsize(os.path.join(gal, n + ".png")) < 20000]
    case("v0o0d.figures_exist", "all six figures are present and not blank",
         not absent and not thin, "%d present, none thin" % len(P.FIGURES),
         "%d present, %d absent, %d thin"
         % (len(P.FIGURES) - len(absent), len(absent), len(thin)), t0)

    t0 = time.time()
    bad = []
    for n in P.FIGURES:
        p = os.path.join(gal, n + ".svg")
        if not os.path.exists(p):
            bad.append(n + " (absent)")
            continue
        try:
            ET.parse(p)
        except Exception as exc:                                 # noqa: BLE001
            bad.append("%s (%s)" % (n, exc))
    case("v0o0d.svgs_parse", "every figure has a well-formed SVG", not bad,
         "all %d parse" % len(P.FIGURES),
         "all %d parse" % len(P.FIGURES) if not bad else ", ".join(bad[:3]), t0)

    t0 = time.time()
    cards = [sidecar(gal, n) for n in P.FIGURES]
    cards = [c for c in cards if c is not None]
    declared = [c for c in cards
                if c.get("renderer_computes_physics") is False and
                c.get("scope") == P.SCOPE]
    case("v0o0d.sidecars_declare_scope_and_no_physics",
         "every sidecar carries the O0D scope and declares no physics computed",
         len(declared) == len(P.FIGURES),
         "all %d" % len(P.FIGURES),
         "%d of %d" % (len(declared), len(P.FIGURES)), t0)

    t0 = time.time()
    drift = []
    for n in P.FIGURES:
        c = sidecar(gal, n)
        if not c:
            continue
        digests = c.get("deterministic_rerender", {}).get("input_digests", {})
        if not digests:
            drift.append("%s: no digests" % n)
        for base, d in digests.items():
            src = os.path.join(rec, base)
            if os.path.exists(src) and O.fnv1a64(src) != d:
                drift.append("%s:%s" % (n, base))
    case("v0o0d.checksums_agree",
         "each figure's recorded input digest still matches its record",
         not drift, "no drift",
         "no drift" if not drift else "drifted: " + ", ".join(drift[:4]), t0)

    t0 = time.time()
    case("v0o0d.scope_statement",
         "the O0D scope statement is present and exact",
         man["scope"] == P.SCOPE, "verbatim", P.SCOPE[:96] + "...", t0)

    passed = sum(1 for r in _R if r["passed"])
    failed = len(_R) - passed
    print("  %d passed, %d failed" % (passed, failed))
    print("  acceptance: 0 Python GR computations, 0 optimizer runs, "
          "0 readiness overclaims, 0 surrogate-as-energy presentations")
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
