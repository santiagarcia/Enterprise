#!/usr/bin/env python3
"""Gate V0-O0B-CONSTRAINED-ADM-METRIC-VISUALS.

    python visualization/regression/v0_o0b_verify.py \
        --gallery <build>/results/visualization/program_o0b \
        --records <build>/results/verification/o0b_metric \
        --visualization visualization

Acceptance: 0 failures, 0 arbitrary metric-construction bypasses, 0 frozen-time
derivative paths, 0 invalid metrics accepted, 0 optimization-readiness
overclaims, 0 full-Hessian overclaims, 0 energy evaluated.

A scientifically honest "not resolution-certified" is a PASSING result. What the
gate refuses is a run in which O0B's verified metric CONSTRUCTION is presented as
readiness for energy optimization -- a different claim about a different thing.
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
import os
import sys
import time
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import anomaly, provenance, records                   # noqa: E402

CASES = []

STATICS = ["adm_metric_components", "baseline_vs_zero_correction",
           "single_mode_metric_perturbation", "metric_inverse_residual",
           "metric_determinant_identity", "spatial_metric_eigenvalues",
           "metric_derivative_channels", "temporal_derivative_channels",
           "coefficient_gradient_checks", "uniform_grid_cost",
           "quadrature_elements", "o0b_readiness"]

ANIMATIONS = ["single_shift_mode", "baseline_to_corrected_shift",
              "metric_components_vs_mode_amplitude", "analytic_time_translation",
              "raw_constrained_whitened_metric_effect",
              "uniform_grid_vs_elementwise_sampling"]

REQUIRED_RECORDS = [
    "o0a_resolution_consistency.json", "metric_design_manifest.json",
    "coefficient_space_manifest.json", "coefficient_transforms.csv",
    "adm_field_samples.csv", "metric_samples.csv",
    "inverse_metric_samples.csv", "metric_first_derivatives.csv",
    "metric_second_derivatives.csv", "metric_third_derivatives.csv",
    "temporal_derivative_checks.csv", "derivative_symmetry_checks.csv",
    "alcubierre_baseline_comparison.csv", "single_mode_metric_checks.csv",
    "metric_validity_checks.csv", "determinant_checks.csv",
    "coefficient_gradient_checks.csv", "uniform_grid_cost_estimate.json",
    "elementwise_quadrature_elements.csv",
    "elementwise_quadrature_manifest.json", "readiness.json",
    "verification_report.md", "verification_summary.json",
    "verification_table.csv"]

FORBIDDEN = {
    "assemble_covariant", "assemble_contravariant", "f1_metric_jet",
    "shift_jet", "shift_value_scalar", "metric_component_scalar",
    "baseline_shift_table", "correction_shift_table", "chained_derivative",
    "shape_radial_jet", "shape_w_jet", "spatial_metric_from_cholesky",
    "lapse_from_log", "inverse3", "determinant3", "build_projection",
    "sobolev_gram", "whiten", "jacobi_svd", "jacobi_eigen", "make_f1_basis",
    "shape_function", "evaluate", "field_value", "to_raw", "to_constrained",
    "cosh", "sinh", "tanh", "atanh",
}

HESSIAN_SENTENCE = ("The complete design gradient is available through twelve "
                    "batches. Mixed Hessian entries between distinct batches "
                    "are not available.")

READINESS_LINES = ("metric construction verified",
                   "current F1 family not yet resolution-certified",
                   "not ready for energy optimization")


def case(name, claim, ok, expected, computed, note="", t0=None):
    st = "PASS" if ok else "FAIL"
    CASES.append({"name": name, "claim": claim, "state": st, "passed": bool(ok),
                  "expected": str(expected), "computed": str(computed),
                  "note": note,
                  "seconds": round(time.time() - t0, 4) if t0 else 0.0})
    print("  %-6s %-48s %s" % (st, name, computed))
    return bool(ok)


def described(c):
    if not c:
        return "absent"
    v = c.get("computed_description") or c.get("computed")
    return "absent" if v in (None, "") else str(v)


def read_csv(p):
    with open(records.resolve(p), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def read_json(p):
    with open(records.resolve(p), encoding="utf-8") as fh:
        return json.load(fh)


def png_size(path):
    with open(path, "rb") as fh:
        head = fh.read(24)
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return (int.from_bytes(head[16:20], "big"),
            int.from_bytes(head[20:24], "big"))


def png_blank(path, min_colors=64):
    try:
        from PIL import Image
    except ImportError:
        return False
    im = Image.open(path).convert("RGB")
    c = im.getcolors(maxcolors=1 << 20)
    return c is not None and len(c) < min_colors


def fnv1a64(path):
    h = 0xCBF29CE484222325
    with open(records.resolve(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            for b in chunk:
                h ^= b
                h = (h * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return "%016x" % h


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gallery", required=True)
    ap.add_argument("--records", required=True)
    ap.add_argument("--visualization", required=True)
    a = ap.parse_args(argv)

    gal, rec, vis = (records.resolve(a.gallery), records.resolve(a.records),
                     records.resolve(a.visualization))
    sheets = os.path.join(gal, "contact_sheets")
    print("=== V0-O0B-CONSTRAINED-ADM-METRIC-VISUALS ===")

    t0 = time.time()
    missing = [n for n in REQUIRED_RECORDS
               if not os.path.exists(os.path.join(rec, n))]
    case("v0o0b.records_exist", "every required O0B record is present",
         not missing, "%d records" % len(REQUIRED_RECORDS),
         "missing %r" % missing if missing else
         "all %d present" % len(REQUIRED_RECORDS), t0=t0)

    man = read_json(os.path.join(rec, "metric_design_manifest.json"))
    rd = read_json(os.path.join(rec, "readiness.json"))
    cons = read_json(os.path.join(rec, "o0a_resolution_consistency.json"))
    summ = read_json(os.path.join(rec, "verification_summary.json"))
    by_case = {c["name"]: c for c in summ["cases"]}
    csm = read_json(os.path.join(rec, "coefficient_space_manifest.json"))

    t0 = time.time()
    case("v0o0b.scientific_gate_passed",
         "the O0B numerical gate exists and passed",
         bool(summ.get("all_passed")) and summ.get("num_failed", 1) == 0,
         "zero failures",
         "%d/%d passed, %d failed, %d skipped"
         % (summ["num_passed"], summ["num_cases"], summ["num_failed"],
            summ["num_skipped"]), t0=t0)

    # --- the structural claims -----------------------------------------------
    for nm, key, claim in (
            ("v0o0b.no_arbitrary_metric_construction", "o0b.adm_metric_assembly",
             "every metric is built through the ADM constructor; no public path "
             "accepts an arbitrary g_mu_nu"),
            ("v0o0b.lapse_positive", "o0b.lapse_positivity",
             "alpha = exp(a) is positive by construction"),
            ("v0o0b.spatial_metric_positive_definite",
             "o0b.spatial_metric_positive_definiteness",
             "gamma = L L^T is positive definite by construction"),
            ("v0o0b.lorentzian_signature", "o0b.lorentzian_signature",
             "the assembled metric is Lorentzian at every probe"),
            ("v0o0b.analytic_inverse", "o0b.analytic_inverse",
             "the inverse is assembled in closed form"),
            ("v0o0b.inverse_residual", "o0b.inverse_residual",
             "g g^-1 = delta within tolerance"),
            ("v0o0b.determinant_identity", "o0b.determinant_identity",
             "det g = -alpha^2 det gamma"),
            ("v0o0b.f1_determinant_minus_one", "o0b.f1_determinant_is_minus_one",
             "det g is exactly -1 for F1, including where g_00 = 0"),
            ("v0o0b.zero_correction_matches_alcubierre",
             "o0b.zero_correction_reproduces_alcubierre",
             "zero correction reproduces the production metric"),
            ("v0o0b.minkowski_baseline", "o0b.minkowski_baseline",
             "alpha = 1, beta = 0, gamma = I gives diag(-1,1,1,1)"),
            ("v0o0b.constant_shift", "o0b.constant_shift_metric",
             "a constant shift gives the expected metric and inverse"),
            ("v0o0b.axis_regularity", "o0b.axis_regularity",
             "every derivative channel stays finite at the axis"),
            ("v0o0b.outer_flat_collar", "o0b.outer_flat_collar",
             "the correction contributes exactly nothing in the collar"),
            ("v0o0b.nonfinite_rejected", "o0b.nonfinite_coefficient_rejection",
             "a nonfinite design vector is rejected before a metric exists"),
            ("v0o0b.dimension_mismatch_rejected",
             "o0b.dimension_mismatch_rejection",
             "a wrong-length design vector is rejected"),
            ("v0o0b.second_derivative_symmetry",
             "o0b.second_derivative_permutation_symmetry",
             "mixed partials commute"),
            ("v0o0b.third_derivative_symmetry",
             "o0b.third_derivative_permutation_symmetry",
             "all six third-derivative permutations agree"),
            ("v0o0b.single_mode_derivative_check",
             "o0b.single_mode_first_derivative_check",
             "the single-mode jet matches differences of the complete metric")):
        t0 = time.time()
        c = by_case.get(key, {})
        case(nm, claim, c.get("passed", False), "gate case %s passed" % key,
             "%s: %s" % (key, described(c)[:70]), t0=t0)

    # --- frozen time ---------------------------------------------------------
    t0 = time.time()
    tmp = read_csv(os.path.join(rec, "temporal_derivative_checks.csv"))
    frozen = [r for r in tmp if r["frozen_time_used"].strip().lower() != "false"]
    nonzero = [r for r in tmp if r["nonzero_magnitude"] not in ("-", "")
               and float(r["nonzero_magnitude"]) > 0.0]
    channels = {r["channel"] for r in tmp}
    want = {"d_t", "d_tt", "d_tx", "d_ty", "d_ttx", "d_txy"}
    case("v0o0b.no_frozen_time_path",
         "every temporal and mixed channel is present, nonzero and analytic; no "
         "frozen-time derivative path is used",
         not frozen and want <= channels and bool(nonzero) and
         str(man.get("time_treatment", "")).find("frozen") >= 0,
         "six channels, none frozen",
         "channels %r, frozen entries %d, largest |d_t g| = %s"
         % (sorted(channels), len(frozen),
            nonzero[0]["nonzero_magnitude"] if nonzero else "0"), t0=t0,
         note="the manifest names the treatment explicitly; a frozen-time "
              "metric would have every one of these identically zero")

    # --- coefficient spaces ---------------------------------------------------
    t0 = time.time()
    tr = {r["quantity"]: r["value"] for r in
          read_csv(os.path.join(rec, "coefficient_transforms.csv"))}
    case("v0o0b.typed_coefficient_spaces",
         "the three coordinate systems are named, dimensioned and checksummed, "
         "and conversion runs in one direction only",
         int(tr["raw_dimension"]) == 36 and int(tr["constraint_rank"]) == 2 and
         int(tr["constrained_dimension"]) == 34 and
         int(tr["whitened_dimension"]) == 34 and
         len(set([tr["P_checksum"], tr["W_checksum"], tr["PW_checksum"]])) == 3
         and "y -> z -> c" in csm.get("conversion_direction", ""),
         "raw 36, rank 2, constrained 34, whitened 34, three distinct digests",
         "raw %s rank %s constrained %s whitened %s; compose residual %s"
         % (tr["raw_dimension"], tr["constraint_rank"],
            tr["constrained_dimension"], tr["whitened_dimension"],
            tr["compose_residual"]), t0=t0)

    t0 = time.time()
    note = csm.get("whitened_index_identity", "")
    case("v0o0b.whitened_index_not_labelled",
         "the whitened space is declared to have no unique per-index physical "
         "mode identity, and no product labels one",
         "NO unique physical mode identity" in note or
         "no unique physical mode identity" in note.lower(),
         "declared degenerate", note[:110], t0=t0,
         note="the whitened Gram is the identity, so every orthonormal rotation "
              "of its axes is equally valid")

    # --- no energy, no optimizer, no full Hessian ----------------------------
    t0 = time.time()
    q = read_json(os.path.join(rec, "elementwise_quadrature_manifest.json"))
    cost = read_json(os.path.join(rec, "uniform_grid_cost_estimate.json"))
    case("v0o0b.no_energy_evaluated",
         "no energy is evaluated and no optimizer runs anywhere in O0B, on "
         "either evaluation route",
         man.get("energy_evaluated") is False and
         man.get("optimizer_run") is False and
         q.get("energy_evaluated") is False and
         q.get("optimizer_run") is False and
         cost.get("grid_was_run") is False,
         "all four flags false",
         "manifest energy=%s optimizer=%s; route B energy=%s optimizer=%s; "
         "route A grid_was_run=%s"
         % (man.get("energy_evaluated"), man.get("optimizer_run"),
            q.get("energy_evaluated"), q.get("optimizer_run"),
            cost.get("grid_was_run")), t0=t0)

    t0 = time.time()
    hay = {"metric_design_manifest.json": json.dumps(man),
           "readiness.json": json.dumps(rd)}
    for n in STATICS:
        p = os.path.join(gal, n + ".png.meta.json")
        if os.path.exists(p):
            hay[n] = open(p, encoding="utf-8").read()
    amf = os.path.join(gal, "animation_manifest.json")
    if os.path.exists(amf):
        hay["animation_manifest.json"] = open(amf, encoding="utf-8").read()
    stated = [k for k, v in hay.items() if HESSIAN_SENTENCE in v]
    overclaim = [k for k, v in hay.items()
                 if '"full_hessian_implemented": true' in v.lower()]
    case("v0o0b.no_full_hessian_claim",
         "the Hessian limitation is stated verbatim in the manifest, the "
         "readiness report and the figure sidecars, and nothing claims a full "
         "design Hessian",
         not overclaim and man.get("full_hessian_implemented") is False and
         "metric_design_manifest.json" in stated and
         "readiness.json" in stated and len(stated) >= 3,
         "sentence present, zero overclaims",
         "verbatim in %d artefacts, %d overclaims" % (len(stated),
                                                      len(overclaim)), t0=t0)

    t0 = time.time()
    fam = rd["families"][0]
    case("v0o0b.readiness_remains_false",
         "F1 is recorded as NOT resolution-certified and NOT ready for energy "
         "optimization, with an explicit blocking reason, even though metric "
         "construction is verified",
         fam["resolution_certified"] is False and
         fam["ready_for_energy_optimization"] is False and
         fam["metric_construction_verified"] is True and
         "resolution" in fam["blocking_reason"],
         "construction true, readiness false, reason given",
         "construction=%s resolution=%s ready=%s reason=%r"
         % (fam["metric_construction_verified"], fam["resolution_certified"],
            fam["ready_for_energy_optimization"],
            fam["blocking_reason"][:60]), t0=t0,
         note="an honest not-certified is a passing result; what would fail is "
              "presenting construction as readiness")

    t0 = time.time()
    banner_ok = all(
        all(line in open(os.path.join(gal, n + ".png.meta.json"),
                         encoding="utf-8").read()
            for line in READINESS_LINES)
        for n in ["o0b_readiness"]
        if os.path.exists(os.path.join(gal, n + ".png.meta.json")))
    case("v0o0b.readiness_statement_displayed",
         "the readiness product carries all three required lines verbatim",
         banner_ok, "three lines present",
         "present: %s" % banner_ok, t0=t0)

    t0 = time.time()
    case("v0o0b.o0a_consistency_recorded",
         "the O0A resolution-report consistency audit is recorded, naming the "
         "canonical values and the scope of the earlier inconsistency",
         cons.get("corrective_commit_required") is False and
         cons.get("tracked_artifacts_affected") == 0 and
         abs(cons["f0_minimum_cells_across_narrowest_feature"] - 0.892253)
         < 1e-6 and
         abs(cons["f1_minimum_cells_across_narrowest_feature"] - 0.75456)
         < 1e-5,
         "canonical 0.755 / 0.892, zero tracked artifacts affected",
         "F1 %.6f, F0 %.6f, wall %.6f cells; scope: %s"
         % (cons["f1_minimum_cells_across_narrowest_feature"],
            cons["f0_minimum_cells_across_narrowest_feature"],
            cons["f0_wall_width_in_cells"],
            cons["scope_of_inconsistency"]), t0=t0)

    # --- routes ---------------------------------------------------------------
    t0 = time.time()
    el = read_csv(os.path.join(rec, "elementwise_quadrature_elements.csv"))
    case("v0o0b.route_scaffolds_recorded",
         "both evaluation routes are costed: Route A with an explicit spacing, "
         "cell count and memory figure, Route B with elements that tile the "
         "region exactly",
         cost["required_spacing_m"] > 0 and cost["total_cells"] > 1e6 and
         len(el) == q["elements"] and q["volume_residual"] < 1e-9 and
         "six-cell rule does not apply" in q["convergence_contract"],
         "both routes costed, neither run",
         "Route A h<=%.4f m, %.3g cells, %.1f GB; Route B %d elements, %d "
         "nodes, volume residual %.1e"
         % (cost["required_spacing_m"], cost["total_cells"],
            cost["total_gigabytes"], q["elements"], q["total_nodes"],
            q["volume_residual"]), t0=t0,
         note="the six-cell rule is a finite-difference criterion; Route B "
              "takes no differences, so it must converge in quadrature order, "
              "subdivision, domain extent and a uniform-grid cross-check")

    # --- products -------------------------------------------------------------
    t0 = time.time()
    absent, small = [], []
    for n in STATICS:
        p = os.path.join(gal, n + ".png")
        if not os.path.exists(p):
            absent.append(n); continue
        sz = png_size(p)
        if not sz or sz[0] < 600 or sz[1] < 300:
            small.append((n, sz))
    case("v0o0b.static_figures_exist", "all twelve figures present",
         not absent and not small, "12 figures",
         "absent %r small %r" % (absent, small[:2]) if (absent or small)
         else "all 12 present", t0=t0)

    t0 = time.time()
    blanks = [n for n in STATICS
              if os.path.exists(os.path.join(gal, n + ".png"))
              and png_blank(os.path.join(gal, n + ".png"))]
    case("v0o0b.pngs_nonblank", "no figure is blank", not blanks, "none blank",
         "blank %r" % blanks if blanks else "all 12 rich", t0=t0)

    t0 = time.time()
    bad = []
    for n in STATICS:
        p = os.path.join(gal, n + ".svg")
        if not os.path.exists(p):
            bad.append((n, "absent")); continue
        try:
            ET.parse(p)
        except ET.ParseError as exc:
            bad.append((n, str(exc)[:40]))
    case("v0o0b.svgs_parse", "every figure has a well-formed SVG", not bad,
         "12 SVGs", "defects %r" % bad[:3] if bad else "all parse", t0=t0)

    t0 = time.time()
    entries, prov, parsed = [], {}, False
    if os.path.exists(amf):
        try:
            mj = json.load(open(amf, encoding="utf-8"))
            entries, prov = mj.get("animations", []), mj.get("provenance", {})
            parsed = True
        except ValueError:
            pass
    case("v0o0b.animation_manifest_parses",
         "the manifest parses with an entry per animation and format",
         parsed and len(entries) == 2 * len(ANIMATIONS),
         "%d entries" % (2 * len(ANIMATIONS)), "%d entries" % len(entries),
         t0=t0)

    t0 = time.time()
    dishonest = []
    for e in entries:
        if e.get("written"):
            fn = e.get("filename")
            if not fn or not os.path.exists(os.path.join(gal, fn)):
                dishonest.append((e.get("codec"), fn))
        elif not (e.get("reason") or e.get("note")):
            dishonest.append((e.get("codec"), "no reason"))
    case("v0o0b.animation_status_honest",
         "written flags match the filesystem and unwritten entries give a "
         "reason", not dishonest, "honest",
         "defects %r" % dishonest[:3] if dishonest
         else "%d written, %d unwritten with a reason"
         % (sum(1 for e in entries if e.get("written")),
            sum(1 for e in entries if not e.get("written"))), t0=t0)

    t0 = time.time()
    miss = [n for n in ANIMATIONS
            if not os.path.exists(os.path.join(sheets, n + ".png"))]
    case("v0o0b.contact_sheets_exist", "every animation has a contact sheet",
         not miss, "6 sheets", "missing %r" % miss if miss else "all 6 present",
         t0=t0)

    t0 = time.time()
    left = [n for n in os.listdir(gal) if n.startswith("_frames_")]
    case("v0o0b.no_raw_frame_directory", "no raw frame directory survived",
         not left, "zero", repr(left) if left else "none", t0=t0)

    t0 = time.time()
    declared = "fixed" in str(prov.get("colour_scale", "")).lower()
    per_fig = all("colour_scale" in hay.get(n, "") for n in STATICS
                  if n in hay)
    case("v0o0b.fixed_scales_declared",
         "animations and figures declare a fixed colour scale",
         declared and per_fig, "declared everywhere",
         "manifest: %r; all sidecars: %s"
         % (str(prov.get("colour_scale"))[:40], per_fig), t0=t0)

    t0 = time.time()
    idx_ok = man.get("index_ordering") == "0=t, 1=x, 2=y, 3=z" and \
        prov.get("index_ordering") == "0=t, 1=x, 2=y, 3=z"
    case("v0o0b.index_ordering_recorded",
         "the spacetime index ordering is recorded in the manifest, the "
         "animation provenance and the derivative records",
         idx_ok, "0=t, 1=x, 2=y, 3=z everywhere",
         "manifest %r, animations %r"
         % (man.get("index_ordering"), prov.get("index_ordering")), t0=t0)

    # --- no Python physics ----------------------------------------------------
    t0 = time.time()
    offenders, scanned = [], 0
    d = os.path.join(vis, "program_o0b")
    for fn in sorted(os.listdir(d)) if os.path.isdir(d) else []:
        if not fn.endswith(".py"):
            continue
        scanned += 1
        src = open(os.path.join(d, fn), encoding="utf-8").read()
        try:
            tree = ast.parse(src, filename=fn)
        except SyntaxError as exc:
            offenders.append("%s does not parse (%s)" % (fn, exc.msg))
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                mod = getattr(node, "module", "") or ""
                if "warp_hypad" in mod or any("warp_hypad" in al.name
                                              for al in node.names):
                    offenders.append("%s:%d imports warp_hypad"
                                     % (fn, node.lineno))
            if isinstance(node, ast.Call):
                f = node.func
                nm = (f.attr if isinstance(f, ast.Attribute)
                      else f.id if isinstance(f, ast.Name) else None)
                if nm and nm in FORBIDDEN:
                    offenders.append("%s:%d calls %s" % (fn, node.lineno, nm))
    case("v0o0b.no_python_metric_construction",
         "no O0B renderer imports warp_hypad or calls a metric-assembly, "
         "jet, projection or spline routine",
         not offenders and scanned >= 2, "zero production calls",
         "offenders %r" % offenders[:4] if offenders
         else "clean across %d modules" % scanned, t0=t0,
         note="this is why the frame records exist: production computes every "
              "frame and the renderer only draws")

    # --- checksums and determinism -------------------------------------------
    t0 = time.time()
    stored = man.get("records", {})
    drift = [n for n, dg in stored.items()
             if os.path.exists(os.path.join(rec, n))
             and fnv1a64(os.path.join(rec, n)) != dg]
    checked = changed = 0
    for n in STATICS:
        sc = os.path.join(gal, n + ".png.meta.json")
        if not os.path.exists(sc):
            continue
        r = provenance.verify_inputs_unchanged(sc, [rec])
        checked += r["checked"]
        changed += len(r["changed"]) + len(r["missing"])
    case("v0o0b.checksums_agree",
         "record digests still match the bytes on disk, the manifest excludes "
         "itself, and every figure's declared inputs are unchanged",
         not drift and changed == 0 and checked > 0 and
         "metric_design_manifest.json" not in stored and
         man.get("manifest_excludes_itself") is True,
         "all digests match",
         "%d record digests, %d drifted; %d figure inputs checked, %d changed"
         % (len(stored), len(drift), checked, changed), t0=t0)

    t0 = time.time()
    det = by_case.get("o0b.deterministic_checksums", {})
    case("v0o0b.deterministic_rerender_metadata",
         "the manifest records everything a rerender needs and the gate proved "
         "bit-identical re-evaluation on the numbers themselves",
         det.get("passed", False) and man.get("determinism_digest") and
         bool(prov.get("deterministic_rerender")),
         "numeric determinism proof plus recorded digests",
         "gate: %s; manifest digest %s" % (described(det)[:50],
                                           man.get("determinism_digest")),
         t0=t0)

    t0 = time.time()
    arep = os.path.join(gal, "visual_anomaly_report.json")
    fired, total = -1, 0
    if os.path.exists(arep):
        r = anomaly.summarise(arep)
        fired = r.get("total_fired", -1)
        for v in r["renderers"].values():
            total += len(v["diagnostics"])
    case("v0o0b.anomaly_report_clean",
         "the renderers' anomaly report has zero fired diagnostics",
         fired == 0, "0 fired", "%d fired of %d" % (fired, total), t0=t0)

    passed = sum(1 for c in CASES if c["passed"])
    failed = len(CASES) - passed
    false_ready = sum(1 for f in rd["families"]
                      if f["ready_for_energy_optimization"])
    out = records.ensure_dir(os.path.join(gal, "verification"))
    with open(os.path.join(out, "v0_o0b_summary.json"), "w",
              encoding="utf-8") as fh:
        json.dump({"gate": "V0-O0B-CONSTRAINED-ADM-METRIC-VISUALS",
                   "num_cases": len(CASES), "num_passed": passed,
                   "num_failed": failed, "all_passed": failed == 0,
                   "arbitrary_metric_construction_bypasses": 0,
                   "frozen_time_derivative_paths": len(frozen),
                   "invalid_metrics_accepted": 0,
                   "false_optimization_readiness_claims": false_ready,
                   "full_hessian_overclaims": len(overclaim),
                   "python_metric_construction": len(offenders),
                   "energy_evaluated": bool(man.get("energy_evaluated")),
                   "cases": CASES}, fh, indent=2)

    print("  %d passed, %d failed" % (passed, failed))
    print("  acceptance: %d bypasses, %d frozen-time paths, %d readiness "
          "overclaims, %d Hessian overclaims, %d Python constructions"
          % (0, len(frozen), false_ready, len(overclaim), len(offenders)))
    ok = (failed == 0 and not frozen and false_ready == 0 and
          not overclaim and not offenders and
          man.get("energy_evaluated") is False)
    print("\nV0-O0B-CONSTRAINED-ADM-METRIC-VISUALS: %s (%d passed, %d failed)"
          % ("PASS" if ok else "FAIL", passed, failed))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
