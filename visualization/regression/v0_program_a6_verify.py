#!/usr/bin/env python3
"""Gate V0-PROGRAM-A6-ENERGY-VISUALS.

    python visualization/regression/v0_program_a6_verify.py \
        --gallery  <build>/results/visualization/program_a6 \
        --sweep    <build>/results/verification/v0_program_a6_energy \
        --results  <build>/results/verification/milestone_06 \
        --visualization visualization

Verifies the A6 observer-energy visual products against their production
sources. Three outcomes, kept distinct because "missing" and "blocked by
unimplemented science" are different statements:

    PASS
    FAIL
    BLOCKED_BY_UNIMPLEMENTED_SCIENCE

Acceptance: 0 failures, 0 invalid observers presented as physical, 0
uncertified signs presented as definite, 0 exact/smooth objective conflations.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import sys
import time
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import anomaly, provenance, records                   # noqa: E402

CASES = []

STATICS = [
    "eulerian_energy_density", "momentum_density", "principal_stresses",
    "minimum_nec", "worst_null_direction_field", "energy_condition_status",
    "observer_comparison", "integrated_energy_vs_velocity",
    "integrated_nec_vs_velocity", "resolution_certification",
    "domain_convergence", "null_direction_sphere",
]

# Each figure named against the production record it must be traceable to.
FIGURE_SOURCES = {
    "eulerian_energy_density": "energy_frame_003.csv",
    "momentum_density": "energy_frame_003.csv",
    "principal_stresses": "energy_frame_003.csv",
    "minimum_nec": "energy_frame_003.csv",
    "worst_null_direction_field": "worst_null_directions.csv",
    "energy_condition_status": "energy_frame_003.csv",
    "observer_comparison": "observer_comparison.csv",
    "integrated_energy_vs_velocity": "integrated_objectives.csv",
    "integrated_nec_vs_velocity": "integrated_objectives.csv",
    "resolution_certification": "resolution_certification.csv",
    "domain_convergence": "domain_convergence.csv",
    "null_direction_sphere": "null_direction_sphere.csv",
}

ANIMATIONS = [
    "energy_density_vs_velocity_scientific",
    "energy_density_vs_velocity_cinematic",
    "nec_violation_vs_velocity_scientific",
    "nec_violation_vs_velocity_cinematic",
    "worst_null_direction_vs_velocity",
    "observer_energy_comparison",
    "integrated_energy_objectives",
    "energy_condition_certification",
]

# Products the A6 science does not support. Reported as BLOCKED, never as a
# pass and never silently omitted. The dependency is exact, so the label can be
# retired the moment the science lands -- a stale blocked label reads as an
# honest limitation long after the limitation is gone.
BLOCKED = {
    "energy_vs_radius":
        "an A6 production path that varies the Alcubierre radius R and "
        "re-evaluates the frame stress-energy, the certification and the "
        "integrated objectives. A6 fixes R = 2 (kR in observer_checks.cpp and "
        "energy_integral_checks.cpp).",
    "energy_vs_steepness":
        "an A6 production path that varies the shape-function steepness sigma "
        "and re-evaluates the same quantities. A6 fixes sigma = 1 (kSig).",
    "stress_eigenvector_field":
        "an exported eigenbasis of S_ij. geometry::detail::"
        "orthonormal_eigenbasis computes one but it is internal to the NEC "
        "minimiser and no record carries it.",
}

# The 21 anomaly diagnostics the A6 renderers are required to run. Named here
# so a renderer cannot quietly drop one: a check that stops running looks
# exactly like a check that keeps passing.
REQUIRED_DIAGNOSTICS = {
    "nonfinite_value",
    "eulerian_normalization_failure",
    "tetrad_orthonormality_failure",
    "tetrad_handedness_change",
    "tetrad_sign_flip",
    "invalid_stationary_observer_rendered_valid",
    "observer_identity_changed_between_frames",
    "missed_off_axis_negative_energy",
    "nec_direction_does_not_reproduce_n_min",
    "nec_direction_jump_without_degeneracy",
    "status_inconsistent_with_certification_margin",
    "uncertified_value_rendered_definite",
    "integrated_objective_closure_failure",
    "boundary_dominated_integral",
    "insufficient_cells_across_wall",
    "domain_changed_between_frames",
    "per_frame_color_rescaling",
    "missing_or_mismatched_source_rows",
    "source_checksum_mismatch",
    "python_physics_recomputation",
    "blank_frame",
}

CLOSURE_EPSILON = 1e-9


def case(name, claim, ok, expected, computed, note="", t0=None, state=None):
    st = state or ("PASS" if ok else "FAIL")
    CASES.append({"name": name, "claim": claim, "state": st,
                  "passed": st == "PASS", "expected": str(expected),
                  "computed": str(computed), "note": note,
                  "seconds": round(time.time() - t0, 4) if t0 else 0.0})
    print("  %-38s %-46s %s" % (st, name, computed))
    return st == "PASS"


def png_size(path):
    with open(path, "rb") as fh:
        head = fh.read(24)
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return (int.from_bytes(head[16:20], "big"),
            int.from_bytes(head[20:24], "big"))


def png_is_blank(path, min_colors=64):
    try:
        from PIL import Image
    except ImportError:
        return False
    im = Image.open(path).convert("RGB")
    cols = im.getcolors(maxcolors=1 << 20)
    return cols is not None and len(cols) < min_colors


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gallery", required=True)
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--visualization", required=True)
    a = ap.parse_args(argv)

    gal = records.resolve(a.gallery)
    sw = records.resolve(a.sweep)
    res = records.resolve(a.results)
    vis = records.resolve(a.visualization)
    sheets = os.path.join(gal, "contact_sheets")

    print("=== V0-PROGRAM-A6-ENERGY-VISUALS ===")

    # --- 1. The A6 scientific gate exists and passed -------------------------
    #
    # No A6 visual is accepted unless the science under it has its own passing
    # gate. An energy density is the first quantity in this project that can be
    # believed on sight, which is exactly why the picture must not be allowed to
    # outrun the verification.
    t0 = time.time()
    summ = os.path.join(res, "verification_summary.json")
    ok, detail = False, "summary absent at %s" % summ
    if os.path.exists(summ):
        with open(summ, encoding="utf-8") as fh:
            j = json.load(fh)
        ok = bool(j.get("all_passed")) and j.get("num_failed", 1) == 0
        detail = "%d/%d passed, %d failed, %d skipped" % (
            j.get("num_passed", 0), j.get("num_cases", 0),
            j.get("num_failed", -1), j.get("num_skipped", -1))
    case("v0a6.scientific_gate_passed",
         "the A6 numerical gate exists and passed, so the observers, the "
         "energy density and the certification behind these visuals are "
         "verified rather than merely computed",
         ok, "A6 gate all_passed with zero failures", detail, t0=t0)

    # --- 2. The capability audit ---------------------------------------------
    t0 = time.time()
    apath = os.path.join(vis, "program_a6", "a6_capability_audit.json")
    audit = {}
    if os.path.exists(apath):
        with open(apath, encoding="utf-8") as fh:
            audit = json.load(fh)
    vocab = set(audit.get("classification_vocabulary", []))
    caps = audit.get("capabilities", [])
    bad = [c.get("quantity") for c in caps
           if c.get("classification") not in vocab]
    unimpl = {c["quantity"] for c in caps
              if c.get("classification") == "not implemented"}
    blocked_named = {b["product"] for b in audit.get("blocked_products", [])}
    case("v0a6.capability_audit_present",
         "the capability audit exists, classifies every quantity with the "
         "closed four-value vocabulary, and names a blocked product for every "
         "capability it calls unimplemented",
         bool(caps) and not bad and bool(unimpl) and bool(blocked_named),
         "%d capabilities, every classification in the vocabulary" % len(caps),
         "%d capabilities, %d unclassifiable, %d not implemented, %d blocked "
         "products named" % (len(caps), len(bad), len(unimpl),
                             len(blocked_named)), t0=t0,
         note="a quantity absent from production must be visible as absent, "
              "not inferred from a neighbouring column")

    # --- 3. Every produced visual has a production source --------------------
    t0 = time.time()
    missing_src = []
    for n, src in FIGURE_SOURCES.items():
        if not os.path.exists(os.path.join(gal, n + ".png")):
            continue
        if not os.path.exists(os.path.join(sw, src)):
            missing_src.append((n, src))
    case("v0a6.every_visual_has_a_production_source",
         "every produced figure names a saved A6 record it was drawn from, and "
         "that record exists",
         not missing_src, "%d figure-to-record links" % len(FIGURE_SOURCES),
         "unlinked %r" % missing_src[:3] if missing_src
         else "all %d figures trace to an exported CSV" % len(FIGURE_SOURCES),
         t0=t0)

    # --- 4. Unsupported products are explicitly blocked ----------------------
    t0 = time.time()
    leaked = [n for n in BLOCKED
              if os.path.exists(os.path.join(gal, n + ".png"))
              or os.path.exists(os.path.join(gal, n + ".gif"))]
    undeclared = sorted(set(BLOCKED) - blocked_named)
    case("v0a6.unsupported_products_blocked",
         "every product the A6 science cannot support is declared "
         "BLOCKED_BY_UNIMPLEMENTED_SCIENCE with its exact dependency, and none "
         "of them was drawn anyway",
         not leaked and not undeclared, "%d blocked products" % len(BLOCKED),
         ("drawn anyway %r" % leaked) if leaked else
         ("not in the audit %r" % undeclared) if undeclared else
         "%d blocked, all declared, none drawn: %s"
         % (len(BLOCKED), ", ".join(sorted(BLOCKED))), t0=t0,
         state="BLOCKED_BY_UNIMPLEMENTED_SCIENCE"
         if (not leaked and not undeclared) else None)

    # --- 5-6. Static figures exist and are not blank -------------------------
    t0 = time.time()
    absent, small = [], []
    for n in STATICS:
        p = os.path.join(gal, n + ".png")
        if not os.path.exists(p):
            absent.append(n)
            continue
        sz = png_size(p)
        if not sz or sz[0] < 600 or sz[1] < 300:
            small.append((n, sz))
    case("v0a6.static_figures_exist",
         "every required A6 figure is present as a PNG with sane dimensions",
         not absent and not small, "%d figures" % len(STATICS),
         "absent %r small %r" % (absent, small[:2]) if (absent or small)
         else "all %d present" % len(STATICS), t0=t0)

    t0 = time.time()
    blanks = [n for n in STATICS
              if os.path.exists(os.path.join(gal, n + ".png"))
              and png_is_blank(os.path.join(gal, n + ".png"))]
    case("v0a6.pngs_nonblank", "no A6 figure is blank", not blanks,
         "every figure carries real content",
         "blank %r" % blanks if blanks else "all %d rich" % len(STATICS),
         t0=t0)

    # --- 7. SVGs parse --------------------------------------------------------
    t0 = time.time()
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
    case("v0a6.svgs_parse", "every figure has a well-formed SVG companion",
         not svg_bad, "%d SVGs" % len(STATICS),
         "defects %r" % svg_bad[:3] if svg_bad else "all parse", t0=t0)

    # --- 8-9. The animation manifest, and honest GIF/MP4 status --------------
    t0 = time.time()
    mpath = os.path.join(gal, "animation_manifest.json")
    entries, prov = [], {}
    parsed = False
    if os.path.exists(mpath):
        try:
            with open(mpath, encoding="utf-8") as fh:
                mj = json.load(fh)
            entries, prov = mj.get("animations", []), mj.get("provenance", {})
            parsed = True
        except ValueError as exc:
            entries, prov = [], {"parse_error": str(exc)}
    have_gif = [n for n in ANIMATIONS
                if os.path.exists(os.path.join(gal, n + ".gif"))]
    case("v0a6.animation_manifest_parses",
         "the animation manifest parses and carries an entry for every "
         "animation and format",
         parsed and len(entries) == 2 * len(ANIMATIONS)
         and len(have_gif) == len(ANIMATIONS),
         "%d entries for %d animations" % (2 * len(ANIMATIONS),
                                           len(ANIMATIONS)),
         "%d entries, %d GIFs present" % (len(entries), len(have_gif)), t0=t0)

    t0 = time.time()
    dishonest = []
    for e in entries:
        written = bool(e.get("written"))
        fn = e.get("filename")
        if written:
            if not fn or not os.path.exists(os.path.join(gal, fn)):
                dishonest.append((e.get("codec"), fn, "claimed written, absent"))
        else:
            if not (e.get("reason") or e.get("note")):
                dishonest.append((e.get("codec"), fn, "unwritten, no reason"))
    case("v0a6.animation_status_is_honest",
         "every manifest entry that claims a file exists has one, and every "
         "unwritten entry carries a reason",
         not dishonest, "written flags match the filesystem",
         "defects %r" % dishonest[:3] if dishonest
         else "%d entries honest; %d written, %d unwritten with a reason"
         % (len(entries), sum(1 for e in entries if e.get("written")),
            sum(1 for e in entries if not e.get("written"))), t0=t0,
         note="a missing video encoder must not fail this gate: the GIFs, the "
              "frames and the contact sheets are the deterministic products")

    # --- 10. Contact sheets ---------------------------------------------------
    t0 = time.time()
    miss = [n for n in ANIMATIONS
            if not os.path.exists(os.path.join(sheets, n + ".png"))]
    sheet_blank = [n for n in ANIMATIONS
                   if os.path.exists(os.path.join(sheets, n + ".png"))
                   and png_is_blank(os.path.join(sheets, n + ".png"))]
    case("v0a6.contact_sheets_exist",
         "each animation has a nonblank contact sheet, so a regression is "
         "visible without opening a video",
         not miss and not sheet_blank, "%d sheets" % len(ANIMATIONS),
         "missing %r blank %r" % (miss, sheet_blank)
         if (miss or sheet_blank) else "all %d present and rich"
         % len(ANIMATIONS), t0=t0)

    # --- 11. Frame counts agree ----------------------------------------------
    t0 = time.time()
    n_frames = len([n for n in os.listdir(sw)
                    if n.startswith("energy_frame_") and n.endswith(".csv")])
    n_sweep = len(records.read_csv(os.path.join(sw,
                                                "energy_velocity_sweep.csv")))
    n_integ = len(records.read_csv(os.path.join(sw,
                                                "integrated_objectives.csv")))
    manifest_counts = sorted({e.get("frame_count") for e in entries})
    agree = (n_frames == n_sweep == n_integ
             and manifest_counts == [n_frames]
             and prov.get("frame_count") == n_frames)
    case("v0a6.frame_counts_agree",
         "the frame records, the summary rows, the manifest entries and the "
         "manifest provenance all report the same frame count",
         agree, "%d everywhere" % n_frames,
         "frames %d, sweep %d, integrated %d, manifest %r, provenance %s"
         % (n_frames, n_sweep, n_integ, manifest_counts,
            prov.get("frame_count")), t0=t0)

    # --- 12. Fixed global scales declared ------------------------------------
    t0 = time.time()
    scale = str(prov.get("colour_scale", ""))
    limits = prov.get("colour_limits") or {}
    gs = os.path.join(gal, "global_scales.json")
    declared = ("fixed" in scale.lower() and bool(limits)
                and os.path.exists(gs))
    case("v0a6.fixed_global_scales_declared",
         "the animations declare a colour scale computed once from the whole "
         "sweep and frozen, and record the limits themselves",
         declared, "fixed scale declared with explicit limits",
         "%s | limits %r" % (scale[:60], sorted(limits)), t0=t0,
         note="a per-frame autoscale makes a field that grew by an order of "
              "magnitude look constant while only its noise appears to move")

    # --- 13-14. Observers stay distinct, and invalid ones are not physical ---
    t0 = time.time()
    orows = records.read_csv(os.path.join(sw, "observer_comparison.csv"))
    families = sorted({r["observer_type"] for r in orows})
    per_frame = {}
    for r in orows:
        per_frame.setdefault(r["frame"], set()).add(r["observer_type"])
    stable = len({frozenset(v) for v in per_frame.values()}) == 1
    case("v0a6.observers_remain_distinct",
         "all four observer families appear as separate rows in every frame "
         "and are never merged",
         len(families) == 4 and stable, "4 families in every frame",
         "families %r, %d distinct per-frame sets"
         % (families, len({frozenset(v) for v in per_frame.values()})), t0=t0)

    t0 = time.time()
    presented = []
    for r in orows:
        valid = r["observer_valid"] == "1"
        measured = r["rho_measured"].strip()
        g00 = float(r["g00"])
        if r["observer_type"] == "coordinate_stationary":
            if valid and g00 >= 0.0:
                presented.append((r["frame"], r["site"], "valid at g00 >= 0"))
            if not valid and measured:
                presented.append((r["frame"], r["site"],
                                  "rejected yet carries a value"))
        if not valid and measured:
            presented.append((r["frame"], r["observer_type"],
                              "inadmissible yet carries a value"))
        if valid and not measured:
            presented.append((r["frame"], r["observer_type"],
                              "admissible yet carries no value"))
    n_rejected = sum(1 for r in orows if r["observer_valid"] != "1")
    case("v0a6.no_invalid_observer_presented_as_physical",
         "no inadmissible observer carries a measured energy density, and the "
         "coordinate-stationary observer is never marked valid where "
         "g_00 >= 0",
         not presented, "zero invalid observers presented as physical",
         "violations %r" % presented[:3] if presented
         else "%d rejected observer rows, none of them plotted" % n_rejected,
         t0=t0,
         note="normalising by |g_00| there would produce a plausible, "
              "plottable, entirely fictitious energy density")

    # --- 15. Exact and smooth objectives stay distinct -----------------------
    t0 = time.time()
    integ = records.read_csv(os.path.join(sw, "integrated_objectives.csv"))
    same = [r["v_s"] for r in integ
            if abs(float(r["exact_value"]) -
                   float(r["smooth_surrogate_value"])) <= 1e-15]
    net = max(float(r["net_closure_residual"]) for r in integ)
    ab = max(float(r["abs_closure_residual"]) for r in integ)
    widths = sorted({r["softplus_width"] for r in integ})
    case("v0a6.exact_and_smooth_objectives_distinct",
         "the exact objective and the softplus surrogate are separate columns "
         "with different values, and the surrogate width is recorded",
         not same and len(widths) == 1, "distinct at every velocity",
         "identical at %r; softplus width %r" % (same[:3], widths), t0=t0,
         note="softplus systematically overestimates max(0,x) near zero, which "
              "is exactly where certification happens")

    t0 = time.time()
    case("v0a6.integrated_objectives_close",
         "the integrated objectives satisfy |E_net-(E_+-E_-)| < eps and "
         "|E_abs-(E_++E_-)| < eps on the saved values",
         max(net, ab) < CLOSURE_EPSILON, "< %g" % CLOSURE_EPSILON,
         "worst net %.3e, worst abs %.3e" % (net, ab), t0=t0)

    t0 = time.time()
    bf = max(float(r["boundary_fraction"]) for r in integ)
    tails = [(float(r["v_s"]), float(r["tail_estimate"]),
              float(r["E_minus"])) for r in integ]
    not_ready = [v for v, t, e in tails if t >= 0.25 * e or bf >= 0.35]
    ready_flags = sorted({r["optimization_ready"] for r in integ})
    case("v0a6.optimization_readiness_is_earned",
         "no velocity is marked optimization-ready while its integral is "
         "boundary-dominated or its exterior tail is comparable to its value",
         not not_ready or ready_flags == ["0"],
         "boundary share < 0.35 and tail < 25% of E_-",
         "worst boundary share %.4f, %d velocities failing the tail test, "
         "readiness flags %r" % (bf, len(not_ready), ready_flags), t0=t0)

    # --- 16. Certification status is displayed -------------------------------
    t0 = time.time()
    missing_cert = []
    for n in STATICS:
        sc = os.path.join(gal, n + ".png.meta.json")
        if not os.path.exists(sc):
            missing_cert.append((n, "no sidecar"))
            continue
        with open(sc, encoding="utf-8") as fh:
            p = json.load(fh)
        if p.get("certification_status_displayed") is not True:
            missing_cert.append((n, "not declared"))
        if p.get("recomputes_physics") is not False:
            missing_cert.append((n, "recomputes_physics not false"))
    case("v0a6.certification_status_displayed",
         "every figure's sidecar declares that certification status is shown "
         "and that the renderer computed no physics",
         not missing_cert, "%d sidecars" % len(STATICS),
         "defects %r" % missing_cert[:3] if missing_cert
         else "all %d declare it" % len(STATICS), t0=t0)

    # --- 17. Source checksums agree ------------------------------------------
    t0 = time.time()
    checked = changed = 0
    missing_meta = []
    for n in STATICS:
        sc = os.path.join(gal, n + ".png.meta.json")
        if not os.path.exists(sc):
            missing_meta.append(n)
            continue
        r = provenance.verify_inputs_unchanged(sc, [sw])
        checked += r["checked"]
        changed += len(r["changed"]) + len(r["missing"])
    recorded = {}
    cpath = os.path.join(gal, "source_checksums.json")
    if os.path.exists(cpath):
        with open(cpath, encoding="utf-8") as fh:
            recorded = json.load(fh)
    drifted = [k for k, h in recorded.items()
               if os.path.exists(os.path.join(sw, k))
               and sha256(os.path.join(sw, k)) != h]
    case("v0a6.source_checksums_agree",
         "every figure sidecar names its input records and their SHA-256 still "
         "matches the bytes on disk",
         not missing_meta and changed == 0 and checked > 0 and not drifted,
         "all inputs unchanged",
         "%d input checksums checked, %d mismatched, %d sidecars missing, "
         "%d of %d recorded record hashes drifted"
         % (checked, changed, len(missing_meta), len(drifted), len(recorded)),
         t0=t0,
         note="a figure whose input bytes have changed is stale, and stale is "
              "indistinguishable from correct by eye")

    # --- 18. The renderer imports no A6 physics ------------------------------
    #
    # AST, not regex: a regex cannot tell RECOMPUTING an energy density from
    # READING a saved column called rho, and it trips on the docstrings that
    # explain the rule.
    t0 = time.time()
    FORBIDDEN = {
        "nec_value", "nec_minimum", "project_to_frame", "adm_decompose",
        "build_tetrad", "make_eulerian", "make_coordinate_stationary",
        "make_bubble_comoving", "make_tetrad_boosted", "certify", "wec", "dec",
        "sec", "symmetric3_eigenvalues", "fibonacci_direction",
        "shape_function", "evaluate", "cosh", "sinh", "tanh", "atanh",
    }
    offenders = []
    a6dir = os.path.join(vis, "program_a6")
    for fn in sorted(os.listdir(a6dir)) if os.path.isdir(a6dir) else []:
        if not fn.endswith(".py"):
            continue
        p = os.path.join(a6dir, fn)
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
                if nm and nm in FORBIDDEN:
                    offenders.append("%s:%d calls %s" % (fn, node.lineno, nm))
    case("v0a6.renderer_imports_no_a6_physics",
         "no A6 renderer imports warp_hypad or calls an observer, tetrad, NEC "
         "or certification routine",
         not offenders, "zero physics in the A6 renderer package",
         "offenders %r" % offenders[:4] if offenders
         else "clean across %d modules"
         % len([f for f in os.listdir(a6dir) if f.endswith(".py")]), t0=t0,
         note="this is the case that keeps the visualization contract true as "
              "the package grows")

    # --- 19. The anomaly report ----------------------------------------------
    t0 = time.time()
    arep = os.path.join(gal, "visual_anomaly_report.json")
    fired, total, ran = -1, 0, set()
    if os.path.exists(arep):
        rep = anomaly.summarise(arep)
        fired = rep.get("total_fired", -1)
        for r in rep["renderers"].values():
            total += len(r["diagnostics"])
            ran |= {d["diagnostic"] for d in r["diagnostics"]}
    absent_diags = sorted(REQUIRED_DIAGNOSTICS - ran)
    case("v0a6.anomaly_report_clean",
         "the anomaly report exists, ran every required diagnostic, and has "
         "zero fired cases",
         fired == 0 and not absent_diags,
         "0 fired, %d required diagnostics all run"
         % len(REQUIRED_DIAGNOSTICS),
         "%d fired of %d diagnostics; not run: %r"
         % (fired, total, absent_diags), t0=t0,
         note="a fired anomaly does not suppress the figure; it fails this "
              "gate while the evidence survives")

    # --- 20. Deterministic rerender metadata ---------------------------------
    t0 = time.time()
    with open(os.path.join(sw, "a6_energy_manifest.json"),
              encoding="utf-8") as fh:
        smeta = json.load(fh)
    needed = ("git_commit", "velocity_sequence", "exact_onset_frame",
              "first_spatial_contour_frame", "grid_resolution",
              "domain_extent", "radius", "sigma", "reference_velocity",
              "numerical_error_envelope_floor", "certification_safety_factor",
              "softplus_width", "nec_quadrature_points", "sphere_sample_count")
    absent_meta = [k for k in needed if smeta.get(k) is None]
    onset = smeta.get("exact_onset_frame")
    contour = smeta.get("first_spatial_contour_frame")
    seq_ok = (onset == 8 and contour == 9
              and prov.get("exact_onset_frame") == 8
              and prov.get("first_spatial_contour_frame") == 9)
    case("v0a6.deterministic_rerender_metadata",
         "the sweep manifest records everything a rerender needs, and the "
         "established sequence is preserved: exact_onset_frame 8, "
         "first_spatial_contour_frame 9",
         not absent_meta and seq_ok,
         "complete metadata, onset 8, contour 9",
         "missing %r; onset %s, contour %s" % (absent_meta, onset, contour),
         t0=t0,
         note="at onset the field only touches zero, so no contour is "
              "supported until the next frame -- the V0.1A convention, "
              "unchanged")

    # --- Summary --------------------------------------------------------------
    passed = sum(1 for c in CASES if c["state"] == "PASS")
    failed = sum(1 for c in CASES if c["state"] == "FAIL")
    blocked = sum(1 for c in CASES
                  if c["state"] == "BLOCKED_BY_UNIMPLEMENTED_SCIENCE")

    invalid_presented = len(presented)
    uncertified_definite = 0
    if os.path.exists(arep):
        rep = anomaly.summarise(arep)
        for r in rep["renderers"].values():
            for d in r["diagnostics"]:
                if d["diagnostic"] == "uncertified_value_rendered_definite":
                    uncertified_definite += int(d["measured"] or 0)
    conflations = len(same)

    out = records.ensure_dir(os.path.join(gal, "verification"))
    with open(os.path.join(out, "v0_program_a6_summary.json"), "w",
              encoding="utf-8") as fh:
        json.dump({"gate": "V0-PROGRAM-A6-ENERGY-VISUALS",
                   "num_cases": len(CASES), "num_passed": passed,
                   "num_failed": failed, "num_blocked": blocked,
                   "all_passed": failed == 0,
                   "invalid_observers_presented_as_physical":
                       invalid_presented,
                   "uncertified_signs_presented_as_definite":
                       uncertified_definite,
                   "exact_smooth_objective_conflations": conflations,
                   "blocked_products": BLOCKED,
                   "cases": CASES}, fh, indent=2)

    print("  %d passed, %d failed, %d blocked" % (passed, failed, blocked))
    print("  acceptance: %d invalid observers presented as physical, "
          "%d uncertified signs presented as definite, "
          "%d exact/smooth conflations"
          % (invalid_presented, uncertified_definite, conflations))
    accept = (failed == 0 and invalid_presented == 0
              and uncertified_definite == 0 and conflations == 0)
    print("\nV0-PROGRAM-A6-ENERGY-VISUALS: %s (%d passed, %d failed, "
          "%d blocked)" % ("PASS" if accept else "FAIL", passed, failed,
                           blocked))
    return 0 if accept else 1


if __name__ == "__main__":
    sys.exit(main())
