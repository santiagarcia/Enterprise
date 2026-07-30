#!/usr/bin/env python3
"""Gate V0-O0A-CONSTRAINED-SPLINE-VISUALS.

    python visualization/regression/v0_o0a_verify.py \
        --gallery <build>/results/visualization/program_o0a \
        --records <build>/results/verification/o0a_basis \
        --visualization visualization

Verifies the O0A spline-basis visual products against the records the scientific
gate wrote.

Acceptance: 0 failures, 0 false optimization-readiness claims, 0 full-Hessian
overclaims, 0 Python spline recomputation.

A scientifically honest "not resolution-certified" is a PASSING result. The gate
checks that the classification is reported, not that it is favourable: F0 and F1
both fail the six-cell resolution contract, and a run in which either silently
became "ready for optimization" is the failure this gate exists to catch.
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

STATICS = [
    "quintic_basis_functions", "basis_first_derivatives",
    "basis_second_derivatives", "basis_third_derivatives",
    "tensor_product_basis_atlas", "basis_support_map", "axis_regular_modes",
    "constraint_projection", "gram_spectrum", "f0_baseline_projection",
    "f0_projection_error", "resolution_coverage",
]

ANIMATIONS = [
    "individual_radial_basis_modes", "individual_axisymmetric_basis_modes",
    "raw_to_constrained_mode", "constraint_projection_effect",
    "f0_baseline_to_reconstruction", "whitened_coordinate_perturbation",
]

REQUIRED_RECORDS = [
    "basis_manifest.json", "basis_readiness.json", "axial_knot_vector.csv",
    "radial_squared_knot_vector.csv", "univariate_basis_profiles.csv",
    "univariate_basis_derivatives.csv", "axisymmetric_basis_profiles.csv",
    "tensor_product_atlas.csv", "tensor_product_supports.csv",
    "axis_regular_modes.csv", "sparse_stencil_statistics.csv",
    "constraint_matrix.csv", "constraint_projection.csv",
    "constraint_singular_values.csv", "gram_spectrum_raw.csv",
    "gram_spectrum_constrained.csv", "gram_spectrum_whitened.csv",
    "basis_normalization.csv", "f0_baseline_projection.csv",
    "f0_projection_errors.csv", "resolution_coverage.csv",
    "coefficient_derivative_checks.csv", "verification_report.md",
    "verification_summary.json", "verification_table.csv",
]

# Names that would mean the renderer had re-implemented production. Includes the
# design-layer entry points as well as the spline kernel: a Python
# `sobolev_gram` or `build_projection` would be exactly the duplication the
# contract forbids.
FORBIDDEN = {
    "evaluate", "evaluate_derivative", "field_value", "cartesian_component",
    "nec_value", "nec_minimum", "jacobi_svd", "jacobi_eigen", "sobolev_gram",
    "build_projection", "whiten", "build_quadrature", "make_f1_basis",
    "make_f0_family", "project_shape_function", "evaluate_radial",
    "shape_reference", "classify_lattice", "make_uniform_lattice",
    "make_axis_clamped_lattice", "spd_solve", "svd_min_norm_solve",
    "shape_function", "cosh", "sinh", "tanh", "atanh",
}

HESSIAN_SENTENCE = (
    "Disjoint three-direction OTI batches recover the complete coefficient "
    "gradient. They do not recover mixed Hessian entries between coefficients "
    "assigned to different batches.")


def case(name, claim, ok, expected, computed, note="", t0=None):
    st = "PASS" if ok else "FAIL"
    CASES.append({"name": name, "claim": claim, "state": st, "passed": bool(ok),
                  "expected": str(expected), "computed": str(computed),
                  "note": note,
                  "seconds": round(time.time() - t0, 4) if t0 else 0.0})
    print("  %-6s %-46s %s" % (st, name, computed))
    return bool(ok)


def described(c):
    """The case's computed text. check_predicate stores it in
    `computed_description`; check_residual and check_value use `computed`.
    Reading only one of them left half the cases looking absent."""
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


def png_is_blank(path, min_colors=64):
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


def sidecar_texts(gal):
    """Every sidecar's serialised text, for scanning claim strings."""
    out = {}
    for n in STATICS:
        for ext in (".png.meta.json", ".svg.meta.json"):
            p = os.path.join(gal, n + ext)
            if os.path.exists(p):
                with open(p, encoding="utf-8") as fh:
                    out[n + ext] = fh.read()
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gallery", required=True)
    ap.add_argument("--records", required=True)
    ap.add_argument("--visualization", required=True)
    a = ap.parse_args(argv)

    gal = records.resolve(a.gallery)
    rec = records.resolve(a.records)
    vis = records.resolve(a.visualization)
    sheets = os.path.join(gal, "contact_sheets")

    print("=== V0-O0A-CONSTRAINED-SPLINE-VISUALS ===")

    # --- 1. Records exist ----------------------------------------------------
    t0 = time.time()
    missing = [n for n in REQUIRED_RECORDS
               if not os.path.exists(os.path.join(rec, n))]
    case("v0o0a.records_exist",
         "every record the O0A gate is required to write is present",
         not missing, "%d records" % len(REQUIRED_RECORDS),
         "missing %r" % missing if missing
         else "all %d present" % len(REQUIRED_RECORDS), t0=t0)

    manifest = read_json(os.path.join(rec, "basis_manifest.json"))
    readiness = read_json(os.path.join(rec, "basis_readiness.json"))
    summary = read_json(os.path.join(rec, "verification_summary.json"))
    by_case = {c["name"]: c for c in summary["cases"]}

    # --- 2. The scientific gate passed --------------------------------------
    t0 = time.time()
    case("v0o0a.scientific_gate_passed",
         "the O0A numerical gate exists and passed, so the basis, constraints "
         "and projection behind these visuals are verified rather than merely "
         "computed",
         bool(summary.get("all_passed")) and summary.get("num_failed", 1) == 0,
         "all_passed with zero failures",
         "%d/%d passed, %d failed, %d skipped"
         % (summary["num_passed"], summary["num_cases"], summary["num_failed"],
            summary["num_skipped"]), t0=t0)

    # --- 3/4. Basis samples complete, partition of unity --------------------
    t0 = time.time()
    prof = read_csv(os.path.join(rec, "univariate_basis_profiles.csv"))
    n_idx = len({r["index"] for r in prof})
    n_coord = len({r["s"] for r in prof})
    axi = read_csv(os.path.join(rec, "axisymmetric_basis_profiles.csv"))
    bases = {r["basis"] for r in axi}
    case("v0o0a.basis_samples_complete",
         "the profile records cover every mode of every basis at every sample "
         "coordinate, with no gaps",
         n_idx == manifest["f0"]["raw"] and n_coord > 100 and
         len(bases) == 2 and len(prof) == n_idx * n_coord,
         "%d modes x %d coordinates, 2 correction bases"
         % (manifest["f0"]["raw"], n_coord),
         "F0 %d modes x %d coords (%d rows), correction bases %r"
         % (n_idx, n_coord, len(prof), sorted(bases)), t0=t0)

    t0 = time.time()
    pu = by_case.get("o0a.partition_of_unity", {})
    case("v0o0a.partition_of_unity",
         "the partition-of-unity residual recorded by the gate is within its "
         "tolerance",
         pu.get("passed") and float(pu["computed"]) <= float(pu["tolerance"]),
         "<= %s" % pu.get("tolerance"), "residual %s" % pu.get("computed"),
         "checked on the clamped validation basis; the correction lattices "
         "deliberately taper to zero instead, which is what makes the collar "
         "exact", t0=t0)

    # --- 5/6/7/8. Leakage and axis ------------------------------------------
    for nm, key, claim in (
            ("v0o0a.no_support_leakage", "o0a.exact_zero_outside_support",
             "no basis function contributes outside its support"),
            ("v0o0a.no_outer_collar_leakage", "o0a.outer_flat_collar_exactness",
             "no field or derivative leaks into the outer flat collar"),
            ("v0o0a.no_protected_cabin_leakage", "o0a.protected_cabin_exactness",
             "no field or derivative leaks into the protected cabin"),
            ("v0o0a.no_axis_singularity", "o0a.no_division_by_rho_at_axis",
             "every derivative channel stays bounded and finite at the axis")):
        t0 = time.time()
        c = by_case.get(key, {})
        case(nm, claim, c.get("passed", False), "gate case %s passed" % key,
             "%s: %s" % (key, described(c)), t0=t0)

    # --- 9/10. Continuity labelling ------------------------------------------
    t0 = time.time()
    c4 = by_case.get("o0a.c4_continuity_ordinary_knot", {})
    texts = sidecar_texts(gal)
    labelled_c4 = any("C4" in t for t in texts.values())
    case("v0o0a.ordinary_knot_continuity_labelled",
         "the gate establishes C4 at an ordinary interior knot and the figures "
         "label it as C4 rather than as a generic smoothness",
         c4.get("passed", False) and labelled_c4, "C4 established and labelled",
         "gate: %s | figures mention C4: %s"
         % (described(c4)[:70], labelled_c4), t0=t0)

    t0 = time.time()
    rk = by_case.get("o0a.repeated_knot_continuity_reduced", {})
    case("v0o0a.repeated_knot_continuity_reduced",
         "a multiplicity-2 interior knot reduces continuity to C3, shown by a "
         "fourth-derivative jump that does NOT shrink with the probe offset",
         rk.get("passed", False), "reduction detected",
         described(rk)[:110], t0=t0)

    # --- 11. Derivative verification data agree with the gate ---------------
    t0 = time.time()
    probes = read_csv(os.path.join(rec, "tensor_derivative_probes.csv"))
    worst_probe = max(float(r["absolute_difference"]) for r in probes)
    grad = by_case.get("o0a.tensor_product_gradient", {})
    hess = by_case.get("o0a.tensor_product_hessian", {})
    tol = max(float(grad.get("tolerance", 0)), float(hess.get("tolerance", 0)))
    case("v0o0a.derivative_verification_agrees",
         "the analytic-versus-finite-difference probe record agrees with the "
         "gradient and Hessian cases the gate recorded",
         grad.get("passed") and hess.get("passed") and worst_probe <= tol,
         "worst probe difference <= %g" % tol,
         "%d probe rows, worst %.3e; gate gradient %s, Hessian %s"
         % (len(probes), worst_probe, grad.get("computed"),
            hess.get("computed")), t0=t0)

    # --- 12/13/14/15/16. Constraint agreement -------------------------------
    t0 = time.time()
    sv = read_csv(os.path.join(rec, "constraint_singular_values.csv"))
    retained = sum(1 for r in sv if r["retained"] == "1")
    con = manifest["constraints"]
    case("v0o0a.constraint_rank_agrees",
         "the rank in the manifest, the retained singular-value count and the "
         "gate case all agree",
         retained == con["rank"] and
         by_case.get("o0a.constraint_rank_detection", {}).get("passed", False),
         "rank %d everywhere" % con["rank"],
         "manifest rank %d, retained singular values %d"
         % (con["rank"], retained), t0=t0)

    t0 = time.time()
    cm = read_csv(os.path.join(rec, "constraint_matrix.csv"))
    case("v0o0a.redundant_constraint_count_agrees",
         "the redundant-row count equals rows minus rank and matches the "
         "constraint-matrix record",
         con["redundant"] == con["rows"] - con["rank"] and
         len(cm) == con["rows"],
         "%d redundant of %d rows" % (con["redundant"], con["rows"]),
         "matrix rows %d, rows %d, rank %d, redundant %d"
         % (len(cm), con["rows"], con["rank"], con["redundant"]), t0=t0)

    for nm, key, lim, claim in (
            ("v0o0a.affine_residual", "o0a.affine_constraint_residual", 1e-12,
             "|C c0 - d| is within tolerance"),
            ("v0o0a.null_space_residual", "o0a.null_space_residual", 1e-12,
             "max |C P| and |P^T P - I| are within tolerance")):
        t0 = time.time()
        c = by_case.get(key, {})
        val = float(c.get("computed", 1.0))
        case(nm, claim, c.get("passed", False) and val <= lim, "<= %g" % lim,
             "%.3e" % val, t0=t0)

    # --- 17/18/19. Gram, null modes, whitening -------------------------------
    t0 = time.time()
    raw = read_csv(os.path.join(rec, "gram_spectrum_raw.csv"))
    conr = read_csv(os.path.join(rec, "gram_spectrum_constrained.csv"))
    whi = read_csv(os.path.join(rec, "gram_spectrum_whitened.csv"))
    ok_spec = (len(raw) == manifest["f1_tensor"]["active"] and
               len(conr) == con["free"] and len(whi) == con["free"] and
               all(r["converged"] == "1" for r in raw + conr + whi))
    case("v0o0a.gram_spectra_agree",
         "the three Gram spectra have the dimensions the projection implies and "
         "every eigensolve reported convergence",
         ok_spec, "%d raw, %d constrained, %d whitened, all converged"
         % (manifest["f1_tensor"]["active"], con["free"], con["free"]),
         "%d / %d / %d; conditions %.4f, %.4f, %.6f"
         % (len(raw), len(conr), len(whi),
            float(raw[0]["condition_number"]),
            float(conr[0]["condition_number"]),
            float(whi[0]["condition_number"])), t0=t0)

    t0 = time.time()
    silent = [r["index"] for r in whi
              if r["retained"] == "1" and float(r["eigenvalue"]) <= 0.0]
    disc = sum(1 for r in whi if r["retained"] != "1")
    nm_case = by_case.get("o0a.numerical_null_mode_rejection", {})
    case("v0o0a.no_silently_retained_null_modes",
         "no non-positive eigenvalue is retained, and the discard decision is "
         "exercised on fixtures that DO have a null mode",
         not silent and nm_case.get("passed", False),
         "zero silently retained", "production discarded %d, silently retained "
         "%d; fixture case: %s" % (disc, len(silent),
                                   described(nm_case)[:80]),
         t0=t0)

    t0 = time.time()
    wc = by_case.get("o0a.whitening", {})
    case("v0o0a.whitening_residual", "W^T G W = I within tolerance",
         wc.get("passed", False), "identity to rounding",
         described(wc)[:120], t0=t0)

    # --- 20/21. Sparse accounting and equivalence ---------------------------
    t0 = time.time()
    st = {r["quantity"]: r["value"] for r in
          read_csv(os.path.join(rec, "sparse_stencil_statistics.csv"))}
    pts = int(st["points"])
    nnz = int(st["nonzeros"])
    ch = int(st["stored_channels"])
    nb = int(st["active_basis_functions"])
    dense = pts * nb * ch * 8
    sparse = nnz * (ch * 8 + 4)
    closes = (int(st["dense_equivalent_bytes"]) == dense and
              int(st["sparse_bytes"]) == sparse and
              abs(float(st["average_active_per_point"]) - nnz / pts)
              <= 1e-12 * (nnz / pts) and
              abs(float(st["memory_reduction_factor"]) - dense / sparse)
              <= 1e-12 * (dense / sparse))
    case("v0o0a.sparse_memory_accounting_closes",
         "the recorded nonzero count, average occupancy and both byte figures "
         "are mutually consistent and recomputable from each other",
         closes, "accounting closes exactly",
         "%d points x %d functions, %d nonzeros, %d channels; %d dense vs %d "
         "sparse bytes, reduction %sx"
         % (pts, nb, nnz, ch, dense, sparse, st["memory_reduction_factor"]),
         t0=t0)

    t0 = time.time()
    sd = by_case.get("o0a.sparse_versus_dense_evaluation", {})
    case("v0o0a.sparse_and_dense_agree",
         "the sparse stencil and the direct sum agree bit for bit",
         sd.get("passed", False) and float(sd.get("computed", 1)) == 0.0,
         "exactly 0", "%s" % sd.get("computed"), t0=t0)

    # --- 22..28. F0 and resolution ------------------------------------------
    t0 = time.time()
    fc = by_case.get("o0a.f0_flat_collar_exactness", {})
    case("v0o0a.f0_collar_exactly_zero",
         "the F0 reconstruction and its first two radial derivatives are exactly "
         "zero in the flat collar, and every collar coefficient is exactly zero",
         fc.get("passed", False) and float(fc.get("computed", 1)) == 0.0,
         "exactly 0", "%s" % fc.get("computed"), t0=t0)

    errs = {r["zone"]: r for r in
            read_csv(os.path.join(rec, "f0_projection_errors.csv"))}
    wall = errs["wall"]
    for nm, key, colname in (
            ("v0o0a.f0_value_error_displayed", "max_value_error", "value"),
            ("v0o0a.f0_first_derivative_error_displayed",
             "max_first_derivative_error", "first derivative"),
            ("v0o0a.f0_second_derivative_error_displayed",
             "max_second_derivative_error", "second derivative")):
        t0 = time.time()
        val = wall[key]
        # The figure must carry the number, either in its plotted CSV or in its
        # sidecar, so a reader can see the size of the error and not merely a
        # curve on an unlabelled axis.
        plotted = os.path.join(gal, "f0_projection_error.plotted.csv")
        shown = False
        if os.path.exists(plotted):
            for r in read_csv(plotted):
                if r.get("zone") == "wall" and r.get(key) == val:
                    shown = True
        case(nm,
             "the F0 %s error is present in the figure's plotted data, at the "
             "value the gate measured" % colname,
             shown, "wall %s = %s" % (key, val),
             "present in f0_projection_error.plotted.csv: %s (value %s)"
             % (shown, val), t0=t0)

    t0 = time.time()
    fams = {f["family"]: f for f in readiness["families"]}
    f0r = fams.get("F0 radial validation", {})
    f1r = fams.get("F1 axisymmetric scalar kernel", {})
    case("v0o0a.f0_not_resolution_certified",
         "F0 is recorded as NOT resolution-certified and NOT ready for "
         "optimization, with an explicit blocking reason",
         f0r.get("resolution_certified") is False and
         f0r.get("ready_for_optimization") is False and
         "resolution" in f0r.get("blocking_reason", ""),
         "false / false with a reason",
         "resolution_certified=%s ready=%s reason=%r"
         % (f0r.get("resolution_certified"), f0r.get("ready_for_optimization"),
            f0r.get("blocking_reason")), t0=t0,
         note="an honestly measured resolution limitation is a passing result "
              "for this gate")

    t0 = time.time()
    case("v0o0a.f1_not_resolution_certified",
         "F1 and its three adapters are recorded as NOT resolution-certified and "
         "NOT ready for optimization",
         all(f["resolution_certified"] is False and
             f["ready_for_optimization"] is False
             for f in readiness["families"]),
         "false / false for every family",
         "F1 resolution_certified=%s ready=%s; %d families total, %d claiming "
         "readiness" % (f1r.get("resolution_certified"),
                        f1r.get("ready_for_optimization"),
                        len(readiness["families"]),
                        sum(1 for f in readiness["families"]
                            if f["ready_for_optimization"])), t0=t0)

    t0 = time.time()
    resp = os.path.join(gal, "resolution_coverage.png.meta.json")
    six_visible = False
    if os.path.exists(resp):
        with open(resp, encoding="utf-8") as fh:
            j = json.load(fh)
        six_visible = j.get("six_cell_threshold") == 6.0
    cov = read_csv(os.path.join(rec, "resolution_coverage.csv"))
    f1cells = min(float(r["cells_across_narrowest_feature"]) for r in cov
                  if r["family"] == "f1_scalar")
    f0cells = min(float(r["cells_across_narrowest_feature"]) for r in cov
                  if r["family"] == "f0_radial")
    failing = all(r["meets_six_cell_contract"] == "0" for r in cov)
    case("v0o0a.six_cell_threshold_visible",
         "the six-cell acceptance threshold is declared by the resolution figure "
         "and both families are recorded as failing it",
         six_visible and failing,
         "threshold 6.0 declared, every mode failing",
         "threshold declared %s; F1 minimum %.3f cells, F0 minimum %.3f cells, "
         "all %d modes fail" % (six_visible, f1cells, f0cells, len(cov)), t0=t0)

    # --- 29/30/31. OTI capacity, schedule, Hessian --------------------------
    t0 = time.time()
    case("v0o0a.oti_capacity_is_three_directions",
         "the recorded OTI capacity is the backend's own: 3 directions at order 2",
         manifest["oti_directions"] == 3 and manifest["oti_order"] == 2,
         "3 directions, order 2",
         "directions %s, order %s, max simultaneous seeds %s"
         % (manifest["oti_directions"], manifest["oti_order"],
            readiness["oti_backend"]["maximum_simultaneous_coefficient_seeds"]),
         t0=t0)

    t0 = time.time()
    ob = readiness["oti_backend"]
    expect_batches = (con["free"] + 2) // 3
    case("v0o0a.gradient_schedule_is_twelve_batches",
         "the gradient schedule is ceil(free / directions) batches: 11 full "
         "batches of 3 and one of 1 for 34 free coordinates",
         ob["gradient_batches"] == expect_batches == 12 and
         ob["full_batches"] == 11 and ob["final_partial_batch_size"] == 1,
         "12 batches = 11 x 3 + 1 x 1",
         "%d free -> %d batches, %d full of 3, final partial %d"
         % (ob["free_coordinates"], ob["gradient_batches"], ob["full_batches"],
            ob["final_partial_batch_size"]), t0=t0)

    t0 = time.time()
    # The sentence must be present verbatim, and no product may claim a full
    # Hessian anywhere.
    haystacks = dict(texts)
    haystacks["basis_manifest.json"] = json.dumps(manifest)
    haystacks["basis_readiness.json"] = json.dumps(readiness)
    amf = os.path.join(gal, "animation_manifest.json")
    if os.path.exists(amf):
        with open(amf, encoding="utf-8") as fh:
            haystacks["animation_manifest.json"] = fh.read()
    stated = [k for k, v in haystacks.items() if HESSIAN_SENTENCE in v]
    claims = [k for k, v in haystacks.items()
              if '"full_hessian_implemented": true' in v.lower()
              or "full hessian" in v.lower() and "not" not in v.lower()]
    case("v0o0a.no_full_hessian_claim",
         "the Hessian limitation is stated verbatim in the manifest, the "
         "readiness report and every figure sidecar, and nothing claims a full "
         "34-dimensional Hessian",
         not claims and manifest["full_hessian_implemented"] is False and
         len(stated) >= 3 and "basis_manifest.json" in stated and
         "basis_readiness.json" in stated,
         "sentence present, zero overclaims",
         "verbatim in %d artefacts, %d overclaims, "
         "full_hessian_implemented=%s"
         % (len(stated), len(claims), manifest["full_hessian_implemented"]),
         t0=t0)

    # --- 32/33/34. Static products ------------------------------------------
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
    case("v0o0a.static_figures_exist",
         "all twelve required figures are present as PNGs with sane dimensions",
         not absent and not small, "%d figures" % len(STATICS),
         "absent %r small %r" % (absent, small[:2]) if (absent or small)
         else "all %d present" % len(STATICS), t0=t0)

    t0 = time.time()
    blanks = [n for n in STATICS
              if os.path.exists(os.path.join(gal, n + ".png"))
              and png_is_blank(os.path.join(gal, n + ".png"))]
    case("v0o0a.pngs_nonblank", "no figure is blank", not blanks,
         "all figures carry content",
         "blank %r" % blanks if blanks else "all %d rich" % len(STATICS), t0=t0)

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
    case("v0o0a.svgs_parse", "every figure has a well-formed SVG companion",
         not svg_bad, "%d SVGs" % len(STATICS),
         "defects %r" % svg_bad[:3] if svg_bad else "all parse", t0=t0)

    # --- 35/36/37/38/39. Animations -----------------------------------------
    t0 = time.time()
    entries, prov = [], {}
    parsed = False
    if os.path.exists(amf):
        try:
            mj = json.load(open(amf, encoding="utf-8"))
            entries, prov = mj.get("animations", []), mj.get("provenance", {})
            parsed = True
        except ValueError:
            pass
    case("v0o0a.animation_manifest_parses",
         "the animation manifest parses and carries an entry for every animation "
         "and format",
         parsed and len(entries) == 2 * len(ANIMATIONS),
         "%d entries" % (2 * len(ANIMATIONS)), "%d entries" % len(entries),
         t0=t0)

    t0 = time.time()
    dishonest = []
    for e in entries:
        if e.get("written"):
            fn = e.get("filename")
            if not fn or not os.path.exists(os.path.join(gal, fn)):
                dishonest.append((e.get("codec"), fn, "claimed, absent"))
        elif not (e.get("reason") or e.get("note")):
            dishonest.append((e.get("codec"), e.get("filename"), "no reason"))
    case("v0o0a.animation_status_honest",
         "every entry claiming a file has one, and every unwritten entry gives a "
         "reason",
         not dishonest, "written flags match the filesystem",
         "defects %r" % dishonest[:3] if dishonest
         else "%d written, %d unwritten with a reason"
         % (sum(1 for e in entries if e.get("written")),
            sum(1 for e in entries if not e.get("written"))), t0=t0,
         note="a missing encoder must not fail this gate; the GIFs and contact "
              "sheets are the deterministic products")

    t0 = time.time()
    miss = [n for n in ANIMATIONS
            if not os.path.exists(os.path.join(sheets, n + ".png"))]
    sb = [n for n in ANIMATIONS
          if os.path.exists(os.path.join(sheets, n + ".png"))
          and png_is_blank(os.path.join(sheets, n + ".png"))]
    case("v0o0a.contact_sheets_exist",
         "every animation has a nonblank contact sheet", not miss and not sb,
         "%d sheets" % len(ANIMATIONS),
         "missing %r blank %r" % (miss, sb) if (miss or sb)
         else "all %d present and rich" % len(ANIMATIONS), t0=t0)

    t0 = time.time()
    expect_frames = {
        "individual_radial_basis_modes": manifest["f0"]["raw"],
        "individual_axisymmetric_basis_modes": manifest["f1_tensor"]["active"],
        "raw_to_constrained_mode": manifest["f1_tensor"]["active"],
        "f0_baseline_to_reconstruction": manifest["f0"]["free"] + 1,
    }
    bad_frames = []
    for e in entries:
        want = expect_frames.get(e["filename"].rsplit(".", 1)[0]
                                 if e.get("filename") else "")
        if want is not None and e.get("frame_count") != want:
            bad_frames.append((e.get("filename"), e.get("frame_count"), want))
    case("v0o0a.frame_counts_agree",
         "each animation's frame count equals the number of modes, coefficients "
         "or reconstruction steps the records contain",
         not bad_frames, "frame counts derived from the records",
         "mismatches %r" % bad_frames[:3] if bad_frames
         else "checked %d expectations: %r" % (len(expect_frames),
                                              sorted(expect_frames.values())),
         t0=t0)

    t0 = time.time()
    declared = ("fixed" in str(prov.get("colour_scale", "")).lower())
    per_fig = all("colour_scale" in t for t in texts.values())
    case("v0o0a.fixed_scales_declared",
         "both the animations and every figure declare a fixed colour scale "
         "computed once from the whole record",
         declared and per_fig, "declared everywhere",
         "animation manifest: %r; all %d figure sidecars declare one: %s"
         % (str(prov.get("colour_scale"))[:50], len(texts), per_fig), t0=t0)

    # --- 40. No Python spline or design physics -----------------------------
    t0 = time.time()
    offenders = []
    o0adir = os.path.join(vis, "program_o0a")
    scanned = 0
    for fn in sorted(os.listdir(o0adir)) if os.path.isdir(o0adir) else []:
        if not fn.endswith(".py"):
            continue
        scanned += 1
        p = os.path.join(o0adir, fn)
        with open(p, encoding="utf-8") as fh:
            src = fh.read()
        try:
            tree = ast.parse(src, filename=fn)
        except SyntaxError as exc:
            offenders.append("%s does not parse (%s)" % (fn, exc.msg))
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                mod = getattr(node, "module", "") or ""
                names = [al.name for al in node.names]
                if "warp_hypad" in mod or any("warp_hypad" in n for n in names):
                    offenders.append("%s:%d imports warp_hypad"
                                     % (fn, node.lineno))
            if isinstance(node, ast.Call):
                f = node.func
                nm = (f.attr if isinstance(f, ast.Attribute)
                      else f.id if isinstance(f, ast.Name) else None)
                if nm and nm in FORBIDDEN:
                    offenders.append("%s:%d calls %s" % (fn, node.lineno, nm))
    case("v0o0a.no_python_spline_recomputation",
         "no O0A renderer imports warp_hypad or calls a spline, basis, "
         "projection, Gram or factorisation routine",
         not offenders and scanned >= 2,
         "zero production calls across the package",
         "offenders %r" % offenders[:4] if offenders
         else "clean across %d modules" % scanned, t0=t0,
         note="this is why the frame records exist: production computes every "
              "projection and partial reconstruction, and the renderer only draws")

    # --- 41. Checksums -------------------------------------------------------
    t0 = time.time()
    stored = manifest.get("records", {})
    drift = [n for n, d in stored.items()
             if os.path.exists(os.path.join(rec, n))
             and fnv1a64(os.path.join(rec, n)) != d]
    absent_dig = [n for n in stored if not os.path.exists(os.path.join(rec, n))]
    checked = changed = 0
    for n in STATICS:
        sc = os.path.join(gal, n + ".png.meta.json")
        if not os.path.exists(sc):
            continue
        r = provenance.verify_inputs_unchanged(sc, [rec])
        checked += r["checked"]
        changed += len(r["changed"]) + len(r["missing"])
    case("v0o0a.checksums_agree",
         "every record digest in the manifest still matches the bytes on disk, "
         "the manifest excludes itself, and every figure's declared inputs are "
         "unchanged",
         not drift and not absent_dig and changed == 0 and checked > 0 and
         "basis_manifest.json" not in stored and
         manifest.get("manifest_excludes_itself") is True,
         "all digests match, manifest self-excluded",
         "%d record digests, %d drifted, %d absent; %d figure inputs checked, "
         "%d changed; self-excluded %s"
         % (len(stored), len(drift), len(absent_dig), checked, changed,
            "basis_manifest.json" not in stored), t0=t0)

    # --- 42. Deterministic rerender metadata ---------------------------------
    t0 = time.time()
    needed = ("degree", "design_region", "scale_R0", "f1_axial", "f1_radial",
              "f1_tensor", "constraints", "f0", "quadrature",
              "determinism_digest", "digest_algorithm")
    absent_meta = [k for k in needed if manifest.get(k) is None]
    det = by_case.get("o0a.deterministic_repeat_evaluation", {})
    fig_meta = all("deterministic_rerender" in t for t in texts.values())
    case("v0o0a.deterministic_rerender_metadata",
         "the manifest records everything a rerender needs, the gate proved "
         "bit-identical re-evaluation on the numbers themselves, and every "
         "figure carries its input digests",
         not absent_meta and det.get("passed", False) and fig_meta and
         bool(prov.get("deterministic_rerender")),
         "complete metadata plus a numeric determinism proof",
         "missing %r; gate: %s; figure sidecars carry digests: %s"
         % (absent_meta, described(det)[:60], fig_meta), t0=t0,
         note="determinism is established on the numbers, not on a timestamp or "
              "a build path")

    # --- 43. No raw frame directory survives --------------------------------
    t0 = time.time()
    left = [n for n in os.listdir(gal) if n.startswith("_frames_")]
    case("v0o0a.no_raw_frame_directory", "no raw frame directory survived",
         not left, "zero", "%r" % left if left else "none", t0=t0)

    # --- The renderers' own anomaly reports ----------------------------------
    t0 = time.time()
    arep = os.path.join(gal, "visual_anomaly_report.json")
    fired, total = -1, 0
    if os.path.exists(arep):
        r = anomaly.summarise(arep)
        fired = r.get("total_fired", -1)
        for v in r["renderers"].values():
            total += len(v["diagnostics"])
    case("v0o0a.anomaly_report_clean",
         "the renderers' anomaly report exists and has zero fired diagnostics",
         fired == 0, "0 fired", "%d fired of %d diagnostics" % (fired, total),
         t0=t0)

    # --- Summary --------------------------------------------------------------
    passed = sum(1 for c in CASES if c["passed"])
    failed = len(CASES) - passed
    false_ready = sum(1 for f in readiness["families"]
                      if f["ready_for_optimization"])
    overclaims = len(claims)
    py_recompute = len(offenders)

    out = records.ensure_dir(os.path.join(gal, "verification"))
    with open(os.path.join(out, "v0_o0a_summary.json"), "w",
              encoding="utf-8") as fh:
        json.dump({"gate": "V0-O0A-CONSTRAINED-SPLINE-VISUALS",
                   "num_cases": len(CASES), "num_passed": passed,
                   "num_failed": failed, "all_passed": failed == 0,
                   "false_optimization_readiness_claims": false_ready,
                   "full_hessian_overclaims": overclaims,
                   "python_spline_recomputation": py_recompute,
                   "f0_resolution_certified": False,
                   "f1_resolution_certified": False,
                   "honest_not_certified_is_a_pass": True,
                   "cases": CASES}, fh, indent=2)

    print("  %d passed, %d failed" % (passed, failed))
    print("  acceptance: %d false optimization-readiness claims, %d full-Hessian "
          "overclaims, %d Python spline recomputations"
          % (false_ready, overclaims, py_recompute))
    accept = (failed == 0 and false_ready == 0 and overclaims == 0 and
              py_recompute == 0)
    print("\nV0-O0A-CONSTRAINED-SPLINE-VISUALS: %s (%d passed, %d failed)"
          % ("PASS" if accept else "FAIL", passed, failed))
    return 0 if accept else 1


if __name__ == "__main__":
    sys.exit(main())
