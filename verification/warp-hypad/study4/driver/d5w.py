"""D5w matched-iterate verdict and gap companion (PROTOCOL-v3 section 8), as pure functions.

Registered table values (not recomputed values) enter G^hi / G^lo (Noether v3 reading 1; builder agrees, logged
before any D5w computation in logs/execution.log).
"""
from __future__ import annotations

import statistics

CONVERGED = ("CONVERGED", "CONVERGED_LS_LIMITED")
REG = {4: dict(s0=84, d0=1.617, pre={20: 2.287, 25: 2.273}),
       6: dict(s0=90, d0=1.146, pre={20: -0.093, 25: -0.084})}
THR = 0.05


def horizon(K: int, runs: list[tuple[str, int]]) -> int:
    """k_c = min(K, terminal iterate of each non-converged run). runs = [(term, k_final)]."""
    k = K
    for term, kf in runs:
        if term not in CONVERGED:
            k = min(k, kf)
    return k


def extended(q: list[float], term: str, K: int) -> tuple[list[float], int]:
    """q[k], k = 0..k_f; converged runs extended as stationary to K. Returns (series, k_f)."""
    kf = len(q) - 1
    if term in CONVERGED and kf < K:
        return q + [q[kf]] * (K - kf), kf
    return list(q), kf


def windows(kc: int, W: int) -> list[tuple[int, int]]:
    return [(a, a + W - 1) for a in range(1, kc + 1, W) if a + W - 1 <= kc]


def rule(series: list[float], pos_label: str, neg_label: str) -> str:
    """>= 2/3 windows >= +0.05 and none <= -0.05 -> pos; mirror -> neg; < 4 windows -> not evaluable."""
    m = len(series)
    if m < 4:
        return "not evaluable (< 4 windows)"
    up = sum(v >= THR for v in series); dn = sum(v <= -THR for v in series)
    if up >= 2 * m / 3 and dn == 0:
        return pos_label
    if dn >= 2 * m / 3 and up == 0:
        return neg_label
    return "not established"


def across_W(a: str, b: str) -> str:
    if "not evaluable" in a or "not evaluable" in b:
        return "not evaluable (< 4 windows)"
    return a if a == b else "not established (window-sensitive)"


def d5w(qA: list[float], termA: str, qB: list[float], termB: str, K: int, density: int,
        b_checks: list[int], b_converged_via_pass: bool) -> dict:
    """qA/qB: in-loop q per accepted study-4 iterate (index = iterate). b_checks: iterates of B's refinement checks."""
    kc = horizon(K, [(termA, len(qA) - 1), (termB, len(qB) - 1)])
    xA, kfA = extended(qA, termA, K); xB, kfB = extended(qB, termB, K)
    reg = REG[density]
    out = dict(k_c=kc, K=K, kfA=kfA, kfB=kfB, termA=termA, termB=termB, dq0_inloop=qA[0] - qB[0],
               dq0_archived=reg["d0"], reconstruction_deviation=abs((qA[0] - qB[0]) - reg["d0"]) > 0.01)
    per = {}
    for W in (20, 25):
        ws = windows(kc, W)
        dw = [statistics.median(xA[a:e + 1]) - statistics.median(xB[a:e + 1]) for a, e in ws]
        n_ext = sum(1 for a, e in ws if e > kfA and termA in CONVERGED or e > kfB and termB in CONVERGED)
        hi = max(reg["pre"][W], reg["d0"]); lo = min(reg["pre"][W], reg["d0"])
        Ghi = [d - hi for d in dw]; Glo = [d - lo for d in dw]
        v = rule(dw, "B lower", "A lower")
        g_grow = rule(Ghi, "gap grows", "_")
        g_shr = rule([-x for x in Glo], "gap shrinks", "_")
        grows = g_grow == "gap grows"; shrinks = g_shr == "gap shrinks"
        single = {}
        for nm, base in (("pre", reg["pre"][W]), ("point", reg["d0"])):
            G = [d - base for d in dw]
            single[nm] = rule(G, "gap grows", "gap shrinks")
        if len(dw) < 4:
            gap = "not evaluable (< 4 windows)"
        elif grows:
            gap = "gap grows"
        elif shrinks:
            gap = "gap shrinks"
        elif any(s in ("gap grows", "gap shrinks") for s in single.values()):
            gap = "not established (baseline-sensitive)"
        else:
            gap = "not established"
        n_after = sum(1 for a, e in ws if b_checks and e > b_checks[0])
        near = [i for i, d in enumerate(dw) if min(abs(d - THR), abs(d + THR), abs(Ghi[i] - THR), abs(Glo[i] + THR)) < 0.001]
        per[W] = dict(windows=ws, delta=dw, G_hi=Ghi, G_lo=Glo, verdict=v, gap=gap, single_baseline=single,
                      m=len(ws), n_ext=n_ext, n_after_chk1=n_after, within_0001_of_threshold=near,
                      maxq_A=[max(xA[a:e + 1]) for a, e in ws], maxq_B=[max(xB[a:e + 1]) for a, e in ws],
                      fracA=[sum(q >= 1.05 for q in xA[a:e + 1]) / W for a, e in ws],
                      fracB=[sum(q >= 1.05 for q in xB[a:e + 1]) / W for a, e in ws])
    out["per_W"] = per
    out["verdict"] = across_W(per[20]["verdict"], per[25]["verdict"])
    out["gap"] = across_W(per[20]["gap"], per[25]["gap"])
    labels = []
    if any(per[W]["m"] and per[W]["n_ext"] >= per[W]["m"] / 2 for W in (20, 25)):
        labels.append("(determined by converged final(s))")
    chk_in = [k for k in b_checks if k <= kc]
    if not b_converged_via_pass:
        if not chk_in:
            labels.append("(no study-4 refinement check)")
        else:
            nW = min((per[W]["n_after_chk1"], per[W]["m"], W) for W in (20, 25))
            if nW[0] < 4:
                labels.append(f"(refinement-limited: {nW[0]} of {nW[1]} windows after the first check)")
    if termA not in CONVERGED or termB not in CONVERGED:
        labels.append("(not converged)")
    out["labels"] = labels
    out["n_check_within_kc"] = len(chk_in)
    refl = any(l.startswith("(no study-4") or l.startswith("(refinement-limited") for l in labels)
    v, g = out["verdict"], out["gap"]
    if v == "B lower":
        base = "B lower at matched cumulative iterates (start gap inherited)"
        if refl:
            lic = base + "; not a formulation effect"
        elif g == "gap grows":
            lic = "B suppresses deepening at matched cumulative iterates (continuation from inherited, unequal starts)"
        elif g == "gap shrinks":
            lic = base + "; gap narrowed against both baselines"
        elif g == "not established (baseline-sensitive)":
            lic = base + "; change in gap not established (baseline-sensitive)"
        else:
            lic = base + "; change not detected"
    elif v == "A lower":
        lic = "A lower at matched cumulative iterates (start gap inherited)"
    else:
        lic = f"no level statement (D5w {v})"
    out["licensed"] = lic + ("" if not labels else " " + " ".join(labels))
    return out
