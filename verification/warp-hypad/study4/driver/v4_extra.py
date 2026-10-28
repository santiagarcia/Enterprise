"""Study-4 v4 additions not covered by analyze_runs.py / hess_analyze.py (PROTOCOL-v4 sections 2.2 NV3, 2.3 NV4,
7 NV6, 8 D3 / D6 / D6x and the NV3 limitation wording).

Usage (Anaconda python -I):
  python -I v4_extra.py STUDY4_DIR eps DIR [DIR ...]       # eps*_min, eps*_all and NV4 label per Hessian directory
  python -I v4_extra.py STUDY4_DIR nv6prep RUN [RUN ...]    # writes hessians/final_RUN/nv6_{plus,minus}.txt, nv6_pred.json
  python -I v4_extra.py STUDY4_DIR nv6eval RUN [RUN ...]    # reads nv6_{zero,plus,minus}.dump (s4 v2dump output)
  python -I v4_extra.py STUDY4_DIR decide                   # D3, D6, D6x from results/decisions.json -> results/v4_decisions.json
"""
from __future__ import annotations

import csv
import io
import json
import sys
from pathlib import Path

import numpy as np

N = 108
U = 2.0 ** -53
NV3_LIMIT = ("τ_H covers node-sum roundoff; a sign change by per-node evaluation error requires a coherent relative "
             "error ≥ ε*_all (reported).")


def meta(d: Path) -> dict:
    m = {}
    for line in (d / "meta.txt").read_text().splitlines():
        k, _, v = line.partition(" ")
        m[k] = v
    return m


def mat(d: Path, f: str) -> np.ndarray:
    return np.fromfile(d / f, dtype="<f8").reshape(N, N)


def spectrum(d: Path, lam: float | None = None) -> dict:
    """H_F spectrum, tau_H (as hess_analyze.py), eps* (NV3) and the NV4 label for one Hessian directory."""
    m = meta(d)
    lam = float(m["lambda"]) if lam is None else lam
    Hf, Hr, Hc, B = (mat(d, f) for f in ("H_R_fwd.bin", "H_R_rev.bin", "H_R_comp.bin", "B_H.bin"))
    y = np.fromfile(d / "y.bin", dtype="<f8")
    gF = np.fromfile(d / "gR_fwd.bin", dtype="<f8") + lam * y
    HF = Hf + lam * np.eye(N)
    mu, Q = np.linalg.eigh(HF)
    eh = max(np.linalg.norm(Hf - Hr, 2), np.linalg.norm(Hf - Hc, 2))
    mumax = float(np.max(np.abs(mu)))
    tau = max(1e-9 * mumax, 10 * eh)
    n_nodes = int(m["nodes"])
    S2 = float(np.linalg.norm(B / ((n_nodes - 1) * U), 2))
    F = float(m["R_fwd"]) + 0.5 * lam * float(y @ y)
    tol = 1e-5 * F / max(float(np.linalg.norm(y)), 1.0)
    neg = mu < -tau
    gn = float(np.linalg.norm(gF))
    c1 = gn <= tol
    nneg = int(neg.sum())
    if nneg:
        lab = f"indefinite ({nneg} negative)"
        if not c1:
            lab += f", non-stationary (‖∇F‖/tol_C1' = {gn / tol:.3g})"
    elif np.any(np.abs(mu) <= tau):
        lab = "singular to τ_H"
    else:
        lab = "PD"
    return dict(design=m.get("tag", d.name), dir=str(d), lam=lam, nodes=n_nodes, mu=mu, Q=Q, gF=gF, y=y, tau=tau,
                mu_min=float(mu[0]), S2=S2, eps_min=abs(float(mu[0])) / S2, eps_all=float(np.min(np.abs(mu))) / S2,
                F=F, tol=tol, gF_over_tol=gn / tol, C1=c1, label_v4=lab, b_F=float(m["R_bound"]), neg=neg)


def cmd_eps(root: Path, dirs: list[str]) -> None:
    rows = []
    for dd in dirs:
        s = spectrum(Path(dd))
        rows.append({k: s[k] for k in ("design", "nodes", "lam", "mu_min", "tau", "S2", "eps_min", "eps_all",
                                       "gF_over_tol", "label_v4")})
    out = root / "results" / "v4_eps.csv"
    new = not out.exists()
    with open(out, "a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        if new:
            w.writeheader()
        for r in rows:
            w.writerow(r)
            print(r)


def cmd_nv6prep(root: Path, runs: list[str]) -> None:
    for r in runs:
        d = root / "hessians" / f"final_{r}"
        s = spectrum(d)
        Pn = s["Q"][:, s["neg"]]
        pg = float(np.linalg.norm(Pn @ (Pn.T @ s["gF"]))) if Pn.shape[1] else 0.0
        w = s["Q"][:, 0]
        mu0 = s["mu_min"]
        alpha = float(np.sqrt(2e3 * s["b_F"] / abs(mu0)))
        for nm, sg in (("plus", 1.0), ("minus", -1.0)):
            np.savetxt(d / f"nv6_{nm}.txt", s["y"] + sg * alpha * w, fmt="%.17g")
        np.savetxt(d / "nv6_zero.txt", s["y"], fmt="%.17g")
        pred = dict(run=r, n_neg=int(s["neg"].sum()), Pneg_gF_over_tol=pg / s["tol"], mu_min=mu0, tau_H=s["tau"],
                    b_F=s["b_F"], alpha=alpha, predicted_dF=0.5 * mu0 * alpha ** 2, minus_1e3_bF=-1e3 * s["b_F"],
                    eps_min=s["eps_min"], eps_all=s["eps_all"])
        json.dump(pred, open(d / "nv6_pred.json", "w"), indent=1)
        print(pred)


def read_dump(p: Path) -> float:
    for line in p.read_text().splitlines():
        if line.startswith("F "):
            return float.fromhex(line.split()[1])
    raise ValueError(p)


def cmd_nv6eval(root: Path, runs: list[str]) -> None:
    res = {}
    for r in runs:
        d = root / "hessians" / f"final_{r}"
        pred = json.load(open(d / "nv6_pred.json"))
        F0, Fp, Fm = (read_dump(d / f"nv6_{k}.dump") for k in ("zero", "plus", "minus"))
        pred.update(F0=F0, dF_plus=Fp - F0, dF_minus=Fm - F0,
                    ratio_plus=(Fp - F0) / pred["predicted_dF"], ratio_minus=(Fm - F0) / pred["predicted_dF"],
                    confirms=bool(Fp - F0 < 0 and Fm - F0 < 0))
        json.dump(pred, open(d / "nv6_result.json", "w"), indent=1)
        res[r] = pred
        print(pred)


def same_point(a: dict, b: dict, root: Path) -> bool:
    ya = np.loadtxt(root / "runs" / a["run"] / "y_final.txt")
    yb = np.loadtxt(root / "runs" / b["run"] / "y_final.txt")
    return (abs(a["q_final"] - b["q_final"]) < 0.01
            and float(np.linalg.norm(ya - yb) / np.linalg.norm(ya)) < 1e-3)


def strict(rec: dict) -> tuple[bool, bool]:
    lab = rec.get("D1b", "")
    return lab.startswith("strict"), lab.endswith("(marginal)")


def conc(rec: dict) -> str:
    q, rE = rec["q_final"], rec["rE_final"]
    return "pocket" if q >= 1.05 else "shallow"


def cmd_decide(root: Path) -> None:
    dec = json.load(open(root / "results" / "decisions.json"))
    R = dec["runs"]
    out: dict = {"D1b_limitation_v4": NV3_LIMIT}

    def conv(r):
        return R.get(r, {}).get("converged") is True

    # D3 per formulation at lambda2 on P1 (study-3 rule)
    d3 = {}
    for form, (r4, r6) in (("A", ("A4L2m", "A6L2m")), ("B", ("B4L2", "B6L2"))):
        a, b = R[r4], R[r6]
        if not (conv(r4) and conv(r6)):
            d3[form] = (f"failed (not converged): {r4} {a.get('state')}, {r6} {b.get('state')}"
                        + (f"; n_outer {a.get('n_outer')}/{b.get('n_outer')}, n_check {a.get('n_check')}/{b.get('n_check')}"
                           if form == "B" else ""))
            continue
        qa, qb = a["q_final"], b["q_final"]
        if qa >= 1.05 and qb >= 1.05 and a["rE_final"] >= 0.10 and b["rE_final"] >= 0.10 and abs(qa - qb) < 0.05:
            v = "concentration at convergence"
        elif qa < 1.01 and qb < 1.01 and abs(qa - qb) < 0.05:
            v = "no concentration at convergence"
        else:
            v = "resolution-dependent"
        if any(strict(x)[1] for x in (a, b)):
            v += " (marginal)"
        d3[form] = f"{v}; D1b {a['D1b']} / {b['D1b']}"
    out["D3_lambda2_P1"] = d3
    # D3 P1-A vs P2-A per density
    d3b = {}
    for dens, (p1, p2) in ((4, ("A4L2m", "A4L2f")), (6, ("A6L2m", "A6L2f"))):
        if R[p1].get("state") == "not run" or R[p2].get("state") == "not run":
            d3b[dens] = "not run"
            continue
        c1, c2 = conv(p1), conv(p2)
        if c1 and c2:
            if same_point(R[p1], R[p2], root):
                d3b[dens] = "same point"
            else:
                s1, m1 = strict(R[p1]); s2, m2 = strict(R[p2])
                if s1 and s2:
                    d3b[dens] = "distinct local minimisers" + (" (marginal)" if (m1 or m2) else "")
                else:
                    d3b[dens] = "distinct stationary points"
        else:
            d3b[dens] = " and ".join(x for x, c in (("P1 not converged", c1), ("P2 not converged", c2)) if not c)
    out["D3_P1A_vs_P2A"] = d3b
    # D6 lambda dependence (needs D3 at both lambda; lambda1 D3 from P4)
    l1 = [R.get(r, {}) for r in ("A4L1f", "A6L1f", "B4L1f", "B6L1f")]
    out["D6"] = ("not evaluable (P4 not run)" if any(x.get("state", "not run") == "not run" for x in l1)
                 else "see D3 at lambda1 (computed below)")
    if out["D6"].startswith("see"):
        d31 = {}
        for form, (r4, r6) in (("A", ("A4L1f", "A6L1f")), ("B", ("B4L1f", "B6L1f"))):
            a, b = R[r4], R[r6]
            if not (conv(r4) and conv(r6)):
                d31[form] = "failed (not converged)"
                continue
            qa, qb = a["q_final"], b["q_final"]
            if qa >= 1.05 and qb >= 1.05 and a["rE_final"] >= 0.10 and b["rE_final"] >= 0.10 and abs(qa - qb) < 0.05:
                d31[form] = "concentration at convergence"
            elif qa < 1.01 and qb < 1.01 and abs(qa - qb) < 0.05:
                d31[form] = "no concentration at convergence"
            else:
                d31[form] = "resolution-dependent"
        out["D3_lambda1_P4"] = d31
        out["D6"] = {f: ("regularisation-dependent" if d31[f] != d3[f].split(";")[0] else "same verdict")
                     if "failed" not in d31[f] and "failed" not in d3[f] else "not evaluable (D3 failed)" for f in d31}
    # D6x
    arms = {"X4L1fromL2": ("pocket", "A4L1f"), "X4L2fromL1": ("shallow", "A4L2f")}
    info = {}
    for arm, (start_char, ref) in arms.items():
        a = R.get(arm, {})
        if a.get("state", "not run") == "not run" or a.get("q_final") is None or a.get("iters", 0) < 1:
            info[arm] = None
            continue
        ch = conc(a)
        info[arm] = dict(char=ch, keeps=(ch == start_char), conv=conv(arm), strict=strict(a), ref=ref,
                         q_final=a["q_final"], ynorm=a["ynorm"])
    if any(v is None for v in info.values()):
        lab = "inconclusive"
    else:
        lab = None
        for arm, v in info.items():
            ref = R.get(v["ref"], {})
            if v["keeps"] and v["conv"] and v["strict"][0]:
                if ref.get("state", "not run") == "not run":
                    continue
                if (conv(v["ref"]) and strict(ref)[0] and conc(ref) != v["char"]
                        and not same_point(R[arm], ref, root)):
                    lab = "multiple minima at λ_k" + (" (marginal)" if (v["strict"][1] or strict(ref)[1]) else "")
                    break
        if lab is None:
            if any(v["keeps"] for v in info.values()):
                lab = "path-dependent (trajectory)"
                det = [a for a, v in info.items() if v["keeps"]]
                if any(R.get(info[a]["ref"], {}).get("state", "not run") == "not run" for a in det):
                    lab += " (reference not run)"
            else:
                lab = "budget-dependent"
                det = list(info)
            if any(not info[a]["conv"] for a in det):
                lab += " (not converged)"
    out["D6x"] = dict(label=lab, arms=info)
    json.dump(out, open(root / "results" / "v4_decisions.json", "w", encoding="utf-8"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    root = Path(sys.argv[1])
    cmd = sys.argv[2]
    {"eps": cmd_eps, "nv6prep": cmd_nv6prep, "nv6eval": cmd_nv6eval}.get(cmd, lambda r, a: cmd_decide(r))(root, sys.argv[3:])
