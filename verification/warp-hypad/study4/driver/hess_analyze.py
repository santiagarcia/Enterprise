"""Study-4 Hessian analysis (PROTOCOL-v3 sections 2.2 H1a-H1c, 2.3): spectra, error measures, labels.

Usage:
  python -I hess_analyze.py STUDY4_DIR design DIR [DIR ...]   # analyse Hessian directories (writes eig.csv, summary.csv)
  python -I hess_analyze.py STUDY4_DIR step0                  # the 10 step-0 designs + checks.csv + labels
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

N, BLK = 108, 36
FIELDS = ("lapse", "axial", "radial")
L1, L2 = 3.75e-8, 3.75e-9
STEP0 = ["A_4_L1", "A_4_L2", "A_6_L1", "A_6_L2", "B_4_L1", "B_4_L2", "B_6_L1", "B_6_L2", "A_4_L2_it84", "A_6_L2_it90"]


def meta(d: Path) -> dict:
    m = {}
    for line in (d / "meta.txt").read_text().splitlines():
        k, _, v = line.partition(" ")
        m[k] = v
    return m


def mat(d: Path, f: str) -> np.ndarray:
    return np.fromfile(d / f, dtype="<f8").reshape(N, N)


def vec(d: Path, f: str) -> np.ndarray:
    return np.fromfile(d / f, dtype="<f8")


def analyse(d: Path, lam: float | None = None) -> dict:
    """All section-2.3 quantities for one Hessian directory. lam overrides meta lambda (cross-lambda reuse)."""
    m = meta(d)
    lam = float(m["lambda"]) if lam is None else lam
    Hf, Hr, Hc, Hl, B = (mat(d, f) for f in ("H_R_fwd.bin", "H_R_rev.bin", "H_R_comp.bin", "H_R_ld.bin", "B_H.bin"))
    y = vec(d, "y.bin")
    g = {k: vec(d, f"gR_{k}.bin") for k in ("fwd", "rev", "comp", "ld", "bound")}
    gF = {k: g[k] + lam * y for k in ("fwd", "rev", "comp", "ld")}
    pen = 0.5 * lam * float(y @ y)
    R = {k: float(m[f"R_{k}"]) for k in ("fwd", "rev", "comp", "ld", "bound")}
    F = {k: R[k] + pen for k in ("fwd", "rev", "comp", "ld")}
    HF = Hf + lam * np.eye(N)
    mu, Q = np.linalg.eigh(HF)
    muR = np.linalg.eigvalsh(Hf)
    e_H = np.linalg.norm(Hf - Hr, 2); e_Hc = np.linalg.norm(Hf - Hc, 2); e_Hld = np.linalg.norm(Hf - Hl, 2)
    BF = np.linalg.norm(B, "fro")
    eh = max(e_H, e_Hc)
    e_g = np.linalg.norm(gF["fwd"] - gF["rev"]); e_gc = np.linalg.norm(gF["fwd"] - gF["comp"]); eg = max(e_g, e_gc)
    e_F = abs(F["fwd"] - F["rev"]); e_Fc = abs(F["fwd"] - F["comp"]); eF = max(e_F, e_Fc)
    ynorm = float(np.linalg.norm(y)); gn = float(np.linalg.norm(gF["fwd"]))
    tol = 1e-5 * F["fwd"] / max(ynorm, 1.0)
    mumax = float(np.max(np.abs(mu)))
    tau = max(1e-9 * mumax, 10 * eh)
    neg = mu < -tau; nz = np.abs(mu) <= tau; pos = mu > tau
    marg = (np.abs(mu) > tau) & (np.abs(mu) < 100 * tau)
    pd = bool(np.all(pos))
    kappa = float(mu[-1] / mu[0]) if pd else float(mumax / np.min(np.abs(mu)))
    shares = (Q.reshape(3, BLK, N) ** 2).sum(axis=1)  # [field, eig]
    out = dict(
        design=m.get("tag", d.name), dir=str(d), nodes=int(m["nodes"]), leaves=int(m["leaves"]), lam=lam,
        ynorm=ynorm, F=F["fwd"], R=R["fwd"],
        H1a_diag_mismatch=int(m["H1a_diag_mismatch"]), H1a_passgrad_mismatch=int(m["H1a_passgrad_mismatch"]),
        H1a_grad_vs_gradient_mismatch=int(m["H1a_grad_vs_gradient_mismatch"]), pair_errors=int(m["pair_errors"]),
        e_H=e_H, e_H_c=e_Hc, ehat_H=eh, e_H_ld=e_Hld, B_H_F=BF, ratio_eHc_eH=(e_Hc / e_H if e_H > 0 else float("inf")),
        reverse_understates=bool(e_Hc > 10 * e_H),
        cons1_eHc_le_BF=bool(e_Hc <= BF),
        cons2_ld=bool(abs(e_Hc - e_Hld) <= 0.1 * e_Hc + 1e-2 * BF),
        e_g=e_g, e_g_c=e_gc, ehat_g=eg, b_g=float(np.linalg.norm(g["bound"])),
        e_F=e_F, e_F_c=e_Fc, ehat_F=eF, b_F=R["bound"],
        tolC1=tol, ehat_g_over_tol=eg / tol, gF=gn, gF_over_tol=gn / tol, gF_over_lam_y=gn / (lam * ynorm),
        mu_min=float(mu[0]), mu_max=float(mu[-1]), max_abs_mu=mumax, tau_H=tau,
        H1b_gate_fail=bool(10 * eh > 1e-6 * mumax),
        n_neg=int(neg.sum()), n_nearzero=int(nz.sum()), n_pos=int(pos.sum()),
        n_marginal=int(marg.sum()), marginal_signs="".join("+" if v > 0 else "-" for v in mu[marg]),
        mu_min_over_tau=float(mu[0] / tau),
        label=("PD" if pd else (f"indefinite ({int(neg.sum())} negative)" if neg.any() else "singular to tau_H")),
        kappa=kappa, kappa_kind=("kappa(H_F)" if pd else "kappa_abs (indefinite)" if neg.any() else "kappa_abs"),
        lmin_HR_over_lam=float(muR[0] / lam), mumax_over_lam=float(mu[-1] / lam),
        n_W0_modified=int(np.sum(~(mu >= lam))),
        secs=float(m["secs"]), threads=int(m["threads"]), ld_digits=int(m["ld_digits"]),
    )
    for fi, fname in enumerate(FIELDS):
        out[f"share_min_eig_{fname}"] = float(shares[fi, 0])
        out[f"share_neg_{fname}"] = float(shares[fi, neg].sum() / max(neg.sum(), 1)) if neg.any() else float("nan")
        out[f"share_nz_{fname}"] = float(shares[fi, nz].sum() / max(nz.sum(), 1)) if nz.any() else float("nan")
    if pd:
        step = np.linalg.solve(HF, gF["fwd"])
        out["newton_step"] = float(np.linalg.norm(step))
        out["decrement_over_F"] = float(0.5 * gF["fwd"] @ step / F["fwd"])
    else:
        out["newton_step"] = float("nan"); out["decrement_over_F"] = float("nan")
    with open(d / "eig.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["i", "mu_HF", "mu_HR", "class", "share_lapse", "share_axial", "share_radial"])
        for i in range(N):
            cl = "neg" if neg[i] else "nearzero" if nz[i] else ("pos_marginal" if marg[i] else "pos")
            w.writerow([i, repr(float(mu[i])), repr(float(muR[i])), cl, *(repr(float(shares[k, i])) for k in range(3))])
    with open(d / "summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(out.keys()); w.writerow(out.values())
    return out


def h4(root: Path, design: str, r: dict) -> dict:
    tr = root.parent / "study3" / "runs" / design / "trace.csv"
    rows = list(csv.DictReader(open(tr)))
    gam = [float(x["gamma"]) for x in rows if x["gamma"] not in ("", "nan")]
    gl = gam[-1] if gam else float("nan")
    inv = 1.0 / gl if gl else float("nan")
    return dict(design=design, gamma_last=gl, inv_gamma=inv, mu_min=r["mu_min"], mu_max=r["mu_max"],
                H4_in_range=bool(r["mu_min"] <= inv <= r["mu_max"]))


def step0(root: Path) -> None:
    res = []
    for dsg in STEP0:
        res.append(analyse(root / "hessians" / dsg))
    keys = list(res[0].keys())
    with open(root / "hessians" / "checks.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); [w.writerow(r) for r in res]
    # labels (Question 1)
    lab = {}
    for lname, lam, finals in (("lambda1", L1, ["A_4_L1", "A_6_L1", "B_4_L1", "B_6_L1"]),
                               ("lambda2", L2, ["A_4_L2", "A_6_L2", "B_4_L2", "B_6_L2"])):
        rr = [r for r in res if r["design"] in finals]
        pdk = [r["kappa"] for r in rr if r["label"] == "PD"]
        if not pdk:
            k_lab = "not evaluable (no PD final)"
        elif all(1e3 <= k <= 1e4 for k in pdk):
            k_lab = "consistent"
        else:
            k_lab = "not supported"
        mx = [r["mumax_over_lam"] for r in rr]
        i5 = "not evaluable" if any(r["mu_max"] <= 0 for r in rr) else (
            "consistent" if all(1e2 <= v <= 1e4 for v in mx) else "not supported")
        lab[lname] = dict(kappa_statement=k_lab, kappa_PD_finals=pdk, kappa_all=[(r["design"], r["label"], r["kappa"]) for r in rr],
                          I5_statement=i5, mumax_over_lam=[(r["design"], v) for r, v in zip(rr, mx)])
    h4s = [h4(root, d, next(r for r in res if r["design"] == d)) for d in STEP0[:8]]
    gate = dict(
        H1a_all_zero=all(r["H1a_diag_mismatch"] == 0 and r["H1a_passgrad_mismatch"] == 0 and r["H1a_grad_vs_gradient_mismatch"] == 0 for r in res),
        pair_errors_zero=all(r["pair_errors"] == 0 for r in res),
        H1b_gate_pass=not any(r["H1b_gate_fail"] for r in res),
        cons1_all=all(r["cons1_eHc_le_BF"] for r in res), cons2_all=all(r["cons2_ld"] for r in res),
        reverse_understates=[r["design"] for r in res if r["reverse_understates"]],
    )
    json.dump(dict(gate=gate, labels=lab, H4=h4s), open(root / "results" / "step0_summary.json", "w"), indent=1)
    print(json.dumps(gate, indent=1))
    print(json.dumps(lab, indent=1))
    for h in h4s:
        print("H4", h)
    cols = ["design", "nodes", "label", "kappa", "mu_min", "mu_max", "tau_H", "n_neg", "n_nearzero", "n_marginal",
            "marginal_signs", "lmin_HR_over_lam", "mumax_over_lam", "e_H", "e_H_c", "e_H_ld", "B_H_F", "ehat_g_over_tol",
            "gF_over_tol", "gF_over_lam_y", "n_W0_modified", "secs"]
    print(",".join(cols))
    for r in res:
        print(",".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c]) for c in cols))


if __name__ == "__main__":
    root = Path(sys.argv[1])
    if sys.argv[2] == "step0":
        step0(root)
    else:
        for dd in sys.argv[3:]:
            r = analyse(Path(dd))
            print(json.dumps(r, indent=1, default=str))
