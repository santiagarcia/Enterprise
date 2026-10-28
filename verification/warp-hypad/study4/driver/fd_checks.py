"""Gates H2/H3 (PROTOCOL-v3 2.2): FD direction files and analysis against the OTI Hessian.

Usage:
  python -I fd_checks.py STUDY4_DIR prepare     # writes hessians/fd/<design>.dirs
  python -I fd_checks.py STUDY4_DIR analyse     # reads hessians/fd/<design>.out, writes hessians/fd/fd_checks.csv
Readings (logged before data): slope rule primary = at least one consecutive step pair on the truncation side with
e(h_k)/e(h_{k+1}) >= 30; strict = two consecutive such pairs. H3 step scale h*max(1,|v.y|).
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

N = 108
HS = [1e-2, 1e-3, 1e-4, 1e-5]
H2 = ["A_4_L2", "A_4_L1"]
H3 = ["A_4_L1", "A_4_L2", "B_4_L2", "A_6_L2", "B_6_L2"]


def Hof(root, d):
    return np.fromfile(root / "hessians" / d / "H_R_fwd.bin", dtype="<f8").reshape(N, N)


def yof(root, d):
    return np.fromfile(root / "hessians" / d / "y.bin", dtype="<f8")


def vdir():
    v = np.random.default_rng(4).standard_normal(N)
    return v / np.linalg.norm(v)


def h2_cols(H):
    j = int(np.argmax(np.abs(np.diag(H))))
    cols = [0, 36, 72]
    return cols + ([j] if j not in cols else []), j


def prepare(root):
    fd = root / "hessians" / "fd"; fd.mkdir(parents=True, exist_ok=True)
    for d in sorted(set(H2) | set(H3)):
        y = yof(root, d); lines = []
        if d in H2:
            cols, _ = h2_cols(Hof(root, d))
            for i in cols:
                e = np.zeros(N); e[i] = 1.0
                for h in HS:
                    lines.append(f"col{i} {h * max(1.0, abs(y[i])):.17g} " + " ".join(repr(float(x)) for x in e))
        if d in H3:
            v = vdir(); s = max(1.0, abs(float(v @ y)))
            for h in HS:
                lines.append(f"v4 {h * s:.17g} " + " ".join(repr(float(x)) for x in v))
        (fd / f"{d}.dirs").write_text("\n".join(lines) + "\n")
        print(d, len(lines), "lines")


def slope_rule(e):
    """e over HS (decreasing h). Truncation side: up to the minimum. Returns (primary, strict, best)."""
    k = int(np.argmin(e)); ratios = [e[i] / e[i + 1] for i in range(k)]
    ok = [r >= 30 for r in ratios]
    strict = any(ok[i] and ok[i + 1] for i in range(len(ok) - 1))
    return any(ok), strict, float(e[k]), ratios


def analyse(root):
    fd = root / "hessians" / "fd"; rows = []
    for d in sorted(set(H2) | set(H3)):
        H = Hof(root, d)
        lines = (fd / f"{d}.out").read_text().splitlines()
        R0 = float(lines[0].split()[1])
        data = {}
        for ln in lines[1:]:
            t = ln.split(); lab, h = t[0], float(t[1]); Rp, Rm = float(t[2]), float(t[3])
            gp = np.array(t[4:4 + N], float); gm = np.array(t[4 + N:4 + 2 * N], float)
            data.setdefault(lab, []).append((h, Rp, Rm, gp, gm))
        fdcols = {}
        for lab, L in data.items():
            if lab.startswith("col"):
                i = int(lab[3:]); ref = H[:, i]; es = []
                for h, Rp, Rm, gp, gm in L:
                    c = (gp - gm) / (2 * h); es.append(np.linalg.norm(c - ref) / np.linalg.norm(ref)); fdcols.setdefault(i, []).append(c)
                p, s, best, rat = slope_rule(np.array(es))
                rows.append(dict(gate="H2", design=d, dir=lab, e=";".join(f"{x:.3e}" for x in es),
                                 ratios=";".join(f"{x:.1f}" for x in rat), primary=p, strict=s, best=best, best_ok=best <= 1e-4))
            else:
                v = vdir(); Hv = H @ v; vHv = float(v @ Hv); es, es2 = [], []
                for h, Rp, Rm, gp, gm in L:
                    es.append(np.linalg.norm((gp - gm) / (2 * h) - Hv) / np.linalg.norm(Hv))
                    es2.append(abs((Rp - 2 * R0 + Rm) / h ** 2 - vHv) / abs(vHv))
                for nm, E in (("Hv", es), ("vHv", es2)):
                    p, s, best, rat = slope_rule(np.array(E))
                    rows.append(dict(gate=f"H3-{nm}", design=d, dir=lab, e=";".join(f"{x:.3e}" for x in E),
                                     ratios=";".join(f"{x:.1f}" for x in rat), primary=p, strict=s, best=best, best_ok=best <= 1e-4))
        # FD symmetry residual (descriptive), at best step per column
        ks = sorted(fdcols)
        for a in ks:
            for b in ks:
                if a < b:
                    rs = [abs(fdcols[a][j][b] - fdcols[b][j][a]) / abs(H[a, b]) if H[a, b] else float("nan") for j in range(len(HS))]
                    rows.append(dict(gate="H2-sym", design=d, dir=f"{a},{b}", e=";".join(f"{x:.3e}" for x in rs), ratios="",
                                     primary="", strict="", best=min(rs), best_ok=""))
    with open(fd / "fd_checks.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); [w.writerow(r) for r in rows]
    for r in rows:
        print(r)
    gated = [r for r in rows if r["gate"] in ("H2", "H3-Hv", "H3-vHv")]
    print("H2/H3 primary reading:", "PASS" if all(r["primary"] and r["best_ok"] for r in gated) else "FAIL")
    print("H2/H3 strict reading :", "PASS" if all(r["strict"] and r["best_ok"] for r in gated) else "FAIL")


if __name__ == "__main__":
    root = Path(sys.argv[1])
    prepare(root) if sys.argv[2] == "prepare" else analyse(root)
