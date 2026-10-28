"""PROTOCOL-v4 s7 (v2, S2): write y + Newton step (-H_F^{-1} grad F) for a converged strict final.

Usage: python -I s2_newton.py HESSDIR OUTFILE   (H_F = H_R_fwd + lambda I; grad F = gR_fwd + lambda y)
"""
import sys
from pathlib import Path

import numpy as np

d = Path(sys.argv[1]); N = 108
lam = float(next(l.split()[1] for l in (d / "meta.txt").read_text().splitlines() if l.startswith("lambda ")))
H = np.fromfile(d / "H_R_fwd.bin", dtype="<f8").reshape(N, N) + lam * np.eye(N)
y = np.fromfile(d / "y.bin", dtype="<f8")
g = np.fromfile(d / "gR_fwd.bin", dtype="<f8") + lam * y
s = -np.linalg.solve(H, g)
np.savetxt(sys.argv[2], y + s, fmt="%.17g")
print(d.name, "newton_step_norm", float(np.linalg.norm(s)), "rel", float(np.linalg.norm(s) / np.linalg.norm(y)))
