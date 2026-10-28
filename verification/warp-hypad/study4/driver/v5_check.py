"""Gate V5e (i) exact check and V5e (ii) norm check (PROTOCOL-v3 section 9).

Usage: python -I v5_check.py STUDY4_DIR
"""
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

d = Path(sys.argv[1]) / "results" / "v5"
lines = (d / "v5e_terms.txt").read_text().split("\n")
neu = float.fromhex(lines[0].split()[1]); fwd = float.fromhex(lines[1].split()[1])
xs = [float.fromhex(s) for s in lines[3:] if s.strip()]
exact = sum((Fraction(x) for x in xs), Fraction(0))
sabs = sum(abs(x) for x in xs)
u = 2.0 ** -53
err_n = abs(Fraction(neu) - exact); err_f = abs(Fraction(fwd) - exact)
ok1 = err_n <= Fraction(2 * u * sabs)
print(f"V5e(i) n={len(xs)} exact={float(exact):.17g} neumaier={neu:.17g} |err|={float(err_n):.3e} "
      f"bound 2u*sum|x|={2*u*sabs:.3e} {'PASS' if ok1 else 'FAIL'} (forward |err|={float(err_f):.3e})")
n = 108
rd = lambda f: np.fromfile(d / "v5c" / f, dtype="<f8").reshape(n, n)
Hf, Hc, Hr, Hl, B, Hx = (rd(f) for f in ("H_R_fwd.bin", "H_R_comp.bin", "H_R_rev.bin", "H_R_ld.bin", "B_H.bin", "H_exact.bin"))
ehc = np.linalg.norm(Hf - Hc, 2); bF = np.linalg.norm(B, "fro")
ok2 = ehc <= bF
print(f"V5e(ii) e_H^c={ehc:.3e} <= |B_H|_F={bF:.3e} {'PASS' if ok2 else 'FAIL'}; e_H(rev)={np.linalg.norm(Hf-Hr,2):.3e}; "
      f"|H_fwd-H_ld|_2={np.linalg.norm(Hf-Hl,2):.3e}; |H_comp-H_exact|_2={np.linalg.norm(Hc-Hx,2):.3e}")
print("V5e-python", "PASS" if ok1 and ok2 else "FAIL")
sys.exit(0 if ok1 and ok2 else 1)
