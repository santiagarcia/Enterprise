"""Summarise wh_probe fd output: OTI vs central FD relative error across the step sweep.

Usage: python analyze_fd.py fd_all.csv > fd_summary.md
Columns of fd_all.csv: tag,objective,kind,i,j,h,oti,fd,relerr,argmin_plus,argmin_minus
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict


def main(path: str) -> None:
    """Print a markdown table: per check, the relative error at each h and the best h."""
    rows: dict[tuple[str, str, str, str, str], list[tuple[float, float, float, float, str, str]]] = defaultdict(list)
    with open(path, newline="") as fh:
        for r in csv.reader(fh):
            if len(r) < 9 or r[2] not in ("grad", "hess", "fwd", "bwd"):
                continue
            tag, obj, kind, i, j, h, oti, fd, rel = r[:9]
            ap = r[9] if len(r) > 9 else ""
            am = r[10] if len(r) > 10 else ""
            o, f = float(oti), float(fd)
            relv = abs(f - o) / max(abs(o), 1e-300)
            rows[(tag, obj, kind, i, j)].append((float(h), o, f, relv, ap, am))
    hs_all = sorted({h for v in rows.values() for (h, *_rest) in v}, reverse=True)
    print("| tag | J | kind | dir | " + " | ".join(f"h={h:g}" for h in hs_all) + " | best rel | best h | OTI value |")
    print("|" + "---|" * (len(hs_all) + 7))
    for key in sorted(rows):
        tag, obj, kind, i, j = key
        vals = {h: rel for (h, _o, _f, rel, _a, _b) in rows[key]}
        best_h = min(vals, key=vals.get)
        oti = rows[key][0][1]
        d = f"{i}{j}" if j else i
        cells = " | ".join(f"{vals[h]:.1e}" if h in vals else "" for h in hs_all)
        switch = any(a and b and a != b for (*_x, a, b) in rows[key])
        flag = " (argmin switches)" if switch and obj == "rho_min" else ""
        print(f"| {tag} | {obj} | {kind}{flag} | {d} | {cells} | {vals[best_h]:.1e} | {best_h:g} | {oti:.6e} |")


if __name__ == "__main__":
    main(sys.argv[1])
