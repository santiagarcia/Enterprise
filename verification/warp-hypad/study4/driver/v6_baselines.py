"""Gate V6 (PROTOCOL-v3 sec 9): recompute registered gap baselines from archived study-3 traces.

Usage: python -I v6_baselines.py WARP_HYPAD_ROOT
q = rhostar_in / b, b = -9.60726526e-4. Iterate k is the trace row with it == k.
"""
import csv, statistics, sys
from pathlib import Path

B = -9.60726526e-4
REG = {4: dict(s0=84, d0=1.617, pre20=2.287, pre25=2.273), 6: dict(s0=90, d0=1.146, pre20=-0.093, pre25=-0.084)}


def q_of(path):
    with open(path) as f:
        return {int(r["it"]): float(r["rhostar_in"]) / B for r in csv.DictReader(f)}


def main(root):
    ok = True
    s3 = Path(root) / "study3" / "runs"
    for d, reg in REG.items():
        qa, qb = q_of(s3 / f"A_{d}_L2/trace.csv"), q_of(s3 / f"B_{d}_L2/trace.csv")
        s0 = reg["s0"]
        d0 = qa[s0] - qb[s0]
        line = f"{d}:1 s0={s0} qA={qa[s0]:.6f} qB={qb[s0]:.6f} D0={d0:+.6f} (reg {reg['d0']:+.3f})"
        ok &= abs(d0 - reg["d0"]) <= 0.001 + 1e-12
        for W in (20, 25):
            its = range(s0 - W + 1, s0 + 1)
            dp = statistics.median(qa[k] for k in its) - statistics.median(qb[k] for k in its)
            r = reg[f"pre{W}"]
            m = abs(dp - r) <= 0.001 + 1e-12
            ok &= m
            line += f" | Dpre(W={W})={dp:+.6f} (reg {r:+.3f}) {'ok' if m else 'MISMATCH'}"
        print(line)
    qa6 = q_of(s3 / "A_6_L2/trace.csv")
    for k, r in ((87, 0.71), (90, 2.02)):
        m = abs(qa6[k] - r) <= 0.005 + 1e-12
        ok &= m
        print(f"A 6:1 in-loop q at iterate {k} = {qa6[k]:.6f} (reg {r:.2f}) {'ok' if m else 'MISMATCH'}")
    print("V6", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
