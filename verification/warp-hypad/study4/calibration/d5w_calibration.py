"""Pre-registration calibration of the study-4 window rule D5w on ARCHIVED study-3 traces only.

Reads study3/runs/*/trace.csv (in-loop q); evaluates no objective. Applies the candidate rule with
window width W from a start iterate s0 over the common iterate horizon of each A/B pair and prints the
per-window median differences and the verdict. Usage: python -I d5w_calibration.py STUDY3_DIR
"""
import csv
import statistics as st
import sys

B_REF = -9.60726526e-4
E0 = 4.634569e-2
THR = 0.05


def load(root: str, run: str) -> dict[int, dict[str, float]]:
    out = {}
    with open(f"{root}/runs/{run}/trace.csv") as fh:
        for r in csv.DictReader(fh):
            out[int(r["it"])] = {"q": float(r["rhostar_in"]) / B_REF, "rE": 1 - float(r["Eneg_v8"]) / E0}
    return out


def plateau(a: dict[int, dict[str, float]], k: int) -> int:
    """First iterate from which A's r_E stays within a 0.01 band for 10 iterates."""
    rE = [a[i]["rE"] for i in range(k + 1)]
    return next(i for i in range(k + 1) if max(rE[i:i + 10]) - min(rE[i:i + 10]) <= 0.01)


def verdict(diffs: list[float], nmin: int = 2) -> str:
    n = len(diffs)
    if n < nmin:
        return f"not evaluable (<{nmin} windows)"
    pos = sum(d >= THR for d in diffs)
    neg = sum(d <= -THR for d in diffs)
    if pos * 3 >= 2 * n and neg == 0:
        return "B lower (suppresses)"
    if neg * 3 >= 2 * n and pos == 0:
        return "A lower (reversed)"
    return "not established"


def main() -> None:
    root = sys.argv[1]
    for W in (25, 20):
        for mode in ("plateau", "fixed40"):
            print(f"== W={W}, start={mode}")
            for lam in ("L1", "L2"):
                per_d = []
                for d in ("4", "6"):
                    a, b = load(root, f"A_{d}_{lam}"), load(root, f"B_{d}_{lam}")
                    k = min(max(a), max(b))
                    s0 = plateau(a, k) if mode == "plateau" else 40
                    diffs, txt = [], []
                    s = s0
                    while s + W - 1 <= k:
                        qa = [a[i]["q"] for i in range(s, s + W)]
                        qb = [b[i]["q"] for i in range(s, s + W)]
                        diffs.append(st.median(qa) - st.median(qb))
                        txt.append(f"[{s}-{s + W - 1}] {st.median(qa):.3f}/{st.median(qb):.3f} ({diffs[-1]:+.3f})")
                        s += W
                    v = verdict(diffs)
                    per_d.append(v)
                    print(f"  {lam} {d}:1 horizon {k} start {s0}: {v}; " + "; ".join(txt))
                both = per_d[0] if per_d[0] == per_d[1] else "not established (densities disagree)"
                print(f"  {lam} lambda-level verdict: {both}")


def registered(root: str) -> None:
    """Final registered form: >= 4 windows, start at the plateau onset, verdict must agree at W = 20 and 25."""
    print("== REGISTERED FORM (nmin 4, W 20 and 25 must agree), start = plateau onset")
    for lam in ("L1", "L2"):
        per_d = []
        for d in ("4", "6"):
            a, b = load(root, f"A_{d}_{lam}"), load(root, f"B_{d}_{lam}")
            k = min(max(a), max(b))
            s0 = plateau(a, k)
            vs = []
            for W in (20, 25):
                diffs = []
                s = s0
                while s + W - 1 <= k:
                    diffs.append(st.median([a[i]["q"] for i in range(s, s + W)])
                                 - st.median([b[i]["q"] for i in range(s, s + W)]))
                    s += W
                vs.append(verdict(diffs, 4))
            v = vs[0] if vs[0] == vs[1] else "not established (window-sensitive)"
            per_d.append(v)
            print(f"  {lam} {d}:1 W20: {vs[0]} | W25: {vs[1]} -> {v}")
        print(f"  {lam} lambda-level: {per_d[0] if per_d[0] == per_d[1] else 'not established (densities disagree)'}")


if __name__ == "__main__":
    main()
    registered(sys.argv[1])
