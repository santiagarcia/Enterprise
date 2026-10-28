"""Protocol-v2 calibration of D5w (Noether S5, S4) on ARCHIVED study-3 traces only.

Reads study3/runs/*/trace.csv (in-loop q); evaluates no objective and uses no study-4 data.
Applies the v2 registered rule with ITERATE-1 ALIGNMENT: a pseudo study-4 start s0 is chosen on the
study-3 lambda1 A/B pair (where A and B are known to be indistinguishable late), windows are
[s0+1, s0+W], [s0+W+1, s0+2W], ... up to the common horizon, complete windows only. Also reports the
v2 change-based companion G_w = Delta_w - Delta_0 with two baselines: "point" (Delta_0 = q_A - q_B at
iterate s0 only) and "pre-window" (registered in v2: Delta_0 = median q_A - median q_B over the W
iterates s0-W+1..s0 of the archived history; at s0 < W, from iterate 0).
A "B lower"/"A lower" verdict or a "gap grows" verdict here is a false positive of the rule.
Usage: python -I d5w_calibration_v2.py STUDY3_DIR
"""
import csv
import statistics as st
import sys

B_REF = -9.60726526e-4
THR = 0.05


def load(root: str, run: str) -> dict[int, float]:
    with open(f"{root}/runs/{run}/trace.csv") as fh:
        return {int(r["it"]): float(r["rhostar_in"]) / B_REF for r in csv.DictReader(fh)}


def verdict(diffs: list[float], pos_lbl: str, neg_lbl: str) -> str:
    n = len(diffs)
    if n < 4:
        return "not evaluable (<4 windows)"
    pos = sum(d >= THR for d in diffs)
    neg = sum(d <= -THR for d in diffs)
    if pos * 3 >= 2 * n and neg == 0:
        return pos_lbl
    if neg * 3 >= 2 * n and pos == 0:
        return neg_lbl
    return "not established"


def windows(a: dict[int, float], b: dict[int, float], s0: int, k: int, W: int) -> list[float]:
    out, s = [], s0 + 1
    while s + W - 1 <= k:
        out.append(st.median([a[i] for i in range(s, s + W)]) - st.median([b[i] for i in range(s, s + W)]))
        s += W
    return out


def main() -> None:
    root = sys.argv[1]
    n_rule = n_gap = n_pgap = n_cases = 0
    for d, registered_s0 in (("4", 84), ("6", 90)):
        a, b = load(root, f"A_{d}_L1"), load(root, f"B_{d}_L1")
        k = min(max(a), max(b))
        for s0 in sorted({registered_s0, 0, 20, 40, 60, 80, 100}):
            d0 = a[s0] - b[s0]
            res, gres, pres = [], [], []
            for W in (20, 25):
                diffs = windows(a, b, s0, k, W)
                lo = max(0, s0 - W + 1)
                p0 = st.median([a[i] for i in range(lo, s0 + 1)]) - st.median([b[i] for i in range(lo, s0 + 1)])
                res.append((verdict(diffs, "B lower", "A lower"), diffs))
                gres.append(verdict([x - d0 for x in diffs], "gap grows", "gap shrinks"))
                pres.append(verdict([x - p0 for x in diffs], "gap grows", "gap shrinks"))
            v = res[0][0] if res[0][0] == res[1][0] else "not established (window-sensitive)"
            g = gres[0] if gres[0] == gres[1] else "not established (window-sensitive)"
            pg = pres[0] if pres[0] == pres[1] else "not established (window-sensitive)"
            tag = " <- registered alignment (B's lambda2 count)" if s0 == registered_s0 else ""
            print(f"L1 {d}:1 s0={s0:3d} horizon={k} Delta_0(point)={d0:+.3f}: D5w {v}; gap(point) {g}; "
                  f"gap(pre-window) {pg}{tag}")
            for W, (vv, diffs) in zip((20, 25), res):
                print(f"    W={W} n={len(diffs)} Delta_w=" + " ".join(f"{x:+.3f}" for x in diffs))
            n_cases += 1
            n_rule += v in ("B lower", "A lower")
            n_gap += g.startswith("gap")
            n_pgap += pg.startswith("gap")
    print(f"SUMMARY: {n_cases} alignments; D5w false positives {n_rule}; gap companion false positives: "
          f"point baseline {n_gap}, pre-window baseline {n_pgap}")


if __name__ == "__main__":
    main()
