"""Protocol-v2 retrospective check of the replacement D4 statistic (Noether R5) on ARCHIVED study-3 data only.

u_w = (rho_min,w - rho*_D) / |rho*_D|, with rho_min,w = last trace row `rhomin_w_post` and rho*_D the assessed
worst value (study3/assess/assess.csv, column 4, first seed). Evaluates no objective; uses no study-4 data.
The 10 % threshold is the study-3 value, fixed before this script was run.
Usage: python -I d4_retro_check.py STUDY3_DIR
"""
import csv
import sys

B_REF = -9.60726526e-4


def main() -> None:
    root = sys.argv[1]
    rs: dict[str, float] = {}
    with open(f"{root}/assess/assess.csv") as fh:
        for r in csv.reader(fh):
            rs.setdefault(r[0], float(r[3]))
    for run in ("A_4_L1", "A_6_L1", "B_4_L1", "B_6_L1", "A_4_L2", "A_6_L2", "B_4_L2", "B_6_L2"):
        with open(f"{root}/runs/{run}/trace.csv") as fh:
            t = list(csv.DictReader(fh))[-1]
        w, d = float(t["rhomin_w_post"]), rs[run]
        u = (w - d) / abs(d)
        print(f"{run} rho*_D={d:.4e} rhomin_w={w:.4e} q={d / B_REF:.3f} u_w={100 * u:.1f}% "
              f"{'FLAG (>10%)' if u > 0.10 else 'no flag'}")


if __name__ == "__main__":
    main()
