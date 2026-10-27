"""Count verification cases per gate from build-*/results/verification/*/verification_summary.json."""
from __future__ import annotations
import glob, json, os, sys

def main(root: str) -> None:
    """Print per-gate pass/fail/skip counts and README-group totals."""
    tot = {}
    for p in sorted(glob.glob(os.path.join(root, "*", "verification_summary.json"))):
        d = json.load(open(p, encoding="utf-8"))
        cs = d.get("cases", [])
        st = ["SKIP" if c.get("skipped") else ("PASS" if c.get("passed") else "FAIL") for c in cs]
        n = len(cs); ps = sum(s.startswith("PASS") for s in st); sk = sum(s.startswith("SKIP") for s in st)
        bl = sum(s.startswith("BLOCK") for s in st); fl = n - ps - sk - bl
        g = os.path.basename(os.path.dirname(p)); tot[g] = (n, ps, fl, sk, bl)
        print(f"{g:28s} cases {n:4d} pass {ps:4d} fail {fl:3d} skip {sk:3d} blocked {bl:3d} all_passed={d.get('all_passed')}")
    a = ["milestone_01","milestone_02","milestone_02t","milestone_03","milestone_04","milestone_05","otilib","program_b1"]
    b = [g for g in tot if g.startswith("program_b")]
    for name, grp in (("README-285 group", a), ("Program B group", b)):
        s = [sum(tot[g][i] for g in grp if g in tot) for i in range(5)]
        print(f"{name}: gates {[g for g in grp if g in tot]} -> cases {s[0]} pass {s[1]} fail {s[2]} skip {s[3]} blocked {s[4]}")
    s = [sum(v[i] for v in tot.values()) for i in range(5)]
    print(f"ALL: gates {len(tot)} cases {s[0]} pass {s[1]} fail {s[2]} skip {s[3]} blocked {s[4]}")

if __name__ == "__main__":
    main(sys.argv[1])
