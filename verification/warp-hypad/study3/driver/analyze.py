"""Study-3 analysis (PROTOCOL.md sections 7-8): iterate selection, decision rules, figures.

Run with the Anaconda interpreter in isolated mode:
    C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I analyze.py select   # writes results/assess_list.txt
    C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I analyze.py report   # tables + figures
All statements concern the Eulerian energy density rho only.
"""

from __future__ import annotations

import csv
import glob
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import SymLogNorm  # noqa: E402

ROOT = "C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad/study3"
COLS = ("tag,yfile,seed,rhoAb,xD,rD,sample_min,prepolish,sub_min,min_est,flagged,d7,starts,mism,chk,fallback,"
        "E8,E8W,negfrac,vtail,L4,nodemin8,refined8,E16,m16,R4,R6,R8,llmin,llmax,bmax,nval,nadm,secs").split(",")
RE_LEVELS = (0.40, 0.50, 0.60)
Y_LEVELS = (10.0, 20.0, 40.0)
LAMS = {"L1": 3.75e-8, "L2": 3.75e-9}


def load_csv(path: str) -> dict[tuple[str, int], dict[str, float | str]]:
    """Return assessment rows keyed by (tag, seed); later rows win."""
    out: dict[tuple[str, int], dict[str, float | str]] = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for row in csv.reader(f):
            d: dict[str, float | str] = {}
            for k, v in zip(COLS, row):
                try:
                    d[k] = float(v)
                except ValueError:
                    d[k] = v
            out[(str(d["tag"]), int(d["seed"]))] = d
    return out


def load_trace(run: str) -> dict[str, np.ndarray]:
    """Load a run's per-iterate trace as float column arrays."""
    with open(f"{ROOT}/runs/{run}/trace.csv") as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def termination(run: str) -> str:
    """Return the termination line of a run, or RUNNING."""
    p = f"{ROOT}/runs/{run}/termination.txt"
    return open(p).read().strip() if os.path.exists(p) else "RUNNING"


def floor_value() -> float:
    """Baseline minimum over D minus Ab, parsed from the finals log."""
    with open(f"{ROOT}/logs/assess_finals.log") as f:
        for line in f:
            if line.startswith("FLOOR"):
                return float(line.split(":")[1].split()[0])
    raise RuntimeError("floor not found")


def runs_list() -> list[str]:
    """All run directories, sorted."""
    return sorted(os.path.basename(p) for p in glob.glob(f"{ROOT}/runs/[AB]_*_L*"))


def first_at(arr: np.ndarray, level: float) -> int | None:
    """Index of the first element >= level, or None."""
    idx = np.nonzero(arr >= level)[0]
    return int(idx[0]) if idx.size else None


def baseline_refs() -> tuple[float, float]:
    """Return (b, E0) from the study-3 baseline assessment (min over seeds and the outside floor)."""
    A = load_csv(f"{ROOT}/assess/assess.csv")
    fl = floor_value()
    b = min(float(A[("baseline", s)]["rhoAb"]) for s in (1, 2) if ("baseline", s) in A)
    return min(b, fl), float(A[("baseline", 1)]["E8"])


def select() -> None:
    """Write the iterate list for assessment priorities 2-5."""
    b, E0 = baseline_refs()
    lines: list[str] = []
    for run in runs_list():
        tr = load_trace(run)
        rE = 1 - tr["Eneg_v8"] / E0
        q = tr["rhostar_in"] / b
        its = tr["it"].astype(int)
        for lv in RE_LEVELS:
            k = first_at(rE, lv)
            if k is not None:
                lines.append(f"{run} {its[k]} D5_rE{lv:.2f}")
        for lv in Y_LEVELS:
            k = first_at(tr["ynorm"], lv)
            if k is not None:
                lines.append(f"{run} {its[k]} D5_y{lv:.0f}")
        ok = rE >= 0.10
        if ok.any():
            k = int(np.argmax(np.where(ok, q, -np.inf)))
            lines.append(f"{run} {its[k]} T_maxq")
        for k in its:
            if k > 0 and k % 50 == 0 and k != its[-1]:
                lines.append(f"{run} {k} snap50")
        for k in its:
            if k > 0 and k % 25 == 0 and k % 50 != 0 and k != its[-1]:
                lines.append(f"{run} {k} snap25")
    # priority order: D5, T, snap50, snap25
    pri = {"D5": 0, "T_": 1, "sn": 2}
    lines.sort(key=lambda s: (pri.get(s.split()[2][:2], 3), "snap25" in s))
    with open(f"{ROOT}/results/assess_list.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"{len(lines)} iterates selected")


def label(q: float, rE: float) -> str:
    """Pre-registered per-design label."""
    if rE >= 0.10 and q >= 1.05:
        return "concentration"
    if q >= 1.01:
        return "marginal"
    return "none"


def report() -> None:
    """Decision tables (results/*.csv, results/decisions.json) and figures."""
    A = load_csv(f"{ROOT}/assess/assess.csv")
    It = load_csv(f"{ROOT}/assess/iters.csv")
    fl = floor_value()
    b, E0 = baseline_refs()
    out: dict[str, object] = {"b": b, "E0": E0, "floor_outside_Ab": fl}
    runs = runs_list()

    def assessed(tag: str, src: dict) -> float | None:
        vals = [float(src[(tag, s)]["rhoAb"]) for s in (1, 2) if (tag, s) in src]
        return min(min(vals), fl) if vals else None

    summ: dict[str, dict[str, object]] = {}
    lines = ["run,term,iters,F,R,pen,ynorm,gF,C1tol,q,q_seed1,q_seed2,r_E,r_E_Dw,E8,E16,E8_E16_rel,R4,R6,R8,"
             "check_worst,underreport_frac,underreport,vtail,L4_rel,llmin,llmax,bmax,valid,flagged,d7,"
             "inloop_q_final,T_label,T_count,T_first,T_last,T_maxq,T_maxq_it,T_maxq_assessed,qrange_last20,"
             "resets,splits_total,nodes_final,last_set_change"]
    base8 = A[("baseline", 1)]
    for run in runs:
        tr = load_trace(run)
        term = termination(run)
        tword = term.split()[0]
        rE_tr = 1 - tr["Eneg_v8"] / E0
        q_tr = tr["rhostar_in"] / b
        ok = rE_tr >= 0.10
        conc = ok & (q_tr >= 1.05)
        its = tr["it"].astype(int)
        if ok.any():
            km = int(np.argmax(np.where(ok, q_tr, -np.inf)))
            maxq, maxq_it = float(q_tr[km]), int(its[km])
        else:
            maxq, maxq_it = float(np.nanmax(q_tr)), int(its[int(np.nanargmax(q_tr))])
        tl = ("concentration along the trajectory" if conc.any() else
              "marginal along the trajectory" if maxq >= 1.01 else "none along the trajectory")
        tmax_ass = assessed(f"{run}@{maxq_it}", It)
        last = slice(max(0, len(its) - 20), None)
        qr = (float(q_tr[last].min()), float(q_tr[last].max()))
        rec: dict[str, object] = {
            "term": tword, "iters": int(its[-1]), "T_label": tl, "T_count": int(conc.sum()),
            "T_first": int(its[conc][0]) if conc.any() else None, "T_last": int(its[conc][-1]) if conc.any() else None,
            "T_maxq": maxq, "T_maxq_it": maxq_it,
            "T_maxq_assessed": (tmax_ass / b) if tmax_ass is not None else None,
            "qrange_last20": qr, "ynorm": float(tr["ynorm"][-1]), "F": float(tr["F"][-1]), "R": float(tr["R"][-1]),
            "pen": float(tr["pen"][-1]),
        }
        if (run, 1) in A:
            a = A[(run, 1)]
            rho = assessed(run, A)
            q = rho / b
            rE = 1 - float(a["E8"]) / E0
            rEw = 1 - float(a["E8W"]) / float(base8["E8W"])
            chk = min(tr["rhomin_w_post"][-1], tr["rhomin_v8"][-1], tr["refined_v8"][-1])
            uf = (chk - rho) / abs(rho)
            rec.update({"q": q, "rE": rE, "rEw": rEw, "underreport_frac": uf, "under": uf > 0.10,
                        "E8": float(a["E8"]), "E16": float(a["E16"]), "R4": float(a["R4"]), "R6": float(a["R6"]),
                        "R8": float(a["R8"]), "label_final": label(q, rE)})
            if tl == "concentration along the trajectory" and rec["label_final"] != "concentration":
                rec["T_label"] = "concentration along the trajectory, not at the final iterate" + (
                    " (converged)" if tword == "CONVERGED" else "")
            qs = [float(A[(run, s)]["rhoAb"]) / b if (run, s) in A else float("nan") for s in (1, 2)]
            lines.append(
                f"{run},{tword},{int(its[-1])},{tr['F'][-1]:.10e},{tr['R'][-1]:.8e},{tr['pen'][-1]:.4e},"
                f"{tr['ynorm'][-1]:.4f},{tr['gF'][-1]:.3e},{1e-5 * tr['F'][-1] / max(tr['ynorm'][-1], 1):.3e},"
                f"{q:.4f},{qs[0]:.4f},{qs[1]:.4f},{rE:.4f},{rEw:.4f},{float(a['E8']):.6e},{float(a['E16']):.6e},"
                f"{(float(a['E8']) - float(a['E16'])) / float(a['E16']):.4f},{float(a['R4']):.8e},"
                f"{float(a['R6']):.8e},{float(a['R8']):.8e},{chk:.6e},{uf:.4f},{int(uf > 0.10)},"
                f"{float(a['vtail']):.4e},{float(a['L4']) / float(base8['L4']):.4f},{float(a['llmin']):.4f},"
                f"{float(a['llmax']):.4f},{float(a['bmax']):.4f},{int(float(a['nval']))}/{int(float(a['nadm']))},"
                f"{int(float(a['flagged']))},{int(float(a['d7']))},{q_tr[-1]:.4f},\"{rec['T_label']}\","
                f"{rec['T_count']},{rec['T_first']},{rec['T_last']},{maxq:.4f},{maxq_it},"
                f"{rec['T_maxq_assessed'] if rec['T_maxq_assessed'] is None else round(rec['T_maxq_assessed'], 4)},"
                f"{qr[0]:.3f}-{qr[1]:.3f},{term.split('resets=')[1].split()[0] if 'resets=' in term else ''},"
                f"{int(tr['splits'].sum())},{int(tr['nodes'][-1])},"
                f"{term.split('last_set_change=')[1].split()[0] if 'last_set_change=' in term else ''}")
        summ[run] = rec
    with open(f"{ROOT}/results/summary.csv", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    # D5 matched progress
    d5_rows = ["lam,density,measure,level,it_A,it_B,qA_assessed,qB_assessed,qA_inloop,qB_inloop,diff_assessed"]
    d5: dict[str, dict[str, object]] = {}
    for L in LAMS:
        verdict_d: dict[str, object] = {}
        for d in ("4", "6"):
            ra, rb = f"A_{d}_{L}", f"B_{d}_{L}"
            if ra not in summ or rb not in summ:
                continue
            ta, tb = load_trace(ra), load_trace(rb)
            hits = 0
            reached = 0
            for meas, levels in (("rE", RE_LEVELS), ("y", Y_LEVELS)):
                for lv in levels:
                    if meas == "rE":
                        ka, kb = first_at(1 - ta["Eneg_v8"] / E0, lv), first_at(1 - tb["Eneg_v8"] / E0, lv)
                    else:
                        ka, kb = first_at(ta["ynorm"], lv), first_at(tb["ynorm"], lv)
                    if ka is None or kb is None:
                        d5_rows.append(f"{L},{d},{meas},{lv},{ka},{kb},,,,,")
                        continue
                    ia, ib = int(ta["it"][ka]), int(tb["it"][kb])
                    qa = assessed(f"{ra}@{ia}", It)
                    qb = assessed(f"{rb}@{ib}", It)
                    qa = qa / b if qa is not None else None
                    qb = qb / b if qb is not None else None
                    qai, qbi = ta["rhostar_in"][ka] / b, tb["rhostar_in"][kb] / b
                    diff = (qa - qb) if (qa is not None and qb is not None) else None
                    if meas == "rE":
                        reached += 1
                        if diff is not None and diff >= 0.05:
                            hits += 1
                    d5_rows.append(f"{L},{d},{meas},{lv},{ia},{ib},{qa},{qb},{qai:.4f},{qbi:.4f},{diff}")
            verdict_d[d] = {"rE_levels_reached_by_both": reached, "hits_diff_ge_0.05": hits}
        if len(verdict_d) == 2:
            evaluable = all(v["rE_levels_reached_by_both"] >= 2 for v in verdict_d.values())
            holds = evaluable and all(v["hits_diff_ge_0.05"] >= 2 for v in verdict_d.values())
            verdict = ("B suppresses deepening at matched E_neg progress" if holds else
                       "not established at matched E_neg progress" if evaluable else
                       "not evaluable at matched progress")
        else:
            verdict = "not evaluable (runs missing)"
        d5[L] = {"per_density": verdict_d, "verdict": verdict}
    with open(f"{ROOT}/results/d5_matched.csv", "w") as f:
        f.write("\n".join(d5_rows) + "\n")
    print("\n".join(d5_rows))
    out["D5"] = d5

    # D2, D3, D4, D6
    def conv(r: str) -> bool:
        return r in summ and summ[r]["term"] == "CONVERGED" and "q" in summ[r]

    d3: dict[str, str] = {}
    for L in LAMS:
        for F in "AB":
            r4, r6 = f"{F}_4_{L}", f"{F}_6_{L}"
            if not (conv(r4) and conv(r6)):
                d3[f"{F}_{L}"] = "failed (not converged)"
                continue
            q4, q6 = float(summ[r4]["q"]), float(summ[r6]["q"])
            e4, e6 = float(summ[r4]["rE"]), float(summ[r6]["rE"])
            if q4 >= 1.05 and q6 >= 1.05 and e4 >= 0.10 and e6 >= 0.10 and abs(q4 - q6) < 0.05:
                d3[f"{F}_{L}"] = "concentration at convergence"
            elif q4 < 1.01 and q6 < 1.01 and abs(q4 - q6) < 0.05:
                d3[f"{F}_{L}"] = "no concentration at convergence"
            else:
                d3[f"{F}_{L}"] = "resolution-dependent"
    out["D3"] = d3
    d2 = {}
    for r in ("B_4_L1", "B_6_L1"):
        if conv(r):
            d2[r] = {"rE": summ[r]["rE"], "dominates": float(summ[r]["rE"]) < 0.40}
    out["D2"] = d2
    d4: dict[str, object] = {}
    for L in LAMS:
        r4, r6 = f"A_4_{L}", f"A_6_{L}"
        if conv(r4) and conv(r6):
            s4, s6 = summ[r4], summ[r6]
            c1 = float(s4["q"]) - float(s6["q"]) >= 0.05 and bool(s4["under"])
            c2 = (abs(float(s4["E8"]) - float(s4["E16"])) / float(s4["E16"]) > 0.02
                  and float(s4["R8"]) > 1.02 * float(s4["R4"]))
            d4[L] = {"exploitable": c1 or c2, "cond1": c1, "cond2": c2}
        else:
            d4[L] = "failed (not converged)"
    out["D4"] = d4
    out["D6"] = {F: ("not evaluable (not converged)" if "failed" in d3[f"{F}_L1"] + d3[f"{F}_L2"] else
                     "regularisation-dependent" if d3[f"{F}_L1"] != d3[f"{F}_L2"] else "consistent")
                 for F in "AB"}
    out["summary"] = summ
    with open(f"{ROOT}/results/decisions.json", "w") as f:
        json.dump(out, f, indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k != "summary"}, indent=1, default=str))

    # in-loop vs assessment agreement at assessed iterates
    agree = ["tag,assessed_rhoD,inloop_rho,rel_diff"]
    for (tag, s), a in It.items():
        run, k = tag.split("@")
        tr = load_trace(run)
        idx = np.nonzero(tr["it"].astype(int) == int(k))[0]
        if idx.size:
            ri = tr["rhostar_in"][idx[0]]
            ra = min(float(a["rhoAb"]), fl)
            agree.append(f"{tag},{ra:.10e},{ri:.10e},{(ri - ra) / abs(ra):.3e}")
    with open(f"{ROOT}/results/inloop_vs_assessment.csv", "w") as f:
        f.write("\n".join(agree) + "\n")

    figures(summ, b, E0, A, It, fl)


def figures(summ: dict, b: float, E0: float, A: dict, It: dict, fl: float) -> None:
    """Trajectory, matched-progress, resolution and meridian figures."""
    os.makedirs(f"{ROOT}/figures", exist_ok=True)
    col = {"4": "tab:blue", "6": "tab:red"}
    for L in LAMS:
        rs = [r for r in summ if r.endswith(L)]
        if not rs:
            continue
        fig, ax = plt.subplots(2, 2, figsize=(14, 9))
        for run in rs:
            tr = load_trace(run)
            ls = "-" if run[0] == "A" else "--"
            c = col[run.split("_")[1]]
            lab = f"{run[0]} {run.split('_')[1]}:1 ({summ[run]['term']})"
            ax[0, 0].plot(tr["it"], 1 - tr["Eneg_v8"] / E0, ls, color=c, label=lab)
            ax[0, 1].plot(tr["it"], tr["rhostar_in"] / b, ls, color=c, label=lab)
            ax[1, 0].plot(tr["it"], tr["ynorm"], ls, color=c, label=lab)
            ax[1, 1].semilogy(tr["it"], tr["gF"] / (1e-5 * tr["F"] / np.maximum(tr["ynorm"], 1)), ls, color=c,
                              label=lab)
        ax[0, 0].set_ylabel("r_E (8:1 validation)")
        ax[0, 1].set_ylabel("q in-loop = rho*_Dw / b")
        ax[0, 1].axhline(1.05, color="k", lw=0.6, ls=":")
        ax[1, 0].set_ylabel("||y||")
        ax[1, 1].set_ylabel("||grad F|| / C1' tolerance")
        ax[1, 1].axhline(1.0, color="k", lw=0.6, ls=":")
        for a_ in ax.flat:
            a_.set_xlabel("accepted iteration")
            a_.grid(alpha=0.3)
            a_.legend(fontsize=8)
        fig.suptitle(f"Study 3, lambda = {LAMS[L]:.3g}, start S1, Eulerian rho only. A solid, B dashed")
        fig.tight_layout()
        fig.savefig(f"{ROOT}/figures/trajectory_{L}.png", dpi=120)
        plt.close(fig)

        fig, ax = plt.subplots(1, 2, figsize=(14, 5))
        for run in rs:
            tr = load_trace(run)
            ls = "-" if run[0] == "A" else "--"
            c = col[run.split("_")[1]]
            lab = f"{run[0]} {run.split('_')[1]}:1"
            ax[0].plot(1 - tr["Eneg_v8"] / E0, tr["rhostar_in"] / b, ls, color=c, lw=0.9, label=lab + " (in-loop)")
            ax[1].plot(tr["ynorm"], tr["rhostar_in"] / b, ls, color=c, lw=0.9, label=lab + " (in-loop)")
            pts = [(int(t.split("@")[1]), a) for (t, s), a in It.items() if t.split("@")[0] == run]
            if pts:
                idx = {int(k): i for i, k in enumerate(tr["it"].astype(int))}
                xe = [1 - float(a["E8"]) / E0 for k, a in pts]
                xy = [tr["ynorm"][idx[k]] for k, a in pts if k in idx]
                qq = [min(float(a["rhoAb"]), fl) / b for k, a in pts]
                mk = "o" if run[0] == "A" else "s"
                ax[0].plot(xe, qq, mk, color=c, mfc="none", ms=6)
                ax[1].plot(xy, qq[: len(xy)], mk, color=c, mfc="none", ms=6)
        for lv in RE_LEVELS:
            ax[0].axvline(lv, color="gray", lw=0.5, ls=":")
        for lv in Y_LEVELS:
            ax[1].axvline(lv, color="gray", lw=0.5, ls=":")
        ax[0].set_xlabel("r_E (progress in validation E_neg)")
        ax[1].set_xlabel("||y|| (design norm)")
        for a_ in ax:
            a_.set_ylabel("q (lines in-loop; markers assessed)")
            a_.axhline(1.05, color="k", lw=0.6, ls=":")
            a_.grid(alpha=0.3)
            a_.legend(fontsize=8)
        fig.suptitle(f"Matched progress (D5), lambda = {LAMS[L]:.3g}. A solid/circles, B dashed/squares")
        fig.tight_layout()
        fig.savefig(f"{ROOT}/figures/matched_{L}.png", dpi=120)
        plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for L, ms in (("L1", "-"), ("L2", ":")):
        for F, mk in (("A", "o"), ("B", "s")):
            xs = [int(d) for d in ("4", "6") if f"{F}_{d}_{L}" in summ and "q" in summ[f"{F}_{d}_{L}"]]
            if xs:
                ax[0].plot(xs, [summ[f"{F}_{d}_{L}"]["q"] for d in map(str, xs)], mk + ms, label=f"{F} {L}")
                ax[1].plot(xs, [summ[f"{F}_{d}_{L}"]["rE"] for d in map(str, xs)], mk + ms, label=f"{F} {L}")
    ax[0].axhline(1.05, color="k", lw=0.6, ls=":")
    ax[0].set_ylabel("final q (assessment)")
    ax[1].set_ylabel("final r_E")
    for a_ in ax:
        a_.set_xlabel("working Gauss order")
        a_.set_xticks([4, 6])
        a_.grid(alpha=0.3)
        a_.legend(fontsize=8)
    fig.suptitle("Final iterates (converged or not; see summary.csv)")
    fig.tight_layout()
    fig.savefig(f"{ROOT}/figures/resolution.png", dpi=120)
    plt.close(fig)

    xs = -2.1 + 0.01 * np.arange(421)
    rr = 0.01 * np.arange(311)
    for L in LAMS:
        for d in ("4", "6"):
            names = ["baseline", f"A_{d}_{L}", f"B_{d}_{L}"]
            if not all(os.path.exists(f"{ROOT}/fields/{n}.txt") for n in names):
                continue
            vals = [np.loadtxt(f"{ROOT}/fields/{n}.txt") for n in names]
            vmax = max(np.abs(v).max() for v in vals)
            norm = SymLogNorm(linthresh=1e-5, vmin=-vmax, vmax=vmax)
            fig, ax = plt.subplots(1, 3, figsize=(17, 5))
            for a_, n, v in zip(ax, names, vals):
                im = a_.pcolormesh(xs, rr, v.T, cmap="RdBu", norm=norm, shading="auto")
                a_.contour(xs, rr, v.T, levels=[0.0], colors="k", linewidths=0.4)
                if (n, 1) in A:
                    mp = f"{ROOT}/assess/assess.csv.{n}.s1.minima.csv"
                    if os.path.exists(mp):
                        m = np.loadtxt(mp, delimiter=",", skiprows=1, ndmin=2)
                        a_.plot(m[:, 0], m[:, 1], "x", color="lime", ms=6, label="distinct minima (NM)")
                    a_.plot(A[(n, 1)]["xD"], A[(n, 1)]["rD"], "*", color="yellow", mec="k", ms=14, label="rho*")
                    a_.set_title(f"{n}: rho*={float(A[(n, 1)]['rhoAb']):.3e}, E_neg={float(A[(n, 1)]['E8']):.3e}",
                                 fontsize=9)
                if n.startswith("B"):
                    with open(f"{ROOT}/runs/{n}/refine_events.csv") as f:
                        er = [r for r in csv.DictReader(f) if r["action"] == "split"]
                    if er:
                        a_.plot([float(r["x"]) for r in er], [float(r["r"]) for r in er], ".", color="m", ms=3,
                                label="B split triggers")
                a_.set_aspect("equal")
                a_.set_xlabel("x")
                a_.set_ylabel("r")
                a_.legend(fontsize=7, loc="upper right")
            fig.colorbar(im, ax=ax, shrink=0.8, label="Eulerian rho (symlog)")
            fig.suptitle(f"Eulerian rho on the active box Ab, {d}:1, lambda {LAMS[L]:.3g}. Black: rho = 0")
            fig.savefig(f"{ROOT}/figures/meridian_{d}_{L}.png", dpi=110)
            plt.close(fig)


if __name__ == "__main__":
    {"select": select, "report": report}[sys.argv[1]]()
