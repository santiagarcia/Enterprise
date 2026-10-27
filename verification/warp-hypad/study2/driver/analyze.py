"""Study-2 analysis: decision-rule table and figures (PROTOCOL.md sections 6-7).

Run with the Anaconda interpreter in isolated mode from the driver directory:
    C:/Users/vfn333/AppData/Local/anaconda3/python.exe -I analyze.py
"""

from __future__ import annotations

import csv
import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import SymLogNorm  # noqa: E402

ROOT = "C:/Users/vfn333/Documents/GitHub/pub-feasibility/warp_hypad/study2"
COLS = ("tag,yfile,seed,rhoD,xD,rD,rhoW,xW,rW,rhoA,xA,rA,rhoI,xI,rI,sample_min,ncells,E8,nf8,nm8,"
        "E16,nf16,nm16,axl,axh,arh,secs").split(",")


def load_assess() -> dict[tuple[str, int], dict[str, float | str]]:
    """Return assessment rows keyed by (tag, seed); later rows win."""
    out: dict[tuple[str, int], dict[str, float | str]] = {}
    with open(f"{ROOT}/assess/assess.csv") as f:
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
    """Load a run's per-iterate trace as column arrays."""
    with open(f"{ROOT}/runs/{run}/trace.csv") as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def termination(run: str) -> str:
    """Return the termination line of a run, or RUNNING/MISSING."""
    p = f"{ROOT}/runs/{run}/termination.txt"
    return open(p).read().strip() if os.path.exists(p) else "MISSING"


def main() -> None:
    """Build the summary table and all figures from archived outputs."""
    A = load_assess()
    base = A[("baseline", 1)]
    b, E0 = float(base["rhoD"]), float(base["E8"])
    floor = float(base["rhoI"])
    runs = sorted(os.path.basename(p) for p in glob.glob(f"{ROOT}/runs/[AB]_*"))
    lines = ["run,termination,iters,q,r_E,rhoD,rhoD_seed2,rhoDw,rho_active,E8,E16,negfrac8,check_worst,"
             "underreport,qrange_last20,ErangeLast20"]
    summ: dict[str, dict[str, float | str]] = {}
    for run in runs:
        if (run, 1) not in A:
            continue
        a = A[(run, 1)]
        a2 = A.get((run, 2))
        tr = load_trace(run)
        rhoD = float(a["rhoD"])
        q = rhoD / b
        rE = 1 - float(a["E8"]) / E0
        # fixed-node check values at the final iterate: working node min, validation node min, repo refine
        chk = min(tr["rhomin_w"][-1], tr["rhomin_v8"][-1], tr["refined_v8"][-1])
        under = (chk - rhoD) > 0.10 * abs(rhoD)
        last = slice(max(0, len(tr["it"]) - 20), None)
        qr = (tr["rhostar_in"][last].min() / b, tr["rhostar_in"][last].max() / b)
        er = (1 - tr["Eneg_v8"][last].max() / E0, 1 - tr["Eneg_v8"][last].min() / E0)
        term = termination(run).split()[0]
        summ[run] = {"q": q, "rE": rE, "term": term, "under": under}
        lines.append(
            f"{run},{term},{int(tr['it'][-1])},{q:.4f},{rE:.4f},{rhoD:.6e},"
            f"{float(a2['rhoD']) if a2 else float('nan'):.6e},{float(a['rhoW']):.6e},{float(a['rhoA']):.6e},"
            f"{float(a['E8']):.6e},{float(a['E16']):.6e},{float(a['nf8']):.4f},{chk:.6e},{int(under)},"
            f"{qr[0]:.3f}-{qr[1]:.3f},{er[0]:.3f}-{er[1]:.3f}")
    with open(f"{ROOT}/results/summary.csv", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"baseline rho*_D {b:.8e}, E0 {E0:.8e}, design-independent floor {floor:.6e} (q_floor {floor / b:.3f})")

    os.makedirs(f"{ROOT}/figures", exist_ok=True)
    starts = sorted({r.split("_")[2] for r in runs})
    dens = sorted({r.split("_")[1] for r in runs})

    # Iteration dependence
    for st in starts:
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
        for run in runs:
            if run.split("_")[2] != st:
                continue
            tr = load_trace(run)
            ls = "-" if run[0] == "A" else "--"
            col = {"4": "tab:blue", "6": "tab:red"}[run.split("_")[1]]
            lab = f"{run[0]} {run.split('_')[1]}:1"
            ax[0].plot(tr["it"], 1 - tr["Eneg_v8"] / E0, ls, color=col, label=lab)
            ax[1].plot(tr["it"], tr["rhostar_in"] / b, ls, color=col, label=lab)
            ax[2].plot(tr["it"], tr["rhomin_w"] / b, ls, color=col, label=lab + " working nodes")
            ax[2].plot(tr["it"], tr["rhomin_v8"] / b, ls, color=col, alpha=0.35)
        ax[0].set_ylabel("r_E = 1 - E_neg/E_neg(baseline)  [8:1 validation]")
        ax[1].set_ylabel("q = rho*/rho*_baseline  [in-loop continuum search, Dw]")
        ax[2].set_ylabel("node rho_min / rho*_baseline (solid-ish: working; faint: 8:1)")
        ax[1].axhline(1.05, color="k", lw=0.6, ls=":")
        ax[1].axhline(floor / b, color="gray", lw=0.6, ls="-.")
        for a_ in ax:
            a_.set_xlabel("accepted iteration")
            a_.legend(fontsize=8)
            a_.grid(alpha=0.3)
        fig.suptitle(f"Iteration dependence, start {st} (Eulerian rho only). A solid, B dashed; "
                     "dotted q=1.05, dash-dot design-independent floor")
        fig.tight_layout()
        fig.savefig(f"{ROOT}/figures/iteration_{st}.png", dpi=130)
        plt.close(fig)

    # Resolution dependence
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for st in starts:
        for F, mk in (("A", "o-"), ("B", "s--")):
            xs, qs, es = [], [], []
            for d in dens:
                r = f"{F}_{d}_{st}"
                if r in summ:
                    xs.append(int(d)); qs.append(summ[r]["q"]); es.append(summ[r]["rE"])
            if xs:
                ax[0].plot(xs, qs, mk, label=f"{F} {st}")
                ax[1].plot(xs, es, mk, label=f"{F} {st}")
    ax[0].axhline(1.0, color="k", lw=0.6)
    ax[0].axhline(1.05, color="k", lw=0.6, ls=":")
    ax[0].set_ylabel("q = rho*_D / rho*_D(baseline)  (final, assessment)")
    ax[1].set_ylabel("r_E (final, 8:1)")
    for a_ in ax:
        a_.set_xlabel("working Gauss order (sub 1)")
        a_.set_xticks([4, 6])
        a_.legend(fontsize=8)
        a_.grid(alpha=0.3)
    fig.suptitle("Resolution dependence of final designs")
    fig.tight_layout()
    fig.savefig(f"{ROOT}/figures/resolution.png", dpi=130)
    plt.close(fig)

    # Meridian maps
    def field(name: str) -> np.ndarray:
        return np.loadtxt(f"{ROOT}/fields/{name}.txt")

    xs = np.linspace(-7, 7, 701)
    rs = np.linspace(0, 7, 351)
    for d in dens:
        for st in starts:
            names = ["baseline", f"A_{d}_{st}", f"B_{d}_{st}"]
            if not all(os.path.exists(f"{ROOT}/fields/{n}.txt") for n in names):
                continue
            fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
            vals = [field(n) for n in names]
            vmax = max(np.abs(v).max() for v in vals)
            norm = SymLogNorm(linthresh=1e-5, vmin=-vmax, vmax=vmax)
            for a_, n, v in zip(ax, names, vals):
                im = a_.pcolormesh(xs, rs, v.T, cmap="RdBu", norm=norm, shading="auto")
                a_.contour(xs, rs, v.T, levels=[0.0], colors="k", linewidths=0.4)
                key = (n, 1)
                if key in A:
                    mp = f"{ROOT}/assess/assess.csv.{n}.minima.csv"
                    if os.path.exists(mp):
                        m = np.loadtxt(mp, delimiter=",", skiprows=1, ndmin=2)
                        a_.plot(m[:8, 0], m[:8, 1], "x", color="lime", ms=6, label="8 deepest local minima")
                    a_.plot(A[key]["xD"], A[key]["rD"], "*", color="yellow", mec="k", ms=14, label="rho*_D")
                    a_.set_title(f"{n}: rho*_D={float(A[key]['rhoD']):.3e}, E_neg={float(A[key]['E8']):.3e}",
                                 fontsize=9)
                if n.startswith("B"):
                    ev = f"{ROOT}/runs/{n}/refine_events.csv"
                    with open(ev) as f:
                        er = [r for r in csv.DictReader(f) if r["action"] == "split"]
                    if er:
                        a_.plot([float(r["x"]) for r in er], [float(r["r"]) for r in er], ".", color="m", ms=3,
                                label="B split triggers")
                a_.add_patch(plt.Rectangle((-3.4, 0), 6.8, 2 * np.sqrt(2.6), fill=False, ec="gray", ls="--", lw=0.8))
                a_.set_xlim(-5, 5)
                a_.set_ylim(0, 5)
                a_.set_aspect("equal")
                a_.set_xlabel("x")
                a_.set_ylabel("r (cylindrical)")
                a_.legend(fontsize=7, loc="upper right")
            fig.colorbar(im, ax=ax, shrink=0.8, label="Eulerian rho (symlog, linthresh 1e-5)")
            fig.suptitle(f"Eulerian energy density on the meridian, {d}:1, start {st}. "
                         "Black: rho=0; dashed: working domain Dw")
            fig.savefig(f"{ROOT}/figures/meridian_{d}_{st}.png", dpi=120)
            plt.close(fig)


if __name__ == "__main__":
    main()
