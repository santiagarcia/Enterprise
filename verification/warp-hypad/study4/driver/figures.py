"""Study-4 figures (Anaconda python -I). Usage: python -I figures.py STUDY4_DIR"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

root = Path(sys.argv[1]); fig_dir = root / "figures"; fig_dir.mkdir(exist_ok=True)
B_REF = -9.60726526e-4


def tr(r):
    p = root / "runs" / r / "trace.csv"
    return list(csv.DictReader(open(p))) if p.exists() else None


# 1. spectra of step-0 designs
designs = [d.name for d in sorted((root / "hessians").iterdir()) if (d / "eig.csv").exists()]
fig, ax = plt.subplots(figsize=(9, 5))
for d in designs:
    mu = np.array([float(r["mu_HF"]) for r in csv.DictReader(open(root / "hessians" / d / "eig.csv"))])
    ax.semilogy(np.arange(len(mu)), np.abs(mu), ".-" if d.startswith("final") else "-", lw=0.8, ms=2, label=d)
    neg = mu < 0
    ax.semilogy(np.arange(len(mu))[neg], np.abs(mu[neg]), "kx", ms=4)
ax.set_xlabel("eigenvalue index (ascending mu)"); ax.set_ylabel("|mu(H_F)|  (x = negative)")
ax.set_title("Measured spectra of H_F (OTI-assembled)"); ax.legend(fontsize=6, ncol=2)
fig.tight_layout(); fig.savefig(fig_dir / "spectra.png", dpi=130); plt.close(fig)

# 2. in-loop q trajectories (P1 pairs and others) and |grad F|/tol
for grp, runs in (("P1_lambda2", ["A4L2m", "B4L2", "A6L2m", "B6L2"]),
                  ("P2_P3_P5", ["A4L2f", "A6L2f", "X4L1fromL2", "X4L2fromL1", "XB4L2fromA"]),
                  ("P4_lambda1", ["A4L1f", "B4L1f", "A6L1f", "B6L1f"])):
    data = {r: tr(r) for r in runs if tr(r)}
    if not data:
        continue
    fig, axs = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for r, t in data.items():
        it = [int(x["it"]) for x in t]
        axs[0].plot(it, [float(x["rhostar_in"]) / B_REF for x in t], label=r, lw=1)
        axs[1].semilogy(it, [float(x["gF"]) / (1e-5 * float(x["F"]) / max(float(x["ynorm"]), 1)) for x in t], label=r, lw=1)
    axs[0].axhline(1.05, color="k", ls=":", lw=0.8); axs[0].set_ylabel("in-loop q")
    for k in range(0, 201, 25):
        axs[0].axvline(k, color="0.85", lw=0.5)
    axs[1].axhline(1, color="k", ls=":", lw=0.8); axs[1].set_ylabel("|grad F| / tol_C1'")
    axs[1].set_xlabel("accepted study-4 iterate"); axs[0].legend(fontsize=7); axs[0].set_title(grp)
    fig.tight_layout(); fig.savefig(fig_dir / f"trajectory_{grp}.png", dpi=130); plt.close(fig)
print("figures written to", fig_dir)
