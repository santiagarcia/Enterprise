#!/usr/bin/env python3
"""Figures for Program B milestone B4 (Maxwell propagation).

Drawn only from files the gate wrote.

    python plot_b4.py <results/verification/program_b4>
"""

import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv(path):
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def savefig(fig, out_dir, name):
    figdir = os.path.join(out_dir, "figures")
    os.makedirs(figdir, exist_ok=True)
    path = os.path.join(figdir, name)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


def plot_dispersion_surfaces(out_dir, rows):
    media = {}
    for r in rows:
        media.setdefault(r["medium"], []).append(r)

    n = len(media)
    cols = min(n, 3)
    rowsn = (n + cols - 1) // cols
    fig, axes = plt.subplots(rowsn, cols, figsize=(4.2 * cols, 4.2 * rowsn),
                             subplot_kw={"projection": "polar"}, squeeze=False)
    flat = [a for row in axes for a in row]

    for ax, (name, rs) in zip(flat, media.items()):
        th = [float(r["theta"]) for r in rs]
        b1 = [float(r["kappa_branch1"]) for r in rs]
        b2 = [float(r["kappa_branch2"]) for r in rs]
        # mirror to the full circle for readability
        th_full = th + [t + 3.141592653589793 for t in th]
        ax.plot(th_full, b1 + b1, lw=1.6, label="branch 1")
        if max(b2) > 0:
            ax.plot(th_full, b2 + b2, lw=1.2, ls="--", label="branch 2")
        ax.set_title(name, fontsize=9, pad=14)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6, loc="lower right")
    for ax in flat[len(media):]:
        ax.axis("off")

    fig.suptitle(r"Constant-frequency dispersion surfaces, $\kappa = k/\omega$"
                 "\ntwo coincident branches = no birefringence = a single "
                 "light cone", y=1.02)
    savefig(fig, out_dir, "b4_dispersion_surfaces.png")


def plot_profile(out_dir, rows):
    x = [float(r["x"]) for r in rows]
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.plot(x, [float(r["eps"]) for r in rows], lw=1.6,
            label=r"$\varepsilon$")
    ax.plot(x, [float(r["mu_inverse"]) for r in rows], lw=1.6,
            label=r"$\mu^{-1}$")
    ax.plot(x, [float(r["zeta"]) for r in rows], lw=1.8,
            label=r"$\zeta$  (dropped by the restricted model)")
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.set_xlabel("x")
    ax.set_ylabel("block value")
    ax.set_title("Layered warp-analogue profile\n"
                 "the constitutive stack propagation is computed through")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    savefig(fig, out_dir, "b4_profile.png")


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    out_dir = sys.argv[1]

    panels = [
        ("b4_dispersion_surfaces.csv", plot_dispersion_surfaces),
        ("b4_profile.csv", plot_profile),
    ]
    drawn = 0
    for name, fn in panels:
        rows = read_csv(os.path.join(out_dir, name))
        if not rows:
            print("skipped: no", name)
            continue
        fn(out_dir, rows)
        drawn += 1

    if drawn == 0:
        print("no figures drawn: run b4_verify first")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
