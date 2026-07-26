#!/usr/bin/env python3
"""Figures for Program B milestone B3 (physical admissibility).

Drawn only from files the gate wrote. A missing CSV skips its panel and says so
rather than inventing data.

    python plot_b3.py <results/verification/program_b3>
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


def plot_dispersion(out_dir, rows):
    w = [float(r["omega"]) for r in rows]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.0))

    for ax, (re_k, im_k, title) in zip(
        axes,
        [("re_eps", "im_eps", r"$\varepsilon(\omega)$"),
         ("re_mu", "im_mu", r"$\mu(\omega)$"),
         ("re_zeta", "im_zeta", r"$\zeta(\omega)$")],
    ):
        ax.plot(w, [float(r[re_k]) for r in rows], lw=1.5, label="real")
        ax.plot(w, [float(r[im_k]) for r in rows], lw=1.5, label="imaginary")
        ax.axhline(0.0, color="0.6", lw=0.8)
        ax.set_xlabel(r"$\omega$")
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(r"Dispersive response, kernel $e^{-i\omega t}$ "
                 r"(passive $\Rightarrow \operatorname{Im}\varepsilon \geq 0$)",
                 y=1.03)
    savefig(fig, out_dir, "b3_dispersion.png")


def plot_passivity(out_dir, rows):
    w = [float(r["omega"]) for r in rows]
    lam = [float(r["min_loss_eigenvalue"]) for r in rows]
    im_eps = [float(r["im_eps"]) for r in rows]
    im_mu = [float(r["im_mu"]) for r in rows]

    fig, ax = plt.subplots(figsize=(8, 4.4))
    ax.plot(w, lam, lw=1.8, color="#C92A2A",
            label=r"$\lambda_{\min}(L)$  (the actual passivity condition)")
    ax.plot(w, im_eps, lw=1.2, ls="--", color="#4C6EF5",
            label=r"$\operatorname{Im}\varepsilon$")
    ax.plot(w, im_mu, lw=1.2, ls=":", color="#2F9E44",
            label=r"$\operatorname{Im}\mu$")
    ax.axhline(0.0, color="0.4", lw=1.0)
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel("value")
    ax.set_title("Passivity: the loss-matrix spectrum, not the diagonal "
                 "imaginary parts\n"
                 r"passive $\Leftrightarrow \lambda_{\min}(L) \geq 0$ for all $\omega > 0$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    savefig(fig, out_dir, "b3_passivity.png")


def plot_histories(out_dir, rows):
    problems = {}
    for r in rows:
        problems.setdefault(r["problem"], []).append(r)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for name, rs in sorted(problems.items()):
        it = [int(r["iteration"]) for r in rs]
        ax1.semilogy(it, [max(float(r["objective"]), 1e-30) for r in rs],
                     lw=1.3, label=name)
        ax2.semilogy(it, [max(float(r["gradient_norm"]), 1e-30) for r in rs],
                     lw=1.3, label=name)
    ax1.set_xlabel("iteration"); ax1.set_ylabel(r"$J_{\rm band}$")
    ax1.set_title("Band objective versus iteration (L-BFGS)")
    ax1.grid(alpha=0.3, which="both")
    ax2.set_xlabel("iteration"); ax2.set_ylabel(r"$\|\nabla J\|$")
    ax2.set_title("Gradient norm versus iteration")
    ax2.grid(alpha=0.3, which="both")
    ax2.legend(fontsize=7)
    savefig(fig, out_dir, "b3_optimizer_histories.png")


def plot_bandwidth(out_dir, rows):
    hw = [float(r["half_width"]) for r in rows]
    J = [max(float(r["J_band"]), 1e-30) for r in rows]
    worst = [max(float(r["worst_mismatch"]), 1e-30) for r in rows]
    centre = [max(float(r["centre_mismatch"]), 1e-30) for r in rows]

    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.loglog(hw, J, "o-", lw=1.6, label=r"$J_{\rm band}$")
    ax.loglog(hw, worst, "s--", lw=1.2, label="worst-frequency mismatch")
    ax.loglog(hw, centre, "^:", lw=1.2, label="centre mismatch")
    ax.set_xlabel("band half-width")
    ax.set_ylabel("mismatch")
    ax.set_title("Bandwidth tradeoff\n"
                 "a one-pole model fitted to a two-pole target")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    savefig(fig, out_dir, "b3_bandwidth_tradeoff.png")


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    out_dir = sys.argv[1]

    panels = [
        ("b3_passivity_sweep.csv", plot_dispersion),
        ("b3_passivity_sweep.csv", plot_passivity),
        ("b3_objective_histories.csv", plot_histories),
        ("b3_bandwidth_tradeoff.csv", plot_bandwidth),
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
        print("no figures drawn: run b3_verify first")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
