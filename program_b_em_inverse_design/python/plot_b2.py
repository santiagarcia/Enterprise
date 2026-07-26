#!/usr/bin/env python3
"""Figures for Program B milestone B2 (abstract constitutive inverse design).

Every figure is drawn from files the gate wrote. Nothing here reruns the
optimizer, and nothing depends on a hand-edited intermediate: if a CSV is
missing the corresponding panel is skipped and reported, rather than fabricated.

    python plot_b2.py <results/verification/program_b2>
"""

import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


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


def plot_histories(out_dir, rows):
    problems = {}
    for r in rows:
        problems.setdefault(r["problem"], []).append(r)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for name, rs in sorted(problems.items()):
        it = [int(r["iteration"]) for r in rs]
        obj = [max(float(r["objective"]), 1e-30) for r in rs]
        gn = [max(float(r["gradient_norm"]), 1e-30) for r in rs]
        ax1.semilogy(it, obj, label=name, lw=1.4)
        ax2.semilogy(it, gn, label=name, lw=1.4)

    ax1.set_xlabel("iteration")
    ax1.set_ylabel(r"$J$")
    ax1.set_title("Objective versus iteration")
    ax1.grid(alpha=0.3, which="both")
    ax2.set_xlabel("iteration")
    ax2.set_ylabel(r"$\|\nabla J\|$")
    ax2.set_title("Gradient norm versus iteration")
    ax2.grid(alpha=0.3, which="both")
    ax2.legend(fontsize=8)
    savefig(fig, out_dir, "b2_optimizer_histories.png")


def plot_recovered(out_dir, rows):
    problems = {}
    for r in rows:
        problems.setdefault(r["problem"], []).append(r)

    n = len(problems)
    fig, axes = plt.subplots(1, n, figsize=(4.0 * n, 3.6), squeeze=False)
    for ax, (name, rs) in zip(axes[0], sorted(problems.items())):
        idx = [int(r["index"]) for r in rs]
        tgt = [float(r["target_q"]) for r in rs]
        rec = [float(r["recovered_q"]) for r in rs]
        ax.plot(idx, tgt, "o", ms=8, mfc="none", label="target")
        ax.plot(idx, rec, "x", ms=6, label="recovered")
        ax.set_xlabel("parameter index")
        ax.set_ylabel("value")
        ax.set_title(name)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("Target versus recovered design parameters", y=1.02)
    savefig(fig, out_dir, "b2_recovered_parameters.png")


def plot_residual_map(out_dir, rows):
    pairs = ["01", "02", "03", "23", "31", "12"]
    idx = {p: i for i, p in enumerate(pairs)}
    exact = [[0.0] * 6 for _ in range(6)]
    restricted = [[0.0] * 6 for _ in range(6)]
    for r in rows:
        i, j = idx[r["pair_A"]], idx[r["pair_B"]]
        exact[i][j] = float(r["exact_residual"])
        restricted[i][j] = float(r["restricted_residual"])

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, data, title in (
        (axes[0], exact, "unrestricted family"),
        (axes[1], restricted, r"restricted family, $\zeta \equiv 0$"),
    ):
        im = ax.imshow(data, cmap="magma")
        ax.set_xticks(range(6), pairs)
        ax.set_yticks(range(6), pairs)
        ax.set_xlabel("index pair (cd)")
        ax.set_ylabel("index pair (ab)")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle(r"$|Z_{\rm design} - Z_{\rm target}|$ by antisymmetric index pair",
                 y=1.02)
    savefig(fig, out_dir, "b2_z_residual_by_pair.png")


def plot_family_residuals(out_dir, rows):
    names = [r["family"] for r in rows]
    npar = [int(r["parameters"]) for r in rows]
    J = [max(float(r["least_squares_J_Z"]), 1e-36) for r in rows]

    fig, ax = plt.subplots(figsize=(7, 4.0))
    bars = ax.bar(range(len(names)), J, color="#4C6EF5")
    ax.set_yscale("log")
    ax.set_xticks(range(len(names)),
                  [f"{n}\n({p} params)" for n, p in zip(names, npar)], fontsize=8)
    ax.set_ylabel(r"least-squares $J_Z$")
    ax.set_title("Cost of restricting the design family\n"
                 "(same magnetoelectric target)")
    ax.grid(alpha=0.3, axis="y", which="both")
    for b, v in zip(bars, J):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.5, f"{v:.2e}",
                ha="center", fontsize=7)
    savefig(fig, out_dir, "b2_family_residuals.png")


def plot_fd_convergence(out_dir, rows):
    h = [float(r["step"]) for r in rows]
    e = [max(float(r["worst_abs_error"]), 1e-20) for r in rows]

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.loglog(h, e, "o-", lw=1.4)
    ax.set_xlabel("finite-difference step $h$")
    ax.set_ylabel(r"worst $|\nabla J_{\rm FD} - \nabla J_{\rm analytic}|$")
    ax.set_title("Finite-difference error versus step size\n"
                 "(objective is exactly quadratic, so the error is roundoff only)")
    ax.grid(alpha=0.3, which="both")
    savefig(fig, out_dir, "b2_fd_convergence.png")


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    out_dir = sys.argv[1]

    panels = [
        ("objective_histories.csv", plot_histories),
        ("recovered_designs.csv", plot_recovered),
        ("z_residual_by_pair.csv", plot_residual_map),
        ("family_residuals.csv", plot_family_residuals),
        ("gradient_fd_convergence.csv", plot_fd_convergence),
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
        print("no figures drawn: run b2_verify first")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
