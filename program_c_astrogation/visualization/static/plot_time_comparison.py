#!/usr/bin/env python3
"""Trajectory-family comparison: clocks, peak velocity and arrival epoch.

    python visualization/static/plot_time_comparison.py \
        results/verification/program_c1/data/trajectory_comparison.csv \
        --metadata results/verification/program_c1/c1_metadata.json \
        --out results/verification/program_c1/figures

READS COLUMNS AND DRAWS THEM. Every time was computed and saved by c1_verify.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import cli, records
from common import render as rnd
from themes import themes


def main(argv=None) -> int:
    p = cli.parser(__doc__.splitlines()[0])
    a = p.parse_args(argv)

    rows = records.read_csv(a.input)
    meta = records.RunMetadata.load(a.metadata)
    out_dir = records.ensure_dir(a.out)
    theme = themes.get(a.theme)
    themes.apply(theme)

    dests = list(dict.fromkeys(records.text_column(rows, "destination")))
    fig, axes = plt.subplots(1, len(dests), figsize=(6.4 * len(dests), 5.6),
                             sharey=False)
    if len(dests) == 1:
        axes = [axes]

    plotted = []
    for ax, dest in zip(axes, dests):
        sub = [r for r in rows if r["destination"] == dest]
        labels, coord, proper, peak = [], [], [], []
        for r in sub:
            model = r["model"]
            g = float(r["acceleration_g"])
            labels.append(model if g == 0 else "%s\n%.2g g" % (model[:9], g))
            coord.append(float(r["coordinate_time_yr"]))
            proper.append(float(r["proper_time_yr"]))
            peak.append(float(r["peak_beta"]))
            plotted.append([dest, r["source_id"], model, g, coord[-1],
                            proper[-1], peak[-1], float(r["arrival_epoch"])])

        idx = range(len(labels))
        w = 0.38
        ax.bar([i - w / 2 for i in idx], coord, w, color=theme.accent,
               label="coordinate time (solar-system frame)")
        ax.bar([i + w / 2 for i in idx], proper, w, color=theme.accent_two,
               label="ship proper time")
        for i in idx:
            ax.text(i, max(coord[i], proper[i]) * 1.02,
                    r"$\beta_{max}$=%.3f" % peak[i], ha="center", fontsize=6.6,
                    color=theme.foreground)
        ax.set_xticks(list(idx))
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_yscale("log")
        ax.set_ylabel("elapsed time [yr, log scale]")
        ax.set_title("%s" % dest, fontsize=9.5)
        ax.grid(True, alpha=theme.grid_alpha, axis="y")
        ax.legend(fontsize=7, facecolor=theme.background, edgecolor=theme.grid,
                  labelcolor=theme.foreground, loc="upper left")

    fig.suptitle("Which clock? Coordinate time and ship proper time are "
                 "different numbers for the same journey", fontsize=11)
    rnd.scope_note(fig, theme,
                   "ConstantCoordinateVelocity is a VERIFICATION MODEL: it needs "
                   "infinite acceleration at both ends. No propulsion energy is "
                   "modelled for any profile, so none is 'more realizable'.")

    prov = meta.figure_provenance(
        color_model="illustrative_only",
        render_parameters=dict(theme.describe(), y_scale="log",
                               bars="saved coordinate and proper times"),
        observer_frame="BarycentricCartesian and the ship rest frame, labelled "
                       "separately on every bar pair")
    for ext in ("png", "svg"):
        pth = os.path.join(out_dir, "time_comparison." + ext)
        fig.savefig(pth, facecolor=fig.get_facecolor())
        records.write_sidecar(pth, prov)
        print("wrote %s" % pth)
    plt.close(fig)

    records.write_plotted_csv(
        os.path.join(out_dir, "time_comparison.plotted.csv"),
        ["destination", "source_id", "model", "acceleration_g",
         "coordinate_time_yr", "proper_time_yr", "peak_beta", "arrival_epoch"],
        plotted)
    return 0


if __name__ == "__main__":
    sys.exit(main())
