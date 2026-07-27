#!/usr/bin/env python3
"""Moving-destination interception, against the naive stale-position aim.

    python visualization/static/plot_destination_interception.py \
        results/verification/program_c1/data/interception.csv \
        --metadata results/verification/program_c1/c1_metadata.json \
        --out results/verification/program_c1/figures

READS COLUMNS AND DRAWS THEM. The four positions -- stale catalogue target,
destination at departure, destination at arrival, interception point -- were all
solved and saved by c1_verify. Nothing here propagates a star.
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

    # Group the saved ledger by destination, preserving file order so the figure
    # is deterministic.
    dests = []
    by_dest = {}
    for r in rows:
        key = (r["destination"], r["source_id"])
        if key not in by_dest:
            by_dest[key] = {}
            dests.append(key)
        by_dest[key][r["quantity"]] = r

    fig, axes = plt.subplots(1, len(dests), figsize=(6.2 * len(dests), 5.8))
    if len(dests) == 1:
        axes = [axes]

    plotted = []
    for ax, key in zip(axes, dests):
        name, sid = key
        g = by_dest[key]

        def vec(q):
            return (float(g[q]["x_pc"]), float(g[q]["y_pc"]), float(g[q]["z_pc"]))

        stale = vec("stale_catalogue_target")
        dep = vec("destination_at_departure")
        arr = vec("destination_at_arrival")
        err_au = float(g["stale_targeting_error"]["scalar"])
        unc_pc = float(g["interception_uncertainty_pc"]["scalar"])
        arrival_epoch = float(g["arrival_epoch"]["scalar"])
        flight = float(g["flight_coordinate_time_yr"]["scalar"])

        # Project onto the plane containing the three saved positions. This is
        # PRESENTATION geometry -- a projection of saved coordinates, which
        # cannot disagree with them -- not a recomputation.
        import math
        n = math.sqrt(sum(v * v for v in arr)) or 1.0
        e1 = tuple(v / n for v in arr)
        raw = tuple(stale[i] - arr[i] for i in range(3))
        dotp = sum(raw[i] * e1[i] for i in range(3))
        perp = tuple(raw[i] - dotp * e1[i] for i in range(3))
        pn = math.sqrt(sum(v * v for v in perp))
        e2 = tuple(v / pn for v in perp) if pn > 1e-18 else (0.0, 0.0, 1.0)

        def proj(v):
            return (sum(v[i] * e1[i] for i in range(3)),
                    sum(v[i] * e2[i] for i in range(3)))

        # Astronomical units on the transverse axis: a parsec-scale transverse
        # axis would render the whole story as a single pixel.
        AU = records.resolve  # (unused; kept out of the numeric path)
        pc_to_au = 206264.806247
        pts = {
            "stale catalogue target (J%s)" % meta.get("reference_epoch", "?"): (stale, theme.uncertainty_color, "x"),
            "destination at departure": (dep, theme.grid, "s"),
            "destination at arrival = interception": (arr, theme.accent, "o"),
        }
        for label, (v, col, mk) in pts.items():
            u, w = proj(v)
            ax.scatter([(u - proj(arr)[0]) * pc_to_au], [(w - proj(arr)[1]) * pc_to_au],
                       s=90, c=col, marker=mk, label=label, zorder=5,
                       edgecolors="none")

        # The uncertainty disc at arrival, to scale, in au.
        circ = plt.Circle((0, 0), unc_pc * pc_to_au, fill=False,
                          edgecolor=theme.accent, ls="--", lw=1.0, alpha=0.8)
        ax.add_patch(circ)
        ax.annotate("interception uncertainty\n%.1f au" % (unc_pc * pc_to_au),
                    xy=(0, unc_pc * pc_to_au), xytext=(6, 8),
                    textcoords="offset points", fontsize=6.8,
                    color=theme.accent)

        # The miss vector: aiming stale versus intercepting.
        su, sw = proj(stale)
        au_ = proj(arr)
        ax.annotate("", xy=(0, 0),
                    xytext=((su - au_[0]) * pc_to_au, (sw - au_[1]) * pc_to_au),
                    arrowprops=dict(arrowstyle="->", color=theme.accent_two,
                                    lw=1.4))
        ax.set_title("%s\nGaia DR3 %s" % (name, sid), fontsize=9.5)
        ax.set_xlabel("radial offset from interception [au]")
        ax.set_ylabel("transverse offset [au]")
        ax.grid(True, alpha=theme.grid_alpha)
        ax.set_aspect("equal", adjustable="datalim")
        ax.legend(fontsize=6.6, facecolor=theme.background,
                  edgecolor=theme.grid, labelcolor=theme.foreground,
                  loc="lower left")
        rnd.stamp(ax, theme, lines=[
            "flight (coord)  %.4f yr" % flight,
            "arrival epoch   J%.4f" % arrival_epoch,
            "stale miss      %.1f au" % err_au,
            "arrival sigma   %.1f au" % (unc_pc * pc_to_au),
        ], loc="upper right")
        plotted.append([name, sid, err_au, unc_pc * pc_to_au, flight,
                        arrival_epoch])

    fig.suptitle("Aiming at a moving destination — the catalogue position is "
                 "not the target", fontsize=11)
    rnd.scope_note(fig, theme,
                   "Special-relativistic kinematics. The star moves during the "
                   "flight; the interception point is where the two worldlines "
                   "meet, not where the catalogue says the star is.")

    prov = meta.figure_provenance(
        color_model="illustrative_only",
        render_parameters=dict(theme.describe(),
                               projection="onto the plane of the three saved "
                                          "positions; a projection of saved "
                                          "coordinates, not a recomputation",
                               transverse_units="au"),
        observer_frame="BarycentricCartesian")

    for ext in ("png", "svg"):
        pth = os.path.join(out_dir, "destination_interception." + ext)
        fig.savefig(pth, facecolor=fig.get_facecolor())
        records.write_sidecar(pth, prov)
        print("wrote %s" % pth)
    plt.close(fig)

    records.write_plotted_csv(
        os.path.join(out_dir, "destination_interception.plotted.csv"),
        ["destination", "source_id", "stale_targeting_error_au",
         "interception_uncertainty_au", "flight_coordinate_time_yr",
         "arrival_epoch"], plotted)
    return 0


if __name__ == "__main__":
    sys.exit(main())
