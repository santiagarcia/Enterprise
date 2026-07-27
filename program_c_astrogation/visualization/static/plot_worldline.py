#!/usr/bin/env python3
"""Worldline figures: clocks, velocity, and the light-cone ledger.

    python visualization/static/plot_worldline.py \
        results/verification/program_c1/data/worldline_proxima_1g.csv \
        --metadata results/verification/program_c1/c1_metadata.json \
        --out results/verification/program_c1/figures

Produces worldline_coordinate_time, worldline_proper_time,
velocity_gamma_acceleration and earth_ship_destination_lightcone.

READS COLUMNS AND DRAWS THEM. Every clock, beta, gamma and signal time was
computed and saved by c1_verify. Nothing here integrates anything.
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

SCOPE = ("Special-relativistic kinematics. No warp metric, no propulsion energy, "
         "no collision hazard. Every clock is named.")


def _save(fig, out_dir, name, prov, plotted=None, header=None):
    written = []
    for ext in ("png", "svg"):
        p = os.path.join(out_dir, "%s.%s" % (name, ext))
        fig.savefig(p, facecolor=fig.get_facecolor())
        records.write_sidecar(p, prov)
        written.append(p)
    plt.close(fig)
    if plotted is not None:
        records.write_plotted_csv(
            os.path.join(out_dir, "%s.plotted.csv" % name), header, plotted)
    for p in written:
        print("wrote %s" % p)


def main(argv=None) -> int:
    p = cli.parser(__doc__.splitlines()[0])
    p.add_argument("--light-cones", default=None,
                   help="path to light_cones.csv; defaults to a sibling of the "
                        "worldline file, resolved absolutely either way")
    a = p.parse_args(argv)

    rows = records.read_csv(a.input)
    meta = records.RunMetadata.load(a.metadata)
    out_dir = records.ensure_dir(a.out)
    theme = themes.get(a.theme)
    themes.apply(theme)

    tau = records.column(rows, "proper_time_yr")
    t = records.column(rows, "coordinate_time_yr")
    beta = records.column(rows, "beta")
    gamma = records.column(rows, "gamma")
    accel = records.column(rows, "accel_pc_yr2")
    remaining = records.column(rows, "remaining_distance_pc")
    xs = records.column(rows, "x_pc")
    ys = records.column(rows, "y_pc")
    zs = records.column(rows, "z_pc")
    # Distance travelled along the (straight) coordinate path. This is a NORM of
    # saved coordinates, not a recomputation of the trajectory.
    dist = [(xs[i] ** 2 + ys[i] ** 2 + zs[i] ** 2) ** 0.5 for i in range(len(xs))]

    dest = meta.destination_name(0)
    sid = meta.destination_source_id(0)
    # Built once, here: figure_provenance() validates the colour model and fills
    # every required field, so a figure cannot be written with a partial sidecar.
    base_prov = meta.figure_provenance(
        color_model="illustrative_only",
        render_parameters=dict(theme.describe(),
                               line_meaning="all curves are saved columns"),
        observer_frame="BarycentricCartesian for coordinate quantities; ship "
                       "proper time is the crew clock along the same worldline")

    # --- 1. Distance against COORDINATE time --------------------------------
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    themes.glow(ax, t, dist, theme, color=theme.accent, base_lw=2.0,
                label="ship, coordinate time")
    ax.plot(t, [dist[-1]] * len(t), color=theme.grid, lw=0.8, ls=":",
            label="interception distance")
    ax.set_xlabel("coordinate time since departure [yr]  (solar-system frame)")
    ax.set_ylabel("distance from the Sun [pc]")
    ax.set_title("Worldline against COORDINATE time — %s (Gaia DR3 %s)"
                 % (dest, sid), fontsize=10)
    ax.grid(True, alpha=theme.grid_alpha)
    ax.legend(fontsize=7.5, facecolor=theme.background, edgecolor=theme.grid,
              labelcolor=theme.foreground)
    rnd.stamp(ax, theme, lines=[
        "clock       coordinate (solar-system barycentric)",
        "elapsed     %.4f yr" % t[-1],
        "ship clock  %.4f yr  (NOT the same clock)" % tau[-1],
        "model       %s" % meta.get("trajectory_model", "?"),
        "catalogue   %s" % meta.get("catalogue_release", "?"),
    ], loc="lower right")
    rnd.scope_note(fig, theme, SCOPE)
    _save(fig, out_dir, "worldline_coordinate_time", base_prov,
          [[t[i], dist[i], remaining[i]] for i in range(len(t))],
          ["coordinate_time_yr", "distance_pc", "remaining_distance_pc"])

    # --- 2. The two clocks, side by side ------------------------------------
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    themes.glow(ax, t, tau, theme, color=theme.accent_two, base_lw=2.0,
                label="ship proper time")
    ax.plot(t, t, color=theme.grid, lw=1.0, ls="--", label="coordinate time")
    ax.fill_between(t, tau, t, color=theme.accent, alpha=0.10,
                    label="time dilation gap")
    ax.set_xlabel("coordinate time since departure [yr]")
    ax.set_ylabel("elapsed time [yr]")
    ax.set_title("Two clocks, not one — %s (Gaia DR3 %s)" % (dest, sid),
                 fontsize=10)
    ax.grid(True, alpha=theme.grid_alpha)
    ax.legend(fontsize=7.5, facecolor=theme.background, edgecolor=theme.grid,
              labelcolor=theme.foreground)
    rnd.stamp(ax, theme, lines=[
        "coordinate  %.4f yr" % t[-1],
        "ship proper %.4f yr" % tau[-1],
        "difference  %.4f yr" % (t[-1] - tau[-1]),
        "ratio       %.6f" % (t[-1] / tau[-1] if tau[-1] else float("nan")),
        'the phrase "travel time" is not used: it does not name a clock',
    ], loc="upper left")
    rnd.scope_note(fig, theme, SCOPE)
    _save(fig, out_dir, "worldline_proper_time", base_prov,
          [[t[i], tau[i], t[i] - tau[i]] for i in range(len(t))],
          ["coordinate_time_yr", "proper_time_yr", "dilation_gap_yr"])

    # --- 3. Velocity, gamma, acceleration -----------------------------------
    fig, axes = plt.subplots(3, 1, figsize=(8.4, 8.0), sharex=True)
    axes[0].plot(tau, beta, color=theme.accent, lw=1.8)
    axes[0].set_ylabel(r"$\beta = v/c$")
    axes[0].set_title("Velocity, Lorentz factor and proper acceleration — %s"
                      % dest, fontsize=10)
    axes[1].plot(tau, gamma, color=theme.accent_two, lw=1.8)
    axes[1].set_ylabel(r"$\gamma$")
    axes[2].plot(tau, accel, color=theme.star_color, lw=1.4)
    axes[2].set_ylabel(r"$|a|$ [pc yr$^{-2}$]")
    axes[2].set_xlabel("ship proper time [yr]")
    for ax_ in axes:
        ax_.grid(True, alpha=theme.grid_alpha)
    axes[0].axhline(1.0, color=theme.uncertainty_color, lw=0.8, ls=":")
    axes[0].text(tau[0], 1.0, " c", color=theme.uncertainty_color, fontsize=7,
                 va="bottom")
    rnd.stamp(axes[0], theme, lines=[
        "peak beta   %.6f" % max(beta),
        "peak gamma  %.6f" % max(gamma),
        "accel       constant in magnitude, reversed at the midpoint",
        "x axis      SHIP proper time",
    ], loc="upper left")
    rnd.scope_note(fig, theme, SCOPE)
    _save(fig, out_dir, "velocity_gamma_acceleration", base_prov,
          [[tau[i], beta[i], gamma[i], accel[i]] for i in range(len(tau))],
          ["proper_time_yr", "beta", "gamma", "accel_pc_yr2"])

    # --- 4. The communication light cones -----------------------------------
    lc_path = a.light_cones or os.path.join(
        os.path.dirname(records.resolve(a.input)), "light_cones.csv")
    lc = records.read_csv(lc_path)
    lct = records.column(lc, "coordinate_time_yr")
    seen = records.column(lc, "earth_time_seen")
    age = records.column(lc, "earth_signal_age_yr")
    to_earth = records.column(lc, "signal_arrives_at_earth")
    epoch0 = float(meta.get("simulation_epoch") or 0.0)

    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    ax.plot(lct, [s - epoch0 for s in seen], color=theme.accent, lw=1.9,
            label="Earth clock the ship can SEE")
    ax.plot(lct, lct, color=theme.grid, lw=1.0, ls="--",
            label="Earth coordinate time")
    ax.plot(lct, [s - epoch0 for s in to_earth], color=theme.accent_two, lw=1.6,
            label="arrival at Earth of a signal sent now")
    ax.fill_between(lct, [s - epoch0 for s in seen], lct,
                    color=theme.uncertainty_color, alpha=0.12,
                    label="information lag")
    ax.set_xlabel("ship coordinate time since departure [yr]")
    ax.set_ylabel("Earth time since departure [yr]")
    ax.set_title("Communication light cones — Earth, ship, %s" % dest,
                 fontsize=10)
    ax.grid(True, alpha=theme.grid_alpha)
    ax.legend(fontsize=7.5, facecolor=theme.background, edgecolor=theme.grid,
              labelcolor=theme.foreground, loc="upper left")
    rnd.stamp(ax, theme, lines=[
        "max signal age    %.4f yr" % max(age),
        "round trip at end %.4f yr" % (2.0 * age[-1]),
        "Earth's CURRENT time is never shown aboard the ship:",
        "  only the latest event whose light has arrived",
    ], loc="lower right")
    rnd.scope_note(fig, theme, SCOPE)
    _save(fig, out_dir, "earth_ship_destination_lightcone", base_prov,
          [[lct[i], seen[i], age[i], to_earth[i]] for i in range(len(lct))],
          ["coordinate_time_yr", "earth_time_seen", "earth_signal_age_yr",
           "signal_arrives_at_earth"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
