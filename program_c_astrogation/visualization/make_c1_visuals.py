#!/usr/bin/env python3
"""Generate every C1 visual from saved solver records.

    python visualization/make_c1_visuals.py \
        --results results/verification/program_c1

Runs each plotting and animation script in turn with explicit paths. This is a
DRIVER: it computes nothing and draws nothing itself, it only sequences the
scripts that read the solver's CSVs.

Run c1_verify first. Every script here fails loudly if its input is absent,
because a visual generated from a stale or missing record is worse than no
visual.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def run(script, *args):
    cmd = [sys.executable, os.path.join(HERE, script)] + list(args)
    print("--- %s" % script)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stdout + r.stderr)
        raise SystemExit("%s failed with %d" % (script, r.returncode))
    for line in r.stdout.splitlines():
        if line.startswith(("wrote", "  ")):
            print(line)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True)
    ap.add_argument("--quick", action="store_true",
                    help="fewer animation frames, for a fast rebuild")
    a = ap.parse_args(argv)

    res = os.path.abspath(a.results)
    data = os.path.join(res, "data")
    meta = os.path.join(res, "c1_metadata.json")
    figs = os.path.join(res, "figures")
    anims = os.path.join(res, "animations")
    sheets = os.path.join(res, "contact_sheets")
    if not os.path.exists(meta):
        raise SystemExit("%s is absent: run c1_verify before generating visuals"
                         % meta)

    run("static/plot_local_neighborhood.py",
        os.path.join(data, "local_neighborhood.csv"), "--metadata", meta,
        "--out", figs)
    run("static/plot_worldline.py",
        os.path.join(data, "worldline_proxima_1g.csv"), "--metadata", meta,
        "--out", figs)
    run("static/plot_destination_interception.py",
        os.path.join(data, "interception.csv"), "--metadata", meta,
        "--out", figs)
    run("static/plot_time_comparison.py",
        os.path.join(data, "trajectory_comparison.csv"), "--metadata", meta,
        "--out", figs)
    run("static/render_apparent_sky.py",
        os.path.join(data, "apparent_sky_earth_frame.csv"), "--metadata", meta,
        "--out", figs, "--name", "sky_earth_frame")
    run("static/render_apparent_sky.py",
        os.path.join(data, "apparent_sky_ship_frame.csv"), "--metadata", meta,
        "--out", figs, "--name", "sky_ship_frame", "--show-shift")

    nf = ["--frames", "48"] if a.quick else []
    run("animation/animate_journey.py",
        os.path.join(data, "worldline_proxima_1g.csv"), "--metadata", meta,
        "--out", anims, "--contact-sheets", sheets, *nf)
    run("animation/animate_star_neighborhood.py",
        os.path.join(data, "local_neighborhood.csv"), "--metadata", meta,
        "--out", anims, "--contact-sheets", sheets, *nf)
    run("animation/animate_relativistic_sky.py",
        os.path.join(data, "apparent_sky_sequence.csv"), "--metadata", meta,
        "--out", anims, "--contact-sheets", sheets)

    # One consolidated manifest beside the per-family ones, since the results
    # layout calls for a single animation_manifest.json.
    import json
    entries, prov = [], {}
    for m in ("journey_manifest.json", "sky_manifest.json",
              "neighborhood_manifest.json"):
        p = os.path.join(anims, m)
        if os.path.exists(p):
            with open(p, encoding="utf-8") as fh:
                j = json.load(fh)
            entries.extend(j["animations"])
            prov = j.get("provenance", prov)
    out = os.path.join(res, "animation_manifest.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump({"animations": entries, "provenance": prov}, fh, indent=2,
                  sort_keys=True)
    print("wrote %s (%d entries)" % (out, len(entries)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
