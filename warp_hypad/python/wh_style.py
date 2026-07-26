"""Shared plotting style for warp_hypad figures.

One module so every figure in the project reads as one system, and so the
palette lives in exactly one place. Values are the validated reference palette
(light mode); the categorical order is fixed and never cycled.

Colour jobs used here:
  categorical  identity of a series (numeric vs analytic, one R value vs another)
  sequential   magnitude, one hue light->dark (determinant, |residual|)
  diverging    polarity about zero, two hues + neutral grey midpoint
               (errors, sensitivities -- anything that can be either sign)

A signed field ALWAYS gets a zero-centred diverging scale. Showing a signed
sensitivity map on a sequential ramp would put the visual midpoint at an
arbitrary value and make "positive" and "negative" regions indistinguishable,
which is the specific thing these maps exist to show.
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

# --- Reference palette, light mode ----------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
INK3 = "#8a8983"
GRID = "#e4e3df"

# Fixed categorical order. Assign in order, never cycle. Charts here use at
# most three, which is the all-pairs-validated prefix.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]

STATUS_GOOD = "#0ca30c"
STATUS_WARN = "#fab219"
STATUS_CRIT = "#d03b3b"

# Sequential: single hue, light -> dark.
SEQ_STEPS = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CMAP_SEQ = LinearSegmentedColormap.from_list("wh_seq", SEQ_STEPS)

# Diverging: blue <-> red with a neutral grey midpoint. Never a hue at the
# midpoint -- the middle must read as "nothing".
CMAP_DIV = LinearSegmentedColormap.from_list(
    "wh_div",
    ["#0d366b", "#256abf", "#3987e5", "#9ec5f4", "#f0efec",
     "#f0a09f", "#e34948", "#c22e2d", "#8f1f1e"],
)

FOOTER_UNITS = ("Geometrized units (G = c = 1). Metric parameters are "
                "dimensionless; coordinates share the grid's length unit.")


def style_axes(ax, logx=False, logy=False):
    ax.set_facecolor(SURFACE)
    ax.grid(True, which="major", color=GRID, linewidth=0.8, zorder=0)
    if logx or logy:
        ax.grid(True, which="minor", color=GRID, linewidth=0.4, alpha=0.6, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK3)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=8, width=0.8)


def title_block(fig, title, subtitle=None, y=0.975):
    fig.suptitle(title, fontsize=13, color=INK, x=0.01, ha="left", y=y)
    if subtitle:
        fig.text(0.01, y - 0.055, subtitle, fontsize=8.5, color=INK2, ha="left")


def footer(fig, text):
    fig.text(0.01, 0.008, text, fontsize=6.5, color=INK3, ha="left",
             va="bottom", linespacing=1.5)


def signed_norm(values):
    """Zero-centred diverging norm; None when the data is single-signed."""
    lo, hi = min(values), max(values)
    if lo >= 0 or hi <= 0:
        return None
    m = max(abs(lo), abs(hi))
    return TwoSlopeNorm(vmin=-m, vcenter=0.0, vmax=m)


def imshow_field(ax, xs, ys, grid, *, signed, label, cbar=True):
    """Draw a 2D field, choosing the colour job from the data's sign structure.

    `signed` only says the field is ALLOWED to change sign. Whether it actually
    does is a property of the data, and it decides the ramp:

      both signs present -> diverging, zero at the neutral grey midpoint
      one sign only      -> sequential in |value|, so zero reads as "nothing"

    Forcing a diverging ramp onto single-signed data is a real error, not a
    cosmetic one: a field spanning [0, 0.4] would render its zeros at one
    saturated end of the ramp, so "no sensitivity here" would look like an
    extreme value. Several of the A3 sensitivity panels are single-signed.
    """
    extent = [min(xs), max(xs), min(ys), max(ys)]
    flat = [v for row in grid for v in row]
    lo, hi = min(flat), max(flat)

    if signed and lo < 0.0 < hi:
        im = ax.imshow(grid, origin="lower", extent=extent, aspect="auto",
                       cmap=CMAP_DIV, norm=signed_norm(flat),
                       interpolation="nearest")
    elif hi <= 0.0:
        # All non-positive: light at zero, dark at the largest magnitude.
        im = ax.imshow(grid, origin="lower", extent=extent, aspect="auto",
                       cmap=CMAP_SEQ.reversed(), vmin=lo, vmax=0.0,
                       interpolation="nearest")
    else:
        # All non-negative: light at zero, dark at the largest magnitude.
        im = ax.imshow(grid, origin="lower", extent=extent, aspect="auto",
                       cmap=CMAP_SEQ, vmin=min(0.0, lo), vmax=hi,
                       interpolation="nearest")
    if cbar:
        cb = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(label, fontsize=7.5, color=INK2)
        cb.ax.tick_params(labelsize=7, colors=INK2)
        cb.outline.set_visible(False)
    return im


def save(fig, path, rect=(0, 0.05, 1, 0.90)):
    fig.tight_layout(rect=rect)
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")
