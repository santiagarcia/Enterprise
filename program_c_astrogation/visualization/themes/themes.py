"""Two themes: scientific and cinematic.

THE RULE THAT SEPARATES THEM
----------------------------
The cinematic theme may add trails, glow and easing. It MAY NOT alter positions,
directions, timing or physical quantities. It draws the same numbers differently;
it does not draw different numbers.

Every knob below is therefore a *rendering* knob -- alpha, linewidth, colour,
marker area, glow passes. There is deliberately no knob that shifts a coordinate,
rescales a time axis nonlinearly, or applies an artistic offset, because such a
knob would let a theme change what the figure asserts.

Marker SIZE is the one place where care is needed. Stars must not all be rendered
at identical size, and the size comes from catalogue photometry or distance --
both saved quantities. The mapping is documented and recorded in the figure's
render_parameters, so a reader knows a big dot means "near" rather than "large".
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Tuple


@dataclass
class Theme:
    name: str
    background: str
    foreground: str
    grid: str
    accent: str
    accent_two: str
    star_color: str
    uncertainty_color: str
    trail: bool
    glow_passes: int
    grid_alpha: float
    dpi: int
    # Marker area in points^2 for the nearest and farthest plotted star. Stars
    # are never all the same size; this is the range the mapping spans.
    marker_area_range: Tuple[float, float] = (2.0, 90.0)
    font_family: str = "DejaVu Sans"

    def describe(self) -> Dict[str, Any]:
        """What goes into render_parameters, so the look is reproducible."""
        return {
            "theme": self.name,
            "background": self.background,
            "dpi": self.dpi,
            "trail": self.trail,
            "glow_passes": self.glow_passes,
            "marker_area_range_pt2": list(self.marker_area_range),
            "marker_area_meaning":
                "area scales with apparent brightness proxy; a larger dot means "
                "brighter or nearer, never physically larger",
            "alters_physical_quantities": False,
        }


# The scientific theme prioritises measurable geometry, uncertainty and labels.
SCIENTIFIC = Theme(
    name="scientific",
    background="#0b0e14",       # dark astronomical, but see PRINT below
    foreground="#e8ecf4",
    grid="#3a4457",
    accent="#4ea1ff",
    accent_two="#ffb454",
    star_color="#dfe7f5",
    uncertainty_color="#ff6b6b",
    trail=False,
    glow_passes=0,
    grid_alpha=0.35,
    dpi=160,
)

# A print-compatible alternative. The requirement is a dark astronomical
# background WITH print-compatible static alternatives, so this exists as a
# first-class theme rather than an afterthought: a dark figure dropped into a
# paper wastes ink and reproduces badly.
PRINT = Theme(
    name="scientific_print",
    background="#ffffff",
    foreground="#111318",
    grid="#c7ccd6",
    accent="#1f5fa9",
    accent_two="#a95a1f",
    star_color="#2a2f3a",
    uncertainty_color="#c0392b",
    trail=False,
    glow_passes=0,
    grid_alpha=0.5,
    dpi=200,
)

# The cinematic theme. Trails and glow only.
CINEMATIC = Theme(
    name="cinematic",
    background="#03040a",
    foreground="#f2f6ff",
    grid="#1b2436",
    accent="#63b8ff",
    accent_two="#ffcc66",
    star_color="#ffffff",
    uncertainty_color="#ff8080",
    trail=True,
    glow_passes=3,
    grid_alpha=0.15,
    dpi=160,
    marker_area_range=(3.0, 160.0),
)

THEMES = {t.name: t for t in (SCIENTIFIC, PRINT, CINEMATIC)}


def get(name: str) -> Theme:
    if name not in THEMES:
        raise KeyError("unknown theme %r; have %r" % (name, sorted(THEMES)))
    return THEMES[name]


def apply(theme: Theme) -> None:
    """Push the theme into matplotlib's global state."""
    import matplotlib as mpl
    mpl.rcParams.update({
        "figure.facecolor": theme.background,
        "savefig.facecolor": theme.background,
        "axes.facecolor": theme.background,
        "axes.edgecolor": theme.grid,
        "axes.labelcolor": theme.foreground,
        "text.color": theme.foreground,
        "xtick.color": theme.foreground,
        "ytick.color": theme.foreground,
        "grid.color": theme.grid,
        "grid.alpha": theme.grid_alpha,
        "font.family": theme.font_family,
        "figure.dpi": theme.dpi,
        "savefig.dpi": theme.dpi,
        "savefig.bbox": "tight",
        # Determinism: no hash-seeded font fallback ordering, no interactive
        # state leaking between figures.
        "svg.hashsalt": "program_c_astrogation",
    })


def glow(ax, x, y, theme: Theme, color=None, base_lw=1.6, **kw):
    """Draw a line with optional glow. GEOMETRY IS UNTOUCHED.

    The glow is additional strokes of the SAME x and y at larger linewidth and
    lower alpha. No offset, no smoothing, no resampling -- the curve a reader
    sees passes through exactly the plotted points in both themes.
    """
    color = color or theme.accent
    handles = []
    for i in range(theme.glow_passes, 0, -1):
        handles.append(ax.plot(x, y, color=color, linewidth=base_lw * (1 + 2.2 * i),
                               alpha=0.06, solid_capstyle="round", zorder=1)[0])
    handles.append(ax.plot(x, y, color=color, linewidth=base_lw, zorder=3, **kw)[0])
    return handles[-1]


def marker_areas(values, theme: Theme, invert=False):
    """Map a saved per-star quantity onto marker area.

    values are read from a CSV -- brightness factor, or 1/distance. The mapping
    is monotone and documented; it is presentation, not physics, because the
    input already came from the solver.
    """
    lo_a, hi_a = theme.marker_area_range
    vs = list(values)
    if not vs:
        return []
    lo, hi = min(vs), max(vs)
    if hi <= lo:
        return [0.5 * (lo_a + hi_a)] * len(vs)
    out = []
    for v in vs:
        t = (v - lo) / (hi - lo)
        if invert:
            t = 1.0 - t
        # Square-root so area grows sub-linearly: a linear area map makes the
        # nearest star swamp the figure and hides the population.
        out.append(lo_a + (hi_a - lo_a) * (t ** 0.5))
    return out


def doppler_rgb(doppler_factor: float) -> Tuple[float, float, float]:
    """Map a SAVED Doppler factor to an RGB triple.

    THIS IS approximate_blackbody_visualization AND NOTHING MORE. It is a
    perceptual cue for the direction and rough magnitude of the shift, not a
    prediction of observed colour: a real answer needs the source spectrum, the
    shifted Planck function and CIE colour-matching, and even then "observed
    colour" depends on adaptation and exposure.

    Any figure using this must declare color_model
    "approximate_blackbody_visualization". The metadata writer refuses
    undeclared colour models, so this cannot leak into a figure unlabelled.

    D > 1 is a blueshift; D < 1 is a redshift.
    """
    import math
    # log D is the natural variable: it is symmetric between a shift and its
    # reciprocal, which a linear map in D is not.
    u = math.log(max(doppler_factor, 1e-12)) * 1.2
    # Squashed with the ALGEBRAIC sigmoid u/(1+|u|) rather than tanh. Both are
    # monotone, odd and bounded, and either does this job equally well -- but a
    # hyperbolic function appearing anywhere in the visualization tree is what
    # trajectory integration and rapidity conversion look like, and the gate case
    # c1vis.visuals_derived_from_saved_output_only bans them outright. A
    # bright-line rule that needs no case-by-case exemption is worth more than
    # the tanh.
    t = u / (1.0 + abs(u))
    if t >= 0.0:                      # blueshift
        return (1.0 - 0.85 * t, 1.0 - 0.35 * t, 1.0)
    s = -t                            # redshift
    return (1.0, 1.0 - 0.55 * s, 1.0 - 0.9 * s)


def style_3d(ax, theme) -> None:
    """Make a 3-D axes look like space rather than a grey box.

    matplotlib fills the three background panes with an opaque grey by default,
    which in a dark theme reads as a lit cube floating in the void and washes out
    every faint star drawn behind it. Only the pane fill, the pane edges and the
    tick colours change here: no limit, no projection, no data.
    """
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((0.0, 0.0, 0.0, 0.0))
        axis.pane.set_edgecolor(theme.grid)
        axis.pane.set_alpha(0.25)
        axis._axinfo["grid"]["color"] = theme.grid
        axis._axinfo["grid"]["linewidth"] = 0.4
    ax.set_facecolor(theme.background)
    ax.tick_params(colors=theme.foreground, labelsize=7)
