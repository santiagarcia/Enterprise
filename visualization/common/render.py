"""Rendering helpers: deterministic frames, GIF/MP4 encoding, contact sheets.

DETERMINISM IS THE POINT. Animation frames must be reproducible, so:

  * the camera path is a closed-form function of the frame index, never
    incremental state that accumulates rounding;
  * axis limits are computed once from the whole dataset and then frozen, so no
    frame autoscales;
  * exposure is fixed across frames unless exposure adaptation is explicitly the
    subject of the animation;
  * labels sit at positions derived from the data, never at random offsets;
  * star identities are stable because rows are keyed by source_id and drawn in a
    fixed, sorted order.

GIF falls back gracefully when no encoder is available. FRAME GENERATION AND PNG
OUTPUT NEVER DO: those are the deterministic products the gate checks, and a
missing video encoder must not be able to weaken them.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from typing import Any, Dict, List, Optional, Sequence

from . import records


# --- Encoder discovery ------------------------------------------------------

def find_ffmpeg() -> Optional[str]:
    """Locate an ffmpeg binary, preferring one already on PATH.

    imageio_ffmpeg ships a self-contained binary, which is why it is an
    acceptable fallback: it adds no system-level dependency.
    """
    exe = shutil.which("ffmpeg")
    if exe:
        return exe
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except Exception:                      # noqa: BLE001 -- optional
        return None


# --- Frame writing ----------------------------------------------------------

def save_frames(fig_factory, n_frames: int, frame_dir: str, *,
                prefix: str = "frame", dpi: int = 120) -> List[str]:
    """Render n_frames PNGs by calling fig_factory(i) -> matplotlib Figure.

    Frames go to a scratch directory that is NOT committed. The requirement is
    explicit: do not commit thousands of raw animation frames. What survives is
    the GIF/MP4 and a small contact sheet.
    """
    import matplotlib.pyplot as plt

    out = records.ensure_dir(frame_dir)
    paths = []
    for i in range(n_frames):
        fig = fig_factory(i)
        p = os.path.join(out, "%s_%05d.png" % (prefix, i))
        fig.savefig(p, dpi=dpi, facecolor=fig.get_facecolor())
        plt.close(fig)
        paths.append(p)
    return paths


def file_sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


# --- Encoding ---------------------------------------------------------------

def write_gif(frames: Sequence[str], out_path: str, fps: int = 20,
              max_frames: int = 48, max_width: int = 480,
              colors: int = 64) -> Optional[str]:
    """A short preview loop, via Pillow. Returns None if it cannot be written.

    A GIF here is a PREVIEW, and the MP4 is the full-quality product. Left at
    full resolution and full frame count a single rotation loop came to 4.3 MB,
    which is not a preview -- it is a large file nobody scrubs, committed forever.
    So: subsampled to max_frames, downscaled to max_width, quantised to a fixed
    palette size.

    Every reduction is deterministic -- a fixed stride, a fixed target width, a
    fixed palette size, no dithering -- so two runs over the same frames produce
    the same bytes.
    """
    try:
        from PIL import Image
    except ImportError:
        return None
    if not frames:
        return None
    stride = max(1, (len(frames) + max_frames - 1) // max_frames)
    picked = frames[::stride]
    try:
        imgs = []
        for p in picked:
            im = Image.open(p).convert("RGB")
            if im.width > max_width:
                h = max(1, round(im.height * max_width / im.width))
                im = im.resize((max_width, h), Image.LANCZOS)
            # Dither.NONE matters for determinism AND for size: error diffusion
            # scatters pixels that would otherwise compress as flat runs.
            imgs.append(im.convert("P", palette=Image.ADAPTIVE, colors=colors,
                                   dither=Image.Dither.NONE))
        imgs[0].save(records.resolve(out_path), save_all=True,
                     append_images=imgs[1:],
                     duration=int(1000 * stride / fps), loop=0, optimize=True)
        return records.resolve(out_path)
    except Exception:                      # noqa: BLE001 -- graceful fallback
        return None


def write_mp4(frames: Sequence[str], out_path: str, fps: int = 24,
              frame_dir: Optional[str] = None) -> Optional[str]:
    """The primary full-quality animation. Returns None if no encoder exists."""
    exe = find_ffmpeg()
    if not exe or not frames:
        return None
    pattern = os.path.join(frame_dir or os.path.dirname(frames[0]),
                           os.path.basename(frames[0])[:-9] + "%05d.png")
    out = records.resolve(out_path)
    cmd = [
        exe, "-y", "-loglevel", "error",
        "-framerate", str(fps), "-i", pattern,
        "-c:v", "libx264", "-pix_fmt", "yuv420p",
        # -crf fixed rather than a bitrate target, so two runs of the same frames
        # encode identically rather than drifting with a rate-control search.
        "-crf", "18", "-preset", "medium",
        # H.264 needs even dimensions.
        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        out,
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
        if r.returncode != 0 or not os.path.exists(out):
            return None
        return out
    except Exception:                      # noqa: BLE001 -- graceful fallback
        return None


# --- Contact sheet ----------------------------------------------------------

def contact_sheet(frames: Sequence[str], out_path: str, *, tiles: int = 9,
                  columns: int = 3, label_prefix: str = "frame") -> Optional[str]:
    """A grid of 6-12 representative frames.

    This is how a regression is spotted without opening a video, and how a
    reviewer checks an animation in a diff. Frames are picked at a fixed even
    stride, never at random, so the same input always yields the same sheet.
    """
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        return None
    if not frames:
        return None
    tiles = max(6, min(12, tiles))
    idx = [round(i * (len(frames) - 1) / (tiles - 1)) for i in range(tiles)] \
        if len(frames) > 1 else [0]
    picked = [frames[i] for i in idx]

    ims = [Image.open(p).convert("RGB") for p in picked]
    w, h = ims[0].size
    scale = min(1.0, 480.0 / max(w, 1))
    tw, th = max(1, int(w * scale)), max(1, int(h * scale))
    rows = (len(ims) + columns - 1) // columns
    sheet = Image.new("RGB", (columns * tw, rows * th), (8, 10, 16))
    draw = ImageDraw.Draw(sheet)
    for k, im in enumerate(ims):
        im = im.resize((tw, th), Image.LANCZOS)
        x, y = (k % columns) * tw, (k // columns) * th
        sheet.paste(im, (x, y))
        draw.text((x + 6, y + 6), "%s %d" % (label_prefix, idx[k]),
                  fill=(235, 240, 250))
    out = records.resolve(out_path)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    sheet.save(out)
    return out


def frame_is_blank(path: str, min_distinct: int = 12) -> bool:
    """True if a frame carries essentially no information.

    A blank frame is the classic silent animation failure: the pipeline runs, the
    file exists, the video plays, and one second of it is empty. Counting
    distinct colours is a cheap, robust detector -- a real plot has hundreds.
    """
    try:
        from PIL import Image
    except ImportError:
        return False
    im = Image.open(path).convert("RGB")
    colors = im.getcolors(maxcolors=1 << 20)
    return colors is not None and len(colors) < min_distinct


# --- Animation manifest -----------------------------------------------------

def animation_entry(name: str, path: Optional[str], *, codec: str,
                    resolution: Sequence[int], fps: int, frames: Sequence[str],
                    duration_s: float, source_csv: str) -> Dict[str, Any]:
    """One record for animation_manifest.json.

    Carries the source data checksum plus first- and final-frame checksums, so a
    regression in either the input or the rendering is visible without decoding
    the video.
    """
    return {
        "filename": os.path.basename(path) if path else None,
        "written": path is not None,
        "codec": codec,
        "resolution": list(resolution),
        "frame_rate": fps,
        "frame_count": len(frames),
        "duration_s": round(duration_s, 6),
        "source_data": os.path.basename(source_csv),
        "source_data_sha256": file_sha256(source_csv),
        "first_frame_sha256": file_sha256(frames[0]) if frames else None,
        "final_frame_sha256": file_sha256(frames[-1]) if frames else None,
        "note": None if path else
                "encoder unavailable; frames and PNG output are unaffected and "
                "remain deterministic",
    }


def write_manifest(path: str, entries: List[Dict[str, Any]],
                   provenance: Dict[str, Any]) -> str:
    out = records.resolve(path)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump({"animations": entries, "provenance": provenance}, fh,
                  indent=2, sort_keys=True)
    return out


# --- Overlays ---------------------------------------------------------------

def stamp(ax, theme, *, lines: Sequence[str], loc: str = "upper left") -> None:
    """The on-frame overlay: timestamp, clocks, beta/gamma, catalogue, target.

    Position is FIXED by `loc`, never chosen to avoid overlap, because a label
    that moves between frames reads as motion in the data.
    """
    x, y, ha, va = {
        "upper left": (0.015, 0.985, "left", "top"),
        "upper right": (0.985, 0.985, "right", "top"),
        "lower left": (0.015, 0.015, "left", "bottom"),
        "lower right": (0.985, 0.015, "right", "bottom"),
    }[loc]
    # Axes3D.text takes (x, y, z, s) and places the label in DATA space; its
    # text2D is the axes-fraction version. Using the wrong one puts the overlay
    # somewhere in the scene, where it moves as the camera rotates -- which is
    # exactly the drifting-label failure the determinism rules forbid.
    draw = getattr(ax, "text2D", ax.text)
    draw(x, y, "\n".join(lines), transform=ax.transAxes, ha=ha, va=va,
         fontsize=7.5, family="monospace", color=theme.foreground,
         linespacing=1.45,
         bbox=dict(facecolor=theme.background, edgecolor=theme.grid,
                   alpha=0.72, boxstyle="round,pad=0.45", linewidth=0.6))


def figure_stamp(fig, theme, *, lines, x=0.012, y=0.5) -> None:
    """The overlay placed in FIGURE coordinates, outside the axes.

    Needed for projections such as Mollweide, where the data fills an ellipse and
    an in-axes box necessarily lands on top of the graticule labels. Position is
    fixed, so it does not move between frames.
    """
    fig.text(x, y, "\n".join(lines), ha="left", va="center", fontsize=7.2,
             family="monospace", color=theme.foreground, linespacing=1.5,
             bbox=dict(facecolor=theme.background, edgecolor=theme.grid,
                       alpha=0.85, boxstyle="round,pad=0.5", linewidth=0.6))


def scope_note(fig, theme, text: str) -> None:
    """The standing caveat, on every C1 visual."""
    fig.text(0.5, 0.005, text, ha="center", va="bottom", fontsize=6.8,
             color=theme.grid, family="monospace")
