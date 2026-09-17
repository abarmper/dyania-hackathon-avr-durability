"""Raster figures -> calibrated step curves.

NOTION Fig 3 (RGB JPEG): two panels sharing one x-axis row; TAVI red, SAVR blue; cumulative incidence.
Wakami Fig 1 (8-bit grayscale PNG): single panel; black curve; light-grey CI band; cumulative incidence.

Axis detection: dark-pixel projections give axis rows/columns; ticks are short dark runs just outside
the axes. Curve tracing: colour (or darkness) mask inside the plot box, 1-px dilation, largest connected
component, per-column median row, monotone enforcement, risers -> step vertices.
"""
from __future__ import annotations

import io
from dataclasses import dataclass

import numpy as np
from pypdf import PdfReader
from scipy import ndimage

from vector_extract import CurveData


@dataclass
class PixelPanel:
    x0: int          # y-axis column
    x1: int          # right end of the x-axis run
    y_axis_row: int  # x-axis row (value 0)
    y_top: int       # top of plot box
    x_ticks: list[int]
    y_ticks: list[int]


def load_image(pdf_path, page_index, xobject):
    reader = PdfReader(str(pdf_path))
    page = reader.pages[page_index]
    obj = page["/Resources"]["/XObject"][xobject].get_object()
    w, h = int(obj["/Width"]), int(obj["/Height"])
    data = obj.get_data()
    flt = str(obj.get("/Filter"))
    if "DCT" in flt:
        from PIL import Image
        return np.asarray(Image.open(io.BytesIO(data)).convert("RGB"))
    cs = str(obj.get("/ColorSpace"))
    nch = 1 if "Gray" in cs else 3
    arr = np.frombuffer(data, np.uint8)
    return arr[: w * h * nch].reshape(h, w, nch) if nch == 3 else arr[: w * h].reshape(h, w)


def _cluster(idx, gap=2):
    groups = []
    for v in idx:
        if groups and v - groups[-1][-1] <= gap:
            groups[-1].append(v)
        else:
            groups.append([v])
    return [float(np.mean(g)) for g in groups]


def _longest_run(v):
    """(start, end) of the longest True run in a 1-D boolean array (end inclusive)."""
    best = (0, -1); start = None
    for i, x in enumerate(list(v) + [False]):
        if x and start is None:
            start = i
        elif not x and start is not None:
            if i - 1 - start > best[1] - best[0]:
                best = (start, i - 1)
            start = None
    return best


def detect_panels(dark: np.ndarray, min_axis_frac=0.25, tick_len=6, gap=3):
    """dark: boolean array. Returns list of PixelPanel (left to right).

    A y-axis is a column whose longest continuous dark run (small gaps closed) spans >= min_axis_frac*H.
    Its x-axis is the row within [bottom - tick_len - 8, bottom] with the longest dark run to the right
    (the run must be >= 10% of the width). Rotated axis labels fail the second test."""
    H, W = dark.shape
    closed = ndimage.binary_closing(dark, structure=np.ones((1, 2 * gap + 1)))
    closed_v = ndimage.binary_closing(dark, structure=np.ones((2 * gap + 1, 1)))
    panels = []
    for cx in range(W):
        a, b = _longest_run(closed_v[:, cx])
        if b - a + 1 < min_axis_frac * H:
            continue
        if panels and cx - panels[-1].x0 <= 3:  # same axis, thicker than 1 px
            continue
        best = None
        for r in range(max(b - tick_len - 8, a), b + 1):
            run = closed[r, cx:]
            L = int(np.argmin(run)) if (~run).any() else W - cx
            if best is None or L > best[1]:
                best = (r, L)
        xr, L = best
        if L < 0.1 * W:
            continue
        end = cx + L - 1
        below = dark[xr + 2: xr + 2 + tick_len + 4, max(cx - 2, 0): end + 3].sum(0)
        xt = [int(round(max(cx - 2, 0) + v)) for v in _cluster(list(np.where(below >= tick_len)[0]), gap=2)]
        left = dark[a: xr + 1, max(cx - tick_len - 6, 0): cx - 1].sum(1)
        yt = [int(round(a + v)) for v in _cluster(list(np.where(left >= tick_len - 1)[0]), gap=2)]
        panels.append(PixelPanel(cx, end, xr, a, xt, yt))
    return panels


def trace_curve(mask: np.ndarray, panel: PixelPanel, rising=True):
    """Largest connected component of mask inside the plot box -> (cols, rows) of the curve."""
    sub = np.zeros_like(mask)
    sub[panel.y_top: panel.y_axis_row, panel.x0 + 1: panel.x1 + 1] = mask[panel.y_top: panel.y_axis_row, panel.x0 + 1: panel.x1 + 1]
    sub = ndimage.binary_dilation(sub, iterations=1)
    lab, n = ndimage.label(sub, structure=np.ones((3, 3)))
    if n == 0:
        raise ValueError("no curve pixels")
    sizes = ndimage.sum(sub, lab, range(1, n + 1))
    comp = lab == (int(np.argmax(sizes)) + 1)
    cols = np.where(comp.any(0))[0]
    rows = np.array([np.median(np.where(comp[:, c])[0]) for c in cols])
    # a rising CIF moves to smaller row indices; enforce monotone
    rows = np.minimum.accumulate(rows) if rising else np.maximum.accumulate(rows)
    return cols, rows


def trace_lowest(mask: np.ndarray, panel: PixelPanel, x_start: int, jump_px=30):
    """Lowest dark pixel per column inside the plot box (for the LOWER of two overlaid cumulative-incidence
    curves, e.g. SVD below all-cause death). Censor marks below the curve are removed by the monotone
    running minimum of the row index; tracing stops at the first upward jump > jump_px (end of the curve)."""
    cols, rows = [], []
    prev = None
    for c in range(x_start, panel.x1 + 1):
        idx = np.where(mask[panel.y_top: panel.y_axis_row, c])[0]
        if len(idx) == 0:
            continue
        low = panel.y_top + int(idx.max())
        if prev is not None and prev - low > jump_px:
            break
        cols.append(c); rows.append(low); prev = low
    cols = np.array(cols); rows = np.minimum.accumulate(np.array(rows, float))
    return cols, rows


def calibrate(panel: PixelPanel, x_max_years: float, y_tick_step_pct: float):
    """x: first tick = 0, last tick = x_max_years (the axis line may start left of the first tick).
    y: lowest tick = 0 (may sit above the x-axis line), tick spacing = y_tick_step_pct."""
    xt = sorted(panel.x_ticks); yt = sorted(panel.y_ticks)
    if len(xt) < 2 or len(yt) < 2:
        raise ValueError(f"need >= 2 ticks on each axis, got x={xt} y={yt}")
    x_origin = xt[0]; px_per_year = (xt[-1] - xt[0]) / x_max_years
    spacing = float(np.median(np.diff(yt))); px_per_pct = spacing / y_tick_step_pct
    # the zero tick may be hidden under the x-axis line: if the lowest detected tick sits about one
    # spacing above the axis, zero is the axis line itself
    y_zero = yt[-1] if (panel.y_axis_row - yt[-1]) < 0.5 * spacing else panel.y_axis_row
    return x_origin, px_per_year, y_zero, px_per_pct


def steps_from_trace(cols, rows, x_origin, px_per_year, px_per_pct, y_zero_row, min_rise_px=1.0):
    t = (cols - x_origin) / px_per_year
    cif = (y_zero_row - rows) / px_per_pct / 100.0
    T, C = [0.0], [0.0]
    last_row = y_zero_row
    for ti, ci, r in zip(t, cif, rows):
        if last_row - r >= min_rise_px:
            T.append(float(max(ti, 0.0))); C.append(float(max(ci, C[-1]))); last_row = r
    return np.array(T), np.clip(np.array(C), 0, 1)


def colour_mask(img, colour):
    r, g, b = img[..., 0].astype(int), img[..., 1].astype(int), img[..., 2].astype(int)
    if colour == "red":
        return (r > 120) & (g < 90) & (b < 90)
    if colour == "blue":
        return (b > 120) & (r < 90) & (g < 110)
    if colour == "black":
        return img.max(-1) < 110 if img.ndim == 3 else img < 60
    raise ValueError(colour)


def extract_notion_fig3(pdf_path, page_index=5, xobject="/Im0", x_max_years=10, y_tick_step_pct=10.0):
    import config as C
    img = load_image(pdf_path, page_index, xobject)
    # axis detection on the MIN channel: saturated curve colours drawn over the axis count as dark
    panels = detect_panels(img.min(-1) < 110)
    if len(panels) < 2:
        raise ValueError(f"expected 2 panels, found {len(panels)}")
    panels = sorted(panels, key=lambda p: p.x0)[:2]
    out = {}
    specs = [s for s in C.NOTION if s.figure == "Fig 3"]
    for spec in specs:
        p = panels[spec.panel_index]
        x_origin, px_year, y_zero, px_pct = calibrate(p, x_max_years, y_tick_step_pct)
        mask = colour_mask(img, spec.colour_hint)
        cols, rows = trace_curve(mask, p, rising=True)
        T, Cif = steps_from_trace(cols, rows, x_origin, px_year, px_pct, y_zero)
        t_end = float((cols.max() - x_origin) / px_year)
        out[spec.id] = CurveData(id=spec.id, figure="Fig 3", panel=spec.panel, name=spec.name, arm=spec.arm, kind="CIF",
                                 t=T, y=1.0 - Cif, censor_times=None, t_end=t_end,
                                 provenance={"page": page_index, "xobject": xobject, "colour": spec.colour_hint,
                                             "panel_px": (p.x0, p.x1, p.y_axis_row, p.y_top), "x_ticks_px": p.x_ticks, "y_ticks_px": p.y_ticks,
                                             "px_per_year": float(px_year), "px_per_pct": float(px_pct), "y_zero_row": int(y_zero), "n_curve_columns": int(len(cols))})
    return out


def extract_wakami_fig1(pdf_path, page_index=3, xobject="/Im0"):
    import config as C
    spec = [s for s in C.WAKAMI if s.figure == "Fig 1"][0]
    img = load_image(pdf_path, page_index, xobject)
    dark = img < 110
    panels = detect_panels(dark)
    p = sorted(panels, key=lambda p: -(p.x1 - p.x0))[0]
    x_max = spec.x_max_years or (len(p.x_ticks) - 1)
    y_step = spec.y_tick_step or 1.0
    x_origin, px_year, y_zero, px_pct = calibrate(p, x_max, y_step)
    mask = img < 130   # anti-aliased 1-px lines have core values ~75-100; text and marks are darker
    # the SVD curve is the lower of the two overlaid CIFs; the dashed death curve + censor marks form the
    # largest connected component, so trace the lowest dark pixel per column instead (min rise 3 px = 0.5 pp)
    cols, rows = trace_lowest(mask, p, x_start=x_origin)
    T, Cif = steps_from_trace(cols, rows, x_origin, px_year, px_pct, y_zero, min_rise_px=3.0)
    t_end = float((cols.max() - x_origin) / px_year)
    return {spec.id: CurveData(id=spec.id, figure="Fig 1", panel="single", name=spec.name, arm="all", kind="CIF", t=T, y=1.0 - Cif,
                               censor_times=None, t_end=t_end,
                               provenance={"page": page_index, "xobject": xobject, "panel_px": (p.x0, p.x1, p.y_axis_row, p.y_top),
                                           "x_ticks_px": p.x_ticks, "y_ticks_px": p.y_ticks, "px_per_year": float(px_year), "px_per_pct": float(px_pct), "y_zero_row": int(y_zero)})}
