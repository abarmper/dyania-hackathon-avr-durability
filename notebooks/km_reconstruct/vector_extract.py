"""Turn PageGraphics (vector figures) into calibrated step curves (CurveData).

Handles the two Kermen 2022 pages:
  page 4  Fig 2A/2B  stroked KM step curves, black tick segments, '%' y-labels, legend line samples
  page 7  Fig 4A     filled CIF polygons, thin tick rectangles, '%' y-labels, legend swatches
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

from pdf_vector import PageGraphics, Polyline, Segment, TextToken, parse_page

BLACK_PREFIX = ("SCN", "scn", "G", "g")


@dataclass
class Axis:
    origin: float          # page coordinate of value 0
    scale: float           # page units per data unit
    ticks: list[float] = field(default_factory=list)
    tick_values: list[float] = field(default_factory=list)

    def to_data(self, v):
        return (np.asarray(v, float) - self.origin) / self.scale

    def to_page(self, d):
        return self.origin + np.asarray(d, float) * self.scale


@dataclass
class Panel:
    letter: str
    x_axis_y: float
    y_axis_x: float
    x_max: float
    y_max: float
    xaxis: Axis | None = None
    yaxis: Axis | None = None

    def contains(self, x, y, pad=3.0):
        return (self.y_axis_x - pad <= x <= self.x_max + pad) and (self.x_axis_y - pad <= y <= self.y_max + pad)


@dataclass
class CurveData:
    id: str
    figure: str
    panel: str
    name: str
    arm: str
    kind: str                      # "KM" | "CIF"
    t: np.ndarray                  # step vertices (data units, years)
    y: np.ndarray                  # survival-type value after the jump at t (fraction); 1-CIF for CIF
    censor_times: np.ndarray | None
    t_end: float
    provenance: dict = field(default_factory=dict)


# ---------------------------------------------------------------- helpers

def is_black(colour) -> bool:
    return colour and colour[0] in BLACK_PREFIX


def colour_key(colour) -> str:
    return str(colour)


def _fit_axis(positions, values) -> Axis:
    positions = np.asarray(positions, float); values = np.asarray(values, float)
    if len(positions) == 1:
        raise ValueError("need >= 2 labelled ticks")
    b, a = np.polyfit(values, positions, 1)      # position = a + b*value
    resid = np.max(np.abs(a + b * values - positions))
    if resid > 0.5:
        raise ValueError(f"axis fit residual {resid:.2f} pt too large")
    return Axis(origin=float(a), scale=float(b), ticks=list(positions), tick_values=list(values))


def percent_tokens(text: list[TextToken]):
    out = []
    for tk in text:
        for run, rx in tk.runs:
            m = re.fullmatch(r"(\d+)%", run)
            if m:
                out.append((float(m.group(1)), rx, tk.y, tk.size))
    return out


def numeric_runs(tk: TextToken):
    """Numbers inside a token, with approximate x (comma decimals allowed)."""
    out = []
    for run, rx in tk.runs:
        for m in re.finditer(r"\d+(?:[.,]\d+)?", run):
            out.append((float(m.group(0).replace(",", ".")), rx + 0.5 * tk.size * m.start()))
    return out


# ---------------------------------------------------------------- panels & axes

def find_panels_from_labels(g: PageGraphics, label_offset=2.5) -> list[Panel]:
    """Group '%' labels into vertical runs (one per panel); return Panel with y-axis calibration."""
    pts = percent_tokens(g.text)
    if not pts:
        raise ValueError("no % labels found")
    # cluster by rounded x of the label's right edge (labels are right-aligned); use y-run continuity instead
    pts.sort(key=lambda p: p[2])
    runs: list[list] = []
    for p in pts:
        placed = False
        for r in runs:
            # same run if x within 12 pt and y gap plausible (< 60 pt)
            if abs(r[-1][1] - p[1]) < 15 and 0 < p[2] - r[-1][2] < 60:
                r.append(p); placed = True; break
        if not placed:
            runs.append([p])
    runs = [r for r in runs if len(r) >= 3]
    panels = []
    for r in sorted(runs, key=lambda r: -r[0][2]):   # top panel first (larger y)
        vals = [p[0] for p in r]; ys = [p[2] + label_offset for p in r]
        yaxis = _fit_axis(ys, vals)
        panels.append(Panel(letter="?", x_axis_y=float(yaxis.to_page(0.0)), y_axis_x=float("nan"),
                            x_max=float("nan"), y_max=float(yaxis.to_page(max(vals))), yaxis=yaxis))
    return panels


def attach_frames_and_letters(g: PageGraphics, panels: list[Panel]):
    """Find the y-axis x and x-extent from long black lines/rects near each panel's x-axis y."""
    long_h = [(min(s.x0, s.x1), max(s.x0, s.x1), s.y0) for s in g.segments if is_black(s.stroke) and s.is_horizontal and s.length > 100 and not s.dash]
    long_v = [(s.x0, min(s.y0, s.y1), max(s.y0, s.y1)) for s in g.segments if is_black(s.stroke) and s.is_vertical and s.length > 60 and not s.dash]
    # rect-drawn axes (page 7): thin long rects
    for r in g.rects:
        if r.mode == "f" and (r.h <= 1.2 and r.w > 100):
            long_h.append((min(r.x0, r.x1), max(r.x0, r.x1), (r.y0 + r.y1) / 2))
        if r.mode == "f" and (r.w <= 1.2 and r.h > 60):
            long_v.append(((r.x0 + r.x1) / 2, min(r.y0, r.y1), max(r.y0, r.y1)))
    letters = sorted([(tk.text, tk.x, tk.y) for tk in g.text if re.fullmatch(r"[A-D]", tk.text) and tk.size >= 9],
                     key=lambda l: -l[2])
    vt, ht = tick_marks(g)
    for p in panels:
        cands = [h for h in long_h if abs(h[2] - p.x_axis_y) < 4]
        if cands:
            h = max(cands, key=lambda h: h[1] - h[0])
            p.y_axis_x = h[0]; p.x_max = h[1]
            p.x_axis_y = h[2]
        else:
            # axes drawn as short rectangles (page 7): y-axis x = right end of the horizontal y-ticks of this panel
            yt = [t for t in ht if p.x_axis_y - 2 <= t[0] <= p.y_max + 2]
            if yt:
                p.y_axis_x = float(np.median([t[2] for t in yt]))
            xt = [t for t in vt if t[2] <= p.x_axis_y + 0.5 and t[1] >= p.x_axis_y - 6 and t[0] >= p.y_axis_x - 3]
            if xt:
                xs = sorted(t[0] for t in xt)
                p.x_max = xs[-1] + (np.median(np.diff(xs)) if len(xs) > 1 else 20.0)
        # y_max: top of this panel's own y-axis (highest label) plus a small margin
        p.y_max = float(p.yaxis.to_page(max(p.yaxis.tick_values))) + 3.0
    # letters by vertical order (top panel gets the first letter)
    for p, l in zip(sorted(panels, key=lambda p: -p.x_axis_y), letters):
        p.letter = l[0]
    return panels


def tick_marks(g: PageGraphics):
    """Short black axis-aligned marks: segments (page 4) or thin rects (page 7). Returns (vertical, horizontal)."""
    v, h = [], []
    for s in g.segments:
        if is_black(s.stroke) and 1.5 <= s.length <= 4.5:
            if s.is_vertical:
                v.append(((s.x0 + s.x1) / 2, min(s.y0, s.y1), max(s.y0, s.y1)))
            elif s.is_horizontal:
                h.append(((s.y0 + s.y1) / 2, min(s.x0, s.x1), max(s.x0, s.x1)))
    for r in g.rects:
        if r.mode == "f" and is_black(r.fill):
            if r.w <= 1.0 and 1.5 <= r.h <= 4.5:
                v.append(((r.x0 + r.x1) / 2, min(r.y0, r.y1), max(r.y0, r.y1)))
            elif r.h <= 1.0 and 1.5 <= r.w <= 4.5:
                h.append(((r.y0 + r.y1) / 2, min(r.x0, r.x1), max(r.x0, r.x1)))
    return v, h


def calibrate_x(g: PageGraphics, p: Panel, x_tick_values: list[float] | None = None, tick_spacing_years: float | None = None):
    """x-axis: ticks just below the x-axis line, labels from numeric tokens on the baseline below."""
    v, _ = tick_marks(g)
    ticks = sorted(set(round(t[0], 2) for t in v if p.y_axis_x - 3 <= t[0] <= p.x_max + 3 and t[2] <= p.x_axis_y + 0.5 and t[1] >= p.x_axis_y - 6))
    # label tokens: numbers whose baseline is 5-16 pt below the axis and x within the panel
    labels = []
    for tk in g.text:
        if p.x_axis_y - 16 <= tk.y <= p.x_axis_y - 4 and tk.size <= 9:
            for val, rx in numeric_runs(tk):
                if p.y_axis_x - 10 <= rx <= p.x_max + 10:
                    labels.append((val, rx))
    ax = None
    spacings = np.diff(ticks) if len(ticks) >= 2 else np.array([])
    regular = len(ticks) >= 5 and np.std(spacings) / np.mean(spacings) < 0.02
    if x_tick_values is not None and len(ticks) == len(x_tick_values):
        ax = _fit_axis(ticks, list(x_tick_values))
    elif regular and tick_spacing_years:
        # evenly spaced ticks starting at the y-axis: value = (x - first tick) / spacing * years_per_tick
        d = float(np.median(spacings))
        ax = Axis(origin=float(ticks[0]), scale=d / tick_spacing_years)
    else:
        # fall back to label matching: nearest tick to each label's centre
        pos, vals = [], []
        for val, rx in labels:
            centre = rx + 0.5 * 8 * len(f"{val:g}") / 2
            cand = [t for t in ticks if abs(t - centre) <= 8]
            if cand:
                pos.append(min(cand, key=lambda t: abs(t - centre))); vals.append(val)
        if x_tick_values is not None and len(pos) < 2 and len(ticks) >= len(x_tick_values):
            pos, vals = ticks[-len(x_tick_values):], list(x_tick_values)
        if len(pos) >= 2:
            ax = _fit_axis(pos, vals)
    if ax is None:
        raise ValueError(f"cannot calibrate x for panel {p.letter}: ticks={ticks} labels={labels}")
    # sanity: the y-axis should sit at (or very near) value 0
    if abs(float(ax.to_data(p.y_axis_x))) > 0.25:
        raise ValueError(f"panel {p.letter}: y-axis at t={float(ax.to_data(p.y_axis_x)):.2f}, expected ~0")
    ax.ticks = ticks
    ax.tick_values = [float(v) for v in ax.to_data(np.array(ticks))]
    p.xaxis = ax
    return ax


# ---------------------------------------------------------------- curves

def step_curves(g: PageGraphics, p: Panel, monotone: str = "down", min_vertices: int = 10) -> list[Polyline]:
    out = []
    for pl in g.polylines:
        if len(pl.points) < min_vertices or is_black(pl.stroke) or pl.dash:
            continue
        xs = np.array([q[0] for q in pl.points]); ys = np.array([q[1] for q in pl.points])
        if not (p.contains(xs.min(), ys.min()) and p.contains(xs.max(), ys.max())):
            continue
        dx = np.diff(xs); dy = np.diff(ys)
        aligned = np.all((np.abs(dx) < 0.05) | (np.abs(dy) < 0.05))
        mono = np.all(dy <= 1e-6) if monotone == "down" else np.all(dy >= -1e-6)
        if aligned and mono and np.all(dx >= -1e-6) and (xs.max() - xs.min()) >= 0.5 * (p.x_max - p.y_axis_x):
            out.append(pl)
    return out


def polyline_to_steps(pl_points, p: Panel):
    """Vertices -> (t, value fraction) right-continuous step: value after any jump at t."""
    pts = sorted(pl_points, key=lambda q: (round(q[0], 4), -q[1]))  # at equal x keep the lowest y last
    t = p.xaxis.to_data(np.array([q[0] for q in pts]))
    v = p.yaxis.to_data(np.array([q[1] for q in pts])) / 100.0
    # collapse equal t: keep last (lowest for KM)
    T, V = [], []
    for ti, vi in zip(t, v):
        if T and abs(ti - T[-1]) < 1e-9:
            V[-1] = vi
        else:
            T.append(float(ti)); V.append(float(vi))
    T = np.array(T); V = np.clip(np.array(V), 0, 1)
    if T[0] > 1e-6:
        T = np.concatenate([[0.0], T]); V = np.concatenate([[1.0], V])
    else:
        T[0] = 0.0; V[0] = 1.0 if V[0] > 0.98 else V[0]
    return T, V


def legend_names(g: PageGraphics, p: Panel, colours: list) -> dict[str, str]:
    """Match legend line samples (horizontal segments 6-16 pt in curve colours) to the label run to their right."""
    names = {}
    for col in colours:
        samples = [s for s in g.segments if colour_key(s.stroke) == col and s.is_horizontal and 6 <= s.length <= 16
                   and not (p.x_axis_y <= s.y0 <= p.y_max and p.y_axis_x <= s.x0 <= p.x_max)]
        best = None
        for s in samples:
            for tk in g.text:
                if abs(tk.y + 2.5 - s.y0) <= 4:
                    for i, (run, rx) in enumerate(tk.runs):
                        if 0 <= rx - max(s.x0, s.x1) <= 12:
                            # take this run and following runs until the next sample position
                            label = run
                            for run2, rx2 in tk.runs[i + 1:]:
                                if any(abs(rx2 - max(s2.x0, s2.x1)) <= 12 and s2 is not s for s2 in samples):
                                    break
                                label += " " + run2
                            best = label; break
                    if best: break
            if best: break
        names[col] = best
    return names


def censor_marks(g: PageGraphics, p: Panel, pl: Polyline) -> np.ndarray:
    """Vertical short segments in the curve colour whose midpoint lies on the step curve."""
    col = colour_key(pl.stroke)
    T, V = polyline_to_steps(pl.points, p)
    xs = np.array([q[0] for q in pl.points]); ys = np.array([q[1] for q in pl.points])
    marks = []
    for s in g.segments:
        if colour_key(s.stroke) == col and s.is_vertical and 1.5 <= s.length <= 4.5 and s.path_id != pl.path_id:
            xm = (s.x0 + s.x1) / 2; ym = (s.y0 + s.y1) / 2
            if not p.contains(xm, ym):
                continue
            # curve y at xm: value after jump (right-continuous) in page units
            idx = np.searchsorted(xs, xm, side="right") - 1
            idx = int(np.clip(idx, 0, len(xs) - 1))
            if abs(ys[idx] - ym) <= 1.5 or abs(ys[max(idx - 1, 0)] - ym) <= 1.5:
                marks.append(xm)
    return np.array(sorted(set(round(float(m), 3) for m in marks)))


def risk_table_from_text(g: PageGraphics, p: Panel, expected_times: list[float]):
    """Numbers on the 'at risk' baseline nearest below the panel, positioned by x -> time."""
    rows = [tk for tk in g.text if "at risk" in tk.text.lower()]
    rows = [tk for tk in rows if tk.y < p.x_axis_y]
    if not rows:
        return None
    row = max(rows, key=lambda tk: tk.y)   # nearest below
    nums = []
    for tk in g.text:
        if abs(tk.y - row.y) <= 1.0:
            for run, rx in tk.runs:
                for m in re.finditer(r"\d+", run):
                    x = rx + 0.5 * tk.size * m.start()
                    nums.append((int(m.group(0)), x))
    # drop huge glued numbers like '248338' by splitting into plausible parts if needed
    fixed = []
    for val, x in nums:
        s = str(val)
        if len(s) == 6 and int(s[:3]) < int(s[3:]) * 3:   # e.g. '248338' -> 248, 338
            fixed.append((int(s[:3]), x)); fixed.append((int(s[3:]), x + 0.5 * 8 * 3))
        else:
            fixed.append((val, x))
    fixed.sort(key=lambda v: v[1])
    if len(fixed) == len(expected_times):
        # positions of merged tokens are approximate; the left-to-right order is reliable
        return {float(t): int(v) for t, (v, _) in zip(sorted(expected_times), fixed)}
    times = p.xaxis.to_data(np.array([x for _, x in fixed]))
    table = {}
    for (val, x), t in zip(fixed, times):
        et = min(expected_times, key=lambda e: abs(e - t))
        if abs(et - t) <= 1.5 and et not in table:
            table[float(et)] = int(val)
    return table


# ---------------------------------------------------------------- polygons -> CIF (page 7)

def polygon_cif_curves(g: PageGraphics, p: Panel, min_vertices=20):
    """Each filled polygon anchored on the x-axis -> its upper boundary as a rising step (CIF fraction)."""
    curves = []
    for poly in g.polygons:
        if len(poly.points) < min_vertices or is_black(poly.fill):
            continue
        pts = np.array(poly.points)
        if not (p.contains(pts[:, 0].min(), pts[:, 1].min()) and p.contains(pts[:, 0].max(), pts[:, 1].max())):
            continue
        top = pts[pts[:, 1] > p.x_axis_y + 0.05]
        if len(top) < 3:
            continue
        order = np.argsort(top[:, 0], kind="stable")
        sx = top[order, 0]; sy = np.maximum.accumulate(top[order, 1])
        t = p.xaxis.to_data(sx); cif = p.yaxis.to_data(sy) / 100.0
        T, C = [0.0], [0.0]
        for ti, ci in zip(t, cif):
            if abs(ti - T[-1]) < 1e-9:
                C[-1] = max(C[-1], ci)
            else:
                T.append(float(ti)); C.append(float(max(ci, C[-1])))
        curves.append((poly, np.array(T), np.clip(np.array(C), 0, 1)))
    return curves


def legend_swatch_names(g: PageGraphics, p: Panel, fills: list) -> dict[str, str]:
    names = {}
    for col in fills:
        sw = [r for r in g.rects if colour_key(r.fill) == col and r.mode == "f" and 5 <= r.w <= 20 and 3 <= r.h <= 12
              and not (p.x_axis_y <= (r.y0 + r.y1) / 2 <= p.y_max and p.y_axis_x <= r.x0 <= p.x_max)]
        best = None
        for r in sw:
            for tk in g.text:
                if abs(tk.y + 2.5 - (r.y0 + r.y1) / 2) <= 5:
                    for i, (run, rx) in enumerate(tk.runs):
                        if 0 <= rx - max(r.x0, r.x1) <= 12:
                            label = run
                            for run2, rx2 in tk.runs[i + 1:]:
                                if any(abs(rx2 - max(r2.x0, r2.x1)) <= 12 and r2 is not r for r2 in g.rects if r2.mode == "f" and 5 <= r2.w <= 20 and 3 <= r2.h <= 12):
                                    break
                                label += " " + run2
                            best = label; break
                    if best: break
            if best: break
        names[col] = best
    return names


# ---------------------------------------------------------------- Kermen drivers

# Expected curve pairs per panel, ordered by end value (lower first). Legend text is used as a cross-check only.
KERMEN_FIG2_NAMES = {
    "A": [("kermen_fig2a_overall_survival", "Overall survival"),
          ("kermen_fig2a_valve_related_survival", "Valve-related survival")],
    "B": [("kermen_fig2b_freedom_stage23", "Freedom from moderate/severe SVD (VARC-3 stage 2/3)"),
          ("kermen_fig2b_freedom_stage3", "Freedom from severe SVD (VARC-3 stage 3)")],
}
# Fig 4A fill colours (CMYK) -> curve id. The three lighter fills on the page are 68% confidence bands.
KERMEN_FIG4A_FILLS = {
    "('k', (1.0, 0.57, 0.0, 0.4))": ("kermen_fig4a_nonvalve_death", "Non valve-related death"),
    "('k', (0.03, 1.0, 0.7, 0.12))": ("kermen_fig4a_valve_related_death", "Valve-related death"),
    "('k', (0.57, 0.06, 0.92, 0.19))": ("kermen_fig4a_explant_svd", "Explantation due to SVD"),
}


def extract_kermen_fig2(pdf_path, page_index=3) -> dict[str, CurveData]:
    g = parse_page(pdf_path, page_index)
    panels = attach_frames_and_letters(g, find_panels_from_labels(g))
    out = {}
    for p in panels:
        calibrate_x(g, p, tick_spacing_years=1.0)
        curves = step_curves(g, p, monotone="down")
        if len(curves) != 2:
            raise ValueError(f"panel {p.letter}: expected 2 step curves, found {len(curves)}")
        cols = sorted(set(colour_key(c.stroke) for c in curves))
        legend = legend_names(g, p, cols)
        rt = risk_table_from_text(g, p, [0, 5, 9.5, 10])
        parsed = []
        for pl in curves:
            T, V = polyline_to_steps(pl.points, p)
            parsed.append((float(V[-1]), pl, T, V))
        parsed.sort(key=lambda z: z[0])              # lower end value first
        for (endv, pl, T, V), (cid, name) in zip(parsed, KERMEN_FIG2_NAMES[p.letter]):
            cm = censor_marks(g, p, pl)
            cd = CurveData(id=cid, figure="Fig 2", panel=p.letter, name=name, arm="all", kind="KM", t=T, y=V,
                           censor_times=p.xaxis.to_data(cm) if len(cm) else None,
                           t_end=float(p.xaxis.to_data(max(q[0] for q in pl.points))),
                           provenance={"page": page_index, "stroke": colour_key(pl.stroke), "n_vertices": len(pl.points),
                                       "x_axis": (p.xaxis.origin, p.xaxis.scale), "y_axis": (p.yaxis.origin, p.yaxis.scale),
                                       "n_censor_marks": int(len(cm)), "legend_text_for_colour": legend.get(colour_key(pl.stroke)),
                                       "risk_table_text": rt, "naming_rule": "end-value order"})
            out[cid] = cd
    return out


def extract_kermen_fig4a(pdf_path, page_index=6) -> dict[str, CurveData]:
    g = parse_page(pdf_path, page_index)
    panels = attach_frames_and_letters(g, find_panels_from_labels(g))
    pA = max(panels, key=lambda p: p.x_axis_y)        # top panel
    pA.letter = "A"
    calibrate_x(g, pA, x_tick_values=[5.0, 10.0])
    curves = polygon_cif_curves(g, pA)
    fills = sorted(set(colour_key(c[0].fill) for c in curves))
    legend = legend_swatch_names(g, pA, fills)
    out = {}
    for poly, T, C in curves:
        key = colour_key(poly.fill)
        if key not in KERMEN_FIG4A_FILLS:
            continue                                  # confidence-band polygons
        cid, name = KERMEN_FIG4A_FILLS[key]
        cd = CurveData(id=cid, figure="Fig 4", panel="A", name=name, arm="all", kind="CIF", t=T, y=1.0 - C,
                       censor_times=None, t_end=float(T[-1]),
                       provenance={"page": page_index, "fill": key, "n_vertices": len(poly.points),
                                   "x_axis": (pA.xaxis.origin, pA.xaxis.scale), "y_axis": (pA.yaxis.origin, pA.yaxis.scale),
                                   "legend_text_for_colour": legend.get(key), "n_polygons_on_panel": len(curves)})
        out[cid] = cd
    if len(out) != 3:
        raise ValueError(f"Fig 4A: expected 3 CIF polygons, matched {len(out)} of {len(curves)}")
    return out


def value_at(cd: CurveData, t: float, side="right") -> float:
    from survival import step_eval
    return float(step_eval(t, cd.t, cd.y, side=side)[0])
