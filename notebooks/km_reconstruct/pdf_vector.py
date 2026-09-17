"""Parse a PDF page's content stream into geometry and text in page coordinates.

Tracks the current transformation matrix (CTM) across `cm`, `q`, `Q`, so that pages
drawn as thousands of one-segment paths each under their own translation (as in the
JTCVS Open figures) come out in a single coordinate system.

Output: PageGraphics with
  segments   -- individual stroked line segments (x0, y0, x1, y1) with stroke state
  polylines  -- stroked subpaths (list of points) with stroke state
  rects      -- `re` rectangles with paint mode ('S', 'f', 'B', 'n')
  polygons   -- filled subpaths (>= 3 points)
  text       -- decoded text tokens with position and effective font size
"""
from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

from pypdf import PdfReader
from pypdf.generic import ArrayObject, ContentStream

logging.getLogger("pypdf").setLevel(logging.ERROR)

Matrix = tuple[float, float, float, float, float, float]
IDENTITY: Matrix = (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)


def mmul(m1: Matrix, m2: Matrix) -> Matrix:
    """PDF matrix product m1 x m2 (apply m1 first, then m2)."""
    a1, b1, c1, d1, e1, f1 = m1
    a2, b2, c2, d2, e2, f2 = m2
    return (
        a1 * a2 + b1 * c2, a1 * b2 + b1 * d2,
        c1 * a2 + d1 * c2, c1 * b2 + d1 * d2,
        e1 * a2 + f1 * c2 + e2, e1 * b2 + f1 * d2 + f2,
    )


def apply(m: Matrix, x: float, y: float) -> tuple[float, float]:
    a, b, c, d, e, f = m
    return (a * x + c * y + e, b * x + d * y + f)


@dataclass
class Segment:
    x0: float; y0: float; x1: float; y1: float
    stroke: tuple; dash: tuple; width: float; path_id: int

    @property
    def length(self) -> float:
        return math.hypot(self.x1 - self.x0, self.y1 - self.y0)

    @property
    def is_vertical(self) -> bool:
        return abs(self.x1 - self.x0) < 0.05

    @property
    def is_horizontal(self) -> bool:
        return abs(self.y1 - self.y0) < 0.05


@dataclass
class Polyline:
    points: list[tuple[float, float]]
    stroke: tuple; dash: tuple; width: float; closed: bool; path_id: int


@dataclass
class Rect:
    x0: float; y0: float; x1: float; y1: float
    fill: tuple; stroke: tuple; mode: str; path_id: int

    @property
    def w(self) -> float:
        return abs(self.x1 - self.x0)

    @property
    def h(self) -> float:
        return abs(self.y1 - self.y0)


@dataclass
class Polygon:
    points: list[tuple[float, float]]
    fill: tuple; path_id: int


@dataclass
class TextToken:
    text: str
    x: float
    y: float
    size: float
    runs: list[tuple[str, float]] = field(default_factory=list)   # (sub-text, approx x)


@dataclass
class PageGraphics:
    segments: list[Segment] = field(default_factory=list)
    polylines: list[Polyline] = field(default_factory=list)
    rects: list[Rect] = field(default_factory=list)
    polygons: list[Polygon] = field(default_factory=list)
    text: list[TextToken] = field(default_factory=list)
    mediabox: tuple[float, float, float, float] = (0, 0, 0, 0)


def _num(v) -> float:
    return float(v)


def _bezier(p0, p1, p2, p3, n=8):
    pts = []
    for i in range(1, n + 1):
        t = i / n
        u = 1 - t
        x = u**3 * p0[0] + 3 * u**2 * t * p1[0] + 3 * u * t**2 * p2[0] + t**3 * p3[0]
        y = u**3 * p0[1] + 3 * u**2 * t * p1[1] + 3 * u * t**2 * p2[1] + t**3 * p3[1]
        pts.append((x, y))
    return pts


def parse_page(pdf_path, page_index: int) -> PageGraphics:
    reader = PdfReader(str(pdf_path))
    page = reader.pages[page_index]
    contents = page.get_contents()
    stream = ContentStream(contents, reader) if not isinstance(contents, ContentStream) else contents

    g = PageGraphics()
    mb = page.mediabox
    g.mediabox = (float(mb.left), float(mb.bottom), float(mb.right), float(mb.top))

    ctm: Matrix = IDENTITY
    stroke: tuple = ("G", (0.0,))
    fill: tuple = ("g", (0.0,))
    dash: tuple = ()
    width: float = 1.0
    stack: list = []
    subpaths: list[list[tuple[float, float]]] = []
    cur: list[tuple[float, float]] = []
    cur_closed = False
    rect_buf: list[tuple] = []
    path_id = 0
    cs_stroke_name = None
    cs_fill_name = None

    def flush_stroke():
        nonlocal path_id
        for sp, closed in [(s, False) for s in subpaths] + ([(cur, cur_closed)] if cur else []):
            if len(sp) >= 2:
                g.polylines.append(Polyline(list(sp), stroke, dash, width, closed, path_id))
                for i in range(len(sp) - 1):
                    g.segments.append(Segment(sp[i][0], sp[i][1], sp[i + 1][0], sp[i + 1][1], stroke, dash, width, path_id))
                if closed and len(sp) >= 3:
                    g.segments.append(Segment(sp[-1][0], sp[-1][1], sp[0][0], sp[0][1], stroke, dash, width, path_id))
        for (x0, y0, x1, y1) in rect_buf:
            g.rects.append(Rect(x0, y0, x1, y1, fill, stroke, "S", path_id))

    def flush_fill():
        for sp in subpaths + ([cur] if cur else []):
            if len(sp) >= 3:
                g.polygons.append(Polygon(list(sp), fill, path_id))
        for (x0, y0, x1, y1) in rect_buf:
            g.rects.append(Rect(x0, y0, x1, y1, fill, stroke, "f", path_id))

    def end_path():
        nonlocal subpaths, cur, cur_closed, rect_buf, path_id
        subpaths = []
        cur = []
        cur_closed = False
        rect_buf = []
        path_id += 1

    for operands, op in stream.operations:
        op = op.decode() if isinstance(op, bytes) else str(op)
        # --- graphics state
        if op == "q":
            stack.append((ctm, stroke, fill, dash, width, cs_stroke_name, cs_fill_name))
        elif op == "Q":
            if stack:
                ctm, stroke, fill, dash, width, cs_stroke_name, cs_fill_name = stack.pop()
        elif op == "cm":
            m = tuple(_num(v) for v in operands[:6])
            ctm = mmul(m, ctm)
        elif op == "w":
            width = _num(operands[0])
        elif op == "d":
            arr = operands[0]
            dash = tuple(_num(v) for v in arr) if isinstance(arr, (list, ArrayObject)) else ()
        # --- colour
        elif op == "G":
            stroke = ("G", tuple(_num(v) for v in operands))
        elif op == "g":
            fill = ("g", tuple(_num(v) for v in operands))
        elif op == "RG":
            stroke = ("RG", tuple(round(_num(v), 3) for v in operands))
        elif op == "rg":
            fill = ("rg", tuple(round(_num(v), 3) for v in operands))
        elif op == "K":
            stroke = ("K", tuple(round(_num(v), 3) for v in operands))
        elif op == "k":
            fill = ("k", tuple(round(_num(v), 3) for v in operands))
        elif op == "CS":
            cs_stroke_name = str(operands[0])
        elif op == "cs":
            cs_fill_name = str(operands[0])
        elif op in ("SC", "SCN"):
            comps = tuple(round(_num(v), 3) for v in operands if not isinstance(v, str) and not hasattr(v, "startswith"))
            stroke = ("SCN", cs_stroke_name, comps)
        elif op in ("sc", "scn"):
            comps = tuple(round(_num(v), 3) for v in operands if not isinstance(v, str) and not hasattr(v, "startswith"))
            fill = ("scn", cs_fill_name, comps)
        # --- path construction
        elif op == "m":
            if cur:
                subpaths.append(cur)
            cur = [apply(ctm, _num(operands[0]), _num(operands[1]))]
            cur_closed = False
        elif op == "l":
            if cur is None:
                cur = []
            cur.append(apply(ctm, _num(operands[0]), _num(operands[1])))
        elif op in ("c", "v", "y"):
            if cur:
                p0 = cur[-1]
                vals = [_num(v) for v in operands]
                if op == "c":
                    p1 = apply(ctm, vals[0], vals[1]); p2 = apply(ctm, vals[2], vals[3]); p3 = apply(ctm, vals[4], vals[5])
                elif op == "v":
                    p1 = p0; p2 = apply(ctm, vals[0], vals[1]); p3 = apply(ctm, vals[2], vals[3])
                else:
                    p1 = apply(ctm, vals[0], vals[1]); p3 = apply(ctm, vals[2], vals[3]); p2 = p3
                cur.extend(_bezier(p0, p1, p2, p3))
        elif op == "h":
            cur_closed = True
        elif op == "re":
            x, y, w, h = (_num(v) for v in operands[:4])
            p = apply(ctm, x, y); q = apply(ctm, x + w, y + h)
            rect_buf.append((p[0], p[1], q[0], q[1]))
        # --- path painting
        elif op in ("S", "s"):
            if op == "s":
                cur_closed = True
            flush_stroke(); end_path()
        elif op in ("f", "F", "f*"):
            flush_fill(); end_path()
        elif op in ("B", "B*", "b", "b*"):
            flush_fill(); flush_stroke(); end_path()
        elif op == "n":
            for (x0, y0, x1, y1) in rect_buf:
                g.rects.append(Rect(x0, y0, x1, y1, fill, stroke, "n", path_id))
            end_path()
        # text is handled by pypdf's extractor below

    # --- text with positions (pypdf handles font decoding)
    tokens: list[TextToken] = []

    def visitor(text, cm, tm, font_dict, font_size):
        if text is None or not text.strip():
            return
        m = mmul(tuple(float(v) for v in tm), tuple(float(v) for v in cm))
        x, y = apply(m, 0.0, 0.0)
        scale = math.hypot(m[0], m[1]) or 1.0
        size = float(font_size or 1.0) * scale
        clean = text.replace("\n", " ").strip()
        # approximate per-run x positions for whitespace-separated sub-tokens
        runs = []
        cx = x
        for part in clean.split(" "):
            if part:
                runs.append((part, cx))
            cx += (len(part) + 1) * 0.5 * size
        tokens.append(TextToken(clean, x, y, size, runs))

    page.extract_text(visitor_text=visitor)
    g.text = tokens
    return g


def summarize(g: PageGraphics) -> str:
    return (f"segments={len(g.segments)} polylines={len(g.polylines)} rects={len(g.rects)} "
            f"polygons={len(g.polygons)} text={len(g.text)} mediabox={g.mediabox}")
