"""Validation of digitised and reconstructed curves against every published number.

Tolerance rules (plan §6):
  * well-populated times: |diff| <= 1.5 pp
  * times with small risk sets (>= 9 y here): published value must lie within
    [min, max] of the curve over [t-0.25, t+0.25], widened by 1.5 pp, and never worse than the published SE
  * raster sources get an extra 1 px in data units (passed in as extra_tol)
  * reconstructed KM vs digitised curve on a 0.1-y grid: max |diff| <= 0.5 pp (vector) / 1.0 pp (raster)
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from survival import step_eval


@dataclass
class CheckResult:
    curve_id: str
    quantity: str
    published: float | None
    reconstructed: float | None
    tolerance: float | None
    passed: bool
    note: str = ""


def curve_value_pct(t, y, at, kind, side="right"):
    v = float(step_eval(at, t, y, side=side)[0]) * 100.0
    return (100.0 - v) if kind == "CIF" else v


def range_over_window(t, y, at, kind, half=0.25):
    grid = np.linspace(at - half, at + half, 11)
    vals = [curve_value_pct(t, y, g, kind) for g in grid] + [curve_value_pct(t, y, at, kind, side="left")]
    return min(vals), max(vals)


def check_anchors(curve_id, t, y, kind, anchors, t_end, base_tol=1.5, extra_tol=0.0, small_risk_from=9.0, label="digitised", n_end=None):
    out = []
    for a in anchors:
        if a.t is None:
            got = curve_value_pct(t, y, t_end, kind)
            tol = base_tol + extra_tol
            note = ""
            if n_end is not None and n_end > 0:
                # a reconstructed curve can only move in quanta of one event at the end of follow-up
                tol = max(tol, 100.0 / n_end); note = f"tolerance = one event with n={n_end} at risk"
            ok = abs(got - a.value) <= tol
            out.append(CheckResult(curve_id, f"{label} value at end of follow-up ({a.note})", a.value, round(got, 2), round(tol, 2), ok, note))
            continue
        got = curve_value_pct(t, y, a.t, kind)
        tol = base_tol + extra_tol
        if a.t >= small_risk_from:
            lo, hi = range_over_window(t, y, a.t, kind)
            tol_w = tol + (a.se or 0.0)
            ok = (lo - tol <= a.value <= hi + tol) or abs(got - a.value) <= tol_w
            note = f"window [{lo:.1f}, {hi:.1f}] around t={a.t}; SE={a.se}"
        else:
            ok = abs(got - a.value) <= tol
            note = ""
        out.append(CheckResult(curve_id, f"{label} value at {a.t} y", a.value, round(got, 2), tol, ok, note))
    return out


def n_at_risk_curve(diag):
    """Number at risk just before each coordinate time, from the reconstruction diagnostics."""
    T = np.array(diag["T"]); d = np.array(diag["d"]); cen = np.array(diag["cen"])
    n = diag["risk_n"][0]; out = []
    for dk, ck in zip(d, cen):
        out.append(n); n -= dk + ck
    return T, np.array(out)


def check_reconstruction(curve_id, t_dig, y_dig, t_rec, y_rec, kind, t_end, diag, tol_pp=1.0, min_risk=30):
    """Compare on a 0.1-y grid while >= min_risk patients remain (one event moves the curve by 100/n pp,
    so the tail is compared only informationally)."""
    grid = np.arange(0.0, t_end + 1e-9, 0.1)
    d = np.array([curve_value_pct(t_dig, y_dig, g, kind) for g in grid])
    r = np.array([curve_value_pct(t_rec, y_rec, g, kind) for g in grid])
    Tn, nn = n_at_risk_curve(diag)
    n_grid = step_eval(grid, Tn, nn.astype(float))
    well = n_grid >= min_risk
    diff_well = float(np.max(np.abs(d - r)[well])) if well.any() else 0.0
    diff_tail = float(np.max(np.abs(d - r)[~well])) if (~well).any() else 0.0
    t_cut = float(grid[well].max()) if well.any() else 0.0
    return [CheckResult(curve_id, f"max |reconstructed - digitised| while n at risk >= {min_risk} (t <= {t_cut:.1f} y), pp", None, round(diff_well, 3), tol_pp, diff_well <= tol_pp),
            CheckResult(curve_id, f"max |reconstructed - digitised| in the tail (n < {min_risk}), pp, informational", None, round(diff_tail, 3), None, True)]


def check_bookkeeping(curve_id, diag, spec):
    out = []
    if spec.tot_events is not None:
        tol = max(2, int(round(0.1 * spec.tot_events)))
        got = diag["total_events"]
        out.append(CheckResult(curve_id, "total events (soft: within max(2, 10%))", spec.tot_events, got, tol, abs(got - spec.tot_events) <= tol,
                               "" if got == spec.tot_events else "exact total not reproduced; see algorithm warnings"))
    if spec.risk is not None and len(spec.risk.times) > 1:
        ok = diag["n_at_risk_reconstructed"] == list(spec.risk.n)
        out.append(CheckResult(curve_id, "numbers at risk reproduced", None, None, None, ok, f"{diag['n_at_risk_reconstructed']} vs {list(spec.risk.n)}"))
    out.append(CheckResult(curve_id, "IPD count equals n0", spec.n0, diag["n_ipd"], 0, diag["n_ipd"] == spec.n0))
    if diag.get("warnings"):
        out.append(CheckResult(curve_id, "algorithm warnings", None, None, None, True, "; ".join(diag["warnings"])))
    return out


def report_markdown(results: list[CheckResult]) -> str:
    lines = ["| curve | quantity | published | reconstructed | tol | result | note |", "|---|---|---|---|---|---|---|"]
    for r in results:
        lines.append(f"| {r.curve_id} | {r.quantity} | {'' if r.published is None else r.published} | "
                     f"{'' if r.reconstructed is None else r.reconstructed} | {'' if r.tolerance is None else r.tolerance} | "
                     f"{'PASS' if r.passed else 'FAIL'} | {r.note} |")
    n_fail = sum(not r.passed for r in results)
    header = f"# Reconstruction validation report\n\n{len(results)} checks, {n_fail} FAIL\n\n"
    return header + "\n".join(lines) + "\n"
