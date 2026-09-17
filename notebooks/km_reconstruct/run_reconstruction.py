"""Driver: extract curves -> Guyot reconstruction -> validation -> outputs.

Usage:
  python run_reconstruction.py --figures kermen_fig2 kermen_fig4a [--strict] [--censor-placement coordinate|placed]

Outputs in output/:
  ipd_<curve_id>.csv          reconstructed individual patient data (time, event)
  targets.json                digitised curves on a 0.1-y grid + anchors + risk tables (simulator targets)
  onset_priors.json           Weibull / log-logistic fits per curve (initialisers only; see plan §7)
  validation_report.md        PASS/FAIL table
  qa_<figure>.png             digitised step, reconstructed KM, anchors
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import config as C  # noqa: E402
from guyot import exact_inversion, reconstruct  # noqa: E402
from survival import fit_loglogistic, fit_weibull, km, step_eval  # noqa: E402
from validate import CheckResult, check_anchors, check_bookkeeping, check_reconstruction, n_at_risk_curve, report_markdown  # noqa: E402
from vector_extract import extract_kermen_fig2, extract_kermen_fig4a  # noqa: E402

FIGURES = {
    "kermen_fig2": lambda: extract_kermen_fig2(C.PDF_KERMEN),
    "kermen_fig4a": lambda: extract_kermen_fig4a(C.PDF_KERMEN),
}


def load_raster_figures():
    try:
        from raster_trace import extract_notion_fig3, extract_wakami_fig1
        FIGURES["notion_fig3"] = lambda: extract_notion_fig3(C.PDF_NOTION)
        FIGURES["wakami_fig1"] = lambda: extract_wakami_fig1(C.PDF_WAKAMI)
    except ImportError:
        pass


def process_curve(cd, spec, censor_placement, results, targets, priors, out_dir, last_interval="observed"):
    is_raster = spec.source.startswith("raster")
    extra_tol = 0.2 if is_raster else 0.0
    results += check_anchors(cd.id, cd.t, cd.y, cd.kind, spec.anchors, cd.t_end, extra_tol=extra_tol, label="digitised")
    if spec.risk is not None and not spec.use_risk_table:
        spec = dataclasses.replace(spec, risk=None)
        results.append(CheckResult(cd.id, "risk table", None, None, None, True, "printed row NOT used (see config notes); Guyot no-numbers-at-risk case"))
    if spec.tot_events is not None and not spec.use_tot_events:
        results.append(CheckResult(cd.id, "published total events (informational; inconsistent with curve, see config notes)", spec.tot_events, None, None, True))
        spec = dataclasses.replace(spec, tot_events=None)
    risk_times = spec.risk.times if spec.risk else [0]
    risk_n = spec.risk.n if spec.risk else [spec.n0]
    obs = cd.censor_times if (last_interval == "observed" and cd.censor_times is not None and len(cd.censor_times)) else None
    if spec.censoring_from:
        src = pd.read_csv(out_dir / f"ipd_{spec.censoring_from}.csv")
        obs = src.loc[src.event == 0, "time"].to_numpy(float)
        results.append(CheckResult(cd.id, "censoring", None, None, None, True, f"{len(obs)} censoring times reused from {spec.censoring_from} (same panel)"))
    times, events, diag = reconstruct(cd.t, cd.y, risk_times, risk_n, tot_events=spec.tot_events, censor_at=censor_placement, obs_censor_times=obs)
    kt, kS = km(times, events)
    results += check_reconstruction(cd.id, cd.t, cd.y, kt, kS, cd.kind, cd.t_end, diag, tol_pp=1.5 if is_raster else 1.0)
    Tn, nn = n_at_risk_curve(diag)
    d_arr = np.array(diag["d"]); n_end = int(nn[d_arr > 0][-1]) if (d_arr > 0).any() else None  # at risk at the last event
    results += check_anchors(cd.id, kt, kS, cd.kind, spec.anchors, cd.t_end, extra_tol=extra_tol, label="reconstructed (Guyot)", n_end=n_end)
    results += check_bookkeeping(cd.id, diag, spec)
    estimator_note = "AJ-as-KM" if cd.kind == "CIF" else "KM"
    df = pd.DataFrame({"time": times, "event": events})
    df["arm"] = spec.arm; df["source"] = spec.id; df["estimator"] = estimator_note; df["stratum"] = spec.stratum
    df.to_csv(out_dir / f"ipd_{cd.id}.csv", index=False)
    # vector-only exact step inversion (second reconstruction, reported alongside)
    if spec.source == "vector":
        xt, xe, xdiag = exact_inversion(cd.t, cd.y, spec.n0, t_end=cd.t_end)
        xkt, xkS = km(xt, xe)
        grid = np.arange(0.0, cd.t_end + 1e-9, 0.1)
        dd = np.max(np.abs(100 * (step_eval(grid, cd.t, cd.y) - step_eval(grid, xkt, xkS))))
        # diagnostic only: single-event drops of < 0.3 pp are at the coordinate precision limit, so the
        # multiplicity choice cascades; the implied risk set at the FIRST drop after 5 y is still informative
        first5 = [n for t, n in xdiag["n_at_risk_trace"] if t >= 5.0][:1]
        results.append(CheckResult(cd.id, "exact step inversion (diagnostic): total events; implied n at risk at first step after 5 y", spec.tot_events,
                                   xdiag["total_events"], None, True, f"implied n at risk after 5 y = {first5[0] if first5 else 'n/a'}; max |diff| vs digitised {dd:.1f} pp"))
        dfx = pd.DataFrame({"time": xt, "event": xe}); dfx["arm"] = spec.arm; dfx["source"] = spec.id; dfx["estimator"] = estimator_note + "_exact"; dfx["stratum"] = spec.stratum
        dfx.to_csv(out_dir / f"ipd_exact_{cd.id}.csv", index=False)
        diag["exact_inversion"] = {k: v for k, v in xdiag.items() if k != "n_at_risk_trace"} | {"n_at_risk_trace": xdiag["n_at_risk_trace"]}
    grid = np.round(np.arange(0.0, cd.t_end + 1e-9, 0.1), 1)
    targets[cd.id] = {
        "figure": cd.figure, "panel": cd.panel, "name": cd.name, "arm": spec.arm, "kind": cd.kind,
        "estimator": spec.estimator, "stratum": spec.stratum, "n0": spec.n0, "t_end": cd.t_end,
        "risk_table": {"times": risk_times, "n": risk_n}, "tot_events": spec.tot_events,
        "anchors": [{"t": a.t, "value": a.value, "se": a.se, "note": a.note} for a in spec.anchors],
        "grid_t": grid.tolist(),
        "grid_value_pct": [float(v) for v in (100 * (1 - step_eval(grid, cd.t, cd.y)) if cd.kind == "CIF" else 100 * step_eval(grid, cd.t, cd.y))],
        "steps_t": cd.t.tolist(), "steps_y": cd.y.tolist(),
        "censor_times": None if cd.censor_times is None else [float(v) for v in cd.censor_times],
        "guyot_case": diag["case"], "guyot_total_events": diag["total_events"], "last_interval_rule": diag.get("last_interval_rule"),
        "exact_inversion": diag.get("exact_inversion"),
        "provenance": {k: (v if isinstance(v, (int, float, str, list, dict, type(None))) else str(v)) for k, v in cd.provenance.items()},
    }
    if events.sum() >= 3:
        w = fit_weibull(times, events); ll = fit_loglogistic(times, events)
        priors[cd.id] = {"stratum": spec.stratum, "estimator": estimator_note, "weibull": w, "loglogistic": ll,
                         "note": "fit to OBSERVED-label incidence; initialiser for calibration only, not the latent onset distribution"}
    return times, events, kt, kS


def qa_plot(figure_name, items, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(items)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 3.6), squeeze=False)
    for ax, (cd, spec, kt, kS) in zip(axes[0], items):
        sign = (lambda v: 100 * (1 - v)) if cd.kind == "CIF" else (lambda v: 100 * v)
        ax.step(cd.t, sign(cd.y), where="post", lw=1.6, label="digitised")
        ax.step(kt, sign(kS), where="post", lw=1.0, ls="--", label="Guyot reconstruction")
        for a in spec.anchors:
            if a.t is not None:
                ax.errorbar([a.t], [a.value], yerr=[a.se or 0], fmt="o", ms=4, capsize=3, label=None)
        if cd.censor_times is not None:
            ax.plot(cd.censor_times, sign(step_eval(cd.censor_times, cd.t, cd.y)), "|", ms=6, alpha=0.6)
        ax.set_title(cd.name[:40], fontsize=8); ax.set_xlabel("years"); ax.set_ylabel("%" if cd.kind == "KM" else "CIF %")
        ax.grid(alpha=0.3); ax.legend(fontsize=6)
    fig.tight_layout(); fig.savefig(out_dir / f"qa_{figure_name}.png", dpi=130); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--figures", nargs="+", default=["kermen_fig2", "kermen_fig4a"])
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--censor-placement", default="coordinate", choices=["coordinate", "placed"])
    ap.add_argument("--last-interval", default="observed", choices=["observed", "guyot"], help="censoring in the last interval: digitised + marks (default) or Guyot span-ratio guess")
    args = ap.parse_args(argv)
    load_raster_figures()
    out_dir = C.OUT_DIR; out_dir.mkdir(exist_ok=True)
    results, targets, priors = [], {}, {}
    if (out_dir / "targets.json").exists():
        targets = json.loads((out_dir / "targets.json").read_text())
    if (out_dir / "onset_priors.json").exists():
        priors = json.loads((out_dir / "onset_priors.json").read_text())
    for fig_name in args.figures:
        curves = FIGURES[fig_name]()
        items = []
        for cid, cd in curves.items():
            spec = C.CURVES_BY_ID.get(cid)
            if spec is None:
                results.append(CheckResult(cid, "config entry", None, None, None, False, "no CurveSpec for this curve id"))
                continue
            times, events, kt, kS = process_curve(cd, spec, args.censor_placement, results, targets, priors, out_dir, last_interval=args.last_interval)
            items.append((cd, spec, kt, kS))
        qa_plot(fig_name, items, out_dir)
    (out_dir / "targets.json").write_text(json.dumps(targets, indent=1))
    (out_dir / "onset_priors.json").write_text(json.dumps(priors, indent=1))
    md = report_markdown([r for r in results if r is not None])
    (out_dir / "validation_report.md").write_text(md)
    print(md)
    n_fail = sum(not r.passed for r in results if r is not None)
    if args.strict and n_fail:
        raise SystemExit(f"{n_fail} validation checks failed")


if __name__ == "__main__":
    main()
