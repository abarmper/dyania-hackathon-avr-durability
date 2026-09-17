"""Sequential calibration of the simulator (plan: physics -> noise -> death -> onset per class -> checks).

Every step records target, simulated value, tolerance and PASS/FAIL in output/calibration_report.md and
the calibrated parameters in output/params_calibrated.json. Curve targets come from
notebooks/km_reconstruct/output/targets.json (reconstructed from the papers' figures, validated there).

Run:  /data/abar/alexenv/bin/python calibrate.py [--n-sim 20000] [--quick]
"""
from __future__ import annotations

import argparse
import copy
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
KM = HERE.parent / "km_reconstruct"
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(KM))

from calibrate_onset import censoring_sampler, observed_curve, risk_weights, target_curve  # noqa: E402
from params import default_params  # noqa: E402
from physics import fit_physics_constants  # noqa: E402
from simulate import _core, simulate_matched  # noqa: E402
from valve_tables import COHORT_PROFILES  # noqa: E402

OUT = HERE / "output"
TARGETS = json.loads((KM / "output" / "targets.json").read_text())
N_WORKERS = max(1, min(os.cpu_count() or 1, 32))


def _pool():
    return ProcessPoolExecutor(max_workers=N_WORKERS, mp_context=mp.get_context("fork"))


def _curve_worker(args):
    """One simulated observed curve for a (paper, class, theta) point; runs in a worker process."""
    p, paper, cls, label, tid, theta, n, seed, grid = args
    target = TARGETS[tid]; censor = censoring_sampler(tid)
    rng = np.random.default_rng(seed)
    tt, st = simulate_matched(paper, {cls: theta}, n, target["t_end"], rng, p, label=label, censor=censor)
    return observed_curve(tt, st, np.asarray(grid), target["estimator"])


def _perigon_worker(args):
    p, n, seed = args
    return perigon_stats(p, n, seed)


def _check_worker(args):
    p, tid, paper, label, n, seed = args
    target = TARGETS[tid]
    return tid, matched_curve(p, paper, label, target, n, seed, censor=censoring_sampler(tid))

# curve targets: target id -> (paper profile, class calibrated, label rule, estimator prefix)
CURVE_STEPS = [
    ("kermen_fig2b_freedom_stage23", "kermen", "perimount", "varc3_single2", 3.0),
    ("wakami_fig1_svd", "wakami", "trifecta", "varc3_single3_or_reint", 2.5),   # n = 110: one event = 0.9 pp
    ("notion_fig3_modsvd_tavi", "notion_tavi", "self_expanding", "varc3_haemo_single2", 3.0),
]
CHECKS = [
    ("kermen_fig2b_freedom_stage3", "kermen", "varc3_single3", 5.0),
    ("notion_fig3_sevsvd_savr", "notion_savr", "varc3_haemo_single3", 3.0),
    ("notion_fig3_sevsvd_tavi", "notion_tavi", "varc3_haemo_single3", 3.0),
    ("kermen_fig4a_explant_svd", "kermen", "reintervention", 3.0),
]


class Report:
    def __init__(self):
        self.rows = []

    def add(self, step, quantity, target, simulated, tol, ok, note=""):
        if step.startswith("noise (pass 1"):      # onset not yet calibrated: informational only
            note = ("informational; " + note).strip("; "); ok = True
        self.rows.append((step, quantity, target, simulated, tol, "PASS" if ok else "FAIL", note))
        print(f"[{step}] {quantity}: target {target} | simulated {simulated} | tol {tol} | {'PASS' if ok else 'FAIL'} {note}")

    def markdown(self):
        n_fail = sum(r[5] == "FAIL" for r in self.rows)
        lines = [f"# Simulator calibration report\n\n{len(self.rows)} checks, {n_fail} FAIL\n",
                 "| step | quantity | target | simulated | tolerance | result | note |", "|---|---|---|---|---|---|---|"]
        for r in self.rows:
            lines.append("| " + " | ".join("" if v is None else str(v) for v in r) + " |")
        return "\n".join(lines) + "\n"


def matched_curve(p, paper, label, target, n, seed, theta=None, censor=None):
    rng = np.random.default_rng(seed)
    grid = np.arange(1.0, target["t_end"] + 1e-9, 1.0)
    if label == "reintervention":
        df, t, vt, obs, ref, stages, times, t_end, ends, latent = _core(p, n, rng, profile=COHORT_PROFILES[paper], theta_by_class=theta,
                                                                         schedule=COHORT_PROFILES[paper].schedule_y,
                                                                         reference_time=COHORT_PROFILES[paper].reference_time_y, t_max=target["t_end"])
        cens = np.minimum(censor(n, rng), target["t_end"]) if censor else np.full(n, target["t_end"])
        tl = np.where(np.isfinite(ends["t_re"]), ends["t_re"], np.inf)
        tt = np.minimum.reduce([tl, ends["t_death"], cens]); st = np.where(tt == tl, 1, np.where(tt == ends["t_death"], 2, 0))
    else:
        tt, st = simulate_matched(paper, theta, n, target["t_end"], rng, p, label=label, censor=censor)
    sim = observed_curve(tt, st, grid, target["estimator"])
    tv = target_curve(target, grid)
    return grid, sim, tv


def risk_ok(target, grid, min_n=30, t_max=10.0):
    """Grid points that count: >= min_n at risk in the paper and <= 10 y (beyond that every risk set is < 30)."""
    rt = np.array(target["risk_table"]["times"], float); rn = np.array(target["risk_table"]["n"], float)
    from survival import step_eval
    return (step_eval(grid, rt, rn) >= min_n) & (grid <= t_max + 1e-9)


def perigon_stats(p, n, seed):
    """Velders-style statistics on a PERIGON-matched cohort (discharge reference, visits 3-6 mo and 1-5 y)."""
    prof = COHORT_PROFILES["perigon"]
    rng = np.random.default_rng(seed)
    df, t, vt, obs, ref, stages, times, t_end, ends, latent = _core(p, n, rng, profile=prof, schedule=prof.schedule_y, reference_time=prof.reference_time_y, t_max=5.5)
    m = obs["mpg"]; lows, highs, mids = [], [], []
    for j in range(1, 5):
        a, b = m[:, j], m[:, j + 1]; ok = ~np.isnan(a) & ~np.isnan(b)
        dec = pd.qcut(a[ok], 10, labels=False, duplicates="drop"); ch = pd.Series(b[ok] - a[ok]).groupby(dec).agg(["mean", "std"])
        lows.append(ch["mean"].iloc[0]); highs.append(ch["mean"].iloc[-1]); mids.append(ch["std"].iloc[len(ch) // 2])
    d5 = m[:, 5] - ref["mpg"]; d5 = d5[~np.isnan(d5)]
    val = {"lowest_decile": float(np.mean(lows)), "highest_decile": float(np.mean(highs)), "middle_sd": float(np.mean(mids)),
           "pi_width": float(np.percentile(d5, 97.5) - np.percentile(d5, 2.5)), "pi": (float(np.percentile(d5, 2.5)), float(np.percentile(d5, 97.5))),
           "change_mean": float(d5.mean())}
    for k, s in (("varc3", "ever_varc3"), ("capodanno", "ever_capodanno"), ("dvir", "ever_dvir")):
        st = stages[k]; lab2 = np.nan_to_num(st >= 2, nan=False); val[s] = float(100 * lab2.any(1).mean())
        if k == "varc3":
            first = lab2.any(1); idx = lab2.argmax(1)
            # persistence is assessed for first labels at visits that HAVE a scheduled next visit (t <= 4.5 y), as in A3 Table 2
            t_first = np.where(first, t[np.arange(len(t)), idx], np.nan)
            first = first & (t_first <= 4.5)
            nxt = np.array([st[i, idx[i] + 1] if idx[i] + 1 < st.shape[1] else np.nan for i in range(len(st))])
            val["persistence"] = float(100 * np.nan_to_num(nxt >= 2, nan=False)[first].mean()) if first.any() else np.nan
            val["next_missing"] = float(100 * np.isnan(nxt)[first].mean()) if first.any() else np.nan
            per_visit = [100 * np.nanmean(np.nan_to_num(st[:, j] >= 2, nan=False)[~np.isnan(st[:, j])]) for j in range(1, 6)]
            val["per_visit_prevalence"] = per_visit
    return val


def step_noise(p, rep, n, seed, label="noise"):
    """Grid over sigma_flow x sigma_cw on a PERIGON-matched cohort against Velders' decile / PI / ever-labelled targets."""
    best = None
    points = [(sf, sc) for sf in (0.045, 0.055, 0.065) for sc in (0.03, 0.045)]
    jobs = []
    for sf, sc in points:
        q = copy.deepcopy(p); q.noise.sigma_flow = sf; q.noise.sigma_cw = sc; jobs.append((q, n, seed))
    with _pool() as ex:
        vals = list(ex.map(_perigon_worker, jobs))
    for (sf, sc), val in zip(points, vals):
        score = (((val["lowest_decile"] - 1.7) / 1.0) ** 2 + ((val["highest_decile"] + 2.5) / 2.0) ** 2 + ((val["middle_sd"] - 3.0) / 0.7) ** 2
                 + ((val["pi_width"] - 17.1) / 3.0) ** 2 + ((val["ever_varc3"] - 3.0) / 1.0) ** 2 + ((val["ever_capodanno"] - 4.6) / 1.2) ** 2
                 + ((np.nan_to_num(val["persistence"], nan=0.0) - 33.0) / 6.0) ** 2)
        print(f"    noise grid sf={sf} sc={sc}: score {score:.2f} ever {val['ever_capodanno']:.1f}/{val['ever_dvir']:.1f}/{val['ever_varc3']:.1f} persist {val['persistence']:.0f} pi_w {val['pi_width']:.1f}")
        if best is None or score < best[0]:
            best = (score, sf, sc, val)
    _, sf, sc, val = best
    p.noise.sigma_flow = sf; p.noise.sigma_cw = sc
    rep.add(label, "sigma_flow / sigma_cw chosen", "", f"{sf} / {sc}", "", True)
    rep.add(label, "lowest-decile mean change (mmHg/interval)", "+1.2..+2.3", round(val["lowest_decile"], 2), "+/-1", 0.2 <= val["lowest_decile"] <= 3.3)
    rep.add(label, "highest-decile mean change", "-1.0..-5.9", round(val["highest_decile"], 2), "+/-1", -6.9 <= val["highest_decile"] <= 0.0)
    rep.add(label, "middle-decile SD of change", "2.4..3.9", round(val["middle_sd"], 2), "+/-0.7", 1.7 <= val["middle_sd"] <= 4.6)
    rep.add(label, "95% PI of 5-y within-patient change", "(-9.6, +7.5), width 17.1", f"({val['pi'][0]:.1f}, {val['pi'][1]:.1f}), width {val['pi_width']:.1f}", "width +/-4", abs(val["pi_width"] - 17.1) <= 4.0)
    rep.add(label, "ever labelled by 5 y, Capodanno / Dvir / VARC-3 (%)", "4.6 / 2.9 / 3.0", f"{val['ever_capodanno']:.1f} / {val['ever_dvir']:.1f} / {val['ever_varc3']:.1f}", "+/-1.5",
            abs(val["ever_capodanno"] - 4.6) <= 1.5 and abs(val["ever_dvir"] - 2.9) <= 1.5 and abs(val["ever_varc3"] - 3.0) <= 1.5)
    rep.add(label, "VARC-3 label persistence at next visit (%)", "~33 (published range 14-50)", round(val["persistence"], 0), "14-50", 14 <= val["persistence"] <= 50,
            f"next visit missing {val['next_missing']:.0f}% (Velders 15-25%); per-visit prevalence " + "/".join(f"{x:.1f}" for x in val["per_visit_prevalence"]) + " (Velders 0.2-1.4)")


def step_death(p, rep, n, seed):
    target = TARGETS["kermen_fig2a_overall_survival"]
    grid = np.arange(1.0, min(target["t_end"], 12.0) + 1e-9, 1.0)
    tv = target_curve(target, grid)
    prof = COHORT_PROFILES["kermen"]
    def sim_surv(a):
        p.death.a = a
        rng = np.random.default_rng(seed)
        df, t, vt, obs, ref, stages, times, t_end, ends, latent = _core(p, n, rng, profile=prof, schedule=prof.schedule_y, reference_time=prof.reference_time_y, t_max=12.5)
        td = ends["t_death"]
        tt = np.minimum(td, 12.5); st = (td <= 12.5).astype(int)
        return observed_curve(tt, st, grid, "KM")
    res = minimize_scalar(lambda la: float(np.sum((sim_surv(np.exp(la)) - tv)[:10] ** 2)), bounds=(np.log(0.002), np.log(0.2)), method="bounded", options={"xatol": 0.02})
    p.death.a = float(np.exp(res.x))
    s = sim_surv(p.death.a)
    diff = float(np.max(np.abs(s - tv)[:10]))
    rep.add("death", "Kermen overall survival, max |sim - target| 1-10 y (pp)", "", round(diff, 2), 3.0, diff <= 3.0, f"a = {p.death.a:.4f}; sim 5/10 y = {s[4]:.1f}/{s[9]:.1f}, target {tv[4]:.1f}/{tv[9]:.1f}")
    # NOTION check: ~63% dead at 10 y at age 79
    rng = np.random.default_rng(seed + 1); prof = COHORT_PROFILES["notion_tavi"]
    df, t, vt, obs, ref, stages, times, t_end, ends, latent = _core(p, n, rng, profile=prof, schedule=prof.schedule_y, reference_time=prof.reference_time_y, t_max=10.0)
    dead10 = float(100 * np.mean(ends["t_death"] <= 10.0))
    rep.add("death", "NOTION all-cause mortality at 10 y (%), age 79", 63.0, round(dead10, 1), 6.0, abs(dead10 - 63.0) <= 6.0, "check, not fitted")


def fit_class(p, tid, paper, cls, label, n, seed, bounds_lm=(-1.0, 1.0), bounds_lk=(-0.7, 0.9), shape_free=True):
    """Fit (log scale multiplier, log shape multiplier) of one class's Weibull onset to a curve target:
    risk-weighted SSE on the yearly grid (<= 10 y, n >= 30), common random numbers; the coarse grid is
    evaluated in parallel processes, then Nelder-Mead refines sequentially with the parameters clipped."""
    from scipy.optimize import minimize
    target = TARGETS[tid]
    grid = np.arange(1.0, target["t_end"] + 1e-9, 1.0)
    tv = target_curve(target, grid); w = risk_weights(target, grid) * risk_ok(target, grid)
    def clipped(theta):
        return (0.0, float(np.clip(theta[0], *bounds_lm)), float(np.clip(theta[1], *bounds_lk)) if shape_free else 0.0)
    def curve(theta):
        return _curve_worker((p, paper, cls, label, tid, clipped(theta), n, seed, grid))
    sse = lambda c: float(np.sum(w * (c - tv) ** 2))
    lks = (-0.5, -0.25, 0.0, 0.25, 0.5) if shape_free else (0.0,)
    pts = [(0.0, 0.0)] + [(lm, lk) for lm in (-0.8, -0.4, 0.0, 0.4, 0.8) for lk in lks]
    with _pool() as ex:
        curves = list(ex.map(_curve_worker, [(p, paper, cls, label, tid, clipped(th), n, seed, grid) for th in pts]))
    direct = curves[0]
    scores = [sse(c) for c in curves[1:]]
    best = (min(scores), pts[1 + int(np.argmin(scores))])
    r = minimize(lambda th: sse(curve(th)), np.array(best[1], float), method="Nelder-Mead", options={"xatol": 0.02, "fatol": 0.05, "maxiter": 60})
    th = tuple(float(x) for x in r.x) if r.fun <= best[0] else best[1]
    th = (float(np.clip(th[0], *bounds_lm)), float(np.clip(th[1], *bounds_lk)) if shape_free else 0.0)
    sim = curve(th)
    return {"theta": th, "sim": sim, "target": tv, "grid": grid, "mask": risk_ok(target, grid),
            "direct_residual": float(np.max(np.abs(direct - tv)[risk_ok(target, grid)]))}


def apply_theta(p, cls, th):
    p.onset.scale[cls] = float(p.onset.scale[cls] * np.exp(th[0])); p.onset.shape[cls] = float(p.onset.shape[cls] * np.exp(th[1]))


def step_onset(p, rep, n, seed, quick=False):
    thetas = {}
    for tid, paper, cls, label, tol in CURVE_STEPS:
        t0 = time.time()
        res = fit_class(p, tid, paper, cls, label, n, seed)
        apply_theta(p, cls, res["theta"]); thetas[cls] = res["theta"]
        diff = float(np.max(np.abs(res["sim"] - res["target"])[res["mask"]]))
        rep.add("onset", f"{tid} ({cls}), max |sim - target| while n>=30, <=10 y (pp)", "", round(diff, 2), tol, diff <= tol,
                f"scale x{np.exp(res['theta'][0]):.2f} -> {p.onset.scale[cls]:.1f} y, shape x{np.exp(res['theta'][1]):.2f} -> {p.onset.shape[cls]:.2f}; prior residual {res['direct_residual']:.1f} pp; {time.time()-t0:.0f}s")
    p.onset.scale["new_savr"] = p.onset.scale["perimount"]; p.onset.shape["new_savr"] = p.onset.shape["perimount"]
    p.onset.scale["new_thv"] = p.onset.scale["self_expanding"]; p.onset.shape["new_thv"] = p.onset.shape["self_expanding"]
    return thetas


def step_notion_savr(p, rep, n, seed):
    """Cross-check: predict NOTION SAVR from the calibrated Perimount and Trifecta classes (age HR transports
    them to 79 y); then free only the porcine/Mitroflow class to close the residual."""
    tid = "notion_fig3_modsvd_savr"; target = TARGETS[tid]; censor = censoring_sampler(tid)
    grid, sim, tv = matched_curve(p, "notion_savr", "varc3_haemo_single2", target, n, seed, censor=censor)
    diff_pred = float(np.max(np.abs(sim - tv)[risk_ok(target, grid)]))
    rep.add("cross-check", "NOTION SAVR >=moderate SVD PREDICTED before freeing the porcine class, max |diff| (pp)", "", round(diff_pred, 2), 5.0, diff_pred <= 5.0,
            f"sim 5/10 y = {sim[4]:.1f}/{sim[9]:.1f}, target {tv[4]:.1f}/{tv[9]:.1f}")
    res = fit_class(p, tid, "notion_savr", "porcine", "varc3_haemo_single2", n, seed, bounds_lm=(-0.8, 0.8), bounds_lk=(-0.5, 0.7))
    apply_theta(p, "porcine", res["theta"])
    diff = float(np.max(np.abs(res["sim"] - res["target"])[res["mask"]]))
    rep.add("onset", f"{tid} (porcine share freed), max |sim - target| (pp)", "", round(diff, 2), 3.0, diff <= 3.0,
            f"porcine scale -> {p.onset.scale['porcine']:.1f} y, shape -> {p.onset.shape['porcine']:.2f}; sim 5/10 y = {res['sim'][4]:.1f}/{res['sim'][9]:.1f}")


def step_balloon(p, rep, n, seed):
    """Balloon-expandable class: scale chosen so that the PARTNER-style 5-y VARC-3 SVD point targets are met
    (Sapien 3 0.7%, XT 1.6%; one class, so the target is their weighted mean at the two cohorts)."""
    def cif5(lm, paper):
        rng = np.random.default_rng(seed)
        tt, st = simulate_matched(paper, {"balloon": (0.0, lm)}, n, 5.0, rng, p, label="varc3_single2")
        return float(observed_curve(tt, st, np.array([5.0]), "AJ")[0])      # already in percent
    # fitted on Sapien XT (PARTNER 2A, 1.6%); Sapien 3 (0.7%) is a check because it sits at the single-echo noise floor
    res = minimize_scalar(lambda lm: (cif5(lm, "partner_xt") - 1.6) ** 2, bounds=(-0.7, 0.7), method="bounded", options={"xatol": 0.05})
    p.onset.scale["balloon"] *= float(np.exp(res.x))
    s3, xt = cif5(0.0, "partner_s3"), cif5(0.0, "partner_xt")
    rep.add("onset", "Sapien XT VARC-3 SVD at 5 y (%), fitted", 1.6, round(xt, 2), 0.8, abs(xt - 1.6) <= 0.8, f"balloon scale x{np.exp(res.x):.2f} -> {p.onset.scale['balloon']:.1f} y")
    rep.add("check", "Sapien 3 VARC-3 SVD at 5 y (%)", 0.7, round(s3, 2), 2.0, abs(s3 - 0.7) <= 2.0,
            "single-echo VARC-3 labels run 3% by 5 y even in PERIGON (A3), so ~1-2% is the noise floor of the definition")


def step_checks(p, rep, n, seed):
    with _pool() as ex:
        results = dict(ex.map(_check_worker, [(p, tid, paper, label, n, seed) for tid, paper, label, tol in CHECKS]))
    for tid, paper, label, tol in CHECKS:
        target = TARGETS[tid]
        grid, sim, tv = results[tid]
        m = risk_ok(target, grid)
        diff = float(np.max(np.abs(sim - tv)[m])) if m.any() else float(np.max(np.abs(sim - tv)))
        at = min(len(grid), 10) - 1
        rep.add("check", f"{tid}: max |sim - target| while n>=30 (pp)", "", round(diff, 2), tol, diff <= tol,
                f"sim at {grid[at]:.0f} y = {sim[at]:.1f}, target {tv[at]:.1f}")
    # reintervention at 10 y in NOTION (2.2-4.3 %)
    rng = np.random.default_rng(seed); prof = COHORT_PROFILES["notion_savr"]
    df, t, vt, obs, ref, stages, times, t_end, ends, latent = _core(p, n, rng, profile=prof, schedule=prof.schedule_y, reference_time=prof.reference_time_y, t_max=10.0)
    r10 = float(100 * np.mean(np.isfinite(ends["t_re"]) & (ends["t_re"] <= 10)))
    rep.add("check", "NOTION SAVR reintervention by 10 y (%), age 79", "2.2-4.3", round(r10, 2), "", 1.0 <= r10 <= 7.0)


def run_all(n_sim=20000, seed=0, quick=False):
    OUT.mkdir(exist_ok=True)
    p = default_params(); rep = Report()
    p.onset.shift = {c: 0.0 for c in p.onset.scale}
    resid = fit_physics_constants(p)
    worst = max(abs(v) for d in resid.values() for v in d.values())
    cs = np.array(list(p.physics.c_size.values())) * 1e4
    rep.add("physics", "c fitted per model x size to the ASE rows (x1e-4; Gorlin = 5.1)", "", f"range {cs.min():.2f}-{cs.max():.2f}", "", True,
            f"a single per-model constant would leave residuals up to {worst:.1f} mmHg at the smallest sizes")
    n = n_sim if not quick else 6000
    T0 = time.time(); print(f"[setup] {N_WORKERS} worker processes")
    step_noise(p, rep, n, seed, label="noise (pass 1, prior onset)")
    step_death(p, rep, n, seed)
    step_onset(p, rep, n, seed, quick)
    # the ever-labelled targets depend weakly on the calibrated Perimount onset: one more pass of each
    step_noise(p, rep, n, seed, label="noise (pass 2)")
    tid, paper, cls, label_, tol = CURVE_STEPS[0]
    res = fit_class(p, tid, paper, cls, label_, n, seed)
    apply_theta(p, cls, res["theta"])
    p.onset.scale["new_savr"] = p.onset.scale["perimount"]; p.onset.shape["new_savr"] = p.onset.shape["perimount"]
    diff = float(np.max(np.abs(res["sim"] - res["target"])[res["mask"]]))
    rep.add("onset (pass 2)", f"{tid} (perimount), max |sim - target| while n>=30 (pp)", "", round(diff, 2), tol, diff <= tol,
            f"scale {p.onset.scale['perimount']:.2f} y, shape {p.onset.shape['perimount']:.2f}; sim 5/10 y = {res['sim'][4]:.1f}/{res['sim'][9]:.1f} vs {res['target'][4]:.1f}/{res['target'][9]:.1f}")
    step_notion_savr(p, rep, n, seed)
    step_balloon(p, rep, n, seed)
    step_checks(p, rep, n, seed)
    rep.add("result", "calibration wall time (s)", "", round(time.time() - T0), "", True, f"{N_WORKERS} processes, n_sim {n}")
    for cls in p.onset.scale:
        rep.add("result", f"onset {cls}: Weibull shape / scale / shift", "", f"{p.onset.shape[cls]:.2f} / {p.onset.scale[cls]:.2f} y / {p.onset.shift[cls]:+.2f} y", "", True)
    rep.add("result", "noise sigma_flow / sigma_cw / sigma_lvot_vti / sigma_lvot_diam", "", f"{p.noise.sigma_flow} / {p.noise.sigma_cw} / {p.noise.sigma_lvot_vti} / {p.noise.sigma_lvot_diam}", "", True)
    rep.add("result", "death a (Gompertz level at age 75)", "", f"{p.death.a:.4f}", "", True)
    p.to_json(OUT / "params_calibrated.json")
    (OUT / "calibration_report.md").write_text(rep.markdown())
    print(rep.markdown())
    return p, rep


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-sim", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    run_all(a.n_sim, a.seed, a.quick)
