"""Simulation-based calibration of the simulator's LATENT onset distribution to the reconstructed
OBSERVED-label curves (design choice D6 in docs/04_synthetic_cohort_spec.md).

Why not fit the IPD directly: the reconstructed curves are echo-detected, visit-timed labels. Using them
as the latent onset distribution would make the simulator's own observed curve later than the paper's by
the onset-to-threshold lag. So the reconstructed curve is the *target*; the latent distribution is
adjusted until the simulated observed curve (same estimator, same risk-set attrition) matches it.

Interface for the real simulator (notebooks/simulator, Idea 1):
    simulate_fn(theta, n, t_max, rng) -> (times, states)   states: 0 censored, 1 SVD label, 2 death
    theta = (shift_years, log_scale_multiplier) applied to the stratum's latent onset distribution.
The built-in `default_simulate` is a deliberately minimal stand-in (latent Weibull onset, fixed
onset-to-threshold lag, annual echo, exponential death, administrative censoring drawn from the
reconstructed IPD) so the pipeline runs end-to-end today.

Usage:
    python calibrate_onset.py [--targets id ...] [--lag 1.5] [--n-sim 4000]
Writes output/calibrated_onset.json: per target, the direct-parametric residual (theta = 0, i.e. the
Weibull fit to the IPD used as the latent onset) and the calibrated theta with its residual.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from survival import aj, km, step_eval  # noqa: E402

OUT = HERE / "output"
SVD_TARGETS = ["kermen_fig2b_freedom_stage23", "kermen_fig2b_freedom_stage3", "kermen_fig4a_explant_svd",
               "notion_fig3_modsvd_tavi", "notion_fig3_modsvd_savr", "notion_fig3_sevsvd_savr", "wakami_fig1_svd"]

# constant all-cause death hazard per year used by the stand-in simulator, from the papers' own numbers
DEATH_RATE = {
    "kermen": -np.log(0.667) / 10.0,       # 66.7% survival at 10 y (Fig 2A)
    "notion": -np.log(0.40) / 10.0,        # ~60% dead at 10 y, age 79 (NOTION 10-y)
    "wakami": -np.log(1 - 0.238) / 7.0,    # 23.8% all-cause death at 7 y (Fig 1 table)
}


def load_targets(ids=None):
    t = json.loads((OUT / "targets.json").read_text())
    p = json.loads((OUT / "onset_priors.json").read_text())
    if ids:
        t = {k: v for k, v in t.items() if k in ids}
    return t, p


def censoring_sampler(target_id):
    """Administrative censoring times from the reconstructed IPD (so simulated risk sets shrink like the paper's)."""
    df = pd.read_csv(OUT / f"ipd_{target_id}.csv")
    c = df.loc[df.event == 0, "time"].to_numpy(float)
    if len(c) == 0:
        c = np.array([df.time.max()])
    return lambda n, rng: rng.choice(c, size=n, replace=True)


def default_simulate(theta, n, t_max, rng, *, weibull, lag, death_rate, censor):
    """Stand-in simulator. Latent onset ~ Weibull(shape, scale*exp(theta[1])) + theta[0] years; the label
    is assigned at the first annual echo at or after onset + lag; death competes; administrative censoring."""
    shift, log_mult = theta
    k, lam = weibull["shape"], weibull["scale"] * float(np.exp(log_mult))
    onset = lam * rng.weibull(k, size=n) + shift
    label_t = np.ceil(np.maximum(onset + lag, 0.0))        # annual echo visits at 1, 2, 3, ...
    label_t = np.where(label_t < 1.0, 1.0, label_t)
    death_t = rng.exponential(1.0 / death_rate, size=n)
    cens_t = np.minimum(censor(n, rng), t_max)
    t = np.minimum.reduce([label_t, death_t, cens_t])
    state = np.where(t == label_t, 1, np.where(t == death_t, 2, 0))
    return t, state


def observed_curve(times, states, grid, estimator):
    """Same estimator as the target: KM with death censored (Kermen Fig 2B) or Aalen-Johansen CIF."""
    if estimator.startswith("KM"):
        tt, S = km(times, (states == 1).astype(int))
        return 100.0 * step_eval(grid, tt, S)                       # freedom from SVD, %
    cif = aj(times, states)
    if 1 not in cif:
        return np.zeros_like(grid)
    tt, C = cif[1]
    return 100.0 * step_eval(grid, tt, C)                          # cumulative incidence, %


def target_curve(target, grid):
    v = np.array(step_eval(grid, np.array(target["grid_t"]), np.array(target["grid_value_pct"]) / 100.0)) * 100.0
    return v


def risk_weights(target, grid):
    rt = np.array(target["risk_table"]["times"], float); rn = np.array(target["risk_table"]["n"], float)
    w = step_eval(grid, rt, rn)
    return w / w.max()


def objective(theta, target, sim, grid, n_sim, seed, estimator):
    """Risk-weighted SSE between the simulated observed curve and the target on the yearly grid.
    One large cohort with common random numbers (same seed for every theta) keeps the objective smooth."""
    tv = target_curve(target, grid); w = risk_weights(target, grid)
    rng = np.random.default_rng(seed)
    t, st = sim(theta, n_sim, target["t_end"], rng)
    return float(np.sum(w * (observed_curve(t, st, grid, estimator) - tv) ** 2))


def calibrate(target, sim, estimator, n_sim=20000, seed=0, grid_step=1.0):
    grid = np.arange(grid_step, target["t_end"] + 1e-9, grid_step)
    f = lambda th: objective(th, target, sim, grid, n_sim, seed, estimator)
    direct = f((0.0, 0.0))
    best = None
    for shift in (-4, -3, -2, -1, 0, 1, 2):
        for lm in (-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6):
            v = f((shift, lm))
            if best is None or v < best[0]:
                best = (v, (shift, lm))
    r = minimize(f, np.array(best[1], float), method="Nelder-Mead", options={"xatol": 0.02, "fatol": 0.1, "maxiter": 200})
    theta = tuple(float(x) for x in r.x) if r.fun <= best[0] else best[1]
    tv = target_curve(target, grid)
    rng = np.random.default_rng(seed + 1)                                # fresh draws for the reported residual
    t, st = sim(theta, n_sim, target["t_end"], rng)
    fit = observed_curve(t, st, grid, estimator)
    t0, st0 = sim((0.0, 0.0), n_sim, target["t_end"], rng)
    fit0 = observed_curve(t0, st0, grid, estimator)
    return {
        "theta": {"shift_years": theta[0], "scale_multiplier": float(np.exp(theta[1]))},
        "objective_direct": direct, "objective_calibrated": float(min(r.fun, best[0])),
        "grid_t": grid.tolist(), "target_pct": tv.tolist(),
        "direct_parametric_pct": fit0.tolist(), "calibrated_pct": fit.tolist(),
        "max_abs_diff_direct_pp": float(np.max(np.abs(fit0 - tv))),
        "max_abs_diff_calibrated_pp": float(np.max(np.abs(fit - tv))),
    }


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="*", default=SVD_TARGETS, help="default: the SVD curves only (onset has no meaning for mortality curves)")
    ap.add_argument("--lag", type=float, default=1.5, help="onset-to-label lag in years used by the stand-in simulator")
    ap.add_argument("--n-sim", type=int, default=20000)
    args = ap.parse_args(argv)
    targets, priors = load_targets(args.targets)
    out = {}
    for tid, target in targets.items():
        if tid not in priors:
            print(f"{tid}: no parametric prior (too few events); skipped"); continue
        w = priors[tid]["weibull"]
        paper = "kermen" if tid.startswith("kermen") else "notion" if tid.startswith("notion") else "wakami"
        estimator = target["estimator"]
        cens = censoring_sampler(tid)
        sim = lambda th, n, tmax, rng, w=w, cens=cens, paper=paper: default_simulate(
            th, n, tmax, rng, weibull=w, lag=args.lag, death_rate=DEATH_RATE[paper], censor=cens)
        res = calibrate(target, sim, estimator, n_sim=args.n_sim)
        res["estimator"] = estimator; res["stratum"] = target["stratum"]; res["prior_weibull"] = w
        out[tid] = res
        print(f"{tid:36s} direct max|diff| {res['max_abs_diff_direct_pp']:5.1f} pp -> calibrated {res['max_abs_diff_calibrated_pp']:5.1f} pp "
              f"(shift {res['theta']['shift_years']:+.2f} y, scale x{res['theta']['scale_multiplier']:.2f})")
    (OUT / "calibrated_onset.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
