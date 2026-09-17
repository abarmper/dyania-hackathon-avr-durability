"""Metrics for competing risks on fixed-horizon eligible subsets, cluster bootstrap by valve, label sweep."""
from __future__ import annotations

import numpy as np
import pandas as pd
from lifelines.utils import concordance_index
from sklearn.metrics import brier_score_loss, roc_auc_score

import analysis_config as C


def _case_control(inst: pd.DataFrame, tau: float, label_time: np.ndarray | None = None):
    """Eligible rows and the case indicator at horizon tau.
    cases: SVD known <= tau (or, for a swapped label, label_time <= landmark + tau: prevalent-or-incident);
    controls: follow-up beyond tau, or a competing event <= tau; rows censored before tau are dropped (LTFU)."""
    tt = inst.time_to.to_numpy(); cause = inst.cause.to_numpy()
    if label_time is None:
        case = (cause == 1) & (tt <= tau + 1e-9)
        elig = case | (tt > tau) | ((cause == 2) & (tt <= tau + 1e-9))
    else:
        lt = label_time - inst.landmark_t.to_numpy()
        case = np.isfinite(lt) & (lt <= tau + 1e-9)
        elig = case | (tt > tau) | ((cause == 2) & (tt <= tau + 1e-9))
    return elig, case


def auc_at(inst, cif_tau, tau, label_time=None):
    elig, case = _case_control(inst, tau, label_time)
    if case[elig].sum() < 5 or (~case[elig]).sum() < 5:
        return np.nan
    return float(roc_auc_score(case[elig], cif_tau[elig]))


def brier_at(inst, cif_tau, tau):
    elig, case = _case_control(inst, tau)
    return float(brier_score_loss(case[elig], np.clip(cif_tau[elig], 1e-6, 1 - 1e-6)))


def c_index_trunc(inst, cif_tau, tau):
    """Wolbers' cause-specific C truncated at tau: T = min(T*, tau), E = SVD within tau; competing events censored at tau.
    lifelines scores higher `predicted` = longer time, so pass -CIF."""
    tt = np.minimum(inst.time_to.to_numpy(), tau); e = ((inst.cause.to_numpy() == 1) & (inst.time_to.to_numpy() <= tau + 1e-9)).astype(int)
    if e.sum() < 5:
        return np.nan
    return float(concordance_index(tt, -cif_tau, e))


def calibration_at(inst, cif_tau, tau, n_bins=10):
    """Decile table (predicted vs observed fraction with SVD known by tau on the eligible subset) and the logistic recalibration slope."""
    elig, case = _case_control(inst, tau)
    p = np.clip(cif_tau[elig], 1e-6, 1 - 1e-6); y = case[elig].astype(int)
    q = pd.qcut(p, n_bins, labels=False, duplicates="drop")
    tab = pd.DataFrame({"bin": q, "pred": p, "obs": y}).groupby("bin").agg(n=("obs", "size"), pred=("pred", "mean"), obs=("obs", "mean")).reset_index()
    from sklearn.linear_model import LogisticRegression
    lr = LogisticRegression(C=1e6).fit(np.log(p / (1 - p))[:, None], y)
    return tab, float(lr.coef_[0, 0]), float(lr.intercept_[0])


def valve_resamples(valve_ids: np.ndarray, B: int, seed: int):
    """Fixed list of valve-level resamples (shared across models so that CIs of differences are paired)."""
    uniq = np.unique(valve_ids); rng = np.random.default_rng(seed)
    return [rng.choice(uniq, size=len(uniq), replace=True) for _ in range(B)]


def bootstrap(fn, inst: pd.DataFrame, resamples, extra=None):
    """fn(rows_index_array) -> float; resample valves with replacement (rows of a resampled valve appear as many times as drawn)."""
    vid = inst.valve_id.to_numpy(); order = np.argsort(vid); vs = vid[order]
    starts = np.searchsorted(vs, np.unique(vs)); ends = np.r_[starts[1:], len(vs)]
    pos = {v: (s, e) for v, s, e in zip(np.unique(vs), starts, ends)}
    vals = []
    for rs in resamples:
        idx = np.concatenate([order[pos[v][0]:pos[v][1]] for v in rs if v in pos])
        vals.append(fn(idx))
    vals = np.array(vals, float)
    return float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


def metrics_table(inst: pd.DataFrame, cif_by_tau: dict, taus, resamples=None, subset_mask=None) -> pd.DataFrame:
    """AUC, C, Brier, calibration slope at each tau with valve-bootstrap CIs; optionally on a subset of instances."""
    rows = []
    m = np.ones(len(inst), bool) if subset_mask is None else subset_mask
    sub = inst[m].reset_index(drop=True)
    for tau in taus:
        cif = cif_by_tau[tau][m]
        auc = auc_at(sub, cif, tau); cidx = c_index_trunc(sub, cif, tau); br = brier_at(sub, cif, tau)
        tab, slope, icpt = calibration_at(sub, cif, tau)
        r = {"tau": tau, "n_eligible": int(_case_control(sub, tau)[0].sum()), "n_cases": int(_case_control(sub, tau)[1].sum()),
             "AUC": auc, "C_trunc": cidx, "Brier": br, "calib_slope": slope, "calib_intercept": icpt}
        if resamples is not None:
            r["AUC_lo"], r["AUC_hi"] = bootstrap(lambda idx: auc_at(sub.iloc[idx].reset_index(drop=True), cif[idx], tau), sub, resamples)
            r["C_lo"], r["C_hi"] = bootstrap(lambda idx: c_index_trunc(sub.iloc[idx].reset_index(drop=True), cif[idx], tau), sub, resamples)
        rows.append(r)
    return pd.DataFrame(rows)


def label_sweep(inst: pd.DataFrame, cif_tau: np.ndarray, tau: float, labels: pd.DataFrame, definitions) -> pd.DataFrame:
    """Same instances, same predictions; only the outcome indicator changes (prevalent-or-incident by landmark + tau)."""
    lab = labels.set_index("valve_id")
    rows = []
    for d in definitions:
        col = C.LABELS[d][0]
        lt = lab[col].reindex(inst.valve_id.to_numpy()).to_numpy(float)
        rows.append({"label": d, "description": C.LABELS[d][1], "AUC": auc_at(inst, cif_tau, tau, label_time=lt),
                     "n_cases": int(_case_control(inst, tau, lt)[1].sum())})
    return pd.DataFrame(rows)


def subgroup_table(inst: pd.DataFrame, cif_tau: np.ndarray, tau: float, pv: pd.DataFrame) -> pd.DataFrame:
    p = pv.set_index("valve_id")
    g = {"approach": p.approach.reindex(inst.valve_id).to_numpy(), "valve_class": p.valve_class.reindex(inst.valve_id).to_numpy(),
         "age_band": pd.cut(p.age_at_implant.reindex(inst.valve_id), [0, 60, 80, 200], labels=["<60", "60-80", ">80"]).astype(str).to_numpy(),
         "size": p.label_size_mm.reindex(inst.valve_id).to_numpy(), "sex": p.sex.reindex(inst.valve_id).to_numpy(),
         "reference_source": p.reference_source.reindex(inst.valve_id).to_numpy()}
    rows = []
    for name, arr in g.items():
        for lev in pd.unique(arr):
            m = arr == lev
            sub = inst[m].reset_index(drop=True)
            if m.sum() < 200:
                continue
            tab, slope, _ = calibration_at(sub, cif_tau[m], tau) if _case_control(sub, tau)[1].sum() >= 10 else (None, np.nan, np.nan)
            rows.append({"factor": name, "level": lev, "n_instances": int(m.sum()), "n_cases": int(_case_control(sub, tau)[1].sum()),
                         "AUC": auc_at(sub, cif_tau[m], tau), "C_trunc": c_index_trunc(sub, cif_tau[m], tau), "calib_slope": slope})
    return pd.DataFrame(rows)
