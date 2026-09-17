"""Landmark features. ONE vectorised function, `history_features`, serves both the training table (looping over the
echo ordinal k) and the online policy loop (called with the history so far), so train and serve cannot diverge."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm

import analysis_config as C
from data import Hist, baseline_frame
from labels import varc3_stage  # simulator's rule (notebooks/simulator/labels.py)


def _stage_haemo(mpg, ar, ref_mpg, ref_ar):
    """Single-echo VARC-3 haemodynamic stage (gradient / AR only) for arrays of any shape."""
    with np.errstate(invalid="ignore"):
        d = mpg - ref_mpg
        st3 = ((mpg >= 30) & (d >= 20)) | ((ar >= 3) & (ref_ar < 3))
        st2 = ((mpg >= 20) & (d >= 10)) | ((ar >= 2) & (ref_ar < 2))
    out = np.where(np.nan_to_num(st3, nan=False), 3.0, np.where(np.nan_to_num(st2, nan=False), 2.0, 0.0))
    return np.where(np.isnan(mpg), np.nan, out)


def history_features(h: Hist, k: int, base: pd.DataFrame, rows: np.ndarray | None = None) -> pd.DataFrame:
    """Features at the k-th post-reference echo (1-based) for the valves in `rows` (default: all with >= k echoes).
    Uses echoes 1..k only. `base` = baseline_frame(pv) aligned with h.valve_ids."""
    if rows is None:
        rows = np.where(h.n_echoes >= k)[0]
    if len(rows) == 0:
        return pd.DataFrame()
    j = k - 1
    b = base.iloc[rows]
    t = h.t[rows, :k]; mpg = np.minimum(h.mpg[rows, :k], C.MPG_CAP); eoa = h.eoa[rows, :k]; dvi = h.dvi[rows, :k]; ar = h.ar[rows, :k]; vt = h.vtype[rows, :k]
    ref_mpg = b.ref_mpg.to_numpy()[:, None]; ref_eoa = b.ref_eoa.to_numpy()[:, None]; ref_dvi = b.ref_dvi.to_numpy()[:, None]; ref_ar = b.ref_ar.to_numpy()[:, None]
    t_ref = 0.08
    tk = t[:, j]; mk = mpg[:, j]
    with np.errstate(invalid="ignore", divide="ignore"):
        # exponentially weighted mean of log MPG over the history (reference included when it was an echo)
        lm = np.log(np.maximum(mpg, 0.5)); w = np.exp(-(tk[:, None] - t) / C.EWMA_TAU_Y); w = np.where(np.isnan(mpg), 0.0, w)
        has_ref_echo = (b.reference_from_table.to_numpy() == 0)
        w_ref = np.where(has_ref_echo, np.exp(-(tk - t_ref) / C.EWMA_TAU_Y), 0.0)
        lm_ref = np.log(np.maximum(ref_mpg[:, 0], 0.5))
        sw = w.sum(1) + w_ref; sw2 = (w ** 2).sum(1) + w_ref ** 2
        mean_log = (np.nansum(w * lm, 1) + w_ref * lm_ref) / sw
        n_eff = sw ** 2 / np.maximum(sw2, 1e-12)
        p_ge20 = norm.cdf((mean_log - np.log(20.0)) / (C.LOG_MPG_SD / np.sqrt(np.maximum(n_eff, 1e-9))))
        # least-squares slope of MPG on time over reference + echoes 1..k (NaN with < 2 post-reference echoes)
        T = np.concatenate([np.full((len(rows), 1), t_ref), t], 1); M = np.concatenate([ref_mpg, mpg], 1)
        ok = ~np.isnan(M)
        n_ok = ok.sum(1); Tm = np.where(ok, T, 0).sum(1) / np.maximum(n_ok, 1); Mm = np.where(ok, M, 0).sum(1) / np.maximum(n_ok, 1)
        cov = np.where(ok, (T - Tm[:, None]) * (M - Mm[:, None]), 0).sum(1); var = np.where(ok, (T - Tm[:, None]) ** 2, 0).sum(1)
        slope = np.where((k >= 2) & (var > 1e-9), cov / np.maximum(var, 1e-9), np.nan)
        prev_m = mpg[:, j - 1] if k >= 2 else ref_mpg[:, 0]; prev_t = t[:, j - 1] if k >= 2 else np.full(len(rows), t_ref)
        slope_last2 = (mk - prev_m) / np.maximum(tk - prev_t, 0.05)
        max_sofar = np.nanmax(mpg, 1)
        rel_eoa = eoa[:, j] / ref_eoa[:, 0] - 1.0; rel_dvi = dvi[:, j] / ref_dvi[:, 0] - 1.0
        # stages on every past echo (single-echo rules) -> unconfirmed flag, prior positives, transient spike
        st_full, st_h = varc3_stage(mpg, eoa, dvi, ar, ref_mpg, ref_eoa, ref_dvi, ref_ar)
        st_h_all = _stage_haemo(mpg, ar, ref_mpg, ref_ar)
        pos_hist = np.nan_to_num(st_h_all >= 2, nan=False)
        prior_pos = pos_hist[:, :j].sum(1); pos_now = pos_hist[:, j]
        # consecutive worsening: how many steps ending at k had a higher MPG than the step before (reference before echo 1)
        M2 = np.concatenate([ref_mpg, mpg], 1); worse = np.nan_to_num(np.diff(M2, axis=1) > 0, nan=False)
        cw = np.zeros(len(rows), int); alive = np.ones(len(rows), bool)
        for jj in range(j, -1, -1):
            alive &= worse[:, jj]; cw += alive
        any_sym_before = (vt[:, :j] == 2).any(1) if j > 0 else np.zeros(len(rows), bool)
    f = pd.DataFrame({
        "valve_id": b.valve_id.to_numpy(), "k": k, "landmark_t": tk, "n_echoes": k,
        "time_since_last": tk - prev_t, "mpg_now": mk, "eoa_now": eoa[:, j], "dvi_now": dvi[:, j], "vmax_now": h.vmax[rows, j],
        "svi_now": h.svi[rows, j], "lvef_now": h.lvef[rows, j], "ar_now": ar[:, j],
        "delta_mpg": mk - ref_mpg[:, 0], "rel_eoa": rel_eoa, "rel_dvi": rel_dvi,
        "ewma_log_mpg": mean_log, "p_true_mpg_ge20": p_ge20, "slope_since_ref": slope, "slope_last2": slope_last2, "max_mpg_sofar": max_sofar,
        "stage_now_full": np.nan_to_num(st_full[:, j], nan=-1.0), "stage_now_haemo": st_h_all[:, j],
        "stage_now_unevaluable": np.isnan(st_full[:, j]).astype(int),
        "unconfirmed_flag": pos_now.astype(int), "prior_positive_count": prior_pos, "consecutive_worsening": cw,
        "transient_spike": ((prior_pos > 0) & ~pos_now).astype(int),
        "symptom_now": (vt[:, j] == 2).astype(int), "any_symptom_before": any_sym_before.astype(int),
        "eoa_missing_now": np.isnan(eoa[:, j]).astype(int), "dvi_missing_now": np.isnan(dvi[:, j]).astype(int),
    })
    return pd.concat([f.reset_index(drop=True), b.drop(columns=["valve_id"]).reset_index(drop=True)], axis=1)


def build_landmark_table(h: Hist, pv: pd.DataFrame) -> pd.DataFrame:
    base = baseline_frame(pv).set_index("valve_id").loc[h.valve_ids].reset_index()
    parts = [history_features(h, k, base) for k in range(1, h.K + 1)]
    df = pd.concat([p for p in parts if len(p)], ignore_index=True)
    return df.sort_values(["valve_id", "k"]).reset_index(drop=True)


def implant_instances(pv: pd.DataFrame) -> pd.DataFrame:
    """One instance per valve at the reference echo (landmark_t = 0.08, no longitudinal history)."""
    b = baseline_frame(pv)
    f = pd.DataFrame({"valve_id": b.valve_id, "k": 0, "landmark_t": 0.08, "n_echoes": 0, "time_since_last": 0.08,
                      "mpg_now": b.ref_mpg, "eoa_now": b.ref_eoa, "dvi_now": b.ref_dvi, "vmax_now": np.nan, "svi_now": np.nan, "lvef_now": np.nan,
                      "ar_now": b.ref_ar, "delta_mpg": 0.0, "rel_eoa": 0.0, "rel_dvi": 0.0, "ewma_log_mpg": np.log(np.maximum(b.ref_mpg, 0.5)),
                      "p_true_mpg_ge20": norm.cdf((np.log(np.maximum(b.ref_mpg, 0.5)) - np.log(20)) / C.LOG_MPG_SD),
                      "slope_since_ref": np.nan, "slope_last2": np.nan, "max_mpg_sofar": b.ref_mpg, "stage_now_full": 0.0, "stage_now_haemo": 0.0,
                      "stage_now_unevaluable": 0, "unconfirmed_flag": 0, "prior_positive_count": 0, "consecutive_worsening": 0, "transient_spike": 0,
                      "symptom_now": 0, "any_symptom_before": 0, "eoa_missing_now": b.ref_eoa.isna().astype(int), "dvi_missing_now": b.ref_dvi.isna().astype(int)})
    df = pd.concat([f, b.drop(columns=["valve_id"])], axis=1)
    return df


def select(df: pd.DataFrame, feature_set: list[str]) -> pd.DataFrame:
    return df[[c for c in feature_set if c in df.columns]]
