"""Latent valve process: onset hazard (proportional hazards on a per-class Weibull), the true orifice-area
trajectory, the regurgitant failure mode, transient thrombosis episodes, and the noise-free VARC-3 stage."""
from __future__ import annotations

import numpy as np
import pandas as pd

from physics import dvi_from_eoa, mpg_from_eoa


def onset_hr(df: pd.DataFrame, p) -> np.ndarray:
    """Multiplicative hazard ratio of latent onset for every valve (biological effects only)."""
    o = p.onset
    is_savr = (df.approach == "SAVR").to_numpy()
    age_ref = df.valve_class.map(o.calibration_age).to_numpy(float)
    hr = np.where(is_savr, o.hr_age_per_year_savr, o.hr_age_per_year_thv) ** (df.age_at_implant.to_numpy() - age_ref)
    hr = hr * o.hr_smoking ** df.smoking.to_numpy() * o.hr_diabetes ** df.diabetes.to_numpy()
    hr = hr * o.hr_bmi_per_unit ** (df.bmi.to_numpy() - o.bmi_centre) * o.hr_ckd ** df.ckd.to_numpy()
    hr = hr * np.where((~is_savr) & (df.anticoag_year1.to_numpy() == 0), o.hr_no_anticoag_year1_tavr, 1.0)
    hr = hr * np.where(is_savr & (df.anticoag_type.to_numpy() == "VKA") & (df.anticoag_longterm.to_numpy() == 1), o.hr_vka_longterm_savr, 1.0)
    hr = hr * np.where(df.ppm_true.to_numpy() == 2, o.hr_severe_ppm_residual, 1.0)
    hr = hr * np.where(df.mpg_ref_true.to_numpy() >= 15.0, getattr(o, "hr_ref_mpg_ge15_residual", 1.0), 1.0)
    hr = hr * np.where((~is_savr) & (df.label_size_mm.to_numpy() <= 23), o.hr_thv_le23_residual, 1.0)
    if "thromb_on" in df:
        hr = hr * np.where(df.thromb_on.to_numpy() == 1, o.hr_prior_thrombosis, 1.0)
    hr = hr * np.where(df.sex.to_numpy() == "F", o.hr_female, 1.0)
    return hr


def draw_thrombosis(df, p, rng):
    """Transient thrombosis / HALT episode in year 1 (switchable)."""
    n = len(df); t = p.thrombosis
    if not t.on:
        df["thromb_on"] = 0; df["thromb_start"] = np.inf; df["thromb_end"] = np.inf; return df
    pr = np.where(df.approach == "SAVR", t.p_year1["SAVR"], t.p_year1["TAVR"]) * np.where(df.anticoag_year1 == 1, t.anticoag_factor, 1.0)
    on = rng.random(n) < pr
    start = rng.uniform(*t.onset_window_y, n); dur = rng.uniform(*t.duration_y, n)
    df["thromb_on"] = on.astype(int)
    df["thromb_start"] = np.where(on, start, np.inf); df["thromb_end"] = np.where(on, start + dur, np.inf)
    return df


def draw_onset(df, p, rng, theta_by_class=None):
    """T_deg = shift + scale * exp(log_mult) * (-ln U / HR)^(1/shape); mode; progression time constant tau."""
    n = len(df); o = p.onset; tr = p.trajectory
    k = df.valve_class.map(o.shape).to_numpy(float); lam = df.valve_class.map(o.scale).to_numpy(float)
    shift = df.valve_class.map(getattr(o, "shift", {})).fillna(0.0).to_numpy(float) if hasattr(o, "shift") else np.zeros(n)
    if theta_by_class:
        for cls, th in theta_by_class.items():
            sh, lm = th[0], th[1]; lk = th[2] if len(th) > 2 else 0.0
            m = (df.valve_class == cls).to_numpy(); lam[m] *= np.exp(lm); k[m] *= np.exp(lk); shift[m] = shift[m] + sh
    hr = onset_hr(df, p)
    U = rng.random(n)
    T = shift + lam * (-np.log(U) / hr) ** (1.0 / k)
    T = np.maximum(T, 0.05)
    p_reg = df.valve_class.map(o.p_regurgitant).to_numpy(float)
    mode = np.where(rng.random(n) < p_reg, "regurgitant", "stenotic").astype(object)
    tau_med = df.valve_class.map(tr.tau_median_y).to_numpy(float) if isinstance(tr.tau_median_y, dict) else np.full(n, float(tr.tau_median_y))
    tau = np.exp(rng.normal(np.log(tau_med), tr.tau_sigma, n))
    lf = df.valve_class.map(lambda c: tr.loss_fraction[c][0]).to_numpy(float); lsd = df.valve_class.map(lambda c: tr.loss_fraction[c][1]).to_numpy(float)
    loss = np.clip(rng.normal(lf, lsd), 0.05, 0.98)
    ar_sev = T + np.exp(rng.normal(np.log(tr.ar_severe_after_median_y), tr.ar_severe_after_sigma, n))
    df["onset_hr"] = hr; df["T_deg"] = T; df["mode"] = mode; df["tau"] = tau; df["loss_fraction"] = loss; df["ar_severe_time"] = ar_sev
    return df


def eoa_star(t, df, p):
    """True EOA at times t (n, m). NaN in t propagates."""
    tr = p.trajectory; th = p.thrombosis
    eoa0 = df.eoa_ref_true.to_numpy()[:, None]
    T = df.T_deg.to_numpy()[:, None]; tau = df.tau.to_numpy()[:, None]
    sten = (df["mode"].to_numpy() == "stenotic")[:, None]
    u = np.maximum(t - T, 0.0)
    L = df.loss_fraction.to_numpy()[:, None] if "loss_fraction" in df else 1.0
    D = np.where(sten, 1.0 - L * (1.0 - np.exp(-u / tau)), np.where(u > 0, tr.regurgitant_eoa_factor, 1.0))
    base = eoa0 * np.maximum(1.0 - tr.drift_per_year * t, 0.3)
    if th.on and "thromb_on" in df:
        s = df.thromb_start.to_numpy()[:, None]; e = df.thromb_end.to_numpy()[:, None]
        base = base * np.where((t >= s) & (t < e), th.eoa_factor, 1.0)
    return np.maximum(base * D, p.physics.eoa_floor)


def ar_star(t, df):
    ar0 = df.ar_ref_true.to_numpy()[:, None]
    reg = (df["mode"].to_numpy() == "regurgitant")[:, None]
    T = df.T_deg.to_numpy()[:, None]; sev = df.ar_severe_time.to_numpy()[:, None]
    g = np.where(reg & (t >= T), np.maximum(ar0, 2), ar0)
    g = np.where(reg & (t >= sev), 3, g)
    return np.where(np.isnan(t), np.nan, g.astype(float))


def flow_at(t, df):
    """Per-valve flow state drifts slowly: Q(t) = Q_ref exp(slope t) (NaN in t propagates)."""
    return df.Q_ref.to_numpy()[:, None] * np.exp(df.flow_slope.to_numpy()[:, None] * np.nan_to_num(t, nan=0.0))


def mpg_star(eoa, df, t=None):
    Q = df.Q_ref.to_numpy()[:, None] if t is None else flow_at(t, df)
    return mpg_from_eoa(eoa, Q, df.c_model.to_numpy()[:, None])


def dvi_star(eoa, df):
    return dvi_from_eoa(eoa, df.a_eff.to_numpy()[:, None])


def latent_stage(t, df, p):
    """Noise-free VARC-3 stage (0/2/3) at times t (n, m), against the true values at the reference time."""
    from labels import varc3_stage
    t_ref = np.full((len(df), 1), p.visits.reference_time_y)
    e_ref = eoa_star(t_ref, df, p); m_ref = mpg_star(e_ref, df, t_ref); d_ref = dvi_star(e_ref, df); a_ref = ar_star(t_ref, df)
    e = eoa_star(t, df, p); m = mpg_star(e, df, t); d = dvi_star(e, df); a = ar_star(t, df)
    st, _ = varc3_stage(m, e, d, a, m_ref, e_ref, d_ref, a_ref)
    return st, (m, e, d, a)


def latent_crossing_times(df, p, t_max=20.0, dt=1 / 12):
    """Continuous-time (monthly) first crossing of stage 2 and stage 3 on the noise-free trajectory."""
    n = len(df); grid = np.arange(dt, t_max + 1e-9, dt)
    t = np.broadcast_to(grid[None, :], (n, len(grid)))
    st, _ = latent_stage(t, df, p)
    def first(mask):
        any_ = mask.any(1); idx = mask.argmax(1)
        return np.where(any_, grid[idx], np.inf)
    return first(st >= 2), first(st >= 3), grid, st
