"""Competing and terminal events: death, endocarditis, loss to follow-up, administrative censoring,
symptom onset (drives extra echoes and reintervention), and reintervention conditional on detection."""
from __future__ import annotations

import numpy as np


def draw_death(df, p, rng, profile=None):
    """Attained-age Gompertz with proportional-hazards multipliers; inverse-transform sampling.
    In the general cohort TAVR patients implanted <= tavr_era_until carry the higher-risk-era multiplier; in a
    paper-matched cohort the profile's own death_hr replaces it (NOTION: both arms 63-64% dead at 10 y)."""
    d = p.death; n = len(df)
    hr = d.hr_lvef_per_point ** (df.lvef.to_numpy() - 60.0)
    hr = hr * np.where(df.dialysis.to_numpy() == 1, d.hr_dialysis, 1.0)
    if profile is None:
        hr = hr * np.where((df.approach.to_numpy() == "TAVR") & (df.implant_year.to_numpy() <= d.tavr_era_until), d.hr_tavr_era, 1.0)
    else:
        hr = hr * float(getattr(profile, "death_hr", 1.0))
    if getattr(d, "frailty_var", 0.0) > 0:
        hr = hr * rng.gamma(1.0 / d.frailty_var, d.frailty_var, n)      # gamma frailty, mean 1
    A = d.a * np.exp(d.b_per_year * (df.age_at_implant.to_numpy() - d.age_ref)) * hr
    E = rng.exponential(1.0, n)
    return np.log1p(d.b_per_year * E / A) / d.b_per_year


def draw_endocarditis(df, p, rng):
    e = p.events; n = len(df)
    t = rng.exponential(1.0 / e.endocarditis_rate, n)
    u = rng.random(n)
    outcome = np.where(u < e.endocarditis_explant_frac, "explant",
                       np.where(u < e.endocarditis_explant_frac + e.endocarditis_death_frac, "death", "treated")).astype(object)
    t_out = np.where(outcome == "explant", t + rng.uniform(0.05, 0.5, n), np.where(outcome == "death", t + rng.uniform(0.02, 0.3, n), np.inf))
    return t, outcome, t_out


def draw_ltfu(df, p, rng):
    return rng.exponential(1.0 / p.events.ltfu_rate, len(df))


def admin_censoring(df, p):
    return np.maximum(p.cohort.censor_year - df.implant_year.to_numpy(), 0.1)


def draw_symptoms(grid, stage_grid, mpg_grid, t_end, p, rng):
    """First valve-symptom time from h(t) = h0 (MPG*/20)^power when latent stage >= 2; inf if none before t_end."""
    v = p.visits; n = stage_grid.shape[0]; dt = grid[1] - grid[0]
    h = np.where(stage_grid >= 2, v.symptom_h0 * (np.maximum(mpg_grid, 1.0) / 20.0) ** v.symptom_power, 0.0)
    H = np.cumsum(h * dt, axis=1)
    E = rng.exponential(1.0, n)[:, None]
    hit = H >= E
    t_sym = np.where(hit.any(1), grid[hit.argmax(1)], np.inf)
    t_bg = rng.exponential(1.0 / v.background_echo_rate, n)
    return np.where(t_sym < t_end, t_sym, np.inf), np.where(t_bg < t_end, t_bg, np.inf)


def draw_reintervention(df, t_detect3, t_detect2c, t_sym, t_detect3_symptomatic, t_end, p, rng, draws=None):
    """Reintervention hazard from the first detection of stage 3 (single echo), or from a confirmed stage 2
    with symptoms; type ViV vs redo by age at reintervention; urgent if stage 3 was first seen at a
    symptom-triggered echo (or symptoms preceded the detection)."""
    e = p.events; n = len(df)
    agef = np.exp(e.reint_age_slope * (df.age_at_implant.to_numpy() - 70.0))
    h3 = np.clip(e.reint_h_stage3 * agef, 0.05, 3.0); h2 = np.clip(e.reint_h_stage2_symptomatic * agef, 0.005, 0.5)
    E2, E3, u_viv = draws if draws is not None else (rng.exponential(1.0, n), rng.exponential(1.0, n), None)
    start2 = np.maximum(t_detect2c, t_sym)
    t2 = np.where(np.isfinite(start2), start2 + E2 / h2, np.inf)
    t3 = np.where(np.isfinite(t_detect3), t_detect3 + E3 / h3, np.inf)
    t_re = np.minimum(t2, t3)
    happens = t_re < t_end
    t_re = np.where(happens, t_re, np.inf)
    age_re = df.age_at_implant.to_numpy() + np.where(happens, t_re, 0.0)
    yr_re = df.implant_year.to_numpy() + np.where(happens, t_re, 0.0)
    p_viv = 1.0 / (1.0 + np.exp(-(age_re - e.viv_logistic[0]) * e.viv_logistic[1]))
    u = rng.random(n) if u_viv is None else u_viv
    rtype = np.where(happens & (yr_re >= e.viv_min_year) & (u < p_viv), "ViV", "redo").astype(object)
    urgent = happens & ((t3 <= t2) & (t_detect3_symptomatic | (t_sym <= t_detect3)))
    indication = np.where(happens, np.where(t3 <= t2, "severe SVD", "moderate SVD with symptoms"), "").astype(object)
    return t_re, rtype, urgent, indication
