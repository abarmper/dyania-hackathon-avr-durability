"""Flow physics linking the latent orifice area to the Doppler quantities clinicians report.

MPG = c_model * (Q / EOA)^2            Q = stroke volume / ejection time (mL/s); Gorlin: c = 1/44.3^2
DVI = EOA / A_eff                      A_eff = effective LVOT area (cm2)
peak gradient = r * MPG; Vmax = sqrt(peak / 4)   (simplified Bernoulli)
c_model is fitted per valve model to the ASE 2024 normal-value rows (EOA, MPG) under the population flow
distribution, which absorbs pressure recovery (supra-annular valves have a lower c).
"""
from __future__ import annotations

import numpy as np

from valve_tables import EOA_TABLE


def flow_rate(svi, bsa, ejection_time):
    """mL/s"""
    return svi * bsa / ejection_time


Q_REF = 280.0
GAMMA = 0.25   # set from params by fit_physics_constants


def mpg_from_eoa(eoa, Q, c, q_ref=None, gamma=None):
    """MPG = c (Q / EOA_eff)^2 with EOA_eff = EOA (Q/q_ref)^gamma (flow dependence of bioprosthetic EOA)."""
    q_ref = Q_REF if q_ref is None else q_ref; gamma = GAMMA if gamma is None else gamma
    eoa_eff = np.maximum(eoa, 1e-3) * (Q / q_ref) ** gamma
    return c * (Q / eoa_eff) ** 2


def eoa_from_mpg(mpg, Q, c, q_ref=None, gamma=None):
    q_ref = Q_REF if q_ref is None else q_ref; gamma = GAMMA if gamma is None else gamma
    eoa_eff = Q * np.sqrt(c / np.maximum(mpg, 1e-6))
    return eoa_eff / (Q / q_ref) ** gamma


def dvi_from_eoa(eoa, a_eff):
    return eoa / a_eff


def vmax_from_mpg(mpg, r):
    return np.sqrt(np.maximum(r * mpg, 0.0) / 4.0)


def ppm_class(eoai, bmi):
    """0 none, 1 moderate, 2 severe (ASE 2024 Table 7 / VARC-3)."""
    eoai = np.asarray(eoai, float); bmi = np.asarray(bmi, float)
    obese = bmi >= 30
    mod = np.where(obese, 0.70, 0.85); sev = np.where(obese, 0.55, 0.65)
    return np.where(eoai <= sev, 2, np.where(eoai <= mod, 1, 0))


def population_Q(p, rng, size):
    """Draw reference flow rates the way baseline.draw_baseline does (BSA ~ N(2.0, 0.2) as a proxy)."""
    b = p.baseline
    svi = np.clip(rng.normal(b.svi[0], b.svi[1], size), b.svi[2], b.svi[3])
    bsa = np.clip(rng.normal(2.0, 0.2, size), 1.3, 2.8)
    et = np.clip(rng.normal(b.ejection_time_s[0], b.ejection_time_s[1], size), 0.22, 0.40)
    return np.clip(svi * bsa / et, *b.q_clip)


def fit_physics_constants(p, n_draws=20000, seed=12345):
    """Monte Carlo fit of c_model: for each model, c = median over sizes of table_MPG / mean[(Q/EOA_eff)^2]
    with Q and true EOA drawn as the simulator draws them. Fills p.physics.c_model; returns per-size
    residuals (implied population mean MPG minus table MPG, mmHg)."""
    global Q_REF, GAMMA
    Q_REF = p.physics.q_ref_ml_s; GAMMA = p.physics.eoa_flow_exponent
    from baseline import MEASUREMENT_CV_EOA
    rng = np.random.default_rng(seed)
    Q = population_Q(p, rng, n_draws)
    resid = {}
    for model, rows in EOA_TABLE.items():
        ratios, keep = {}, []
        for size, r in rows.items():
            sd_true = np.sqrt(max(r.eoa_sd ** 2 - (MEASUREMENT_CV_EOA * r.eoa) ** 2, 0.06 ** 2))
            sd_true = min(sd_true, 0.30)
            eoa = np.clip(rng.normal(r.eoa, sd_true, n_draws), 0.7, 3.2)
            ratios[size] = float(np.mean(mpg_from_eoa(eoa, Q, 1.0)))
            keep.append(r.mpg / ratios[size])
            p.physics.c_size[f"{model}|{size}"] = float(r.mpg / ratios[size])   # exact to the ASE row
        c = float(np.median(keep))
        p.physics.c_model[model] = c
        # residual if a single per-model constant were used (informational: per-size constants are exact)
        resid[model] = {size: round(c * ratios[size] - r.mpg, 2) for size, r in rows.items()}
    return resid
