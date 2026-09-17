"""Echo visit process: schedule, jitter, informative missingness, symptom-triggered echoes, and the
observation model (four independent error sources with the correct cross-channel structure)."""
from __future__ import annotations

import numpy as np

from latent import ar_star, dvi_star, eoa_star, mpg_star

VTYPE_REF, VTYPE_SCHED, VTYPE_SYMPTOM, VTYPE_OTHER = 0, 1, 2, 3


def miss_prob(year_index, p):
    v = p.visits
    k = np.asarray(year_index, float)
    base = np.array(v.miss_by_year)
    out = np.where(k <= len(base), base[np.minimum(k.astype(int), len(base)) - 1], np.nan)
    # logistic saturation beyond the observed years
    sat = v.miss_saturation; last = base[-1]
    extra = last + (sat - last) * (1 - np.exp(-0.5 * (k - len(base))))
    return np.where(np.isnan(out), extra, out)


def schedule_matrix(n, schedule, rng, p):
    sched = np.asarray(schedule, float)
    return sched[None, :] + rng.normal(0.0, p.visits.jitter_sd_y, (n, len(sched)))


def visit_times(df, t_end, t_sym, t_bg, p, rng, schedule=None, reference_time=None, stage_at=None, miss_scale=1.0):
    """Returns (t, vtype) arrays (n, S) sorted per row with NaN padding.
    stage_at(t) -> latent stage array, used for the 'feeling well' missingness multiplier.
    miss_scale scales the scheduled-visit missingness (trial cohorts < 1)."""
    v = p.visits; n = len(df)
    sched = list(range(1, v.max_years + 1)) if schedule is None else list(schedule)
    t_ref = v.reference_time_y if reference_time is None else reference_time
    ts = schedule_matrix(n, sched, rng, p)
    pm = miss_scale * np.broadcast_to(miss_prob(np.arange(1, len(sched) + 1), p)[None, :], ts.shape).copy()
    if stage_at is not None:
        st = stage_at(ts)
        symptomatic = ts >= t_sym[:, None]
        well = (st == 0) & ~symptomatic
        pm = pm * np.where(well, v.miss_mult_well, 1.0) * np.where(symptomatic, v.miss_mult_symptomatic, 1.0)
    keep = (rng.random(ts.shape) >= np.clip(pm, 0, 0.95)) & (ts < t_end[:, None]) & (ts > t_ref + 0.2)
    ts = np.where(keep, ts, np.nan)
    ref = np.where(df.reference_present.to_numpy() == 1, t_ref, np.nan)[:, None]
    sym = np.where(np.isfinite(t_sym) & (t_sym + v.symptom_echo_delay_y < t_end), t_sym + v.symptom_echo_delay_y, np.nan)[:, None]
    bg = np.where(np.isfinite(t_bg) & (t_bg > t_ref + 0.2), t_bg, np.nan)[:, None]
    t = np.concatenate([ref, ts, sym, bg], axis=1)
    vt = np.concatenate([np.full((n, 1), VTYPE_REF), np.full(ts.shape, VTYPE_SCHED), np.full((n, 1), VTYPE_SYMPTOM), np.full((n, 1), VTYPE_OTHER)], axis=1).astype(float)
    order = np.argsort(np.where(np.isnan(t), np.inf, t), axis=1)
    t = np.take_along_axis(t, order, 1); vt = np.take_along_axis(vt, order, 1)
    vt = np.where(np.isnan(t), np.nan, vt)
    return t, vt


NOISE_SOURCES = ("f", "f_ref", "eL", "eA", "eD", "vmax_eps", "lvef_eps", "ar_flip", "ar_sign", "pvl_flip", "pvl_sign", "eoa_miss", "dvi_miss", "lvef_miss")


def draw_noise_bank(n, K, p, rng):
    """Pre-drawn observation noise, one column per echo ordinal (common random numbers across surveillance policies).
    Same distributions as `observe` draws internally."""
    nz = p.noise
    return {"f": np.exp(rng.normal(0, nz.sigma_flow, (n, K))), "f_ref": np.exp(rng.normal(nz.flow_ref_bias, nz.sigma_flow_ref_extra, (n, K))),
            "eL": np.exp(rng.normal(0, nz.sigma_lvot_vti, (n, K))), "eA": np.exp(rng.normal(0, nz.sigma_cw, (n, K))), "eD": np.exp(rng.normal(0, nz.sigma_lvot_diam, (n, K))),
            "vmax_eps": np.exp(rng.normal(0, 0.03, (n, K))), "lvef_eps": rng.normal(0, 4.0, (n, K)),
            "ar_flip": rng.random((n, K)) < nz.ar_misclass, "ar_sign": rng.choice([-1, 1], (n, K)), "pvl_flip": rng.random((n, K)) < nz.ar_misclass, "pvl_sign": rng.choice([-1, 1], (n, K)),
            "eoa_miss": rng.random((n, K)) < nz.eoa_missing, "dvi_miss": rng.random((n, K)) < nz.dvi_missing, "lvef_miss": rng.random((n, K)) < nz.lvef_missing}


def observe(t, df, p, rng, vt=None, noise=None):
    """Observed channels at visit times t (n, S). Returns dict of (n, S) arrays incl. the true values.
    noise: optional dict of pre-drawn (n, S) arrays (see draw_noise_bank); otherwise drawn from rng in a fixed order."""
    nz = p.noise; n, S = t.shape
    e = eoa_star(t, df, p); m = mpg_star(e, df, t); d = dvi_star(e, df); a = ar_star(t, df)
    g = (lambda k, draw: noise[k]) if noise is not None else (lambda k, draw: draw())
    f = g("f", lambda: np.exp(rng.normal(0, nz.sigma_flow, (n, S))))
    if vt is not None:   # the reference (30-day) echo is taken in a more variable, slightly hyperdynamic post-operative state
        is_ref = vt == VTYPE_REF
        f = np.where(is_ref, f * g("f_ref", lambda: np.exp(rng.normal(nz.flow_ref_bias, nz.sigma_flow_ref_extra, (n, S)))), f)
    eL = g("eL", lambda: np.exp(rng.normal(0, nz.sigma_lvot_vti, (n, S))))
    eA = g("eA", lambda: np.exp(rng.normal(0, nz.sigma_cw, (n, S))))
    eD = g("eD", lambda: np.exp(rng.normal(0, nz.sigma_lvot_diam, (n, S))))
    bias = df.lvot_bias.to_numpy()[:, None]
    mpg_obs = m * f ** 2 * eA ** 2
    eoa_obs = e * eL * eD ** 2 * bias / eA
    dvi_obs = d * eL * bias / eA
    svi_obs = df.svi_true.to_numpy()[:, None] * f * eD ** 2 * bias
    vmax_obs = np.sqrt(np.maximum(df.peak_ratio.to_numpy()[:, None] * mpg_obs, 0) / 4.0) * g("vmax_eps", lambda: np.exp(rng.normal(0, 0.03, (n, S))))
    lvef_obs = np.clip(df.lvef.to_numpy()[:, None] + g("lvef_eps", lambda: rng.normal(0, 4.0, (n, S))), 15, 80)
    flip = g("ar_flip", lambda: rng.random((n, S)) < nz.ar_misclass)
    ar_obs = np.clip(a + np.where(flip, g("ar_sign", lambda: rng.choice([-1, 1], (n, S))), 0), 0, 3)
    pvl_true = df.pvl_ref_true.to_numpy()[:, None].astype(float)
    flip2 = g("pvl_flip", lambda: rng.random((n, S)) < nz.ar_misclass)
    pvl_obs = np.clip(pvl_true + np.where(flip2, g("pvl_sign", lambda: rng.choice([-1, 1], (n, S))), 0), 0, 3)
    # channel missingness given an echo
    eoa_obs = np.where(g("eoa_miss", lambda: rng.random((n, S)) < nz.eoa_missing), np.nan, eoa_obs)
    dvi_obs = np.where(g("dvi_miss", lambda: rng.random((n, S)) < nz.dvi_missing), np.nan, dvi_obs)
    lvef_obs = np.where(g("lvef_miss", lambda: rng.random((n, S)) < nz.lvef_missing), np.nan, lvef_obs)
    nan = np.isnan(t)
    out = {"mpg": np.round(mpg_obs, 1), "eoa": np.round(eoa_obs, 2), "dvi": np.round(dvi_obs, 2), "vmax": np.round(vmax_obs, 2),
           "svi": np.round(svi_obs, 1), "lvef": np.round(lvef_obs, 0), "ar": ar_obs, "pvl": pvl_obs,
           "mpg_true": m, "eoa_true": e, "dvi_true": d, "ar_true": a}
    for k in out:
        out[k] = np.where(nan, np.nan, out[k])
    return out
