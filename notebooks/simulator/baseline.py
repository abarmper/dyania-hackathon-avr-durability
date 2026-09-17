"""Draw the patients_valves table (one row per implanted valve) and the latent reference haemodynamics."""
from __future__ import annotations

import numpy as np
import pandas as pd

from physics import dvi_from_eoa, flow_rate, mpg_from_eoa, ppm_class
from valve_tables import (A_EFF, A_EFF_SAVR_DIAM_CM, EOA_TABLE, MODEL_CLASS, MODEL_PROBS_SAVR, MODEL_PROBS_TAVR,
                          SIZE_PROBS, CLASS_PARENT)

MEASUREMENT_CV_EOA = 0.15   # the ASE rows are measured values (single echo, incl. flow variation); true between-valve SD is smaller


def trunc_normal(rng, mean, sd, lo, hi, size):
    x = rng.normal(mean, sd, size)
    bad = (x < lo) | (x > hi)
    while bad.any():
        x[bad] = rng.normal(mean, sd, int(bad.sum()))
        bad = (x < lo) | (x > hi)
    return x


def choice(rng, probs: dict, size):
    keys = list(probs)
    pr = np.array([probs[k] for k in keys], float); pr = pr / pr.sum()
    return np.array(keys, dtype=object)[rng.choice(len(keys), size=size, p=pr)]


def _size_probs_for(model, approach):
    sizes = sorted(EOA_TABLE[model].keys())
    if approach == "SAVR":
        base = {s: SIZE_PROBS["SAVR"].get(s, 0.02) for s in sizes}
    elif MODEL_CLASS[model] == "balloon":
        base = {s: SIZE_PROBS["TAVR_balloon"].get(s, 0.03) for s in sizes}
    else:
        base = {s: SIZE_PROBS["TAVR_self_expanding"].get(s, 0.05) for s in sizes}
    return base


def draw_sizes(rng, models, approach, bsa, profile=None):
    """Labelled size per valve; larger BSA tilts towards larger sizes (weights x exp(1.5*(bsa-2)*rank))."""
    n = len(models); out = np.zeros(n, int)
    for m in np.unique(models):
        idx = np.where(models == m)[0]
        probs = dict(profile.sizes) if profile is not None else _size_probs_for(m, approach[idx[0]])
        probs = {s: p for s, p in probs.items() if s in EOA_TABLE[m]} or _size_probs_for(m, approach[idx[0]])
        sizes = np.array(sorted(probs)); base = np.array([probs[s] for s in sizes], float)
        rank = np.arange(len(sizes)) - (len(sizes) - 1) / 2
        w = base[None, :] * np.exp(1.5 * (bsa[idx, None] - 2.0) * rank[None, :])
        w = w / w.sum(1, keepdims=True)
        cum = np.cumsum(w, 1); u = rng.random(len(idx))[:, None]
        out[idx] = sizes[np.minimum((u > cum).sum(1), len(sizes) - 1)]
    return out


def draw_baseline(n, p, rng, profile=None) -> pd.DataFrame:
    b = p.baseline
    if profile is None:
        approach = np.where(rng.random(n) < p.cohort.savr_frac, "SAVR", "TAVR").astype(object)
        implant_year = rng.uniform(p.cohort.implant_year_min, p.cohort.implant_year_max, n)
    else:
        approach = np.full(n, profile.approach, dtype=object)
        implant_year = rng.uniform(2008.0, 2016.0, n)
    is_savr = approach == "SAVR"

    model = np.empty(n, dtype=object)
    if profile is None:
        model[is_savr] = choice(rng, MODEL_PROBS_SAVR, int(is_savr.sum()))
        grp = choice(rng, MODEL_PROBS_TAVR, int((~is_savr).sum())); yrs = implant_year[~is_savr]
        model[~is_savr] = np.where(grp == "balloon", np.where(yrs < 2016.0, "Sapien XT", "Sapien 3"),
                                   np.where(grp == "self_expanding", np.where(yrs < 2015.0, "CoreValve", "Evolut R"), "New THV"))
    else:
        model[:] = choice(rng, profile.models, n)
    vclass = np.array([MODEL_CLASS[m] for m in model], dtype=object)
    is_new = np.isin(model, ["Inspiris Resilia", "New SAVR valve", "New THV"])

    if profile is None:
        age = np.where(is_savr, trunc_normal(rng, *b.age_savr, n), trunc_normal(rng, *b.age_tavr, n))
    else:
        age = trunc_normal(rng, profile.age_mean, profile.age_sd, 40.0, 95.0, n)
    female_frac = b.female_frac if profile is None else profile.female_frac
    sex = np.where(rng.random(n) < female_frac, "F", "M").astype(object)
    height = np.where(sex == "M", rng.normal(*b.height_cm["M"], n), rng.normal(*b.height_cm["F"], n))
    bmi = trunc_normal(rng, *b.bmi, n)
    if profile is not None and profile.name == "wakami":       # Japanese cohort: BSA 1.54 +/- 0.18
        height = height - 12.0; bmi = trunc_normal(rng, 23.0, 3.5, 16.0, 40.0, n)
    weight = bmi * (height / 100.0) ** 2
    bsa = np.sqrt(height * weight / 3600.0)

    diabetes = (rng.random(n) < (b.diabetes if profile is None else profile.diabetes)).astype(int)
    ckd = (rng.random(n) < (b.ckd if profile is None else profile.ckd)).astype(int)
    dialysis = (ckd == 1) & (rng.random(n) < b.dialysis_given_ckd)
    smoking = (rng.random(n) < (b.smoking if profile is None else profile.smoking)).astype(int)
    af_p = np.where(is_savr, b.af["SAVR"], b.af["TAVR"]) if profile is None else profile.af
    af = (rng.random(n) < af_p).astype(int)
    ac = np.where(af == 1, rng.random(n) < b.anticoag_if_af, rng.random(n) < b.anticoag_if_no_af)
    ac_type = np.where(ac, np.where(implant_year < b.doac_from_year, "VKA", "DOAC"), "none").astype(object)
    postop = is_savr & ~ac & (rng.random(n) < b.postop_vka_savr)
    anticoag_year1 = (ac | postop).astype(int)
    anticoag_longterm = ac.astype(int)

    lvef = trunc_normal(rng, *b.lvef, n)
    svi = trunc_normal(rng, b.svi[0], b.svi[1], b.svi[2], b.svi[3], n) + b.svi_per_lvef_point * (lvef - 59.0)
    svi = np.clip(svi, 22.0, 70.0)
    et = trunc_normal(rng, b.ejection_time_s[0], b.ejection_time_s[1], 0.22, 0.40, n)
    Q = np.clip(flow_rate(svi, bsa, et), *b.q_clip)

    size = draw_sizes(rng, model, approach, bsa, profile)
    rows = [EOA_TABLE[m][s] for m, s in zip(model, size)]
    eoa_mean = np.array([r.eoa for r in rows]); eoa_sd = np.array([r.eoa_sd for r in rows])
    sd_true = np.sqrt(np.maximum(eoa_sd ** 2 - (MEASUREMENT_CV_EOA * eoa_mean) ** 2, 0.06 ** 2))
    sd_true = np.minimum(sd_true, 0.30)                      # a few ASE rows have SDs > 0.5 (small n); cap the true spread
    eoa_ref = trunc_normal(rng, 0.0, 1.0, -2.5, 2.5, n) * sd_true + eoa_mean
    eoa_ref = np.clip(eoa_ref, 0.7, 3.2)
    a_eff = np.where(is_savr, np.pi * (rng.normal(*A_EFF_SAVR_DIAM_CM, n) / 2.0) ** 2,
                     np.array([A_EFF.get(m, 3.4) for m in model]) * np.exp(rng.normal(0, 0.08, n)))
    lvot_bias = np.exp(rng.normal(0.0, b.lvot_bias_sd, n))
    c = np.array([p.physics.c_size.get(f"{m}|{s}", p.physics.c_model[m]) for m, s in zip(model, size)])
    peak_ratio = rng.normal(*p.physics.peak_over_mean, n)
    flow_slope = rng.normal(0.0, p.trajectory.flow_slope_sd_per_year, n)
    mpg_ref = mpg_from_eoa(eoa_ref, Q, c)
    dvi_ref = dvi_from_eoa(eoa_ref, a_eff)
    eoai = eoa_ref / bsa
    ppm = ppm_class(eoai, bmi)

    ar_ref = choice(rng, {0: b.ar_ref_probs[0], 1: b.ar_ref_probs[1], 2: b.ar_ref_probs[2]}, n).astype(int)
    early_thv = np.isin(model, ["Sapien XT", "CoreValve"]) | ((~is_savr) & (implant_year < 2016.0))
    pvl = np.zeros(n, int)
    for mask, key in ((is_savr, "SAVR"), ((~is_savr) & early_thv, "TAVR_early"), ((~is_savr) & ~early_thv, "TAVR_late")):
        k = int(mask.sum())
        if k:
            pr = b.pvl_ref_probs[key]; pvl[mask] = choice(rng, {0: pr[0], 1: pr[1], 2: pr[2]}, k).astype(int)
    ref_present = (rng.random(n) >= b.reference_missing_frac) if profile is None else np.ones(n, bool)

    df = pd.DataFrame({
        "valve_id": np.arange(1, n + 1), "approach": approach, "valve_model": model, "valve_class": vclass,
        "is_new_model": is_new.astype(int), "implant_year": implant_year, "age_at_implant": age, "sex": sex,
        "height_cm": height, "weight_kg": weight, "bsa": bsa, "bmi": bmi, "diabetes": diabetes, "ckd": ckd,
        "dialysis": dialysis.astype(int), "smoking": smoking, "af": af, "anticoag_discharge": ac.astype(int),
        "anticoag_type": ac_type, "anticoag_year1": anticoag_year1, "anticoag_longterm": anticoag_longterm,
        "lvef": lvef, "svi_true": svi, "ejection_time_s": et, "Q_ref": Q, "label_size_mm": size,
        "eoa_ref_true": eoa_ref, "a_eff": a_eff, "lvot_bias": lvot_bias, "c_model": c, "peak_ratio": peak_ratio, "flow_slope": flow_slope,
        "mpg_ref_true": mpg_ref, "dvi_ref_true": dvi_ref, "eoai_true": eoai, "ppm_true": ppm,
        "ar_ref_true": ar_ref, "pvl_ref_true": pvl, "reference_present": ref_present.astype(int),
        "table_eoa": eoa_mean, "table_mpg": np.array([r.mpg for r in rows]),
    })
    return df
