"""Orchestrator: baseline -> latent -> events -> visits -> observation -> labels -> reintervention -> tables.

simulate_cohort(p, n, seed)             -> dict of DataFrames (patients_valves, echo_visits, events, labels, latent_truth)
simulate_matched(paper, theta, n, ...)  -> (times, states) for km_reconstruct/calibrate_onset.calibrate
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from baseline import draw_baseline                                   # noqa: E402
from events import admin_censoring, draw_death, draw_endocarditis, draw_ltfu, draw_reintervention, draw_symptoms  # noqa: E402
from labels import capodanno, dvir, first_times, varc3_stage         # noqa: E402
from latent import draw_onset, draw_thrombosis, latent_crossing_times, latent_stage  # noqa: E402
from physics import fit_physics_constants                            # noqa: E402
from valve_tables import COHORT_PROFILES, MISS_SCALE                 # noqa: E402
from visits import VTYPE_SYMPTOM, observe, visit_times               # noqa: E402

MEASUREMENTS = ["mpg", "vmax", "eoa", "dvi", "svi", "lvef", "ar", "pvl"]
MEASUREMENT_NAMES = {"mpg": "av_mean_gradient_mmHg", "vmax": "av_peak_velocity_m_s", "eoa": "av_area_cm2", "dvi": "doppler_velocity_index",
                     "svi": "stroke_volume_index_ml_m2", "lvef": "lvef_pct", "ar": "intraprosthetic_ar_grade", "pvl": "paravalvular_leak_grade"}


def pre_visits(p, n, rng, profile=None, theta_by_class=None, t_max=None):
    """Every draw that precedes the echo visits (baseline, thrombosis, onset, death, endocarditis, LTFU, symptoms)."""
    if not p.physics.c_model:
        fit_physics_constants(p)
    df = draw_baseline(n, p, rng, profile)
    df = draw_thrombosis(df, p, rng)
    df = draw_onset(df, p, rng, theta_by_class)
    t_death = draw_death(df, p, rng, profile=profile)
    t_endo, endo_out, t_endo_out = draw_endocarditis(df, p, rng)
    t_death = np.where((endo_out == "death") & (t_endo_out < t_death), t_endo_out, t_death)
    t_ltfu = draw_ltfu(df, p, rng)
    t_admin = admin_censoring(df, p) if t_max is None else np.full(n, float(t_max))
    t_endo_explant = np.where(endo_out == "explant", t_endo_out, np.inf)
    t_end0 = np.minimum.reduce([t_death, t_ltfu, t_admin, t_endo_explant])
    t2_lat, t3_lat, grid, st_grid = latent_crossing_times(df, p, t_max=float(p.visits.max_years) + 1.0)
    _, (mpg_grid, _, _, _) = latent_stage(np.broadcast_to(grid[None, :], (n, len(grid))), df, p)
    t_sym, t_bg = draw_symptoms(grid, st_grid, mpg_grid, t_end0, p, rng)
    return {"df": df, "t_death": t_death, "t_endo": t_endo, "endo_outcome": endo_out, "t_endo_out": t_endo_out, "t_ltfu": t_ltfu,
            "t_admin": t_admin, "t_end0": t_end0, "t2_latent": t2_lat, "t3_latent": t3_lat, "t_sym": t_sym, "t_bg": t_bg, "profile": profile}


def post_visits(pre, t, vt, obs, p, rng, reference_time=None, reint_draws=None):
    """Reference values, labels, reintervention (conditional on detection) and truncation, given the visit matrices."""
    df = pre["df"]; n = len(df); t_end0 = pre["t_end0"]; t_endo = pre["t_endo"]; t_sym = pre["t_sym"]
    ref_present = df.reference_present.to_numpy() == 1
    is_ref = (vt == 0)
    def at_ref(arr):
        v = np.where(is_ref, arr, np.nan); return np.nanmax(np.where(np.isnan(v), -np.inf, v), axis=1)
    with np.errstate(invalid="ignore"):
        ref_mpg = np.where(ref_present, at_ref(obs["mpg"]), df.table_mpg.to_numpy())
        ref_eoa = np.where(ref_present, at_ref(obs["eoa"]), df.table_eoa.to_numpy())
        ref_dvi = np.where(ref_present, at_ref(obs["dvi"]), df.table_eoa.to_numpy() / df.a_eff.to_numpy())
        ref_ar = np.where(ref_present, at_ref(obs["ar"]), 0.0)
    ref_eoa = np.where(np.isfinite(ref_eoa), ref_eoa, df.table_eoa.to_numpy())
    ref_dvi = np.where(np.isfinite(ref_dvi), ref_dvi, df.table_eoa.to_numpy() / df.a_eff.to_numpy())
    ref_mpg = np.where(np.isfinite(ref_mpg), ref_mpg, df.table_mpg.to_numpy())
    ref = {"mpg": ref_mpg, "eoa": ref_eoa, "dvi": ref_dvi, "ar": ref_ar}
    post = ~is_ref
    R = lambda k: ref[k][:, None]
    def masked(arr): return np.where(post, arr, np.nan)
    st_full, st_h = varc3_stage(masked(obs["mpg"]), masked(obs["eoa"]), masked(obs["dvi"]), masked(obs["ar"]), R("mpg"), R("eoa"), R("dvi"), R("ar"))
    st_cap = capodanno(masked(obs["mpg"]), R("mpg"))
    st_dvir = dvir(masked(obs["mpg"]), masked(obs["eoa"]), masked(obs["dvi"]), R("mpg"), R("eoa"), R("dvi"))
    st_or, _ = latent_stage(np.where(post, t, np.nan), df, p)
    stages = {"varc3": st_full, "varc3_haemo": st_h, "capodanno": st_cap, "dvir": st_dvir, "oracle_varc3": st_or}
    endo_mask = (t >= t_endo[:, None])
    for k in stages:
        stages[k] = np.where(endo_mask, np.nan, stages[k])
    times = first_times(stages, np.where(np.isnan(t), np.inf, t))
    t_det3 = times["varc3"]["single3"]
    det3_idx = np.nan_to_num(stages["varc3"] >= 3, nan=False).argmax(1)
    det3_sym = np.isfinite(t_det3) & (vt[np.arange(n), det3_idx] == VTYPE_SYMPTOM)
    t_re, rtype, urgent, indication = draw_reintervention(df, t_det3, times["varc3"]["conf2"], t_sym, det3_sym, t_end0, p, rng, draws=reint_draws)
    t_end = np.minimum(t_end0, t_re)
    after = t > t_end[:, None]
    t = np.where(after, np.nan, t); vt = np.where(after, np.nan, vt)
    for k in obs: obs[k] = np.where(after, np.nan, obs[k])
    for k in stages: stages[k] = np.where(after, np.nan, stages[k])
    times = first_times(stages, np.where(np.isnan(t), np.inf, t))
    ends = {"t_death": pre["t_death"], "t_endo": t_endo, "endo_outcome": pre["endo_outcome"], "t_endo_out": pre["t_endo_out"], "t_ltfu": pre["t_ltfu"],
            "t_admin": pre["t_admin"], "t_sym": t_sym, "t_bg": pre["t_bg"], "t_re": t_re, "re_type": rtype, "re_urgent": urgent, "re_indication": indication, "t_end": t_end}
    latent = {"t2_latent": pre["t2_latent"], "t3_latent": pre["t3_latent"]}
    return df, t, vt, obs, ref, stages, times, t_end, ends, latent


def _core(p, n, rng, profile=None, theta_by_class=None, schedule=None, reference_time=None, t_max=None):
    """Everything up to labels. Returns (df, t, vt, obs, ref, stages, times, t_end, ends, latent)."""
    pre = pre_visits(p, n, rng, profile, theta_by_class, t_max)
    df = pre["df"]
    stage_at = lambda ts: latent_stage(ts, df, p)[0]
    miss_scale = MISS_SCALE.get(profile.name, 1.0) if profile is not None else 1.0
    t, vt = visit_times(df, pre["t_end0"], pre["t_sym"], pre["t_bg"], p, rng, schedule=schedule, reference_time=reference_time, stage_at=stage_at, miss_scale=miss_scale)
    obs = observe(t, df, p, rng, vt=vt)
    return post_visits(pre, t, vt, obs, p, rng, reference_time=reference_time)


def simulate_cohort(p, n=None, seed=None, policy="guideline"):
    n = n or p.cohort.n; seed = p.cohort.seed if seed is None else seed
    rng = np.random.default_rng(seed)
    df, t, vt, obs, ref, stages, times, t_end, ends, latent = _core(p, n, rng)
    vid = df.valve_id.to_numpy()
    # ---- patients_valves
    pv = df.drop(columns=["eoa_ref_true", "mpg_ref_true", "dvi_ref_true", "ppm_true", "ar_ref_true", "pvl_ref_true", "svi_true",
                          "Q_ref", "lvot_bias", "c_model", "peak_ratio", "onset_hr", "T_deg", "mode", "tau", "ar_severe_time",
                          "thromb_on", "thromb_start", "thromb_end", "eoai_true", "a_eff", "ejection_time_s", "loss_fraction", "flow_slope"]).copy()
    pv["reference_source"] = np.where(df.reference_present == 1, "30d_echo", "model_size_table")
    pv["ref_mean_gradient_mmHg"] = np.round(ref["mpg"], 1); pv["ref_av_area_cm2"] = np.round(ref["eoa"], 2)
    pv["ref_dvi"] = np.round(ref["dvi"], 2); pv["ref_ar_grade"] = ref["ar"]
    pv["ref_eoai_cm2_m2"] = np.round(ref["eoa"] / df.bsa, 2)
    from physics import ppm_class
    pv["ppm_class_ref"] = ppm_class(pv["ref_eoai_cm2_m2"], df.bmi)
    pv["implant_date"] = pd.to_datetime((df.implant_year - 1970) * 365.25, unit="D").dt.date
    pv["followup_end_years"] = np.round(t_end, 3)
    # ---- echo_visits (long)
    n_, S = t.shape
    rows = []
    for k in MEASUREMENTS:
        m = ~np.isnan(t) & ~np.isnan(obs[k])
        rows.append(pd.DataFrame({"valve_id": np.repeat(vid, S)[m.ravel()], "visit_time_years": np.round(t, 3).ravel()[m.ravel()],
                                  "visit_type": vt.ravel()[m.ravel()], "measurement": MEASUREMENT_NAMES[k], "result": obs[k].ravel()[m.ravel()]}))
    ev = pd.concat(rows, ignore_index=True)
    ev["visit_type"] = ev.visit_type.map({0: "reference", 1: "scheduled", 2: "symptom_triggered", 3: "other_indication"})
    ev["visit_date"] = pd.to_datetime((df.implant_year.to_numpy()[ev.valve_id - 1] + ev.visit_time_years - 1970) * 365.25, unit="D").dt.date
    ev = ev.sort_values(["valve_id", "visit_time_years", "measurement"]).reset_index(drop=True)
    # ---- events
    evs = []
    def add(mask, etype, tt, detail):
        for i in np.where(mask)[0]:
            evs.append((vid[i], etype, round(float(tt[i]), 3), detail[i] if not isinstance(detail, str) else detail))
    add(ends["t_death"] <= t_end + 1e-9, "death", ends["t_death"], np.where(ends["endo_outcome"] == "death", "endocarditis", "other").astype(object))
    add(ends["t_endo"] < t_end, "endocarditis", ends["t_endo"], ends["endo_outcome"])
    add(np.isfinite(ends["t_re"]), "reintervention", ends["t_re"],
        np.array([f"{a}|{b}|{'urgent' if u else 'elective'}" for a, b, u in zip(ends["re_type"], ends["re_indication"], ends["re_urgent"])], dtype=object))
    add((ends["endo_outcome"] == "explant") & np.isclose(t_end, ends["t_endo_out"]), "reintervention", t_end, "redo|endocarditis|urgent")
    add((ends["t_ltfu"] <= t_end + 1e-9) & (ends["t_ltfu"] < ends["t_death"]) & ~np.isfinite(ends["t_re"]), "lost_to_follow_up", ends["t_ltfu"], "")
    add((ends["t_admin"] <= t_end + 1e-9) & (ends["t_admin"] < ends["t_death"]) & (ends["t_admin"] < ends["t_ltfu"]) & ~np.isfinite(ends["t_re"]), "administrative_censoring", ends["t_admin"], "2026-09-15")
    events = pd.DataFrame(evs, columns=["valve_id", "event_type", "time_years", "detail"]).sort_values(["valve_id", "time_years"]).reset_index(drop=True)
    # ---- labels (per valve)
    lab = pd.DataFrame({"valve_id": vid, "followup_end_years": t_end,
                        "death": (ends["t_death"] <= t_end + 1e-9).astype(int), "reintervention": np.isfinite(ends["t_re"]).astype(int),
                        "endocarditis": (ends["t_endo"] < t_end).astype(int)})
    for name, d in times.items():
        for k, v in d.items():
            lab[f"{name}_{k}"] = np.where(np.isfinite(v), v, np.nan)
    lab["primary_event"] = np.isfinite(times["varc3"]["conf2"]).astype(int)
    lab["primary_time"] = np.where(lab.primary_event == 1, times["varc3"]["conf2"], t_end)
    lab["primary_competing"] = np.where(lab.primary_event == 1, 0, np.where(lab.death == 1, 2, np.where(lab.reintervention == 1, 3, 0)))
    # ---- latent truth (never for training)
    lt = df[["valve_id", "T_deg", "mode", "tau", "loss_fraction", "flow_slope", "onset_hr", "eoa_ref_true", "mpg_ref_true", "dvi_ref_true", "thromb_on", "thromb_start", "thromb_end"]].copy()
    lt["stage2_crossing_years"] = np.where(np.isfinite(latent["t2_latent"]), latent["t2_latent"], np.nan)
    lt["stage3_crossing_years"] = np.where(np.isfinite(latent["t3_latent"]), latent["t3_latent"], np.nan)
    lt["symptom_onset_years"] = np.where(np.isfinite(ends["t_sym"]), ends["t_sym"], np.nan)
    lt["death_time_years"] = ends["t_death"]
    return {"patients_valves": pv, "echo_visits": ev, "events": events, "labels": lab, "latent_truth": lt,
            "_wide": {"t": t, "vt": vt, "obs": obs, "ref": ref, "stages": stages, "df": df, "ends": ends}}


def simulate_matched(paper, theta_by_class, n, t_max, rng, p, label="varc3_single2", censor=None):
    """Cohort matched to a published paper; returns (times, states) with states 0 censored, 1 label, 2 death.
    label: 'varc3_single2' (Kermen: VARC-3 full, stage >= 2, single echo), 'varc3_haemo_single2' (NOTION),
    'varc3_haemo_single3', 'varc3_single3_or_reint' (Wakami)."""
    prof = COHORT_PROFILES[paper]
    p_local = p
    df, t, vt, obs, ref, stages, times, t_end, ends, latent = _core(p_local, n, rng, profile=prof, theta_by_class=theta_by_class,
                                                                     schedule=prof.schedule_y, reference_time=prof.reference_time_y, t_max=t_max)
    if label == "varc3_single2": tl = times["varc3"]["single2"]
    elif label == "varc3_haemo_single2": tl = times["varc3_haemo"]["single2"]
    elif label == "varc3_haemo_single3": tl = times["varc3_haemo"]["single3"]
    elif label == "varc3_single3": tl = times["varc3"]["single3"]
    elif label == "varc3_single3_or_reint": tl = np.minimum(times["varc3"]["single3"], np.where(np.isfinite(ends["t_re"]), ends["t_re"], np.inf))
    else: raise ValueError(label)
    t_death = ends["t_death"]
    cens = np.full(n, float(t_max)) if censor is None else np.minimum(censor(n, rng), float(t_max))
    cens = np.minimum(cens, np.where(np.isfinite(ends["t_re"]), ends["t_re"], np.inf))
    tt = np.minimum.reduce([tl, t_death, cens])
    state = np.where(tt == tl, 1, np.where(tt == t_death, 2, 0))
    return tt, state
