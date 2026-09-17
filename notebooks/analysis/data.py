"""Load the synthetic tables and shape them for modelling: one row per echo, NaN-padded per-valve history arrays."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

import analysis_config as C

MEAS = {"av_mean_gradient_mmHg": "mpg", "av_peak_velocity_m_s": "vmax", "av_area_cm2": "eoa", "doppler_velocity_index": "dvi",
        "stroke_volume_index_ml_m2": "svi", "lvef_pct": "lvef", "intraprosthetic_ar_grade": "ar", "paravalvular_leak_grade": "pvl"}
VTYPE = {"reference": 0, "scheduled": 1, "symptom_triggered": 2, "other_indication": 3}


def load_tables(data_dir=None) -> dict:
    d = C.DATA if data_dir is None else data_dir
    out = {}
    for name in ("patients_valves", "echo_visits", "events", "labels", "latent_truth"):
        pq = d / f"{name}.parquet"
        out[name] = pd.read_parquet(pq) if pq.exists() else pd.read_csv(d / f"{name}.csv")
    return out


def echo_wide(ev: pd.DataFrame) -> pd.DataFrame:
    """One row per (valve, echo) with the eight channels as columns."""
    w = ev.pivot_table(index=["valve_id", "visit_time_years", "visit_type"], columns="measurement", values="result", aggfunc="first").reset_index()
    w = w.rename(columns=MEAS)
    for c in MEAS.values():
        if c not in w:
            w[c] = np.nan
    w["vtype"] = w.visit_type.map(VTYPE).astype(int)
    w = w.rename(columns={"visit_time_years": "t"}).sort_values(["valve_id", "t"]).reset_index(drop=True)
    return w[["valve_id", "t", "vtype"] + list(MEAS.values())]


@dataclass
class Hist:
    """Post-reference echo history per valve, NaN-padded (n_valves, K). Row order = `valve_ids`."""
    valve_ids: np.ndarray
    t: np.ndarray
    vtype: np.ndarray
    mpg: np.ndarray
    vmax: np.ndarray
    eoa: np.ndarray
    dvi: np.ndarray
    svi: np.ndarray
    lvef: np.ndarray
    ar: np.ndarray
    n_echoes: np.ndarray      # (n_valves,)

    @property
    def K(self):
        return self.t.shape[1]


def history_arrays(wide: pd.DataFrame, valve_ids: np.ndarray) -> Hist:
    post = wide[wide.vtype != VTYPE["reference"]].copy()
    post["k"] = post.groupby("valve_id").cumcount()
    K = int(post.k.max()) + 1 if len(post) else 1
    n = len(valve_ids); pos = pd.Series(np.arange(n), index=valve_ids)
    r = pos.loc[post.valve_id].to_numpy(); c = post.k.to_numpy()
    def mat(col, fill=np.nan, dtype=float):
        m = np.full((n, K), fill, dtype=dtype); m[r, c] = post[col].to_numpy(); return m
    counts = np.zeros(n, int); np.add.at(counts, r, 1)
    return Hist(valve_ids=valve_ids, t=mat("t"), vtype=mat("vtype", fill=-1, dtype=int), mpg=mat("mpg"), vmax=mat("vmax"),
                eoa=mat("eoa"), dvi=mat("dvi"), svi=mat("svi"), lvef=mat("lvef"), ar=mat("ar"), n_echoes=counts)


def baseline_frame(pv: pd.DataFrame) -> pd.DataFrame:
    """Baseline and reference-echo columns with the names used by the feature sets."""
    b = pd.DataFrame({
        "valve_id": pv.valve_id.to_numpy(),
        "age_at_implant": pv.age_at_implant, "female": (pv.sex == "F").astype(int), "bsa": pv.bsa, "bmi": pv.bmi,
        "diabetes": pv.diabetes, "ckd": pv.ckd, "dialysis": pv.dialysis, "smoking": pv.smoking, "af": pv.af,
        "anticoag_type": pv.anticoag_type.map(C.CODES["anticoag_type"]).astype(float), "approach": pv.approach.map(C.CODES["approach"]).astype(float),
        "valve_model": pv.valve_model.map(C.CODES["valve_model"]).astype(float), "valve_class": pv.valve_class.map(C.CODES["valve_class"]).astype(float),
        "is_new_model": pv.is_new_model, "label_size_mm": pv.label_size_mm, "implant_year": pv.implant_year,
        "reference_from_table": (pv.reference_source == "model_size_table").astype(int),
        "ref_mpg": pv.ref_mean_gradient_mmHg, "ref_eoa": pv.ref_av_area_cm2, "ref_eoai": pv.ref_eoai_cm2_m2, "ref_dvi": pv.ref_dvi,
        "ref_ar": pv.ref_ar_grade, "ppm_class_ref": pv.ppm_class_ref, "expected_eoa": pv.table_eoa,
        "eoa_residual": pv.ref_av_area_cm2 - pv.table_eoa,
        "followup_end_years": pv.followup_end_years,
    })
    # flow-corrected reference gradient: the simulator's SVi at the reference is not stored in patients_valves;
    # use the ASE expected gradient for the model x size as the flow-neutral comparator instead
    b["ref_mpg_flow_corrected"] = pv.ref_mean_gradient_mmHg / pv.table_mpg
    return b.reset_index(drop=True)
