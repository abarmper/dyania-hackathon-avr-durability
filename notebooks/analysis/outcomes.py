"""Outcomes per landmark instance (C0 of docs/05): the event is known at the CONFIRMING echo; competing causes;
eligibility at fixed horizons; person-period expansion for the discrete-time cause-specific hazard model."""
from __future__ import annotations

import numpy as np
import pandas as pd

import analysis_config as C


def reintervention_times(events: pd.DataFrame) -> pd.DataFrame:
    """Per valve: time of reintervention and whether its indication was SVD (vs endocarditis)."""
    r = events[events.event_type == "reintervention"].copy()
    r["svd_indication"] = ~r.detail.str.contains("endocarditis")
    r = r.sort_values("time_years").drop_duplicates("valve_id")
    return r.set_index("valve_id")[["time_years", "svd_indication"]]


def known_times(labels: pd.DataFrame, definition: str = C.PRIMARY_LABEL, events: pd.DataFrame | None = None, svd_reint_as_event: bool = True) -> pd.DataFrame:
    """Per valve: T_known (time the label is known, or end of follow-up) and cause:
    1 = SVD label known (or, if svd_reint_as_event, a reintervention for SVD that truncated follow-up before the confirming echo:
        the 284 such valves in the cohort all carry an SVD indication), 2 = death or endocarditis (labels are voided after
        endocarditis), 0 = censored (lost to follow-up, administrative). The pure-label coding is the sensitivity analysis."""
    col = C.LABELS[definition][0]
    t_ev = labels[col].to_numpy(float)
    t_end = labels.followup_end_years.to_numpy(float)
    ev = np.isfinite(t_ev) & (t_ev <= t_end + 1e-9)
    T = np.where(ev, t_ev, t_end)
    cause = np.where(ev, 1, 0)
    if events is not None and svd_reint_as_event:
        r = reintervention_times(events).reindex(labels.valve_id.to_numpy())
        t_re = r.time_years.to_numpy(float); svd_re = np.where(r.svd_indication.isna(), False, r.svd_indication).astype(bool)
        re_event = svd_re & np.isfinite(t_re) & ~ev
        T = np.where(re_event, t_re, T); cause = np.where(re_event, 1, cause)
        ev = ev | re_event
    death = labels.death.to_numpy() == 1; endo = labels.endocarditis.to_numpy() == 1
    reint_other = (labels.reintervention.to_numpy() == 1) & ~ev
    cause = np.where(ev, 1, np.where(death | endo | reint_other, 2, 0))
    return pd.DataFrame({"valve_id": labels.valve_id.to_numpy(), "T_known": T, "cause": cause, "t_end": t_end})


def make_instances(landmarks: pd.DataFrame, labels: pd.DataFrame, definition: str = C.PRIMARY_LABEL, events: pd.DataFrame | None = None,
                   svd_reint_as_event: bool = True) -> pd.DataFrame:
    """Join outcomes to landmark rows; keep the risk set (T_known > landmark_t); add time_to and eligibility flags."""
    kt = known_times(labels, definition, events, svd_reint_as_event).set_index("valve_id")
    d = landmarks.copy()
    d["T_known"] = kt.T_known.loc[d.valve_id].to_numpy(); d["cause"] = kt.cause.loc[d.valve_id].to_numpy()
    d = d[d.T_known > d.landmark_t + 1e-9].copy()
    d["time_to"] = d.T_known - d.landmark_t
    for hzn in sorted(set(C.HORIZONS_Y) | set(C.LANDMARK_TAUS_Y)):
        # eligible at horizon: outcome within the horizon is known (event or competing event before it, or follow-up beyond it)
        d[f"elig_{int(hzn)}"] = ((d.cause > 0) & (d.time_to <= hzn)) | (d.time_to > hzn)
        d[f"y_{int(hzn)}"] = ((d.cause == 1) & (d.time_to <= hzn)).astype(int)
    d["currently_negative"] = (d.unconfirmed_flag == 0).astype(int)
    return d.reset_index(drop=True)


def person_periods(inst: pd.DataFrame, feature_cols: list[str], horizon_y: float = C.LANDMARK_HORIZON_Y, period_y: float = C.PERIOD_Y):
    """Expand each instance into periods 1..J, J = min(ceil(time_to / period), horizon / period).
    y = 0 in every period except the last one if the event/competing event happens within the horizon."""
    J_max = int(round(horizon_y / period_y))
    tt = inst.time_to.to_numpy(); cause = inst.cause.to_numpy()
    J = np.minimum(np.ceil(tt / period_y - 1e-9).astype(int), J_max); J = np.maximum(J, 1)
    rep = np.repeat(np.arange(len(inst)), J)
    period = np.concatenate([np.arange(1, j + 1) for j in J])
    last = np.concatenate([np.r_[np.zeros(j - 1, bool), True] for j in J])
    within = (tt <= horizon_y + 1e-9)
    y = np.where(last & within[rep] & (cause[rep] > 0), cause[rep], 0)
    X = inst.iloc[rep][feature_cols].reset_index(drop=True)
    X["period"] = period
    return X, y, rep
