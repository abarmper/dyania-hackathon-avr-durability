"""SVD label definitions applied per echo against the reference echo; single-echo and confirmed variants;
first-event times per valve. All functions take (n, S) arrays with NaN for absent visits / channels."""
from __future__ import annotations

import numpy as np


def _ge(a, b):
    with np.errstate(invalid="ignore"):
        return np.nan_to_num(a >= b, nan=False) if isinstance(a, np.ndarray) else a >= b


def varc3_stage(mpg, eoa, dvi, ar, ref_mpg, ref_eoa, ref_dvi, ref_ar):
    """Returns (stage_full, stage_haemodynamic). stage_full is NaN where the flow-independent criterion
    cannot be evaluated (EOA and DVI both missing) and the gradient criterion is met."""
    with np.errstate(invalid="ignore"):
        d = mpg - ref_mpg
        eoa_fall = ref_eoa - eoa; eoa_pct = eoa_fall / ref_eoa
        dvi_fall = ref_dvi - dvi; dvi_pct = dvi_fall / ref_dvi
        haemo2 = (mpg >= 20) & (d >= 10); haemo3 = (mpg >= 30) & (d >= 20)
        fi2 = (eoa_fall >= 0.3) | (eoa_pct >= 0.25) | (dvi_fall >= 0.1) | (dvi_pct >= 0.2)
        fi3 = (eoa_fall >= 0.6) | (eoa_pct >= 0.5) | (dvi_fall >= 0.2) | (dvi_pct >= 0.4)
        new_mod_ar = (ar >= 2) & (ref_ar < 2); new_sev_ar = (ar >= 3) & (ref_ar < 3)
        both_missing = np.isnan(eoa) & np.isnan(dvi)
    haemo2 = np.nan_to_num(haemo2, nan=False); haemo3 = np.nan_to_num(haemo3, nan=False)
    fi2 = np.nan_to_num(fi2, nan=False); fi3 = np.nan_to_num(fi3, nan=False)
    new_mod_ar = np.nan_to_num(new_mod_ar, nan=False); new_sev_ar = np.nan_to_num(new_sev_ar, nan=False)
    st3 = (haemo3 & fi3) | new_sev_ar
    st2 = (haemo2 & fi2) | new_mod_ar
    stage = np.where(st3, 3.0, np.where(st2, 2.0, 0.0))
    na = both_missing & haemo2 & ~new_mod_ar & ~np.isnan(mpg)
    stage = np.where(na, np.nan, stage)
    stage = np.where(np.isnan(mpg), np.nan, stage)
    st3h = haemo3 | new_sev_ar; st2h = haemo2 | new_mod_ar
    stage_h = np.where(np.isnan(mpg), np.nan, np.where(st3h, 3.0, np.where(st2h, 2.0, 0.0)))
    return stage, stage_h


def capodanno(mpg, ref_mpg):
    with np.errstate(invalid="ignore"):
        d = mpg - ref_mpg
        mod = (d >= 10) & (mpg >= 20); sev = (d >= 20) & (mpg >= 40)
    out = np.where(np.nan_to_num(sev, nan=False), 3.0, np.where(np.nan_to_num(mod, nan=False), 2.0, 0.0))
    return np.where(np.isnan(mpg), np.nan, out)


def dvir(mpg, eoa, dvi, ref_mpg, ref_eoa, ref_dvi):
    with np.errstate(invalid="ignore"):
        ok = (mpg - ref_mpg > 10) & (eoa < ref_eoa) & (dvi < ref_dvi)
    out = np.where(np.nan_to_num(ok, nan=False), 2.0, 0.0)
    return np.where(np.isnan(mpg) | np.isnan(eoa) | np.isnan(dvi), np.nan, out)


def next_available(stage, t):
    """For each visit, the stage at the next visit with a non-NaN stage (NaN if none). Visits sorted by t per row."""
    n, S = stage.shape
    nxt = np.full_like(stage, np.nan)
    for j in range(S - 2, -1, -1):
        cand = stage[:, j + 1]
        nxt[:, j] = np.where(np.isnan(cand), nxt[:, j + 1] if j + 1 < S - 1 else np.nan, cand)
    return nxt


def confirmed(stage, t, level):
    """stage >= level at this visit AND at the next available visit."""
    nxt = next_available(stage, t)
    with np.errstate(invalid="ignore"):
        return np.nan_to_num(stage >= level, nan=False) & np.nan_to_num(nxt >= level, nan=False)


def first_time(mask, t):
    any_ = mask.any(1); idx = mask.argmax(1)
    return np.where(any_, t[np.arange(len(t)), idx], np.inf)


def first_times(stage_arrays: dict, t):
    """stage_arrays: name -> (n, S) stage array. Returns dict name -> {single2, single3, conf2, conf3, conf2_confirming}."""
    out = {}
    for name, st in stage_arrays.items():
        with np.errstate(invalid="ignore"):
            s2 = np.nan_to_num(st >= 2, nan=False); s3 = np.nan_to_num(st >= 3, nan=False)
        c2 = confirmed(st, t, 2); c3 = confirmed(st, t, 3)
        nxt_t = np.full_like(t, np.nan)
        for j in range(t.shape[1] - 1):
            nxt_t[:, j] = t[:, j + 1]
        out[name] = {"single2": first_time(s2, t), "single3": first_time(s3, t), "conf2": first_time(c2, t),
                     "conf3": first_time(c3, t), "conf2_confirming": first_time(c2, np.nan_to_num(nxt_t, nan=np.inf))}
    return out
