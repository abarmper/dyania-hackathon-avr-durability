"""C10: secondary endpoints, descriptive. Time from first confirmed stage 2 to stage 3 / reintervention (KM by class),
urgent share, annualised gradient slope (mixed model), BVF and reintervention cumulative incidence (Aalen-Johansen)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from lifelines import AalenJohansenFitter, KaplanMeierFitter

import analysis_config as C
from data import echo_wide, load_tables

OUT = C.OUT; TAB = OUT / "tables"; FIG = OUT / "figures"


def run():
    T = load_tables(); pv, lab, ev = T["patients_valves"], T["labels"], T["events"]
    # ---- stage 2 -> stage 3 or reintervention
    t2 = lab.varc3_conf2_confirming.to_numpy(float); t3 = lab.varc3_single3.to_numpy(float)
    re = ev[ev.event_type == "reintervention"].set_index("valve_id").time_years.reindex(lab.valve_id).to_numpy(float)
    end = lab.followup_end_years.to_numpy()
    has2 = np.isfinite(t2)
    t_next = np.nanmin(np.vstack([np.where(np.isfinite(t3) & (t3 > t2), t3, np.inf), np.where(np.isfinite(re) & (re > t2), re, np.inf)]), 0)
    dur = np.where(np.isfinite(t_next), t_next - t2, end - t2); evt = np.isfinite(t_next)
    rows = []
    for cls in ["all"] + sorted(pv.valve_class.unique()):
        m = has2 & ((pv.valve_class.to_numpy() == cls) if cls != "all" else True)
        if m.sum() < 20: continue
        km = KaplanMeierFitter().fit(dur[m], evt[m])
        med = km.median_survival_time_
        rows.append({"valve_class": cls, "n_confirmed_stage2": int(m.sum()), "n_progressed_or_reintervened": int(evt[m].sum()),
                     "median_years_to_stage3_or_reintervention": float(med) if np.isfinite(med) else np.nan,
                     "share_within_2y": float(1 - km.predict(2.0))})
    pd.DataFrame(rows).to_csv(TAB / "secondary_stage2_to_stage3.csv", index=False)
    # ---- reintervention: urgent share, type
    r = ev[ev.event_type == "reintervention"].copy(); parts = r.detail.str.split("|", expand=True); r["type"] = parts[0]; r["indication"] = parts[1]; r["urgency"] = parts[2]
    (r.groupby(["indication", "type", "urgency"]).size().rename("n").reset_index()).to_csv(TAB / "secondary_reintervention_breakdown.csv", index=False)
    # ---- BVF (severe SVD confirmed or SVD reintervention) and reintervention CIF with death competing
    t_bvf = np.nanmin(np.vstack([np.where(np.isfinite(lab.varc3_conf3.to_numpy(float)), lab.varc3_conf3.to_numpy(float), np.inf), np.where(np.isfinite(re), re, np.inf)]), 0)
    ev_bvf = np.where(np.isfinite(t_bvf), 1, np.where(lab.death == 1, 2, 0)); tt = np.where(np.isfinite(t_bvf), t_bvf, end)
    aj = AalenJohansenFitter(calculate_variance=False).fit(tt, ev_bvf, event_of_interest=1)
    cif = aj.cumulative_density_
    at = {h: float(cif[cif.index <= h].iloc[-1, 0]) if (cif.index <= h).any() else 0.0 for h in (5, 8, 10)}
    ev_re = np.where(np.isfinite(re), 1, np.where(lab.death == 1, 2, 0)); tr = np.where(np.isfinite(re), re, end)
    aj2 = AalenJohansenFitter(calculate_variance=False).fit(tr, ev_re, event_of_interest=1); c2 = aj2.cumulative_density_
    at2 = {h: float(c2[c2.index <= h].iloc[-1, 0]) if (c2.index <= h).any() else 0.0 for h in (5, 8, 10)}
    pd.DataFrame([{"endpoint": "BVF (confirmed severe SVD or reintervention)", **{f"cif_{h}y": v for h, v in at.items()}},
                  {"endpoint": "any reintervention", **{f"cif_{h}y": v for h, v in at2.items()}}]).to_csv(TAB / "secondary_bvf_reintervention_cif.csv", index=False)
    # ---- annualised gradient slope: per-valve OLS on post-reference echoes (mixed model summary of the slope distribution)
    w = echo_wide(ev if False else T["echo_visits"]); w = w[w.vtype != 0]
    g = w.groupby("valve_id")
    def slope(d):
        if len(d) < 3: return np.nan
        x = d.t.to_numpy(); y = np.minimum(d.mpg.to_numpy(), 100); ok = ~np.isnan(y)
        return np.polyfit(x[ok], y[ok], 1)[0] if ok.sum() >= 3 else np.nan
    s = g.apply(slope).rename("slope_mmHg_per_y").reset_index().merge(pv[["valve_id", "approach", "valve_class"]], on="valve_id")
    try:
        import statsmodels.formula.api as smf
        d = w.merge(pv[["valve_id", "approach"]], on="valve_id"); d["mpg_c"] = np.minimum(d.mpg, 100)
        mm = smf.mixedlm("mpg_c ~ t * approach", d.dropna(subset=["mpg_c"]), groups=d.dropna(subset=["mpg_c"]).valve_id, re_formula="~t").fit(method="lbfgs", maxiter=200)
        fixed = mm.params.to_dict()
    except Exception as e:
        fixed = {"error": str(e)[:200]}
    summ = s.groupby("approach").slope_mmHg_per_y.describe(percentiles=[.1, .5, .9]).reset_index()
    summ.to_csv(TAB / "secondary_gradient_slope.csv", index=False)
    pd.Series(fixed).to_csv(TAB / "secondary_gradient_mixedlm_fixed_effects.csv")
    print("[secondary] stage2->3:\n" + pd.DataFrame(rows).round(2).to_string(index=False)); print("[secondary] BVF/reint CIF:", at, at2)
    print("[secondary] slope by approach:\n" + summ.round(2).to_string(index=False)); print("[secondary] mixed model fixed effects:", {k: round(v, 3) for k, v in fixed.items()} if "error" not in fixed else fixed)


if __name__ == "__main__":
    run()
