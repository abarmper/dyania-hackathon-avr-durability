"""C9: analyses only a simulated cohort with an oracle can deliver. Writes tables/oracle_*.csv and figures."""
from __future__ import annotations

import json
import time

import joblib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed, parallel_config
from threadpoolctl import threadpool_limits

import analysis_config as C
from data import echo_wide, history_arrays, load_tables
from evaluate import auc_at, c_index_trunc, valve_resamples, bootstrap
from features import build_landmark_table
from models import DiscreteTimeHazard
from outcomes import make_instances

OUT = C.OUT; TAB = OUT / "tables"; FIG = OUT / "figures"; MOD = OUT / "models"


def _split_parts():
    T = load_tables(); L = pd.read_parquet(OUT / "landmark_table.parquet")
    sp = {k: set(v) for k, v in json.load(open(OUT / "split.json")).items()}
    part = lambda d, k: d[d.valve_id.isin(sp[k])].reset_index(drop=True)
    return T, L, part


def ceiling(T, L, part, tau=5.0):
    """Achievable-AUC ceiling: the same features trained on the LATENT stage-2 crossing, evaluated against the oracle label and
    the confirmed label on the test instances; decomposes the gap into label noise vs model/feature limits."""
    lat = T["latent_truth"].set_index("valve_id")
    lab = T["labels"].copy()
    lab["latent_cross"] = lat.stage2_crossing_years.reindex(lab.valve_id).to_numpy()
    C.LABELS["latent_cross"] = ("latent_cross", "latent stage-2 crossing (noise-free, continuous time)")
    inst_lat = make_instances(L, lab, "latent_cross", T["events"], svd_reint_as_event=False)
    inst_obs = make_instances(L, T["labels"], C.PRIMARY_LABEL, T["events"])
    with threadpool_limits(C.OMP_THREADS):
        m_lat = DiscreteTimeHazard(C.FULL).fit(part(inst_lat, "train"), part(inst_lat, "val"))
        m_obs = joblib.load(MOD / "dth_full.joblib")
    te_lat = part(inst_lat, "test"); te_obs = part(inst_obs, "test")
    oracle_t = T["labels"].set_index("valve_id")["oracle_varc3_single2"]
    rows = []
    for mname, m in (("trained_on_latent_crossing", m_lat), ("trained_on_confirmed_label", m_obs)):
        for ename, te, lt in (("latent_crossing", te_lat, None), ("confirmed_label", te_obs, None), ("oracle_label", te_obs, oracle_t.reindex(te_obs.valve_id).to_numpy(float))):
            cif = m.predict_cif(te, [tau])[:, 0]
            rows.append({"model": mname, "evaluated_against": ename, "AUC_5y": auc_at(te, cif, tau, label_time=lt), "C_5y": c_index_trunc(te, cif, tau) if lt is None else np.nan})
    d = pd.DataFrame(rows); d.to_csv(TAB / "oracle_ceiling.csv", index=False)
    return d


def definition_audit(T):
    """Sensitivity, false-positive share and lead time of each definition against the latent stage-2 crossing; kappa between definitions."""
    lab = T["labels"]; lat = T["latent_truth"].set_index("valve_id").reindex(lab.valve_id)
    t2 = lat.stage2_crossing_years.to_numpy(float); t_end = lab.followup_end_years.to_numpy()
    truth = np.isfinite(t2) & (t2 <= t_end)
    rows = []; flags = {}
    for d in ("varc3_single", "varc3_conf", "varc3_haemo_conf", "capodanno_conf", "dvir_conf", "oracle"):
        col = C.LABELS[d][0]; tl = lab[col].to_numpy(float); pos = np.isfinite(tl) & (tl <= t_end)
        flags[d] = pos
        tp = pos & truth; fp = pos & ~truth
        lead = (tl - t2)[tp]
        rows.append({"definition": d, "n_positive": int(pos.sum()), "sensitivity_for_latent_stage2": float(tp.sum() / max(truth.sum(), 1)),
                     "false_positive_share": float(fp.sum() / max(pos.sum(), 1)), "ppv": float(tp.sum() / max(pos.sum(), 1)),
                     "median_lag_after_crossing_y": float(np.median(lead)), "share_labelled_before_crossing": float((lead < 0).mean())})
    aud = pd.DataFrame(rows); aud.to_csv(TAB / "oracle_definition_audit.csv", index=False)
    names = list(flags); kap = pd.DataFrame(index=names, columns=names, dtype=float)
    for a in names:
        for b in names:
            x, y = flags[a], flags[b]; po = (x == y).mean(); pe = x.mean() * y.mean() + (1 - x.mean()) * (1 - y.mean())
            kap.loc[a, b] = (po - pe) / (1 - pe) if pe < 1 else 1.0
    kap.to_csv(TAB / "oracle_definition_kappa.csv")
    return aud, kap


def reference_value(T, L, part, tau=5.0):
    """Value of the 30-day reference echo: model performance and detection lag by reference source."""
    inst = make_instances(L, T["labels"], C.PRIMARY_LABEL, T["events"]); te = part(inst, "test")
    m = joblib.load(MOD / "dth_full.joblib"); cif = m.predict_cif(te, [tau])[:, 0]
    lat = T["latent_truth"].set_index("valve_id"); lab = T["labels"].set_index("valve_id")
    rows = []
    for src, mask in (("30d_echo", te.reference_from_table.to_numpy() == 0), ("model_size_table", te.reference_from_table.to_numpy() == 1)):
        sub = te[mask].reset_index(drop=True)
        vids = T["patients_valves"].valve_id[(T["patients_valves"].reference_source == src)].to_numpy()
        t2 = lat.stage2_crossing_years.reindex(vids).to_numpy(float); conf = lab.varc3_conf2_confirming.reindex(vids).to_numpy(float)
        ok = np.isfinite(t2) & np.isfinite(conf)
        rows.append({"reference_source": src, "n_instances": int(mask.sum()), "AUC_5y": auc_at(sub, cif[mask], tau), "C_5y": c_index_trunc(sub, cif[mask], tau),
                     "median_confirmed_lag_y": float(np.median(np.maximum(conf[ok] - t2[ok], 0))), "detection_rate_of_latent_stage2": float(ok.sum() / max(np.isfinite(t2).sum(), 1))})
    d = pd.DataFrame(rows); d.to_csv(TAB / "oracle_reference_value.csv", index=False)
    return d


def _regen_eval(kind, factor, seed=C.SEED + 3, n=10000, tau=5.0):
    """Regenerate a cohort under a mechanism shift with the calibrated parameters, rebuild the landmark table, score dth_full."""
    from threadpoolctl import threadpool_limits
    with threadpool_limits(4):
        from params import SimParams
        from physics import fit_physics_constants
        from simulate import simulate_cohort
        p = SimParams.from_json(C.PARAMS_CALIBRATED); fit_physics_constants(p)
        if kind == "noise":
            for a in ("sigma_flow", "sigma_cw", "sigma_lvot_vti", "sigma_lvot_diam"): setattr(p.noise, a, getattr(p.noise, a) * factor)
        elif kind == "missingness":
            p.visits.miss_by_year = tuple(min(0.9, x * factor) for x in p.visits.miss_by_year); p.visits.miss_saturation = min(0.9, p.visits.miss_saturation * factor)
        elif kind == "progression":
            p.trajectory.tau_median_y = {k: v * factor for k, v in p.trajectory.tau_median_y.items()}
        out = simulate_cohort(p, n=n, seed=seed)
        w = echo_wide(out["echo_visits"]); h = history_arrays(w, out["patients_valves"].valve_id.to_numpy())
        L = build_landmark_table(h, out["patients_valves"]); inst = make_instances(L, out["labels"], C.PRIMARY_LABEL, out["events"])
        m = joblib.load(MOD / "dth_full.joblib"); cif = m.predict_cif(inst, [tau])[:, 0]
        return {"shift": f"{kind} x{factor:g}", "n_instances": int(len(inst)), "AUC_5y": auc_at(inst, cif, tau), "C_5y": c_index_trunc(inst, cif, tau),
                "event_rate_5y": float(inst.y_5.mean())}


def mechanism_shifts():
    tasks = [("none", 1.0), ("noise", 0.5), ("noise", 1.5), ("missingness", 2.0), ("progression", 0.7)]
    with parallel_config(backend="loky", n_jobs=len(tasks), inner_max_num_threads=4):
        rows = Parallel()(delayed(_regen_eval)(k, f) for k, f in tasks)
    d = pd.DataFrame(rows); d.to_csv(TAB / "oracle_mechanism_shifts.csv", index=False)
    return d


def _sample_size_one(n_valves, rep, T, L, part, tau=5.0):
    from threadpoolctl import threadpool_limits
    with threadpool_limits(4):
        inst = make_instances(L, T["labels"], C.PRIMARY_LABEL, T["events"])
        rng = np.random.default_rng(C.SEED + 100 * rep + n_valves)
        vids = rng.choice(np.unique(inst.valve_id), size=n_valves, replace=False)
        sub = inst[inst.valve_id.isin(vids)].reset_index(drop=True)
        u = np.unique(sub.valve_id); rng.shuffle(u); n_te = max(int(0.3 * len(u)), 50)
        te = sub[sub.valve_id.isin(u[:n_te])].reset_index(drop=True); tr = sub[sub.valve_id.isin(u[n_te:])].reset_index(drop=True)
        m = DiscreteTimeHazard(C.FULL, max_iter=200).fit(tr)
        cif = m.predict_cif(te, [tau])[:, 0]
        rs = valve_resamples(te.valve_id.to_numpy(), 100, rep)
        auc = auc_at(te, cif, tau); lo, hi = bootstrap(lambda idx: auc_at(te.iloc[idx].reset_index(drop=True), cif[idx], tau), te, rs)
        return {"n_valves": n_valves, "rep": rep, "AUC_5y": auc, "ci_width": hi - lo, "n_test_cases": int(te.y_5.sum())}


def sample_size(T, L, part, sizes=(500, 1000, 2000, 5000), reps=10):
    tasks = [(n, r) for n in sizes for r in range(reps)]
    with parallel_config(backend="loky", n_jobs=min(len(tasks), C.N_JOBS), inner_max_num_threads=4):
        rows = Parallel()(delayed(_sample_size_one)(n, r, T, L, part) for n, r in tasks)
    d = pd.DataFrame(rows); d.to_csv(TAB / "oracle_sample_size_raw.csv", index=False)
    s = d.groupby("n_valves").agg(AUC_mean=("AUC_5y", "mean"), AUC_sd=("AUC_5y", "std"), ci_width_mean=("ci_width", "mean"), test_cases=("n_test_cases", "mean")).reset_index()
    s.to_csv(TAB / "oracle_sample_size.csv", index=False)
    return s


def run():
    t0 = time.time(); T, L, part = _split_parts()
    ce = ceiling(T, L, part); print("[oracle] ceiling:\n" + ce.round(3).to_string(index=False), flush=True)
    aud, kap = definition_audit(T); print("[oracle] definition audit:\n" + aud.round(3).to_string(index=False), flush=True)
    rv = reference_value(T, L, part); print("[oracle] reference value:\n" + rv.round(3).to_string(index=False), flush=True)
    ms = mechanism_shifts(); print("[oracle] mechanism shifts:\n" + ms.round(3).to_string(index=False), flush=True)
    ss = sample_size(T, L, part); print("[oracle] sample size:\n" + ss.round(3).to_string(index=False), flush=True)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 3.5)); ax.errorbar(ss.n_valves, ss.AUC_mean, yerr=ss.ci_width_mean / 2, fmt="o-"); ax.set_xscale("log")
    ax.set_xlabel("valves in the study"); ax.set_ylabel("5-y AUC (mean, half CI width)"); ax.axhline(0.75, ls="--", color="k", lw=1); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(FIG / "oracle_sample_size.png", dpi=130); plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 3.5)); ax.barh(aud.definition, aud.sensitivity_for_latent_stage2, color="C0", label="sensitivity for latent stage 2")
    ax.barh(aud.definition, -aud.false_positive_share, color="C3", label="false-positive share (negative axis)"); ax.axvline(0, color="k", lw=1); ax.legend(fontsize=8); ax.set_xlim(-0.5, 1)
    fig.tight_layout(); fig.savefig(FIG / "oracle_definition_audit.png", dpi=130); plt.close(fig)
    print(f"[oracle] done in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    run()
