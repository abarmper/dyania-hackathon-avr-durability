"""End-to-end pipeline. Usage: python run_all.py [features|models|evaluate|explain|oracle|policy|secondary|all]

Writes notebooks/analysis/output/: landmark_table.parquet, instances.parquet, models/*.joblib, metrics.json, tables/*.md, figures/*.png,
results_summary.md. Every stage can be re-run alone; it loads what the previous stage wrote."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import analysis_config as C
from data import baseline_frame, echo_wide, history_arrays, load_tables
from evaluate import label_sweep, metrics_table, subgroup_table, valve_resamples
from features import build_landmark_table, implant_instances
from models import CoxCauseSpecific, DiscreteTimeHazard, RiskFactorCoxB4, TimeSinceImplantOnly
from outcomes import known_times, make_instances

OUT = C.OUT; TAB = OUT / "tables"; FIG = OUT / "figures"; MOD = OUT / "models"
for d in (OUT, TAB, FIG, MOD):
    d.mkdir(parents=True, exist_ok=True)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def split_valves(valve_ids, seed=C.SEED):
    rng = np.random.default_rng(seed); u = np.unique(valve_ids); rng.shuffle(u)
    n_test = int(len(u) * C.TEST_FRAC); n_val = int(len(u) * C.VAL_FRAC)
    return {"test": set(u[:n_test]), "val": set(u[n_test:n_test + n_val]), "train": set(u[n_test + n_val:])}


def stage_features():
    T = load_tables()
    w = echo_wide(T["echo_visits"]); h = history_arrays(w, T["patients_valves"].valve_id.to_numpy())
    L = build_landmark_table(h, T["patients_valves"]); I0 = implant_instances(T["patients_valves"])
    L.to_parquet(OUT / "landmark_table.parquet", index=False); I0.to_parquet(OUT / "implant_table.parquet", index=False)
    log(f"features: landmark {L.shape}, implant {I0.shape}")
    return T, L, I0


def stage_models(T=None, L=None, I0=None):
    if T is None:
        T = load_tables(); L = pd.read_parquet(OUT / "landmark_table.parquet"); I0 = pd.read_parquet(OUT / "implant_table.parquet")
    inst = make_instances(L, T["labels"], C.PRIMARY_LABEL, T["events"])
    inst0 = make_instances(I0, T["labels"], C.PRIMARY_LABEL, T["events"])
    sp = split_valves(T["patients_valves"].valve_id.to_numpy())
    part = lambda d, k: d[d.valve_id.isin(sp[k])].reset_index(drop=True)
    inst.to_parquet(OUT / "instances.parquet", index=False); inst0.to_parquet(OUT / "instances_implant.parquet", index=False)
    json.dump({k: sorted(int(x) for x in v) for k, v in sp.items()}, open(OUT / "split.json", "w"))
    models = {}
    with threadpool_limits(C.OMP_THREADS):
        t0 = time.time()
        models["dth_full"] = DiscreteTimeHazard(C.FULL).fit(part(inst, "train"), part(inst, "val"))
        log(f"dth_full fit: {models['dth_full'].n_train_rows_} rows, {models['dth_full'].n_iter_} iters, {time.time()-t0:.0f}s")
        t0 = time.time(); models["dth_registry"] = DiscreteTimeHazard(C.REGISTRY_ONLY).fit(part(inst, "train"), part(inst, "val")); log(f"dth_registry fit {time.time()-t0:.0f}s")
        t0 = time.time(); models["dth_nosymptom"] = DiscreteTimeHazard([c for c in C.FULL if c not in ("symptom_now", "any_symptom_before")]).fit(part(inst, "train"), part(inst, "val")); log(f"dth_nosymptom fit {time.time()-t0:.0f}s")
        t0 = time.time(); models["dth_time_only"] = TimeSinceImplantOnly().fit(part(inst, "train"), part(inst, "val")); log(f"time-only fit {time.time()-t0:.0f}s")
        t0 = time.time(); models["dth_implant"] = DiscreteTimeHazard(C.BASELINE + C.REFERENCE_ECHO, horizon_y=C.IMPLANT_HORIZON_Y).fit(part(inst0, "train"), part(inst0, "val")); log(f"implant-time fit {time.time()-t0:.0f}s")
        t0 = time.time(); models["cox_full"] = CoxCauseSpecific(C.BASELINE + C.REFERENCE_ECHO).fit(part(inst0, "train")); log(f"cox fit {time.time()-t0:.0f}s")
        models["cox_b4"] = RiskFactorCoxB4().fit(part(inst0, "train"))
    for k, m in models.items():
        joblib.dump(m, MOD / f"{k}.joblib")
    models["cox_full"].hr_table().to_csv(TAB / "cox_hr_table.csv", index=False)
    return T, inst, inst0, sp, models


def stage_evaluate(T=None, inst=None, inst0=None, sp=None, models=None):
    if models is None:
        T = load_tables(); inst = pd.read_parquet(OUT / "instances.parquet"); inst0 = pd.read_parquet(OUT / "instances_implant.parquet")
        sp = {k: set(v) for k, v in json.load(open(OUT / "split.json")).items()}
        models = {p.stem: joblib.load(p) for p in MOD.glob("*.joblib")}
    part = lambda d, k: d[d.valve_id.isin(sp[k])].reset_index(drop=True)
    te = part(inst, "test"); te0 = part(inst0, "test")
    rs = valve_resamples(te.valve_id.to_numpy(), C.BOOT_B, C.SEED + 1); rs0 = valve_resamples(te0.valve_id.to_numpy(), C.BOOT_B, C.SEED + 2)
    results = {"n_test_instances": int(len(te)), "n_test_valves": int(te.valve_id.nunique()), "n_test_implant": int(len(te0))}
    tables = []
    with threadpool_limits(C.OMP_THREADS):
        # landmark models at 3 and 5 y
        for name in ("dth_full", "dth_registry", "dth_nosymptom", "dth_time_only"):
            m = models[name]; cif = {tau: m.predict_cif(te, [tau])[:, 0] for tau in C.LANDMARK_TAUS_Y}
            for subset, mask in (("all", None), ("currently_negative", te.currently_negative.to_numpy() == 1)):
                t = metrics_table(te, cif, C.LANDMARK_TAUS_Y, rs, mask); t.insert(0, "subset", subset); t.insert(0, "model", name); tables.append(t)
            if name == "dth_full":
                sw = label_sweep(te, cif[5.0], 5.0, T["labels"], ["varc3_conf", "varc3_single", "varc3_haemo_conf", "capodanno_conf", "dvir_conf", "oracle"])
                sw.to_csv(TAB / "label_sweep_5y.csv", index=False); results["label_sweep_5y"] = sw.to_dict("records")
                sg = subgroup_table(te, cif[5.0], 5.0, T["patients_valves"]); sg.to_csv(TAB / "subgroups_5y.csv", index=False)
                results["subgroups_5y"] = sg.to_dict("records")
                tab, slope, icpt = __import__("evaluate").calibration_at(te, cif[5.0], 5.0); tab.to_csv(TAB / "calibration_5y_dth_full.csv", index=False)
        # implant-time models at 5/8/10 y
        for name in ("dth_implant", "cox_full", "cox_b4"):
            m = models[name]
            cif = {tau: m.predict_cif(te0, [tau])[:, 0] for tau in C.HORIZONS_Y}
            t = metrics_table(te0, cif, C.HORIZONS_Y, rs0); t.insert(0, "subset", "implant_time"); t.insert(0, "model", name); tables.append(t)
    metrics = pd.concat(tables, ignore_index=True); metrics.to_csv(TAB / "metrics.csv", index=False)
    results["metrics"] = metrics.round(4).to_dict("records")
    # usefulness thresholds (docs/03 section 0) on the primary model, confirmed label, 5 y, all instances
    row = metrics[(metrics.model == "dth_full") & (metrics.subset == "all") & (metrics.tau == 5.0)].iloc[0]
    results["usefulness"] = {"AUC_5y": float(row.AUC), "AUC_ge_0.75": bool(row.AUC >= 0.75), "C_5y": float(row.C_trunc), "C_ge_0.70": bool(row.C_trunc >= 0.70),
                             "calib_slope_5y": float(row.calib_slope), "slope_in_0.8_1.2": bool(0.8 <= row.calib_slope <= 1.2)}
    json.dump(results, open(OUT / "metrics.json", "w"), indent=1, default=float)
    _figures_evaluate(metrics, results)
    log("evaluate: " + json.dumps(results["usefulness"]))
    return results


def _figures_evaluate(metrics, results):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    sw = pd.DataFrame(results["label_sweep_5y"])
    fig, ax = plt.subplots(figsize=(7, 3.6)); ax.barh(sw.label, sw.AUC, color="C0"); ax.axvline(0.75, ls="--", color="k", lw=1, label="usefulness threshold 0.75")
    ax.set_xlabel("time-dependent AUC at 5 y (same instances, same predictions)"); ax.set_xlim(0.5, 1.0); ax.legend(fontsize=8); ax.set_title("Discrimination vs label definition (dth_full)", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "auc_vs_label_5y.png", dpi=130); plt.close(fig)
    tab = pd.read_csv(TAB / "calibration_5y_dth_full.csv")
    fig, ax = plt.subplots(figsize=(4.2, 4)); ax.plot([0, tab.pred.max() * 1.1], [0, tab.pred.max() * 1.1], "k--", lw=1); ax.plot(tab.pred, tab.obs, "o-", color="C1")
    ax.set_xlabel("predicted 5-y CIF (decile mean)"); ax.set_ylabel("observed fraction with confirmed SVD by 5 y"); ax.set_title("Calibration, landmark model, 5 y", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "calibration_5y.png", dpi=130); plt.close(fig)
    m = metrics[(metrics.subset.isin(["all", "implant_time"])) & (metrics.tau == 5.0)]
    fig, ax = plt.subplots(figsize=(7, 3.6)); ax.barh(m.model, m.AUC, xerr=[m.AUC - m.AUC_lo, m.AUC_hi - m.AUC], color="C2"); ax.axvline(0.75, ls="--", color="k", lw=1)
    ax.set_xlabel("AUC at 5 y (95% valve-bootstrap CI)"); ax.set_xlim(0.5, 1.0); fig.tight_layout(); fig.savefig(FIG / "auc_models_5y.png", dpi=130); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(); ap.add_argument("stage", nargs="?", default="all")
    a = ap.parse_args(argv)
    t0 = time.time()
    if a.stage in ("features", "all"):
        T, L, I0 = stage_features()
    if a.stage in ("models", "all"):
        T, inst, inst0, sp, models = stage_models(*(locals().get(k) for k in ("T", "L", "I0"))) if a.stage == "all" else stage_models()
    if a.stage in ("evaluate", "all"):
        stage_evaluate(*( (T, inst, inst0, sp, models) if a.stage == "all" else (None,) * 5))
    if a.stage == "explain":
        import explain; explain.run()
    if a.stage == "oracle":
        import oracle_analyses; oracle_analyses.run()
    if a.stage == "policy":
        import policy; policy.run()
    if a.stage == "secondary":
        import secondary; secondary.run()
    log(f"done ({a.stage}) in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
