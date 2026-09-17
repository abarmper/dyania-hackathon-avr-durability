"""Explainability (C8) and feature selection (C5): SHAP on the boosted landmark model, permutation importance,
greedy backward elimination to a parsimonious model, per-inference output records."""
from __future__ import annotations

import json
import time

import joblib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed, parallel_config
from sklearn.inspection import permutation_importance
from threadpoolctl import threadpool_limits

import analysis_config as C
from evaluate import auc_at
from models import DiscreteTimeHazard

OUT = C.OUT; TAB = OUT / "tables"; FIG = OUT / "figures"; MOD = OUT / "models"


def _load():
    inst = pd.read_parquet(OUT / "instances.parquet")
    sp = {k: set(v) for k, v in json.load(open(OUT / "split.json")).items()}
    part = lambda k: inst[inst.valve_id.isin(sp[k])].reset_index(drop=True)
    return part("train"), part("val"), part("test"), joblib.load(MOD / "dth_full.joblib")


def _period_rows(model, inst, period):
    X = inst[model.feature_cols].to_numpy(float)
    return np.hstack([X, np.full((len(inst), 1), float(period))])


def shap_report(model: DiscreteTimeHazard, te: pd.DataFrame, n_background=3000, seed=C.SEED):
    """Global importance and per-instance contributions for the SVD hazard in the 5th period (~2.5 y after the landmark),
    a representative horizon; the sign convention is 'towards SVD'."""
    import shap
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(te), size=min(n_background, len(te)), replace=False)
    X = _period_rows(model, te.iloc[idx], period=5)
    names = model.feature_cols + ["period"]
    expl = shap.TreeExplainer(model.clf)
    sv = expl.shap_values(X)
    sv = sv[:, :, 1] if isinstance(sv, np.ndarray) and sv.ndim == 3 else (sv[1] if isinstance(sv, list) else sv)   # class 1 = SVD
    imp = pd.DataFrame({"feature": names, "mean_abs_shap": np.abs(sv).mean(0)}).sort_values("mean_abs_shap", ascending=False)
    imp.to_csv(TAB / "shap_importance.csv", index=False)
    fig = plt.figure(figsize=(7, 7)); shap.summary_plot(sv, X, feature_names=names, max_display=20, show=False); plt.title("SHAP, SVD hazard, period 5 (landmark model)", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "shap_summary.png", dpi=130, bbox_inches="tight"); plt.close(fig)
    for f in ("slope_since_ref", "ewma_log_mpg", "ref_eoai", "delta_mpg"):
        if f in names:
            j = names.index(f); fig, ax = plt.subplots(figsize=(4.5, 3.4)); ax.scatter(X[:, j], sv[:, j], s=4, alpha=0.4)
            ax.set_xlabel(f); ax.set_ylabel("SHAP (SVD hazard)"); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(FIG / f"shap_dependence_{f}.png", dpi=130); plt.close(fig)
    # per-instance top-3 for three example patients at the 5-y CIF terciles of the sampled instances
    cif5 = model.predict_cif(te.iloc[idx], [5.0])[:, 0]
    order = np.argsort(cif5); examples = [order[len(order) // 10], order[len(order) // 2], order[-len(order) // 20]]
    records = []
    for e in examples:
        contrib = pd.Series(sv[e], index=names).drop("period"); top = contrib.abs().sort_values(ascending=False).head(3).index.tolist()
        tier = "low" if cif5[e] < 0.05 else ("moderate" if cif5[e] < 0.15 else "high")
        records.append({"valve_id": int(te.iloc[idx[e]].valve_id), "landmark_t": float(te.iloc[idx[e]].landmark_t), "cif_5y": float(cif5[e]), "risk_tier": tier,
                        "top3": [{"feature": f, "value": float(X[e, names.index(f)]), "shap": float(contrib[f])} for f in top],
                        "unconfirmed_echo_flag": int(te.iloc[idx[e]].unconfirmed_flag)})
        fig = plt.figure(figsize=(7, 4)); shap.plots._waterfall.waterfall_legacy(expl.expected_value[1] if np.ndim(expl.expected_value) else expl.expected_value, sv[e], feature_names=names, max_display=10, show=False)
        fig.tight_layout(); fig.savefig(FIG / f"shap_waterfall_{tier}.png", dpi=130, bbox_inches="tight"); plt.close(fig)
    json.dump(records, open(TAB / "example_inference_records.json", "w"), indent=1)
    return imp, records


def permutation_report(model, te, tau=5.0, n_repeats=5, seed=C.SEED):
    """Permutation importance of the 5-y AUC (valve-level shuffling is not needed: rows are permuted within the test set)."""
    cols = model.feature_cols
    base = auc_at(te, model.predict_cif(te, [tau])[:, 0], tau)
    rng = np.random.default_rng(seed)
    def one(c):
        vals = []
        for r in range(n_repeats):
            d = te.copy(); d[c] = rng.permutation(d[c].to_numpy())
            vals.append(base - auc_at(d, model.predict_cif(d, [tau])[:, 0], tau))
        return c, float(np.mean(vals))
    with parallel_config(backend="loky", n_jobs=min(16, C.N_JOBS), inner_max_num_threads=4):
        res = Parallel()(delayed(one)(c) for c in cols)
    imp = pd.DataFrame(res, columns=["feature", "auc_drop"]).sort_values("auc_drop", ascending=False)
    imp.to_csv(TAB / "permutation_importance_5y.csv", index=False)
    return imp, base


def parsimonious_model(tr, va, te, ranking: list[str], full_auc: float, tau=5.0, max_features=8, keep_frac=0.95):
    """Greedy backward elimination along the permutation ranking (least important first) until <= max_features,
    stopping earlier if the AUC would fall below keep_frac x full AUC."""
    cols = list(ranking); history = []
    while len(cols) > max_features:
        cand = cols[-1]; trial = [c for c in cols if c != cand]
        m = DiscreteTimeHazard(trial).fit(tr, va); auc = auc_at(te, m.predict_cif(te, [tau])[:, 0], tau)
        history.append({"n_features": len(trial), "dropped": cand, "AUC": auc})
        if auc < keep_frac * full_auc:
            break
        cols = trial
    m = DiscreteTimeHazard(cols).fit(tr, va); auc = auc_at(te, m.predict_cif(te, [tau])[:, 0], tau)
    joblib.dump(m, MOD / "dth_parsimonious.joblib")
    return cols, auc, pd.DataFrame(history)


def run():
    t0 = time.time(); tr, va, te, model = _load()
    with threadpool_limits(C.OMP_THREADS):
        imp, records = shap_report(model, te)
        print(f"[explain] SHAP top: {imp.head(8).feature.tolist()} ({time.time()-t0:.0f}s)", flush=True)
        pimp, base = permutation_report(model, te)
        print(f"[explain] permutation top: {pimp.head(8).feature.tolist()} base AUC {base:.3f} ({time.time()-t0:.0f}s)", flush=True)
        ranking = pimp.feature.tolist()
        cols, auc, hist = parsimonious_model(tr, va, te, ranking, base)
    hist.to_csv(TAB / "backward_elimination.csv", index=False)
    json.dump({"parsimonious_features": cols, "AUC_5y": auc, "full_AUC_5y": base, "shap_top10": imp.head(10).to_dict("records")}, open(TAB / "parsimonious_model.json", "w"), indent=1)
    print(f"[explain] parsimonious {len(cols)} features AUC {auc:.3f} vs full {base:.3f}: {cols} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    run()
