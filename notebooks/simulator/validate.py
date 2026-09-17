"""Acceptance tests for a generated cohort (docs/04 section 5) -> output/validation_report.md + figures/.

Run:  /data/abar/alexenv/bin/python validate.py [--params output/params_calibrated.json] [--n 10000]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
KM = HERE.parent / "km_reconstruct"
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(KM))

from calibrate import CHECKS, CURVE_STEPS, Report, TARGETS, matched_curve, perigon_stats, risk_ok  # noqa: E402
from calibrate_onset import censoring_sampler, observed_curve  # noqa: E402
from params import SimParams, default_params  # noqa: E402
from physics import fit_physics_constants  # noqa: E402
from simulate import _core, simulate_cohort  # noqa: E402
from survival import aj, km, step_eval  # noqa: E402
from valve_tables import COHORT_PROFILES, EOA_TABLE  # noqa: E402

OUT = HERE / "output"; FIG = OUT / "figures"


def check_reference_echo(p, rep, n, seed):
    rng = np.random.default_rng(seed); prof = COHORT_PROFILES["kermen"]
    df, t, vt, obs, ref, stages, times, t_end, ends, latent = _core(p, n, rng, profile=prof, schedule=prof.schedule_y, reference_time=prof.reference_time_y, t_max=12)
    ppm = pd.Series(df.ppm_true).value_counts(normalize=True)
    rep.add("reference echo", "Kermen Magna Ease: discharge MPG (mmHg)", "12.6 +/- 3.0", f"{np.nanmean(ref['mpg']):.1f} +/- {np.nanstd(ref['mpg']):.1f}", "mean +/-2", abs(np.nanmean(ref["mpg"]) - 12.6) <= 2.0)
    rep.add("reference echo", "Kermen: EOA (cm2)", "1.6 +/- 0.3", f"{df.eoa_ref_true.mean():.2f} +/- {df.eoa_ref_true.std():.2f}", "mean +/-0.15", abs(df.eoa_ref_true.mean() - 1.6) <= 0.15)
    rep.add("reference echo", "Kermen: EOAi (cm2/m2)", "0.87 +/- 0.17", f"{df.eoai_true.mean():.2f} +/- {df.eoai_true.std():.2f}", "mean +/-0.08", abs(df.eoai_true.mean() - 0.87) <= 0.08)
    rep.add("reference echo", "Kermen: moderate / severe PPM (%)", "51 / 2", f"{100*ppm.get(1,0):.0f} / {100*ppm.get(2,0):.0f}", "+/-10 / +/-4", abs(100 * ppm.get(1, 0) - 51) <= 10 and abs(100 * ppm.get(2, 0) - 2) <= 4)
    rng = np.random.default_rng(seed + 1); prof = COHORT_PROFILES["perigon"]
    df, t, vt, obs, ref, stages, times, t_end, ends, latent = _core(p, n, rng, profile=prof, schedule=prof.schedule_y, reference_time=prof.reference_time_y, t_max=5.5)
    rep.add("reference echo", "Velders (Avalus ~ Perimount-class): MPG / EOA / DVI at discharge", "13.1 +/- 4.7 / 1.54 +/- 0.36 / 0.49 +/- 0.10",
            f"{np.nanmean(ref['mpg']):.1f} +/- {np.nanstd(ref['mpg']):.1f} / {np.nanmean(ref['eoa']):.2f} +/- {np.nanstd(ref['eoa']):.2f} / {np.nanmean(ref['dvi']):.2f} +/- {np.nanstd(ref['dvi']):.2f}",
            "MPG +/-2, EOA +/-0.2, DVI +/-0.08", abs(np.nanmean(ref["mpg"]) - 13.1) <= 2 and abs(np.nanmean(ref["eoa"]) - 1.54) <= 0.2 and abs(np.nanmean(ref["dvi"]) - 0.49) <= 0.08)


def check_models_vs_ase(pv, rep):
    rows = []
    for model, tab in EOA_TABLE.items():
        sub = pv[pv.valve_model == model]
        if len(sub) < 30 or model == "Perimount (classic)":
            continue
        # compare per size where n >= 20
        for size, r in tab.items():
            s2 = sub[sub.label_size_mm == size]
            if len(s2) >= 20:
                rows.append((model, size, r.mpg, s2.ref_mean_gradient_mmHg.mean(), r.eoa, s2.ref_av_area_cm2.mean()))
    d = pd.DataFrame(rows, columns=["model", "size", "ase_mpg", "sim_mpg", "ase_eoa", "sim_eoa"])
    worst = float((d.sim_mpg - d.ase_mpg).abs().max()) if len(d) else 0.0
    rep.add("reference echo", "per model x size mean gradient vs ASE 2024 rows (mmHg), worst |diff|", "", round(worst, 2), 3.5, worst <= 3.5,
            f"{len(d)} model-size cells with n>=20")
    return d


def check_velders(p, rep, n, seed):
    val = perigon_stats(p, n, seed)
    rep.add("label flicker (A3)", "ever labelled by 5 y, Capodanno / Dvir / VARC-3 (%)", "4.6 / 2.9 / 3.0", f"{val['ever_capodanno']:.1f} / {val['ever_dvir']:.1f} / {val['ever_varc3']:.1f}", "+/-1.5",
            abs(val["ever_capodanno"] - 4.6) <= 1.5 and abs(val["ever_dvir"] - 2.9) <= 1.5 and abs(val["ever_varc3"] - 3.0) <= 1.5)
    rep.add("label flicker (A3)", "VARC-3 persistence at the next visit (%)", "~33 (published range 14-50)", round(val["persistence"], 0), "14-50", 14 <= val["persistence"] <= 50,
            f"next visit missing {val['next_missing']:.0f}%")
    rep.add("label flicker (A3)", "regression to the mean: lowest / highest decile change (mmHg per interval)", "+1.2..+2.3 / -1.0..-5.9",
            f"{val['lowest_decile']:+.1f} / {val['highest_decile']:+.1f}", "", 0.2 <= val["lowest_decile"] <= 3.3 and -6.9 <= val["highest_decile"] <= 0.0)
    rep.add("label flicker (A3)", "95% PI of within-patient 5-y change (mmHg)", "(-9.6, +7.5)", f"({val['pi'][0]:.1f}, {val['pi'][1]:.1f})", "width +/-4", abs(val["pi_width"] - 17.1) <= 4)
    rep.add("label flicker (A3)", "per-visit VARC-3 prevalence, years 1-5 (%)", "1.0 / 1.1 / 1.4 / 1.4 / 0.5", "/".join(f"{x:.1f}" for x in val["per_visit_prevalence"]), "0.2-3", all(0.2 <= x <= 3.0 for x in val["per_visit_prevalence"]))
    return val


def check_curves(p, rep, n, seed):
    curves = {}
    for tid, paper, cls, label, tol in CURVE_STEPS:
        target = TARGETS[tid]; grid, sim, tv = matched_curve(p, paper, label, target, n, seed, censor=censoring_sampler(tid))
        m = risk_ok(target, grid); diff = float(np.max(np.abs(sim - tv)[m]))
        rep.add("curve targets", f"{tid}: max |sim - target| while n>=30 (pp)", "", round(diff, 2), tol, diff <= tol, f"at {min(len(grid),10)} y sim {sim[min(len(grid),10)-1]:.1f} vs {tv[min(len(grid),10)-1]:.1f}")
        curves[tid] = (grid, sim, tv)
    tid = "notion_fig3_modsvd_savr"; target = TARGETS[tid]
    grid, sim, tv = matched_curve(p, "notion_savr", "varc3_haemo_single2", target, n, seed, censor=censoring_sampler(tid))
    diff = float(np.max(np.abs(sim - tv)[risk_ok(target, grid)]))
    rep.add("curve targets", f"{tid}: max |sim - target| (pp)", "", round(diff, 2), 3.0, diff <= 3.0); curves[tid] = (grid, sim, tv)
    for tid, paper, label, tol in CHECKS:
        target = TARGETS[tid]; grid, sim, tv = matched_curve(p, paper, label, target, n, seed, censor=censoring_sampler(tid))
        m = risk_ok(target, grid); diff = float(np.max(np.abs(sim - tv)[m])) if m.any() else float(np.max(np.abs(sim - tv)))
        rep.add("curve checks", f"{tid}: max |sim - target| while n>=30 (pp)", "", round(diff, 2), tol, diff <= tol); curves[tid] = (grid, sim, tv)
    target = TARGETS["kermen_fig2a_overall_survival"]; grid, sim, tv = matched_curve(p, "kermen", "death_only", target, n, seed) if False else (None, None, None)
    return curves


def check_death_reint(p, rep, n, seed):
    rng = np.random.default_rng(seed); prof = COHORT_PROFILES["kermen"]
    df, t, vt, obs, ref, stages, times, t_end, ends, latent = _core(p, n, rng, profile=prof, schedule=prof.schedule_y, reference_time=prof.reference_time_y, t_max=12.5)
    target = TARGETS["kermen_fig2a_overall_survival"]; grid = np.arange(1.0, 11.0)
    tt = np.minimum(ends["t_death"], 12.5); st = (ends["t_death"] <= 12.5).astype(int)
    s = observed_curve(tt, st, grid, "KM"); tv = 100 * step_eval(grid, np.array(target["grid_t"]), np.array(target["grid_value_pct"]) / 100)
    diff = float(np.max(np.abs(s - tv)))
    rep.add("death", "Kermen overall survival 1-10 y, max |sim - target| (pp)", "", round(diff, 2), 3.0, diff <= 3.0, f"10 y sim {s[-1]:.1f} vs {tv[-1]:.1f}")
    det3 = np.isfinite(times["varc3"]["single3"]); re = np.isfinite(ends["t_re"])
    share = float(100 * np.mean(re[det3])) if det3.any() else np.nan
    rep.add("reintervention", "Kermen: share of detected stage-3 valves reintervened (%)", "9/11 = 82", round(share, 0), ">= 55", share >= 55)
    rng = np.random.default_rng(seed + 1); prof = COHORT_PROFILES["notion_tavi"]
    df, t, vt, obs, ref, stages, times, t_end, ends, latent = _core(p, n, rng, profile=prof, schedule=prof.schedule_y, reference_time=prof.reference_time_y, t_max=10.0)
    dead10 = float(100 * np.mean(ends["t_death"] <= 10.0)); r10 = float(100 * np.mean(np.isfinite(ends["t_re"]) & (ends["t_re"] <= 10)))
    rep.add("death", "NOTION TAVI all-cause mortality at 10 y (%)", 62.7, round(dead10, 1), 6.0, abs(dead10 - 62.7) <= 6.0)
    rep.add("reintervention", "NOTION TAVI reintervention by 10 y (%)", 4.3, round(r10, 2), "1-7", 1.0 <= r10 <= 7.0)


def check_emergent_hr(out, rep):
    from lifelines import CoxPHFitter
    pv, lab = out["patients_valves"], out["labels"]
    d = pv.merge(lab[["valve_id", "primary_event", "primary_time"]], on="valve_id")
    d = d[d.primary_time > 0.1].copy()
    d["severe_ppm"] = (d.ppm_class_ref == 2).astype(int); d["mpg15"] = (d.ref_mean_gradient_mmHg >= 15).astype(int)
    d["thv_le23"] = ((d.approach == "TAVR") & (d.label_size_mm <= 23)).astype(int); d["female"] = (d.sex == "F").astype(int)
    d["tavr"] = (d.approach == "TAVR").astype(int)
    cols = ["age_at_implant", "severe_ppm", "mpg15", "thv_le23", "female", "smoking", "diabetes", "ckd", "bmi", "tavr"]
    cph = CoxPHFitter().fit(d[cols + ["primary_time", "primary_event"]], "primary_time", "primary_event")
    hr = np.exp(cph.params_)
    checks = [("age_at_implant", "age per year", 0.951, (0.93, 0.98)), ("severe_ppm", "severe PPM (total effect)", 1.85, (1.3, 2.8)),
              ("mpg15", "reference MPG >= 15 mmHg", 1.30, (1.05, 1.8)), ("thv_le23", "THV <= 23 mm", 2.07, (1.3, 3.2)),
              ("female", "female sex (null)", 1.0, (0.85, 1.18)), ("smoking", "smoking", 2.58, (1.9, 3.4)), ("diabetes", "diabetes", 1.33, (1.05, 1.7))]
    for c, name, tgt, (lo, hi) in checks:
        rep.add("emergent HR (Cox on confirmed VARC-3 >=2)", name, tgt, round(float(hr[c]), 2), f"{lo}-{hi}", lo <= hr[c] <= hi)
    return hr


def descriptives(out, rep):
    pv, ev, events, lab = out["patients_valves"], out["echo_visits"], out["events"], out["labels"]
    n = len(pv); vis = ev.drop_duplicates(["valve_id", "visit_time_years"])
    t = lab.primary_time.to_numpy(); st = np.where(lab.primary_event == 1, 1, np.where(lab.death == 1, 2, 0))
    cif = aj(t, st); c1 = cif.get(1, (np.array([0.0]), np.array([0.0])))
    at = 100 * step_eval(np.array([5.0, 8.0, 10.0]), c1[0], c1[1])
    rep.add("dataset", "valves / echoes / events", "", f"{n} / {len(vis)} / {len(events)}", "", True)
    rep.add("dataset", "primary endpoint (confirmed VARC-3 >=2) cumulative incidence at 5 / 8 / 10 y (%)", "", "/".join(f"{x:.1f}" for x in at), "", True,
            f"{int(lab.primary_event.sum())} events; single-echo labels {int(np.isfinite(lab.varc3_single2).sum())}")
    first_single = np.isfinite(lab.varc3_single2); conf = np.isfinite(lab.varc3_conf2)
    rep.add("dataset", "share of single-echo VARC-3 labels never confirmed (%)", "Velders: 59-65 absent at next visit", round(100 * (1 - conf.sum() / max(first_single.sum(), 1)), 0), "", True)
    rep.add("dataset", "echoes per valve (mean) / follow-up years (mean)", "", f"{len(vis)/n:.1f} / {pv.followup_end_years.mean():.1f}", "", True)
    rep.add("dataset", "events: " + ", ".join(f"{k} {v}" for k, v in events.event_type.value_counts().items()), "", "", "", True)
    urgent = events[(events.event_type == "reintervention")].detail.str.contains("urgent").mean() if (events.event_type == "reintervention").any() else np.nan
    rep.add("dataset", "reinterventions urgent share (%)", "", round(100 * urgent, 0), "", True)
    rep.add("dataset", "valve models: " + ", ".join(f"{k} {v}" for k, v in pv.valve_model.value_counts().items()), "", "", "", True, "provided notes: ~10 models, 3 THV generations, sizes 19-29 (docs/01)")


def figures(curves, out, val, model_tab):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    FIG.mkdir(parents=True, exist_ok=True)
    ids = list(curves); nc = 3; nr = int(np.ceil(len(ids) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(4.2 * nc, 3.2 * nr), squeeze=False)
    for ax, tid in zip(axes.ravel(), ids):
        g, s, tv = curves[tid]; ax.plot(g, tv, "k-", lw=2, label="published (reconstructed)"); ax.plot(g, s, "o--", color="C1", ms=4, label="simulated")
        ax.set_title(tid, fontsize=8); ax.set_xlabel("years"); ax.grid(alpha=.3); ax.legend(fontsize=7)
    for ax in axes.ravel()[len(ids):]: ax.axis("off")
    fig.tight_layout(); fig.savefig(FIG / "curves_sim_vs_target.png", dpi=130); plt.close(fig)
    if len(model_tab):
        fig, ax = plt.subplots(figsize=(8, 3.5)); x = np.arange(len(model_tab))
        ax.bar(x - 0.2, model_tab.ase_mpg, 0.4, label="ASE 2024 table"); ax.bar(x + 0.2, model_tab.sim_mpg, 0.4, label="simulated reference")
        ax.set_xticks(x); ax.set_xticklabels([f"{m[:10]} {s}" for m, s in zip(model_tab.model, model_tab["size"])], rotation=70, fontsize=7); ax.set_ylabel("mean gradient (mmHg)"); ax.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(FIG / "reference_mpg_vs_ase.png", dpi=130); plt.close(fig)
    # label-flicker heatmap (Velders Fig 3 style) for valves with at least one single-echo VARC-3 label
    w = out["_wide"]; st = w["stages"]["varc3"]; t = w["t"]
    lab2 = np.nan_to_num(st >= 2, nan=False); rows = np.where(lab2.any(1))[0][:60]
    if len(rows):
        yrs = np.arange(1, 11); M = np.full((len(rows), len(yrs)), np.nan)
        for i, r in enumerate(rows):
            for j, y in enumerate(yrs):
                k = np.where(~np.isnan(t[r]) & (np.abs(t[r] - y) < 0.5))[0]
                if len(k): M[i, j] = st[r, k[0]]
        fig, ax = plt.subplots(figsize=(6, 7)); im = ax.imshow(np.where(np.isnan(M), -1, M), aspect="auto", cmap="viridis", vmin=-1, vmax=3)
        ax.set_xticks(range(len(yrs))); ax.set_xticklabels(yrs); ax.set_xlabel("year"); ax.set_ylabel("valves with >=1 VARC-3 label"); ax.set_title("VARC-3 stage per annual echo (-1 = no echo)")
        fig.colorbar(im, ax=ax, fraction=0.03); fig.tight_layout(); fig.savefig(FIG / "label_flicker_heatmap.png", dpi=130); plt.close(fig)


def run_all(p, n=10000, seed=1, n_matched=20000):
    OUT.mkdir(exist_ok=True); rep = Report()
    if not p.physics.c_model: fit_physics_constants(p)
    out = simulate_cohort(p, n=n, seed=seed)
    check_reference_echo(p, rep, n_matched, seed)
    model_tab = check_models_vs_ase(out["patients_valves"], rep)
    val = check_velders(p, rep, n_matched, seed)
    curves = check_curves(p, rep, n_matched, seed)
    check_death_reint(p, rep, n_matched, seed)
    try:
        check_emergent_hr(out, rep)
    except Exception as e:  # lifelines failures should not kill the report
        rep.add("emergent HR", "Cox fit", "", "", "", False, str(e)[:120])
    descriptives(out, rep)
    figures(curves, out, val, model_tab)
    (OUT / "validation_report.md").write_text(rep.markdown()); print(rep.markdown())
    return out, rep


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default=str(OUT / "params_calibrated.json"))
    ap.add_argument("--n", type=int, default=10000); ap.add_argument("--seed", type=int, default=1); ap.add_argument("--n-matched", type=int, default=20000)
    a = ap.parse_args()
    p = SimParams.from_json(a.params) if Path(a.params).exists() else default_params()
    run_all(p, a.n, a.seed, a.n_matched)
