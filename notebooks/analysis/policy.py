"""Idea 2: surveillance policies on the simulator with feedback (C12-C15 of docs/05).

A sequential visit generator replaces the simulator's cohort-wide schedule: for every echo ordinal k it asks the policy
for each valve's next echo time, applies attendance, forces symptom-triggered and confirmation echoes, observes with a
PRE-DRAWN noise bank (common random numbers across policies), computes the online features with the SAME function as the
training table, and finally hands the visit matrices to the simulator's post_visits (labels, reintervention conditional on
detection, truncation). Metrics compare detection delay (vs the latent stage-2 crossing), echoes per patient-decade,
stage-3-before-stage-2 share, urgent reinterventions and false-positive-driven reinterventions."""
from __future__ import annotations

import json
import time
from dataclasses import dataclass

import joblib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed, parallel_config

import analysis_config as C
from data import Hist, baseline_frame
from features import history_features
from params import SimParams
from physics import fit_physics_constants
from simulate import MEASUREMENT_NAMES, post_visits, pre_visits
from visits import VTYPE_OTHER, VTYPE_REF, VTYPE_SCHED, VTYPE_SYMPTOM, draw_noise_bank, miss_prob, observe
from latent import latent_stage

OUT = C.OUT; TAB = OUT / "tables"; FIG = OUT / "figures"; MOD = OUT / "models"
VTYPE_CONFIRM = 4
K_MAX = 48
CONFIRM_DELAY = (0.25, 0.5)        # years after any single positive echo (VARC-3 / EAPCI: repeat at 3-6 months)
POST_CONFIRM_INTERVAL = 0.5        # after a confirmed label every policy echoes 6-monthly
MISSED_RESCHEDULE = 1.0            # a missed scheduled echo is re-offered a year later
STAGE3_CAP = 0.02                  # safety constraint of the risk-adapted policy


# ----------------------------------------------------------------------------- policies
class Policy:
    name = "base"
    confirmation = True

    def next_time(self, t_now, feats, idx):
        """t_now: echo times of the attending valves; feats: their online features (None for grid policies at start);
        idx: their row indices in the cohort (None = all valves)."""
        raise NotImplementedError


class Guideline(Policy):
    """30 d, 1 y, then annually; no confirmation rule (the pure standard of care)."""
    name = "guideline"; confirmation = False

    def next_time(self, t_now, feats, idx):
        return np.floor(t_now + 1e-6) + 1.0


class GuidelinePlusConfirmation(Guideline):
    name = "guideline_confirm"; confirmation = True


class Deescalation(Policy):
    """Annual to 5 y, then every 2 y while nothing has been seen (no positive echo so far and gradient < 20)."""
    name = "deescalation"

    def next_time(self, t_now, feats, idx):
        base = np.floor(t_now + 1e-6)
        low = (feats.prior_positive_count.to_numpy() == 0) & (feats.unconfirmed_flag.to_numpy() == 0) & (feats.mpg_now.to_numpy() < 20)
        return np.where((base >= 5) & low, base + 2.0, base + 1.0)


class TimeSinceImplantOnly(Policy):
    """Null policy: interval from the time-only hazard model with the same tau rule (no echo information)."""
    name = "time_only"

    def __init__(self, model, tau, candidates=C.INTERVAL_CANDIDATES_Y):
        self.model = model; self.tau = tau; self.cand = np.array(candidates)

    def next_time(self, t_now, feats, idx):
        cif = self.model.predict_cif(feats, self.cand)
        ok = cif <= self.tau
        choice = np.where(ok.any(1), self.cand[np.maximum((ok * np.arange(1, len(self.cand) + 1)).max(1) - 1, 0)], self.cand[0])
        return t_now + choice


class RiskAdapted(Policy):
    """Largest candidate interval whose predicted confirmed-stage->=2 risk is <= tau and predicted stage-3 risk <= STAGE3_CAP."""
    name = "risk_adapted"

    def __init__(self, model, model3, tau, candidates=C.INTERVAL_CANDIDATES_Y, cap3=STAGE3_CAP):
        self.model = model; self.model3 = model3; self.tau = tau; self.cand = np.array(candidates); self.cap3 = cap3
        self.name = f"risk_adapted_tau{tau:g}"

    def next_time(self, t_now, feats, idx):
        cif = self.model.predict_cif(feats, self.cand)
        ok = cif <= self.tau
        if self.model3 is not None:
            ok &= self.model3.predict_cif(feats, self.cand) <= self.cap3
        best = np.where(ok.any(1), self.cand[np.maximum((ok * np.arange(1, len(self.cand) + 1)).max(1) - 1, 0)], self.cand[0])
        return t_now + best


class Oracle(Policy):
    """Ceiling: no surveillance echo until the latent stage-2 crossing, then 3-monthly until confirmed."""
    name = "oracle"

    def __init__(self, t2_latent):
        self.t2 = t2_latent

    def next_time(self, t_now, feats, idx):
        t2 = self.t2 if idx is None else self.t2[idx]
        return np.where(np.isfinite(t2), np.where(t_now < t2, t2, t_now + 0.25), np.inf)


# ----------------------------------------------------------------------------- generator
@dataclass
class Bank:
    noise: dict
    attend: np.ndarray
    jitter: np.ndarray
    confirm_delay: np.ndarray
    reint: tuple


def make_bank(n, p, seed, K=K_MAX):
    rng = np.random.default_rng(seed)
    noise = draw_noise_bank(n, K, p, rng)
    return Bank(noise=noise, attend=rng.random((n, K)), jitter=rng.normal(0.0, p.visits.jitter_sd_y, (n, K)),
                confirm_delay=rng.uniform(*CONFIRM_DELAY, (n, K)), reint=(rng.exponential(1.0, n), rng.exponential(1.0, n), rng.random(n)))


def _hist_from(cols, k, valve_ids):
    """Build a Hist from the growing (n, K) arrays for the first k post-reference echoes."""
    return Hist(valve_ids=valve_ids, t=cols["t"][:, :k], vtype=cols["vt"][:, :k], mpg=cols["mpg"][:, :k], vmax=cols["vmax"][:, :k], eoa=cols["eoa"][:, :k],
                dvi=cols["dvi"][:, :k], svi=cols["svi"][:, :k], lvef=cols["lvef"][:, :k], ar=cols["ar"][:, :k], n_echoes=np.full(len(valve_ids), k))


def run_policy(policy: Policy, pre: dict, p: SimParams, bank: Bank, pv_like: pd.DataFrame, reference_time=0.08):
    """Sequential visits under `policy`; returns the same tuple as simulate._core (df, t, vt, obs, ref, stages, times, t_end, ends, latent)."""
    df = pre["df"]; n = len(df); t_end0 = pre["t_end0"]; t_sym = pre["t_sym"]; t_bg = pre["t_bg"]
    ref_present = df.reference_present.to_numpy() == 1
    # column 0 = reference echo (NaN when absent), then post-reference echoes in order
    K = K_MAX + 1
    T = np.full((n, K), np.nan); VT = np.full((n, K), np.nan)
    T[:, 0] = np.where(ref_present, reference_time, np.nan); VT[:, 0] = np.where(ref_present, VTYPE_REF, np.nan)
    obs_cols = {k: np.full((n, K), np.nan) for k in ("mpg", "eoa", "dvi", "vmax", "svi", "lvef", "ar", "pvl", "mpg_true", "eoa_true", "dvi_true", "ar_true")}
    o0 = observe(T[:, :1], df, p, None, vt=VT[:, :1], noise={k: v[:, :1] for k, v in bank.noise.items()})
    for k in obs_cols: obs_cols[k][:, 0] = o0[k][:, 0]
    base = baseline_frame(pv_like)   # reference values as the policy sees them (observed or table)
    base["ref_mpg"] = np.where(ref_present, o0["mpg"][:, 0], df.table_mpg.to_numpy()); base["ref_eoa"] = np.where(ref_present & ~np.isnan(o0["eoa"][:, 0]), o0["eoa"][:, 0], df.table_eoa.to_numpy())
    base["ref_dvi"] = np.where(ref_present & ~np.isnan(o0["dvi"][:, 0]), o0["dvi"][:, 0], df.table_eoa.to_numpy() / df.a_eff.to_numpy()); base["ref_ar"] = np.where(ref_present, o0["ar"][:, 0], 0.0)
    base["ref_eoai"] = base.ref_eoa / df.bsa.to_numpy(); base["eoa_residual"] = base.ref_eoa - df.table_eoa.to_numpy(); base["ref_mpg_flow_corrected"] = base.ref_mpg / df.table_mpg.to_numpy()
    hist_cols = {k: np.full((n, K_MAX), np.nan) for k in ("t", "mpg", "eoa", "dvi", "vmax", "svi", "lvef", "ar")}; hist_cols["vt"] = np.full((n, K_MAX), -1, int)
    # per-valve state
    t_last = np.full(n, reference_time); next_due = policy.next_time(np.full(n, reference_time), None, None) if not _needs_feats(policy) else np.full(n, reference_time + 1.0)
    next_due = np.asarray(next_due, float) + bank.jitter[:, 0]
    sym_done = np.zeros(n, bool); bg_done = np.zeros(n, bool); confirm_due = np.full(n, np.inf); confirmed = np.zeros(n, bool); consec_pos = np.zeros(n, int)
    active = np.ones(n, bool); k = 0
    while active.any() and k < K_MAX:
        col = k + 1
        t_sym_echo = np.where(~sym_done & np.isfinite(t_sym), t_sym + p.visits.symptom_echo_delay_y, np.inf)
        t_bg_echo = np.where(~bg_done & np.isfinite(t_bg) & (t_bg > reference_time + 0.2), t_bg, np.inf)
        cand = np.minimum.reduce([next_due, t_sym_echo, t_bg_echo, confirm_due])
        kind = np.select([cand == t_sym_echo, cand == confirm_due, cand == t_bg_echo], [VTYPE_SYMPTOM, VTYPE_CONFIRM, VTYPE_OTHER], VTYPE_SCHED)
        due = active & np.isfinite(cand) & (cand < t_end0) & (cand > reference_time + 0.2)
        active &= np.isfinite(cand) & (cand < t_end0)
        # attendance for policy-scheduled echoes (year-indexed missingness, well / symptomatic multipliers); forced echoes always happen
        yr = np.maximum(1, np.round(np.where(np.isfinite(cand), cand, 1.0))).astype(int)
        pm = miss_prob(yr, p)
        st_lat = latent_stage(np.where(np.isfinite(cand), cand, np.nan)[:, None], df, p)[0][:, 0]
        well = (st_lat == 0) & (cand < np.where(np.isfinite(t_sym), t_sym, np.inf)); symptomatic = cand >= np.where(np.isfinite(t_sym), t_sym, np.inf)
        pm = np.clip(pm * np.where(well, p.visits.miss_mult_well, 1.0) * np.where(symptomatic, p.visits.miss_mult_symptomatic, 1.0), 0, 0.95)
        attend = due & ((kind != VTYPE_SCHED) | (bank.attend[:, k] >= pm))
        missed = due & ~attend
        if not attend.any():
            next_due = np.where(missed & (kind == VTYPE_SCHED), cand + MISSED_RESCHEDULE, np.where(missed, cand + 0.05, next_due))
            k += 1
            continue
        # observe the attended echoes
        tcol = np.where(attend, cand, np.nan)
        o = observe(tcol[:, None], df, p, None, vt=np.where(attend, kind, np.nan)[:, None].astype(float), noise={kk: v[:, col:col + 1] for kk, v in bank.noise.items()})
        T[:, col] = tcol; VT[:, col] = np.where(attend, kind, np.nan)
        for kk in obs_cols: obs_cols[kk][:, col] = o[kk][:, 0]
        for kk in ("mpg", "eoa", "dvi", "vmax", "svi", "lvef", "ar"): hist_cols[kk][:, k] = o[kk][:, 0]
        hist_cols["t"][:, k] = tcol; hist_cols["vt"][:, k] = np.where(attend, kind, -1).astype(int)
        # bookkeeping of forced echoes
        sym_done |= attend & (kind == VTYPE_SYMPTOM); bg_done |= attend & (kind == VTYPE_OTHER)
        confirm_due = np.where(attend & (kind == VTYPE_CONFIRM), np.inf, confirm_due)
        # online features for the attended valves (history = echoes 1..k+1 for valves with that many, NaN-padded otherwise)
        feats = _online_features(hist_cols, k + 1, base, df.valve_id.to_numpy(), attend)
        # confirmation logic on the observed single-echo stage (haemodynamic or full VARC-3 >= 2)
        pos = np.zeros(n, bool); pos[attend] = feats.unconfirmed_flag.to_numpy() == 1
        consec_pos = np.where(attend, np.where(pos, consec_pos + 1, 0), consec_pos)
        newly_confirmed = attend & (consec_pos >= 2)
        confirmed |= newly_confirmed
        if policy.confirmation:
            need = attend & pos & ~confirmed & (kind != VTYPE_CONFIRM)
            confirm_due = np.where(need, cand + bank.confirm_delay[:, k], confirm_due)
        # next policy-scheduled echo
        nxt = np.full(n, np.inf)
        nxt[attend] = np.asarray(policy.next_time(cand[attend], feats, np.where(attend)[0]), float)
        nxt = np.where(attend & confirmed, np.minimum(nxt, cand + POST_CONFIRM_INTERVAL), nxt)
        nxt = nxt + np.where(attend, bank.jitter[:, col] * (kind == VTYPE_SCHED), 0.0)
        next_due = np.where(attend, np.maximum(nxt, cand + 0.2), np.where(missed & (kind == VTYPE_SCHED), cand + MISSED_RESCHEDULE, next_due))
        next_due = np.where(missed & (kind != VTYPE_SCHED), cand + 0.05, next_due)   # forced echo blocked by bounds: retry shortly
        t_last = np.where(attend, cand, t_last)
        k += 1
    # sort columns by time, NaN last, then hand to the simulator's post-visit half
    order = np.argsort(np.where(np.isnan(T), np.inf, T), axis=1)
    T = np.take_along_axis(T, order, 1); VT = np.take_along_axis(VT, order, 1)
    for kk in obs_cols: obs_cols[kk] = np.take_along_axis(obs_cols[kk], order, 1)
    return post_visits(pre, T, VT, obs_cols, p, np.random.default_rng(0), reference_time=reference_time, reint_draws=bank.reint)


def _needs_feats(policy):
    return isinstance(policy, (RiskAdapted, TimeSinceImplantOnly, Deescalation))


def _online_features(hist_cols, k, base, valve_ids, rows_mask):
    rows = np.where(rows_mask)[0]
    h = _hist_from(hist_cols, k, valve_ids)
    # valves attending at ordinal k have exactly k post-reference echoes recorded in columns 0..k-1 (NaN where missed);
    # history_features tolerates NaN columns (missed echoes) because every reduction is NaN-aware
    f = history_features(h, k, base, rows)
    f = f.copy(); f["n_echoes"] = np.sum(~np.isnan(hist_cols["t"][rows, :k]), 1)
    return f


# ----------------------------------------------------------------------------- metrics
def policy_metrics(res, pre, name):
    np.seterr(invalid="ignore")
    df, t, vt, obs, ref, stages, times, t_end, ends, latent = res
    n = len(df); fu = t_end
    echoes = np.sum(~np.isnan(t) & (vt != VTYPE_REF), 1)
    t2 = pre["t2_latent"]; t3 = pre["t3_latent"]
    conf = times["varc3"]["conf2_confirming"]; first = times["varc3"]["single2"]
    has_cross = np.isfinite(t2) & (t2 < t_end)
    detected = has_cross & np.isfinite(conf)
    delay = np.maximum(conf - t2, 0.0)[detected]; delay_first = np.maximum(first - t2, 0.0)[has_cross & np.isfinite(first)]
    stage3 = np.isfinite(t3) & (t3 < t_end)
    seen_as_2 = stage3 & np.isfinite(first) & (first < t3)
    re = np.isfinite(ends["t_re"]); urgent = re & ends["re_urgent"]
    fp_re = re & (ends["re_indication"] != "") & ~(np.isfinite(t2) & (t2 <= ends["t_re"]))
    return {"policy": name, "echoes_per_patient_decade": float(echoes.sum() / fu.sum() * 10), "mean_followup_y": float(fu.mean()),
            "latent_stage2_within_followup": int(has_cross.sum()), "detected_confirmed": int(detected.sum()), "detection_rate": float(detected.mean() / max(has_cross.mean(), 1e-9)),
            "confirmed_delay_mean_y": float(delay.mean()) if len(delay) else np.nan, "confirmed_delay_p90_y": float(np.percentile(delay, 90)) if len(delay) else np.nan,
            "confirmed_delay_median_y": float(np.median(delay)) if len(delay) else np.nan, "first_flag_delay_mean_y": float(delay_first.mean()) if len(delay_first) else np.nan,
            "stage3_cases": int(stage3.sum()), "stage3_seen_while_stage2_share": float(seen_as_2.sum() / max(stage3.sum(), 1)),
            "reinterventions": int(re.sum()), "urgent_share": float(urgent.sum() / max(re.sum(), 1)), "false_positive_reinterventions": int(fp_re.sum()),
            "deaths_per_10k_illustrative": float(10000 / n * (urgent.sum() * 0.226 + (re.sum() - urgent.sum()) * 0.014))}


def _load_models():
    m = joblib.load(MOD / "dth_policy.joblib") if (MOD / "dth_policy.joblib").exists() else joblib.load(MOD / "dth_full.joblib")
    m3 = joblib.load(MOD / "dth_stage3.joblib") if (MOD / "dth_stage3.joblib").exists() else None
    m0 = joblib.load(MOD / "dth_time_only.joblib")
    return m, m3, m0


def _fit_policy_models():
    """Spacing-robust landmark model for the policy and a stage-3 model for the safety constraint (fitted on the training valves)."""
    from models import DiscreteTimeHazard
    from outcomes import make_instances
    from data import load_tables
    T = load_tables(); L = pd.read_parquet(OUT / "landmark_table.parquet"); sp = {k: set(v) for k, v in json.load(open(OUT / "split.json")).items()}
    part = lambda d, k: d[d.valve_id.isin(sp[k])].reset_index(drop=True)
    inst = make_instances(L, T["labels"], C.PRIMARY_LABEL, T["events"]); inst3 = make_instances(L, T["labels"], "varc3_single3", T["events"])
    m = DiscreteTimeHazard(C.SPACING_ROBUST).fit(part(inst, "train"), part(inst, "val")); joblib.dump(m, MOD / "dth_policy.joblib")
    m3 = DiscreteTimeHazard(C.SPACING_ROBUST).fit(part(inst3, "train"), part(inst3, "val")); joblib.dump(m3, MOD / "dth_stage3.joblib")
    return m, m3


def _one_policy(name, tau, n, seed):
    """Worker: regenerate the pre-visit draws (same seed for every policy), run the policy, return metrics."""
    from threadpoolctl import threadpool_limits
    with threadpool_limits(1):
        p = SimParams.from_json(C.PARAMS_CALIBRATED); fit_physics_constants(p)
        rng = np.random.default_rng(seed); pre = pre_visits(p, n, rng)
        bank = make_bank(n, p, seed + 7)
        m, m3, m0 = _load_models()
        pol = {"guideline": Guideline(), "guideline_confirm": GuidelinePlusConfirmation(), "deescalation": Deescalation(),
               "time_only": TimeSinceImplantOnly(m0, tau if tau else 0.05), "oracle": Oracle(pre["t2_latent"]),
               "risk_adapted": RiskAdapted(m, m3, tau if tau else 0.05)}[name]
        pv_like = _pv_like(pre["df"])
        res = run_policy(pol, pre, p, bank, pv_like)
        met = policy_metrics(res, pre, pol.name); met["tau"] = tau
        return met


def _pv_like(df):
    """A patients_valves-like frame (the columns baseline_frame needs) from the simulator's internal df."""
    pv = df.copy()
    pv["reference_source"] = np.where(df.reference_present == 1, "30d_echo", "model_size_table")
    pv["ref_mean_gradient_mmHg"] = df.table_mpg; pv["ref_av_area_cm2"] = df.table_eoa; pv["ref_dvi"] = df.table_eoa / df.a_eff
    pv["ref_ar_grade"] = 0.0; pv["ref_eoai_cm2_m2"] = df.table_eoa / df.bsa; pv["ppm_class_ref"] = df.ppm_true; pv["followup_end_years"] = np.nan
    return pv


def run(n=10000, seed=C.SEED + 11, taus=(0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15)):
    t0 = time.time()
    if not (MOD / "dth_policy.joblib").exists():
        _fit_policy_models(); print(f"[policy] policy models fitted ({time.time()-t0:.0f}s)", flush=True)
    tasks = [("guideline", None), ("guideline_confirm", None), ("deescalation", None), ("oracle", None), ("time_only", 0.05)] + [("risk_adapted", tau) for tau in taus]
    with parallel_config(backend="loky", n_jobs=min(len(tasks), C.N_JOBS), inner_max_num_threads=1):
        rows = Parallel()(delayed(_one_policy)(name, tau, n, seed) for name, tau in tasks)
    res = pd.DataFrame(rows); res.to_csv(TAB / "policy_metrics.csv", index=False)
    g = res[res.policy == "guideline"].iloc[0]
    res["echo_reduction_vs_guideline_pct"] = 100 * (1 - res.echoes_per_patient_decade / g.echoes_per_patient_decade)
    res["extra_delay_mean_months"] = 12 * (res.confirmed_delay_mean_y - g.confirmed_delay_mean_y)
    res["extra_delay_p90_months"] = 12 * (res.confirmed_delay_p90_y - g.confirmed_delay_p90_y)
    res["passes_section0"] = (res.echo_reduction_vs_guideline_pct >= 25) & (res.extra_delay_mean_months <= 3) & (res.extra_delay_p90_months <= 3)
    # the fair baseline for the confirmation-including policies is the guideline WITH the confirmation echo
    gc = res[res.policy == "guideline_confirm"].iloc[0]
    res["echo_reduction_vs_guideline_confirm_pct"] = 100 * (1 - res.echoes_per_patient_decade / gc.echoes_per_patient_decade)
    res["extra_delay_vs_confirm_mean_months"] = 12 * (res.confirmed_delay_mean_y - gc.confirmed_delay_mean_y)
    res["extra_delay_vs_confirm_p90_months"] = 12 * (res.confirmed_delay_p90_y - gc.confirmed_delay_p90_y)
    res["passes_section0_vs_confirm"] = (res.echo_reduction_vs_guideline_confirm_pct >= 25) & (res.extra_delay_vs_confirm_mean_months <= 3) & (res.extra_delay_vs_confirm_p90_months <= 3)
    res.to_csv(TAB / "policy_metrics.csv", index=False)
    json.dump(res.to_dict("records"), open(TAB / "policy_metrics.json", "w"), indent=1, default=float)
    _figure(res)
    print(res[["policy", "tau", "echoes_per_patient_decade", "confirmed_delay_mean_y", "confirmed_delay_p90_y", "detection_rate", "urgent_share", "false_positive_reinterventions", "echo_reduction_vs_guideline_pct", "extra_delay_mean_months", "passes_section0", "echo_reduction_vs_guideline_confirm_pct", "extra_delay_vs_confirm_mean_months", "extra_delay_vs_confirm_p90_months", "passes_section0_vs_confirm"]].round(3).to_string(index=False), flush=True)
    print(f"[policy] done in {time.time()-t0:.0f}s", flush=True)
    return res


def _figure(res):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ra = res[res.policy.str.startswith("risk_adapted")].sort_values("tau")
    ax.plot(ra.echoes_per_patient_decade, 12 * ra.confirmed_delay_mean_y, "o-", color="C1", label="risk-adapted (tau sweep)")
    for _, r in ra.iterrows(): ax.annotate(f"{r.tau:g}", (r.echoes_per_patient_decade, 12 * r.confirmed_delay_mean_y), fontsize=7, xytext=(3, 3), textcoords="offset points")
    for name, mk, col in (("guideline", "s", "k"), ("guideline_confirm", "D", "C0"), ("deescalation", "^", "C2"), ("time_only", "v", "C3"), ("oracle", "*", "C4")):
        r = res[res.policy == name]
        if len(r): ax.plot(r.echoes_per_patient_decade, 12 * r.confirmed_delay_mean_y, mk, color=col, ms=9, label=name)
    ax.set_xlabel("echoes per patient-decade"); ax.set_ylabel("mean confirmed-detection delay (months)"); ax.grid(alpha=.3); ax.legend(fontsize=8)
    ax.set_title("Surveillance policies: echoes vs delay (common random numbers)", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "policy_frontier.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    run()
