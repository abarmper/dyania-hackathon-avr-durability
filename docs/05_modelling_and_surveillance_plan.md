# Modelling (Idea 1) and surveillance policy (Idea 2): approach and design choices

Companion to `03_ideas_brainstorm.md` (§0 endpoints, Ideas 1–2) and `04_synthetic_cohort_spec.md` (the cohort). The cohort exists (`data/synthetic/`, 10,000 valves, 50,849 echoes, 1,641 confirmed primary events of which 597 / 1,003 / 1,322 by 5 / 8 / 10 y, 4,009 deaths, 284 competing reinterventions, 1,007 "new model" valves, 31% without a 30-day reference echo). This document fixes *how* it is modelled and how a risk-adapted echo schedule is tested, and lists the design choices still open. Recommended options come first and are marked **(rec.)**; each has a one-line cost. Section 6 collects the decisions.

---

## 0. What the models must deliver (from §0) and where it goes

| Requirement (§0) | Value / definition | Produced by | Fills |
|---|---|---|---|
| Primary endpoint | confirmed VARC-3 stage ≥2 (two consecutive echoes), competing death and non-SVD reintervention; cumulative incidence at 5 / 8 / 10 y | implant-time model + landmark model (§2) | protocol §3, approach §1 |
| Usefulness threshold, statistical | time-dependent AUC ≥0.75 and C ≥0.70 at 5 y; calibration slope 0.8–1.2 at 5 / 8 / 10 y, under the confirmed label | evaluation (C7) | protocol §5 metrics |
| Usefulness threshold, decision-level | ≥25% fewer echoes per patient-decade with detection delay ≤3 months longer than annual surveillance, or same echo count with earlier detection | Idea 2 (§3) | protocol §3, approach §5–§6 |
| Secondary objective 1 | SAVR vs TAVR and valve-model trajectories | subgroup metrics, class-level CIF, held-out new-valve test (C6) | protocol §5 subgroups |
| Secondary objective 2 | risk-adapted vs annual surveillance | Idea 2 | approach §6 |
| Secondary objective 3 | reintervention window: time from first moderate deterioration to severe / failure; urgent share | C10 + Idea 2 metrics | protocol §3 secondary endpoints |
| Secondary objective 4 | reliability of the SVD label across definitions | C3 label sweep + C9 audit against truth | protocol §4 ground truth, approach §7 |
| Per-inference output | CIF at 5 / 8 / 10 y; risk tier; top-3 features; next-echo interval; unconfirmed-echo flag | C2 (CIF), C13 (tier and interval), C8 (top-3), C4 (flag) | approach §5 |
| Subgroups | approach, valve class, age band (<60, 60–80, >80), size and EOAi, sex | C7 | protocol §5 |
| Comparators | time since implant alone; published risk-factor Cox (`B4` Table 2 set); STS-PROM for death only; gestalt not testable | C2 comparators | protocol §5 comparator |
| Fairness | sex simulated with no direct effect | C7 | protocol §6 |
| Reporting | TRIPOD+AI (`D8`) | C8 | protocol §5, approach §4 |

---

## 1. The data as the models see it, and the rules that keep it honest

**Prediction instances.** One instance per (valve, echo) at every post-reference echo (the *landmark*), plus one instance per valve at the reference echo (the implant-time model). Features use only echoes at or before the landmark. `latent_truth` is never an input; it is used only for the oracle label, the achievable-AUC ceiling and the detection-delay metrics. The reference echo appears legitimately on both sides (anchor of the predictor deltas and of the label's change-from-baseline); that is shared conditioning, not leakage.

**C0. Outcome timing (fixed, not a choice).** `labels.primary_time` is the *first* echo of the confirmed pair, but the label only becomes known at the *confirming* echo (`varc3_conf2_confirming`, median 0.94 y later, 90th percentile 2.1 y). A landmark model must therefore use the confirming time: the risk set at landmark L is the valves with `varc3_conf2_confirming > L`, and the event time is the confirming echo ("when will we know"). Using the first-of-pair time would either drop the hardest instances from the risk set or create negative times to event. The implant-time 5 / 8 / 10-y cumulative incidence keeps the first-of-pair convention (the trial convention), and the report states the difference in one sentence.

**Two instance sets are always reported.** All landmark instances, and the *currently negative* ones (the echo itself does not meet a single-echo stage ≥2 rule). Eleven percent of landmark echoes are already single-echo positive; predicting "confirmed at the next echo" for them is near-tautological and inflates the AUC. Surveillance decisions are taken on the currently-negative set.

**The symptom-triggered flag.** In the simulator, symptoms only occur once the latent stage is ≥2, so `visit_type == symptom_triggered` is stronger than in reality. It is a legitimate real-world field (indication for echo), so it stays, but every model is reported with and without it.

---

## 2. Idea 1: design choices

### C1. Prediction framing
- **(a, rec.) Landmark dynamic prediction plus an implant-time model.** Each echo is a decision point; the model predicts the cumulative incidence over the next 3 and 5 years from everything known at that echo. The implant-time model gives the 5 / 8 / 10-y CIF the protocol asks for. Cost: a wide landmark table (~44,000 instances) and online feature computation; half a day.
- (b) Implant-time model only. Mirrors published risk-factor models and is a required comparator, but cannot personalise the echo interval after year 1.
- (c) Joint longitudinal–survival model (`D1`). The principled dynamic method and the source of the information-gain scheduling rule; no mature Python implementation and not needed to demonstrate the point. Cited as the method for the real study.

### C2. Competing-risk model
- **(a, rec.) Discrete-time cause-specific hazard model.** Person-period expansion of each landmark instance into 6-month periods; a multinomial `HistGradientBoostingClassifier` (scikit-learn, installed) over {SVD label, death or other reintervention, none} with the period index as a feature. Censoring is handled by the expansion (loss to follow-up is ~1%/y and stated; there is no administrative censoring before 5.7 y because implants end in 2020 and follow-up in 2026), so no inverse-probability weights are needed. One fit yields P(event in (L, L + Δ]) for any Δ: the 5 / 8 / 10-y CIF for Idea 1 and the interval risk Idea 2 needs. Native handling of missing values and categorical valve model. About 250,000 rows, seconds to fit. Cost: the expansion code and a clear explanation for clinicians.
- (b) Cause-specific Cox + Aalen–Johansen (lifelines, installed). The interpretable statistical model, reported with a hazard-ratio table per `D5`. Built alongside (a).
- (c) Fixed-horizon multinomial classifier on the eligible subset (5 y: all valves; 8 y: 7,937; 10 y: 6,172). Simplest fallback if (a) misbehaves.
- (d) Fine–Gray. No maintained Python implementation; the subdistribution view is reported through (b) with the cumulative incidence from Aalen–Johansen.
- (e) Random survival forest and DeepHit. scikit-survival is not installed and `docs/03` already rules out deep survival models at this scale; both are cited as next steps at registry volume.
- Note on interpretation: deaths (2,256 by 5 y) outnumber SVD events (597 by 5 y), so the multinomial is partly a mortality model; SVD discrimination is reported cause-specifically (C7). Death and SVD are independent in the simulator (no SVD-to-mortality pathway), which is stated wherever mortality appears.

### C3. Which label to train and evaluate under
Train on the confirmed VARC-3 label (the primary endpoint). Evaluate every model against single-echo VARC-3, the NOTION haemodynamic variant, Capodanno, Dvir and the oracle, with the horizon and the eligible set held fixed and only the outcome indicator swapped, so that risk sets do not change between labels. The headline figure is discrimination as a function of label noise: single-echo → confirmed → oracle. Choice: (a, rec.) all five plus oracle; (b) confirmed and oracle only if time is short.

### C4. Feature engineering (the edge)
Every derived feature is framed as an ASE 2024 (`A5`) or VARC-3 recommendation with its citation, not as knowledge of the generator, and a **registry-only** feature set (what TVT and STS actually hold) is shipped as the deployment-realism comparator.

| Group | Features | Rationale |
|---|---|---|
| Baseline | age, sex, BSA, BMI, diabetes, CKD/dialysis, smoking, AF and anticoagulation type, approach, valve model with class parent, `is_new_model`, labelled size, implant year, `reference_source` | `docs/04` §7 audit; every field exists on TVT v3.0 / STS v4.20.2 |
| Reference echo | MPG, EOA, EOAi, DVI, AR grade, PPM class; **expected EOA for model × size** from the ASE table and the **residual** (measured − expected); flow-corrected gradient (MPG normalised by stroke-volume index) | ASE 2024 Table 5 compares measured EOA with the reference EOA for the model and size; DVI is the flow-independent index |
| Longitudinal, at the landmark | last MPG and change from reference; **log-scale exponentially weighted mean of MPG** (Doppler error is multiplicative, ~14% per echo) and the **posterior probability that the true gradient is ≥20 mmHg**; per-year slope since reference (NaN with fewer than two post-reference echoes); slope over the last two echoes; max MPG so far; relative change in EOA and DVI; **current single-echo stage and count of prior positive flags** (this is the unconfirmed-echo flag of the output format); consecutive-worsening count; time since last echo; number of echoes; symptom-triggered flag (reported with/without); transient-spike flag (rise then fall) | VARC-3 requires serial echoes; the two-echo rule and Velders' flicker motivate the smoothed and count features |
| Missingness | indicators for missing EOA/DVI channels and for the table-based reference; native NaN handling in the boosted model | Velders Table E1; the table-reference arm has a different noise structure (persistent LVOT bias no longer cancels in deltas) |

Choice: **(a, rec.)** full set plus registry-only set, both reported; (b) registry-only set alone (faster, loses the edge).

### C5. Feature selection
- **(a, rec.)** Domain prefilter (the `docs/04` §7 audit already removed variables without a quantified effect) → permutation importance on the held-out fold → greedy backward elimination to ≤8 features that keep ≥95% of the full model's AUC. The result is a hand-computable "parsimonious model" for the protocol. Cost: ~30 min.
- (b) Stability selection (bootstrap × L1-penalised Cox). Half a day for a slide nobody asks for; dropped.

### C6. Valve-model pooling (patch P4)
- **(a, rec.)** `valve_model` and `valve_class` as native categorical features of the boosted model; the deliverable is the **held-out new-valve test**: train without the 1,007 `is_new_model == 1` valves, evaluate on them with the model name masked to its class. Shows what the system does for a device with no history.
- (b) Empirical-Bayes target encoding with shrinkage to the class. Unnecessary with native categoricals.

### C7. Validation and metrics
- **(a, rec.)** Split by valve (70 / 15 / 15) with `GroupKFold` for tuning; **cluster bootstrap by valve** for every confidence interval (44,000 instances come from 10,000 valves). Temporal split (implants ≤2015 vs ≥2016 at the 5-y horizon) as a sensitivity analysis only, because it is confounded with device generation and the VKA-to-DOAC switch. **Case-mix transport**: evaluate on one freshly simulated paper-matched cohort (NOTION-like, age 79) — called robustness, not external validation, because the generative mechanism is the same. **Mechanism shifts only a simulator affords**: echo noise ×0.5 and ×1.5, missingness ×2, progression speed ×0.7.
- (b) Temporal split only. Cheaper, weaker.
- Metrics: cause-specific time-dependent AUC at 5 / 8 / 10 y (competing events count as controls, Blanche 2013); cause-specific C (competing events ranked as lower-risk non-events; on the eligible subset this is a truncated Harrell's C and is named as such rather than "Uno's C"); Brier score; decile calibration against Aalen–Johansen on the same eligible subset with a logistic recalibration slope; all with bootstrap CIs; subgroups by approach, class, age band, size, sex.
- Fairness: sex has **no direct effect** in the simulator but acts through labelled size and PPM (severe PPM 8% in women vs 12% in men), so the claim to make is "calibration and discrimination within sex", not "no association".

### C8. Explainability
- **(a, rec.)** SHAP `TreeExplainer` on the boosted model: global importance, dependence plots for the slope, the smoothed gradient and EOAi, a per-patient top-3 waterfall (the output format's "top-3 contributing features"), and the unconfirmed-echo flag surfaced next to the prediction. Requires one `pip install shap` (alternative: `xgboost`, whose `pred_contribs` gives exact TreeSHAP). Plus the Cox hazard-ratio table with CIs and the TRIPOD+AI checklist (`D8`).
- (b) scikit-learn permutation importance, partial dependence, and "feature → population median" contrasts for the per-patient top-3, if no package may be installed. Same story, less polished.

### C9. Oracle-only analyses (what no registry can do)
- **Achievable-AUC ceiling.** Train the same features on the *latent* stage-2 crossing and evaluate on the oracle label. The gap between the observed-label AUC and this ceiling decomposes into label noise, model error and feature insufficiency, and answers whether 0.75 is reachable at all.
- **Definition audit against truth.** Sensitivity of the confirmed VARC-3 label for a latent stage-2 crossing within follow-up (≈60% in the cohort), false-positive share (≈6% of confirmed labels precede any latent crossing), lead time, and kappa between the definitions.
- **Value of the 30-day reference echo.** AUC and detection delay in the 31% table-reference arm vs the echo-reference arm: a concrete registry recommendation.
- **Measurement-noise elasticity.** AUC and delay under echo noise ×0.5 / ×1.5: "improving echo reproducibility is worth X AUC".
- **Simulation-based sample size for the real study.** Width of the 5-y AUC confidence interval at n = 500 / 1,000 / 2,000 / 5,000 valves; the protocol's sample-size section.
- Choice: **(a, rec.)** all five (≈2 h total); (b) the first two only.

### C10. Secondary endpoints (descriptive)
Time from first confirmed stage 2 to stage 3 or reintervention (Kaplan–Meier by valve class) and the urgent share (objective 3); annualised gradient slope per patient from a mixed model (`statsmodels.MixedLM`, patch P2); BVF and reintervention cumulative incidence (P1). No dedicated model.

### C11. Uncertainty
Cluster-bootstrap confidence intervals on every metric (rec.). Per-prediction intervals from a bagged ensemble are skipped.

---

## 3. Idea 2: design choices

### C12. What quantity schedules the next echo
- **(a, rec.)** The predicted probability of a *confirmed stage ≥2* label within the candidate interval, from the landmark model, with the predicted probability of stage 3 within the interval as a safety constraint. This is the quantity the §0 criterion measures (delay to confirmed detection).
- (b) Stage-3 risk alone. Only 354 confirmed / 1,427 single-echo stage-3 events to learn from, and not what the criterion measures.

### C13. Policies compared (every policy includes the confirmation rule)
Any single positive echo triggers a repeat at 3–6 months regardless of risk (VARC-3 / EAPCI). Without it, sparse schedules are penalised by their own spacing on the confirmed-delay metric, since the two-echo rule alone costs a median 0.94 y under annual spacing.
1. Guideline: 30 days, 1 year, then annually (`A2`).
2. Guideline plus confirmation echo. Isolates what the confirmation rule itself buys.
3. De-escalation: annual to 5 y, then every two years if low risk (`D10`). The policy to beat: in years 1–5 a risk model will say "36 months" for nearly everyone, so most savings are de-escalation in disguise unless the increment is shown.
4. Risk-adapted: the next interval is the largest of {6, 12, 18, 24, 36} months whose predicted risk is ≤τ; τ is swept to draw the echoes-vs-delay frontier; the clinical risk tiers (low / moderate / high → 36 / 12 / 6 months) are read off τ.
5. Time-since-implant-only policy (no echo information): the null.
6. Oracle policy (echo at the latent crossing, then confirmation): the ceiling; the distance to it is the model's headroom.
The information-gain rule (`D1`) is cited, not built.

### C14. Simulator integration
- **(a, rec.)** Fork the simulator core after the pre-visit draws. Baseline, thrombosis, onset, death, endocarditis, loss to follow-up and symptom times are all drawn before the visits, so common random numbers across policies come for free. Replace the cohort-wide visit generation by a per-valve sequential generator (loop over visit index; observe one column at a time; update the online features; call the policy; apply the year-specific missingness; force the symptom-triggered echo). Reuse the labels, the reintervention module (conditional on *detected* status, so policies have consequences) and the truncation at reintervention. About 150 lines. Note: the simulator's `policy` argument is currently unused and the schedule is one cohort-wide list, so no shortcut exists.
- (b) Post-hoc thinning of the guideline visits. Cheap but ignores the feedback (detection → reintervention → end of follow-up); rejected.
- Distribution shift: the model is trained at annual spacing; under the policy spacing is 6–36 months, so slope, time-since-last-echo and worsening counts leave their training support. Options: restrict the policy model to spacing-robust features (per-year slope, smoothed gradient, max so far, current stage) **(rec.)**, or train it on a cohort generated with randomised intervals (off-policy training, only a simulator affords it) as a robustness run.

### C15. Metrics and the pass criterion
Echoes per patient-decade; first-flag delay and confirmed delay measured from `latent_truth.stage2_crossing_years` (clipped at zero because ~6% of confirmed labels precede the latent crossing), reported as mean and 90th percentile; share of stage-3 cases first seen while still stage 2; urgent-reintervention share; harm: reinterventions triggered by false-positive labels; expected deaths per 10,000 from the 1.4% elective vs 22.6% emergency redo mortality in `B3` (illustrative, since death is independent of SVD in the simulator); everything by approach and valve class. Pass criterion: the §0 rule (≥25% fewer echoes with ≤3 months extra delay, on the mean and the 90th percentile). The guideline baseline is already computable from the generated cohort and is the slide-2 number: about 6.1 echoes per patient-decade, median confirmed delay 1.4 y, urgent share 52%.

---

## 4. What gives the team an edge

1. Physics-informed and noise-aware features cited to ASE 2024 and VARC-3, with a registry-only comparator that shows their increment.
2. The label-noise decomposition: single-echo → confirmed → oracle → achievable ceiling.
3. The definition audit against truth (sensitivity, false-positive share, lead time of each consensus definition).
4. The value of the 30-day reference echo as a registry recommendation.
5. Robustness under mechanism shifts (noise, missingness, progression speed) and case-mix transport.
6. A feedback-aware surveillance experiment with a confirmation rule, a null policy and an oracle frontier, so the model's increment over de-escalation is visible.
7. A simulation-based sample size for the real study.
8. Provenance of every number and TRIPOD+AI reporting.

---

## 5. Implementation order (after the decisions; ≈12 h)

| Step | Content | Time |
|---|---|---|
| 1 | Landmark table from `echo_visits` + online features (`notebooks/analysis/features.py`) | 1.5 h |
| 2 | Outcomes per instance with the confirming time; eligible sets at 5 / 8 / 10 y (`outcomes.py`) | 1 h |
| 3 | Discrete-time hazard model, cause-specific Cox + Aalen–Johansen, the two comparators (`models.py`) | 2 h |
| 4 | Metrics, calibration, cluster bootstrap, label sweep (`evaluate.py`) — **minimum viable line** | 2 h |
| 5 | Explainability (`explain.py`) | 1 h |
| 6 | Sequential visit generator and the six policies (`policy.py`) | 3 h |
| 7 | Oracle-only analyses | 1.5 h |
| 8 | Figures, `02_features_and_models.ipynb`, `03_surveillance_policy.ipynb`, write-up into the templates | 2 h |

Environment: everything except SHAP runs on the installed stack (scikit-learn 1.8, lifelines 0.30, statsmodels); one `pip install shap` unlocks C8(a). scikit-survival is optional and only after a dry run against scikit-learn 1.8.

---

## 6. Decisions

All recommended options were adopted and implemented in `notebooks/analysis/` (results in `notebooks/analysis/output/results_summary.md`). Deviations found necessary during implementation: SVD-driven reinterventions that truncate follow-up before the confirming echo count as events (the 284 such valves all carry an SVD indication); categorical features are ordinal codes so that SHAP explains the model's own splits; the confirmation echo is added to every policy except the pure guideline, and the policy comparison is reported against both the pure guideline and the guideline-plus-confirmation baseline; cohort regenerations use the calibrated parameter file.

| # | Choice | Recommended | Decision |
|---|---|---|---|
| C0 | Outcome timing | confirming echo for landmark models; first-of-pair for the implant-time CIF (fixed) | — |
| C1 | Framing | landmark dynamic + implant-time model | implemented as recommended (Sept 17, `notebooks/analysis/`) |
| C2 | Competing-risk model | discrete-time cause-specific hazard (boosted multinomial) + Cox/AJ | implemented as recommended (Sept 17, `notebooks/analysis/`) |
| C3 | Labels | train on confirmed; evaluate on all five + oracle | implemented as recommended (Sept 17, `notebooks/analysis/`) |
| C4 | Features | full physics-informed set + registry-only comparator | implemented as recommended (Sept 17, `notebooks/analysis/`) |
| C5 | Feature selection | permutation importance + backward elimination to ≤8 | implemented as recommended (Sept 17, `notebooks/analysis/`) |
| C6 | Valve-model pooling | native categoricals + held-out new-valve test | implemented as recommended (Sept 17, `notebooks/analysis/`) |
| C7 | Validation | valve-level split, cluster bootstrap, temporal sensitivity, case-mix + mechanism shifts | implemented as recommended (Sept 17, `notebooks/analysis/`) |
| C8 | Explainability | SHAP (needs one install) | implemented as recommended (Sept 17, `notebooks/analysis/`) |
| C9 | Oracle-only analyses | all five | implemented as recommended (Sept 17, `notebooks/analysis/`) |
| C10 | Secondary endpoints | descriptive | implemented as recommended (Sept 17, `notebooks/analysis/`) |
| C11 | Uncertainty | cluster bootstrap CIs | implemented as recommended (Sept 17, `notebooks/analysis/`) |
| C12 | Scheduling quantity | confirmed stage ≥2 risk within the interval, stage 3 as constraint | implemented as recommended (Sept 17, `notebooks/analysis/`) |
| C13 | Policies | six, all with the confirmation rule | implemented as recommended (Sept 17, `notebooks/analysis/`) |
| C14 | Simulator integration | sequential visit generator forked after the pre-visit draws; spacing-robust policy features | implemented as recommended (Sept 17, `notebooks/analysis/`) |
| C15 | Metrics | §0 criterion on mean and P90; harm metrics; guideline baseline as slide-2 number | implemented as recommended (Sept 17, `notebooks/analysis/`) |
