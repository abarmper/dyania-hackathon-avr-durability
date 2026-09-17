# Modelling (Idea 1) and surveillance policy (Idea 2)

Implements the recommended options of `docs/05_modelling_and_surveillance_plan.md` on the synthetic cohort in `data/synthetic/`.

```bash
cd notebooks/analysis
/data/abar/alexenv/bin/python run_all.py all         # features -> models -> evaluation, ~6 min
/data/abar/alexenv/bin/python run_all.py explain     # SHAP, permutation importance, parsimonious model, example inference records
/data/abar/alexenv/bin/python run_all.py policy      # sequential surveillance simulation, 13 policy runs in parallel
/data/abar/alexenv/bin/python run_all.py oracle      # achievable ceiling, definition audit, reference-echo value, mechanism shifts, sample size
/data/abar/alexenv/bin/python run_all.py secondary   # stage 2 -> 3, reintervention breakdown, BVF / reintervention CIF, gradient slope
/data/abar/alexenv/bin/python summarise.py           # output/results_summary.md (numbers for the templates)
/data/abar/alexenv/bin/python -m pytest tests -q
```

Outputs go to `output/` (tables as CSV, figures as PNG, models as joblib, `metrics.json`, `results_summary.md`); the notebooks
`../02_features_and_models.ipynb` and `../03_surveillance_policy.ipynb` display them.

## What is modelled

- **Instances**: one per post-reference echo (the landmark; 44,000) plus one per valve at the reference echo (implant-time model).
- **Outcome (C0)**: confirmed VARC-3 stage ≥2 SVD, *known at the confirming echo*; competing: death, endocarditis; reinterventions
  for SVD that truncated follow-up before the confirming echo count as events (all 284 carry an SVD indication); loss to follow-up and
  administrative censoring are censoring. The risk set at landmark L is `T_known > L`.
- **Primary model**: discrete-time cause-specific hazard — person-period expansion (6-month periods, 5-y horizon for landmark models,
  10-y for the implant-time model), multinomial `HistGradientBoostingClassifier` (scikit-learn, OpenMP-threaded, native NaN handling,
  early stopping on a valve-grouped validation split); CIF = Σ S_{j−1} h_{1j}.
- **Comparators**: registry-only feature set, no-symptom-flag set, time-since-implant-only hazard, implant-time boosted model,
  cause-specific Cox + Aalen–Johansen (lifelines; HR table), published risk-factor Cox (B4 Table 2 set).
- **Features** (`features.py`, one vectorised function used both for the training table and online in the policy loop): baseline,
  reference echo (incl. ASE expected EOA and residual), and longitudinal (last values and change from reference, log-scale EWMA of the
  gradient with the posterior probability that the true gradient is ≥20 mmHg, slopes, max so far, current single-echo stage, count of
  prior positive flags = the unconfirmed-echo flag, consecutive worsening, symptom-triggered flag, channel-missing indicators).
- **Evaluation** (`evaluate.py`): cause-specific time-dependent AUC on fixed-horizon eligible subsets (competing events are controls;
  LTFU before the horizon dropped and counted), truncated Wolbers C, Brier, decile calibration with logistic recalibration slope;
  valve-level bootstrap (500, paired resamples across models); reported on all instances and on currently-negative instances;
  label sweep with fixed instances (prevalent-or-incident cases); subgroups; fairness by sex.
- **Explainability** (`explain.py`): SHAP TreeExplainer on the boosted model (ordinal codes, so SHAP explains the model's own splits),
  permutation importance, greedy backward elimination to ≤8 features keeping ≥95% AUC, per-inference records (CIF, tier, top-3, flag).
- **Oracle-only analyses** (`oracle_analyses.py`): achievable-AUC ceiling, definition audit against the latent crossing, value of the
  30-day reference echo, mechanism shifts (noise ×0.5/×1.5, missingness ×2, progression ×0.7; regenerated with the calibrated
  parameters), simulation-based sample size.
- **Surveillance policies** (`policy.py`): sequential visit generator on top of the simulator's `pre_visits`/`post_visits` halves
  (behaviour-preserving refactor, regression-tested); pre-drawn noise bank by echo ordinal and pre-drawn reintervention variates
  (common random numbers across policies); confirmation echo 3–6 months after any positive single echo; forced symptom-triggered and
  other-indication echoes; year-indexed missingness with the well/symptomatic multipliers; policies: guideline, guideline+confirmation,
  de-escalation, time-only null, risk-adapted (τ sweep, stage-3 safety cap), oracle ceiling. Metrics: echoes per patient-decade,
  confirmed and first-flag delay vs the latent crossing (mean, P90), stage-3-seen-as-stage-2 share, urgent share, false-positive
  reinterventions, illustrative deaths per 10,000.

## Parallelism

Main-process boosted fits use OpenMP (`ANALYSIS_OMP_THREADS`, default min(cpu, 64)); independent runs use loky processes
(`ANALYSIS_N_JOBS`, default min(cpu, 48)) with inner thread limits to avoid oversubscription: permutation importance and backward
elimination, mechanism-shift regenerations, the sample-size curve, and the 13 policy runs. The person-period table is < 1 M rows, so
the GPUs are not needed; `xgboost` is installed as an optional backend but not used by default.

## Caveats stated in the results

The symptom-triggered echo flag is stronger than in reality (symptoms only occur at latent stage ≥2 in the simulator); death and SVD
are independent in the simulator, so mortality translations are illustrative; the physics-informed features are framed as ASE 2024 /
VARC-3 recommendations and their increment over the registry-only set is reported as measured.
