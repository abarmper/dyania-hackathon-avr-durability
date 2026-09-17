# Synthetic bioprosthetic-AVR cohort simulator

Generates the prototype cohort for Idea 1 (`docs/03_ideas_brainstorm.md`) from the specification in
`docs/04_synthetic_cohort_spec.md`: one row per implanted valve, a long-format table of echo measurements,
an events table, per-valve label times, and a separate latent-truth table that is never used for training.

```bash
/data/abar/alexenv/bin/python cli.py calibrate            # ~12 min: physics -> noise -> death -> onset per valve class -> checks
/data/abar/alexenv/bin/python cli.py validate             # acceptance tests -> output/validation_report.md + output/figures/
/data/abar/alexenv/bin/python cli.py generate --n 10000   # -> data/synthetic/ (CSV + parquet, data_dictionary.md)
/data/abar/alexenv/bin/python -m pytest tests -q
```

## Why the approach is valid

1. **One latent state, everything else derived.** The true effective orifice area `EOA*(t)` is the only quantity
   that degenerates. Mean gradient, peak velocity and DVI follow from flow physics
   (`MPG = c·(Q/EOA_eff)²`, `EOA_eff = EOA·(Q/Q_ref)^0.35`, `DVI = EOA/A_LVOT`), with the constant `c` fitted per
   valve model and size to the normal-value tables of the 2024 ASE prosthetic-valve guideline (`references/A5`,
   Appendix Tables A1/A2/A4). Label flicker, regression to the mean and disagreement between SVD definitions
   therefore *emerge* from measurement error; nothing is injected.
2. **Covariates act on hazards, never on times.** Onset comes from a per-class Weibull with proportional hazards,
   `T = λ·(−ln U / HR(x))^{1/k}`. Only biological effects enter `HR(x)` (smoking 2.58, diabetes 1.33, BMI 1.08
   per unit, CKD 1.3, age 0.951 per year for surgical valves, anticoagulation and thrombosis effects). Effects with
   a mechanical pathway (PPM, high baseline gradient, small THV) act through the reference orifice area and are
   checked afterwards as *emergent* Cox hazard ratios on the simulated labels, so nothing is double counted.
3. **Calibrate to what the papers measured.** Each curve target is matched by simulating a cohort that mirrors the
   paper's baseline table (age, sizes, PPM, valve class, reference-echo timing, visit schedule, completeness,
   mortality) and applying the paper's **single-echo** label rule and estimator (KM with death censored for
   Kermen Fig 2B, Aalen–Johansen elsewhere). The confirmed two-echo primary endpoint is an *output*, never a target.
4. **Over-identification.** Kermen (Perimount-class, 70 y) fixes that class, Wakami fixes Trifecta-class, NOTION
   TAVI fixes self-expanding THV; NOTION SAVR at 79 y is then *predicted* (Perimount + Trifecta transported to 79 y by
   the age hazard) before the porcine/Mitroflow share is freed. Kermen stage 3 and NOTION severe identify the
   progression speed given onset. The report prints every one of these checks, pass or fail.
5. **Competing risks and informative processes.** Death depends on age (which also protects against SVD),
   missingness rises with time and is higher when the patient feels well, symptoms trigger extra echoes, and
   reintervention depends on *detected* status, so a surveillance policy has consequences (Idea 2).
6. **Provenance.** Every parameter carries a source tag or an `ASSUMPTION` flag (`params.py`, rendered into
   `data/synthetic/data_dictionary.md`); assumptions are the sensitivity-analysis list.

## Modules

| file | role |
|---|---|
| `params.py` | all parameters with source / assumption tags; `SimParams` (JSON round trip) |
| `valve_tables.py` | ASE 2024 normal values by model × size, valve classes, size and model mixes, published point targets, paper cohort profiles |
| `physics.py` | flow–gradient–DVI relations; per model × size constants fitted by Monte Carlo to the ASE rows |
| `baseline.py` | `patients_valves` draw (general cohort or paper-matched) |
| `latent.py` | onset hazard, `EOA*(t)`, regurgitant mode, transient thrombosis, noise-free VARC-3 stage, crossing times |
| `events.py` | death (attained-age Gompertz, LVEF, dialysis, frailty), endocarditis, loss to follow-up, symptoms, reintervention |
| `visits.py` | schedule, jitter, informative missingness, symptom-triggered echoes, four-source observation noise, channel missingness |
| `labels.py` | Capodanno, Dvir, VARC-3 (full and NOTION haemodynamic), single-echo and confirmed, first-event times |
| `simulate.py` | orchestrator; `simulate_matched()` implements the `calibrate_onset` interface |
| `calibrate.py` | sequential calibration → `output/params_calibrated.json`, `output/calibration_report.md` |
| `validate.py` | acceptance tests → `output/validation_report.md`, `output/figures/` |
| `cli.py` | `calibrate` / `validate` / `generate`; writes the data dictionary |

## Generative model in one paragraph each

**Baseline.** 55% SAVR (Perimount Magna Ease 45%, Inspiris 10%, Trifecta 25%, Mitroflow/Epic/Mosaic 15%, new 5%) and
45% TAVR (Sapien XT/3 55% by era, CoreValve/Evolut R 40%, new 5%); age N(70, 9) / N(78, 6); height and weight drawn
jointly by sex (BSA 2.0 ± 0.2, BMI 29 ± 5); diabetes 27%, CKD 11%, smoking 15%, AF 30% generating anticoagulation
(VKA before 2015, DOAC after); LVEF N(59, 10); stroke-volume index N(42, 4); labelled size conditional on approach and
BSA; reference EOA from the model × size table (measured SD deflated by the 17% measurement CV); 30% of valves lack a
30-day reference echo and use the table instead (D1); implant dates 2010–2020, censored 2026-09-15.

**Latent process.** Onset `T_deg` ~ PH-Weibull per class (calibrated). Stenotic mode:
`EOA*(t) = EOA_ref·(1 − 0.011 t)·exp(−(t − T_deg)+/τ)` with `τ` lognormal (median 7.5 y surgical, 6 y Trifecta, 10 y THV;
σ 0.4). Regurgitant mode (25% surgical, 60% Trifecta, 5% THV): intraprosthetic AR jumps to moderate at onset and to
severe ~2 y later. Transient thrombosis/HALT in year 1 (TAVR 8%, SAVR 3%, halved on anticoagulation) reduces `EOA*`
by 25% for 3–9 months and raises later onset hazard ×1.5. A per-valve slow drift of flow state (log Q, SD 0.05/yr)
gives the between-patient trend heterogeneity seen in Velders' 5-year change interval.

**Observation.** Four independent per-visit error sources with the correct cross-channel structure: flow (MPG ∝ f²;
EOA, DVI unchanged), LVOT VTI sampling (EOA and DVI × e_L), CW alignment (MPG × e_A², EOA and DVI ÷ e_A), LVOT
diameter (EOA × e_D²); plus a per-valve persistent LVOT bias that cancels in changes from reference. Sigmas
calibrated to Velders' decile regression-to-the-mean table, 5-year change interval, ever-labelled rates and label
persistence. AR/PVL grade misclassified ±1 with p = 0.03; EOA missing in 8% and DVI in 5% of echoes.

**Visits and events.** Reference at 30 d, then 1 y and annually (ESC/EACTS 2025: baseline within 3 months; 2021: within 30 days) with jitter; scheduled-visit missingness
9/15/22/35% at years 1–4 (Velders) saturating at 45%, × 1.3 when well, × 0.6 when symptomatic; symptom-triggered
echoes; death by attained-age Gompertz (slope 8.5%/y, LVEF 0.974 per point, dialysis ×2, gamma frailty, higher-risk
TAVR era ×1.8) calibrated to Kermen's survival curve; endocarditis 0.75%/y; reintervention only after *detected*
stage 3 (1.5/y at age 70, × exp(−0.15·(age − 70))) or confirmed stage 2 with symptoms; ViV vs redo by age; urgent if
the first stage-3 detection came from a symptom-triggered echo. Follow-up of the valve ends at reintervention (D7).

**Labels.** Per echo vs the reference: Capodanno, Dvir, VARC-3 full (NA when EOA and DVI are both missing) and the
NOTION haemodynamic variant; each single-echo and confirmed (same or higher stage at the next available echo).
Primary endpoint = confirmed VARC-3 stage ≥ 2, time = first echo of the confirmed pair. Oracle = the same rules on
the noise-free channels at the same visits, plus continuous-time crossing times for the surveillance-delay metric.

## Calibration and validation results (Sept 16)

`output/calibration_report.md`: 40 checks, 2 marginal misses (NOTION severe SVD at 10 y, SAVR 13.4 vs 9.9 and
TAVI 4.6 vs 1.4, i.e. 3–4 pp against a 3 pp tolerance). Curve fits: Kermen stage 2/3 within 1.3 pp, Wakami 2.1,
NOTION TAVI ≥moderate 2.9; the zero-parameter prediction of NOTION SAVR lands at 21.6% vs 20.7% at 10 y. Wall
time 4 min on 32 processes (`N_WORKERS` in `calibrate.py`).

`output/validation_report.md` (37 checks, 7 misses): all Velders label-flicker statistics pass (ever-labelled
5.2 / 3.6 / 4.2% vs 4.6 / 2.9 / 3.0; persistence 20% within the published 14–50%; 5-y change interval
(−9.6, +10.3) vs (−9.6, +7.5)); mortality, explant and reintervention pass; Kermen moderate PPM 40% vs 51% and one
model × size gradient cell are marginal. Three emergent hazard ratios fall short of their published values
(reference MPG ≥15: 0.91 vs 1.30; THV ≤23 mm: 1.07 vs 2.07; smoking: 1.79 vs 2.58): effects that act through the
reference gradient are attenuated on observed labels by regression to the mean of the 30-day reference, and label
noise dilutes every effect. This is stated as a finding in `docs/04` §5 rather than tuned away.

Generated cohort (`data/synthetic/`, seed 1): 10,000 valves, 50,849 echoes (397,711 measurement rows), 4,009
deaths, 1,589 reinterventions (52% urgent), 572 endocarditis; confirmed VARC-3 stage ≥2 incidence 6.2 / 11.3 /
16.8% at 5 / 8 / 10 y. `echo_visits.csv` (23 MB) is not committed; use the parquet copy or regenerate.

The notebook `notebooks/01_generate_and_validate.ipynb` displays both reports with the figures.

## Known limitations (stated in the protocol)

HALT gradient effect unsourced; NOTION SAVR valve mix assumed; population age-mortality slope vs Kermen's adjusted
within-cohort HR; single-echo VARC-3 labels have a noise floor of ~1–2% by 5 y, so the Sapien 3 point (0.7%) is a
check, not a target; no inter-centre variation beyond a per-valve bias; PVL as null covariate only; progression speed
identified from few stage-3 events; symptom hazard shape assumed; labels applied without the morphology criteria of
ASE Table 6.
