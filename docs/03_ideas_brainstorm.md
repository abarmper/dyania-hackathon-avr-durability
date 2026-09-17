# Ideas that can be done by Sept 17 and will produce results

Written after re-reading all 23 references in `references/` and re-profiling the notes data, then checked line by line against every endpoint, objective and output the challenge templates ask for (`protocol/study_protocol.md` §1, §3, §4, §5; `model/approach.md` §3, §5, §6; `presentation/slides.md` slides 2, 5, 8, 9; the rubric).

Three ideas, deliberately few, that chain into one story: **define the label properly → model it → turn it into a surveillance decision**, anchored by a real-data extraction demo. A fourth is optional. Section 0 fixes the endpoints first, because the ideas only make sense against them.

Facts that constrain what is feasible:

- The notes file has **50 patients with a follow-up gradient, but 0 patients with gradients in two different follow-up years** and only 1 with both an index-year and a later gradient. The real data cannot give a longitudinal series per patient; it can give a cross-sectional "latest echo" snapshot and ~15 documented failure events.
- MIMIC-IV-ECHO needs credentialing that takes days to weeks. Not a prototype option.
- Events are rare and late (NOTION: severe SVD 1.5–10% at 10 years) and the label flickers (Velders: 59–65% of first SVD classifications are absent at the next visit). Any prototype must confront both.

So the prototype cohort has to be **simulated, calibrated to published numbers, and shown to reproduce them**. Done transparently, that is a result nobody else will have.

---

## 0. The endpoints we commit to (copy into the protocol)

**Primary objective.** In patients with a bioprosthetic aortic valve (SAVR or TAVR), can routinely collected post-implant echocardiographic and EHR data predict an individual's time to haemodynamic structural valve deterioration accurately enough to personalise the echo surveillance interval?

**Secondary objectives.**
1. Differentiate SAVR vs TAVR, and valve-model-specific, failure trajectories.
2. Optimise the echo surveillance interval: risk-adapted versus the guideline's fixed annual echo (ESC/EACTS 2025, verbatim: "within 3 months after valve implantation, again at 1 year, and annually thereafter"; `A2` 2021 said within 30 days).
3. Predict the reintervention window: expected time from first moderate deterioration to severe deterioration or valve failure, and the fraction of reinterventions that would be urgent rather than elective.
4. Quantify the reliability of the SVD label itself: agreement and persistence across the three consensus definitions (a methodological objective the field currently lacks).

**Primary endpoint.** Time from the 3-month reference echo to **VARC-3 ≥moderate (stage 2 or 3) haemodynamic SVD, confirmed on two consecutive echoes**, with death and non-SVD reintervention as competing risks. VARC-3 stage 2: mean gradient ≥20 mmHg *and* increase ≥10 mmHg from reference, or new ≥moderate intraprosthetic regurgitation; stage 3: ≥30 mmHg *and* increase ≥20, or new severe regurgitation (verbatim in `B1`). Reported as cumulative incidence at 5, 8 and 10 years. Rationale for stage ≥2 rather than failure: a surveillance model has to fire *before* failure, and moderate deterioration is the consensus trigger for closer follow-up (`A4`); the two-echo confirmation is what VARC-3 and the 2017 EAPCI statement already recommend (`A1`, point 4) and what Velders showed is necessary (`A3`). Sensitivity analyses: Capodanno 2017 and Dvir 2018 definitions; single-echo versus confirmed.

**Performance threshold for clinical usefulness.** Statistical: time-dependent AUC ≥0.75 and Uno's C ≥0.70 at 5 years, calibration slope 0.8–1.2 at 5/8/10 years, evaluated under the confirmed label. Decision-level, which is the one that matters: a risk-adapted schedule must deliver **≥25% fewer echoes per patient-decade with detection delay not more than 3 months longer than annual surveillance**, or the same echo count with earlier detection.

**Secondary endpoints.**

| Endpoint | Measurement | Timeframe |
|---|---|---|
| Bioprosthetic valve failure (VARC-3) | severe SVD, or aortic valve reintervention after a dysfunction diagnosis, or valve-related death | cumulative incidence at 5/8/10 y |
| Time to aortic valve reintervention | valve-in-valve TAVR or redo SAVR, classified elective vs urgent | 10 y |
| Mean gradient progression | annualised change in mean transprosthetic gradient from the 3-month reference, from a mixed model | whole follow-up |
| Non-structural dysfunction and other BVD | new or worsening ≥moderate paravalvular leak; clinical valve thrombosis; endocarditis | adjudicated and reported descriptively as competing events, not predicted |

Paravalvular leak and thrombosis are **not SVD** under every consensus definition, so they enter the study as adjudication rules for the label and as competing events, and thrombosis-related exposure (HALT, anticoagulation status) enters as a **predictor** (no anticoagulation HR 3.35, `B4` Table 2). Predicting them is out of scope for the prototype and said so.

**Ground truth hierarchy.** (1) Explant or autopsy pathology; (2) reintervention record with a documented dysfunction diagnosis; (3) VARC-3 echo criteria on two consecutive echoes against the 3-month reference; (4) clinical adjudication for discordant cases. Endocarditis, thrombosis and paravalvular leak are adjudicated out of the SVD label at every level.

**Subgroups.** Implant approach (SAVR/TAVR), valve model class, age at implant (<60, 60–80, >80, per `B5`), labelled valve size and indexed orifice area, sex.

**Comparators.** Time since implant alone; a "published risk-factor" Cox model using the predictor set in `B4` Table 2 as the nomogram proxy; STS-PROM as a covariate for the competing risk of death only, since it predicts mortality, not durability; cardiologist gestalt is not testable on simulated data and is planned for the real study.

**Model output per inference.** Cumulative incidence of the primary endpoint at 5/8/10 years; risk tier (low/moderate/high); top-3 contributing features; recommended next echo interval; a flag when the input echo is unconfirmed (single threshold crossing) so the clinician knows the label is provisional.

---

## Coverage check: which idea delivers which requirement

| Requirement (where asked) | Idea 1 simulator + model | Idea 2 surveillance policy | Idea 3 real-note extraction |
|---|---|---|---|
| Primary endpoint, confirmed VARC-3 stage ≥2 at 5/8/10 y | **yes, core** | uses it | cross-sectional check on ~50 patients |
| BVF and time to reintervention (secondary) | yes **after patch P1** | urgent vs elective share | real ViV/redo events, ~15 |
| Gradient progression (secondary) | yes, patch P2 makes it a reported output | – | latest gradient only |
| PVL, thrombosis, endocarditis (secondary) | endocarditis as competing event, patch P7; PVL/thrombosis not simulated | – | **yes**, extracted and tabulated |
| SAVR vs TAVR and valve-model trajectories (objective 1, subgroups) | yes **after patch P4** | policy results by stratum | valve model × outcome table |
| Surveillance interval optimisation (objective 2, slide 8) | – | **yes, core** | – |
| Reintervention timing window (objective 3) | needs P1 | **yes** with P1 | interval index→event on real cases |
| "How late are failing valves found today" (slide 2) | – | **yes**: detection delay under annual echo | last normal echo → failure on real cases |
| Label reliability across definitions (objective 4) | **yes, core** | – | flicker-prone borderline cases |
| Ground-truth hierarchy (protocol §4) | latent truth = oracle | – | **yes**: op reports, pathology, reintervention |
| Metrics: td-AUC, C, calibration (slide 9, §5) | **yes** | echoes vs delay | agreement with manual annotation |
| SHAP / feature importance (slide 9, output format) | yes **after patch P5** | – | – |
| Engineered features: gradient rate, indexed EOA, PPM flag (approach §3) | yes **after patch P3** | – | PPM flag where BSA is quoted |
| Missingness incl. informative (approach §3, data plan §4) | **yes** | – | real missingness rates |
| Fairness (protocol §6) | sex as null covariate, patch P6 | – | sex extractable for 89 patients |
| Clinical integration and safeguards (approach §6) | – | interval + referral trigger | – |
| Comparators (protocol §5) | **yes** | annual schedule is the comparator | – |

Gaps found in the first draft and now patched: reintervention was not simulated (P1); gradient rate, indexed EOA and PPM flag were implicit, not outputs (P2, P3); valve-model strata were not explicit (P4); no SHAP (P5); no fairness check (P6); endocarditis and non-structural events not handled (P7); no usefulness threshold (now in §0).

---

## Idea 1 — A literature-calibrated simulator that reproduces the label-flicker problem, then a competing-risk model trained under three label definitions

**What it is.** A generative model of each patient's *latent* transprosthetic mean gradient trajectory plus *observed* noisy annual echoes, with death, endocarditis and reintervention as competing events, and informative missingness. Run the three published SVD definitions on the observed data and show the simulator reproduces the inconsistency Velders found. Then fit survival models under (a) single-threshold VARC-3 labels, (b) two-consecutive-echo confirmed labels, and (c) the latent truth as an oracle, and report how much discrimination and calibration are lost to label noise and recovered by confirmation and by using the trajectory.

**Calibration inputs, all cited.**

| Component | Value | Source |
|---|---|---|
| Baseline (3-month) mean gradient | N(13.1, 4.7) mmHg; Magna Ease 12.6 ± 3.0; Trifecta 9.5 ± 4.1 (no SVD) vs 13.9 ± 3.8 (SVD) | `A3`; `A4`; `C1` |
| Baseline indexed EOA and PPM (patch P3) | EOAi 0.87 ± 0.17 cm²/m²; moderate PPM 51%, severe 2% (Magna Ease); EOAi 0.77 vs 0.96 in SVD vs no SVD (Trifecta); baseline gradient generated from EOAi so the two are consistent | `A4` Table 1; `C1` Table 2 |
| Per-echo measurement noise | SD ≈ 3 mmHg (5-year within-patient change 95% PI −9.6 to +7.5; per-interval SD 2–5) | `A3` Table E4 |
| Population gradient drift | 14 → 20 mmHg mean by year 15; 12.6 → 15.0 at 6.6 y | `B5`; `A4` |
| Degeneration onset time | SAVR age <70: N(10, 5) y; ≥70: N(17, 5) y; TAVR bimodal 20% N(4, 1.5) + 80% N(11.5, 3.5) | `B6` Monte Carlo inputs |
| Valve-model strata (patch P4) | Perimount-class: explant for SVD 1.9% at 10 y; Trifecta-class: SVD 4.8% at 5 y, 6.6% at 7 y; balloon-expandable THV: SVD 1.6% (XT) / 0.7% (S3) at 5 y; self-expanding THV: 2.2% at 5 y, ≥moderate 15.4% at 10 y | `B5`; `C1`; `B4` Table 1; `B1` |
| Age effect | HR 0.907 per year for severe SVD; explant for SVD at 10/15/20 y: <60 y 5.6/20/45%, 60–80 y 1.5/5.1/8.1%, >80 y 0% | `A4`; `B5` |
| Other predictors | PPM HR 7.7 (Trifecta) / severe PPM 1.85; smoking 2.58; diabetes 1.33; BMI 1.08 per unit; CKD 1.10; no anticoagulation 3.35; baseline gradient ≥15 mmHg 1.30; THV ≤23 mm 2.07 | `C1`; `B4` Table 2 |
| Reintervention given deterioration (patch P1) | explant for SVD 8.0% at 10 y as cumulative incidence; 9 of 11 stage-3 patients reintervened; NOTION reintervention 2.2–4.3% at 10 y at age 79; index-to-reintervention 6.8 ± 4.4 y (redo) and 9.8 ± 4.2 y (ViV); reintervention probability falls with age | `A4`; `B1`; `B7`; `B5` |
| Death | ~63% at 10 y at age 79; 19% at 5 y, 33% at 10 y at age 70; death before explant 46% at 10 y, 76% at 20 y | `B1`; `A4`; `B5` |
| Endocarditis (patch P7) | 7.2–7.4% at 10 y | `B1` |
| Echo missingness | 9% at 1 y, 15% at 2 y, 22% at 3 y, 35% at 4 y; more likely when the patient is well | `A3` Table E1 |

**Validation targets the simulator must hit before anything is trained.** This is the credibility step and is itself a reportable result.
- Ever-labelled SVD by 5 years: 4.6% / 2.9% / 3.0% under Capodanno / Dvir / VARC-3 (`A3`).
- Persistence of an SVD label at the next visit: 14–50%, mostly ~33% (`A3` Table 2).
- Regression to the mean by decile: lowest decile rises 1.2–2.3 mmHg per interval, highest falls 1.0–5.9 (`A3` Table E4). This should *emerge* from the noise model.
- ≥moderate SVD 15–21% and severe 1.5–10% at 10 years in a 79-year-old cohort (`B1`).
- Freedom from stage 2/3: 98.5% at 5 y, 61% at 10 y at age 70 (`A4`).

**Models and outputs.** Cause-specific Cox and Fine-Gray (`D5`) on baseline features; random survival forest (`D4`); a landmark model at 1, 3 and 5 years using last gradient, slope from the reference, and baseline covariates; a hierarchical version with partial pooling across the four valve strata, tested by holding out one small stratum as the "new valve with little history" case (patch P4). Comparators: time since implant alone; the `B4` risk-factor Cox. Metrics: time-dependent AUC and Uno's C (`D6`) at 5/8/10 years, calibration plots, all under the three labels and the oracle, overall and by stratum. Secondary endpoint outputs: annualised gradient slope per patient (patch P2); BVF and reintervention cumulative incidence (P1). Feature contributions (patch P5): permutation importance for the forest, plus SHAP on a gradient-boosted 5-year-horizon classifier to produce the "top-3 features" per patient the output format asks for. Fairness (patch P6): sex simulated with no true effect; report calibration by sex and confirm the model does not manufacture an effect. Report to TRIPOD+AI (`D8`).

**Expected headline.** "Under single-echo VARC-3 labels the C-index is X with poor calibration because ~60% of 'events' are noise; requiring two consecutive echoes raises it to Y; using the trajectory to Z, close to the oracle." The shape of that result is guaranteed by construction, which is the point.

**Effort.** One person, one day, Python with `lifelines`, `scikit-survival`, `shap`. Simulator ~200 lines with the patches.

---

## Idea 2 — Turn the prediction into a surveillance interval and count echoes saved versus events missed

**What it is.** On the same simulated cohort, compare surveillance policies head to head:
1. **Guideline**: reference echo (30 days, within the 2025 three-month window), 1 year, then annually (ESC/EACTS 2025; `A2`).
2. **Risk-adapted**: after each echo, schedule the next when the model's predicted cumulative risk of stage-3 SVD reaches a threshold (e.g. 5%), bounded to [6 months, 3 years]. A simplified form of the Rizopoulos information-gain rule (`D1`), which sets the next measurement where conditional survival drops to a constant κ.
3. **De-escalation only**: annual to 5 years, then 2-yearly if low risk, the precedent from endovascular aneurysm surveillance (`D10`).

**Metrics.** Echoes per patient over 10 years; delay from latent stage-2 onset to confirmed detection; fraction of stage-3 cases detected while still stage 2; fraction of reinterventions that are urgent (stage 3 reached before any echo caught stage 2), using the reintervention sub-model from P1. Translate urgent share with the elective-vs-emergency redo mortality figures in `B3` (1.4% vs 22.6%) into expected deaths per 10,000, labelled illustrative. Report by SAVR/TAVR and valve stratum.

**This is also the number for slide 2.** The detection delay under the annual policy, on a simulator that reproduces the published flicker, is a defensible estimate of "how late failing valves are identified today", which no reference we found actually reports.

**Effort.** Half a day once Idea 1 exists; the policies are loops over the simulated trajectories.

---

## Idea 3 — Extraction and label construction on the 117 real notes

**What it is.** Run a language model over the 215 notes to extract, per patient: procedure type, valve model and label size, approach, native morphology, baseline post-implant gradient if quoted, latest gradient, regurgitation grade, **paravalvular leak grade, thrombosis or HALT mention, endocarditis**, ejection fraction, BSA where quoted (for the PPM flag), sex, any valve-in-valve or redo event with year and the surgeon's explant findings (e.g. "cusp tear", "prosthetic thickening/calcification"). Hand-check ~30 notes for agreement. Apply the three SVD definitions where both a reference and a latest gradient exist; stage the ~15 known failure patients through the ground-truth hierarchy (pathology → reintervention record → echo criteria).

**Why it works.** Fine-tuned local transformers reach 95.7% accuracy on echo reports and 93.6% on the prosthetic aortic valve field (`E1`); open-weight LLMs reach >96% on prosthetic valve detection (JAMIA 2025, in `references/README.md`). 36 progress notes already quote a prosthetic valve size inside echo text; 50 patients have a follow-up gradient; earlier profiling found paravalvular mentions in 18 patients, endocarditis in 9, thrombosis in 4.

**What it produces.** A structured table of ~100 patients by valve model, approach and outcome; a figure of valve model × failure; the non-structural events table the secondary endpoints ask for; the real-data missingness rates for the data plan; the interval from last documented normal echo to failure for the cases where it exists (slide 2); and the flicker-prone group, patients whose single latest gradient is 20–29 mmHg with no reference echo, who are unclassifiable under all three definitions. That last group is the real-world argument for Idea 1.

**Honest limits to state.** Cross-sectional; failure-enriched sample; year-only dates; brand names sometimes eaten by de-identification. A pipeline demo and case series, not a validation set.

**Effort.** Half a day for a second team member, in parallel with Idea 1.

---

## Idea 4 (optional) — Reconstruct published Kaplan-Meier curves to replace the parametric onset assumptions

Use the Guyot algorithm (`D9`) on the curves in `A4` (numbers at risk 338 / 247 / 34 / 20 at 0 / 5 / 9.5 / 10 years are printed) and `B1`, then fit the onset-time distribution to those instead of the `B6` normals. Upgrades "assumed" to "reconstructed from published data" in the data plan. Needs manual digitising, ~30 minutes per curve. Only if a third person is free; otherwise the `B4` Table 1 point estimates are enough to sanity-check the simulator.

---

## What not to spend time on

- Training anything on the 117 notes. Fifty snapshots and fifteen events cannot support a model.
- Deep survival models (`D3`, `D7`) on a few thousand simulated patients. Cite as the next step at registry scale.
- Applying for MIMIC for the hackathon. Do apply, because it is the real validation cohort, but frame it as pre-registered external validation.
- Predicting paravalvular leak or thrombosis. They are not SVD; adjudicate them out and report them.

## Suggested split of the remaining time

| When | Person A | Person B |
|---|---|---|
| Sept 16 afternoon | Idea 1: simulator with P1, P3, P4, P7; hit the validation targets | Idea 3: extraction on notes |
| Sept 16 evening | Idea 1: models under three labels, metrics, SHAP, subgroups | Idea 3: hand-check, definitions, non-structural events, figures |
| Sept 17 morning | Idea 2: policies, echoes-vs-delay figure, urgent share | Protocol §1–§5 from §0 above and `docs/02` |
| Sept 17 midday | Notebook cleanup, README | Slides |

The three ideas produce, respectively: a results table (label noise vs model performance, overall and by valve stratum, with the secondary endpoints), a figure (echoes vs detection delay under two policies, plus the slide-2 lateness number), and a real-data table (valve model × outcome and non-structural events from free text). That covers every endpoint in §0 and every output the templates ask for.
