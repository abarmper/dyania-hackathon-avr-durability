# What we will model: the synthetic cohort, its schema, and where each feature exists in real data

Companion to `03_ideas_brainstorm.md` (Idea 1, then 2, optionally 4). This fixes *what* is simulated, *why*, and *which real dataset carries the same field*, so the data plan can say honestly: "every feature in the prototype exists in a named registry or EHR table; the values are simulated, the schema is not."

Sources verified in this session: the STS/ACC TVT Registry v3.0 TAVR Data Collection Form (public PDF, field numbers quoted below), the STS Adult Cardiac Surgery Database v4.20.2 Data Specifications (public PDF, short names quoted), the MIMIC-IV-ECHO v1.0.1 page, and the 23 papers in `references/`. Counts "in provided notes" are from re-profiling `data/notes_deidentified.xlsx` (117 patients) and the labs/meds files (17 patients) today.

---

## 1. Anchoring decision: which real schema do we imitate?

No single real dataset has everything. Each has a different piece:

| Real source | Has | Lacks | Access |
|---|---|---|---|
| **TVT Registry v3.0 (TAVR)** | baseline, discharge, 30-day and 1-year echo (AV mean gradient, AV area, LVEF, central and paravalvular regurgitation), device, comorbidities, discharge antithrombotics, reintervention with indication, valve thrombosis, death | follow-up stops at 1 year (later only via CMS linkage); no explicit SVD staging; TAVR only | research proposal, analysis at their centre |
| **STS ACSD v4.20.2 (SAVR)** | prosthesis implant with device identifier (`AorticImplant`, `VSAoImUDI`, `VSAVPr`), explant with device identifier and year (`ValExpUDI`, `ValExp2YrImplantKn`), pre/post LVEF, aortic gradient (`VDGradA`), comorbidities, discharge anticoagulant and statin (`DCDirOralAnticoag`, `MedLipid`) | no post-discharge echo; reintervention seen only as a new admission | member PUF program |
| **MIMIC-IV + MIMIC-IV-ECHO** | serial structured echo 2008–2022 in long format (`subject_id, measurement, result, unit`), 180–230 variables incl. valvular function and Doppler haemodynamics; ICD-coded AVR/TAVR; labs; meds; death | valve model and size (free text only, notes not yet released); exact gradient variable names unverified until credentialed | credentialed, days to weeks |
| **Trial CRFs (PERIGON, NOTION, PARTNER)** | core-lab MPG, EOA, DVI at discharge, 3–6 mo, then yearly; adjudicated events | proprietary | data-sharing request, months |

**Decision (recommended):** the prototype's schema is the **union**, organised as three tables that mirror how the real sources are shaped:

1. `patients_valves` — one row per implanted valve, baseline fields named after TVT DCF and STS ACSD fields.
2. `echo_visits` — one row per echo per valve, **long format identical in spirit to MIMIC-IV-ECHO** (`valve_id, visit_time, measurement, result`), with the VARC-3 variables (mean gradient, EOA, DVI, regurgitation grades) as the measurement names.
3. `events` — death, reintervention (type and indication, as in the TVT follow-up form), endocarditis.

This lets the same notebook code run unchanged on MIMIC-IV-ECHO later, and lets the protocol name the TVT and STS fields as the deployment source for every baseline variable.

---

## 2. What is modeled, and why

The object being modeled is **one prosthetic valve followed from implant to reintervention, death or censoring**, with a hidden state (true haemodynamic function) observed through noisy, irregularly timed echoes. Everything else exists to make that realistic enough to test a surveillance model.

### 2.1 Baseline patient and valve (`patients_valves`)

| Variable | How generated | Why it is in | Source of distribution |
|---|---|---|---|
| `approach` SAVR / TAVR | Bernoulli, default 55/45 | secondary objective 1; different onset distributions | era mix 2010–2020 |
| `valve_model` (nested in `valve_class`) | Named models with a class parent for partial pooling. SAVR: Perimount/Magna Ease 45%, Inspiris-RESILIA 10% (Perimount-class); Trifecta 25% (Trifecta-class); Mitroflow/Epic/porcine 15% (porcine-early-failure class); one "new valve" model 5% with almost no history. TAVR: Sapien XT 10%, Sapien 3 45% (balloon-expandable); CoreValve/Evolut 40% (self-expanding); new valve 5% | exact model where a published rate exists, class otherwise; the hierarchical model borrows strength for sparse models | `B5`, `A4`, `C1`, `B4` Table 1, `B1`, `A1` Table 1 (Mitroflow) |
| `age_at_implant` | SAVR N(70, 9) truncated 40–90; TAVR N(78, 6) truncated 60–95 | strongest known predictor; competing death | `A3` Table 1; `B1`; `B5` |
| `sex` | 45% female | fairness check (no true effect simulated) | `A4`; `A3` |
| `height_cm`, `weight_kg` → `bsa`, `bmi` | height and weight drawn jointly by sex (male 175 ± 7 cm, female 162 ± 7 cm; weight to give BMI 29 ± 5), then BSA by Mosteller and BMI computed, never drawn independently (targets BSA 2.0 ± 0.2, `A3`) | indexed EOA, PPM, BMI hazard; avoids impossible BSA/BMI pairs | `A3` Table 1; `B4`; TVT Height 6000 / Weight 6005 |
| `diabetes`, `ckd` (eGFR <60 or dialysis), `smoking` | prevalences 27%, 11%, current smoker 15% | the three comorbidities with a quantified SVD hazard: diabetes 1.33, CKD 1.10–1.34, smoking 2.58 | `B4` Table 2; `A4` Table 3; `A3` Table 1 |
| `af` | 30% at implant (higher after SAVR), drives anticoagulation | comorbidity; the indication behind anticoagulation, so the model faces real confounding | `B1` (new-onset AF 52–74%); `A4` |
| `anticoag_discharge`, `anticoag_type` (VKA / DOAC / none) | generated **through** AF and era: 85% of AF patients anticoagulated (VKA before 2015, DOAC after), 10% of non-AF; plus a short post-SAVR course in some centres | effect is type- and time-specific, not one arrow: absence of anticoagulation in the **first year after TAVR** raises early HVD (thrombosis pathway, HR up to 3.35); **long-term VKA after SAVR** modestly raises late SVD (calcification pathway, Salaun's late-HVD factor); DOAC neutral | `B4` Table 2 and text (refs 84–86, GALILEO); `B4` Salaun early-vs-late factors; `A2` HALT 12.4% vs 32.4% |
| `vmax_ref`, `svi_ref` | peak velocity derived from the latent gradient (peak ≈ 1.6–1.8 × mean, own noise); stroke-volume index from LVEF and LV size, low-flow if <35 mL/m² | Vmax is what many real reports quote; SVi is needed to interpret flow-dependent gradients | `B5`; `B3`; TVT Low Flow 13700 |
| `ar_ref`, `pvl_ref` | intraprosthetic AR grade at reference (none 85%, trace/mild 14%, ≥moderate 1%); paravalvular grade for TAVR (none 50%, mild 40%, ≥moderate 10% for early-generation, lower for Sapien 3) and SAVR (≥mild 5%) | baseline for "new or worsening" regurgitation in all three definitions; PVL is a covariate only, never an SVD outcome | `B1` (PVL 18–53% TAVI, 5.2% SAVR; intraprosthetic 5.8 / 2.2%); `A3` (0.2% ≥moderate at 5 y) |
| `label_size_mm` | conditional on approach and BSA: SAVR {19:15%, 21:35%, 23:30%, 25:15%, 27:5%}; TAVR {23:25%, 26:40%, 29:30%, 34:5%} | PPM; THV ≤23 mm hazard | `A4` Table 1; `B1` |
| `eoa_ref`, `eoai_ref`, `ppm` | EOA from a per-class-and-size normal table with SD 0.25 cm²; EOAi = EOA/BSA; PPM moderate if EOAi ≤0.85 (0.70 if BMI ≥30), severe if ≤0.65 (0.55) | VARC-3 PPM; the mismatch → SVD pathway | `B1` thresholds; `A4` EOAi 0.87 ± 0.17 |
| `mpg_ref` | derived from EOA and stroke volume through a Gorlin-type relation, then jittered: target mean 13.1 ± 4.7 mmHg (SAVR), lower for TAVR; TAVR small-annulus higher | primary reference for all three definitions; baseline gradient ≥15 mmHg hazard | `A3`; `A4`; `C1`; `B1` |
| `lvef` | N(59, 10) | flow correction; low-flow gradients | `A3` |
| `implant_year` | uniform 2010–2020; administrative censoring at 2026-09 | realistic right-censoring pattern, temporal validation split | design |

Not simulated at baseline, stated as such: race/ethnicity (absent from every source we can calibrate to; fairness across race is deferred to the real cohort), STS-PROM (a mortality score; enters only via age and comorbidity in the death model), frailty and NYHA class (only as the symptom mechanism below), annulus CT geometry, commissural alignment.

### 2.2 Latent valve function over time (the thing that actually degenerates)

Model the **true effective orifice area** `EOA*(t)` of the valve as the latent state; derive gradient and DVI from it. Reason: degeneration is physically a loss of orifice area (calcification, stiffening) or a tear (regurgitation); gradient is a flow-dependent consequence, which is exactly why gradients flicker.

- **Onset time** `T_deg` per valve, drawn from the class- and age-specific distribution: SAVR age <70 N(10, 5) y, ≥70 N(17, 5) y; TAVR 20% N(4, 1.5) + 80% N(11.5, 3.5); Trifecta-class shifted earlier to match 4.8% at 5 y and 6.6% at 7 y; multiplied by the published hazard ratios (PPM 1.85–7.7 depending on class, smoking 2.58, diabetes 1.33, BMI 1.08/unit, CKD 1.10, no anticoagulation 3.35, baseline gradient ≥15 mmHg 1.30, THV ≤23 mm 2.07, age per Kermen 0.907/yr). Sources: `B6`, `C1`, `B4` Table 2, `A4`, `B5`.
- **Trajectory**: before `T_deg` a slow drift (population mean gradient 14 → 20 mmHg by year 15, `B5`); after `T_deg` an accelerating decline in `EOA*` so that stage 2 is reached on average ~1.5–2.5 years after onset and stage 3 ~2–4 years later (Kermen: stage 2 at mean 7.6 y, stage 3 at 7.8 y, i.e. they can be close; Johnston: gradient rises "exponentially faster" in eventual explants).
- **Regurgitant mode** for a minority (~25% of degenerations, Trifecta-class higher because of cusp tears, `C1`): new intraprosthetic AR that jumps grades rather than a gradient rise.
- **Derived true values** at any time: `MPG*(t) = c · (SV / EOA*(t))²` with stroke volume from BSA and LVEF; `DVI*(t) = EOA*(t) / LVOT_area`.

### 2.3 Observation process (what an echo actually returns)

Each echo returns noisy copies: `MPG_obs = MPG* + ε`, ε ~ N(0, 3²) mmHg (from Velders' within-patient 5-year change interval −9.6 to +7.5 and per-interval SDs 2–5, `A3`); `EOA_obs` and `DVI_obs` with proportional noise of ~12% and ~10% (NOTION and Velders both note EOA and DVI are noisier than gradient); regurgitation grade with a one-grade misclassification probability of ~10%. Noise is independent between visits: this single assumption is what produces regression to the mean and the label flicker, and the simulator is validated by checking it reproduces Velders' decile table.

**Echo timing (design choice D4 below):** guideline schedule (discharge, 30 days, 1 year, annual) with jitter; missingness rising with time (9% at 1 y to 35% at 4 y, `A3` Table E1) and **higher for patients who feel well**; plus **symptom-triggered echoes** when the latent state is stage ≥2 and the patient is symptomatic. Symptom onset is modeled as a hazard that rises with `MPG*` (most stage-2 disease is subclinical: Salaun found 31% degeneration, mostly silent).

### 2.4 Competing and terminal events (`events`)

- **Death**: Gompertz hazard by age (HR 1.033 per year, `A4` Table E2) with an LVEF multiplier (HR 0.974 per percentage point, `A4` Table E2) and a dialysis multiplier; no other comorbidity multipliers, because none is quantified in our sources; calibrated to ~63% at 10 y at age 79 (`B1`), 33% at 10 y at age 70 (`A4`), 46% dead before explant at 10 y (`B5`). Death depends on age, which is also protective for SVD, so censoring is informative, exactly as `A1` warns.
- **Reintervention**: hazard conditional on *detected* stage 3 (or symptomatic stage 2), scaled down with age and frailty; calibrated to 8% cumulative explant for SVD at 10 y in a 70-year-old cohort (`A4`), 2–4% at 10 y at 79 (`B1`), 9 of 11 stage-3 patients reintervened (`A4`). Type: ViV vs redo drawn by age (`B7`: ViV recipients 79 vs redo 71 years). **Urgent** if the first detection of stage 3 came from a symptom-triggered echo. The valve is the unit of analysis; follow-up of that valve ends at reintervention.
- **Endocarditis**: constant hazard giving ~7% at 10 y (`B1`); removes the valve from the SVD label (adjudicated out) and usually triggers redo.
- **Transient dysfunction (thrombosis / HALT), switchable module, default on.** A reversible episode with hazard concentrated in the first year (TAVR ~8%, SAVR ~3%, halved on anticoagulation), producing a +8 to +15 mmHg gradient bump for 3–9 months that then resolves. It is *not* SVD and is adjudicated out of the oracle. Its purpose is to test the label's specificity: a single-echo rule will call these SVD; a two-consecutive-echo rule mostly will not. Sources: `A2` (HALT 12.4% on anticoagulation vs 32.4% on antiplatelets at 3 months), `B4` (SLT 6–15%, up to 31% at 1 year; resolves with 3–6 months of anticoagulation), `A3` (HALT named as a transient cause of inconsistent classification).
- Not simulated: paravalvular leak *progression* (PVL enters only as a baseline covariate), pannus as a distinct process.

### 2.5 Labels computed from the observed data

Three definitions applied exactly as published, against the reference echo: Capodanno 2017 (MPG ≥20 and rise ≥10, or ≥40 / rise ≥20), Dvir 2018 (rise >10 with EOA and DVI decrease), VARC-3 (MPG ≥20 and rise ≥10 with EOA fall ≥0.3 cm² or ≥25% or DVI fall ≥0.1 or ≥20%; stage 3 at ≥30 / rise ≥20 / EOA fall ≥0.6 or ≥50%), each as single-echo and as two-consecutive-echo confirmed. The oracle label is the same rule applied to the noise-free `MPG*`, `EOA*`, `DVI*`.

---

## 3. Feature-by-feature: where each simulated variable exists in real data

Counts under "provided notes" are patients out of 117 whose free text contains the item; labs/meds are out of 17.

| Simulated variable | Provided notes / labs / meds | TVT Registry v3.0 field | STS ACSD v4.20.2 | MIMIC-IV / -ECHO | Trials |
|---|---|---|---|---|---|
| approach (SAVR/TAVR) | 56 / 57 | whole form | `VSAVSurgType`, `VSAVPr` | ICD-10-PCS 02RF0*/02RF3* | yes |
| valve class / model | 97 named | device fields (model list), sheath 13509 | `AorticImplant`, `VSAoImUDI` (UDI) | free text only | yes |
| label size | 107 | device size | implant record | no | yes |
| age, sex | 117 (redacted) / 89 | demographics | demographics | `patients` | yes |
| height, weight, BMI, BSA | 68 BMI/weight, 1 BSA | Height 6000, Weight 6005 | `CalculatedBMI` | `omr` | yes |
| diabetes | 43 | Diabetes Mellitus (history) | `Diabetes` | ICD; HbA1c 11/17 in labs | yes |
| CKD / dialysis / creatinine | 43; creatinine 17/17 | Currently on Dialysis 13880; Creatinine 6050 | `Dialysis`, `CreatLst` | `labevents` | yes |
| smoking | 62 | Tobacco Use 4625 | `TobaccoUse` | ICD / notes | yes |
| hypertension, dyslipidaemia (*not simulated*, see §7) | 85 / 80; statin 16/17 | history section | `MedLipid` | ICD, `prescriptions` | yes |
| atrial fibrillation | 46 | AFib class 13179 | `AFib` | ICD | yes |
| anticoagulant at discharge | 34; warfarin 6/17, DOAC 3/17 | Discharge meds: Warfarin, Factor Xa, thrombin inhibitor (10205 block) | `DCDirOralAnticoag`, `DCOthAnticoag` | `prescriptions` | yes |
| antiplatelet | 76; 16/17 | Aspirin, P2Y12 (discharge) | discharge meds | `prescriptions` | yes |
| LVEF | 79; 17/17 labs | LVEF 13305 (pre), 13690 (FU) | `PreLVEF`, `PostLVEF`, `HDEF` | `lvef` variants in ECHO | yes |
| pre-op AVA, indexed AVA | 32 | AV Area 13481 | — | ECHO (name unverified) | yes |
| **reference post-implant mean gradient** | 44 intra-op TEE mentions; ~5 with numeric | AV Gradient (mean) 14303 post-implant; 13675 discharge; **13676 at 30 d / 1 y** | `VDGradA` (pre-op) | ECHO serial | core lab |
| follow-up mean gradient | 62 patients with a gradient, 50 in progress notes | 13676 (30 d, 1 y) | — | ECHO serial 2008–2022 | core lab |
| EOA / DVI at follow-up | DVI 51 | AV Area 13669 (FU) | — | ECHO (LVOT VTI / AV VTI likely) | core lab |
| intraprosthetic AR, PVL grade | 19 | Central 14500, Paravalvular 14504 (FU) | — | ECHO valvular function | core lab |
| PPM flag | 17 mention | derivable from 13669 + BSA | derivable | derivable | yes |
| valve thrombosis / HALT | 2 | Valve Thrombosis Noted 13694, Leaflet Dysfunction 13695 | — | notes only | CT substudies |
| endocarditis | 9 | Endocarditis (event) | readmission | ICD | adjudicated |
| reintervention, type, indication | 11 (ViV/redo) | Aortic Valve Re-intervention (F-U) 14398, type 14354, indication 14355/14399 | `ValExpUDI`, `ReadmAortInt*` | ICD procedures | adjudicated |
| death | notes (few) | Date of Death 14315/14388 | `DischMortStat` | `patients.dod` | yes |
| symptoms / NYHA / KCCQ | 19 NYHA | NYHA 13688, KCCQ-12 | — | notes | yes |

Two things this table makes visible:

1. **The reference echo is the scarce item in real EHR data.** The provided notes have gradients for 62 patients but a *post-implant reference* value for only a handful (intra-op TEE seating comments, "pk/mn 14.5/9" in a few TEE reports). TVT captures it explicitly at discharge and 30 days; STS does not capture it post-discharge at all; MIMIC will have it only if the patient had an echo at that hospital in that window. The simulator must therefore include *missing reference* as a first-class scenario, and the model needs a fallback: the valve-class-and-size normal reference table (the JAHA review's "predicted EOA from the normal reference value for the model and size").
2. **Valve model and size are in the notes (97 and 107 of 117) and in both registries, but not in structured MIMIC.** That is the argument for the extraction step (Idea 3) being part of the deployment pipeline, not a side demo.

---

## 4. Design choices you need to make

Each changes what gets built. My recommendation is first.

**D1. Reference echo for "change from baseline".**
- (a) **30-day / discharge echo** as reference, matching the TVT fields (13675 discharge, 13676 at 30 days) and Velders' main analysis. *Recommended*: it is what real registries hold, and it makes the prototype directly portable.
- (b) 3-month echo, the strict VARC-3 wording ("1 to 3 months"). Cleaner biologically (gradients settle), but almost never exists in routine EHR.
- (c) Both, with (b) as a sensitivity analysis. Costs one extra simulated visit and one extra results column.
- Whichever you pick, include a **missing-reference arm** (say 30% of valves) handled by the valve-class-size reference table, because that is the real-world case.

**D2. How many haemodynamic channels to simulate.**
- (a) **Mean gradient, EOA and DVI, physically linked through a latent orifice area** (§2.2). *Recommended*: it is the only way to evaluate all three definitions and PPM, and it is what makes the flicker emerge naturally rather than being injected. Roughly 60 extra lines.
- (b) Mean gradient only, using the "haemodynamic" VARC-3 variant NOTION used. Half a day faster, but objective 4 collapses to one definition and Dvir cannot be evaluated.

**D3. What "truth" means for the oracle label.**
- (a) **The noise-free trajectory passed through the same VARC-3 rule.** *Recommended*: keeps truth and label on the same scale, so the gap between them is purely measurement noise plus visit timing.
- (b) The latent onset time `T_deg` itself. Cleaner mathematically, but then part of the "label error" is really the biological lag from onset to threshold, which muddies the message.

**D4. Symptom-triggered echoes.**
- (a) **Include them** (echo when symptomatic and latent stage ≥2). *Recommended*: without them the annual-surveillance arm looks worse than reality and the "urgent reintervention" metric in Idea 2 is meaningless. Adds a symptom hazard sub-model.
- (b) Guideline schedule only. Simpler, but the surveillance comparison becomes a comparison against a straw man.

**D5. Cohort composition and size.**
- (a) **10,000 valves, SAVR:TAVR 55:45, implant years 2010–2020, censor at Sept 2026**, ages as in §2.1. *Recommended*: enough events (~500 stage-2 by 10 years) for stable time-dependent AUC by stratum, and a realistic mix of short and long follow-up. Runs in seconds.
- (b) Two separate cohorts (SAVR at 70, TAVR at 79) reported side by side. Closer to how trials are published; loses the pooled valve-model-aware model.

**D6. Where the onset-time distributions come from.** *Decided Sept 16: (b), done and validated; simulator built and calibrated the same day (see §5 results).*
- (a) Parametric (the `B6` normals and `C1`/`B5` rates). Never fitted to data; describes "time to valve failure", not our label.
- (b) **Reconstruct the published curves automatically** (Guyot 2012, `D9`) and use them as calibration targets. *Chosen and implemented* in `notebooks/km_reconstruct/` (see its README): Kermen Fig 2A/2B/4A (vector PDF, zero manual input), NOTION Fig 3 (raster, risk rows typed from the figure) and Wakami Fig 1 (raster, risk row and 7 events from the figure/text). 112 validation checks against every published number pass; the primary-label curve (Kermen freedom from VARC-3 stage 2/3 SVD) reproduces its risk row 338/247/34/20 and its 49 events exactly.
- **How the simulator uses it (this is the point that makes (b) correct rather than just better-sourced).** The reconstructed curves are *observed-label* incidences (echo-detected at visits, death competing for the CIFs). They are not the latent onset distribution `T_deg`. `calibrate_onset.py` simulates a cohort matched on n, stratum, follow-up and censoring, computes the same estimator (KM with death censored for Kermen Fig 2B, Aalen–Johansen for the rest) on the *simulated observed labels*, and adjusts the latent distribution (location shift, scale multiplier) until the simulated observed curve matches `output/targets.json`. With the stand-in simulator, using the fitted curves directly as `T_deg` misses the targets by 2–25 percentage points; after calibration the residuals are 0.8–3.9 points for the ≥moderate-SVD curves (Kermen stage 2/3: 25.1 → 3.9; NOTION SAVR: 9.7 → 0.8; NOTION TAVI: 10.0 → 2.7; NOTION severe SAVR: 5.2 → 1.7; Wakami: 2.4 → 1.3). Strata without a curve (Sapien 3, Evolut, newer valves) keep the parametric form scaled to the `B4` five-year point estimates.
- One documented inconsistency: the risk row printed under Kermen Fig 2 is the all-cause survival risk set, which is incompatible with the severe-SVD curve's step sizes and its 11 reported events; that curve is reconstructed with the stage-2/3 curve's censoring and its published total is reported as informational.

**D7. Unit of analysis after reintervention.**
- (a) **The valve.** Follow-up ends at reintervention; the new valve is not followed. Matches NOTION's censoring of echo at reintervention and keeps the label clean. *Recommended.*
- (b) The patient, with the second valve continuing. More realistic for lifetime management, but doubles the simulator's complexity for no gain in the hackathon.

Defaults I will take unless you say otherwise: D1(a) with a missing-reference arm, D2(a), D3(a), D4(a), D5(a), **D6(b) (decided, implemented)**, D7(a).

---

## 5. Validation targets, restated as acceptance tests

The simulator is accepted only when, on a cohort matched to each paper's age and valve mix, it produces:

| Target | Value | Source | Reconstructed curve available? (digitised value) |
|---|---|---|---|
| Ever-labelled SVD by 5 y (Capodanno / Dvir / VARC-3), single echo | 4.6 / 2.9 / 3.0% | `A3` | no (point estimates only) |
| Label persistence at the next visit | 14–50%, typically ~33% | `A3` Table 2 | no |
| Regression to the mean by decile of MPG | lowest +1.2 to +2.3; highest −1.0 to −5.9 mmHg per interval | `A3` Table E4 | no |
| Within-patient 5-y change in MPG, 95% interval | −9.6 to +7.5 mmHg | `A3` | no |
| ≥moderate / severe SVD at 10 y, age 79, TAVI vs SAVR | 15.4 / 1.5% vs 20.8 / 10.0% | `B1` | **yes, whole curves** (15.4 / 1.4 vs 20.7 / 9.9), risk rows exact |
| Freedom from stage 2/3 at 5 and 10 y, age 70, Perimount-class | 98.5%, 61% | `A4` | **yes, whole curve + censor marks** (98.5, 60.9; 49 events exact) |
| Freedom from stage 3 at 5 and 10 y | 99.6%, 88.3% | `A4` | yes (100.0, 87.2); see D6 note |
| Explant for SVD at 10 y, competing-risk | 8.0% | `A4` | yes (8.4) |
| Trifecta-class SVD at 5 and 7 y | 4.8%, 6.6% | `C1` | **yes, whole curve** (4.8, 6.3; 7 events exact) |
| Death by 10 y at age 79 / 70 | ~63% / 33% | `B1`, `A4` | Kermen yes (overall survival 66.8 at 10 y; valve-related / non-valve death CIFs 7.4 / 25.4) |

For the rows marked "yes" the acceptance test is no longer two points but the whole curve on a yearly grid, risk-weighted, computed with the paper's own estimator on the simulated observed labels (`notebooks/simulator/calibrate.py`, `validate.py`; targets in `notebooks/km_reconstruct/output/targets.json`).

**Result (Sept 16, simulator built and calibrated; `notebooks/simulator/output/calibration_report.md` and `validation_report.md`).** Every parameter carries a source or an assumption flag (`params.py`). Calibration is sequential: physics constants per valve model and size from the ASE 2024 normal-value tables (`A5`), observation noise on a PERIGON-matched cohort against `A3`, mortality on Kermen's survival curve, then the latent onset Weibull (scale, shape) per valve class on the reconstructed curves with each paper's own single-echo label rule and estimator. Runs in ~4 min on 32 processes.

| Target | Published | Simulated | Verdict |
|---|---|---|---|
| Kermen freedom from stage 2/3 at 5 / 10 y (curve, KM) | 97.0 / 63.8 (digitised) | 96.9 / 62.9; max deviation 1.3–2.1 pp | pass |
| Wakami Trifecta SVD (curve, AJ) | 4.8 / 6.3 at 5 / 7 y | max deviation 2.0 pp | pass |
| NOTION TAVI ≥moderate SVD at 10 y (curve, AJ) | 15.4 | 13.2; max deviation 2.9–3.0 pp | pass |
| NOTION SAVR ≥moderate SVD at 10 y, **predicted with no free parameter** (Perimount + Trifecta transported to 79 y by the age HR) | 20.7 | 21.6 at 10 y (7.7 vs 10.7 at 5 y); 3.1 pp | pass (cross-check) |
| NOTION SAVR after freeing the porcine/Mitroflow share | 10.7 / 20.7 | 9.7 / 22.9 | pass |
| Kermen stage 3 at 10 y (check on progression speed) | 87.2 | 85.5 | pass |
| NOTION severe SVD at 10 y, SAVR / TAVI (checks) | 9.9 / 1.4 | 13.4 / 4.6; 3.5–4.1 pp | marginal miss, accepted |
| Kermen explant for SVD at 10 y (AJ) | 8.4 | 7.6 | pass |
| Kermen overall survival (curve) / NOTION mortality at 10 y | 66.8 / 63 | 68.3 / 61.6 | pass |
| Sapien XT / Sapien 3 VARC-3 SVD at 5 y | 1.6 / 0.7 | 1.6 / 2.1 (single-echo noise floor ~1–2%) | pass / check |
| Velders ever-labelled by 5 y, Capodanno / Dvir / VARC-3 | 4.6 / 2.9 / 3.0 | 5.2 / 3.6 / 4.2 | pass |
| Velders label persistence at the next visit | ~33 (14–50) | 20 | pass (low end) |
| Velders regression to the mean, lowest / highest decile | +1.2…+2.3 / −1.0…−5.9 | +0.8 / −1.2 | pass |
| Velders 95% PI of 5-y within-patient change | (−9.6, +7.5) | (−9.6, +10.3) | pass |
| Kermen reference echo: EOA / EOAi / moderate PPM | 1.6 / 0.87 / 51% | 1.58 / 0.81 / 40% | pass / pass / marginal |
| Emergent Cox HR on the confirmed label: age per year / severe PPM / diabetes / female (null) | 0.951 / 1.85 / 1.33 / 1.0 | 0.98 / 1.70 / 1.30 / 0.90 | pass |
| Emergent Cox HR: reference MPG ≥15 / THV ≤23 mm / smoking | 1.30 / 2.07 / 2.58 | 0.91 / 1.07 / 1.79 | **short**; see note |

Note on the last row (a finding, not only a residual): effects that act through the reference gradient are attenuated, and can even flip sign, on *observed* labels because a high measured reference gradient is partly a high flow state at the 30-day echo and regresses to the mean, so later increases are smaller. Label noise dilutes every effect (smoking 3.5 latent → 1.8 observed). The protocol should state that published risk-factor hazard ratios on single-echo labels are not the biological effects, and the prototype's model comparison under the three labels and the oracle makes this visible.

Generated cohort (`data/synthetic/`, seed 1): 10,000 valves, 50,849 echoes, 4,009 deaths, 1,589 reinterventions (52% urgent), 572 endocarditis; confirmed VARC-3 stage ≥2 cumulative incidence 6.2 / 11.3 / 16.8% at 5 / 8 / 10 y (1,641 events; 2,504 single-echo labels of which 34% never confirmed).

If a target is missed by more than its published confidence interval, the offending module is tuned before any model is trained, and the tuning is logged in the notebook. That log is part of the deliverable.

---

## 6. Revisions after the "2.1 modifications" clinical review (Sept 16)

A collaborator reviewed §2.1 variable by variable. Outcome of that review:

**Accepted and changed above**
- BSA and BMI are now derived from jointly drawn height and weight, not drawn independently. My original draw could produce impossible body-size pairs. Correct criticism.
- Atrial fibrillation is an explicit variable that *generates* anticoagulation, so the model has to handle confounding by indication instead of seeing anticoagulation as a free-standing flag. Coronary disease was accepted at first and then removed in the audit in §7: it has no quantified effect on SVD or, in our sources, on post-AVR death.
- Anticoagulation is kept but restructured: type (VKA / DOAC / none) and time-specific effects, protective early after TAVR through the thrombosis pathway, modestly harmful long-term for VKA after SAVR through the calcification pathway. See the rejection note below for why it is not removed.
- Peak velocity, stroke-volume index, baseline intraprosthetic AR grade and baseline paravalvular grade are added as reference-echo channels. PVL is a covariate only.
- Valve strata are named after real models with a class parent, so "exact model where available, class otherwise" is exactly what the hierarchical model does.
- A transient thrombosis/HALT module is added (default on). Previously "not simulated"; the review's point that a thrombosis-driven gradient rise must not be labelled SVD is only testable if such rises exist in the data.
- Features for the landmark model now explicitly include time since last echo and the number of consecutive worsening echoes, alongside the per-year deltas of gradient, EOA, DVI and Vmax.

**Rejected, with reasons**
- *Prediction horizon of 3 years.* Almost no events exist that early: freedom from VARC-3 stage 2/3 is 98.5% at 5 years (`A4`), per-visit SVD is 0.2–2% in the first 5 years (`A3` Table E2), severe SVD does not appear before year 5 in Perimount-class valves. A 3-year model would be estimated on a few dozen events out of 10,000 and would be uninformative. Horizons stay 5, 8 and 10 years; 3 years is reported descriptively only.
- *"Time to first incident SVD" as the primary outcome.* Taken literally, the first threshold crossing is the label that flickers: 59–65% of first classifications are absent at the next visit (`A3`). The primary outcome stays **confirmed on two consecutive echoes**, which VARC-3 and the 2017 consensus already ask for (`A1` point 4). The review's own "number of consecutive worsening measurements" feature is the same idea applied on the predictor side.
- *Remove anticoagulation as a predictor because the relationship is "arbitrary".* It is not arbitrary; it is one of the largest published effects: absence of anticoagulation at discharge was an independent predictor of haemodynamic deterioration in a 1,521-patient TAVR registry (HR 3.35), anticoagulation halved dysfunction at 3 years in another (OR 0.54), and Salaun found warfarin use among the factors for *late* deterioration after SAVR (`B4` Table 2 and text). What is wrong is a single naive arrow; the fix is the type- and time-specific structure above, generated through the AF indication. It is also in every real schema we anchor to (TVT discharge medications, STS `DCDirOralAnticoag`).
- *Add calcium / phosphate / PTH, haemoglobin, LV mass, heart-failure hospitalisation to the simulator.* No published SVD hazard exists for any of them in our references, so simulating them can only inject invented structure. They belong in the real-study feature list, not the simulator. Heart-failure hospitalisation is additionally outcome-adjacent (it is a VARC-3 stage-1 failure criterion), so using it as a predictor risks the leakage the review itself warns about.
- *"Add time since implantation as a predictor."* In a time-to-event model it is the time axis, not a covariate; in the landmark model the landmark time already plays that role. Nothing to add.
- *Leakage framing.* The rule "use only information available at prediction time" is right and is how the landmark design works. One correction to the example: the reference echo legitimately appears on both sides, as the anchor of the predictor deltas and of the label's change-from-baseline. That is shared conditioning, not leakage. Leakage would be using the year-4 echo that defines a year-4 event to predict it.

---

## 7. Variable audit: effect on the target, and whether clinicians actually collect it

Two admission tests for every variable in the simulator. (1) *Effect*: is there a quantified effect on SVD in our references, or an explicit non-predictor role (label component, physiological link, confounder of a quantified effect, competing-event driver)? (2) *Collected*: is it acquired under the ESC 2021 follow-up recommendation (`A2`: baseline echo with transprosthetic gradients within 30 days, at 1 year, annually), required by the VARC-3 haemodynamic definitions (`A1`, `B1`), or a named field on the TVT v3.0 or STS v4.20.2 forms? A variable that fails (1) or (2) is out of the simulator, whatever its biological appeal. "Assumption" means a modelling choice with no published number; each is listed so it can be varied in a sensitivity run.

| Variable | Effect on SVD (direction, size, source) | Collected under | Verdict |
|---|---|---|---|
| age at implant | strong protective per year: severe SVD HR 0.907/yr (`A4` Table 3); explant for SVD 5.6/20/45% at 10/15/20 y under 60 vs 1.5/5.1/8.1% at 60–80 (`B5`) | universal | **keep** |
| approach, valve model | TAVI vs SAVR severe SVD HR 0.2 at 10 y (`B1`); model-specific rates: Sapien XT 1.6% vs Sapien 3 0.7% at 5 y, self-expanding 2.2% vs SAVR 4.4% (`B4` Table 1); Trifecta 4.8% at 5 y (`C1`); Perimount explant 1.9% at 10 y (`B5`); Mitroflow accelerated (`A1` Table 1) | STS `AorticImplant` / UDI; TVT device fields | **keep** |
| label size | TAVR ≤23 mm HR 2.07 (`B4`); for SAVR, size is *not* independent once the post-implant gradient is in the model (`B5`; `C1` 19 mm NS) → acts only through EOA/PPM/gradient | STS implant record; TVT device size | **keep**, effect routed through gradient for SAVR |
| height, weight → BSA, BMI | BMI HR 1.08 per kg/m² (`B4`); BSA needed for EOAi/PPM | TVT Height 6000 / Weight 6005; STS `CalculatedBMI` | **keep** |
| sex | contradictory: male HR 2.17 in one study, female a late-HVD factor in Salaun (`B4`) → simulated with **no effect**; used for the fairness check | universal | keep as null covariate |
| diabetes | HR 1.33 (`B4` Table 2, from Salaun early HVD) | TVT history; STS `Diabetes` | **keep** |
| CKD / dialysis | HR 1.10 (`B4` Table 2), 1.29–1.34 point estimates in `A4`; renal insufficiency an early-HVD factor in Salaun | TVT Dialysis 13880 / Creatinine 6050; STS `Dialysis`, `CreatLst` | **keep** at HR ≈1.3 for eGFR <60 |
| smoking | HR 2.58 (`B4` Table 2); only clinical univariable predictor in the PET study (`B3`) | TVT Tobacco Use 4625; STS `TobaccoUse` | **keep** |
| hypertension | one estimate only, severe SVD HR 6.1 with CI 1.1–33 on 11 events (`A4`); absent from the pooled predictor table (`B4`) | TVT history | **drop** from simulator; prevalence 76% carries little information anyway |
| dyslipidaemia (binary) | no quantified hazard for the binary in our set; `A4` NS; `B5` lipid-lowering drugs associated with explant but flagged as confounded; metabolic-syndrome papers (Briand, Mahjoub) paywalled and unverified; the quantified "dysmetabolic" effects are biomarkers (PCSK9, Lp-PLA2, HOMA) nobody collects routinely | STS `MedLipid` (statin) only | **drop** from simulator; list metabolic syndrome as a real-study candidate |
| coronary disease | no SVD effect in any source; no verified post-AVR mortality HR in our set | TVT prior PCI/CABG; STS coronary fields | **drop** |
| atrial fibrillation | no direct SVD effect; it is the *indication* for anticoagulation, which has a quantified effect, so it must exist for the confounding to be real | TVT AFib class 13179; STS `AFib` | keep as confounder |
| anticoagulation (type, timing) | absence at discharge → HVD HR 3.35 after TAVR (`B4` Table 2, 1,521-patient registry); anticoagulation → BVD OR 0.54 at 3 y (`B4` ref 86); warfarin a *late*-HVD factor after SAVR (`B4`, Salaun); HALT 12.4% vs 32.4% with vs without (`A2`) | TVT discharge medications block 10205; STS `DCDirOralAnticoag`, `DCOthAnticoag` | **keep**, type- and time-specific |
| implant year | no effect of its own; carries device generation (XT → S3) and defines censoring | universal | keep as design variable |
| reference mean gradient | the central quantified predictor: ≥15 mmHg HR 1.30 (`B4`, Salaun); continuous in `B5` Table 2, a 10 mmHg higher post-op peak gradient more than doubles 20-y explant risk under 60; 13.9 vs 9.5 mmHg in SVD vs no SVD (`C1`) | ESC 2021 baseline echo (`A2`); VARC-3; TVT 13675 / 13676 | **keep** |
| reference EOA, EOAi, PPM | severe PPM HR 1.85 after SAVR, 1.70 after TAVR (`B4` Table 2); PPM HR 7.7 in Trifecta (`C1`); Flameng 2010; `B5` shows gradient dominates size | VARC-3 requires EOA; TVT AV Area 13495/13669; BSA from height/weight | **keep** |
| reference DVI | no separate hazard (collinear with EOA); required by the VARC-3 and Dvir definitions, missing 5–8% at baseline (`A3` Table E1) | VARC-3; ASE prosthetic-valve echo | keep as label channel |
| peak velocity (Vmax) | peak gradient effect is quantified (`B5` Table 2); peak = 4·Vmax², so this is the same effect as the gradient, not an extra one | TVT AV Peak Velocity 13703; standard report item | keep as **derived** channel, no separate hazard (avoids double counting) |
| stroke-volume index, LVEF | no SVD hazard; they set the flow that turns EOA into gradient, which is the mechanism of gradient flicker and of low-flow low-gradient states | TVT Low Flow SVI <35 field 13700; LVEF 13305/13690; STS `PreLVEF`/`PostLVEF` | keep as **physiological link**, LVEF also in the death model (HR 0.974/%) |
| reference intraprosthetic AR grade | label component ("new or worsening"); ≥mild transprosthetic regurgitation is an early-HVD factor in Salaun, direction only, no HR available | VARC-3; TVT Central Regurgitation 14499/14500 | keep as label anchor; **hazard set to null**, flagged as assumption |
| reference paravalvular grade (TAVR) | none: PVL is non-structural, not associated with mortality in NOTION (`B1`), no SVD link | VARC-3; TVT Paravalvular 14503/14504 | optional null covariate, for the "PVL ≠ SVD" adjudication demo only; first to drop if time is short |
| calcium, phosphate, PTH, haemoglobin, LV mass | no SVD hazard in any source; not part of guideline follow-up | not on TVT/STS forms except Hb (13775), which has no link | **drop** |
| heart-failure hospitalisation | it is a VARC-3 stage-1 *failure criterion*, so it is outcome, not predictor; `A2` says HF after valve surgery should trigger a search for SVD | events, not baseline | **drop** as predictor |
| CT anatomy, commissural alignment | misalignment HR 1.03 per 10° (`B4`), but CT is not routine follow-up and the effect is small | not routine | **drop** (secondary model only) |
| transient thrombosis / HALT module | HALT prevalence 12.4% vs 32.4% at 3 months (`A2`), SLT 6–15% up to 31% at 1 y, resolves with 3–6 months anticoagulation, associated with later SVD (`B4`) | the *episode* needs CT; its *gradient effect* appears on routine echo; TVT "Valve Thrombosis Noted" 13694 | keep as latent process; magnitude of the bump (+8–15 mmHg) is an **assumption** |
| endocarditis | competing event, 7% at 10 y (`B1`), removes the valve from the SVD label | TVT event; ICD | keep |
| death | age HR 1.033/yr, LVEF 0.974/% (`A4` Table E2); dialysis multiplier | universal | keep |
| symptom onset (drives symptom-triggered echo and reintervention) | no published hazard as a function of gradient; shape is an **assumption** calibrated only by the fraction of deterioration that is subclinical (most of 31% in Salaun; `B4`) | NYHA 13688, KCCQ on TVT | keep, flagged for sensitivity analysis |
| echo missingness rising with time, and higher when well | rates from `A3` Table E1; the "higher when well" dependence is an **assumption** | n/a | keep, flagged |

Net effect of the audit: three variables the review asked to keep or add are removed from the simulator (hypertension, dyslipidaemia, coronary disease), five optional ones it proposed never enter (calcium/phosphate/PTH, haemoglobin, LV mass, heart-failure admission, CT), and every remaining variable is one of: a quantified predictor, a label component, a physiological link between them, a confounder of a quantified effect, or a driver of a competing event. Four modelling assumptions are named and will each get a sensitivity run: the AR-grade hazard (null), the thrombosis bump size, the symptom hazard shape, and the well-patient missingness dependence.
