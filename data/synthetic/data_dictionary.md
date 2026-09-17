# Synthetic cohort data dictionary

Every column exists in a named real-world source (docs/04 section 3); values are simulated, the schema is not.


## `patients_valves`

| column(s) | meaning / real-world field |
|---|---|
| `valve_id` | synthetic key (one row per implanted valve; D7: the valve is the unit of analysis) |
| `approach` | SAVR / TAVR (TVT DCF procedure; STS ACSD `VSAVPr`) |
| `valve_model` | named prosthesis (STS `AorticImplant` / UDI `VSAoImUDI`; TVT device fields) |
| `valve_class` | calibration stratum (Perimount-class, Trifecta-class, porcine/Mitroflow, balloon-expandable, self-expanding, new) |
| `is_new_model` | 1 if the model has little published history (partial pooling target) |
| `implant_year / implant_date` | index procedure (STS/TVT procedure date; shifted dates in the provided notes) |
| `age_at_implant, sex, height_cm, weight_kg, bsa, bmi` | TVT Height 6000 / Weight 6005; STS `CalculatedBMI`; BSA by Mosteller |
| `diabetes, ckd, dialysis, smoking, af` | TVT history fields (Dialysis 13880, Tobacco 4625, AFib 13179); STS `Diabetes`, `Dialysis`, `TobaccoUse`, `AFib` |
| `anticoag_discharge, anticoag_type, anticoag_year1, anticoag_longterm` | TVT discharge medications block 10205; STS `DCDirOralAnticoag` |
| `lvef` | TVT LVEF 13305/13690; STS `PreLVEF`/`PostLVEF`; MIMIC-IV-ECHO measurement |
| `label_size_mm` | labelled prosthesis size (STS implant record; TVT device size) |
| `reference_source` | 30d_echo (TVT 13675/13676; baseline echo, ESC/EACTS 2025 window within 3 months) or model_size_table (D1 missing-reference arm) |
| `ref_mean_gradient_mmHg, ref_av_area_cm2, ref_dvi, ref_ar_grade, ref_eoai_cm2_m2, ppm_class_ref` | reference (30-day) echo values; TVT AV Mean Gradient 13675, AV Area 13495, Central AR 14499; PPM per ASE 2024 Table 7 |
| `table_eoa, table_mpg` | ASE 2024 normal values for the model x size (reference when the 30-day echo is missing) |
| `followup_end_years` | end of valve follow-up (death, reintervention, lost, administrative censoring) |

## `echo_visits`

| column(s) | meaning / real-world field |
|---|---|
| `valve_id, visit_time_years, visit_date` | one row per measurement per echo (long format as MIMIC-IV-ECHO `structured_measurement`) |
| `visit_type` | reference | scheduled | symptom_triggered | other_indication |
| `measurement` | av_mean_gradient_mmHg (TVT 13674-13676), av_peak_velocity_m_s (13703), av_area_cm2 (13481/13495/13669), doppler_velocity_index, stroke_volume_index_ml_m2 (Low Flow SVI 13700), lvef_pct, intraprosthetic_ar_grade (14499/14500), paravalvular_leak_grade (14503/14504) |
| `result` | value as reported (rounded as in clinical reports; EOA/DVI missing in a share of echoes as in PERIGON) |

## `events`

| column(s) | meaning / real-world field |
|---|---|
| `event_type` | death | endocarditis | reintervention | lost_to_follow_up | administrative_censoring (TVT follow-up form; STS `ValExpUDI`) |
| `detail` | reintervention: type|indication|urgency (ViV / redo; severe SVD, moderate SVD with symptoms, endocarditis); death: cause class |

## `labels`

| column(s) | meaning / real-world field |
|---|---|
| `<definition>_single2 / _single3` | first echo meeting the definition (Capodanno, Dvir, VARC-3 full, VARC-3 haemodynamic, oracle) at stage >=2 / stage 3, years |
| `<definition>_conf2 / _conf3` | first echo of a pair of consecutive echoes both meeting the stage (persistence-confirmed) |
| `primary_event, primary_time, primary_competing` | PRIMARY ENDPOINT: confirmed VARC-3 stage >=2 (time = first echo of the confirmed pair); competing 2 = death, 3 = reintervention |

## `latent_truth`

| column(s) | meaning / real-world field |
|---|---|
| `*` | NOT FOR TRAINING. Latent onset T_deg, mode, progression tau, thrombosis episodes, noise-free crossing times, symptom onset; used only for the oracle and the surveillance-delay metrics |

## Parameter provenance

| parameter | source | assumption |
|---|---|---|
| `cohort.*` | docs/04 D5: 10,000 valves, SAVR:TAVR 55:45, implants 2010-2020, censor Sept 2026 |  |
| `baseline.age_*` | A3 Table 1 (70.2 +/- 9.0), B1 (79.1 +/- 4.8), B5 |  |
| `baseline.height/bmi` | A3 Table 1: BSA 2.0 +/- 0.2, BMI 29.4 +/- 5.4 |  |
| `baseline.diabetes/ckd/smoking` | A3 Table 1: 27% / 11% renal / smoking prevalence ASSUMED 15% | yes |
| `baseline.af, anticoag` | B1 new-onset AF 52-74%; A4 5.6% at implant; anticoagulation via AF ASSUMED shares | yes |
| `baseline.lvef` | A3 Table 1: 59 +/- 10 |  |
| `baseline.svi, ejection_time` | physiological ranges; low-flow < 35 mL/m2 ~16% (ASSUMPTION) | yes |
| `baseline.reference_missing_frac` | docs/04 D1(a) missing-reference arm 30% (ASSUMPTION; ~5/117 numeric TEE gradients in notes) | yes |
| `baseline.ar_ref/pvl_ref` | B1 (PVL 18-53% TAVI, 5.2% SAVR; intraprosthetic 5.8 / 2.2%), A3 (0.2% >= moderate at 5 y) |  |
| `physics.c_model, A_eff` | A5 Appendix Tables A1/A2/A4 (fitted); Gorlin 44.3 constant |  |
| `physics.peak_over_mean` | A5 Table A4 (Perimount 19.9/11.5, Mosaic 23.8/13.7) |  |
| `onset.shape/scale` | CALIBRATED to km_reconstruct/output/targets.json (Kermen, Wakami, NOTION) and B4 point targets |  |
| `onset.hr_age_per_year_savr` | A4 Table 3: 0.951 (stage 2/3) |  |
| `onset.hr_age_per_year_thv` | ASSUMPTION (B4 Table 2: inconclusive) | yes |
| `onset.hr_smoking/diabetes/bmi/ckd` | B4 Table 2: 2.58 / 1.33 / 1.08 per unit / 1.10 (A4 1.29-1.34) |  |
| `onset.hr_no_anticoag_year1_tavr` | biological residual of B4 3.35 (rest through thrombosis module) | yes |
| `onset.hr_vka_longterm_savr` | ASSUMPTION from Salaun late-HVD factor (B4) | yes |
| `onset.hr_severe_ppm_residual / hr_thv_le23_residual` | residuals; totals 1.85 / 2.07 (B4) are emergent checks | yes |
| `onset.p_regurgitant` | C1: 5/7 Trifecta failures were cusp tears; others ASSUMED | yes |
| `trajectory.drift_per_year` | B5: mean gradient 14 -> 20 mmHg by year 15 |  |
| `trajectory.tau_*` | ASSUMPTION; checked against Kermen stage 3 and NOTION severe curves | yes |
| `thrombosis.p_year1, anticoag_factor, duration` | A2 (HALT 12.4% vs 32.4%), B4 (SLT 6-15%, resolves in 3-6 months) |  |
| `thrombosis.eoa_factor` | ASSUMPTION (no published gradient effect) | yes |
| `noise.sigma_*` | CALIBRATED to A3 Table E4 (decile regression to the mean) and 5-y change PI (-9.6, +7.5) |  |
| `noise.eoa_missing/dvi_missing` | A3 Table E1 |  |
| `visits.miss_by_year` | A3 Table E1 (9.3 / 14.5 / 22.4 / 35.0 %) |  |
| `visits.miss_saturation, miss_mult_*` | ASSUMPTION | yes |
| `visits.symptom_*` | ASSUMPTION (most stage-2 disease silent) | yes |
| `death.a` | CALIBRATED to A4 Fig 2A overall survival |  |
| `death.b_per_year` | population Gompertz slope ~9%/yr (ASSUMPTION; A4 adjusted within-cohort HR 1.033) | yes |
| `death.hr_lvef_per_point` | A4 Table E2: 0.974 |  |
| `death.hr_dialysis, hr_tavr_era` | ASSUMPTION (NOTION 63% dead at 10 y is the check) | yes |
| `events.endocarditis_rate` | B1: 7.2-7.4% at 10 y |  |
| `events.reint_*` | CALIBRATED to A4 explant CIF 8.4% at 10 y, 9/11 stage-3 reintervened; B1 2.2-4.3% |  |
| `events.viv_logistic` | B7: ViV 79.0 +/- 5.3 vs redo 70.6 +/- 8.8 y |  |
| `events.ltfu_rate` | A4: 2% lost over 6.6 y |  |
