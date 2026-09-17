# Simulator calibration report

37 checks, 7 FAIL

| step | quantity | target | simulated | tolerance | result | note |
|---|---|---|---|---|---|---|
| reference echo | Kermen Magna Ease: discharge MPG (mmHg) | 12.6 +/- 3.0 | 14.3 +/- 5.9 | mean +/-2 | PASS |  |
| reference echo | Kermen: EOA (cm2) | 1.6 +/- 0.3 | 1.58 +/- 0.23 | mean +/-0.15 | PASS |  |
| reference echo | Kermen: EOAi (cm2/m2) | 0.87 +/- 0.17 | 0.81 +/- 0.14 | mean +/-0.08 | PASS |  |
| reference echo | Kermen: moderate / severe PPM (%) | 51 / 2 | 40 / 5 | +/-10 / +/-4 | FAIL |  |
| reference echo | Velders (Avalus ~ Perimount-class): MPG / EOA / DVI at discharge | 13.1 +/- 4.7 / 1.54 +/- 0.36 / 0.49 +/- 0.10 | 13.9 +/- 5.7 / 1.67 +/- 0.36 / 0.54 +/- 0.13 | MPG +/-2, EOA +/-0.2, DVI +/-0.08 | PASS |  |
| reference echo | per model x size mean gradient vs ASE 2024 rows (mmHg), worst |diff| |  | 3.82 | 3.5 | FAIL | 52 model-size cells with n>=20 |
| label flicker (A3) | ever labelled by 5 y, Capodanno / Dvir / VARC-3 (%) | 4.6 / 2.9 / 3.0 | 5.2 / 3.6 / 4.2 | +/-1.5 | PASS |  |
| label flicker (A3) | VARC-3 persistence at the next visit (%) | ~33 (published range 14-50) | 20.0 | 14-50 | PASS | next visit missing 6% |
| label flicker (A3) | regression to the mean: lowest / highest decile change (mmHg per interval) | +1.2..+2.3 / -1.0..-5.9 | +0.8 / -1.2 |  | PASS |  |
| label flicker (A3) | 95% PI of within-patient 5-y change (mmHg) | (-9.6, +7.5) | (-9.6, 10.3) | width +/-4 | PASS |  |
| label flicker (A3) | per-visit VARC-3 prevalence, years 1-5 (%) | 1.0 / 1.1 / 1.4 / 1.4 / 0.5 | 0.6/0.6/0.5/0.8/1.6 | 0.2-3 | PASS |  |
| curve targets | kermen_fig2b_freedom_stage23: max |sim - target| while n>=30 (pp) |  | 2.06 | 3.0 | PASS | at 10 y sim 62.0 vs 63.8 |
| curve targets | wakami_fig1_svd: max |sim - target| while n>=30 (pp) |  | 2.01 | 2.5 | PASS | at 7 y sim 9.1 vs 6.3 |
| curve targets | notion_fig3_modsvd_tavi: max |sim - target| while n>=30 (pp) |  | 3.0 | 3.0 | PASS | at 10 y sim 13.2 vs 15.4 |
| curve targets | notion_fig3_modsvd_savr: max |sim - target| (pp) |  | 2.66 | 3.0 | PASS |  |
| curve checks | kermen_fig2b_freedom_stage3: max |sim - target| while n>=30 (pp) |  | 2.18 | 5.0 | PASS |  |
| curve checks | notion_fig3_sevsvd_savr: max |sim - target| while n>=30 (pp) |  | 4.06 | 3.0 | FAIL |  |
| curve checks | notion_fig3_sevsvd_tavi: max |sim - target| while n>=30 (pp) |  | 3.83 | 3.0 | FAIL |  |
| curve checks | kermen_fig4a_explant_svd: max |sim - target| while n>=30 (pp) |  | 0.96 | 3.0 | PASS |  |
| death | Kermen overall survival 1-10 y, max |sim - target| (pp) |  | 2.13 | 3.0 | PASS | 10 y sim 68.3 vs 66.8 |
| reintervention | Kermen: share of detected stage-3 valves reintervened (%) | 9/11 = 82 | 74.0 | >= 55 | PASS |  |
| death | NOTION TAVI all-cause mortality at 10 y (%) | 62.7 | 61.6 | 6.0 | PASS |  |
| reintervention | NOTION TAVI reintervention by 10 y (%) | 4.3 | 1.21 | 1-7 | PASS |  |
| emergent HR (Cox on confirmed VARC-3 >=2) | age per year | 0.951 | 0.98 | 0.93-0.98 | PASS |  |
| emergent HR (Cox on confirmed VARC-3 >=2) | severe PPM (total effect) | 1.85 | 1.7 | 1.3-2.8 | PASS |  |
| emergent HR (Cox on confirmed VARC-3 >=2) | reference MPG >= 15 mmHg | 1.3 | 0.91 | 1.05-1.8 | FAIL |  |
| emergent HR (Cox on confirmed VARC-3 >=2) | THV <= 23 mm | 2.07 | 1.07 | 1.3-3.2 | FAIL |  |
| emergent HR (Cox on confirmed VARC-3 >=2) | female sex (null) | 1.0 | 0.9 | 0.85-1.18 | PASS |  |
| emergent HR (Cox on confirmed VARC-3 >=2) | smoking | 2.58 | 1.79 | 1.9-3.4 | FAIL |  |
| emergent HR (Cox on confirmed VARC-3 >=2) | diabetes | 1.33 | 1.3 | 1.05-1.7 | PASS |  |
| dataset | valves / echoes / events |  | 10000 / 50849 / 10572 |  | PASS |  |
| dataset | primary endpoint (confirmed VARC-3 >=2) cumulative incidence at 5 / 8 / 10 y (%) |  | 6.2/11.3/16.8 |  | PASS | 1641 events; single-echo labels 2504 |
| dataset | share of single-echo VARC-3 labels never confirmed (%) | Velders: 59-65 absent at next visit | 34.0 |  | PASS |  |
| dataset | echoes per valve (mean) / follow-up years (mean) |  | 5.1 / 7.2 |  | PASS |  |
| dataset | events: death 4009, administrative_censoring 3715, reintervention 1589, lost_to_follow_up 687, endocarditis 572 |  |  |  | PASS |  |
| dataset | reinterventions urgent share (%) |  | 52.0 |  | PASS |  |
| dataset | valve models: Perimount Magna Ease 2435, Trifecta 1414, Sapien XT 1397, Sapien 3 1115, Evolut R 984, CoreValve 846, Inspiris Resilia 552, Mitroflow 318, Epic 265, New SAVR valve 238, Mosaic 219, New THV 217 |  |  |  | PASS | provided notes: ~10 models, 3 THV generations, sizes 19-29 (docs/01) |
