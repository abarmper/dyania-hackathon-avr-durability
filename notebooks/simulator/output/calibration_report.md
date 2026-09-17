# Simulator calibration report

40 checks, 2 FAIL

| step | quantity | target | simulated | tolerance | result | note |
|---|---|---|---|---|---|---|
| physics | c fitted per model x size to the ASE rows (x1e-4; Gorlin = 5.1) |  | range 1.88-5.45 |  | PASS | a single per-model constant would leave residuals up to 11.2 mmHg at the smallest sizes |
| noise (pass 1, prior onset) | sigma_flow / sigma_cw chosen |  | 0.055 / 0.03 |  | PASS | informational |
| noise (pass 1, prior onset) | lowest-decile mean change (mmHg/interval) | +1.2..+2.3 | 0.79 | +/-1 | PASS | informational |
| noise (pass 1, prior onset) | highest-decile mean change | -1.0..-5.9 | -0.79 | +/-1 | PASS | informational |
| noise (pass 1, prior onset) | middle-decile SD of change | 2.4..3.9 | 2.49 | +/-0.7 | PASS | informational |
| noise (pass 1, prior onset) | 95% PI of 5-y within-patient change | (-9.6, +7.5), width 17.1 | (-8.5, 9.5), width 18.0 | width +/-4 | PASS | informational |
| noise (pass 1, prior onset) | ever labelled by 5 y, Capodanno / Dvir / VARC-3 (%) | 4.6 / 2.9 / 3.0 | 4.7 / 3.2 / 5.8 | +/-1.5 | PASS | informational |
| noise (pass 1, prior onset) | VARC-3 label persistence at next visit (%) | ~33 (published range 14-50) | 43.0 | 14-50 | PASS | informational; next visit missing 9% (Velders 15-25%); per-visit prevalence 0.6/0.6/0.6/1.0/2.6 (Velders 0.2-1.4) |
| death | Kermen overall survival, max |sim - target| 1-10 y (pp) |  | 1.99 | 3.0 | PASS | a = 0.0315; sim 5/10 y = 83.1/67.6, target 81.1/66.8 |
| death | NOTION all-cause mortality at 10 y (%), age 79 | 63.0 | 62.3 | 6.0 | PASS | check, not fitted |
| onset | kermen_fig2b_freedom_stage23 (perimount), max |sim - target| while n>=30, <=10 y (pp) |  | 1.1 | 3.0 | PASS | scale x1.00 -> 10.5 y, shape x1.40 -> 5.60; prior residual 8.6 pp; 50s |
| onset | wakami_fig1_svd (trifecta), max |sim - target| while n>=30, <=10 y (pp) |  | 2.0 | 2.5 | PASS | scale x1.44 -> 31.6 y, shape x0.58 -> 0.92; prior residual 3.8 pp; 62s |
| onset | notion_fig3_modsvd_tavi (self_expanding), max |sim - target| while n>=30, <=10 y (pp) |  | 2.77 | 3.0 | PASS | scale x0.44 -> 7.0 y, shape x1.00 -> 2.50; prior residual 8.8 pp; 36s |
| noise (pass 2) | sigma_flow / sigma_cw chosen |  | 0.055 / 0.045 |  | PASS |  |
| noise (pass 2) | lowest-decile mean change (mmHg/interval) | +1.2..+2.3 | 0.96 | +/-1 | PASS |  |
| noise (pass 2) | highest-decile mean change | -1.0..-5.9 | -1.37 | +/-1 | PASS |  |
| noise (pass 2) | middle-decile SD of change | 2.4..3.9 | 2.77 | +/-0.7 | PASS |  |
| noise (pass 2) | 95% PI of 5-y within-patient change | (-9.6, +7.5), width 17.1 | (-9.1, 9.3), width 18.4 | width +/-4 | PASS |  |
| noise (pass 2) | ever labelled by 5 y, Capodanno / Dvir / VARC-3 (%) | 4.6 / 2.9 / 3.0 | 4.3 / 2.9 / 4.0 | +/-1.5 | PASS |  |
| noise (pass 2) | VARC-3 label persistence at next visit (%) | ~33 (published range 14-50) | 20.0 | 14-50 | PASS | next visit missing 7% (Velders 15-25%); per-visit prevalence 0.6/0.6/0.5/0.7/1.4 (Velders 0.2-1.4) |
| onset (pass 2) | kermen_fig2b_freedom_stage23 (perimount), max |sim - target| while n>=30 (pp) |  | 1.42 | 3.0 | PASS | scale 10.49 y, shape 5.60; sim 5/10 y = 97.0/62.9 vs 97.0/63.8 |
| cross-check | NOTION SAVR >=moderate SVD PREDICTED before freeing the porcine class, max |diff| (pp) |  | 3.06 | 5.0 | PASS | sim 5/10 y = 7.7/21.6, target 10.7/20.7 |
| onset | notion_fig3_modsvd_savr (porcine share freed), max |sim - target| (pp) |  | 2.88 | 3.0 | PASS | porcine scale -> 12.0 y, shape -> 1.21; sim 5/10 y = 10.2/23.5 |
| onset | Sapien XT VARC-3 SVD at 5 y (%), fitted | 1.6 | 1.67 | 0.8 | PASS | balloon scale x1.76 -> 52.7 y |
| check | Sapien 3 VARC-3 SVD at 5 y (%) | 0.7 | 2.34 | 2.0 | PASS | single-echo VARC-3 labels run 3% by 5 y even in PERIGON (A3), so ~1-2% is the noise floor of the definition |
| check | kermen_fig2b_freedom_stage3: max |sim - target| while n>=30 (pp) |  | 2.27 | 5.0 | PASS | sim at 10 y = 85.3, target 87.2 |
| check | notion_fig3_sevsvd_savr: max |sim - target| while n>=30 (pp) |  | 4.18 | 3.0 | FAIL | sim at 10 y = 14.1, target 9.9 |
| check | notion_fig3_sevsvd_tavi: max |sim - target| while n>=30 (pp) |  | 3.66 | 3.0 | FAIL | sim at 10 y = 5.1, target 1.4 |
| check | kermen_fig4a_explant_svd: max |sim - target| while n>=30 (pp) |  | 1.03 | 3.0 | PASS | sim at 10 y = 8.1, target 8.4 |
| check | NOTION SAVR reintervention by 10 y (%), age 79 | 2.2-4.3 | 4.86 |  | PASS |  |
| result | calibration wall time (s) |  | 256 |  | PASS | 32 processes, n_sim 20000 |
| result | onset perimount: Weibull shape / scale / shift |  | 5.60 / 10.49 y / +0.00 y |  | PASS |  |
| result | onset trifecta: Weibull shape / scale / shift |  | 0.92 / 31.64 y / +0.00 y |  | PASS |  |
| result | onset porcine: Weibull shape / scale / shift |  | 1.21 / 12.00 y / +0.00 y |  | PASS |  |
| result | onset new_savr: Weibull shape / scale / shift |  | 5.60 / 10.49 y / +0.00 y |  | PASS |  |
| result | onset balloon: Weibull shape / scale / shift |  | 2.00 / 52.67 y / +0.00 y |  | PASS |  |
| result | onset self_expanding: Weibull shape / scale / shift |  | 2.50 / 6.98 y / +0.00 y |  | PASS |  |
| result | onset new_thv: Weibull shape / scale / shift |  | 2.50 / 6.98 y / +0.00 y |  | PASS |  |
| result | noise sigma_flow / sigma_cw / sigma_lvot_vti / sigma_lvot_diam |  | 0.055 / 0.045 / 0.08 / 0.05 |  | PASS |  |
| result | death a (Gompertz level at age 75) |  | 0.0315 |  | PASS |  |
