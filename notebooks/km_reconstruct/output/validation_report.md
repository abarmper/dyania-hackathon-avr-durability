# Reconstruction validation report

112 checks, 0 FAIL

| curve | quantity | published | reconstructed | tol | result | note |
|---|---|---|---|---|---|---|
| kermen_fig2a_overall_survival | digitised value at 5 y | 80.9 | 81.12 | 1.5 | PASS |  |
| kermen_fig2a_overall_survival | digitised value at 10 y | 66.7 | 66.81 | 1.5 | PASS | window [62.4, 70.0] around t=10; SE=4.4 |
| kermen_fig2a_overall_survival | max |reconstructed - digitised| while n at risk >= 30 (t <= 9.8 y), pp |  | 0.302 | 1.0 | PASS |  |
| kermen_fig2a_overall_survival | max |reconstructed - digitised| in the tail (n < 30), pp, informational |  | 1.664 |  | PASS |  |
| kermen_fig2a_overall_survival | reconstructed (Guyot) value at 5 y | 80.9 | 81.19 | 1.5 | PASS |  |
| kermen_fig2a_overall_survival | reconstructed (Guyot) value at 10 y | 66.7 | 67.32 | 1.5 | PASS | window [64.0, 70.2] around t=10; SE=4.4 |
| kermen_fig2a_overall_survival | total events (soft: within max(2, 10%)) | 88 | 84 | 9 | PASS | exact total not reproduced; see algorithm warnings |
| kermen_fig2a_overall_survival | numbers at risk reproduced |  |  |  | PASS | [338, 247, 34, 20] vs [338, 247, 34, 20] |
| kermen_fig2a_overall_survival | IPD count equals n0 | 338 | 338 | 0 | PASS |  |
| kermen_fig2a_overall_survival | algorithm warnings |  |  |  | PASS | cannot reach total events 88; got 84 with no censoring left |
| kermen_fig2a_overall_survival | exact step inversion (diagnostic): total events; implied n at risk at first step after 5 y | 88 | 73 |  | PASS | implied n at risk after 5 y = 56; max |diff| vs digitised 27.3 pp |
| kermen_fig2a_valve_related_survival | digitised value at 5 y | 92.5 | 93.02 | 1.5 | PASS |  |
| kermen_fig2a_valve_related_survival | digitised value at 10 y | 86.0 | 92.5 | 1.5 | PASS | window [86.4, 92.5] around t=10; SE=6.1 |
| kermen_fig2a_valve_related_survival | max |reconstructed - digitised| while n at risk >= 30 (t <= 9.9 y), pp |  | 0.177 | 1.0 | PASS |  |
| kermen_fig2a_valve_related_survival | max |reconstructed - digitised| in the tail (n < 30), pp, informational |  | 1.61 |  | PASS |  |
| kermen_fig2a_valve_related_survival | reconstructed (Guyot) value at 5 y | 92.5 | 93.05 | 1.5 | PASS |  |
| kermen_fig2a_valve_related_survival | reconstructed (Guyot) value at 10 y | 86.0 | 92.63 | 1.5 | PASS | window [88.0, 92.6] around t=10; SE=6.1 |
| kermen_fig2a_valve_related_survival | total events (soft: within max(2, 10%)) | 25 | 23 | 2 | PASS | exact total not reproduced; see algorithm warnings |
| kermen_fig2a_valve_related_survival | numbers at risk reproduced |  |  |  | PASS | [338, 247, 34, 20] vs [338, 247, 34, 20] |
| kermen_fig2a_valve_related_survival | IPD count equals n0 | 338 | 338 | 0 | PASS |  |
| kermen_fig2a_valve_related_survival | algorithm warnings |  |  |  | PASS | cannot reach total events 25; got 23 with no censoring left |
| kermen_fig2a_valve_related_survival | exact step inversion (diagnostic): total events; implied n at risk at first step after 5 y | 25 | 23 |  | PASS | implied n at risk after 5 y = 168; max |diff| vs digitised 4.1 pp |
| kermen_fig2b_freedom_stage23 | digitised value at 5 y | 98.5 | 97.03 | 1.5 | PASS |  |
| kermen_fig2b_freedom_stage23 | digitised value at 10 y | 60.9 | 63.77 | 1.5 | PASS | window [56.5, 73.6] around t=10; SE=7.0 |
| kermen_fig2b_freedom_stage23 | digitised value at end of follow-up (end of follow-up) | 26.7 | 26.43 | 1.5 | PASS |  |
| kermen_fig2b_freedom_stage23 | max |reconstructed - digitised| while n at risk >= 30 (t <= 9.8 y), pp |  | 0.448 | 1.0 | PASS |  |
| kermen_fig2b_freedom_stage23 | max |reconstructed - digitised| in the tail (n < 30), pp, informational |  | 3.211 |  | PASS |  |
| kermen_fig2b_freedom_stage23 | reconstructed (Guyot) value at 5 y | 98.5 | 97.04 | 1.5 | PASS |  |
| kermen_fig2b_freedom_stage23 | reconstructed (Guyot) value at 10 y | 60.9 | 64.7 | 1.5 | PASS | window [58.0, 73.8] around t=10; SE=7.0 |
| kermen_fig2b_freedom_stage23 | reconstructed (Guyot) value at end of follow-up (end of follow-up) | 26.7 | 29.65 | 25.0 | PASS | tolerance = one event with n=4 at risk |
| kermen_fig2b_freedom_stage23 | total events (soft: within max(2, 10%)) | 49 | 49 | 5 | PASS |  |
| kermen_fig2b_freedom_stage23 | numbers at risk reproduced |  |  |  | PASS | [338, 247, 34, 20] vs [338, 247, 34, 20] |
| kermen_fig2b_freedom_stage23 | IPD count equals n0 | 338 | 338 | 0 | PASS |  |
| kermen_fig2b_freedom_stage23 | exact step inversion (diagnostic): total events; implied n at risk at first step after 5 y | 49 | 42 |  | PASS | implied n at risk after 5 y = 94; max |diff| vs digitised 9.9 pp |
| kermen_fig2b_freedom_stage3 | digitised value at 5 y | 99.6 | 100.0 | 1.5 | PASS |  |
| kermen_fig2b_freedom_stage3 | digitised value at 10 y | 88.3 | 87.21 | 1.5 | PASS | window [87.2, 91.5] around t=10; SE=5.0 |
| kermen_fig2b_freedom_stage3 | digitised value at end of follow-up (end of follow-up) | 58.9 | 57.92 | 1.5 | PASS |  |
| kermen_fig2b_freedom_stage3 | risk table |  |  |  | PASS | printed row NOT used (see config notes); Guyot no-numbers-at-risk case |
| kermen_fig2b_freedom_stage3 | published total events (informational; inconsistent with curve, see config notes) | 11 |  |  | PASS |  |
| kermen_fig2b_freedom_stage3 | censoring |  |  |  | PASS | 289 censoring times reused from kermen_fig2b_freedom_stage23 (same panel) |
| kermen_fig2b_freedom_stage3 | max |reconstructed - digitised| while n at risk >= 30 (t <= 12.0 y), pp |  | 0.488 | 1.0 | PASS |  |
| kermen_fig2b_freedom_stage3 | max |reconstructed - digitised| in the tail (n < 30), pp, informational |  | 0.0 |  | PASS |  |
| kermen_fig2b_freedom_stage3 | reconstructed (Guyot) value at 5 y | 99.6 | 100.0 | 1.5 | PASS |  |
| kermen_fig2b_freedom_stage3 | reconstructed (Guyot) value at 10 y | 88.3 | 87.7 | 1.5 | PASS | window [87.7, 91.4] around t=10; SE=5.0 |
| kermen_fig2b_freedom_stage3 | reconstructed (Guyot) value at end of follow-up (end of follow-up) | 58.9 | 57.63 | 2.86 | PASS | tolerance = one event with n=35 at risk |
| kermen_fig2b_freedom_stage3 | IPD count equals n0 | 338 | 338 | 0 | PASS |  |
| kermen_fig2b_freedom_stage3 | exact step inversion (diagnostic): total events; implied n at risk at first step after 5 y |  | 12 |  | PASS | implied n at risk after 5 y = 115; max |diff| vs digitised 3.3 pp |
| kermen_fig4a_nonvalve_death | digitised value at 5 y | 12.0 | 11.61 | 1.5 | PASS |  |
| kermen_fig4a_nonvalve_death | digitised value at 10 y | 25.4 | 25.92 | 1.5 | PASS | window [22.8, 25.9] around t=10; SE=4.1 |
| kermen_fig4a_nonvalve_death | max |reconstructed - digitised| while n at risk >= 30 (t <= 9.9 y), pp |  | 0.42 | 1.0 | PASS |  |
| kermen_fig4a_nonvalve_death | max |reconstructed - digitised| in the tail (n < 30), pp, informational |  | 1.674 |  | PASS |  |
| kermen_fig4a_nonvalve_death | reconstructed (Guyot) value at 5 y | 12.0 | 11.76 | 1.5 | PASS |  |
| kermen_fig4a_nonvalve_death | reconstructed (Guyot) value at 10 y | 25.4 | 26.49 | 1.5 | PASS | window [23.1, 26.5] around t=10; SE=4.1 |
| kermen_fig4a_nonvalve_death | total events (soft: within max(2, 10%)) | 63 | 57 | 6 | PASS | exact total not reproduced; see algorithm warnings |
| kermen_fig4a_nonvalve_death | numbers at risk reproduced |  |  |  | PASS | [338, 248, 34, 20] vs [338, 248, 34, 20] |
| kermen_fig4a_nonvalve_death | IPD count equals n0 | 338 | 338 | 0 | PASS |  |
| kermen_fig4a_nonvalve_death | algorithm warnings |  |  |  | PASS | cannot reach total events 63; got 57 with no censoring left |
| kermen_fig4a_nonvalve_death | exact step inversion (diagnostic): total events; implied n at risk at first step after 5 y | 63 | 59 |  | PASS | implied n at risk after 5 y = 111; max |diff| vs digitised 15.3 pp |
| kermen_fig4a_valve_related_death | digitised value at 5 y | 7.1 | 7.49 | 1.5 | PASS |  |
| kermen_fig4a_valve_related_death | digitised value at 10 y | 7.4 | 7.89 | 1.5 | PASS | window [7.9, 11.7] around t=10; SE=1.5 |
| kermen_fig4a_valve_related_death | max |reconstructed - digitised| while n at risk >= 30 (t <= 9.9 y), pp |  | 0.176 | 1.0 | PASS |  |
| kermen_fig4a_valve_related_death | max |reconstructed - digitised| in the tail (n < 30), pp, informational |  | 1.107 |  | PASS |  |
| kermen_fig4a_valve_related_death | reconstructed (Guyot) value at 5 y | 7.1 | 7.57 | 1.5 | PASS |  |
| kermen_fig4a_valve_related_death | reconstructed (Guyot) value at 10 y | 7.4 | 7.99 | 1.5 | PASS | window [8.0, 12.8] around t=10; SE=1.5 |
| kermen_fig4a_valve_related_death | total events (soft: within max(2, 10%)) | 25 | 25 | 2 | PASS |  |
| kermen_fig4a_valve_related_death | numbers at risk reproduced |  |  |  | PASS | [338, 248, 34, 20] vs [338, 248, 34, 20] |
| kermen_fig4a_valve_related_death | IPD count equals n0 | 338 | 338 | 0 | PASS |  |
| kermen_fig4a_valve_related_death | exact step inversion (diagnostic): total events; implied n at risk at first step after 5 y | 25 | 22 |  | PASS | implied n at risk after 5 y = 119; max |diff| vs digitised 6.9 pp |
| kermen_fig4a_explant_svd | digitised value at 5 y | 0.0 | 0.51 | 1.5 | PASS |  |
| kermen_fig4a_explant_svd | digitised value at 10 y | 8.0 | 8.39 | 1.5 | PASS | window [5.4, 8.4] around t=10; SE=3.4 |
| kermen_fig4a_explant_svd | max |reconstructed - digitised| while n at risk >= 30 (t <= 9.9 y), pp |  | 0.325 | 1.0 | PASS |  |
| kermen_fig4a_explant_svd | max |reconstructed - digitised| in the tail (n < 30), pp, informational |  | 1.213 |  | PASS |  |
| kermen_fig4a_explant_svd | reconstructed (Guyot) value at 5 y | 0.0 | 0.59 | 1.5 | PASS |  |
| kermen_fig4a_explant_svd | reconstructed (Guyot) value at 10 y | 8.0 | 9.61 | 1.5 | PASS | window [5.1, 9.6] around t=10; SE=3.4 |
| kermen_fig4a_explant_svd | total events (soft: within max(2, 10%)) | 9 | 8 | 2 | PASS | exact total not reproduced; see algorithm warnings |
| kermen_fig4a_explant_svd | numbers at risk reproduced |  |  |  | PASS | [338, 248, 34, 20] vs [338, 248, 34, 20] |
| kermen_fig4a_explant_svd | IPD count equals n0 | 338 | 338 | 0 | PASS |  |
| kermen_fig4a_explant_svd | algorithm warnings |  |  |  | PASS | cannot reach total events 9; got 8 with no censoring left |
| kermen_fig4a_explant_svd | exact step inversion (diagnostic): total events; implied n at risk at first step after 5 y | 9 | 11 |  | PASS | implied n at risk after 5 y = 286; max |diff| vs digitised 0.4 pp |
| notion_fig3_modsvd_tavi | digitised value at 10 y | 15.4 | 15.38 | 1.7 | PASS | window [14.4, 15.4] around t=10; SE=None |
| notion_fig3_modsvd_tavi | max |reconstructed - digitised| while n at risk >= 30 (t <= 10.0 y), pp |  | 0.997 | 1.5 | PASS |  |
| notion_fig3_modsvd_tavi | max |reconstructed - digitised| in the tail (n < 30), pp, informational |  | 0.0 |  | PASS |  |
| notion_fig3_modsvd_tavi | reconstructed (Guyot) value at 10 y | 15.4 | 14.96 | 1.7 | PASS | window [15.0, 15.0] around t=10; SE=None |
| notion_fig3_modsvd_tavi | numbers at risk reproduced |  |  |  | PASS | [134, 131, 128, 117, 109, 96, 82, 71, 56, 44, 30] vs [134, 131, 128, 117, 109, 96, 82, 71, 56, 44, 30] |
| notion_fig3_modsvd_tavi | IPD count equals n0 | 134 | 134 | 0 | PASS |  |
| notion_fig3_modsvd_savr | digitised value at 10 y | 20.8 | 20.66 | 1.7 | PASS | window [19.9, 20.7] around t=10; SE=None |
| notion_fig3_modsvd_savr | max |reconstructed - digitised| while n at risk >= 30 (t <= 10.0 y), pp |  | 1.233 | 1.5 | PASS |  |
| notion_fig3_modsvd_savr | max |reconstructed - digitised| in the tail (n < 30), pp, informational |  | 0.0 |  | PASS |  |
| notion_fig3_modsvd_savr | reconstructed (Guyot) value at 10 y | 20.8 | 21.89 | 1.7 | PASS | window [19.4, 21.9] around t=10; SE=None |
| notion_fig3_modsvd_savr | numbers at risk reproduced |  |  |  | PASS | [123, 122, 116, 107, 96, 84, 69, 61, 48, 41, 32] vs [123, 122, 116, 107, 96, 84, 69, 61, 48, 41, 32] |
| notion_fig3_modsvd_savr | IPD count equals n0 | 123 | 123 | 0 | PASS |  |
| notion_fig3_sevsvd_tavi | digitised value at 10 y | 1.5 | 1.42 | 1.7 | PASS | window [1.4, 1.4] around t=10; SE=None |
| notion_fig3_sevsvd_tavi | max |reconstructed - digitised| while n at risk >= 30 (t <= 10.0 y), pp |  | 0.64 | 1.5 | PASS |  |
| notion_fig3_sevsvd_tavi | max |reconstructed - digitised| in the tail (n < 30), pp, informational |  | 0.0 |  | PASS |  |
| notion_fig3_sevsvd_tavi | reconstructed (Guyot) value at 10 y | 1.5 | 0.78 | 1.7 | PASS | window [0.8, 0.8] around t=10; SE=None |
| notion_fig3_sevsvd_tavi | numbers at risk reproduced |  |  |  | PASS | [134, 132, 129, 118, 109, 96, 82, 73, 62, 51, 40] vs [134, 132, 129, 118, 109, 96, 82, 73, 62, 51, 40] |
| notion_fig3_sevsvd_tavi | IPD count equals n0 | 134 | 134 | 0 | PASS |  |
| notion_fig3_sevsvd_savr | digitised value at 10 y | 10.0 | 9.91 | 1.7 | PASS | window [9.9, 9.9] around t=10; SE=None |
| notion_fig3_sevsvd_savr | max |reconstructed - digitised| while n at risk >= 30 (t <= 10.0 y), pp |  | 0.769 | 1.5 | PASS |  |
| notion_fig3_sevsvd_savr | max |reconstructed - digitised| in the tail (n < 30), pp, informational |  | 0.0 |  | PASS |  |
| notion_fig3_sevsvd_savr | reconstructed (Guyot) value at 10 y | 10.0 | 10.68 | 1.7 | PASS | window [10.7, 10.7] around t=10; SE=None |
| notion_fig3_sevsvd_savr | numbers at risk reproduced |  |  |  | PASS | [123, 122, 119, 110, 100, 91, 79, 70, 58, 50, 39] vs [123, 122, 119, 110, 100, 91, 79, 70, 58, 50, 39] |
| notion_fig3_sevsvd_savr | IPD count equals n0 | 123 | 123 | 0 | PASS |  |
| wakami_fig1_svd | digitised value at 5 y | 4.8 | 4.84 | 1.7 | PASS |  |
| wakami_fig1_svd | digitised value at 7 y | 6.6 | 6.29 | 1.7 | PASS |  |
| wakami_fig1_svd | max |reconstructed - digitised| while n at risk >= 30 (t <= 6.6 y), pp |  | 0.925 | 1.5 | PASS |  |
| wakami_fig1_svd | max |reconstructed - digitised| in the tail (n < 30), pp, informational |  | 1.8 |  | PASS |  |
| wakami_fig1_svd | reconstructed (Guyot) value at 5 y | 4.8 | 5.37 | 1.7 | PASS |  |
| wakami_fig1_svd | reconstructed (Guyot) value at 7 y | 6.6 | 5.37 | 1.7 | PASS |  |
| wakami_fig1_svd | total events (soft: within max(2, 10%)) | 7 | 7 | 2 | PASS |  |
| wakami_fig1_svd | numbers at risk reproduced |  |  |  | PASS | [110, 100, 84, 58, 24] vs [110, 100, 84, 58, 24] |
| wakami_fig1_svd | IPD count equals n0 | 110 | 110 | 0 | PASS |  |
