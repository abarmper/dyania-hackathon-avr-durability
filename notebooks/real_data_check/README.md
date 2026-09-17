# Real-data check on the provided notes

A small face-validity check of the label rule and both models on the de-identified notes in `data/notes_deidentified.xlsx`. It is **not a validation**: the notes have no serial prosthetic echoes, year-only dates and redacted ages, and the sample is enriched for failures.

```bash
/data/abar/alexenv/bin/python real_data_check.py     # uses the trained models in notebooks/analysis/output/models/
```

**Inputs.** Deleted notes are dropped. `curated_cohort.csv` holds the hand adjudication, one row per patient; the verbatim note excerpt behind each row is kept locally in `curated_cohort_with_quotes.csv`, which, like the provided `data/*.xlsx`, is not committed. It records the index valve, implant year, the follow-up echo of the *index* valve (the echo before any valve-in-valve procedure) and documented structural dysfunction. Baseline covariates come from regex with negation and family-history filtering (`real_data_check.py`).

**Cohorts.**
- **Primary:** 26 patients with a known index valve, implant year and documented valve status at least a year after implant; 11 have documented structural dysfunction (8 reinterventions or planned redo, 3 moderate or worse dysfunction).
- **Label only:** 3 more patients with an echo but a redacted implant year.
- **Uncertain:** 3 patients, used in a sensitivity analysis. Excluded as non-structural: endocarditis (2 patients) and paravalvular leak (1).

## Results

**1. What a notes-based pipeline can recover (all 117 patients).**

| Item | Patients |
|---|---|
| Valve model named | 108 (92%) |
| Labelled size | 113 (97%) |
| BMI | 55 (47%) |
| Numeric age | 0 (redacted everywhere) |
| A prosthetic-valve echo with a mean gradient | 42 (36%) |
| Prosthetic mean gradients in two or more different years | 2 (2%) |

**2. VARC-3 single-echo rule with the model × size reference (the missing-reference arm) vs documented dysfunction, 25 patients.**

| Rule | True pos. | False neg. | False pos. | True neg. | Sensitivity | Specificity |
|---|---|---|---|---|---|---|
| Haemodynamic VARC-3 (gradient or new regurgitation) | 6 | 1 | 2 | 16 | 86% | 89% |
| Full VARC-3 (adds DVI fall when reported) | 6 | 1 | 2 | 16 | 86% | 89% |

- **False negative:** Patient_042, moderate prosthetic stenosis documented at a mean gradient of 11 mmHg with LVEF 40%. A low-flow state hides the gradient rise.
- **False positives:** Patient_086 and Patient_106, stable gradients of 21 mmHg in 25 mm and 23 mm valves. The ASE table under-estimates their individual baseline. Patient_106's own prior echo (19 mmHg) as reference makes the change 2 mmHg and the label negative. This is the error the simulation predicted for the missing-reference arm, and it is why the output flags a table-based reference.

**3. Models.**

| Check | n (events) | Result |
|---|---|---|
| Landmark model at the follow-up echo: 5-y risk vs documented dysfunction | 20 (5) | AUC 0.97 (0.88–1.00) |
| Same, surgical valves ≥5 y after implant only | 9 (5) | AUC 0.90; the single-echo rule alone gives 0.78 |
| Implant-time boosted model: 10-y risk vs time to documented dysfunction | 26 (11) | Harrell's C 0.59 (0.36–0.82) |
| Cause-specific Cox, full baseline | 26 (11) | C 0.60 (0.38–0.84) |
| Published risk-factor Cox (B4 set) | 26 (11) | C 0.51 (0.19–0.86) |
| Implant-time boosted model, including uncertain cases | 29 (14) | C 0.64 (0.42–0.82) |

## How to read this

- The landmark AUC is **partly circular**: documented dysfunction is often diagnosed from the same echo the model scores. It shows the model ranks real failing valves above well-functioning ones, and slightly better than the binary rule. It does not show prediction ahead of time.
- The implant-time check is **uninformative**: every confidence interval crosses 0.5. Age, the strongest predictor, is redacted in all notes and was imputed; times are year-level; the sample is small and failure-enriched; and all 11 failures are older surgical valves.
- What it does show: the pipeline runs end to end on real clinical text; valve model and size are recoverable; the label rule with a table reference behaves as simulated, including its false positives. A real validation needs structured serial echoes, exact dates and ages, which only 2 of 117 patients approach here. That is the data the protocol requests.
- Comorbidities are taken from all notes, including later history sections, because operative reports rarely list them. This is a mild look-ahead, acceptable for chronic conditions in a face-validity check.
