# Reconstructing the published SVD curves (design choice D6)

Automatic reconstruction of individual-patient time-to-event data from the published Kaplan–Meier /
cumulative-incidence figures of three cohorts, with every published number used as a validation check.
The result (`output/targets.json`) is the calibration target for the synthetic cohort's structural
valve deterioration (SVD) onset module; see `docs/04_synthetic_cohort_spec.md` §4 (D6) and §5.

Run with the project environment:

```bash
/data/abar/alexenv/bin/python run_reconstruction.py --figures kermen_fig2 kermen_fig4a notion_fig3 wakami_fig1 --strict
/data/abar/alexenv/bin/python -m pytest tests -q
/data/abar/alexenv/bin/python calibrate_onset.py
```

## Sources and what was automated

| Source (references/) | Figure | Format | Curves | Manual input | Guyot case |
|---|---|---|---|---|---|
| A4 Kermen 2022, Magna Ease, n = 338 | Fig 2A/2B (KM), Fig 4A (competing-risk CIF) | vector PDF paths | overall survival, valve-related survival, freedom from VARC-3 stage 2/3 SVD (**primary label**), freedom from stage 3, explant for SVD, valve-related death, non-valve death | none (risk row and legend parsed from the PDF text) | all information |
| B1 NOTION 10-year, TAVI n = 134 / SAVR n = 123 at risk at baseline | Fig 3 (≥ moderate and severe SVD, CIF) | RGB JPEG inside the PDF | 4 curves (2 panels × 2 arms) | the two "Patients at risk" rows per panel, typed from `output/qa_notion_fig3_risk_table.png`; y-tick step | risk table, no total |
| C1 Wakami 2022, Trifecta, n = 110 | Fig 1 (SVD CIF, death competing) | 8-bit grayscale image | 1 curve | "Number at risk" row (110/100/84/58/24 at 0/20/40/60/80 months) and 7 events, from the figure/text | all information |

Pipeline: `pdf_vector.py` (content stream → segments, polygons, text with the current transformation
matrix) → `vector_extract.py` (axis calibration from tick marks, step-curve classification, legend
naming, censor "+" marks, risk-row parsing) or `raster_trace.py` (axis and tick detection from pixel
projections, colour / darkness masks, connected components or lowest-pixel tracing) → `guyot.py`
(Guyot et al. 2012 steps 1–8, four information cases, plus digitised censor marks for the last interval)
→ `validate.py` (every published anchor, risk row, event total, IPD count, reconstructed-vs-digitised
curve) → `run_reconstruction.py` (outputs below).

## Validation summary (`output/validation_report.md`, 112 checks, 0 FAIL)

| Curve | Published | Digitised | Guyot reconstruction | Risk row / events |
|---|---|---|---|---|
| Kermen overall survival 5 / 10 y | 80.9 / 66.7 | 81.1 / 66.8 | 81.1 / 66.8 | 338/247/34/20 exact; 84 events (88 published; 4 fall after the last risk time) |
| Kermen valve-related survival 5 / 10 y | 92.5 / 86.0 | 92.5 / 86.0 | 92.5 / 86.0 | exact; 23 (25) |
| Kermen freedom from stage 2/3 SVD 5 / 10 / end | 98.5 / 60.9 / 26.7 | 98.5 / 60.9 / 26.4 | 98.5 / 60.9 / 29.7 | exact; **49 (49)** |
| Kermen freedom from stage 3 SVD 5 / 10 / end | 99.6 / 88.3 / 58.9 | 100.0 / 87.2 / 57.9 | 100.0 / 87.2 / 57.6 | see note 1 |
| Kermen explant for SVD CIF 5 / 10 y | 0.0 / 8.0 | 0.0 / 8.4 | 0.0 / 8.4 | exact; 8 (9) |
| Kermen valve-related death CIF 5 / 10 y | 7.1 / 7.4 | 7.1 / 7.4 | 7.1 / 7.4 | exact; 25 (25) |
| Kermen non-valve death CIF 5 / 10 y | 12.0 / 25.4 | 12.0 / 25.4 | 12.0 / 25.4 | exact; 57 (63) |
| NOTION ≥ moderate SVD at 10 y, TAVI / SAVR | 15.4 / 20.8 | 15.4 / 20.7 | 15.0 / 21.9 | both rows exact |
| NOTION severe SVD at 10 y, TAVI / SAVR | 1.5 / 10.0 | 1.4 / 9.9 | 0.8 / 10.7 | both rows exact |
| Wakami SVD CIF 5 / 7 y | 4.8 / 6.6 | 4.8 / 6.3 | 5.4 / 5.4 | exact; **7 (7)** |

Tolerances: 1.5 percentage points at well-populated times (+0.2 for rasters); at ≥ 9 y the published
value must fall inside the curve's range over ±0.25 y widened by 1.5 points, or within the published
standard error; reconstructed vs digitised curve ≤ 1.0 (vector) / 1.5 (raster) points on a 0.1-y grid
while ≥ 30 patients remain, informational beyond; end-of-follow-up values within one event of the
reconstructed risk set; event totals within max(2, 10%).

Note 1 (documented inconsistency in the paper): the risk row printed under Kermen Fig 2 is the all-cause
survival risk set. With 247 at risk at 5 y the first stage-3 drop (0.9 points) would be two to three
events and the curve would imply at least 14 events, but the text reports 11; an exact step inversion
of the vector coordinates implies about 115 echo-followed patients at risk at 5 y. The stage-3 curve is
therefore reconstructed with the stage-2/3 curve's censoring times (same panel, same population) and the
published total is reported as informational. The digitised curve, which is what the simulator is
calibrated to, is unaffected.

## Outputs

| File | Content |
|---|---|
| `output/targets.json` | per curve: digitised step function on a 0.1-y grid, estimator (KM with death censored or Aalen–Johansen), n0, risk table, anchors, censor-mark times, provenance (page, calibration) |
| `output/ipd_<curve>.csv` | reconstructed pseudo individual-patient data (time, event) — the Guyot output |
| `output/ipd_exact_<curve>.csv` | vector-only exact step inversion, diagnostic only (single-event drops below 0.3 points are at the coordinate precision limit) |
| `output/onset_priors.json` | Weibull and log-logistic fits to the reconstructed IPD: **initialisers for calibration, not the latent onset distribution** (they describe observed, visit-timed labels) |
| `output/validation_report.md` | the 112 checks |
| `output/qa_<figure>.png` | digitised curve, Guyot reconstruction, published anchors with SE, censor marks |
| `output/calibrated_onset.json` | `calibrate_onset.py` result: for each SVD target, the residual of the direct parametric variant (latent onset = fit to the IPD) and of the calibrated variant (location shift and scale multiplier), with the stand-in simulator |

## What the simulator must do with this

Do not sample latent onset times from `onset_priors.json`. The curves are echo-detected labels: using
them as the latent distribution would make the simulator's own observed curve later than the paper's by
the onset-to-threshold lag. `calibrate_onset.py` shows the hook: simulate a cohort matched on n, stratum,
follow-up and censoring, compute the same estimator on the simulated observed labels, and adjust the
latent distribution until the simulated observed curve matches `targets.json`. Replace
`default_simulate` with the real simulator's observed-label generator (same signature).

Stratum mapping: Kermen → SAVR Perimount-class, age ~70; NOTION SAVR → SAVR mixed, age 79;
NOTION TAVI → self-expanding TAVR, age 79; Wakami → Trifecta-class. Cumulative-incidence targets are
fed to Guyot as 1 − CIF ("AJ-as-KM"); their pseudo-IPD conflate death with censoring and are only used
for curve-to-curve comparison.

No PDF content is copied into the repository; the CSVs and JSON are derived data.
