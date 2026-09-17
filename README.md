# Dyania Health Hackathon 2026 — Team Submission Repository

**Challenge:** Build a study — using machine learning, a statistical model, or whatever approach you prefer — proposing a protocol to predict aortic valve durability in patients with a bioprosthetic aortic valve replacement.
**Event:** September 15–17, 2026 (3 days)
**Team size:** 2–3 ML engineers
**Team:** CardioNTUA
**Members:** Athanasia Karagiannopoulou — ML Engineer, Anna Panagiotakopoulou — ML Engineer, Evangelie Sintou — Clinician, Alexandros Barmperis — ML Engineer

---

## Problem Statement

Bioprosthetic valves deteriorate, and patients are followed with the same annual echo whatever their risk (ESC/EACTS 2025). Single-echo SVD labels are unreliable (59–65% vanish at the next visit, Velders 2024), so deterioration is confirmed late: in our calibrated simulation, annual surveillance confirms moderate deterioration a mean 1.5 years (P90 2.9 years) after it begins. Late recognition turns elective reinterventions into urgent ones (redo mortality 22.6% emergency vs 1.4% elective), while low-risk patients receive echoes they do not need, and valve-in-valve volume is rising fast.

## Our Approach

A study protocol for **dynamic, competing-risks prediction of confirmed VARC-3 moderate SVD at every follow-up echo**, turned into a **risk-adapted echo schedule**. Target data: echo reporting databases, STS ACSD and STS/ACC TVT registries, EHR and claims linkage. The prototype is built on a synthetic cohort of 10,000 valves that we generated and calibrated to curves reconstructed from three published cohorts (Kermen 2022, NOTION 10-year, Wakami 2022), physics-anchored to ASE 2024 valve tables, reproducing label flicker. The model is a discrete-time cause-specific hazard model (gradient boosting over 6-month periods) with SHAP explanations. A surveillance simulation with feedback compares schedules on the same patients.

| Headline (synthetic test set) | Value |
|---|---|
| 5-y AUC / C / calibration slope on echoes not already abnormal | 0.875 (0.851–0.899) / 0.855 / 1.04 |
| Time-since-implant comparator | AUC 0.651 |
| Implant-time model vs published risk-factor Cox, AUC at 10 y | 0.823 vs 0.771 |
| Risk-adapted schedule (1% rule) vs guideline plus confirmation echo | 44% fewer echoes; +0.7 months mean, +1.8 months P90 detection delay |
| Confirmation echo added to annual surveillance | mean delay 1.5 → 1.1 y |
| Real notes, 26 patients / 11 documented failures (face validity only) | label rule sensitivity 86%, specificity 89%; landmark model AUC 0.97 at the follow-up echo (partly circular); implant-time C 0.59, uninformative |

## Key Design Decisions

| Decision | Rationale |
|---|---|
| Primary endpoint = VARC-3 stage ≥2 **confirmed on two consecutive echoes**; dynamic models date the event at the confirming echo | single-echo labels flicker; the label only exists when confirmed, so this timing avoids leaking future echoes into predictions |
| Landmark discrete-time **cause-specific hazard** model (boosted multinomial) | handles censoring and competing death (40% die in follow-up) and yields risk over any interval, which the scheduling rule needs |
| **Simulated cohort calibrated to reconstructed published curves**, with a noise-free "oracle" | no public dataset has serial prosthetic-valve echoes with valve model; the oracle lets us audit the label definitions themselves (confirmed VARC-3: 60% sensitivity, 1.5-y median lag) |
| Features only if quantified in the literature and present on TVT/STS forms; physics features from ASE 2024; registry-only comparator | avoids invented structure; showed the engineered features add only 0.003 AUC, so the gain comes from the framing |
| Every surveillance policy includes a **confirmation echo** at 3–6 months and is compared on identical simulated patients | separates what the model buys from what the confirmation rule buys |

---

## Reproduce

```bash
# Python env with numpy, pandas, scipy, scikit-learn, lifelines, statsmodels, shap, matplotlib, pyarrow
cd notebooks/km_reconstruct && python run_reconstruction.py --figures kermen_fig2 kermen_fig4a notion_fig3 wakami_fig1 --strict
cd ../simulator && python cli.py calibrate && python cli.py validate && python cli.py generate --n 10000 --seed 1
cd ../analysis && python run_all.py all && python run_all.py explain && python run_all.py policy \
  && python run_all.py oracle && python run_all.py secondary && python summarise.py
```

Results: `notebooks/analysis/output/results_summary.md`; notebooks `notebooks/01_*`, `02_*`, `03_*`.

## Repository Structure

```
.
├── README.md                     # this file
├── protocol/study_protocol.md    # full study design (main deliverable)
├── model/approach.md             # ML methodology, validation, outputs, limitations
├── data/data_plan.md             # data sources, preprocessing, availability; data/synthetic/ = generated cohort
├── presentation/                 # slides.pptx + slides.pdf, built by build_deck.py with the charts and schematics in figures/
├── report/main.pdf               # 13-page write-up of the whole submission (LaTeX source alongside)
├── docs/                         # problem analysis, literature, ideas, cohort spec, modelling design
├── references/                   # curated open-access papers with README
└── notebooks/
    ├── km_reconstruct/           # individual-patient data reconstructed from published curves
    ├── simulator/                # calibrated synthetic cohort generator
    ├── analysis/                 # models, evaluation, SHAP, oracle analyses, surveillance policies
    └── real_data_check/          # face-validity check of the label rule and models on the provided notes
```

## Submission Checklist

- [x] `README.md` — team overview, problem framing, key design decisions
- [x] `protocol/study_protocol.md` — complete study protocol
- [x] `model/approach.md` — ML methodology
- [x] `data/data_plan.md` — data plan
- [x] `presentation/slides.pdf` — slide deck (also `slides.pptx`)
- [x] `notebooks/` — proof-of-concept

> ⚠️ No real patient data is used for modelling; `data/synthetic/` is simulated. Submission: branch `team/CardioNTUA`, pull request to `dyaniahealth/dyania-hackathon-avr-durability`, left open.
