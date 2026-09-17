# Presentation Slide Content

Outline of the 10 slides. The deck is built by `build_deck.py` (one slide spec, rendered to `slides.pptx` and to `slides.pdf`; charts and schematics in `figures/`).

---

## Slide 1 — Title

- **Knowing when to look again: dynamic SVD risk to personalise echo surveillance after bioprosthetic AVR**
- Team CardioNTUA — Athanasia Karagiannopoulou (ML Engineer), Anna Panagiotakopoulou (ML Engineer), Evangelie Sintou (Clinician), Alexandros Barmperis (ML Engineer)
- Dyania Health Hackathon, 17 September 2026

---

## Slide 2 — The Clinical Problem

- Bioprosthetic valves deteriorate (structural valve deterioration, SVD); severe SVD at 10 y: 10.0% SAVR vs 1.5% TAVI (NOTION), with the SAVR vs TAVR durability question still open.
- The label is noisy: 59–65% of first SVD classifications are gone at the next echo (Velders 2024).
- **The gap:** under annual surveillance, moderate deterioration is confirmed a mean **1.5 y (P90 2.9 y)** after it begins (our calibrated simulation).

## Slide 3 — Current Workflow Failure

- Same echo every year for everyone (ESC/EACTS 2025), regardless of valve, age or trajectory.
- Bottlenecks: missed visits (35% by year 4), echo variability, no confirmation rule, reference echo often absent.
- Cost: urgent redo surgery mortality 22.6% vs 1.4% elective; unnecessary echoes for low-risk patients; valve-in-valve volume rising.

## Slide 4 — Hypothesis

- A risk model updated at every echo can lengthen intervals for low-risk valves and shorten them for rising trajectories.
- Target: ≥25% fewer echoes with ≤3 months extra detection delay, plus earlier confirmation for high-risk valves.

---

## Slide 5 — Study Overview

- **Population:** adults with a first bioprosthetic SAVR or TAVR, reference echo or model × size norm, ≥1 follow-up echo; exclude mechanical, Ross, multi-valve, index ViV.
- **Primary endpoint:** VARC-3 stage ≥2 confirmed on two consecutive echoes, competing death/endocarditis; incidence at 5/8/10 y.
- **Secondary:** bioprosthetic valve failure, reintervention (elective vs urgent), gradient progression, surveillance performance.
- **Ground truth:** pathology > SVD reintervention > confirmed VARC-3 > heart-team adjudication; thrombosis, PVL, PPM, endocarditis adjudicated out.

## Slide 6 — Data Sources

- **Deployment:** echo reporting database, STS ACSD, STS/ACC TVT Registry, EHR with NLP, CMS/NDI linkage; target ≥5,000 development and ≥2,000 external valves.
- **Prototype:** 10,000-valve synthetic cohort calibrated to curves reconstructed from Kermen 2022, NOTION 10-y and Wakami 2022 (112 checks pass), physics from ASE 2024 tables; reproduces label flicker.
- Figure: `curves_sim_vs_target.png` (simulated vs published curves).

## Slide 7 — Validation Plan

- Valve-level split and bootstrap; the event dated at the confirming echo (no leakage); temporal and external validation in the real study.
- Metrics: cause-specific time-dependent AUC, C-index, calibration slope, plus decision metrics (echoes, delay).
- Comparators: time since implant, published risk-factor Cox, registry-only features; guideline, confirmation and de-escalation schedules.

---

## Slide 8 — Model

- Discrete-time cause-specific hazard model: gradient boosting over 6-month periods after each echo → risk over any interval.
- Features: change from reference, smoothed log gradient, slope, current and prior abnormal echoes, EOA/DVI change, expected EOA for model × size, PPM.
- Output: 5/8/10-y risk, tier, top-3 SHAP features, next-echo interval (longest of 6–36 months with risk ≤1%), unconfirmed-echo flag.

## Slide 9 — Results (synthetic test set)

- Echoes not already abnormal: **AUC 0.875 (0.851–0.899), C 0.855, calibration slope 1.04**; time since implant 0.651.
- Implant-time model 0.865/0.850/0.823 at 5/8/10 y vs published risk-factor Cox 0.756/0.756/0.771.
- **Surveillance:** 44% fewer echoes than guideline-plus-confirmation, +0.7 months mean delay; the confirmation echo alone cuts mean delay 1.5 → 1.1 y.
- Oracle audit: confirmed VARC-3 detects 60% of true moderate deterioration, median lag 1.5 y.
- Real notes (26 patients, 11 failures): label rule 86% sensitive / 89% specific, false positives from the missing reference echo; model ranks failing valves first (AUC 0.97, same-echo). Face validity only.
- Figures: `policy_frontier.png`, `shap_summary.png`, `calibration_5y.png`.

---

## Slide 10 — Impact and Next Steps

- Enables: risk-adapted echo scheduling with a confirmation rule, and an audit of SVD definitions against truth that no registry can do.
- Validation needs: linked echo database + STS/TVT at ≥3 centres; external registry.
- Next 3 months: run the pipeline on one institution's echo database (same schema), silent-mode prospective evaluation.
- Readiness: code and schema are registry-shaped; real-world performance is unproven until external validation.
