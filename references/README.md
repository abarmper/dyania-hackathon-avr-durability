# References

Curated literature for the AVR durability challenge. Deliberately small and directed rather than exhaustive: 23 items, of which 21 are downloaded PDFs and 2 are full text where no open-access PDF exists.

Every entry below was verified against Europe PMC or the publisher landing page — no DOI here is reconstructed from memory. Items marked **paywalled** could not be downloaded from this machine; the DOI is given so a team member with institutional access can fetch them.

Naming: `A` = definitions and guidelines, `B` = durability evidence, `C` = risk factors, `D` = methods, `E` = data extraction.

---

## A. Definitions and guidelines (what "SVD" legally means for our label)

| File | Citation | Why it matters |
|---|---|---|
| `A1_capodanno2017_eapci_esc_eacts_svd_definitions.pdf` | Capodanno D, et al. Standardized definitions of structural deterioration and valve failure in assessing long-term durability of transcatheter and surgical aortic bioprosthetic valves. *Eur Heart J* 2017;38(45):3382–90. doi:10.1093/eurheartj/ehx303 | The original EAPCI/ESC/EACTS staging. Defines morphological vs. moderate vs. severe haemodynamic SVD. One of the three definitions we must test our label against. |
| (not downloaded; see below) | **Current guideline.** Praz F, Borger MA, et al. 2025 ESC/EACTS Guidelines for the management of valvular heart disease. *Eur Heart J* 2025;46(44):4635–4736. doi:10.1093/eurheartj/ehaf194. Also in *Eur J Cardiothorac Surg* 2025;67(8):ezaf276, doi:10.1093/ejcts/ezaf276. | **Supersedes `A2`.** §14.2, verbatim: TTE "should be performed in patients receiving a BHV within 3 months after valve implantation, again at 1 year, and annually thereafter, or sooner if new cardiovascular symptoms occur" (2021: within 30 days). Table 12 gives moderate/severe haemodynamic valve deterioration criteria identical to VARC-3, against the echo 1–3 months post-procedure. §14.4.4.1: HALT on CT in 10–30% of aortic BHVs, routine CT not indicated. Valve choice: mechanical preferred <60 years, biological >65 years (aortic). TAVI recommended ≥70 years with tricuspid AS if anatomy suitable; SAVR <70 at low surgical risk. Lifelong low-dose aspirin after TAVI without an OAC indication. Not redistributable (publisher copyright), so cite by DOI. |
| `A2_vahanian2021_esc_eacts_vhd_guidelines.pdf` | *Superseded by the 2025 edition above; kept for the HALT figures.* Vahanian A, et al. 2021 ESC/EACTS Guidelines for the management of valvular heart disease. *Eur Heart J* 2022;43(7):561–632. doi:10.1093/eurheartj/ehab395 | **Contains the exact standard-of-care our proposal would replace.** Verbatim: echocardiography "should be performed within 30 days after valve implantation (i.e. baseline), at 1 year, and annually thereafter." Also gives HALT prevalence of 12.4% on anticoagulation vs 32.4% on antiplatelets at 3 months. |
| `A3_velders2024_svd_definitions_lack_consistency.pdf` | Velders BJJ, et al. Current definitions of hemodynamic structural valve deterioration after bioprosthetic aortic valve replacement lack consistency. *JTCVS Open* 2024;19:68–90. doi:10.1016/j.xjon.2024.02.023 | **The most important methodological paper we found.** In 1,118 patients with one valve model and core-lab echoes, the three definitions flagged 4.6%, 2.9% and 3.0% of patients. After a first SVD classification, 59–65% had no SVD at the next visit. Conclusion: current definitions are "strong negative predictors but inconsistent positive discriminators." This means a single threshold crossing is a noisy label, and it is the basis for our persistence-confirmed label design. |
| `A4_velders2022_varc3_durability_bovine_pericardial.pdf` | **Correction of authorship:** Kermen S, Strella J, Aupart A, Espitalier F, Aupart M, Bernard A, Bourguignon T. Durability of a bovine pericardial aortic bioprosthesis based on VARC-3 echocardiographic criteria. *JTCVS Open* 2022;11:72–80. doi:10.1016/j.xjon.2022.05.008 (filename kept for stability). | Magna Ease, n = 338, mean follow-up 6.6 y, annual echo, VARC-3 staging. **The only published curve of freedom from VARC-3 stage 2/3 SVD** (98.5% at 5 y, 60.9% at 10 y; 49 events), plus competing-risk explant and death curves. Reconstructed in `notebooks/km_reconstruct/` and used as the primary calibration target of the simulator. |
| `A5_zoghbi2024_ase_prosthetic_valve_guidelines.pdf` | Zoghbi WA, et al. Guidelines for the evaluation of prosthetic valve function with cardiovascular imaging: a report from the American Society of Echocardiography. *J Am Soc Echocardiogr* 2024;37(1):2–63. doi:10.1016/j.echo.2023.10.004 (free PDF hosted by ASE at asecho.org) | Appendix Tables A1/A2/A4 give **normal EOA, mean gradient and DVI by valve model and labelled size** (Sapien family, CoreValve/Evolut R, Trifecta, Epic, Perimount, Mitroflow, Mosaic, Inspiris); Table 5 the stenosis thresholds, Table 6 the SVD criteria, Table 7 the PPM cut-offs. Anchors the simulator's physics module (`notebooks/simulator/valve_tables.py`). |

**Paywalled but canonical — cite, do not download:**
- VARC-3: Généreux P, et al. Valve Academic Research Consortium 3: updated endpoint definitions for aortic valve clinical research. *Eur Heart J* 2021;42(19):1825–57. doi:10.1093/eurheartj/ehaa799 (also *J Am Coll Cardiol* 2021;77(21):2717–46, doi:10.1016/j.jacc.2021.02.038). The VARC-3 thresholds are quoted verbatim inside `B1` if you need them without the original.
- Dvir D, et al. Standardized definition of structural valve degeneration for surgical and transcatheter bioprosthetic aortic valves. *Circulation* 2018;137(4):388–99. doi:10.1161/CIRCULATIONAHA.117.030729
- Capodanno D, et al. Standardized definitions for bioprosthetic valve dysfunction following aortic or mitral valve replacement: JACC state-of-the-art review. *J Am Coll Cardiol* 2022;80(5):545–61. doi:10.1016/j.jacc.2022.06.002
- Otto CM, et al. 2020 ACC/AHA guideline for the management of patients with valvular heart disease. *Circulation* 2021;143(5):e72–227. doi:10.1161/CIR.0000000000000923
- Salaun E, et al. Rate, timing, correlates, and outcomes of hemodynamic valve deterioration after bioprosthetic surgical aortic valve replacement. *Circulation* 2018;138(10):971–85. doi:10.1161/CIRCULATIONAHA.118.035150 (source of the "early vs late deterioration factors" quoted in `B4`; not downloaded)
- Hahn RT, et al. Comprehensive echocardiographic assessment of normal transcatheter valve function. *JACC Cardiovasc Imaging* 2019;12(1):25–34. doi:10.1016/j.jcmg.2018.04.010 (its tables are reproduced in `A5` Appendix Tables A1–A2)

### VARC-3 thresholds, quoted from `B1` (verified)

- **Moderate SVD:** mean transprosthetic gradient ≥20 mmHg **and** an increase ≥10 mmHg from the 3-month echo, **or** new ≥moderate intraprosthetic regurgitation.
- **Severe SVD:** mean gradient ≥30 mmHg **and** an increase ≥20 mmHg from the 3-month echo, **or** new severe intraprosthetic regurgitation.
- **Bioprosthetic valve failure (BVF):** valve-related death, **or** severe haemodynamic SVD, **or** reintervention following a diagnosis of bioprosthetic valve dysfunction.
- **Patient–prosthesis mismatch**, indexed effective orifice area, body mass index <30 kg/m²: moderate 0.65–0.85 cm²/m², severe ≤0.65 cm²/m². If body mass index ≥30: moderate 0.55–0.70, severe ≤0.55.

Note the baseline is the **3-month** echo, not the discharge echo. Our label construction must respect this.

---

## B. Durability evidence and event rates

| File | Citation | Key numbers |
|---|---|---|
| `B1_thyregod2024_notion_10year_tavi_vs_savr.pdf` | Thyregod HGH, et al. Transcatheter or surgical aortic valve implantation: 10-year outcomes of the NOTION trial. *Eur Heart J* 2024;45(13):1116–24. doi:10.1093/eurheartj/ehae043 | The longest randomised TAVI-vs-SAVR follow-up. At 10 years: ≥moderate SVD 15.4% TAVI vs 20.8% SAVR (HR 0.7, CI 0.4–1.3); **severe SVD 1.5% vs 10.0% (HR 0.2, CI 0.04–0.7)**; BVF ~10% vs 13.8% (HR 0.7, CI 0.4–1.5); reintervention 4.3% vs 2.2%; severe non-structural dysfunction (mostly mismatch) 20.5% vs 43.0%; endocarditis 7.2% vs 7.4%. Also quotes PARTNER 2 at 5 years: SVD 2.20% TAVR vs 4.38% SAVR. |
| `B2_kostyunin2020_degeneration_mechanisms_update.pdf` | Kostyunin AE, et al. Degeneration of bioprosthetic heart valves: update 2020. *J Am Heart Assoc* 2020;9(19):e018506. doi:10.1161/JAHA.120.018506 | Mechanism review: calcification, lipid infiltration, immune response to alpha-Gal, mechanical fatigue, thrombosis. Source for biologically motivated features. |
| `B3_salaun2019_detection_prediction_jacc_review.pdf` | **Correction:** this is Cartlidge TRG, Doris MK, Sellers SL, … Newby DE, Dweck MR. Detection and prediction of bioprosthetic aortic valve degeneration. *J Am Coll Cardiol* 2019;73(10):1107–19. doi:10.1016/j.jacc.2018.12.056 (filename kept for stability). | Prospective 18F-fluoride PET-CT study, not a review. In 71 patients with apparently normal bioprostheses, PET uptake predicted deterioration (annualised change in peak velocity 0.30 vs 0.01 m/s/yr) and was the only independent predictor on multivariable analysis; every patient who developed new dysfunction had uptake at baseline. Gives the "subclinical degeneration is detectable years before echo" argument and the emergency-vs-elective redo mortality figures (22.6% vs 1.4%). |
| `B4_jaha2025_tavr_durability_svd_review.pdf` | Transcatheter aortic valve durability: focus on structural valve deterioration. *J Am Heart Assoc* 2025;14:e041505. doi:10.1161/JAHA.125.041505 | Most recent synthesis of TAVR durability, including where the evidence still runs out beyond 8–10 years. |
| `B5_johnston2015_perimount_12569_implants.txt` | Johnston DR, et al. Long-term durability of bioprosthetic aortic valves: implications from 12,569 implants. *Ann Thorac Surg* 2015;99(4):1239–47. doi:10.1016/j.athoracsur.2014.10.070 | Largest single-centre durability series (Cleveland Clinic), and the same institution as our sample data. Full text only, no open-access PDF. |
| `B6_bvf_treatment_volume_projection_2025.pdf` | Predicting treatment of bioprosthetic aortic valve failure in the United States: a proposed model. *Struct Heart* 2025;9:100339. doi:10.1016/j.shj.2024.100339 | Projects valve-in-valve volume to roughly 42,000 procedures by 2035. Good for the "why this matters now" slide. |
| `B7_redo_savr_vs_viv_patient_profiles_2026.pdf` | Reintervention for failed aortic bioprostheses: distinct patient profiles for redo surgery and valve-in-valve TAVR. *J Clin Med* 2026;15(2):474. doi:10.3390/jcm15020474 | Shows that which reintervention a patient gets depends on their profile. Relevant to how our risk output maps to an action. |

**Paywalled but worth citing:** Bourguignon 2015 Perimount 20-year (doi:10.1016/j.athoracsur.2014.09.030); Rodriguez-Gabella 2017 JACC durability review (doi:10.1016/j.jacc.2017.07.715); Salaun 2018 *Circulation* haemodynamic valve deterioration in 1,387 patients, 31% degeneration mostly subclinical (doi:10.1161/CIRCULATIONAHA.118.035150); Pibarot 2020 PARTNER-2 SVD (doi:10.1016/j.jacc.2020.08.049); Blackman 2019 UK TAVI durability (doi:10.1016/j.jacc.2018.10.078); Sénage 2014 Mitroflow early SVD (doi:10.1161/CIRCULATIONAHA.114.010400); Yongue 2021 Trifecta durability in 2,298 valves (doi:10.1016/j.athoracsur.2020.07.040).

---

## C. Risk factors

| File | Citation | Contribution |
|---|---|---|
| `C1_ppm_trifecta_early_svd_2022.pdf` | Impact of postoperative patient-prosthesis mismatch as a risk factor for early structural valve deterioration after aortic valve replacement with Trifecta bioprosthesis. *J Cardiothorac Surg* 2022;17:172. doi:10.1186/s13019-022-01918-3 | Open-access evidence for the mismatch-to-SVD pathway in a valve model known for early failure. Directly supports indexed effective orifice area as a feature. |

**Paywalled, but these are the canonical risk-factor papers:** Flameng 2010 *Circulation*, mismatch predicts SVD (doi:10.1161/CIRCULATIONAHA.109.901272); Briand 2006 *Circulation*, metabolic syndrome and faster degeneration (doi:10.1161/CIRCULATIONAHA.105.000422); Makkar 2015 *NEJM* subclinical leaflet thrombosis (doi:10.1056/NEJMoa1509233); Chakravarty 2017 *Lancet* (doi:10.1016/S0140-6736(17)30757-2).

Consensus predictor set assembled from the above: younger age at implant, patient–prosthesis mismatch / small indexed orifice area, specific valve models, renal failure, diabetes and metabolic syndrome, dyslipidaemia, calcium-phosphate disturbance, elevated early post-implant gradient, and absence of early anticoagulation.

---

## D. Methods

| File | Citation | Role in our design |
|---|---|---|
| `D1_rizopoulos_personalized_screening_intervals_joint_models.pdf` | Rizopoulos D, et al. Personalized screening intervals for biomarkers using joint models for longitudinal and survival data. *Biostatistics* 2016;17(1):149–64. arXiv:1503.06448 | **The method that turns a risk score into a surveillance interval.** Schedules the next measurement by information gain from a joint longitudinal-survival model. This is the backbone of the "when should the next echo be" output. |
| `D2_dynamic_prediction_savr_mortality_2023.pdf` | Pollack J, et al. Dynamic prediction modeling of postoperative mortality among patients undergoing surgical aortic valve replacement in a statewide cohort over a 12-year period. *JTCVS Open* 2023;15:94–112. doi:10.1016/j.xjon.2023.07.011 | **Nuance:** "dynamic" here means updating model coefficients over calendar time (calibration regression, dynamic logistic state-space model) to counter drift, on 44,546 isolated SAVR patients, not prediction from serial biomarkers. Still useful: it shows static risk models decay in this population and that yearly re-calibration is feasible, which is a deployment argument for our surveillance model. |
| `D3_lee2018_deephit_competing_risks.pdf` | Lee C, Zame W, Yoon J, van der Schaar M. DeepHit: a deep learning approach to survival analysis with competing risks. *AAAI* 2018. | Competing-risk deep survival baseline. Its successor Dynamic-DeepHit (*IEEE TBME* 2020, doi:10.1109/TBME.2019.2909027, paywalled) handles longitudinal inputs and is the natural advanced model here. |
| `D4_ishwaran2008_random_survival_forests.pdf` | Ishwaran H, et al. Random survival forests. *Ann Appl Stat* 2008;2(3):841–60. arXiv:0811.1645 | Robust nonlinear survival baseline for small event counts. |
| `D5_austin_fine2017_finegray_reporting.pdf` | Austin PC, Fine JP. Practical recommendations for reporting Fine-Gray model analyses for competing risk data. *Stat Med* 2017;36(27):4391–400. doi:10.1002/sim.7501 | Death is a competing risk for SVD in an elderly cohort. This tells us how to report it correctly. |
| `D6_uno2011_c_statistic_censored.txt` | Uno H, et al. On the C-statistics for evaluating overall adequacy of risk prediction procedures with censored survival data. *Stat Med* 2011;30(10):1105–17. doi:10.1002/sim.4154 | Justifies the discrimination metric. Full text only. |
| `D7_deep_learning_survival_analysis_review.pdf` | Wiegrebe S, et al. Deep learning for survival analysis: a review. arXiv:2305.14961 | Survey to justify model choice against simpler baselines. |
| `D8_collins2024_tripod_ai_reporting.pdf` | Collins GS, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning. *BMJ* 2024;385:e078378. doi:10.1136/bmj-2023-078378 | The reporting checklist our protocol should explicitly follow. Cheap way to score on study-design quality. |
| `D9_guyot2012_reconstruct_ipd_from_km_curves.pdf` | Guyot P, et al. Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan-Meier survival curves. *BMC Med Res Methodol* 2012;12:9. doi:10.1186/1471-2288-12-9 | **Our route to a defensible prototype cohort with no data access.** Reconstruct patient-level survival data from published durability curves by valve model, then calibrate a synthetic cohort to it. |
| `D10_evar_individualized_surveillance_2024.pdf` | Individualizing surveillance after endovascular aortic repair using a modular imaging algorithm. *Diagnostics* 2024;14(9):930. doi:10.3390/diagnostics14090930 | Precedent from an adjacent field where risk-stratified imaging surveillance replaced fixed intervals. Useful analogy when the panel asks whether this is realistic. |

**Paywalled:** Riley RD, et al. Calculating the sample size required for developing a clinical prediction model. *BMJ* 2020;368:m441. doi:10.1136/bmj.m441 — needed for the sample-size section; the formulae are widely reproduced.

---

## E. Getting structured data out of free text

| File | Citation | Contribution |
|---|---|---|
| `E1_transformer_structuring_cath_echo_reports_2026.pdf` | Converting unstructured cardiac catheterization and echocardiography reports into structured data using transformer-based language models. *JAMIA Open* 2026. doi:10.1093/jamiaopen/ooag036 | Open-access evidence that the extraction step in our pipeline is solved and validated. |

**Paywalled but directly on point:** A comparative analysis of privacy-preserving large language models for automated echocardiography report analysis. *JAMIA* 2025;32(7):1120. doi:10.1093/jamia/ocaf056 — reports **>96% accuracy for prosthetic valve detection** from echo reports (100% and 99.9% for the two best open-weight models with chain-of-thought prompting), while severity grading was much weaker at 54–86%. This is the single best citation for "we can reliably identify prosthetic valves in free text, locally, without sending data to a vendor." Also: Automated structured data extraction from intraoperative echocardiography reports using large language models. *Br J Anaesth* 2025. doi:10.1016/j.bja.2025.01.028

---

## How these PDFs were obtained

Publisher sites, PubMed Central's web interface, Europe PMC and figshare all block scripted downloads from this machine. Everything here came from the NIH PMC Open Access mirror on AWS (`pmc-oa-opendata` S3 bucket), arXiv, AAAI, or institutional repositories. The helper script that resolves a DOI or PMC ID to a PDF through that mirror lives in the session scratchpad; re-create it with `fetch_ref.py` logic if more papers are needed.

Two items are `.txt` rather than `.pdf` because the article is in PMC's open-access subset without a PDF rendition.

## Licensing

Every file here is redistributable. The PMC-sourced items carry Creative Commons licences, verified from the PMC metadata: nine are CC BY, one CC BY-NC, and five CC BY-NC-ND. The no-derivatives clause restricts modification, not redistribution, so committing them to a public fork is permitted with attribution — which the citations above provide. The remaining items come from arXiv, the AAAI proceedings, or institutional repositories (Liège), all of which permit redistribution.

The folder is about 33 MB. If the upstream maintainers object to binaries in the pull request, delete the PDFs and keep this README — every entry has a DOI and the retrieval method is documented above, so the set can be rebuilt in a few minutes.
