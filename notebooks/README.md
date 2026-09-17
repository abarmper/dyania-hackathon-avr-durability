# Notebooks — Proof-of-Concept Implementation

## What is here (team addition)

- `km_reconstruct/` — reconstructs individual-patient data from the published SVD curves (Kermen 2022, NOTION 10-year, Wakami 2022) with the Guyot algorithm; 112 validation checks; produces the simulator's calibration targets. See its README.
- `simulator/` — the literature-calibrated synthetic cohort generator (physics-based latent valve process, four-source echo noise, competing events, three SVD definitions, oracle); `cli.py calibrate | validate | generate`. See its README.
- `01_generate_and_validate.ipynb` — thin notebook that loads the calibrated parameters, generates the cohort, and displays the calibration and validation reports with figures.
- `analysis/` — Idea 1 models (landmark discrete-time cause-specific hazard, Cox, comparators), evaluation with valve-level bootstrap, SHAP, oracle-only analyses, and the Idea 2 surveillance-policy simulation; `run_all.py all|explain|policy|oracle|secondary`, `summarise.py`. See its README.
- `02_features_and_models.ipynb`, `03_surveillance_policy.ipynb` — thin notebooks displaying the analysis outputs.
- Data: `../data/synthetic/` (CSV + parquet, `data_dictionary.md`).

Environment: `/data/abar/alexenv/bin/python` (numpy, pandas, scipy, matplotlib, lifelines, pyarrow, pytest); `simulator/requirements.txt`.

> This folder is **optional** but evaluated positively if present.

Place your Jupyter notebooks here. A strong submission includes:

- Data loading and exploratory analysis
- Feature engineering pipeline
- Model training and evaluation
- SHAP / feature importance visualisation

## Suggested Notebook Structure

```
notebooks/
├── 01_eda.ipynb                 # Exploratory data analysis on your chosen dataset
├── 02_feature_engineering.ipynb # Feature construction (gradient progression rate, EOA index, PPM flag, etc.)
├── 03_model_training.ipynb      # Model training, cross-validation, hyperparameter tuning
└── 04_evaluation.ipynb          # Metrics, calibration, SHAP plots, subgroup analysis
```

You can combine these into a single notebook if preferred — the split is just for readability.

## Environment

Document your dependencies here so the panel can reproduce your results:

```
python >= 3.10
pandas
numpy
scikit-learn
xgboost       # or your chosen framework
lifelines     # or scikit-survival, for time-to-event modelling
shap
matplotlib
jupyter
```

Or include a `requirements.txt` / `environment.yml` in this folder.
