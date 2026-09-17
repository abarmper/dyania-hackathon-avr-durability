"""Small real-data check on the provided de-identified notes (data/notes_deidentified.xlsx).

What the notes allow: an index procedure note and, for some patients, one later note with a prosthetic-valve echo
summary. There are no serial prosthetic echoes, dates are year-only and ages are redacted, and the sample is enriched
for failures. So this is a FACE-VALIDITY check, not a validation:

  1. extraction yield on all 117 patients (what a registry-style pipeline can recover from notes);
  2. label-rule check: the VARC-3 single-echo rule with the model x size reference (our missing-reference arm)
     against clinician-documented structural valve dysfunction;
  3. model check: the implant-time model and the published risk-factor Cox on the hand-adjudicated cohort
     (Harrell's C against time to documented dysfunction), and the landmark model scored at the follow-up echo.

Hand adjudication with source quotes: curated_cohort.csv. Baseline covariates: regex on the notes (negation-aware).
Run: /data/abar/alexenv/bin/python real_data_check.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT / "notebooks" / "analysis"))
import analysis_config as C  # noqa: E402  (puts the simulator on sys.path too)
from data import Hist, baseline_frame  # noqa: E402
from features import history_features, implant_instances  # noqa: E402
from lifelines.utils import concordance_index  # noqa: E402
from physics import ppm_class  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from valve_tables import A_EFF, EOA_TABLE, MODEL_CLASS  # noqa: E402

OUT = HERE / "output"; OUT.mkdir(exist_ok=True)
MOD = ROOT / "notebooks" / "analysis" / "output" / "models"
NEG = r"(?:no|denies|denied|without|negative for|not)\s+(?:\w+\s+){0,2}"
PATTERNS = {
    "diabetes": r"diabet|\bDM\b|\bT2DM\b|\bIDDM\b|\bNIDDM\b",
    "af": r"atrial fibrillation|\bafib\b|\ba-?fib\b|\bPAF\b",
    "ckd": r"chronic kidney|\bCKD\b|renal insufficiency|dialysis|\bESRD\b",
    "dialysis": r"dialysis|\bESRD\b|hemodialysis",
    "smoking": r"smok|tobacco",
}


def load_notes():
    n = pd.read_excel(ROOT / "data" / "notes_deidentified.xlsx")
    n = n[n["Signed Status"] != "Deleted"].copy()
    n["year"] = n["Service Date"].astype(int); n["text"] = n.Notes.fillna("").astype(str).str.replace("\n", " ")
    return n.sort_values(["Profile Key", "year"])


FAMILY = r"^\W{0,3}(?:mellitus\W{0,3})?(?:brother|sister|mother|father|son|daughter|grand\w+|aunt|uncle|family)"


def flag(text, pat):
    """Keyword present, not negated in the preceding words, and not a family-history entry ("Diabetes Brother")."""
    for m in re.finditer(pat, text, flags=re.I):
        pre = text[max(0, m.start() - 40): m.start()]; post = text[m.end(): m.end() + 25]
        if re.search(NEG + r"$", pre, flags=re.I) or re.search(FAMILY, post, flags=re.I):
            continue
        return 1
    return 0


def smoking_flag(text):
    """Ever smoking from EHR social-history fields: 'Smoking status: Former/Current' or explicit smoker history."""
    st = [m.group(1).lower() for m in re.finditer(r"smoking status\W{0,3}(former|current|never|every|some)", text, flags=re.I)]
    if any(x in ("former", "current", "every", "some") for x in st):
        return 1
    if st:
        return 0
    if re.search(r"never smok|non-?smoker|never used tobacco", text, flags=re.I):
        return 0
    return int(bool(re.search(r"former smoker|current smoker|tobacco (?:abuse|use disorder)|smoking addiction|\bsmoker\b|pack[- ]years?", text, flags=re.I)))


def baseline_from_notes(notes: pd.DataFrame, patient: str, implant_year) -> dict:
    g = notes[notes["Profile Key"] == patient]
    idx = g[g.year <= implant_year] if pd.notna(implant_year) and (g.year <= implant_year).any() else g.head(1)
    txt = " ".join(idx.text); alltxt = " ".join(g.text)
    # chronic comorbidities are mostly recorded in later history sections, not in operative reports: use all notes (limitation)
    # explicit "[AGE] male/female/gentleman/woman" descriptors first; pronouns (which also refer to physicians) only as fallback
    exp_m = len(re.findall(r"\[AGE\]\s*(?:y/?o\s*|year[- ]old\s*)?(?:male|man|gentleman|M\b)|\bMr\.", alltxt, flags=re.I))
    exp_f = len(re.findall(r"\[AGE\]\s*(?:y/?o\s*|year[- ]old\s*)?(?:female|woman|lady|F\b)|\b(?:Mrs|Ms)\.", alltxt, flags=re.I))
    if exp_m + exp_f > 0:
        he, she = exp_m, exp_f
    else:
        he = len(re.findall(r"\b(he|his|him)\b", alltxt, flags=re.I)); she = len(re.findall(r"\b(she|her|hers)\b", alltxt, flags=re.I))
    age = re.search(r"\b([4-9]\d)[- ]?(?:year[- ]old|yo\b|y/o|y\.o\.)", txt, flags=re.I)
    bmi = re.search(r"\bBMI[^0-9]{0,12}(\d{2}(?:\.\d)?)", alltxt, flags=re.I)
    ac = "none"
    if re.search(r"warfarin|coumadin", txt, flags=re.I): ac = "VKA"
    if re.search(r"apixaban|eliquis|rivaroxaban|xarelto|dabigatran|edoxaban", txt, flags=re.I): ac = "DOAC"
    out = {k: flag(alltxt, p) for k, p in PATTERNS.items() if k != "smoking"}
    out["smoking"] = smoking_flag(alltxt)
    out.update({"sex": "F" if she > he else "M", "age_extracted": float(age.group(1)) if age else np.nan,
                "bmi_extracted": float(bmi.group(1)) if bmi and 15 <= float(bmi.group(1)) <= 60 else np.nan, "anticoag_type": ac})
    return out


def extraction_yield(notes):
    rows = []
    for pk, g in notes.groupby("Profile Key"):
        t = " ".join(g.text)
        rows.append({"patient": pk,
                     "valve_model_named": bool(re.search(r"magna|perimount|carpentier|inspiris|trifecta|mitroflow|biocor|(?<!epic care )\bepic\b(?! care)|mosaic|hancock|freestyle|sapien|\bS3\b|corevalve|evolut", t, flags=re.I)),
                     "size_mm": bool(re.search(r"(?:#\s?|size\s?#?\s?)(19|2[0-9]|3[0-4])\b|\b(19|2[0-9]|3[0-4])\s?-?\s?mm", t, flags=re.I)),
                     "sex": bool(re.search(r"\b(he|she|his|her|male|female)\b", t, flags=re.I)),
                     "numeric_age": bool(re.search(r"\b[4-9]\d[- ]?(?:year[- ]old|yo\b|y/o)", t, flags=re.I)),
                     "bmi": bool(re.search(r"\bBMI[^0-9]{0,12}\d{2}", t, flags=re.I)),
                     "diabetes_mentioned": bool(re.search(PATTERNS["diabetes"], t, flags=re.I)),
                     "af_mentioned": bool(re.search(PATTERNS["af"], t, flags=re.I)),
                     "prosthetic_echo_with_gradient": bool(re.search(r"prosthe\w+.{0,300}?mean gradient", t, flags=re.I)),
                     "notes_in_2plus_years": g.year.nunique() >= 2,
                     "prosthetic_gradient_in_2plus_years": sum(bool(re.search(r"prosthe\w+.{0,300}?mean gradient", " ".join(gy.text), flags=re.I)) for _, gy in g.groupby("year")) >= 2})
    d = pd.DataFrame(rows)
    y = d.drop(columns="patient").mean().rename("share_of_117").to_frame(); y["n"] = d.drop(columns="patient").sum()
    return d, y


def reference(model, size, approach):
    tab = EOA_TABLE[model]; s = min(tab, key=lambda k: abs(k - size)); r = tab[s]
    a_eff = A_EFF.get(model, np.pi * 1.0 ** 2) if approach == "TAVR" else np.pi * 1.0 ** 2
    return {"size_used": s, "table_mpg": r.mpg, "table_eoa": r.eoa, "table_dvi": r.eoa / a_eff}


def label_check(cc):
    rows = []
    for _, r in cc[cc.fu_mpg.notna() | cc.fu_ar.ge(2)].iterrows():
        ref = reference(r.valve_model, r.size_mm, r.approach)
        mpg = r.fu_mpg; d = mpg - ref["table_mpg"] if pd.notna(mpg) else np.nan
        haemo2 = bool(pd.notna(mpg) and mpg >= 20 and d >= 10) or bool(pd.notna(r.fu_ar) and r.fu_ar >= 2)
        haemo3 = bool(pd.notna(mpg) and mpg >= 30 and d >= 20) or bool(pd.notna(r.fu_ar) and r.fu_ar >= 3)
        dvi_fall = (ref["table_dvi"] - r.fu_dvi) if pd.notna(r.fu_dvi) else np.nan
        full2 = haemo2 and (pd.isna(r.fu_dvi) or (dvi_fall >= 0.1 or dvi_fall / ref["table_dvi"] >= 0.2) or (pd.notna(r.fu_ar) and r.fu_ar >= 2))
        rows.append({"patient": r.patient, "group": r.group, "valve": f"{r.valve_model} {int(r.size_mm)}", "fu_mpg": mpg, "table_mpg": round(ref["table_mpg"], 1),
                     "delta_vs_table": None if pd.isna(d) else round(d, 1), "fu_dvi": r.fu_dvi, "table_dvi": round(ref["table_dvi"], 2), "fu_ar": r.fu_ar,
                     "varc3_haemo_stage": 3 if haemo3 else (2 if haemo2 else 0), "varc3_full_positive": int(full2), "documented_svd": int(r.event)})
    return pd.DataFrame(rows)


def two_by_two(df, col):
    tp = int(((df[col] > 0) & (df.documented_svd == 1)).sum()); fn = int(((df[col] == 0) & (df.documented_svd == 1)).sum())
    fp = int(((df[col] > 0) & (df.documented_svd == 0)).sum()); tn = int(((df[col] == 0) & (df.documented_svd == 0)).sum())
    return {"rule": col, "TP": tp, "FN": fn, "FP": fp, "TN": tn, "sensitivity": tp / max(tp + fn, 1), "specificity": tn / max(tn + fp, 1)}


def pv_like(cc, notes):
    rows = []
    for i, r in cc.reset_index(drop=True).iterrows():
        b = baseline_from_notes(notes, r.patient, r.implant_year)
        ref = reference(r.valve_model, r.size_mm, r.approach)
        age = b["age_extracted"] if pd.notna(b["age_extracted"]) else (70.0 if r.approach == "SAVR" else 78.0)
        bmi = b["bmi_extracted"] if pd.notna(b["bmi_extracted"]) else 28.5
        bsa = 2.0 if b["sex"] == "M" else 1.75
        eoai = ref["table_eoa"] / bsa
        rows.append({"valve_id": i + 1, "patient": r.patient, "approach": r.approach, "valve_model": r.valve_model, "valve_class": MODEL_CLASS[r.valve_model],
                     "is_new_model": int(r.valve_model in ("Inspiris Resilia",)), "implant_year": float(r.implant_year) if pd.notna(r.implant_year) else np.nan,
                     "age_at_implant": age, "age_imputed": int(pd.isna(b["age_extracted"])), "sex": b["sex"], "bsa": bsa, "bmi": bmi, "bmi_imputed": int(pd.isna(b["bmi_extracted"])),
                     "diabetes": b["diabetes"], "ckd": b["ckd"], "dialysis": b["dialysis"], "smoking": b["smoking"], "af": b["af"], "anticoag_type": b["anticoag_type"],
                     "label_size_mm": ref["size_used"], "reference_source": "model_size_table", "ref_mean_gradient_mmHg": ref["table_mpg"], "ref_av_area_cm2": ref["table_eoa"],
                     "ref_eoai_cm2_m2": eoai, "ref_dvi": ref["table_dvi"], "ref_ar_grade": 0.0, "ppm_class_ref": int(ppm_class(eoai, bmi)),
                     "table_eoa": ref["table_eoa"], "table_mpg": ref["table_mpg"], "followup_end_years": np.nan})
    return pd.DataFrame(rows)


def harrell(time, score, event):
    return float(concordance_index(time, -score, event))


def boot_ci(fn, n, B=2000, seed=0):
    rng = np.random.default_rng(seed); vals = []
    for _ in range(B):
        idx = rng.integers(0, n, n)
        try:
            v = fn(idx)
            if np.isfinite(v): vals.append(v)
        except Exception:
            pass
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))) if len(vals) > 50 else (np.nan, np.nan)


def model_check(cc, notes):
    res = {}
    rows = []
    for cohort_name, sub in (("primary", cc[cc.group == "primary"]), ("primary_plus_uncertain", cc[cc.group.isin(["primary", "uncertain"])])):
        sub = sub.reset_index(drop=True)
        pv = pv_like(sub, notes)
        I0 = implant_instances(pv)
        t = np.where(sub.event == 1, sub.event_year - sub.implant_year, sub.fu_year - sub.implant_year).astype(float)
        t = np.maximum(t, 0.5); e = sub.event.to_numpy().astype(int)
        dth = joblib.load(MOD / "dth_implant.joblib"); b4 = joblib.load(MOD / "cox_b4.joblib"); cox = joblib.load(MOD / "cox_full.joblib")
        scores = {"implant-time boosted model (10-y CIF)": dth.predict_cif(I0, [10.0])[:, 0],
                  "cause-specific Cox, full baseline (10-y CIF)": cox.predict_cif(I0, [10.0])[:, 0],
                  "published risk-factor Cox, B4 set (10-y CIF)": b4.predict_cif(I0, [10.0])[:, 0]}
        for strat, mask in (("all", np.ones(len(sub), bool)), ("SAVR only", (sub.approach == "SAVR").to_numpy())):
            for name, s in scores.items():
                tt, ss, ee = t[mask], s[mask], e[mask]
                if ee.sum() < 3 or (1 - ee).sum() < 3:
                    continue
                c = harrell(tt, ss, ee); lo, hi = boot_ci(lambda idx: harrell(tt[idx], ss[idx], ee[idx]), len(tt))
                rows.append({"cohort": cohort_name, "stratum": strat, "model": name, "n": int(len(tt)), "events": int(ee.sum()), "harrell_C": c, "C_lo": lo, "C_hi": hi})
        pv.assign(cif10_boosted=scores["implant-time boosted model (10-y CIF)"], cif10_b4=scores["published risk-factor Cox, B4 set (10-y CIF)"],
                  time_years=t, event=e).to_csv(OUT / f"scored_{cohort_name}.csv", index=False)
    res["implant_time"] = pd.DataFrame(rows)
    # landmark model at the follow-up echo (single echo, reference from table); outcome = documented SVD (same or later note)
    sub = cc[(cc.group == "primary") & cc.fu_mpg.notna()].reset_index(drop=True)
    pv = pv_like(sub, notes); base = baseline_frame(pv)
    lt = (sub.fu_year - sub.implant_year).astype(float).clip(lower=0.5).to_numpy()
    arr = lambda v: np.asarray(v, float)[:, None]
    h = Hist(valve_ids=pv.valve_id.to_numpy(), t=arr(lt), vtype=np.ones((len(sub), 1), int), mpg=arr(sub.fu_mpg), vmax=np.full((len(sub), 1), np.nan),
             eoa=np.full((len(sub), 1), np.nan), dvi=arr(sub.fu_dvi), svi=np.full((len(sub), 1), np.nan), lvef=np.full((len(sub), 1), np.nan),
             ar=arr(sub.fu_ar.fillna(0)), n_echoes=np.ones(len(sub), int))
    feats = history_features(h, 1, base)
    m = joblib.load(MOD / "dth_full.joblib"); cif5 = m.predict_cif(feats, [5.0])[:, 0]
    y = sub.event.to_numpy().astype(int)
    auc = float(roc_auc_score(y, cif5)); lo, hi = boot_ci(lambda idx: roc_auc_score(y[idx], cif5[idx]) if 0 < y[idx].sum() < len(idx) else np.nan, len(y))
    lm = sub[["patient", "approach", "valve_model", "size_mm", "fu_mpg", "fu_dvi", "fu_ar", "event"]].assign(years_since_implant=lt, cif_5y=np.round(cif5, 3),
                                                                                                         unconfirmed_flag=feats.unconfirmed_flag.to_numpy())
    res["landmark"] = {"n": int(len(y)), "events": int(y.sum()), "AUC": auc, "AUC_lo": lo, "AUC_hi": hi}
    like = ((sub.approach == "SAVR") & (lt >= 5)).to_numpy()   # all documented failures are surgical valves >= 7 y: compare like with like
    ys, cs = y[like], cif5[like]
    res["landmark_savr_ge5y"] = {"n": int(like.sum()), "events": int(ys.sum()), "AUC": float(roc_auc_score(ys, cs)) if 0 < ys.sum() < len(ys) else None}
    # does the model add anything beyond the single-echo rule at that echo?
    res["landmark_vs_rule"] = {"AUC_rule_unconfirmed_flag": float(roc_auc_score(y, feats.unconfirmed_flag.to_numpy())),
                               "AUC_rule_unconfirmed_flag_savr_ge5y": float(roc_auc_score(ys, feats.unconfirmed_flag.to_numpy()[like])) if 0 < ys.sum() < len(ys) else None}
    res["landmark_table"] = lm.sort_values("cif_5y", ascending=False)
    return res


def main():
    notes = load_notes()
    per_patient, yld = extraction_yield(notes)
    yld.to_csv(OUT / "extraction_yield.csv"); per_patient.to_csv(OUT / "extraction_per_patient.csv", index=False)
    cc = pd.read_csv(HERE / "curated_cohort.csv")
    lc = label_check(cc[cc.group.isin(["primary", "label_only"])]); lc.to_csv(OUT / "label_check.csv", index=False)
    lab_summary = pd.DataFrame([two_by_two(lc, "varc3_haemo_stage"), two_by_two(lc, "varc3_full_positive")]); lab_summary.to_csv(OUT / "label_check_summary.csv", index=False)
    mc = model_check(cc, notes)
    mc["implant_time"].to_csv(OUT / "model_check_implant_time.csv", index=False); mc["landmark_table"].to_csv(OUT / "model_check_landmark.csv", index=False)
    base = pv_like(cc[cc.group == "primary"], notes)
    summary = {"patients_in_notes": int(notes["Profile Key"].nunique()), "primary_cohort": int((cc.group == "primary").sum()), "primary_events": int(cc[cc.group == "primary"].event.sum()),
               "uncertain": int((cc.group == "uncertain").sum()), "label_only": int((cc.group == "label_only").sum()),
               "age_imputed_share": float(base.age_imputed.mean()), "bmi_imputed_share": float(base.bmi_imputed.mean()),
               "baseline_prevalence": base[["diabetes", "af", "ckd", "smoking"]].mean().round(2).to_dict(), "sex_female_share": float((base.sex == "F").mean()),
               "landmark": mc["landmark"], "landmark_savr_ge5y": mc["landmark_savr_ge5y"], "landmark_vs_rule": mc["landmark_vs_rule"]}
    json.dump(summary, open(OUT / "summary.json", "w"), indent=1)
    pd.set_option("display.width", 220)
    print("EXTRACTION YIELD\n", yld.round(2).to_string())
    print("\nLABEL CHECK\n", lc.to_string(index=False), "\n", lab_summary.round(2).to_string(index=False))
    print("\nIMPLANT-TIME MODEL CHECK\n", mc["implant_time"].round(3).to_string(index=False))
    print("\nLANDMARK MODEL AT THE FOLLOW-UP ECHO\n", mc["landmark_table"].to_string(index=False), "\n", json.dumps(mc["landmark"]))
    print("\nSUMMARY\n", json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
