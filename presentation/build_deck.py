"""Build the presentation: slides.pptx (editable) and slides.pdf (same content, rendered with matplotlib).

One declarative slide spec (SLIDES) in inch coordinates on a 13.333 x 7.5 in canvas; two renderers, so the PDF
is a faithful preview of the PPTX. Charts and schematics are regenerated into figures/ from
notebooks/analysis/output/tables/ and from the simulator outputs.

Run:  /data/abar/alexenv/bin/python build_deck.py
"""
from __future__ import annotations

import subprocess
import textwrap
from pathlib import Path

import matplotlib
import matplotlib.font_manager as fm
import numpy as np
import pandas as pd
from PIL import Image

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import patches  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIG = HERE / "figures"; FIG.mkdir(exist_ok=True)
ICO = FIG / "icons"
TAB = ROOT / "notebooks" / "analysis" / "output" / "tables"
W, H = 13.333, 7.5

C = dict(navy="1E2761", panel="243170", teal="028090", light="F4F7FC", border="B9C6E0", body="3A3A3A", muted="6B7590",
         pale="CADCFC", white="FFFFFF", orange="E07A1F", grey="9AA3B8", sep="D7DEEC")
FONT = {"T": "Cambria", "B": "Calibri"}          # PowerPoint names
MPL_FONT = {"T": "Caladea", "B": "Carlito"}      # metric-compatible substitutes for the PDF

TEAM = "CardioNTUA"
MEMBERS = [("Athanasia Karagiannopoulou", "ML Engineer"), ("Anna Panagiotakopoulou", "ML Engineer"),
           ("Evangelie Sintou", "Clinician"), ("Alexandros Barmperis", "ML Engineer")]


# ----------------------------------------------------------------------------------------------- primitives
def rect(x, y, w, h, fill=None, line=None, r=0.0):
    return dict(k="rect", x=x, y=y, w=w, h=h, fill=fill, line=line, r=r)


def circle(cx, cy, d, fill):
    return dict(k="circle", x=cx - d / 2, y=cy - d / 2, w=d, h=d, fill=fill)


def text(x, y, w, h, s, size=14, font="B", bold=False, color="body", align="left", valign="top", italic=False, ls=1.12):
    return dict(k="text", x=x, y=y, w=w, h=h, s=s, size=size, font=font, bold=bold, color=color, align=align, valign=valign, italic=italic, ls=ls)


def image(x, y, w, h, path):
    return dict(k="image", x=x, y=y, w=w, h=h, path=str(path))


def line(x1, y1, x2, y2, color="sep", width=1.0):
    return dict(k="line", x1=x1, y1=y1, x2=x2, y2=y2, color=color, width=width)


def footnote(s, color="muted"):
    return text(0.6, 7.0, 12.1, 0.4, s, size=10, color=color)


def title(s, sub=None, dark=False):
    out = [text(0.6, 0.42, 12.0, 0.8, s, size=32, font="T", bold=True, color="white" if dark else "navy")]
    if sub:
        out.append(text(0.6, 1.15, 12.0, 0.45, sub, size=15, color="pale" if dark else "muted", italic=True))
    return out


def lead(x, y, w, head, body, size=14, color="body", head_color="navy"):
    """A short bold lead line followed by one or two plain lines. No box."""
    return [text(x, y, w, 0.4, head, size=size, bold=True, color=head_color),
            text(x, y + 0.38, w, 1.2, body, size=size, color=color)]


# ----------------------------------------------------------------------------------------------- charts
def _style():
    for name in ("Carlito", "Caladea"):
        for p in subprocess.run(["fc-list", f":family={name}", "file"], capture_output=True, text=True).stdout.split("\n"):
            p = p.strip().rstrip(":")
            if p.endswith((".ttf", ".otf")):
                fm.fontManager.addfont(p)
    plt.rcParams.update({"font.family": "Carlito", "font.size": 11, "axes.edgecolor": "#" + C["muted"], "axes.labelcolor": "#" + C["body"],
                         "xtick.color": "#" + C["body"], "ytick.color": "#" + C["body"], "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": "#E4E8F0", "grid.linewidth": 0.8, "axes.axisbelow": True, "figure.facecolor": "white"})


def hx(k): return "#" + C[k]


def chart_auc():
    m = pd.read_csv(TAB / "metrics.csv"); m = m[m.tau == 5.0]
    rows = [("dth_full", "all", "Dynamic model, all echoes", 1), ("dth_full", "currently_negative", "Dynamic model, echoes not yet abnormal", 1),
            ("dth_registry", "currently_negative", "Registry-only features (not yet abnormal)", 0), ("dth_time_only", "currently_negative", "Time since implant only (not yet abnormal)", 0),
            ("dth_implant", "implant_time", "Implant-time boosted model", 1), ("cox_full", "implant_time", "Cause-specific Cox, full baseline", 0),
            ("cox_b4", "implant_time", "Published risk-factor Cox", 0)]
    fig, ax = plt.subplots(figsize=(6.2, 3.3), dpi=200)
    ys = np.arange(len(rows))[::-1]
    for yy, (mod, sub, lab, ours) in zip(ys, rows):
        r = m[(m.model == mod) & (m.subset == sub)].iloc[0]
        ax.barh(yy, r.AUC - 0.5, left=0.5, height=0.55, color=hx("navy") if ours else hx("grey"))
        ax.plot([r.AUC_lo, r.AUC_hi], [yy, yy], color=hx("body"), lw=1.2)
        ax.text(r.AUC_hi + 0.008, yy, f"{r.AUC:.3f}", va="center", ha="left", fontsize=9.5, color=hx("body"))
    ax.set_yticks(ys); ax.set_yticklabels([r[2] for r in rows], fontsize=9.5)
    ax.axvline(0.75, color=hx("teal"), ls="--", lw=1.2); ax.text(0.752, ys[-1] - 0.62, "pre-specified 0.75", color=hx("teal"), fontsize=8.5, va="center")
    ax.set_xlim(0.5, 1.0); ax.set_ylim(-0.9, len(rows) - 0.4); ax.set_xlabel("AUC at 5 years (95% valve-bootstrap CI)")
    ax.grid(axis="y", visible=False)
    ax.plot([], [], color=hx("navy"), lw=6, label="proposed"); ax.plot([], [], color=hx("grey"), lw=6, label="comparator")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "chart_auc.png"); plt.close(fig)


def chart_policy():
    d = pd.read_csv(TAB / "policy_metrics.csv")
    d["delay_m"] = d.confirmed_delay_mean_y * 12
    fig, ax = plt.subplots(figsize=(5.9, 3.3), dpi=200)
    ra = d[d.policy.str.startswith("risk_adapted")].sort_values("tau")
    ax.plot(ra.echoes_per_patient_decade, ra.delay_m, color=hx("orange"), lw=2, marker="o", ms=6, zorder=3, label="risk-adapted schedule (threshold sweep)")
    for _, r in ra.iterrows():
        if r.tau in (0.005, 0.01, 0.02, 0.05):
            ax.annotate(f"{r.tau*100:g}%", (r.echoes_per_patient_decade, r.delay_m), xytext=(6, -11 if r.tau == 0.01 else 5), textcoords="offset points", fontsize=8.5, color=hx("body"))
    fixed = {"guideline": "guideline", "guideline_confirm": "guideline + confirmation", "deescalation": "de-escalation", "time_only": "time since implant only"}
    for key, lab in fixed.items():
        r = d[d.policy == key].iloc[0]
        ax.plot(r.echoes_per_patient_decade, r.delay_m, "s", color=hx("navy"), ms=8, zorder=3, label="fixed schedules" if key == "guideline" else None)
        dx, dy = {"guideline": (-8, 6), "guideline_confirm": (-8, -14), "deescalation": (-8, 8), "time_only": (8, -2)}[key]
        ax.annotate(lab, (r.echoes_per_patient_decade, r.delay_m), xytext=(dx, dy), textcoords="offset points", fontsize=8.5, ha="right" if dx < 0 else "left", color=hx("body"))
    o = d[d.policy == "oracle"].iloc[0]
    ax.plot(o.echoes_per_patient_decade, o.delay_m, "*", color=hx("teal"), ms=13, zorder=3, label="oracle (sees the truth)")
    ax.annotate("oracle", (o.echoes_per_patient_decade, o.delay_m), xytext=(8, -3), textcoords="offset points", fontsize=8.5, color=hx("body"))
    gc = d[d.policy == "guideline_confirm"].iloc[0]
    ax.axhline(gc.delay_m + 3, color=hx("muted"), ls=":", lw=1.2); ax.text(2.0, gc.delay_m + 3.15, "+3-month bound vs guideline + confirmation", fontsize=8, color=hx("muted"))
    r1 = ra[ra.tau == 0.01].iloc[0]
    ax.annotate("1% rule: 44% fewer echoes,\n+0.7 months mean delay", (r1.echoes_per_patient_decade, r1.delay_m), xytext=(4.9, 11.9), fontsize=8.5, color=hx("orange"),
                arrowprops=dict(arrowstyle="-", color=hx("orange"), lw=1))
    ax.set_xlabel("echoes per patient-decade"); ax.set_ylabel("mean delay to confirmed detection (months)")
    ax.set_xlim(1.5, 9.2); ax.set_ylim(10.5, 19.2); ax.legend(loc="upper left", frameon=False, fontsize=8.5)
    fig.tight_layout(); fig.savefig(FIG / "chart_policy.png"); plt.close(fig)


def chart_shap():
    s = pd.read_csv(TAB / "shap_importance.csv").head(8)
    names = {"approach": "Approach (SAVR / TAVR)", "valve_model": "Valve model", "period": "Prediction horizon", "delta_mpg": "Gradient change from reference",
             "landmark_t": "Time since implant", "age_at_implant": "Age at implant", "vmax_now": "Peak velocity", "stage_now_haemo": "Current haemodynamic stage",
             "max_mpg_sofar": "Highest gradient so far", "stage_now_full": "Current VARC-3 stage", "mpg_now": "Current mean gradient"}
    fig, ax = plt.subplots(figsize=(5.2, 2.6), dpi=200)
    ys = np.arange(len(s))[::-1]
    ax.barh(ys, s.mean_abs_shap, height=0.6, color=hx("navy"))
    for yy, v in zip(ys, s.mean_abs_shap):
        ax.text(v + 0.006, yy, f"{v:.2f}", va="center", fontsize=8.5, color=hx("body"))
    ax.set_yticks(ys); ax.set_yticklabels([names.get(f, f) for f in s.feature], fontsize=9)
    ax.set_xlabel("mean |SHAP| on the SVD hazard (log-odds)", fontsize=9); ax.grid(axis="y", visible=False); ax.set_xlim(0, 0.56)
    ax.set_title("What drives the prediction (test set)", fontsize=10, color=hx("navy"), loc="left")
    fig.tight_layout(); fig.savefig(FIG / "chart_shap.png"); plt.close(fig)


def chart_curves():
    from PIL import ImageDraw, ImageFont
    src = Image.open(ROOT / "notebooks" / "simulator" / "output" / "figures" / "curves_sim_vs_target.png")
    w, h = src.size
    crop = src.crop((0, 0, w, int(h * 0.335))).convert("RGB"); d = ImageDraw.Draw(crop)
    d.rectangle((0, 0, w, 34), fill="white")
    font_path = [ln.split(":")[0] for ln in subprocess.run(["fc-list", ":family=Carlito", "file", "style"], capture_output=True, text=True).stdout.split("\n") if "Bold" in ln and "Italic" not in ln]
    fnt = ImageFont.truetype(font_path[0], 22) if font_path else ImageFont.load_default()
    titles = ["Kermen 2022: freedom from VARC-3 stage 2/3 SVD (%)", "Wakami 2022: Trifecta SVD (%)", "NOTION 10-y: TAVI moderate SVD (%)"]
    for i, t in enumerate(titles):
        d.text((w * (2 * i + 1) / 6, 18), t, fill=hx("navy"), font=fnt, anchor="mm")
    crop.save(FIG / "chart_curves.png")


# ----------------------------------------------------------------------------------------------- schematics
def _clean(ax, xlim, ylim=None, xlabel=None):
    ax.grid(False); ax.set_xlim(*xlim)
    if ylim: ax.set_ylim(*ylim)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.set_yticks([])
    if xlabel: ax.set_xlabel(xlabel, fontsize=9.5)


def fig_flicker():
    """One valve: true gradient, noisy annual echoes, the VARC-3 threshold, a false alarm and the late confirmation."""
    t = np.linspace(0, 10, 400); ref = 11.0
    true = ref + 0.25 * t + 24 / (1 + np.exp(-(t - 7.5) / 0.8))
    thr = max(20.0, ref + 10)
    echo_t = np.array([0.1] + list(range(1, 11)), float)
    echo_v = np.array([ref, 9.6, 12.4, 22.0, 13.5, 12.6, 17.5, 21.5, 26.0, 31.8, 34.5])
    fig, ax = plt.subplots(figsize=(6.6, 3.5), dpi=200)
    ax.plot(t, true, color=hx("navy"), lw=2, label="true mean gradient")
    ax.axhline(thr, color=hx("muted"), ls="--", lw=1.1); ax.text(0.15, thr - 1.7, "VARC-3 moderate threshold: ≥ 20 mmHg and ≥ 10 above reference", fontsize=8.5, color=hx("muted"))
    ax.plot(echo_t, echo_v, "o", color=hx("teal"), ms=6.5, zorder=4, label="annual echo reading")
    for k, lab, dx, dy, ha in ((3, "isolated reading above the\nthreshold: reverses next year", 0, 12, "center"), (7, "first reading\nabove threshold", -10, 8, "right"), (8, "confirmed", 8, 4, "left")):
        ax.plot(echo_t[k], echo_v[k], "o", ms=11, mfc="none", mec=hx("orange"), mew=1.6, zorder=5)
        ax.annotate(lab, (echo_t[k], echo_v[k]), xytext=(dx, dy), textcoords="offset points", fontsize=8.5, color=hx("orange"), ha=ha)
    tc = t[np.argmax(true >= thr)]
    ax.axvline(tc, color=hx("navy"), ls=":", lw=1.2); ax.text(tc - 0.12, 39, "true crossing", rotation=90, fontsize=8.5, color=hx("navy"), va="top", ha="right")
    ax.annotate("", xy=(8, 6.0), xytext=(tc, 6.0), arrowprops=dict(arrowstyle="<->", color=hx("body"), lw=1)); ax.text((tc + 8) / 2, 7.0, "delay to\nconfirmation", ha="center", fontsize=8.5, color=hx("body"))
    ax.set_xlabel("years after implant", fontsize=9.5); ax.set_ylabel("mean gradient (mmHg)", fontsize=9.5)
    ax.set_xlim(-0.2, 10.4); ax.set_ylim(3, 40); ax.grid(False)
    ax.legend(loc="upper left", frameon=False, fontsize=8.5)
    fig.tight_layout(); fig.savefig(FIG / "fig_flicker.png"); plt.close(fig)


def fig_schedule():
    """The guideline calendar and what happens to it in practice."""
    fig, ax = plt.subplots(figsize=(6.6, 2.7), dpi=200)
    yrs = [0.25] + list(range(1, 9))
    for row, ylab in ((1.0, "Guideline\n(ESC/EACTS 2025)"), (0.0, "In practice")):
        ax.plot([0, 8.6], [row, row], color=hx("navy"), lw=2, solid_capstyle="round", zorder=1)
        ax.text(-0.35, row, ylab, ha="right", va="center", fontsize=9.5, color=hx("navy"), fontweight="bold")
    for y in yrs:
        ax.plot(y, 1.0, "o", color=hx("teal"), ms=9, zorder=3)
    ax.text(8.75, 1.0, "…", fontsize=13, color=hx("navy"), va="center")
    for y in yrs:
        if y == 0.25: ax.plot(y, 0, "o", ms=9, mfc="white", mec=hx("grey"), mew=1.5, zorder=3)
        elif y == 4: ax.plot(y, 0, "o", ms=9, mfc="white", mec=hx("grey"), mew=1.5, ls="", zorder=3)
        elif y == 6: ax.plot(y, 0, "o", color=hx("orange"), ms=9, zorder=3)
        else: ax.plot(y, 0, "o", color=hx("teal"), ms=9, zorder=3)
    ax.text(8.75, 0, "…", fontsize=13, color=hx("navy"), va="center")
    for y, lab in zip(yrs, ["≤3 mo", "1 y", "2 y", "3 y", "4 y", "5 y", "6 y", "7 y", "8 y"]):
        ax.text(y, 1.22, lab, ha="center", fontsize=8.5, color=hx("body"))
    ax.annotate("no reference echo\nin ~30% of patients", (0.25, 0), xytext=(0.25, -0.75), fontsize=8.5, ha="center", color=hx("muted"), arrowprops=dict(arrowstyle="-", color=hx("grey"), lw=0.8))
    ax.annotate("visit missed\n(35% by year 4)", (4, 0), xytext=(4, -0.75), fontsize=8.5, ha="center", color=hx("muted"), arrowprops=dict(arrowstyle="-", color=hx("grey"), lw=0.8))
    ax.annotate("abnormal reading:\nno confirmation step", (6, 0), xytext=(6.6, -0.75), fontsize=8.5, ha="center", color=hx("orange"), arrowprops=dict(arrowstyle="-", color=hx("orange"), lw=0.8))
    ax.text(4.3, 0.55, "same interval whatever the valve, the age or the trend", ha="center", fontsize=8.5, color=hx("muted"), style="italic")
    _clean(ax, (-2.6, 9.2), (-1.15, 1.5)); ax.set_xticks([]); ax.spines["bottom"].set_visible(False)
    fig.tight_layout(); fig.savefig(FIG / "fig_schedule.png"); plt.close(fig)


def fig_adaptive():
    """Two valves under a risk-adapted schedule: intervals lengthen when stable and shorten when the gradient rises."""
    fig, axes = plt.subplots(2, 1, figsize=(6.6, 3.2), dpi=200, sharex=True)
    t = np.linspace(0, 10, 300)
    stable = 11 + 0.2 * t
    rising = 11 + 0.2 * t + np.where(t > 4, 2.2 * (np.exp(0.5 * (t - 4)) - 1), 0)
    for ax, g, echoes, lab, extra in ((axes[0], stable, [0.25, 1, 3, 6, 9], "Stable valve: intervals lengthen", None),
                                       (axes[1], rising, [0.25, 1, 3, 5, 5.5, 6.5, 7], "Rising gradient: intervals shorten, then a confirmation echo", 6.0)):
        ax.plot(t, g, color=hx("navy"), lw=1.8)
        ev = np.interp(echoes, t, g)
        ax.plot(echoes, ev, "o", color=hx("teal"), ms=7, zorder=3)
        for a, b in zip(echoes[:-1], echoes[1:]):
            ax.annotate("", xy=(b, 4), xytext=(a, 4), arrowprops=dict(arrowstyle="<->", color=hx("grey"), lw=0.9, shrinkA=0, shrinkB=0))
            ax.text((a + b) / 2, 5.2, f"{12*(b-a):.0f} mo", ha="center", fontsize=7.5, color=hx("muted"))
        if extra:
            ax.annotate("first abnormal echo\n→ confirmation echo 6 months later", (extra, np.interp(extra, t, g)), xytext=(0.3, 24), fontsize=8, color=hx("orange"),
                        arrowprops=dict(arrowstyle="-", color=hx("orange"), lw=0.8))
        ax.text(0.05, 0.92, lab, transform=ax.transAxes, fontsize=9.5, color=hx("navy"), fontweight="bold", va="top")
        ax.set_ylim(1, 42); ax.set_xlim(-0.2, 10.3); ax.grid(False); ax.set_yticks([10, 20, 30]); ax.tick_params(labelsize=8)
        ax.axhline(20, color=hx("muted"), ls="--", lw=0.8)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
    axes[1].set_xlabel("years after implant", fontsize=9.5); axes[0].set_ylabel("mmHg", fontsize=8.5); axes[1].set_ylabel("mmHg", fontsize=8.5)
    fig.tight_layout(h_pad=0.4); fig.savefig(FIG / "fig_adaptive.png"); plt.close(fig)


def fig_landmark():
    """Where a prediction is made, what it may use, and when the outcome is dated."""
    fig, ax = plt.subplots(figsize=(6.6, 2.6), dpi=200)
    ax.plot([0, 9.6], [0, 0], color=hx("navy"), lw=2, zorder=1)
    ax.axvspan(4, 9, ymin=0.35, ymax=0.62, color=hx("light"), zorder=0)
    ax.text(6.5, 0.42, "prediction window: risk over the next 0.5–5 y", ha="center", fontsize=8.5, color=hx("navy"))
    pts = {0: ("implant", "navy"), 0.25: ("reference echo", "teal"), 1: ("", "teal"), 2: ("", "teal"), 3: ("", "teal"), 4: ("landmark echo", "orange"),
           6: ("first abnormal echo", "grey"), 7: ("confirming echo", "navy"), 8.6: ("death (competing)", "muted")}
    for x, (lab, col) in pts.items():
        if x == 8.6: ax.plot(x, 0, "x", color=hx(col), ms=9, mew=2, zorder=3)
        elif x == 0: ax.plot(x, 0, "|", color=hx(col), ms=14, mew=2, zorder=3)
        else: ax.plot(x, 0, "o", color=hx(col), ms=8 if x not in (4, 7) else 10, zorder=3)
        if lab: ax.text(x, -0.22 if x not in (0.25, 6) else -0.40, lab, ha="center", fontsize=8.5, color=hx(col) if col != "grey" else hx("muted"))
    ax.annotate("prediction made here,\nusing only echoes up to this point", (4, 0), xytext=(2.2, 0.95), fontsize=8.5, color=hx("orange"), ha="center",
                arrowprops=dict(arrowstyle="-", color=hx("orange"), lw=0.9))
    ax.annotate("outcome dated here:\nthe label exists only when confirmed", (7, 0), xytext=(8.0, 0.95), fontsize=8.5, color=hx("navy"), ha="center",
                arrowprops=dict(arrowstyle="-", color=hx("navy"), lw=0.9))
    ax.annotate("", xy=(3.95, 0.18), xytext=(0.0, 0.18), arrowprops=dict(arrowstyle="-|>", color=hx("teal"), lw=1)); ax.text(2, 0.24, "features: history so far", ha="center", fontsize=8, color=hx("teal"))
    _clean(ax, (-0.4, 10), (-0.6, 1.3), "years after implant"); ax.set_xticks(range(0, 10)); ax.tick_params(labelsize=8)
    fig.tight_layout(); fig.savefig(FIG / "fig_landmark.png"); plt.close(fig)


def fig_model():
    """Flow of the discrete-time cause-specific hazard model, from features to the next-echo decision."""
    fig = plt.figure(figsize=(6.6, 3.3), dpi=200); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 6.6); ax.set_ylim(0, 3.3); ax.axis("off"); ax.grid(False)
    def box(x, y, w, h, s, fill="light", fc="navy", fs=8.5, bold=False):
        ax.add_patch(patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.08", facecolor=hx(fill), edgecolor="none"))
        ax.text(x + w / 2, y + h / 2, s, ha="center", va="center", fontsize=fs, color=hx(fc), fontweight="bold" if bold else "normal", linespacing=1.15)
    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=hx("teal"), lw=1.3))
    box(0.15, 2.25, 1.75, 0.85, "Baseline\nage, comorbidities,\nvalve model, size, mismatch")
    box(0.15, 1.25, 1.75, 0.85, "Reference echo\n(or model × size norm,\nflagged)")
    box(0.15, 0.25, 1.75, 0.85, "Trajectory\nΔ gradient, slope, prior\nabnormal echoes, EOA / DVI")
    box(2.35, 1.05, 1.6, 1.25, "Gradient-boosted\nmultinomial classifier\n\nper 6-month period:\nhazard of SVD,\nhazard of competing event", fill="navy", fc="white", fs=8)
    for y in (2.67, 1.67, 0.67): arrow(1.9, y, 2.35, 1.67 if y == 1.67 else (2.2 if y > 1.67 else 1.15))
    # cumulative incidence inset
    ins = fig.add_axes([0.63, 0.42, 0.2, 0.42]); tt = np.linspace(0, 5, 100); cif = 1 - np.exp(-0.035 * tt ** 1.6)
    ins.plot(tt, cif, color=hx("navy"), lw=1.8); ins.fill_between(tt, 0, cif, color=hx("light"))
    ins.axhline(0.01, color=hx("orange"), ls="--", lw=1); ins.axvline(2.0, color=hx("orange"), ls=":", lw=1)
    ins.set_xlim(0, 5); ins.set_ylim(0, 0.16); ins.set_xticks([0, 1, 2, 3, 4, 5]); ins.set_yticks([0, 0.05, 0.10, 0.15]); ins.tick_params(labelsize=6.5); ins.grid(False)
    ins.set_title("cumulative incidence of\nconfirmed SVD after this echo", fontsize=7.5, color=hx("navy")); ins.set_xlabel("years", fontsize=7)
    ins.text(2.1, 0.013, "1% rule", fontsize=6.5, color=hx("orange"))
    arrow(3.95, 1.67, 4.15, 1.67)
    box(4.25, 0.25, 2.2, 0.85, "Next echo: longest of 6–36 months\nwith risk ≤ 1%; risk tier, top drivers;\nany abnormal echo → confirmation echo", fill="light", fs=7.8)
    arrow(5.35, 1.35, 5.35, 1.12)
    fig.savefig(FIG / "fig_model.png"); plt.close(fig)


def chart_calibration():
    d = pd.read_csv(TAB / "calibration_5y_dth_full.csv")
    fig, ax = plt.subplots(figsize=(3.4, 2.3), dpi=200)
    ax.plot([0, 1], [0, 1], color=hx("grey"), lw=1, ls="--")
    ax.plot(d.pred, d.obs, "o-", color=hx("navy"), ms=5, lw=1.5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1]); ax.tick_params(labelsize=8)
    ax.set_xlabel("predicted 5-y risk (decile mean)", fontsize=8.5); ax.set_ylabel("observed", fontsize=8.5)
    ax.set_title("Calibration at 5 y, slope 1.04", fontsize=9.5, color=hx("navy"), loc="left")
    fig.tight_layout(); fig.savefig(FIG / "chart_calibration.png"); plt.close(fig)


def chart_samplesize():
    d = pd.read_csv(TAB / "oracle_sample_size.csv")
    fig, ax = plt.subplots(figsize=(3.4, 2.3), dpi=200)
    ax.plot(d.n_valves, d.ci_width_mean, "o-", color=hx("navy"), ms=5, lw=1.5)
    for n, w in zip(d.n_valves, d.ci_width_mean):
        ax.annotate(f"{w:.3f}", (n, w), xytext=(6, 4), textcoords="offset points", fontsize=8, color=hx("body"))
    ax.set_xscale("log"); ax.set_xticks([500, 1000, 2000, 5000]); ax.set_xticklabels(["500", "1,000", "2,000", "5,000"]); ax.tick_params(labelsize=8)
    ax.set_ylim(0, 0.12); ax.set_xlim(380, 8000); ax.minorticks_off()
    ax.set_xlabel("valves in the cohort (simulation)", fontsize=8.5); ax.set_ylabel("95% CI width, 5-y AUC", fontsize=8.5)
    ax.set_title("Sample size: 5,000 valves for ±0.01", fontsize=9.5, color=hx("navy"), loc="left")
    fig.tight_layout(); fig.savefig(FIG / "chart_samplesize.png"); plt.close(fig)


def make_figures():
    _style(); chart_auc(); chart_policy(); chart_shap(); chart_curves(); chart_calibration(); chart_samplesize()
    fig_flicker(); fig_schedule(); fig_adaptive(); fig_landmark(); fig_model()


# ----------------------------------------------------------------------------------------------- slide spec
def slide_title():
    names = "   ·   ".join(f"{n} ({r})" for n, r in MEMBERS)
    return dict(bg="navy", items=[
        circle(12.85, 1.05, 6.5, "panel"), circle(0.25, 7.55, 5.5, "panel"),
        image(6.07, 1.0, 1.2, 1.2, ICO / "icon_s1_1.png"),
        text(0.8, 2.5, 11.7, 1.0, "Knowing When to Look Again", size=40, font="T", bold=True, color="white", align="center", valign="middle"),
        text(1.4, 3.5, 10.5, 0.8, "Dynamic risk of structural valve deterioration to personalise echo surveillance after bioprosthetic aortic valve replacement",
             size=18, color="pale", align="center", valign="middle"),
        text(1.4, 4.75, 10.5, 0.5, f"Team {TEAM}", size=16, bold=True, color="white", align="center", valign="middle"),
        text(0.9, 5.25, 11.5, 0.8, names, size=13.5, color="pale", align="center", valign="middle"),
        text(1.4, 6.65, 10.5, 0.4, "Dyania Health Hackathon  ·  AVR durability challenge  ·  17 September 2026", size=13, color="pale", align="center"),
    ], notes="""Good [morning/afternoon] everyone, thank you for having us. We are team CardioNTUA: Evangelie Sintou, our clinician, and Athanasia Karagiannopoulou, Anna Panagiotakopoulou and Alexandros Barmperis, ML engineers.
We are going to talk about a question that is not asked enough after a patient's aortic valve replacement: once we have chosen a bioprosthetic valve, how do we know when to look again?
Our project is a model that predicts a patient's risk of structural valve deterioration at every echo and uses it to turn follow-up into a personalised surveillance plan instead of a one-size-fits-all calendar.""")


def slide_problem():
    return dict(bg="white", items=[
        *title("The clinical problem", "Deterioration is slow, silent, and hard to confirm"),
        image(0.5, 1.75, 6.9, 3.7, FIG / "fig_flicker.png"),
        text(0.6, 5.5, 6.7, 0.7, "Schematic of one valve: the true gradient (line), its annual echo readings (dots) and the VARC-3 threshold. Measurement noise creates false alarms; a label needs two consecutive abnormal echoes and arrives late.",
             size=10.5, color="muted", italic=True),
        *lead(7.8, 1.8, 4.9, "Durability differs by valve", "Severe SVD at 10 years: 10.0% after SAVR vs 1.5% after TAVI (NOTION); the SAVR–TAVR question is still open."),
        *lead(7.8, 3.4, 4.9, "The label is noisy", "59–65% of first SVD calls are reversed at the next echo (Velders 2024)."),
        *lead(7.8, 4.75, 4.9, "Confirmation is late", "Under annual surveillance, moderate deterioration is confirmed a mean 1.5 years after onset; 2.9 years for one patient in ten (our calibrated simulation)."),
        footnote("Praz 2025 Eur Heart J (ESC/EACTS)  ·  Thyregod 2024 Eur Heart J (NOTION 10-year)  ·  Velders 2024 JTCVS Open  ·  Kostyunin 2020 JAHA"),
    ], notes="""SVD is the permanent, intrinsic deterioration of a bioprosthetic valve: fibrosis, calcification, tearing. It is why durability is finite in a way it is not for mechanical valves.
The figure shows the difficulty. The valve's true gradient rises slowly; each annual echo reads it with noise, so an isolated reading can cross the threshold and reverse the next year. Velders found that 59 to 65 percent of first SVD calls do not hold at the next echo. A label therefore needs two consecutive abnormal echoes, and in our calibrated simulation that confirmation comes a mean 1.5 years after the true crossing, almost 3 years for one patient in ten.
Durability also differs by valve: NOTION shows a real gap between SAVR and TAVI at ten years. A noisy signal, confirmed slowly. That is the problem.""")


def slide_workflow():
    return dict(bg="white", items=[
        *title("Current workflow", "One calendar for every patient"),
        image(0.5, 1.7, 7.0, 2.9, FIG / "fig_schedule.png"),
        text(0.6, 4.65, 6.8, 0.5, "Guideline schedule: echo within 3 months, at 1 year, then every year (ESC/EACTS 2025); ACC/AHA similar.", size=10.5, color="muted", italic=True),
        *lead(7.9, 1.8, 4.8, "Same echo for everyone", "The interval does not depend on the valve model, the patient's age or the trend of the readings."),
        *lead(7.9, 3.35, 4.8, "Timing errors are expensive", "Urgent redo surgery: 22.6% mortality vs 1.4% when elective (Cartlidge 2019)."),
        *lead(7.9, 4.75, 4.8, "Effort in the wrong places", "Low-risk patients are over-tested while valve-in-valve volume heads towards 42,000 a year by 2035 (Genereux 2025)."),
        footnote("Praz 2025 Eur Heart J (ESC/EACTS)  ·  Otto 2021 Circulation (ACC/AHA)  ·  Cartlidge 2019 JACC  ·  Genereux 2025 Struct Heart"),
    ], notes="""If the signal is this noisy and this slow to confirm, you would hope the surveillance built around it was adaptive. It is not.
The guideline calls for the same echo, within three months, at one year and then every year, for essentially every patient, regardless of valve, age, or how the valve is trending. In practice the calendar leaks: about 35 percent of patients have missed their visit by year four, about 30 percent have no reference echo to measure change against, and there is no rule for confirming an abnormal reading before acting on it.
The cost is not abstract: urgent redo surgery carries 22.6 percent mortality against 1.4 percent electively. At the same time we over-test low-risk patients while valve-in-valve volumes climb. We spend surveillance effort in the wrong places.""")


def slide_hypothesis():
    return dict(bg="white", items=[
        *title("Hypothesis", "Let the interval follow the risk"),
        text(0.6, 1.75, 12.1, 0.9, "A risk estimate updated at every echo can lengthen the surveillance interval for stable valves and shorten it for valves that begin to deteriorate, without losing time to detection.",
             size=17, color="navy"),
        image(0.5, 2.85, 7.0, 3.3, FIG / "fig_adaptive.png"),
        text(7.9, 2.9, 4.8, 0.45, "Targets, fixed before any analysis", size=14, bold=True, color="navy"),
        text(7.9, 3.35, 4.8, 2.6,
             "≥ 25% fewer echoes per patient-decade\n\n≤ 3 months added detection delay, on average and for the slowest-detected tenth\n\nEarlier confirmation for the high-risk valves where timing matters most\n\nAny single abnormal echo still triggers a confirmation echo; symptoms always override the schedule.",
             size=14),
        text(0.6, 6.2, 6.8, 0.6, "Schematic: two valves under a risk-adapted schedule. Dots are echoes; arrows show the chosen interval.", size=10.5, color="muted", italic=True),
    ], notes="""Clinicians already exercise judgement at the point of valve choice. Our hypothesis is that the same personalised thinking should extend into follow-up: a risk model, updated at every echo, that lengthens intervals for valves that are stable and shortens them for valves that are starting to deteriorate.
The figure shows the idea: a stable valve is seen at 1, 3, 6 and 9 years; a valve whose gradient starts to rise is seen more often, and the first abnormal reading triggers a confirmation echo.
We fixed the targets before running anything: at least 25 percent fewer echoes, no more than three months of added detection delay, on average and for the slowest-detected tenth, and earlier confirmation for the high-risk valves.
[Hand over to the engineer for the study design.]""")


def slide_study():
    col = lambda x, head, body: [text(x, 1.8, 3.8, 0.45, head, size=14.5, bold=True, color="navy"), text(x, 2.3, 3.8, 4.5, body, size=13.5)]
    return dict(bg="white", items=[
        *title("Study design", "Who is in, what counts as deterioration, and what the truth is"),
        *col(0.6, "Population",
             "Adults with a first bioprosthetic SAVR or TAVR; valve model and size recorded; reference echo within 3 months (or a model × size norm, flagged); at least one follow-up echo.\n\n"
             "Excluded: mechanical, homograft and Ross valves; index valve-in-valve; multiple prostheses.\n\n"
             "≥ 5,000 valves from ≥ 3 centres; ≥ 2,000 external."),
        *col(4.75, "Endpoints",
             "Primary: VARC-3 moderate SVD (stage ≥ 2) confirmed on two consecutive echoes; death, endocarditis and non-SVD reintervention compete; incidence at 5, 8 and 10 years.\n\n"
             "Secondary: valve failure; reintervention, elective vs urgent; gradient progression; echoes used and detection delay under each schedule."),
        *col(8.9, "Ground truth, highest level wins",
             "1.  Explant or autopsy pathology\n2.  Reintervention for SVD\n3.  VARC-3 criteria on two consecutive echoes\n4.  Blinded heart-team adjudication\n\n"
             "Endocarditis, thrombosis, paravalvular leak and mismatch are adjudicated out."),
        line(4.5, 1.9, 4.5, 6.6), line(8.65, 1.9, 8.65, 6.6),
        footnote("Généreux 2021 Eur Heart J (VARC-3)  ·  Capodanno 2017 Eur Heart J  ·  Dvir 2018 Circulation  ·  Praz 2025 (Table 12 adopts the same criteria)"),
    ], notes="""The population is any adult with a first bioprosthetic valve, surgical or transcatheter, with the valve model recorded and at least one follow-up echo. We exclude what would confound gradients: mechanical valves, valve-in-valve as the index procedure, multiple prostheses.
The primary endpoint is VARC-3 moderate SVD confirmed on two consecutive echoes, with death as a competing risk. Confirmation is what removes the label noise from the second slide.
Ground truth is a hierarchy: pathology, then a reintervention for SVD, then the echo criteria, then adjudication. And we fixed the usefulness bar before running anything: AUC 0.75, calibration slope within 0.8 to 1.2, and the surveillance target from the previous slide.""")


def slide_data():
    rows = [("Echo reporting database", "serial gradients, EOA, DVI, LVEF"), ("STS ACSD, STS/ACC TVT Registry", "valve model, size, comorbidities"),
            ("EHR notes, local NLP", "atrial fibrillation, anticoagulation, symptoms"), ("CMS claims, National Death Index", "reintervention and death elsewhere")]
    items = [*title("Data", "Registry-shaped for deployment; a calibrated simulation for the prototype"),
             text(0.6, 1.8, 5.6, 0.45, "Deployment: routinely collected, linked by one key", size=14, bold=True, color="navy")]
    y = 2.3
    for src, what in rows:
        items += [text(0.6, y, 2.6, 0.65, src, size=12.5, bold=True, color="body"), text(3.25, y, 3.0, 0.65, what, size=12.5), line(0.6, y + 0.68, 6.2, y + 0.68)]
        y += 0.78
    items += [text(0.6, 5.5, 5.6, 1.1, "Needs exact dates and ages (a HIPAA limited data set). In the provided notes, a prosthetic echo with a gradient was found for only 36% of patients.", size=12, color="muted"),
              text(6.7, 1.8, 6.1, 0.45, "Prototype: no public serial-echo dataset, so we built one", size=14, bold=True, color="navy")]
    bx = [(6.7, "Published curves\nKermen, NOTION, Wakami"), (8.25, "Patient-level data\nreconstructed (Guyot)"),
          (9.8, "Simulator calibrated\nASE 2024 physics,\nVelders 2024 noise"), (11.35, "10,000 valves,\n50,849 echoes,\ntrue state kept aside")]
    for i, (x, s) in enumerate(bx):
        items += [rect(x, 2.3, 1.4, 1.05, fill="navy" if i == 3 else "light", r=0.06),
                  text(x + 0.05, 2.33, 1.3, 1.0, s, size=9.5, color="white" if i == 3 else "body", valign="middle", align="center")]
        if i < 3:
            items.append(text(x + 1.4, 2.55, 0.15, 0.55, "→", size=13, bold=True, color="teal", align="center", valign="middle"))
    items += [image(6.7, 3.5, 6.1, 2.4, FIG / "chart_curves.png"),
              text(6.7, 5.95, 6.1, 0.8, "Simulated (orange) vs reconstructed published curves (black): within 1–3 points. NOTION's surgical arm predicted with no free parameter: 21.6% vs 20.7% at 10 y.", size=10.5, color="muted", italic=True),
              footnote("Guyot 2012 BMC Med Res Methodol  ·  Kermen 2022 JTCVS Open  ·  Thyregod 2024  ·  Wakami 2022 J Cardiothorac Surg  ·  Zoghbi 2024 JASE")]
    return dict(bg="white", items=items, notes="""For deployment the data already exists: the echo reporting database gives the serial measurements, the STS and TVT registries give the valve and the comorbidities, the EHR fills the gaps through local language-model extraction, and claims linkage catches reinterventions and deaths elsewhere. The one dependency is a linkage key.
For the prototype there is no public dataset with serial prosthetic-valve echoes, so we built one and made it accountable. We reconstructed patient-level data from three published durability curves with the Guyot algorithm, 112 checks pass. We built a simulator where the only thing that degenerates is the valve's true orifice area; gradients follow from flow physics anchored to the ASE normal-value tables, and echo noise is calibrated to the Velders label-flicker statistics.
We calibrated it to the reconstructed curves, and as a test we predicted the NOTION surgical arm with no free parameter: 21.6 percent against the published 20.7. The result is 10,000 valves and 50,000 echoes, with the true deterioration kept aside so we can audit the labels themselves.""")


def slide_validation():
    return dict(bg="white", items=[
        *title("Validation", "No peeking at the future, and death is a competing risk"),
        image(0.5, 1.6, 7.0, 2.75, FIG / "fig_landmark.png"),
        text(0.6, 4.3, 6.8, 0.45, "Every follow-up echo is a landmark: features use only the history up to it; the outcome is dated at the confirming echo.", size=10.5, color="muted", italic=True),
        image(0.5, 4.8, 3.4, 2.15, FIG / "chart_calibration.png"),
        image(4.0, 4.8, 3.4, 2.15, FIG / "chart_samplesize.png"),
        *lead(7.9, 1.8, 4.8, "Design", "Split by valve, 70 / 15 / 15, bootstrap by valve. Death, endocarditis and other reintervention as competing risks. Temporal and external validation in the real study."),
        *lead(7.9, 3.55, 4.8, "Metrics", "Time-dependent AUC and C-index; calibration slope; echoes per patient-decade and delay to confirmed detection. Reported separately on echoes not yet abnormal."),
        *lead(7.9, 5.3, 4.8, "Comparators", "Time since implant alone; the published risk-factor Cox model; for surveillance, the guideline schedule with and without a confirmation echo."),
        footnote("Uno 2011 Stat Med  ·  Austin & Fine 2017 Stat Med  ·  Collins 2024 BMJ (TRIPOD+AI)  ·  Riley 2020 BMJ"),
    ], notes="""Validation is where surveillance models usually cheat, so three rules. Split by valve, never by echo. Date the outcome at the confirming echo, so the model never sees a future echo; the figure shows where a prediction is made and what it may use. And treat death as a competing risk, because forty percent of this cohort dies during follow-up.
We report cause-specific AUC and C-index, but also calibration, because a schedule needs the absolute risk right. And we report everything separately on echoes that are not yet abnormal, because that is where the decision to wait is taken.
The comparators are the honest ones: time since implant, which is what the guideline assumes; the published risk-factor Cox model; and for surveillance, the guideline schedule with and without a confirmation echo. The simulation also gave us the sample size: 5,000 valves bring the AUC interval down to two points.""")


def slide_model():
    return dict(bg="white", items=[
        *title("Model", "A hazard at every echo, turned into a date for the next one"),
        image(0.5, 1.7, 7.0, 3.5, FIG / "fig_model.png"),
        text(0.6, 5.3, 6.8, 1.2, "Discrete-time cause-specific hazard model: survival × hazard summed over 6-month periods gives the cumulative incidence over any horizon, so one model answers \"what is the risk if we wait 6, 12, 24 or 36 months?\"",
             size=13, color="body"),
        image(7.7, 1.75, 5.2, 2.6, FIG / "chart_shap.png"),
        *lead(7.9, 4.55, 5.0, "What the clinician sees", "Risk at 1–5 years, a tier, the top drivers and the next echo interval. Test-set example: \"9.9 y after TAVR: 5-y risk 5.1%, moderate; drivers TAVR\u00a0↓, gradient +4.9\u00a0mmHg\u00a0↑, time since\u00a0implant\u00a0↑.\""),
        footnote("Rizopoulos 2016 Biostatistics (personalised screening intervals)  ·  Lee 2018 AAAI (DeepHit) and Ishwaran 2008 (random survival forests) as next steps at registry scale"),
    ], notes="""The model is a discrete-time cause-specific hazard model. Every echo is a landmark. We cut the time after it into six-month periods and a gradient-boosted classifier gives, for each period, the hazard of confirmed SVD and of a competing event. Summing survival times hazard gives the cumulative incidence over any horizon, so the same model answers the scheduling question directly: what is the risk if we wait six, twelve or thirty-six months?
Features come in three groups: the patient and valve at baseline, the reference echo, and the trajectory since then: change from reference, a smoothed gradient, slope, prior abnormal echoes. The SHAP chart shows what drives the prediction: approach, valve model, the horizon, the change in gradient, time since implant.
The output is a result card: the risk, a tier, the top three drivers, and the next echo interval. Two safety rules sit outside the model: any single abnormal echo triggers a confirmation echo, and symptoms always override the schedule. The model schedules surveillance; it never recommends reintervention.""")


def slide_results():
    return dict(bg="white", items=[
        *title("Results on the synthetic test set", "1,401 held-out valves; every pre-specified bar cleared"),
        image(0.5, 1.6, 6.2, 3.05, FIG / "chart_auc.png"),
        image(6.9, 1.6, 5.9, 3.05, FIG / "chart_policy.png"),
        *lead(0.6, 4.75, 5.9, "Prediction", "AUC 0.875 (0.85–0.90) on echoes not yet abnormal, against 0.65 for time since implant; calibration slope 1.04. At implant, 0.865 vs 0.756 for the published risk-factor Cox."),
        *lead(6.75, 4.75, 6.0, "Surveillance", "The confirmation echo alone cuts the mean delay from 1.5 to 1.1 years. The 1% rule then uses 44% fewer echoes than guideline + confirmation at +0.7 months mean delay, inside the 3-month bound."),
        text(0.6, 6.38, 12.1, 0.6, "Audit against the simulator's truth: even the confirmed VARC-3 label detects 60% of true deterioration, median lag 1.5 years. On 26 real de-identified patients the label rule was 86% sensitive and 89% specific (face validity only).",
             size=12.5, color="body"),
        footnote("All numbers: notebooks/analysis/output/results_summary.md  ·  real-data check: notebooks/real_data_check/  ·  full report: report/main.pdf"),
    ], notes="""On held-out valves the dynamic model reaches an AUC of 0.875 on echoes that are not yet abnormal, against 0.65 for time since implant, which is what the guideline implicitly uses. Calibration slope 1.04, so the absolute risks are usable. At implant time the boosted model beats the published risk-factor Cox at every horizon.
The right chart is the surveillance experiment, six schedules on identical simulated patients. Two findings. The confirmation echo alone cuts the mean delay from 1.5 to 1.1 years. And the risk-adapted schedule with a one percent rule uses 44 percent fewer echoes than guideline-plus-confirmation, at plus 0.7 months mean delay, inside our three-month bound.
Two honest numbers. Because the simulator keeps the truth, we audited the label: even the confirmed VARC-3 label detects only 60 percent of true deterioration, with a median lag of 1.5 years. And on the 26 real de-identified patients we could adjudicate, the label rule behaved as simulated, 86 percent sensitive, 89 percent specific, but that is face validity, not validation.""")


def slide_impact():
    blocks = [(0.9, 1.9, "What it enables", "A surveillance interval per valve, recomputed at every echo, with a confirmation rule built in; and an audit of the SVD definitions against the truth."),
              (7.0, 1.9, "What validation needs", "A linked echo database with STS/TVT at three or more centres, with exact dates and ages; then external validation in a second health system."),
              (0.9, 4.1, "Next three months", "Run the same pipeline on one institution's echo database (the schema is registry-shaped), then evaluate it in silent mode."),
              (7.0, 4.1, "Readiness", "Code and reports are public and reproducible. Real-world performance is unproven until external validation; the system schedules surveillance and never recommends reintervention.")]
    items = [circle(12.9, 0.9, 5.5, "panel"), *title("Impact and next steps", "What changes, and what it will take to prove it", dark=True)]
    for x, y, h, b in blocks:
        items += [text(x, y, 5.5, 0.45, h, size=16, bold=True, color="teal"), text(x, y + 0.5, 5.5, 1.5, b, size=14.5, color="white")]
    items += [text(0.9, 6.3, 11.5, 0.45, f"Team {TEAM}  ·  github.com/abarmper/Dyania-Docathon  ·  report/main.pdf", size=13, color="pale", align="center"),
              footnote("Hofmann 2024 Diagnostics (risk-adapted imaging surveillance precedent)  ·  Genereux 2025 Struct Heart  ·  Collins 2024 BMJ (TRIPOD+AI)", color="pale")]
    return dict(bg="navy", items=items, notes="""What this enables is a surveillance interval per valve, recomputed at every echo, with a confirmation rule built in, and something no registry can do: an audit of the SVD definitions themselves against the truth.
To prove it we need a linked echo database with the STS and TVT registries at a few centres, with exact dates and ages, and then an external system. In the next three months we would run this same pipeline on one institution's echo database, the schema is already registry-shaped, and then evaluate it in silent mode.
Where we stand: the code, the data dictionary and every report are public and reproducible; the real-world performance is unproven until external validation. Thank you, we are happy to take questions.""")


SLIDES = [slide_title, slide_problem, slide_workflow, slide_hypothesis, slide_study, slide_data, slide_validation, slide_model, slide_results, slide_impact]


# ----------------------------------------------------------------------------------------------- pptx renderer
def fit_box(it):
    iw, ih = Image.open(it["path"]).size
    scale = min(it["w"] / iw, it["h"] / ih); w, h = iw * scale, ih * scale
    return it["x"] + (it["w"] - w) / 2, it["y"] + (it["h"] - h) / 2, w, h


def to_pptx(path):
    from pptx import Presentation
    from pptx.dml.color import RGBColor
    from pptx.enum.shapes import MSO_CONNECTOR, MSO_SHAPE
    from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
    from pptx.util import Inches, Pt

    rgb = lambda k: RGBColor.from_string(C[k])
    prs = Presentation(); prs.slide_width = Inches(W); prs.slide_height = Inches(H)
    blank = prs.slide_layouts[6]
    for build in SLIDES:
        sp = build(); s = prs.slides.add_slide(blank)
        s.background.fill.solid(); s.background.fill.fore_color.rgb = rgb(sp["bg"])
        for it in sp["items"]:
            k = it["k"]
            if k in ("rect", "circle"):
                shape = MSO_SHAPE.OVAL if k == "circle" else (MSO_SHAPE.ROUNDED_RECTANGLE if it["r"] else MSO_SHAPE.RECTANGLE)
                sh = s.shapes.add_shape(shape, Inches(it["x"]), Inches(it["y"]), Inches(it["w"]), Inches(it["h"]))
                if k == "rect" and it["r"]:
                    sh.adjustments[0] = min(0.5, it["r"] / min(it["w"], it["h"]))
                sh.shadow.inherit = False
                if it.get("fill"): sh.fill.solid(); sh.fill.fore_color.rgb = rgb(it["fill"])
                else: sh.fill.background()
                if it.get("line"): sh.line.color.rgb = rgb(it["line"]); sh.line.width = Pt(1)
                else: sh.line.fill.background()
            elif k == "text":
                tb = s.shapes.add_textbox(Inches(it["x"]), Inches(it["y"]), Inches(it["w"]), Inches(it["h"]))
                tf = tb.text_frame; tf.word_wrap = True
                tf.margin_left = tf.margin_right = Inches(0.05); tf.margin_top = tf.margin_bottom = Inches(0.03)
                tf.vertical_anchor = {"top": MSO_ANCHOR.TOP, "middle": MSO_ANCHOR.MIDDLE, "bottom": MSO_ANCHOR.BOTTOM}[it["valign"]]
                for i, para in enumerate(it["s"].split("\n")):
                    p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
                    p.alignment = {"left": PP_ALIGN.LEFT, "center": PP_ALIGN.CENTER, "right": PP_ALIGN.RIGHT}[it["align"]]
                    p.line_spacing = it["ls"]; p.space_after = Pt(3)
                    r = p.add_run(); r.text = para
                    f = r.font; f.name = FONT[it["font"]]; f.size = Pt(it["size"]); f.bold = it["bold"]; f.italic = it["italic"]; f.color.rgb = rgb(it["color"])
            elif k == "image":
                x, y, w, h = fit_box(it)
                s.shapes.add_picture(it["path"], Inches(x), Inches(y), Inches(w), Inches(h))
            elif k == "line":
                cn = s.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(it["x1"]), Inches(it["y1"]), Inches(it["x2"]), Inches(it["y2"]))
                cn.line.color.rgb = rgb(it["color"]); cn.line.width = Pt(it["width"])
        s.notes_slide.notes_text_frame.text = sp["notes"]
    prs.save(path)


# ----------------------------------------------------------------------------------------------- pdf renderer
CHAR_W = {"B": 0.50, "T": 0.53}   # average glyph width as a fraction of the font size (Carlito / Caladea)


def _wrap(it):
    cw = CHAR_W[it["font"]] * (1.06 if it["bold"] else 1.0)
    n = max(4, int((it["w"] - 0.1) * 72 / (it["size"] * cw)))
    lines = []
    for para in it["s"].split("\n"):
        lines += textwrap.wrap(para, n, break_long_words=False) or [""]
    return lines


def to_pdf(path):
    _style()
    with PdfPages(path) as pdf:
        for build in SLIDES:
            sp = build()
            fig = plt.figure(figsize=(W, H)); fig.patch.set_facecolor(hx(sp["bg"]))
            ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off"); ax.grid(False)
            for it in sp["items"]:
                k = it["k"]
                if k == "rect":
                    kw = dict(facecolor=hx(it["fill"]) if it["fill"] else "none", edgecolor=hx(it["line"]) if it["line"] else "none", linewidth=1)
                    if it["r"]:
                        ax.add_patch(patches.FancyBboxPatch((it["x"] + it["r"], it["y"] + it["r"]), it["w"] - 2 * it["r"], it["h"] - 2 * it["r"],
                                                            boxstyle=f"round,pad={it['r']},rounding_size={it['r']}", **kw))
                    else:
                        ax.add_patch(patches.Rectangle((it["x"], it["y"]), it["w"], it["h"], **kw))
                elif k == "circle":
                    ax.add_patch(patches.Circle((it["x"] + it["w"] / 2, it["y"] + it["h"] / 2), it["w"] / 2, facecolor=hx(it["fill"]), edgecolor="none"))
                elif k == "line":
                    ax.plot([it["x1"], it["x2"]], [it["y1"], it["y2"]], color=hx(it["color"]), lw=it["width"], solid_capstyle="butt")
                elif k == "image":
                    x, y, w, h = fit_box(it)
                    ax.imshow(np.asarray(Image.open(it["path"]).convert("RGBA")), extent=(x, x + w, y + h, y), aspect="auto", interpolation="lanczos")
                elif k == "text":
                    lines = _wrap(it); lh = it["size"] * it["ls"] * 1.18 / 72; total = lh * len(lines)
                    y0 = {"top": it["y"] + 0.05, "middle": it["y"] + (it["h"] - total) / 2, "bottom": it["y"] + it["h"] - total - 0.05}[it["valign"]]
                    x0 = {"left": it["x"] + 0.05, "center": it["x"] + it["w"] / 2, "right": it["x"] + it["w"] - 0.05}[it["align"]]
                    for i, ln in enumerate(lines):
                        ax.text(x0, y0 + i * lh, ln, fontsize=it["size"], fontfamily=MPL_FONT[it["font"]], fontweight="bold" if it["bold"] else "normal",
                                fontstyle="italic" if it["italic"] else "normal", color=hx(it["color"]), ha=it["align"], va="top")
            pdf.savefig(fig); plt.close(fig)


def word_counts():
    for i, build in enumerate(SLIDES, 1):
        sp = build(); n = sum(len(it["s"].split()) for it in sp["items"] if it["k"] == "text" and it["y"] < 7.0)
        print(f"  slide {i}: {n} words on slide (excluding footnote)")


if __name__ == "__main__":
    make_figures()
    to_pptx(HERE / "slides.pptx")
    to_pdf(HERE / "slides.pdf")
    word_counts()
    print("wrote slides.pptx and slides.pdf")
