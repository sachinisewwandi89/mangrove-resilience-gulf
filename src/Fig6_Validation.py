# Figure 6 - Model validation and covariate effects.
#
# Run from C:\mangrove-resilience   (python Fig6_Validation.py)
#
# INPUT   tables/tableS1_variogram.csv      (empirical variogram)
#         tables/stats_robustness.json      (repeated blocked CV, ALE)
# OUTPUT  figures/fig6_validation.png
#
# a  Empirical variogram of climate-adjusted recovery, with the block size.
# b  Spatially blocked cross-validated R2, with predictor-set ablation.
# c  Accumulated local effects, ordered by effect range.
#
# This replaces the old composite. The previous Figure 6 was assembled from
# figS1_variogram.png plus fig14_statistical_robustness.png, which together
# supplied five panels for a three-panel caption, and one of them repeated
# Figure 4a. Every value below is read from the saved tables, so the figure
# cannot disagree with Section 3.4. Key values are printed at the end.

import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
vg = pd.read_csv(ROOT / "tables" / "tableS1_variogram.csv")
st = json.loads((ROOT / "tables" / "stats_robustness.json").read_text())

# --- house style, shared with Figs 3, 4 and 5 -----------------------------
ORANGE, BLUE, GREEN = "#D55E00", "#0072B2", "#009E73"   # Okabe-Ito, CVD-safe
INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#e6e6e6"
plt.rcParams.update({
    "font.size": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "legend.fontsize": 6.3,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.labelcolor": INK, "savefig.dpi": 400, "savefig.bbox": "tight", "savefig.pad_inches": 0.12,
    "font.family": "DejaVu Sans"})

BLOCK_KM = 27.0          # 0.25 degrees


def tidy(ax, grid="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis=grid, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)


def tag(ax, t, sub):
    ax.text(-0.19, 1.10, t, transform=ax.transAxes, fontsize=9,
            fontweight="bold", ha="left", va="bottom")
    ax.set_title(sub, fontsize=6.8, color=MUTED, loc="left", pad=6)


fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.45),
                         gridspec_kw={"width_ratios": [1, 1, 1.28]})

# ============================================== a. empirical variogram ====
ax = axes[0]
ax.plot(vg.mid_km, vg.pct_of_sill, "-o", color=BLUE, lw=1.6, ms=2.8,
        mec=BLUE, mew=0.0, zorder=3)
ax.axhline(100, color=MUTED, lw=0.8, zorder=2)
ax.axvline(BLOCK_KM, color=ORANGE, lw=1.0, ls="--", zorder=2)
ax.set_xscale("log")

near = float(vg.pct_of_sill.iloc[0])                    # 0-1 km
at_block = float(vg.loc[vg.to_km == 27, "pct_of_sill"].iloc[0])   # 20-27 km

ax.text(0.04, 0.88, f"{near:.0f}% of the variance\nlies within 1 km",
        transform=ax.transAxes, fontsize=6.3, color=INK, ha="left",
        va="center", linespacing=1.3)
ax.text(BLOCK_KM * 1.25, 62, f"block size\n{BLOCK_KM:.0f} km", fontsize=6.3,
        color=ORANGE, ha="left", va="bottom", linespacing=1.25)
ax.set_xlabel("Separation distance (km, log scale)")
ax.set_ylabel("Semivariance (% of total variance)")
ax.set_ylim(55, 130)
tidy(ax)
tag(ax, "a", "Dependence extends past the block")

# ==================================== b. blocked CV and ablation ==========
ax = axes[1]
sets = st["cv"]["sets"]
order = sorted(sets, key=lambda k: sets[k]["r2_mean"], reverse=True)
LABEL = {"all predictors": "all predictors",
         "projectable set": "projectable set",
         "without rmse": "without fit error",
         "without climate": "without climate",
         "intrinsic only": "stand covariates",
         "without pre_nbr": "without pre-NBR",
         "environment only": "surroundings only"}
COLOR = {"all predictors": BLUE, "projectable set": GREEN}

full = sets["all predictors"]["r2_mean"]
ax.axvline(full, color=BLUE, lw=0.8, ls=":", zorder=2)
ax.axvline(0, color=MUTED, lw=0.8, zorder=2)

for i, k in enumerate(order):
    y = len(order) - 1 - i
    s = sets[k]
    col = COLOR.get(k, ORANGE if s["delta"] < -0.05 else MUTED)
    ax.plot([s["fold_min"], s["fold_max"]], [y, y], color=col, lw=0.9,
            alpha=0.35, zorder=3, solid_capstyle="butt")
    ax.plot([s["r2_mean"] - s["r2_sd"], s["r2_mean"] + s["r2_sd"]], [y, y],
            color=col, lw=2.6, zorder=4, solid_capstyle="butt")
    ax.plot([s["r2_mean"]], [y], "o", ms=4.2, color=col, mec="white",
            mew=0.7, zorder=5)
    ax.text(0.86, y, f"{s['r2_mean']:.2f}", transform=ax.get_yaxis_transform(),
            fontsize=6.5, color=INK, va="center", ha="right")

ax.set_yticks(range(len(order)))
ax.set_yticklabels([LABEL[k] for k in reversed(order)], fontsize=6.4)
ax.set_ylim(-0.7, len(order) - 0.3 + 0.55)
ax.set_xlim(-0.15, 0.90)
ax.set_xlabel("Spatially blocked R$^2$")
ax.text(0.98, 0.985, f"random folds R$^2$ = {st['cv']['r2_random']:.2f}",
        transform=ax.transAxes, fontsize=6.3, color=MUTED, ha="right",
        va="top")
tidy(ax, grid="x")
tag(ax, "b", "Pre-disturbance state carries the model")

# ============================================ c. accumulated local effects =
ax = axes[2]
NAME = {"pre_nbr": "Pre-disturbance NBR",
        "climate_shift": "Climate shift",
        "dist_to_channel_m": "Distance to channel",
        "magnitude": "Magnitude of loss",
        "rain_recovery_mm": "Rainfall 2016-24",
        "rain_2015_mm": "Rainfall 2015"}
STYLE = {"pre_nbr": (BLUE, "-", 2.0),
         "climate_shift": (GREEN, "-", 1.3),
         "dist_to_channel_m": (ORANGE, "-", 1.3),
         "magnitude": (MUTED, "-", 1.1),
         "rain_recovery_mm": (GREEN, "--", 1.1),
         "rain_2015_mm": (ORANGE, "--", 1.1)}
ale = st["ale"]
keys = sorted(NAME, key=lambda k: ale[k]["range"], reverse=True)
ax.axhline(0, color=MUTED, lw=0.8, zorder=2)
for k in keys:
    x = np.asarray(ale[k]["x"], float)
    y = np.asarray(ale[k]["ale"], float)
    xs = (x - x.min()) / (x.max() - x.min())
    col, ls, lw = STYLE[k]
    ax.plot(xs, y, color=col, ls=ls, lw=lw, zorder=3,
            label=f"{NAME[k]} ({ale[k]['range']:.3f})")
ax.set_xlabel("Covariate, scaled to observed range")
ax.set_ylabel("Accumulated local effect (NBR)")
ax.set_xlim(0, 1)
ax.set_ylim(-0.135, 0.145)
ax.legend(frameon=False, loc="upper left", handlelength=1.5,
          handletextpad=0.5, labelspacing=0.26, borderaxespad=0.15,
          fontsize=5.9)
tidy(ax)
tag(ax, "c", "Effect range in brackets")

fig.subplots_adjust(wspace=0.72)
out = ROOT / "figures" / "fig6_validation.png"
try:
    fig.savefig(out)
except OSError:
    import time
    alt = out.with_name(f"fig6_validation_{time.strftime('%H%M%S')}.png")
    fig.savefig(alt)
    print(f"\n!! could not overwrite {out.name} - close it in any viewer")
    print(f"!! saved as {alt.name} instead\n")
    out = alt
print("saved", out)

print("")
print("CHECK THESE AGAINST THE MANUSCRIPT (Section 3.4)")
print(f"  semivariance at 0-1 km      {near:.1f}%     draft says 61%")
print(f"  semivariance at 20-27 km    {at_block:.1f}%     draft says 84%")
print(f"  largest lag reached         {vg.pct_of_sill.max():.0f}%      (so it does not plateau)")
print(f"  blocked R2, all predictors  {full:.3f}    text: 0.365")
print(f"  blocked R2, projectable     {sets['projectable set']['r2_mean']:.3f}    text: 0.349")
print(f"  random-fold R2              {st['cv']['r2_random']:.3f}    text: 0.759")
print(f"  ALE range, pre-NBR          {ale['pre_nbr']['range']:.3f}    text: 0.197")
