"""
Figure 3 - the detection floor, and the two reference-sample confirmations.

"""
 
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.stats import mannwhitneyu, fisher_exact

ROOT = Path(__file__).resolve().parent
REF = pd.read_csv(ROOT / "tables" / "reference_point_geometry.csv")
PIX = pd.read_csv(ROOT / "data" / "pixels_with_covariates.csv")
 
# house style, shared with the other figures 
BLUE, ORANGE, GREEN = "#0072B2", "#D55E00", "#009E73"      # Okabe-Ito, CVD-safe
INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#e6e6e6"
plt.rcParams.update({
    "font.size": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "legend.fontsize": 6.5,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.labelcolor": INK, "savefig.dpi": 400, "savefig.bbox": "tight",
    "font.family": "DejaVu Sans"})
 
NBR_CAN, NBR_BG, THR = 0.80, 0.00, 0.15
CM = THR / (NBR_CAN - NBR_BG)          # 0.1875, invariant to cell size
 
 
def tidy(ax, grid="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis=grid, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
 
 
def tag(ax, t, sub):
    ax.text(-0.17, 1.13, t, transform=ax.transAxes, fontsize=9,
            fontweight="bold", ha="left", va="bottom")
    ax.set_title(sub, fontsize=6.8, color=MUTED, loc="left", pad=7)
 
 
def ecdf(v):
    v = np.sort(np.asarray(v, float))
    return v, np.arange(1, len(v) + 1) / len(v)
 
 
det = REF[(REF.map_class == 2) & (REF.ref_class == 2)]
mis = REF[(REF.map_class != 2) & (REF.ref_class == 2)]
GRP = (("detected by the map", det, BLUE, "-"),
       ("missed by the map", mis, ORANGE, "--"))

# the tests, computed here rather than written in by hand 
# The mixture model predicts the DIRECTION of both differences - missed
# dieback should stand in sparser canopy and lie nearer a stand margin - so
# both are one-sided Mann-Whitney tests, matching the methods section. The
# panel subtitles are formatted from these values
P_CANOPY = mannwhitneyu(mis.canopy, det.canopy, alternative="less")[1]
P_EDGE = mannwhitneyu(mis.edge_m, det.edge_m, alternative="less")[1]
P_MARGIN = fisher_exact([[(mis.edge_m == 0).sum(), (mis.edge_m > 0).sum()],
                         [(det.edge_m == 0).sum(), (det.edge_m > 0).sum()]])[1]


SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


def pfmt(p):
    """P = 0.007 down to 0.001, then one significant figure x 10^-n."""
    if p >= 1e-3:
        return ("P = %.3f" % p).rstrip("0").rstrip(".")
    e = int(np.floor(np.log10(p)))
    return "P = %.0f × 10⁻%s" % (p / 10 ** e, str(-e).translate(SUP))

fig, axes = plt.subplots(2, 2, figsize=(7.09, 4.75))
 
# a. the floor 
ax = axes[0, 0]
w = np.linspace(0, 12, 400)
for res, col, ls in ((30, ORANGE, "-"), (20, BLUE, "--"), (10, GREEN, ":")):
    ax.plot(w, (w / res) * (NBR_CAN - NBR_BG), color=col, lw=1.8, ls=ls,
            label=f"{res} m cell", zorder=3)
    wc = CM * res
    if wc <= w.max():
        ax.plot([wc], [THR], "o", ms=5, color=col, mec="white", mew=0.8, zorder=4)
        dx, dy, ha = {30: (0.5, 0.105, "left"),
                      20: (0.5, 0.060, "left"),
                      10: (-0.5, 0.030, "right")}[res]
        ax.annotate(f"{wc:.1f} m", (wc, THR), (wc + dx, THR + dy),
                    fontsize=6.5, ha=ha, color=col,
                    arrowprops=dict(arrowstyle="-", color=col, lw=0.5))
ax.axhline(THR, color=INK, lw=1.0)
ax.text(11.8, THR - 0.024, "detection threshold 0.15", fontsize=6.5,
        ha="right", color=INK)
ax.set_xlabel("Width of fully killed strip (m)")
ax.set_ylabel("Resulting change in NBR")
ax.set_xlim(0, 12); ax.set_ylim(0, 0.42)
ax.legend(frameon=False, loc="upper left")
tidy(ax); tag(ax, "a", "Narrower loss cannot register, however complete")
 
# b. canopy fraction ECDF 
ax = axes[0, 1]
for nm, g, col, ls in GRP:
    x, y = ecdf(g.canopy)
    ax.step(x, y, where="post", color=col, lw=1.8, ls=ls,
            label=f"{nm}  (n={len(g)})", zorder=3)
for nm, g, col, off in (("detected", det, BLUE, (0.055, 0.30)),
                        ("missed", mis, ORANGE, (-0.055, 0.70))):
    m = float(np.median(g.canopy))
    ax.plot([m], [0.5], "o", ms=5, color=col, mec="white", mew=0.8, zorder=5)
    ax.annotate(f"median {m:.2f}\nneeds {CM/m:.0%} mortality",
                (m, 0.5), (m + off[0], off[1]), fontsize=6.5, color=col,
                ha="left" if off[0] > 0 else "right",
                va="center", linespacing=1.35,
                arrowprops=dict(arrowstyle="-", color=col, lw=0.5))
ax.set_xlabel("Pre-disturbance canopy fraction")
ax.set_ylabel("Cumulative share of points")
ax.set_xlim(0.30, 1.05); ax.set_ylim(0, 1.02)
ax.legend(frameon=False, loc="upper left", borderpad=0.1)
tidy(ax); tag(ax, "b", f"Missed dieback stood in sparser canopy  ({pfmt(P_CANOPY)})")
 
# c. distance to boundary 
ax = axes[1, 0]
XMAX = float(np.ceil(np.percentile(pd.concat([det.edge_m, mis.edge_m]), 95) / 50) * 50)
for nm, g, col, ls in GRP:
    x, y = ecdf(np.clip(g.edge_m, 0, XMAX))
    x = np.r_[0, x]; y = np.r_[0, y]
    ax.step(x, y, where="post", color=col, lw=1.8, ls=ls,
            label=f"{nm}  (n={len(g)})", zorder=3)
    ax.plot([0], [(g.edge_m == 0).mean()], "o", ms=5, color=col,
            mec="white", mew=0.8, zorder=5)

ax.annotate(f"missed by the map  (n = {len(mis)})\n"
            f"median {mis.edge_m.median():.0f} m, "
            f"{(mis.edge_m == 0).mean():.0%} on a margin",
            xy=(0, (mis.edge_m == 0).mean()), xycoords="data",
            xytext=(0.30, 0.46), textcoords="axes fraction",
            fontsize=6.5, color=ORANGE, va="center", ha="left", linespacing=1.35,
            arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.5,
                            shrinkA=0, shrinkB=3))
ax.annotate(f"detected by the map  (n = {len(det)})\n"
            f"median {det.edge_m.median():.0f} m, "
            f"{(det.edge_m == 0).mean():.0%} on a margin",
            xy=(0, (det.edge_m == 0).mean()), xycoords="data",
            xytext=(0.30, 0.17), textcoords="axes fraction",
            fontsize=6.5, color=BLUE, va="center", ha="left", linespacing=1.35,
            arrowprops=dict(arrowstyle="-", color=BLUE, lw=0.5,
                            shrinkA=0, shrinkB=3))
ax.set_xlabel(f"Distance to nearest stand boundary (m), truncated at {XMAX:.0f}")
ax.set_ylabel("Cumulative share of points")
ax.set_xlim(-XMAX * 0.02, XMAX); ax.set_ylim(0, 1.02)
tidy(ax); tag(ax, "c", f"Missed dieback lay nearer the margin  ({pfmt(P_EDGE)})")
 
# d. the noise floor
ax = axes[1, 1]
ax.hist(PIX.rmse, bins=np.linspace(0, 0.22, 60), color=GREEN, alpha=0.85, lw=0)
ax.axvline(THR, color=ORANGE, lw=1.6)
ax.axvline(PIX.rmse.median(), color=INK, lw=1.0, ls=":")
ax.text(THR - 0.004, ax.get_ylim()[1] * 0.95, "detection\nthreshold 0.15",
        fontsize=6.5, ha="right", va="top", color=ORANGE)
ax.text(PIX.rmse.median() + 0.004, ax.get_ylim()[1] * 0.62,
        f"median fit error\n{PIX.rmse.median():.3f}"
        f"\n({PIX.rmse.median()/THR:.0%} of threshold)",
        fontsize=6.5, color=INK, va="top")
ax.set_xlabel("Segmentation fit error (NBR)")
ax.set_ylabel("Pixels")
ax.set_xlim(0, 0.22)
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
tidy(ax); tag(ax, "d", "Residual noise left by harmonisation")
 
fig.subplots_adjust(wspace=0.34, hspace=0.62)
out = ROOT / "figures" / "fig3_detection_floor.png"
fig.savefig(out)
print("saved", out)

print(f"  panel b  canopy fraction, detected {det.canopy.median():.2f} (n={len(det)})"
      f"  |  missed {mis.canopy.median():.2f} (n={len(mis)})  |  {pfmt(P_CANOPY)}")
print(f"           required mortality {CM/det.canopy.median():.0%} against"
      f" {CM/mis.canopy.median():.0%}")
print(f"  panel c  distance to boundary, detected {det.edge_m.median():.0f} m"
      f"  |  missed {mis.edge_m.median():.0f} m  |  {pfmt(P_EDGE)}")
print(f"           on a margin, detected {(det.edge_m==0).mean():.0%}"
      f"  |  missed {(mis.edge_m==0).mean():.0%}"
      f"  |  Fisher exact {pfmt(P_MARGIN)}")
print(f"panel d  median fit error {PIX.rmse.median():.3f}"
      f"  = {PIX.rmse.median()/THR:.0%} of threshold")
print(f"  median census canopy fraction = {(PIX.pre_nbr/NBR_CAN).median():.2f}")
 