"""
Figure 4 - the magnitude, type and direction of disagreement.
 
"""
 
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
 
ROOT = Path(__file__).resolve().parent
ACC = json.load(open(ROOT / "tables" / "accuracy_summary.json"))
 
BLUE, ORANGE, GREEN = "#0072B2", "#D55E00", "#009E73"      # Okabe-Ito, CVD-safe
INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#e6e6e6"
plt.rcParams.update({
    "font.size": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "legend.fontsize": 6.5,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.labelcolor": INK, "savefig.dpi": 400, "savefig.bbox": "tight",
    "font.family": "DejaVu Sans"})
 
 
def tidy(ax, grid="x"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis=grid, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
 
 
def tag(ax, t, sub):
    ax.text(-0.02, 1.14, t, transform=ax.transAxes, fontsize=9,
            fontweight="bold", ha="right", va="bottom")
    ax.set_title(sub, fontsize=6.8, color=MUTED, loc="left", pad=7)
 
 
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.09, 2.85),
                               gridspec_kw={"width_ratios": [1.0, 1.15]})
 
# a. components of disagreement 
FRAMES = (("Sample-based", ACC["disagreement_sample"]),
          ("Area-weighted", ACC["disagreement"]))
PARTS = (("Quantity", "quantity", BLUE),
         ("Exchange", "exchange", ORANGE),
         ("Shift", "shift", GREEN))
 
for i, (fname, dd) in enumerate(FRAMES):
    left = 0.0
    for pname, key, col in PARTS:
        v = float(dd[key]) * 100
        axA.barh(i, v, left=left, height=0.52, color=col, lw=0, zorder=3)
        # 2px surface gap between stacked segments
        axA.barh(i, 0.18, left=left + v - 0.09, height=0.52,
                 color="white", lw=0, zorder=4)
        if v > 2.5:                       # direct label, so colour is not alone
            axA.text(left + v / 2, i, f"{v:.2f}", ha="center", va="center",
                     fontsize=6.5, color="white", fontweight="bold", zorder=5)
        left += v
    axA.text(left + 0.7, i, f"total {float(dd['total'])*100:.2f}", va="center",
             fontsize=6.5, color=INK)
    axA.annotate(f"shift {float(dd['shift'])*100:.2f}",
                 xy=(left - float(dd['shift']) * 100 / 2, i - 0.30),
                 xytext=(left - 6.0, i - 0.46), fontsize=6.3, color=GREEN,
                 ha="center", va="bottom",
                 arrowprops=dict(arrowstyle="-", color=GREEN, lw=0.5,
                                 shrinkA=0, shrinkB=1))
 
axA.set_yticks(range(len(FRAMES)))
axA.set_yticklabels([f for f, _ in FRAMES])
axA.set_ylim(-0.6, len(FRAMES) - 0.1)
axA.set_xlim(0, 40)
axA.set_xlabel("Disagreement (percentage points of the map)")
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in PARTS]
axA.legend(handles, [p for p, _, _ in PARTS], frameon=False, ncol=3,
           loc="lower left", bbox_to_anchor=(0.0, -0.02), handlelength=1.1,
           columnspacing=1.1, handletextpad=0.5)
tidy(axA); tag(axA, "a", "Misplacement dominates, and it is mutual")
 
# b. direction of confusion
BT = ACC["bradley_terry"]
STD = BT.get("standardised")
if STD is None:
    raise SystemExit("accuracy_summary.json has no 'standardised' Bradley-Terry ")
CLASSES = ["Undisturbed", "2015-16 dieback", "Other-year change"]
K = len(CLASSES)
 
SERIES = (("standardised (result)", STD, INK, "o", 5.2, True),
          ("unstandardised, sample", BT.get("unstandardised_sample"), BLUE, "s", 4.4, False),
          ("unstandardised, area", BT.get("unstandardised_area"), ORANGE, "^", 4.8, False))
OFFS = (0.0, 0.23, -0.23)
 
axB.axvline(0, color=MUTED, lw=0.8, zorder=2)
for (nm, res, col, mk, ms, primary), off in zip(SERIES, OFFS):
    if res is None:
        continue
    sc = np.asarray(res["score"], float)
    se = np.asarray(res["se"], float)
    ys = np.arange(K) + off
    anchored = se == 0                     # the reference class, fixed by design
    if primary:
        axB.errorbar(sc[~anchored], ys[~anchored],
                     xerr=1.96 * se[~anchored], fmt="none",
                     ecolor=col, elinewidth=1.4, capsize=2.2, zorder=4)
        axB.plot(sc[~anchored], ys[~anchored], mk, ms=ms, color=col,
                 mec="white", mew=0.8, zorder=5, label=nm)
        axB.plot(sc[anchored], ys[anchored], mk, ms=ms, mfc="white",
                 mec=col, mew=1.2, zorder=5)
    else:
        axB.plot(sc, ys, mk, ms=ms, mfc="white", mec=col, mew=1.2,
                 zorder=4, label=nm, ls="none")
 
axB.set_yticks(range(K))
axB.set_yticklabels(CLASSES)
axB.set_ylim(-0.75, K + 0.45)
axB.invert_yaxis()
axB.set_xlabel("Bradley–Terry strength (log odds)")
axB.set_xlim(-3.9, 2.3)
 
# short note beside the anchored class, no leader line needed
axB.text(0.22, -0.44, "reference class, fixed at 0", fontsize=6.3,
         color=MUTED, ha="left", va="center")
axB.annotate("the two unstandardised fits\ndisagree about the sign",
             xy=(float(np.asarray(BT["unstandardised_area"]["score"])[1]), 1 - 0.23),
             xytext=(-3.75, 0.42), fontsize=6.3, color=ORANGE,
             va="center", ha="left", linespacing=1.3,
             arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.5,
                             shrinkA=0, shrinkB=4))
axB.legend(frameon=False, loc="lower left", bbox_to_anchor=(-0.01, -0.02),
           handletextpad=0.4, labelspacing=0.3)
tidy(axB); tag(axB, "b", "Standardising removes the design, and the ambiguity")
 
fig.subplots_adjust(wspace=0.42)
out = ROOT / "figures" / "fig4_disagreement.png"
try:
    fig.savefig(out)
except OSError:
    import time
    alt = out.with_name(f"fig4_disagreement_{time.strftime('%H%M%S')}.png")
    fig.savefig(alt)
    print(f"\n!! could not overwrite {out.name} ")
    print(f"!! saved as {alt.name} instead\n")
    out = alt
print("saved", out)
 
print("\npanel a")
for fname, dd in FRAMES:
    print(f"  {fname:<15} total {dd['total']*100:5.2f}"
          f"  quantity {dd['quantity']*100:5.2f}"
          f"  allocation {dd['allocation']*100:5.2f}"
          f"  exchange {dd['exchange']*100:5.2f}"
          f"  shift {dd['shift']*100:4.2f}"
          f"   exchange = {dd['exchange']/dd['allocation']:.1%} of allocation")
print("\npanel b")
for nm, res, *_ in SERIES:
    if res is None:
        continue
    print(f"  {nm:<24}" + "  ".join(
        f"{c} {v:+.3f}" for c, v in zip(CLASSES, res["score"])))