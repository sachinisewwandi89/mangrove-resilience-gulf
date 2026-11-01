
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent
import pandas as pd

C_DIEBACK = "#D55E00"
C_REF     = "#0072B2"
C_INK     = "#0E3B43"
C_MUTED   = "#5A6B6E"
C_GRID    = "#E7ECEE"
W2 = 7.20

# font 
from matplotlib import font_manager as _fm


def pick_font(*names):
    have = {f.name for f in _fm.fontManager.ttflist}
    for n in names:
        if n in have:
            return n
    return "DejaVu Sans"


FONT = pick_font("Segoe UI", "Lato", "DejaVu Sans")

plt.rcParams.update({
    "font.family": FONT, "savefig.transparent": True,
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.edgecolor": "#C3D0D3", "axes.linewidth": 0.8,
    "xtick.color": C_MUTED, "ytick.color": C_MUTED,
    "text.color": C_INK, "axes.labelcolor": C_INK,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
})

# titles


traj = pd.read_csv(ROOT / "data" / "phase3_trajectories.csv")

fig, ax = plt.subplots(figsize=(W2, 3.4))
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.grid(True, color=C_GRID, lw=0.5, zorder=0)
ax.set_axisbelow(True)

# Decline window, 2014 to 2016 — the corrected timing.
ax.axvspan(2014, 2016, color="#f0f0f0", zorder=0)
ax.annotate("decline 2014–2016", xy=(2015, 0.985),
            xycoords=("data", "axes fraction"),
            ha="center", va="top", fontsize=7, color=C_MUTED)

ax.plot(traj.year, traj.fit_reference, color=C_REF, lw=2.0, ls=(0, (5, 2)),
        label="Undisturbed reference", zorder=4)
ax.plot(traj.year, traj.fit_dieback, color=C_DIEBACK, lw=2.0,
        label="2015–16 dieback", zorder=3)

# Markers on the three decline vertices so the two-stage shape is visible
# rather than inferred: 2014 -> 2015 carries 59% of the loss.
d = traj.set_index("year")["fit_dieback"]
ax.plot([2014, 2015, 2016], [d[2014], d[2015], d[2016]], "o",
        ms=4.5, color=C_DIEBACK, mec="white", mew=0.9, zorder=6)

for yr, dx, dy in ((2014, -22, 6), (2015, 7, -3), (2016, 9, -10)):
    ax.annotate(f"{d[yr]:.2f}", xy=(yr, d[yr]), xytext=(dx, dy),
                textcoords="offset points", fontsize=6.5, color=C_MUTED)

# The two gaps that carry the recovery argument.
for yr, val, lab in ((2014, 0.112, "+0.112"), (2025, 0.098, "+0.098")):
    lo = traj.set_index("year")["fit_reference"][yr]
    hi = traj.set_index("year")["fit_dieback"][yr]
    ax.annotate("", xy=(yr, hi), xytext=(yr, lo),
                arrowprops=dict(arrowstyle="<->", color=C_MUTED, lw=0.7))
    ax.annotate(lab, xy=(yr, (hi + lo) / 2), xytext=(4, 0),
                textcoords="offset points", fontsize=6.5, color=C_MUTED,
                va="center")

xe = traj.year.max()
ax.annotate("Dieback", xy=(xe, traj.fit_dieback.iloc[-1]),
            xytext=(4, 0), textcoords="offset points",
            fontsize=7, color=C_INK, va="center")
ax.annotate("Undisturbed", xy=(xe, traj.fit_reference.iloc[-1]),
            xytext=(4, 0), textcoords="offset points",
            fontsize=7, color=C_INK, va="center")

ax.set_xlabel("Year")
ax.set_ylabel("Leaf density of the canopy")
ax.set_xlim(1987, 2026)
ax.set_ylim(0.26, 0.75)
ax.set_xticks([1990, 1995, 2000, 2005, 2010, 2015, 2020, 2025])
ax.legend(frameon=False, loc="lower left")

# Title block 
if TITLES:
    ax.set_title(
        "Canopy trajectories, 1987–2025: dieback stands versus undisturbed reference",
        fontsize=9.5, fontweight="bold", color=C_INK, loc="left", pad=22)
ax.annotate(
    "Gulf of Carpentaria  ·  how dense the leafy canopy was, each year  ·  "
    "29,453 pixels that died, 765,138 that did not",
    xy=(0, 1.015), xycoords="axes fraction", ha="left", va="bottom",
    fontsize=7.5, color=C_MUTED)

fig.savefig(ROOT / "figures" / "fig1_trajectories.png")
print("ok")
