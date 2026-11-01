# Figure 7 - Recovery capacity, observed and predicted.


from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter

ROOT = Path(__file__).resolve().parent
seg = pd.read_csv(ROOT / "tables" / "table5_recovery_capacity.csv")
cell = pd.read_csv(ROOT / "tables" / "table7_predicted_capacity.csv")
vul = pd.read_csv(ROOT / "tables" / "table8b_vulnerability_with_aoa.csv")

SUBSET = "inside AOA"    # or "all pixels" for the unrestricted totals
SPLIT = 138.5            # western / southeastern division

# house style, shared with Figs 3 to 6 
ORANGE, BLUE, GREEN = "#D55E00", "#0072B2", "#009E73"   # Okabe-Ito, CVD-safe
INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#e6e6e6"
plt.rcParams.update({
    "font.size": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "legend.fontsize": 6.3,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.labelcolor": INK, "savefig.dpi": 400, "savefig.bbox": "tight", "savefig.pad_inches": 0.12,
    "font.family": "DejaVu Sans"})

# blue (predicted above baseline) - grey - orange (predicted below baseline)
CMAP = LinearSegmentedColormap.from_list(
    "deficit", [BLUE, "#9ecae8", "#efefef", "#f0b27a", ORANGE])


def tidy(ax, grid="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis=grid, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)


def tag(ax, t, sub, dx=-0.16):
    ax.text(dx, 1.10, t, transform=ax.transAxes, fontsize=9,
            fontweight="bold", ha="left", va="bottom")
    ax.set_title(sub, fontsize=6.8, color=MUTED, loc="left", pad=6)


fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.05),
                         gridspec_kw={"width_ratios": [1, 1.42]})

# a. observed, ranked segments 
ax = axes[0]
s = seg.sort_values("frac_unrecovered", ascending=False).reset_index(drop=True)
x = np.arange(1, len(s) + 1)
west = (s.region == "west").values

for sel, col, lab in ((west, ORANGE, "Western Gulf"),
                      (~west, BLUE, "Southeastern Gulf")):
    ax.vlines(x[sel], s.unrec_lo[sel], s.unrec_hi[sel], color=col, lw=0.8,
              alpha=0.30, zorder=2)
    ax.plot(x[sel], s.frac_unrecovered[sel], "o", ms=2.6, color=col,
            mew=0, zorder=3, label=f"{lab}  (n = {sel.sum()})")

for sel, col, va in ((west, ORANGE, "bottom"), (~west, BLUE, "top")):
    m = float(s.frac_unrecovered[sel].median())
    ax.axhline(m, color=col, lw=0.9, ls="--", zorder=4)
    ax.text(2, m, f"median {m:.2f}", fontsize=6.3, color=col,
            ha="left", va=va)

ax.set_xlabel("Coastal segment, ranked")
ax.set_ylabel("Fraction of dieback area unrecovered")
ax.set_xlim(0, len(s) + 1)
ax.set_ylim(-0.03, 1.05)
ax.legend(frameon=False, loc="upper right", handlelength=0.9,
          handletextpad=0.5, labelspacing=0.3, markerscale=1.6)
tidy(ax)
tag(ax, "a", "Observed: the west fills the bad end")

# b. predicted deficit across the Gulf 
ax = axes[1]
norm = TwoSlopeNorm(vmin=-0.15, vcenter=0.0, vmax=0.15)
sc = ax.scatter(cell.lon, cell.lat, c=cell.pred_deficit.clip(-0.15, 0.15),
                cmap=CMAP, norm=norm, s=5.5, marker="s", lw=0, zorder=3)
ax.axvline(SPLIT, color=MUTED, lw=0.8, ls="--", zorder=2)
ax.set_xlim(135.0, 141.5)
ax.set_ylim(-18.0, -14.0)
ax.set_aspect(1.0 / np.cos(np.radians(16.0)))
ax.set_xlabel("Longitude (\u00b0E)")
ax.set_ylabel("Latitude (\u00b0S)")
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v):.1f}"))
ax.grid(color=GRID, lw=0.5, zorder=0)
ax.set_axisbelow(True)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

v = vul[vul.subset == SUBSET].set_index("region")
ANN = [("Western Gulf", ORANGE, "west"), ("Southeastern Gulf", BLUE, "east")]
yy = 0.97
for name, col, reg in ANN:
    ax.text(0.38, yy, name, transform=ax.transAxes, fontsize=6.6, color=col,
            fontweight="bold", ha="left", va="top")
    ax.text(0.38, yy - 0.055,
            f"{v.loc[reg,'vulnerable_pct']:.1f}% vulnerable, "
            f"{v.loc[reg,'vulnerable_ha']:,.0f} ha",
            transform=ax.transAxes, fontsize=6.3, color=INK, ha="left",
            va="top")
    yy -= 0.135
tag(ax, "b", "Predicted: deficit after a similar event", dx=-0.13)

cax = ax.inset_axes([0.05, 0.225, 0.36, 0.030])
cb = fig.colorbar(sc, cax=cax, orientation="horizontal",
                  ticks=[-0.15, 0.0, 0.15])
cb.set_label("Predicted recovery deficit (NBR)", fontsize=6.3, labelpad=2.5)
cb.ax.set_xticklabels(["\u2264 \u22120.15", "0", "\u2265 0.15"], fontsize=6.0)
cb.outline.set_visible(False)
cb.ax.tick_params(length=2, width=0.5, pad=1.5)

fig.subplots_adjust(wspace=0.34)
out = ROOT / "figures" / "fig7_capacity.png"
try:
    fig.savefig(out)
except OSError:
    import time
    alt = out.with_name(f"fig7_capacity_{time.strftime('%H%M%S')}.png")
    fig.savefig(alt)
    print(f"\n!! could not overwrite {out.name} - close it in any viewer")
    print(f"!! saved as {alt.name} instead\n")
    out = alt
print("saved", out)

print("")
print(f"  west vulnerable    {v.loc['west','vulnerable_pct']:.1f}%  "
      f"{v.loc['west','vulnerable_ha']:,.0f} ha     text: 46.6%, 13,511 ha")
print(f"  southeast          {v.loc['east','vulnerable_pct']:.1f}%   "
      f"{v.loc['east','vulnerable_ha']:,.0f} ha      text: 7.3%, 1,728 ha")
print(f"  ratio west:east    {v.loc['west','vulnerable_pct']/v.loc['east','vulnerable_pct']:.1f}x")
print(f"  observed median unrecovered, west {seg[seg.region=='west'].frac_unrecovered.median():.3f}")
print(f"  observed median unrecovered, east {seg[seg.region=='east'].frac_unrecovered.median():.3f}")
print(f"  segments: {len(seg)}   display cells: {len(cell)}")
