# Figure 5 - Recovery a decade after the event.


from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

ROOT = Path(__file__).resolve().parent
d = pd.read_csv(ROOT / "data" / "pixels_with_covariates.csv")

# house style, shared with Figs 3 and 4 
ORANGE, BLUE, GREEN = "#D55E00", "#0072B2", "#009E73"   # Okabe-Ito, CVD-safe
INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#e6e6e6"
plt.rcParams.update({
    "font.size": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "legend.fontsize": 6.5,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.labelcolor": INK, "savefig.dpi": 400, "savefig.bbox": "tight", "savefig.pad_inches": 0.12,
    "font.family": "DejaVu Sans"})

THR = 0.80          # recovery threshold, 80% of pre-disturbance NBR
SPLIT = 138.5       # western / southeastern division


def tidy(ax, grid="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis=grid, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)


def tag(ax, t, sub):
    ax.text(-0.16, 1.10, t, transform=ax.transAxes, fontsize=9,
            fontweight="bold", ha="left", va="bottom")
    ax.set_title(sub, fontsize=6.8, color=MUTED, loc="left", pad=6)


# the three outcome groups 
never = d.unrecovered.astype(bool).values                 # never reached 0.80
declined = (~never) & (d.recovery_ratio < THR).values     # reached, then fell
holding = (~never) & (d.recovery_ratio >= THR).values     # reached, still above
f_never, f_decl, f_hold = never.mean(), declined.mean(), holding.mean()

fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.45))

#  a. Kaplan-Meier 
ax = axes[0]
# Recovery is observed at years_to_recovery; pixels that never recovered are
# right-censored at the end of the record.
t_obs = d.years_to_recovery.values
tmax = int(np.nanmax(t_obs))
dur = np.where(np.isnan(t_obs), tmax, t_obs)
event = ~np.isnan(t_obs)

times = np.arange(0, tmax + 1)
S, s = [], 1.0
for t in times:
    at_risk = (dur >= t).sum()
    failed = ((dur == t) & event).sum()
    if at_risk > 0 and t > 0:
        s *= 1 - failed / at_risk
    S.append(s)
S = np.array(S)

ax.step(times, S, where="post", color=ORANGE, lw=2.0, zorder=3)
ax.axhline(f_never, color=INK, lw=1.0, ls="--", zorder=2)
ax.text(0.3, f_never + 0.04, f"{f_never:.1%} never recover",
        fontsize=6.5, ha="left", color=INK)
med = times[np.argmax(S <= 0.5)] if (S <= 0.5).any() else np.nan
if not np.isnan(med):
    ax.plot([med, med], [0, 0.5], color=MUTED, lw=0.7, ls=":", zorder=2)
    ax.plot([0, med], [0.5, 0.5], color=MUTED, lw=0.7, ls=":", zorder=2)
    ax.text(med - 0.25, 0.545, f"median {med:.0f} yr", fontsize=6.5,
            color=MUTED, ha="right")
ax.set_xlabel("Years since disturbance")
ax.set_ylabel("Proportion not yet recovered")
ax.set_xlim(0, tmax); ax.set_ylim(0, 1.02)
tidy(ax); tag(ax, "a", "Recovery is slow, and incomplete")

# b. recovery ratio, 3 groups 
ax = axes[1]
# The axis is truncated at -0.2 and 1.8, but the values outside it are NOT
# clipped into the end bins: doing so builds two artificial spikes in the
# tails, where the neighbouring bins are small, and those spikes read as
# real modes. The out-of-range counts are stated under the axis instead.
LO, HI = -0.2, 1.8
ratio = d.recovery_ratio_adj
inrange = ratio[(ratio >= LO) & (ratio <= HI)]
n_lo, n_hi = int((ratio < LO).sum()), int((ratio > HI).sum())
bins = np.linspace(LO, HI, 61)
ax.hist(inrange, bins=bins, color="#c9ccd1", lw=0, zorder=2)
ax.axvline(THR, color=INK, lw=1.2, zorder=4)
ax.set_xlabel("Climate-adjusted recovery ratio")
ax.set_ylabel("Pixels")
ax.set_xlim(LO, HI)
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))

# headroom for the outcome strip, drawn inside the axes
top = ax.get_ylim()[1]
ax.set_ylim(0, top * 1.42)
ys, h = top * 1.30, top * 0.085
left = -0.2
for frac, col, lab in ((f_never, ORANGE, "never recovered"),
                       (f_decl, BLUE, "recovered, declined again"),
                       (f_hold, GREEN, "recovered, above threshold")):
    wdt = frac * 2.0
    ax.barh(ys, wdt, left=left, height=h, color=col, lw=0, zorder=5,
            label=f"{lab}  {frac:.1%}")
    ax.barh(ys, 0.014, left=left + wdt - 0.007, height=h, color="white",
            lw=0, zorder=6)                      # 2px surface gap
    if wdt > 0.38:                               # direct label where it fits
        ax.text(left + wdt / 2, ys, f"{frac:.1%}", ha="center", va="center",
                fontsize=6.5, color="white", fontweight="bold", zorder=7)
    left += wdt
ax.legend(loc="upper left", bbox_to_anchor=(0.015, 0.875), frameon=False,
          fontsize=5.9, handlelength=0.9, handletextpad=0.45,
          labelspacing=0.28, borderaxespad=0.0)
ax.text(THR - 0.04, top * 0.02, "80% threshold", fontsize=6.5, ha="right",
        va="bottom", color=INK, rotation=90)
# State the truncation inside the panel, in the empty upper right, rather
# than folding the out-of-range pixels back into the end bins.
ax.text(0.985, 0.70, f"axis truncated\n{n_lo:,} pixels below {LO:g}\n"
                     f"{n_hi:,} above {HI:g}",
        transform=ax.transAxes, fontsize=6.0, color=MUTED, ha="right",
        va="top", linespacing=1.4)
tidy(ax); tag(ax, "b", "A third sits below threshold")

# c. western vs southeastern 
ax = axes[2]
west = (d.longitude < SPLIT).values
for lab, sel, col, ls in (("Western Gulf", west, ORANGE, "-"),
                          ("Southeastern Gulf", ~west, BLUE, "--")):
    v = np.sort(d.recovery_ratio_adj.values[sel])
    y = np.arange(1, len(v) + 1) / len(v)
    ax.step(v, y, where="post", color=col, lw=1.8, ls=ls, zorder=3,
            label=f"{lab}  (n = {sel.sum():,})")
    m = float(np.median(v))
    ax.plot([m], [0.5], "o", ms=5, color=col, mec="white", mew=0.8, zorder=5)
    ax.annotate(f"median {m:.2f}", (m, 0.5),
                (0.04 if col == ORANGE else 0.52, 0.80 if col == ORANGE else 0.28),
                textcoords="axes fraction", fontsize=6.5, color=col,
                va="center", ha="left",
                arrowprops=dict(arrowstyle="-", color=col, lw=0.5,
                                shrinkA=0, shrinkB=3))
ax.axvline(THR, color=INK, lw=1.0, zorder=2)
ax.set_xlabel("Climate-adjusted recovery ratio")
ax.set_ylabel("Cumulative share of pixels")
ax.set_xlim(-0.2, 1.8); ax.set_ylim(0, 1.02)
ax.legend(frameon=False, loc="upper left", handlelength=1.2, fontsize=6.0,
          handletextpad=0.45, labelspacing=0.3, borderaxespad=0.2)
tidy(ax); tag(ax, "c", "The west lags the southeast")

fig.subplots_adjust(wspace=0.46)
out = ROOT / "figures" / "fig5_recovery.png"
try:
    fig.savefig(out)
except OSError:
    import time
    alt = out.with_name(f"fig5_recovery_{time.strftime('%H%M%S')}.png")
    fig.savefig(alt)
    print(f"\n!! could not overwrite {out.name} - close it in any viewer")
    print(f"!! saved as {alt.name} instead\n")
    out = alt
print("saved", out)

print("")
print(f"  never recovered            {f_never:.1%}      text: 26.7%")
print(f"  recovered, declined again  {f_decl:.1%}       text: 9.6%")
print(f"  below threshold now        {f_never + f_decl:.1%}      text: 36.3%")
print(f"  median years to recovery   {med:.0f}          text: 8")
print(f"  west median (adjusted)     {d.recovery_ratio_adj[west].median():.3f}     text: 0.592")
print(f"  southeast median           {d.recovery_ratio_adj[~west].median():.3f}     text: 0.929")
print("")
print("PANEL b AXIS TRUNCATION")
print(f"  below {LO:g}                 {n_lo:,} pixels ({n_lo/len(d):.2%})")
print(f"  above {HI:g}                  {n_hi:,} pixels ({n_hi/len(d):.2%})")
print(f"  plotted                    {len(inrange):,} pixels")
