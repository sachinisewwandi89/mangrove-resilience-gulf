# Study area map 

import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle
from scipy.ndimage import maximum_filter

ROOT = Path(__file__).resolve().parent
if not (ROOT / "data").exists():
    ROOT = ROOT.parent
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

ORANGE, BLUE, GREEN, AMBER = "#D55E00", "#0072B2", "#009E73", "#E69F00"
INK, MUTED = "#0E3B43", "#5A6B6E"
OCEAN, LAND, COAST = "#e9f1f6", "#f6f3ee", "#b9c3c9"

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
    "font.size": 11, "text.color": INK, "axes.labelcolor": INK,
    "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.16,
    "font.family": FONT, "savefig.transparent": True})

EXT = (134.95, 141.62, -18.05, -13.95)
SPLIT = 138.5

land = json.loads((ROOT / "gulf_coastline.json").read_text())["polygons"]
LANDPOLY = [np.asarray(p) for p in land]

und = pd.read_csv(ROOT / "data" / "undisturbed_census.csv",
                  usecols=["longitude", "latitude"])
dieback = pd.read_csv(ROOT / "data" / "pixels_with_covariates.csv",
                      usecols=["longitude", "latitude"])

# thicken the ribbon 
# A 30 m pixel is 0.013 mm on a slide 700 km across, so the mangrove band has
# to be drawn wider than it is. Both layers get the same treatment, so their
# relative extent is preserved.
CELL = 0.004          # about 440 m
def ribbon(df, grow):
    gx = np.floor((df.longitude.values - EXT[0]) / CELL).astype(int)
    gy = np.floor((df.latitude.values - EXT[2]) / CELL).astype(int)
    nx = int((EXT[1] - EXT[0]) / CELL) + 1
    ny = int((EXT[3] - EXT[2]) / CELL) + 1
    ok = (gx >= 0) & (gx < nx) & (gy >= 0) & (gy < ny)
    grid = np.zeros((ny, nx), bool)
    grid[gy[ok], gx[ok]] = True
    grid = maximum_filter(grid, size=grow)
    yy, xx = np.nonzero(grid)
    return EXT[0] + (xx + 0.5) * CELL, EXT[2] + (yy + 0.5) * CELL

ux, uy = ribbon(und, 3)
dx, dy = ribbon(dieback, 3)

fig, ax = plt.subplots(figsize=(9.4, 5.6))
ax.set_facecolor(OCEAN)
ax.add_collection(PolyCollection(LANDPOLY, facecolors=LAND, edgecolors=COAST,
                                 linewidths=0.8, zorder=2))
ax.scatter(ux, uy, s=1.6, c=GREEN, lw=0, zorder=3, alpha=0.85)
ax.scatter(dx, dy, s=2.4, c=ORANGE, lw=0, zorder=4)

ax.plot([SPLIT, SPLIT], [EXT[2], EXT[3]], color=MUTED, lw=1.3, ls="--", zorder=5)
ax.text(SPLIT - 0.08, -17.92, "west", fontsize=10, color=MUTED, ha="right")
ax.text(SPLIT + 0.08, -17.92, "southeast", fontsize=10, color=MUTED, ha="left")

ax.text(139.20, -14.60, "GULF OF\nCARPENTARIA", fontsize=10.5, color="#7f98a6",
        ha="center", va="center", linespacing=1.5, zorder=3,
        fontweight="bold")
ax.text(135.72, -16.30, "NORTHERN\nTERRITORY", fontsize=9.5, color="#a69c8c",
        ha="center", va="center", linespacing=1.5, zorder=3)
ax.text(140.45, -17.98, "QUEENSLAND", fontsize=9.5, color="#a69c8c",
        ha="center", va="bottom", zorder=3)

# legend 
items = [(GREEN, f"Mangrove, undisturbed   {len(und):,} pixels"),
         (ORANGE, f"Died back in 2015\u201316   {len(dieback):,} pixels")]
for i, (col, lab) in enumerate(items):
    y = 0.225 - i * 0.070
    ax.plot([0.030, 0.062], [y, y], color=col, lw=5, solid_capstyle="round",
            transform=ax.transAxes, zorder=8, clip_on=False)
    ax.text(0.078, y, lab, transform=ax.transAxes, fontsize=10.5, color=INK,
            ha="left", va="center", zorder=8)

# scale bar
KM = 100.0
dlon = KM / (111.32 * np.cos(np.radians(16.0)))
x0, y0 = 139.05, -16.20
ax.plot([x0, x0 + dlon], [y0, y0], color=INK, lw=2.6, solid_capstyle="butt",
        zorder=8)
ax.text(x0 + dlon / 2, y0 + 0.10, "100 km", fontsize=10, color=INK,
        ha="center", va="bottom", zorder=8)

# north
ax.annotate("", xy=(135.30, -14.22), xytext=(135.30, -14.62),
            arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.4), zorder=8)
ax.text(135.30, -14.14, "N", fontsize=11, fontweight="bold", color=INK,
        ha="center", va="bottom", zorder=8)

ax.set_xlim(EXT[0], EXT[1])
ax.set_ylim(EXT[2], EXT[3])
ax.set_aspect(1 / np.cos(np.radians(16)))
ax.set_xticks(range(135, 142))
ax.set_yticks(range(-18, -13))
ax.set_xticklabels([f"{v}\u00b0E" for v in range(135, 142)], fontsize=10,
                   color=MUTED)
ax.set_yticklabels([f"{abs(v)}\u00b0S" for v in range(-18, -13)], fontsize=10,
                   color=MUTED)
ax.tick_params(length=3, color=MUTED)
# Name the axes. Degree ticks on their own leave the reader to infer what
# the numbers are, and say nothing about the coordinate system.
ax.set_xlabel("Longitude", fontsize=10.5, color=MUTED, labelpad=5)
ax.set_ylabel("Latitude", fontsize=10.5, color=MUTED, labelpad=5)
for sp in ax.spines.values():
    sp.set_edgecolor(COAST)
    sp.set_linewidth(0.9)

# inset
aus = json.loads((ROOT / "aus_outline.json").read_text())
polys = []
for part in aus["coordinates"]:
    polys.append(np.asarray(part[0]))
iax = ax.inset_axes([0.760, 0.620, 0.235, 0.345])
iax.set_facecolor(OCEAN)
iax.add_collection(PolyCollection(polys, facecolors=LAND, edgecolors=COAST,
                                  linewidths=0.6))
iax.add_patch(Rectangle((EXT[0], EXT[2]), EXT[1] - EXT[0], EXT[3] - EXT[2],
                        fill=False, edgecolor=ORANGE, lw=1.6, zorder=5))
iax.set_xlim(112, 155)
iax.set_ylim(-44, -9)
iax.set_aspect(1 / np.cos(np.radians(26)))
iax.set_xticks([])
iax.set_yticks([])
for sp in iax.spines.values():
    sp.set_edgecolor(COAST)


ax.text(0, 1.012, "135.0\u2013141.5\u00b0E, 14.0\u201318.0\u00b0S   \u00b7   "
        "Landsat 5, 7, 8 and 9, 1987\u20132025   \u00b7   30 m pixels, "
        "dry-season composites", transform=ax.transAxes, fontsize=10.5,
        color=MUTED, va="bottom")
ax.text(0, -0.168, "Geographic coordinates, WGS 84 (EPSG:4326).  Mangrove "
        "extent from Global Mangrove Watch v3 (Bunting et al. 2022).\n"
        "The band is drawn wider than 30 m so that it is visible at this "
        "scale.", transform=ax.transAxes,
        fontsize=9, color=MUTED, va="top", linespacing=1.5)

fig.savefig(OUT / "fig0_study_area.png")
plt.close(fig)
print("saved", OUT / "fig0_study_area.png")
print(f"  undisturbed {len(und):,}   dieback {len(dieback):,}   "
      f"total {len(und) + len(dieback):,} pixels")
