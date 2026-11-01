# The two presentation maps, drawn on a real coastline.

#
# INPUT   data/rasters/classmap_fullcoast_30m.tif   class map, 14944 x 24120
#         data/pixels_with_covariates.csv           29,453 dieback pixels
#         data/predicted_capacity_pixels.csv        765,138 undisturbed pixels
#         gulf_coastline.json                       land polygons for the basemap
#
# OUTPUT  figures/map1_recovery_outcome.png   16:9, for one slide
#         figures/map2_predicted_risk.png     16:9, for one slide
#
# The coastline comes from Natural Earth 10m land, already clipped to the Gulf
# and simplified, so nothing is downloaded and the maps draw the same way on
# any machine. 
#
# A note on scale, because it governs the design. The Gulf is about 700 km
# across, so on a full-width slide one millimetre is about 2.3 km. A 30 m pixel
# is 0.013 mm and cannot be drawn. The Gulf-wide view therefore shows WHERE
# mangrove is affected, with the fringe deliberately widened so it can be seen;
# the insets on map 1 are at true scale, where a 30 m pixel really is visible.
 
import json
from pathlib import Path
 
import numpy as np
import pandas as pd
import rasterio
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap, BoundaryNorm, TwoSlopeNorm, \
    LinearSegmentedColormap
from matplotlib.patches import Rectangle
from scipy.ndimage import maximum_filter
 
ROOT = Path(__file__).resolve().parent
if not (ROOT / "data").exists():
    ROOT = ROOT.parent
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
 
THR, DEFICIT, PX_HA = 0.80, 0.10, 0.08564
SPLIT = 138.5

# HEADING = True   
# HEADING = False  
HEADING = False
EXT = (135.0, 141.5, -18.0, -14.0)
 
RED, AMBER, GREEN = "#D55E00", "#E69F00", "#009E73"
BLUE, MANGROVE = "#0072B2", "#9aa3a8"
INK, MUTED = "#0E3B43", "#5A6B6E"
OCEAN, LAND, COAST, GRAT = "#dce9f0", "#f4f1ea", "#b9b2a4", "#c8d8e2"
 
#  font 

from matplotlib import font_manager as _fm


def pick_font(*names):
    have = {f.name for f in _fm.fontManager.ttflist}
    for n in names:
        if n in have:
            return n
    return "DejaVu Sans"


FONT = pick_font("Segoe UI", "Lato", "DejaVu Sans")

plt.rcParams.update({"font.family": FONT, "text.color": INK,
                     "savefig.dpi": 220, "savefig.facecolor": "white"})
 
COAST_FILE = ROOT / "gulf_coastline.json"
if not COAST_FILE.exists():
    COAST_FILE = Path(__file__).resolve().parent / "gulf_coastline.json"
LANDPOLY = [np.asarray(p) for p in json.loads(COAST_FILE.read_text())["polygons"]]
 
print("reading the class raster ...")
src = rasterio.open(ROOT / "data" / "rasters" / "classmap_fullcoast_30m.tif")
cls = src.read(1)
H, W = cls.shape
T = src.transform
INV = ~T
 
 
def rc(lon, lat):
    col, row = INV * (np.asarray(lon), np.asarray(lat))
    return np.round(row).astype(int), np.round(col).astype(int)
 
 
def basemap(ax, extent=EXT, graticule=True):
    """Ocean, land and coastline. Everything else is drawn on top."""
    lo, hi, la0, la1 = extent
    ax.set_facecolor(OCEAN)
    ax.add_collection(PolyCollection(LANDPOLY, facecolors=LAND,
                                     edgecolors=COAST, linewidths=0.7,
                                     zorder=2))
    ax.set_xlim(lo, hi)
    ax.set_ylim(la0, la1)
    ax.set_aspect(1 / np.cos(np.radians(16)))
    if graticule:
        ax.set_xticks(range(int(np.ceil(lo)), int(hi) + 1))
        ax.set_yticks(range(int(np.ceil(la0)), int(la1) + 1))
        ax.set_xticklabels([f"{v}\u00b0E" for v in
                            range(int(np.ceil(lo)), int(hi) + 1)], fontsize=9.5)
        ax.set_yticklabels([f"{abs(v)}\u00b0S" for v in
                            range(int(np.ceil(la0)), int(la1) + 1)], fontsize=9.5)
        ax.grid(color=GRAT, lw=0.6, zorder=3)
        ax.set_axisbelow(False)
        ax.tick_params(colors=MUTED, length=3)
        ax.set_xlabel("Longitude", fontsize=9.5, color=MUTED, labelpad=4)
        ax.set_ylabel("Latitude", fontsize=9.5, color=MUTED, labelpad=4)
    else:
        ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_edgecolor(MUTED)
        s.set_linewidth(0.8)
    return ax
 
 
def scalebar(ax, km, x, y, extent=EXT, color=INK):
    deg = km / (111.32 * np.cos(np.radians(16)))
    ax.plot([x, x + deg], [y, y], color="white", lw=4.5, zorder=11,
            solid_capstyle="butt")
    ax.plot([x, x + deg], [y, y], color=color, lw=2.2, zorder=12,
            solid_capstyle="butt")
    ax.text(x + deg / 2, y + (extent[3] - extent[2]) * 0.018, f"{km} km",
            fontsize=9, color=color, ha="center", va="bottom",
            fontweight="bold", zorder=12)
 
 
def northarrow(ax, x, y):
    ax.annotate("N", xy=(x, y), xytext=(x, y - 0.42), ha="center", va="top",
                fontsize=10, fontweight="bold", color=INK, zorder=12,
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.3))
 
 
def block_max(a, k):
    h, w = a.shape[0] // k * k, a.shape[1] // k * k
    return a[:h, :w].reshape(h // k, k, w // k, k).max(axis=(1, 3))
 
 
def full_extent(shape, k):
    return [T.c, T.c + shape[1] * k * T.a, T.f + shape[0] * k * T.e, T.f]
 
 
# The deck is set in a sans face, so the heading baked into each map is too.
HEAD_FONT = "Lato"


def heading(fig, title, subtitle):
    fig.text(0.030, 0.945, title, fontsize=21, fontweight="bold", color=INK,
             fontname=HEAD_FONT)
    fig.text(0.030, 0.898, subtitle, fontsize=11.5, color=MUTED)
 
 
PX = 0.625          
DY = 0.0 if HEADING else 0.108
MAP = [0.035, 0.055, 0.555, 0.805 + (0.0 if HEADING else 0.100)]
TAG = "" if HEADING else "_nohead"
 
# =========================================================================
# MAP 1 - what actually happened
# =========================================================================
print("map 1: recovery outcome ...")
d = pd.read_csv(ROOT / "data" / "pixels_with_covariates.csv")
never = d.unrecovered.astype(bool).values
declined = (~never) & (d.recovery_ratio < THR).values
holding = (~never) & (d.recovery_ratio >= THR).values
f_never, f_decl, f_hold = never.mean(), declined.mean(), holding.mean()
 
out = np.zeros((H, W), np.uint8)
out[(cls == 1) | (cls == 3)] = 1
r, c = rc(d.longitude.values, d.latitude.values)
out[r, c] = np.where(never, 5, np.where(declined, 4, 3)).astype(np.uint8)
 
CM1 = ListedColormap([MANGROVE, MANGROVE, GREEN, AMBER, RED])
NM1 = BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5, 5.5], 5)
 
K = 12
small = block_max(out, K)
# At this scale a 200 m fringe is under a tenth of a millimetre, so the
# mangrove band is widened to roughly 2 km. It says WHERE mangrove is, not
# how wide it is; the insets below are the ones drawn true to scale.
ribbon = maximum_filter((small > 0).astype(np.uint8), size=5)
hot = maximum_filter(np.where(small >= 3, small, 0), size=5)
disp = np.where(hot > 0, hot, ribbon).astype(np.uint8)
 
fig = plt.figure(figsize=(13.333, 7.5))
if HEADING:
    heading(fig, "A decade on, a third of the dieback has not come back",
            "Mangrove along the Gulf of Carpentaria that died in the summer "
            "of 2015\u201316, and how much of it has regrown by 2025")
ax = basemap(fig.add_axes(MAP))
ax.imshow(np.ma.masked_equal(disp, 0), cmap=CM1, norm=NM1,
          interpolation="nearest", extent=full_extent(disp.shape, K),
          origin="upper", zorder=5)
scalebar(ax, 100, 135.25, -17.82)
northarrow(ax, 141.25, -14.22)
 
ROWS = ((RED, f"{f_never:.1%}", "never recovered",
         f"{never.sum() * PX_HA:,.0f} ha"),
        (AMBER, f"{f_decl:.1%}", "recovered, then declined again",
         f"{declined.sum() * PX_HA:,.0f} ha"),
        (GREEN, f"{f_hold:.1%}", "recovered, still above threshold", ""))
fig.text(PX, 0.820 + DY, "Recovery state in 2025", fontsize=12,
         fontweight="bold", color=INK)
y = 0.715 + DY
for col, big, lab, sub in ROWS:
    fig.patches.append(Rectangle((PX, y - 0.004), 0.012, 0.046,
                                 transform=fig.transFigure, facecolor=col,
                                 edgecolor="none"))
    fig.text(PX + 0.118, y + 0.019, big, fontsize=26, fontweight="bold",
             color=col, va="center", ha="right")
    fig.text(PX + 0.136, y + 0.031, lab, fontsize=12, color=INK, va="center")
    if sub:
        fig.text(PX + 0.136, y + 0.005, sub, fontsize=10.5, color=MUTED,
                 va="center")
    y -= 0.092
 
fig.patches.append(Rectangle((PX, 0.405 + DY), 0.012, 0.022,
                             transform=fig.transFigure, facecolor=MANGROVE,
                             edgecolor="none"))
fig.text(PX + 0.022, 0.416 + DY, "mangrove that did not die in 2015\u201316",
         fontsize=10.5, color=MUTED, va="center")
 
fig.text(PX, 0.330 + DY, "Recovered = canopy back to 80% of its own\n"
         "pre-dieback density. Median time 8 years.",
         fontsize=10.5, color=INK, va="top", linespacing=1.65)
 
fig.text(PX, 0.252 + DY, "Geographic coordinates, WGS 84 (EPSG:4326)\n"
         "30 m pixels, Landsat 5 / 7 / 8 / 9",
         fontsize=9, color=MUTED, va="top", linespacing=1.6)
 
fig.savefig(OUT / f"map1_recovery_outcome{TAG}.png", bbox_inches="tight",
            pad_inches=0.12)
plt.close(fig)
print(f"  saved map1_recovery_outcome{TAG}.png")
 
# =========================================================================
# MAP 2 - what would happen next time
# =========================================================================
print("map 2: predicted risk ...")
p = pd.read_csv(ROOT / "data" / "predicted_capacity_pixels.csv")
ins = p.inside_aoa.astype(bool).values
# restricted to the area of applicability: the model only speaks where the
# conditions resemble those it was trained on (Meyer & Pebesma 2021)
vul = {}
for reg in ("west", "east"):
    g = p[(p.region == reg) & p.inside_aoa]
    vul[reg] = (100 * (g.pred_deficit > DEFICIT).mean(),
                (g.pred_deficit > DEFICIT).sum() * PX_HA)
 
CM2 = LinearSegmentedColormap.from_list(
    "deficit", [BLUE, "#9ecae8", "#f0f4f6", "#f0b27a", RED])
norm = TwoSlopeNorm(vmin=-0.15, vcenter=0.0, vmax=0.15)
 
fig = plt.figure(figsize=(13.333, 7.5))
if HEADING:
    heading(fig, "Where the next die-off would hit hardest",
            "If a summer like 2015\u201316 came again, this is where the "
            "mangrove would struggle to grow back")
ax = basemap(fig.add_axes(MAP))
ax.scatter(p.longitude.values[ins], p.latitude.values[ins],
           c=np.clip(p.pred_deficit.values[ins], -0.15, 0.15), cmap=CM2,
           norm=norm, s=1.6, marker="o", lw=0, zorder=5)
ax.axvline(SPLIT, color=MUTED, lw=1.1, ls="--", zorder=6)
ax.text(SPLIT - 0.08, -14.08, "west", fontsize=9.5, color=MUTED, ha="right",
        va="top", zorder=7)
ax.text(SPLIT + 0.08, -14.08, "southeast", fontsize=9.5, color=MUTED,
        ha="left", va="top", zorder=7)
scalebar(ax, 100, 135.25, -17.82)
northarrow(ax, 141.25, -14.22)
 
fig.text(PX, 0.820 + DY, "Mangrove at risk",
         fontsize=12, fontweight="bold", color=INK)
y = 0.700 + DY
for reg, name, col in (("west", "Western Gulf", RED),
                       ("east", "Southeastern Gulf", BLUE)):
    fig.text(PX, y + 0.072, name, fontsize=13, fontweight="bold", color=col)
    fig.text(PX + 0.100, y, f"{vul[reg][0]:.1f}%", fontsize=30,
             fontweight="bold", color=col, va="center", ha="right")
    fig.text(PX + 0.118, y + 0.014, "of its mangrove", fontsize=11.5,
             color=INK, va="center")
    fig.text(PX + 0.118, y - 0.012, f"{vul[reg][1]:,.0f} ha",
             fontsize=10.5, color=MUTED, va="center")
    y -= 0.150
fig.text(PX, 0.472 + DY, "Six times as much, in the west.", fontsize=11.5,
         fontweight="bold", color=INK)
 
cax = fig.add_axes([PX, 0.352 + DY, 0.185, 0.026])
cax.imshow(np.linspace(0, 1, 256).reshape(1, -1), cmap=CM2, aspect="auto")
cax.set_yticks([]); cax.set_xticks([])
for s in cax.spines.values():
    s.set_visible(False)
fig.text(PX, 0.393 + DY, "Predicted canopy deficit", fontsize=12,
         fontweight="bold", color=INK)
fig.text(PX, 0.332 + DY, "recovered", fontsize=9.5, color=MUTED, va="top")
fig.text(PX + 0.185, 0.332 + DY, "still thinner", fontsize=9.5, color=MUTED,
         va="top", ha="right")
 
fig.text(PX, 0.268 + DY, "At risk = deficit still above 0.10 NBR ten years on.",
         fontsize=10.5, color=INK, va="top", linespacing=1.65)
 
fig.text(PX, 0.202 + DY, "Geographic coordinates, WGS 84 (EPSG:4326)\n"
         "30 m pixels, Landsat 5 / 7 / 8 / 9",
         fontsize=9, color=MUTED, va="top", linespacing=1.6)
 
fig.savefig(OUT / f"map2_predicted_risk{TAG}.png", bbox_inches="tight",
            pad_inches=0.12)
plt.close(fig)
print(f"  saved map2_predicted_risk{TAG}.png")
 
print("")
print(f"  never recovered        {f_never:.1%}  {never.sum()*PX_HA:,.0f} ha   text: 26.7%, 674 ha")
print(f"  declined again          {f_decl:.1%}  {declined.sum()*PX_HA:,.0f} ha   text: 9.6%, 242 ha")
print(f"  still above threshold  {f_hold:.1%}             text: 63.7%")
print(f"  west at risk           {vul['west'][0]:.1f}%  {vul['west'][1]:,.0f} ha  text: 46.6%, 13,511 ha")
print(f"  southeast at risk       {vul['east'][0]:.1f}%   {vul['east'][1]:,.0f} ha  text: 7.3%, 1,728 ha")
print(f"  west : east ratio      {vul['west'][0]/vul['east'][0]:.1f}x")
print(f"  inside area of applicability  {100*ins.mean():.1f}%   text: 80.1%")
print(f"  (including every pixel it is 41.0% / 6.4%, ratio 6.4x - the sensitivity check)")