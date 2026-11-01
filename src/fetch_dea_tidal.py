"""
fetch_dea_tidal.py
------------------

"""

from pathlib import Path
import json
import sys
import urllib.parse

import numpy as np
import pandas as pd

try:
    import rasterio
    from rasterio.errors import RasterioIOError
    import requests
    from pyproj import Transformer
except ImportError as e:
    sys.exit(f"missing package: {e.name}\n"
             "  conda install -c conda-forge rasterio pyproj requests -y")

# COGs over HTTP
import os
os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif,.TIF")
os.environ.setdefault("GDAL_HTTP_MAX_RETRY", "3")
os.environ.setdefault("GDAL_HTTP_RETRY_DELAY", "2")

# configuration 
BBOX  = [135.0, -18.0, 141.5, -14.0]
STAC  = "https://explorer.dea.ga.gov.au/stac/search"

INTERTIDAL = "ga_s2ls_intertidal_cyear_3"
INTERTIDAL_YEAR = "2024"          # each year is a 3-year rolling composite


WATER_CANDIDATES = ["ga_ls_wo_fq_myear_3", "ga_ls_wo_fq_cyear_3"]
WATER_YEAR = "2023"               # used only for the cyear fallback
WATER_THRESHOLD = 0.20            # water seen >=20% of clear observations

TARGETS = ["data/phase3_pixels.csv", "data/undisturbed_census.csv"]

ROOT = Path(__file__).resolve().parent
OUT  = ROOT / "data" / "covariates"
OUT.mkdir(parents=True, exist_ok=True)


# points to sample 
pts = []
for t in TARGETS:
    p = ROOT / t
    if not p.exists():
        sys.exit(f"missing {p}")
    pts.append(pd.read_csv(p, usecols=["longitude", "latitude"]))
pts = pd.concat(pts, ignore_index=True).drop_duplicates().reset_index(drop=True)
print(f"{len(pts):,} pixels to sample\n")

to3577 = Transformer.from_crs("EPSG:4326", "EPSG:3577", always_xy=True)
X, Y = to3577.transform(pts.longitude.values, pts.latitude.values)


# STAC search 
def stac_items(collection, datetime_range):
    params = {
        "collections": collection,
        "bbox": json.dumps(BBOX),
        "datetime": datetime_range,
        "limit": 200,
    }
    url = f"{STAC}?{urllib.parse.urlencode(params)}"
    items, seen = [], 0
    while url and seen < 20:
        r = requests.get(url, timeout=90)
        if r.status_code != 200:
            return []
        j = r.json()
        items += j.get("features", [])
        url = next((l["href"] for l in j.get("links", [])
                    if l.get("rel") == "next"), None)
        seen += 1
    return items


def asset_href(item, *names):
    a = item.get("assets", {})
    for n in names:
        if n in a and a[n].get("href"):
            return a[n]["href"]
    # fall back to any asset whose key contains the name
    for n in names:
        for k, v in a.items():
            if n in k and v.get("href"):
                return v["href"]
    return None


def sample_band(items, band_names, label):
   
    out = np.full(len(pts), np.nan)
    used = 0
    for k, it in enumerate(items, 1):
        href = asset_href(it, *band_names)
        if not href:
            continue
        href = href.replace("s3://dea-public-data/",
                            "https://dea-public-data.s3.ap-southeast-2"
                            ".amazonaws.com/")
        try:
            with rasterio.open(href) as src:
                left, bottom, right, top = src.bounds
                inside = ((X >= left) & (X <= right) &
                          (Y >= bottom) & (Y <= top) & np.isnan(out))
                if not inside.any():
                    continue
                idx = np.flatnonzero(inside)
                vals = np.array([v[0] for v in src.sample(
                    zip(X[idx], Y[idx]))], dtype="float64")
                nod = src.nodata
                if nod is not None:
                    vals[vals == nod] = np.nan
                out[idx] = vals
                used += 1
        except RasterioIOError as e:
            print(f"    tile {k} unreadable: {str(e)[:80]}")
            continue
        if k % 10 == 0:
            print(f"    {k}/{len(items)} tiles, "
                  f"{np.isfinite(out).mean():.1%} covered")
    print(f"  {label}: {used} tiles used, "
          f"{np.isfinite(out).mean():.1%} of pixels have a value")
    return out


def write(name, values, unit):
    s = pd.DataFrame({"longitude": pts.longitude,
                      "latitude": pts.latitude,
                      "value": values}).dropna()
    if s.empty:
        print(f"  {name}: nothing to write")
        return
    s.to_csv(OUT / f"{name}.csv", index=False)
    w = s[s.longitude < 138.5].value.median()
    e = s[s.longitude >= 138.5].value.median()
    print(f"  wrote {name}.csv  {len(s):,} rows  "
          f"west {w:.2f}, east {e:.2f} {unit} (difference {e-w:+.2f})")


# 1. DEA Intertidal: elevation and exposure 
print("DEA Intertidal", INTERTIDAL, INTERTIDAL_YEAR)
items = stac_items(INTERTIDAL, "2016-01-01/2024-12-31")
print(f"  {len(items)} tiles found")
if items:
    elev = sample_band(items, ["elevation"], "elevation")
    write("tidal_elevation_m", elev, "m")
    expo = sample_band(items, ["exposure"], "exposure")
    write("tidal_exposure_pct", expo, "%")
else:
    print("  no tiles returned ")


# --- 2. Water frequency -> distance to the nearest permanent water ----
print("\nDEA Water Observations frequency")
witems, wname = [], None
for cand in WATER_CANDIDATES:
    rng = ("1986-01-01/2024-12-31" if "myear" in cand
           else f"{WATER_YEAR}-01-01/{WATER_YEAR}-12-31")
    witems = stac_items(cand, rng)
    if witems:
        wname = cand
        break
print(f"  {wname}: {len(witems)} tiles")

if witems:
    # Sample water frequency at the mangrove pixels themselves, and also
    # on a 200 m lattice across the study area, so that channels lying
    # OUTSIDE the mangrove mask are found. Distance is then measured to
    # the nearest lattice node that is water.
    step_m = 200.0
    gx = np.arange(X.min() - 5000, X.max() + 5000, step_m)
    gy = np.arange(Y.min() - 5000, Y.max() + 5000, step_m)
    GX, GY = np.meshgrid(gx, gy)
    GX, GY = GX.ravel(), GY.ravel()
    print(f"  lattice of {len(GX):,} nodes at {step_m:.0f} m")

    keepX, keepY = X, Y
    X, Y = GX, GY                      # sample_band works on X, Y
    n_pts = len(pts)
    pts_backup = pts
    pts = pd.DataFrame({"longitude": np.zeros(len(GX)),
                        "latitude": np.zeros(len(GX))})
    freq = sample_band(witems, ["frequency", "count_wet"], "water frequency")
    X, Y, pts = keepX, keepY, pts_backup

    is_water = np.isfinite(freq) & (freq >= WATER_THRESHOLD)
    print(f"  {is_water.sum():,} lattice nodes are water "
          f"(frequency >= {WATER_THRESHOLD})")

    if is_water.sum() > 10:
        from scipy.spatial import cKDTree
        tree = cKDTree(np.column_stack([GX[is_water], GY[is_water]]))
        dist, _ = tree.query(np.column_stack([X, Y]), k=1)
        # the lattice is 200 m, so resolution below that is not real
        write("dist_to_channel_m", dist, "m")
        print("  NOTE: measured against a 200 m lattice, so distances "
              "below ~200 m are not resolved. ")
    else:
        print("  too few water nodes ")

print("\nDone")
