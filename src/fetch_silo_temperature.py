"""
fetch_silo_temperature.py
"""

from pathlib import Path
import io
import os
import sys
import time

import numpy as np
import pandas as pd

try:
    import requests
except ImportError:
    sys.exit("conda install -c conda-forge requests -y")

# configuration 
EMAIL = os.environ.get("SILO_EMAIL", "")                     

GRID        = 0.05             # SILO grid spacing, degrees
YEAR_EVENT  = 2015
YEARS_REC   = (2016, 2024)
PAUSE_S     = 0.4              

TARGETS = ["data/phase3_pixels.csv", "data/undisturbed_census.csv"]
                                            

URL = "https://www.longpaddock.qld.gov.au/cgi-bin/silo/DataDrillDataset.php"

ROOT = Path(__file__).resolve().parent
OUT  = ROOT / "data" / "covariates"
CACHE = ROOT / "data" / "silo_points"
OUT.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)

if not EMAIL or "@" not in EMAIL:
    sys.exit("Set the SILO_EMAIL environment variable to your email "
             "address before running this script. SILO uses it as the "
             "API username.")



pts = []
for t in TARGETS:
    p = ROOT / t
    if not p.exists():
        sys.exit(f"missing {p}")
    df = pd.read_csv(p, usecols=["longitude", "latitude"])
    pts.append(df)
pts = pd.concat(pts, ignore_index=True)

cells = pd.DataFrame({
    "lon": (np.round(pts.longitude / GRID) * GRID).round(3),
    "lat": (np.round(pts.latitude  / GRID) * GRID).round(3),
}).drop_duplicates().sort_values(["lon", "lat"]).reset_index(drop=True)

print(f"{len(pts):,} pixels fall in {len(cells):,} SILO grid cells")
print(f"about {len(cells)*PAUSE_S/60:.0f} minutes of requests\n")



def fetch_point(lat, lon):
    """Daily series for one grid cell. Cached, so re-running is cheap."""
    key = CACHE / f"{lon:.3f}_{lat:.3f}.txt"
    if key.exists() and key.stat().st_size > 500:
        return key.read_text()
    r = requests.get(URL, timeout=90, params={
        "format":   "alldata",
        "lat":      f"{lat:.4f}",
        "lon":      f"{lon:.4f}",
        "start":    f"{YEAR_EVENT}0101",
        "finish":   f"{YEARS_REC[1]}1231",
        "username": EMAIL,
        "password": "apirequest",
    })
    r.raise_for_status()
    key.write_text(r.text)
    return r.text


def parse(text):
    """
    The 'alldata' response is a header block followed by a whitespace
    table. The header row is found rather than assumed, because a change
    in column order would otherwise be silently wrong.
    """
    lines = text.splitlines()
    hdr = None
    for i, ln in enumerate(lines):
        parts = ln.split()
        if "Date" in parts and any(p.startswith("T.Max") or p == "T.Max"
                                   or p.lower().startswith("max") for p in parts):
            hdr = i
            break
    if hdr is None:
        raise ValueError("could not find the header row")

    cols = lines[hdr].split()
    body = [ln for ln in lines[hdr + 1:] if ln.strip() and ln.split()[0][:1].isdigit()]
    df = pd.read_csv(io.StringIO("\n".join(body)), sep=r"\s+",
                     names=cols, engine="python")

    def pick(*cands):
        for c in cands:
            for col in df.columns:
                if col.lower().replace(".", "") == c.lower().replace(".", ""):
                    return col
        return None

    dcol = pick("Date", "Date2", "YYYY-MM-DD")
    xcol = pick("T.Max", "Tmax", "max_temp")
    ncol = pick("T.Min", "Tmin", "min_temp")
    if dcol is None or xcol is None:
        raise ValueError(f"columns not recognised: {list(df.columns)}")

    d = pd.to_datetime(df[dcol].astype(str), errors="coerce", format="mixed")
    out = pd.DataFrame({"year": d.dt.year,
                        "tmax": pd.to_numeric(df[xcol], errors="coerce")})
    if ncol is not None:
        out["tmin"] = pd.to_numeric(df[ncol], errors="coerce")
    return out.dropna(subset=["year"])


# --- probe first, so a schema change is caught on request one ---------
print("probe request ...")
probe = fetch_point(cells.lat[0], cells.lon[0])
print("---- first 14 lines of the response ----")
print("\n".join(probe.splitlines()[:14]))
print("----------------------------------------")
try:
    p = parse(probe)
    print(f"parsed {len(p):,} daily rows, years "
          f"{int(p.year.min())}-{int(p.year.max())}, "
          f"tmax {p.tmax.mean():.1f} C\n")
except Exception as e:
    sys.exit(f"\nPARSE FAILED: {e}\n")
             


#  all cells 
rows = []
for i, c in cells.iterrows():
    try:
        p = parse(fetch_point(c.lat, c.lon))
    except Exception as e:
        print(f"  {c.lon:.2f},{c.lat:.2f} failed: {e}")
        continue
    rec = p[p.year.between(*YEARS_REC)]
    rows.append({
        "longitude": c.lon, "latitude": c.lat,
        "tmax_2015":    p.loc[p.year == YEAR_EVENT, "tmax"].mean(),
        "tmax_recovery": rec.tmax.mean(),
        "tmin_recovery": rec.tmin.mean() if "tmin" in rec else np.nan,
    })
    if (i + 1) % 25 == 0:
        print(f"  {i+1:,} / {len(cells):,}")
    time.sleep(PAUSE_S)

t = pd.DataFrame(rows)
print(f"\n{len(t):,} cells retrieved")

for col, name in (("tmax_2015", "tmax_2015_C"),
                  ("tmax_recovery", "tmax_recovery_C"),
                  ("tmin_recovery", "tmin_recovery_C")):
    s = t[["longitude", "latitude", col]].dropna()
    if s.empty:
        continue
    s.columns = ["longitude", "latitude", "value"]
    s.to_csv(OUT / f"{name}.csv", index=False)
    w = s[s.longitude < 138.5].value.median()
    e = s[s.longitude >= 138.5].value.median()
    print(f"  {name}.csv  {len(s):,} cells  "
          f"west {w:.1f} C, east {e:.1f} C, difference {e-w:+.1f} C")
    print("done")


