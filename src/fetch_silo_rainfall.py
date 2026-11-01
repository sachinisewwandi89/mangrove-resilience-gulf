"""
fetch_silo_rainfall.py
----------------------
Builds the external rainfall covariates for notebook 03, section 1b.

"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

try:
    import requests
    import xarray as xr
except ImportError as e:
    sys.exit(f"missing package: {e.name}\n"
             "  conda install -c conda-forge xarray netcdf4 requests -y")

# configuration 
YEAR_EVENT     = 2015
YEARS_RECOVERY = list(range(2016, 2025))      # 2016-2024

# Gulf study area, padded by ~0.2 deg so the nearest-neighbour join in the
# notebook never falls off the edge of the grid.
LON0, LON1 = 134.8, 141.7
LAT0, LAT1 = -18.2, -13.8

BASE = ("https://s3-ap-southeast-2.amazonaws.com/silo-open-data"
        "/Official/annual/monthly_rain/{year}.monthly_rain.nc")

ROOT = Path(__file__).resolve().parent
RAW  = ROOT / "data" / "silo_raw"
OUT  = ROOT / "data" / "covariates"
RAW.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)


# download 
def fetch(year):
    dst = RAW / f"{year}.monthly_rain.nc"
    if dst.exists() and dst.stat().st_size > 1_000_000:
        print(f"  {year}  already downloaded")
        return dst
    url = BASE.format(year=year)
    print(f"  {year}  downloading ...", end="", flush=True)
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        tmp = dst.with_suffix(".part")
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
        tmp.replace(dst)
    print(f" {dst.stat().st_size/1e6:.0f} MB")
    return dst


#  read one year, return annual total on the Gulf subset 
def annual_total(path):
    with xr.open_dataset(path) as ds:
        # SILO names the variable monthly_rain, but detect it rather than
        # assume: a silently wrong variable would poison the covariate.
        names = [v for v in ds.data_vars
                 if ds[v].ndim == 3 and "rain" in v.lower()]
        if not names:
            names = [v for v in ds.data_vars if ds[v].ndim == 3]
        if not names:
            raise RuntimeError(f"no 3-D variable found in {path.name}")
        da = ds[names[0]]

        ynm = "lat" if "lat" in da.dims else "latitude"
        xnm = "lon" if "lon" in da.dims else "longitude"

        
        y0, y1 = (LAT0, LAT1) if da[ynm][0] < da[ynm][-1] else (LAT1, LAT0)
        sub = da.sel({xnm: slice(LON0, LON1), ynm: slice(y0, y1)})

        if sub.size == 0:
            raise RuntimeError(
                f"{path.name}: the Gulf window is empty. Check the "
                f"coordinate ranges: {xnm} {float(da[xnm].min()):.1f} to "
                f"{float(da[xnm].max()):.1f}, {ynm} "
                f"{float(da[ynm].min()):.1f} to {float(da[ynm].max()):.1f}")

        total = sub.sum(dim=[d for d in sub.dims if d not in (ynm, xnm)],
                        skipna=False)
        return total.rename("value"), xnm, ynm


def to_frame(da, xnm, ynm):
    df = da.to_dataframe().reset_index()
    df = df.rename(columns={xnm: "longitude", ynm: "latitude"})
    df = df[["longitude", "latitude", "value"]].dropna()
    return df


# main 
print("SILO monthly rainfall, Gulf of Carpentaria")
print(f"event year {YEAR_EVENT}, recovery {YEARS_RECOVERY[0]}-"
      f"{YEARS_RECOVERY[-1]}\n")

print("download:")
paths = {y: fetch(y) for y in [YEAR_EVENT] + YEARS_RECOVERY}

print("\nprocessing:")
ev, xnm, ynm = annual_total(paths[YEAR_EVENT])
df_event = to_frame(ev, xnm, ynm)
df_event.to_csv(OUT / "rain_2015_mm.csv", index=False)
print(f"  rain_2015_mm.csv       {len(df_event):,} cells, "
      f"{df_event.value.min():.0f}-{df_event.value.max():.0f} mm")

stack = []
for y in YEARS_RECOVERY:
    tot, _, _ = annual_total(paths[y])
    stack.append(tot)
rec = xr.concat(stack, dim="year").mean("year")
df_rec = to_frame(rec, xnm, ynm)
df_rec.to_csv(OUT / "rain_recovery_mm.csv", index=False)
print(f"  rain_recovery_mm.csv   {len(df_rec):,} cells, "
      f"{df_rec.value.min():.0f}-{df_rec.value.max():.0f} mm")

# the check that matters 
# The east-west recovery gradient splits at 138.5 E. If rainfall is to
# explain it, rainfall has to differ across that line.
for name, df in (("2015", df_event), ("2016-2024 mean", df_rec)):
    w = df[df.longitude < 138.5].value.median()
    e = df[df.longitude >= 138.5].value.median()
    print(f"\n{name}: west median {w:.0f} mm, east median {e:.0f} mm, "
          f"difference {e - w:+.0f} mm")
