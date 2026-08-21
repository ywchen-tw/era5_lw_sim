"""Winter (Nov-Jan) cloud-thickness statistics from Cloudnet at Ny-Alesund.

Computes phase-resolved SINGLE-LAYER geometric cloud thickness from the
ACTRIS Cloudnet target-classification product at the AWIPEV observatory,
Ny-Alesund (78.92 N, 11.93 E), for the two winters Nov 2016 - Jan 2017 and
Nov 2017 - Jan 2018 — the season-matched observational reference for the
CRE study's cloud-geometry comparison (cre_sim.py census-layers panel c).
The published climatology (Nomokonova et al. 2019, ACP 19, 4105-4126,
doi:10.5194/acp-19-4105-2019) pools ALL seasons of June 2016 - July 2017;
this script applies the same definitions to the winter months only:

  cloud bins           target_classification in {1, 3, 4, 5, 6, 7}
                       (cloud droplets and/or ice hydrometeors)
  liquid-precip screen profiles containing {2, 3, 6, 7} (drizzle/rain or
                       melting ice) are EXCLUDED, as in the paper's Fig. 10
  single layer         exactly one contiguous cloudy run of >= 3 bins
  thickness            (n_bins x bin spacing): upper border of top bin minus
                       lower border of bottom bin (their definition)
  phase                droplets {1, 3, 5, 7} and/or ice {4, 5, 6, 7}

Data: data/nyalesund/classification/YYYYMMDD_ny-alesund_classification.nc,
downloaded from https://cloudnet.fmi.fi (legacy products, CC-BY 4.0; the
underlying dataset of the paper is doi:10.1594/PANGAEA.898556).

Usage: conda run -n era5 python src/nya_cloudnet_winter.py
"""

from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import numpy as np
import xarray as xr

REPO = Path(__file__).resolve().parent.parent
CLS_DIR = REPO / "data" / "nyalesund" / "classification"
OUT = REPO / "derived" / "nyalesund" / "nya_winter_thickness.json"

CLOUD_CATS = (1, 3, 4, 5, 6, 7)
DROPLET_CATS = (1, 3, 5, 7)
ICE_CATS = (4, 5, 6, 7)
LIQ_PRECIP_CATS = (2, 3, 6, 7)
MIN_BINS = 3                       # single-layer definition of the paper


def main() -> int:
    files = sorted(glob.glob(str(CLS_DIR / "*_classification.nc")))
    if not files:
        sys.exit(f"no classification files under {CLS_DIR}")
    thick = {"liquid": [], "mixed": [], "ice": []}
    n_prof = n_cloudy = n_single = n_precip_excl = 0
    for f in files:
        ds = xr.open_dataset(f)
        tc = ds["target_classification"].values          # (time, height)
        dz = float(np.median(np.diff(ds["height"].values)))
        ds.close()
        cloudy = np.isin(tc, CLOUD_CATS)
        has_liq_precip = np.isin(tc, LIQ_PRECIP_CATS).any(axis=1)
        n_prof += tc.shape[0]
        for it in range(tc.shape[0]):
            col = cloudy[it]
            if not col.any():
                continue
            n_cloudy += 1
            if has_liq_precip[it]:
                n_precip_excl += 1
                continue
            # contiguous cloudy runs (any clear bin splits, 20 m resolution)
            idx = np.where(col)[0]
            splits = np.where(np.diff(idx) > 1)[0]
            runs = np.split(idx, splits + 1)
            if len(runs) != 1 or runs[0].size < MIN_BINS:
                continue
            n_single += 1
            g = runs[0]
            cats = tc[it, g]
            has_drop = np.isin(cats, DROPLET_CATS).any()
            has_ice = np.isin(cats, ICE_CATS).any()
            phase = ("mixed" if has_drop and has_ice
                     else "liquid" if has_drop else "ice")
            thick[phase].append(g.size * dz / 1000.0)    # km, border-to-border

    stats = {}
    print(f"{len(files)} days, {n_prof} profiles, {n_cloudy} cloudy, "
          f"{n_precip_excl} excluded (liquid precip), "
          f"{n_single} single-layer")
    print(f"{'phase':7s} {'n':>7s} {'median':>7s} {'mean':>7s} {'p99':>7s} km")
    for ph, v in thick.items():
        v = np.array(v)
        stats[ph] = {"n": int(v.size), "median_km": float(np.median(v)),
                     "mean_km": float(np.mean(v)),
                     "p99_km": float(np.percentile(v, 99))}
        print(f"{ph:7s} {v.size:7d} {np.median(v):7.2f} {np.mean(v):7.2f} "
              f"{np.percentile(v, 99):7.2f}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(
        {"site": "Ny-Alesund AWIPEV (78.92N, 11.93E)",
         "period": "Nov 2016 - Jan 2017 + Nov 2017 - Jan 2018",
         "source": "ACTRIS Cloudnet legacy classification "
                   "(cloudnet.fmi.fi, CC-BY 4.0); method of Nomokonova "
                   "et al. 2019 (doi:10.5194/acp-19-4105-2019)",
         "n_profiles": n_prof, "n_cloudy": n_cloudy,
         "n_liquid_precip_excluded": n_precip_excl,
         "n_single_layer": n_single, "thickness": stats}, indent=1))
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
