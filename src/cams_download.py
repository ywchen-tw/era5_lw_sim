#!/usr/bin/env python
"""Download CAMS EGG4 monthly-mean CO2 and CH4 profiles (for the LW simulation).

Fetches the CAMS global greenhouse-gas reanalysis (EGG4) monthly means on
pressure levels from the Atmosphere Data Store (ADS) into
data/cams/cams_egg4_co2_ch4_YYYYMM.nc. EGG4 covers 2003-2020.

The ADS shares the ECMWF login with the CDS, but has its own API endpoint and
its own licence: log in at https://ads.atmosphere.copernicus.eu, accept the
CAMS licence on the EGG4 dataset page once, and the ~/.cdsapirc token works.

If ADS access is not set up yet, run the LW simulation with constant gases
instead: lrt_sim.py prep --fallback-constants (CO2 415 ppm, CH4 1.9 ppm).

Examples:
    python src/cams_download.py --year 2020 --month 1
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reanlib.config import REPO_ROOT, load_config

ADS_URL = "https://ads.atmosphere.copernicus.eu/api"
EGG4_LEVELS = ["1", "2", "3", "5", "7", "10", "20", "30", "50", "70", "100",
               "150", "200", "250", "300", "400", "500", "600", "700", "800",
               "850", "900", "925", "950", "1000"]

_ADS_HELP = """\
CAMS download failed. The ADS needs one-time setup even with working CDS
credentials:
  1. log in at https://ads.atmosphere.copernicus.eu with your ECMWF account
  2. open the "CAMS global greenhouse gas reanalysis (EGG4) monthly averaged
     fields" dataset page and accept its licence under "Terms of use"
  3. your ~/.cdsapirc token is reused automatically (the script overrides the
     endpoint URL)
Until then, use: lrt_sim.py prep --fallback-constants
Original error:
"""


def cams_path(year: int, month: int) -> Path:
    return REPO_ROOT / "data" / "cams" / f"cams_egg4_co2_ch4_{year:04d}{month:02d}.nc"


def read_cdsapirc_key() -> str:
    rc = Path.home() / ".cdsapirc"
    if not rc.exists():
        sys.exit("no ~/.cdsapirc found — set up CDS credentials first "
                 "(see README Setup)")
    for line in rc.read_text().splitlines():
        if line.strip().lower().startswith("key"):
            return line.split(":", 1)[1].strip()
    sys.exit("no 'key:' line found in ~/.cdsapirc")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--year", type=int, required=True, help="2003-2020 (EGG4 span)")
    parser.add_argument("--month", type=int, required=True)
    parser.add_argument("--area", type=float, nargs=4, default=None,
                        metavar=("N", "W", "S", "E"))
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    if not 2003 <= args.year <= 2020:
        sys.exit(f"EGG4 covers 2003-2020; {args.year} is outside")
    target = cams_path(args.year, args.month)
    if target.exists() and not args.force:
        print(f"{target} exists, skipping (use --force to re-download)")
        return 0

    cfg = load_config()
    area = list(args.area) if args.area is not None else cfg["area"]
    request = {
        "variable": ["carbon_dioxide", "methane"],
        "pressure_level": EGG4_LEVELS,
        "year": [f"{args.year:04d}"],
        "month": [f"{args.month:02d}"],
        "product_type": "monthly_mean",
        "data_format": "netcdf",
        "area": area,
    }

    import cdsapi
    try:
        client = cdsapi.Client(url=ADS_URL, key=read_cdsapirc_key(),
                               quiet=True, progress=False)
        target.parent.mkdir(parents=True, exist_ok=True)
        part = target.with_suffix(".nc.part")
        client.retrieve("cams-global-ghg-reanalysis-egg4-monthly",
                        request).download(str(part))
        os.replace(part, target)
    except Exception as exc:
        sys.exit(_ADS_HELP + f"  {exc}")
    print(f"wrote {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
