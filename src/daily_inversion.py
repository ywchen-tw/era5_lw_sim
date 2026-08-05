#!/usr/bin/env python
"""Compute temperature-inversion metrics from downloaded ERA5 files.

Reads data/YYYY/MM/DD/era5_{plev,sfc}_YYYYMMDD.nc and writes
derived/YYYY/MM/DD/era5_inversion_YYYYMMDD.nc with SBI (profile scan),
dt_850_2m and dt_925_1000 metrics. See reanlib/inversion.py for definitions
and references.

Examples:
    python src/daily_inversion.py --year 2025 --month 1 --days 1
    python src/daily_inversion.py --year 2025 --month 1 --days 1-7 --check
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from era5_download import parse_days
from reanlib.config import inversion_path, load_config, plev_path, sfc_path
from reanlib.inversion import compute_inversion_dataset
from reanlib.io_era5 import open_era5


def check_report(ds) -> str:
    """Physical-plausibility summary of one day's metrics."""
    lines = []
    found = ds["sbi_found"].values.astype(bool)
    lines.append(f"grid points x times : {found.size}")
    lines.append(f"SBI found fraction  : {found.mean():.1%}")
    for name in ("sbi_strength", "sbi_depth_p", "sbi_depth_z",
                 "dt_850_2m", "dt_925_1000"):
        v = ds[name].values
        v = v[np.isfinite(v)]
        if v.size:
            lines.append(f"{name:<20}: min {v.min():7.2f}  median {np.median(v):7.2f}  "
                         f"max {v.max():7.2f} {ds[name].attrs['units']}")
        else:
            lines.append(f"{name:<20}: all NaN")
    sp_hpa = ds["sp"].values / 100.0
    n_low = int((sp_hpa < 1000).sum())
    n_masked = int(np.isnan(ds["dt_925_1000"].values).sum())
    lines.append(f"sp < 1000 hPa points: {n_low} (dt_925_1000 NaN: {n_masked})")
    top_ok = np.all(np.isnan(ds["sbi_top_p"].values)
                    | (ds["sbi_top_p"].values <= sp_hpa + 1e-6))
    lines.append(f"sbi_top_p <= sp     : {'OK' if top_ok else 'VIOLATED'}")
    return "\n".join("  " + s for s in lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--month", type=int, required=True)
    parser.add_argument("--days", nargs="+", required=True,
                        help="day numbers and/or A-B ranges, e.g. 1 2 5-7")
    parser.add_argument("--hours", type=int, nargs="+", default=None,
                        help="subset of UTC hours (default: all in the file)")
    parser.add_argument("--source", choices=["era5", "merra2"], default=None,
                        help="data source (default from config.yaml)")
    parser.add_argument("--config", default=None, help="path to config.yaml")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--check", action="store_true",
                        help="print physical-plausibility summary per day")
    args = parser.parse_args(argv)

    cfg = load_config(args.config, source=args.source)
    for day in parse_days(args.days, args.year, args.month):
        date = dt.date(args.year, args.month, day)
        target = inversion_path(cfg, date)
        if target.exists() and not args.overwrite:
            print(f"{date}: exists, skipping (use --overwrite to redo)")
            if args.check:
                print(check_report(open_era5(target)))
            continue

        missing = [p for p in (plev_path(cfg, date), sfc_path(cfg, date))
                   if not p.exists()]
        if missing:
            names = ", ".join(str(m) for m in missing)
            sys.exit(f"{date}: missing input {names}\n"
                     f"fix with: python src/era5_download.py --year {date.year} "
                     f"--month {date.month} --days {date.day}")

        ds_plev = open_era5(plev_path(cfg, date))
        ds_sfc = open_era5(sfc_path(cfg, date))
        if args.hours:
            keep = [t for t in ds_plev["valid_time"].values
                    if t.astype("datetime64[s]").item().hour in set(args.hours)]
            ds_plev = ds_plev.sel(valid_time=keep)
            ds_sfc = ds_sfc.sel(valid_time=keep, method=None)

        out = compute_inversion_dataset(ds_plev, ds_sfc, cfg)
        out.attrs["source_files"] = (f"{plev_path(cfg, date).name}, "
                                     f"{sfc_path(cfg, date).name}")
        target.parent.mkdir(parents=True, exist_ok=True)
        encoding = {v: {"zlib": True, "complevel": 4} for v in out.data_vars}
        out.to_netcdf(target, encoding=encoding)
        print(f"{date}: wrote {target} ({out.sizes['valid_time']} time steps)")
        if args.check:
            print(check_report(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
