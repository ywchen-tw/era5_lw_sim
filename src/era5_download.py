#!/usr/bin/env python
"""Download ERA5 pressure-level and single-level data for given days/hours.

One CDS request per day per dataset (all hours bundled). Files land in
data/YYYY/MM/DD/era5_{plev,sfc}_YYYYMMDD.nc. Re-running skips complete files
and re-requests the union of hours if a file is missing some.

Examples:
    python src/era5_download.py --year 2025 --month 1 --days 1
    python src/era5_download.py --year 2025 --month 1 --days 1 2 5-7 --hours 0 12
    python src/era5_download.py --year 2025 --month 1 --days 15 --datasets plev --dry-run
"""

from __future__ import annotations

import argparse
import calendar
import concurrent.futures
import datetime as dt
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reanlib.config import load_config, plev_path, sfc_path, source_area
from reanlib.io_era5 import open_era5, require_cds_credentials

DATASET_IDS = {
    "plev": "reanalysis-era5-pressure-levels",
    "sfc": "reanalysis-era5-single-levels",
}


def parse_days(tokens: list[str], year: int, month: int) -> list[int]:
    """Expand day tokens like ['1', '2', '5-7'] into sorted unique day numbers."""
    ndays = calendar.monthrange(year, month)[1]
    days: set[int] = set()
    for tok in tokens:
        if "-" in tok:
            a, b = tok.split("-", 1)
            days.update(range(int(a), int(b) + 1))
        else:
            days.add(int(tok))
    bad = [d for d in days if not 1 <= d <= ndays]
    if bad:
        sys.exit(f"invalid day(s) {bad} for {year}-{month:02d} (month has {ndays} days)")
    return sorted(days)


def build_request(cfg: dict, kind: str, date: dt.date, hours: list[int],
                  area: list[float]) -> dict:
    request = {
        "product_type": ["reanalysis"],
        "year": [f"{date.year:04d}"],
        "month": [f"{date.month:02d}"],
        "day": [f"{date.day:02d}"],
        "time": [f"{h:02d}:00" for h in hours],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": area,
    }
    if kind == "plev":
        request["variable"] = cfg["download"]["plev_variables"]
        request["pressure_level"] = [str(p) for p in cfg["download"]["pressure_levels"]]
    else:
        request["variable"] = cfg["download"]["sfc_variables"]
    return request


def hours_in_file(path: Path, date: dt.date) -> set[int]:
    """Hours of `date` already present in an existing file; empty set if unreadable."""
    try:
        ds = open_era5(path)
    except Exception:
        return set()
    times = ds["valid_time"].values
    ds.close()
    hours = set()
    for t in times:
        ts = t.astype("datetime64[s]").item()
        if ts.date() == date:
            hours.add(ts.hour)
    return hours


def normalize_download(path: Path) -> None:
    """Merge a zipped multi-stream CDS delivery into a single netCDF, in place.

    The CDS ignores download_format "unarchived" when a request mixes GRIB
    streams (e.g. instantaneous t2m/skt/sp with accumulated fluxes) and
    returns a zip of one netCDF per stream with identical coordinates.
    """
    import tempfile
    import zipfile

    import xarray as xr

    if not zipfile.is_zipfile(path):
        return
    with tempfile.TemporaryDirectory(dir=path.parent) as tmp:
        with zipfile.ZipFile(path) as zf:
            members = zf.namelist()
            zf.extractall(tmp)
        parts = [xr.open_dataset(Path(tmp) / m) for m in members]
        merged = xr.merge(parts, compat="no_conflicts", combine_attrs="drop_conflicts")
        tmp_nc = path.with_suffix(".nc.merge")
        merged.to_netcdf(tmp_nc, encoding={v: {"zlib": True, "complevel": 4}
                                           for v in merged.data_vars})
        for p in parts:
            p.close()
        os.replace(tmp_nc, path)


def download_one(kind: str, cfg: dict, date: dt.date, hours: list[int],
                 area: list[float], target: Path, force: bool, dry_run: bool) -> str:
    """Returns 'downloaded' | 'merged' | 'skipped' | 'dry-run' | 'FAILED: ...'.

    Creates its own CDS client so calls are safe to run from worker threads.
    """
    status = "downloaded"
    want = set(hours)
    if target.exists() and not force:
        have = hours_in_file(target, date)
        if want <= have:
            return "skipped"
        hours = sorted(want | have)
        status = "merged"

    request = build_request(cfg, kind, date, hours, area)
    if dry_run:
        print(f"  request for {target.name}:")
        print("  " + json.dumps(request, indent=2).replace("\n", "\n  "))
        return "dry-run"

    import cdsapi
    try:
        client = cdsapi.Client(quiet=True, progress=False)
        target.parent.mkdir(parents=True, exist_ok=True)
        part = target.with_suffix(".nc.part")
        client.retrieve(DATASET_IDS[kind], request).download(str(part))
        normalize_download(part)
        os.replace(part, target)
    except Exception as exc:
        return f"FAILED: {exc}"
    return status


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--month", type=int, required=True)
    parser.add_argument("--days", nargs="+", required=True,
                        help="day numbers and/or A-B ranges, e.g. 1 2 5-7")
    parser.add_argument("--hours", type=int, nargs="+", default=None,
                        help="UTC hours 0-23 (default from config)")
    parser.add_argument("--area", type=float, nargs=4, default=None,
                        metavar=("N", "W", "S", "E"), help="override config area")
    parser.add_argument("--datasets", nargs="+", choices=["plev", "sfc"],
                        default=["plev", "sfc"])
    parser.add_argument("--config", default=None, help="path to config.yaml")
    parser.add_argument("--jobs", type=int, default=4,
                        help="concurrent CDS requests (default 4); queue waits "
                             "overlap so this speeds up multi-day downloads")
    parser.add_argument("--force", action="store_true", help="re-download even if complete")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the CDS requests without downloading")
    args = parser.parse_args(argv)

    cfg = load_config(args.config, source="era5")  # this downloader is ERA5-only
    hours = args.hours if args.hours is not None else cfg["download"]["default_hours"]
    bad = [h for h in hours if not 0 <= h <= 23]
    if bad:
        sys.exit(f"invalid hour(s): {bad}")
    area = list(args.area) if args.area is not None else source_area(cfg)
    days = parse_days(args.days, args.year, args.month)

    if not args.dry_run:
        require_cds_credentials()

    path_fn = {"plev": plev_path, "sfc": sfc_path}
    tasks = [(dt.date(args.year, args.month, day), kind)
             for day in days for kind in args.datasets]

    if args.dry_run:
        for date, kind in tasks:
            print(f"{date} {kind}:")
            download_one(kind, cfg, date, hours, area, path_fn[kind](cfg, date),
                         args.force, True)
        return 0

    failures = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {pool.submit(download_one, kind, cfg, date, hours, area,
                               path_fn[kind](cfg, date), args.force, False): (date, kind)
                   for date, kind in tasks}
        for fut in concurrent.futures.as_completed(futures):
            date, kind = futures[fut]
            result = fut.result()
            print(f"{date} {kind}: {result}", flush=True)
            if result.startswith("FAILED"):
                failures += 1
    done = len(tasks) - failures
    print(f"{done}/{len(tasks)} requests completed" +
          (f", {failures} FAILED (rerun the same command to retry)" if failures else ""))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
