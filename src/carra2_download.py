#!/usr/bin/env python
"""Download CARRA-2 pressure-level and single-level analyses from the CDS.

CARRA-2 (Copernicus pan-Arctic Regional Reanalysis, HARMONIE-AROME at 2.5 km)
gathers every level type in one catalogue entry, `reanalysis-pan-carra`, and
supports geographic subsetting — so one request per day per level type,
trimmed to the configured area, mirrors the ERA5 downloader. Deliveries are
normalized on write (reanlib/io_carra2.normalize_carra2) to the ERA5 variable
names, keeping the native polar-stereographic y/x grid with 2-D latitude and
longitude; downstream stages read it with --source carra2.

Note this is the *sub-daily* entry: `reanalysis-pan-carra-means` holds only
daily and monthly aggregates, which cannot feed the inversion metrics.

Analyses exist 3-hourly (00, 03, ... 21 UTC). CARRA-2 publishes no specific
humidity on pressure levels, so `q` is derived from relative humidity; and no
ozone at all, which is why stage 7 does not accept --source carra2.

Files land in data/carra2/YYYY/MM/DD/carra2_{plev,sfc}_YYYYMMDD.nc.
Re-running skips complete files and re-requests the union of hours if a file
is missing some.

Examples:
    python src/carra2_download.py --year 2020 --month 1 --days 1
    python src/carra2_download.py --year 2020 --month 1 --days 1-31 --jobs 4
    python src/carra2_download.py --year 2020 --month 1 --days 1 --dry-run
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import json
import os
import sys
import tempfile
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from era5_download import parse_days
from reanlib.config import load_config, plev_path, sfc_path, source_area
from reanlib.io_carra2 import hours_in_file, normalize_carra2
from reanlib.io_era5 import require_cds_credentials

LEVEL_TYPES = {"plev": "pressure_levels", "sfc": "single_levels"}

#: CARRA-2 analyses run 3-hourly; other hours exist only in the forecast stream.
ANALYSIS_STEP_H = 3


def build_request(cfg: dict, kind: str, date: dt.date, hours: list[int],
                  area: list[float]) -> dict:
    c2 = cfg["carra2"]
    request = {
        "level_type": LEVEL_TYPES[kind],
        "product_type": "analysis",
        "variable": c2["plev_variables"] if kind == "plev" else c2["sfc_variables"],
        "year": [f"{date.year:04d}"],
        "month": [f"{date.month:02d}"],
        "day": [f"{date.day:02d}"],
        "time": [f"{h:02d}:00" for h in hours],
        "data_format": "netcdf",
        "area": area,
    }
    if kind == "plev":
        request["level_location"] = [str(p) for p in c2["pressure_levels"]]
    return request


def open_delivery(path: Path):
    """Open a CDS delivery as one dataset, merging a zipped multi-part reply.

    The CDS returns a zip of one netCDF per GRIB stream when a request spans
    several, exactly as it does for ERA5.
    """
    import xarray as xr

    if not zipfile.is_zipfile(path):
        return xr.open_dataset(path)
    with tempfile.TemporaryDirectory(dir=path.parent) as tmp:
        with zipfile.ZipFile(path) as zf:
            members = zf.namelist()
            zf.extractall(tmp)
        parts = [xr.open_dataset(Path(tmp) / m) for m in members]
        merged = xr.merge(parts, compat="no_conflicts",
                          combine_attrs="drop_conflicts").load()
        for p in parts:
            p.close()
    return merged


def write_normalized(raw: Path, target: Path, kind: str, cfg: dict,
                     area: list[float]) -> None:
    """Normalize a raw delivery and write the daily file, replacing it atomically."""
    ds = open_delivery(raw)
    try:
        out = normalize_carra2(ds, kind, cfg, area=area)
        out = out.sortby("valid_time")
        tmp_nc = target.with_suffix(".nc.norm")
        encoding = {v: {"zlib": True, "complevel": 4} for v in out.data_vars}
        out.to_netcdf(tmp_nc, encoding=encoding)
    finally:
        ds.close()
    os.replace(tmp_nc, target)


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
        raw = target.with_suffix(".nc.part")
        client.retrieve(cfg["carra2"]["dataset"], request).download(str(raw))
        write_normalized(raw, target, kind, cfg, area)
        raw.unlink(missing_ok=True)
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
                        help=f"UTC hours, multiples of {ANALYSIS_STEP_H} "
                             "(default from config)")
    parser.add_argument("--area", type=float, nargs=4, default=None,
                        metavar=("N", "W", "S", "E"), help="override config area")
    parser.add_argument("--datasets", nargs="+", choices=["plev", "sfc"],
                        default=["plev", "sfc"])
    parser.add_argument("--config", default=None, help="path to config.yaml")
    parser.add_argument("--jobs", type=int, default=4,
                        help="concurrent CDS requests (default 4); queue waits "
                             "overlap so this speeds up multi-day downloads")
    parser.add_argument("--force", action="store_true",
                        help="re-download even if complete")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the CDS requests without downloading")
    args = parser.parse_args(argv)

    cfg = load_config(args.config, source="carra2")  # this downloader is CARRA-2-only
    hours = args.hours if args.hours is not None else cfg["carra2"]["default_hours"]
    bad = [h for h in hours if not 0 <= h <= 23]
    if bad:
        sys.exit(f"invalid hour(s): {bad}")
    off_cadence = [h for h in hours if h % ANALYSIS_STEP_H]
    if off_cadence:
        sys.exit(f"CARRA-2 analyses are {ANALYSIS_STEP_H}-hourly; "
                 f"hour(s) {off_cadence} are not available "
                 "(they exist only in the forecast stream)")
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
