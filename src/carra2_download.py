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


def build_request(cfg: dict, kind: str, dates: list[dt.date], hours: list[int],
                  area: list[float]) -> dict:
    """One CDS request covering every date in ``dates``.

    The CDS expands year/month/day as a cartesian product, so a request must
    not span months — ``main`` only ever chunks within the month it was asked
    for, which keeps the product exact.
    """
    months = {(d.year, d.month) for d in dates}
    if len(months) != 1:
        raise ValueError(f"a request must stay within one month, got {sorted(months)}")
    c2 = cfg["carra2"]
    request = {
        "level_type": LEVEL_TYPES[kind],
        "product_type": "analysis",
        "variable": c2["plev_variables"] if kind == "plev" else c2["sfc_variables"],
        "year": [f"{dates[0].year:04d}"],
        "month": [f"{dates[0].month:02d}"],
        "day": [f"{d.day:02d}" for d in dates],
        "time": [f"{h:02d}:00" for h in hours],
        "data_format": "netcdf",
        "area": area,
    }
    if kind == "plev":
        request["level_location"] = [str(p) for p in c2["pressure_levels"]]
    return request


def chunked(items: list, size: int) -> list[list]:
    """Split into consecutive chunks of at most ``size`` items."""
    size = max(1, int(size))
    return [items[i:i + size] for i in range(0, len(items), size)]


def chunk_days_for(cfg: dict, kind: str, override: int | None = None) -> int:
    """Days per request for one level type.

    ``carra2.chunk_days`` may be a single number or a per-level-type mapping,
    because the two differ by more than an order of magnitude in CDS "cost".
    The archive rejects a request above a cost limit of 12000; at the default
    variables/levels/area a profile day costs 960 and a surface day 72, so
    plev tops out at 12 days per request and sfc comfortably takes a month.
    Changing ``plev_variables``, ``pressure_levels`` or ``area`` changes the
    cost proportionally — the rejection is immediate and says so.
    """
    if override is not None:
        return int(override)
    value = cfg["carra2"].get("chunk_days", 12)
    if isinstance(value, dict):
        return int(value.get(kind, 12))
    return int(value)


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


def _split_static(ds):
    """Separate the time-varying variables from the grid description.

    ``domain_mask`` and the grid-mapping variable have no time dimension, so
    they must be held aside when concatenating slabs along ``valid_time``.
    """
    static = [v for v in ds.data_vars if "valid_time" not in ds[v].dims]
    return ds.drop_vars(static), ds[static]


def _merge_existing(out, target: Path):
    """Union of a fresh slab with the hours already in ``target``.

    Existing hours are preserved rather than replaced, so re-running for a
    different hour set never silently drops previously downloaded times. If
    the file on disk holds a different variable set (e.g. it predates a
    change to ``carra2.plev_variables``) it is replaced outright, since
    concatenating mismatched variables is not meaningful.
    """
    import numpy as np
    import xarray as xr

    if not target.exists():
        return out, "downloaded"
    try:
        with xr.open_dataset(target) as handle:
            old = handle.load()
    except Exception:
        return out, "downloaded"          # unreadable: replace it
    if set(old.data_vars) != set(out.data_vars):
        return out, "replaced (variable set changed)"
    # The horizontal grid must match exactly. If it does not — an older file
    # written before the projection constants were corrected, say — xarray
    # would outer-join the two coordinate sets and silently fill the result
    # with NaN, so replace the file instead of concatenating.
    for coord in ("x", "y"):
        if coord in old.coords and coord in out.coords:
            a, b = old[coord].values, out[coord].values
            if a.shape != b.shape or not np.allclose(a, b):
                return out, "replaced (grid changed)"

    old_t, _ = _split_static(old)
    new_t, static = _split_static(out)
    # join="exact": any residual coordinate mismatch is an error, never a
    # silently NaN-padded union
    combined = xr.concat([old_t, new_t], dim="valid_time", join="exact")
    _, keep = np.unique(combined["valid_time"].values, return_index=True)
    combined = combined.isel(valid_time=np.sort(keep)).sortby("valid_time")
    # static vars were dropped from `combined`, so there is nothing to reconcile
    merged = combined.merge(static, compat="override")
    merged.attrs = out.attrs
    status = "merged" if merged.sizes["valid_time"] > out.sizes["valid_time"] else "downloaded"
    return merged, status


def write_normalized_chunk(raw: Path, kind: str, cfg: dict, area: list[float],
                           dates: list[dt.date]) -> dict[dt.date, str]:
    """Normalize a multi-day delivery and split it into one file per day."""
    import numpy as np
    import pandas as pd

    from reanlib.io_carra2 import _standardize_names

    ds = _standardize_names(open_delivery(raw))
    path_fn = plev_path if kind == "plev" else sfc_path
    stamps = pd.DatetimeIndex(ds["valid_time"].values)
    results: dict[dt.date, str] = {}
    try:
        for date in dates:
            idx = np.flatnonzero(stamps.date == date)
            if idx.size == 0:
                results[date] = "FAILED: no time steps for this day in the delivery"
                continue
            # one day at a time: a whole plev chunk is ~63 GB decompressed on
            # the full CARRA-2 mesh, so it must never be materialized at once
            day = normalize_carra2(ds.isel(valid_time=idx), kind, cfg, area=area)
            target = path_fn(cfg, date)
            target.parent.mkdir(parents=True, exist_ok=True)
            out, status = _merge_existing(day.load(), target)
            tmp_nc = target.with_suffix(".nc.norm")
            out.to_netcdf(tmp_nc, encoding={v: {"zlib": True, "complevel": 4}
                                            for v in out.data_vars})
            os.replace(tmp_nc, target)
            results[date] = status
    finally:
        ds.close()
    return results


def download_chunk(kind: str, cfg: dict, dates: list[dt.date], hours: list[int],
                   area: list[float], force: bool, dry_run: bool) -> dict[dt.date, str]:
    """Fetch one chunk of days in a single CDS request, then split it per day.

    Batching matters far more than request size when the archive is busy: the
    CDS queues per request, so a month costs 62 queue positions one day at a
    time and 2 as whole-month chunks. Returns a per-date status.

    Runs in its own process (see ``main``), so the CDS client it builds shares
    no socket state with any other worker.
    """
    path_fn = plev_path if kind == "plev" else sfc_path
    want = set(hours)
    todo = [d for d in dates
            if force or not (want <= hours_in_file(path_fn(cfg, d), d))]
    skipped = {d: "skipped" for d in dates if d not in todo}
    if not todo:
        return skipped

    request = build_request(cfg, kind, todo, hours, area)
    if dry_run:
        span = f"{todo[0]:%Y-%m-%d}..{todo[-1]:%Y-%m-%d}" if len(todo) > 1 else f"{todo[0]}"
        print(f"  request for {kind} {span} ({len(todo)} day(s)):")
        print("  " + json.dumps(request, indent=2).replace("\n", "\n  "))
        return {**skipped, **{d: "dry-run" for d in todo}}

    import cdsapi

    raw = path_fn(cfg, todo[0]).with_suffix(f".{kind}-chunk.part")
    try:
        client = cdsapi.Client(quiet=True, progress=False)
        raw.parent.mkdir(parents=True, exist_ok=True)
        client.retrieve(cfg["carra2"]["dataset"], request).download(str(raw))
        results = write_normalized_chunk(raw, kind, cfg, area, todo)
        raw.unlink(missing_ok=True)
    except Exception as exc:
        # the raw delivery is deliberately left in place: a normalization bug
        # is then fixable offline instead of costing another queue position
        hint = ""
        if "cost limits exceeded" in str(exc) or "too large" in str(exc):
            hint = (f" -- lower carra2.chunk_days (or --chunk-days) below "
                    f"{len(todo)} for {kind}, or trim variables/levels")
        return {**skipped, **{d: f"FAILED: {exc}{hint}" for d in todo}}
    return {**skipped, **results}


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
    parser.add_argument("--chunk-days", type=int, default=None,
                        help="days bundled into ONE CDS request (default from "
                             "config). The archive queues per request, not per "
                             "byte, so bundling is what shortens a month: 31 "
                             "gives 2 queue positions for January, 1 gives 62")
    parser.add_argument("--jobs", type=int, default=4,
                        help="concurrent CDS requests (default 4); queue waits "
                             "overlap so this speeds up multi-chunk downloads")
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

    dates = [dt.date(args.year, args.month, day) for day in days]
    tasks, plan = [], []
    for kind in args.datasets:
        n = chunk_days_for(cfg, kind, args.chunk_days)
        groups = chunked(dates, n)
        tasks += [(chunk, kind) for chunk in groups]
        plan.append(f"{kind} {len(groups)} x <={n}d")
    print(f"{len(dates)} day(s) -> {len(tasks)} CDS request(s)  "
          f"({', '.join(plan)})", flush=True)

    if args.dry_run:
        for chunk, kind in tasks:
            download_chunk(kind, cfg, chunk, hours, area, args.force, True)
        return 0

    # separate processes, not threads: the CDS client keeps socket state that
    # has been seen to break when several of them poll concurrently in one
    # interpreter ([Errno 9] Bad file descriptor on every status check)
    failures = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {pool.submit(download_chunk, kind, cfg, chunk, hours, area,
                               args.force, False): (chunk, kind)
                   for chunk, kind in tasks}
        for fut in concurrent.futures.as_completed(futures):
            chunk, kind = futures[fut]
            try:
                results = fut.result()
            except Exception as exc:                      # worker died outright
                results = {d: f"FAILED: {exc}" for d in chunk}
            for date in chunk:
                result = results.get(date, "FAILED: no result returned")
                print(f"{date} {kind}: {result}", flush=True)
                if result.startswith("FAILED"):
                    failures += 1
    total = len(dates) * len(args.datasets)
    print(f"{total - failures}/{total} day-files written" +
          (f", {failures} FAILED (rerun the same command to retry)" if failures else ""))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
