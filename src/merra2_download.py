#!/usr/bin/env python
"""Download MERRA-2 daily Arctic subsets from NASA GES DISC.

Collections (config `merra2:` block): M2I3NPASM (3-hourly INSTANTANEOUS
pressure-level state, 42 levels) and M2I1NXASM (hourly INSTANTANEOUS
single-level state) — both snapshots, so surface and profile describe the
same instant. Granules are found with earthaccess (Earthdata login in
~/.netrc) and subset server-side through the GES DISC OPeNDAP endpoints, so
only the configured area (default 80-90N) is transferred instead of the
~1.5 GB/day global files. Subsets are normalized on write to the ERA5 file
conventions (reanlib/io_merra2.normalize_merra2): variables t/q/o3/clwc/ciwc
and t2m/skt/sp, coords valid_time/latitude/longitude/pressure_level — every
downstream stage then reads them with --source merra2, unchanged.

Files land in data/merra2/YYYY/MM/DD/merra2_{plev,sfc}_YYYYMMDD.nc.
Re-running skips complete files and re-fetches the union of hours if a file
is missing some.

Examples:
    python src/merra2_download.py --year 2020 --month 1 --days 1
    python src/merra2_download.py --year 2020 --month 1 --days 1-31 --jobs 4
    python src/merra2_download.py --year 2020 --month 1 --days 1 --dry-run
    python src/merra2_download.py --year 2020 --month 1 --days 1 --full-granule
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from era5_download import parse_days
from reanlib.config import load_config, plev_path, sfc_path
from reanlib.io_merra2 import hours_in_file, normalize_merra2, require_earthdata_credentials


def find_granule(kind: str, cfg: dict, date: dt.date):
    """Locate the day's granule via CMR (handles MERRA-2 stream numbering)."""
    import earthaccess

    m2 = cfg["merra2"]
    results = earthaccess.search_data(
        short_name=m2["collections"][kind], version=m2["version"],
        temporal=(date.isoformat(), (date + dt.timedelta(days=1)).isoformat()))
    tag = f".{date:%Y%m%d}."
    matches = [g for g in results if tag in g["umm"]["GranuleUR"]]
    if not matches:
        raise FileNotFoundError(
            f"no {m2['collections'][kind]} granule found for {date} "
            f"(GES DISC latency is a few weeks; check the date)")
    return matches[0]


def opendap_url(granule) -> str | None:
    """The granule's OPeNDAP endpoint, from its related URLs if advertised."""
    for entry in granule["umm"].get("RelatedUrls", []):
        url = entry.get("URL", "")
        if "opendap" in url.lower() and (entry.get("Subtype", "") == "OPENDAP DATA"
                                         or entry.get("Type", "") == "USE SERVICE API"):
            return url
    for url in granule.data_links(access="external"):
        if "gesdisc.eosdis.nasa.gov/data/" in url:  # goldsmr4/5 archive -> Hyrax
            return url.replace("/data/", "/opendap/")
    return None


def subset_via_opendap(url: str, kind: str, cfg: dict, date: dt.date,
                       hours: list[int]):
    """Open the granule over OPeNDAP and pull only the wanted slab."""
    import earthaccess
    import xarray as xr

    store = xr.backends.PydapDataStore.open(
        url, session=earthaccess.get_requests_https_session())
    ds = xr.open_dataset(store)
    return subset_local(ds, kind, cfg, date, hours)


def subset_local(ds, kind: str, cfg: dict, date: dt.date, hours: list[int]):
    north, west, south, east = cfg["area"]
    m2 = cfg["merra2"]
    variables = m2["plev_variables"] if kind == "plev" else m2["sfc_variables"]
    sel = ds[list(variables)].sel(lat=slice(south, north), lon=slice(west, east))
    keep = [t for t in sel["time"].values
            if t.astype("datetime64[s]").item().hour in set(hours)]
    return sel.sel(time=keep).load()


def download_one(kind: str, cfg: dict, date: dt.date, hours: list[int],
                 target: Path, force: bool, dry_run: bool,
                 full_granule: bool) -> str:
    """Returns 'downloaded' | 'merged' | 'skipped' | 'dry-run' | 'FAILED: ...'."""
    import earthaccess

    m2 = cfg["merra2"]
    status = "downloaded"
    want = set(hours)
    if target.exists() and not force:
        have = hours_in_file(target, date)
        if want <= have:
            return "skipped"
        want |= have
        status = "merged"
    hours = sorted(want)

    if dry_run:
        print(f"  {m2['collections'][kind]} v{m2['version']} {date} "
              f"hours {hours} area {cfg['area']} -> {target}")
        return "dry-run"

    try:
        granule = find_granule(kind, cfg, date)
        subset = None
        if not full_granule:
            url = opendap_url(granule)
            if url is not None:
                try:
                    subset = subset_via_opendap(url, kind, cfg, date, hours)
                except Exception as exc:
                    print(f"  {date} {kind}: OPeNDAP subset failed ({exc}); "
                          f"falling back to full granule", flush=True)
        if subset is None:
            import xarray as xr
            with tempfile.TemporaryDirectory() as tmp:
                files = earthaccess.download([granule], tmp, threads=1)
                with xr.open_dataset(files[0]) as full:
                    subset = subset_local(full, kind, cfg, date, hours)

        out = normalize_merra2(subset, kind)
        target.parent.mkdir(parents=True, exist_ok=True)
        part = target.with_suffix(".nc.part")
        out.to_netcdf(part, encoding={v: {"zlib": True, "complevel": 4}
                                      for v in out.data_vars})
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
                        help="UTC hours, multiples of 3 (default from config)")
    parser.add_argument("--datasets", nargs="+", choices=["plev", "sfc"],
                        default=["plev", "sfc"])
    parser.add_argument("--config", default=None, help="path to config.yaml")
    parser.add_argument("--jobs", type=int, default=4,
                        help="concurrent downloads (default 4)")
    parser.add_argument("--force", action="store_true",
                        help="re-download even if complete")
    parser.add_argument("--full-granule", action="store_true",
                        help="download whole granules and subset locally "
                             "instead of OPeNDAP server-side subsetting")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the planned requests without downloading")
    args = parser.parse_args(argv)

    cfg = load_config(args.config, source="merra2")  # this downloader is MERRA-2-only
    hours = args.hours if args.hours is not None else cfg["merra2"]["default_hours"]
    bad = [h for h in hours if not 0 <= h <= 23 or h % 3]
    if bad:
        sys.exit(f"invalid hour(s) {bad}: the M2I3NPASM pressure-level "
                 f"collection is 3-hourly (00, 03, ..., 21 UTC)")
    days = parse_days(args.days, args.year, args.month)

    path_fn = {"plev": plev_path, "sfc": sfc_path}
    tasks = [(dt.date(args.year, args.month, day), kind)
             for day in days for kind in args.datasets]

    if args.dry_run:
        for date, kind in tasks:
            print(f"{date} {kind}:")
            download_one(kind, cfg, date, hours, path_fn[kind](cfg, date),
                         args.force, True, args.full_granule)
        return 0

    strategy = require_earthdata_credentials()
    import earthaccess
    auth = earthaccess.login(strategy=strategy)
    if not auth.authenticated:
        sys.exit("Earthdata login failed — check the urs.earthdata.nasa.gov "
                 "entry in ~/.netrc (or EARTHDATA_USERNAME/PASSWORD)")

    failures = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {pool.submit(download_one, kind, cfg, date, hours,
                               path_fn[kind](cfg, date), args.force, False,
                               args.full_granule): (date, kind)
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
