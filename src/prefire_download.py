#!/usr/bin/env python
"""Download PREFIRE TIRS L1B radiance granules and SRF files (stage 8 input).

Run under the ``era5`` conda env (needs ``earthaccess``; Earthdata login is
read from ``~/.netrc``). Granules of ``PREFIRE_SATx_1B-RAD`` (R01, NASA ASDC)
that intersect the configured Arctic domain land in
``data/prefire/YYYY/MM/``; the TIRS spectral-response-function files (v13,
Zenodo record 16638853 — the version the R01 granules were calibrated with,
see their ``SRF_NEdR_version`` attribute) land in ``data/prefire/srf/``.
Both downloads are idempotent: existing files are skipped.

Examples:
    python src/prefire_download.py --year 2025 --month 1 --days 1 --sat 1
    python src/prefire_download.py --srf-only
"""

from __future__ import annotations

import argparse
import calendar
import datetime as dt
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reanlib.config import load_config

SRF_ZENODO_URL = ("https://zenodo.org/api/records/16638853/files/"
                  "PREFIRE_SRFs-v13_2024-09-15.zip/content")
SRF_ZIP_NAME = "PREFIRE_SRFs-v13_2024-09-15.zip"
SRF_FILES = {1: "PREFIRE_TIRS1_SRF_v13_2024-09-15.nc",
             2: "PREFIRE_TIRS2_SRF_v13_2024-09-15.nc"}


def parse_days(tokens: "list[str]", year: int, month: int) -> "list[int]":
    """'1 2 5-7' -> [1, 2, 5, 6, 7], validated against the month length."""
    ndays = calendar.monthrange(year, month)[1]
    days: "set[int]" = set()
    for tok in tokens:
        if "-" in tok:
            a, b = tok.split("-")
            days.update(range(int(a), int(b) + 1))
        else:
            days.add(int(tok))
    bad = [d for d in days if not 1 <= d <= ndays]
    if bad:
        sys.exit(f"invalid day(s) {bad} for {year}-{month:02d} "
                 f"(month has {ndays} days)")
    return sorted(days)


def prefire_dir(cfg: dict) -> Path:
    return Path(cfg["paths"]["data"]) / "prefire"


def srf_path(cfg: dict, sat: int) -> Path:
    return prefire_dir(cfg) / "srf" / SRF_FILES[sat]


def granule_dir(cfg: dict, year: int, month: int) -> Path:
    return prefire_dir(cfg) / f"{year:04d}" / f"{month:02d}"


def download_srf(cfg: dict) -> None:
    """Fetch + unzip the v13 SRF pair from Zenodo (open access, no auth)."""
    out = prefire_dir(cfg) / "srf"
    if all(srf_path(cfg, s).exists() for s in (1, 2)):
        print(f"SRF files present in {out} — skipping")
        return
    import urllib.request

    out.mkdir(parents=True, exist_ok=True)
    zip_file = out / SRF_ZIP_NAME
    if not zip_file.exists():
        print(f"downloading {SRF_ZIP_NAME} from Zenodo ...")
        urllib.request.urlretrieve(SRF_ZENODO_URL, zip_file)
    with zipfile.ZipFile(zip_file) as zf:
        for member in zf.namelist():
            name = Path(member).name
            if name in SRF_FILES.values():
                target = out / name
                target.write_bytes(zf.read(member))
                print(f"  extracted {target}")


def download_granules(cfg: dict, year: int, month: int, days: "list[int]",
                      sat: int) -> "list[Path]":
    import earthaccess

    auth = earthaccess.login(strategy="netrc")
    if not auth.authenticated:
        sys.exit("Earthdata login failed — check the urs.earthdata.nasa.gov "
                 "entry in ~/.netrc")

    north, west, south, east = cfg["area"]
    out = granule_dir(cfg, year, month)
    out.mkdir(parents=True, exist_ok=True)

    t0 = dt.datetime(year, month, min(days))
    t1 = dt.datetime(year, month, max(days)) + dt.timedelta(days=1)
    granules = earthaccess.search_data(
        short_name=f"PREFIRE_SAT{sat}_1B-RAD", version="R01",
        temporal=(t0.isoformat(), t1.isoformat()),
        bounding_box=(west, south, east, north))
    wanted_days = set(days)

    def overlaps_wanted(g) -> bool:
        rng = g["umm"]["TemporalExtent"]["RangeDateTime"]
        gs = dt.datetime.fromisoformat(rng["BeginningDateTime"].rstrip("Z"))
        ge = dt.datetime.fromisoformat(rng["EndingDateTime"].rstrip("Z"))
        d = gs.date()
        while d <= ge.date():
            if (d.year, d.month) == (year, month) and d.day in wanted_days:
                return True
            d += dt.timedelta(days=1)
        return False

    granules = [g for g in granules if overlaps_wanted(g)]
    missing = [g for g in granules
               if not (out / g["umm"]["GranuleUR"]).exists()]
    print(f"SAT{sat} {year}-{month:02d} days {days}: {len(granules)} granules "
          f"intersect the domain, {len(missing)} to download -> {out}")
    if missing:
        earthaccess.download(missing, str(out))
    return [out / g["umm"]["GranuleUR"] for g in granules]


def main(argv: "list[str] | None" = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--year", type=int)
    parser.add_argument("--month", type=int)
    parser.add_argument("--days", nargs="+", default=None,
                        help="days of month, e.g. 1 2 5-7")
    parser.add_argument("--sat", type=int, default=1, choices=[1, 2],
                        help="PREFIRE satellite / TIRS instrument (default 1)")
    parser.add_argument("--srf-only", action="store_true",
                        help="only fetch the SRF files")
    parser.add_argument("--config", default=None)
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    download_srf(cfg)
    if args.srf_only:
        return 0
    if args.year is None or args.month is None or args.days is None:
        sys.exit("--year, --month and --days are required unless --srf-only")
    days = parse_days(args.days, args.year, args.month)
    files = download_granules(cfg, args.year, args.month, days, args.sat)
    print(f"{len(files)} granule file(s) ready")
    return 0


if __name__ == "__main__":
    sys.exit(main())
