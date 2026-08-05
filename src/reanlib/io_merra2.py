"""MERRA-2 (NASA GES DISC) support: credential checks and normalization.

MERRA-2 granules are normalized *on write* by src/merra2_download.py using
``normalize_merra2``: variables and coordinates are renamed to the same
conventions as the ERA5 files (see reanlib/io_era5.open_era5), so every
analysis stage opens either source's daily files unchanged.

Units already agree between the two sources: T/TS/T2M in K, PS in Pa,
QV/O3/QL/QI as kg/kg mass mixing ratios, lev in hPa. One behavioral
difference remains in the data: ERA5 extrapolates pressure levels below the
surface, while MERRA-2 stores fill values there (NaN after decoding); the
inversion scan treats non-finite levels as below-ground.
"""

from __future__ import annotations

import netrc
import os
import sys
from pathlib import Path

import numpy as np
import xarray as xr

# MERRA-2's canonical fill value is 1e15 (undefined, e.g. below-ground levels).
# OPeNDAP transport can drop the _FillValue attribute, so xarray's automatic
# mask-and-scale cannot be relied on; normalize_merra2 masks explicitly.
FILL_THRESHOLD = 1.0e14

EARTHDATA_HOST = "urs.earthdata.nasa.gov"

_CREDENTIAL_HELP = f"""\
NASA Earthdata credentials not found.

MERRA-2 is distributed by the NASA GES DISC and needs an Earthdata login.
Register (free) at https://urs.earthdata.nasa.gov, then EITHER add to
~/.netrc (recommended):

    machine {EARTHDATA_HOST} login <username> password <password>

and `chmod 600 ~/.netrc`, OR export environment variables (e.g. in ~/.zshrc):

    export EARTHDATA_USERNAME=<username>
    export EARTHDATA_PASSWORD=<password>

Also authorize the "NASA GESDISC DATA ARCHIVE" application once at
https://urs.earthdata.nasa.gov (Applications -> Authorized Apps) —
otherwise the first data request is rejected.
"""

# MERRA-2 short name -> pipeline (ERA5-convention) name
PLEV_RENAMES = {"T": "t", "QV": "q", "O3": "o3", "QL": "clwc", "QI": "ciwc"}
SFC_RENAMES = {"T2M": "t2m", "TS": "skt", "PS": "sp"}
COORD_RENAMES = {"time": "valid_time", "lat": "latitude", "lon": "longitude",
                 "lev": "pressure_level"}


def require_earthdata_credentials() -> str:
    """The earthaccess login strategy to use ('netrc' or 'environment');
    exits with setup instructions if neither is configured."""
    try:
        if netrc.netrc().authenticators(EARTHDATA_HOST):
            return "netrc"
    except (FileNotFoundError, netrc.NetrcParseError):
        pass
    if os.environ.get("EARTHDATA_USERNAME") and os.environ.get("EARTHDATA_PASSWORD"):
        return "environment"
    sys.exit(_CREDENTIAL_HELP)


def normalize_merra2(ds: xr.Dataset, kind: str) -> xr.Dataset:
    """Rename a MERRA-2 subset to pipeline conventions (kind: 'plev' or 'sfc').

    Returns a dataset satisfying the same guarantees as io_era5.open_era5:
    time dim 'valid_time'; 'latitude' descending; vertical dim (plev only)
    'pressure_level' in hPa sorted descending (index 0 = 1000 hPa).
    """
    var_renames = PLEV_RENAMES if kind == "plev" else SFC_RENAMES
    missing = [v for v in var_renames if v not in ds]
    if missing:
        raise KeyError(f"MERRA-2 {kind} granule is missing variable(s) {missing}")
    ds = ds[list(var_renames)]
    ds = ds.rename({**var_renames,
                    **{k: v for k, v in COORD_RENAMES.items() if k in ds.dims}})
    for name in ds.data_vars:
        ds[name] = ds[name].where(np.abs(ds[name]) < FILL_THRESHOLD)
    ds = ds.sortby("latitude", ascending=False)
    if "pressure_level" in ds.dims:
        ds["pressure_level"] = ds["pressure_level"].astype(float)
        ds = ds.sortby("pressure_level", ascending=False)
    ds.attrs = {
        "source": "MERRA-2 (GMAO, GEOS 5.12.4), NASA GES DISC",
        "expver_values": "n/a (MERRA-2)",
    }
    return ds


def hours_in_file(path: Path, date) -> set[int]:
    """Hours already present in a normalized daily file (empty if unreadable)."""
    if not Path(path).exists():
        return set()
    try:
        with xr.open_dataset(path) as ds:
            times = ds["valid_time"].values
    except Exception:
        return set()
    import pandas as pd

    idx = pd.DatetimeIndex(times)
    return {int(h) for d, h in zip(idx.date, idx.hour) if d == date}
