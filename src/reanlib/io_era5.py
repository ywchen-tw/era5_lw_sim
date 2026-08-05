"""Opening and normalizing ERA5 netCDF files from the CDS, plus credential checks."""

from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

import xarray as xr

CDSAPIRC = Path.home() / ".cdsapirc"

_CREDENTIAL_HELP = f"""\
CDS API credentials not found.

Create {CDSAPIRC} with exactly two lines:

    url: https://cds.climate.copernicus.eu/api
    key: <your-personal-access-token>

Your Personal Access Token is shown at https://cds.climate.copernicus.eu/profile
after logging in (see https://cds.climate.copernicus.eu/how-to-api).

Also make sure you have accepted the licence terms of the ERA5 datasets on the
CDS website (each dataset's download page has a "Terms of use" section) —
otherwise the first request will be rejected.
"""


def require_cds_credentials() -> None:
    """Exit with setup instructions unless CDS API credentials are configured."""
    if CDSAPIRC.exists():
        return
    if os.environ.get("CDSAPI_URL") and os.environ.get("CDSAPI_KEY"):
        return
    sys.exit(_CREDENTIAL_HELP)


def open_era5(path: str | Path) -> xr.Dataset:
    """Open a CDS-delivered ERA5 netCDF file and normalize its conventions.

    Handles both the current CDS output (dims valid_time/pressure_level, expver as
    a per-time string coordinate) and legacy deliveries (dims time/level, expver as
    an extra dimension splitting ERA5 and ERA5T into NaN-padded slabs).

    Guarantees on the returned dataset:
      - time dim is 'valid_time'; vertical dim (if any) is 'pressure_level' in hPa
      - pressure_level is sorted descending (index 0 = lowest level, e.g. 1000 hPa)
      - no 'expver' or 'number' dims/vars; expver values seen are recorded in
        ds.attrs['expver_values']
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"{path} does not exist — download it first (src/era5_download.py "
            "for ERA5, src/merra2_download.py for MERRA-2)"
        )
    ds = xr.open_dataset(path)

    renames = {}
    if "time" in ds.dims and "valid_time" not in ds.dims:
        renames["time"] = "valid_time"
    for legacy in ("level", "isobaricInhPa"):
        if legacy in ds.dims and "pressure_level" not in ds.dims:
            renames[legacy] = "pressure_level"
    if renames:
        ds = ds.rename(renames)

    expver_values: list[str] = []
    if "expver" in ds.dims:
        # Legacy ERA5/ERA5T split: two slabs, each NaN where the other has data.
        expver_values = [str(v) for v in ds["expver"].values]
        first = ds.isel(expver=0, drop=True)
        merged = first
        for i in range(1, ds.sizes["expver"]):
            merged = merged.combine_first(ds.isel(expver=i, drop=True))
        ds = merged
    elif "expver" in ds.coords or "expver" in ds.variables:
        vals = ds["expver"].values
        expver_values = sorted({str(v) for v in vals.ravel()}) if vals.ndim else [str(vals)]
        ds = ds.drop_vars("expver")
    if any(v.strip("0") == "5" for v in expver_values):
        warnings.warn(
            f"{path.name} contains ERA5T (expver 0005) data, which is preliminary "
            "and may be revised in the final ERA5 release.",
            stacklevel=2,
        )
    if expver_values:
        ds.attrs["expver_values"] = ",".join(expver_values)
    else:
        ds.attrs.setdefault("expver_values", "unknown")  # MERRA-2 files carry their own

    if "number" in ds.coords or "number" in ds.variables:
        ds = ds.drop_vars("number")

    if "pressure_level" in ds.dims:
        ds = ds.sortby("pressure_level", ascending=False)

    return ds
