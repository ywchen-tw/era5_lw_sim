"""CARRA-2 (Copernicus pan-Arctic Regional Reanalysis) support.

CARRA-2 granules are normalized *on write* by src/carra2_download.py using
``normalize_carra2``: variables are renamed to the same conventions as the
ERA5 files, so every analysis stage reads the same names from either source.

Unlike ERA5 and MERRA-2, CARRA-2 is a **regional** reanalysis on a 2.5 km
north polar-stereographic grid, so the horizontal dims stay ``y``/``x`` with
2-D ``latitude(y, x)`` / ``longitude(y, x)`` coordinates. Stages handle both
layouts through ``reanlib.grid`` rather than assuming 1-D lat/lon.

Two further differences from the global sources are handled here:

* **No specific humidity on pressure levels.** CARRA-2 publishes relative
  humidity only, so ``q`` is derived with ``specific_humidity_from_rh``.
  CARRA documents its relative humidity against saturation over *water*
  (config ``carra2.rh_over``).
* **The profile top is 50 hPa**, against 1 hPa for ERA5 and 0.1 hPa for
  MERRA-2. That is far above the inversion scan's 500 hPa limit and so does
  not affect stages 2-6, but any radiative-transfer use would have to splice
  a standard atmosphere from 50 hPa upward instead of from 1 hPa.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

#: Ratio of the gas constants for dry air and water vapour.
EPSILON = 0.621981

# CARRA-2's grid, recovered from a delivered file rather than assumed: these
# are the only parameters for which projecting the 2-D latitude/longitude
# yields a *regular* mesh, and they reproduce the documented 2.5 km spacing
# exactly (dx = dy = 2500.000 m, residual non-regularity < 0.2 m, against
# 320 km for a straight-up-Greenwich guess). The CDS delivery carries no x/y
# coordinates of its own, so they are reconstructed with these.
PROJ_CENTRAL_LON = -30.0        # straight vertical longitude from the pole
PROJ_STANDARD_PARALLEL = 90.0
PROJ_EARTH_RADIUS = 6371229.0   # m
PROJ_GRID_SPACING = 2500.0      # m, for the regularity check

# Variables recognized in each delivery, in pipeline (ERA5-convention)
# spelling. CARRA-2's GRIB short names already match ERA5's for these fields;
# the CDS may instead use the long CDS names, which VAR_ALIASES maps. 'r'
# (relative humidity) is converted to 'q' during normalization, not renamed.
PLEV_VARS = ("t", "r", "clwc", "ciwc", "cc")
SFC_VARS = ("t2m", "skt", "sp", "lsm")
# Everything the inversion metrics actually need is mandatory; the rest is
# accepted when present and skipped when the config does not request it, so
# `carra2.plev_variables` can be trimmed to shrink a download without
# breaking normalization. Stages 1-6 read only t and q (from r), plus
# t2m/skt/sp -- the cloud fields and lsm are for future radiative-transfer work.
OPTIONAL_VARS = frozenset({"clwc", "ciwc", "cc", "lsm"})
KIND_VARS = {"plev": PLEV_VARS, "sfc": SFC_VARS}

# CDS deliveries have used several coordinate spellings; map them all.
COORD_ALIASES = {
    "time": "valid_time", "forecast_reference_time": "valid_time",
    "isobaricInhPa": "pressure_level", "level": "pressure_level",
    "plev": "pressure_level", "lat": "latitude", "lon": "longitude",
    "projection_x_coordinate": "x", "projection_y_coordinate": "y",
}

# Variables the CDS may name differently from the GRIB short name.
VAR_ALIASES = {
    "2t": "t2m", "t2m": "t2m", "air_temperature": "t",
    "skt": "skt", "skin_temperature": "skt",
    "sp": "sp", "surface_pressure": "sp",
    "r": "r", "relative_humidity": "r",
    "t": "t", "temperature": "t",
    "clwc": "clwc", "specific_cloud_liquid_water_content": "clwc",
    "ciwc": "ciwc", "specific_cloud_ice_water_content": "ciwc",
    "cc": "cc", "cloud_cover": "cc",
    "lsm": "lsm", "land_sea_mask": "lsm",
}


def saturation_vapour_pressure(t_k, over: str = "water"):
    """Saturation vapour pressure [Pa] from temperature [K].

    Alduchov & Eskridge (1996) improved Magnus coefficients, which stay
    accurate to the Arctic winter temperatures this pipeline works at
    (better than 0.4 % down to -80 C over ice).

    ``over``: 'water', 'ice', or 'mixed' — the IFS-style blend that uses ice
    below -23 C, water above 0 C, and a quadratic ramp between.

    Only true ufuncs are used, so numpy arrays and xarray DataArrays both work
    and DataArray inputs keep their dims (the caller broadcasts a level
    coordinate against a 4-D field).
    """
    tc = t_k - 273.15
    e_water = 610.94 * np.exp(17.625 * tc / (tc + 243.04))
    e_ice = 611.21 * np.exp(22.587 * tc / (tc + 273.86))
    if over == "water":
        return e_water
    if over == "ice":
        return e_ice
    if over == "mixed":
        ramp = (t_k - 250.16) / (273.16 - 250.16)
        alpha = np.minimum(np.maximum(ramp, 0.0), 1.0) ** 2
        return alpha * e_water + (1.0 - alpha) * e_ice
    raise ValueError(f"unknown saturation reference {over!r} "
                     "(expected 'water', 'ice' or 'mixed')")


def specific_humidity_from_rh(rh_percent, t_k, p_pa, over: str = "water"):
    """Specific humidity [kg/kg] from relative humidity [%], T [K], p [Pa].

    Relative humidity is taken as the vapour-pressure ratio e/e_sat (the GRIB
    and ECMWF convention), so e = RH * e_sat and
    q = eps*e / (p - (1 - eps)*e).
    """
    e = np.maximum(rh_percent, 0.0) / 100.0 * saturation_vapour_pressure(t_k, over)
    e = np.minimum(e, 0.99 * p_pa)
    return EPSILON * e / (p_pa - (1.0 - EPSILON) * e)


#: Scalar GRIB coordinates that carry no information once the level type is
#: known (analysis step is always 0; the rest just name the level surface).
DROP_SCALAR_COORDS = ("step", "surface", "heightAboveGround", "number",
                      "depthBelowLandLayer", "entireAtmosphere")


def _standardize_names(ds: xr.Dataset) -> xr.Dataset:
    """Rename delivered coordinates and variables to pipeline spellings.

    Idempotent, so it can run again on an already-standardized dataset.

    The CDS delivers CARRA-2 through cfgrib, which gives a ``time`` dimension
    *and* a ``valid_time`` coordinate along it (identical for analyses, where
    step = 0). Renaming ``time`` to ``valid_time`` would collide with that
    existing coordinate, so the dimension is swapped onto it instead.
    """
    if "valid_time" in ds.coords and "time" in ds.dims:
        ds = ds.swap_dims({"time": "valid_time"})
        if "time" in ds.coords:
            ds = ds.drop_vars("time")
    taken = set(ds.dims) | set(ds.coords)
    renames = {k: v for k, v in COORD_ALIASES.items()
               if (k in ds.dims or k in ds.coords) and v not in taken}
    renames.update({k: v for k, v in VAR_ALIASES.items()
                    if k in ds.data_vars and k != v and v not in ds.data_vars})
    ds = ds.rename(renames) if renames else ds
    drop = [c for c in DROP_SCALAR_COORDS if c in ds.coords]
    return ds.drop_vars(drop) if drop else ds


def _projection_coords(ds: xr.Dataset) -> xr.Dataset:
    """Ensure 1-D projection coordinates x/y plus a CF grid-mapping variable.

    The CDS netCDF converter is documented as experimental and has not always
    carried x/y; when they are missing they are reconstructed by projecting
    the 2-D latitude/longitude onto the CARRA-2 polar stereographic.
    """
    from .grid import GRID_MAPPING_VAR

    if "x" not in ds.coords or "y" not in ds.coords:
        import pyproj

        proj = pyproj.Proj(proj="stere", lat_0=90.0,
                           lat_ts=PROJ_STANDARD_PARALLEL,
                           lon_0=PROJ_CENTRAL_LON, R=PROJ_EARTH_RADIUS)
        lat = ds["latitude"].transpose("y", "x").values
        lon = ds["longitude"].transpose("y", "x").values
        xx, yy = proj(lon, lat)
        # The grid must come out regular: x depending only on the column and y
        # only on the row. If it does not, the projection constants above no
        # longer match the delivery and the x/y written here — and every map
        # drawn from them — would be silently wrong, so fail loudly instead.
        spread = max(np.nanmax(np.nanstd(xx, axis=0)),
                     np.nanmax(np.nanstd(yy, axis=1)))
        if spread > 0.05 * PROJ_GRID_SPACING:
            raise ValueError(
                f"CARRA-2 grid is not regular under lat_ts="
                f"{PROJ_STANDARD_PARALLEL}, lon_0={PROJ_CENTRAL_LON} "
                f"(spread {spread:.1f} m vs {PROJ_GRID_SPACING:.0f} m spacing) "
                "— the delivered projection differs from the one assumed in "
                "reanlib/io_carra2.py")
        # regular by the check above, so any row/column gives the coordinates
        ds = ds.assign_coords(x=("x", xx[xx.shape[0] // 2, :]),
                              y=("y", yy[:, yy.shape[1] // 2]))
    ds["x"].attrs.update(units="m", standard_name="projection_x_coordinate")
    ds["y"].attrs.update(units="m", standard_name="projection_y_coordinate")

    ds[GRID_MAPPING_VAR] = xr.DataArray(
        np.int32(0),
        attrs={
            "grid_mapping_name": "polar_stereographic",
            "straight_vertical_longitude_from_pole": PROJ_CENTRAL_LON,
            "standard_parallel": PROJ_STANDARD_PARALLEL,
            "latitude_of_projection_origin": 90.0,
            "earth_radius": PROJ_EARTH_RADIUS,
        },
    )
    for name in ds.data_vars:
        if name != GRID_MAPPING_VAR and {"y", "x"} <= set(ds[name].dims):
            ds[name].attrs["grid_mapping"] = GRID_MAPPING_VAR
    return ds


def normalize_carra2(ds: xr.Dataset, kind: str, cfg: dict | None = None,
                     *, area: list[float] | None = None) -> xr.Dataset:
    """Rename a CARRA-2 delivery to pipeline conventions (kind: plev/sfc).

    Returns a dataset satisfying the same guarantees as io_era5.open_era5,
    except that the horizontal dims are ``y``/``x`` with 2-D latitude and
    longitude coordinates:

      - time dim 'valid_time'; vertical dim (plev only) 'pressure_level' in
        hPa sorted descending (index 0 = 1000 hPa)
      - plev files carry 'q', derived from the delivered relative humidity
      - a CF 'polar_stereographic' grid-mapping variable and x/y coordinates

    ``area`` [N, W, S, E] crops to the cells whose centres fall inside it.
    The CDS masks rather than crops: it returns the entire 2869x2869
    pan-Arctic mesh carrying values only inside the requested area, so this
    discards empty canvas, never data. It happens before any arithmetic —
    see the comment at the call site.
    """
    rh_over = (cfg or {}).get("carra2", {}).get("rh_over", "water")
    ds = _standardize_names(ds)

    wanted = KIND_VARS[kind]
    missing = [v for v in wanted if v not in ds and v not in OPTIONAL_VARS]
    if missing:
        raise KeyError(f"CARRA-2 {kind} delivery is missing variable(s) {missing}; "
                       f"got {sorted(ds.data_vars)}")
    ds = ds[[v for v in wanted if v in ds]]

    for coord in ("latitude", "longitude"):
        if coord not in ds.coords:
            raise KeyError(f"CARRA-2 delivery has no {coord} coordinate; "
                           f"got {sorted(ds.coords)}")
    if ds["latitude"].ndim != 2:
        raise ValueError("expected 2-D latitude/longitude on the CARRA-2 "
                         f"projected grid, got {ds['latitude'].ndim}-D")

    # Clip FIRST, before any arithmetic. The CDS masks rather than crops: it
    # returns the whole 2869x2869 pan-Arctic mesh with values only inside the
    # requested area (~2 % of cells for 85-90N). Deriving q on the full mesh
    # would materialize ~63 GB for a 12-day chunk; on the clipped box it is
    # ~1.5 GB. isel keeps the selection lazy, so nothing large is ever read.
    if area is not None:
        ds = clip_to_area(ds, area)

    if "pressure_level" in ds.dims:
        ds["pressure_level"] = ds["pressure_level"].astype(float)
        ds = ds.sortby("pressure_level", ascending=False)
        # q from relative humidity: CARRA-2 has no specific humidity field
        p_pa = ds["pressure_level"] * 100.0
        ds["q"] = specific_humidity_from_rh(ds["r"], ds["t"], p_pa, rh_over)
        ds["q"].attrs = {
            "units": "kg kg**-1",
            "long_name": "specific humidity",
            "comment": (f"derived from CARRA-2 relative humidity (saturation over "
                        f"{rh_over}); CARRA-2 publishes no specific humidity on "
                        "pressure levels"),
        }
        ds = ds.drop_vars("r")
        # cloud cover is delivered in %, ERA5's cc is a 0-1 fraction
        if "cc" in ds and float(np.nanmax(ds["cc"].values)) > 1.5:
            ds["cc"] = ds["cc"] / 100.0
            ds["cc"].attrs["units"] = "1"

    ds = _projection_coords(ds)

    ds.attrs = {
        "source": ("CARRA-2 (Copernicus pan-Arctic Regional Reanalysis, "
                   "HARMONIE-AROME 2.5 km), via CDS"),
        "expver_values": "n/a (CARRA-2)",
        "grid": "north polar stereographic, 2.5 km",
        "rh_over": rh_over,
    }
    return ds


def clip_to_area(ds: xr.Dataset, area: list[float]) -> xr.Dataset:
    """Trim a projected grid to the y/x span whose cells fall in [N, W, S, E].

    The CDS subsets CARRA-2 in projection space, so the delivered box always
    over-covers the requested latitude band. Keeping a rectangle of y/x (the
    smallest one containing every in-band cell) preserves the regular
    projection coordinates; cells outside the band are masked instead of
    dropped so the array stays rectangular.
    """
    north, west, south, east = area
    lat = ds["latitude"].transpose("y", "x").values
    lon = ((ds["longitude"].transpose("y", "x").values + 180.0) % 360.0) - 180.0
    inside = (lat >= south) & (lat <= north)
    if west <= east and (west, east) != (-180.0, 180.0):
        inside &= (lon >= west) & (lon <= east)
    elif west > east:  # box crossing the date line
        inside &= (lon >= west) | (lon <= east)
    if not inside.any():
        raise ValueError(f"no CARRA-2 cells inside area {area}")

    rows = np.flatnonzero(inside.any(axis=1))
    cols = np.flatnonzero(inside.any(axis=0))
    ds = ds.isel(y=slice(int(rows[0]), int(rows[-1]) + 1),
                 x=slice(int(cols[0]), int(cols[-1]) + 1))
    keep = xr.DataArray(inside[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1],
                        dims=("y", "x"))
    from .grid import GRID_MAPPING_VAR

    for name in ds.data_vars:
        if name != GRID_MAPPING_VAR and {"y", "x"} <= set(ds[name].dims):
            ds[name] = ds[name].where(keep)
    ds["domain_mask"] = keep.astype(np.int8)
    ds["domain_mask"].attrs = {
        "units": "1",
        "long_name": "1 where the cell centre falls inside the requested area",
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
