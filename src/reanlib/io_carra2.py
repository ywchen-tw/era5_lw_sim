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
  humidity only, so ``q`` is derived with ``specific_humidity_from_rh``
  (config ``carra2.rh_over``, default water). NOTE the pressure-level r was
  shown empirically to follow the model's temperature-dependent saturation
  (ice-like for T below about -25 C, water-like above; see PLAN_TODO), while
  the CARRA docs define only the 2 m RH (over water) and are silent for
  pressure levels. With the over-water default, derived q above ~750 hPa is
  biased high by construction, up to the e_w/e_i factor (~1.5 at 300 hPa).
* **The profile top is 50 hPa**, against 1 hPa for ERA5 and 0.1 hPa for
  MERRA-2. That is far above the inversion scan's 500 hPa limit and so does
  not affect stages 2-6; stage 7 splices the standard atmosphere from 50 hPa
  upward instead of from 1 hPa.
* **Surface radiation is forecast-stream only** (accumulated J m-2 from each
  3-hourly cycle start). ``normalize_carra2_rad`` differences consecutive
  lead times into 1-h mean W m-2 windows stamped at the window END, matching
  the ERA5 convention that stage 7 reads.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from .humidity import specific_humidity_from_rh

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
SFC_VARS = ("t2m", "skt", "sp", "lsm", "siconc", "tcc")
CLOUD_VARS = ("clwc", "ciwc", "cc")   # stage-7 top-up of a trimmed plev file
RAD_VARS = ("str", "strd")            # forecast stream, accumulated J m-2
# Everything the inversion metrics actually need is mandatory; the rest is
# accepted when present and skipped when the config does not request it, so
# `carra2.plev_variables` can be trimmed to shrink a download without
# breaking normalization. Stages 1-6 read only t and q (from r), plus
# t2m/skt/sp -- the cloud fields and lsm are for radiative-transfer work
# (stage 7), which fetches them separately via the `cloud` kind.
OPTIONAL_VARS = frozenset({"clwc", "ciwc", "cc", "lsm", "siconc", "tcc"})
KIND_VARS = {"plev": PLEV_VARS, "sfc": SFC_VARS, "cloud": CLOUD_VARS}

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
    # CARRA-2's GRIB short name for per-level cloud cover is `ccl`, not
    # ERA5's `cc` (seen in the first real cloud delivery)
    "cc": "cc", "cloud_cover": "cc", "ccl": "cc",
    "lsm": "lsm", "land_sea_mask": "lsm",
    "str": "str", "surface_net_thermal_radiation": "str",
    "strd": "strd", "thermal_surface_radiation_downwards": "strd",
    # classification fields; short names guessed until the first real
    # delivery confirms them (cf. ccl above) — the unknown-variable guard in
    # normalize_carra2 turns a wrong guess into a loud failure, raw retained
    "siconc": "siconc", "ci": "siconc", "icec": "siconc",
    "sea_ice_area_fraction": "siconc",
    "tcc": "tcc", "tcdc": "tcc", "total_cloud_cover": "tcc",
}


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
    if missing or not any(v in ds for v in wanted):
        raise KeyError(f"CARRA-2 {kind} delivery is missing variable(s) "
                       f"{missing or list(wanted)}; got {sorted(ds.data_vars)}")
    # a delivered data variable that maps to nothing we recognize means a
    # GRIB short name VAR_ALIASES does not know (cloud cover arrived as `ccl`
    # once already) — dropping it silently would discard paid-for data, so
    # fail loudly; the raw delivery is retained by the caller on exception
    unknown = sorted(set(ds.data_vars) - set(wanted))
    if unknown:
        raise KeyError(f"CARRA-2 {kind} delivery contains unrecognized "
                       f"variable(s) {unknown} — add them to "
                       "io_carra2.VAR_ALIASES")
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

    # fraction fields CARRA-2 delivers in percent (cc and tcc observed at
    # 0-100; siconc has arrived as 0-1 but gets the same guard)
    for name in ("cc", "tcc", "siconc"):
        if name in ds and float(np.nanmax(ds[name].values)) > 1.5:
            ds[name] = ds[name] / 100.0
            ds[name].attrs["units"] = "1"

    if "pressure_level" in ds.dims:
        ds["pressure_level"] = ds["pressure_level"].astype(float)
        ds = ds.sortby("pressure_level", ascending=False)
        # q from relative humidity: CARRA-2 has no specific humidity field.
        # (the `cloud` kind carries no r, so nothing to derive there)
        if "r" in ds:
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

    ds = _projection_coords(ds)

    ds.attrs = {
        "source": ("CARRA-2 (Copernicus pan-Arctic Regional Reanalysis, "
                   "HARMONIE-AROME 2.5 km), via CDS"),
        "expver_values": "n/a (CARRA-2)",
        "grid": "north polar stereographic, 2.5 km",
        "rh_over": rh_over,
    }
    return ds


def normalize_carra2_rad(ds: xr.Dataset, cfg: dict | None = None,
                         *, area: list[float] | None = None) -> xr.Dataset:
    """Normalize a CARRA-2 forecast-radiation delivery to 1-h flux windows.

    CARRA-2 has no radiation in the analysis stream: ``str``/``strd`` exist
    only as forecasts, accumulated J m-2 from each cycle start. Consecutive
    lead times are differenced into per-window mean fluxes and stamped at the
    window END (``valid_time = cycle + lead``), which is ERA5's convention
    ("accumulation over the hour ending at valid_time") — so stage 7 can read
    either source the same way. Output variables are ``lwdn``/``lwup`` in
    W m-2 (already divided by the window length; LWup = LWdn - net), with
    ``cycle_time`` and ``lead_h`` kept as provenance coordinates.

    Duplicate window ends (two cycles covering the same hour) keep the
    shorter lead, which sits closer to its initializing analysis.
    """
    # the analysis-oriented time->valid_time swap in _standardize_names must
    # not run here: in the forecast layout valid_time is 2-D (time, step)
    renames = {k: v for k, v in VAR_ALIASES.items()
               if k in ds.data_vars and k != v and v not in ds.data_vars}
    renames.update({k: v for k, v in
                    {"forecast_reference_time": "time", "forecast_period": "step",
                     "lat": "latitude", "lon": "longitude"}.items()
                    if (k in ds.coords or k in ds.dims) and v not in ds.dims
                    and v not in ds.coords})
    ds = ds.rename(renames) if renames else ds
    missing = [v for v in RAD_VARS if v not in ds]
    if missing:
        raise KeyError(f"CARRA-2 rad delivery is missing variable(s) {missing}; "
                       f"got {sorted(ds.data_vars)}")
    ds = ds[list(RAD_VARS)]
    for dim in ("time", "step"):     # single-cycle/-lead deliveries drop dims
        if dim not in ds.dims:
            ds = ds.expand_dims(dim)

    # clip before any arithmetic, exactly as normalize_carra2 does
    if area is not None:
        ds = clip_to_area(ds, area)

    step = ds["step"].values
    step_h = (step / np.timedelta64(1, "h") if np.issubdtype(step.dtype, np.timedelta64)
              else step.astype(float))
    order = np.argsort(step_h)
    ds = ds.isel(step=order)
    step_h = step_h[order]
    if step_h[0] <= 0:
        raise ValueError(f"lead times must be positive, got {step_h}")
    window_h = np.diff(np.concatenate([[0.0], step_h]))

    cycles = ds["time"].values
    ny, nx = ds.sizes["y"], ds.sizes["x"]
    flux = {}
    for v in RAD_VARS:
        acc = ds[v].transpose("time", "step", "y", "x").values.astype(float)
        w = np.empty_like(acc)
        w[:, 0] = acc[:, 0] / (window_h[0] * 3600.0)
        if step_h.size > 1:
            w[:, 1:] = np.diff(acc, axis=1) / (window_h[1:, None, None] * 3600.0)
        flux[v] = w.reshape(-1, ny, nx)

    vt = (cycles[:, None] + (step_h[None, :] * 3600.0e9).astype("timedelta64[ns]")
          ).reshape(-1)
    cyc_flat = np.repeat(cycles, step_h.size)
    lead_flat = np.tile(step_h, cycles.size)
    win_flat = np.tile(window_h, cycles.size)

    # sort by window end; on ties (overlapping cycles) keep the shortest lead
    order = np.lexsort((lead_flat, vt))
    _, first = np.unique(vt[order], return_index=True)
    keep = order[first]

    lwdn = flux["strd"][keep]
    lwup = lwdn - flux["str"][keep]
    out = xr.Dataset(
        {
            "lwdn": (("valid_time", "y", "x"), lwdn),
            "lwup": (("valid_time", "y", "x"), lwup),
        },
        coords={
            "valid_time": vt[keep],
            "cycle_time": (("valid_time",), cyc_flat[keep]),
            "lead_h": (("valid_time",), lead_flat[keep]),
            "window_h": (("valid_time",), win_flat[keep]),
            **{c: ds[c] for c in ("latitude", "longitude", "y", "x")
               if c in ds.coords},
        },
    )
    if "domain_mask" in ds:
        out["domain_mask"] = ds["domain_mask"]
    for name, long in (("lwdn", "surface LW down"), ("lwup", "surface LW up")):
        out[name].attrs = {
            "units": "W m**-2",
            "long_name": f"{long}, mean over the window ending at valid_time",
            "comment": ("differenced CARRA-2 forecast accumulations (str/strd, "
                        "J m-2 from cycle start); LWup = LWdn - net"),
        }
    out = _projection_coords(out)
    out.attrs = {
        "source": ("CARRA-2 forecast-stream surface thermal radiation "
                   "(reanalysis-pan-carra, product_type forecast), via CDS"),
        "expver_values": "n/a (CARRA-2)",
        "grid": "north polar stereographic, 2.5 km",
        "note": ("1-h (window_h) mean fluxes stamped at the window END, ERA5 "
                 "convention; on duplicate window ends the shorter lead wins"),
    }
    return out


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
