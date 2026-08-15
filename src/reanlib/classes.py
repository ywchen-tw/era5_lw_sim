"""Surface-type and sky classification for stratified inversion statistics.

Classes are computed per (time, cell) at the analysis instants and stored as
int8 variables in the daily inversion files, so every downstream stage
stratifies identically and re-runs are reproducible.

Surface type (``surface_class``), decided in this order:
  0 land          land fraction >= 0.5
  1 open_water    sea-ice fraction <  0.05
  2 marginal_ice  0.05 <= sea-ice fraction <= 0.95
  3 pack_ice      sea-ice fraction >  0.95
 -1 unknown       classification inputs missing/NaN at this cell

Sky (``sky_class``), from total cloud cover:
  0 clear         tcc <= 0.05
  1 partial       between the thresholds (radiatively ambiguous; reported
                  but usually excluded from clear/cloudy contrasts)
  2 cloudy        tcc >= 0.95
 -1 unknown

The clear threshold is 0.05 rather than stage 7's 0.01 because MERRA-2's
CLDTOT almost never reaches 0.01 (trace fractions everywhere); 0.05 keeps
the definition satisfiable by every source. Stage 7's stricter screens are
unchanged — these classes are for stratified statistics, not for picking
radiative-transfer columns.

Where the inputs come from, per source:
  ERA5     sfc file: lsm, siconc, tcc at the analysis instants
  CARRA-2  sfc file: lsm, siconc, tcc (fields added with the 80-90N domain)
  MERRA-2  const file (FRLAND+FRLANDICE), ocn file (FRSEAICE) and rad file
           (CLDTOT) — the M2T1NX* collections are 1-h means stamped HH:30,
           so the two windows bracketing each instant are averaged, exactly
           as reanlib.fluxes does.
"""

from __future__ import annotations

import datetime as dt

import numpy as np
import xarray as xr

from .config import const_path, ocn_path, rad_path, sfc_path
from .io_era5 import open_era5

LAND_MIN_FRACTION = 0.5
SIC_OPEN_MAX = 0.05
SIC_PACK_MIN = 0.95
TCC_CLEAR_MAX = 0.05
TCC_CLOUDY_MIN = 0.95

SURFACE_LABELS = ("land", "open_water", "marginal_ice", "pack_ice")
SKY_LABELS = ("clear", "partial", "cloudy")


def surface_class(lsm: np.ndarray, siconc: np.ndarray) -> np.ndarray:
    """int8 class array from land fraction and sea-ice fraction (0..1)."""
    lsm = np.asarray(lsm, dtype=float)
    sic = np.asarray(siconc, dtype=float)
    out = np.full(np.broadcast(lsm, sic).shape, -1, dtype=np.int8)
    land = np.isfinite(lsm) & (lsm >= LAND_MIN_FRACTION)
    out[land] = 0
    sea = np.isfinite(lsm) & (lsm < LAND_MIN_FRACTION) & np.isfinite(sic)
    out[sea & (sic < SIC_OPEN_MAX)] = 1
    out[sea & (sic >= SIC_OPEN_MAX) & (sic <= SIC_PACK_MIN)] = 2
    out[sea & (sic > SIC_PACK_MIN)] = 3
    return out


def sky_class(tcc: np.ndarray) -> np.ndarray:
    """int8 class array from total cloud cover (0..1)."""
    tcc = np.asarray(tcc, dtype=float)
    out = np.full(tcc.shape, -1, dtype=np.int8)
    ok = np.isfinite(tcc)
    out[ok] = 1
    out[ok & (tcc <= TCC_CLEAR_MAX)] = 0
    out[ok & (tcc >= TCC_CLOUDY_MIN)] = 2
    return out


def _bracket_mean(da: xr.DataArray, when: np.datetime64) -> np.ndarray:
    """Mean of the HH:30-stamped windows bracketing `when` (MERRA-2 means)."""
    stamps = [when - np.timedelta64(30, "m"), when + np.timedelta64(30, "m")]
    have = [s for s in stamps if s in da["valid_time"].values]
    if not have:
        raise KeyError(f"no stamps bracketing {when} in the collection")
    return da.sel(valid_time=have).mean("valid_time").values


def load_class_inputs(cfg: dict, date: dt.date, times: np.ndarray):
    """(lsm, siconc, tcc) arrays shaped (ntime, ny, nx) at the instants.

    Raises KeyError/FileNotFoundError with a download hint when the source's
    classification fields have not been fetched yet.
    """
    src = cfg["source"]
    if src in ("era5", "carra2"):
        sfc = open_era5(sfc_path(cfg, date))
        try:
            sel = sfc.sel(valid_time=times)
            missing = [v for v in ("lsm", "siconc", "tcc") if v not in sel]
            if missing:
                dl = ("era5_download.py --datasets sfc --force" if src == "era5"
                      else "carra2_download.py --datasets sfc --force")
                raise KeyError(
                    f"{src} sfc file lacks {missing} — re-download with "
                    f"the classification variables ({dl})")
            lsm2 = sel["lsm"].values
            if lsm2.ndim == 2:                     # static field, no time dim
                lsm2 = np.broadcast_to(lsm2, sel["siconc"].shape)
            return lsm2, sel["siconc"].values, sel["tcc"].values
        finally:
            sfc.close()

    if src == "merra2":
        const = open_era5(const_path(cfg))
        land = (const["frland"] + const["frlandice"]).squeeze().values
        const.close()
        ocn = open_era5(ocn_path(cfg, date))
        rad = open_era5(rad_path(cfg, date))
        try:
            sic = np.stack([_bracket_mean(ocn["siconc"], np.datetime64(t))
                            for t in times])
            tcc = np.stack([_bracket_mean(rad["cldtot"], np.datetime64(t))
                            for t in times])
        finally:
            ocn.close()
            rad.close()
        lsm2 = np.broadcast_to(land, sic.shape)
        return lsm2, sic, tcc

    raise ValueError(f"no classification inputs defined for source {src!r}")


def class_dataset(cfg: dict, date: dt.date, times: np.ndarray,
                  dims, coords) -> xr.Dataset:
    """surface_class / sky_class as an int8 dataset on the analysis grid."""
    lsm2, sic, tcc = load_class_inputs(cfg, date, times)
    out = xr.Dataset(
        {
            "surface_class": (dims, surface_class(lsm2, sic)),
            "sky_class": (dims, sky_class(tcc)),
        },
        coords=coords,
    )
    out["surface_class"].attrs = {
        "flag_values": np.arange(len(SURFACE_LABELS), dtype=np.int8),
        "flag_meanings": " ".join(SURFACE_LABELS),
        "comment": (f"land: land fraction >= {LAND_MIN_FRACTION}; else by "
                    f"sea-ice fraction: open < {SIC_OPEN_MAX}, pack > "
                    f"{SIC_PACK_MIN}, marginal between; -1 unknown"),
    }
    out["sky_class"].attrs = {
        "flag_values": np.arange(len(SKY_LABELS), dtype=np.int8),
        "flag_meanings": " ".join(SKY_LABELS),
        "comment": (f"total cloud cover: clear <= {TCC_CLEAR_MAX}, cloudy >= "
                    f"{TCC_CLOUDY_MIN}, partial between; -1 unknown"),
    }
    return out
