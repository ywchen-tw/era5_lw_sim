"""Configuration loading and canonical file paths for the reanalysis pipeline.

Every path helper takes the loaded ``cfg`` dict, whose ``source`` key
("era5", "merra2" or "carra2") selects the per-source subtree and filename stem:

    data/<source>/YYYY/MM/DD/<source>_{plev,sfc}_YYYYMMDD.nc
    derived/<source>/YYYY/MM/DD/<source>_inversion_YYYYMMDD.nc
    derived/<source>/YYYY/MM/<source>_inversion_monthly_YYYYMM.nc  (etc.)
    figures/<source>/
"""

from __future__ import annotations

import copy
import datetime as dt
from pathlib import Path

import yaml

# this file lives at <repo>/src/reanlib/config.py
REPO_ROOT = Path(__file__).resolve().parent.parent.parent

SOURCES = ("era5", "merra2", "carra2")
SOURCE_LABELS = {"era5": "ERA5", "merra2": "MERRA-2", "carra2": "CARRA-2"}


def source_label(cfg: dict) -> str:
    """Display name of the configured data source (e.g. for figure titles)."""
    return SOURCE_LABELS[cfg["source"]]

DEFAULTS: dict = {
    "source": "era5",
    "area": [90, -180, 80, 180],  # N, W, S, E
    "paths": {"data": "data", "derived": "derived", "figures": "figures"},
    "download": {
        "default_hours": [0, 6, 12, 18],
        "state_cadence_h": 6,
        "plev_variables": [
            "fraction_of_cloud_cover",
            "ozone_mass_mixing_ratio",
            "relative_humidity",
            "specific_cloud_ice_water_content",
            "specific_cloud_liquid_water_content",
            "specific_humidity",
            "temperature",
        ],
        "sfc_variables": ["2m_temperature", "skin_temperature", "surface_pressure",
                          "surface_sensible_heat_flux", "surface_latent_heat_flux",
                          "surface_net_solar_radiation", "surface_net_thermal_radiation",
                          "surface_solar_radiation_downwards",
                          "surface_thermal_radiation_downwards"],
        "pressure_levels": [1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175,
                            200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650,
                            700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000],
    },
    "merra2": {
        "collections": {"plev": "M2I3NPASM", "sfc": "M2I1NXASM",
                        "rad": "M2T1NXRAD"},
        "version": "5.12.4",
        "plev_variables": ["T", "QV", "O3", "QL", "QI"],
        "sfc_variables": ["T2M", "TS", "PS"],
        "rad_variables": ["LWGAB", "LWGEM", "LWGABCLR", "EMIS", "CLDTOT"],
        "default_hours": [0, 6, 12, 18],
        "state_cadence_h": 3,
    },
    "carra2": {
        # the sub-daily pan-Arctic entry; reanalysis-pan-carra-means holds only
        # daily/monthly aggregates and cannot feed the inversion metrics
        "dataset": "reanalysis-pan-carra",
        # CARRA-2 keeps its own, smaller domain: cells scale as
        # tan^2((90-lat)/2), so 85-90N is a quarter of the 80-90N cap
        # (199 k vs 796 k cells per field, ~318 MB vs 1.27 GB of profiles per
        # day). The global `area` stays 80-90N for ERA5 and MERRA-2.
        "area": [90, -180, 85, 180],
        # CARRA-2 carries no specific humidity on pressure levels, so q is
        # derived from relative_humidity (see reanlib/io_carra2.py)
        "plev_variables": ["temperature", "relative_humidity",
                           "specific_cloud_liquid_water_content",
                           "specific_cloud_ice_water_content", "cloud_cover"],
        "sfc_variables": ["2m_temperature", "skin_temperature", "surface_pressure",
                          "land_sea_mask"],
        # the 20 CARRA-2 pressure levels; the top is 50 hPa, not ERA5's 1 hPa
        "pressure_levels": [50, 70, 100, 150, 200, 250, 300, 400, 500, 600, 700,
                            750, 800, 825, 850, 875, 900, 925, 950, 1000],
        "default_hours": [0, 6, 12, 18],
        "state_cadence_h": 3,
        # saturation reference for the RH -> q conversion. CARRA documents its
        # relative humidity against saturation over water; "ice" and "mixed"
        # (water above 0 C, ice below -23 C, blended between) are available for
        # sensitivity tests.
        "rh_over": "water",
    },
    "sbi": {"top_limit_hpa": 500, "max_embedded_levels": 1, "min_strength_k": 0.5},
    "masking": {"mask_fixed_below_ground": True},
}


def _merge(base: dict, overlay: dict) -> dict:
    out = copy.deepcopy(base)
    for key, val in overlay.items():
        if isinstance(val, dict) and isinstance(out.get(key), dict):
            out[key] = _merge(out[key], val)
        else:
            out[key] = copy.deepcopy(val)
    return out


def load_config(path: str | Path | None = None, source: str | None = None) -> dict:
    """Built-in defaults overlaid with config.yaml (repo root, or explicit path).

    ``source`` (e.g. a --source CLI argument) overrides the YAML ``source:`` key.
    """
    if path is None:
        path = REPO_ROOT / "config.yaml"
    path = Path(path)
    cfg = copy.deepcopy(DEFAULTS)
    if path.exists():
        with open(path) as f:
            overlay = yaml.safe_load(f) or {}
        cfg = _merge(cfg, overlay)
    if source is not None:
        cfg["source"] = source
    if cfg["source"] not in SOURCES:
        raise ValueError(f"unknown source {cfg['source']!r}; expected one of {SOURCES}")
    return cfg


def source_block(cfg: dict) -> dict:
    """The config block holding the current source's download settings.

    ERA5's lives under ``download:`` for historical reasons; the other sources
    use a block named after themselves.
    """
    source = cfg["source"]
    return cfg["download"] if source == "era5" else cfg[source]


def source_area(cfg: dict) -> list[float]:
    """Download domain [N, W, S, E] for the current source.

    Falls back to the global ``area:`` unless the source's own block overrides
    it. CARRA-2 does: at 2.5 km the 80-90N cap is ~796 k cells per field, so
    its domain is set separately from the global one the coarser sources use.
    """
    return list(source_block(cfg).get("area", cfg["area"]))


def state_cadence_h(cfg: dict) -> int:
    """Hour spacing of the analysis states used for satellite collocation.

    ERA5 is downloaded at synoptic hours (6-hourly by default); MERRA-2's plev
    collection (M2I3NPASM) and CARRA-2's analyses are natively 3-hourly,
    halving the worst-case state-time offset. Override per source with
    ``state_cadence_h`` in the ``download:`` / ``merra2:`` / ``carra2:``
    config blocks.
    """
    return int(source_block(cfg).get("state_cadence_h",
                                     6 if cfg["source"] == "era5" else 3))


def _root(cfg: dict, kind: str) -> Path:
    root = Path(cfg["paths"][kind])
    if not root.is_absolute():
        root = REPO_ROOT / root
    return root / cfg["source"]


def _day_dir(cfg: dict, kind: str, date: dt.date) -> Path:
    return _root(cfg, kind) / f"{date.year:04d}" / f"{date.month:02d}" / f"{date.day:02d}"


def _month_dir(cfg: dict, kind: str, year: int, month: int) -> Path:
    return _root(cfg, kind) / f"{year:04d}" / f"{month:02d}"


def plev_path(cfg: dict, date: dt.date) -> Path:
    return _day_dir(cfg, "data", date) / f"{cfg['source']}_plev_{date:%Y%m%d}.nc"


def sfc_path(cfg: dict, date: dt.date) -> Path:
    return _day_dir(cfg, "data", date) / f"{cfg['source']}_sfc_{date:%Y%m%d}.nc"


def rad_path(cfg: dict, date: dt.date) -> Path:
    """Surface-radiation daily file (MERRA-2 only; ERA5 radiation lives in
    the sfc file as strd/str accumulations)."""
    if cfg["source"] != "merra2":
        raise ValueError("rad_path is MERRA-2-only; ERA5 radiation is in sfc_path")
    return _day_dir(cfg, "data", date) / f"{cfg['source']}_rad_{date:%Y%m%d}.nc"


def inversion_path(cfg: dict, date: dt.date) -> Path:
    return _day_dir(cfg, "derived", date) / f"{cfg['source']}_inversion_{date:%Y%m%d}.nc"


def monthly_path(cfg: dict, year: int, month: int) -> Path:
    return (_month_dir(cfg, "derived", year, month)
            / f"{cfg['source']}_inversion_monthly_{year:04d}{month:02d}.nc")


def analysis_path(cfg: dict, year: int, month: int) -> Path:
    return (_month_dir(cfg, "derived", year, month)
            / f"{cfg['source']}_profile_analysis_{year:04d}{month:02d}.nc")


def pairs_path(cfg: dict, year: int, month: int) -> Path:
    return (_month_dir(cfg, "derived", year, month)
            / f"{cfg['source']}_mosaic_pairs_{year:04d}{month:02d}.nc")


def figures_dir(cfg: dict) -> Path:
    root = Path(cfg["paths"]["figures"])
    if not root.is_absolute():
        root = REPO_ROOT / root
    return root / cfg["source"]
