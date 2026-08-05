"""Configuration loading and canonical file paths for the reanalysis pipeline.

Every path helper takes the loaded ``cfg`` dict, whose ``source`` key
("era5" or "merra2") selects the per-source subtree and filename stem:

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

SOURCES = ("era5", "merra2")
SOURCE_LABELS = {"era5": "ERA5", "merra2": "MERRA-2"}


def source_label(cfg: dict) -> str:
    """Display name of the configured data source (e.g. for figure titles)."""
    return SOURCE_LABELS[cfg["source"]]

DEFAULTS: dict = {
    "source": "era5",
    "area": [90, -180, 80, 180],  # N, W, S, E
    "paths": {"data": "data", "derived": "derived", "figures": "figures"},
    "download": {
        "default_hours": [0, 6, 12, 18],
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
        "collections": {"plev": "M2I3NPASM", "sfc": "M2I1NXASM"},
        "version": "5.12.4",
        "plev_variables": ["T", "QV", "O3", "QL", "QI"],
        "sfc_variables": ["T2M", "TS", "PS"],
        "default_hours": [0, 6, 12, 18],
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
