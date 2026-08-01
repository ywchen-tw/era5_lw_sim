"""Configuration loading and canonical file paths for the ERA5 pipeline."""

from __future__ import annotations

import copy
import datetime as dt
from pathlib import Path

import yaml

# this file lives at <repo>/src/era5lib/config.py
REPO_ROOT = Path(__file__).resolve().parent.parent.parent

DEFAULTS: dict = {
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


def load_config(path: str | Path | None = None) -> dict:
    """Built-in defaults overlaid with config.yaml (repo root, or explicit path)."""
    if path is None:
        path = REPO_ROOT / "config.yaml"
    path = Path(path)
    cfg = copy.deepcopy(DEFAULTS)
    if path.exists():
        with open(path) as f:
            overlay = yaml.safe_load(f) or {}
        cfg = _merge(cfg, overlay)
    return cfg


def _day_dir(cfg: dict, kind: str, date: dt.date) -> Path:
    root = Path(cfg["paths"][kind])
    if not root.is_absolute():
        root = REPO_ROOT / root
    return root / f"{date.year:04d}" / f"{date.month:02d}" / f"{date.day:02d}"


def plev_path(cfg: dict, date: dt.date) -> Path:
    return _day_dir(cfg, "data", date) / f"era5_plev_{date:%Y%m%d}.nc"


def sfc_path(cfg: dict, date: dt.date) -> Path:
    return _day_dir(cfg, "data", date) / f"era5_sfc_{date:%Y%m%d}.nc"


def inversion_path(cfg: dict, date: dt.date) -> Path:
    return _day_dir(cfg, "derived", date) / f"era5_inversion_{date:%Y%m%d}.nc"


def figures_dir(cfg: dict) -> Path:
    root = Path(cfg["paths"]["figures"])
    if not root.is_absolute():
        root = REPO_ROOT / root
    return root
