"""Temperature-inversion strength metrics on reanalysis pressure-level profiles.

Three metrics are computed at every grid point and time:

1. Surface-based inversion (SBI), profile scan
   The temperature profile is scanned upward from the surface (base = 2 m
   temperature at surface pressure). The inversion top is the last running
   temperature maximum before the profile stops increasing; up to
   ``max_embedded_levels`` consecutive non-increasing levels are tolerated
   inside the inversion (the ~25 hPa level spacing makes one tolerated level
   roughly analogous to the 100 m embedded-layer rule used with radiosondes).
   Strength = T(top) - T(2 m). An inversion is only counted when the strength
   is >= ``min_strength_k`` (default 0.5 K, the criterion used with this
   algorithm at Arctic radiosonde stations, doi:10.1175/JAMC-D-21-0054.1; a
   0.3 K variant appears in the Ny-Alesund high-resolution radiosonde
   climatology, Atmos. Res. 2021).
   References: Kahl (1990), Int. J. Climatol. 10, 537-548;
   Serreze, Kahl & Schnell (1992), J. Climate 5, 615-629;
   Zhang, Seidel, Golaz, Deser & Tomas (2011), J. Climate 24, 5167-5186.

2. dt_850_2m = T(850 hPa) - T(2 m)
   Fixed-level estimate standard in Arctic climate-model studies; 850 hPa is
   chosen to sit above the boundary layer, and the 2 m temperature is used
   instead of T(1000 hPa) because winter surface pressure deviates from
   1000 hPa. References: Medeiros, Deser, Tomas & Kay (2011), J. Climate 24,
   doi:10.1175/2011JCLI3968.1; Pavelsky, Boe, Hall & Fetzer (2011),
   Clim. Dyn., doi:10.1007/s00382-010-0756-8.

3. dt_925_1000 = T(925 hPa) - T(1000 hPa)
   Simplest fixed-level metric, requiring only pressure-level data; used with
   ERA5 in e.g. arXiv:2011.11127. Related two-level stability metrics: LTS
   (Klein & Hartmann 1993, J. Climate 6, 1587-1606) and EIS (Wood &
   Bretherton 2006, J. Climate 19, 6425-6432).

Pressure levels with p > surface pressure are below ground (ERA5 extrapolates
them; MERRA-2 stores fill values, NaN after decoding) and are excluded from
the SBI scan — as are levels with non-finite temperature; the fixed-level
metrics are masked there when ``masking.mask_fixed_below_ground`` is enabled.
"""

from __future__ import annotations

import numpy as np
import xarray as xr

RD = 287.04    # dry-air gas constant [J kg-1 K-1]
G0 = 9.80665   # standard gravity [m s-2]

REFERENCES = {
    "sbi": ("Kahl (1990) Int. J. Climatol. 10, 537-548; "
            "Serreze, Kahl & Schnell (1992) J. Climate 5, 615-629; "
            "Zhang et al. (2011) J. Climate 24, 5167-5186"),
    "dt_850_2m": ("Medeiros, Deser, Tomas & Kay (2011) J. Climate 24, "
                  "doi:10.1175/2011JCLI3968.1; Pavelsky et al. (2011) "
                  "Clim. Dyn., doi:10.1007/s00382-010-0756-8"),
    "dt_925_1000": ("e.g. arXiv:2011.11127; cf. LTS: Klein & Hartmann (1993) "
                    "J. Climate 6, 1587-1606; EIS: Wood & Bretherton (2006) "
                    "J. Climate 19, 6425-6432"),
    "min_strength": ("0.5 K minimum: doi:10.1175/JAMC-D-21-0054.1 "
                     "(three Arctic stations, Kahl/Serreze algorithm); "
                     "0.3 K variant: Ny-Alesund radiosonde climatology, "
                     "Atmos. Res. 2021"),
}

SOURCE_CITATIONS = {
    "era5": "ERA5 (Hersbach et al. 2020, QJRMS 146, 1999-2049), via CDS",
    "merra2": ("MERRA-2 (Gelaro et al. 2017, J. Climate 30, 5419-5454; "
               "GMAO GEOS 5.12.4), via NASA GES DISC"),
}


def sbi_scan(T: np.ndarray, p_hpa: np.ndarray, t2m: np.ndarray, sp_hpa: np.ndarray,
             q: np.ndarray | None = None, *, top_limit_hpa: float = 500.0,
             max_embedded_levels: int = 1,
             min_strength_k: float = 0.0) -> dict[str, np.ndarray]:
    """Vectorized surface-based-inversion scan.

    Parameters
    ----------
    T : (ntime, nlev, nlat, nlon) temperature [K], levels ordered surface-first
        (p_hpa decreasing: index 0 = e.g. 1000 hPa)
    p_hpa : (nlev,) pressure levels [hPa], strictly decreasing
    t2m, sp_hpa : (ntime, nlat, nlon) 2 m temperature [K], surface pressure [hPa]
    q : optional specific humidity like T, for virtual-temperature layer
        thicknesses; omitted -> dry heights (small error, < a few % in depth_z)

    Returns dict of (ntime, nlat, nlon) arrays:
    strength [K], top_p [hPa], depth_p [hPa], depth_z [m, hypsometric,
    approximate], found [bool]. Non-found points are NaN (found False).
    """
    if not np.all(np.diff(p_hpa) < 0):
        raise ValueError("p_hpa must be strictly decreasing (surface-first)")
    T = np.asarray(T, dtype=np.float64)
    t2m = np.asarray(t2m, dtype=np.float64)
    sp_hpa = np.asarray(sp_hpa, dtype=np.float64)
    nlev = T.shape[1]
    grid_shape = t2m.shape

    runmax = t2m.copy()                                  # running profile maximum
    top_p = np.full(grid_shape, np.nan)
    top_z = np.full(grid_shape, np.nan)
    top_found = np.zeros(grid_shape, dtype=bool)
    embed = np.zeros(grid_shape, dtype=np.int16)         # consecutive dips
    terminated = np.zeros(grid_shape, dtype=bool)

    # hydrostatic state for approximate heights above the surface
    z_prev = np.zeros(grid_shape)
    p_prev = sp_hpa.copy()
    tv_prev = t2m.copy()

    for k in range(nlev):
        pk = float(p_hpa[k])
        if pk < top_limit_hpa:
            break
        tk = T[:, k]
        # above ground and finite (MERRA-2 fills below-ground levels with NaN)
        above = (pk <= sp_hpa) & np.isfinite(tk)
        tvk = tk * (1.0 + 0.61 * q[:, k]) if q is not None else tk
        dz = (RD / G0) * 0.5 * (tv_prev + tvk) * np.log(p_prev / pk)
        zk = np.where(above, z_prev + dz, 0.0)

        active = above & ~terminated
        warmer = active & (tk > runmax)
        runmax = np.where(warmer, tk, runmax)
        top_p = np.where(warmer, pk, top_p)
        top_z = np.where(warmer, zk, top_z)
        top_found |= warmer
        embed = np.where(warmer, 0, embed)
        dip = active & ~warmer
        embed = np.where(dip, embed + 1, embed).astype(np.int16)
        terminated |= embed > max_embedded_levels

        z_prev = np.where(above, zk, z_prev)
        p_prev = np.where(above, pk, p_prev)
        tv_prev = np.where(above, tvk, tv_prev)

    strength = runmax - t2m
    found = top_found & (strength >= min_strength_k)
    nan = np.full(grid_shape, np.nan)
    return {
        "strength": np.where(found, strength, nan),
        "top_p": np.where(found, top_p, nan),
        "depth_p": np.where(found, sp_hpa - top_p, nan),
        "depth_z": np.where(found, top_z, nan),
        "found": found,
    }


def column_heights(T: np.ndarray, p_hpa: np.ndarray, t2m: float, sp_hpa: float,
                   q: np.ndarray | None = None) -> np.ndarray:
    """Approximate hypsometric heights [m above surface] of one column's levels.

    T, p_hpa, q: (nlev,) surface-first (p decreasing). Below-ground levels
    (p > sp_hpa) return NaN.
    """
    z = np.full(p_hpa.shape, np.nan)
    z_prev, p_prev = 0.0, float(sp_hpa)
    tv_prev = float(t2m)
    for k in range(p_hpa.size):
        pk = float(p_hpa[k])
        if pk > sp_hpa or not np.isfinite(T[k]):
            continue
        tvk = float(T[k]) * (1.0 + 0.61 * float(q[k])) if q is not None else float(T[k])
        z[k] = z_prev + (RD / G0) * 0.5 * (tv_prev + tvk) * np.log(p_prev / pk)
        z_prev, p_prev, tv_prev = z[k], pk, tvk
    return z


def fixed_level_metrics(ds_plev: xr.Dataset, t2m: xr.DataArray,
                        sp_hpa: xr.DataArray, cfg: dict) -> xr.Dataset:
    """dt_850_2m and dt_925_1000 with optional below-ground masking."""
    t = ds_plev["t"]
    dt_850_2m = t.sel(pressure_level=850) - t2m
    dt_925_1000 = t.sel(pressure_level=925) - t.sel(pressure_level=1000)
    if cfg["masking"]["mask_fixed_below_ground"]:
        dt_850_2m = dt_850_2m.where(sp_hpa >= 850)
        dt_925_1000 = dt_925_1000.where(sp_hpa >= 1000)
    out = xr.Dataset({"dt_850_2m": dt_850_2m, "dt_925_1000": dt_925_1000})
    return out.drop_vars("pressure_level", errors="ignore")


def compute_inversion_dataset(ds_plev: xr.Dataset, ds_sfc: xr.Dataset,
                              cfg: dict) -> xr.Dataset:
    """Run all three metrics; returns one dataset on the (time, lat, lon) grid."""
    for coord in ("latitude", "longitude"):
        if not np.array_equal(ds_plev[coord].values, ds_sfc[coord].values):
            raise ValueError(f"{coord} grids of plev and sfc files differ")
    common = np.intersect1d(ds_plev["valid_time"].values, ds_sfc["valid_time"].values)
    if common.size == 0:
        raise ValueError("plev and sfc files share no valid_time steps")
    ds_plev = ds_plev.sel(valid_time=common)
    ds_sfc = ds_sfc.sel(valid_time=common)

    dims = ("valid_time", "pressure_level", "latitude", "longitude")
    ds_plev = ds_plev.transpose(*dims)
    t2m, skt = ds_sfc["t2m"], ds_sfc["skt"]
    sp_hpa = ds_sfc["sp"] / 100.0

    sbi_cfg = cfg["sbi"]
    sbi = sbi_scan(
        ds_plev["t"].values,
        ds_plev["pressure_level"].values.astype(float),
        t2m.values, sp_hpa.values,
        q=ds_plev["q"].values if "q" in ds_plev else None,
        top_limit_hpa=float(sbi_cfg["top_limit_hpa"]),
        max_embedded_levels=int(sbi_cfg["max_embedded_levels"]),
        min_strength_k=float(sbi_cfg["min_strength_k"]),
    )

    hdims = ("valid_time", "latitude", "longitude")
    coords = {d: ds_sfc[d] for d in hdims}
    out = xr.Dataset(
        {
            "sbi_strength": (hdims, sbi["strength"].astype(np.float32)),
            "sbi_top_p": (hdims, sbi["top_p"].astype(np.float32)),
            "sbi_depth_p": (hdims, sbi["depth_p"].astype(np.float32)),
            "sbi_depth_z": (hdims, sbi["depth_z"].astype(np.float32)),
            "sbi_found": (hdims, sbi["found"].astype(np.int8)),
        },
        coords=coords,
    )
    out = out.merge(fixed_level_metrics(ds_plev, t2m, sp_hpa, cfg).astype(np.float32))
    out["t2m"] = t2m.astype(np.float32)
    out["skt"] = skt.astype(np.float32)
    out["sp"] = ds_sfc["sp"].astype(np.float32)

    meta = {
        "sbi_strength": ("K", "SBI strength: T(inversion top) - T(2 m)", "sbi"),
        "sbi_top_p": ("hPa", "pressure of the SBI top", "sbi"),
        "sbi_depth_p": ("hPa", "SBI depth: surface pressure - top pressure", "sbi"),
        "sbi_depth_z": ("m", "SBI depth above surface (hypsometric, approximate)", "sbi"),
        "sbi_found": ("1", "1 where a surface-based inversion was detected", "sbi"),
        "dt_850_2m": ("K", "T(850 hPa) - T(2 m)", "dt_850_2m"),
        "dt_925_1000": ("K", "T(925 hPa) - T(1000 hPa)", "dt_925_1000"),
        "t2m": ("K", "2 m temperature (from the single-level file)", ""),
        "skt": ("K", "skin temperature (from the single-level file)", ""),
        "sp": ("Pa", "surface pressure (from the single-level file)", ""),
    }
    for name, (units, desc, refkey) in meta.items():
        out[name].attrs["units"] = units
        out[name].attrs["long_name"] = desc
        if refkey:
            out[name].attrs["references"] = REFERENCES[refkey]

    src = cfg.get("source", "era5")
    out.attrs = {
        "title": f"{src.upper()} Arctic temperature-inversion strength metrics",
        "sbi_parameters": (f"top_limit_hpa={sbi_cfg['top_limit_hpa']}, "
                           f"max_embedded_levels={sbi_cfg['max_embedded_levels']}, "
                           f"min_strength_k={sbi_cfg['min_strength_k']}"),
        "sbi_min_strength_reference": REFERENCES["min_strength"],
        "expver_plev": ds_plev.attrs.get("expver_values", "unknown"),
        "expver_sfc": ds_sfc.attrs.get("expver_values", "unknown"),
        "source": SOURCE_CITATIONS.get(src, src),
    }
    return out
