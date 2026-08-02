#!/usr/bin/env python
"""PREFIRE-TIRS brightness-temperature simulation from ERA5 columns (stage 8).

Collocates PREFIRE 1B-RAD footprints (data/prefire/, see
era5_prefire_download.py) with ERA5 analysis columns, simulates the TIRS
channel spectrum with libRadtran (thermal-source nadir-ish radiance,
``reptran fine``), converts to channel brightness temperature with the
mission SRF files, and cross-checks at band level with RRTMG-LW
(``return_spectral_olr``). Finite-difference Jacobians (skt, T/q per level,
cloud water path / r_eff / cloud-top, emissivity) are the input for a later
cloud-property retrieval (validation against collocated EarthCARE products
is a separate future stage).

Channel BT convention: TIRS channel values are SRF-weighted mean spectral
radiances (W m-2 sr-1 um-1); BT inverts the SRF file's own blackbody
radiance lookup ``rad(T_grid, channel, scene)``, so simulated and observed
BT are on the identical channel radiometric scale. Each of the 8 cross-track
scenes has its own wavelength registration — channel BT is always computed
per scene with that scene's SRF.

Subcommands / envs:
  collocate (era5 env)     footprints -> (ERA5 cell, hour) columns, pick test set
  prep      (era5 env)     per-column uvspec profile/cloud files + manifest
  run       (er3t_env)     baseline spectral radiance -> channel BT
  rrtmg     (era5 env)     RRTMG-LW 16-band flux-equivalent BT cross-check
  jacobian  (er3t_env lrt / era5 rrtmg)  perturbation runs -> K matrices
  figure    (era5 env)     BT spectra vs obs + Jacobian heatmaps

Examples:
    conda activate era5     && python src/era5_prefire_bt.py collocate --year 2025 --month 1 --sat 1
    conda activate era5     && python src/era5_prefire_bt.py prep --year 2025 --month 1 --sat 1
    conda activate er3t_env && python src/era5_prefire_bt.py run --year 2025 --month 1 --sat 1
    conda activate era5     && python src/era5_prefire_bt.py rrtmg --year 2025 --month 1 --sat 1
    conda activate er3t_env && python src/era5_prefire_bt.py jacobian --year 2025 --month 1 --sat 1
    conda activate era5     && python src/era5_prefire_bt.py figure --year 2025 --month 1 --sat 1
"""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import glob
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from era5lib.config import REPO_ROOT, figures_dir, load_config, plev_path, sfc_path
from era5lib.io_era5 import open_era5
from era5_lrt_sim import (EMISSIVITY, REFF_ICE_UM, REFF_LIQ_UM,
                          build_profile_files, load_cams_profiles,
                          write_cloud_file)
from era5_prefire_download import srf_path

G0 = 9.80665
TRAPZ = getattr(np, "trapezoid", getattr(np, "trapz", None))
N_CHANNEL, N_SCENE = 63, 8

# 2025 well-mixed GHG levels (CAMS EGG4 ends in 2020; flat profiles keep the
# stage-7 gas plumbing unchanged — a ~10 ppm CO2 offset moves 15-um BT < 0.2 K)
CO2_2025_PPM = 424.0
CH4_2025_PPM = 1.93

# perturbation sizes for the finite-difference Jacobian
DT_K = 1.0            # skt and per-level air temperature
DQ_FRAC = 0.05        # per-level fractional water-vapour change
DWP_FRAC = 0.10       # fractional liquid/ice water path change
DREFF_UM = 1.0        # effective-radius change
DEMIS = -0.01         # emissivity change
DTAU_FRAC = 0.05      # fractional optical-thickness change (cotscan)


# ---------------------------------------------------------------- paths

def pbt_dir(cfg: dict, year: int, month: int) -> Path:
    root = Path(cfg["paths"]["derived"])
    if not root.is_absolute():
        root = REPO_ROOT / root
    return root / f"{year:04d}" / f"{month:02d}" / "prefire_bt"


def collocation_path(cfg, year, month, sat) -> Path:
    return pbt_dir(cfg, year, month) / f"collocation_{year:04d}{month:02d}_sat{sat}.json"


def manifest_path(cfg, year, month, sat) -> Path:
    return pbt_dir(cfg, year, month) / f"manifest_{year:04d}{month:02d}_sat{sat}.json"


def results_path(cfg, year, month, sat) -> Path:
    return pbt_dir(cfg, year, month) / f"results_{year:04d}{month:02d}_sat{sat}.json"


def jacobian_path(cfg, year, month, sat, label, simulator) -> Path:
    return (pbt_dir(cfg, year, month) /
            f"jacobian_{simulator}_{label}_sat{sat}.nc")


# ---------------------------------------------------------------- SRF helpers

def load_srf(path: str) -> dict:
    """TIRS SRF file -> numpy dict (wavelengths um, radiances W m-2 sr-1 um-1)."""
    import netCDF4 as nc

    ds = nc.Dataset(path)
    out = {
        "wavelen": ds["wavelen"][:].filled(np.nan),
        "srf_normed": ds["srf_normed"][:].filled(0.0),        # (spec, ch, scene)
        "wl1": ds["channel_wavelen1"][:].filled(np.nan),      # (ch, scene)
        "wl2": ds["channel_wavelen2"][:].filled(np.nan),
        "mean_wl": ds["channel_mean_wavelen"][:].filled(np.nan),
        "bitflags": ds["detector_bitflags"][:].filled(1),
        "T_grid": ds["T_grid"][:].filled(np.nan),
        "rad_lut": ds["rad"][:].filled(np.nan),               # (T, ch, scene)
        "nedr": ds["NEDR"][:].filled(np.nan),                 # (ch, scene)
    }
    ds.close()
    return out


def good_channels(srf: dict, scene: int) -> np.ndarray:
    return np.where(srf["bitflags"][:, scene] == 0)[0]


def sim_wvl_range_nm(srf: dict, margin_um: float = 0.5) -> "tuple[float, float]":
    """Wavelength span (nm) covering every usable channel of every scene."""
    good = srf["bitflags"] == 0
    lo = float(np.nanmin(np.where(good, srf["wl1"], np.nan))) - margin_um
    hi = float(np.nanmax(np.where(good, srf["wl2"], np.nan))) + margin_um
    return max(lo, 2.6) * 1000.0, hi * 1000.0   # reptran thermal starts 2.5 um


def channel_bt(srf: dict, scene: int, wl_um: np.ndarray,
               rad_per_um: np.ndarray) -> "tuple[np.ndarray, np.ndarray]":
    """SRF-convolve a spectrum -> (channel radiance, channel BT), NaN if unusable.

    rad_per_um: spectral radiance W m-2 sr-1 um-1 on wl_um (ascending um).
    Channel radiance = sum(srf_normed * L) (srf_normed columns sum to 1);
    BT inverts the SRF file's blackbody channel-radiance lookup.
    """
    rad_c = np.full(N_CHANNEL, np.nan)
    bt_c = np.full(N_CHANNEL, np.nan)
    L = np.interp(srf["wavelen"], wl_um, rad_per_um,
                  left=np.nan, right=np.nan)
    cov = np.isfinite(L)
    for c in good_channels(srf, scene):
        w = srf["srf_normed"][:, c, scene]
        wsum = float(np.sum(w[cov]))
        # SRFs carry per-mil-level far tails across the whole spectrograph
        # range; renormalize over the simulated span but skip a channel whose
        # response is genuinely not covered
        if wsum < 0.995:
            continue
        rad_c[c] = float(np.sum(w[cov] * L[cov]) / wsum)
        lut = srf["rad_lut"][:, c, scene]
        bt_c[c] = float(np.interp(rad_c[c], lut, srf["T_grid"]))
    return rad_c, bt_c


def obs_channel_bt_to_rad(srf: dict, scene: int, bt: np.ndarray) -> np.ndarray:
    """Observed channel BT -> channel radiance via the same blackbody lookup."""
    rad = np.full(N_CHANNEL, np.nan)
    ok = np.isfinite(bt)
    for c in np.where(ok)[0]:
        rad[c] = float(np.interp(bt[c], srf["T_grid"], srf["rad_lut"][:, c, scene]))
    return rad


# ---------------------------------------------------------------- collocate

def granule_files(cfg, year, month, sat) -> "list[str]":
    pat = str(Path(cfg["paths"]["data"]) / "prefire" / f"{year:04d}" /
              f"{month:02d}" / f"PREFIRE_SAT{sat}_1B-RAD_*.nc")
    return sorted(glob.glob(pat))


def cmd_collocate(args) -> int:
    import netCDF4 as nc

    cfg = load_config(args.config)
    files = granule_files(cfg, args.year, args.month, args.sat)
    if not files:
        sys.exit("no granules — run era5_prefire_download.py first")
    south = float(cfg["area"][2])
    out_dir = pbt_dir(cfg, args.year, args.month)
    out_dir.mkdir(parents=True, exist_ok=True)

    era5_cache = {}

    def era5_day(date):
        if date not in era5_cache:
            p = plev_path(cfg, date)
            era5_cache[date] = open_era5(p) if p.exists() else None
        return era5_cache[date]

    groups = {}          # (etime iso, iy, ix) -> footprint list
    n_seen = n_dom = n_noera5 = 0
    for fn in files:
        ds = nc.Dataset(fn)
        g, b = ds["Geometry"], ds["BT"]
        lat = g["latitude"][:].filled(np.nan)
        lon = g["longitude"][:].filled(np.nan)
        vza = g["viewing_zenith_angle"][:].filled(np.nan)
        tv = g["time_UTC_values"][:].filled(0)
        bt = b["spectral_BT"][:].filled(np.nan)
        btq = b["BT_quality_flag"][:].filled(9)
        bt[btq != 0] = np.nan
        n_seen += lat.size

        rows = np.where(np.nanmax(lat, axis=1) >= south)[0]
        for ia in rows:
            t = dt.datetime(*[int(v) for v in tv[ia, :5]])
            hour6 = int(round((t.hour + t.minute / 60.0) / 6.0)) * 6
            et = (dt.datetime(t.year, t.month, t.day)
                  + dt.timedelta(hours=hour6))
            plev = era5_day(et.date())
            if plev is None:
                n_noera5 += 1
                continue
            lats = plev["latitude"].values
            lons = plev["longitude"].values
            for sc in range(N_SCENE):
                if not np.isfinite(lat[ia, sc]) or lat[ia, sc] < south:
                    continue
                if not np.any(np.isfinite(bt[ia, sc])):
                    continue
                n_dom += 1
                iy = int(np.abs(lats - lat[ia, sc]).argmin())
                dlon = (lons - lon[ia, sc] + 180.0) % 360.0 - 180.0
                ix = int(np.abs(dlon).argmin())
                key = (et.isoformat(), iy, ix)
                groups.setdefault(key, []).append({
                    "granule": os.path.basename(fn), "atrack": int(ia),
                    "scene": int(sc), "time": t.isoformat(),
                    "lat": float(lat[ia, sc]), "lon": float(lon[ia, sc]),
                    "vza": float(vza[ia, sc]),
                    "bt": [None if not np.isfinite(v) else round(float(v), 3)
                           for v in bt[ia, sc]],
                })
        ds.close()
        print(f"  {os.path.basename(fn)}: cumulative {len(groups)} columns")

    print(f"{n_seen} footprints scanned, {n_dom} in-domain with good BT, "
          f"{n_noera5} skipped (no ERA5 day file), {len(groups)} unique "
          f"(cell, hour) columns")

    # classify sky per column from the ERA5 column, then pick the test set
    columns = []
    for (et_iso, iy, ix), fps in sorted(groups.items()):
        et = dt.datetime.fromisoformat(et_iso)
        plev = era5_day(et.date())
        p5 = plev.sel(valid_time=et).transpose("pressure_level", "latitude",
                                               "longitude")
        p = p5["pressure_level"].values.astype(float)
        cc = p5["cc"].values[:, iy, ix].astype(float)
        clwc = p5["clwc"].values[:, iy, ix].astype(float)
        ciwc = p5["ciwc"].values[:, iy, ix].astype(float)
        p_pa_asc = p[::-1] * 100.0
        lwp = float(TRAPZ(clwc[::-1], x=p_pa_asc) / G0 * 1e3)
        iwp = float(TRAPZ(ciwc[::-1], x=p_pa_asc) / G0 * 1e3)
        cc_max = float(np.nanmax(cc))
        twp = lwp + iwp
        if cc_max <= 0.01 and twp <= 1.0:
            sky = "clear"
        elif cc_max >= 0.99 and twp > 1.0:
            sky = "overcast"
        else:
            sky = "partial"
        columns.append({
            "label": f"{et:%Y%m%dT%H}Z_{iy:02d}_{ix:03d}",
            "etime": et_iso, "iy": iy, "ix": ix,
            "lat": float(p5["latitude"].values[iy]),
            "lon": float(p5["longitude"].values[ix]),
            "sky": sky, "cc_max": cc_max, "lwp_g": lwp, "iwp_g": iwp,
            "n_footprints": len(fps),
            "mean_dt_h": float(np.mean(
                [abs((dt.datetime.fromisoformat(f["time"]) - et)
                     .total_seconds()) / 3600.0 for f in fps])),
            "footprints": fps,
        })

    n_sky = {s: sum(1 for c in columns if c["sky"] == s)
             for s in ("clear", "partial", "overcast")}
    print(f"sky classes: {n_sky}")

    # test set: clear columns with the most footprints (best obs averaging);
    # overcast ranked by water path, alternating ice- and liquid-dominated so
    # both phases appear if available. Partial columns are excluded — BT
    # blending across a broken scene is nonlinear and not a clean test.
    clear = sorted([c for c in columns if c["sky"] == "clear"],
                   key=lambda c: -c["n_footprints"])[:args.n_clear]
    over = sorted([c for c in columns if c["sky"] == "overcast"],
                  key=lambda c: -(c["lwp_g"] + c["iwp_g"]))
    ice = [c for c in over if c["iwp_g"] >= c["lwp_g"]]
    liq = [c for c in over if c["lwp_g"] > c["iwp_g"]]
    cloudy = []
    while len(cloudy) < args.n_cloudy and (ice or liq):
        for pool in (ice, liq):
            if pool and len(cloudy) < args.n_cloudy:
                cloudy.append(pool.pop(0))
    selected = {c["label"] for c in clear + cloudy}
    for c in columns:
        c["selected"] = c["label"] in selected
        if not c["selected"]:
            c["footprints"] = []      # keep the file small
    print("selected test columns:")
    for c in clear + cloudy:
        print(f"  {c['label']} {c['sky']:9s} lat {c['lat']:6.2f} "
              f"lwp {c['lwp_g']:6.1f} iwp {c['iwp_g']:6.1f} g/m2 "
              f"footprints {c['n_footprints']}")

    path = collocation_path(cfg, args.year, args.month, args.sat)
    path.write_text(json.dumps({
        "year": args.year, "month": args.month, "sat": args.sat,
        "n_columns": len(columns), "sky_counts": n_sky,
        "columns": columns}, indent=1))
    print(f"wrote {path}")
    return 0


# ---------------------------------------------------------------- prep

def cmd_prep(args) -> int:
    cfg = load_config(args.config)
    cpath = collocation_path(cfg, args.year, args.month, args.sat)
    if not cpath.exists():
        sys.exit(f"missing {cpath} — run collocate first")
    coll = json.loads(cpath.read_text())
    out_dir = pbt_dir(cfg, args.year, args.month)

    cams = load_cams_profiles(args.year, args.month)
    if cams is None:
        cams = {"p_hpa": np.array([1.0, 1100.0]),
                "x_co2": np.full(2, CO2_2025_PPM * 1e-6),
                "x_ch4": np.full(2, CH4_2025_PPM * 1e-6),
                "source": (f"constants CO2 {CO2_2025_PPM} ppm / CH4 "
                           f"{CH4_2025_PPM} ppm (CAMS EGG4 ends 2020)")}

    spath = srf_path(cfg, args.sat)
    if not spath.exists():
        sys.exit(f"missing {spath} — run era5_prefire_download.py --srf-only")
    srf = load_srf(str(spath))
    wvl_lo_nm, wvl_hi_nm = sim_wvl_range_nm(srf)
    print(f"simulation range {wvl_lo_nm / 1000.0:.2f}-{wvl_hi_nm / 1000.0:.2f} um "
          f"(usable TIRS{args.sat} channels + margin)")

    era5_cache = {}
    columns = []
    for col in [c for c in coll["columns"] if c["selected"]]:
        et = dt.datetime.fromisoformat(col["etime"])
        if et.date() not in era5_cache:
            era5_cache[et.date()] = (open_era5(plev_path(cfg, et.date())),
                                     open_era5(sfc_path(cfg, et.date())))
        plev, sfc = era5_cache[et.date()]
        iy, ix = col["iy"], col["ix"]
        p5 = plev.sel(valid_time=et).transpose("pressure_level", "latitude",
                                               "longitude")
        s5 = sfc.sel(valid_time=et)
        p = p5["pressure_level"].values.astype(float)
        sp_hpa = float(s5["sp"].values[iy, ix]) / 100.0
        t2m = float(s5["t2m"].values[iy, ix])
        skt = float(s5["skt"].values[iy, ix])
        above = p <= sp_hpa
        cols = {v: p5[v].values[above, iy, ix].astype(float)
                for v in ("t", "q", "o3", "cc", "clwc", "ciwc")}

        label = col["label"]
        atm_file, ch4_file, keep, z_keep = build_profile_files(
            out_dir, label, p[above], cols["t"], cols["q"], cols["o3"],
            t2m, sp_hpa, cams)
        wc_file = write_cloud_file(out_dir / f"wc_{label}.dat", z_keep,
                                   cols["clwc"][keep], p[above][keep],
                                   cols["t"][keep], REFF_LIQ_UM)
        ic_file = write_cloud_file(out_dir / f"ic_{label}.dat", z_keep,
                                   cols["ciwc"][keep], p[above][keep],
                                   cols["t"][keep], REFF_ICE_UM)

        # per-scene mean observed BT (scenes have distinct wavelength grids)
        scenes = {}
        for f in col["footprints"]:
            scenes.setdefault(f["scene"], []).append(
                [np.nan if v is None else v for v in f["bt"]])
        obs_bt = {str(sc): np.nanmean(np.array(v, dtype=float), axis=0)
                  .round(3).tolist() for sc, v in sorted(scenes.items())}
        vza = float(np.mean([f["vza"] for f in col["footprints"]]))

        columns.append({
            **{k: col[k] for k in ("label", "etime", "iy", "ix", "lat", "lon",
                                   "sky", "cc_max", "lwp_g", "iwp_g",
                                   "n_footprints", "mean_dt_h")},
            "skt": skt, "t2m": t2m, "sp_hpa": sp_hpa, "vza": vza,
            "obs_bt": obs_bt,
            "atm_file": atm_file, "ch4_file": ch4_file,
            "wc_file": wc_file, "ic_file": ic_file,
        })
        print(f"  {label}: {col['sky']}, skt {skt:.1f} K, vza {vza:.1f} deg, "
              f"scenes {sorted(scenes)}")

    manifest = {
        "year": args.year, "month": args.month, "sat": args.sat,
        "srf_file": str(spath), "emissivity": EMISSIVITY,
        "wavelength_range_nm": [round(wvl_lo_nm), round(wvl_hi_nm)],
        "reff_um": {"liquid": REFF_LIQ_UM, "ice": REFF_ICE_UM},
        "gas_source": cams["source"],
        "columns": columns,
    }
    mpath = manifest_path(cfg, args.year, args.month, args.sat)
    mpath.write_text(json.dumps(manifest, indent=1))
    print(f"wrote {mpath} ({len(columns)} columns)")
    return 0


# ---------------------------------------------------------------- run (er3t_env)

def resolve_mol_abs_param(args) -> str:
    """coarse on the local Mac, fine on Linux/CURC.

    reptran fine loads a large thermal table per uvspec worker; running it
    locally with parallel workers has exhausted memory and shut the machine
    down. Refuse it on Darwin unless explicitly forced.
    """
    import platform as _platform

    local = _platform.system() == "Darwin"
    choice = args.mol_abs_param or ("coarse" if local else "fine")
    if local and choice != "coarse" and not args.force_local:
        sys.exit(f"refusing 'reptran {choice}' on this machine — it has "
                 "caused an out-of-memory shutdown. Run on CURC "
                 "(slurm/curc_prefire_bt.sh) or pass --force-local to "
                 "override.")
    return choice


def radiance_init(tmp: Path, tag: str, manifest: dict, col: dict,
                  atm_file: str, wc_file, ic_file, skt: float,
                  emissivity: float, lrt_cfg_base, er3t, extra_opts=None):
    """One thermal-source spectral radiance uvspec job (TOA, sensor vza)."""
    lrt_cfg = copy.deepcopy(lrt_cfg_base)
    lrt_cfg["atmosphere_file"] = atm_file
    extra = {
        "source": "thermal",
        "sur_temperature": f"{skt:.2f}",
        "wavelength_add": "{:.0f} {:.0f}".format(
            *manifest["wavelength_range_nm"]),
        "mol_file": "CH4 " + col["ch4_file"],
    }
    if wc_file:
        extra["wc_file 1D"] = wc_file
        extra["wc_properties"] = "mie interpolate"
    if ic_file:
        extra["ic_file 1D"] = ic_file
        extra["ic_properties"] = "yang2013 interpolate"
    if extra_opts:
        extra.update(extra_opts)
    init = er3t.rtm.lrt.lrt_init_mono_rad(
        input_file=str(tmp / f"input_{tag}.txt"),
        output_file=str(tmp / f"output_{tag}.txt"),
        date=dt.datetime.fromisoformat(col["etime"]),
        surface_albedo=1.0 - emissivity,
        solar_zenith_angle=80.0,
        sensor_zenith_angle=float(col["vza"]),
        output_altitude="toa",
        input_dict_extra=extra,
        mute_list=["wavelength", "spline", "source solar",
                   "slit_function_file"],
        lrt_cfg=lrt_cfg, cld_cfg=None, aer_cfg=None)
    return init


def read_spectrum(output_file: str) -> "tuple[np.ndarray, np.ndarray]":
    """uvspec 'output_user lambda uu' rows -> (wl um asc, radiance per um).

    With ``source thermal`` + a band parameterization uvspec reports radiance
    per WAVENUMBER, W m-2 sr-1 (cm-1)-1 (verified against the Planck curve);
    convert to per-micron with |dnu/dlambda| = 1e4 / lambda_um^2.
    """
    data = np.atleast_2d(np.loadtxt(output_file))
    wl_um = data[:, 0] / 1000.0
    rad_per_um = data[:, 1] * 1e4 / wl_um**2
    order = np.argsort(wl_um)
    return wl_um[order], rad_per_um[order]


def variants_of(col: dict) -> "list[str]":
    """Simulated variants: clear always, cloudy when ERA5 holds condensate."""
    v = ["clear"]
    if col.get("wc_file") or col.get("ic_file"):
        v.append("cloudy")
    return v


def cmd_run(args) -> int:
    import platform

    import er3t

    cfg = load_config(args.config)
    mpath = manifest_path(cfg, args.year, args.month, args.sat)
    if not mpath.exists():
        sys.exit(f"missing {mpath} — run prep first")
    manifest = json.loads(mpath.read_text())
    tmp = pbt_dir(cfg, args.year, args.month) / "tmp"
    tmp.mkdir(exist_ok=True)
    srf = load_srf(manifest["srf_file"])
    mol_abs = resolve_mol_abs_param(args)

    lrt_cfg_base = er3t.rtm.lrt.get_lrt_cfg()
    lrt_cfg_base["mol_abs_param"] = f"reptran {mol_abs}"
    lrt_cfg_base["number_of_streams"] = args.streams
    lrt_cfg_base["solar_file"] = None

    inits, to_run = [], []
    for col in manifest["columns"]:
        for vname in variants_of(col):
            tag = f"{col['label']}_{vname}"
            wc = col.get("wc_file") if vname == "cloudy" else None
            ic = col.get("ic_file") if vname == "cloudy" else None
            init = radiance_init(tmp, tag, manifest, col, col["atm_file"],
                                 wc, ic, col["skt"], manifest["emissivity"],
                                 lrt_cfg_base, er3t)
            inits.append((col, vname, init))
            if args.overwrite or not (os.path.exists(init.output_file)
                                      and os.path.getsize(init.output_file) > 50):
                to_run.append(init)

    if to_run:
        workers = args.workers or max((os.cpu_count() or 2) - 2, 1)
        print(f"running {len(to_run)} uvspec job(s) (reptran "
              f"{mol_abs}, {args.streams} streams, "
              f"{workers} workers) ...")
        t0 = time.time()
        er3t.rtm.lrt.lrt_run_mp(to_run, Ncpu=min(workers, len(to_run)))
        print(f"uvspec done in {time.time() - t0:.0f} s "
              f"({(time.time() - t0) / len(to_run):.0f} s/job)")
    else:
        print("all uvspec outputs present (use --overwrite to rerun)")

    results = {}
    for col, vname, init in inits:
        wl_um, rad = read_spectrum(init.output_file)
        np.savez(tmp / f"spectrum_{col['label']}_{vname}.npz",
                 wl_um=wl_um, rad_per_um=rad)
        entry = {"spectrum": f"tmp/spectrum_{col['label']}_{vname}.npz",
                 "bt": {}, "rad": {}}
        for sc_str in col["obs_bt"]:
            rad_c, bt_c = channel_bt(srf, int(sc_str), wl_um, rad)
            entry["bt"][sc_str] = [None if not np.isfinite(v)
                                   else round(float(v), 3) for v in bt_c]
            entry["rad"][sc_str] = [None if not np.isfinite(v)
                                    else float(v) for v in rad_c]
        results.setdefault(col["label"], {})[vname] = entry

        # sanity line: window (10-12 um) BT vs skt for the matching variant
        sc0 = sorted(col["obs_bt"])[0]
        win = [c for c in good_channels(srf, int(sc0))
               if 10.0 <= srf["mean_wl"][c, int(sc0)] <= 12.0]
        bt0 = np.array([np.nan if results[col["label"]][vname]["bt"][sc0][c]
                        is None else results[col["label"]][vname]["bt"][sc0][c]
                        for c in win])
        print(f"  {col['label']:24s} {vname:6s} window BT "
              f"{np.nanmean(bt0):6.1f} K (skt {col['skt']:.1f}, "
              f"sky {col['sky']})")

    rpath = results_path(cfg, args.year, args.month, args.sat)
    merged = json.loads(rpath.read_text()) if rpath.exists() else {}
    for k, v in results.items():
        for vname, entry in v.items():
            # merge at the variant level so rrtmg entries survive a rerun
            merged.setdefault(k, {}).setdefault(vname, {}).update(entry)
    rpath.write_text(json.dumps(merged, indent=1))
    print(f"wrote {rpath}")
    return 0


# ---------------------------------------------------------------- rrtmg (era5)

RRTMG_BANDS_CM = [(10, 350), (350, 500), (500, 630), (630, 700), (700, 820),
                  (820, 980), (980, 1080), (1080, 1180), (1180, 1390),
                  (1390, 1480), (1480, 1800), (1800, 2080), (2080, 2250),
                  (2250, 2380), (2380, 2600), (2600, 3250)]

H_PLANCK, C_LIGHT, K_BOLTZ = 6.62607015e-34, 2.99792458e8, 1.380649e-23


def planck_band_flux(T: float, nu1_cm: float, nu2_cm: float) -> float:
    """pi * integral of B_nu over [nu1, nu2] cm-1 -> W m-2."""
    nu = np.linspace(nu1_cm, nu2_cm, 200) * 100.0       # m-1
    x = H_PLANCK * C_LIGHT * nu / (K_BOLTZ * T)
    b = 2.0 * H_PLANCK * C_LIGHT**2 * nu**3 / np.expm1(x)   # W m-2 sr-1 / m-1
    return float(np.pi * TRAPZ(b, nu))                  # W m-2


def band_flux_to_bt(flux: float, nu1_cm: float, nu2_cm: float) -> float:
    """Flux-equivalent brightness temperature of an RRTMG band (bisection)."""
    lo, hi = 100.0, 400.0
    if not np.isfinite(flux) or flux <= 0:
        return np.nan
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if planck_band_flux(mid, nu1_cm, nu2_cm) < flux:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def cmd_rrtmg(args) -> int:
    try:
        import climlab
        from climlab.domain.axis import Axis
    except ImportError:
        sys.exit("climlab not importable — run under the `era5` conda env")
    from era5_rrtmg_sim import (_BOTTOM_AIR_T, DGE_FROM_REFF, N2O_PPM,
                                cloud_layer_paths, layer_mean,
                                patch_interface_temperature, q_from_vmr)

    cfg = load_config(args.config)
    mpath = manifest_path(cfg, args.year, args.month, args.sat)
    manifest = json.loads(mpath.read_text())
    patch_interface_temperature()

    results = {}
    for col in manifest["columns"]:
        rows = np.loadtxt(col["atm_file"], comments="#")
        z, p, T, air, o3, o2, h2o, co2, _no2 = rows.T
        ch4 = np.atleast_2d(np.loadtxt(col["ch4_file"], comments="#"))[:, 1]
        for vname in variants_of(col):
            state = climlab.column_state(lev=Axis(axis_type="lev", bounds=p))
            state.Ts[:] = col["skt"]
            state.Tatm[:] = layer_mean(T)
            _BOTTOM_AIR_T["value"] = T[-1]
            absorber_vmr = {
                "O3": layer_mean(o3 / air), "CO2": layer_mean(co2 / air),
                "CH4": layer_mean(ch4 / air), "O2": layer_mean(o2 / air),
                "N2O": N2O_PPM * 1e-6,
                "CCL4": 0.0, "CFC11": 0.0, "CFC12": 0.0, "CFC22": 0.0,
            }
            kwargs = {"icld": 0}
            if vname == "cloudy":
                clwp = (cloud_layer_paths(col["wc_file"], z)
                        if col.get("wc_file") else np.zeros(p.size - 1))
                ciwp = (cloud_layer_paths(col["ic_file"], z)
                        if col.get("ic_file") else np.zeros(p.size - 1))
                reff = manifest["reff_um"]
                kwargs = {
                    "icld": 1, "inflglw": 2, "liqflglw": 1, "iceflglw": 3,
                    "cldfrac": ((clwp + ciwp) > 0).astype(float),
                    "clwp": clwp, "ciwp": ciwp,
                    "r_liq": np.full(clwp.size, reff["liquid"]),
                    "r_ice": np.full(ciwp.size, DGE_FROM_REFF * reff["ice"]),
                }
            rad = climlab.radiation.RRTMG_LW(
                state=state,
                specific_humidity=q_from_vmr(layer_mean(h2o / air)),
                emissivity=manifest["emissivity"],
                absorber_vmr=absorber_vmr,
                return_spectral_olr=True, **kwargs)
            rad.compute_diagnostics()
            olr_band = np.squeeze(np.asarray(rad.OLR_spectral))    # (16,)
            band_bt = [band_flux_to_bt(f, *RRTMG_BANDS_CM[i])
                       for i, f in enumerate(olr_band)]
            results.setdefault(col["label"], {})[vname] = {
                "olr_band": [float(f) for f in olr_band],
                "band_bt": [None if not np.isfinite(b) else round(b, 3)
                            for b in band_bt],
                "olr": float(np.squeeze(np.asarray(rad.OLR))),
            }
            print(f"  {col['label']:24s} {vname:6s} OLR "
                  f"{float(np.squeeze(np.asarray(rad.OLR))):7.2f} W/m2")

    rpath = results_path(cfg, args.year, args.month, args.sat)
    merged = json.loads(rpath.read_text()) if rpath.exists() else {}
    for label, per_variant in results.items():
        for vname, entry in per_variant.items():
            merged.setdefault(label, {}).setdefault(vname, {})["rrtmg"] = entry
    rpath.write_text(json.dumps(merged, indent=1))
    print(f"wrote {rpath} (rrtmg band entries)")
    return 0


# ---------------------------------------------------------------- jacobian

def era5_row_index(atm_file: str) -> np.ndarray:
    """Row indices (TOA-first file) belonging to the ERA5 part of the profile.

    The file is afglsw rows on top, ERA5 rows below; the ERA5 block starts at
    the first pressure >= 0.9 hPa jump discontinuity — simpler: rows below
    the afglsw/ERA5 splice are contiguous at the end, identified by pressure
    increasing monotonically past the last afglsw pressure. We tag them as
    all rows with p >= p of the first row whose spacing drops under 2 km.
    """
    rows = np.loadtxt(atm_file, comments="#")
    z = rows[:, 0]
    dz = -np.diff(z)
    i0 = int(np.argmax(dz < 2.0))       # first fine-spaced gap = ERA5 block
    return np.arange(i0, rows.shape[0])


def build_state_list(col: dict, args) -> "list[dict]":
    """Ordered perturbation states for one column (finite forward differences)."""
    states = [{"name": "skt", "kind": "skt", "dx": DT_K, "units": "K"}]
    rows = np.loadtxt(col["atm_file"], comments="#")
    idx = era5_row_index(col["atm_file"])
    for i in idx:
        states.append({"name": f"t_{rows[i, 1]:.0f}hPa", "kind": "t_level",
                       "row": int(i), "p_hpa": float(rows[i, 1]),
                       "dx": DT_K, "units": "K"})
    for i in idx:
        states.append({"name": f"q_{rows[i, 1]:.0f}hPa", "kind": "q_level",
                       "row": int(i), "p_hpa": float(rows[i, 1]),
                       "dx": DQ_FRAC, "units": "frac"})
    if col.get("wc_file"):
        states += [{"name": "ln_lwp", "kind": "wp", "file": "wc_file",
                    "dx": DWP_FRAC, "units": "frac"},
                   {"name": "reff_liq", "kind": "reff", "file": "wc_file",
                    "dx": DREFF_UM, "units": "um"}]
    if col.get("ic_file"):
        states += [{"name": "ln_iwp", "kind": "wp", "file": "ic_file",
                    "dx": DWP_FRAC, "units": "frac"},
                   {"name": "reff_ice", "kind": "reff", "file": "ic_file",
                    "dx": DREFF_UM, "units": "um"}]
    if col.get("wc_file") or col.get("ic_file"):
        states.append({"name": "cth", "kind": "cth", "dx": 1.0,
                       "units": "layer"})
    states.append({"name": "emissivity", "kind": "emis", "dx": DEMIS,
                   "units": "1"})
    if args.states:
        keep = set(args.states)
        states = [s for s in states
                  if s["kind"] in keep or s["name"] in keep]
    return states


def perturbed_files(col: dict, state: dict, jdir: Path) -> dict:
    """Write perturbed copies of the column's input files; return overrides."""
    over = {"atm_file": col["atm_file"], "wc_file": col.get("wc_file"),
            "ic_file": col.get("ic_file"), "skt": col["skt"],
            "emissivity": None}
    jdir.mkdir(parents=True, exist_ok=True)

    def shift_cloud(path_in, path_out):
        rows = np.atleast_2d(np.loadtxt(path_in, comments="#"))
        wc = rows[:, 1].copy()
        nz = np.where(wc > 0)[0]
        wc_new = np.zeros_like(wc)
        for i in nz:                    # rows are TOA-first: up = index - 1
            wc_new[max(i - 1, 0)] += wc[i]
        rows[:, 1] = wc_new
        np.savetxt(path_out, rows, fmt="%9.3f %11.5e %7.2f",
                   header="z(km)     wc(g/m^3)  reff(um)")
        return str(path_out)

    kind = state["kind"]
    if kind == "skt":
        over["skt"] = col["skt"] + state["dx"]
    elif kind == "emis":
        over["emissivity"] = state["dx"]        # applied as ems + dx by caller
    elif kind in ("t_level", "q_level"):
        rows = np.loadtxt(col["atm_file"], comments="#")
        if kind == "t_level":
            rows[state["row"], 2] += state["dx"]
        else:
            rows[state["row"], 6] *= (1.0 + state["dx"])
        out = jdir / f"atm_{state['name']}.dat"
        np.savetxt(out, rows,
                   fmt="%11.3f %11.5f %11.3f %12.6e %12.6e %12.6e %12.6e "
                       "%12.6e %12.6e",
                   header=("perturbed ERA5 profile\n     z(km)      p(mb)   "
                           "     T(K)    air(cm-3)    o3(cm-3)     o2(cm-3)"
                           "    h2o(cm-3)    co2(cm-3)     no2(cm-3)"))
        over["atm_file"] = str(out)
    elif kind in ("wp", "reff"):
        src = col[state["file"]]
        rows = np.atleast_2d(np.loadtxt(src, comments="#"))
        if kind == "wp":
            rows[:, 1] *= (1.0 + state["dx"])
        else:
            rows[rows[:, 1] > 0, 2] += state["dx"]
        out = jdir / f"{Path(src).stem}_{state['name']}.dat"
        np.savetxt(out, rows, fmt="%9.3f %11.5e %7.2f",
                   header="z(km)     wc(g/m^3)  reff(um)")
        over[state["file"]] = str(out)
    elif kind == "cth":
        for key in ("wc_file", "ic_file"):
            if col.get(key):
                out = jdir / f"{Path(col[key]).stem}_cth.dat"
                over[key] = shift_cloud(col[key], out)
    return over


def cmd_jacobian(args) -> int:
    cfg = load_config(args.config)
    mpath = manifest_path(cfg, args.year, args.month, args.sat)
    manifest = json.loads(mpath.read_text())
    columns = [c for c in manifest["columns"]
               if not args.columns or c["label"] in args.columns]
    if args.simulator == "lrt":
        return jacobian_lrt(cfg, args, manifest, columns)
    return jacobian_rrtmg(cfg, args, manifest, columns)


def jacobian_variant(col: dict) -> str:
    """Jacobians are computed around the column's actual sky state."""
    return "cloudy" if (col.get("wc_file") or col.get("ic_file")) else "clear"


def save_jacobian(cfg, args, manifest, col, states, bt0, bt_pert,
                  extra_coords) -> None:
    """bt0 (channel,), bt_pert (state, channel) -> netCDF K matrix."""
    import xarray as xr

    k = (bt_pert - bt0[None, :]) / np.array([s["dx"] for s in states])[:, None]
    ds = xr.Dataset(
        {"bt0": (("channel",), bt0),
         "K": (("state", "channel"), k)},
        coords={
            "channel": np.arange(bt0.size),
            "state": [s["name"] for s in states],
            "state_kind": ("state", [s["kind"] for s in states]),
            "state_dx": ("state", [s["dx"] for s in states]),
            "state_units": ("state", [s["units"] for s in states]),
            "state_p_hpa": ("state", [s.get("p_hpa", np.nan) for s in states]),
            **extra_coords,
        },
        attrs={"label": col["label"], "sky": col["sky"],
               "variant": jacobian_variant(col),
               "simulator": args.simulator, "sat": manifest["sat"],
               "skt": col["skt"], "lat": col["lat"], "lon": col["lon"],
               "note": ("K = dBT/dx, forward finite differences; fractional "
                        "states (q, wp) are per unit fractional change")})
    out = jacobian_path(cfg, args.year, args.month, manifest["sat"],
                        col["label"], args.simulator)
    ds.to_netcdf(out)
    print(f"  wrote {out}")


def jacobian_lrt(cfg, args, manifest, columns) -> int:
    import er3t

    srf = load_srf(manifest["srf_file"])
    tmp = pbt_dir(cfg, args.year, args.month) / "tmp"
    lrt_cfg_base = er3t.rtm.lrt.get_lrt_cfg()
    lrt_cfg_base["mol_abs_param"] = f"reptran {resolve_mol_abs_param(args)}"
    lrt_cfg_base["number_of_streams"] = args.streams
    lrt_cfg_base["solar_file"] = None

    for col in columns:
        vname = jacobian_variant(col)
        states = build_state_list(col, args)
        jdir = tmp / "jac" / col["label"]
        sc = int(sorted(col["obs_bt"])[0])      # reference scene for K
        print(f"{col['label']} ({vname}): {len(states)} states, scene {sc}")

        inits, to_run = [], []
        jobs = [("base", None)] + [(s["name"], s) for s in states]
        for name, state in jobs:
            over = (perturbed_files(col, state, jdir) if state else
                    {"atm_file": col["atm_file"],
                     "wc_file": col.get("wc_file"),
                     "ic_file": col.get("ic_file"), "skt": col["skt"],
                     "emissivity": None})
            ems = manifest["emissivity"] + (over["emissivity"] or 0.0)
            wc = over["wc_file"] if vname == "cloudy" else None
            ic = over["ic_file"] if vname == "cloudy" else None
            init = radiance_init(jdir, name, manifest, col, over["atm_file"],
                                 wc, ic, over["skt"], ems, lrt_cfg_base, er3t)
            inits.append((name, init))
            if args.overwrite or not (os.path.exists(init.output_file)
                                      and os.path.getsize(init.output_file) > 50):
                to_run.append(init)

        if to_run:
            workers = args.workers or max((os.cpu_count() or 2) - 2, 1)
            print(f"  running {len(to_run)} uvspec job(s) ...")
            t0 = time.time()
            er3t.rtm.lrt.lrt_run_mp(to_run, Ncpu=min(workers, len(to_run)))
            print(f"  done in {time.time() - t0:.0f} s")

        bt = {}
        for name, init in inits:
            wl_um, rad = read_spectrum(init.output_file)
            _, bt_c = channel_bt(srf, sc, wl_um, rad)
            bt[name] = bt_c
        bt_pert = np.array([bt[s["name"]] for s in states])
        save_jacobian(cfg, args, manifest, col, states, bt["base"], bt_pert,
                      {"channel_wavelen":
                       ("channel", srf["mean_wl"][:, sc]),
                       "channel_nedr": ("channel", srf["nedr"][:, sc])})
    return 0


def jacobian_rrtmg(cfg, args, manifest, columns) -> int:
    try:
        import climlab
        from climlab.domain.axis import Axis
    except ImportError:
        sys.exit("climlab not importable — run under the `era5` conda env")
    from era5_rrtmg_sim import (_BOTTOM_AIR_T, DGE_FROM_REFF, N2O_PPM,
                                cloud_layer_paths, layer_mean,
                                patch_interface_temperature, q_from_vmr)

    patch_interface_temperature()
    tmp = pbt_dir(cfg, args.year, args.month) / "tmp"

    def band_bt_of(atm_file, wc_file, ic_file, skt, ems, cloudy, reff):
        rows = np.loadtxt(atm_file, comments="#")
        z, p, T, air, o3, o2, h2o, co2, _no2 = rows.T
        # CH4 file is never perturbed; reuse the baseline
        ch4 = np.atleast_2d(np.loadtxt(ch4_file, comments="#"))[:, 1]
        state = climlab.column_state(lev=Axis(axis_type="lev", bounds=p))
        state.Ts[:] = skt
        state.Tatm[:] = layer_mean(T)
        _BOTTOM_AIR_T["value"] = T[-1]
        kwargs = {"icld": 0}
        if cloudy:
            clwp = (cloud_layer_paths(wc_file, z) if wc_file
                    else np.zeros(p.size - 1))
            ciwp = (cloud_layer_paths(ic_file, z) if ic_file
                    else np.zeros(p.size - 1))
            r_liq, r_ice = reff
            kwargs = {"icld": 1, "inflglw": 2, "liqflglw": 1, "iceflglw": 3,
                      "cldfrac": ((clwp + ciwp) > 0).astype(float),
                      "clwp": clwp, "ciwp": ciwp,
                      "r_liq": np.full(clwp.size, r_liq),
                      "r_ice": np.full(ciwp.size, DGE_FROM_REFF * r_ice)}
        rad = climlab.radiation.RRTMG_LW(
            state=state, specific_humidity=q_from_vmr(layer_mean(h2o / air)),
            emissivity=ems,
            absorber_vmr={"O3": layer_mean(o3 / air),
                          "CO2": layer_mean(co2 / air),
                          "CH4": layer_mean(ch4 / air),
                          "O2": layer_mean(o2 / air), "N2O": N2O_PPM * 1e-6,
                          "CCL4": 0.0, "CFC11": 0.0, "CFC12": 0.0,
                          "CFC22": 0.0},
            return_spectral_olr=True, **kwargs)
        rad.compute_diagnostics()
        olr_band = np.squeeze(np.asarray(rad.OLR_spectral))
        return np.array([band_flux_to_bt(f, *RRTMG_BANDS_CM[i])
                         for i, f in enumerate(olr_band)])

    for col in columns:
        vname = jacobian_variant(col)
        cloudy = vname == "cloudy"
        states = build_state_list(col, args)
        jdir = tmp / "jac" / col["label"]
        ch4_file = col["ch4_file"]
        reff0 = (manifest["reff_um"]["liquid"], manifest["reff_um"]["ice"])
        print(f"{col['label']} ({vname}): {len(states)} states [rrtmg]")

        t0 = time.time()
        bt0 = band_bt_of(col["atm_file"], col.get("wc_file"),
                         col.get("ic_file"), col["skt"],
                         manifest["emissivity"], cloudy, reff0)
        bt_pert = []
        for state in states:
            over = perturbed_files(col, state, jdir)
            ems = manifest["emissivity"] + (over["emissivity"] or 0.0)
            reff = list(reff0)
            # r_eff perturbations live in the cloud-file 3rd column for
            # uvspec; RRTMG takes them as arguments instead
            if state["kind"] == "reff":
                if state["file"] == "wc_file":
                    reff[0] += state["dx"]
                else:
                    reff[1] += state["dx"]
                over[state["file"]] = col[state["file"]]
            bt_pert.append(band_bt_of(over["atm_file"], over["wc_file"],
                                      over["ic_file"], over["skt"], ems,
                                      cloudy, tuple(reff)))
        print(f"  {len(states) + 1} RRTMG columns in {time.time() - t0:.1f} s")
        centers = [1e4 / (0.5 * (b[0] + b[1])) for b in RRTMG_BANDS_CM]
        save_jacobian(cfg, args, manifest, col, states, bt0,
                      np.array(bt_pert),
                      {"channel_wavelen": ("channel", np.array(centers))})
    return 0


# ---------------------------------------------------------------- cotscan

def cmd_cotscan(args) -> int:
    """TOA channel BT + dBT/dtau + dBT/dreff vs synthetic-cloud COT.

    Reproduces the ARCSIX BT-vs-COT sensitivity figure with PREFIRE
    channels: a single-layer cloud of given phase / r_eff / top / base is
    inserted into one collocated ERA5 column, its 550-nm optical thickness
    swept over a log grid via ``{wc,ic}_modify tau set``, and each state is
    one spectral uvspec run convolved to channel BT. Jacobians are forward
    finite differences (tau +5 %, r_eff +1 um at fixed tau).
    """
    import er3t

    cfg = load_config(args.config)
    manifest = json.loads(manifest_path(cfg, args.year, args.month,
                                        args.sat).read_text())
    cols = manifest["columns"]
    if args.column:
        col = next(c for c in cols if c["label"] == args.column)
    else:
        col = next((c for c in cols if c["sky"] == "clear"), cols[0])
    sc = int(sorted(col["obs_bt"])[0])
    srf = load_srf(manifest["srf_file"])
    tmp = pbt_dir(cfg, args.year, args.month) / "tmp" / "cotscan"
    tmp.mkdir(parents=True, exist_ok=True)
    mol_abs = resolve_mol_abs_param(args)

    cots = np.logspace(np.log10(args.cot_min), np.log10(args.cot_max),
                       args.n_cot)

    def cloud_file(reff: float) -> str:
        # TOA-first rows; the CBH row's content fills the layer up to CTH.
        # The water content is arbitrary: "tau set" rescales the layer.
        path = tmp / f"cld_{args.phase}_reff{reff:05.2f}.dat"
        path.write_text(
            "#   z(km)     wc(g/m^3)  reff(um)\n"
            f"{args.cth:9.3f} {0.0:11.5e} {reff:7.2f}\n"
            f"{args.cbh:9.3f} {0.01:11.5e} {reff:7.2f}\n")
        return str(path)

    lrt_cfg_base = er3t.rtm.lrt.get_lrt_cfg()
    lrt_cfg_base["mol_abs_param"] = f"reptran {mol_abs}"
    lrt_cfg_base["number_of_streams"] = args.streams
    lrt_cfg_base["solar_file"] = None
    mod_key = "ic_modify" if args.phase == "ice" else "wc_modify"

    inits, to_run = [], []
    for i, cot in enumerate(cots):
        for kind, tau, reff in (
                ("base", cot, args.cer),
                ("dtau", cot * (1.0 + DTAU_FRAC), args.cer),
                ("dreff", cot, args.cer + DREFF_UM)):
            tag = f"{i:02d}_{kind}"
            cld = cloud_file(reff)
            wc = cld if args.phase == "liquid" else None
            ic = cld if args.phase == "ice" else None
            init = radiance_init(
                tmp, tag, manifest, col, col["atm_file"], wc, ic,
                col["skt"], manifest["emissivity"], lrt_cfg_base, er3t,
                extra_opts={mod_key: f"tau set {tau:.5f}"})
            inits.append((i, kind, init))
            if args.overwrite or not (os.path.exists(init.output_file)
                                      and os.path.getsize(init.output_file) > 50):
                to_run.append(init)

    if to_run:
        workers = args.workers or max((os.cpu_count() or 2) - 2, 1)
        print(f"running {len(to_run)} uvspec job(s) (reptran {mol_abs}, "
              f"{args.phase} cloud, scene {sc}) ...")
        t0 = time.time()
        er3t.rtm.lrt.lrt_run_mp(to_run, Ncpu=min(workers, len(to_run)))
        print(f"uvspec done in {time.time() - t0:.0f} s")

    bt = {k: np.full((cots.size, N_CHANNEL), np.nan)
          for k in ("base", "dtau", "dreff")}
    for i, kind, init in inits:
        wl_um, rad = read_spectrum(init.output_file)
        _, bt[kind][i] = channel_bt(srf, sc, wl_um, rad)
    k_tau = (bt["dtau"] - bt["base"]) / (DTAU_FRAC * cots)[:, None]
    k_reff = (bt["dreff"] - bt["base"]) / DREFF_UM

    stem = (f"prefire_cotscan_{args.phase}_cer{args.cer:g}"
            f"_cth{args.cth:g}_{col['label']}_sat{args.sat}")
    npz = pbt_dir(cfg, args.year, args.month) / f"{stem}.npz"
    np.savez(npz, cots=cots, bt=bt["base"], k_tau=k_tau, k_reff=k_reff,
             mean_wl=srf["mean_wl"][:, sc], scene=sc)
    print(f"wrote {npz}")

    # figure: BT / dBT_dtau / dBT_dreff vs COT for a channel subset
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from era5lib.plotstyle import apply_agu_style, panel_label
    apply_agu_style()

    good = good_channels(srf, sc)
    wl_good = srf["mean_wl"][good, sc]
    targets = args.wavelengths or [4.8, 8.4, 11.4, 17.3, 22.4, 26.5]
    chs = []
    for t in targets:
        c = int(good[np.argmin(np.abs(wl_good - t))])
        if c not in chs:
            chs.append(c)
    colors = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7",
              "#56B4E9", "#000000"]

    fig, axes = plt.subplots(1, 3, figsize=(13.4, 3.9))
    for j, c in enumerate(chs):
        wl_c = srf["mean_wl"][c, sc]
        lab = f"{wl_c:.1f} $\\mu$m ({1e4 / wl_c:.0f} cm$^{{-1}}$)"
        kw = {"color": colors[j % len(colors)], "lw": 1.5}
        axes[0].plot(cots, bt["base"][:, c], label=lab, **kw)
        axes[1].plot(cots, k_tau[:, c], **kw)
        axes[2].plot(cots, k_reff[:, c], **kw)
    for k, (ax, ylab) in enumerate(zip(axes, (
            "simulated TOA BT (K)",
            "dBT/d$\\tau$ (K per unit $\\tau$)",
            "dBT/dr$_{\\rm eff}$ (K $\\mu$m$^{-1}$)"))):
        ax.set_xscale("log")
        ax.set_xlabel("cloud optical thickness (550 nm)")
        ax.set_ylabel(ylab)
        if k > 0:
            ax.axhline(0, color="k", lw=0.8, ls="--")
        panel_label(ax, chr(97 + k))
    axes[0].legend(fontsize=6.4, loc="best")
    fig.suptitle(
        f"TIRS{args.sat} scene {sc} — {args.phase} cloud, "
        f"r$_{{\\rm eff}}$={args.cer:g} $\\mu$m, CTH={args.cth:g} km, "
        f"CBH={args.cbh:g} km — atm/sfc: {col['label']} "
        f"(skt {col['skt']:.1f} K)", y=1.02, fontsize=9)
    fig.tight_layout()
    fpath = figures_dir(cfg) / f"{stem}.png"
    fig.savefig(fpath, dpi=300, bbox_inches="tight")
    print(f"wrote {fpath}")
    return 0


# ---------------------------------------------------------------- figure

def cmd_figure(args) -> int:
    import matplotlib.pyplot as plt

    from era5lib.plotstyle import apply_agu_style, panel_label

    cfg = load_config(args.config)
    manifest = json.loads(manifest_path(cfg, args.year, args.month,
                                        args.sat).read_text())
    rpath = results_path(cfg, args.year, args.month, args.sat)
    if not rpath.exists():
        sys.exit(f"missing {rpath} — run `run` (and optionally `rrtmg`) first")
    results = json.loads(rpath.read_text())
    srf = load_srf(manifest["srf_file"])
    apply_agu_style()

    cols = manifest["columns"]
    n = len(cols)
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.9 * ncol, 3.6 * nrow),
                             squeeze=False)
    C_OBS, C_LIB, C_RRT = "#000000", "#0072B2", "#009E73"

    for k, col in enumerate(cols):
        ax = axes[k // ncol][k % ncol]
        vname = jacobian_variant(col)
        res = results.get(col["label"], {}).get(vname, {})
        for i, (sc_str, obs) in enumerate(sorted(col["obs_bt"].items())):
            sc = int(sc_str)
            wl = srf["mean_wl"][:, sc]
            obs = np.array([np.nan if v is None else v for v in obs])
            ax.plot(wl, obs, "o", ms=3.4, mfc="none", color=C_OBS,
                    label="PREFIRE obs" if i == 0 else None)
            bt = res.get("bt", {}).get(sc_str)
            if bt:
                bt = np.array([np.nan if v is None else v for v in bt])
                ax.plot(wl, bt, "s", ms=2.8, color=C_LIB,
                        label="libRadtran sim" if i == 0 else None)
        rrt = res.get("rrtmg")
        if rrt:
            for i, ((nu1, nu2), b) in enumerate(zip(RRTMG_BANDS_CM,
                                                    rrt["band_bt"])):
                if b is None:
                    continue
                ax.plot([1e4 / nu2, 1e4 / nu1], [b, b], "-", lw=1.4,
                        color=C_RRT,
                        label="RRTMG band (flux-equiv.)" if i == 0 else None)
        ax.axhline(col["skt"], color="0.6", lw=0.8, ls=":",
                   label="ERA5 skt")
        ax.set_xlim(3.5, 30.0)
        ax.set_xlabel("wavelength (µm)")
        ax.set_ylabel("brightness temperature (K)")
        panel_label(ax, chr(97 + k))
        ax.set_title(f"{col['label']} — {col['sky']}", fontsize=8)
        if k == 0:
            ax.legend(loc="lower right", fontsize=6.5)
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.tight_layout()
    fdir = figures_dir(cfg)
    fpath = fdir / (f"prefire_bt_{args.year:04d}{args.month:02d}"
                    f"_sat{args.sat}.png")
    fig.savefig(fpath, dpi=300, bbox_inches="tight")
    print(f"wrote {fpath}")

    # Jacobian heatmaps, one figure per column that has a K file
    import xarray as xr
    for col in cols:
        for sim in ("lrt", "rrtmg"):
            jpath = jacobian_path(cfg, args.year, args.month, args.sat,
                                  col["label"], sim)
            if not jpath.exists():
                continue
            ds = xr.open_dataset(jpath)
            fig2 = jacobian_figure(ds, col, sim)
            fp2 = fdir / (f"prefire_jacobian_{sim}_{col['label']}"
                          f"_sat{args.sat}.png")
            fig2.savefig(fp2, dpi=300, bbox_inches="tight")
            print(f"wrote {fp2}")
    return 0


def jacobian_figure(ds, col, sim):
    import matplotlib.pyplot as plt

    from era5lib.plotstyle import panel_label

    kinds = ds["state_kind"].values.astype(str)
    wl = ds["channel_wavelen"].values
    ok = np.isfinite(ds["bt0"].values)
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.2),
                             gridspec_kw={"width_ratios": [3, 3, 1.6]})
    for ax, kind, title, cmap in ((axes[0], "t_level",
                                   "K$_T$ (K BT / K)", "RdBu_r"),
                                  (axes[1], "q_level",
                                   "K$_q$ (K BT / frac q)", "RdBu_r")):
        sel = kinds == kind
        if not sel.any():
            ax.axis("off")
            continue
        K = ds["K"].values[sel][:, ok]
        p = ds["state_p_hpa"].values[sel]
        vmax = np.nanmax(np.abs(K)) or 1.0
        m = ax.pcolormesh(wl[ok], p, K, cmap=cmap, vmin=-vmax, vmax=vmax,
                          shading="nearest")
        ax.invert_yaxis()
        ax.set_xlabel("channel wavelength (µm)")
        ax.set_ylabel("perturbed level (hPa)")
        ax.set_title(title, fontsize=9)
        plt.colorbar(m, ax=ax, pad=0.02)
    scalars = [i for i, k in enumerate(kinds)
               if k in ("skt", "wp", "reff", "cth", "emis")]
    ax = axes[2]
    for i in scalars:
        ax.plot(wl[ok], ds["K"].values[i, ok], lw=1.1,
                label=str(ds["state"].values[i]))
    ax.axhline(0, color="0.7", lw=0.6)
    ax.set_xlabel("channel wavelength (µm)")
    ax.set_ylabel("dBT/dx (K per unit)")
    ax.set_title("scalar states", fontsize=9)
    ax.legend(fontsize=6)
    for i, ax in enumerate(axes):
        panel_label(ax, chr(97 + i))
    fig.suptitle(f"{col['label']} ({col['sky']}) — {sim} Jacobian",
                 y=1.02, fontsize=10)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------- main

def add_common(sp):
    sp.add_argument("--year", type=int, required=True)
    sp.add_argument("--month", type=int, required=True)
    sp.add_argument("--sat", type=int, default=1, choices=[1, 2])
    sp.add_argument("--config", default=None)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("collocate", help="footprints -> ERA5 columns")
    add_common(sp)
    sp.add_argument("--n-clear", type=int, default=3)
    sp.add_argument("--n-cloudy", type=int, default=3)
    sp.set_defaults(func=cmd_collocate)

    sp = sub.add_parser("prep", help="profile files + manifest")
    add_common(sp)
    sp.set_defaults(func=cmd_prep)

    sp = sub.add_parser("run", help="baseline spectral radiance (er3t_env)")
    add_common(sp)
    sp.add_argument("--streams", type=int, default=16)
    sp.add_argument("--mol-abs-param", default=None,
                    choices=["coarse", "medium", "fine"],
                    help="default: coarse on Darwin, fine on Linux/CURC")
    sp.add_argument("--force-local", action="store_true",
                    help="override the Darwin OOM guard for medium/fine")
    sp.add_argument("--workers", type=int, default=None)
    sp.add_argument("--overwrite", action="store_true")
    sp.set_defaults(func=cmd_run)

    sp = sub.add_parser("rrtmg", help="RRTMG band BT cross-check (era5 env)")
    add_common(sp)
    sp.set_defaults(func=cmd_rrtmg)

    sp = sub.add_parser("jacobian", help="finite-difference K matrices")
    add_common(sp)
    sp.add_argument("--simulator", default="rrtmg", choices=["rrtmg", "lrt"])
    sp.add_argument("--columns", nargs="*", default=None,
                    help="restrict to these column labels")
    sp.add_argument("--states", nargs="*", default=None,
                    help="restrict to these state kinds/names")
    sp.add_argument("--streams", type=int, default=16)
    sp.add_argument("--mol-abs-param", default=None,
                    choices=["coarse", "medium", "fine"],
                    help="default: coarse on Darwin, fine on Linux/CURC")
    sp.add_argument("--force-local", action="store_true",
                    help="override the Darwin OOM guard for medium/fine")
    sp.add_argument("--workers", type=int, default=None)
    sp.add_argument("--overwrite", action="store_true")
    sp.set_defaults(func=cmd_jacobian)

    sp = sub.add_parser("cotscan",
                        help="BT + Jacobians vs synthetic-cloud COT (er3t_env)")
    add_common(sp)
    sp.add_argument("--column", default=None,
                    help="background column label (default: first clear)")
    sp.add_argument("--phase", default="ice", choices=["ice", "liquid"])
    sp.add_argument("--cer", type=float, default=20.0,
                    help="cloud effective radius (um)")
    sp.add_argument("--cth", type=float, default=10.0,
                    help="cloud-top height (km)")
    sp.add_argument("--cbh", type=float, default=9.9,
                    help="cloud-base height (km)")
    sp.add_argument("--cot-min", type=float, default=0.08)
    sp.add_argument("--cot-max", type=float, default=11.0)
    sp.add_argument("--n-cot", type=int, default=14)
    sp.add_argument("--wavelengths", type=float, nargs="*", default=None,
                    help="channel wavelengths to plot (um)")
    sp.add_argument("--streams", type=int, default=16)
    sp.add_argument("--mol-abs-param", default=None,
                    choices=["coarse", "medium", "fine"],
                    help="default: coarse on Darwin, fine on Linux/CURC")
    sp.add_argument("--force-local", action="store_true",
                    help="override the Darwin OOM guard for medium/fine")
    sp.add_argument("--workers", type=int, default=None)
    sp.add_argument("--overwrite", action="store_true")
    sp.set_defaults(func=cmd_cotscan)

    sp = sub.add_parser("figure", help="BT spectra + Jacobian heatmaps")
    add_common(sp)
    sp.set_defaults(func=cmd_figure)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
