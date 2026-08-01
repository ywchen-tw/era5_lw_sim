#!/usr/bin/env python
"""Broadband LW (thermal) flux simulation from ERA5 profiles via libRadtran/er3t.

Run under the ``er3t_env`` conda env (needs er3t + libRadtran; see README).
Three subcommands, all operating on one time snapshot:

  prep     select N cloud-free pixels spanning the SBI-strength range, build
           libRadtran atmosphere/CH4 profile files from ERA5 (+ CAMS EGG4
           CO2/CH4), write a manifest with the matching ERA5 fluxes
  run      run uvspec (source thermal, integrated 4-100 um flux) per profile
  compare  table + figure of simulated vs ERA5 surface LW fluxes; becomes a
           three-way comparison when era5_rrtmg_sim.py (run under the era5
           env) has added RRTMG-LW rows to the results CSV

Physics choices (see plan/README): surface emission from ERA5 skin temperature
(uvspec ``sur_temperature``) with emissivity 0.99; ERA5 LWdn = strd/3600 and
LWup = (strd-str)/3600 (1-h accumulations, W m-2); simulated range 4-100 um vs
ERA5's RRTMG-LW 3.08-1000 um — the missing far-IR tail (~1% of sigma*T^4 at
Arctic winter temperatures) is estimated analytically in ``compare``.

Examples:
    python src/era5_lrt_sim.py prep    --year 2020 --month 1 --day 1 --hour 12
    python src/era5_lrt_sim.py run     --year 2020 --month 1 --day 1 --hour 12
    python src/era5_lrt_sim.py compare --year 2020 --month 1 --day 1 --hour 12
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from era5lib.config import REPO_ROOT, figures_dir, inversion_path, load_config, plev_path, sfc_path
from era5lib.inversion import column_heights
from era5lib.io_era5 import open_era5

K_B = 1.380649e-23        # J K-1
G0 = 9.80665              # m s-2
SIGMA = 5.670374419e-8    # W m-2 K-4
M_DRY, M_H2O, M_O3, M_CO2, M_CH4 = 28.9647, 18.0153, 47.9982, 44.0095, 16.0425
EMISSIVITY = 0.99
# reptran thermal covers 2500-100000 nm; start at ERA5's RRTMG-LW lower edge
# (3.08 um). The 3.08-4 um band adds only ~0.05 W/m2 at Arctic winter
# temperatures; the remaining mismatch vs ERA5 is the far-IR >100 um tail.
WVL_RANGE_NM = (3080.0, 100000.0)
CO2_CONST_PPM = 415.0
CH4_CONST_PPM = 1.9
# effective radii for cloudy runs (not in ERA5 pressure-level output; ERA5's
# radiation diagnoses them internally — a documented uncertainty of the setup)
REFF_LIQ_UM = 10.0
REFF_ICE_UM = 25.0
R_DRY = 287.04            # J kg-1 K-1

LIBRADTRAN_DIR = os.environ.get(
    "LIBRADTRAN_V2_DIR", "/Users/yuch8913/programming/er3t/libRadtran-2.0.6")


def lw_sim_dir(cfg: dict, date: dt.date) -> Path:
    return inversion_path(cfg, date).parent / "lw_sim"


def sky_tag(sky: str) -> str:
    return "" if sky == "clear" else "_" + sky


def manifest_path(cfg: dict, date: dt.date, hour: int, sky: str = "clear") -> Path:
    return lw_sim_dir(cfg, date) / f"manifest_{date:%Y%m%d}T{hour:02d}Z{sky_tag(sky)}.json"


def results_path(cfg: dict, date: dt.date, hour: int, sky: str = "clear") -> Path:
    return lw_sim_dir(cfg, date) / f"results_{date:%Y%m%d}T{hour:02d}Z{sky_tag(sky)}.csv"


# ---------------------------------------------------------------- gas helpers

def vmr_from_q(q: np.ndarray) -> np.ndarray:
    """Specific humidity (kg/kg) -> H2O volume mixing ratio."""
    w = q / (1.0 - q)
    ratio = w * (M_DRY / M_H2O)
    return ratio / (1.0 + ratio)


def load_cams_profiles(year: int, month: int):
    """Arctic-mean CAMS EGG4 CO2/CH4 vmr profiles vs pressure, or None."""
    import xarray as xr

    from era5_cams_download import cams_path
    path = cams_path(year, month)
    if not path.exists():
        return None
    ds = xr.open_dataset(path).squeeze()
    w = np.cos(np.deg2rad(ds["latitude"]))
    co2 = ds["co2"].weighted(w).mean(("latitude", "longitude")).values
    ch4 = ds["ch4"].weighted(w).mean(("latitude", "longitude")).values
    p = ds["pressure_level"].values.astype(float)
    order = np.argsort(p)   # ascending pressure for np.interp
    return {
        "p_hpa": p[order],
        "x_co2": co2[order] * (M_DRY / M_CO2),
        "x_ch4": ch4[order] * (M_DRY / M_CH4),
        "source": str(path.name),
    }


def gas_vmr_at(p_hpa: np.ndarray, cams: "dict | None", gas: str) -> np.ndarray:
    """CO2/CH4 vmr interpolated in log-p; constants when CAMS is unavailable."""
    const = {"co2": CO2_CONST_PPM * 1e-6, "ch4": CH4_CONST_PPM * 1e-6}[gas]
    if cams is None:
        return np.full(p_hpa.shape, const)
    return np.interp(np.log(p_hpa), np.log(cams["p_hpa"]), cams["x_" + gas])


# ---------------------------------------------------------------- prep

def write_cloud_file(path: Path, z_km: np.ndarray, wc_kgkg: np.ndarray,
                     p_hpa: np.ndarray, T: np.ndarray, reff_um: float) -> "str | None":
    """Write a 1D libRadtran cloud file (z, water content g/m3, r_eff um).

    Rows are TOA-first; the content at a level applies to the layer beneath it
    (libRadtran convention). Returns None when the column holds no condensate.
    """
    rho = (p_hpa * 100.0) / (R_DRY * T)                  # kg m-3
    wc_gm3 = np.maximum(wc_kgkg, 0.0) * rho * 1000.0
    wc_gm3[wc_gm3 < 1e-5] = 0.0
    if not np.any(wc_gm3 > 0):
        return None
    rows = np.column_stack([z_km, wc_gm3, np.full(z_km.shape, reff_um)])[::-1]
    lines = ["{:9.3f} {:11.5e} {:7.2f}".format(*r) for r in rows]
    path.write_text("#   z(km)     wc(g/m^3)  reff(um)\n" + "\n".join(lines) + "\n")
    return str(path)


def pick_cloudy_pixels(clear_ok, lwp_g, iwp_g, n):
    """Pick up to n overcast pixels spanning phase and water path.

    clear_ok: mask of eligible (overcast, full-column) pixels; lwp/iwp in g/m2.
    Categories: liquid- and ice-dominated, mixed-phase, thin and thick total
    path; empty categories are topped up with the largest remaining paths.
    """
    twp = lwp_g + iwp_g
    with np.errstate(invalid="ignore", divide="ignore"):
        liq_frac = np.where(twp > 0, lwp_g / twp, np.nan)
    picks = []

    def grab(tag, mask, score):
        s = np.where(mask & clear_ok, score, np.nan)
        if np.isfinite(s).any():
            iy, ix = np.unravel_index(np.nanargmax(s), s.shape)
            if all(p[1] != (iy, ix) for p in picks):
                picks.append((tag, (iy, ix)))

    grab("liquid", liq_frac > 0.8, lwp_g)
    grab("ice", liq_frac < 0.2, iwp_g)
    grab("mixed", (liq_frac >= 0.4) & (liq_frac <= 0.6), twp)
    grab("thin", (twp > 5) & (twp < 20), -twp)
    grab("thick", np.isfinite(liq_frac), twp)
    while len(picks) < n:          # top up with the next-largest total path
        s = twp.astype(float).copy()
        s[~clear_ok] = np.nan
        for _, (iy, ix) in picks:
            s[iy, ix] = np.nan
        if not np.isfinite(s).any():
            break
        iy, ix = np.unravel_index(np.nanargmax(s), s.shape)
        picks.append((f"extra{len(picks)}", (iy, ix)))
    return picks[:n]


def sample_pixels(mask: np.ndarray, n: int, seed: int):
    """Seeded random sample of eligible pixels, labeled px0000, px0001, ..."""
    iy, ix = np.where(mask)
    rng = np.random.default_rng(seed)
    sel = rng.choice(iy.size, size=min(n, iy.size), replace=False)
    return [("px%04d" % k, (int(iy[s]), int(ix[s]))) for k, s in enumerate(sel)]


def read_afglsw_above(p_top_hpa: float, z_top_km: float) -> np.ndarray:
    """Rows (9 col) of the subarctic-winter standard atmosphere above ERA5's top."""
    path = Path(LIBRADTRAN_DIR) / "data" / "atmmod" / "afglsw.dat"
    rows = np.loadtxt(str(path), comments="#")
    keep = (rows[:, 1] < 0.9 * p_top_hpa) & (rows[:, 0] > z_top_km + 2.0)
    return rows[keep]   # afglsw is TOA-first already


def build_profile_files(out_dir: Path, label: str, p_hpa, T, q, o3_mmr,
                        t2m: float, sp_hpa: float, cams) -> "tuple[str, str]":
    """Write the 9-column atmosphere file and the CH4 mol_file for one pixel.

    p_hpa/T/q/o3_mmr are surface-first (descending pressure) above-ground level
    arrays; the surface point (sp_hpa, t2m) is appended as the lowest row.
    """
    z_km = column_heights(T, p_hpa, t2m, sp_hpa, q) / 1000.0

    # drop levels within 10 m of the surface: at 3-decimal km precision they
    # would duplicate the surface row's height and break uvspec's z grid
    keep = z_km >= 0.01
    p_hpa, T, q, o3_mmr, z_km = (a[keep] for a in (p_hpa, T, q, o3_mmr, z_km))

    # surface-first arrays including the surface row
    p_all = np.concatenate([[sp_hpa], p_hpa])
    t_all = np.concatenate([[t2m], T])
    q_all = np.concatenate([[q[0]], q])
    o3_all = np.concatenate([[o3_mmr[0]], o3_mmr])
    z_all = np.concatenate([[0.0], z_km])

    air = (p_all * 100.0) / (K_B * t_all) * 1e-6           # cm-3
    n_h2o = vmr_from_q(q_all) * air
    n_o3 = o3_all * (M_DRY / M_O3) * air
    n_o2 = 0.2095 * air
    n_co2 = gas_vmr_at(p_all, cams, "co2") * air
    n_ch4 = gas_vmr_at(p_all, cams, "ch4") * air
    n_no2 = np.zeros_like(air)

    era5_rows = np.column_stack(
        [z_all, p_all, t_all, air, n_o3, n_o2, n_h2o, n_co2, n_no2])[::-1]
    upper = read_afglsw_above(float(p_all.min()), float(z_all.max()))
    rows = np.vstack([upper, era5_rows])

    atm_file = out_dir / f"atm_profile_{label}.dat"
    header = ("# ERA5-derived atmospheric profile (upper levels: afglsw)\n"
              "#      z(km)      p(mb)        T(K)    air(cm-3)    o3(cm-3)"
              "     o2(cm-3)    h2o(cm-3)    co2(cm-3)     no2(cm-3)\n")
    lines = ["{:11.3f} {:11.5f} {:11.3f} {:12.6e} {:12.6e} {:12.6e} "
             "{:12.6e} {:12.6e} {:12.6e}".format(*row) for row in rows]
    atm_file.write_text(header + "\n".join(lines) + "\n")

    # CH4 mol_file on the same height grid (upper rows: vmr at ERA5 top x afglsw air)
    x_ch4_top = float(gas_vmr_at(np.array([p_all.min()]), cams, "ch4")[0])
    ch4_rows = np.vstack([
        np.column_stack([upper[:, 0], x_ch4_top * upper[:, 3]]),
        np.column_stack([era5_rows[:, 0], n_ch4[::-1]]),
    ])
    ch4_file = out_dir / f"ch4_profile_{label}.dat"
    ch4_lines = ["{:11.3f} {:12.6e}".format(*row) for row in ch4_rows]
    ch4_file.write_text("# CH4 profile (ERA5 grid, CAMS/constant vmr)\n"
                        "#      z(km)      ch4(cm-3)\n" + "\n".join(ch4_lines) + "\n")
    # keep: which input levels survived the near-surface filter; z_km matches it
    return str(atm_file), str(ch4_file), keep, z_km


def cmd_prep(args) -> int:
    cfg = load_config(args.config)
    date = dt.date(args.year, args.month, args.day)
    when = np.datetime64(f"{date:%Y-%m-%d}T{args.hour:02d}:00")

    plev = open_era5(plev_path(cfg, date)).sel(valid_time=when)
    sfc = open_era5(sfc_path(cfg, date)).sel(valid_time=when)
    inv = open_era5(inversion_path(cfg, date)).sel(valid_time=when)
    plev = plev.transpose("pressure_level", "latitude", "longitude")

    p = plev["pressure_level"].values.astype(float)          # descending
    sp_hpa = sfc["sp"].values / 100.0

    # column condensate paths and cloud-fraction maximum
    cc_max = plev["cc"].max("pressure_level").values
    p_pa_asc = p[::-1] * 100.0
    lwp_g = np.trapz(plev["clwc"].values[::-1], x=p_pa_asc, axis=0) / G0 * 1000.0
    iwp_g = np.trapz(plev["ciwc"].values[::-1], x=p_pa_asc, axis=0) / G0 * 1000.0

    strength = np.nan_to_num(inv["sbi_strength"].values)
    found = inv["sbi_found"].values.astype(bool)
    lat2d, lon2d = np.meshgrid(inv["latitude"].values, inv["longitude"].values,
                               indexing="ij")

    if args.sky == "clear":
        # no cloud fraction, negligible condensate path, full column
        clear = ((cc_max <= 0.01) & (lwp_g + iwp_g <= 1.0) & (sp_hpa >= 1000.0))
        print(f"clear-sky pixels: {clear.sum()} / {clear.size} "
              f"({clear.mean():.1%} of the domain)")
        if clear.sum() < args.n:
            sys.exit("not enough clear pixels — relax thresholds or pick another time")

        if args.n > 5:      # statistics mode: random sample of the population
            picks = sample_pixels(clear, args.n, args.seed)
        else:
            # pick pixels across the clear-sky SBI-strength distribution
            picks = []
            s_found = strength.copy()
            s_found[~clear | ~found] = np.nan
            def nearest_to(target):
                return np.unravel_index(np.nanargmin(np.abs(s_found - target)),
                                        s_found.shape)
            picks.append(("strongest",
                          np.unravel_index(np.nanargmax(s_found), s_found.shape)))
            finite = s_found[np.isfinite(s_found)]
            for tag, pct in (("p75", 75), ("median", 50), ("p25", 25)):
                picks.append((tag, nearest_to(np.percentile(finite, pct))))
            none_mask = clear & ~found
            if none_mask.any():
                cand = np.where(none_mask)
                picks.append(("no_sbi", (cand[0][0], cand[1][0])))
            else:
                picks.append(("weakest",
                              np.unravel_index(np.nanargmin(s_found),
                                               s_found.shape)))
            picks = picks[:args.n]
    else:
        # cloudy: near-overcast columns so the plane-parallel (cloud fraction 1)
        # simulation matches ERA5's all-sky flux as closely as possible
        overcast = (cc_max >= 0.99) & (sp_hpa >= 1000.0) & (lwp_g + iwp_g > 1.0)
        print(f"overcast pixels: {overcast.sum()} / {overcast.size} "
              f"({overcast.mean():.1%} of the domain)")
        if overcast.sum() < args.n:
            sys.exit("not enough overcast pixels — relax thresholds")
        picks = (sample_pixels(overcast, args.n, args.seed) if args.n > 5
                 else pick_cloudy_pixels(overcast, lwp_g, iwp_g, args.n))

    out_dir = lw_sim_dir(cfg, date)
    out_dir.mkdir(parents=True, exist_ok=True)
    cams = None if args.fallback_constants else load_cams_profiles(args.year, args.month)
    if cams is None and not args.fallback_constants:
        sys.exit("no CAMS file found — run src/era5_cams_download.py first, or "
                 "pass --fallback-constants for CO2 415 ppm / CH4 1.9 ppm")

    profiles = []
    for label, (iy, ix) in picks:
        above = p <= sp_hpa[iy, ix]
        skt = float(sfc["skt"].values[iy, ix])
        t2m = float(sfc["t2m"].values[iy, ix])
        p_col = p[above]
        t_col = plev["t"].values[above, iy, ix].astype(float)
        atm_file, ch4_file, keep, z_col = build_profile_files(
            out_dir, label, p_col, t_col,
            plev["q"].values[above, iy, ix].astype(float),
            plev["o3"].values[above, iy, ix].astype(float),
            t2m, float(sp_hpa[iy, ix]), cams)
        strd = float(sfc["strd"].values[iy, ix]) / 3600.0
        str_net = float(sfc["str"].values[iy, ix]) / 3600.0
        prof = {
            "label": label,
            "lat": float(lat2d[iy, ix]), "lon": float(lon2d[iy, ix]),
            "sbi_strength": float(strength[iy, ix]),
            "sbi_found": int(found[iy, ix]),
            "skt": skt, "t2m": t2m,
            "era5_lwdn": strd, "era5_lwup": strd - str_net,
            "atm_file": atm_file, "ch4_file": ch4_file,
        }
        cloud_note = ""
        if args.sky == "cloudy":
            prof["lwp_g"] = float(lwp_g[iy, ix])
            prof["iwp_g"] = float(iwp_g[iy, ix])
            prof["cc_max"] = float(cc_max[iy, ix])
            prof["wc_file"] = write_cloud_file(
                out_dir / f"wc_{label}.dat", z_col,
                plev["clwc"].values[above, iy, ix].astype(float)[keep],
                p_col[keep], t_col[keep], REFF_LIQ_UM)
            prof["ic_file"] = write_cloud_file(
                out_dir / f"ic_{label}.dat", z_col,
                plev["ciwc"].values[above, iy, ix].astype(float)[keep],
                p_col[keep], t_col[keep], REFF_ICE_UM)
            cloud_note = (f"  LWP {prof['lwp_g']:6.1f}  IWP {prof['iwp_g']:6.1f}"
                          f" g/m2")
        profiles.append(prof)
        if len(picks) <= 12:
            print(f"  {label:<10} ({lat2d[iy, ix]:.2f}N, {lon2d[iy, ix]:.2f}E)  "
                  f"SBI {strength[iy, ix]:5.2f} K  skt {skt:6.2f} K  "
                  f"ERA5 LWdn {strd:6.1f}  LWup {strd - str_net:6.1f} W/m2"
                  + cloud_note)
        elif len(profiles) % 100 == 0:
            print(f"  prepared {len(profiles)}/{len(picks)} profiles ...")

    manifest = {
        "snapshot": f"{date:%Y-%m-%d}T{args.hour:02d}:00Z",
        "sky": args.sky,
        "emissivity": EMISSIVITY,
        "wavelength_range_nm": list(WVL_RANGE_NM),
        "reff_um": {"liquid": REFF_LIQ_UM, "ice": REFF_ICE_UM},
        "gas_source": cams["source"] if cams else
                      f"constants: CO2 {CO2_CONST_PPM} ppm, CH4 {CH4_CONST_PPM} ppm",
        "notes": ("ERA5 fluxes are 1-h accumulations ending at the snapshot time; "
                  "sim range vs ERA5 RRTMG-LW 3.08-1000 um; cloudy runs assume "
                  "overcast plane-parallel clouds with fixed effective radii"),
        "profiles": profiles,
    }
    mpath = manifest_path(cfg, date, args.hour, args.sky)
    mpath.write_text(json.dumps(manifest, indent=2))
    print(f"wrote {mpath}")
    return 0


# ---------------------------------------------------------------- run

def cmd_run(args) -> int:
    import copy
    import er3t

    cfg = load_config(args.config)
    date = dt.date(args.year, args.month, args.day)
    mpath = manifest_path(cfg, date, args.hour, args.sky)
    if not mpath.exists():
        sys.exit(f"missing {mpath} — run the prep subcommand first "
                 f"(with --sky {args.sky})")
    manifest = json.loads(mpath.read_text())

    if args.streams is None:
        args.streams = 4 if platform.system() == "Darwin" else 8
    tmp_dir = lw_sim_dir(cfg, date) / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    lrt_cfg_base = er3t.rtm.lrt.get_lrt_cfg()
    lrt_cfg_base["mol_abs_param"] = "reptran " + args.mol_abs_param
    lrt_cfg_base["number_of_streams"] = args.streams
    mute_list = ["albedo", "wavelength", "spline", "source solar",
                 "slit_function_file"]

    tag = sky_tag(args.sky)
    inits, to_run = [], []
    for prof in manifest["profiles"]:
        label = prof["label"]
        in_path = str(tmp_dir / f"input_{label}{tag}.txt")
        out_path = str(tmp_dir / f"output_{label}{tag}.txt")
        lrt_cfg = copy.deepcopy(lrt_cfg_base)
        lrt_cfg["atmosphere_file"] = prof["atm_file"]
        input_dict_extra = {
            "source": "thermal",
            "albedo_add": f"{1.0 - manifest['emissivity']:.4f}",
            "sur_temperature": f"{prof['skt']:.2f}",
            "wavelength_add": "{:.0f} {:.0f}".format(*manifest["wavelength_range_nm"]),
            "output_process": "integrate",
            "mol_file": "CH4 " + prof["ch4_file"],
        }
        # cloudy pixels: separate liquid/ice 1D cloud files with
        # phase-appropriate optical properties (Hu & Stamnes for droplets,
        # Fu for ice crystals — both built into libRadtran, thermal-capable)
        if prof.get("wc_file"):
            input_dict_extra["wc_file 1D"] = prof["wc_file"]
            input_dict_extra["wc_properties"] = "hu interpolate"
        if prof.get("ic_file"):
            input_dict_extra["ic_file 1D"] = prof["ic_file"]
            input_dict_extra["ic_properties"] = "fu interpolate"
        init = er3t.rtm.lrt.lrt_init_mono_flx(
            input_file=in_path,
            output_file=out_path,
            date=dt.datetime(args.year, args.month, args.day, args.hour),
            solar_zenith_angle=80.0,     # irrelevant for source thermal
            Nx=1,
            output_altitude=[0, "toa"],
            input_dict_extra=input_dict_extra,
            mute_list=mute_list,
            lrt_cfg=lrt_cfg,
            cld_cfg=None,
            aer_cfg=None,
        )
        inits.append((label, init))
        if args.overwrite or not (os.path.exists(out_path)
                                  and os.path.getsize(out_path) > 50):
            to_run.append(init)

    if to_run:
        workers = args.workers or max((os.cpu_count() or 2) - 2, 1)
        print(f"running {len(to_run)} uvspec thermal job(s) "
              f"(reptran {args.mol_abs_param}, {args.streams} streams, "
              f"{workers} workers) ...")
        er3t.rtm.lrt.lrt_run_mp(to_run, Ncpu=min(workers, len(to_run)))
    else:
        print("all uvspec outputs present; reading existing results "
              "(use --overwrite to rerun)")

    rows = []
    for label, init in inits:
        data = er3t.rtm.lrt.lrt_read_uvspec_flx([init])
        # er3t divides uvspec output by 1000 (mW solar-file convention);
        # thermal-source output is natively W m-2, so undo that here
        f_dn = np.squeeze(data.f_down) * 1000.0     # [surface, toa]
        f_up = np.squeeze(data.f_up) * 1000.0
        rows.append({
            "label": label,
            "simulator": "libradtran-2.0.6 uvspec (disort)",
            "mol_abs_param": "reptran " + args.mol_abs_param,
            "streams": args.streams,
            "sim_lwdn_sfc": float(f_dn[0]),
            "sim_lwup_sfc": float(f_up[0]),
            "sim_olr_toa": float(f_up[1]),
        })
        if len(inits) <= 12:
            print(f"  {label:<10} LWdn {f_dn[0]:7.2f}  LWup {f_up[0]:7.2f}  "
                  f"OLR {f_up[1]:7.2f} W/m2")

    import pandas as pd
    rpath = results_path(cfg, date, args.hour, args.sky)
    out = pd.DataFrame(rows)
    if rpath.exists():          # keep rows from other simulators (era5_rrtmg_sim.py)
        old = pd.read_csv(rpath)
        other = old[~old["simulator"].astype(str).str.startswith("libradtran")]
        if len(other):
            out = pd.concat([out, other], ignore_index=True)
    out.to_csv(rpath, index=False)
    print(f"wrote {rpath}")
    return 0


# ---------------------------------------------------------------- compare

def planck_band_fraction(T: float, wvl1_um: float, wvl2_um: float) -> float:
    """Fraction of sigma*T^4 emitted within [wvl1, wvl2] um."""
    h, c, kb = 6.62607015e-34, 2.99792458e8, K_B
    wvl = np.logspace(np.log10(0.5e-6), np.log10(3000e-6), 20000)
    planck = (2 * h * c**2 / wvl**5) / (np.expm1(h * c / (wvl * kb * T)))
    band = (wvl >= wvl1_um * 1e-6) & (wvl <= wvl2_um * 1e-6)
    return float(np.trapz(planck[band], wvl[band]) / np.trapz(planck, wvl))


def cmd_compare(args) -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    from era5lib.plotstyle import apply_agu_style, panel_label

    cfg = load_config(args.config)
    date = dt.date(args.year, args.month, args.day)
    manifest = json.loads(manifest_path(cfg, date, args.hour, args.sky).read_text())
    res = pd.read_csv(results_path(cfg, date, args.hour, args.sky))
    sims = res["simulator"].astype(str)
    prof = pd.DataFrame(manifest["profiles"])
    df = prof.merge(res[sims.str.startswith("libradtran")], on="label") \
             .sort_values("sbi_strength", ascending=False)
    rrt = res[sims.str.startswith("rrtmg")]     # era5_rrtmg_sim.py cross-check
    if len(rrt):
        df = df.merge(
            rrt[["label", "sim_lwdn_sfc", "sim_lwup_sfc", "sim_olr_toa"]].rename(
                columns={"sim_lwdn_sfc": "rrtmg_lwdn_sfc",
                         "sim_lwup_sfc": "rrtmg_lwup_sfc",
                         "sim_olr_toa": "rrtmg_olr_toa"}),
            on="label", how="left")
    has_rrtmg = "rrtmg_lwdn_sfc" in df.columns

    w1, w2 = np.array(manifest["wavelength_range_nm"]) / 1000.0
    df["tail_up"] = [(1 - planck_band_fraction(t, w1, w2)) * SIGMA * t**4
                     * manifest["emissivity"] for t in df["skt"]]
    t_eff_dn = (df["era5_lwdn"] / SIGMA) ** 0.25
    df["tail_dn"] = [(1 - planck_band_fraction(t, w1, w2)) * SIGMA * t**4
                     for t in t_eff_dn]

    if len(df) > 12:
        return compare_stats(cfg, args, manifest, df, w1, w2)

    print(f"\nERA5 vs libRadtran {manifest.get('sky', 'clear')}-sky LW fluxes — "
          f"{manifest['snapshot']}"
          f"  (emissivity {manifest['emissivity']}, {w1:.0f}-{w2:.0f} um)")
    print(f"{'label':<10}{'SBI(K)':>7}{'skt(K)':>8} | "
          f"{'LWdn sim':>9}{'ERA5':>7}{'bias':>7}{'tail~':>6} | "
          f"{'LWup sim':>9}{'ERA5':>7}{'bias':>7}{'tail~':>6}"
          + (f" | {'RRTMGdn':>8}{'RRTMGup':>8}" if has_rrtmg else ""))
    for _, r in df.iterrows():
        print(f"{r['label']:<10}{r['sbi_strength']:>7.2f}{r['skt']:>8.2f} | "
              f"{r['sim_lwdn_sfc']:>9.1f}{r['era5_lwdn']:>7.1f}"
              f"{r['sim_lwdn_sfc'] - r['era5_lwdn']:>+7.1f}{r['tail_dn']:>6.1f} | "
              f"{r['sim_lwup_sfc']:>9.1f}{r['era5_lwup']:>7.1f}"
              f"{r['sim_lwup_sfc'] - r['era5_lwup']:>+7.1f}{r['tail_up']:>6.1f}"
              + (f" | {r['rrtmg_lwdn_sfc']:>8.1f}{r['rrtmg_lwup_sfc']:>8.1f}"
                 if has_rrtmg else ""))
    print("tail~ = analytic estimate of flux outside the simulated range "
          "(ERA5 includes it; expected sim low bias)\n")

    # surface blackbody estimate for LW up: emission + reflected downwelling
    eps = manifest["emissivity"]
    df["bb_up"] = eps * SIGMA * df["skt"] ** 4 + (1.0 - eps) * df["era5_lwdn"]

    apply_agu_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.9), layout="constrained")
    x = np.arange(len(df))
    bb_label = ("$\\epsilon\\sigma T_{skin}^{4} + (1{-}\\epsilon)\\,"
                "LW\\!\\downarrow$ (blackbody est.)")
    for ax, kind, letter in ((axes[0], "lwdn", "a"), (axes[1], "lwup", "b")):
        era5 = df["era5_" + kind].values
        sim = df["sim_" + kind + "_sfc"].values
        tail = df["tail_" + kind.replace("lw", "")].values
        for xi, (e, s, t) in enumerate(zip(era5, sim, tail)):
            ax.plot([xi, xi], [e, s], color="#bbbbbb", lw=1, zorder=1)
        if kind == "lwup":
            ax.plot(x, df["bb_up"].values, "D", ms=7, mfc="none",
                    mec="#555555", mew=1.2, label=bb_label, zorder=2)
        ax.plot(x, era5, "o", ms=8, color="#0072B2", label="ERA5", zorder=3)
        ax.plot(x, sim, "s", ms=7, color="#D55E00", label="libRadtran", zorder=3)
        ax.plot(x, sim + tail, "s", ms=7, mfc="none", mec="#D55E00", mew=1.2,
                label="libRadtran + far-IR tail est.", zorder=3)
        if has_rrtmg:
            ax.plot(x, df["rrtmg_" + kind + "_sfc"].values, "^", ms=7,
                    color="#009E73", label="RRTMG-LW (climlab)", zorder=3)
        ax.set_xticks(x)
        if "lwp_g" in df.columns:
            ax.set_xticklabels([f"{l}\n{lw:.0f}/{iw:.0f}" for l, lw, iw in
                                zip(df["label"], df["lwp_g"], df["iwp_g"])],
                               fontsize=8)
            ax.set_xlabel("profile (LWP/IWP, g m$^{-2}$)")
        else:
            ax.set_xticklabels([f"{l}\n{s:.1f} K" for l, s in
                                zip(df["label"], df["sbi_strength"])], fontsize=8)
            ax.set_xlabel("profile (SBI strength)")
        ax.set_ylabel("LW$\\downarrow$ surface (W m$^{-2}$)" if kind == "lwdn"
                      else "LW$\\uparrow$ surface (W m$^{-2}$)")
        ax.grid(color="#e3e3e3", lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        panel_label(ax, letter, x=-0.14, y=1.06)
    # one shared legend below the panels so it never overlaps data points
    handles, labels = axes[1].get_legend_handles_labels()
    want = ["ERA5", "libRadtran", "libRadtran + far-IR tail est.",
            "RRTMG-LW (climlab)", bb_label]
    order = [labels.index(w) for w in want if w in labels]
    fig.legend([handles[i] for i in order], [labels[i] for i in order],
               loc="lower center", bbox_to_anchor=(0.5, -0.1),
               ncol=3 if has_rrtmg else 2, frameon=False, fontsize=9)
    tail_eq = (
        "far-IR tail $= \\left[1 - F_{%.2f\\mathrm{-}%.0f\\,\\mu m}(T)\\right]"
        "\\,\\epsilon\\,\\sigma T^{4}$,   "
        "$F_{band}(T) = \\int_{band} B_{\\lambda}(T)\\,d\\lambda \\; / \\;"
        "\\sigma T^{4}$\n"
        "$T = T_{skin},\\ \\epsilon{=}%.2f$ for LW$\\uparrow$;   "
        "$T = (LW\\!\\downarrow_{ERA5}/\\sigma)^{1/4},\\ \\epsilon{=}1$ "
        "for LW$\\downarrow$" % (w1, w2, eps))
    sky = manifest.get("sky", "clear")
    sims_txt = "libRadtran & RRTMG-LW" if has_rrtmg else "libRadtran"
    fig.suptitle(f"{sky.capitalize()}-sky broadband LW: {sims_txt} "
                 f"(ERA5 profiles) vs ERA5 — {manifest['snapshot']}\n"
                 f"$\\epsilon$={manifest['emissivity']}, sur_temperature = ERA5 skt; "
                 f"sim {w1:.2f}–{w2:.0f} µm vs ERA5 3.08–1000 µm; "
                 "ERA5 fluxes are 1-h means\n" + tail_eq, y=1.26)

    outdir = figures_dir(cfg)
    outdir.mkdir(parents=True, exist_ok=True)
    fpath = outdir / f"lw_sim_{date:%Y%m%d}T{args.hour:02d}Z{sky_tag(sky)}.png"
    fig.savefig(fpath, bbox_inches="tight")
    print(f"wrote {fpath}")
    return 0


def compare_stats(cfg, args, manifest, df, w1, w2):
    """Statistical comparison for large pixel samples (n > 12).

    With era5_rrtmg_sim.py results present the figure becomes a 2x2 three-way
    comparison; panel (d) pits RRTMG (full 3.08-1000 um band) against
    libRadtran + the analytic far-IR tail, a direct check of that estimate.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    from era5lib.plotstyle import apply_agu_style, panel_label

    sky = manifest.get("sky", "clear")
    date = dt.date(args.year, args.month, args.day)
    has_rrtmg = ("rrtmg_lwdn_sfc" in df.columns
                 and df["rrtmg_lwdn_sfc"].notna().all())
    print(f"\nERA5 vs simulated {sky}-sky LW statistics — {manifest['snapshot']}"
          f"  (n = {len(df)}, emissivity {manifest['emissivity']}, "
          f"libRadtran {w1:.2f}-{w2:.0f} um"
          + (", RRTMG-LW 3.08-1000 um" if has_rrtmg else "") + ")")
    stats = {}
    for kind in ("lwdn", "lwup"):
        e = df["era5_" + kind].values
        s = df["sim_" + kind + "_sfc"].values
        t = df["tail_" + kind.replace("lw", "")].values
        stats[kind] = (e, s, t)
        print(f"  {kind} libRadtran: r={np.corrcoef(e, s)[0, 1]:+.3f}  "
              f"bias={(s - e).mean():+.2f}  "
              f"rmse={np.sqrt(((s - e) ** 2).mean()):.2f}  "
              f"bias+tail={(s + t - e).mean():+.2f} W/m2")
        if has_rrtmg:
            g = df["rrtmg_" + kind + "_sfc"].values
            print(f"  {kind} RRTMG-LW:   r={np.corrcoef(e, g)[0, 1]:+.3f}  "
                  f"bias={(g - e).mean():+.2f}  "
                  f"rmse={np.sqrt(((g - e) ** 2).mean()):.2f}  "
                  f"RRTMG-(lib+tail)={(g - s - t).mean():+.2f} W/m2")

    apply_agu_style()
    if has_rrtmg:
        fig, ax4 = plt.subplots(2, 2, figsize=(9.8, 9.4), layout="constrained")
        axes, lab_x = ax4.ravel(), -0.14
    else:
        fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.6), layout="constrained")
        lab_x = -0.18
    for i, kind in enumerate(("lwdn", "lwup")):
        ax = axes[i]
        e, s, t = stats[kind]
        vals = [e, s]
        if has_rrtmg:
            g = df["rrtmg_" + kind + "_sfc"].values
            vals.append(g)
        lim = (min(v.min() for v in vals) - 3, max(v.max() for v in vals) + 3)
        ax.plot(lim, lim, color="#bbbbbb", lw=1)
        if has_rrtmg:
            ax.scatter(e, s, s=9, color="#D55E00", alpha=0.3,
                       edgecolors="none", label="libRadtran")
            ax.scatter(e, g, s=9, color="#009E73", alpha=0.3, marker="^",
                       edgecolors="none", label="RRTMG-LW")
            txt = (f"n = {len(e)}\n"
                   f"libRadtran: r = {np.corrcoef(e, s)[0, 1]:+.3f}, "
                   f"bias = {(s - e).mean():+.2f} "
                   f"(+tail: {(s + t - e).mean():+.2f})\n"
                   f"RRTMG-LW: r = {np.corrcoef(e, g)[0, 1]:+.3f}, "
                   f"bias = {(g - e).mean():+.2f} W m$^{{-2}}$")
            ax.legend(frameon=False, fontsize=8, loc="lower right")
        else:
            ax.scatter(e, s, s=9, color="#0072B2", alpha=0.3, edgecolors="none")
            txt = (f"n = {len(e)}\nr = {np.corrcoef(e, s)[0, 1]:+.3f}\n"
                   f"bias = {(s - e).mean():+.2f}\n"
                   f"bias+tail = {(s + t - e).mean():+.2f} W m$^{{-2}}$")
        ax.text(0.04, 0.96, txt, transform=ax.transAxes, va="top",
                fontsize=8 if has_rrtmg else 9,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#cccccc"))
        arrow = "\\downarrow" if kind == "lwdn" else "\\uparrow"
        ax.set_xlabel(f"ERA5 LW${arrow}$ (W m$^{{-2}}$)")
        ax.set_ylabel(("simulated" if has_rrtmg else "libRadtran")
                      + f" LW${arrow}$ (W m$^{{-2}}$)")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect("equal")
        ax.grid(color="#e3e3e3", lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        panel_label(ax, "ab"[i], x=lab_x, y=1.06)

    ax = axes[2]
    if sky == "cloudy":
        xvals = df["lwp_g"].values + df["iwp_g"].values
        ax.set_xlabel("LWP + IWP (g m$^{-2}$)")
        ax.set_xscale("log")
    else:
        xvals = df["sbi_strength"].values
        ax.set_xlabel("SBI strength (K)")
    ax.axhline(0, color="#bbbbbb", lw=1)
    for kind, color, lab in (("lwdn", "#0072B2", "LW$\\downarrow$"),
                             ("lwup", "#D55E00", "LW$\\uparrow$")):
        e, s, _ = stats[kind]
        ax.scatter(xvals, s - e, s=9, color=color, alpha=0.3, edgecolors="none",
                   label=lab + (" libRadtran" if has_rrtmg else ""))
        if has_rrtmg:
            g = df["rrtmg_" + kind + "_sfc"].values
            ax.scatter(xvals, g - e, s=9, color=color, alpha=0.3, marker="^",
                       edgecolors="none", label=lab + " RRTMG")
    ax.set_ylabel("bias sim − ERA5 (W m$^{-2}$)")
    ax.legend(frameon=False, fontsize=8 if has_rrtmg else 9,
              ncol=2 if has_rrtmg else 1)
    ax.grid(color="#e3e3e3", lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    panel_label(ax, "c", x=lab_x, y=1.06)

    if has_rrtmg:
        # head-to-head on identical inputs: isolates RT method + spectral
        # coverage, and directly validates the analytic far-IR tail estimate
        ax = axes[3]
        allv = []
        for kind, color, lab in (("lwdn", "#0072B2", "LW$\\downarrow$"),
                                 ("lwup", "#D55E00", "LW$\\uparrow$")):
            e, s, t = stats[kind]
            g = df["rrtmg_" + kind + "_sfc"].values
            diff = g - s - t
            ax.scatter(s + t, g, s=9, color=color, alpha=0.3, edgecolors="none",
                       label=(f"{lab}  $\\Delta$ = {diff.mean():+.2f}"
                              f" $\\pm$ {diff.std():.2f}"))
            allv += [s + t, g]
        lim = (min(v.min() for v in allv) - 3, max(v.max() for v in allv) + 3)
        ax.plot(lim, lim, color="#bbbbbb", lw=1, zorder=0)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect("equal")
        ax.set_xlabel("libRadtran + far-IR tail (W m$^{-2}$)")
        ax.set_ylabel("RRTMG-LW (W m$^{-2}$)")
        ax.legend(frameon=False, fontsize=8, loc="lower right",
                  title="$\\Delta$ = RRTMG $-$ (lib+tail), W m$^{-2}$",
                  title_fontsize=8)
        ax.grid(color="#e3e3e3", lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        panel_label(ax, "d", x=lab_x, y=1.06)

    sims_txt = "libRadtran & RRTMG-LW" if has_rrtmg else "libRadtran"
    fig.suptitle(f"{sky.capitalize()}-sky broadband LW statistics: {sims_txt} "
                 f"(ERA5 profiles) vs ERA5 — {manifest['snapshot']}, "
                 f"n = {len(df)} pixels", y=1.03 if has_rrtmg else 1.06)
    outdir = figures_dir(cfg)
    outdir.mkdir(parents=True, exist_ok=True)
    fpath = outdir / f"lw_sim_stats_{date:%Y%m%d}T{args.hour:02d}Z{sky_tag(sky)}.png"
    fig.savefig(fpath, bbox_inches="tight")
    print(f"wrote {fpath}")
    return 0


# ---------------------------------------------------------------- CLI

def main(argv: "list[str] | None" = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--year", type=int, required=True)
    common.add_argument("--month", type=int, required=True)
    common.add_argument("--day", type=int, required=True)
    common.add_argument("--hour", type=int, required=True)
    common.add_argument("--sky", default="clear", choices=["clear", "cloudy"],
                        help="clear: cloud-free pixels; cloudy: near-overcast "
                             "pixels with ERA5 clwc/ciwc as wc/ic cloud files")
    common.add_argument("--config", default=None)

    p_prep = sub.add_parser("prep", parents=[common],
                            help="select clear pixels, build profile files")
    p_prep.add_argument("--n", type=int, default=5,
                        help="pixels to simulate; > 5 switches to a seeded "
                             "random sample of the eligible population")
    p_prep.add_argument("--seed", type=int, default=0)
    p_prep.add_argument("--fallback-constants", action="store_true",
                        help=f"skip CAMS; CO2 {CO2_CONST_PPM} ppm, "
                             f"CH4 {CH4_CONST_PPM} ppm")
    p_prep.set_defaults(func=cmd_prep)

    p_run = sub.add_parser("run", parents=[common], help="run uvspec thermal jobs")
    p_run.add_argument("--streams", type=int, default=None,
                       help="disort streams (default: 4 on macOS, 8 on Linux)")
    p_run.add_argument("--mol-abs-param", default="coarse",
                       choices=["coarse", "medium", "fine"])
    p_run.add_argument("--workers", type=int, default=None)
    p_run.add_argument("--overwrite", action="store_true")
    p_run.set_defaults(func=cmd_run)

    p_cmp = sub.add_parser("compare", parents=[common],
                           help="table + figure vs ERA5 fluxes")
    p_cmp.set_defaults(func=cmd_compare)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
