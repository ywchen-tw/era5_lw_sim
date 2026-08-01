#!/usr/bin/env python
"""MOSAiC case study: one clear-sky and one cloudy ERA5 pixel, end to end.

Picks the two ERA5 pixels matched to MOSAiC radiosondes (Jozef et al. 2023)
that best represent clear and overcast conditions in the study month, then
walks them through the whole pipeline in one figure each:

  (a) ERA5 temperature profile vs the observed MOSAiC structure
  (b) how the inversion strengths are computed (ERA5 SBI scan vs obs criterion)
  (c) the radiative-transfer input profile (ERA5 levels + afglsw splice, clouds)
  (d) surface LW fluxes: libRadtran & RRTMG-LW vs ERA5 and the MOSAiC obs

Subcommands / envs:
  prep    (era5 env)     select cases, build atm/cloud files + case manifest
  run     (er3t_env)     uvspec thermal jobs for the two cases
  figure  (era5 env)     RRTMG-LW inline + the two case-study figures

Examples:
    conda activate era5     && python src/era5_case_study.py prep
    conda activate er3t_env && python src/era5_case_study.py run
    conda activate era5     && python src/era5_case_study.py figure
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import platform
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from era5lib.config import REPO_ROOT, figures_dir, inversion_path, load_config, plev_path, sfc_path
from era5lib.inversion import column_heights
from era5lib.io_era5 import open_era5
from era5_lrt_sim import (EMISSIVITY, REFF_ICE_UM, REFF_LIQ_UM, SIGMA,
                         WVL_RANGE_NM, build_profile_files, load_cams_profiles,
                         write_cloud_file)

TRAPZ = getattr(np, "trapezoid", getattr(np, "trapz", None))
G0 = 9.80665
K_B = 1.380649e-23
CASE_DIR = REPO_ROOT / "derived" / "case_study"
MANIFEST = CASE_DIR / "case_manifest.json"
RESULTS = CASE_DIR / "case_results.json"
# MOSAiC level-2 radiosonde profiles (Maturilli et al. 2021, CC-BY-4.0)
SOUNDING_TAB = REPO_ROOT / "data" / "mosaic" / "soundings" / "PS122_2_radiosonde_202001.tab"
SOUNDING_URL = "https://doi.pangaea.de/10.1594/PANGAEA.928659?format=textfile"
_SOUNDING_CACHE = {}


def load_sounding_profile(sounding_time, z_max_m=5000.0):
    """T(z) of the MOSAiC level-2 radiosonde nearest the given launch time."""
    import pandas as pd
    if not SOUNDING_TAB.exists():
        SOUNDING_TAB.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading MOSAiC radiosonde level-2 data (~60 MB) -> "
              f"{SOUNDING_TAB.name} ...")
        try:
            import urllib.request
            urllib.request.urlretrieve(SOUNDING_URL, SOUNDING_TAB)
        except Exception as e:
            print(f"  download failed ({e}); figure falls back to the derived "
                  "inversion layers")
            return None
    if "df" not in _SOUNDING_CACHE:
        with open(SOUNDING_TAB) as f:
            skip = next(i for i, line in enumerate(f) if line.startswith("*/")) + 1
        df = pd.read_csv(SOUNDING_TAB, sep="\t", skiprows=skip,
                         usecols=["Date/Time", "Altitude [m]", "TTT [°C]"])
        df["Date/Time"] = pd.to_datetime(df["Date/Time"], format="ISO8601")
        _SOUNDING_CACHE["df"] = df
    df = _SOUNDING_CACHE["df"]
    st = pd.Timestamp(sounding_time)
    m = (df["Date/Time"] - st).abs() <= pd.Timedelta("90min")
    sel = df[m & (df["Altitude [m]"] <= z_max_m)].sort_values("Altitude [m]")
    if sel.empty:
        return None
    return {"z_m": [float(v) for v in sel["Altitude [m]"].values[::4]],
            "t_k": [float(v) + 273.15 for v in sel["TTT [°C]"].values[::4]]}


def era5_match_time(sounding_time):
    """Same 6-h rounding as era5_mosaic_compare: nearest of 00/06/12/18."""
    import pandas as pd
    return (pd.Timestamp(sounding_time) + pd.Timedelta(hours=3)).floor("6h")


def planck_band_fraction(T, wvl1_um, wvl2_um):
    h, c = 6.62607015e-34, 2.99792458e8
    wvl = np.logspace(np.log10(0.5e-6), np.log10(3000e-6), 20000)
    planck = (2 * h * c**2 / wvl**5) / np.expm1(h * c / (wvl * K_B * T))
    band = (wvl >= wvl1_um * 1e-6) & (wvl <= wvl2_um * 1e-6)
    return float(TRAPZ(planck[band], wvl[band]) / TRAPZ(planck, wvl))


# ---------------------------------------------------------------- prep

def cmd_prep(args) -> int:
    import pandas as pd
    import xarray as xr

    cfg = load_config(args.config)
    pairs = xr.open_dataset(
        REPO_ROOT / "derived" / f"{args.year}" / f"{args.month:02d}" /
        f"era5_mosaic_pairs_{args.year}{args.month:02d}.nc")
    mos = xr.open_dataset(REPO_ROOT / "data" / "mosaic" / "MOSAiC_Atm_Properties.nc")

    cache = {}

    def day_data(date):
        if date not in cache:
            cache[date] = (open_era5(plev_path(cfg, date)),
                           open_era5(sfc_path(cfg, date)),
                           open_era5(inversion_path(cfg, date)))
        return cache[date]

    # score every matched sounding for clear/cloudy suitability
    rows = []
    for k in range(pairs.time.size):
        st = pd.Timestamp(pairs.time.values[k])
        et = era5_match_time(st)
        if et.month != args.month or et.year != args.year:
            continue
        date = et.date()
        if not plev_path(cfg, date).exists():
            continue
        plev, sfc, invd = day_data(date)
        p5 = plev.sel(valid_time=et).transpose("pressure_level", "latitude", "longitude")
        iy = int(np.abs(p5["latitude"].values - float(pairs["lat"][k])).argmin())
        ix = int(np.abs(p5["longitude"].values - float(pairs["lon"][k])).argmin())
        p_pa_asc = p5["pressure_level"].values[::-1].astype(float) * 100.0
        cc_max = float(p5["cc"].values[:, iy, ix].max())
        lwp = float(TRAPZ(p5["clwc"].values[::-1, iy, ix], x=p_pa_asc) / G0 * 1e3)
        iwp = float(TRAPZ(p5["ciwc"].values[::-1, iy, ix], x=p_pa_asc) / G0 * 1e3)
        ob = mos.sel(time=st)
        rows.append({"k": k, "sounding_time": st, "era5_time": et,
                     "iy": iy, "ix": ix, "cc_max": cc_max, "lwp": lwp,
                     "iwp": iwp, "obs_cc": float(ob["cc"]),
                     "obs_cbh": float(ob["cbh"]),
                     "era5_sbi": float(pairs["era5_strength"][k]),
                     "obs_sbi": float(pairs["obs_strength"][k])})
    df = pd.DataFrame(rows)

    clear = df[(df.cc_max <= 0.01) & (df.lwp + df.iwp <= 1.0) & (df.obs_cc <= 0.2)]
    if clear.empty:
        clear = df[(df.cc_max <= 0.05) & (df.lwp + df.iwp <= 1.0) & (df.obs_cc <= 0.3)]
        print("note: relaxed clear thresholds (cc_max<=0.05, obs_cc<=0.3)")
    cloudy = df[(df.cc_max >= 0.99) & (df.lwp + df.iwp > 5.0) & (df.obs_cc >= 0.8)]
    if clear.empty or cloudy.empty:
        sys.exit(f"no suitable case: clear n={len(clear)}, cloudy n={len(cloudy)}")
    pick_clear = clear.loc[clear.era5_sbi.idxmax()]     # strongest SBI: best demo
    pick_cloudy = cloudy.loc[(cloudy.lwp + cloudy.iwp).idxmax()]
    print(f"clear candidates {len(clear)}, cloudy candidates {len(cloudy)}")

    CASE_DIR.mkdir(parents=True, exist_ok=True)
    cams = load_cams_profiles(args.year, args.month)
    cases = []
    for label, pick in (("clear", pick_clear), ("cloudy", pick_cloudy)):
        k, et, st = int(pick.k), pick.era5_time, pick.sounding_time
        date = et.date()
        plev, sfc, invd = day_data(date)
        p5 = plev.sel(valid_time=et).transpose("pressure_level", "latitude", "longitude")
        s5, i5 = sfc.sel(valid_time=et), invd.sel(valid_time=et)
        iy, ix = int(pick.iy), int(pick.ix)

        p = p5["pressure_level"].values.astype(float)
        sp_hpa = float(s5["sp"].values[iy, ix]) / 100.0
        t2m = float(s5["t2m"].values[iy, ix])
        skt = float(s5["skt"].values[iy, ix])
        above = p <= sp_hpa
        cols = {v: p5[v].values[above, iy, ix].astype(float)
                for v in ("t", "q", "o3", "clwc", "ciwc")}
        z_km = column_heights(cols["t"], p[above], t2m, sp_hpa, cols["q"]) / 1e3

        atm_file, ch4_file, keep, z_keep = build_profile_files(
            CASE_DIR, label, p[above], cols["t"], cols["q"], cols["o3"],
            t2m, sp_hpa, cams)
        wc_file = ic_file = None
        if label == "cloudy":
            wc_file = write_cloud_file(CASE_DIR / f"wc_{label}.dat", z_keep,
                                       cols["clwc"][keep], p[above][keep],
                                       cols["t"][keep], REFF_LIQ_UM)
            ic_file = write_cloud_file(CASE_DIR / f"ic_{label}.dat", z_keep,
                                       cols["ciwc"][keep], p[above][keep],
                                       cols["t"][keep], REFF_ICE_UM)

        ob = mos.sel(time=st)
        inv_obs = []
        for j in range(int(mos.sizes["inversion_number"])):
            alt = float(ob["inv_alt"][j])
            if np.isfinite(alt):
                inv_obs.append({"alt_m": alt,
                                "dz_m": float(ob["inv_dz"][j]),
                                "t_k": float(ob["inv_t"][j]) + 273.15,
                                "dt_k": float(ob["inv_dt"][j])})
        strd = float(s5["strd"].values[iy, ix]) / 3600.0
        str_net = float(s5["str"].values[iy, ix]) / 3600.0
        cases.append({
            "label": label,
            "sounding_time": st.isoformat(), "era5_time": et.isoformat(),
            "pixel_lat": float(p5["latitude"].values[iy]),
            "pixel_lon": float(p5["longitude"].values[ix]),
            "obs_lat": float(pairs["lat"][k]), "obs_lon": float(pairs["lon"][k]),
            "match_km": float(pairs["match_km"][k]),
            "match_dt_h": float(pairs["match_dt_h"][k]),
            "skt": skt, "t2m": t2m, "sp_hpa": sp_hpa,
            "era5_lwdn": strd, "era5_lwup": strd - str_net,
            "cc_max": float(pick.cc_max), "lwp_g": float(pick.lwp),
            "iwp_g": float(pick.iwp),
            "profile": {kk: list(vv) for kk, vv in
                        [("p_hpa", p[above]), ("t", cols["t"]), ("q", cols["q"]),
                         ("clwc", cols["clwc"]), ("ciwc", cols["ciwc"]),
                         ("z_km", z_km)]},
            "era5_sbi": {"found": bool(i5["sbi_found"].values[iy, ix]),
                         "strength": float(i5["sbi_strength"].values[iy, ix]),
                         "top_p": float(i5["sbi_top_p"].values[iy, ix]),
                         "depth_z": float(i5["sbi_depth_z"].values[iy, ix]),
                         "dt_850_2m": float(i5["dt_850_2m"].values[iy, ix]),
                         "dt_925_1000": float(i5["dt_925_1000"].values[iy, ix])},
            "obs": {"t_2m": float(ob["t_2m"]) + 273.15,
                    "t_10m": float(ob["t_10m"]) + 273.15,
                    "t_h": float(ob["t_h"]) + 273.15, "h_m": float(ob["h"]),
                    "cc": float(ob["cc"]), "cbh_m": float(ob["cbh"]),
                    "lwdn": float(ob["lwdn"]), "lwup": float(ob["lwup"]),
                    "inversions": inv_obs,
                    "sbi_found": bool(pairs["obs_sbi_found"][k]),
                    "sbi_strength": float(pairs["obs_strength"][k]),
                    "sbi_base_m": float(pairs["obs_base"][k]),
                    "sbi_depth_m": float(pairs["obs_depth"][k])},
            "obs_profile": load_sounding_profile(st),
            "atm_file": atm_file, "ch4_file": ch4_file,
            "wc_file": wc_file, "ic_file": ic_file,
        })
        print(f"{label:>6}: sounding {st:%Y-%m-%d %H:%M} -> ERA5 {et:%Y-%m-%d %H}Z"
              f"  ({cases[-1]['pixel_lat']:.2f}N, {cases[-1]['pixel_lon']:.2f}E)"
              f"  cc_max {pick.cc_max:.2f}  LWP+IWP {pick.lwp + pick.iwp:.1f} g/m2"
              f"  ERA5 SBI {pick.era5_sbi:.1f} K  obs SBI {pick.obs_sbi:.1f} K")

    manifest = {"emissivity": EMISSIVITY,
                "n_plev_total": int(p.size),
                "wavelength_range_nm": list(WVL_RANGE_NM),
                "reff_um": {"liquid": REFF_LIQ_UM, "ice": REFF_ICE_UM},
                "gas_source": cams["source"] if cams else "constants",
                "cases": cases}
    MANIFEST.write_text(json.dumps(manifest, indent=2))
    print(f"wrote {MANIFEST}")
    return 0


# ---------------------------------------------------------------- run (er3t_env)

def cmd_run(args) -> int:
    import copy
    import er3t

    manifest = json.loads(MANIFEST.read_text())
    tmp = CASE_DIR / "tmp"
    tmp.mkdir(exist_ok=True)
    streams = 4 if platform.system() == "Darwin" else 8

    lrt_cfg_base = er3t.rtm.lrt.get_lrt_cfg()
    lrt_cfg_base["mol_abs_param"] = "reptran coarse"
    lrt_cfg_base["number_of_streams"] = streams
    mute = ["albedo", "wavelength", "spline", "source solar", "slit_function_file"]

    inits = []
    for case in manifest["cases"]:
        label = case["label"]
        et = dt.datetime.fromisoformat(case["era5_time"])
        lrt_cfg = copy.deepcopy(lrt_cfg_base)
        lrt_cfg["atmosphere_file"] = case["atm_file"]
        extra = {
            "source": "thermal",
            "albedo_add": f"{1.0 - manifest['emissivity']:.4f}",
            "sur_temperature": f"{case['skt']:.2f}",
            "wavelength_add": "{:.0f} {:.0f}".format(*manifest["wavelength_range_nm"]),
            "output_process": "integrate",
            "mol_file": "CH4 " + case["ch4_file"],
        }
        if case.get("wc_file"):
            extra["wc_file 1D"] = case["wc_file"]
            extra["wc_properties"] = "hu interpolate"
        if case.get("ic_file"):
            extra["ic_file 1D"] = case["ic_file"]
            extra["ic_properties"] = "fu interpolate"
        init = er3t.rtm.lrt.lrt_init_mono_flx(
            input_file=str(tmp / f"input_{label}.txt"),
            output_file=str(tmp / f"output_{label}.txt"),
            date=et, solar_zenith_angle=80.0, Nx=1,
            output_altitude=[0, "toa"], input_dict_extra=extra,
            mute_list=mute, lrt_cfg=lrt_cfg, cld_cfg=None, aer_cfg=None)
        inits.append((label, init))

    er3t.rtm.lrt.lrt_run_mp([i for _, i in inits], Ncpu=2)
    results = {}
    for label, init in inits:
        data = er3t.rtm.lrt.lrt_read_uvspec_flx([init])
        f_dn = np.squeeze(data.f_down) * 1000.0     # er3t /1000 undo (thermal is W)
        f_up = np.squeeze(data.f_up) * 1000.0
        results[label] = {"lib_lwdn": float(f_dn[0]), "lib_lwup": float(f_up[0]),
                          "lib_olr": float(f_up[1])}
        print(f"  {label:<7} LWdn {f_dn[0]:7.2f}  LWup {f_up[0]:7.2f} W/m2")
    RESULTS.write_text(json.dumps(results, indent=2))
    print(f"wrote {RESULTS}")
    return 0


# ---------------------------------------------------------------- figure (era5 env)

def draw_case(cfg, manifest, case, res, rrtmg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from era5lib.plotstyle import apply_agu_style, panel_label

    C_ERA5, C_OBS, C_LIB, C_RRT = "#0072B2", "#000000", "#D55E00", "#009E73"
    prof = {k: np.array(v) for k, v in case["profile"].items()}
    obs, sbi = case["obs"], case["era5_sbi"]
    label = case["label"]
    w1, w2 = np.array(manifest["wavelength_range_nm"]) / 1000.0
    eps = manifest["emissivity"]

    apply_agu_style()
    fig, axes = plt.subplots(1, 4, figsize=(14.2, 6.2), layout="constrained")

    # (a) ERA5 profile vs MOSAiC structure -----------------------------------
    ax = axes[0]
    z_m = np.concatenate([[2.0], prof["z_km"] * 1e3])
    t_all = np.concatenate([[case["t2m"]], prof["t"]])
    ax.plot(t_all, z_m, "-o", ms=3.5, color=C_ERA5, lw=1.2,
            label="ERA5 (37 plev + 2 m)")
    ax.plot(case["skt"], 0.0, "s", ms=6, color=C_ERA5, mfc="white",
            label="ERA5 skin T")
    op = case.get("obs_profile")
    if op:
        ax.plot(op["t_k"], op["z_m"], "-", color=C_OBS, lw=1.0,
                label="MOSAiC radiosonde (level 2)")
    else:       # fallback when the level-2 sounding file is unavailable
        for j, iv in enumerate(obs["inversions"]):
            ax.plot([iv["t_k"], iv["t_k"] + iv["dt_k"]],
                    [iv["alt_m"], iv["alt_m"] + iv["dz_m"]], "-", lw=3,
                    color=C_OBS, alpha=0.45,
                    label="MOSAiC inversion layers" if j == 0 else None)
    ax.plot([obs["t_2m"], obs["t_10m"]], [2, 10], "*", ms=8, color=C_OBS,
            ls="none", label="MOSAiC 2 m / 10 m")
    if np.isfinite(obs["h_m"]) and np.isfinite(obs["t_h"]):
        ax.plot(obs["t_h"], obs["h_m"], "P", ms=6, color=C_OBS, mfc="white",
                label="MOSAiC BL top")
    if label == "cloudy":
        if np.isfinite(obs["cbh_m"]) and obs["cbh_m"] < 2400:
            ax.axhline(obs["cbh_m"], color="#CC79A7", lw=1.1, ls="-.",
                       label=f"obs cloud base ({obs['cbh_m']:.0f} m)")
        elif np.isfinite(obs["cbh_m"]):
            ax.text(0.97, 0.62, f"obs cloud base: {obs['cbh_m'] / 1e3:.1f} km\n"
                    "(above panel)", transform=ax.transAxes, ha="right",
                    fontsize=6.6, color="#666666")
    ax.set_ylim(-40, 2500)
    tvis = [t_all[z_m < 2600], [case["skt"], obs["t_2m"], obs["t_10m"]]]
    if op:
        opz, opt = np.array(op["z_m"]), np.array(op["t_k"])
        tvis.append(opt[opz < 2600])
    else:
        tvis += [[iv["t_k"], iv["t_k"] + iv["dt_k"]] for iv in obs["inversions"]
                 if iv["alt_m"] < 2600]
    tvis = np.concatenate([np.atleast_1d(v) for v in tvis])
    tvis = tvis[np.isfinite(tvis)]
    ax.set_xlim(tvis.min() - 4, tvis.max() + 4)
    if label == "cloudy":
        # ERA5 cloud water as a scaled fill along the left edge (peak in legend)
        x0, x1 = ax.get_xlim()
        files = []
        for f, cc, lab in ((case.get("wc_file"), "#56B4E9", "LWC"),
                           (case.get("ic_file"), "#CC79A7", "IWC")):
            if f:
                cw = np.atleast_2d(np.loadtxt(f, comments="#"))
                if (cw[:, 1] > 0).any():
                    files.append((cw, cc, lab))
        wmax = max(float(cw[:, 1].max()) for cw, _, _ in files) if files else 0.0
        for cw, cc, lab in files:
            ax.fill_betweenx(cw[:, 0] * 1e3, x0,
                             x0 + cw[:, 1] / wmax * 0.28 * (x1 - x0),
                             color=cc, alpha=0.3, lw=0.8, zorder=0,
                             label=f"ERA5 {lab} (peak "
                                   f"{cw[:, 1].max() * 1e3:.0f} mg m$^{{-3}}$)")
        ax.set_xlim(x0, x1)
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel("height (m)")
    ax.legend(frameon=False, fontsize=6.6, ncol=2, loc="lower left",
              bbox_to_anchor=(0.0, 1.0), borderaxespad=0.2)
    panel_label(ax, "a", x=-0.22, y=1.15)
    ax.set_title("ERA5 profile vs MOSAiC", fontsize=9, pad=44)

    # (b) inversion-strength construction ------------------------------------
    ax = axes[1]
    ztop = 900.0
    if sbi["found"]:
        ztop = max(ztop, 1.4 * sbi["depth_z"])
    if obs["sbi_found"] and np.isfinite(obs["sbi_depth_m"]):
        ztop = max(ztop, 1.4 * (obs["sbi_base_m"] + obs["sbi_depth_m"]))
    if op:
        ax.plot(opt, opz, "-", color=C_OBS, lw=1.0, alpha=0.55)
    else:
        for iv in obs["inversions"]:
            ax.plot([iv["t_k"], iv["t_k"] + iv["dt_k"]],
                    [iv["alt_m"], iv["alt_m"] + iv["dz_m"]], "-", lw=3,
                    color=C_OBS, alpha=0.45)
    ax.plot(t_all, z_m, "-o", ms=3.5, color=C_ERA5, lw=1.2)
    lines = []
    if sbi["found"]:
        i_top = int(np.argmin(np.abs(prof["p_hpa"] - sbi["top_p"])))
        t_base, t_top = case["t2m"], prof["t"][i_top]
        z_top = prof["z_km"][i_top] * 1e3
        ax.plot([t_base, t_top], [2, z_top], "o", ms=7, mfc="none",
                mec=C_ERA5, mew=1.6)
        ax.annotate("", xy=(t_top, z_top), xytext=(t_base, z_top),
                    arrowprops=dict(arrowstyle="<->", color=C_ERA5, lw=1.2))
        ax.axvline(t_base, color=C_ERA5, lw=0.7, ls=":")
        ax.text(0.5 * (t_base + t_top), z_top + 0.02 * ztop,
                f"ERA5 $\\Delta T$ = {sbi['strength']:.1f} K", color=C_ERA5,
                fontsize=7, ha="center", va="bottom", fontweight="bold")
        lines.append(f"ERA5 SBI: $\\Delta T = T_{{top}} - T_{{2m}}$ = "
                     f"{sbi['strength']:.1f} K (crit. $\\geq$ 0.5 K),")
        lines.append(f"   top {sbi['top_p']:.0f} hPa, "
                     f"depth {sbi['depth_z']:.0f} m")
    else:
        lines.append("ERA5 SBI: none detected")
    if obs["sbi_found"] and op and np.isfinite(obs["sbi_depth_m"]):
        zb = obs["sbi_base_m"]
        zt = zb + obs["sbi_depth_m"]
        seg = (opz >= zb) & (opz <= zt)
        ax.plot(opt[seg], opz[seg], "-", color=C_OBS, lw=2.6)
        tb_o, tt_o = np.interp([zb, zt], opz, opt)
        ax.annotate("", xy=(tt_o, zt), xytext=(tb_o, zt),
                    arrowprops=dict(arrowstyle="<->", color=C_OBS, lw=1.2))
        ax.axvline(tb_o, color=C_OBS, lw=0.7, ls=":")
        ax.text(0.5 * (tb_o + tt_o), zt + 0.02 * ztop,
                f"obs $\\Delta T$ = {obs['sbi_strength']:.1f} K", color=C_OBS,
                fontsize=7, ha="center", va="bottom", fontweight="bold")
    if obs["sbi_found"]:
        lines.append(f"obs SBI: {obs['sbi_strength']:.1f} K over "
                     f"{obs['sbi_depth_m']:.0f} m")
        lines.append("   (crit. 0.65 $^{\\circ}$C/100 m over $\\geq$ 25 m)")
    else:
        lines.append("obs SBI: none (base > 100 m or below criterion)")
    lines.append(f"fixed metrics: T$_{{850}}$$-$T$_{{2m}}$ = "
                 f"{sbi['dt_850_2m']:+.1f} K,")
    lines.append(f"   T$_{{925}}$$-$T$_{{1000}}$ = {sbi['dt_925_1000']:+.1f} K")
    ax.text(0.0, -0.19, "\n".join(lines), transform=ax.transAxes, va="top",
            fontsize=7.2,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#cccccc"))
    ax.set_ylim(-25, ztop)
    tvis = [t_all[z_m < ztop * 1.1]]
    if op:
        tvis.append(opt[opz < ztop * 1.1])
    else:
        tvis.append(np.array([t for iv in obs["inversions"]
                              if iv["alt_m"] < ztop
                              for t in (iv["t_k"], iv["t_k"] + iv["dt_k"])]))
    tvis = np.concatenate(tvis)
    tvis = tvis[np.isfinite(tvis)]
    ax.set_xlim(tvis.min() - 3, tvis.max() + 3)
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel("height (m)")
    panel_label(ax, "b", x=-0.22, y=1.15)
    ax.set_title("inversion strength (zoom of a)", fontsize=9, pad=44)
    # mark the zoom window in panel (a)
    from matplotlib.patches import Rectangle
    xb = ax.get_xlim()
    axes[0].add_patch(Rectangle((xb[0], -25), xb[1] - xb[0], ztop + 25,
                                fill=False, ec="#888888", ls="--", lw=0.9))
    axes[0].text(xb[1], ztop + 45, "(b)", fontsize=7, color="#666666",
                 ha="right", va="bottom")

    # (c) simulation input profile -------------------------------------------
    ax = axes[2]
    rows = np.loadtxt(case["atm_file"], comments="#")      # TOA-first
    z_a, p_a, t_a = rows[:, 0], rows[:, 1], rows[:, 2]
    is_era5 = p_a >= prof["p_hpa"].min() - 1e-6
    ax.plot(t_a[~is_era5], p_a[~is_era5], "^-", ms=3.5, lw=0.8, color="#999999",
            label="afglsw splice (> ERA5 top)")
    ax.plot(t_a[is_era5], p_a[is_era5], "o-", ms=3, lw=1.0, color=C_ERA5,
            label="ERA5 levels + surface row")
    ax.set_yscale("log")
    ax.set_ylim(1100, 0.05)
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel("pressure (hPa)")
    n_total = int(manifest.get("n_plev_total", 37))
    n_above = prof["p_hpa"].size
    n_used = int(is_era5.sum()) - 1        # minus the added 2 m surface row
    setup = (f"levels: {n_total} plev; {n_total - n_above} below ground,\n"
             f"   {n_above - n_used} near-sfc dropped "
             f"$\\rightarrow$ {n_used} used + 2 m row\n"
             f"$\\epsilon$ = {eps}, sur_temperature = skt = {case['skt']:.1f} K\n"
             f"{w1:.2f}$-${w2:.0f} $\\mu$m, reptran coarse, DISORT\n"
             f"CO2/CH4: CAMS EGG4; O3, T, q: ERA5\n"
             f"pixel 0.25$^{{\\circ}}$ vs sonde: "
             f"{case['match_km']:.0f} km, {abs(case['match_dt_h']):.1f} h")
    if label == "cloudy":
        setup += (f"\nLWP {case['lwp_g']:.0f} / IWP {case['iwp_g']:.0f} "
                  f"g m$^{{-2}}$,\n"
                  f"   r$_{{eff}}$ {manifest['reff_um']['liquid']:.0f}/"
                  f"{manifest['reff_um']['ice']:.0f} $\\mu$m (Hu/Fu)")
    ax.text(0.0, -0.19, setup, transform=ax.transAxes, va="top", fontsize=7.2,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#cccccc"))
    ax.legend(frameon=False, fontsize=6.6, loc="lower left",
              bbox_to_anchor=(0.0, 1.0), borderaxespad=0.2)
    panel_label(ax, "c", x=-0.22, y=1.15)
    ax.set_title("simulation input", fontsize=9, pad=44)

    # (d) LW fluxes: sims vs ERA5 vs MOSAiC obs ------------------------------
    ax = axes[3]
    tail_up = (1 - planck_band_fraction(case["skt"], w1, w2)) * eps * SIGMA \
        * case["skt"] ** 4
    t_eff = (case["era5_lwdn"] / SIGMA) ** 0.25
    tail_dn = (1 - planck_band_fraction(t_eff, w1, w2)) * SIGMA * t_eff ** 4
    series = [
        ("MOSAiC obs", C_OBS, "*", 11, [obs["lwdn"], obs["lwup"]]),
        ("ERA5 (1-h mean)", C_ERA5, "o", 7, [case["era5_lwdn"], case["era5_lwup"]]),
        ("libRadtran", C_LIB, "s", 6, [res["lib_lwdn"], res["lib_lwup"]]),
        ("libRadtran + tail", C_LIB, "s", 6,
         [res["lib_lwdn"] + tail_dn, res["lib_lwup"] + tail_up]),
        ("RRTMG-LW", C_RRT, "^", 7, [rrtmg["sim_lwdn_sfc"], rrtmg["sim_lwup_sfc"]]),
    ]
    for i, (name, color, marker, ms, vals) in enumerate(series):
        x = np.array([0.0, 1.0]) + (i - 2) * 0.09
        open_face = "tail" in name
        ax.plot(x, vals, marker, ms=ms, color=color, ls="none",
                mfc="white" if open_face else color, mew=1.3, label=name)
        for xi, v in zip(x, vals):
            if np.isfinite(v):
                ax.annotate(f"{v:.0f}", (xi, v), textcoords="offset points",
                            xytext=(0, 7), ha="center", fontsize=6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["LW$\\downarrow$", "LW$\\uparrow$"])
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel("surface flux (W m$^{-2}$)")
    ax.legend(frameon=False, fontsize=6.6, loc="best")
    panel_label(ax, "d", x=-0.22, y=1.15)
    ax.set_title("surface LW fluxes", fontsize=9, pad=44)

    for ax in axes:
        ax.grid(color="#e3e3e3", lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    st = dt.datetime.fromisoformat(case["sounding_time"])
    et = dt.datetime.fromisoformat(case["era5_time"])
    fig.suptitle(
        f"MOSAiC case study — {label} sky: sonde {st:%Y-%m-%d %H:%M} UTC vs "
        f"ERA5 {et:%Y-%m-%d %H} UTC at ({case['pixel_lat']:.2f}$^\\circ$N, "
        f"{case['pixel_lon']:.2f}$^\\circ$E)", y=1.08)
    out = figures_dir(cfg) / f"case_study_{label}_{et:%Y%m%d}T{et:%H}Z.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")


def cmd_figure(args) -> int:
    import climlab
    from climlab.domain.axis import Axis
    import era5_rrtmg_sim as R

    cfg = load_config(args.config)
    manifest = json.loads(MANIFEST.read_text())
    results = json.loads(RESULTS.read_text())
    R.patch_interface_temperature()
    for case in manifest["cases"]:
        rrtmg = R.run_pixel(case, manifest, climlab, Axis)
        print(f"{case['label']:>6}: RRTMG LWdn {rrtmg['sim_lwdn_sfc']:7.2f}  "
              f"LWup {rrtmg['sim_lwup_sfc']:7.2f} W/m2")
        draw_case(cfg, manifest, case, results[case["label"]], rrtmg)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    for name, fn in (("prep", cmd_prep), ("run", cmd_run), ("figure", cmd_figure)):
        p = sub.add_parser(name)
        p.add_argument("--year", type=int, default=2020)
        p.add_argument("--month", type=int, default=1)
        p.add_argument("--config", default=None)
        p.set_defaults(func=fn)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
