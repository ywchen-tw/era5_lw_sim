"""Surface LW cloud radiative effect vs inversion strength (CRE study).

Simulates the surface longwave cloud radiative effect of every sampled
overcast column twice — all-sky (profile condensate as plane-parallel
liquid/ice clouds, fixed effective radii) and clear-sky (same column, clouds
removed) — with both libRadtran (uvspec/DISORT via er3t) and RRTMG-LW
(climlab), on the 75-90N MERRA-2 study domain (config_cre.yaml), and
stratifies CRE by inversion strength for 100% sea-ice vs 100% open-ocean
cells.

    CRE_dn  = LWdn_all - LWdn_clr          (surface downwelling CRE)
    CRE_net = (LWdn-LWup)_all - (LWdn-LWup)_clr = emissivity * CRE_dn
              (the skin temperature is held fixed between the paired runs)

Selection (per analysis instant, 00/06/12/18 UTC): sea cells in polar night
(geometric noon sun below the horizon) with bracket-mean FRSEAICE >= 0.99
(ice100) or <= 0.01 (ocean100), bracket-mean CLDTOT >= 0.99 and condensate
path > 1 g/m2 (the stage-7 plane-parallel screen), stratified-sampled across
fixed bins of dt_850_skt = T(850 hPa) - T(skin) — the study's inversion-
strength coordinate — so rare strata (e.g. inverted open ocean) are
populated. MERRA-2's own fluxes (LWGAB, LWGABCLR, EMIS) ride along as the
model's internal CRE estimate.

Subcommands (env in parentheses):
    select     (era5)     sample columns, write profiles + cloud files + manifest
    run-lrt    (er3t_env) paired uvspec thermal runs        -> results_libradtran.csv
    run-rrtmg  (era5)     paired climlab RRTMG-LW runs      -> results_rrtmg.csv
    analyze    (era5)     stats + figure stratified by inversion strength

Gases are held at the stage-7 constants (CO2 415 ppm, CH4 1.9 ppm) for all
columns: CRE is an all-sky-minus-clear difference of the same column, so the
gas description cancels almost exactly and month-to-month CAMS files are not
needed. libRadtran integrates 3.08-100 um; the far-IR >100 um tail largely
cancels in the CRE difference, and RRTMG (3.08-1000 um, no tail) quantifies
the residual.
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
from lrt_sim import (EMISSIVITY, REFF_ICE_UM, REFF_LIQ_UM, TRAPZ,
                     WVL_RANGE_NM, afglsw_o3_mmr, build_profile_files,
                     write_cloud_file)
from reanlib.classes import load_class_inputs
from reanlib.config import (area_tag, inversion_path, load_config, plev_path,
                            rad_path, sfc_path)
from reanlib.io_era5 import open_era5

G0 = 9.80665
HOURS = (0, 6, 12, 18)
SIC_ICE100_MIN = 0.99
SIC_OCEAN100_MAX = 0.01
TCC_OVERCAST_MIN = 0.99
TWP_MIN_G = 1.0
CLASSES = ("ice100", "ocean100")
# dt_850_skt bin edges (K): T(850) - T(skin). Negative = lapse (typical open
# ocean under cold-air advection), positive = inversion (typical pack ice).
DT_BIN_EDGES = (-np.inf, -4.0, 0.0, 4.0, 8.0, 12.0, np.inf)
DT_BIN_LABELS = ("<-4", "-4..0", "0..4", "4..8", "8..12", ">12")
CLASS_COLOURS = {"ice100": "#0072B2", "ocean100": "#D55E00"}


def cre_dir(cfg: dict) -> Path:
    from reanlib.config import REPO_ROOT
    root = Path(cfg["paths"]["derived"])
    if not root.is_absolute():
        root = REPO_ROOT / root
    return root / cfg["source"] / f"cre_sim{area_tag(cfg)}"


def season_tag(args) -> str:
    return f"{args.start.replace('-', '')[:6]}-{args.end.replace('-', '')[:6]}"


def manifest_file(cfg: dict, args) -> Path:
    return cre_dir(cfg) / f"manifest_{season_tag(args)}.json"


def solar_declination_deg(when: dt.datetime) -> float:
    """Spencer (1971) declination, degrees."""
    doy = when.timetuple().tm_yday
    g = 2.0 * np.pi / 365.0 * (doy - 1 + (when.hour - 12) / 24.0)
    d = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g)
         - 0.006758 * np.cos(2 * g) + 0.000907 * np.sin(2 * g)
         - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    return float(np.rad2deg(d))


def date_range(start: str, end: str):
    d = dt.date.fromisoformat(start)
    stop = dt.date.fromisoformat(end)
    while d <= stop:
        yield d
        d += dt.timedelta(days=1)


def bracket_mean_2d(da, when: np.datetime64) -> np.ndarray:
    """Mean of the two HH:30 windows bracketing an analysis instant."""
    stamps = [when - np.timedelta64(30, "m"), when + np.timedelta64(30, "m")]
    have = [s for s in stamps if s in da["valid_time"].values]
    return da.sel(valid_time=have).mean("valid_time").values


def dt_bin_index(x: float) -> int:
    return int(np.digitize([x], DT_BIN_EDGES[1:-1])[0])


# ---------------------------------------------------------------- select

def cmd_select(args) -> int:
    cfg = load_config(args.config)
    out_root = cre_dir(cfg)
    rng = np.random.default_rng(args.seed)

    profiles, counts = [], {}
    for date in date_range(args.start, args.end):
        pp = plev_path(cfg, date)
        if not pp.exists():
            print(f"{date}: no data, skipped")
            continue
        plev = open_era5(pp)
        sfc = open_era5(sfc_path(cfg, date))
        inv = open_era5(inversion_path(cfg, date))
        rad = open_era5(rad_path(cfg, date))
        times = np.array([np.datetime64(f"{date}T{h:02d}:00") for h in HOURS])
        lsm, sic_all, tcc_all = load_class_inputs(cfg, date, times)

        p = plev["pressure_level"].values.astype(float)      # descending hPa
        p_pa_asc = p[::-1] * 100.0
        lat = plev["latitude"].values
        lon = plev["longitude"].values
        lat2 = np.broadcast_to(lat[:, None], lsm.shape[1:])
        decl = solar_declination_deg(dt.datetime(date.year, date.month,
                                                 date.day, 12))
        night = (np.sin(np.deg2rad(lat2)) * np.sin(np.deg2rad(decl))
                 + np.cos(np.deg2rad(lat2)) * np.cos(np.deg2rad(decl))) <= 0.0

        prof_dir = out_root / "profiles" / f"{date:%Y%m%d}"
        for k, (hour, when) in enumerate(zip(HOURS, times)):
            pl = plev.sel(valid_time=when).transpose("pressure_level",
                                                     "latitude", "longitude")
            sf = sfc.sel(valid_time=when)
            iv = inv.sel(valid_time=when)
            sic, tcc = sic_all[k], tcc_all[k]

            clwc = np.nan_to_num(pl["clwc"].values)
            ciwc = np.nan_to_num(pl["ciwc"].values)
            lwp = np.abs(TRAPZ(clwc[::-1], x=p_pa_asc, axis=0)) / G0 * 1e3
            iwp = np.abs(TRAPZ(ciwc[::-1], x=p_pa_asc, axis=0)) / G0 * 1e3
            dt_skt = iv["dt_850_skt"].values
            sea = (lsm[k] < 0.5) & np.isfinite(sic)
            base = (sea & night & (tcc >= TCC_OVERCAST_MIN)
                    & (lwp + iwp > TWP_MIN_G) & np.isfinite(dt_skt))
            masks = {"ice100": base & (sic >= SIC_ICE100_MIN),
                     "ocean100": base & (sic <= SIC_OCEAN100_MAX)}

            lw_dn = bracket_mean_2d(rad["lwgab"], when)
            lw_dn_clr = bracket_mean_2d(rad["lwgabclr"], when)
            emis2d = bracket_mean_2d(rad["emis"], when)

            for cname, cmask in masks.items():
                for b in range(len(DT_BIN_LABELS)):
                    lo, hi = DT_BIN_EDGES[b], DT_BIN_EDGES[b + 1]
                    iys, ixs = np.where(cmask & (dt_skt > lo) & (dt_skt <= hi))
                    if iys.size == 0:
                        continue
                    take = rng.choice(iys.size, size=min(args.per_bin,
                                                         iys.size),
                                      replace=False)
                    for s in take:
                        iy, ix = int(iys[s]), int(ixs[s])
                        pid = f"{date:%Y%m%d}T{hour:02d}_{cname[:3]}_{iy:02d}_{ix:03d}"
                        prof = build_column(pl, sf, iv, p, iy, ix, pid,
                                            prof_dir)
                        if prof is None:
                            continue
                        prof.update({
                            "class": cname, "date": f"{date:%Y-%m-%d}",
                            "hour": hour,
                            "lat": float(lat[iy]), "lon": float(lon[ix]),
                            "sic": float(sic[iy, ix]),
                            "cldtot": float(tcc[iy, ix]),
                            "dt_bin": b,
                            "merra2_lwdn": float(lw_dn[iy, ix]
                                                 / emis2d[iy, ix]),
                            "merra2_lwdn_clr": float(lw_dn_clr[iy, ix]
                                                     / emis2d[iy, ix]),
                            "merra2_emis": float(emis2d[iy, ix]),
                        })
                        profiles.append(prof)
                        counts[cname, b] = counts.get((cname, b), 0) + 1
        for ds in (plev, sfc, inv, rad):
            ds.close()
        if date.day in (1, 15) or date == dt.date.fromisoformat(args.end):
            print(f"{date}: {len(profiles)} columns so far")

    manifest = {
        "season": season_tag(args),
        "source": cfg["source"], "area_tag": area_tag(cfg).lstrip("_"),
        "emissivity": EMISSIVITY,
        "wavelength_range_nm": list(WVL_RANGE_NM),
        "reff_um": {"liquid": REFF_LIQ_UM, "ice": REFF_ICE_UM},
        "gas_source": "constants: CO2 415 ppm, CH4 1.9 ppm (stage-7 values; "
                      "cancels in the all-minus-clear CRE difference)",
        "selection": {
            "sic_ice100_min": SIC_ICE100_MIN,
            "sic_ocean100_max": SIC_OCEAN100_MAX,
            "tcc_overcast_min": TCC_OVERCAST_MIN, "twp_min_g": TWP_MIN_G,
            "polar_night": "geometric noon sun below horizon",
            "strat_metric": "dt_850_skt",
            "dt_bin_edges": [float(e) for e in DT_BIN_EDGES[1:-1]],
            "per_bin_per_snapshot": args.per_bin, "seed": args.seed,
        },
        "profiles": profiles,
    }
    out_root.mkdir(parents=True, exist_ok=True)
    mpath = manifest_file(cfg, args)
    mpath.write_text(json.dumps(manifest, indent=1))
    print(f"\nwrote {mpath}  ({len(profiles)} columns)")
    print(f"{'bin (K)':>10s}  " + "  ".join(f"{c:>8s}" for c in CLASSES))
    for b, lab in enumerate(DT_BIN_LABELS):
        print(f"{lab:>10s}  " + "  ".join(
            f"{counts.get((c, b), 0):8d}" for c in CLASSES))
    return 0


def build_column(pl, sf, iv, p, iy, ix, pid, prof_dir):
    """Profile + cloud files for one cell; returns the profile dict or None."""
    sp_hpa = float(sf["sp"].values[iy, ix]) / 100.0
    t_col_all = pl["t"].values[:, iy, ix].astype(float)
    above = (p <= sp_hpa) & np.isfinite(t_col_all)
    if above.sum() < 20:                     # degenerate column
        return None
    skt = float(sf["skt"].values[iy, ix])
    t2m = float(sf["t2m"].values[iy, ix])
    p_col = p[above]
    t_col = t_col_all[above]
    q_col = pl["q"].values[above, iy, ix].astype(float)
    o3_col = pl["o3"].values[above, iy, ix].astype(float)
    o3_col = np.where(np.isfinite(o3_col), o3_col, afglsw_o3_mmr(p_col))
    clwc = np.nan_to_num(pl["clwc"].values[above, iy, ix].astype(float))
    ciwc = np.nan_to_num(pl["ciwc"].values[above, iy, ix].astype(float))

    prof_dir.mkdir(parents=True, exist_ok=True)
    atm_file, ch4_file, keep, z_col = build_profile_files(
        prof_dir, pid, p_col, t_col, q_col, o3_col, t2m, sp_hpa, None)

    lwp = float(np.abs(TRAPZ(clwc[::-1], x=p_col[::-1] * 100.0)) / G0 * 1e3)
    iwp = float(np.abs(TRAPZ(ciwc[::-1], x=p_col[::-1] * 100.0)) / G0 * 1e3)
    wc_file = write_cloud_file(prof_dir / f"wc_{pid}.dat", z_col,
                               clwc[keep], p_col[keep], t_col[keep],
                               REFF_LIQ_UM)
    ic_file = write_cloud_file(prof_dir / f"ic_{pid}.dat", z_col,
                               ciwc[keep], p_col[keep], t_col[keep],
                               REFF_ICE_UM)
    if wc_file is None and ic_file is None:   # condensate below file threshold
        return None

    # cloud base/top from the condensate profile (mechanism diagnostics)
    cld = (clwc + ciwc) > 1e-7
    icld = np.where(cld)[0]
    p_base, p_top = float(p_col[icld].max()), float(p_col[icld].min())
    t_base = float(t_col[icld[np.argmax(p_col[icld])]])
    t_top = float(t_col[icld[np.argmin(p_col[icld])]])

    return {
        "label": pid, "sp_hpa": sp_hpa, "skt": skt, "t2m": t2m,
        "dt_850_skt": float(iv["dt_850_skt"].values[iy, ix]),
        "dt_850_2m": float(iv["dt_850_2m"].values[iy, ix]),
        "sbi_strength": float(np.nan_to_num(
            iv["sbi_strength"].values[iy, ix])),
        "sbi_found": int(iv["sbi_found"].values[iy, ix]),
        "lwp_g": lwp, "iwp_g": iwp,
        "cld_base_hpa": p_base, "cld_top_hpa": p_top,
        "cld_base_t": t_base, "cld_top_t": t_top,
        "atm_file": atm_file, "ch4_file": ch4_file,
        "wc_file": wc_file, "ic_file": ic_file,
    }


# ---------------------------------------------------------------- run-lrt

def cmd_run_lrt(args) -> int:
    import copy
    import inspect

    import er3t
    import pandas as pd

    cfg = load_config(args.config)
    manifest = json.loads(manifest_file(cfg, args).read_text())
    tmp_dir = cre_dir(cfg) / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    rpath = cre_dir(cfg) / f"results_libradtran_{season_tag(args)}.csv"

    if args.streams is None:
        args.streams = 4 if platform.system() == "Darwin" else 8
    lrt_cfg_base = er3t.rtm.lrt.get_lrt_cfg()
    lrt_cfg_base["mol_abs_param"] = "reptran " + args.mol_abs_param
    lrt_cfg_base["number_of_streams"] = args.streams
    mute_list = ["albedo", "wavelength", "spline", "source solar",
                 "slit_function_file"]
    accepted = inspect.signature(
        er3t.rtm.lrt.lrt_init_mono_flx.__init__).parameters

    done = set()
    if rpath.exists() and not args.overwrite:
        done = set(pd.read_csv(rpath)["label"])

    inits, meta = [], []
    for prof in manifest["profiles"]:
        if prof["label"] in done:
            continue
        when = dt.datetime.fromisoformat(prof["date"]) \
            + dt.timedelta(hours=prof["hour"])
        for mode in ("all", "clr"):
            in_path = str(tmp_dir / f"in_{prof['label']}_{mode}.txt")
            out_path = str(tmp_dir / f"out_{prof['label']}_{mode}.txt")
            lrt_cfg = copy.deepcopy(lrt_cfg_base)
            lrt_cfg["atmosphere_file"] = prof["atm_file"]
            extra = {
                "source": "thermal",
                "albedo_add": f"{1.0 - manifest['emissivity']:.4f}",
                "sur_temperature": f"{prof['skt']:.2f}",
                "wavelength_add": "{:.0f} {:.0f}".format(
                    *manifest["wavelength_range_nm"]),
                "output_process": "integrate",
                "mol_file": "CH4 " + prof["ch4_file"],
            }
            if mode == "all":
                if prof.get("wc_file"):
                    extra["wc_file 1D"] = prof["wc_file"]
                    extra["wc_properties"] = "hu interpolate"
                if prof.get("ic_file"):
                    extra["ic_file 1D"] = prof["ic_file"]
                    extra["ic_properties"] = "fu interpolate"
            kwargs = dict(input_file=in_path, output_file=out_path,
                          date=when, solar_zenith_angle=80.0, Nx=1,
                          output_altitude=[0, "toa"],
                          input_dict_extra=extra, mute_list=mute_list,
                          lrt_cfg=lrt_cfg, cld_cfg=None, aer_cfg=None)
            init = er3t.rtm.lrt.lrt_init_mono_flx(
                **{k: v for k, v in kwargs.items() if k in accepted})
            inits.append(init)
            meta.append((prof["label"], mode, out_path))

    workers = args.workers or max((os.cpu_count() or 2) - 2, 1)
    print(f"{len(inits)} uvspec jobs to run ({len(meta) // 2} columns x 2), "
          f"reptran {args.mol_abs_param}, {args.streams} streams, "
          f"{workers} workers")
    chunk = 200
    for i0 in range(0, len(inits), chunk):
        batch = [ini for ini in inits[i0:i0 + chunk]
                 if args.overwrite or not (
                     os.path.exists(ini.output_file)
                     and os.path.getsize(ini.output_file) > 50)]
        if batch:
            er3t.rtm.lrt.lrt_run_mp(batch, Ncpu=min(workers, len(batch)))
        rows = {}
        for (label, mode, out_path), ini in zip(meta[i0:i0 + chunk],
                                                inits[i0:i0 + chunk]):
            data = er3t.rtm.lrt.lrt_read_uvspec_flx([ini])
            f_dn = np.squeeze(data.f_down) * 1000.0     # [surface, toa]
            f_up = np.squeeze(data.f_up) * 1000.0
            r = rows.setdefault(label, {"label": label})
            r[f"lwdn_{mode}"] = float(f_dn[0])
            r[f"lwup_{mode}"] = float(f_up[0])
            r[f"olr_{mode}"] = float(f_up[1])
        df = pd.DataFrame(list(rows.values()))
        header = not rpath.exists()
        df.to_csv(rpath, mode="a", header=header, index=False)
        print(f"  {min(i0 + chunk, len(inits))}/{len(inits)} jobs done")
    print(f"wrote {rpath}")
    return 0


# ---------------------------------------------------------------- run-rrtmg

def cmd_run_rrtmg(args) -> int:
    import pandas as pd

    try:
        import climlab
        from climlab.domain.axis import Axis
    except ImportError:
        sys.exit("climlab not importable — run under the `era5` conda env")
    from rrtmg_sim import patch_interface_temperature, run_pixel
    patch_interface_temperature()

    cfg = load_config(args.config)
    manifest = json.loads(manifest_file(cfg, args).read_text())
    rpath = cre_dir(cfg) / f"results_rrtmg_{season_tag(args)}.csv"
    done = set()
    if rpath.exists() and not args.overwrite:
        done = set(pd.read_csv(rpath)["label"])

    rows = []
    todo = [p for p in manifest["profiles"] if p["label"] not in done]
    print(f"RRTMG-LW on {len(todo)} columns x 2 runs")
    for i, prof in enumerate(todo):
        out_all = run_pixel(prof, manifest, climlab, Axis)
        p_clr = {k: v for k, v in prof.items()
                 if k not in ("wc_file", "ic_file")}
        out_clr = run_pixel(p_clr, manifest, climlab, Axis)
        rows.append({
            "label": prof["label"],
            "lwdn_all": out_all["sim_lwdn_sfc"],
            "lwup_all": out_all["sim_lwup_sfc"],
            "olr_all": out_all["sim_olr_toa"],
            "lwdn_clr": out_clr["sim_lwdn_sfc"],
            "lwup_clr": out_clr["sim_lwup_sfc"],
            "olr_clr": out_clr["sim_olr_toa"],
        })
        if (i + 1) % 200 == 0 or i + 1 == len(todo):
            df = pd.DataFrame(rows)
            df.to_csv(rpath, mode="a", header=not rpath.exists(), index=False)
            rows = []
            print(f"  {i + 1}/{len(todo)} columns done")
    print(f"wrote {rpath}")
    return 0


# ---------------------------------------------------------------- analyze

def load_merged(cfg, args, keep_files: bool = False):
    """Manifest + both result CSVs merged into one DataFrame with CRE columns."""
    import pandas as pd

    manifest = json.loads(manifest_file(cfg, args).read_text())
    prof = pd.DataFrame([{k: v for k, v in p.items()
                          if keep_files or not str(k).endswith("_file")}
                         for p in manifest["profiles"]])
    tag = season_tag(args)
    frames = {"libradtran": f"results_libradtran_{tag}.csv",
              "rrtmg": f"results_rrtmg_{tag}.csv"}
    for code, fname in frames.items():
        fp = cre_dir(cfg) / fname
        if not fp.exists():
            sys.exit(f"missing {fp} — run the run-{code.replace('libradtran', 'lrt')} "
                     f"subcommand first")
        r = pd.read_csv(fp).drop_duplicates("label", keep="last")
        r[f"cre_dn_{code}"] = r["lwdn_all"] - r["lwdn_clr"]
        r[f"cre_net_{code}"] = ((r["lwdn_all"] - r["lwup_all"])
                                - (r["lwdn_clr"] - r["lwup_clr"]))
        prof = prof.merge(
            r[["label", f"cre_dn_{code}", f"cre_net_{code}",
               "lwdn_all", "lwdn_clr"]].rename(columns={
                   "lwdn_all": f"lwdn_all_{code}",
                   "lwdn_clr": f"lwdn_clr_{code}"}),
            on="label", how="inner")
    prof["cre_dn_merra2"] = prof["merra2_lwdn"] - prof["merra2_lwdn_clr"]
    # MERRA-2 net CRE: (LWGAB - LWGABCLR) - the emitted flux cancels (same
    # skin state in the model's all-sky and clear-sky diagnostic calls)
    prof["cre_net_merra2"] = prof["cre_dn_merra2"] * prof["merra2_emis"]
    prof["twp_g"] = prof["lwp_g"] + prof["iwp_g"]
    return manifest, prof


def cmd_analyze(args) -> int:
    import matplotlib.pyplot as plt

    from reanlib.plotstyle import apply_agu_style, panel_label

    cfg = load_config(args.config)
    manifest, prof = load_merged(cfg, args)
    tag = season_tag(args)
    q = f"cre_{args.flux}"                                   # plotted quantity
    qlab = {"dn": "CRE$_{dn}$", "net": "CRE$_{net}$"}[args.flux]
    print(f"{len(prof)} columns with both simulators (plotting {q})")

    # ---------------- stats table
    stats = {}
    lines = [f"{'class':9s} {'bin (K)':>8s} {'n':>5s} "
             f"{'CRE_dn lrt':>12s} {'CRE_dn rrtmg':>13s} {'CRE_dn M2':>10s} "
             f"{'TWP med':>8s} {'LWP med':>8s}"]
    for cname in CLASSES:
        for b, lab in enumerate(DT_BIN_LABELS):
            d = prof[(prof["class"] == cname) & (prof["dt_bin"] == b)]
            if not len(d):
                continue
            e = {
                "n": int(len(d)),
                "cre_dn_lrt_med": float(d["cre_dn_libradtran"].median()),
                "cre_dn_lrt_q25": float(d["cre_dn_libradtran"].quantile(.25)),
                "cre_dn_lrt_q75": float(d["cre_dn_libradtran"].quantile(.75)),
                "cre_dn_rrtmg_med": float(d["cre_dn_rrtmg"].median()),
                "cre_dn_merra2_med": float(d["cre_dn_merra2"].median()),
                "cre_net_lrt_med": float(d["cre_net_libradtran"].median()),
                "twp_med": float(d["twp_g"].median()),
                "lwp_med": float(d["lwp_g"].median()),
                "dt_850_skt_med": float(d["dt_850_skt"].median()),
            }
            stats[f"{cname}|{lab}"] = e
            lines.append(
                f"{cname:9s} {lab:>8s} {e['n']:5d} "
                f"{e['cre_dn_lrt_med']:12.1f} {e['cre_dn_rrtmg_med']:13.1f} "
                f"{e['cre_dn_merra2_med']:10.1f} "
                f"{e['twp_med']:8.1f} {e['lwp_med']:8.1f}")
    print("\n".join(lines))
    spath = cre_dir(cfg) / f"cre_stats_{tag}.json"
    spath.write_text(json.dumps(
        {"season": tag, "strat_metric": "dt_850_skt",
         "units": "W m-2 (CRE_dn); K (bins)", "bins": list(DT_BIN_LABELS),
         "stats": stats}, indent=1))
    print(f"wrote {spath}")

    # ---------------- figure
    apply_agu_style()
    from reanlib.config import REPO_ROOT, figures_dir
    fdir = figures_dir(cfg)
    fdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 6.2))
    (ax_a, ax_b), (ax_c, ax_d) = axes

    # (a) binned median CRE vs dt_850_skt
    centers = {}
    for cname in CLASSES:
        col = CLASS_COLOURS[cname]
        for code, ls, mk in (("libradtran", "-", "o"), ("rrtmg", "--", "s"),
                             ("merra2", ":", "^")):
            xs, ys, lo, hi = [], [], [], []
            for b in range(len(DT_BIN_LABELS)):
                d = prof[(prof["class"] == cname) & (prof["dt_bin"] == b)]
                if len(d) < 5:
                    continue
                xs.append(d["dt_850_skt"].median())
                v = d[f"{q}_{code}"]
                ys.append(v.median())
                lo.append(v.quantile(.25))
                hi.append(v.quantile(.75))
            if code == "libradtran":
                ax_a.fill_between(xs, lo, hi, color=col, alpha=0.15, lw=0)
                centers[cname] = xs
            ax_a.plot(xs, ys, ls, marker=mk, ms=3.5, color=col, lw=1.2,
                      label=f"{cname} {code}" if code != "merra2"
                      else f"{cname} MERRA-2")
    ax_a.axvline(0, color="0.75", lw=0.6, zorder=0)
    ax_a.set_xlabel("T(850 hPa) $-$ T(skin) (K)")
    ax_a.set_ylabel(f"surface {qlab} (W m$^{{-2}}$)")
    ax_a.legend(fontsize=6, ncol=2, frameon=False)
    panel_label(ax_a, "a", outside=True)

    # (b) CRE distributions per class (libRadtran)
    for i, cname in enumerate(CLASSES):
        d = prof[prof["class"] == cname][f"{q}_libradtran"]
        ax_b.hist(d, bins=np.arange(-5, 105, 2.5), density=True,
                  histtype="stepfilled", alpha=0.45,
                  color=CLASS_COLOURS[cname],
                  label=f"{cname} (n={len(d)}, med {d.median():.1f})")
    ax_b.set_xlabel(f"surface {qlab}, libRadtran (W m$^{{-2}}$)")
    ax_b.set_ylabel("density")
    ax_b.legend(fontsize=6, frameon=False)
    panel_label(ax_b, "b", outside=True)

    # (c) CRE vs total water path, colored by inversion strength
    sc = None
    for cname, mk in (("ice100", "o"), ("ocean100", "^")):
        d = prof[prof["class"] == cname]
        sc = ax_c.scatter(d["twp_g"], d[f"{q}_libradtran"], s=4, marker=mk,
                          c=d["dt_850_skt"], cmap="coolwarm", vmin=-12,
                          vmax=12, lw=0, alpha=0.7)
    ax_c.set_xscale("log")
    ax_c.set_xlabel("LWP + IWP (g m$^{-2}$)")
    ax_c.set_ylabel(f"surface {qlab} (W m$^{{-2}}$)")
    cb = plt.colorbar(sc, ax=ax_c, pad=0.02)
    cb.set_label("T(850) $-$ T(skin) (K)", fontsize=7)
    ax_c.text(0.03, 0.95, "circles ice100, triangles ocean100",
              transform=ax_c.transAxes, fontsize=6, va="top")
    panel_label(ax_c, "c", outside=True)

    # (d) code comparison
    for cname in CLASSES:
        d = prof[prof["class"] == cname]
        ax_d.scatter(d[f"{q}_libradtran"], d[f"{q}_rrtmg"], s=4,
                     color=CLASS_COLOURS[cname], lw=0, alpha=0.6,
                     label=cname)
    lim = [-5, float(prof[[f"{q}_libradtran", f"{q}_rrtmg"]].max().max())
           + 5]
    ax_d.plot(lim, lim, color="0.6", lw=0.7)
    bias = (prof[f"{q}_rrtmg"] - prof[f"{q}_libradtran"]).mean()
    ax_d.text(0.03, 0.95, f"RRTMG $-$ libRadtran: {bias:+.1f} W m$^{{-2}}$",
              transform=ax_d.transAxes, fontsize=6, va="top")
    ax_d.set_xlabel(f"{qlab} libRadtran (W m$^{{-2}}$)")
    ax_d.set_ylabel(f"{qlab} RRTMG (W m$^{{-2}}$)")
    ax_d.legend(fontsize=6, frameon=False)
    panel_label(ax_d, "d", outside=True)

    fig.suptitle(f"MERRA-2 {manifest['area_tag']}  {tag}  overcast "
                 f"polar-night columns — surface LW {qlab} vs inversion "
                 f"strength", fontsize=8, y=0.995)
    fig.tight_layout()
    suffix = "" if args.flux == "dn" else f"_{args.flux}"
    fpath = fdir / f"cre_inversion{suffix}_{tag}.png"
    fig.savefig(fpath, dpi=300)
    print(f"wrote {fpath}")
    return 0


# ------------------------------------------------------------- analyze-cloud

# geometric-optics visible extinction optical depth tau = 3 W / (2 rho reff)
# with the fixed simulation radii; rho_ice = 0.917 g/cm3 (Fu's D_ge nuances
# are ignored — this is a descriptive coordinate, not an RT input)
TAU_PER_G_LIQ = 3.0 / (2.0 * 1.0 * REFF_LIQ_UM)      # per g/m2 of LWP
TAU_PER_G_ICE = 3.0 / (2.0 * 0.917 * REFF_ICE_UM)    # per g/m2 of IWP


# ---- multi-layer segmentation, calibrated against the Cloudnet radar/lidar
# climatology at Ny-Alesund (Nomokonova et al. 2019, ACP 19, 4105-4126,
# doi:10.5194/acp-19-4105-2019; 94 GHz radar + ceilometer + MWR). Observed
# SINGLE-LAYER geometric thickness there:
#     liquid       median 0.24 km, <1% above 0.8 km
#     mixed-phase  median 1.1  km, <1% above 3.0 km (range to 8.5 km)
#     ice          median 1.5  km, <1% above 4.2 km (range to 9.5 km)
# i.e. multi-km cloud thickness is a ~1% tail, not p90. Deep Arctic systems
# are typically a thin liquid(-mixed) layer seeding ice/snow below it
# (seeder-feeder; Arctic multi-layer seeding statistics: ACP 26, 3049, 2026,
# doi:10.5194/acp-26-3049-2026); Cloudnet counts that falling ice as
# precipitating ice, NOT cloud thickness. The segmentation mirrors this:
# ice-only condensate contiguous below a liquid layer's base is VIRGA —
# excluded from the layer's thickness, kept as a separate attribute (it
# still emits: the per-layer RT validation includes it). Reanalyses also
# misplace Arctic cloud water vertically (Graham et al. 2019, J. Climate 32,
# doi:10.1175/JCLI-D-18-0643.1; Yeo et al. 2022, Atmos. Res. 270,
# doi:10.1016/j.atmosres.2022.106080), and 0.5-degree grid-mean condensate
# cannot reproduce radar layering exactly — the liquid anchor minimizes
# that damage because liquid, unlike ice, is vertically concentrated.
GAP_SPLIT_KM = 1.0
# ANY clear level splits (Cloudnet splits at any clear 30 m radar bin; one
# MERRA-2 level is already 200-900 m of clear air)
GAP_SPLIT_LEVELS = 1
TAU_SIG_MIN = 0.05
# a level only counts as CLOUD above this water content (~ condensate mixing
# ratio 1e-6 kg/kg, the common model cloud-boundary criterion). Sensitivity
# scan (Nov 2019-Jan 2020 sample): thickness p90 moves only 4.5 -> 4.1 km
# from 1e-5 to 1e-3 g/m3 — the deep tail is dense cyclone cloud, not wisps.
CLOUD_CONTENT_MIN_GM3 = 1e-3
# liquid content that anchors a cloud layer (typical Arctic liquid layers
# carry 0.05-0.3 g/m3; 0.01 keeps weak but genuine liquid bases)
LIQ_CORE_MIN_GM3 = 0.01
# split a contiguous run at an interior content minimum this far below both
# flanking peaks (grid-mean condensate rarely reaches exactly zero between
# genuinely separate decks, so the gap rule alone cannot see the boundary);
# 0.3 = a 3x dip separates decks
VALLEY_FRAC = 0.3
# tau-effective band (kept as secondary diagnostics dz_eff/z_top_eff)
TAU_EFF_BAND = (0.05, 0.95)
# LW absorption optical depth ~ half the visible extinction tau (thermal
# single-scattering albedo is small); the layer's emission depth dz_emit is
# the distance from its bottommost condensate to cumulative tau_abs = 1
LW_ABS_PER_VIS = 0.5
# Ny-Alesund Cloudnet single-layer thickness benchmarks (km): median and
# ~99th percentile per phase (Nomokonova et al. 2019, Sect. 4/Fig. 6) —
# printed next to the model numbers by analyze-cloud as a calibration check.
# Caveats when reading that check: (1) the obs EXCLUDE liquid-precipitating
# profiles, so the model's deep drizzling marine decks inflate the "liquid"
# row (they would not enter the obs statistic at all) — the mixed-phase row
# is the meaningful anchor; (2) 0.5-degree grid-mean condensate smears radar
# layering, so a residual high bias in the ice row is expected and
# irreducible (see the reference block above).
NYA_OBS_THICKNESS = {"liquid": (0.24, 0.8), "mixed": (1.1, 3.0),
                     "ice": (1.5, 4.2)}
# SEASON-MATCHED (Nov-Jan) benchmarks, computed by src/nya_cloudnet_winter.py
# from the ACTRIS Cloudnet legacy target classification at the same site
# (cloudnet.fmi.fi, CC-BY 4.0; 182 days of Nov 2016 - Jan 2017 + Nov 2017 -
# Jan 2018; 120,264 single-layer profiles; same definitions as the paper,
# liquid-precipitating profiles excluded). WINTER clouds are much deeper
# than the annual climatology: mixed 1.64 km median / 6.0 km p99, ice
# 1.49 / 7.9 km — multi-km winter cloud IS observed at this site, so
# MERRA-2's deep tail is realistic; its residual biases are the ice-phase
# median (~2x) and the drizzling liquid decks (not in the obs statistic).
NYA_OBS_WINTER = {"liquid": (0.16, 0.54), "mixed": (1.64, 6.02),
                  "ice": (1.49, 7.88)}


def read_cloud_ladder(wc_file, ic_file):
    """Ascending-z ladder (km) with liquid/ice water content (g/m3).

    The wc/ic files the RT consumed share one z ladder (zeros below the
    1e-5 g/m3 write floor), so the geometry describes the simulated cloud.
    """
    z = wc = ic = None
    for f, kind in ((wc_file, "wc"), (ic_file, "ic")):
        if not isinstance(f, str):           # None/NaN: no file of this phase
            continue
        rows = np.atleast_2d(np.loadtxt(f, comments="#"))[::-1]  # ascending z
        if z is None:
            z = rows[:, 0]
            wc = np.zeros(z.size)
            ic = np.zeros(z.size)
        if kind == "wc":
            wc = rows[:, 1]
        else:
            ic = rows[:, 1]
    return z, wc, ic


def _valley_split(g, c):
    """Recursively split an index run at content minima far below both
    flanking peaks (VALLEY_FRAC); the valley level joins the lower deck."""
    if g.size < 3:
        return [g]
    vals = c[g]
    i_min = int(np.argmin(vals[1:-1])) + 1
    if vals[i_min] <= VALLEY_FRAC * min(vals[:i_min].max(),
                                        vals[i_min + 1:].max()):
        return _valley_split(g[:i_min + 1], c) + _valley_split(g[i_min + 1:],
                                                               c)
    return [g]


def segment_layers(z, wc, ic, content_min=CLOUD_CONTENT_MIN_GM3):
    """Cloud layers (ascending): gap rule + valley splitting, with geometry,
    water paths, tau, and tau-effective boundaries per layer."""
    dz = np.empty_like(z)                     # cell-centred level thickness
    dz[1:-1] = (z[2:] - z[:-2]) / 2.0
    dz[0], dz[-1] = z[1] - z[0], z[-1] - z[-2]
    wc = np.where(wc >= content_min, wc, 0.0)
    ic = np.where(ic >= content_min, ic, 0.0)
    c = wc + ic
    idx = np.where(c > 0)[0]
    if not idx.size:
        return []
    groups = [[int(idx[0])]]
    for i_prev, i in zip(idx[:-1], idx[1:]):
        n_clear = i - i_prev - 1
        if n_clear >= GAP_SPLIT_LEVELS or (
                n_clear >= 1 and z[i] - z[i_prev] >= GAP_SPLIT_KM):
            groups.append([int(i)])
        else:
            groups[-1].append(int(i))
    groups = [h for g in groups for h in _valley_split(np.asarray(g), c)]
    layers = []
    for g in groups:
        # phase-aware virga trim (see the reference block above): the CLOUD
        # part of a liquid-anchored layer spans lowest to highest liquid-core
        # level; contiguous ice-only condensate below it is virga, above it
        # seeding/anvil ice. A layer with no liquid core is an ice cloud and
        # keeps its full extent (obs ice layers: median 1.5 km).
        liq = g[wc[g] >= LIQ_CORE_MIN_GM3]
        if liq.size:
            cloud = g[(g >= liq[0]) & (g <= liq[-1])]
            virga = g[g < liq[0]]
            above = g[g > liq[-1]]
        else:
            cloud, virga, above = g, g[:0], g[:0]

        def paths(sel):
            return (float(np.sum(wc[sel] * dz[sel]) * 1000.0),   # g/m2
                    float(np.sum(ic[sel] * dz[sel]) * 1000.0))
        lwp_c, iwp_c = paths(cloud)
        lwp_t, iwp_t = paths(g)
        phase = ("ice" if not liq.size else
                 "liquid" if iwp_c <= 0.1 * (lwp_c + iwp_c) else "mixed")

        # tau-effective boundaries within the cloud part (secondary metric)
        tau_lvl = (TAU_PER_G_LIQ * wc[cloud] + TAU_PER_G_ICE * ic[cloud]) \
            * dz[cloud] * 1000.0
        cum = np.cumsum(tau_lvl)
        i_lo = int(np.searchsorted(cum, TAU_EFF_BAND[0] * cum[-1]))
        i_hi = int(np.searchsorted(cum, TAU_EFF_BAND[1] * cum[-1]))

        # emission depth: from the layer's bottommost condensate (virga
        # included — it emits) up to cumulative LW-absorption tau = 1
        abs_lvl = LW_ABS_PER_VIS * (TAU_PER_G_LIQ * wc[g]
                                    + TAU_PER_G_ICE * ic[g]) * dz[g] * 1000.0
        cum_abs = np.cumsum(abs_lvl)
        i_em = int(np.searchsorted(cum_abs, 1.0)) if cum_abs[-1] >= 1.0 \
            else g.size - 1
        layers.append({
            "z_base": float(z[cloud[0]]), "z_top": float(z[cloud[-1]]),
            "dz": float(z[cloud[-1]] - z[cloud[0]]),
            "z_top_eff": float(z[cloud[i_hi]]),
            "dz_eff": float(z[cloud[i_hi]] - z[cloud[i_lo]]),
            "phase": phase,
            "virga_dz": float(z[cloud[0]] - z[g[0]]),
            "virga_iwp": paths(virga)[1],
            "ice_above_iwp": paths(above)[1],
            "dz_emit": float(z[g[i_em]] - z[g[0]]),
            "lwp": lwp_t, "iwp": iwp_t,
            # significance uses ALL condensate assigned to the layer (a
            # strong virga under a weak liquid base still counts)
            "tau": TAU_PER_G_LIQ * lwp_t + TAU_PER_G_ICE * iwp_t,
            "idx": g,
        })
    return layers


def significant(layers, tau_min=TAU_SIG_MIN):
    """Layers that count as CLOUD: tau >= TAU_SIG_MIN (= 0.05, i.e. LWP
    >= 0.33 g/m2 pure liquid or IWP >= 0.76 g/m2 pure ice at the fixed
    radii). Wisps below that are never counted as cloud layers."""
    return [l for l in layers if l["tau"] >= tau_min]


def layer_stats(prof, content_min=CLOUD_CONTENT_MIN_GM3,
                tau_min=TAU_SIG_MIN):
    """Per-column layering summary (lowest SIGNIFICANT layer + counts)."""
    import pandas as pd

    rows = []
    for wcf, icf in zip(prof["wc_file"], prof["ic_file"]):
        z, wc, ic = read_cloud_ladder(wcf, icf)
        sig = significant(segment_layers(z, wc, ic, content_min),
                          tau_min)
        if sig:
            low = sig[0]
            # geometric CLOUD-part thickness (virga excluded) — the
            # Cloudnet-comparable metric; dz_eff stays available in the
            # layer dicts as a secondary diagnostic
            rows.append({"n_sig": len(sig),
                         "z_top_col": sig[-1]["z_top"],
                         "dz_top": sig[-1]["dz"],
                         "z_top_low": low["z_top"],
                         "dz_low": low["dz"], "tau_low": low["tau"],
                         "z_base_low": low["z_base"],
                         "phase_low": low["phase"],
                         "virga_dz_low": low["virga_dz"],
                         "dz_emit_low": low["dz_emit"]})
        else:   # only sub-threshold wisps: no cloud layer at all
            rows.append({"n_sig": 0, "z_top_col": np.nan, "dz_top": np.nan,
                         "z_top_low": np.nan, "dz_low": np.nan,
                         "tau_low": np.nan, "z_base_low": np.nan,
                         "phase_low": None, "virga_dz_low": np.nan,
                         "dz_emit_low": np.nan})
    return pd.DataFrame(rows, index=prof.index)


def cmd_analyze_cloud(args) -> int:
    import matplotlib.pyplot as plt

    from reanlib.config import figures_dir
    from reanlib.plotstyle import apply_agu_style, panel_label

    cfg = load_config(args.config)
    manifest, prof = load_merged(cfg, args, keep_files=True)
    tag = season_tag(args)
    q = f"cre_{args.flux}"
    qlab = {"dn": "CRE$_{dn}$", "net": "CRE$_{net}$"}[args.flux]
    prof["tau"] = (TAU_PER_G_LIQ * prof["lwp_g"]
                   + TAU_PER_G_ICE * prof["iwp_g"])
    prof = prof.join(layer_stats(prof, args.content_min, args.tau_sig))
    print(f"{len(prof)} columns (plotting {q}); tau median "
          f"ice {prof[prof['class'] == 'ice100']['tau'].median():.1f}, "
          f"ocean {prof[prof['class'] == 'ocean100']['tau'].median():.1f}")
    print(f"  cloud thresholds: level >= {args.content_min:g} g/m3, "
          f"layer tau >= {args.tau_sig:g} (LWP >= "
          f"{args.tau_sig / TAU_PER_G_LIQ:.2f} / IWP >= "
          f"{args.tau_sig / TAU_PER_G_ICE:.2f} g/m2); "
          f"{int((prof['n_sig'] == 0).sum())} columns with no qualifying "
          f"cloud layer excluded from layer panels")
    for cname in CLASSES:
        d = prof[prof["class"] == cname]
        print(f"  {cname}: {float((d['n_sig'] >= 2).mean()):.1%} multi-layer "
              f"(2: {int((d['n_sig'] == 2).sum())}, "
              f"3+: {int((d['n_sig'] >= 3).sum())}); "
              f"lowest-layer thickness p50/p90/p99 = "
              f"{d['dz_low'].quantile(.5):.1f}/{d['dz_low'].quantile(.9):.1f}"
              f"/{d['dz_low'].quantile(.99):.1f} km; median virga depth "
              f"{d['virga_dz_low'].median():.1f} km")
    # NOTE on representativeness (full-season census-layers, 702k columns,
    # cos-lat weighted — supersedes an earlier 5-snapshot audit whose tail
    # fractions were unstable): area-true lowest-layer thickness p50/p90/p99
    # = 0.2/3.0/5.8 km over ice100 (10.6% of area > 3 km) and 1.6/4.4/5.4 km
    # over ocean100 (30.3% > 3 km). The stratified sample reproduces these
    # distributions closely (sample p90 2.8/3.9 km), so thickness statistics
    # here ARE representative; only stability-related occupancy reflects the
    # stratified weights. The >3 km frequency is a genuine feature of
    # MERRA-2's overcast winter columns — ~2x the Ny-Alesund radar
    # climatology per phase even after virga trimming (see census-layers
    # panel c) — consistent with grid-mean smearing plus the strict-overcast
    # conditioning (radar statistics include all cloudy scenes; overcast
    # instants are biased toward synoptic systems).
    # calibration check against the Ny-Alesund Cloudnet single-layer
    # thickness climatology (Nomokonova et al. 2019; see constants block)
    single = prof[prof["n_sig"] == 1]
    print("  single-layer thickness vs Ny-Alesund Cloudnet "
          "(model p50/p99 vs obs p50/~p99, km):")
    for ph, (o50, o99) in NYA_OBS_THICKNESS.items():
        d = single[single["phase_low"] == ph]["dz_low"]
        if len(d):
            print(f"    {ph:6s} (n={len(d):4d}): "
                  f"{d.quantile(.5):.2f}/{d.quantile(.99):.2f} "
                  f"vs {o50:.2f}/{o99:.1f}")

    apply_agu_style()
    fig, axes = plt.subplots(2, 3, figsize=(10.4, 7.2))
    ylab = f"surface {qlab} (W m$^{{-2}}$)"
    ybins = np.arange(0.0, 96.0, 2.5)

    # (a)/(b) LWP and IWP separately; symlog keeps the zero-path columns
    # (pure-ice or pure-liquid clouds) visible instead of dropping them.
    # linscale 0.5 narrows the 0-1 linear segment to half a decade so the
    # log decades carry the visual weight.
    panels = [
        (axes[0, 0], "lwp_g", "LWP (g m$^{-2}$)", "a",
         ("symlog", {"linthresh": 1.0, "linscale": 0.5})),
        (axes[0, 1], "iwp_g", "IWP (g m$^{-2}$)", "b",
         ("symlog", {"linthresh": 1.0, "linscale": 0.5})),
        (axes[0, 2], "tau", r"cloud optical depth $\tau$ (10/25 $\mu$m)",
         "c", ("log", {})),
        (axes[1, 0], "z_top_low",
         "lowest-layer top height (km)", "d", None),
        (axes[1, 1], "dz_low",
         "lowest-layer thickness (km)", "e", None),
    ]
    for ax, xcol, xlabel, _letter, xscale in panels:
        for cname in CLASSES:                   # ice first, ocean on top
            d = prof[prof["class"] == cname]
            ax.scatter(d[xcol], d[f"{q}_libradtran"], s=3.5,
                       color=CLASS_COLOURS[cname], lw=0, alpha=0.5,
                       label=cname)
        if xscale:
            ax.set_xscale(xscale[0], **xscale[1])
            if xscale[0] == "symlog":
                ax.set_xlim(left=0.0)     # no padding below the zero column
            ax.grid(which="major", axis="x", color="0.90", lw=0.6)
            ax.set_axisbelow(True)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylab)
    axes[0, 0].legend(fontsize=6, frameon=False, loc="upper left")

    # (f) binned median CRE vs tau, per class (with IQR band)
    ax = axes[1, 2]
    edges = np.logspace(np.log10(0.2), np.log10(60), 13)
    for cname in CLASSES:
        d = prof[prof["class"] == cname]
        xs, ys, lo, hi = [], [], [], []
        for e0, e1 in zip(edges[:-1], edges[1:]):
            m = d[(d["tau"] > e0) & (d["tau"] <= e1)]
            if len(m) < 10:
                continue
            xs.append(float(m["tau"].median()))
            v = m[f"{q}_libradtran"]
            ys.append(v.median())
            lo.append(v.quantile(.25))
            hi.append(v.quantile(.75))
        ax.fill_between(xs, lo, hi, color=CLASS_COLOURS[cname], alpha=0.18,
                        lw=0)
        ax.plot(xs, ys, "-o", ms=3, lw=1.2, color=CLASS_COLOURS[cname],
                label=cname)
    ax.set_xscale("log")
    ax.set_xlabel(r"cloud optical depth $\tau$ (10/25 $\mu$m)")
    ax.set_ylabel(f"median {qlab} (W m$^{{-2}}$)")
    ax.legend(fontsize=6, frameon=False, loc="upper left")
    panel_label(ax, "f", outside=True)

    fig.suptitle(f"MERRA-2 {manifest['area_tag']}  {tag}  overcast "
                 f"polar-night columns — surface LW {qlab} vs cloud "
                 f"properties", fontsize=8, y=0.995)
    # rect reserves the right 4% of the canvas for the last column's
    # y-marginal (inset axes are invisible to tight_layout and would be
    # clipped by the figure edge otherwise)
    fig.tight_layout(w_pad=4.0, h_pad=3.4, rect=(0, 0, 0.96, 0.95))

    # marginal per-class histograms on the scatter panels (a)-(e); added
    # after tight_layout so the main-axes geometry is final. Panel labels
    # move to the top marginal so they stay above everything. Bins are built
    # from the FINAL axis limits (so bars end exactly at the frame), and bars
    # show the fraction of the class per bin — mass, not density — so heights
    # stay comparable across the mixed linear/log bin widths.
    for ax, xcol, _xlabel, letter, xscale in panels:
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        if xscale and xscale[0] == "symlog":
            xbins = np.concatenate([np.linspace(0.0, 1.0, 5)[:-1],
                                    np.logspace(0.0, np.log10(xlim[1]), 21)])
        elif xscale and xscale[0] == "log":
            xbins = np.logspace(*np.log10(xlim), 25)
        else:
            xbins = np.linspace(0.0, xlim[1], 33)
        axt = ax.inset_axes([0.0, 1.02, 1.0, 0.15], sharex=ax)
        axr = ax.inset_axes([1.02, 0.0, 0.15, 1.0], sharey=ax)
        for cname in CLASSES:
            d = prof[prof["class"] == cname].dropna(subset=[xcol])
            w = np.full(len(d), 1.0 / len(d))
            axt.hist(d[xcol], bins=xbins, weights=w,
                     histtype="stepfilled", alpha=0.45,
                     color=CLASS_COLOURS[cname])
            axr.hist(d[f"{q}_libradtran"], bins=ybins, weights=w,
                     orientation="horizontal", histtype="stepfilled",
                     alpha=0.45, color=CLASS_COLOURS[cname])
        for a in (axt, axr):
            a.tick_params(which="both", left=False, bottom=False,
                          labelleft=False, labelbottom=False)
            for s in a.spines.values():
                s.set_visible(False)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        panel_label(axt, letter, outside=True)
    fdir = figures_dir(cfg)
    fdir.mkdir(parents=True, exist_ok=True)
    suffix = "" if args.flux == "dn" else f"_{args.flux}"
    fpath = fdir / f"cre_cloud{suffix}{thr_tag(args)}_{tag}.png"
    fig.savefig(fpath, dpi=300)
    print(f"wrote {fpath}")

    # ---------------- layering diagnostics figure
    fig2, (axg, axh, axi) = plt.subplots(1, 3, figsize=(9.8, 3.4))
    nmax, width = 4, 0.38
    for j, cname in enumerate(CLASSES):
        d = prof[(prof["class"] == cname) & (prof["n_sig"] > 0)]
        frac = [float((d["n_sig"] == k).mean()) for k in range(1, nmax)]
        frac.append(float((d["n_sig"] >= nmax).mean()))
        axg.bar(np.arange(1, nmax + 1) + (j - 0.5) * width, frac, width=width,
                color=CLASS_COLOURS[cname], label=cname)
    axg.set_xticks(range(1, nmax + 1),
                   [str(k) for k in range(1, nmax)] + [f"{nmax}+"])
    axg.set_xlabel("significant cloud layers per column")
    axg.set_ylabel("fraction of columns")
    axg.legend(fontsize=6, frameon=False)
    panel_label(axg, "a", outside=True)

    edges = np.logspace(np.log10(0.2), np.log10(60), 13)
    for cname in CLASSES:
        for sel, ls, sub in ((prof["n_sig"] == 1, "-", "single"),
                             (prof["n_sig"] >= 2, "--", "multi")):
            d = prof[(prof["class"] == cname) & sel]
            xs, ys = [], []
            for e0, e1 in zip(edges[:-1], edges[1:]):
                m = d[(d["tau"] > e0) & (d["tau"] <= e1)]
                if len(m) < 10:
                    continue
                xs.append(float(m["tau"].median()))
                ys.append(float(m[f"{q}_libradtran"].median()))
            axh.plot(xs, ys, ls, marker="o", ms=2.5, lw=1.1,
                     color=CLASS_COLOURS[cname],
                     label=f"{cname} {sub} (n={len(d)})")
    axh.set_xscale("log")
    axh.grid(which="major", axis="x", color="0.90", lw=0.6)
    axh.set_axisbelow(True)
    axh.set_xlabel(r"cloud optical depth $\tau$ (10/25 $\mu$m)")
    axh.set_ylabel(f"median {qlab} (W m$^{{-2}}$)")
    axh.legend(fontsize=5.5, frameon=False, loc="lower right")
    panel_label(axh, "b", outside=True)

    for cname in CLASSES:
        d = prof[prof["class"] == cname]
        axi.scatter(d["z_top_low"], d["z_top_col"], s=3.5,
                    color=CLASS_COLOURS[cname], lw=0, alpha=0.5, label=cname)
    top = float(prof["z_top_col"].max()) + 0.5
    axi.plot([0, top], [0, top], color="0.6", lw=0.7)
    axi.set_xlabel("lowest significant layer top (km)")
    axi.set_ylabel("whole-column cloud top (km)")
    axi.legend(fontsize=6, frameon=False, loc="lower right")
    panel_label(axi, "c", outside=True)

    fig2.suptitle(f"MERRA-2 {manifest['area_tag']}  {tag}  cloud layering "
                  f"(split: any clear level or 3x content dip; "
                  f"liquid-anchored, virga excluded; significant: "
                  f"$\\tau \\geq$ {args.tau_sig:g})", fontsize=8, y=0.99)
    fig2.tight_layout(rect=(0, 0, 1, 0.90))
    fpath2 = fdir / f"cre_layers{suffix}{thr_tag(args)}_{tag}.png"
    fig2.savefig(fpath2, dpi=300)
    print(f"wrote {fpath2}")
    return 0


def thr_tag(args) -> str:
    """Filename tag for non-default cloud thresholds (sensitivity runs)."""
    t = ""
    if args.content_min != CLOUD_CONTENT_MIN_GM3:
        t += f"_cmin{args.content_min:g}"
    if args.tau_sig != TAU_SIG_MIN:
        t += f"_tau{args.tau_sig:g}"
    return t


# ------------------------------------------------------------ analyze-drivers

def cloud_base_temperature(prof):
    """T at the lowest significant layer's base, from each atm profile file."""
    t_cb = []
    for atm, zb in zip(prof["atm_file"], prof["z_base_low"]):
        if not np.isfinite(zb):
            t_cb.append(np.nan)
            continue
        rows = np.loadtxt(atm, comments="#")          # TOA-first: z descending
        t_cb.append(float(np.interp(zb, rows[::-1, 0], rows[::-1, 2])))
    return np.array(t_cb)


def cmd_analyze_drivers(args) -> int:
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize

    from reanlib.config import figures_dir
    from reanlib.plotstyle import apply_agu_style, panel_label

    cfg = load_config(args.config)
    manifest, prof = load_merged(cfg, args, keep_files=True)
    prof["tau"] = (TAU_PER_G_LIQ * prof["lwp_g"]
                   + TAU_PER_G_ICE * prof["iwp_g"])
    prof = prof.join(layer_stats(prof, args.content_min, args.tau_sig))
    prof["t_cb"] = cloud_base_temperature(prof)
    prof["dt_cb"] = prof["t_cb"] - prof["skt"]
    tag = season_tag(args)
    q = f"cre_{args.flux}"
    qlab = {"dn": "CRE$_{dn}$", "net": "CRE$_{net}$"}[args.flux]
    d_ok = prof.dropna(subset=["z_top_low"])
    print(f"{len(d_ok)} columns with a qualifying cloud layer; dT_cb "
          f"median ice "
          f"{d_ok[d_ok['class'] == 'ice100']['dt_cb'].median():+.1f} K, "
          f"ocean {d_ok[d_ok['class'] == 'ocean100']['dt_cb'].median():+.1f} K")

    apply_agu_style()
    fdir = figures_dir(cfg)
    fdir.mkdir(parents=True, exist_ok=True)
    suffix = "" if args.flux == "dn" else f"_{args.flux}"

    # ---------------- figure 1: (top height x thickness) hexbin heat maps,
    # once for the LOWEST significant layer (CRE-dominant; validation:
    # median 84% of full CRE standalone) and once for the TOPMOST (the layer
    # a satellite IR/lidar view sees first). Single-layer columns appear
    # identically in both.
    norms = {"n": None, "tau": LogNorm(0.05, 60.0),
             "cre": Normalize(0.0, 85.0)}
    cmaps = {"n": "viridis", "tau": "magma", "cre": "plasma"}
    titles = {"n": "columns per hexbin", "tau": r"median $\tau$",
              "cre": f"median {qlab} (W m$^{{-2}}$)"}
    for lay_tag, xcol, ycol, which in (
            ("", "z_top_low", "dz_low", "lowest"),
            ("_top", "z_top_col", "dz_top", "topmost")):
        x_max = float(np.ceil(d_ok[xcol].max()) + 0.5)
        y_max = float(np.ceil(d_ok[ycol].max()) + 0.5)
        extent = (0.0, x_max, 0.0, y_max)
        gridsize = (22, 13)                  # hexagons ~0.5 km across
        fig, axes = plt.subplots(2, 3, figsize=(9.8, 6.4), sharex=True,
                                 sharey=True)
        letters = iter("abcdef")
        hbs = {}
        for i, cname in enumerate(CLASSES):
            d = d_ok[d_ok["class"] == cname]
            for j, key in enumerate(("n", "tau", "cre")):
                ax = axes[i, j]
                kw = dict(gridsize=gridsize, extent=extent, cmap=cmaps[key],
                          norm=norms[key], linewidths=0.2)
                if key == "n":
                    hb = ax.hexbin(d[xcol], d[ycol], mincnt=1, **kw)
                else:
                    vals = d["tau"] if key == "tau" else d[f"{q}_libradtran"]
                    # mincnt for C-hexbins counts members ABOVE the
                    # threshold, so 2 keeps hexagons with >= 3 columns
                    hb = ax.hexbin(d[xcol], d[ycol], C=vals,
                                   reduce_C_function=np.median, mincnt=2,
                                   **kw)
                hbs[cname, key] = hb
                lim = min(x_max, y_max)
                ax.plot([0, lim], [0, lim], color="0.5", lw=0.7, ls="--")
                if i == 0:
                    ax.set_title(titles[key], fontsize=8)
                if i == 1:
                    ax.set_xlabel(f"{which}-layer top height (km)")
                if j == 0:
                    ax.set_ylabel(f"{cname}\nthickness (km)")
                panel_label(ax, next(letters), outside=True)
        # shared count norm across both classes, one colorbar per column
        n_max = max(float(hbs[c, "n"].get_array().max()) for c in CLASSES)
        for cname in CLASSES:
            hbs[cname, "n"].set_norm(LogNorm(1, n_max))
        for j, key in enumerate(("n", "tau", "cre")):
            fig.colorbar(hbs[CLASSES[1], key], ax=[axes[0, j], axes[1, j]],
                         shrink=0.85, pad=0.03)
        fig.suptitle(f"MERRA-2 {manifest['area_tag']}  {tag}  {which}-layer "
                     f"geometry (dashed: thickness = top height)",
                     fontsize=8, y=0.98)
        fpath = fdir / f"cre_geom{lay_tag}{suffix}{thr_tag(args)}_{tag}.png"
        fig.savefig(fpath, dpi=300)
        print(f"wrote {fpath}")

    # ---------------- figure 2: cloud-base - surface temperature driver
    fig2, (axa, axb, axc) = plt.subplots(1, 3, figsize=(9.8, 3.4))
    dt_edges = np.arange(-16.0, 17.0, 2.0)

    def driver_panel(ax, dd, note):
        for cname in CLASSES:
            d = dd[dd["class"] == cname]
            ax.scatter(d["dt_cb"], d[f"{q}_libradtran"], s=3, lw=0,
                       alpha=0.2, color=CLASS_COLOURS[cname])
            xs, ys, lo, hi = [], [], [], []
            for e0, e1 in zip(dt_edges[:-1], dt_edges[1:]):
                m = d[(d["dt_cb"] > e0) & (d["dt_cb"] <= e1)]
                if len(m) < 10:
                    continue
                xs.append(float(m["dt_cb"].median()))
                v = m[f"{q}_libradtran"]
                ys.append(v.median())
                lo.append(v.quantile(.25))
                hi.append(v.quantile(.75))
            ax.fill_between(xs, lo, hi, color=CLASS_COLOURS[cname],
                            alpha=0.18, lw=0)
            ok = np.isfinite(d["dt_cb"]) & np.isfinite(d[f"{q}_libradtran"])
            r = np.corrcoef(d["dt_cb"][ok], d[f"{q}_libradtran"][ok])[0, 1]
            ax.plot(xs, ys, "-o", ms=3, lw=1.2, color=CLASS_COLOURS[cname],
                    label=f"{cname} (r = {r:+.2f})")
        ax.axvline(0, color="0.75", lw=0.6, zorder=0)
        ax.set_xlabel("T(cloud base) $-$ T(skin) (K)")
        ax.set_ylabel(f"surface {qlab} (W m$^{{-2}}$)")
        ax.legend(fontsize=6, frameon=False, loc="lower right")
        ax.text(0.03, 0.97, note, transform=ax.transAxes, fontsize=6,
                va="top")

    driver_panel(axa, d_ok, "all columns")
    driver_panel(axb, d_ok[d_ok["tau"] >= 5.0],
                 r"$\tau \geq 5$ (emissivity-saturated)")
    panel_label(axa, "a", outside=True)
    panel_label(axb, "b", outside=True)

    sc = None
    for cname, mk in (("ice100", "o"), ("ocean100", "^")):
        d = d_ok[d_ok["class"] == cname]
        sc = axc.scatter(d["dt_850_skt"], d["dt_cb"], s=4, marker=mk,
                         c=d[f"{q}_libradtran"], cmap="plasma", vmin=0,
                         vmax=85, lw=0, alpha=0.7)
        ok = np.isfinite(d["dt_cb"])
        r = np.corrcoef(d["dt_850_skt"][ok], d["dt_cb"][ok])[0, 1]
        print(f"  {cname}: r(dt_850_skt, dt_cb) = {r:+.2f}")
    lim = [-18, 18]
    axc.plot(lim, lim, color="0.6", lw=0.7)
    axc.set_xlim(lim)
    axc.set_xlabel("T(850 hPa) $-$ T(skin) (K)")
    axc.set_ylabel("T(cloud base) $-$ T(skin) (K)")
    cb = plt.colorbar(sc, ax=axc, pad=0.02)
    cb.set_label(f"{qlab} (W m$^{{-2}}$)", fontsize=7)
    axc.text(0.03, 0.97, "circles ice100, triangles ocean100",
             transform=axc.transAxes, fontsize=6, va="top")
    panel_label(axc, "c", outside=True)

    fig2.suptitle(f"MERRA-2 {manifest['area_tag']}  {tag}  cloud-base vs "
                  f"surface temperature as the {qlab} driver", fontsize=8,
                  y=0.99)
    fig2.tight_layout(rect=(0, 0, 1, 0.92))
    fpath2 = fdir / f"cre_drivers{suffix}{thr_tag(args)}_{tag}.png"
    fig2.savefig(fpath2, dpi=300)
    print(f"wrote {fpath2}")
    return 0


# ------------------------------------------------------------- census-layers

def cmd_census_layers(args) -> int:
    """AREA-TRUE layer statistics over the FULL eligible population.

    The analysis sample is stratified by dt_850_skt and therefore
    over-represents rare stability regimes (deep frontal cloud over pack
    ice); this census segments EVERY eligible overcast polar-night column of
    every snapshot (no RT), weights by cos(lat), and renders area-weighted
    geometry/thickness figures comparable to radar climatologies
    (Nomokonova et al. 2019 — see the reference block at the constants).
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.colors import LogNorm, Normalize

    from reanlib.classes import load_class_inputs
    from reanlib.config import figures_dir
    from reanlib.inversion import column_heights
    from reanlib.plotstyle import apply_agu_style, panel_label

    R_DRY = 287.04
    cfg = load_config(args.config)
    tag = season_tag(args)
    out_root = cre_dir(cfg)
    out_root.mkdir(parents=True, exist_ok=True)
    cpath = out_root / f"layer_census{thr_tag(args)}_{tag}.csv.gz"

    if cpath.exists() and not args.overwrite:
        df = pd.read_csv(cpath)
        print(f"loaded {len(df)} columns from {cpath}")
    else:
        rows = []
        for date in date_range(args.start, args.end):
            pp = plev_path(cfg, date)
            if not pp.exists():
                continue
            plev = open_era5(pp)
            sfc = open_era5(sfc_path(cfg, date))
            times = np.array([np.datetime64(f"{date}T{h:02d}:00")
                              for h in HOURS])
            lsm, sic_a, tcc_a = load_class_inputs(cfg, date, times)
            lat = plev["latitude"].values
            lat2 = np.broadcast_to(lat[:, None], lsm.shape[1:])
            decl = solar_declination_deg(dt.datetime(date.year, date.month,
                                                     date.day, 12))
            night = (np.sin(np.deg2rad(lat2)) * np.sin(np.deg2rad(decl))
                     + np.cos(np.deg2rad(lat2))
                     * np.cos(np.deg2rad(decl))) <= 0.0
            for k, when in enumerate(times):
                pl = plev.sel(valid_time=when).transpose(
                    "pressure_level", "latitude", "longitude")
                sf = sfc.sel(valid_time=when)
                p = pl["pressure_level"].values.astype(float)
                clwc = np.nan_to_num(pl["clwc"].values)
                ciwc = np.nan_to_num(pl["ciwc"].values)
                lwp = np.abs(TRAPZ(clwc[::-1], x=p[::-1] * 100.0,
                                   axis=0)) / G0 * 1e3
                iwp = np.abs(TRAPZ(ciwc[::-1], x=p[::-1] * 100.0,
                                   axis=0)) / G0 * 1e3
                sic, tcc = sic_a[k], tcc_a[k]
                sea = (lsm[k] < 0.5) & np.isfinite(sic)
                base = (sea & night & (tcc >= TCC_OVERCAST_MIN)
                        & (lwp + iwp > TWP_MIN_G))
                t3 = pl["t"].values
                q3 = pl["q"].values
                sp2 = sf["sp"].values / 100.0
                t2m2 = sf["t2m"].values
                for cname, cmask in (
                        ("ice100", base & (sic >= SIC_ICE100_MIN)),
                        ("ocean100", base & (sic <= SIC_OCEAN100_MAX))):
                    for iy, ix in zip(*np.where(cmask)):
                        above = (p <= sp2[iy, ix]) & np.isfinite(t3[:, iy, ix])
                        if above.sum() < 20:
                            continue
                        pc = p[above]
                        tc = t3[above, iy, ix]
                        z = column_heights(tc, pc, float(t2m2[iy, ix]),
                                           float(sp2[iy, ix]),
                                           q3[above, iy, ix]) / 1000.0
                        keep = z >= 0.01
                        rho = pc[keep] * 100.0 / (R_DRY * tc[keep])
                        wc = np.maximum(clwc[above, iy, ix][keep], 0) \
                            * rho * 1000.0
                        ic = np.maximum(ciwc[above, iy, ix][keep], 0) \
                            * rho * 1000.0
                        wc[wc < 1e-5] = 0.0        # cloud-file write floor
                        ic[ic < 1e-5] = 0.0
                        if not (wc + ic).any():
                            continue
                        sig = significant(
                            segment_layers(z[keep], wc, ic,
                                           args.content_min), args.tau_sig)
                        if not sig:
                            continue
                        low = sig[0]
                        rows.append({
                            "date": f"{date:%Y-%m-%d}", "hour": HOURS[k],
                            "class": cname, "lat": float(lat[iy]),
                            "n_sig": len(sig),
                            "dz_low": low["dz"], "z_top_low": low["z_top"],
                            "phase_low": low["phase"],
                            "tau_low": low["tau"],
                            "virga_dz_low": low["virga_dz"],
                            "dz_top": sig[-1]["dz"],
                            "z_top_col": sig[-1]["z_top"],
                            "lwp_g": float(lwp[iy, ix]),
                            "iwp_g": float(iwp[iy, ix]),
                        })
            plev.close()
            sfc.close()
            if date.day in (1, 15):
                print(f"{date}: {len(rows)} columns so far", flush=True)
        df = pd.DataFrame(rows)
        df.to_csv(cpath, index=False)
        print(f"wrote {cpath}  ({len(df)} columns)")

    df["w"] = np.cos(np.deg2rad(df["lat"]))
    df["tau"] = TAU_PER_G_LIQ * df["lwp_g"] + TAU_PER_G_ICE * df["iwp_g"]
    df["month"] = df["date"].str[:7]

    def wq(v, w, q):
        o = np.argsort(v.values)
        v2, w2 = v.values[o], w.values[o]
        return float(np.interp(q, np.cumsum(w2) / w2.sum(), v2))

    print("\nAREA-TRUE lowest-layer statistics (cos-lat weighted):")
    for cname in CLASSES:
        d = df[df["class"] == cname]
        f3 = float(d["w"][d["dz_low"] > 3.0].sum() / d["w"].sum())
        fm = float(d["w"][d["n_sig"] >= 2].sum() / d["w"].sum())
        print(f"  {cname:9s} n={len(d):7d}: thickness p50/p90/p99 = "
              f"{wq(d['dz_low'], d['w'], .5):.1f}/"
              f"{wq(d['dz_low'], d['w'], .9):.1f}/"
              f"{wq(d['dz_low'], d['w'], .99):.1f} km | >3 km {f3:.1%} | "
              f"multi-layer {fm:.1%}")
        for mo, dm in d.groupby("month"):
            print(f"      {mo}: p50/p90 = {wq(dm['dz_low'], dm['w'], .5):.1f}"
                  f"/{wq(dm['dz_low'], dm['w'], .9):.1f} km")
    single = df[df["n_sig"] == 1]
    print("  single-layer thickness, area-weighted model p50/p99 vs "
          "Cloudnet Nov-Jan (km):")
    for ph, (w50, w99) in NYA_OBS_WINTER.items():
        d = single[single["phase_low"] == ph]
        if len(d):
            print(f"    {ph:6s} (n={len(d):6d}): "
                  f"{wq(d['dz_low'], d['w'], .5):.2f}/"
                  f"{wq(d['dz_low'], d['w'], .99):.2f} vs "
                  f"{w50:.2f}/{w99:.1f}")

    # ---------------- figures
    apply_agu_style()
    fdir = figures_dir(cfg)
    fdir.mkdir(parents=True, exist_ok=True)

    # A: area-true geometry hexbins (lowest layer): area fraction, median
    # column tau, multi-layer fraction
    x_max = float(np.ceil(df["z_top_low"].max()) + 0.5)
    y_max = float(np.ceil(df["dz_low"].max()) + 0.5)
    fig, axes = plt.subplots(2, 3, figsize=(9.8, 6.4), sharex=True,
                             sharey=True)
    letters = iter("abcdef")
    specs = [("area fraction per hexbin", "viridis", None),
             (r"median column $\tau$", "magma", LogNorm(0.05, 60.0)),
             ("multi-layer fraction", "cividis", Normalize(0, 1))]
    hbs = {}
    for i, cname in enumerate(CLASSES):
        d = df[df["class"] == cname]
        wfrac = d["w"] / d["w"].sum()
        for j, (title, cmap, norm) in enumerate(specs):
            ax = axes[i, j]
            kw = dict(gridsize=(22, 13), extent=(0, x_max, 0, y_max),
                      cmap=cmap, linewidths=0.2)
            if j == 0:
                hb = ax.hexbin(d["z_top_low"], d["dz_low"], C=wfrac,
                               reduce_C_function=np.sum, mincnt=1, **kw)
            elif j == 1:
                hb = ax.hexbin(d["z_top_low"], d["dz_low"], C=d["tau"],
                               reduce_C_function=np.median, mincnt=2,
                               norm=norm, **kw)
            else:
                hb = ax.hexbin(d["z_top_low"], d["dz_low"],
                               C=(d["n_sig"] >= 2).astype(float),
                               reduce_C_function=np.mean, mincnt=2,
                               norm=norm, **kw)
            hbs[cname, j] = hb
            lim = min(x_max, y_max)
            ax.plot([0, lim], [0, lim], color="0.5", lw=0.7, ls="--")
            if i == 0:
                ax.set_title(title, fontsize=8)
            if i == 1:
                ax.set_xlabel("lowest-layer top height (km)")
            if j == 0:
                ax.set_ylabel(f"{cname}\nthickness (km)")
            panel_label(ax, next(letters), outside=True)
    a_max = max(float(hbs[c, 0].get_array().max()) for c in CLASSES)
    for cname in CLASSES:
        hbs[cname, 0].set_norm(LogNorm(1e-5, a_max))
    for j in range(3):
        fig.colorbar(hbs[CLASSES[1], j], ax=[axes[0, j], axes[1, j]],
                     shrink=0.85, pad=0.03)
    fig.suptitle(f"MERRA-2 {manifest_area(cfg)}  {tag}  AREA-TRUE "
                 f"lowest-layer geometry, full eligible population",
                 fontsize=8, y=0.98)
    fp_a = fdir / f"cre_geom_true{thr_tag(args)}_{tag}.png"
    fig.savefig(fp_a, dpi=300)
    print(f"wrote {fp_a}")

    # B: area-weighted thickness/top distributions, stratified sample
    # overlaid as outline for contrast
    try:
        _, sample = load_merged(cfg, args, keep_files=True)
        sample = sample.join(layer_stats(sample, args.content_min,
                                         args.tau_sig))
        sample = sample.dropna(subset=["dz_low"])
    except SystemExit:
        sample = None
    fig2, (axa, axb, axc) = plt.subplots(1, 3, figsize=(9.8, 3.4))
    for xcol, ax, bins, xlabel in (
            ("dz_low", axa, np.arange(0, 8.25, 0.25),
             "lowest-layer thickness (km)"),
            ("z_top_low", axb, np.arange(0, 13.5, 0.5),
             "lowest-layer top height (km)")):
        for cname in CLASSES:
            d = df[df["class"] == cname]
            ax.hist(d[xcol], bins=bins, weights=d["w"] / d["w"].sum(),
                    histtype="stepfilled", alpha=0.45,
                    color=CLASS_COLOURS[cname],
                    label=f"{cname} area-true")
            if sample is not None:
                s = sample[sample["class"] == cname]
                ax.hist(s[xcol], bins=bins,
                        weights=np.full(len(s), 1.0 / len(s)),
                        histtype="step", lw=1.0, ls="--",
                        color=CLASS_COLOURS[cname],
                        label=f"{cname} strat. sample")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("fraction of area / of sample")
    axa.legend(fontsize=5.5, frameon=False)
    panel_label(axa, "a", outside=True)
    panel_label(axb, "b", outside=True)

    # (c) phase-resolved single-layer thickness vs the season-matched
    # Nov-Jan Cloudnet statistics (nya_cloudnet_winter.py)
    xpos = np.arange(len(NYA_OBS_WINTER))
    for k, (ph, (w50, w99)) in enumerate(NYA_OBS_WINTER.items()):
        d = single[single["phase_low"] == ph]
        if not len(d):
            continue
        m50 = wq(d["dz_low"], d["w"], .5)
        m99 = wq(d["dz_low"], d["w"], .99)
        axc.bar(k - 0.18, m50, width=0.34, color="#0072B2",
                label="MERRA-2 area-true (Nov–Jan)" if k == 0 else None)
        axc.bar(k + 0.18, w50, width=0.34, color="0.35",
                label="Cloudnet Nov–Jan 2016/17+17/18" if k == 0 else None)
        axc.plot([k - 0.18], [m99], "v", ms=5, color="#0072B2")
        axc.plot([k + 0.18], [w99], "v", ms=5, color="0.2")
    axc.set_xticks(xpos, list(NYA_OBS_WINTER))
    axc.set_ylabel("single-layer thickness (km)")
    axc.text(0.02, 0.97, "bars: median, triangles: ~99th pct\n"
             "Ny-Ålesund AWIPEV (78.9°N, 11.9°E);\n"
             "obs exclude liquid-precipitating profiles",
             transform=axc.transAxes, fontsize=5.5, va="top")
    # legend below the frame (the panel interior stays clear for the bars)
    axc.legend(fontsize=5.5, frameon=False, loc="upper center",
               bbox_to_anchor=(0.5, -0.10), ncol=1)
    panel_label(axc, "c", outside=True)

    fig2.suptitle(f"MERRA-2 {manifest_area(cfg)}  {tag}  area-true cloud "
                  f"geometry vs stratified sample and winter Cloudnet "
                  f"observations at Ny-Ålesund", fontsize=8, y=0.99)
    # bottom margin reserves room for panel (c)'s below-frame legend
    fig2.tight_layout(rect=(0, 0.12, 1, 0.92))
    fpb = fdir / f"cre_thick_true{thr_tag(args)}_{tag}.png"
    fig2.savefig(fpb, dpi=300)
    print(f"wrote {fpb}")
    return 0


def manifest_area(cfg) -> str:
    return area_tag(cfg).lstrip("_") or "full-domain"


# ------------------------------------------------------------ validate-layers

def write_layer_file(path, z, val, idx, reff_um):
    """Cloud file holding ONE layer's content (zeros elsewhere), TOA-first."""
    out = np.zeros_like(val)
    out[idx] = val[idx]
    if not np.any(out > 0):
        return None
    rows = np.column_stack([z, out, np.full(z.size, reff_um)])[::-1]
    lines = ["{:9.3f} {:11.5e} {:7.2f}".format(*r) for r in rows]
    path.write_text("#   z(km)     wc(g/m^3)  reff(um)\n" + "\n".join(lines)
                    + "\n")
    return str(path)


def cmd_validate_layers(args) -> int:
    """Per-layer CRE attribution on a subsample of multi-layer columns.

    Each significant layer is simulated ALONE (other layers removed); its
    standalone CRE_dn against the column's existing clear-sky run measures
    how much of the full CRE each layer carries (er3t_env).
    """
    import copy
    import inspect

    import er3t
    import pandas as pd

    cfg = load_config(args.config)
    manifest, prof = load_merged(cfg, args, keep_files=True)
    prof = prof.join(layer_stats(prof, args.content_min, args.tau_sig))
    tag = season_tag(args)
    multi = prof[prof["n_sig"] >= 2]
    take = multi.sample(n=min(args.n, len(multi)), random_state=args.seed)
    print(f"validating {len(take)} of {len(multi)} multi-layer columns")

    tmp_dir = cre_dir(cfg) / "tmp_layers"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    if args.streams is None:
        args.streams = 4 if platform.system() == "Darwin" else 8
    lrt_cfg_base = er3t.rtm.lrt.get_lrt_cfg()
    lrt_cfg_base["mol_abs_param"] = "reptran " + args.mol_abs_param
    lrt_cfg_base["number_of_streams"] = args.streams
    mute_list = ["albedo", "wavelength", "spline", "source solar",
                 "slit_function_file"]
    accepted = inspect.signature(
        er3t.rtm.lrt.lrt_init_mono_flx.__init__).parameters

    inits, meta = [], []
    for _, row in take.iterrows():
        z, wc, ic = read_cloud_ladder(row["wc_file"], row["ic_file"])
        sig = significant(segment_layers(z, wc, ic, args.content_min),
                          args.tau_sig)
        if len(sig) < 2:
            continue
        when = dt.datetime.fromisoformat(row["date"]) \
            + dt.timedelta(hours=int(row["hour"]))
        for k, lay in enumerate(sig):
            stem = f"{row['label']}_L{k}"
            wcf = write_layer_file(tmp_dir / f"wc_{stem}.dat", z, wc,
                                   lay["idx"], REFF_LIQ_UM)
            icf = write_layer_file(tmp_dir / f"ic_{stem}.dat", z, ic,
                                   lay["idx"], REFF_ICE_UM)
            lrt_cfg = copy.deepcopy(lrt_cfg_base)
            lrt_cfg["atmosphere_file"] = row["atm_file"]
            extra = {
                "source": "thermal",
                "albedo_add": f"{1.0 - manifest['emissivity']:.4f}",
                "sur_temperature": f"{row['skt']:.2f}",
                "wavelength_add": "{:.0f} {:.0f}".format(
                    *manifest["wavelength_range_nm"]),
                "output_process": "integrate",
                "mol_file": "CH4 " + row["ch4_file"],
            }
            if wcf:
                extra["wc_file 1D"] = wcf
                extra["wc_properties"] = "hu interpolate"
            if icf:
                extra["ic_file 1D"] = icf
                extra["ic_properties"] = "fu interpolate"
            kwargs = dict(input_file=str(tmp_dir / f"in_{stem}.txt"),
                          output_file=str(tmp_dir / f"out_{stem}.txt"),
                          date=when, solar_zenith_angle=80.0, Nx=1,
                          output_altitude=[0, "toa"],
                          input_dict_extra=extra, mute_list=mute_list,
                          lrt_cfg=lrt_cfg, cld_cfg=None, aer_cfg=None)
            inits.append(er3t.rtm.lrt.lrt_init_mono_flx(
                **{kk: v for kk, v in kwargs.items() if kk in accepted}))
            meta.append((row["label"], k, len(sig), lay))

    workers = args.workers or max((os.cpu_count() or 2) - 2, 1)
    to_run = [ini for ini in inits
              if args.overwrite or not (os.path.exists(ini.output_file)
                                        and os.path.getsize(ini.output_file)
                                        > 50)]
    print(f"{len(to_run)} single-layer uvspec jobs "
          f"({len(inits)} layers total), {workers} workers")
    if to_run:
        er3t.rtm.lrt.lrt_run_mp(to_run, Ncpu=min(workers, len(to_run)))

    by_label = prof.set_index("label")
    rows = []
    for (label, k, n_sig, lay), ini in zip(meta, inits):
        data = er3t.rtm.lrt.lrt_read_uvspec_flx([ini])
        lwdn = float(np.squeeze(data.f_down)[0] * 1000.0)
        r = by_label.loc[label]
        cre_layer = lwdn - float(r["lwdn_clr_libradtran"])
        rows.append({
            "label": label, "class": r["class"], "layer": k, "n_sig": n_sig,
            "z_top": lay["z_top"], "tau": lay["tau"],
            "cre_dn_layer": cre_layer,
            "cre_dn_full": float(r["cre_dn_libradtran"]),
            "share": cre_layer / float(r["cre_dn_libradtran"]),
        })
    df = pd.DataFrame(rows)
    cpath = cre_dir(cfg) / f"cre_layer_validation_{tag}.csv"
    df.to_csv(cpath, index=False)
    print(f"wrote {cpath}")

    low = df[df["layer"] == 0]
    upper = df[df["layer"] > 0]
    tot = df.groupby("label")["share"].sum()
    print(f"lowest-layer standalone CRE / full CRE: "
          f"median {low['share'].median():.2f} "
          f"(IQR {low['share'].quantile(.25):.2f}-"
          f"{low['share'].quantile(.75):.2f})")
    print(f"upper layers standalone: median {upper['share'].median():.2f}")
    print(f"sum of standalone shares (overlap non-additivity): "
          f"median {tot.median():.2f}")
    return 0


# ---------------------------------------------------------------- CLI

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--config", default="config_cre.yaml")
        p.add_argument("--start", default="2019-11-01")
        p.add_argument("--end", default="2020-01-31")

    p = sub.add_parser("select", help="sample columns, write profiles+manifest")
    common(p)
    p.add_argument("--per-bin", type=int, default=2,
                   help="columns per (class, dt bin, snapshot) (default 2)")
    p.add_argument("--seed", type=int, default=42)
    p.set_defaults(fn=cmd_select)

    p = sub.add_parser("run-lrt", help="paired uvspec runs (er3t_env)")
    common(p)
    p.add_argument("--streams", type=int, default=None)
    p.add_argument("--mol-abs-param", default="coarse",
                   choices=["coarse", "medium", "fine"])
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--overwrite", action="store_true")
    p.set_defaults(fn=cmd_run_lrt)

    p = sub.add_parser("run-rrtmg", help="paired RRTMG-LW runs (era5 env)")
    common(p)
    p.add_argument("--overwrite", action="store_true")
    p.set_defaults(fn=cmd_run_rrtmg)

    p = sub.add_parser("analyze", help="stats + stratified figure")
    p.add_argument("--flux", default="dn", choices=["dn", "net"],
                   help="plot downwelling (dn) or net surface CRE "
                        "(net = eps*dn for the sims; LWGAB-LWGABCLR for "
                        "MERRA-2)")
    common(p)
    p.set_defaults(fn=cmd_analyze)

    p = sub.add_parser("analyze-cloud",
                       help="CRE vs cloud properties (LWP/IWP/tau/geometry)")
    p.add_argument("--flux", default="dn", choices=["dn", "net"])
    common(p)
    p.add_argument("--content-min", type=float,
                   default=CLOUD_CONTENT_MIN_GM3,
                   help="level cloud-content threshold (g/m3)")
    p.add_argument("--tau-sig", type=float, default=TAU_SIG_MIN,
                   help="layer significance threshold (tau)")
    p.set_defaults(fn=cmd_analyze_cloud)

    p = sub.add_parser("analyze-drivers",
                       help="geometry heat maps + cloud-base-temperature "
                            "driver analysis")
    p.add_argument("--flux", default="dn", choices=["dn", "net"])
    common(p)
    p.add_argument("--content-min", type=float,
                   default=CLOUD_CONTENT_MIN_GM3,
                   help="level cloud-content threshold (g/m3)")
    p.add_argument("--tau-sig", type=float, default=TAU_SIG_MIN,
                   help="layer significance threshold (tau)")
    p.set_defaults(fn=cmd_analyze_drivers)

    p = sub.add_parser("census-layers",
                       help="area-true layer census over the FULL eligible "
                            "population (no RT) + area-weighted figures")
    p.add_argument("--overwrite", action="store_true",
                   help="rebuild the census table even if cached")
    common(p)
    p.add_argument("--content-min", type=float,
                   default=CLOUD_CONTENT_MIN_GM3,
                   help="level cloud-content threshold (g/m3)")
    p.add_argument("--tau-sig", type=float, default=TAU_SIG_MIN,
                   help="layer significance threshold (tau)")
    p.set_defaults(fn=cmd_census_layers)

    p = sub.add_parser("validate-layers",
                       help="per-layer CRE attribution on a multi-layer "
                            "subsample (er3t_env)")
    p.add_argument("--n", type=int, default=100)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--streams", type=int, default=None)
    p.add_argument("--mol-abs-param", default="coarse",
                   choices=["coarse", "medium", "fine"])
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--overwrite", action="store_true")
    common(p)
    p.add_argument("--content-min", type=float,
                   default=CLOUD_CONTENT_MIN_GM3,
                   help="level cloud-content threshold (g/m3)")
    p.add_argument("--tau-sig", type=float, default=TAU_SIG_MIN,
                   help="layer significance threshold (tau)")
    p.set_defaults(fn=cmd_validate_layers)

    args = parser.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
