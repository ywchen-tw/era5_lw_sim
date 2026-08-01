#!/usr/bin/env python
"""LW flux simulation at every ERA5 column matched to a MOSAiC sounding.

Extends the two-pixel case study to the full month: each matched sounding's
(ERA5 pixel, 6-h time) column is simulated with libRadtran and RRTMG-LW —
once clear-sky and, wherever ERA5 holds condensate, once overcast — and an
all-sky flux is blended with the random-overlap effective cloud fraction
f = 1 - prod(1 - cc_k). Comparisons run against BOTH the MOSAiC surface
radiometers (Jozef et al. 2023) and ERA5's own flux product.

Duplicate handling: the simulation unit is the (pixel, valid_time) key.
Soundings sharing a key (possible with sub-6-h launch cadence) are simulated
once; their observed fluxes are averaged before statistics so one simulated
column is never counted twice.

Subcommands / envs:
  prep    (era5 env)     build per-column profile/cloud files + manifest
  run     (er3t_env)     uvspec thermal jobs (clear + overcast variants)
  figure  (era5 env)     RRTMG-LW inline + drift time series & scatters

Examples:
    conda activate era5     && python src/era5_mosaic_flux.py prep
    conda activate er3t_env && python src/era5_mosaic_flux.py run
    conda activate era5     && python src/era5_mosaic_flux.py figure
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from era5lib.config import REPO_ROOT, figures_dir, load_config, plev_path, sfc_path
from era5lib.io_era5 import open_era5
from era5_lrt_sim import (EMISSIVITY, REFF_ICE_UM, REFF_LIQ_UM, SIGMA,
                          WVL_RANGE_NM, build_profile_files,
                          load_cams_profiles, write_cloud_file)
from era5_case_study import TRAPZ, planck_band_fraction

G0 = 9.80665
MOSAIC_FILE = REPO_ROOT / "data" / "mosaic" / "MOSAiC_Atm_Properties.nc"


def mflux_dir(cfg: dict, year: int, month: int) -> Path:
    root = Path(cfg["paths"]["derived"])
    if not root.is_absolute():
        root = REPO_ROOT / root
    return root / f"{year:04d}" / f"{month:02d}" / "mosaic_flux"


def manifest_path(cfg, year, month) -> Path:
    return mflux_dir(cfg, year, month) / f"manifest_{year:04d}{month:02d}.json"


def results_path(cfg, year, month) -> Path:
    return mflux_dir(cfg, year, month) / f"results_{year:04d}{month:02d}.json"


# ---------------------------------------------------------------- prep

def cmd_prep(args) -> int:
    import pandas as pd
    import xarray as xr

    cfg = load_config(args.config)
    pairs = xr.open_dataset(
        REPO_ROOT / "derived" / f"{args.year}" / f"{args.month:02d}" /
        f"era5_mosaic_pairs_{args.year}{args.month:02d}.nc")
    mos = xr.open_dataset(MOSAIC_FILE)
    out_dir = mflux_dir(cfg, args.year, args.month)
    out_dir.mkdir(parents=True, exist_ok=True)
    cams = load_cams_profiles(args.year, args.month)
    if cams is None:
        sys.exit("no CAMS file — run src/era5_cams_download.py first")

    # group soundings by (ERA5 time, pixel): one simulation per key
    cache = {}

    def day_data(date):
        if date not in cache:
            cache[date] = (open_era5(plev_path(cfg, date)),
                           open_era5(sfc_path(cfg, date)))
        return cache[date]

    groups = {}
    for k in range(pairs.time.size):
        if not np.isfinite(pairs["match_dt_h"].values[k]):
            continue
        st = pd.Timestamp(pairs.time.values[k])
        et = (st + pd.Timedelta(hours=3)).floor("6h")
        date = et.date()
        if not plev_path(cfg, date).exists():
            continue
        plev, _ = day_data(date)
        iy = int(np.abs(plev["latitude"].values - float(pairs["lat"][k])).argmin())
        ix = int(np.abs(plev["longitude"].values - float(pairs["lon"][k])).argmin())
        groups.setdefault((et, iy, ix), []).append(st)

    print(f"{sum(len(v) for v in groups.values())} soundings -> "
          f"{len(groups)} unique (pixel, time) columns")

    columns = []
    for n, ((et, iy, ix), sts) in enumerate(sorted(groups.items())):
        plev, sfc = day_data(et.date())
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

        label = f"{et:%Y%m%dT%H}Z_{iy:02d}_{ix:03d}"
        atm_file, ch4_file, keep, z_keep = build_profile_files(
            out_dir, label, p[above], cols["t"], cols["q"], cols["o3"],
            t2m, sp_hpa, cams)
        wc_file = write_cloud_file(out_dir / f"wc_{label}.dat", z_keep,
                                   cols["clwc"][keep], p[above][keep],
                                   cols["t"][keep], REFF_LIQ_UM)
        ic_file = write_cloud_file(out_dir / f"ic_{label}.dat", z_keep,
                                   cols["ciwc"][keep], p[above][keep],
                                   cols["t"][keep], REFF_ICE_UM)

        p_pa_asc = p[above][::-1] * 100.0
        lwp = float(TRAPZ(cols["clwc"][::-1], x=p_pa_asc) / G0 * 1e3)
        iwp = float(TRAPZ(cols["ciwc"][::-1], x=p_pa_asc) / G0 * 1e3)
        cc_max = float(np.nanmax(cols["cc"]))
        f_eff = float(1.0 - np.prod(1.0 - np.clip(cols["cc"], 0.0, 1.0)))
        twp = lwp + iwp
        if cc_max <= 0.01 and twp <= 1.0:
            sky = "clear"
        elif cc_max >= 0.99 and twp > 1.0:
            sky = "overcast"
        else:
            sky = "partial"

        obs_dn = [float(mos["lwdn"].sel(time=st)) for st in sts]
        obs_up = [float(mos["lwup"].sel(time=st)) for st in sts]
        strd = float(s5["strd"].values[iy, ix]) / 3600.0
        str_net = float(s5["str"].values[iy, ix]) / 3600.0
        lat = float(p5["latitude"].values[iy])
        lon = float(p5["longitude"].values[ix])
        columns.append({
            "label": label, "etime": et.isoformat(),
            "sounding_times": [st.isoformat() for st in sts],
            "n_soundings": len(sts),
            "lat": lat, "lon": lon, "iy": iy, "ix": ix,
            "skt": skt, "t2m": t2m, "sp_hpa": sp_hpa,
            "era5_lwdn": strd, "era5_lwup": strd - str_net,
            "obs_lwdn": float(np.nanmean(obs_dn)),
            "obs_lwup": float(np.nanmean(obs_up)),
            "cc_max": cc_max, "f_eff": f_eff,
            "lwp_g": lwp, "iwp_g": iwp, "sky": sky,
            "mean_dt_h": float(np.mean([abs((et - st).total_seconds()) / 3600.0
                                        for st in sts])),
            "atm_file": atm_file, "ch4_file": ch4_file,
            "wc_file": wc_file, "ic_file": ic_file,
        })
        if (n + 1) % 25 == 0:
            print(f"  prepared {n + 1}/{len(groups)} columns ...")

    n_sky = {s: sum(1 for c in columns if c["sky"] == s)
             for s in ("clear", "partial", "overcast")}
    print(f"sky classes: {n_sky}")
    manifest = {"year": args.year, "month": args.month,
                "emissivity": EMISSIVITY,
                "wavelength_range_nm": list(WVL_RANGE_NM),
                "reff_um": {"liquid": REFF_LIQ_UM, "ice": REFF_ICE_UM},
                "gas_source": cams["source"],
                "columns": columns}
    mpath = manifest_path(cfg, args.year, args.month)
    mpath.write_text(json.dumps(manifest, indent=2))
    print(f"wrote {mpath}")
    return 0


# ---------------------------------------------------------------- run (er3t_env)

def cmd_run(args) -> int:
    import copy
    import er3t

    cfg = load_config(args.config)
    mpath = manifest_path(cfg, args.year, args.month)
    if not mpath.exists():
        sys.exit(f"missing {mpath} — run prep first")
    manifest = json.loads(mpath.read_text())
    tmp = mflux_dir(cfg, args.year, args.month) / "tmp"
    tmp.mkdir(exist_ok=True)
    streams = 4 if platform.system() == "Darwin" else 8

    lrt_cfg_base = er3t.rtm.lrt.get_lrt_cfg()
    lrt_cfg_base["mol_abs_param"] = "reptran coarse"
    lrt_cfg_base["number_of_streams"] = streams
    mute = ["albedo", "wavelength", "spline", "source solar",
            "slit_function_file"]

    inits, to_run = [], []
    for col in manifest["columns"]:
        variants = [("clear", None, None)]
        if col.get("wc_file") or col.get("ic_file"):
            variants.append(("cloudy", col.get("wc_file"), col.get("ic_file")))
        for vname, wc, ic in variants:
            label = f"{col['label']}_{vname}"
            out_path = str(tmp / f"output_{label}.txt")
            lrt_cfg = copy.deepcopy(lrt_cfg_base)
            lrt_cfg["atmosphere_file"] = col["atm_file"]
            extra = {
                "source": "thermal",
                "albedo_add": f"{1.0 - manifest['emissivity']:.4f}",
                "sur_temperature": f"{col['skt']:.2f}",
                "wavelength_add": "{:.0f} {:.0f}".format(
                    *manifest["wavelength_range_nm"]),
                "output_process": "integrate",
                "mol_file": "CH4 " + col["ch4_file"],
            }
            if wc:
                extra["wc_file 1D"] = wc
                extra["wc_properties"] = "hu interpolate"
            if ic:
                extra["ic_file 1D"] = ic
                extra["ic_properties"] = "fu interpolate"
            init = er3t.rtm.lrt.lrt_init_mono_flx(
                input_file=str(tmp / f"input_{label}.txt"),
                output_file=out_path,
                date=dt.datetime.fromisoformat(col["etime"]),
                solar_zenith_angle=80.0, Nx=1,
                output_altitude=[0, "toa"], input_dict_extra=extra,
                mute_list=mute, lrt_cfg=lrt_cfg, cld_cfg=None, aer_cfg=None)
            inits.append((col["label"], vname, init))
            import os
            if args.overwrite or not (os.path.exists(out_path)
                                      and os.path.getsize(out_path) > 50):
                to_run.append(init)

    if to_run:
        import os
        workers = args.workers or max((os.cpu_count() or 2) - 2, 1)
        print(f"running {len(to_run)} uvspec job(s) "
              f"({streams} streams, {workers} workers) ...")
        t0 = time.time()
        er3t.rtm.lrt.lrt_run_mp(to_run, Ncpu=min(workers, len(to_run)))
        print(f"uvspec done in {time.time() - t0:.0f} s")
    else:
        print("all uvspec outputs present (use --overwrite to rerun)")

    results = {}
    for label, vname, init in inits:
        data = er3t.rtm.lrt.lrt_read_uvspec_flx([init])
        f_dn = np.squeeze(data.f_down) * 1000.0     # undo er3t /1000 (thermal)
        f_up = np.squeeze(data.f_up) * 1000.0
        results.setdefault(label, {})[vname] = {
            "lwdn": float(f_dn[0]), "lwup": float(f_up[0]),
            "olr": float(f_up[1])}
    rpath = results_path(cfg, args.year, args.month)
    rpath.write_text(json.dumps(results, indent=2))
    print(f"wrote {rpath}")
    return 0


# ---------------------------------------------------------------- figure (era5)

SKY_COLOR = {"clear": "#0072B2", "partial": "#E69F00", "overcast": "#CC79A7"}


def blend(f_eff, clear, cloudy):
    if cloudy is None:
        return clear
    return (1.0 - f_eff) * clear + f_eff * cloudy


def cmd_figure(args) -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    import climlab
    from climlab.domain.axis import Axis
    import era5_rrtmg_sim as R

    from era5lib.plotstyle import apply_agu_style, panel_label

    cfg = load_config(args.config)
    manifest = json.loads(manifest_path(cfg, args.year, args.month).read_text())
    results = json.loads(results_path(cfg, args.year, args.month).read_text())
    R.patch_interface_temperature()
    w1, w2 = np.array(manifest["wavelength_range_nm"]) / 1000.0
    eps = manifest["emissivity"]

    rows = []
    t0 = time.time()
    for k, col in enumerate(manifest["columns"]):
        lib = results[col["label"]]
        # RRTMG: clear variant, plus cloudy variant when condensate exists
        clear_prof = dict(col, wc_file=None, ic_file=None)
        rc = R.run_pixel(clear_prof, manifest, climlab, Axis)
        rr = {"clear": {"lwdn": rc["sim_lwdn_sfc"], "lwup": rc["sim_lwup_sfc"]}}
        if col.get("wc_file") or col.get("ic_file"):
            r2 = R.run_pixel(col, manifest, climlab, Axis)
            rr["cloudy"] = {"lwdn": r2["sim_lwdn_sfc"],
                            "lwup": r2["sim_lwup_sfc"]}
        tail_up = (1 - planck_band_fraction(col["skt"], w1, w2)) * eps \
            * SIGMA * col["skt"] ** 4
        t_eff = (col["era5_lwdn"] / SIGMA) ** 0.25
        tail_dn = (1 - planck_band_fraction(t_eff, w1, w2)) * SIGMA * t_eff ** 4
        f = col["f_eff"]
        row = {
            "time": dt.datetime.fromisoformat(col["sounding_times"][0]),
            "sky": col["sky"], "f_eff": f,
            "obs_lwdn": col["obs_lwdn"], "obs_lwup": col["obs_lwup"],
            "era5_lwdn": col["era5_lwdn"], "era5_lwup": col["era5_lwup"],
            "lib_lwdn": blend(f, lib["clear"]["lwdn"],
                              lib.get("cloudy", {}).get("lwdn")) + tail_dn,
            "lib_lwup": blend(f, lib["clear"]["lwup"],
                              lib.get("cloudy", {}).get("lwup")) + tail_up,
            "rrt_lwdn": blend(f, rr["clear"]["lwdn"],
                              rr.get("cloudy", {}).get("lwdn")),
            "rrt_lwup": blend(f, rr["clear"]["lwup"],
                              rr.get("cloudy", {}).get("lwup")),
        }
        rows.append(row)
        if (k + 1) % 25 == 0:
            print(f"  RRTMG + blend {k + 1}/{len(manifest['columns'])} ...")
    import pandas as pd
    df = pd.DataFrame(rows).sort_values("time")
    print(f"assembled {len(df)} columns in {time.time() - t0:.0f} s")

    # ------------------------------------------------ statistics
    def stat_lines(kind):
        obs = df["obs_" + kind].values
        out = []
        for name, sim in (("libRadtran+tail", df["lib_" + kind].values),
                          ("RRTMG-LW", df["rrt_" + kind].values),
                          ("ERA5", df["era5_" + kind].values)):
            m = np.isfinite(obs) & np.isfinite(sim)
            r = np.corrcoef(obs[m], sim[m])[0, 1]
            bias = (sim[m] - obs[m]).mean()
            rmse = np.sqrt(((sim[m] - obs[m]) ** 2).mean())
            out.append((name, m.sum(), r, bias, rmse))
        return out

    print(f"\nvs MOSAiC radiometers ({args.year}-{args.month:02d}, "
          f"n columns = {len(df)}):")
    for kind in ("lwdn", "lwup"):
        for name, n, r, bias, rmse in stat_lines(kind):
            print(f"  {kind} {name:<15} (n={n:3d}): r={r:+.3f}  "
                  f"bias={bias:+.2f}  rmse={rmse:.2f} W/m2")

    # ------------------------------------------------ figure
    apply_agu_style()
    fig = plt.figure(figsize=(13.6, 10.6), layout="constrained")
    gs = fig.add_gridspec(3, 2, height_ratios=[0.75, 0.75, 1.1])

    series = (("obs", "#000000", "MOSAiC obs", "*", 6),
              ("era5", "#0072B2", "ERA5 (1-h mean)", "o", 3),
              ("lib", "#D55E00", "libRadtran + tail", "s", 3),
              ("rrt", "#009E73", "RRTMG-LW", "^", 3))
    for i, kind in enumerate(("lwdn", "lwup")):
        ax = fig.add_subplot(gs[i, :])
        for pre, color, lab, mk, ms in series:
            ax.plot(df["time"], df[pre + "_" + kind], mk + "-", color=color,
                    lw=0.9, ms=ms, label=lab if i == 0 else None, alpha=0.85)
        arrow = "\\downarrow" if kind == "lwdn" else "\\uparrow"
        ax.set_ylabel(f"LW${arrow}$ (W m$^{{-2}}$)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax.grid(color="#e3e3e3", lw=0.6)
        ax.set_axisbelow(True)
        if i == 0:
            ax.legend(frameon=False, ncol=4, fontsize=8, loc="upper left")
        panel_label(ax, "ab"[i], x=-0.05, y=1.06)

    for i, kind in enumerate(("lwdn", "lwup")):
        ax = fig.add_subplot(gs[2, i])
        obs = df["obs_" + kind].values
        vals = [obs[np.isfinite(obs)]]
        for pre, mk in (("lib", "o"), ("rrt", "^")):
            sim = df[pre + "_" + kind].values
            for sky, color in SKY_COLOR.items():
                m = (df["sky"] == sky).values & np.isfinite(obs) & np.isfinite(sim)
                ax.scatter(obs[m], sim[m], s=13, marker=mk, color=color,
                           alpha=0.65, edgecolors="none")
            vals.append(sim[np.isfinite(sim)])
        sim = df["era5_" + kind].values
        m = np.isfinite(obs) & np.isfinite(sim)
        ax.scatter(obs[m], sim[m], s=16, marker="+", color="#555555",
                   alpha=0.7, lw=0.9)
        vals.append(sim[m])
        lim = (min(v.min() for v in vals) - 4, max(v.max() for v in vals) + 4)
        ax.plot(lim, lim, color="#bbbbbb", lw=1, zorder=0)
        txt = "\n".join(f"{name}: r={r:+.2f}, bias={bias:+.1f}, rmse={rmse:.1f}"
                        for name, n, r, bias, rmse in stat_lines(kind))
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=7,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#cccccc"))
        arrow = "\\downarrow" if kind == "lwdn" else "\\uparrow"
        ax.set_xlabel(f"MOSAiC LW${arrow}$ (W m$^{{-2}}$)")
        ax.set_ylabel(f"simulated / ERA5 LW${arrow}$ (W m$^{{-2}}$)")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect("equal")
        ax.grid(color="#e3e3e3", lw=0.6)
        ax.set_axisbelow(True)
        panel_label(ax, "cd"[i], x=-0.14, y=1.04)
        if i == 1:
            from matplotlib.lines import Line2D
            handles = [Line2D([], [], ls="none", marker="o", color=c,
                              label=f"{s} (circle: libRadtran)")
                       for s, c in SKY_COLOR.items()]
            handles += [Line2D([], [], ls="none", marker="^", color="#888888",
                               label="triangle: RRTMG"),
                        Line2D([], [], ls="none", marker="+", color="#555555",
                               label="ERA5 product")]
            ax.legend(handles=handles, frameon=False, fontsize=6.4,
                      loc="lower right")

    fig.suptitle(
        f"Surface LW along the MOSAiC drift — {args.year}-{args.month:02d}: "
        f"simulations (ERA5 columns, all-sky blend by f = 1$-\\Pi$(1$-$cc)) "
        "vs MOSAiC radiometers and ERA5")
    out = figures_dir(cfg) / f"mosaic_flux_{args.year:04d}{args.month:02d}.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    for name, fn in (("prep", cmd_prep), ("run", cmd_run), ("figure", cmd_figure)):
        p = sub.add_parser(name)
        p.add_argument("--year", type=int, default=2020)
        p.add_argument("--month", type=int, default=1)
        p.add_argument("--config", default=None)
        if name == "run":
            p.add_argument("--workers", type=int, default=None)
            p.add_argument("--overwrite", action="store_true")
        p.set_defaults(func=fn)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
