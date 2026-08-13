#!/usr/bin/env python
"""Sounding-by-sounding comparison of reanalysis inversion metrics with MOSAiC.

Matches each MOSAiC radiosonde (Jozef et al. 2023, ESSD 15, 4983-5010;
doi:10.1594/PANGAEA.957760) launched in the chosen month against the nearest
reanalysis grid point and 6-hourly time step of the derived inversion files, and
compares surface-based inversion presence, strength, and depth, plus the 2 m
temperature.

Definitions are not identical and the comparison keeps them explicit:
- MOSAiC (obs): temperature inversions detected where dT/dz > 0.65 C/100 m
  over >= 25 m (5 m resolution). A sounding counts as having an SBI when its
  lowest inversion base (inv_alt) is <= --sbi-base-max (default 100 m); its
  strength/depth are that layer's inv_dt / inv_dz.
- Reanalysis: Kahl/Serreze-style scan on pressure levels, strength >= 0.5 K
  (see reanlib/inversion.py).

Outputs derived/<source>/YYYY/MM/<source>_mosaic_pairs_YYYYMM.nc and one figure.

Examples:
    python src/mosaic_compare.py --year 2020 --month 1
"""

from __future__ import annotations

import argparse
import calendar
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reanlib.config import (REPO_ROOT, SOURCE_LABELS, SOURCES, figures_dir,
                            inversion_path, load_config, pairs_path)
from reanlib.grid import GridIndex
from reanlib.io_era5 import open_era5
from reanlib.mapping import polar_panel
from reanlib.mosaic import sounding_profile
from reanlib.plotstyle import apply_agu_style, panel_label

MOSAIC_FILE = REPO_ROOT / "data" / "mosaic" / "MOSAiC_Atm_Properties.nc"
C_OBS = "#0072B2"
C_ERA = "#D55E00"



def load_mosaic(year: int, month: int, sbi_base_max: float) -> xr.Dataset:
    """One row per sounding in the month: obs SBI presence/strength/depth."""
    if not MOSAIC_FILE.exists():
        sys.exit(f"missing {MOSAIC_FILE}\nfetch with: curl -L -o {MOSAIC_FILE} "
                 "https://download.pangaea.de/dataset/957760/files/MOSAiC_Atm_Properties.nc")
    ds = xr.open_dataset(MOSAIC_FILE)
    t = ds["time"].values
    sel = np.array([(ts.astype("datetime64[s]").item().year == year
                     and ts.astype("datetime64[s]").item().month == month)
                    for ts in t])
    ds = ds.isel(time=sel)

    base = ds["inv_alt"].values           # (inversion_number, time), m
    dt_k = ds["inv_dt"].values            # K (degC differences)
    dz_m = ds["inv_dz"].values
    # lowest inversion layer of each sounding
    low = np.nanargmin(np.where(np.isfinite(base), base, np.inf), axis=0)
    idx = (low, np.arange(base.shape[1]))
    lowest_base = base[idx]
    obs_found = np.isfinite(lowest_base) & (lowest_base <= sbi_base_max)
    return xr.Dataset(
        {
            "lat": ds["lat"], "lon": ds["lon"],
            "obs_sbi_found": (("time",), obs_found),
            "obs_strength": (("time",), np.where(obs_found, dt_k[idx], np.nan)),
            "obs_depth": (("time",), np.where(obs_found, dz_m[idx], np.nan)),
            "obs_base": (("time",), lowest_base),
            "obs_t2m": ds["t_2m"] + 273.15,
        },
        coords={"time": ds["time"]},
    )


def obs_fixed_metrics(obs: xr.Dataset) -> xr.Dataset:
    """Observed counterparts of the reanalysis fixed-level metrics, per sounding.

    T(850 hPa) - T(2 m) uses the tower 2 m temperature and the level-2
    radiosonde profile interpolated in log-p; T(925) - T(1000) is NaN when
    the sounding's surface pressure is below 1000 hPa (level underground) —
    the same masking the reanalysis metric applies.
    """
    dt850, dt925 = [], []
    for ts, t2m in zip(obs["time"].values, obs["obs_t2m"].values):
        sel = sounding_profile(ts)
        if sel is None:
            dt850.append(np.nan)
            dt925.append(np.nan)
            continue
        p = sel["PPPP [hPa]"].values
        t = sel["TTT [°C]"].values + 273.15
        good = np.isfinite(p) & np.isfinite(t)
        order = np.argsort(p[good])
        p, t = p[good][order], t[good][order]

        def t_at(level):
            if p.size > 1 and p[0] <= level <= p[-1]:
                return float(np.interp(np.log(level), np.log(p), t))
            return np.nan

        dt850.append(t_at(850.0) - t2m if np.isfinite(t2m) else np.nan)
        dt925.append(t_at(925.0) - t_at(1000.0))
    return obs.assign(obs_dt_850_2m=(("time",), np.array(dt850)),
                      obs_dt_925_1000=(("time",), np.array(dt925)))


def match_reanalysis(cfg: dict, obs: xr.Dataset, max_dt_h: float = 3.0) -> xr.Dataset:
    """Nearest-in-space-and-time reanalysis values for each sounding."""
    fields = {k: [] for k in ("rean_sbi_found", "rean_strength", "rean_depth",
                              "rean_dt_850_2m", "rean_dt_925_1000",
                              "rean_t2m", "match_dt_h", "match_km")}
    cache: dict = {}
    index = None
    for ts, la, lo in zip(obs["time"].values, obs["lat"].values, obs["lon"].values):
        t_near = (ts + np.timedelta64(3, "h")).astype("datetime64[6h]").astype("datetime64[ns]")
        date = t_near.astype("datetime64[D]").item()
        if date not in cache:
            p = inversion_path(cfg, date)
            cache.clear()
            cache[date] = open_era5(p) if p.exists() else None
            # one KD-tree per grid, not per sounding: on CARRA-2's 2.5 km mesh
            # building it is the expensive part of the match
            index = GridIndex(cache[date]) if cache[date] is not None else None
        inv = cache[date]
        dt_h = abs((t_near - ts) / np.timedelta64(1, "h"))
        if inv is None or dt_h > max_dt_h or not np.isfinite(la):
            for k in fields:
                fields[k].append(np.nan)
            continue
        (iy, ix), km = index.query(la, lo)
        col = inv.sel(valid_time=t_near).isel(**index.isel(iy, ix))
        fields["rean_sbi_found"].append(float(col["sbi_found"]))
        fields["rean_strength"].append(float(col["sbi_strength"]))
        fields["rean_depth"].append(float(col["sbi_depth_z"]))
        fields["rean_dt_850_2m"].append(float(col["dt_850_2m"]))
        fields["rean_dt_925_1000"].append(float(col["dt_925_1000"]))
        fields["rean_t2m"].append(float(col["t2m"]))
        fields["match_dt_h"].append(dt_h)
        fields["match_km"].append(km)
    return obs.assign({k: (("time",), np.array(v, dtype=float))
                       for k, v in fields.items()})


def source_label(pairs: xr.Dataset) -> str:
    return pairs.attrs.get("source_label", "ERA5")


def stats_block(pairs: xr.Dataset) -> str:
    lab = source_label(pairs)
    ok = np.isfinite(pairs["match_dt_h"].values)
    obs_f = pairs["obs_sbi_found"].values.astype(bool) & ok
    era_f = (pairs["rean_sbi_found"].values == 1) & ok
    both = obs_f & era_f
    km = pairs["match_km"].values[ok]
    hrs = pairs["match_dt_h"].values[ok]
    lines = [f"matched soundings          : {int(ok.sum())} / {ok.size}",
             # max distance is the tell-tale for a sounding outside the
             # downloaded domain: it still matches, to the nearest edge cell
             f"match offset (median/max)  : {np.median(hrs):.1f} h / "
             f"{np.median(km):.0f} km  (max {np.max(hrs):.1f} h / "
             f"{np.max(km):.0f} km)",
             f"obs SBI frequency          : {obs_f.sum() / ok.sum():.1%}",
             f"{lab} SBI frequency".ljust(27) + f": {era_f.sum() / ok.sum():.1%}",
             f"detection agreement        : "
             f"{((obs_f == era_f) & ok).sum() / ok.sum():.1%}"]
    for name, o, e, unit in (
            ("strength", pairs["obs_strength"].values, pairs["rean_strength"].values, "K"),
            ("depth", pairs["obs_depth"].values, pairs["rean_depth"].values, "m")):
        m = both & np.isfinite(o) & np.isfinite(e)
        if m.sum() > 2:
            r = np.corrcoef(o[m], e[m])[0, 1]
            bias = (e[m] - o[m]).mean()
            rmse = np.sqrt(((e[m] - o[m]) ** 2).mean())
            lines.append(f"{name:<10} (both detect, n={m.sum():3d}): r={r:+.2f}  "
                         f"bias={bias:+.2f} {unit}  rmse={rmse:.2f} {unit}")
    # fixed-level metrics vs the SAME metric derived from the soundings
    for name, evar, ovar in (("dt_850_2m", "rean_dt_850_2m", "obs_dt_850_2m"),
                             ("dt_925_1000", "rean_dt_925_1000",
                              "obs_dt_925_1000")):
        if ovar not in pairs:
            lines.append(f"{name:<11}: no obs counterpart (level-2 soundings "
                         "unavailable; rerun with --overwrite)")
            continue
        o = pairs[ovar].values
        e = pairs[evar].values
        m = ok & np.isfinite(o) & np.isfinite(e)
        if m.sum() > 2:
            r = np.corrcoef(o[m], e[m])[0, 1]
            lines.append(f"{name:<11} vs obs {name} (n={m.sum():3d}): "
                         f"r={r:+.2f}  bias={(e[m] - o[m]).mean():+.2f} K  "
                         f"rmse={np.sqrt(((e[m] - o[m]) ** 2).mean()):.2f} K")
    o, e = pairs["obs_t2m"].values, pairs["rean_t2m"].values
    m = ok & np.isfinite(o) & np.isfinite(e)
    r = np.corrcoef(o[m], e[m])[0, 1]
    lines.append(f"t2m        (n={m.sum():3d})           : r={r:+.2f}  "
                 f"bias={(e[m] - o[m]).mean():+.2f} K  "
                 f"rmse={np.sqrt(((e[m] - o[m])**2).mean()):.2f} K")
    return "\n".join("  " + s for s in lines)


C_DT850 = "#009E73"
C_DT925 = "#CC79A7"


def fig_compare(pairs: xr.Dataset, year: int, month: int, outdir: Path) -> Path:
    import cartopy.crs as ccrs

    lab = source_label(pairs)
    fig = plt.figure(figsize=(13.6, 10.2), layout="constrained")
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.85])
    t = pairs["time"].values
    obs_s = np.where(pairs["obs_sbi_found"], pairs["obs_strength"], 0.0)
    era_s = np.nan_to_num(pairs["rean_strength"].values)

    # (a) drift track on the domain, colored by observed SBI strength
    ax = fig.add_subplot(gs[0, 0], projection=ccrs.NorthPolarStereo())
    ax.set_facecolor("#f2f2f2")
    ax.coastlines(resolution="50m", lw=0.7, color="#444444")
    ax.gridlines(draw_labels=False, lw=0.4, color="#bbbbbb",
                 xlocs=range(-180, 181, 45), ylocs=[80, 85])
    ax.set_extent([-180, 180, 80, 90], crs=ccrs.PlateCarree())
    import matplotlib.path as mpath
    theta = np.linspace(0, 2 * np.pi, 200)
    ax.set_boundary(mpath.Path(np.column_stack([np.sin(theta), np.cos(theta)]) * 0.5 + 0.5),
                    transform=ax.transAxes)
    sc = ax.scatter(pairs["lon"], pairs["lat"], c=obs_s, s=14, cmap="viridis",
                    vmin=0, vmax=max(np.nanmax(obs_s), 1.0),
                    transform=ccrs.PlateCarree(), zorder=5)
    cb = plt.colorbar(sc, ax=ax, orientation="horizontal", pad=0.05, fraction=0.045)
    cb.set_label("observed SBI strength (K)", fontsize=10)
    ax.set_title("Polarstern drift track")
    panel_label(ax, "a", x=-0.02, y=1.05)

    # (b) time series along the drift: obs SBI vs the three reanalysis metrics
    ax = fig.add_subplot(gs[0, 1:])
    ax.plot(t, obs_s, "o-", color=C_OBS, lw=1.4, ms=3.5, label="MOSAiC SBI")
    ax.plot(t, era_s, "s-", color=C_ERA, lw=1.2, ms=3.0, label=f"{lab} SBI")
    ax.plot(t, pairs["rean_dt_850_2m"], "-", color=C_DT850, lw=1.0,
            label=f"{lab} T850 − T2m")
    ax.plot(t, pairs["rean_dt_925_1000"], "-", color=C_DT925, lw=1.0,
            label=f"{lab} T925 − T1000")
    # observed counterparts of the fixed-level metrics: dashed, same colors
    if "obs_dt_850_2m" in pairs:
        ax.plot(t, pairs["obs_dt_850_2m"], "--", color=C_DT850, lw=1.0,
                label="MOSAiC T850 − T2m")
        ax.plot(t, pairs["obs_dt_925_1000"], "--", color=C_DT925, lw=1.0,
                label="MOSAiC T925 − T1000")
    ax.axhline(0, color="#bbbbbb", lw=0.8)
    ax.set_ylabel("inversion strength (K)  [SBI: 0 = none]")
    ax.legend(frameon=False, ncol=3, fontsize=8, loc="upper left")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.grid(color="#e3e3e3", lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    panel_label(ax, "b", x=-0.06, y=1.06)

    # (c-e) scatter: each reanalysis metric vs the SAME metric from the observations
    obs_found = pairs["obs_sbi_found"].values.astype(bool)
    nanvec = np.full(pairs["time"].shape, np.nan)
    allsky = np.ones(pairs["time"].shape, dtype=bool)
    panels = (
        (pairs["obs_strength"].values, pairs["rean_strength"].values,
         "MOSAiC SBI strength", f"{lab} SBI strength", C_OBS,
         obs_found & (pairs["rean_sbi_found"].values == 1)),
        (pairs["obs_dt_850_2m"].values if "obs_dt_850_2m" in pairs else nanvec,
         pairs["rean_dt_850_2m"].values,
         "MOSAiC T850 − T2m", f"{lab} T850 − T2m", C_DT850, allsky),
        (pairs["obs_dt_925_1000"].values if "obs_dt_925_1000" in pairs else nanvec,
         pairs["rean_dt_925_1000"].values,
         "MOSAiC T925 − T1000", f"{lab} T925 − T1000", C_DT925, allsky),
    )
    for i, (ox, ex, xlabel, ylabel, color, m) in enumerate(panels):
        ax = fig.add_subplot(gs[1, i])
        m = m & np.isfinite(ox) & np.isfinite(ex)
        panel_label(ax, "cde"[i], x=-0.16, y=1.06)
        if m.sum() < 3:
            ax.text(0.5, 0.5, "level-2 soundings\nunavailable", ha="center",
                    va="center", transform=ax.transAxes, color="#888888")
            ax.set_xlabel(f"{xlabel} (K)")
            ax.set_ylabel(f"{ylabel} (K)")
            continue
        lim = (min(np.nanmin(ox[m]), np.nanmin(ex[m]), 0.0) - 1,
               max(np.nanmax(ox[m]), np.nanmax(ex[m])) + 1)
        ax.plot(lim, lim, color="#bbbbbb", lw=1)
        ax.scatter(ox[m], ex[m], s=14, color=color, alpha=0.6, edgecolors="none")
        r = np.corrcoef(ox[m], ex[m])[0, 1]
        bias = (ex[m] - ox[m]).mean()
        ax.text(0.04, 0.96, f"n = {m.sum()}\nr = {r:+.2f}\nbias = {bias:+.2f} K",
                transform=ax.transAxes, va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#cccccc"))
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel(f"{xlabel} (K)")
        ax.set_ylabel(f"{ylabel} (K)")
        ax.set_aspect("equal")
        ax.grid(color="#e3e3e3", lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    fig.suptitle(f"{lab} vs. MOSAiC radiosondes — {calendar.month_name[month]} {year}\n"
                 "SBI: gradient criterion (Jozef et al. 2023) vs Kahl/Serreze "
                 "scan; fixed metrics from level-2 profiles (Maturilli et al. "
                 "2021) + tower 2 m", y=1.04)
    path = outdir / f"mosaic_compare_{year:04d}{month:02d}.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--month", type=int, required=True)
    parser.add_argument("--sbi-base-max", type=float, default=100.0,
                        help="max obs inversion-base altitude (m) to count as "
                             "surface-based (default 100)")
    parser.add_argument("--source", choices=list(SOURCES), default=None,
                        help="data source (default from config.yaml)")
    parser.add_argument("--config", default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args(argv)

    apply_agu_style()
    cfg = load_config(args.config, source=args.source)
    target = pairs_path(cfg, args.year, args.month)
    if target.exists() and not args.overwrite:
        print(f"{target} exists, loading it (use --overwrite to recompute)")
        pairs = xr.open_dataset(target)
    else:
        obs = load_mosaic(args.year, args.month, args.sbi_base_max)
        if obs.sizes["time"] == 0:
            sys.exit(f"no MOSAiC soundings in {args.year}-{args.month:02d}")
        print(f"{obs.sizes['time']} MOSAiC soundings in {args.year}-{args.month:02d}, "
              f"matching against {SOURCE_LABELS[cfg['source']]} ...")
        obs = obs_fixed_metrics(obs)
        pairs = match_reanalysis(cfg, obs)
        lab = SOURCE_LABELS[cfg["source"]]
        pairs.attrs = {
            "title": f"{lab} vs MOSAiC radiosonde inversion comparison",
            "source_label": lab,
            "obs_source": "Jozef et al. 2023, doi:10.1594/PANGAEA.957760 (CC-BY-4.0)",
            "obs_profile_source": ("level-2 radiosondes, Maturilli et al. 2021, "
                                   "doi:10.1594/PANGAEA.928656 (CC-BY-4.0)"),
            "obs_criterion": ("dT/dz > 0.65 C/100m over >= 25 m; SBI = lowest "
                              f"inversion base <= {args.sbi_base_max:g} m"),
            "rean_criterion": "Kahl/Serreze pressure-level scan, strength >= 0.5 K",
        }
        target.parent.mkdir(parents=True, exist_ok=True)
        pairs.to_netcdf(target)
        print(f"wrote {target}")

    print(stats_block(pairs))
    if not args.no_figures:
        outdir = figures_dir(cfg)
        outdir.mkdir(parents=True, exist_ok=True)
        print(f"wrote {fig_compare(pairs, args.year, args.month, outdir)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
