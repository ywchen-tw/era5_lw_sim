#!/usr/bin/env python
"""Monthly inversion-strength analysis from the daily derived files.

Aggregates derived/YYYY/MM/DD/era5_inversion_YYYYMMDD.nc over a month into
derived/YYYY/MM/era5_inversion_monthly_YYYYMM.nc and renders three figures:
monthly maps, strength distributions, and a domain-mean time series. All
domain statistics are area-weighted by cos(latitude).

Monthly per-grid-point statistics:
  sbi_frequency        fraction of time steps with a detected SBI
  sbi_strength_cond    mean strength over detected cases only ("intensity",
                       Zhang et al. 2011)
  sbi_strength_uncond  mean strength counting not-found as 0 K
  sbi_depth_p_mean, sbi_depth_z_mean, sbi_top_p_mean   (detected cases only)
  dt_850_2m_mean, dt_925_1000_mean   means over unmasked time steps

Examples:
    python src/monthly_stats.py --year 2020 --month 1
    python src/monthly_stats.py --year 2020 --month 1 --days 1-7 --skip-missing
"""

from __future__ import annotations

import argparse
import calendar
import datetime as dt
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from era5_download import parse_days
from reanlib.config import (SOURCES, figures_dir, inversion_path, load_config,
                            monthly_path, source_label)
from reanlib.grid import GRID_MAPPING_VAR, area_weights, domain_mask, grid_template
from reanlib.grid import hdims as grid_hdims
from reanlib.grid import horizontal_coords, weighted_mean
from reanlib.io_era5 import open_era5
from reanlib.mapping import grid_kwargs, polar_panel
from reanlib.plotstyle import apply_agu_style, panel_label

STRENGTH_BINS = np.arange(0.0, 25.5, 0.5)    # K
DEPTH_BINS = np.arange(0.0, 2050.0, 50.0)    # m


def aggregate(files: list[Path], label: str = "ERA5") -> xr.Dataset:
    """Stream the daily files, accumulating monthly statistics."""
    sums: dict[str, np.ndarray] = {}
    hist_strength = np.zeros(STRENGTH_BINS.size - 1)
    hist_sd = np.zeros((STRENGTH_BINS.size - 1, DEPTH_BINS.size - 1))
    ts_time, ts_found, ts_cond, ts_uncond = [], [], [], []
    template = w2d = valid = None
    n_time = 0

    for path in files:
        ds = open_era5(path)
        if template is None:
            # keep the grid description of the first file: coordinates, the CF
            # grid mapping (CARRA-2) and the domain mask all come from it
            template = grid_template(ds)
            w2d = area_weights(ds)
            valid = domain_mask(ds)
            shape = w2d.shape
            for key in ("found", "strength_found", "depth_p", "depth_z", "top_p",
                        "strength_uncond", "dt850", "dt925"):
                sums[key] = np.zeros(shape)
            for key in ("n_dt850", "n_dt925"):
                sums[key] = np.zeros(shape)

        found = ds["sbi_found"].values.astype(bool)
        strength = ds["sbi_strength"].values.astype(float)
        depth_z = ds["sbi_depth_z"].values.astype(float)
        n_time += found.shape[0]

        sums["found"] += found.sum(axis=0)
        sums["strength_found"] += np.nansum(strength, axis=0)
        sums["depth_p"] += np.nansum(ds["sbi_depth_p"].values, axis=0)
        sums["depth_z"] += np.nansum(depth_z, axis=0)
        sums["top_p"] += np.nansum(ds["sbi_top_p"].values, axis=0)
        sums["strength_uncond"] += np.nan_to_num(strength).sum(axis=0)
        for key, var in (("dt850", "dt_850_2m"), ("dt925", "dt_925_1000")):
            vals = ds[var].values.astype(float)
            sums[key] += np.nansum(vals, axis=0)
            sums["n_" + key] += np.isfinite(vals).sum(axis=0)

        for it in range(found.shape[0]):
            s, f = strength[it], found[it]
            wsum = w2d.sum()
            ts_time.append(ds["valid_time"].values[it])
            ts_found.append(float((f * w2d).sum() / wsum))
            wf = (w2d * f).sum()
            ts_cond.append(float(np.nansum(s * w2d) / wf) if wf > 0 else np.nan)
            ts_uncond.append(float(np.nan_to_num(s * w2d).sum() / wsum))
            ok = f & np.isfinite(s)
            hist_strength += np.histogram(s[ok], STRENGTH_BINS, weights=w2d[ok])[0]
            ok2 = ok & np.isfinite(depth_z[it])
            hist_sd += np.histogram2d(s[ok2], depth_z[it][ok2],
                                      [STRENGTH_BINS, DEPTH_BINS],
                                      weights=w2d[ok2])[0]
        expver = ds.attrs.get("expver_plev", "unknown")
        ds.close()

    with np.errstate(invalid="ignore", divide="ignore"):
        freq = sums["found"] / n_time
        cond = np.where(sums["found"] > 0, sums["strength_found"] / sums["found"], np.nan)
        depth_p = np.where(sums["found"] > 0, sums["depth_p"] / sums["found"], np.nan)
        depth_z = np.where(sums["found"] > 0, sums["depth_z"] / sums["found"], np.nan)
        top_p = np.where(sums["found"] > 0, sums["top_p"] / sums["found"], np.nan)
        dt850 = np.where(sums["n_dt850"] > 0, sums["dt850"] / sums["n_dt850"], np.nan)
        dt925 = np.where(sums["n_dt925"] > 0, sums["dt925"] / sums["n_dt925"], np.nan)

    # cells outside the requested area exist only because a projected bounding
    # box is rectangular; blank them rather than reporting them as "no SBI"
    freq, cond, depth_p, depth_z, top_p, dt850, dt925 = (
        np.where(valid, f, np.nan)
        for f in (freq, cond, depth_p, depth_z, top_p, dt850, dt925))
    strength_uncond = np.where(valid, sums["strength_uncond"] / n_time, np.nan)

    hdims = grid_hdims(template)
    out = xr.Dataset(
        {
            "sbi_frequency": (hdims, freq.astype(np.float32)),
            "sbi_strength_cond": (hdims, cond.astype(np.float32)),
            "sbi_strength_uncond": (hdims, strength_uncond.astype(np.float32)),
            "sbi_depth_p_mean": (hdims, depth_p.astype(np.float32)),
            "sbi_depth_z_mean": (hdims, depth_z.astype(np.float32)),
            "sbi_top_p_mean": (hdims, top_p.astype(np.float32)),
            "dt_850_2m_mean": (hdims, dt850.astype(np.float32)),
            "dt_925_1000_mean": (hdims, dt925.astype(np.float32)),
            "hist_strength": (("strength_bin",), hist_strength),
            "hist_strength_depth": (("strength_bin", "depth_bin"), hist_sd),
            "ts_sbi_frequency": (("time",), np.array(ts_found, dtype=np.float32)),
            "ts_strength_cond": (("time",), np.array(ts_cond, dtype=np.float32)),
            "ts_strength_uncond": (("time",), np.array(ts_uncond, dtype=np.float32)),
        },
        coords={
            **horizontal_coords(template),
            "strength_bin": 0.5 * (STRENGTH_BINS[:-1] + STRENGTH_BINS[1:]),
            "depth_bin": 0.5 * (DEPTH_BINS[:-1] + DEPTH_BINS[1:]),
            "time": np.array(ts_time),
        },
    )
    for name in (GRID_MAPPING_VAR, "domain_mask"):
        if name in template.variables:
            out[name] = template[name]
    units = {"sbi_frequency": "1", "sbi_strength_cond": "K", "sbi_strength_uncond": "K",
             "sbi_depth_p_mean": "hPa", "sbi_depth_z_mean": "m", "sbi_top_p_mean": "hPa",
             "dt_850_2m_mean": "K", "dt_925_1000_mean": "K",
             "ts_sbi_frequency": "1", "ts_strength_cond": "K", "ts_strength_uncond": "K"}
    for name, unit in units.items():
        out[name].attrs["units"] = unit
    out.attrs = {
        "title": f"Monthly {label} Arctic temperature-inversion statistics",
        "source_label": label,
        "n_time_steps": n_time,
        "n_days": len(files),
        "note": ("area statistics weighted by relative cell area (cos(latitude) "
                 "on a regular grid, (1+sin(latitude))^2 on the CARRA-2 polar "
                 "stereographic); *_cond over detected SBIs only (intensity, "
                 "Zhang et al. 2011); *_uncond counts not-found as 0 K"),
        "expver_last_file": expver,
    }
    return out


def fig_maps(out: xr.Dataset, year: int, month: int, outdir: Path) -> Path:
    import cartopy.crs as ccrs

    gkw = grid_kwargs(out)
    panels = [
        ("sbi_frequency", "SBI frequency  (fraction of time)", "seq", dict(vmax=1.0)),
        ("sbi_strength_cond", "SBI strength, detected cases  (K)", "seq", {}),
        ("sbi_depth_z_mean", "SBI depth, detected cases  (m)", "seq", {}),
        ("dt_850_2m_mean", "T(850 hPa) − T(2 m)  (K)", "div", {}),
        ("dt_925_1000_mean", "T(925 hPa) − T(1000 hPa)  (K)", "div", {}),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(14.5, 11.4), layout="constrained",
                             subplot_kw={"projection": ccrs.NorthPolarStereo()})
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.95))
    axes = axes.ravel()
    for i, (var, label, kind, kw) in enumerate(panels):
        polar_panel(axes[i], field=out[var].values, kind=kind,
                    cbar_label=label, **gkw, **kw)
        panel_label(axes[i], "abcdefgh"[i], x=-0.02, y=1.05)
    axes[-1].set_visible(False)
    fig.suptitle(f"{out.attrs.get('source_label', 'ERA5')} monthly temperature-inversion statistics — "
                 f"{calendar.month_name[month]} {year} — gray: no data / below ground",
                 y=0.99)
    path = outdir / f"monthly_maps_{year:04d}{month:02d}.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_distributions(out: xr.Dataset, year: int, month: int, outdir: Path) -> Path:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.4), layout="constrained")

    s_bin = out["strength_bin"].values
    hist = out["hist_strength"].values
    pdf = hist / hist.sum() / (s_bin[1] - s_bin[0])
    ax1.bar(s_bin, pdf, width=s_bin[1] - s_bin[0], color="#0072B2",
            edgecolor="white", linewidth=0.3)
    med = s_bin[np.searchsorted(np.cumsum(hist) / hist.sum(), 0.5)]
    ax1.axvline(med, color="#D55E00", lw=1.5, ls="--")
    ax1.annotate(f"median {med:.1f} K", (med, ax1.get_ylim()[1] * 0.95),
                 xytext=(6, 0), textcoords="offset points", color="#D55E00")
    ax1.set_xlabel("SBI strength (K)")
    ax1.set_ylabel("probability density (K$^{-1}$)")
    ax1.set_title("detected SBIs, area-weighted")
    ax1.grid(color="#e3e3e3", lw=0.6)
    ax1.set_axisbelow(True)
    for side in ("top", "right"):
        ax1.spines[side].set_visible(False)
    panel_label(ax1, "a", x=-0.12, y=1.06)

    h2 = out["hist_strength_depth"].values
    d_bin = out["depth_bin"].values
    mesh = ax2.pcolormesh(s_bin, d_bin, (h2 / h2.sum()).T, cmap="viridis",
                          shading="nearest")
    ax2.set_xlabel("SBI strength (K)")
    ax2.set_ylabel("SBI depth (m)")
    ax2.set_title("strength vs. depth, detected SBIs")
    cb = plt.colorbar(mesh, ax=ax2, fraction=0.05)
    cb.set_label("fraction")
    panel_label(ax2, "b", x=-0.16, y=1.06)

    fig.suptitle(f"SBI strength distributions — {calendar.month_name[month]} {year}")
    path = outdir / f"monthly_distributions_{year:04d}{month:02d}.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_timeseries(out: xr.Dataset, year: int, month: int, outdir: Path) -> Path:
    t = out["time"].values
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.0, 5.6), sharex=True,
                                   layout="constrained")
    ax1.plot(t, out["ts_strength_cond"], color="#0072B2", lw=1.5,
             label="detected cases (conditional)")
    ax1.plot(t, out["ts_strength_uncond"], color="#D55E00", lw=1.5,
             label="all points (not-found = 0)")
    ax1.set_ylabel("mean SBI strength (K)")
    ax1.legend(frameon=False, loc="upper right", ncol=2)
    panel_label(ax1, "a", x=-0.08, y=1.1)

    ax2.plot(t, out["ts_sbi_frequency"], color="#555555", lw=1.5)
    ax2.set_ylabel("SBI frequency")
    ax2.set_ylim(0, 1)
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    panel_label(ax2, "b", x=-0.08, y=1.1)

    for ax in (ax1, ax2):
        ax.grid(color="#e3e3e3", lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.suptitle(f"Domain-mean (area-weighted) inversion metrics — "
                 f"{calendar.month_name[month]} {year}")
    path = outdir / f"monthly_timeseries_{year:04d}{month:02d}.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--month", type=int, required=True)
    parser.add_argument("--days", nargs="+", default=None,
                        help="subset of days (default: whole month)")
    parser.add_argument("--skip-missing", action="store_true",
                        help="aggregate whatever daily files exist")
    parser.add_argument("--source", choices=list(SOURCES), default=None,
                        help="data source (default from config.yaml)")
    parser.add_argument("--config", default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args(argv)

    apply_agu_style()
    cfg = load_config(args.config, source=args.source)
    ndays = calendar.monthrange(args.year, args.month)[1]
    tokens = args.days if args.days else [f"1-{ndays}"]
    days = parse_days(tokens, args.year, args.month)

    files, missing = [], []
    for day in days:
        p = inversion_path(cfg, dt.date(args.year, args.month, day))
        (files if p.exists() else missing).append(p)
    if missing and not args.skip_missing:
        sys.exit(f"missing {len(missing)} daily file(s), first: {missing[0]}\n"
                 f"fix with: python src/daily_inversion.py --year {args.year} "
                 f"--month {args.month} --days 1-{ndays}\n"
                 f"or rerun with --skip-missing to aggregate the "
                 f"{len(files)} available day(s)")
    if not files:
        sys.exit("no daily derived files found for the requested days")
    if missing:
        print(f"note: skipping {len(missing)} missing day(s), using {len(files)}")

    target = monthly_path(cfg, args.year, args.month)
    if target.exists() and not args.overwrite:
        print(f"{target} exists, loading it (use --overwrite to recompute)")
        out = open_era5(target)
    else:
        out = aggregate(files, source_label(cfg))
        target.parent.mkdir(parents=True, exist_ok=True)
        out.to_netcdf(target, encoding={v: {"zlib": True, "complevel": 4}
                                        for v in out.data_vars})
        print(f"wrote {target}")

    w = area_weights(out)
    freq = out["sbi_frequency"].values
    cond = out["sbi_strength_cond"].values
    print(f"  domain mean SBI frequency      : {weighted_mean(freq, w):.1%}")
    print(f"  domain mean conditional strength: {weighted_mean(cond, w):.2f} K")

    if not args.no_figures:
        outdir = figures_dir(cfg)
        outdir.mkdir(parents=True, exist_ok=True)
        for fn in (fig_maps, fig_distributions, fig_timeseries):
            print(f"wrote {fn(out, args.year, args.month, outdir)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
