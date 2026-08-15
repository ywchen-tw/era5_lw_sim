#!/usr/bin/env python
"""Level-by-level profile comparison of a reanalysis against MOSAiC radiosondes.

Stage 6b. Where `mosaic_compare.py` scores the derived inversion metrics, this
compares the *state* itself — temperature and humidity at standard pressure
levels — so a bias can be attributed to a layer rather than to a metric.

For every level-2 sounding in the month (Maturilli et al. 2021) the ascent is
interpolated in log-p onto the report levels and matched to the nearest
analysis time (one per sounding). The nearest reanalysis cell is found PER
LEVEL, at the balloon's own drifted position — the level-2 files carry
latitude/longitude for every sample, and the drift (median 16.5 km by
300 hPa) exceeds every source's grid spacing aloft, so matching everything
to the launch point would misplace exactly the layers where the sources
disagree on humidity. Three quantities are compared per level:

  T    [K]      directly
  q    [g/kg]   observed q derived from the sonde's RH, T and p
  RH   [%]      recomputed from q for EVERY source, including the observations

RH is deliberately not taken from each product as published: ERA5 reports RH
over a mixed phase below 0 C while radiosondes and CARRA report it over water,
and at Arctic winter temperatures that convention alone is worth tens of
percent. Recomputing all of them from q with one saturation formula
(reanlib/humidity.py) makes the columns comparable; the cost is that the
comparison tests q, not each product's own RH diagnostic.

Below-ground levels (p > surface pressure) are excluded, as are levels outside
the sounding's own pressure range — the observed profile is never
extrapolated.

Outputs derived/<source>/YYYY/MM/<source>_mosaic_profiles_YYYYMM.nc and, with
--report, a three-source figure plus a printed table.

Examples:
    python src/mosaic_profiles.py --year 2020 --month 1 --source era5
    python src/mosaic_profiles.py --year 2020 --month 1 --report
"""

from __future__ import annotations

import argparse
import calendar
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reanlib.config import (SOURCE_LABELS, SOURCES, figures_dir, inversion_path,
                            load_config, plev_path)
from reanlib.grid import GridIndex
from reanlib.humidity import relative_humidity_from_q, specific_humidity_from_rh
from reanlib.io_era5 import open_era5
from reanlib.mosaic import load_sounding_month
from reanlib.plotstyle import apply_agu_style, panel_label

#: Pressure levels present in all three sources, covering the inversion layer
#: and the free troposphere above it.
REPORT_LEVELS = np.array([1000, 950, 925, 900, 875, 850, 825, 800, 750, 700,
                          600, 500, 400, 300], dtype=float)

#: Colour per source, Okabe-Ito (CVD-safe).
SOURCE_COLOURS = {"era5": "#D55E00", "merra2": "#0072B2", "carra2": "#009E73"}


def soundings_in_month(year: int, month: int) -> list[dict]:
    """One entry per ascent: launch time, position, and the profile itself."""
    df = load_sounding_month(year, month)
    if df is None:
        sys.exit(f"no MOSAiC level-2 soundings registered for {year}-{month:02d}")
    need = {"Event", "PPPP [hPa]", "TTT [°C]", "RH [%]", "Latitude", "Longitude"}
    missing = need - set(df.columns)
    if missing:
        sys.exit(f"sounding file lacks column(s) {sorted(missing)}; "
                 "reanlib/mosaic.load_sounding_month must read them")
    out = []
    for event, g in df.groupby("Event", sort=True):
        g = g.dropna(subset=["PPPP [hPa]", "TTT [°C]"])
        if len(g) < 20:
            continue
        first = g.iloc[0]
        out.append({
            "event": event,
            "time": g["Date/Time"].min(),
            "lat": float(first["Latitude"]),
            "lon": float(first["Longitude"]),
            "p": g["PPPP [hPa]"].to_numpy(dtype=float),
            "t": g["TTT [°C]"].to_numpy(dtype=float) + 273.15,
            "rh": g["RH [%]"].to_numpy(dtype=float),
            "lat_prof": g["Latitude"].to_numpy(dtype=float),
            "lon_prof": g["Longitude"].to_numpy(dtype=float),
        })
    return out


def observed_on_levels(snd: dict, levels: np.ndarray):
    """Interpolate one ascent onto ``levels`` in log-p. No extrapolation."""
    p, t, rh = snd["p"], snd["t"], snd["rh"]
    good = np.isfinite(p) & np.isfinite(t)
    p, t, rh = p[good], t[good], rh[good]
    order = np.argsort(p)                       # np.interp needs ascending x
    p, t, rh = p[order], t[order], rh[order]
    inside = (levels >= p.min()) & (levels <= p.max())
    lp, lt, lrh = np.log(p), t, rh
    t_out = np.where(inside, np.interp(np.log(levels), lp, lt), np.nan)
    rh_out = np.where(inside, np.interp(np.log(levels), lp, lrh), np.nan)
    with np.errstate(invalid="ignore"):
        q_out = specific_humidity_from_rh(rh_out, t_out, levels * 100.0)
    return t_out, q_out, rh_out


def positions_on_levels(snd: dict, levels: np.ndarray):
    """Balloon position interpolated onto ``levels`` in log-p.

    Interpolates the components of the position's 3-D unit vector rather than
    latitude/longitude directly, so the result stays correct next to the pole
    and across the date line. Outside the ascent's pressure range np.interp
    clamps to the endpoint positions, and a sounding without usable positions
    falls back to its launch point — in both cases the observed values there
    are NaN anyway, so no pair is scored at a made-up position.
    """
    p, lat, lon = snd["p"], snd["lat_prof"], snd["lon_prof"]
    good = np.isfinite(p) & np.isfinite(lat) & np.isfinite(lon)
    if good.sum() < 2:
        return (np.full(levels.size, snd["lat"]),
                np.full(levels.size, snd["lon"]))
    p, lat, lon = p[good], lat[good], lon[good]
    order = np.argsort(p)
    latr, lonr = np.radians(lat[order]), np.radians(lon[order])
    xyz = np.stack([np.cos(latr) * np.cos(lonr),
                    np.cos(latr) * np.sin(lonr),
                    np.sin(latr)])
    lp = np.log(p[order])
    comp = np.stack([np.interp(np.log(levels), lp, c) for c in xyz])
    comp /= np.linalg.norm(comp, axis=0)
    lat_out = np.degrees(np.arcsin(np.clip(comp[2], -1.0, 1.0)))
    lon_out = np.degrees(np.arctan2(comp[1], comp[0]))
    return lat_out, lon_out


def build_pairs(cfg: dict, year: int, month: int, levels: np.ndarray,
                max_dt_h: float = 3.0) -> xr.Dataset:
    """Pair every sounding with its matched reanalysis profile."""
    snds = soundings_in_month(year, month)
    nl = levels.size
    fields = {k: np.full((len(snds), nl), np.nan) for k in
              ("obs_t", "obs_q", "obs_rh", "rean_t", "rean_q", "rean_rh",
               "match_km")}
    meta = {k: np.full(len(snds), np.nan) for k in ("lat", "lon", "match_dt_h")}
    times = np.empty(len(snds), dtype="datetime64[ns]")

    cache: dict = {}
    index = None
    for i, snd in enumerate(snds):
        ts = np.datetime64(snd["time"])
        times[i] = ts
        meta["lat"][i], meta["lon"][i] = snd["lat"], snd["lon"]
        obs = observed_on_levels(snd, levels)
        fields["obs_t"][i], fields["obs_q"][i], fields["obs_rh"][i] = obs
        lat_lev, lon_lev = positions_on_levels(snd, levels)

        t_near = ((ts + np.timedelta64(3, "h")).astype("datetime64[6h]")
                  .astype("datetime64[ns]"))
        date = t_near.astype("datetime64[D]").item()
        if date not in cache:
            cache.clear()
            pp, ip = plev_path(cfg, date), inversion_path(cfg, date)
            cache[date] = ((open_era5(pp), open_era5(ip))
                           if pp.exists() and ip.exists() else None)
            index = GridIndex(cache[date][0]) if cache[date] else None
        if cache[date] is None:
            continue
        plev, inv = cache[date]
        dt_h = abs((t_near - ts) / np.timedelta64(1, "h"))
        if dt_h > max_dt_h:
            continue
        try:
            snap = plev.sel(valid_time=t_near)
            sp_map = inv["sp"].sel(valid_time=t_near)
        except KeyError:
            continue
        meta["match_dt_h"][i] = dt_h

        # one analysis time per sounding (the balloon reaches 300 hPa ~28 min
        # after launch, small against the 6 h cadence), but the cell is
        # re-queried per level at the balloon's drifted position; columns are
        # cached because adjacent levels usually land in the same cell
        p_src = snap["pressure_level"].values.astype(float)
        cols: dict[tuple[int, int], tuple] = {}
        for j, lev in enumerate(levels):
            k = np.flatnonzero(np.isclose(p_src, lev))
            if k.size == 0:
                continue
            (iy, ix), km = index.query(lat_lev[j], lon_lev[j])
            if (iy, ix) not in cols:
                sel = index.isel(iy, ix)
                cell = snap.isel(**sel)
                cols[iy, ix] = (
                    cell["t"].values.astype(float),
                    cell["q"].values.astype(float) if "q" in cell else None,
                    float(sp_map.isel(**sel)) / 100.0,
                )
            t_src, q_src, sp_hpa = cols[iy, ix]
            if lev > sp_hpa:                     # below ground at this cell
                continue
            fields["match_km"][i, j] = km
            fields["rean_t"][i, j] = t_src[k[0]]
            if q_src is not None:
                fields["rean_q"][i, j] = q_src[k[0]]

    with np.errstate(invalid="ignore", divide="ignore"):
        fields["rean_rh"] = relative_humidity_from_q(
            fields["rean_q"], fields["rean_t"], levels[None, :] * 100.0)

    ds = xr.Dataset(
        {k: (("sounding", "level"), v) for k, v in fields.items()}
        | {k: (("sounding",), v) for k, v in meta.items()},
        coords={"sounding": np.arange(len(snds)), "level": levels,
                "time": (("sounding",), times)},
    )
    ds["obs_q"] *= 1000.0          # g/kg reads better than kg/kg at these values
    ds["rean_q"] *= 1000.0
    for name, unit in (("obs_t", "K"), ("rean_t", "K"), ("obs_q", "g kg**-1"),
                       ("rean_q", "g kg**-1"), ("obs_rh", "%"), ("rean_rh", "%"),
                       ("match_km", "km"), ("match_dt_h", "h")):
        ds[name].attrs["units"] = unit
    ds.attrs = {
        "title": f"MOSAiC radiosonde vs {SOURCE_LABELS[cfg['source']]} profiles",
        "source_label": SOURCE_LABELS[cfg["source"]],
        "n_soundings": len(snds),
        "rh_convention": ("RH recomputed from q over water for observations "
                          "and reanalysis alike (reanlib/humidity.py)"),
        "collocation": ("nearest cell PER LEVEL at the balloon's drifted "
                        "position (level-2 per-sample lat/lon, interpolated "
                        "on unit vectors in log-p); one analysis time per "
                        "sounding, nearest within 3 h of launch"),
    }
    return ds


def level_stats(ds: xr.Dataset, var: str) -> dict[str, np.ndarray]:
    """Per-level n, bias, rmse, spread (σ of differences) and correlation."""
    o = ds[f"obs_{var}"].values
    r = ds[f"rean_{var}"].values
    ok = np.isfinite(o) & np.isfinite(r)
    n = ok.sum(axis=0)
    bias = np.full(ds.sizes["level"], np.nan)
    rmse = np.full(ds.sizes["level"], np.nan)
    sigma = np.full(ds.sizes["level"], np.nan)
    corr = np.full(ds.sizes["level"], np.nan)
    for j in range(ds.sizes["level"]):
        m = ok[:, j]
        if m.sum() < 3:
            continue
        d = r[m, j] - o[m, j]
        bias[j] = d.mean()
        rmse[j] = np.sqrt((d ** 2).mean())
        sigma[j] = d.std()                       # rmse² = bias² + σ²
        if o[m, j].std() > 0 and r[m, j].std() > 0:
            corr[j] = np.corrcoef(o[m, j], r[m, j])[0, 1]
    return {"n": n, "bias": bias, "rmse": rmse, "sigma": sigma, "r": corr}


VARS = (("t", "temperature", "K"), ("q", "specific humidity", "g/kg"),
        ("rh", "relative humidity", "%"))


def print_table(per_source: dict[str, xr.Dataset]) -> str:
    lines = []
    for var, long, unit in VARS:
        lines.append(f"\n{long} ({unit}) — bias / rmse / r, by level")
        head = "  level  " + "".join(f"{SOURCE_LABELS[s]:>26}" for s in per_source)
        lines.append(head)
        levels = next(iter(per_source.values()))["level"].values
        stats = {s: level_stats(d, var) for s, d in per_source.items()}
        for j, lev in enumerate(levels):
            row = f"  {lev:>5.0f}  "
            for s in per_source:
                st = stats[s]
                row += (f"{st['bias'][j]:>+8.2f}{st['rmse'][j]:>8.2f}"
                        f"{st['r'][j]:>+7.2f}   " if np.isfinite(st["bias"][j])
                        else f"{'—':>25}   ")
            lines.append(row)
        n = {s: level_stats(d, var)["n"] for s, d in per_source.items()}
        lines.append("  n     " + "".join(f"{int(np.nanmax(v)):>26}" for v in n.values()))
    return "\n".join(lines)


def fig_report(per_source: dict[str, xr.Dataset], year: int, month: int,
               outdir: Path) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(13.0, 8.6), layout="constrained",
                             sharey=True)
    levels = next(iter(per_source.values()))["level"].values
    for col, (var, long, unit) in enumerate(VARS):
        for row, metric in enumerate(("bias", "rmse")):
            ax = axes[row, col]
            for src, ds in per_source.items():
                st = level_stats(ds, var)
                if metric == "bias":             # ±1σ sounding-to-sounding
                    ax.fill_betweenx(levels, st["bias"] - st["sigma"],
                                     st["bias"] + st["sigma"], lw=0,
                                     color=SOURCE_COLOURS[src], alpha=0.13)
                ax.plot(st[metric], levels, "-o", ms=3, lw=1.4,
                        color=SOURCE_COLOURS[src], label=SOURCE_LABELS[src])
            if metric == "bias":
                ax.axvline(0, color="#888888", lw=0.8, zorder=0)
            ax.set_xlabel(f"{metric} ({unit})")
            if col == 0:
                ax.set_ylabel("pressure (hPa)")
            ax.set_ylim(1010, 290)
            ax.grid(color="#e3e3e3", lw=0.5)
            ax.set_axisbelow(True)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            panel_label(ax, "abcdef"[row * 3 + col], x=-0.02, y=1.04)
            if row == 0:
                ax.set_title(long, fontsize=11)
    # figure-level legend below every panel: inside panel (a) it sat directly
    # on the 950-1000 hPa lines, which is where the largest biases are
    from matplotlib.patches import Patch

    handles, labels = axes[0, 0].get_legend_handles_labels()
    handles.append(Patch(facecolor="#888888", alpha=0.25, lw=0))
    labels.append("±1σ of differences")
    fig.legend(handles, labels, loc="outside lower center", ncol=len(labels),
               frameon=False, fontsize=10)
    fig.suptitle(f"Reanalysis vs MOSAiC radiosondes, {calendar.month_name[month]} "
                 f"{year} — profile bias (top) and rmse (bottom)", y=1.02)
    path = outdir / f"mosaic_profiles_{year:04d}{month:02d}.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def profiles_path(cfg: dict, year: int, month: int) -> Path:
    from reanlib.config import _month_dir

    return (_month_dir(cfg, "derived", year, month)
            / f"{cfg['source']}_mosaic_profiles_{year:04d}{month:02d}.nc")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--month", type=int, required=True)
    parser.add_argument("--source", choices=list(SOURCES), default=None)
    parser.add_argument("--report", action="store_true",
                        help="combine every source that has a file into one "
                             "figure and table")
    parser.add_argument("--config", default=None)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)

    apply_agu_style()
    if not args.report:
        cfg = load_config(args.config, source=args.source)
        target = profiles_path(cfg, args.year, args.month)
        if target.exists() and not args.overwrite:
            print(f"{target} exists (use --overwrite to recompute)")
            return 0
        ds = build_pairs(cfg, args.year, args.month, REPORT_LEVELS)
        target.parent.mkdir(parents=True, exist_ok=True)
        ds.to_netcdf(target)
        matched = int(np.isfinite(ds["match_dt_h"].values).sum())
        km = ds["match_km"].values
        print(f"wrote {target}")
        print(f"  {matched} / {ds.sizes['sounding']} soundings matched")
        with np.errstate(all="ignore"):
            med = np.nanmedian(km, axis=0)
        print("  median match km by level: "
              + "  ".join(f"{int(l)}:{m:.1f}" for l, m in
                          zip(ds["level"].values, med) if np.isfinite(m)))
        print(print_table({cfg["source"]: ds}))
        return 0

    per_source = {}
    for src in SOURCES:
        cfg = load_config(args.config, source=src)
        p = profiles_path(cfg, args.year, args.month)
        if p.exists():
            per_source[src] = xr.open_dataset(p)
    if not per_source:
        sys.exit("no per-source profile files yet — run without --report first")
    print(f"sources: {', '.join(SOURCE_LABELS[s] for s in per_source)}")
    print(print_table(per_source))
    outdir = figures_dir(load_config(args.config, source=next(iter(per_source))))
    outdir = outdir.parent
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"\nwrote {fig_report(per_source, args.year, args.month, outdir)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
