#!/usr/bin/env python
"""Plot temperature profiles illustrating the inversion-strength metrics.

For selected grid points of one time snapshot, draws the ERA5 temperature
profile (above-ground levels only) with the detected surface-based inversion
layer shaded and all three metric values annotated, so the derived numbers can
be checked visually against the raw profile.

Examples:
    python src/plot_profiles.py --year 2025 --month 1 --day 1 --hour 12
    python src/plot_profiles.py --year 2025 --month 1 --day 1 --hour 12 --lat 85 --lon -150
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reanlib.config import figures_dir, inversion_path, load_config, plev_path, sfc_path, source_label
from reanlib.inversion import column_heights
from reanlib.io_era5 import open_era5
from reanlib.plotstyle import apply_agu_style, panel_label

# Okabe-Ito CVD-safe hues; identity is also carried by marker shape + labels.
C_PROFILE = "#333333"
C_T2M = "#0072B2"       # blue, circle
C_SKT = "#D55E00"       # vermillion, diamond
C_SHADE = "#0072B2"     # inversion layer tint (low alpha)
C_GRID = "#d9d9d9"


def pick_points(inv: "xr.Dataset", how: list[str]) -> list[dict]:
    """Select illustration points from one time snapshot of the derived dataset."""
    import xarray as xr  # noqa: F401  (type hint only)

    strength = inv["sbi_strength"]
    found = inv["sbi_found"].astype(bool)
    points = []
    for tag in how:
        if tag == "strongest":
            target = strength.fillna(-np.inf).argmax(("latitude", "longitude"))
        elif tag == "weakest":
            # clearest no-inversion column: most negative T850-T2m among not-found
            cand = inv["dt_850_2m"].where(~found)
            if np.isfinite(cand).any():
                target = cand.fillna(np.inf).argmin(("latitude", "longitude"))
            else:  # every column has an inversion: take the weakest one
                target = strength.fillna(np.inf).argmin(("latitude", "longitude"))
        elif tag == "median":
            med = strength.median()
            target = abs(strength - med).fillna(np.inf).argmin(("latitude", "longitude"))
        else:
            raise ValueError(f"unknown point selector {tag!r}")
        lat = float(inv["latitude"][target["latitude"]])
        lon = float(inv["longitude"][target["longitude"]])
        points.append({"lat": lat, "lon": lon, "label": tag})
    return points


def plot_profile_panel(ax, ds_plev, inv, lat: float, lon: float, label: str,
                       src_lab: str = "ERA5") -> None:
    col = ds_plev.sel(latitude=lat, longitude=lon, method="nearest")
    ivn = inv.sel(latitude=lat, longitude=lon, method="nearest")
    sp_hpa = float(ivn["sp"]) / 100.0
    t2m, skt = float(ivn["t2m"]), float(ivn["skt"])

    p = col["pressure_level"].values.astype(float)      # descending
    T = col["t"].values.astype(float)
    q = col["q"].values.astype(float) if "q" in col else None
    z = column_heights(T, p, t2m, sp_hpa, q)
    above = p <= sp_hpa
    show = above & (p >= 400)

    # detected SBI layer
    top_p = float(ivn["sbi_top_p"])
    if np.isfinite(top_p):
        ax.axhspan(top_p, sp_hpa, color=C_SHADE, alpha=0.12, zorder=0)
        ax.axhline(top_p, color=C_SHADE, lw=1, alpha=0.5, zorder=1)

    tc = T - 273.15
    ax.plot(tc[show], p[show], "-", color=C_PROFILE, lw=2, zorder=3,
            label=f"{src_lab} T profile")
    for ref_p, m in ((925, "s"), (850, "s")):
        if ref_p <= sp_hpa:
            k = int(np.argmin(np.abs(p - ref_p)))
            ax.plot(tc[k], p[k], m, ms=5, mfc="white", mec=C_PROFILE, mew=1.2,
                    zorder=4)
            ax.annotate(f"{ref_p}", (tc[k], p[k]), textcoords="offset points",
                        xytext=(-6, -3), ha="right", fontsize=8, color="#666666")
    ax.plot(t2m - 273.15, sp_hpa, "o", ms=8, color=C_T2M, zorder=5, label="T 2 m")
    ax.plot(skt - 273.15, sp_hpa, "D", ms=7, color=C_SKT, zorder=5, label="T skin")

    ax.set_yscale("log")
    ax.set_ylim(1010 if sp_hpa <= 1000 else sp_hpa + 10, 400)
    yticks = [1000, 925, 850, 700, 600, 500, 400]
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(t) for t in yticks])
    ax.minorticks_off()
    ax.grid(color=C_GRID, lw=0.6, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    # right-hand axis: approximate height of the same pressure ticks, this column
    axr = ax.twinx()
    axr.set_yscale("log")
    axr.set_ylim(ax.get_ylim())
    axr.set_yticks(yticks)
    zlab = []
    for tp in yticks:
        if tp <= sp_hpa:
            k = int(np.argmin(np.abs(p - tp)))
            zlab.append(f"{z[k] / 1000:.1f}" if np.isfinite(z[k]) else "")
        else:
            zlab.append("")
    axr.set_yticklabels(zlab, fontsize=8, color="#666666")
    axr.minorticks_off()
    axr.spines[["top", "left", "right"]].set_visible(False)
    axr.tick_params(length=0)
    axr.set_ylabel("approx. height (km)", fontsize=8, color="#666666")

    def fmt(v, unit=""):
        return f"{float(v):.1f}{unit}" if np.isfinite(float(v)) else "none"

    text = (f"SBI strength: {fmt(ivn['sbi_strength'], ' K')}\n"
            f"SBI top: {fmt(ivn['sbi_top_p'], ' hPa')}"
            f"  depth: {fmt(ivn['sbi_depth_p'], ' hPa')} / {fmt(ivn['sbi_depth_z'], ' m')}\n"
            f"T850 − T2m: {fmt(ivn['dt_850_2m'], ' K')}\n"
            f"T925 − T1000: {fmt(ivn['dt_925_1000'], ' K')}")
    ax.text(0.03, 0.97, text, transform=ax.transAxes, va="top", ha="left",
            fontsize=8, family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#cccccc", alpha=0.9))

    lon_txt = f"{abs(lon):g}°{'E' if lon >= 0 else 'W'}"
    ax.set_title(f"{label}  ({lat:g}°N, {lon_txt})", fontsize=10)
    ax.set_xlabel("temperature (°C)")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--month", type=int, required=True)
    parser.add_argument("--day", type=int, required=True)
    parser.add_argument("--hour", type=int, required=True)
    parser.add_argument("--points", nargs="+", default=["strongest", "median", "weakest"],
                        choices=["strongest", "median", "weakest"])
    parser.add_argument("--lat", type=float, default=None,
                        help="plot a single explicit point instead of --points")
    parser.add_argument("--lon", type=float, default=None)
    parser.add_argument("--source", choices=["era5", "merra2"], default=None,
                        help="data source (default from config.yaml)")
    parser.add_argument("--config", default=None)
    parser.add_argument("--outdir", default=None, help="default: figures/ from config")
    args = parser.parse_args(argv)

    apply_agu_style()
    cfg = load_config(args.config, source=args.source)
    date = dt.date(args.year, args.month, args.day)
    when = np.datetime64(f"{date:%Y-%m-%d}T{args.hour:02d}:00")

    inv_file = inversion_path(cfg, date)
    if not inv_file.exists():
        sys.exit(f"missing {inv_file}\nfix with: python src/daily_inversion.py "
                 f"--year {date.year} --month {date.month} --days {date.day}")
    inv = open_era5(inv_file).sel(valid_time=when)
    ds_plev = open_era5(plev_path(cfg, date)).sel(valid_time=when)

    if (args.lat is None) != (args.lon is None):
        sys.exit("--lat and --lon must be given together")
    if args.lat is not None:
        points = [{"lat": args.lat, "lon": args.lon, "label": "selected"}]
    else:
        points = pick_points(inv, args.points)

    n = len(points)
    fig, axes = plt.subplots(1, n, figsize=(4.6 * n, 5.2), squeeze=False)
    for i, (ax, pt) in enumerate(zip(axes[0], points)):
        plot_profile_panel(ax, ds_plev, inv, pt["lat"], pt["lon"], pt["label"],
                           src_lab=source_label(cfg))
        panel_label(ax, "abcdefgh"[i], x=-0.16, y=1.08)
    axes[0][0].set_ylabel("pressure (hPa)")
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=9)
    fig.suptitle(f"{source_label(cfg)} temperature profiles & inversion metrics — "
                 f"{date:%Y-%m-%d} {args.hour:02d} UTC\n"
                 f"shaded: detected surface-based inversion layer", fontsize=11)
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))

    outdir = Path(args.outdir) if args.outdir else figures_dir(cfg)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"profiles_{date:%Y%m%d}T{args.hour:02d}Z.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
