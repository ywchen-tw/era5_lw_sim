#!/usr/bin/env python
"""Polar-stereographic maps of inversion strength for one time snapshot.

Examples:
    python src/plot_maps.py --year 2025 --month 1 --day 1 --hour 12
    python src/plot_maps.py --year 2025 --month 1 --day 1 --hour 12 --metrics sbi
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
from reanlib.config import (SOURCES, figures_dir, inversion_path, load_config,
                            source_label)
from reanlib.io_era5 import open_era5
from reanlib.grid import domain_mask, latlon2d
from reanlib.mapping import grid_kwargs, polar_panel
from reanlib.plotstyle import apply_agu_style, panel_label

# metric key -> (variable, title, colormap kind)
METRICS = {
    "sbi": ("sbi_strength", "SBI strength  T(top) − T(2 m)", "seq"),
    "dt850": ("dt_850_2m", "T(850 hPa) − T(2 m)", "div"),
    "dt925": ("dt_925_1000", "T(925 hPa) − T(1000 hPa)", "div"),
}


def plot_metric_panel(ax, inv, key: str, south_lat: float, gkw: dict) -> None:
    var, title, kind = METRICS[key]
    field = inv[var].values
    if key == "sbi":
        # not-found = no inversion = 0 K strength, but only inside the domain:
        # a projected bounding box has corners outside the requested area, and
        # those carry no data rather than a zero-strength inversion
        field = np.where(domain_mask(inv), np.nan_to_num(field), np.nan)
    polar_panel(ax, field=field, kind=kind, cbar_label=f"{title}  (K)",
                south_lat=south_lat, **gkw)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--month", type=int, required=True)
    parser.add_argument("--day", type=int, required=True)
    parser.add_argument("--hour", type=int, required=True)
    parser.add_argument("--metrics", nargs="+", default=["sbi", "dt850", "dt925"],
                        choices=list(METRICS))
    parser.add_argument("--source", choices=list(SOURCES), default=None,
                        help="data source (default from config.yaml)")
    parser.add_argument("--config", default=None)
    parser.add_argument("--outdir", default=None, help="default: figures/ from config")
    args = parser.parse_args(argv)

    import cartopy.crs as ccrs

    apply_agu_style()
    cfg = load_config(args.config, source=args.source)
    date = dt.date(args.year, args.month, args.day)
    when = np.datetime64(f"{date:%Y-%m-%d}T{args.hour:02d}:00")

    inv_file = inversion_path(cfg, date)
    if not inv_file.exists():
        sys.exit(f"missing {inv_file}\nfix with: python src/daily_inversion.py "
                 f"--year {date.year} --month {date.month} --days {date.day}")
    inv = open_era5(inv_file).sel(valid_time=when)
    gkw = grid_kwargs(inv)
    # boundary exactly at the southernmost cell so no background ring shows
    south_lat = float(np.nanmin(latlon2d(inv)[0]))

    n = len(args.metrics)
    fig, axes = plt.subplots(1, n, figsize=(4.8 * n, 5.9), squeeze=False,
                             subplot_kw={"projection": ccrs.NorthPolarStereo()},
                             layout="constrained")
    # keep the axes region below the title so it clears the 180° labels
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.92))
    for i, (ax, key) in enumerate(zip(axes[0], args.metrics)):
        plot_metric_panel(ax, inv, key, south_lat, gkw)
        panel_label(ax, "abcdefgh"[i], x=-0.02, y=1.05)
    masked_note = (" — gray: level below ground"
                   if any(METRICS[k][2] == "div" for k in args.metrics) else "")
    fig.suptitle(f"{source_label(cfg)} temperature-inversion strength — "
                 f"{date:%Y-%m-%d} {args.hour:02d} UTC{masked_note}", y=0.99)

    outdir = Path(args.outdir) if args.outdir else figures_dir(cfg)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"maps_{date:%Y%m%d}T{args.hour:02d}Z.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
