#!/usr/bin/env python
"""Profile PCA and surface-temperature analysis of the monthly inversion data.

Three analyses over one month (samples = grid point x time):

1. PCA/EOF of lower-tropospheric temperature profiles (1000 hPa up to --top,
   default 400 hPa). Only samples with surface pressure >= 1000 hPa are used
   so every level is above ground (this excludes Greenland/terrain). The
   covariance matrix is accumulated from raw moments in one pass, so memory
   stays flat. EOFs are oriented so the lowest-level loading is positive.
2. Surface temperature: monthly mean and variability of t2m.
3. Correlations of SBI strength with t2m and with the leading PC scores —
   per-grid-point temporal correlations (maps) and domain-wide 2-D histograms.
   Strength is the unconditional series (not-found = 0 K) so the time series
   is gap-free; the correlation therefore mixes occurrence and intensity.

Outputs derived/YYYY/MM/era5_profile_analysis_YYYYMM.nc and two figures.

Examples:
    python src/era5_profile_analysis.py --year 2020 --month 1
"""

from __future__ import annotations

import argparse
import calendar
import datetime as dt
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from era5_download import parse_days
from era5lib.config import REPO_ROOT, figures_dir, inversion_path, load_config, plev_path
from era5lib.io_era5 import open_era5
from era5lib.mapping import polar_panel
from era5lib.plotstyle import apply_agu_style, panel_label

T2M_BINS = np.arange(-50.0, 1.0, 1.0)        # degC
STRENGTH_BINS = np.arange(0.0, 25.5, 0.5)    # K
PC_BINS = np.arange(-40.0, 40.5, 1.0)        # K (PC scores carry K units)
N_MODES = 6


def analysis_path(cfg: dict, year: int, month: int) -> Path:
    root = Path(cfg["paths"]["derived"])
    if not root.is_absolute():
        root = REPO_ROOT / root
    return (root / f"{year:04d}" / f"{month:02d}"
            / f"era5_profile_analysis_{year:04d}{month:02d}.nc")


class CorrAccum:
    """Streaming Pearson correlation between two (time-varying) fields."""

    def __init__(self, shape):
        self.n = np.zeros(shape)
        self.sx = np.zeros(shape)
        self.sy = np.zeros(shape)
        self.sxx = np.zeros(shape)
        self.syy = np.zeros(shape)
        self.sxy = np.zeros(shape)

    def add(self, x, y, valid):
        """x, y, valid: (ntime, *shape) — reduced over the leading time axis."""
        v = valid.astype(float)
        x = np.where(valid, x, 0.0)
        y = np.where(valid, y, 0.0)
        self.n += v.sum(axis=0)
        self.sx += x.sum(axis=0)
        self.sy += y.sum(axis=0)
        self.sxx += (x * x).sum(axis=0)
        self.syy += (y * y).sum(axis=0)
        self.sxy += (x * y).sum(axis=0)

    def corr(self, min_n: int = 30):
        with np.errstate(invalid="ignore", divide="ignore"):
            cov = self.n * self.sxy - self.sx * self.sy
            var = ((self.n * self.sxx - self.sx**2)
                   * (self.n * self.syy - self.sy**2))
            r = cov / np.sqrt(var)
        return np.where(self.n >= min_n, r, np.nan)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--month", type=int, required=True)
    parser.add_argument("--days", nargs="+", default=None)
    parser.add_argument("--top", type=float, default=400.0,
                        help="top pressure level of the PCA profiles (hPa)")
    parser.add_argument("--config", default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args(argv)

    apply_agu_style()
    cfg = load_config(args.config)
    ndays = calendar.monthrange(args.year, args.month)[1]
    tokens = args.days if args.days else [f"1-{ndays}"]
    days = parse_days(tokens, args.year, args.month)
    dates = [dt.date(args.year, args.month, d) for d in days]
    missing = [d for d in dates if not (plev_path(cfg, d).exists()
                                        and inversion_path(cfg, d).exists())]
    if missing:
        sys.exit(f"missing plev/derived files for {len(missing)} day(s), "
                 f"first: {missing[0]} — run era5_download.py / era5_inversion.py")

    target = analysis_path(cfg, args.year, args.month)
    if target.exists() and not args.overwrite:
        print(f"{target} exists, loading it (use --overwrite to recompute)")
        out = xr.open_dataset(target)
    else:
        out = compute(cfg, dates, args.top)
        target.parent.mkdir(parents=True, exist_ok=True)
        out.to_netcdf(target)
        print(f"wrote {target}")

    ev = out["explained_variance_ratio"].values
    print(f"  PCA explained variance: PC1 {ev[0]:.1%}, PC2 {ev[1]:.1%}, PC3 {ev[2]:.1%}")
    for name in ("r_strength_t2m", "r_strength_pc1", "r_strength_pc2"):
        w = np.cos(np.deg2rad(out["latitude"].values))[:, None]
        r = out[name].values
        ok = np.isfinite(r)
        print(f"  domain mean {name}: {np.nansum(r * w * ok) / (w * ok).sum():+.2f}")

    if not args.no_figures:
        outdir = figures_dir(cfg)
        outdir.mkdir(parents=True, exist_ok=True)
        print(f"wrote {fig_pca(out, args.year, args.month, outdir)}")
        print(f"wrote {fig_surface_t(out, args.year, args.month, outdir)}")
    return 0


def compute(cfg: dict, dates: list[dt.date], top_hpa: float) -> xr.Dataset:
    # ---- pass 1: raw moments of the profile vectors -> mean + covariance
    sum_x = m2 = None
    n_samples = 0
    lat = lon = levels = None
    for date in dates:
        plev = open_era5(plev_path(cfg, date))
        inv = open_era5(inversion_path(cfg, date))
        if lat is None:
            lat = plev["latitude"].values
            lon = plev["longitude"].values
            p = plev["pressure_level"].values.astype(float)
            levels = p[(p <= 1000.0) & (p >= top_hpa)]
            nlev = levels.size
            sum_x = np.zeros(nlev)
            m2 = np.zeros((nlev, nlev))
        T = (plev["t"].sel(pressure_level=levels)
             .transpose("valid_time", "pressure_level", "latitude", "longitude")
             .values.astype(np.float64))
        sp_hpa = inv["sp"].values / 100.0
        valid = sp_hpa >= 1000.0                       # (ntime, nlat, nlon)
        x = T.transpose(0, 2, 3, 1)[valid]             # (nsamples, nlev)
        sum_x += x.sum(axis=0)
        m2 += x.T @ x
        n_samples += x.shape[0]
        plev.close()
        inv.close()

    mean_profile = sum_x / n_samples
    cov = m2 / n_samples - np.outer(mean_profile, mean_profile)
    eigval, eigvec = np.linalg.eigh(cov)
    order = np.argsort(eigval)[::-1]
    eigval, eigvec = eigval[order], eigvec[:, order]
    eofs = eigvec[:, :N_MODES].T                       # (mode, nlev)
    signs = np.where(eofs[:, 0] >= 0, 1.0, -1.0)       # lowest level positive
    eofs = eofs * signs[:, None]
    evr = eigval[:N_MODES] / eigval.sum()

    # ---- pass 2: PC scores, correlations, histograms
    shape = (lat.size, lon.size)
    acc_t2m = CorrAccum(shape)
    acc_pc1 = CorrAccum(shape)
    acc_pc2 = CorrAccum(shape)
    h_t2m = np.zeros((T2M_BINS.size - 1, STRENGTH_BINS.size - 1))
    h_pc1 = np.zeros((PC_BINS.size - 1, STRENGTH_BINS.size - 1))
    sum_t2m = np.zeros(shape)
    sum_t2m2 = np.zeros(shape)
    n_time = 0
    w2d = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, lon.size))

    for date in dates:
        plev = open_era5(plev_path(cfg, date))
        inv = open_era5(inversion_path(cfg, date))
        T = (plev["t"].sel(pressure_level=levels)
             .transpose("valid_time", "pressure_level", "latitude", "longitude")
             .values.astype(np.float64))
        t2m = inv["t2m"].values.astype(np.float64)
        sp_hpa = inv["sp"].values / 100.0
        strength = np.nan_to_num(inv["sbi_strength"].values.astype(np.float64))
        valid = sp_hpa >= 1000.0
        anom = T.transpose(0, 2, 3, 1) - mean_profile
        pc1 = anom @ eofs[0]
        pc2 = anom @ eofs[1]
        n_time += T.shape[0]

        all_valid = np.ones_like(valid)
        acc_t2m.add(strength, t2m, all_valid)
        acc_pc1.add(strength, pc1, valid)
        acc_pc2.add(strength, pc2, valid)
        sum_t2m += t2m.sum(axis=0)
        sum_t2m2 += (t2m**2).sum(axis=0)

        for it in range(T.shape[0]):
            w = w2d.ravel()
            h_t2m += np.histogram2d((t2m[it] - 273.15).ravel(), strength[it].ravel(),
                                    [T2M_BINS, STRENGTH_BINS], weights=w)[0]
            v = valid[it].ravel()
            h_pc1 += np.histogram2d(pc1[it].ravel()[v], strength[it].ravel()[v],
                                    [PC_BINS, STRENGTH_BINS], weights=w[v])[0]
        plev.close()
        inv.close()

    t2m_mean = sum_t2m / n_time
    t2m_std = np.sqrt(np.maximum(sum_t2m2 / n_time - t2m_mean**2, 0.0))

    hdims = ("latitude", "longitude")
    out = xr.Dataset(
        {
            "mean_profile": (("level",), mean_profile),
            "eofs": (("mode", "level"), eofs),
            "explained_variance_ratio": (("mode",), evr),
            "r_strength_t2m": (hdims, acc_t2m.corr()),
            "r_strength_pc1": (hdims, acc_pc1.corr()),
            "r_strength_pc2": (hdims, acc_pc2.corr()),
            "t2m_mean": (hdims, t2m_mean.astype(np.float32)),
            "t2m_std": (hdims, t2m_std.astype(np.float32)),
            "hist_t2m_strength": (("t2m_bin", "strength_bin"), h_t2m),
            "hist_pc1_strength": (("pc_bin", "strength_bin"), h_pc1),
        },
        coords={
            "latitude": lat, "longitude": lon, "level": levels,
            "mode": np.arange(1, N_MODES + 1),
            "t2m_bin": 0.5 * (T2M_BINS[:-1] + T2M_BINS[1:]),
            "strength_bin": 0.5 * (STRENGTH_BINS[:-1] + STRENGTH_BINS[1:]),
            "pc_bin": 0.5 * (PC_BINS[:-1] + PC_BINS[1:]),
        },
    )
    out.attrs = {
        "title": "Monthly ERA5 profile PCA / surface-T / correlation analysis",
        "n_pca_samples": n_samples,
        "n_time_steps": n_time,
        "pca_levels": f"1000-{top_hpa:g} hPa, samples restricted to sp >= 1000 hPa",
        "strength_series": "unconditional SBI strength (not-found = 0 K)",
        "note": "EOFs oriented so the lowest-level loading is positive; "
                "PC scores in K (projection of T anomaly profiles onto EOFs)",
    }
    return out


def fig_pca(out: xr.Dataset, year: int, month: int, outdir: Path) -> Path:
    import cartopy.crs as ccrs

    fig = plt.figure(figsize=(13.2, 5.4), layout="constrained")
    gs = fig.add_gridspec(1, 4, width_ratios=[1.1, 0.9, 1.25, 1.25])
    colors = ["#0072B2", "#D55E00", "#555555"]

    ax = fig.add_subplot(gs[0])
    for m in range(3):
        ev = float(out["explained_variance_ratio"][m])
        ax.plot(out["eofs"][m], out["level"], color=colors[m], lw=2,
                label=f"EOF{m + 1} ({ev:.0%})")
    ax.axvline(0, color="#bbbbbb", lw=0.8)
    ax.set_yscale("log")
    ax.set_ylim(1010, float(out["level"].min()) - 10)
    yticks = [1000, 925, 850, 700, 600, 500, 400]
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(t) for t in yticks])
    ax.minorticks_off()
    ax.set_xlabel("EOF loading")
    ax.set_ylabel("pressure (hPa)")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(color="#e3e3e3", lw=0.6)
    ax.set_axisbelow(True)
    panel_label(ax, "a", x=-0.22, y=1.05)

    ax = fig.add_subplot(gs[1])
    evr = out["explained_variance_ratio"].values
    ax.bar(np.arange(1, evr.size + 1), evr * 100, color="#0072B2")
    ax.set_xlabel("mode")
    ax.set_ylabel("explained variance (%)")
    ax.grid(color="#e3e3e3", lw=0.6, axis="y")
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    panel_label(ax, "b", x=-0.24, y=1.05)

    lat, lon = out["latitude"].values, out["longitude"].values
    for i, (var, label) in enumerate((
            ("r_strength_pc1", "corr(SBI strength, PC1)"),
            ("r_strength_pc2", "corr(SBI strength, PC2)"))):
        axm = fig.add_subplot(gs[2 + i], projection=ccrs.NorthPolarStereo())
        polar_panel(axm, lat, lon, out[var].values, kind="div",
                    vmin=-1.0, vmax=1.0, cbar_label=label, cbar_label_size=10)
        panel_label(axm, "cd"[i], x=-0.02, y=1.08)

    fig.suptitle(f"Temperature-profile PCA and its relation to SBI strength — "
                 f"{calendar.month_name[month]} {year}", y=1.06)
    path = outdir / f"profile_pca_{year:04d}{month:02d}.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def fig_surface_t(out: xr.Dataset, year: int, month: int, outdir: Path) -> Path:
    import cartopy.crs as ccrs

    lat, lon = out["latitude"].values, out["longitude"].values
    fig = plt.figure(figsize=(13.2, 11.0), layout="constrained")
    gs = fig.add_gridspec(2, 2)

    t2m_c = out["t2m_mean"].values - 273.15
    ax = fig.add_subplot(gs[0, 0], projection=ccrs.NorthPolarStereo())
    polar_panel(ax, lat, lon, t2m_c, kind="seq", cmap=plt.get_cmap("magma"),
                vmin=float(np.nanpercentile(t2m_c, 1)),
                vmax=float(np.nanpercentile(t2m_c, 99)),
                cbar_label="monthly mean T 2 m  (°C)", cbar_label_size=10)
    panel_label(ax, "a", x=-0.02, y=1.05)

    ax = fig.add_subplot(gs[0, 1], projection=ccrs.NorthPolarStereo())
    polar_panel(ax, lat, lon, out["t2m_std"].values, kind="seq",
                cbar_label="T 2 m std. dev.  (K)", cbar_label_size=10)
    panel_label(ax, "b", x=-0.02, y=1.05)

    ax = fig.add_subplot(gs[1, 0], projection=ccrs.NorthPolarStereo())
    polar_panel(ax, lat, lon, out["r_strength_t2m"].values, kind="div",
                vmin=-1.0, vmax=1.0,
                cbar_label="corr(SBI strength, T 2 m)", cbar_label_size=10)
    panel_label(ax, "c", x=-0.02, y=1.05)

    ax = fig.add_subplot(gs[1, 1])
    h = out["hist_t2m_strength"].values
    frac = h / h.sum()
    mesh = ax.pcolormesh(out["t2m_bin"], out["strength_bin"], frac.T,
                         cmap="viridis", shading="nearest")
    ax.set_xlabel("T 2 m (°C)")
    ax.set_ylabel("SBI strength (K)")
    ax.set_title("all samples (not-found = 0 K), area-weighted")
    cb = plt.colorbar(mesh, ax=ax, fraction=0.05)
    cb.set_label("fraction")
    panel_label(ax, "d", x=-0.14, y=1.06)

    fig.suptitle(f"Surface temperature and its relation to SBI strength — "
                 f"{calendar.month_name[month]} {year}", y=1.03)
    path = outdir / f"surface_t_{year:04d}{month:02d}.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


if __name__ == "__main__":
    sys.exit(main())
