#!/usr/bin/env python
"""State-time test: does averaging two hourly snapshots close the LW gap?

ERA5 surface radiation (strd/str) is accumulated over the hour ENDING at
valid_time, while the stage-7 simulations are instantaneous snapshots of the
analysis at valid_time. This script tests the stage-7 diagnosis of the
clear-sky LWdn scatter by comparing three estimators of the 11-12Z ERA5
accumulation at the same pixels:

    inst 12Z         the end-of-window snapshot (what stage 7 uses)
    inst 11Z         the start-of-window snapshot
    (11Z + 12Z)/2    trapezoidal average over the accumulation window

Inputs are the stage-7 manifests/results for both hours; the second hour is
produced with  lrt_sim.py prep --hour 11 --pixels-from <12Z manifest>
so the pixel set is identical. Pixels that are no longer cloud-free at the
earlier hour are excluded from the headline statistics (their ERA5
accumulation contains cloudy radiation) and reported separately.

Runs under either conda env (numpy/pandas/matplotlib only).

Example:
    python src/statetime_test.py --year 2020 --month 1 --day 1 --hour 12
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reanlib.config import figures_dir, load_config
from lrt_sim import SIGMA, manifest_path, planck_band_fraction, results_path

CLEAR_CC_MAX = 0.01     # same eligibility thresholds as lrt_sim.py prep
CLEAR_PATH_G = 1.0


def load_hour(cfg, date: dt.date, hour: int, sky: str):
    """Manifest profiles + libRadtran/RRTMG results for one snapshot, by label."""
    import pandas as pd

    manifest = json.loads(manifest_path(cfg, date, hour, sky).read_text())
    res = pd.read_csv(results_path(cfg, date, hour, sky))
    sims = res["simulator"].astype(str)
    df = pd.DataFrame(manifest["profiles"]).merge(
        res[sims.str.startswith("libradtran")]
        [["label", "sim_lwdn_sfc", "sim_lwup_sfc"]], on="label")
    rrt = res[sims.str.startswith("rrtmg")]
    if len(rrt):
        df = df.merge(
            rrt[["label", "sim_lwdn_sfc", "sim_lwup_sfc"]].rename(
                columns={"sim_lwdn_sfc": "rrtmg_lwdn_sfc",
                         "sim_lwup_sfc": "rrtmg_lwup_sfc"}),
            on="label", how="left")
    return manifest, df


def stat_line(err: np.ndarray, e: np.ndarray, s: np.ndarray) -> str:
    return (f"bias={err.mean():+6.2f}  rmse={np.sqrt((err ** 2).mean()):6.2f}  "
            f"r={np.corrcoef(e, s)[0, 1]:+.3f}")


def main(argv: "list[str] | None" = None) -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from reanlib.plotstyle import apply_agu_style, panel_label

    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--month", type=int, required=True)
    parser.add_argument("--day", type=int, required=True)
    parser.add_argument("--hour", type=int, default=12,
                        help="end of the ERA5 accumulation window (default 12)")
    parser.add_argument("--prev-hour", type=int, default=None,
                        help="start-of-window snapshot (default: hour - 1)")
    parser.add_argument("--sky", default="clear", choices=["clear", "cloudy"])
    parser.add_argument("--config", default=None)
    args = parser.parse_args(argv)
    h1 = args.prev_hour if args.prev_hour is not None else args.hour - 1

    cfg = load_config(args.config)
    date = dt.date(args.year, args.month, args.day)
    man2, df2 = load_hour(cfg, date, args.hour, args.sky)      # end of window
    man1, df1 = load_hour(cfg, date, h1, args.sky)             # start of window
    df = df2.merge(df1, on="label", suffixes=("_h2", "_h1"))
    has_rrtmg = ("rrtmg_lwdn_sfc_h2" in df.columns
                 and "rrtmg_lwdn_sfc_h1" in df.columns)

    # ERA5 truth: the accumulation over [h1, h2] ending at h2
    era5 = {"lwdn": df["era5_lwdn_h2"].values, "lwup": df["era5_lwup_h2"].values}

    # far-IR tail estimates for libRadtran (RRTMG covers the full band):
    # LWup from each hour's skin temperature, LWdn from the (shared) ERA5 flux
    w1, w2 = np.array(man2["wavelength_range_nm"]) / 1000.0
    eps = man2["emissivity"]
    tails = {}
    for h in ("h1", "h2"):
        tails["lwup", h] = np.array(
            [(1 - planck_band_fraction(t, w1, w2)) * SIGMA * t ** 4 * eps
             for t in df["skt_" + h]])
    t_eff = (era5["lwdn"] / SIGMA) ** 0.25
    tails["lwdn", "h1"] = tails["lwdn", "h2"] = np.array(
        [(1 - planck_band_fraction(t, w1, w2)) * SIGMA * t ** 4 for t in t_eff])

    # estimators per component: {name: (lib+tail, rrtmg)}
    est = {}
    for kind in ("lwdn", "lwup"):
        lib2 = df[f"sim_{kind}_sfc_h2"].values + tails[kind, "h2"]
        lib1 = df[f"sim_{kind}_sfc_h1"].values + tails[kind, "h1"]
        rrt2 = df[f"rrtmg_{kind}_sfc_h2"].values if has_rrtmg else None
        rrt1 = df[f"rrtmg_{kind}_sfc_h1"].values if has_rrtmg else None
        est[kind] = {
            f"inst {args.hour:02d}Z": (lib2, rrt2),
            f"inst {h1:02d}Z": (lib1, rrt1),
            "average": ((lib1 + lib2) / 2.0,
                        (rrt1 + rrt2) / 2.0 if has_rrtmg else None),
        }

    # headline subset: pixels still cloud-free at the earlier hour.
    # (cc/lwp/iwp columns are only suffixed when both manifests record them —
    # manifests written before the state-time work lack these fields.)
    def h1col(name):
        return df[name + "_h1"] if name + "_h1" in df.columns else df[name]

    still = ((h1col("cc_max").values <= CLEAR_CC_MAX)
             & (h1col("lwp_g").values + h1col("iwp_g").values <= CLEAR_PATH_G))
    n_all, n_clear = len(df), int(still.sum())

    print(f"\nState-time test — ERA5 {h1:02d}-{args.hour:02d}Z accumulation vs "
          f"snapshot estimators, {date} ({args.sky} sky)")
    print(f"pixels: {n_all} total; {n_clear} still cloud-free at {h1:02d}Z "
          f"(headline stats), {n_all - n_clear} excluded (cloud in the window)")
    for kind in ("lwdn", "lwup"):
        e = era5[kind][still]
        print(f"  LW{'dn' if kind == 'lwdn' else 'up'}:")
        for name, (lib, rrt) in est[kind].items():
            s = lib[still]
            line = f"    {name:<9} libRadtran+tail: {stat_line(s - e, e, s)}"
            if has_rrtmg:
                g = rrt[still]
                line += f"   RRTMG: {stat_line(g - e, e, g)}"
            print(line)
    # the excluded pixels: ERA5 saw cloud inside the window, the clear-sky
    # snapshot cannot — expect ERA5 > sim in LWdn
    if n_all > n_clear:
        e = era5["lwdn"][~still]
        s = est["lwdn"][f"inst {args.hour:02d}Z"][0][~still]
        print(f"  excluded (cloudy at {h1:02d}Z) LWdn inst {args.hour:02d}Z: "
              f"{stat_line(s - e, e, s)}  (ERA5 > sim expected)")

    # cross-window check: window timing vs state identity. If the snapshot/
    # accumulation TIMING drove the mismatch, each snapshot would fit the
    # window it starts best; if instead one analysis STATE is off-trajectory,
    # that snapshot fits no window (and the other fits them all).
    def tail_dn(e):
        return np.array([(1 - planck_band_fraction(t, w1, w2)) * SIGMA * t ** 4
                         for t in (e / SIGMA) ** 0.25])

    lib1_raw = df["sim_lwdn_sfc_h1"].values
    lib2_raw = df["sim_lwdn_sfc_h2"].values
    rows = [
        (f"sim {h1:02d}Z vs {h1 - 1:02d}-{h1:02d}Z accum",
         lib1_raw, df["era5_lwdn_h1"].values, still),
        (f"sim {h1:02d}Z vs {h1:02d}-{args.hour:02d}Z accum",
         lib1_raw, era5["lwdn"], still),
        (f"sim {args.hour:02d}Z vs {h1:02d}-{args.hour:02d}Z accum",
         lib2_raw, era5["lwdn"], still),
    ]
    h3 = args.hour + 1
    try:
        man3 = json.loads(manifest_path(cfg, date, h3, args.sky).read_text())
        p3 = {p["label"]: p for p in man3["profiles"]}
        e3 = np.array([p3[l]["era5_lwdn"] for l in df["label"]])
        clear3 = np.array([(p3[l]["cc_max"] <= CLEAR_CC_MAX
                            and p3[l]["lwp_g"] + p3[l]["iwp_g"] <= CLEAR_PATH_G)
                           for l in df["label"]])
        rows += [
            (f"sim {h1:02d}Z vs {args.hour:02d}-{h3:02d}Z accum",
             lib1_raw, e3, still),
            (f"sim {args.hour:02d}Z vs {args.hour:02d}-{h3:02d}Z accum",
             lib2_raw, e3, clear3),
        ]
    except FileNotFoundError:
        print(f"  (no {h3:02d}Z manifest — prep it with --pixels-from for the "
              "next-window cross-check)")
    print("  cross-window (LWdn, libRadtran+tail, pixels clear at the "
          "snapshot hours):")
    for tag, sim, truth, mask in rows:
        s = sim[mask] + tail_dn(truth[mask])
        e = truth[mask]
        print(f"    {tag:<28} {stat_line(s - e, e, s)}  n={mask.sum()}")

    # ---------------------------------------------------------------- figure
    apply_agu_style()
    fig, ax4 = plt.subplots(2, 2, figsize=(9.8, 9.4), layout="constrained")
    axes = ax4.ravel()
    C_INST, C_AVG, C_PREV = "#999999", "#0072B2", "#D55E00"
    name2 = f"inst {args.hour:02d}Z"
    name1 = f"inst {h1:02d}Z"

    for i, kind in enumerate(("lwdn", "lwup")):
        ax = axes[i]
        e = era5[kind][still]
        s2 = est[kind][name2][0][still]
        sa = est[kind]["average"][0][still]
        lim = (min(e.min(), s2.min(), sa.min()) - 3,
               max(e.max(), s2.max(), sa.max()) + 3)
        ax.plot(lim, lim, color="#bbbbbb", lw=1, zorder=0)
        ax.scatter(e, s2, s=9, color=C_INST, alpha=0.35, edgecolors="none",
                   label=name2)
        ax.scatter(e, sa, s=9, color=C_AVG, alpha=0.35, edgecolors="none",
                   label=f"({h1:02d}Z+{args.hour:02d}Z)/2")
        txt = "\n".join(
            f"{lab}: bias {(s - e).mean():+.2f}, rmse "
            f"{np.sqrt(((s - e) ** 2).mean()):.2f}, r {np.corrcoef(e, s)[0, 1]:+.3f}"
            for lab, s in ((name2, s2), ("average", sa)))
        ax.text(0.04, 0.96, f"n = {len(e)}\n" + txt, transform=ax.transAxes,
                va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#cccccc"))
        arrow = "\\downarrow" if kind == "lwdn" else "\\uparrow"
        ax.set_xlabel(f"ERA5 LW${arrow}$ {h1:02d}–{args.hour:02d}Z accum. "
                      "(W m$^{-2}$)")
        ax.set_ylabel(f"libRadtran+tail LW${arrow}$ (W m$^{{-2}}$)")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect("equal")
        ax.legend(frameon=False, fontsize=8, loc="lower right")
        ax.grid(color="#e3e3e3", lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        panel_label(ax, "ab"[i], x=-0.14, y=1.06)

    # (c) LWdn bias distributions for the three estimators
    ax = axes[2]
    e = era5["lwdn"][still]
    errs = [(name1, est["lwdn"][name1][0][still] - e, C_PREV, False),
            (name2, est["lwdn"][name2][0][still] - e, C_INST, False),
            (f"({h1:02d}Z+{args.hour:02d}Z)/2",
             est["lwdn"]["average"][0][still] - e, C_AVG, True)]
    span = max(abs(np.concatenate([x for _, x, _, _ in errs])).max(), 1.0)
    bins = np.linspace(-span, span, 41)
    for lab, err, color, fill in errs:
        ax.hist(err, bins=bins, histtype="stepfilled" if fill else "step",
                color=color, alpha=0.45 if fill else 1.0, lw=1.2,
                label=f"{lab}  ({err.mean():+.1f} $\\pm$ {err.std():.1f})")
    ax.axvline(0, color="#555555", lw=0.8)
    ax.set_xlabel("LW$\\downarrow$ bias, sim $-$ ERA5 (W m$^{-2}$)")
    ax.set_ylabel("pixels")
    ax.legend(frameon=False, fontsize=8, title="bias $\\pm$ std (W m$^{-2}$)",
              title_fontsize=8)
    ax.grid(color="#e3e3e3", lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    panel_label(ax, "c", x=-0.14, y=1.06)

    # (d) does the within-window drift explain the snapshot error pixel-by-pixel?
    ax = axes[3]
    lib2 = est["lwdn"][name2][0][still]
    lib1 = est["lwdn"][name1][0][still]
    x = (lib1 - lib2) / 2.0                 # correction the average applies
    y = era5["lwdn"][still] - lib2          # gap the snapshot leaves
    lim = (min(x.min(), y.min()) - 1, max(x.max(), y.max()) + 1)
    ax.plot(lim, lim, color="#bbbbbb", lw=1, zorder=0,
            label="1:1 (ERA5 = window mean)")
    ax.plot(lim, (2 * lim[0], 2 * lim[1]), color="#D55E00", lw=1, ls="--",
            zorder=0, label=f"2:1 (ERA5 = inst {h1:02d}Z)")
    ax.scatter(x, y, s=9, color=C_AVG, alpha=0.35, edgecolors="none")
    slope = np.polyfit(x, y, 1)[0]
    ax.text(0.04, 0.96,
            f"r = {np.corrcoef(x, y)[0, 1]:+.3f}, slope = {slope:.2f}\n"
            "slope 1 $\\Rightarrow$ accumulation = window mean\n"
            f"slope 2 $\\Rightarrow$ accumulation = {h1:02d}Z snapshot",
            transform=ax.transAxes, va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#cccccc"))
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_xlabel(f"applied correction $({h1:02d}Z - {args.hour:02d}Z)/2$ "
                  "(W m$^{-2}$)")
    ax.set_ylabel(f"snapshot gap, ERA5 $-$ inst {args.hour:02d}Z (W m$^{-2}$)")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.grid(color="#e3e3e3", lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    panel_label(ax, "d", x=-0.14, y=1.06)

    fig.suptitle(
        f"State-time test, {args.sky} sky — ERA5 {h1:02d}–{args.hour:02d}Z "
        f"accumulation vs analysis snapshots, {date} "
        f"(n = {n_clear} pixels cloud-free at both hours; libRadtran + far-IR "
        "tail)", y=1.03)
    outdir = figures_dir(cfg)
    outdir.mkdir(parents=True, exist_ok=True)
    fpath = (outdir / f"lw_statetime_{date:%Y%m%d}T{args.hour:02d}Z"
             f"{'' if args.sky == 'clear' else '_' + args.sky}.png")
    fig.savefig(fpath, bbox_inches="tight")
    print(f"wrote {fpath}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
