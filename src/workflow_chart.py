#!/usr/bin/env python
"""Draw the pipeline workflow charts -> docs/workflow.png (+ _prefire.png).

Two charts: the core inversion + LW-closure pipeline, and a separate one for
the PREFIRE brightness-temperature / Jacobian stage. They live under docs/
(not figures/) so they can be committed while all generated figures stay
untracked.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reanlib.config import REPO_ROOT
from reanlib.plotstyle import apply_agu_style

C_EXT = "#ffffff"     # external services
C_DATA = "#ececec"    # local data stores
C_PROC = "#dbe9f6"    # pipeline scripts
C_FIG = "#e3f2e1"     # figure outputs
EDGE = "#555555"


def box(ax, x, y, w, h, text, fc, fontsize=7.2, bold_first=True, ec=EDGE,
        ls="-"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.35",
                                fc=fc, ec=ec, lw=1.0, ls=ls))
    lines = text.split("\n")
    if bold_first and len(lines) > 1:
        ax.text(x + w / 2, y + h - 0.4, lines[0], ha="center", va="top",
                fontsize=fontsize, fontweight="bold")
        ax.text(x + w / 2, y + h - 0.6 - 1.55, "\n".join(lines[1:]),
                ha="center", va="top", fontsize=fontsize - 0.6)
    else:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=fontsize,
                fontweight="bold" if bold_first else "normal")
    return (x, y, w, h)


def arrow(ax, a, b, label=None, ls="-", con="arc3,rad=0.0", lab_dy=0.5):
    xa = a[0] + a[2]
    ya = a[1] + a[3] / 2
    xb = b[0]
    yb = b[1] + b[3] / 2
    ax.add_patch(FancyArrowPatch((xa, ya), (xb, yb), arrowstyle="-|>",
                                 mutation_scale=9, color="#777777", lw=1.0,
                                 ls=ls, connectionstyle=con,
                                 shrinkA=2, shrinkB=2))
    if label:
        ax.text((xa + xb) / 2, (ya + yb) / 2 + lab_dy, label, ha="center",
                va="bottom", fontsize=6.0, color="#555555", style="italic")


def varrow(ax, a, b):
    """Vertical arrow from bottom of a to top of b (same column)."""
    xa = a[0] + a[2] / 2
    ax.add_patch(FancyArrowPatch((xa, a[1]), (b[0] + b[2] / 2, b[1] + b[3]),
                                 arrowstyle="-|>", mutation_scale=9,
                                 color="#777777", lw=1.0,
                                 shrinkA=2, shrinkB=2))


def toparc(ax, a, b, rad, x_frac=0.4, start_frac=0.5):
    """Arc from the top of box a over intervening boxes to the top of box b."""
    ax.add_patch(FancyArrowPatch((a[0] + a[2] * start_frac, a[1] + a[3]),
                                 (b[0] + b[2] * x_frac, b[1] + b[3]),
                                 arrowstyle="-|>", mutation_scale=9,
                                 color="#777777", lw=1.0,
                                 connectionstyle=f"arc3,rad={rad}",
                                 shrinkA=2, shrinkB=2))


def elbow(ax, pts, ls="-"):
    """Right-angle connector through the given waypoints, arrowhead at the end."""
    xs = [p[0] for p in pts[:-1]]
    ys = [p[1] for p in pts[:-1]]
    ax.plot(xs, ys, color="#777777", lw=1.0, ls=ls, solid_capstyle="round",
            zorder=1)
    ax.add_patch(FancyArrowPatch(pts[-2], pts[-1], arrowstyle="-|>",
                                 mutation_scale=9, color="#777777", lw=1.0,
                                 ls=ls, shrinkA=0, shrinkB=2))


def legend_row(ax, y, entries):
    for x, fc, lab in entries:
        ax.add_patch(FancyBboxPatch((x, y), 2.6, 2.2,
                                    boxstyle="round,pad=0.25", fc=fc, ec=EDGE,
                                    lw=0.8))
        ax.text(x + 3.4, y + 1.1, lab, va="center", fontsize=7)


def draw_main() -> "plt.Figure":
    fig, ax = plt.subplots(figsize=(13.0, 7.0))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 67)
    ax.axis("off")

    ax.text(50, 66.4, "ERA5 Arctic temperature-inversion & LW-closure pipeline",
            ha="center", va="top", fontsize=12, fontweight="bold")

    # column 1: external sources ------------------------------------------
    s1 = box(ax, 1, 50, 14.5, 7, "Copernicus CDS\nERA5 pressure levels\n+ single levels", C_EXT)
    s2 = box(ax, 1, 40, 14.5, 6, "Copernicus ADS\nCAMS EGG4\nCO2 / CH4", C_EXT)
    s3 = box(ax, 1, 20, 14.5, 7.5, "PANGAEA\nMOSAiC obs:\nJozef et al. 2023,\nMaturilli et al. 2021", C_EXT)

    # column 2: local data ------------------------------------------------
    d1 = box(ax, 21.5, 50, 15, 7, "data/<src>/YYYY/MM/DD/\n<src>_{plev,sfc}_*.nc\nvia era5_download.py or\nmerra2_download.py", C_DATA)
    d2 = box(ax, 21.5, 40, 15, 6, "data/cams/\nCO2/CH4 profiles via\ncams_download.py", C_DATA)
    d3 = box(ax, 21.5, 20, 15, 7.5, "data/mosaic/\nAtm_Properties.nc\n+ soundings/ (auto)", C_DATA)
    arrow(ax, s1, d1)
    arrow(ax, s2, d2)
    arrow(ax, s3, d3)

    # column 3: inversion analysis ----------------------------------------
    a1 = box(ax, 42, 50, 17, 7, "daily_inversion.py\nSBI scan ($\\Delta$T $\\geq$ 0.5 K),\nT$_{850}$$-$T$_{2m}$, T$_{925}$$-$T$_{1000}$", C_PROC)
    a2 = box(ax, 42, 43, 17, 4.5, "plot_profiles.py\nplot_maps.py", C_PROC)
    a3 = box(ax, 42, 36.5, 17, 4.5, "monthly_stats.py\nmonthly climatology", C_PROC)
    a4 = box(ax, 42, 30, 17, 4.5, "profile_analysis.py\nPCA, T$_{2m}$ correlations", C_PROC)
    a5 = box(ax, 42, 20, 17, 6, "mosaic_compare.py\n123 matched soundings,\nbias statistics", C_PROC)
    arrow(ax, d1, a1)
    varrow(ax, a1, a2)
    varrow(ax, a2, a3)
    varrow(ax, a3, a4)
    varrow(ax, a4, a5)
    arrow(ax, d3, a5)

    # column 4: radiative closure -----------------------------------------
    p1 = box(ax, 64, 50, 16, 6.5, "lrt_sim.py prep\npixel selection,\natm + cloud files", C_PROC)
    p2 = box(ax, 64, 41.5, 16, 5.5, "lrt_sim.py run\nuvspec thermal\n(libRadtran, DISORT)", C_PROC)
    p3 = box(ax, 64, 33.5, 16, 5.5, "rrtmg_sim.py\nRRTMG-LW (climlab),\nfull 3.08$-$1000 $\\mu$m", C_PROC)
    c1 = box(ax, 85, 41.5, 13.5, 5.5, "results CSV\n(simulator column)", C_DATA)
    c2 = box(ax, 85, 33.5, 13.5, 5.5, "lrt_sim.py\ncompare\n3-way vs reference", C_PROC)
    cs = box(ax, 64, 20, 16, 6.5, "case_study.py\nmosaic_flux.py\nMOSAiC pixels: walkthrough\n+ full-drift LW closure", C_PROC)
    fg = box(ax, 85, 20, 13.5, 6.5, "figures/*.png", C_FIG)
    # ERA5 data over the top; CAMS through the corridor between columns
    toparc(ax, d1, p1, rad=-0.18, x_frac=0.4)
    elbow(ax, [(36.5, 42.0), (61.6, 42.0), (61.6, 53.25), (64, 53.25)])
    varrow(ax, p1, p2)
    varrow(ax, p2, p3)
    arrow(ax, p2, c1)
    arrow(ax, p3, c1, con="arc3,rad=-0.2")
    varrow(ax, c1, c2)
    varrow(ax, c2, fg)
    arrow(ax, a5, cs)
    arrow(ax, cs, fg)

    # legend ---------------------------------------------------------------
    legend_row(ax, 12.5, ((1, C_EXT, "external service"),
                          (17, C_DATA, "local data"),
                          (31, C_PROC, "pipeline script"),
                          (46, C_FIG, "figure output")))
    ax.text(1, 16.8, "Stages run left to right; every stage is an idempotent CLI "
            "(python src/<script>.py --help). PREFIRE BT simulation: see "
            "docs/workflow_prefire.png.", fontsize=6.6, color="#555555")
    ax.set_ylim(11, 67)
    return fig


def draw_prefire() -> "plt.Figure":
    fig, ax = plt.subplots(figsize=(12.0, 6.2))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 60)
    ax.axis("off")

    ax.text(50, 59.4, "PREFIRE brightness-temperature simulation & Jacobians",
            ha="center", va="top", fontsize=12, fontweight="bold")

    # sources --------------------------------------------------------------
    s1 = box(ax, 1, 44, 15, 7, "NASA ASDC\nPREFIRE SAT1/2\n1B-RAD (R01)", C_EXT)
    s2 = box(ax, 1, 33, 15, 7, "Zenodo\nTIRS SRF v13\n(BT lookup + NEdR)", C_EXT)
    s3 = box(ax, 1, 20, 15, 7, "data/YYYY/MM/DD/\nERA5 columns\n(main pipeline)", C_DATA)

    # data ------------------------------------------------------------------
    d1 = box(ax, 21.5, 38.5, 16, 7.5, "data/prefire/\nYYYY/MM/ granules + srf/\nvia prefire_download.py", C_DATA)
    arrow(ax, s1, d1)
    arrow(ax, s2, d1)

    # prefire_bt.py subcommand stack -----------------------------------
    p1 = box(ax, 43, 44, 17, 8, "collocate\nfootprints $\\rightarrow$ (cell, hour)\ncolumns; pick clear +\novercast test set", C_PROC)
    p2 = box(ax, 43, 33.5, 17, 7, "prep\natm + cloud files,\nper-scene obs BT manifest", C_PROC)
    p3 = box(ax, 43, 23, 17, 7, "run\nuvspec thermal radiance,\nSRF $\\rightarrow$ channel BT", C_PROC)
    p4 = box(ax, 43, 12, 17, 7.5, "jacobian\nperturb skt, T(z), q(z),\ncloud, emissivity", C_PROC)
    ax.text(51.5, 53.6, "prefire_bt.py", ha="center", fontsize=7.4,
            fontweight="bold", style="italic", color="#333333")
    arrow(ax, d1, p1)
    elbow(ax, [(16, 23.5), (40.5, 23.5), (40.5, 48), (43, 48)])
    varrow(ax, p1, p2)
    varrow(ax, p2, p3)
    varrow(ax, p3, p4)

    # cross-check + outputs --------------------------------------------------
    r1 = box(ax, 66, 34, 15, 6.5, "rrtmg\n16-band flux-equiv.\nBT cross-check", C_PROC)
    o1 = box(ax, 66, 23.5, 15, 6.5, "results JSON\nchannel BT + spectra\n(sim vs obs)", C_DATA)
    o2 = box(ax, 66, 12, 15, 6.5, "jacobian_*.nc\nK matrices per column\n(+ NEdR for S$_e$)", C_DATA)
    arrow(ax, p2, r1)
    varrow(ax, r1, o1)
    arrow(ax, p3, o1)
    arrow(ax, p4, o2)

    f1 = box(ax, 86, 18, 13, 12, "figure\nBT spectra vs obs,\nK heatmaps\n$\\rightarrow$ figures/\nprefire_*.png", C_FIG)
    arrow(ax, o1, f1)
    arrow(ax, o2, f1)

    # planned follow-on ------------------------------------------------------
    fut = box(ax, 86, 5, 13, 8, "planned:\ncloud-property\nretrieval (OE) +\nEarthCARE validation", C_EXT, ls="--", ec="#999999")
    arrow(ax, o2, fut, ls="--")

    legend_row(ax, 2.5, ((1, C_EXT, "external service"),
                         (17, C_DATA, "local data"),
                         (31, C_PROC, "pipeline script"),
                         (46, C_FIG, "figure output")))
    ax.text(1, 7.0, "collocate/prep run on the ERA5 side; run/jacobian are "
            "libRadtran jobs (fine spectral grid on the cluster).",
            fontsize=6.6, color="#555555")
    ax.set_ylim(1.5, 60)
    return fig


def main() -> int:
    apply_agu_style()
    out_dir = REPO_ROOT / "docs"
    out_dir.mkdir(exist_ok=True)
    for fname, draw in (("workflow.png", draw_main),
                        ("workflow_prefire.png", draw_prefire)):
        fig = draw()
        out = out_dir / fname
        fig.savefig(out, bbox_inches="tight")
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
