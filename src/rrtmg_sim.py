#!/usr/bin/env python
"""RRTMG-LW cross-check of the libRadtran LW simulation (stage 7b).

Run under the ``era5`` conda env (needs climlab + climlab-rrtmg; installed via
conda-forge). Reads the manifest written by ``lrt_sim.py prep`` and reruns
every pixel through climlab's RRTMG_LW — the same radiation-scheme family ERA5
uses, covering the full 3.08-1000 um range, so no far-IR tail correction is
needed. Rows are appended to the shared results CSV with simulator
``rrtmg-lw ...``; ``lrt_sim.py compare`` (under er3t_env) then produces a
three-way libRadtran / RRTMG / ERA5 comparison.

Input identity with the libRadtran runs is guaranteed by parsing the very
profile files uvspec consumed: the 9-column atmosphere file (incl. the afglsw
upper-atmosphere splice and CAMS CO2), the CH4 mol_file, and the wc/ic 1D
cloud files (mirroring libRadtran's layer convention: a row's water content
fills the layer up to the level above it). Cloud optics are matched by
parameterization family: liquid Hu & Stamnes (liqflglw=1, reff 10 um) as
"wc_properties hu", ice Fu 1996 (iceflglw=3, dge = 1.0315 reff) as
"ic_properties fu". Known input differences vs the uvspec runs: N2O is a
0.32 ppm constant here vs libRadtran's built-in default profile, and RRTMG
takes layer means of the level profiles (both sub-W/m2 effects).

Examples:
    python src/rrtmg_sim.py --year 2020 --month 1 --day 1 --hour 12
    python src/rrtmg_sim.py --year 2020 --month 1 --day 1 --hour 12 --sky cloudy
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reanlib.config import load_config, source_label
from lrt_sim import M_DRY, M_H2O, REF_KEY_COMPAT, manifest_path, results_path

N2O_PPM = 0.32          # libRadtran default profile is ~this near the surface
DGE_FROM_REFF = 1.0315  # Fu 1996 generalized effective size from reff

_BOTTOM_AIR_T = {"value": np.nan}   # per-pixel t2m for the patched tlev


def patch_interface_temperature():
    """Bottom interface temperature = surface AIR temperature, not skin.

    climlab's interface_temperature() sets the lowest interface to Ts (skin
    temperature), which leaks the radiative skin temperature into the lowest
    air layer's emission — several W/m2 of LWdn error under strong surface
    inversions where skt is much colder than t2m. libRadtran (and ERA5's IFS)
    keep the air column at t2m and use skt only for the surface emission, so
    match that here; Ts still drives surface emission via tsfc.
    """
    from climlab.radiation.rrtm import utils as rrtm_utils
    from scipy.interpolate import interp1d

    def interface_temperature_air_bottom(Ts, Tatm, **kwargs):
        lev = Tatm.domain.axes["lev"].points
        bounds = Tatm.domain.axes["lev"].bounds
        t_mid = interp1d(lev, Tatm, axis=-1)(bounds[1:-1])
        t_toa = np.asarray(Tatm)[..., 0, np.newaxis]
        t_bot = np.full_like(t_toa, _BOTTOM_AIR_T["value"])
        return np.concatenate((t_toa, t_mid, t_bot), axis=-1)

    rrtm_utils.interface_temperature = interface_temperature_air_bottom


def q_from_vmr(x: np.ndarray) -> np.ndarray:
    """H2O volume mixing ratio -> specific humidity (kg/kg)."""
    w = x / (1.0 - x) * (M_H2O / M_DRY)
    return w / (1.0 + w)


def layer_mean(a: np.ndarray) -> np.ndarray:
    return 0.5 * (a[:-1] + a[1:])


def cloud_layer_paths(cloud_file: str, z_lev: np.ndarray) -> np.ndarray:
    """Per-layer in-cloud water path (g/m2) on the atmosphere layer grid.

    Mirrors libRadtran's 1D cloud-file convention (default since v1.4): the
    value on a row fills the layer from that altitude up to the row above;
    the top row's value is unused. z_lev is the TOA-first level grid of the
    atmosphere file; layer i sits between z_lev[i] and z_lev[i+1].
    """
    rows = np.atleast_2d(np.loadtxt(cloud_file, comments="#"))
    wp = np.zeros(z_lev.size - 1)
    for j in range(1, rows.shape[0]):
        z_hi, z_lo, wc = rows[j - 1, 0], rows[j, 0], rows[j, 1]
        if wc <= 0.0:
            continue
        i = int(np.argmin(np.abs(z_lev - z_hi)))
        wp[i] += wc * (z_hi - z_lo) * 1000.0
    return wp


def run_pixel(prof: dict, manifest: dict, climlab, Axis) -> dict:
    """One RRTMG_LW column on the exact grid/gases of the uvspec profile files."""
    rows = np.loadtxt(prof["atm_file"], comments="#")      # TOA-first
    z, p, T, air, o3, o2, h2o, co2, _no2 = rows.T
    ch4 = np.atleast_2d(np.loadtxt(prof["ch4_file"], comments="#"))[:, 1]

    state = climlab.column_state(lev=Axis(axis_type="lev", bounds=p))
    state.Ts[:] = prof["skt"]
    state.Tatm[:] = layer_mean(T)
    _BOTTOM_AIR_T["value"] = T[-1]      # surface-row air temperature (t2m)

    absorber_vmr = {
        "O3": layer_mean(o3 / air),
        "CO2": layer_mean(co2 / air),
        "CH4": layer_mean(ch4 / air),
        "O2": layer_mean(o2 / air),
        "N2O": N2O_PPM * 1e-6,
        "CCL4": 0.0, "CFC11": 0.0, "CFC12": 0.0, "CFC22": 0.0,
    }
    kwargs = {"icld": 0}
    clwp = cloud_layer_paths(prof["wc_file"], z) if prof.get("wc_file") \
        else np.zeros(p.size - 1)
    ciwp = cloud_layer_paths(prof["ic_file"], z) if prof.get("ic_file") \
        else np.zeros(p.size - 1)
    if (clwp + ciwp).any():
        # reff from the manifest (fixed values; see lrt_sim.py). Mind
        # that changing effective radius moves LWdn by up to ~20 W/m2 for
        # optically thin clouds but has almost no effect on thick,
        # emissivity-saturated ones (absorption ~ 1/reff until tau >> 1).
        reff = manifest["reff_um"]
        kwargs = {
            "icld": 1, "inflglw": 2, "liqflglw": 1, "iceflglw": 3,
            "cldfrac": ((clwp + ciwp) > 0).astype(float),
            "clwp": clwp, "ciwp": ciwp,
            "r_liq": np.full(clwp.size, reff["liquid"]),
            "r_ice": np.full(ciwp.size, DGE_FROM_REFF * reff["ice"]),
        }

    rad = climlab.radiation.RRTMG_LW(
        state=state,
        specific_humidity=q_from_vmr(layer_mean(h2o / air)),
        emissivity=manifest["emissivity"],
        absorber_vmr=absorber_vmr,
        **kwargs)
    rad.compute_diagnostics()
    return {
        "sim_lwdn_sfc": float(rad.LW_flux_down[-1]),
        "sim_lwup_sfc": float(rad.LW_flux_up[-1]),
        "sim_olr_toa": float(np.squeeze(rad.OLR)),
    }


def main(argv: "list[str] | None" = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--month", type=int, required=True)
    parser.add_argument("--day", type=int, required=True)
    parser.add_argument("--hour", type=int, required=True)
    parser.add_argument("--sky", default="clear", choices=["clear", "cloudy"])
    parser.add_argument("--source", choices=["era5", "merra2", "carra2"],
                        default=None,
                        help="data source (default from config.yaml)")
    parser.add_argument("--config", default=None)
    args = parser.parse_args(argv)

    try:
        import climlab
        from climlab.domain.axis import Axis
    except ImportError:
        sys.exit("climlab not importable — run under the `era5` conda env "
                 "(conda install -n era5 -c conda-forge climlab climlab-rrtmg)")
    import pandas as pd

    patch_interface_temperature()

    cfg = load_config(args.config, source=args.source)
    date = dt.date(args.year, args.month, args.day)
    mpath = manifest_path(cfg, date, args.hour, args.sky)
    if not mpath.exists():
        sys.exit(f"missing {mpath} — run lrt_sim.py prep first "
                 f"(with --sky {args.sky})")
    manifest = json.loads(mpath.read_text())
    profiles = manifest["profiles"]
    for p in profiles:      # pre-rename manifests used era5_lw* keys
        for old, new in REF_KEY_COMPAT.items():
            if old in p:
                p[new] = p.pop(old)
    simulator = f"rrtmg-lw (climlab {climlab.__version__})"

    print(f"RRTMG-LW on {len(profiles)} {manifest.get('sky', 'clear')}-sky "
          f"profiles — {manifest['snapshot']} (emissivity "
          f"{manifest['emissivity']}, full 3.08-1000 um band)")
    t0 = time.time()
    rows = []
    for k, prof in enumerate(profiles):
        out = run_pixel(prof, manifest, climlab, Axis)
        rows.append({"label": prof["label"], "simulator": simulator, **out})
        if len(profiles) <= 12:
            print(f"  {prof['label']:<10} LWdn {out['sim_lwdn_sfc']:7.2f}  "
                  f"LWup {out['sim_lwup_sfc']:7.2f}  "
                  f"OLR {out['sim_olr_toa']:7.2f} W/m2")
        elif (k + 1) % 100 == 0:
            print(f"  {k + 1}/{len(profiles)} columns "
                  f"({(time.time() - t0) / (k + 1) * 1000.0:.0f} ms/column) ...")

    new = pd.DataFrame(rows)
    ref_dn = np.array([p["ref_lwdn"] for p in profiles])
    ref_up = np.array([p["ref_lwup"] for p in profiles])
    print(f"done in {time.time() - t0:.1f} s — RRTMG bias vs {source_label(cfg)}: "
          f"LWdn {(new['sim_lwdn_sfc'] - ref_dn).mean():+.2f}, "
          f"LWup {(new['sim_lwup_sfc'] - ref_up).mean():+.2f} W/m2")

    rpath = results_path(cfg, date, args.hour, args.sky)
    if rpath.exists():
        old = pd.read_csv(rpath)
        old = old[~old["simulator"].astype(str).str.startswith("rrtmg")]
        new = pd.concat([old, new], ignore_index=True)
    new.to_csv(rpath, index=False)
    print(f"wrote {rpath} ({len(new)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
