# Plan & TODO

Status ledger for the pipeline. Stages run left to right in
`docs/workflow.png` / `docs/workflow_prefire.png`; details and usage live in
the README. Update this file when a stage lands or a plan changes.

## Done

- [x] **Stages 1–6** — ERA5 download (plev + sfc, idempotent), SBI +
      fixed-level inversion metrics, profile/map figures, monthly
      climatology, profile PCA, MOSAiC comparison (incl. matched-metric
      obs counterparts from the level-2 soundings).
- [x] **Stage 7 — LW flux closure** (`era5_lrt_sim.py`): libRadtran thermal
      fluxes vs ERA5 strd/str; n=500 clear + cloudy statistics.
      Clear-sky LW↓ mismatch diagnosed as forecast-accumulation vs analysis
      state-time, not radiative transfer.
- [x] **Stage 7b — RRTMG-LW cross-check** (`era5_rrtmg_sim.py`): validates
      the far-IR tail estimate (RRTMG − (lib+tail) = +0.02 ± 0.07 W/m²);
      bottom-interface skin-temperature patch.
- [x] **MOSAiC case study + full-drift closure** (`era5_case_study.py`,
      `era5_mosaic_flux.py`): 123 matched columns, Jan 2020; LW↑ warm-skin
      bias +11.5 W/m² shared by sims and ERA5.
- [x] **Stage 8 — PREFIRE BT simulation + Jacobians**
      (`era5_prefire_download.py`, `era5_prefire_bt.py`): collocation
      (2025-01-01 SAT1 test set: 3 clear + 3 overcast), channel BT on the
      mission SRF/blackbody-lookup scale, RRTMG band cross-check,
      finite-difference K matrices (skt, T/q per level, cloud, emissivity),
      cotscan BT-vs-COT sensitivity figure.

## Next (in rough priority order)

- [ ] **CURC fine-grid stage-8 run** — local `reptran coarse` (~15 cm⁻¹)
      under-resolves far-IR channels (3–9 cm⁻¹ wide); science-grade BT and
      K need `reptran fine` on CURC (`slurm/curc_prefire_bt.sh`; NEVER
      fine/medium locally — OOM). Sync `derived/YYYY/MM/prefire_bt/` first.
- [ ] **Extend PREFIRE collocation** — more January 2025 days (ERA5
      download for the extra days; optionally hourly `--hours` at overpass
      times to shrink the ≤3 h state-time offset), then statistics over
      many clear columns instead of 3.
- [ ] **SAT2 / TIRS2 set** — rerun collocate→figure with `--sat 2`
      (second wavelength registration); compare the two instruments on
      shared scenes.
- [ ] **Cloud-property retrieval (OE)** — optimal estimation on the K
      netCDFs (channel NEdR is already stored for S_e); start with
      ln IWP / r_eff / CTH on overcast columns; averaging kernels + DOF
      figures.
- [ ] **EarthCARE validation of retrieved cloud properties** — collocate
      PREFIRE footprints with EarthCARE cloud products (MSI M-CLD,
      CPR/ATLID synergy ACM-CAP). Prerequisite: ESA EO account.
- [ ] **Hourly ERA5 (11+12 UTC) state-time test** — definitive check of
      the stage-7 clear-sky LW↓ diagnosis: simulate the 11–12Z
      accumulation window instead of the 12Z instant.

## Backlog / ideas

- [ ] Spectral far-IR snow/ice surface emissivity for stage 8 (currently
      constant ε = 0.99; the far-IR dirty-window K_skt gradient is where it
      matters — this is PREFIRE's own science target).
- [ ] Post-2020 greenhouse-gas source: CAMS EGG4 ends Dec 2020; stage 8
      uses flat 2025 constants (424 ppm CO₂ / 1.93 ppm CH₄). Consider CAMS
      operational GHG forecasts or NOAA marine-boundary-layer values.
- [ ] Partial-cloud BT columns — currently excluded from the stage-8 test
      set (BT does not blend linearly over a broken scene); needs
      radiance-space blending across independent clear/cloudy runs.
- [ ] IFS-style diagnosed r_eff (Martin et al. liquid, Sun–Rikus ice) as an
      alternative to fixed 10/25 µm in the flux stages.
- [ ] In-cloud condensate scaling (clwc/cc) for partial-cloud flux columns.
- [ ] TOA comparison: download `ttr` and compare simulated OLR against
      ERA5 at TOA (stage 7 currently closes at the surface only).
- [ ] RRTMG (and BT sims) on ERA5's 137 model levels instead of 37 pressure
      levels.
- [ ] Elevated-inversion (EI) metric alongside the surface-based SBI scan.
- [ ] Add `geopotential` to `download.plev_variables` for exact heights
      (SBI depth is currently hypsometric).

## Standing constraints

- Run stages under the right conda env (era5 vs er3t_env; see README).
- `reptran fine`/`medium` never on the local Mac — CURC only.
- Commits are made locally on request; pushing is done by the user.
