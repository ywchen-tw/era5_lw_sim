# Plan & TODO

Status ledger for the pipeline. Stages run left to right in
`docs/workflow.png` / `docs/workflow_prefire.png`; details and usage live in
the README. Update this file when a stage lands or a plan changes.

## Done

- [x] **Stages 1–6** — ERA5 download (plev + sfc, idempotent), SBI +
      fixed-level inversion metrics, profile/map figures, monthly
      climatology, profile PCA, MOSAiC comparison (incl. matched-metric
      obs counterparts from the level-2 soundings).
- [x] **Stage 7 — LW flux closure** (`lrt_sim.py`): libRadtran thermal
      fluxes vs ERA5 strd/str; n=500 clear + cloudy statistics.
      Clear-sky LW↓ mismatch diagnosed as forecast-accumulation vs analysis
      state-time, not radiative transfer.
- [x] **Stage 7b — RRTMG-LW cross-check** (`rrtmg_sim.py`): validates
      the far-IR tail estimate (RRTMG − (lib+tail) = +0.02 ± 0.07 W/m²);
      bottom-interface skin-temperature patch.
- [x] **MOSAiC case study + full-drift closure** (`case_study.py`,
      `mosaic_flux.py`): 123 matched columns, Jan 2020; LW↑ warm-skin
      bias +11.5 W/m² shared by sims and ERA5.
- [x] **Stage 8 — PREFIRE BT simulation + Jacobians**
      (`prefire_download.py`, `prefire_bt.py`): collocation
      (2025-01-01 SAT1 test set: 3 clear + 3 overcast), channel BT on the
      mission SRF/blackbody-lookup scale, RRTMG band cross-check,
      finite-difference K matrices (skt, T/q per level, cloud, emissivity),
      cotscan BT-vs-COT sensitivity figure.
- [x] **MERRA-2 as a second data source (stages 1–6)** — `source` option
      (config `source:` + `--source` on stages 2–6), per-source trees
      (`data/<source>/…`, `derived/<source>/…`, `figures/<source>/`),
      `merra2_download.py` (earthaccess + GES DISC OPeNDAP subsetting,
      M2I3NPASM 3-hourly instantaneous plev + M2I1NXASM hourly instantaneous
      sfc, normalize-on-write to ERA5 conventions via
      `reanlib/io_merra2.py`). Renamed source-agnostic code: `era5lib` →
      `reanlib`, stage scripts dropped the `era5_` prefix
      (`daily_inversion.py`, `monthly_stats.py`, …); MOSAiC pairs variables
      `era5_*` → `rean_*`. RT stages (7/7b/7c/8) remain ERA5-only.
- [x] **Stage 7c — hourly state-time test** (`statetime_test.py`,
      `lrt_sim.py prep --pixels-from`): 11Z + 13Z added for
      2020-01-01, the 12Z 500-pixel set re-simulated at 11Z. Verdict: the
      11Z analysis matches ALL adjacent hourly accumulations (rmse ~3 W/m²,
      r ~0.86–0.91), the (11+12)/2 average only halves the 12Z error, and
      the 12Z analysis fits no window (r ~0.15 even for 12–13Z, which it
      starts) — clear-sky LW↓ scatter is the synoptic-time (12Z) analysis
      increment displacing the state off the flux-producing trajectory,
      not accumulation-window timing (that is only ~3 W/m² rmse). LW↑
      insensitive (0.5 W/m²) throughout.

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

## Backlog / ideas

- [x] **MERRA-2 stage 7, clear-sky** — `lrt_sim.py`/`rrtmg_sim.py --source
      merra2`; `M2T1NXRAD` fetched as the `rad` dataset (1-h means stamped
      HH:30; the two windows bracketing the instant are averaged;
      LW↓ = LWGAB/EMIS, LW↑ = LWGEM + (1−EMIS)·LW↓). The three duplicated
      flux-conversion sites now go through
      `reanlib/fluxes.load_surface_lw()`; manifest keys renamed
      `era5_lw*` → `ref_lw*` (readers accept both).
- [ ] **MERRA-2 clear-sky LW↓ offset** — first run (2020-01-01 12Z, 5
      pixels): sim ≈ +6–7 W/m² above MERRA-2's flux (libRadtran and RRTMG
      agree; LWGABCLR ≈ all-sky there, so not cloud). Candidates: GEOS
      Chou–Suarez LW scheme vs RRTMG-family physics, IAU
      analysis-vs-trajectory state (ERA5 stage-7c analogue; MERRA-2 plev is
      3-hourly, so the 11Z trick becomes 09/12/15Z), and the 42-level
      pressure-level truncation of the model state. LW↑ closes to −0.7 W/m².
- [ ] **MERRA-2 stage 7 cloudy** — needs 3-D cloud fraction, which only
      exists time-averaged (`M2T3NPCLD`, stamps 01:30/04:30/…) — new
      state-time semantics; the stage-7c framing changes shape
      (instantaneous state at HH vs mean flux centered HH:30).
- [ ] **MERRA-2 stage 8 (PREFIRE)** — generalize the hard-coded 6-h snap in
      `prefire_bt.py` collocation to the source cadence; MERRA-2's 3-hourly
      plev state would halve the ≤3 h state-time offset listed above.
- [ ] Generalize the stage-7c finding: repeat the state-time test on more
      days and around other hours — especially 00Z (the other radiosonde
      synoptic time; expect the same off-trajectory jump) and the 18–19Z
      window (fluxes valid 19Z come from the fresh 18Z forecast, step 1) —
      to separate the synoptic-increment interpretation from a fortuitous
      trajectory crossing at 11Z on 2020-01-01.
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
- Data source: `--source {era5,merra2}` on stages 2–6 (default from
  config.yaml `source:`); each downloader is pinned to its own source.
  MERRA-2 needs Earthdata credentials in `~/.netrc`; ERA5 needs
  `~/.cdsapirc`.
- `reptran fine`/`medium` never on the local Mac — CURC only.
- Commits are made locally on request; pushing is done by the user.
