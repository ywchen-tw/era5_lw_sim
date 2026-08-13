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
      `era5_*` → `rean_*`. Stages 7/7b/8 have since gained MERRA-2 support
      (see below); only 7c remains ERA5-only.
- [x] **CARRA-2 as a third data source (stages 1–6)** — `--source carra2`,
      `carra2_download.py` (CDS `reanalysis-pan-carra`, the SUB-daily entry;
      `-means` is aggregates only), `reanlib/io_carra2.py` (normalize on
      write: ERA5 names, `q` from relative humidity via Alduchov & Eskridge
      Magnus, cloud cover to fraction, x/y + CF grid mapping, `domain_mask`).
      CARRA-2 is *regional*, so the native 2.5 km polar-stereographic grid is
      kept (dims `y`/`x`, 2-D lat/lon) rather than regridded — stages 2–6 were
      generalized through the new `reanlib/grid.py` (`hdims`, `area_weights`
      with the (1+sin φ)² polar-stereographic weight, `GridIndex` KD-tree
      nearest-neighbour, `grid_template`, `horizontal_coords`,
      `projection_crs`) and `mapping.grid_kwargs`, which draws a projected
      grid in its own CRS so no date-line seam appears at the pole. ERA5 and
      MERRA-2 outputs verified numerically unchanged (Jan 2020: 61.6 % SBI
      frequency / 6.43 K, EOF1 73.6 % / EOF2 14.4 %, MOSAiC 78.9 % agreement,
      r = +0.27, biases −2.00 K / +315 m / +2.95 K — all as documented).

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
- [ ] **Averaging kernels + DOF from the existing K files** — no new
      simulations needed: per column compute
      A = (Kᵀ S_e⁻¹ K + S_a⁻¹)⁻¹ Kᵀ S_e⁻¹ K from the `jacobian_*.nc`
      matrices (channel NEdR already stored → S_e; S_a from a monthly
      reanalysis profile covariance), report DOF = trace(A) and AK-row
      figures per state kind. This is the natural first step of the OE
      item below and quantifies what the TIRS channels can actually
      resolve per sky class. Note the forward sim itself needs no AK —
      the spectral weighting is fully handled by the per-scene SRF
      convolution + blackbody-lookup BT, and the vertical weighting
      functions ARE the K_T rows (closure check: K_skt + ΣK_T ≈ 1).
- [ ] **Cloud-property retrieval (OE)** — optimal estimation on the K
      netCDFs (channel NEdR is already stored for S_e); start with
      ln IWP / r_eff / CTH on overcast columns; averaging kernels + DOF
      figures.
- [ ] **EarthCARE validation of retrieved cloud properties** — collocate
      PREFIRE footprints with EarthCARE cloud products (MSI M-CLD,
      CPR/ATLID synergy ACM-CAP). Prerequisite: ESA EO account.
      When comparing against any retrieved L2 product (PREFIRE 2B-ATM or
      EarthCARE), apply the product's averaging kernels to the reanalysis
      profile first (Rodgers & Connor 2003: x̂ = x_a + A(x_model − x_a))
      so vertical-resolution smoothing isn't misread as state bias.

## Backlog / ideas

- [ ] **CARRA-2 stage 7 (LW closure)** — blocked on two source gaps, both
      documented in the README: CARRA-2 publishes **no ozone** (needs a
      climatological profile) and its profile top is **50 hPa**, so the
      standard-atmosphere splice starts there instead of at 1 hPa. Its
      surface radiation is also **forecast-stream only**
      (`product_type: forecast` + `leadtime_hour`, accumulated from the cycle
      start), so the reference flux needs a second request and leadtime
      differencing rather than ERA5's single sfc file. At 2.5 km a
      lower-resolution pixel sample would be the sane starting point.
- [ ] **CARRA-2 vs ERA5 at matched resolution** — the obvious science use of
      the new source: does a 2.5 km regional model produce systematically
      stronger/shallower SBIs than 0.25° ERA5 over the same MOSAiC soundings?
      Stage 6 already runs on both; the comparison needs one shared month
      downloaded (mind the volume note in the README) and a joint figure.
- [ ] **CARRA-2 first real-data validation** — the stages were verified on a
      synthetic delivery through the real normalizer plus the unchanged
      ERA5/MERRA-2 regressions; the first genuine CDS delivery should be
      checked for coordinate spelling (the netCDF converter is documented as
      experimental — GRIB is the native format) and for whether below-ground
      pressure levels arrive as fill values or extrapolated.

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
      Chou–Suarez LW scheme vs RRTMG-family physics and the 42-level
      pressure-level truncation of the model state (no stage-7c-style
      hourly state-time test planned for MERRA-2). LW↑ closes to −0.7 W/m².
- [x] **MERRA-2 stage 7 cloudy (overcast)** — screened with the rad file's
      CLDTOT (1-h means bracketing the instant) ≥ 0.99 + condensate;
      per-level cloud *fraction* (partial-cloud work, `M2T3NPCLD`) remains
      future.
- [x] **MERRA-2 stage 8 (PREFIRE)** — `prefire_bt.py --source merra2`; the
      collocation snap is now the source state cadence (6 h ERA5 / 3 h
      MERRA-2, `--cadence` overrides), sky classification falls back to the
      stage-7 CLDTOT/condensate screens (no per-level cc in M2I3NPASM), and
      `--ic-properties baum` substitutes for the locally-absent yang2013
      tables. 2025-01-01 SAT1 coarse: 15Z clear columns close to +1.0 K
      bias / 2.3 K rmse vs PREFIRE (12Z clear +6.2 K — synoptic time);
      overcast −7…−11 K (cloud placement).
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
- Data source: `--source {era5,merra2,carra2}` on stages 2–6 (default from
  config.yaml `source:`); each downloader is pinned to its own source.
  MERRA-2 needs Earthdata credentials in `~/.netrc`; ERA5 and CARRA-2 need
  `~/.cdsapirc` (and each CDS dataset's licence accepted once, separately).
- Never assume 1-D `latitude`/`longitude`: CARRA-2 is on a projected `y`/`x`
  grid. Horizontal reductions, nearest-cell lookups and map panels go through
  `reanlib/grid.py` and `mapping.grid_kwargs`.
- `reptran fine`/`medium` never on the local Mac — CURC only.
- Commits are made locally on request; pushing is done by the user.
