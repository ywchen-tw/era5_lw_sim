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

- [x] **Stage 6b — per-level profile comparison** (`mosaic_profiles.py`,
      `reanlib/humidity.py`): T, q and RH at 14 levels common to all three
      sources against the same 123 soundings, with RH recomputed from q for
      every source *including* the observations so saturation conventions
      cannot masquerade as moisture differences. Collocation is per level at
      the balloon's drifted position, and the bias panels carry ±1σ bands
      (see the two backlog items marked done below). Result: the ~+3 K
      warm-surface bias of both global reanalyses is confined to the 2 m
      diagnostic — at 1000 hPa ERA5 is +1.01 K and MERRA-2 already −0.51 K,
      and from 850 hPa upward all three are within ±0.4 K — which is why both
      underestimate SBI strength while CARRA-2 (cold at 2 m and 1000 hPa)
      overestimates it. Humidity ranks the sources in the reverse order of
      resolution (q rmse over 1000–700 hPa: 0.07–0.10 / 0.09–0.16 /
      0.12–0.21 g/kg).
      NOTE the soundings were GTS-transmitted and assimilated, so agreement
      aloft is not independent validation; cross-source comparison is on
      firmer ground than absolute skill.

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

- [x] **Collocate each level at its own balloon position** (stage 6b) — each
      report level is now matched at the balloon's drifted position, taken
      from the level-2 per-sample Latitude/Longitude and interpolated on 3-D
      unit vectors in log-p (pole- and date-line-safe); one analysis time per
      sounding is kept (300 hPa is reached ~28 min after launch, small vs the
      6 h cadence), and `match_km` is per (sounding, level). Match distance is
      now bounded by half a cell at EVERY level — CARRA-2 max 1.7 km even at
      300 hPa, where launch-position matching could be off by up to 43 km.
      Outcome: no bias or rmse moved by more than 0.04 K / 0.01 g/kg /
      0.8 % RH (mostly a small rmse reduction aloft, the expected sign), so
      the launch-position shortcut was adequate after all — and CARRA-2's
      moist bias aloft survives proper collocation unchanged (+37 % RH at
      300 hPa), eliminating drift as its explanation. Stage 6 (`mosaic_compare.py`)
      still matches at the launch position, which is fine: its metrics live
      below ~800 hPa where drift ≤ 3 km.

- [x] **Spread on the stage-6b profile figure** — bias panels now carry a
      shaded ±1σ band of the per-sounding differences (σ from the same pairs;
      rmse² = bias² + σ²), with the band explained in the figure legend. The
      bands make the moist-bias story legible at a glance: CARRA-2's upper-
      level RH band sits entirely on the moist side of zero (consistent
      offset) while ERA5/MERRA-2 bands straddle it (scatter).

- [x] **CARRA-2 moist bias aloft — RESOLVED: saturation-convention artifact,
      not model moisture.** Settled against the published `r` field read
      straight from the retained raw chunks at the same drifted-balloon
      matched cells as stage 6b. Three results: (1) ingest validated end to
      end — published r equals the stage-6b RH recomputed from our derived q
      to 0.0001 % RH, so the +37 % at 300 hPa is a property of the delivered
      field under an over-water reading; (2) read as over-ICE, the upper-air
      bias collapses (600 hPa +22.2 → +1.1 %, 500 +28.6 → +0.04 %, 300
      +36.7 → +0.15 %) while below 800 hPa the over-water reading stays the
      good one (850 hPa: +2.0 % water vs −12.2 % ice); (3) binning
      r_published/RH_sonde by CARRA temperature tracks the theoretical
      e_w/e_i(T) curve for T ≲ −25 °C (1.56 vs 1.64 at −52 °C, 1.45 vs 1.50
      at −42 °C) and sits near 1.09 for T ≳ −20 °C. So CARRA-2's
      pressure-level r follows the model's temperature-dependent saturation
      (ice-like in cold air), unlike its 2 m RH which the CARRA docs define
      over water; the docs are silent on the plev convention. The earlier
      "ice rejected" test failed because it applied ice at ALL temperatures,
      which wrecks the (water-correct) moist lower levels that dominate the
      column q. CONSEQUENCES: `rh_over: water` stays at ingest (the
      empirical transition is under-constrained between −30 and −20 °C, and
      an invented ramp would bake a guess into the data); stage-6b upper-air
      CARRA-2 q/RH biases are conversion artifacts, not model moisture; the
      lower-troposphere humidity ranking (where the water reading is right
      and nearly all vapour lives) is unaffected. DECISION OPEN: implement an
      empirical ice-below-243 K / water-above-253 K ramp in the ingest and
      re-normalize the month from the retained chunks, if upper-air CARRA-2
      humidity is ever used quantitatively.

- [ ] **CARRA-2 stage 7 (LW closure)** — blocked on two source gaps, both
      documented in the README: CARRA-2 publishes **no ozone** (needs a
      climatological profile) and its profile top is **50 hPa**, so the
      standard-atmosphere splice starts there instead of at 1 hPa. Its
      surface radiation is also **forecast-stream only**
      (`product_type: forecast` + `leadtime_hour`, accumulated from the cycle
      start), so the reference flux needs a second request and leadtime
      differencing rather than ERA5's single sfc file. At 2.5 km a
      lower-resolution pixel sample would be the sane starting point.
      The forecast request has now been TESTED and works: `product_type:
      forecast` + `leadtime_hour: [1, 2]` returns `str` and `strd` — ERA5's
      own short names — dimensioned (step, y, x) with `time` the cycle start
      (12 UTC) and `valid_time` 13:00/14:00, so an hourly flux is the
      difference of consecutive steps. Also confirmed on an 80-90N probe:
      CARRA-2 EXTRAPOLATES below-ground pressure levels rather than filling
      them (244,300 of 624,693 band cells have sp < 1000 hPa and every one
      has a finite t(1000)), so it behaves like ERA5, not MERRA-2.
- [x] **Three-source comparison over matched DOMAINS** — the
      sounding-by-sounding comparison was already like-for-like; the
      matched-domain *climatology* now exists too, via `monthly_stats.py
      --area 90 -180 85 180 --hours 0 6 12 18` (see README tables). On
      identical 85-90N cells and hours: SBI frequency 60.7 / 50.2 / 92.0 %
      (ERA5 / MERRA-2 / CARRA-2), conditional strength 6.22 / 6.04 / 9.78 K,
      T850-T2m +3.81 / +2.47 / +7.89 K — CARRA-2's stronger-inversion
      character survives domain matching intact, consistent with the
      sounding result (over-detects, overestimates) while both global
      sources underestimate. The CARRA-2 run under `--area` reproduces its
      full-domain numbers exactly (mask identity check). NOTE the raw
      retained chunks were verified 85-90N-masked (0 finite cells below
      85N), so extending CARRA-2 southward still needs a re-download.
- [x] **Skin-referenced inversion strength `dt_850_skt`** — T(850) − T(skin)
      added to the fixed-level metrics (daily files recomputed for Jan 2020,
      all existing variables verified bit-identical; monthly + matched-domain
      aggregations and the maps figure carry it). Motivation: the ~+3 K
      warm-surface bias lives in the 2 m *diagnostic* (stage 6b), so the 2 m-
      and skin-referenced strengths bracket the surface-coupling uncertainty.
      Result (85-90N): the skin−2m offset differs in SIGN between the global
      sources — ERA5 skin 0.3 K colder than its 2 m, MERRA-2 skin 0.6 K
      warmer, CARRA-2 within 0.03 K — so the cross-source spread in
      surface-referenced strength widens from 5.4 to 6.1 K with the skin.
      No MOSAiC scoring for this metric: the radiosondes measure air, not
      skin (the sounding files have no skin-temperature counterpart).

- [ ] **ERA5-MERRA-2 SBI detection gap is a latitude story** — on their
      shared 80-90N domain they differ by 17.8 points in SBI frequency
      (61.6 vs 43.8 %) despite near-identical surface biases. Band-split
      shows ERA5 nearly flat with latitude (61.8 % at 80-85N, 60.7 % at
      85-90N) while MERRA-2 climbs poleward (41.2 % → 50.2 %): the
      disagreement is 20.6 points in the outer band — which contains the
      Atlantic-sector ice edge — and half that over the central pack. So the
      candidate explanations are the detection criterion interacting with
      vertical level spacing (37 vs 42 levels) AND the surface state near
      open water / thinner ice; a map of the frequency *difference* field
      would separate sea-ice-margin structure from uniform criterion effects.
- [x] **Analysis-stage `--area` for cross-source domain means** —
      `monthly_stats.py --area N W S E` (CDS order) masks the statistics to
      a sub-box at analysis time via the new `grid.box_mask` /
      `grid.area_tag`; outputs are tagged (`*_85-90N.nc`, figures likewise)
      so full-domain files are never overwritten, and the restricted mask is
      stored as the file's own `domain_mask` so downstream readers reproduce
      the same domain. A `--hours` filter landed with it, which also defuses
      the stage-7c footgun (extra 11Z/13Z hours in a daily file skewing a
      re-aggregation). `profile_analysis` (PCA) does not have the knob yet —
      add it there if a matched-domain EOF comparison is ever needed.
- [ ] **Widen CARRA-2 to 80-90N if the spatial story matters** — the current
      85-90N domain is all pack ice and drops N Greenland / Ellesmere /
      Svalbard / Franz Josef Land / Severnaya Zemlya, i.e. exactly the
      terrain where 2.5 km has most reason to beat 0.25 deg. 4x the volume
      (41 GB/month vs 10 GB). Trimming `plev_variables` to temperature +
      relative_humidity would pay for most of that (stages 1-6 read no other
      profile field), at the cost of re-downloading if stage 7 is ever
      ported — though that is blocked on ozone regardless.
- [x] **CARRA-2 first real-data validation** — January 2020 downloaded and
      run end to end through stages 1–6. The delivery exposed five bugs the
      synthetic test could not (it built its grid from the same assumptions
      the code made): the cfgrib `time`/`valid_time` collision, a guessed
      projection (true lon_0 = −30, recovered by requiring a regular grid),
      clipping that had to precede the q conversion or materialize 63 GB,
      plev/sfc coordinate drift silently NaN-ing every fixed-level metric via
      xarray alignment, and a CRS built without the earth radius. All fixed
      and committed; ERA5/MERRA-2 regressions unchanged throughout.
      Key delivery facts: the CDS honours `area` by MASKING a full
      2869×2869 canvas, not cropping it (155,588 populated cells for
      85–90°N); the netCDF converter yields `t`/`r` GRIB short names; and
      below-ground pressure levels arrive populated, not as fill values.

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
