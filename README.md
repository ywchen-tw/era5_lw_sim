# ERA5 Arctic Temperature-Inversion Pipeline

Automated pipeline that downloads ERA5 reanalysis data over the Arctic
(80–90°N by default), computes temperature-inversion strength with three
metrics, aggregates monthly climatologies, analyzes profile variability (PCA)
and surface-temperature relationships, validates against MOSAiC radiosondes,
and simulates clear-sky broadband LW fluxes with libRadtran for comparison
with ERA5 radiation — with AGU-style figures at every stage.

## Workflow at a glance

![Pipeline workflow](docs/workflow.png)

Regenerate with `python src/era5_workflow_chart.py`.

## Repository layout

```
era5_analysis/
├── config.yaml                  # user-editable defaults (area, variables, SBI params, paths)
├── environment.yml              # conda env "era5"
├── src/
│   ├── era5_download.py         # stage 1: CDS downloads (parallel, idempotent)
│   ├── era5_inversion.py        # stage 2: daily inversion metrics
│   ├── era5_plot_profiles.py    # stage 3a: profile illustration figures
│   ├── era5_plot_maps.py        # stage 3b: polar snapshot maps
│   ├── era5_monthly.py          # stage 4: monthly statistics + figures
│   ├── era5_profile_analysis.py # stage 5: profile PCA, surface-T, correlations
│   ├── era5_mosaic_compare.py   # stage 6: ERA5 vs MOSAiC radiosondes
│   ├── era5_cams_download.py    # CAMS EGG4 CO2/CH4 profiles (for stage 7)
│   ├── era5_lrt_sim.py          # stage 7: libRadtran LW fluxes vs ERA5 (er3t_env!)
│   ├── era5_rrtmg_sim.py        # stage 7b: RRTMG-LW cross-check via climlab (era5 env)
│   ├── era5_case_study.py       # MOSAiC clear/cloudy single-pixel walkthrough figures
│   ├── era5_mosaic_flux.py      # LW simulation at every MOSAiC-matched column
│   └── era5lib/                 # shared code: config, CDS I/O, science, maps, style
├── slurm/                       # CURC job templates (stage 7)
├── data/YYYY/MM/DD/             # raw ERA5: era5_{plev,sfc}_YYYYMMDD.nc
├── data/mosaic/                 # MOSAiC observations (PANGAEA download)
├── data/cams/                   # CAMS EGG4 greenhouse-gas profiles
├── derived/YYYY/MM/DD/          # daily metrics + lw_sim/ (profiles, manifest, results)
├── derived/YYYY/MM/             # monthly products (stats, PCA, MOSAiC pairs)
└── figures/
```

## Setup

```bash
conda env create -f environment.yml
conda activate era5
```

CDS API credentials are required for downloading. Create `~/.cdsapirc`:

```
url: https://cds.climate.copernicus.eu/api
key: <your-personal-access-token>
```

The token is shown at <https://cds.climate.copernicus.eu/profile> (see
<https://cds.climate.copernicus.eu/how-to-api>). You must also accept the
licence terms on each ERA5 dataset's download page once, or requests are
rejected.

## Usage

```bash
# 1. download pressure-level + single-level data (one request per day per dataset)
python src/era5_download.py --year 2020 --month 1 --days 1-31 --jobs 6

# 2. compute daily inversion metrics
python src/era5_inversion.py --year 2020 --month 1 --days 1-31 --check

# 3a. profile illustrations for one snapshot (strongest / median / weakest points)
python src/era5_plot_profiles.py --year 2020 --month 1 --day 1 --hour 12

# 3b. polar-stereographic inversion-strength maps for one snapshot
python src/era5_plot_maps.py --year 2020 --month 1 --day 1 --hour 12

# 4. monthly aggregation: stats netCDF + map/distribution/time-series figures
python src/era5_monthly.py --year 2020 --month 1

# 5. profile PCA, surface-T statistics, strength correlations (maps + figures)
python src/era5_profile_analysis.py --year 2020 --month 1

# 6. sounding-by-sounding comparison against MOSAiC radiosondes
python src/era5_mosaic_compare.py --year 2020 --month 1
```

- All stages are idempotent: existing complete files are skipped
  (`--force`/`--overwrite` to redo). Stable defaults (area, variables, hours,
  SBI parameters, paths) live in `config.yaml`; per-run selection is CLI
  arguments.
- `era5_download.py --jobs N` (default 4) runs CDS requests concurrently —
  queue waits overlap, so month-scale downloads speed up several-fold.
  `--dry-run` prints the CDS requests without downloading. Requests mixing
  instantaneous and accumulated variables are delivered by the CDS as zips of
  per-stream netCDFs; the script merges them into one netCDF automatically.
- Stage 6 needs `data/mosaic/MOSAiC_Atm_Properties.nc`; fetch it once with
  `curl -L -o data/mosaic/MOSAiC_Atm_Properties.nc
  https://download.pangaea.de/dataset/957760/files/MOSAiC_Atm_Properties.nc`.
- All-domain statistics are cos(latitude) area-weighted. Figures follow AGU
  style (Arial, ≥8 pt, 300 dpi, (a)/(b)/(c) panel labels) via
  `src/era5lib/plotstyle.py`.

## Metric definitions & references

| Variable | Definition |
|---|---|
| `sbi_strength` | T(inversion top) − T(2 m), profile scan from the surface |
| `sbi_top_p`, `sbi_depth_p`, `sbi_depth_z` | SBI top pressure and depth (hPa / approx. m) |
| `dt_850_2m` | T(850 hPa) − T(2 m) |
| `dt_925_1000` | T(925 hPa) − T(1000 hPa) |

**Surface-based inversion (SBI), profile scan.** The temperature profile is
scanned upward from the surface (base = 2 m temperature at surface pressure);
the inversion top is the last running temperature maximum before the profile
stops increasing, with up to `sbi.max_embedded_levels` consecutive
non-increasing levels tolerated inside the inversion — at the ~25 hPa level
spacing of ERA5, one tolerated level is roughly analogous to the 100 m
embedded-layer rule used with radiosondes. Strength = T(top) − T(2 m).
Below-ground pressure levels (p > surface pressure; ERA5 extrapolates them)
are excluded from the scan.

An inversion is only counted when its strength is **≥ 0.5 K**
(`sbi.min_strength_k`); weaker cases are reported as "no SBI" (NaN,
`sbi_found = 0`). The 0.5 K minimum follows the criterion used with the same
first-derivative algorithm at three Arctic radiosonde stations
([J. Appl. Meteor. Climatol., 2022, doi:10.1175/JAMC-D-21-0054.1](https://journals.ametsoc.org/view/journals/apme/61/4/JAMC-D-21-0054.1.xml));
a 0.3 K variant is used in the Ny-Ålesund high-resolution radiosonde
climatology
([Atmos. Res., 2021](https://www.sciencedirect.com/science/article/abs/pii/S016980952100082X))
to suppress inversions within measurement uncertainty. No established Arctic
SBI study uses a 1 K minimum, so 0.5 K is the default here.

- Kahl, J. D. (1990): Characteristics of the low-level temperature inversion
  along the Alaskan Arctic coast. *Int. J. Climatol.* **10**, 537–548.
- Serreze, M. C., J. D. Kahl and R. C. Schnell (1992): Low-level temperature
  inversions of the Eurasian Arctic and comparisons with Soviet drifting
  station data. *J. Climate* **5**, 615–629.
- Zhang, Y., D. J. Seidel, J.-C. Golaz, C. Deser and R. A. Tomas (2011):
  Climatological characteristics of Arctic and Antarctic surface-based
  inversions. *J. Climate* **24**, 5167–5186.

**T(850 hPa) − T(2 m).** Fixed-level estimate standard in Arctic
climate-model studies; 850 hPa sits above the boundary layer, and the 2 m
temperature is preferred over T(1000 hPa) because wintertime surface pressure
deviates from 1000 hPa.

- Medeiros, B., C. Deser, R. A. Tomas and J. E. Kay (2011): Arctic inversion
  strength in climate models. *J. Climate* **24**, doi:10.1175/2011JCLI3968.1.
- Pavelsky, T. M., J. Boé, A. Hall and E. J. Fetzer (2011): Atmospheric
  inversion strength over polar oceans in winter regulated by sea ice.
  *Clim. Dyn.*, doi:10.1007/s00382-010-0756-8.

**T(925 hPa) − T(1000 hPa).** Simplest fixed-level metric, requiring only
pressure-level data (~600 m minus ~30 m a.s.l.); masked by default where the
surface pressure is below 1000 hPa. Used with ERA5 in e.g. arXiv:2011.11127;
related two-level stability metrics are LTS (Klein & Hartmann 1993,
*J. Climate* **6**, 1587–1606) and EIS (Wood & Bretherton 2006, *J. Climate*
**19**, 6425–6432).

## January 2020 case study (current results)

January 2020 was chosen as the study month: deep polar night with the MOSAiC
expedition drifting at 85–88.6°N inside the domain (final ERA5, `expver 0001`).

- **Monthly climatology** (`monthly_*_202001.png`): domain-mean SBI frequency
  61.6%, conditional strength 6.4 K (median 6.8 K), depths mostly 400–800 m.
  Spatial pattern matches the winter radiosonde climatologies: near-permanent
  strong inversions over Greenland/CAA, ~60–85% over the pack ice, minimum
  over the Atlantic sector.
- **Profile PCA** (`profile_pca_202001.png`, `surface_t_202001.png`): EOF1
  (73% of 1000–400 hPa temperature variance) is a whole-column warm/cold
  mode; EOF2 (15%) is the shallow inversion mode. SBI strength correlates
  positively with both PCs and negatively with T2m over the pack ice — the
  strength–surface-temperature anti-correlation of Zhang et al. (2011).
- **MOSAiC validation** (`mosaic_compare_202001.png`): 123 January soundings
  matched (median offsets 1.1 h / 7 km). Detection agreement 78.9%; ERA5 SBI
  strength bias −2.0 K (r = +0.27), depth bias +315 m, T2m warm bias +2.95 K
  (r = +0.79) — consistent with the documented winter ERA5 warm-surface /
  weak-inversion bias. The fixed-level metrics are compared against the same
  metrics derived from the level-2 radiosonde profiles (Maturilli et al.
  2021) + tower 2 m temperature: T925−T1000 r = +0.80 (bias −1.0 K) and
  T850−T2m r = +0.60 (bias −3.0 K — almost entirely the +2.95 K T2m warm
  bias). Both track their observed counterparts far better than either SBI
  definition tracks the other (r = +0.27), showing the SBI disagreement is
  mostly criterion/vertical-resolution, not thermodynamic error. (For
  reference vs obs SBI strength the proxies read ~4 K low.)
- **LW closure test** (`lw_sim_20200101T12Z*.png`, `lw_sim_stats_*.png`):
  libRadtran + RRTMG-LW broadband fluxes vs ERA5 at 12 UTC. With n = 500
  random pixels per sky: clear-sky LW↑ closes to **−0.16 W/m²** population
  bias after the far-IR tail correction (libRadtran, r = 0.997) and
  **−0.14 W/m²** for RRTMG-LW (rmse 0.57); the two codes agree on LW↑ to
  +0.02 ± 0.07 W/m², directly validating the analytic tail estimate.
  Clear-sky LW↓ scatters vs ERA5 (r ≈ 0.1, rmse 10.7/12.6 W/m²) for **both**
  codes even though they agree with each other in the mean (+0.11 W/m²) —
  diagnosed as a state-time mismatch, not radiative transfer: ERA5 `strd` is
  an 11–12 UTC accumulation produced inside the 06 UTC forecast, while the
  simulations use the 12 UTC analysis profile. Clear-sky LW↓ recomputed from
  the 06 vs 12 UTC analyses at the same (still-clear) pixels shifts by
  −10.7 ± 11.7 W/m² while LW↑ moves only ±0.9 — reproducing the observed
  LW↓/LW↑ error asymmetry. Cloud-edge proximity and `strd` spatial-gradient
  tests came back null (|r| ≈ 0.1), ruling out radiation-grid smearing as the
  dominant term. Cloudy (overcast, ERA5 clwc/ciwc as wc/ic files): LW↓
  r = 0.947 (bias +1.3, rmse 6.2 W/m², widest for thin clouds), LW↑ r = 0.983
  (bias+tail −1.3 W/m²); RRTMG's cloudy LW↓ runs ~+5.8 W/m² above
  libRadtran+tail from cloud-optics parameterization differences. Extending
  the simulated range from 4 µm down to ERA5's 3.08 µm edge was verified to
  matter only ~0.05 W/m².

## LW flux simulation (stage 7)

Broadband clear-sky thermal fluxes simulated with libRadtran (uvspec/DISORT,
via er3t) from ERA5 profiles + skin temperature, compared against ERA5 `strd`
and `LW_up = strd − str`. **Runs under the `er3t_env` conda env** (needs er3t
importable and `$LIBRADTRAN_V2_DIR` set), not the `era5` env.

```bash
conda activate er3t
# one-time: CAMS EGG4 CO2/CH4 monthly profiles (ADS licence required; see --help)
python src/era5_cams_download.py --year 2020 --month 1

python src/era5_lrt_sim.py prep    --year 2020 --month 1 --day 1 --hour 12   # 5 clear pixels across SBI range
python src/era5_lrt_sim.py run     --year 2020 --month 1 --day 1 --hour 12   # uvspec thermal jobs
python src/era5_lrt_sim.py compare --year 2020 --month 1 --day 1 --hour 12   # table + figure

# cloudy mode: 5 near-overcast pixels (liquid/ice/mixed/thin/thick when
# available); ERA5 clwc/ciwc become 1D wc/ic cloud files with fixed effective
# radii (liquid 10 um via Hu & Stamnes, ice 25 um via Fu)
python src/era5_lrt_sim.py prep    --year 2020 --month 1 --day 1 --hour 12 --sky cloudy
python src/era5_lrt_sim.py run     --year 2020 --month 1 --day 1 --hour 12 --sky cloudy
python src/era5_lrt_sim.py compare --year 2020 --month 1 --day 1 --hour 12 --sky cloudy
```

Cost: ~1.4 s per clear pixel and ~4 s per cloudy pixel serially on the local
Mac (reptran coarse, 4 streams); jobs parallelize across `--workers`.
Cloudy caveats: plane-parallel overcast assumption (pixels are screened for
cloud fraction ≥ 0.99), fixed effective radii (ERA5 diagnoses its own
internally), and ERA5's cloud overlap scheme — expect larger spread than
clear-sky, especially for optically thin clouds. The fixed radii matter most
exactly there: LW absorption scales ~1/r_eff, so thin clouds (τ ~ 1) shift
by up to ~20 W/m² across plausible radii while emissivity-saturated thick
clouds (liquid LWP ≳ 30–40 g/m²) are insensitive (< 0.1 W/m²).

Setup: surface emission from ERA5 `skt` (uvspec `sur_temperature`) with
emissivity 0.99; ERA5 T/q/O3 profiles + CAMS CO2/CH4 (or
`--fallback-constants`); subarctic-winter standard atmosphere above 1 hPa;
integrated flux over 3.08–100 µm (reptran thermal supports 2.5–100 µm; the
lower edge matches ERA5's RRTMG-LW, and the 3.08–4 µm band was verified to
contribute only ~0.05 W/m² at Arctic winter temperatures). Resolution knobs
for local vs cluster:
`run --streams {4,8,...} --mol-abs-param {coarse,medium,fine}` (defaults:
4/coarse on macOS, 8/coarse on Linux); `slurm/curc_lrt_sim.sh` is a CURC
template. Caveats printed with every comparison: ERA5's RRTMG-LW extends to
1000 µm (its first band, 10–350 cm⁻¹, covers the far-IR explicitly), so the
simulation misses the far-IR >100 µm tail (~1.5–2 W/m² at Arctic winter
temperatures; an analytic estimate is reported), and ERA5 fluxes are
1-h accumulations while the simulation is instantaneous — residual cloud
within the accumulation hour shows up as ERA5 > simulation in LWdn.

The `compare` figure annotates the far-IR tail equation
(tail = [1 − F_band(T)]·εσT⁴ with F_band the fractional blackbody emissive
power in the simulated band; T = T_skin, ε = 0.99 for LW↑ and
T = (LW↓_ERA5/σ)^¼, ε = 1 for LW↓) and shows the surface blackbody estimate
εσT_skin⁴ + (1−ε)·LW↓ alongside the upwelling fluxes as a closure check.

## RRTMG-LW cross-check (stage 7b)

`src/era5_rrtmg_sim.py` reruns every pixel of an existing manifest through
RRTMG-LW — the same radiation-scheme family ERA5's IFS uses, covering the
full 3.08–1000 µm range (no far-IR tail correction needed) — via
[climlab](https://climlab.readthedocs.io) + climlab-rrtmg (conda-forge).
**Runs under the `era5` env** (climlab 0.9.2 / climlab-rrtmg 0.4.2 are
installed there), not `er3t_env`:

```bash
conda activate era5
python src/era5_rrtmg_sim.py --year 2020 --month 1 --day 1 --hour 12
python src/era5_rrtmg_sim.py --year 2020 --month 1 --day 1 --hour 12 --sky cloudy
# then regenerate the (now three-way) comparison under er3t_env:
conda activate er3t_env
python src/era5_lrt_sim.py compare --year 2020 --month 1 --day 1 --hour 12
```

Input identity with the uvspec runs is guaranteed by parsing the same profile
files (atmosphere `.dat` incl. the afglsw splice, CH4 mol_file, wc/ic cloud
files with libRadtran's layer convention). Cloud optics are matched by family:
liquid Hu & Stamnes (`liqflglw=1`), ice Fu 1996 (`iceflglw=3`,
dge = 1.0315·reff). Rows land in the same results CSV under the reserved
`simulator` column; `compare` then adds RRTMG points and a head-to-head
panel (d) — RRTMG vs libRadtran+tail — that isolates RT-scheme + spectral-
coverage differences on identical inputs. Cost: ~3 ms per column (500 pixels
in under 2 s).

Implementation note: climlab's wrapper sets the bottom interface temperature
to the *skin* temperature, which leaks skt into the lowest air layer's
emission and adds several W/m² of spurious LW↓ error under strong surface
inversions; `era5_rrtmg_sim.py` patches it to the surface *air* temperature
(t2m), matching libRadtran and ERA5's own use of RRTMG. Remaining documented
input differences: N2O constant 0.32 ppm (libRadtran uses its default
profile) and per-layer means vs level profiles — a ×4 sublayer-refinement
test showed the residual clear-sky LW↓ code-to-code spread (±6.6 W/m²,
correlated with column water vapor, r = +0.43) is band-model/continuum
physics, not discretization.

Future extension (recorded, not implemented): RRTMG comparison directly on
ERA5 137-model-level profiles to quantify the 37-pressure-level truncation.

## MOSAiC case study (single-pixel walkthrough)

`src/era5_case_study.py` picks the best clear-sky and overcast ERA5 pixels
among the MOSAiC-matched soundings and produces one 4-panel figure per case
(`figures/case_study_{clear,cloudy}_*.png`): (a) ERA5 profile vs the full
MOSAiC level-2 radiosonde profile (Maturilli et al. 2021,
doi:10.1594/PANGAEA.928659; auto-downloaded to `data/mosaic/soundings/`)
plus tower/BL-top points and cloud water, (b) a zoom showing how the ERA5
and observed SBI strengths are constructed (ΔT arrows, criteria below the
panel), (c) the radiative-transfer input (all 37 levels + surface row +
afglsw splice, level bookkeeping, solver settings, collocation), and
(d) surface LW↓/LW↑ from libRadtran (± far-IR tail) and RRTMG-LW against
both the ERA5 flux and the MOSAiC surface radiometer observation.

```bash
conda activate era5     && python src/era5_case_study.py prep
conda activate er3t_env && python src/era5_case_study.py run
conda activate era5     && python src/era5_case_study.py figure
```

Jan 2020 picks — clear: 2020-01-21 18 UTC (87.50°N, 95.75°E; ERA5 SBI 7.1 K
vs obs 8.9 K; both simulators within 1.5 W/m² of ERA5 LW↑ after tail, and
both ~13 W/m² below ERA5/obs LW↓, the documented clear-sky LW↓ state-time
mismatch). Cloudy: 2020-01-08 06 UTC (87.00°N, 115.25°E; IWP 47 g/m²;
simulators within ~7 W/m² of ERA5 for both components — but ERA5 places ice
near the surface while the observed lowest cloud base was 5.5 km, and the
observed LW↓ is ~40 W/m² below ERA5, a reminder that closure against ERA5 is
not closure against reality).

## MOSAiC drift flux simulation (all matched columns)

`src/era5_mosaic_flux.py` (prep/run/figure, same env split as the case
study) extends the walkthrough to **every** ERA5 column matched to a MOSAiC
sounding — 123 unique (pixel, 6-h time) columns in January 2020. Each column
is simulated twice (clear-sky, and overcast plane-parallel wherever ERA5
holds condensate) with both libRadtran and RRTMG-LW; an all-sky flux is
blended with the random-overlap effective cloud fraction f = 1 − Π(1 − cc).
Soundings sharing a (pixel, time) key are simulated once with their observed
fluxes averaged (no pseudo-replication; 123 soundings → 123 unique keys this
month). `figures/mosaic_flux_YYYYMM.png` shows the month-long drift time
series of LW↓/LW↑ and scatters vs the MOSAiC surface radiometers, with the
ERA5 flux product as a third contender. Cost: ~1 min of uvspec locally.

Jan 2020 verdict vs the radiometers (n = 118): LW↑ is biased **+11.5 W/m²
for simulations and ERA5 alike** (r ≈ 0.79) — the documented ERA5 warm skin
bias, inherited by the simulations through skt. LW↓: libRadtran r = 0.72 /
RRTMG r = 0.70 vs ERA5's own product r = 0.67 (all ~ +10 W/m² high, rmse
24 W/m²) — the offline simulations track the radiometer slightly better
than ERA5's flux product, and the shared overestimate points at ERA5's
cloud state along the drift, not radiative transfer.

## PREFIRE brightness-temperature simulation + Jacobians (stage 8)

![PREFIRE workflow](docs/workflow_prefire.png)

`src/era5_prefire_download.py` (era5 env) fetches PREFIRE TIRS L1B spectral
radiance granules (`PREFIRE_SATx_1B-RAD` R01, NASA ASDC via `earthaccess`;
Earthdata login in `~/.netrc`) and the TIRS spectral-response files (v13,
Zenodo record 16638853 — the version the R01 calibration used). PREFIRE
observes to ~83.8°N, so the overlap with the 80–90°N domain is the
80–83.8°N band; viewing zenith angles are 5–16°, and each of the 8
cross-track scenes has its own wavelength registration.

`src/era5_prefire_bt.py` runs the stage:

```
conda activate era5     && python src/era5_prefire_bt.py collocate --year 2025 --month 1 --sat 1
conda activate era5     && python src/era5_prefire_bt.py prep      --year 2025 --month 1 --sat 1
conda activate er3t_env && python src/era5_prefire_bt.py run       --year 2025 --month 1 --sat 1
conda activate era5     && python src/era5_prefire_bt.py rrtmg     --year 2025 --month 1 --sat 1
conda activate er3t_env && python src/era5_prefire_bt.py jacobian  --year 2025 --month 1 --sat 1 --simulator lrt
conda activate era5     && python src/era5_prefire_bt.py figure    --year 2025 --month 1 --sat 1
```

- **collocate** maps every good-quality footprint to its ERA5 0.25° cell and
  nearest analysis hour (≤ 3 h offset — the same state-time caveat as stage
  7) and picks a test set of clear and single-class overcast columns;
  partially cloudy columns are excluded because BT does not blend linearly
  across a broken scene.
- **run** simulates one thermal-source spectral radiance per column with
  libRadtran at the footprint viewing angle (`mie` liquid / `yang2013` ice —
  radiance-grade optics, unlike the flux-oriented Hu & Stamnes / Fu of stage
  7), convolves with the scene's SRF and inverts the SRF file's own
  blackbody channel-radiance lookup, so simulated and observed BT share one
  radiometric scale. uvspec's thermal band output is radiance per
  wavenumber; the reader converts accordingly (verified against the Planck
  curve).
- **rrtmg** is a 16-band sanity check: `climlab` RRTMG-LW spectral OLR
  converted to flux-equivalent band BT. It agrees with band-aggregated
  libRadtran to ±1–2 K (flux-equivalent vs nadir-radiance BT differ by
  limb darkening, so this is a consistency check, not validation).
- **jacobian** builds finite-difference K matrices around each column's
  actual sky state: skt and per-level T (+1 K), per-level q (+5 %),
  ln LWP/IWP (+10 %), r_eff (+1 µm), cloud-top (one layer up), emissivity
  (−0.01). `--simulator rrtmg` sweeps all states in ~2 s/column (structure
  prototyping); `--simulator lrt` gives the channel-resolved K for the
  planned cloud-property retrieval (future validation target: collocated
  EarthCARE cloud products, which need an ESA EO account).

A `cotscan` subcommand (er3t_env) reproduces the ARCSIX-style
BT-vs-cloud-optical-thickness sensitivity figure with PREFIRE channels: a
synthetic single-layer cloud (`--phase/--cer/--cth/--cbh`) is inserted into
one collocated column, its 550-nm optical thickness swept on a log grid via
`ic_modify tau set`, and each state is one spectral run convolved to channel
BT — panels show BT(τ), dBT/dτ and dBT/dr_eff
(`figures/prefire_cotscan_*.png`).

**Resolution note:** `reptran fine` must NOT run on a laptop — parallel
fine-grid workers exhaust memory (this has crashed a machine); the CLI
refuses fine/medium on Darwin and defaults to `coarse` locally. Production
fine-grid runs belong on CURC (`slurm/curc_prefire_bt.sh`). The local coarse
grid (~15 cm⁻¹) under-resolves far-IR channels whose widths shrink to
~3–9 cm⁻¹, so science numbers should come from the fine run.

First results (2025-01-01, SAT1, coarse, 3 clear + 3 overcast columns):
clear-sky sim − obs ≈ **+5 K** in window and far-IR channels — consistent
in direction with the +11.5 W/m² warm-skin LW↑ bias found along the MOSAiC
drift; overcast columns scatter −15…+2 K (ERA5 cloud placement). Clear-sky
Jacobians close (ΣK_T + K_skt ≈ 1), K_skt falls from ~0.98 (11 µm window)
to 0.08 (26.5 µm) across the far-IR dirty window, and opaque ice cloud
gives ∂BT/∂lnIWP ≈ −7 K per 100 % IWP with the surface fully masked
(K_skt = 0) — the expected information content for a retrieval.

## Notes on surface fluxes

The single-levels download includes turbulent fluxes (`sshf`, `slhf`) and
radiation (`ssrd`, `ssr`, `strd`, `str`) for future analyses. All are
accumulations over the hour ending at `valid_time` in J m⁻² (divide by 3600
for mean W m⁻²), positive downward. ERA5 has no explicit upward radiation
variables; recover them from downward minus net:

```
SW_up = ssrd - ssr        LW_up = strd - str
```

## Notes

- Recent months come from ERA5T (preliminary, `expver 0005`); the pipeline
  warns when a file contains ERA5T data and records `expver` in the outputs.
- `sbi_depth_z` is a hypsometric estimate from T and q; adding `geopotential`
  to `download.plev_variables` in `config.yaml` would allow exact heights.

## Data sources

- ERA5 (Hersbach et al. 2020, *QJRMS* **146**, 1999–2049), Copernicus Climate
  Data Store.
- MOSAiC lower-atmospheric properties: Jozef, G. C., et al. (2023), *ESSD*
  **15**, doi:10.5194/essd-15-4983-2023; dataset doi:10.1594/PANGAEA.957760
  (CC-BY-4.0).
