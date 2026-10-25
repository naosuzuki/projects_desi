# Handover: quasar eigenspectra (new thread)

Date: 2026-10-24. From: the "projects_desi" session (galaxy/star/quasar template design, pass B running).
Purpose: let a fresh session start the **quasar template work** without re-deriving the context.
Read this, then `~/github/papers/26_sdssdesi_spectra/ms.tex` §3.7 (quasar plan) and §3.9 (processing recipe).

## 1. The project in one paragraph
We build eigenspectra for a spectral classification/redshifting tool from **SDSS DR8 + DR17 and DESI DR2 (loa)**
spectra, separately for **galaxies, quasars, stars**. Common preprocessing (ms.tex §3.9): (1) spectrophotometric
scale against **fiber** photometry, (2) Milky Way extinction (Fitzpatrick99, R_V=3.1, E(B−V)×0.86), (3) rest frame on a
log10 λ grid (Δ=1e-4 = 69 km/s, integer-pixel shift), (4) absolute magnitudes from **whole-galaxy** flux with flat ΛCDM
h=0.7 Ωm=0.3, K from the object's own spectrum, (5) per-spectrum measurements, (6) stacking in a grid of physical
properties, then a single weighted EM-PCA (lines + continuum together — Nao's decision) on stacks + high-S/N
individuals. All documentation goes into **ms.tex** (Overleaf repo `~/github/papers/26_sdssdesi_spectra`, branch main;
Nao commits). Scripts: `~/github/projects_desi/programs_quasars/` (README.md there; common pieces imported from `programs_templates/`).

## 2. Data — everything is on this Mac (Naos-Mac-mini, M4 Pro 14 cores, 64 GB, 20-core GPU/Metal)
| data | where | notes |
|---|---|---|
| DESI DR2 loa healpix coadds, main dark+bright | `/Volumes/exdisk2/desiredux/DR2/spectro/redux/loa/healpix/main/{dark,bright}/...` (33,539) + 969 overflow on the exFAT disk, mounted as `/Volumes/exdisk1` (formerly `/Volumes/My Book`) | 34,508 files, 24 TB, checksum-verified; **read-only, never write there** |
| `zall-pix-loa.fits` | `~/desiredux/` (51 GB) | 64.1 M rows |
| SDSS DR17 spPlate (v5_13_2) | `~/sdssredux/v5_13_2/PLATE/` | 4,131 files, all of spAll |
| SDSS DR8 spPlate | `~/sdssredux/DR8/{26,103,104}/PLATE/` | 3,975 files, all of specObj-dr8; `sdss_db.py` searches these |
| `spAll-v5_13_2.fits`, `specObj-dr8.fits` | `~/sdssredux/` | |
| catalog cache (Parquet) | `~/data/sdssdesi_products/`: `zall_loa.parquet` (4.7 GB), `spall_dr17.parquet`, `specobj_dr8.parquet`, `zlookup_main.parquet`, `zlookup_compact.npz` | persistent; rebuild notes in programs_templates/README.md |
| Gaia DR3 × SDSS quasars | `~/github/projects_gaia/csvfiles/gaiadr3_sdssdr{8,17}_quasar*.csv`, XP spectra `~/github/projects_gaia/data/gaiadr3_xpspec_sdssdr{8,17}_quasar.fits` | for absolute flux checks of bright quasars |
| redrock + templates | `~/github/redrock` (0.22.2), `~/github/redrock-templates` (`RR_TEMPLATE_DIR`) | QSO-LOZ/HIZ v1.1 = reference |
| sdsslib | `~/github/sdsslib/pylib` (`PYLIB`) | SDSSspec reader, log-rebin, 2D stack format |

Env (in `~/.zshrc`): `SPECTRO_REDUX=~/sdssredux/`, `DESI_REDUX=~/desiredux/`, `PYLIB`, `DESILIB`, `RR_TEMPLATE_DIR`.
Python: /opt/miniconda3 python3.13; numpy 2.1, astropy 7.2, fitsio, speclite, extinction, dustmaps (SFD maps work),
healpy, redrock, torch 2.12 (MPS works), psutil.

## 3. What is running right now — do not disturb
**Pass B** (per-spectrum measurements over all 34,508 DESI coadds) started 2026-10-24 11:53, two detached workers
(`02_passB_measure.py --all --shard {0,1}/2`, nice 10, OPENBLAS_NUM_THREADS=1), ~30 files/min, **ETA Sunday
2026-10-25 ~07:00**. Output: `~/data/sdssdesi_passB/{dark,bright}/hpx<N>.parquet` (one per coadd; resumable).
Monitor + resource guard: `programs_monitor/monitor.py` (all steps, ETA) → **http://127.0.0.1:8765/monitor.html**; logs in `~/data/sdssdesi_status/`.
`caffeinate` holds keep the Mac awake until the workers exit. Resources used: ~2.5 cores, ~10 GB; ~33 GB free.
A separate Prospector job (HSC IGM project) also runs on the 10 performance cores.
**Rules for the new thread while pass B runs:** no reboot; nothing heavy on the external disks (they are the
bottleneck); keep ≥ 16 GB memory free (the guard pauses pass B below 8 GB); catalog-level work is fine.

## 4. What pass B gives quasars, and what it does NOT
Every spectrum (stars, galaxies, **quasars**) gets: TARGETID, Z/ZWARN/SPECTYPE/DELTACHI2 (from zall), S/N per camera
and total (`SNR_ALL`, `SNR_SDSSPIX`), synthetic fiber g/r/z vs catalog (`CAL_DG`, `CAL_DR`, `CAL_A` scale, `CAL_B_EST`
tilt estimate, `CAL_OK`, `CAL_MODE`), Dn4000, fiber fraction, M_R/M_G/M_U/M_NUV with K, and the **galaxy** emission
lines (Hα … [O II]). It does **not** measure quasar lines (Lyα, C IV, C III], Mg II, Hβ/[O III] at z<0.9 are covered
by the galaxy list only), BAL flags, or anything in the Lyα forest. Columns: see `02_passB_measure.py` (`measure_chunk`).

## 5. The quasar plan (ms.tex §3.7) — summary and open decisions
Sample: DESI 2.72 M quasars after cuts (median S/N≈3; 913k with S/N>5, 94k with S/N>20, at all z); SDSS DR17 667k
(233k >5), DR8 83k (67k >5). S/N does not fall with z. Rest coverage: 2400–6550 Å (z=0.5) … 600–1640 Å (z=5).
Plan: single weighted EM-PCA (lines+continuum) on rest 600–9800 Å with the **Lyα forest masked** (rest <1216 Å) and
a mean IGM transmission (Inoue+14 / Kamble+20) applied at fit time; **BALs excluded** from training; systemic z from
Mg II/[O III] for z<2.5 before de-redshifting; **host-galaxy component** (galaxy eigenspectra) in the fit at z<1;
stacking cells (Δz=0.1, luminosity, continuum-slope tertile); training set = all non-BAL S/N>5 quasars (~1.2 M) +
stacks; validation vs Redrock QSO v1.1 + QuasarNET on held-out spectra and DESI EDR VI truth tables; BOSS DR17
quasar plates need the **Margala et al. 2016** blue/red flux correction.
Decisions written as recommendations, not yet confirmed by Nao: (a) re-derive systemic z before de-redshifting;
(b) apply the Margala correction in Step 1 for DR17 quasars.

## 6. Suggested first tasks for the quasar thread (catalog-level, safe while pass B runs)
1. **BAL catalog**: DESI DR2 QSO VAC / BAL catalog (`BAL_PROB`) is not on disk; find the public DR2 (or DR1) QSO
   catalog and download (small). SDSS: DR16Q (Lyke+2020) has BAL flags and systemic redshifts — download (~1 GB).
2. **Quasar line list for a pass-B extension**: Lyα 1215.67, N V 1240, Si IV 1397, C IV 1549, C III] 1909, Mg II 2799,
   Hβ/[O III], Hα; plan a `--qso` mode of `02_passB_measure.py` (or a small pass over the ~2.7 M quasar spectra only,
   using the pass-B table to pick files/rows) measuring these lines, continuum slope α_λ, and a forest flag.
3. **Systemic-redshift recipe**: Mg II (z<2.5), [O III] (z<0.9); compare to pipeline Z; estimate C IV blueshift vs L.
4. **Luminosity**: M_i(z=2) à la Richards+06 from FLUX_R/FLUX_Z with the quasar K-correction, for the cells.
5. **Margala correction**: locate the published correction (per plate/fiber, BOSS) and plan its application to DR17.
6. **Gaia XP check**: synthetic photometry of the brightest SDSS quasars vs XP spectra (projects_gaia) → absolute
   flux check of the BOSS correction.
7. When pass B finishes (Sunday): pull the quasar rows from `~/data/sdssdesi_passB`, confirm S/N distribution vs the
   census in ms.tex Table (qsocensus), select S/N>5 non-BAL, and build the first stacks.

## 7. Conventions
- Document in ms.tex (Nao commits); scripts in `projects_desi/programs_templates/` (numbered); caches in `~/data/`.
- Never write to the external disks. Scratchpad files are lost between sessions — keep anything useful in the repo
  or `~/data/`.
- Session is filed under the sidebar group **projects_desi**.

## 8. Status of the quasar thread (2026-10-24, evening) — what was done, what is next
Pass B still running (checked 2026-10-24: both shards alive, monitor up, ~33 GB free); nothing below touched the
external disks beyond one test coadd.

### 8.1 Catalogs located (nothing downloaded yet except the 6 KB Richards table; sizes from HTTP HEAD / index pages)
| catalog | URL | size | what it gives us |
|---|---|---|---|
| SDSS DR16Q (Lyke+20) | https://data.sdss.org/sas/dr16/eboss/qso/DR16Q/DR16Q_v4.fits | 1.72 GB, 750,414 rows, 183 cols | `BAL_PROB`, `BI_CIV`, `AI_CIV`, `BI_SIIV`, `AI_SIIV`, `Z_PCA`, `Z_MGII`, `Z_CIV`, `Z_HBETA`, `Z_LYA` (+ `ZWARN_*`), `M_I`, PLATE/MJD/FIBERID → BAL flags + systemic z for the DR17 (eBOSS) quasars and the DR8 ones re-observed by BOSS |
| DESI DR1 Lyα redshift/BAL VAC ("zlya") | https://data.desi.lbl.gov/public/dr1/vac/dr1/zlya/v1.0/qso_cat_dr1_main_dark_healpix_zlya-v0.fits | 281 MB | main dark, z>1.6 (ZWARN 0 or 4): `BI_CIV`, `AI_CIV`, `NCIV_2000/450`, `VMIN/VMAX_CIV_*`, `ZMASK` (z refit with BAL masked). **No `BAL_PROB`** — the DR1/DR2 catalogs give AI/BI, not a probability |
| DESI DR1 QSO catalog (afterburner) | https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/QSO/iron/QSO_cat_iron_cumulative_v0.fits | 816 MB | QuasarNET / Mg II afterburner redshifts and the final QSO selection for DR1 (column list not on the public doc page; read the header after download) |
| DESI EDR BAL VAC (Filbert+23) | https://data.desi.lbl.gov/public/edr/vac/edr/bal/fuji/v1.0/QSO_cat_fuji_sv{1,3}_dark_bal_mask_v1.0.fits | 45 + 60 MB | SV BAL catalog (balfinder columns, AI/BI); useful for validating our own BAL proxy |
| DESI DR1 BH-mass VAC (Mg II, 0.6<z<1.6) | https://data.desi.lbl.gov/public/dr1/vac/dr1/qmassiron/v1.7/VAC_BHmass_338_v1.7.fits | 109 MB | Mg II line fits → check of our Mg II measurements |
| Richards+06 i-band K(z) to z=2 | `data/quasar/richards06_table4_kcorr.dat` (saved; from CDS J/AJ/131/2766 table4, 0≤z≤5.49 step 0.01) | 6 KB | M_i(z=2) for the luminosity axis of the stacking cells |
DESI **DR2 is not public** (https://data.desi.lbl.gov/public/dr2/ is empty), so there is no DR2 BAL catalog to fetch:
BAL flags for the DR2 training set come from (i) DR1 zlya for the overlap (TARGETID match; DR1 ⊂ DR2 for main dark),
(ii) DR16Q for SDSS, (iii) our own AI_CIV proxy from pass Q (below) + sigma clipping for the rest.

### 8.2 Margala correction — resolved: nothing to apply to DR17
The SDSS DR14+ pipeline (v5_10_0 onward; our DR17 spPlates are v5_13_2, header `RUN2D = v5_13_2`) applies the
offset-fiber throughput corrections (Margala+16 recipe) to the individual frames (`spXYthrucorr` files) **before**
coaddition (SDSS DR17 spectrophotometry page, "In DR14, these corrections are applied to individual frames ... The
frames are co-added after this correction is applied"; "Spectrophotometry algorithms are unchanged in DR17").
DR8 is SDSS-I/II (no offset fibers). The standalone correction site http://darkmatter.ps.uci.edu/tpcorr/ no longer
resolves. → ms.tex §3.7 item 6 should say that DR17 quasar spectra already carry the correction and that the Gaia XP
check (task 6) is the test, not a step; decision (b) of §5 is moot. Open: Lan+21 (A&A 655, A53) note the DR14+
correction ignores the wavelength dependence of seeing (σ ∝ λ^−0.2); residual slope errors of a few % remain possible.

### 8.3 Pass Q — quasar measurements (`programs_templates/05_passQ_measure.py`, written and tested)
Reuses pass B (import): for the QSO rows of each coadd (SPECTYPE=QSO, ZCAT_PRIMARY; ZWARN kept as a column) it runs
the full pass-B chain (same columns: CAL_*, SNR_*, MABS_*, galaxy lines incl. narrow Hβ/[O III]/Hα, Dn4000), then
rebuilds the calibrated de-reddened spectrum on a 600–9900 Å rest grid and measures: power-law slope `ALPHA_LAM`
(χ² in linear flux over the Vanden Berk line-free windows, α on a 0.02 grid, `CONT_LOGF3000`, `CONT_CHI2NU`);
broad lines O VI 1034, Lyα+N V, O I 1304, Si IV, C IV, O III] 1663, Al III+C III], Mg II, broad Hβ, broad Hα (Gaussian
σ grid 600–5000 km/s on a local linear continuum: `_FLUX/_EW/_SIGMA/_VOFF/_CONT/_SNR`), narrow [O III] 5007 and
[O I] 6300 (σ 100–600); [O II] 3727 and [O III] 4959/5007 fluxes come from the pass-B list (`OII_SUM_*`, `OIII_5007_*`); line redshifts
`Z_MGII, Z_OIII, Z_CIV, Z_HBETA` from the centroid offsets; forest `FOREST_NPIX`, `FOREST_TRANS`; BAL proxy
`BAL_AI_CIV/_SIIV` (AI-type index, troughs ≥450 km/s, 5-pixel smoothed flux >1σ below 0.9×continuum),
`BAL_NTROUGH_*`, `BAL_VMIN/VMAX_*`; `REST_LAMMIN/MAX`. Output `~/data/sdssdesi_passQ/{dark,bright}/hpx<N>.parquet`.
Test (coadd-main-dark-11274, 112 quasars, 1 s + 0.5 s read): α_λ median −1.60 ± 0.01 for S/N>5 (Vanden Berk+01:
−1.56); C IV `VOFF` median −313 km/s, Mg II −28, [O III] −21, Lyα −293 (as expected); Z_MGII − Z median −2.4e-4.
Caveats: the BAL proxy still fires on ~50% of z>2.2 spectra at S/N<2 (1 of 6 at S/N>3) → use it only at S/N>3 and
prefer catalog AI/BI; FOREST_TRANS uses the power law extrapolated below 1275 Å (indicative only).
Cost: the file read dominates (same as pass B) → full run ≈ 26 h with one worker, ≈ 13 h with two; **do not start it
while pass B runs** (disk bottleneck). Suggested: start Sunday after pass B, same detached/nice recipe
(`05_passQ_measure.py --all --shard i/2`, `OPENBLAS_NUM_THREADS=1`), monitor dir `~/data/sdssdesi_passQ`.

### 8.4 Next steps (in order)
1. Nao: approve the downloads in 8.1 (DR16Q 1.7 GB; DR1 zlya 281 MB; DR1 QSO cat 816 MB; EDR BAL 105 MB; BH-mass
   109 MB) into `~/data/sdssdesi_products/quasar/` (catalog-level, safe during pass B).
2. Sunday, after pass B: run pass Q over all coadds (13 h, two workers); then SDSS pass (spPlates, internal disk).
3. Luminosity axis: M_i(z=2) = m_i(dereddened) − DM(z) − K(z) with Richards table; DESI has no i → use
   synthetic i from the calibrated spectrum (pass Q can add `SYN_I`), or FLUX_R/FLUX_Z interpolation. Decide.
4. Systemic z: Z_MGII (z<2.5) / Z_OIII (z<0.9) from pass Q vs DR16Q Z_MGII / BH-mass VAC on overlaps; C IV blueshift
   vs luminosity (Baldwin/blueshift trend) as a check — then decision (a) of §5.
5. BAL exclusion: zlya AI_CIV>0 (DR1 overlap) ∪ DR16Q BAL_PROB>0.5 ∪ pass-Q AI proxy at S/N>3; sigma clipping backstop.
6. Gaia XP check of the DR17 flux calibration (task 6 of §6).
7. Edit ms.tex §3.7 item 6 (Margala) as in 8.2; add the pass-Q line list to §3.9.

### 8.5 Proposed quasar grid and cell counts (2026-10-24; catalog estimate, S/N>5, DESI+DR17+DR8 = 1,242,186 quasars)
Inputs: DESI S/N predicted from log S/N = −0.449 + 0.846 log FIBERFLUX_R + 0.430 log TSNR2_QSO (fit on 5,946 quasars
of the pass-B test coadds, rms 0.13 dex); SDSS SN_MEDIAN(_ALL); M_i(z=2) = m_i − DM − K_Richards(z), m_i from LS r,z
interpolated in log λ (DESI) or PSFMAG_i (SDSS). Table: `data/quasar/grid_cells_qso_z_Mi_snr5.csv`.
Axes: z Δ=0.1 (0–3.5), 0.25 (3.5–4), 0.5 (4–5) → 39 slices; M_i(z=2) 1-mag bins −31…−19; continuum-slope tertiles
within each (z, L) cell; BAL class for z ≥ 1.5 only: non-BAL / weak (AI_CIV>0, BI=0) / strong (BI_CIV>0), the two
BAL classes split into 3 trough-velocity bins (VMIN_CIV: 0–5000, 5000–15000, >15000 km/s) and not by slope.
Counts with floor 300 per stack: (z, L) cells 165 (holding 1,236,713 of the quasars); non-BAL × tertiles 131 cells =
393 stacks; BAL (fractions assumed: weak 7%, strong 8%) 18 weak and 21 strong (z, L) cells × 3 velocity bins = 117
stacks; **total ≈ 510 stacks** (≈ 480 without the velocity split). The −28…−24 range dominates at 1<z<2.5; the
z>3.5 slices have 2–4 luminosity bins each; 4.5–5 has none ≥ 300 (individuals only, or merge with 4–4.5).
Tune after pass Q (real S/N and slopes) and after the BAL catalogs are in (real fractions and velocities).

### 8.6 Nao's choice (2026-10-24 evening): classify by emission-line strength, not slope/BAL
Primary line by redshift (line ±5% inside 3600–9800 Å): Hα (z<0.40), Hβ (0.40–0.90), Mg II (0.90–2.30), C IV (2.30–5).
Rest EW bins, factor 2 wide around the Vanden Berk+01 composite: Hα <100/100–200/200–400/>400 Å; Hβ <25/25–50/
50–100/>100; Mg II <15/15–30/30–60/>60; C IV <12/12–25/25–50/>50. Secondary lines where in the window: Si IV
(z>1.8; <4/4–8/8–16/>16). Oxygen axis (Nao, 2026-10-24): [O III] 5007 narrow (z<0.9; <5/5–15/>15 Å), [O II] 3727
(z<1.6; <2/2–5/>5 Å, from pass B), O I 1304 (1.8<z<5; <2/2–5/>5 Å), O III] 1663 (1.2<z<4.9; <1/1–3/>3 Å), O VI 1034
(z>2.6; <5/5–15/>15 Å); all are in pass Q now (test coadd: O VI EW med 16 Å at n=12, O I and O III] ≈ 0–2 Å, measurable
only at S/N>3 per spectrum → bin on them for stacks/high-S/N individuals, not for the S/N~3 majority).
Estimated stacks (EW log-normal σ=0.3 dex assumed, floor 300): z × M_i × EW(4) = 490; adding [O III](3) at z<0.9 = 606;
z × EW only (no luminosity axis) = 152. Real EW distributions come from pass Q → re-bin then. BALs: not an axis; flag
and keep (strong-BI quasars can be clipped at stack time if they distort a cell).

### 8.7 The quasar grid — consolidated plan (2026-10-24, after Nao's choices in 8.6)
Same construction as the galaxy grid (ms.tex Table 3): physical axes, a floor per stack, weighted mean + outlier pass.
| axis | bins | by redshift |
|---|---|---|
| 1. redshift | Δz=0.1 (0–3.5), 0.25 (3.5–4), 0.5 (4–5): 39 slices | all |
| 2. luminosity M_i(z=2) | 1 mag, −31…−19 (Richards+06 K; DESI i from LS r,z) | all |
| 3. primary broad line EW | 4 bins, factor 2 around the composite | Hα z<0.4 (<100/200/400); Hβ 0.4–0.9 (<25/50/100); Mg II 0.9–2.3 (<15/30/60); C IV 2.3–5 (<12/25/50) |
| 4. oxygen line EW | 3 bins | [O III] 5007 z<0.9 (<5/15); [O II] 3727 0.9–1.6 (<2/5); O I 1304 1.8–2.6 (<2/5); O VI 1034 z>2.6 (<5/15); O III] 1663 as alternative 1.2–2.6 |
| recorded, not binned | continuum slope α_λ, Si IV EW, C IV blueshift, BAL AI/BI + trough velocities, host fraction (z<1), forest transmission | — |
Estimated stacks (S/N>5 sample of 1.24 M; EW fractions assumed 16/34/34/16 % and 25/50/25 %), floor 300:
z × L = 165; × primary EW = 490 (1.21 M quasars inside); × oxygen = 986 (floor 500: 155 / 420 / 754).
By range: z<0.4: 46 → 84 with oxygen; 0.4–0.9: 70 → 148; 0.9–1.6: 108 → 248; 1.6–2.3: 108 → 242; 2.3–3.5: 136 → 248;
3.5–5: 22 → 16 (oxygen axis dropped there: too few). Target ≈ 500 stacks with the oxygen axis used only where every
sub-cell keeps ≥300 (else collapse it) — same rule as the galaxy grid's "measured on 3 axes".
Procedure: pass Q (Sunday) → real EW distributions → fix bin edges per line (keep the factor-2 scheme, shift the
center to the measured median per z slice) → cell table (`data/quasar/cells_qso.csv`) → pass C stacks per z slice
(2D FITS, sdsslib format) → EM-PCA on stacks + S/N>5 individuals. Cell assignment of low-S/N spectra uses the
primary line only (oxygen EW at S/N<3 is noise → those spectra go to the collapsed cell).

### 8.8 Grid v2 (Nao, 2026-10-24 night): 8 EW bins, 8 oxygen bins, 8 BAL bins — supersedes 8.7's bin counts
Axes 1–2 unchanged (39 z slices; M_i(z=2) 1-mag). All EW ladders are factor-2 (log2) bins like the galaxy grid.
| axis | 8 bins (rest EW, Å) | line by z |
|---|---|---|
| 3. primary line | Hα: <12.5, 12.5–25, 25–50, 50–100, 100–200, 200–400, 400–800, >800 | z<0.4 |
| | Hβ: <3, 3–6, 6–12, 12–25, 25–50, 50–100, 100–200, >200 | 0.4–0.9 |
| | Mg II: <2, 2–4, 4–8, 8–15, 15–30, 30–60, 60–120, >120 | 0.9–2.3 |
| | C IV: <1.5, 1.5–3, 3–6, 6–12, 12–25, 25–50, 50–100, >100 | 2.3–5 |
| 4. oxygen line | one ladder for all: <0.5, 0.5–1, 1–2, 2–4, 4–8, 8–16, 16–32, >32 | [O III] 5007 z<0.9; [O II] 3727 0.9–1.6; O I 1304 1.8–2.6 (O III] 1663 alt. 1.2–2.6); O VI 1034 z>2.6 |
| 5. BAL | 1 non-BAL (AI_CIV=0, no Mg II/Al III trough); weak HiBAL (AI>0, BI=0) × VMIN_CIV 0–5000 / 5000–15000 / >15000 km/s (3); strong HiBAL (BI>0) × same 3 velocity bins (3); LoBAL (Mg II or Al III trough, any C IV) (1) | 8 bins for z≥1.5; z<1.5 only non-BAL / LoBAL (2) |
Estimated stacks (floor 300; EW fractions log-normal σ=0.3 dex primary → 0.1/2.1/13.6/34.1/34.1/13.6/2.1/0.1 %,
σ=0.5 dex oxygen → 3.6/7.9/15.9/22.6/22.6/15.9/7.9/3.6 %; BAL 84/2.1/2.8/2.1/2.4/3.2/2.4/1.0 %):
z×L 165 → ×EW8 552 (1.18 M quasars inside) → ×O8 1168 → ×BAL8 1072 (z×L×EW8×BAL8 without oxygen: 522).
By range (×EW8 / ×EW8×O8 / all four): z<0.4 50/92/92; 0.4–0.9 82/188/188; 0.9–1.6 132/364/364; 1.6–2.3 126/308/264;
2.3–3.5 142/212/160; 3.5–5 20/4/4. The outer EW bins (<2 % each) and the BAL bins only populate at 1<z<3 in the
M_i −28…−25 bins; elsewhere they collapse to the neighbouring bin (galaxy-grid rule). Nominal grid 39×12×8×8×8 =
239,616 cells; ≈1,100 stacks expected.

### 8.9 Grid v3 (Nao, 2026-10-24 night): + continuum slope axis — supersedes 8.8
Axis 6. continuum slope α_λ (F_λ ∝ λ^α, pass Q `ALPHA_LAM`; the quasar equivalent of the galaxy color axis): 8 bins
of 0.25: < −2.25, −2.25…−2.0, −2.0…−1.75, −1.75…−1.5, −1.5…−1.25, −1.25…−1.0, −1.0…−0.75, > −0.75
(assumed N(−1.6, 0.4) → 5/11/20/25/21/12/5/2 %; red tail = dust-reddened quasars). Alternative: observed g−r / r−z color
in 0.25-mag bins as for galaxies — equivalent but z-dependent; slope is the rest-frame quantity, preferred.
Full grid: 39 z × 12 M_i × 8 EW × 8 slope × 8 oxygen × 8 BAL = 1.9 M nominal cells.
Estimated stacks at floor 300 (same assumed fractions as 8.8): z×L 165 → ×EW8 552 → ×slope8 1204 → ×O8 336 →
×BAL8 268 (all six). Without the oxygen axis: 1090. By range (+EW / +slope / +O / all six): z<0.4 50/104/24/24;
0.4–0.9 82/196/44/44; 0.9–1.6 132/358/160/144; 1.6–2.3 126/322/104/56; 2.3–3.5 142/216/4/0; 3.5–5 20/8/0/0.
Reading: the sample supports three fine axes at a time (z, L, EW, slope → 1,204 stacks), not six; the product of six
8-bin axes only populates at 0.9<z<2.3 in the −27…−25 bins. Rule, as for galaxies ("measured on 3 axes"): a cell is
defined on all six axes; when a cell is below the floor it is merged along the axis with the least spread in the
stack (BAL first, then oxygen, then slope) until ≥300. Expected: ≈1,200 stacks of which ≈270 are full six-axis cells.

### 8.10 Bin reduction options (2026-10-24 night; same assumed fractions, floor 300)
Number of stacks with ≥300 quasars (S/N>5) for each redshift slicing (rows) and binning of the four physical axes
(columns; bins of primary-line EW / slope / oxygen EW / BAL class, all on top of z × M_i):
| redshift slices | z × M_i only | + EW 8, slope 8, oxygen 8, BAL 8 | + EW 6, slope 6, oxygen 4, BAL 4 | + EW 6, slope 4, oxygen 3, BAL 4 | + EW 6, slope 6 (no oxygen, no BAL) |
|---|---|---|---|---|---|
| 39 (Δz 0.1 to 3.5, 0.25, 0.5) | 165 | 268 | 756 | 1104 | 1244 |
| 24 (Δz 0.1 to 1, 0.2 to 3, 0.5 to 5) | 104 | 576 | 960 | 1200 | 1040 |
| 18 (Δz 0.2 to 3, 0.5 to 5) | 86 | 692 | 1052 | 1200 | 942 |
Reduced ladders: EW6 = merge the two outer bins on each side (<¼, ¼–½, ½–1, 1–2, 2–4, >4 × composite);
slope6 = < −2.0, −2.0…−1.75, −1.75…−1.5, −1.5…−1.25, −1.25…−1.0, > −1.0 (slope4: < −1.9, −1.9…−1.6, −1.6…−1.3, > −1.3);
O4 = <1, 1–4, 4–16, >16 Å (O3 = <2, 2–8, >8); BAL4 = non-BAL / weak HiBAL / strong HiBAL / LoBAL (no velocity split).
Recommendation: z24 with EW6·slope6·O4·BAL4 → ≈960 stacks, all six axes kept, full cells at 0.4<z<2.3; z>3.5 by
individuals + z×L×EW stacks only. Δz=0.2 at 1<z<3 is 3–7 % in (1+z), one pixel shift of 69 km/s is 2.3e-4 in (1+z),
so a slice's rest coverage differs by ≤ 100 Å at the edges — harmless for a rest-frame PCA.

### 8.11 Steps, as settled with Nao (2026-10-24 night)
Redshift rule (Nao): de-redshift with the **peak of a narrow/low-ionization line when available, not the whole-spectrum
(Redrock) redshift**: [O III] 5007 peak (z<0.9) → Mg II peak (0.9–2.3) → O I 1304 peak when detected ≥3σ (z>2.3) →
pipeline z only as fallback. Implemented in pass Q as `Z_SYS`, `Z_SYS_SRC` (1/2/3/0) from the new `_VPEAK` columns
(peak of the smoothed continuum-subtracted line, parabolic refinement). Test coadd: 68 Mg II, 15 [O III], 7 O I,
22 fallback; Z_SYS − Z_pipe median −34 km/s, MAD 410 km/s; peak − centroid: [O III] 10 km/s, Mg II 190, C IV 265.
1. Catalogs (now): DR16Q, DR1 zlya, DR1 QSO cat, EDR BAL, BH-mass VACs → `~/data/sdssdesi_products/quasar/`.
2. Pass Q over all DESI coadds (Sunday, after pass B; 2 workers, ~13 h), then the SDSS spPlates (internal disk).
   Per quasar: pass-B columns + slope, line EWs/peaks (incl. oxygen lines), Z_SYS, forest, BAL proxy.
3. Master quasar table: pass Q + catalog BAL class (DR16Q BAL_PROB/BI/AI, zlya AI/BI; proxy where neither) + Z_SYS;
   S/N>5 selection; 10 % hold-out.
4. Bin edges from the measured distributions (factor-2 EW ladders centred on the per-slice median; slope 0.25 bins);
   cell labels on all six axes (z, M_i, EW, slope, oxygen, BAL); fine bins only where the line is detected ≥3σ.
5. Stacks (pass C): rest frame with Z_SYS (integer-pixel shift), forest pixels <1216 Å masked, normalization window
   per slice, weighted mean + one outlier pass, merge under-filled cells (BAL → oxygen → slope) to floor 300;
   24 z slices; ≈1,000 stacks; 2D FITS per slice (sdsslib format) + cell table.
6. Weighted EM-PCA on stacks + S/N>5 individuals (also de-redshifted with Z_SYS), lines+continuum, rest 600–9800 Å.
7. Fit model: quasar components × mean IGM transmission + galaxy components (z<1) + velocity width.
8. Validation: hold-out and BAL subsets vs Redrock QSO v1.1 / QuasarNET / DR16Q Z_PCA / EDR VI truth; line-confusion
   rates; C IV blueshift distribution vs luminosity as a physics check.
9. ms.tex: §3.7 (Margala rewording, line list, grid, Z_SYS rule), §3.9 (pass Q), table of stacks.

### 8.12 Own BAL measurement (Nao, 2026-10-24 night: measure BALs ourselves; catalogs only as cross-check)
Pass Q now measures, on the smoothed (5 px ≈ 350 km/s) flux normalized by the power-law continuum, requiring the
smoothed flux to be ≥2σ below 0.9 C: **AI** (Hall+02; troughs ≥450 km/s, −25,000<v<0) for C IV, Si IV, Mg II, Al III;
**BI** (Weymann+91; ≥2000 km/s, −25,000<v<−3000) for C IV and Si IV; NTROUGH, VMIN/VMAX of the deepest trough.
`BAL_CLASS`: 0 non-BAL, 1 weak HiBAL (AI>0, BI=0), 2 strong HiBAL (BI>0), 3 LoBAL (Mg II/Al III trough together with
a C IV/Si IV trough; at z<1.5, where C IV is not covered, only AI_MGII>1000 km/s), −1 no BAL region covered.
Test coadd (112 quasars): 11 % flagged (9 weak, 1 strong BI=4800 km/s, 1–3 LoBAL), i.e. the expected 10–15 %; the
1σ version flagged 65 % (noise), the 2σ test fixed it. Known limits: (i) the global power law is a poor continuum
under Mg II (Fe II pseudo-continuum) and at the blue edge → LoBAL and Al III are the least reliable; (ii) at S/N<2
only BI ≳ 1000 km/s troughs are detectable. Planned iteration: after the version-0 PCA, redo the trough search with
the non-BAL eigenspectra as the continuum (as balfinder / Guo & Martini 2019 do), which also makes the low-S/N
classification possible. DR16Q and DR1 zlya remain useful as a cross-check (optional download).

### 8.13 Downloaded 2026-10-24 20:37 (Nao approved the first two catalogs), under ~/sdssredux
| file | size | content |
|---|---|---|
| `~/sdssredux/DR16Q/DR16Q_v4.fits` | 1.72 GB, 750,414 rows, 183 cols | BAL_PROB, BI/AI_CIV, BI/AI_SIIV, Z_PCA, Z_MGII, Z_CIV, Z_HBETA, Z_HALPHA, Z_CIII, Z_LYA (+ZWARN_*), PLATE/MJD/FIBERID |
| `~/sdssredux/DESI_DR1_zlya/qso_cat_dr1_main_dark_healpix_zlya-v0.fits` | 281 MB, 856,767 rows, 24 cols | TARGETID, Z, ZMASK, BI_CIV, AI_CIV (+ERR, NCIV, VMIN/VMAX); sha256 verified |

### 8.14 Step 1 (spectrophotometric scale to the fiber flux) for quasars — unchanged, and why only the scale
Pass Q runs pass B's `measure_chunk`, so every quasar gets Step 1 exactly as galaxies do: spectrum × MEAN_PSF_TO_FIBER_SPECFLUX,
synthetic g/r vs FIBERFLUX_G/R (DECam or BASS/MzLS by PHOTSYS), scale CAL_A applied, tilt CAL_B_EST stored but not
applied; then Step 2 (F99) and Step 3 (rest frame with Z_SYS). SDSS: same against FIBER2MAG/FIBERMAG.
Test coadd: CAL_A median 1.00, scatter 0.22 (stars: 0.05) — the extra scatter is quasar variability between the LS
imaging epoch (2014–2019) and the DESI spectrum (2021–2024), typically 0.1–0.2 mag. Consequences: (i) the scale ties
the spectrum to the imaging-epoch flux, harmless for the shape (every spectrum is renormalized before stacking/PCA)
but it means M_i is the imaging-epoch luminosity; (ii) a per-object color tilt from imaging would be wrong for a
variable source (quasars are bluer when brighter), so the tilt stays unapplied for quasars — the per-fiber tilt is to
be solved from stars per exposure/petal, as for galaxies. The BOSS offset-fiber problem is already corrected in the
DR17 pipeline (8.2); the Gaia XP check of bright quasars remains the absolute test.

### 8.15 The quasar steps, complete (2026-10-24 night; supersedes 8.11 as the reference list)
A. Catalogs (done): DR16Q and DR1 zlya under ~/sdssredux (8.13) — cross-checks only.
B. Per-spectrum pass (pass Q, `05_passQ_measure.py`; DESI Sunday after pass B, then SDSS spPlates):
   1. spectrophotometric scale to the fiber flux (×MEAN_PSF_TO_FIBER_SPECFLUX, g/r vs FIBERFLUX; scale only; SDSS vs FIBER2MAG);
   2. Milky Way extinction (F99, R_V=3.1, 0.86×E(B−V)_SFD) on spectrum and photometry;
   3. observed-frame measurements: S/N per camera/total, synthetic mags, Dn4000-type indices;
   4. first rest frame with the pipeline z (log grid 1e-4, integer shift, F×(1+z)): power-law slope α_λ, broad lines
      (Lyα+NV, OVI, OI 1304, SiIV, CIV, OIII] 1663, AlIII+CIII], MgII, Hβ, Hα) and narrow lines ([OIII], [OI], [OII] via
      pass B): flux, EW, σ, centroid and peak velocities; forest coverage/transmission; AI/BI troughs for CIV, SiIV,
      MgII, AlIII → BAL_CLASS;
   5. systemic redshift Z_SYS = peak of [OIII] 5007 (z<0.9) → MgII (0.9–2.3) → OI 1304 (z>2.3, ≥3σ) → pipeline z;
   6. absolute magnitude: flat ΛCDM h=0.7 Ωm=0.3; M_i(z=2) = i_dered − DM(z) − K_Richards(z) (i from LS r,z or
      synthetic i — decide); pass-B M_R/G/U/NUV kept.
C. Master table + training set: pass Q ∪ catalogs; ZWARN=0, S/N>5 (≈1.2 M incl. SDSS); 10 % hold-out.
D. Cell labels (six axes): z (24 slices), M_i (1 mag), primary-line EW (6 log2 bins; Hα/Hβ/MgII/CIV by z), slope
   (6 × 0.25), oxygen EW (4; [OIII]/[OII]/OI/OVI by z), BAL class (4). Edges from measured per-slice medians; fine
   bins only for lines detected ≥3σ.
E. Stacks (pass C): re-shift every spectrum to Z_SYS; forest (<1216 Å rest) masked; per-slice normalization window;
   weighted mean Σwf/Σw (w = ivar·s²) + median/scatter/N; one χ² outlier pass; merge under-filled cells
   (BAL → oxygen → slope) to floor 300 → ≈1,000 stacks; 2D FITS per z slice (sdsslib) + cell table.
F. Eigenspectra: weighted EM-PCA on stacks + S/N>5 individuals (at Z_SYS), lines+continuum, rest 600–9800 Å;
   number of components from the hold-out residual.
G. BAL iteration: trough search redone with the non-BAL eigenspectra as continuum (fixes MgII/Fe II and low S/N);
   re-label the BAL axis, rebuild affected stacks, PCA once more.
H. Fit model: quasar eigenspectra × mean IGM transmission (Inoue+14/Kamble+20) + galaxy eigenspectra (z<1)
   + velocity width; forest down-weighted.
I. Validation: hold-out and BAL subsets vs Redrock QSO v1.1, QuasarNET, DR16Q Z_PCA/Z_MGII, zlya ZMASK, EDR VI truth;
   catastrophic-failure and line-confusion rates; CIV blueshift vs luminosity; Gaia XP absolute-flux check of the
   bright SDSS quasars.
J. ms.tex: §3.7 (Margala rewording, line list, grid, Z_SYS rule, own BAL measurement), §3.9 (pass Q), stack table.

### 8.16 ms.tex updated (2026-10-24 night; Nao commits in ~/github/papers/26_sdssdesi_spectra)
§3.7 (sec:qsotemplates): item 2 (BALs: own AI/BI measurement, BAL class as a grid axis, catalogs as cross-check),
item 3 (systemic z from the line peak: [O III] / Mg II / O I 1304, pipeline z only as fallback), item 6 (Margala
correction already in DR14+ pipeline; Gaia XP as the check); the "Plan" paragraph replaced by "Processing steps"
(steps 1–6 incl. fiber-flux scale, F99 extinction, rest frame, M_i(z=2) with flat ΛCDM + Richards K), the training
set, stacking/merge rule, EM-PCA, fit model, validation) and a "Grid" paragraph + Table tab:qsogrid (6 axes, 24×12×6×6×4×4
= 165,888 nominal cells, 104 populated (z,L) cells, ≈960 expected). New bibitems: Hall2002, Lyke2020, Richards2006,
VandenBerk2001, Weymann1991, DESI2025zlya. Not compiled here (no pdflatex on this Mac) — braces/environments/
citations/refs checked by script; please compile on Overleaf. Backup of the pre-edit file: scratchpad
ms.tex.bak_20261024 (session-only).

### 8.17 Scripts moved to `programs_quasars/`; SDSS first (Nao, 2026-10-24 night)
`programs_quasars/qsolib.py` (shared measurements), `01_passQ_sdss.py` (SDSS spPlates, **running since 2026-10-24
~21:30**, log `~/data/sdssdesi_status/passQ_sdss.log`, ~3 h), `02_passQ_desi.py` (DESI; after pass B). README there.
SDSS findings: the pipeline calibrates spectra to PSF flux → Step 1 compares with PSFMAG (DR17; CAL_A median 1.00,
MAD 0.24 = variability), DR8 has no PSF mags in specObj → no scale (CAL_MODE 0; fiber offsets stored). All 3,933 DR17
and 2,699 DR8 plate-MJDs with primary quasars are on disk. Slopes at S/N>5: DR17 −1.32, DR8 −1.44, DESI −1.60 (one
coadd) — check with full samples; if DR17 stays bluer, the BOSS offset-fiber correction (8.2) needs the Gaia XP test.
Star work branched to its own thread: `programs_stars/HANDOVER_stars.md`.

### 8.18 Monitor (2026-10-24 21:03): `programs_monitor/monitor.py` replaces 04_monitor.py on port 8765
One page for all tracks and steps (pass B, pass Q SDSS/DESI, cell tables, stacks, EM-PCA, star steps) with gauge,
rate, ETA/finish time, per-step worker CPU/RSS, system CPU/memory/swap/disk gauges and a 1-h history; same guard and
status dir. New steps report through `~/data/sdssdesi_status/step_<name>.json` (README there). Pass B ETA slipped to
Sunday ~23:00 under the evening load (Prospector + pass B + SDSS pass); SDSS pass reniced to 19.

### 8.19 Products directory renamed (2026-10-24 22:22): `~/data/sdssdesi_products/` (was `sdssdesi_cache`)
Nao: these are products, not caches (the zall Parquet is a 51 GB read, the master table an hour, the stacks a day).
A symlink `~/data/sdssdesi_cache -> sdssdesi_products` keeps older scripts working (the star thread's scripts still
use the old name; the running pass-B workers loaded the lookup at start and a restart by the autoscaler resolves
through the symlink). Scripts and docs of the templates / quasars tracks now use the new name. Everything under
`~/data` is on the internal SSD and survives reboots; only the session scratchpad (/private/tmp) is lost.
Quasar products: `~/data/sdssdesi_products/quasar/` (master tables, checks, EW quantiles, later cell tables).
Cleanup (Nao, 2026-10-24): pass-Q test directories deleted. **To do after pass B finishes:** delete
`~/data/sdssdesi_passB_test`, `~/data/sdssdesi_passB_test2` and `~/data/sdssdesi_products/zlookup_main.parquet`
(rebuildable from zall_loa.parquet; pass B reads only zlookup_compact.npz).

### 8.20 SDSS master table v1 and checks (2026-10-24 23:05; `03_qso_table_sdss.py`, full report in
`data/quasar/qso_checks_sdss_v1_airwave.txt`, quantiles `data/quasar/qso_ew_quantiles_sdss_v1.csv`)
Rows 1,140,616 (DR17 978k, DR8 163k); after cuts (ZWARN=0, plate good, 0<z<5) 748,305; S/N>5: 302,883 (DR17 236k,
DR8 67k; ms.tex census said 233k + 67k ✓). DR16Q match 727k.
**Bug found and fixed in pass Q**: the optical lines (Hα, Hβ, [O III], [O II], [O I], …) were at AIR wavelengths while
SDSS/DESI spectra are vacuum → our [O III] peak sat +76 km/s (= 1.4 Å at 5007) above DR16Q Z_PCA with MAD 74 km/s.
qsolib now overrides pass B's LINES with vacuum values and uses vacuum Hβ/Hα/[O III]/[O I]; SDSS pass rerunning
(started 23:05, ~1 h at nice 19), table to be rebuilt. **The same air-wavelength list is in `02_passB_measure.py`
(galaxy lines of pass B, running now)** — for the galaxy thread: narrow lines fitted 80 km/s off centre with the
σ=80 km/s Gaussian lose flux; fix the LINES table (vacuum) and either re-measure lines in the cell-assignment pass
or rerun pass B's line step. Decision for the galaxy thread / Nao.
Findings (S/N>5, before the fix; offsets will shift by the 76–84 km/s correction):
- Z_SYS: [O III] peak vs DR16Q Z_PCA median +76, MAD 74 km/s, 0.2 % > 1000 km/s (excellent); Mg II peak vs DR16Q
  Z_MGII +212 / MAD 315 (DR16Q uses 2799.12? we use 2798.75 = 39 km/s of it; the rest is the redvsblue "peak" vs
  ours), vs Z_PCA +150 / 353, 11 % > 1000; **O I 1304 peak is poor** (MAD 600–900 km/s, 15–30 % > 1000 km/s) → drop
  O I from the Z_SYS chain (keep the pipeline z at z>2.3, or C III]/Mg II where in window) — Nao to confirm.
  Pipeline Z vs DR16Q Z_PCA: −7 / 90 km/s.
- BAL (DR17, 1.57<z<5, S/N>5, 120k): our class>0 catches 96 % of DR16Q BI_CIV>0 (67 % as class 2); false positives
  14.5 % among BAL_PROB<0.1; our BAL fraction 29 % (S/N 5–10) / 39 % (S/N>10) vs DR16Q BAL_PROB≥0.5 31 % / 33 %;
  agreement 80 %. Usable as a grid axis now; the eigenspectra-continuum iteration should cut the false positives.
- Slope: DR17 and DR8 agree at z>1.2 (α_λ ≈ −1.30 to −1.40, both); at 0.3<z<1 DR17 is bluer by 0.3–0.5 (different
  programs/hosts, not calibration). The DESI test coadd gave −1.60 → compare with the full DESI pass.
- Measured EW medians (S/N>5, ≥3σ): Hα 80–176 Å (rising with z = luminosity), Hβ 35–38, Mg II 26, C IV 22–25
  (composite 24 ✓); [O III] 12–16, [O II] 2.4–3.2, O I 1304 1.6–2.1, O VI 10–14; slope 16/50/84 % = −1.7/−1.35/−0.9
  at z>0.9 → the slope bins −2.0…−1.0 (6 × 0.25) are well placed; the factor-2 EW ladders centred on these medians.
  **Oxygen-axis gap**: at 1.6<z<1.8 neither [O II] (obs >9700 Å) nor O I 1304 (obs <3650 Å) is in the window
  (167 detections) → use O III] 1663 there (obs 4300–4700 Å).
- BAL fraction (our class) 0.22–0.39 at z>1.6; 0 at z<1.5 except LoBAL ≈ 0.3–6 % at z<0.6 (host/Fe II false
  positives likely — check).

### 8.21 SDSS pass Q v2 (vacuum wavelengths) and table rebuilt 2026-10-25 00:03 — `data/quasar/qso_checks_sdss_v2_vacuum.txt`
[O III] peak vs DR16Q Z_PCA: median −5 km/s, MAD 77 (was +76): the air/vacuum bug is confirmed fixed. Mg II unchanged
(+212 / 315 km/s vs DR16Q Z_MGII; UV lines were already vacuum). Master table: `~/data/sdssdesi_products/quasar/
qso_master_sdss.parquet` (1,140,616 rows, 168 cols, 506 MB). Grid quantiles: `data/quasar/qso_ew_quantiles_sdss.csv`.
Machine 2026-10-25 02:10: Prospector finished (load 1); pass B restarted with 4 workers at 02:11 (resume.log), but
**exdisk2 is on a 480 Mb/s (USB 2) link** (system_profiler; disk read 30–35 MB/s): pass B is disk-bound at ≈13–15
files/min whatever the worker count (ETA Mon 00:40). Re-plugging exdisk2 on a USB 3 port (5 Gb/s) is the only real
acceleration (the autoscaler then allows 4 workers; ≈100 files/min → hours instead of a day). Pass Q DESI shares that
disk → it stays queued behind pass B; the quasar track continues on the SDSS half (cells, stacks) meanwhile.

### 8.22 Cells on the SDSS half (2026-10-25 17:15; `04_cells.py --survey sdss`; tables `data/quasar/qso_stack_cells_sdss.csv`,
`qso_primary_ew_medians_sdss.csv`; per-object labels `~/data/sdssdesi_products/quasar/qso_cells_sdss.parquet`)
Training 302,840 (QCUT, S/N>5, M_i). Six-axis cells occupied 15,385; ≥300 as full cells: 152 (82.5k quasars). After the
merge rule (BAL → oxygen → slope): **537 stacks holding 265,598 quasars (87.7 %)**; 37k quasars in 817 pool cells still
below 300 (mostly z<0.4, z>3.5 and the bright/faint M_i edges). Stages: 152 full, 18 BAL-collapsed, 209 +oxygen,
158 +slope. Fine bins available: EW 95 %, slope 100 %, oxygen 50 % (detected ≥3σ), BAL known 94 %.
Per slice: 7–10 stacks at z<0.6, 18–26 at 0.6–1.0, 43–49 at 1.0–2.4 (12–25 full six-axis cells each), 10–36 at 2.4–3.5,
4 at 3.5–4, none above 4. DESI (≈940k at S/N>5) will roughly triple this → the ≈1,000–1,200 expected in 8.10.
Coordination note for all threads: `docs/COORDINATION.md` (2026-10-25); messages sent to the Galaxy and Stellar sessions.

### 8.23 Internal disk nearly full (2026-10-25 18:00): 60 GB free of 1.8 TB (97 %)
Still to land tonight: pass B +≈7 GB, pass Q DESI +≈5 GB, pass S +≈3 GB, EXP_FIBERMAP dump ≈1.5 GB → fits, but the
stacks (pass C, quasar stacks) and any rest-frame spectrum store will not. Candidates to free (Nao decides; none deleted):
`~/github/projects_desi/programs/zall-pix-loa.csv` 10 GB (zall extract, superseded by products/zall_loa.parquet);
`~/github/projects_desi/csvfiles/` 7.6 GB (desidr2star float32 csvs — duplicates of projects_gaia/csvfiles 12 GB?);
`~/desiredux/zall-pix-loa.fits` 47 GB (the source; products/zall_loa.parquet holds the used columns — keep unless desperate);
`~/data/exdisk2_keep_tmp` 55 GB, `~/data/snana` 63 GB, `~/data/HSC` 167 GB (other projects); products/zlookup_main.parquet
1.7 GB (approved, deleted after pass B). The 1–2 TB SSD for rest-frame spectra (programs_templates/README.md) is now needed.
`02_passQ_desi.py --with-stars --expfibermap` (Stellar hook `write_star_file` + EXP_FIBERMAP dump, all rows) tested on one
coadd: 112 QSO, 12 STAR, 784 exposure rows (44 KB); `after_passB.sh` relaunched with both flags.
