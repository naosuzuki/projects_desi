# programs_stars — stellar templates (pass S) runbook

Steps (Nao, 2026-10-24 night): 1 catalog with Gaia parallax → 2 spectrophotometric calibration (2b joint fit) → 3 RV and rest frame → 3b dust map /
Milky Way extinction → 4 parallax distance → absolute flux → 5 HR grid, test that the stars of a cell are alike (report shared cells) →
6 stack the alike cells → 7 EM-PCA (normalize to Gaia-G flux = 1, factor stored). SDSS first (internal disk), DESI after pass B. Details: `docs/HANDOVER_stars.md` §3.

Design and status: `docs/HANDOVER_stars.md` (steps, census, HR-cell counts, decisions); paper: ms.tex §3.5 ("Stars: combining Gaia
and DESI"). Imports the common pieces from `../programs_templates/02_passB_measure.py` (filters, rebin, extinction, Step 1).

| script | step | what | status |
|---|---|---|---|
| `step1_desi_catalog.py` | 1 | main-survey STAR rows of zall-pix-loa with SUBTYPE and Gaia columns → `~/data/sdssdesi_cache/zall_stars_main.parquet` (12.18 M) + census | done |
| `step1_desi_hr.py` | 1 | Gaia×DESI quality sample, M_G0 / (BP−RP)_0 / [M/H], HR-cell counts → `gaia_desi_stars_hr.parquet` (3.92 M) | done |
| `step1_sdss_catalog.py` | 1 | spAll DR17 + specObj DR8 primary stars joined with the Gaia crossmatch (nearest per plate/mjd/fiber), ELODIE, PSF mags, quality flag → `sdss_stars_gaia.parquet` (1.06 M rows, 780 k with parallax) | done |
| `step2_3_pass_sdss.py` | 2–3 | spPlates, catalog stars with a Gaia parallax: scale to PSF g,r (quasar thread's `sdss_scale`), Gaia G/BP/RP check, RV frame (no rebin), indices, distance/extinction/absolute-flux columns → `~/data/sdssdesi_passS_sdss/{dr17,dr8}/spPlate-*.parquet` (0.15 s/plate, ≈25 min for 8,084 plates) | **running** 2026-10-24 22:00 |
| `step2_3_pass_desi.py` | 2–3 | DESI coadds, STAR rows: pass-B chain + star lookup + starlib → `~/data/sdssdesi_passS/{dark,bright}/hpx<N>.parquet`; run after pass B, 2 workers (or inside the pass-Q walk) | tested on 1 coadd |
| `step2b_spcal.py` | 2b | **joint spectrophotometric fit** (Nao, 2026-10-24 night) on the pass tables: per star a gray scale A_MAG and a tilt (mag/dex in log λ), per band a zero point and a color term, Gaia G/BP/RP fixed as the reference system, SDSS PSF u g r i z (u dropped: 85 % extrapolated) / LS g r z (DECam S, BASS+MzLS N) constrain the tilt and get their zero points from the stars; alternating weighted LSQ, 3σ clipping → `~/data/sdssdesi_cache/stars/spcal_<survey>.parquet`, `zp_<survey>.csv`. SDSS filters = Doi+2010 (speclite `sdss2010`, airmass 1.3). step 6 applies A_FLUX and the tilt where TILT_ERR < 0.1 mag/dex. DESI test coadd: scale +0.076 ± 0.05, tilt −0.25 mag/dex (= the +0.045 BP−RP offset), tilt scatter beyond errors 0.22 mag/dex; LS zero points vs Gaia g −0.07, r −0.06, z −0.02 (one coadd) | run on SDSS: 714 k stars, ZP g −0.042 r −0.025 i −0.037 z −0.020 vs Gaia, tilt scatter beyond errors 0.24 mag/dex |
| `step3b_ebv3d.py` | 3b | 3D E(B−V) at the parallax distance: Edenhofer+24 (< 1.25 kpc) + Bayestar19, maps in `~/data/GAIA`; SDSS catalog columns EBV_3D/EBV_B19/EBV_E24/EBV_3D_SRC; DESI → `desi_stars_ebv3d.parquet` (10.3 M). Units calibrated on the data: E(B−V)_SF11 = 0.86 E_B19 (matches 0.86 SFD at > 3 kpc to 0.2 %), = 1.00 E_E24 (matches Bayestar in the overlap to 0.6 %). Nearby stars (< 300 pc) have 56 % of the SFD column in front | done (1.5 min) |
| `step5_hr_cells_check.py` | 5 | HR grid on the pass tables; per-cell similarity test on in-cell residuals (plane fit in M_G0, color; tilt-corrected colors): index spread relative to the median cell, bimodality, absolute-flux scatter, Hα emitters, [M/H] spread → `~/data/sdssdesi_cache/stars/cells_<tag>.csv`, members; report of the cells shared by different groups with the split axis. SDSS: 216 cells, 38 homogeneous, 178 mixed (docs/HANDOVER_stars.md §5) | run on SDSS 2026-10-24 |
| `step6_stack.py` | 6 | streaming two-pass stacks of the alike cells (Gaia-G scale, RV shift, 3D dereddening, × (d/10 pc)²), mixed cells split on their axis → `~/data/sdssdesi_stars/stacks_<tag>.fits` (FLUX/ERR/SCATTER/NSPEC + CELLS table, sdsslib log-grid header) | run on SDSS 2026-10-24: 316 stacks from 216 cells, 91,726 stars, 11 min; `~/data/sdssdesi_stars/stacks_sdss_0.5x0.1_gspphot.fits` + quick-look PNG (docs/HANDOVER_stars.md §6) |
| `check_pass.py` | — | Gaia / PSF / fiber checks, RV, indices per subtype on any pass output | written |
| (monitor) | the stars track with its gauges and ETAs is on the main monitor http://127.0.0.1:8765/monitor.html (`programs_monitor/monitor.py`), sections separated per track; the separate 8767 page was retired 2026-10-24 night | |
| `run_steps.sh` | — | the steps in order: `./run_steps.sh 1 \| 2-3sdss \| 3b \| 2-3desi \| 5 \| 6 \| monitor` | |
| `starlib.py` | lib | Gaia synthetic photometry (Vega, Riello+21 zero points), RV frame, Lick/Cenarro indices, Reid+95 ratios, colors, distance/extinction/absolute-flux columns | |

Run: `./run_steps.sh <step>` (see the table); everything but the DESI pass is internal-disk work and safe while pass B runs.
Test 2026-10-24 night (coadd-main-bright-0, 381 stars, 10 s + 9 s read):
- Gaia check (catalog − synthetic, Vega): G +0.06 (σ68 0.03), BP +0.075, RP +0.025, BP−RP +0.045, flat in S/N and SUBTYPE →
  the pipeline flux is ≈6 % brighter than Gaia G (same sense as the +0.03–0.05 offset vs the LS fiber fluxes in pass B) and the
  spectrum is ≈0.045 mag bluer than Gaia BP−RP; the band-weight fraction in the flat padding is 1–2 %. One coadd: confirm at scale,
  check the BP Vega zero point and the blue-end calibration before using BP−RP as the tilt.
- The Gaia G residual and the LS fiber r residual of the same star correlate at 0.80 (n=340): the per-star scale error is in the spectrum, and Gaia G (σ 0.03) measures it better than the fiber flux (σ 0.05).
- Step 1 (pass B) stars: g +0.042, r +0.019 (σ68 0.04–0.05), z after scale −0.007 ± 0.046, CAL_OK 0.90 — as in the 20-coadd test.
- RV median +10 km/s, σ68 48 km/s, none beyond ±300; integer-shift residual σ68 23 km/s.
- Indices behave: A stars HδA 10.9 / Hβ 7.4 / Hα 5.4 Å; G 8.6 Ca K, 4.0 G band; K 10.2 Ca K, 5.4 Mg b, 4.3 Na D; M dwarfs TiO5 0.58,
  CaH2 0.56, Na D 11 Å, 7 of 378 stars with Hα in emission (>3σ, active M dwarfs). Coverage 97–99 % of good stars.

Rules: never write to the external disks; pass S over the coadds only after pass B (disk-bound), preferably inside the pass-Q walk.
