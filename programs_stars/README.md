# programs_stars — stellar templates (pass S) runbook

Steps (Nao, 2026-10-24 night): 1 catalog with Gaia parallax → 2 spectrophotometric calibration → 3 RV and rest frame → 3b dust map /
Milky Way extinction → 4 parallax distance → absolute flux → 5 HR grid, test that the stars of a cell are alike (report shared cells) →
6 stack the alike cells (**goal for now**) → 7 EM-PCA. SDSS first (internal disk), DESI after pass B. Details: `docs/HANDOVER_stars.md` §3.

Design and status: `docs/HANDOVER_stars.md` (steps, census, HR-cell counts, decisions); paper: ms.tex §3.5 ("Stars: combining Gaia
and DESI"). Imports the common pieces from `../programs_templates/02_passB_measure.py` (filters, rebin, extinction, Step 1).

| script | what | status |
|---|---|---|
| `00_build_star_lookup.py` | main-survey STAR rows of zall-pix-loa with SUBTYPE and Gaia columns → `~/data/sdssdesi_cache/zall_stars_main.parquet` (12.18 M rows, 2 min, internal disk) + census | done 2026-10-24 |
| `01_hr_cells.py` | Gaia×DESI quality sample, M_G0 / (BP−RP)_0 / [M/H], HR-cell counts → `~/data/sdssdesi_cache/gaia_desi_stars_hr.parquet` (3.92 M) | done 2026-10-24 |
| `04_sdss_star_catalog.py` | **step 1 (SDSS)**: primary stars of spAll DR17 + specObj DR8 joined with the Gaia DR3 crossmatch of projects_gaia (nearest match per plate/mjd/fiber), ELODIE parameters, PSF mags (DR8 from the crossmatch), quality flag, M_G0/(BP−RP)_0 → `~/data/sdssdesi_cache/sdss_stars_gaia.parquet` (1.06 M rows, 780 k with parallax, 25 s) | done 2026-10-24 |
| `05_passS_sdss.py` | **steps 2–3 (SDSS)**: spPlate rows of the catalog with a Gaia parallax; scale to PSF g,r (quasar thread's `sdss_scale`), Gaia G/BP/RP check, RV frame (no rebin: spPlates are on the log grid), indices, distance/extinction/absolute-flux columns → `~/data/sdssdesi_passS_sdss/{dr17,dr8}/spPlate-PPPP-MMMMM.parquet` | to write |
| `06_hr_cells_check.py` | **step 5**: HR grid on the pass-S tables; per-cell spectral homogeneity (index scatter vs errors, bimodality, Hα emitters, [Fe/H] spread); report of cells shared by different stellar groups | to write |
| `07_stack_stars.py` | **step 6**: RV-shifted, dereddened, 10-pc stacks of the homogeneous cells; mixed cells split first | to write |
| `starlib.py` | measurements on the RV-corrected 3600–9824 Å grid: synthetic Gaia G/BP/RP (Vega, Riello+21 zero points) vs catalog, Lick/Cenarro-type indices (Ca K, HδA, G4300, HγA, Hβ, Mg b, Fe5270/5335, He II 4686, He I 5876, Na D, Hα, Ca II triplet), Reid+95 ratios TiO5/CaH2/CaH3, colors C_BLUE/C_RED raw and dereddened | written, tested |
| `02_passS_measure.py` | DESI coadds, STAR rows (SPECTYPE=STAR, ZCAT_PRIMARY): pass-B chain + star lookup columns (SUBTYPE, Gaia, parallax, PM) + starlib → `~/data/sdssdesi_passS/{dark,bright}/hpx<N>.parquet` (197 columns); same CLI as pass B/Q; to be called from the pass-Q DESI walk so each coadd is read once | written, tested on 1 coadd |
| `03_passS_check.py` | Gaia G/BP/RP check, Step-1 fiber-flux check, RV, indices per SUBTYPE | written |
| `08_empca_stars.py` | weighted EM-PCA per temperature class on stacks + S/N>20 individuals; validation vs Redrock STAR v0.1 | to write |

Run (catalog level, safe while pass B runs):
```
cd ~/github/projects_desi/programs_stars
python3 00_build_star_lookup.py      # rebuilds the cache if missing, else prints the census
python3 01_hr_cells.py
```
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
