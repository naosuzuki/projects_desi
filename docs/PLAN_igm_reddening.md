# Plan: IGM reddening — the dust content of the intergalactic medium (DESI quasars × LS/DES/HSC)

Date: 2026-11-05, reframed 2026-11-06. **Science goal (Nao): study IGM reddening** — the dust budget of the
low-z intergalactic medium — via three mutually calibrating routes, with the Ménard+10 reproduction as the
validation backbone:

* **Route 1 — halo term**: E(B−V)(r_p) profile around galaxies (Ménard+10 redone with DESI); integrate over the
  galaxy population → Ω_dust in halos (lower bound on IGM dust; M10 got ~5×10⁻⁶, MF12: ~half of cosmic dust is
  outside galaxies).
* **Route 2 — large-scale tail**: the same cross-correlation at Mpc–30 Mpc scales → dust tracing large-scale
  structure beyond individual halos (McCleary+25 saw signal to 30 h⁻¹Mpc; we add wavelength detail).
* **Route 3 — mean opacity (the new DESI measurement)**: ⟨dE(B−V)/dz⟩ from the evolution of mean quasar
  slope/color with z_q after controlling intrinsic trends (matched-L, matched-EW stacks; dust's λ-signature
  separates reddening from gray/calibration drifts). Total dust column vs z; difference from routes 1+2 = the
  homogeneous/diffuse component.

Imaging datasets (Nao, 2026-11-06, final): **(1) DESI Legacy Surveys** (primary — area + every DESI object has
forced photometry), **(2) DES latest release**, **(3) HSC DR4 (internal — insider access; deeper than public
PDR3)**. Euclid dropped.

Backbone reference: **Ménard, Scranton, Fukugita & Richards (2010, MNRAS 405, 1025)** — circumgalactic dust
reddening and magnification of background quasars by foreground galaxies — redone with DESI spectroscopic samples
and the three imaging surveys instead of SDSS.

## 1. The reference measurement (Ménard+10)

Cross-correlated ~85,000 SDSS quasars with ~24 M SDSS photometric galaxies (17 < i < 21, z_med ≈ 0.36, ~3,800 deg²):

* **Reddening**: mean quasar color excess around galaxies → ⟨E(B−V)⟩(θ) ∝ θ^−0.86 from ~20 kpc to several Mpc;
  ~10⁻³ mag at 100 kpc. Wavelength dependence of the excess in (g−i, g−r, …) → crude extinction-curve constraint,
  consistent with SMC-type dust.
* **Magnification** ("brightening"): simultaneous achromatic ⟨δm⟩(θ) → galaxy–mass correlation; combined with the
  reddening → dust-to-mass ratio of halos.
* Physics: dust mass in the halo ≈ 5×10⁷ M☉ per L* galaxy (comparable to the disk's), Ω_dust ≈ 5×10⁻⁶.

Follow-ups: Peek, Ménard & Corrales 2015 (passive galaxies as "standard crayons", 10 → 0.5 mmag from 30 kpc to 1 Mpc);
Ménard & Fukugita 2012 (Mg II absorbers).

**Citation sweep (phase 0, done 2026-11-06).** The galaxy–dust correlation HAS been redone post-SDSS, but only with
*background-galaxy photometry*, never with quasars:

* **McCleary, Huff, Bartlett & Hensley 2025 (arXiv:2503.04098, ApJ)** — DES Y3 redMaGiC backgrounds ×
  WISExSuperCOSMOS / low-z LRG foregrounds; A_V(r) = C (r/100 h⁻¹kpc)^α with α = −0.87 ± 0.04,
  C = (3.1–4.6)×10⁻³ — "excellent agreement" with M10. The closest prior work.
* **Genc, Wright & Hildebrandt 2025 (2508.15490) + Genc+26 (2609.34615)** — KiDS DR4/DR5 magnitude shifts + GGL:
  M_dust = 7.2(+0.5/−0.6)×10⁷ M☉ at log M_halo = 12.3; linear M*–M_dust relation.
* **Crenshaw, McQuinn & Werk 2026 (2606.24159, ApJL)** — Rubin ComCam DP1 (4.6 deg²), faint foregrounds:
  A_V = 0.13(r/20 kpc)^−1.8 inside 120 kpc (steeper inner slope); red-sample dust/stellar mass ≈ 2 %.
* DESI-adjacent: **Napolitano+25 (2412.15383)** — E(B−V) of Mg II absorbers in 50k DESI spectra (updates MF12);
  **Lan+25 (2511.03195)** — DESI ELG–QSO pair absorption stacking, CGM dust mass ~ ISM. *Ting-Wen Lan is the
  nearest person to this space inside DESI — coordinate.* **Guha+26 (2601.02480)** — 1,345 GOTOQs, median
  E(B−V) = 0.087 at D ≲ 20 kpc. **Zhou+24 (2409.05140)** — DESI stellar reddening map (SFD deviations up to
  80 mmag): essential systematics reference. **Chiang & Ménard 2019 (1808.03294)** — extragalactic imprints in
  SFD: both citation and systematic.

**What remains open (our paper):** (a) the quasar-based estimator — unrepeated since SDSS; (b) Legacy Surveys
~14,000 deg², larger than DES/KiDS/DP1 combined; (c) spectroscopic foregrounds (BGS/LRG) → physical-radius and
host-property splits; (d) the stacked-spectrum extinction curve and ALPHA_LAM slope profiles — completely
unexploited. Frame against McCleary+25 and Crenshaw+26; simulations to compare: McKinnon+16/17 (1505.04792,
1606.02714).

## 2. Survey choice (settled 2026-11-05)

| | DESI Legacy Surveys DR11 | DES (latest release) | HSC DR4 (internal) |
|---|---|---|---|
| area | 22,731 deg² (3.9 B sources; released 2026) | 5,000 deg² | ~1,100+ deg² Wide |
| bands | g,r,z (+i south) + WISE W1–W4 forced | g,r,i,z,Y | g,r,i,z,y |
| depth (r) | 23.4–23.9 | ~24.1 | ~26 |
| calibration | ~5 mmag | ~2 mmag (FGCM) | ~1 % |
| spec quasars inside | **DESI: millions** | few (OzDES) | some (SDSS/DESI overlap) |

Final dataset set (Nao, 2026-11-06): **(1) DESI Legacy Surveys** (primary), **(2) DES latest release**,
**(3) HSC DR4 — internal data, insider access** (deeper than public PDR3; check collaboration data-rights policy
before publishing HSC-based numbers). Euclid dropped. LS carries the headline measurement; DES and HSC provide
calibration cross-checks and the depth axis.

**Choice: Legacy Surveys + DESI spectroscopy.** The signal is mmag-level and statistics-limited → area wins; the DESI
spectroscopic quasars remove the photometric-selection systematics that limited Ménard+10 (their quasars were
photometric, their galaxies photo-z); WISE forced photometry gives a reddening-free long-wavelength anchor; and both
foreground (BGS/LRG) and background (QSO) samples have spectroscopic redshifts → exact foreground/background split.
DES stays as a calibration cross-check in the overlap; HSC only if a deep-magnification extension is wanted later.

**The novel step beyond a reproduction**: we have DESI *spectra* of every background quasar. Stacking quasar spectra
in bins of impact parameter to foreground galaxies and ratioing to the unassociated stack measures the **extinction
curve directly and continuously in λ** — no 3–5 band photometry can match that. Pass Q already measures the per-quasar
continuum slope `ALPHA_LAM`, an immediate reddening proxy.

## 3. What we have, and where

Machines (updated 2026-11-06): the pilot and catalog-level work run on **astrofs4** (Nao; internal FSU host,
reached by ssh to the gateway `suzuki@astro4.physics.fsu.edu`, then ssh to `astrofs4`; use a ProxyJump +
ControlMaster connection from the laptop). Anything touching DESI DR2 spectra runs on **Naos-Mac-mini** (DR2 loa
coadds on exdisk2, pass Q products, SDSS spPlates). The FSU laptop is for drafting and small catalog tests.

| ingredient | source | where |
|---|---|---|
| background quasars | DESI DR2: 2.72 M spec quasars; pass Q gives `ALPHA_LAM`, line EWs, `BAL_CLASS`, `Z_SYS`, calibrated spectra | Mac mini (`~/data/sdssdesi_passQ`, master tables in `~/data/sdssdesi_products/quasar/`) |
| | SDSS master table 1.14 M (`qso_master_sdss.parquet`) — overlap cross-check | Mac mini |
| | public fallback for the laptop: DESI DR1 QSO catalog (iron, 816 MB) — ~1.5 M quasars | download |
| foreground galaxies | DESI BGS (z_med≈0.2) + LRG (0.4–1.0) — spec-z; DR1 **LSS catalogs with randoms** (clustering-ready, fiber-assignment encoded in the randoms) | download (few GB) / zall on Mac mini |
| | deeper photometric sample: LS DR9/DR10 sweeps cut to r < 21.5 (or NOIRLab Data Lab SQL pull of only the needed columns — much smaller than full sweeps) | phase 3 download |
| photometry of all DESI objects | LS DR9 `FLUX_G/R/Z/W1/W2` + `MW_TRANSMISSION_*` are in the targeting columns (zall / LSS / QSO catalogs) | already attached to every object |
| dust maps | SFD (dustmaps works), CSFD / Planck GNILC as alternatives | Mac mini / small downloads |

## 4. Estimators

1. **Photometric reddening** ⟨δc⟩(r_p) for c = g−r, r−z, g−z, z−W1: mean background-source color in annuli of
   projected physical separation r_p around foreground galaxies, minus the same around the randoms. Colors
   MW-dereddened per object (`MW_TRANSMISSION_*`). Convert the color-excess vector to E(B−V) and a dust-law slope
   by fitting F99/SMC/LMC curves at the galaxy redshift (the dust is at z_gal — Ménard's "(1+z) rest-frame" factor).
   Two background-source classes (Nao, 2026-11-06):
   * **1a — quasars** (Ménard+10's estimator, unrepeated since SDSS): DESI spectroscopic QSOs, z_q > z_gal + 0.2.
   * **1b — LRGs as "standard crayons"** (Peek+15 / McCleary+25 style, but with spec-z instead of redMaGiC
     photo-z): DESI LRGs (0.4 < z < 1.0) behind BGS foregrounds; passive SEDs → small intrinsic color scatter →
     per-object E(B−V) sensitivity ~10× better than quasars; direct comparison with McCleary+25's DES result.
2. **Magnification** ⟨δm⟩(r_p): the achromatic part of the same measurement (all bands shift together). Number-count
   magnification (Scranton+05 style) needs the quasar selection function → optional, phase 3+.
3. **Spectroscopic reddening (new)**:
   a. ⟨δ ALPHA_LAM⟩(r_p) from pass Q — per-quasar continuum slope vs impact parameter; immediate, already measured.
   b. Stacked-spectrum ratio: median calibrated quasar spectrum in r_p bins ÷ reference stack (matched in z, L, EW to
      kill intrinsic trends) → extinction curve A(λ)/A(V) of circumgalactic dust, continuous over rest 1216–9800 Å of
      the quasars, at the foreground (1+z_gal) scaling. The headline plot.
   c. Later: per-quasar E(B−V) from the eigenspectra × extinction-law fit (fit model of HANDOVER_quasars 8.15 H), once
      the quasar PCA is final.
4. **Mean IGM opacity (route 3, new)**: ⟨dE(B−V)/dz⟩ from the evolution of mean quasar slope (ALPHA_LAM) and
   matched stacked spectra vs z_q over 1.5 < z_q < 3, controlling intrinsic trends by matching in L, EW;
   dust's λ-signature separates reddening from gray evolution/calibration drift. Total cosmic dust column vs z;
   subtracting routes 1+2 isolates the diffuse/homogeneous IGM component. Hardest systematics; phase 2+.

Binning: r_p = 20 kpc – 10 Mpc, 4 bins/dex (physical, at z_gal); split foreground by class (BGS vs LRG), stellar
mass / M_r, and z slices. Require z_q > z_gal + 0.2 and r_p computed at z_gal.

## 5. Systematics and null tests

* **Galactic dust residuals** (SFD errors correlate with everything on the sky): the key null test is **stars as
  background sources** — LS point sources behind the same galaxies must show zero excess reddening. Repeat with CSFD
  and Planck GNILC corrections; difference maps bound the residual.
* **Imaging systematics** (depth, seeing, stellar density modulate both galaxy selection and photometry): use the LSS
  randoms + LS maskbits; cross-correlate the color maps with the survey-property maps; jackknife over ~100 sky patches
  for errors.
* **Fiber assignment** imprints angular selection on the spectroscopic galaxy sample → use the LSS randoms (they carry
  it); cross-check the profile with the photometric galaxy sample (phase 3), whose positions are selection-free.
* **Intrinsic quasar population trends**: cut `BAL_CLASS > 0`; verify the color/slope distribution of quasars does not
  vary with galaxy density after the cuts (it shouldn't — intrinsic properties don't know about z≈0.3 foregrounds).
* **Variability**: LS colors are imaging-epoch, DESI slopes are spectrum-epoch; agreement of estimators 1 and 3a is
  itself a systematics check (variability adds noise, not bias, to either).
* **Redshift leakage / physical association**: spec-z on both sides makes this exact (unlike Ménard's photo-z galaxy
  sample) — verify with the z_q > z_gal + 0.2 cut varied.

## 6. Feasibility estimate (first concrete task)

Back-of-envelope to be done properly in `01_feasibility.py`: ~2 M usable quasars over ~14,000 deg² (~140/deg²) ×
~10 M BGS+LRG galaxies → ~10⁷–10⁸ pairs inside 1 Mpc (3.7′ at z=0.3). Per-quasar color error ~0.03–0.05 mag at r≈21
→ σ⟨δc⟩ ≈ 0.04/√N_pairs ≈ 10⁻⁵–10⁻⁴ per r_p bin vs a 10⁻³–10⁻² mag signal at 30–100 kpc: S/N per bin of tens,
comfortably better than Ménard+10, now systematics-limited. The notebook should produce the predicted S/N(r_p) curve
per estimator and per foreground split, and decide the binning.

## 7. Phases

* **Phase 0 (now)**: ADS citation sweep of Ménard+10 (confirm nobody has redone it; collect the follow-up
  measurements to compare to); feasibility notebook (§6); download DR1 LSS catalogs (BGS_ANY, LRG + randoms) and the
  DR1 QSO catalog into `~/data/igmdust/` on astrofs4 [sizes measured 2026-11-05: QSO cat 0.82 GB +
  BGS_ANY/LRG clustering 0.86 GB + randoms 1.9–3.9 GB per set ⇒ 5.5 GB minimum, 9.4 GB with 2 randoms/cap;
  awaiting Nao's approval and astrofs4 access].
* **Phase 1 — pilot (astrofs4, DR1 public data only)**: estimator 1 on DR1 quasars × DR1 BGS/LRG with targeting
  photometry; jackknife errors; star null test; SFD vs CSFD. **Gate: signal detected with the expected amplitude and
  θ-slope → proceed.**
* **Phase 2 (Mac mini)**: scale to DR2 (2.72 M quasars, pass Q table); estimator 3a (slope profile); estimator 3b
  (stacked-spectrum extinction curve) using the pass-Q calibrated rest spectra machinery; magnification ⟨δm⟩.
* **Phase 3 (optional, downloads)**: LS photometric foreground sample (Data Lab pull, r<21.5, photo-z) → selection-free
  cross-check + z to ~0.8; DES overlap cross-check; number-count magnification.
* **Phase 4 — physics & paper**: Σ_dust(r_p) profile, halo dust mass vs stellar mass (BGS splits), extinction-curve
  type (SMC/MW, R_V), Ω_dust; comparison table vs Ménard+10 / Peek+15. Paper repo: [decide — standalone letter
  (reproduction + spectroscopic extinction curve) vs section of a larger paper].

## 8. Conventions

Scripts in `programs_igmdust/` (numbered, same style as programs_quasars; shared pieces imported from qsolib/pass Q).
Products in `~/data/igmdust/` on whichever machine runs the step. Document in a new ms (phase 4 decision). This file
is the plan of record; status sections get appended below as work proceeds (HANDOVER style).

## 9. Open decisions for Nao

1. ~~Project short name~~ — "igmdust" (settled).
2. ~~Approve phase-0 downloads~~ — DONE 2026-11-06: 9.4 GB DR1 set (QSO_cat_iron_cumulative_v0 sha256-verified +
   BGS_ANY/LRG clustering dat + 2 randoms/cap, sizes byte-verified) in
   `/data/astrofs2_1/suzuki/data/DESI/catalog/dr1/`. LS sweeps on hold for Nao's go (→ `DESI/catalog/ls/`):
   **DR11 is the program photometry** (dr11/{north,south}/sweep/11.0/ MAIN sweeps only: 1,600 S files ≈3.0 TB
   + 526 N files ≈0.1 TB ≈ **3.1 TB**; the -lightcurves/-extra/-photo-z row-matched trees bring the full
   directory to ~10.5 TB and are NOT needed — photo-z maybe in phase 3; latest reductions, i-band + full
   DES/DELVE depth south, Gaia DR3) — fits the 49 TB volume; DR9 sweeps (370 GB) optional as the
   targeting-photometry cross-check. Run downloads on astro4 (gateway; mounts the volume,
   no SLURM), NEVER as an astrofs4 job — see
   machine notes in memory: astrofs4 wraps every ssh login/command in srun (system-level, no root; admin fix
   pending), so jobs there can lock everyone out.
3. ~~Survey set~~ — settled 2026-11-06: LS primary + DES latest + HSC DR4 internal; Euclid dropped.
4. DR1-public-first (astrofs4, randoms available, reproducible) then DR2 — proceeding DR1-first unless Nao says
   otherwise.
5. Foreground priority: BGS first (closest analog of Ménard's z≈0.36 sample) with LRG as the second bin?
6. Standalone paper vs section; and whether the spectroscopic extinction curve (estimator 3b) and the mean-opacity
   measurement (estimator 4) are in the first paper or split into a second.
7. HSC DR4 internal-data publication rights — confirm with the collaboration before HSC numbers go in a paper.
