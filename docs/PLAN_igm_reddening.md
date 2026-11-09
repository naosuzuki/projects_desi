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

## 10. Status 2026-11-08 — why the Legacy QSO results are weaker than Ménard+10, and the fix

Nao's question: Legacy is deeper than SDSS and we have 28× Ménard's quasars — why are the first results *less*
powerful than M10? Measured answer (per-pair scatter = `sigma × sqrt(npair)` from `paircount_bgs_mq_qso.npz`):

1. **Depth is irrelevant; the noise is intrinsic.** Per-pair QSO color scatter is 0.33–0.36 (g−r), 0.35–0.45
   (r−z), 0.85–0.93 (z−W1) mag — quasar SED diversity + variability, NOT photometric error. The §6 feasibility
   assumption (0.03–0.05 mag) was ~8× optimistic; the promised "S/N of tens per bin" is really 2–4σ per bin.
   Same floor applied to M10 — deeper imaging cannot beat it; only pair counts and intrinsically-narrow
   backgrounds (crayons) can.
2. **Inside ~500 kpc (shot-noise regime): the foreground side is the deficit.** BGS ≈ 370/deg² vs M10's 24M
   photometric galaxies ≈ 6,300/deg² — 17× sparser. Net pairs: only ~8× M10 at 144 kpc (5.8×10⁵ vs ~7×10⁴),
   nowhere near the naive 34×, and the estimator systematics (anchor, SFD imprint, g-band covariance) spent
   that margin.
3. **Beyond ~1 Mpc (saturation regime): only N_QSO matters, and we win 5× — if systematics allow.** The
   effective per-pair scatter grows to 3.8 (g−r) by 7.5 Mpc because each quasar sits in the annuli of
   hundreds of foregrounds (pairs stop being independent); the error floor is σ_c/√N_QSO. With 2.37M vs 85k
   QSOs our floor is √28 ≈ 5.3× below M10's. The 2 Mpc offset anchor was subtracting the signal there; the
   fix wave (15 Mpc anchor + CSFD + cov-ML + star null) is exactly the repair. (McCleary verdict on DES shows
   the recipe works.)

**Revision (promotes the phase-3 photometric-foreground item to active):** the Ménard-replica run —
foreground = every 16 < r < 21 galaxy in the DR11 sweeps (TYPE≠PSF/DUP, NOBS_GRZ≥1, MASKBITS clean; ~17×
BGS pair counts inside 500 kpc, z_eff = 0.30 for physical binning), randoms = official DR11 imaging randoms
(2 files, cut with the same footprint rules — LSS randoms would imprint DESI fiber-assignment holes the
photometric sample doesn't have), backgrounds = merged QSOs with CSFD, 10 kpc – 45 Mpc, 4 shards, cov-ML.
Scripts: `11_sweepgal.py` (sweep scan + randoms cut), `02_paircount.py --fg photgal`,
`00d_download_dr11_randoms.sh` (astro4), `12_submit_menard_astrofs4.sh` → `13_slurm_menard.sbatch`
(astrofs4 SLURM chain scan → count0+count1 → finish; 2×10 cores/20G so logins stay possible).
Expected: ~4× smaller errors than BGS×QSO inside 500 kpc, M10-beating errors at every radius; BGS/LRG
spec-z runs remain the physical-radius/host-split measurement, the replica is the raw-power measurement.

**Forecast (`14_forecast.py`, empirical noise model from the measured pair counts, 2026-11-08):**

| route | per-bin S/N at 27–450 kpc | profile S/N vs Ménard law | tail amplitude δC/C |
|---|---|---|---|
| BGS × QSO (measured) | 2.8–3.7σ | 8.5σ | 79% |
| M10 shot-only (flattered: no systematics charged) | 1.4–2.0σ | 7.1σ | ~21% |
| **replica ×17 pairs (forecast)** | **3.6–14.4σ** | **25σ** | 78% (systematics-floored) |
| DES crayons (measured, iteration 4) | — | — | **14.5%** |
| **Legacy crayons ×2.8 area (lcr, running)** | — | — | **8.7%** (spec-z: no f_abs) |

Verdict: the replica beats M10 ~7× per bin in the halo term; the tail via QSO colors stays
systematics-floored until the fix wave + star null land (best case √28 ≈ 5× below M10's floor);
the tail *amplitude* crown comes from the crayon route — DES already matches McCleary's ~20%, and
the Legacy spec-z crayons forecast 8.7%, 2.3× better than McCleary+26, with no photo-z absorption
systematic.

### Fix-wave results (landed 2026-11-08 evening)

bqw + lqw (CSFD, anchor 15 Mpc, cov-ML, 45 Mpc) finished in ~7 h wall (6 BGS shards + LRG):

- **LRG tail detected**: cov-ML E(B−V) = +1.24±0.25e−3 @ 671 kpc (5.0σ), +2.71±0.84e−4
  @ 2.06 Mpc (3.2σ), +1.24±0.50e−4 @ 3.6 Mpc (2.5σ); positive every bin to 19 Mpc;
  **5.0σ combined over the contamination-safe 1.2–19 Mpc**, tracking M10 r^−0.86 in slope,
  amplitude ~2–3× (LRG bias). The anchor-move prediction (§10) confirmed: the hidden 2.8σ
  anchor-zone excess resolved into the profile.
- Contamination χ² decays outward: 157/1 @ 219 kpc → 40/1 @ 671 kpc → ≤4.8/1 from 1.2 Mpc.
  671 kpc bin flagged. Joint fit still prefers δμ<0 at 1.2–2.1 Mpc (−5σ!) — unphysical sign,
  so quoted tail = conservative lower bookmark until snl arbitrates.
- **BGS**: 23 kpc +4.83±2.0e−3 (cov-ML) — headline survives CSFD+anchor+ML. Tail still
  noise-dominated (as budgeted). 144 kpc dip → 125 kpc −2.6σ: reduced by CSFD, not gone;
  flagged for snl.
- ms.tex + av_profile_legacy.pdf updated (Overleaf re-upload pending).

### Star null verdict (snl landed 2026-11-09 02:55, paircount_bgsstars_qso.npz)
- **Quasar-side offset confirmed**: star-crossed QSO colors sit at a *flat* −1.2e−3
  E(B−V)-equiv from 125 kpc all the way to 34 Mpc (−1.32±0.30e−3 @ 671 kpc,
  −1.04±0.24e−3 @ 34 Mpc). Stars cannot redden quasars → flat large-scale offsets are
  photometric zero-point/selection artifacts, exactly what the ≥15 Mpc anchor removes.
  **Anchor method validated; anchored fits stay the results of record.**
- Null subtraction is a shape test only: star and BGS samples carry *different* flat
  offsets (Δ≈1.3e−3; ML --null plateau = +1.0–1.45e−3 over 219 kpc–34 Mpc), as expected
  from different sky weighting (PSF-star density tracks galactic latitude → extinction
  residuals). So --null cannot replace the anchor.
- **125 kpc BGS dip survives the null**: the star profile is featureless there
  (−1.16±0.80e−3, same as neighbors), so subtraction leaves the dip intact relative to
  the plateau (−0.47e−3 vs +1.25e−3 at 219 kpc). Not a star-proximity artifact; stays an
  open feature of BGS × QSO colors.
- **23 kpc headline strengthens** under null subtraction: +8.44±2.0e−3 (4.2σ, vs
  +4.83e−3 anchored) — the star template is *negative* inward (−9.2e−3 @ 23 kpc,
  dominated by z−W1 = −5.1e−2: WISE blending/scattered light near stars), so if anything
  the anchored headline is conservative.
- LRG δμ<0 @ 1.2–2.1 Mpc: null colors are featureless there → not a star-proximity
  *color* artifact. Magnitude channel (04_dustmag --fg bgsstars, 2026-11-09): star
  δμ = +0.000±0.011 and −0.001±0.011 at 1.2/2.1 Mpc — no systematic; generic
  quasar-side magnitude systematics cannot supply the LRG −0.021±0.004, though the
  null's ±0.011 still admits ~half that amplitude. ms.tex item (iv) added
  (§legacyresults) + DES quasar-tail and anchor-diagnostic passages updated
  (**Overleaf re-upload pending**).

### HSC LRG star null (hls landed 2026-11-09 13:18, paircount_hsclrgstarscsfd_qso.npz)
- **HSC anchor validated in the E(B−V) channel**: star-crossed QSO E(B−V) is zero
  within 1σ at every bin from 71 kpc to 34 Mpc (|E(B−V)| ≤ 1.5e−3, typical err
  2–4e−3; e.g. −1.49±4.2e−3 @ 125 kpc, +0.24±1.3e−3 @ 34 Mpc). Cleaner than the
  Legacy null — no flat offset to speak of, so the ≥15 Mpc anchor costs nothing
  in HSC.
- **Coherent quasar-side r−i zero point**: d(r−i) sits at +1.1–1.3e−2 from 125 kpc
  to ~6 Mpc, declining to +3.0e−3 at 34 Mpc (per-bin ~2σ but coherent across 10
  bins). It does not project onto the extinction vector (E(B−V) column stays
  null), so dust results are immune; worth remembering for any future per-color
  HSC analysis.
- **Star-proximity artifacts confined to <~50 kpc**: +3.97±2.2e−2 @ 13 kpc,
  +2.43±1.6e−2 @ 23 kpc, +6.3±6.8e−3 @ 41 kpc — positive (opposite sign to
  Legacy's −9.2e−3 @ 23 kpc, which was WISE-driven; HSC has no WISE bands, this
  is optical blending/scattered light). HSC crayon bins inside ~50 kpc get a
  star-contamination caveat; beyond that the profiles are clean.
- Run context: inside job 84 (astrofs4, 2c/16G); randoms pass took 42,384 s
  (~11.8 h) for 6.47M fg × 268.8e9 pairs. Job 84 moved straight on to hlr
  (HSC LRG × quiescent crayons) at 13:24.
- HSC DR4 is collaboration-internal — these numbers stay in the plan, not in
  ms.tex, until publication rights are confirmed.

### Replica logistics (2026-11-09 morning)
- SLURM auto-killed queued jobs 79/80/81 (DependencyNeverSatisfied on job 78's exit
  status) although the scan output is good: photgal_dr11.fits, 135,952,763 galaxies.
- A pgw shard needs **~23 GB RSS** (measured on astro4) → can never run inside the
  20 GB astrofs4 cap (jobs 89 shards OOM-killed in <1 min). All four shards now run
  on astro4: 2+3 in parallel since 10:49; runner_pgw01.sh (nohup, ~/runner on astro4)
  waits for them, runs 0+1, then 02b merge + 04 anchor 15 + 10_ml + 10_ml --null.
- astrofs4 steady state: job 84 (HSC LRG, 2c/16G) + job 87 (HSC BGS wave, 10c/18G);
  memory, not cores, binds — 20.5 GB headroom fits no pgw shard.

## 11. Backgrounds v2 (Nao, 2026-11-08): multi-survey quasar photometry + DR2 spec-z

Nao's direction: integrate SDSS + Pan-STARRS + Legacy photometry, DESI DR2 for quasar spec-z. Division of labor
per the §10 power budget:

* **Quasar side (do it)** — beyond ~1 Mpc only N_QSO matters, so (a) rebuild the merged background on **DESI DR2**
  spec quasars (2.72M, vs the DR1 piece in the current 2.84M merge; 05_qso_merge → 06_sweepmatch rerun; DR2 QSO
  catalog source: Mac mini pass-Q master table or public DR2 release — confirm which), plus **SDSS DR8 + DR17
  spec-z quasars** (DR8 ⊂ DR17 specObjAll, cumulative: one DR17 pull of CLASS=QSO, ZWARNING=0 covers both; dedupe
  1" against DR16Q/DESI/Quaia, z priority DESI > SDSS > Quaia); (b) attach **SDSS DR16 and
  Pan-STARRS DR2 mean photometry** per quasar (PS1 via MAST CasJobs upload-match on astro4): multi-epoch averaging
  cuts the variability part of the 0.33 mag per-pair scatter, and u→y→W1 colors extend the extinction-curve lever
  arm. Keep per-survey color systems separate in the fit (never mix zero-points into one "color"); each survey
  contributes its own Δcolor vector with its own R coefficients, combined in the cov-ML solve.
* **Galaxy side (don't)** — SDSS/PS1 galaxy catalogs are shallower than Legacy everywhere Legacy exists, with worse
  star-galaxy separation at r≈21; the §10 photgal run already supplies the dense foreground. No integration.
* Sequence: after the fix wave + replica land (they validate the estimator), not before — one variable at a time.
