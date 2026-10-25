# programs_quasars — quasar pass (pass Q) and templates

Design: `docs/HANDOVER_quasars.md` §8 (decisions by Nao, 2026-10-24); paper: ms.tex §3.7.
Imports the common pieces from `../programs_templates/02_passB_measure.py` (filters, rebin, extinction, galaxy lines).

| script | what |
|---|---|
| `qsolib.py` | shared measurements on the 600–9900 Å rest grid: power-law slope, broad + oxygen lines (flux, EW, σ, centroid/peak velocity), narrow [O III]/[O I], forest, AI/BI troughs for C IV/Si IV/Mg II/Al III → BAL_CLASS, line redshifts, Z_SYS (peak of [O III] → Mg II → O I 1304 → pipeline z) |
| `01_passQ_sdss.py` | **runs first** (internal disk): DR17 v5_13_2 and DR8 spPlates, CLASS=QSO & SPECPRIMARY rows; Step 1 scale to PSFMAG g,r (DR17; DR8 no scale — specObj has no PSF mags), F99 with plugmap SFD_EBV, rest frame by integer shift (spPlate grid is already log10 1e-4), M_I(z=2) with Richards+06 K; output `~/data/sdssdesi_passQ_sdss/{dr17,dr8}/spPlate-PPPP-MMMMM.parquet` |
| `02_passQ_desi.py` | DESI DR2 loa coadds, QSO rows; pass-B chain then quasar measurements; output `~/data/sdssdesi_passQ/{dark,bright}/hpx<N>.parquet`. Run only after pass B (external disks) |

Status 2026-10-24 night: SDSS run started (`~/data/sdssdesi_status/passQ_sdss.log`, ~1.3 s/plate → ~3 h for 4131 + 3953
plates; 978 k DR17 + 163 k DR8 primary QSO rows before cuts). DESI driver tested on one coadd (identical to the earlier
`programs_templates/05_passQ_measure.py`, now removed). Test findings (24 SDSS plates): Step 1 vs PSF mags gives
CAL_A median 1.00 (vs fiber mags −0.70 mag BOSS / −0.30 SDSS-I/II: the pipeline calibrates to PSF flux); S/N>5 slopes
α_λ = −1.32 (DR17), −1.44 (DR8) vs −1.60 (DESI test coadd) — to be examined with the full samples (Gaia XP check).

Run:
```
cd ~/github/projects_desi/programs_quasars
OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 01_passQ_sdss.py --survey all --all > ~/data/sdssdesi_status/passQ_sdss.log 2>&1 &
# after pass B:
for i in 0 1; do OPENBLAS_NUM_THREADS=1 nohup nice -n 10 python3 02_passQ_desi.py --all --shard $i/2 > ~/data/sdssdesi_status/passQ_desi$i.log 2>&1 & done
```
Next: cell table + stacks (pass C) scripts `03_cells.py`, `04_stack.py`; EM-PCA `05_empca.py`.
