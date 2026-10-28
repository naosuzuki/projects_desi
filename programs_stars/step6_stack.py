#!/usr/bin/env python3
"""Step 6 (Nao, 2026-10-24 night): stack the stars of each HR cell whose stars are alike -> the empirical stellar library.

Inputs: the cell table and member list of step5_hr_cells_check.py (cells_<tag>.csv, members_<tag>.parquet), the pass
tables of steps 2-3 (per-star calibration, RV shift, distance, extinction), and the spectra (SDSS spPlates; DESI coadds).
Per star, before it enters its cell:
  step 2  flux x A_FLUX and, where TILT_OK, x 10^(-0.4 TILT log10(lambda/5500)) from the joint fit of step2b_spcal.py
          (fallback: 10^(-0.4 GAIA_DG), i.e. the Gaia-G scale alone; CAL_A when no Gaia G)
  step 3  integer-pixel shift by the pipeline RV onto the common log grid (starlib.SWAVE, 3550-9880 A)
  step 3b dereddening with EBV_USED (F99, R_V = 3.1; the 3D value where step3b_ebv3d.py has run)
  step 4  x ABS_SCALE = (d / 10 pc)^2  -> F_lambda at 10 pc
Mixed cells (HOMOGENEOUS = False) are split along SPLIT_AXIS before stacking: FEH (two halves at the cell median [M/H]),
HALPHA (emitters / non-emitters), BALMER (Hbeta EW above / below the cell median, e.g. DA / non-DA white dwarfs);
CONT / COLOR / INDEX are flagged but not split (the member list keeps the flag).
Stacking is streaming and two-pass (the spectra are never all in memory): pass A accumulates per cell and pixel the
weighted sum, weight, weighted square and count with per-star weights w = ivar / median(ivar) (every star counts once;
pixel quality within a star); pass B repeats with pixels more than --nsig (4) sigma from the pass-A mean rejected.
Output (OUTDIR): stacks_<tag>.fits with HDUs FLUX [ncell, npix] (erg/s/cm2/A at 10 pc), ERR (scatter / sqrt N), SCATTER,
NSPEC, and a CELLS table (cell id, HR coordinates, [M/H], N, split label, homogeneity flags); members_<tag>.csv.
Header: COEFF0 = log10 of the first pixel, COEFF1 = 1e-4 (sdsslib 2D convention).
DESI (2026-10-27): ~10 M stars, spectra on exdisk2 (~1.4 s of disk per coadd). --maxper N keeps per stack at most N stars,
the highest-S/N ones, chosen greedily to come from as few coadds as possible (candidates = the top 3N by S/N of each stack;
files with the most candidates of unfilled stacks first); --onepass reads each coadd once and keeps the selected spectra
in memory (float32) for passes A and B.
Usage:
  python3 step6_stack.py --survey sdss --tag sdss_0.5x0.1_gspphot [--nmin 100] [--snmin 5] [--maxfiles N]
  python3 step6_stack.py --survey desi --tag desi_0.5x0.1_gspphot --maxper 300 --onepass
"""
import os, sys, glob, argparse, time
import numpy as np, pandas as pd, fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import starlib
passB = starlib.passB
CACHE = os.path.expanduser('~/data/sdssdesi_cache/stars')
OUTDIR = os.path.expanduser('~/data/sdssdesi_stars')
PASSDIR = {'sdss': os.path.expanduser('~/data/sdssdesi_passS_sdss'), 'desi': os.path.expanduser('~/data/sdssdesi_passS')}
DLOG = starlib.DLOG; SWAVE = starlib.SWAVE; NPIX = len(SWAVE); SRID_MIN = starlib.SRID_MIN
A1_F99 = passB.extinction.fitzpatrick99(SWAVE, 1.0, 3.1)        # A(lambda) for A_V = 1 on the common grid
X5500 = np.log10(SWAVE / 5500.0)                                  # tilt coordinate of step 2b


def members(tag, nmin):
    cells = pd.read_csv(os.path.join(CACHE, f'cells_{tag}.csv'))
    mem = pd.read_parquet(os.path.join(CACHE, f'members_{tag}.parquet'))
    cells = cells[cells.N >= nmin]
    mem = mem[mem.CELL.isin(cells.CELL)].copy()
    return cells, mem


def split_cells(cells, mem, df):
    """apply SPLIT_AXIS to the mixed cells: returns mem with a STACK column (cell id + suffix)."""
    mem = mem.merge(df[['SID', 'HALPHA_ST_EW', 'HALPHA_ST_EW_ERR', 'HBETA_LK_EW']], on='SID', how='left')
    mem['STACK'] = mem.CELL
    for _, c in cells[~cells.HOMOGENEOUS].iterrows():
        m = mem.CELL == c.CELL
        if c.SPLIT_AXIS == 'FEH':
            med = mem.loc[m, 'FEH'].median()
            mem.loc[m, 'STACK'] = c.CELL + np.where(mem.loc[m, 'FEH'] < med, '_fehlo', '_fehhi')
        elif c.SPLIT_AXIS == 'HALPHA':
            em = (mem.loc[m, 'HALPHA_ST_EW'] < -1) & (mem.loc[m, 'HALPHA_ST_EW'] / mem.loc[m, 'HALPHA_ST_EW_ERR'] < -3)
            mem.loc[m, 'STACK'] = c.CELL + np.where(em, '_haem', '_haabs')
        elif c.SPLIT_AXIS == 'BALMER':
            med = mem.loc[m, 'HBETA_LK_EW'].median()
            mem.loc[m, 'STACK'] = c.CELL + np.where(mem.loc[m, 'HBETA_LK_EW'] >= med, '_balmhi', '_balmlo')
    return mem


class Accum:
    def __init__(self, stacks):
        self.idx = {s: i for i, s in enumerate(stacks)}; n = len(stacks)
        self.s0 = np.zeros((n, NPIX)); self.s1 = np.zeros((n, NPIX)); self.s2 = np.zeros((n, NPIX)); self.n = np.zeros((n, NPIX), np.int32)
        self.nstar = np.zeros(n, int)

    def add(self, stack, f, w):
        i = self.idx[stack]; g = w > 0
        self.s0[i] += np.where(g, w, 0); self.s1[i] += np.where(g, w * f, 0); self.s2[i] += np.where(g, w * f * f, 0); self.n[i] += g
        self.nstar[i] += 1

    def mean(self):
        with np.errstate(all='ignore'):
            m = np.where(self.s0 > 0, self.s1 / self.s0, np.nan)
            var = np.where(self.s0 > 0, self.s2 / self.s0 - m**2, np.nan)
        return m, np.sqrt(np.clip(var, 0, None))


def star_spectrum(row, f, iv, mask, rid0):
    """pipeline spectrum of one star (flux unit 1e-17, log grid from rid0) -> F_lambda at 10 pc on SWAVE, with weights."""
    good = (iv > 0) & (mask == 0)
    if np.isfinite(row.A_FLUX):
        scale = passB.FLUX_UNIT * row.A_FLUX                                             # step 2b joint fit
    else:
        scale = passB.FLUX_UNIT * (10**(-0.4 * row.GAIA_DG) if np.isfinite(row.GAIA_DG) else (row.CAL_A if np.isfinite(row.CAL_A) else 1.0))
    tilt = 10**(-0.4 * row.TILT * X5500) if (bool(row.TILT_OK) and np.isfinite(row.TILT)) else 1.0    # mag per dex -> flux factor on the grid
    A = 3.1 * row.EBV_USED * A1_F99                               # on the rest grid (RV shift < 0.4 %: negligible for the curve)
    z = np.clip(row.Z if np.isfinite(row.Z) else 0.0, -0.004, 0.004); shift = int(np.round(np.log10(1 + z) / DLOG))
    r0 = rid0 + np.arange(len(f)) - shift - SRID_MIN; ok = (r0 >= 0) & (r0 < NPIX) & good
    out = np.zeros(NPIX); w = np.zeros(NPIX)
    out[r0[ok]] = f[ok] * scale; w[r0[ok]] = iv[ok] / scale**2
    out *= tilt * 10**(0.4 * A) * row.ABS_SCALE; w *= 10**(-0.8 * A) / (tilt * row.ABS_SCALE)**2
    med = np.median(w[w > 0]) if (w > 0).any() else 1.0
    return out, w / med                                           # every star counts once; pixel quality inside the star


def iter_sdss(df, maxfiles):
    groups = df.groupby('FILE')
    for k, (base, g) in enumerate(groups):
        if maxfiles and k >= maxfiles: break
        _, plate, mjd = base.split('-')
        path = None
        for d in passB_qs_dirs(g.SURVEY.iloc[0]):
            cand = os.path.join(d, plate, base + '.fits')
            if os.path.exists(cand): path = cand; break
        if path is None: print('missing', base); continue
        F = fitsio.FITS(path); h = F[0].read_header(); rid0 = int(round(float(h['COEFF0']) / DLOG))
        fib = g.FIBERID.to_numpy() - 1
        flux = F[0].read()[fib]; ivar = F[1].read()[fib]; mask = F[2].read()[fib]; F.close()
        for j, (_, row) in enumerate(g.iterrows()):
            yield row, flux[j].astype(float), ivar[j].astype(float), mask[j], rid0


def passB_qs_dirs(survey):
    S = os.path.expanduser('~/sdssredux')
    return [f'{S}/v5_13_2'] if survey == 'dr17' else [f'{S}/DR8/26', f'{S}/DR8/103', f'{S}/DR8/104']


def iter_desi(df, maxfiles):
    groups = df.groupby('FILE')
    for k, (base, g) in enumerate(groups):
        if maxfiles and k >= maxfiles: break
        prog = 'dark' if 'dark' in base else 'bright'; hpx = base.split('-')[-1].replace('.fits', '')
        p = passB.local_path(f'/Volumes/exdisk2/desiredux/DR2/spectro/redux/loa/healpix/main/{prog}/{int(hpx)//100}/{hpx}/{base}')
        if p is None: print('missing', base); continue
        F = fitsio.FITS(p); tid = F['FIBERMAP'].read(columns=['TARGETID'])['TARGETID']
        rows = np.array([np.where(tid == t)[0][0] for t in g.TARGETID])
        w, flux, ivar, mask, _ = passB.combine_cameras(F); F.close()
        rid, M, _ = passB.rebin_matrix(w)
        for j, (_, row) in enumerate(g.iterrows()):
            r = rows[j]
            f = flux[r].astype(float); iv = ivar[r].astype(float); mk = mask[r]
            _, rf, ri = passB.rebin(w, f[None, :], iv[None, :], mk[None, :])
            yield row, rf[0] , ri[0], (ri[0] <= 0).astype(np.int8), int(rid[0])


def select_per_stack(df, maxper, pool=3, batch=200):
    """at most maxper stars per STACK, highest S/N, concentrated in few files (greedy over files)."""
    cand = df.sort_values('SNR_ALL', ascending=False).groupby('STACK').head(pool * maxper)
    need = cand.STACK.value_counts().clip(upper=maxper).to_dict()
    chosen = []
    while len(cand):
        r = cand[cand.STACK.map(need) > 0]
        if r.empty: break
        top = r.FILE.value_counts().index[:batch]
        take = r[r.FILE.isin(top)].sort_values('SNR_ALL', ascending=False)
        take = take[take.groupby('STACK').cumcount() < take.STACK.map(need)]
        for st, n in take.STACK.value_counts().items(): need[st] -= n
        chosen.append(take); cand = cand[~cand.FILE.isin(top)]
    return pd.concat(chosen, ignore_index=True)


def _progress(name, done, total, status='done', note=''):
    """progress file of the monitors (programs_monitor/monitor.py)."""
    try:
        import json; S = os.path.expanduser('~/data/sdssdesi_status'); os.makedirs(S, exist_ok=True)
        json.dump({'done': int(done), 'total': int(total), 'status': status, 'note': note}, open(f'{S}/step_star_{name}.json', 'w'))
    except Exception:
        pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--survey', default='sdss', choices=['sdss', 'desi']); ap.add_argument('--tag', required=True)
    ap.add_argument('--passdir'); ap.add_argument('--nmin', type=int, default=100); ap.add_argument('--snmin', type=float, default=5.0)
    ap.add_argument('--nsig', type=float, default=4.0); ap.add_argument('--maxfiles', type=int, default=0); ap.add_argument('--outdir', default=OUTDIR)
    ap.add_argument('--maxper', type=int, default=0, help='at most N stars per stack (highest S/N, few files); 0 = all')
    ap.add_argument('--onepass', action='store_true', help='read each spectrum once, keep it in memory for passes A and B')
    a = ap.parse_args(); t0 = time.time(); os.makedirs(a.outdir, exist_ok=True)
    cells, mem = members(a.tag, a.nmin)
    files = sorted(glob.glob(os.path.join(a.passdir or PASSDIR[a.survey], '*', '*.parquet')))
    need = ['FILE', 'Z', 'GAIA_DG', 'CAL_A', 'EBV_USED', 'ABS_SCALE', 'SNR_ALL', 'HALPHA_ST_EW', 'HALPHA_ST_EW_ERR', 'HBETA_LK_EW']
    need += ['PLATE', 'MJD', 'FIBERID', 'SURVEY'] if a.survey == 'sdss' else ['TARGETID', 'PROGRAM']
    need += ['EBV_SFD']
    import pyarrow.parquet as pq
    files = [f for f in files if 'FILE' in pq.ParquetFile(f).schema_arrow.names]          # skip the empty per-file tables
    df = pd.concat([pd.read_parquet(f, columns=[c for c in need if c in pq.ParquetFile(f).schema_arrow.names]) for f in files], ignore_index=True)
    # step 3b at stack time: take the 3D E(B-V) from the catalogs (step3b_ebv3d.py) when it exists, else 0.86 x SFD
    if a.survey == 'sdss':
        import pyarrow.parquet as pq
        cp = os.path.expanduser('~/data/sdssdesi_cache/sdss_stars_gaia.parquet'); have = pq.ParquetFile(cp).schema_arrow.names
        cat = pd.read_parquet(cp, columns=['PLATE', 'MJD', 'FIBERID'] + (['EBV_3D'] if 'EBV_3D' in have else []))
        if 'EBV_3D' in cat:
            df = df.drop(columns=[c for c in ['EBV_3D'] if c in df]).merge(cat, on=['PLATE', 'MJD', 'FIBERID'], how='left')
    else:
        p3 = os.path.expanduser('~/data/sdssdesi_cache/desi_stars_ebv3d.parquet')
        if os.path.exists(p3):
            df = df.drop(columns=[c for c in ['EBV_3D'] if c in df]).merge(pd.read_parquet(p3, columns=['TARGETID', 'EBV_3D']), on='TARGETID', how='left')
    if 'EBV_3D' in df:
        e3 = df.EBV_3D.to_numpy(float); df['EBV_USED'] = np.where(np.isfinite(e3), e3, 0.86 * df.EBV_SFD.to_numpy(float))
        print(f'3D E(B-V) used for {np.isfinite(e3).sum():,} of {len(df):,} stars (SFD fallback for the rest)', flush=True)
    df['SID'] = (df.PLATE.astype(str) + '-' + df.MJD.astype(str) + '-' + df.FIBERID.astype(str)) if a.survey == 'sdss' else starlib.desi_sid(df)
    df = df[df.SID.isin(mem.SID) & (df.SNR_ALL > a.snmin) & np.isfinite(df.ABS_SCALE)]
    # step 2b: joint spectrophotometric solution (scale A_FLUX, tilt) when it exists; else the Gaia-G scale alone
    sp = os.path.join(CACHE, f'spcal_{a.survey}.parquet')
    if os.path.exists(sp):
        s2 = pd.read_parquet(sp, columns=['SID', 'A_FLUX', 'TILT', 'TILT_OK'])
        df = df.merge(s2, on='SID', how='left')
        print(f'step-2b solution for {df.A_FLUX.notna().sum():,} of {len(df):,} stars; tilt applied where TILT_OK ({int(df.TILT_OK.fillna(False).sum()):,})', flush=True)
    else:
        df['A_FLUX'] = np.nan; df['TILT'] = np.nan; df['TILT_OK'] = False
    mem = split_cells(cells, mem, df)
    df = df.merge(mem[['SID', 'CELL', 'STACK']], on='SID', how='inner')
    stacks = sorted(df.STACK.unique())
    if a.maxper:
        n0, f0 = len(df), df.FILE.nunique(); df = select_per_stack(df, a.maxper)
        print(f'--maxper {a.maxper}: {len(df):,} of {n0:,} stars from {df.FILE.nunique():,} of {f0:,} files', flush=True)
    print(f'{a.survey}: {len(cells)} cells >= {a.nmin}, {len(stacks)} stacks after splits, {len(df):,} stars, {df.FILE.nunique()} files; {time.time()-t0:.0f}s', flush=True)
    it = iter_sdss if a.survey == 'sdss' else iter_desi
    nfile = df.FILE.nunique()
    def spectra():
        last = None; k = 0
        for row, f, iv, mk, rid0 in it(df, a.maxfiles):
            if row.FILE != last:
                last = row.FILE; k += 1
                if k % 200 == 0:
                    print(f'  {k:,}/{nfile:,} files, {time.time()-t0:.0f}s', flush=True); _progress('stack', 1, 2, 'running', f'DESI: {k:,}/{nfile:,} coadds read')
            spec, w = star_spectrum(row, f, iv, mk, rid0)
            yield row.STACK, spec, w
    if a.onepass:
        mem_spec = [(st, sp.astype(np.float32), w.astype(np.float32)) for st, sp, w in spectra()]
        print(f'read once: {len(mem_spec):,} spectra in memory, {time.time()-t0:.0f}s', flush=True)
        source = lambda: iter(mem_spec)
    else:
        source = spectra
    accA = Accum(stacks)
    for st, spec, w in source():
        accA.add(st, spec, w)
    meanA, sigA = accA.mean()
    print(f'pass A done: {accA.nstar.sum():,} stars, {time.time()-t0:.0f}s', flush=True)
    accB = Accum(stacks)
    for st, spec, w in source():
        i = accB.idx[st]
        with np.errstate(all='ignore'):
            sig = np.where(np.isfinite(sigA[i]) & (sigA[i] > 0), sigA[i], np.inf)
            reject = np.abs(spec - meanA[i]) > a.nsig * sig
        accB.add(st, spec, np.where(reject, 0.0, w))
    mean, scatter = accB.mean()
    with np.errstate(all='ignore'):
        err = scatter / np.sqrt(np.where(accB.n > 0, accB.n, 1))
    # cell table for the output
    rows = []
    for s in stacks:
        base = s.split('_')[0] + '_' + s.split('_')[1] + '_' + s.split('_')[2]; c = cells[cells.CELL == base].iloc[0]
        g = df[df.STACK == s]
        rows.append({'STACK': s, 'CELL': base, 'SPLIT': s[len(base):].lstrip('_') or 'none', 'N': int(accB.nstar[accB.idx[s]]), 'MG0_MED': c.MG0_MED, 'COL0_MED': c.COL0_MED,
                     'FEH_MED': mem[mem.STACK == s].FEH.median(), 'SNR_MED': g.SNR_ALL.median(), 'HOMOGENEOUS': bool(c.HOMOGENEOUS), 'CAUSES': c.CAUSES if isinstance(c.CAUSES, str) else ''})
    tab = pd.DataFrame(rows)
    out = os.path.join(a.outdir, f'stacks_{a.tag}.fits')
    if os.path.exists(out): os.remove(out)
    hdr = {'COEFF0': SRID_MIN * DLOG, 'COEFF1': DLOG, 'BUNIT': 'erg/s/cm2/A at 10 pc', 'SURVEY': a.survey, 'TAG': a.tag, 'NSIG': a.nsig}
    with fitsio.FITS(out, 'rw') as F:
        F.write(mean.astype(np.float32), header=hdr, extname='FLUX'); F.write(err.astype(np.float32), extname='ERR')
        F.write(scatter.astype(np.float32), extname='SCATTER'); F.write(accB.n, extname='NSPEC')
        F.write(tab.to_records(index=False), extname='CELLS')
    df[['SID', 'CELL', 'STACK', 'FILE', 'SNR_ALL']].to_csv(os.path.join(a.outdir, f'members_{a.tag}.csv'), index=False)
    print(f'written {out}: {len(stacks)} stacks x {NPIX} pixels; median N per stack {int(np.median(accB.nstar))}; {time.time()-t0:.0f}s')
    print(tab.head(10).to_string())
    _progress('stack', 1 if a.survey == 'sdss' else 2, 2, 'waiting' if a.survey == 'sdss' else 'done', f'{a.survey}: {len(stacks)} stacks' + ('; DESI after pass B' if a.survey == 'sdss' else ''))


if __name__ == '__main__':
    main()
