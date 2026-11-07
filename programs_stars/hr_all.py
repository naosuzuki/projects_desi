#!/usr/bin/env python3
"""Stellar version 8 (Nao 2026-11-06): the 2D HR groups of hr_build.py from every catalog star that has a DESI spectrum on disk, Gaia single stars only, with the He I white dwarfs as their own group.
Versions 6 and 7 used only the v2 store, which holds the members of the v2 stacks (cells of at least 100 stars, at most 300 stars per stack): 18 % of the bright and dark catalog stars. No exdisk2 read here either: the shared v2 walk kept the native spectra of every
DESI main star with SPECTYPE STAR in /Volumes/exdisk3/sdssdesi_native/desi (programs_quasars/native_hook.py), and the step 2-3 pass tables have the calibration of every one.
Steps (python3 hr_all.py <step>):
  select  every catalog star (hr_groups.load) in the step 2-3 pass tables with S/N > 5 and a distance; one spectrum per star (the higher S/N one when a star has a bright and a dark coadd); columns as the step-6 selection (step 2b A_FLUX, TILT; 3D E(B-V) or 0.86 SFD).
  store   rest-frame store from the native spectra, one resampling at the pipeline RV (star_hook.StarHook, version 2, as the v2 store): OUT/store/<coadd>.npz. Resumable, NPROC readers.
  lines   equivalent widths of He I 4471, 5876, 6678, H beta and H alpha of every star on the WD locus (S/N from the pass tables): OUT/wd_line_ew.parquet.
  stack   groups of hr_groups.assign; the Gaia single-star flags of version 7 (hr_build.single_ids); WD_HE = stars on the WD locus with S/N > 10, EW(5876) > 1.5 A at > 4 sigma, EW(4471) > 1 A and |EW(H beta)| < 3 A (the He-I-like test of hr_compare.py);
          BINARY and OTHER left out. A stack is a 0.5 mag x 0.1 colour box of one group with at least NMIN stars, no cap; T0 as version 7. Two passes over the store (the spectra do not fit in memory).
  pca     hr_build.pca: sqrt(N) stack weights, ten vectors, and up to NIND stars per stack with S/N > SN read from the store.
Outputs in OUT: selection_all.parquet, store/, wd_line_ew.parquet, selection_desi_hr2d.parquet, stacks_desi_hr2d.fits, members_desi_hr2d.csv, tell_coeff_desi_hr2d.parquet, individuals_desi_hr2d.npz, empca_desi_<group>.fits.
"""
import os, sys, time, glob, argparse
import multiprocessing as mp
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path[:0] = [HERE, os.path.join(os.path.dirname(HERE), 'programs_telluric')]
import hr_groups as H, hr_build as B, step6_stack as S6, star_hook, tell_hook, tellapply, starlib
OUT = '/Volumes/exdisk3/sdssdesi_stars/hr2d_all'; STORE = os.path.join(OUT, 'store'); NATIVE = '/Volumes/exdisk3/sdssdesi_native/desi'; SELALL = os.path.join(OUT, 'selection_all.parquet'); TAG = 'desi_hr2d'
EBV3D = os.path.expanduser('~/data/sdssdesi_cache/desi_stars_ebv3d.parquet'); NMIN = 30; SNMIN = 5.0; DROP = ('OTHER', 'BINARY'); NIND = 5; SN = 20.0; NPROC = int(os.environ.get('HR_NPROC', 1))     # one reader: exdisk3 is a hard disk, and three readers seeking against each other read 40 MB/s
LINES = dict(He4471=(4460, 4485), He5876=(5865, 5890), He6678=(6670, 6690), Hb=(4840, 4885), Ha=(6545, 6585))


def native_files():
    return {e.name[:-4] + '.fits' for e in os.scandir(NATIVE) if e.name.endswith('.npz') and '.tmp' not in e.name}


def store_files():
    return {e.name[:-4] + '.fits' for e in os.scandir(STORE) if e.name.endswith('.npz') and '.tmp' not in e.name} if os.path.isdir(STORE) else set()


def select():
    import pyarrow.parquet as pq
    t0 = time.time(); os.makedirs(OUT, exist_ok=True); need = ['FILE', 'Z', 'GAIA_DG', 'CAL_A', 'ABS_SCALE', 'SNR_ALL', 'TARGETID', 'PROGRAM', 'EBV_SFD', 'HALPHA_ST_EW', 'HALPHA_ST_EW_ERR', 'HBETA_LK_EW']
    files = [f for f in sorted(glob.glob(os.path.join(S6.PASSDIR['desi'], '*', '*.parquet'))) if 'FILE' in pq.ParquetFile(f).schema_arrow.names]
    df = pd.concat([pd.read_parquet(f, columns=[c for c in need if c in pq.ParquetFile(f).schema_arrow.names]) for f in files], ignore_index=True); h = H.load(); n0 = len(df); df = df[df.TARGETID.isin(h.TARGETID)]
    print(f'{n0:,} pass-table rows, {len(df):,} of catalog stars ({df.TARGETID.nunique():,} stars); {time.time()-t0:.0f}s', flush=True)
    k = (df.SNR_ALL > SNMIN) & np.isfinite(df.ABS_SCALE); print(f'S/N > {SNMIN:g} and a distance: {int(k.sum()):,} rows of {df[k].TARGETID.nunique():,} stars (dropped {int((~k).sum()):,} rows)', flush=True)
    df = df[k].sort_values('SNR_ALL', ascending=False); n1 = len(df); df = df.drop_duplicates('TARGETID'); print(f'one spectrum per star: {len(df):,} (dropped {n1 - len(df):,} second coadds)', flush=True)
    e3 = pd.read_parquet(EBV3D, columns=['TARGETID', 'EBV_3D']).drop_duplicates('TARGETID'); df = df.merge(e3, on='TARGETID', how='left'); e = df.EBV_3D.to_numpy(float); df['EBV_USED'] = np.where(np.isfinite(e), e, 0.86 * df.EBV_SFD.to_numpy(float))
    df['SID'] = starlib.desi_sid(df); sp = pd.read_parquet(os.path.join(S6.CACHE, 'spcal_desi.parquet'), columns=['SID', 'A_FLUX', 'TILT', 'TILT_OK']).drop_duplicates('SID'); df = df.merge(sp, on='SID', how='left'); df['TILT_OK'] = df.TILT_OK.fillna(False).astype(bool)
    nat = native_files(); n2 = len(df); df = df[df.FILE.isin(nat)]
    print(f'3D E(B-V) for {int(np.isfinite(e).sum()):,}; step-2b solution for {int(df.A_FLUX.notna().sum()):,}, tilt for {int(df.TILT_OK.sum()):,}; coadd in the native store for {len(df):,} of {n2:,}; {df.FILE.nunique():,} coadds', flush=True)
    df.reset_index(drop=True).to_parquet(SELALL); print(f'written {SELALL}; {time.time()-t0:.0f}s', flush=True)


class Hook(star_hook.StarHook):
    def gauge(self):                     # the monitor's gauge belongs to the v2 store
        pass


def _store_worker(files):
    SH = Hook(sel=SELALL, store=STORE, sel300=SELALL, archive='', version=2); n = 0; t0 = time.time()
    for f in files:
        if not SH.needs(f): continue
        z = np.load(os.path.join(NATIVE, f.replace('.fits', '.npz'))); fl = z['FLUX']; iv = z['IVAR']; w = float(z['W0']) + float(z['DW']) * np.arange(fl.shape[1])
        SH.process(f, z['TARGETID'], w, fl, iv, (iv <= 0).astype(np.int8)); n += 1
        if n % 500 == 0: print(f'  [{os.getpid()}] {n:,}/{len(files):,} coadds, {time.time()-t0:.0f}s', flush=True)
    return n


def store():
    t0 = time.time(); os.makedirs(STORE, exist_ok=True); files = sorted(pd.read_parquet(SELALL, columns=['FILE']).FILE.unique()); todo = [f for f in files if f not in store_files()]
    print(f'{len(files):,} coadds hold selected stars, {len(todo):,} to convert, {NPROC} readers', flush=True)
    with mp.get_context('spawn').Pool(NPROC) as P:
        n = sum(P.map(_store_worker, [todo[i::NPROC] for i in range(NPROC)]))
    left = [f for f in files if f not in store_files()]; print(f'done: {n:,} converted; {len(files) - len(left):,}/{len(files):,} in the store; {time.time()-t0:.0f}s', flush=True)


def ew(w, f, iv, lo, hi):
    i = (w >= lo) & (w <= hi) & (iv > 0); c = (((w >= lo - 25) & (w < lo)) | ((w > hi) & (w <= hi + 25))) & (iv > 0)
    if i.sum() < 5 or c.sum() < 8: return np.nan, np.nan
    C = np.median(f[c])
    if not C > 0: return np.nan, np.nan
    dw = np.gradient(w)[i]; return np.sum((1 - f[i] / C) * dw), np.sqrt(np.sum(dw**2 / iv[i])) / C


def groups():
    h = H.load(); g, dw, _ = H.assign(h, H.ms_ridge(h)); h['GROUP'] = g; h['DW'] = dw; return h


def lines():
    t0 = time.time(); h = groups(); d = pd.read_parquet(SELALL, columns=['SID', 'TARGETID', 'FILE', 'SNR_ALL']).merge(h[['TARGETID', 'GROUP', 'MG0', 'BPRP0', 'DW']], on='TARGETID'); d = d[d.GROUP.str.startswith('WD') & d.FILE.isin(store_files())]; rows = []
    print(f'{len(d):,} stars on the WD locus in the store', flush=True)
    for f, g in d.groupby('FILE'):
        z = np.load(os.path.join(STORE, f.replace('.fits', '') + '.npz')); pos = {int(t): i for i, t in enumerate(z['TARGETID'])}; F_ = z['FLUX']; I_ = z['IVAR']; w = 10**(starlib.DLOG * (int(z['RID0']) + np.arange(F_.shape[1])))
        for r in g.itertuples():
            i = pos.get(int(r.TARGETID))
            if i is None: continue
            o = dict(SID=r.SID, TARGETID=r.TARGETID, GROUP=r.GROUP, MG0=r.MG0, BPRP0=r.BPRP0, DW=r.DW, SNR=r.SNR_ALL)
            for n, (a, b) in LINES.items(): o[n], o[n + '_s'] = ew(w, F_[i].astype(float), I_[i].astype(float), a, b)
            rows.append(o)
    t = pd.DataFrame(rows); t.to_parquet(os.path.join(OUT, 'wd_line_ew.parquet')); s = he_like(t); print(f'{len(t):,} measured; He-I-like {int(s.sum()):,}: ' + str(t[s].GROUP.value_counts().to_dict()) + f'; {time.time()-t0:.0f}s', flush=True)


def he_like(t):
    return (t.SNR > 10) & (t.He5876 > 1.5) & (t.He5876 / t.He5876_s > 4) & (t.He4471 > 1) & (t.Hb.abs() < 3)


class T0Two(B.DirectT0):
    """DirectT0 for two passes: run_stack writes the table before the passes when it does not hold the spectra in memory, and the passes apply T0 twice; write once afterwards, one row per star."""
    last = None

    def __init__(self, fam, tag, df):
        super().__init__(fam, tag, df); T0Two.last = self

    def write(self, outdir):
        pass

    def flush(self, outdir):
        seen = set(); self.rows = [r for r in self.rows if not (r[0] in seen or seen.add(r[0]))]; B.DirectT0.write(self, outdir)


def gaia_distance(d):
    """ABS_SCALE from the Gaia DR3 parallax (+0.017 mas) that gives MG0 in the catalog, not from the fibermap: where the fibermap has PARALLAX = 0 the pass tables hold 0.017 mas, 59 kpc (5,836 stars of the first v8 run, 69 stacks spoiled).
    Spectra need a flux scale (A_FLUX or GAIA_DG, as every star of version 7 had), and a Gaia source gives one spectrum (two TARGETIDs can share one)."""
    g = pd.concat([pd.read_csv(B.GCSV.format(i=i), usecols=['targetid', 'source_id', 'parallax']) for i in range(1, 6)], ignore_index=True).drop_duplicates('targetid').rename(columns=dict(targetid='TARGETID', source_id='SOURCE_ID', parallax='PLX'))
    d = d.merge(g, on='TARGETID', how='left'); new = (100.0 / (d.PLX + 0.017))**2; r = (d.ABS_SCALE / new).to_numpy()
    print(f'  pass-table ABS_SCALE / Gaia DR3: within 1 % for {int((np.abs(r - 1) < 0.01).sum()):,}, off for {int((np.abs(r - 1) >= 0.01).sum()):,} (59 kpc placeholder {int((d.ABS_SCALE > 3.4e7).sum()):,}), no DR3 parallax {int(new.isna().sum()):,}', flush=True)
    d['ABS_SCALE'] = new; n0 = len(d); cal = d.A_FLUX.notna() | d.GAIA_DG.notna(); d = d[cal & np.isfinite(d.ABS_SCALE) & (d.PLX > 0)]; n1 = len(d)
    d = d.sort_values('SNR_ALL', ascending=False).drop_duplicates('SOURCE_ID'); print(f'  no flux scale or parallax: {n0 - n1:,} dropped; second target of a Gaia source: {n1 - len(d):,} dropped', flush=True); return d


def stack():
    t0 = time.time(); h = groups(); n_all = h.GROUP.value_counts(); h = h[h.TARGETID.isin(B.single_ids(h))]
    t = pd.read_parquet(os.path.join(OUT, 'wd_line_ew.parquet')); he = set(t.TARGETID[he_like(t)]); k = h.TARGETID.isin(he) & h.GROUP.str.startswith('WD'); h.loc[k, 'GROUP'] = 'WD_HE'
    print(pd.DataFrame(dict(catalog=n_all, single=h.GROUP.value_counts())).fillna(0).astype(int).to_string(), flush=True)
    feh = pd.read_parquet(H.CAT, columns=['targetid', 'mh_gspphot']).drop_duplicates('targetid').rename(columns=dict(targetid='TARGETID', mh_gspphot='FEH'))
    d = pd.read_parquet(SELALL).merge(h[['TARGETID', 'GROUP', 'MG0', 'BPRP0']], on='TARGETID').merge(feh, on='TARGETID', how='left'); n0 = len(d); d = d[~d.GROUP.isin(DROP) & d.FILE.isin(store_files())]; d = gaia_distance(d)
    d['CELL'] = np.floor(d.MG0 / 0.5).astype(int).astype(str) + '_' + np.floor(d.BPRP0 / 0.1).astype(int).astype(str) + '_' + d.GROUP.str.replace('_', ''); n = d.CELL.map(d.CELL.value_counts())
    print(f'{n0:,} single stars selected, {len(d):,} in the store and a kept group, {int((n < NMIN).sum()):,} in boxes below {NMIN}', flush=True); d = d[n >= NMIN].copy(); d['STACK'] = d.CELL
    c = d.groupby('CELL').agg(N=('SID', 'size'), MG0_MED=('MG0', 'median'), COL0_MED=('BPRP0', 'median')).reset_index(); c['HOMOGENEOUS'] = True; c['CAUSES'] = ''; print(d.GROUP.value_counts().to_dict(), f'{len(c)} stacks, largest {int(c.N.max()):,} stars', flush=True)
    d.reset_index(drop=True).to_parquet(os.path.join(OUT, f'selection_{TAG}.parquet')); S6.OUTDIR = OUT; tell_hook.StarTelluric = T0Two
    a = argparse.Namespace(survey='desi', tag=TAG, onepass=False, telluric=True, maxfiles=0, nsig=4.0, outdir=OUT)
    S6.run_stack(a, c, None, d, sorted(d.STACK.unique()), lambda x, m: S6.iter_store(x, m, STORE), t0); T0Two.last.flush(OUT); print(f'stack done, {time.time()-t0:.0f}s', flush=True)
    gauge('stack', f'v8 (hr_all.py, DESI only): {len(c)} stacks, {len(d):,} stars')


def individuals(d):
    """up to NIND stars per stack with S/N > SN, calibrated as in step 6 and divided by T0 (as hr_build.individuals does with the v2 cache); cached in OUT."""
    cp = os.path.join(OUT, f'individuals_{TAG}.npz'); s = pd.read_parquet(os.path.join(OUT, f'selection_{TAG}.parquet')); s = s[s.SNR_ALL > SN].sort_values('SNR_ALL', ascending=False).groupby('STACK').head(NIND); key = (s.SID.astype(str) + ':' + s.ABS_SCALE.map('{:.6g}'.format)).to_numpy()     # the scale is in the key: a new distance makes a new cache
    if os.path.exists(cp):
        z = np.load(cp)
        if 'KEY' in z.files and set(z['KEY'].tolist()) == set(key): print(f'  individuals from {cp}', flush=True); return z['FLUX'].astype(float), z['W'].astype(float), z['IDS'].astype(str), z['STACKS'].astype(str)
    b = tellapply.basis('desi'); zero = np.zeros(len(b.T)); F, W, I, S = [], [], [], []
    for row, f, iv, mk, rid0 in S6.iter_store(s, 0, STORE):
        x, w = S6.star_spectrum(row, f, iv, mk, rid0); C = b.factor_of(b.grid(starlib.SWAVE * 10**float(tell_hook.logz(row.Z_TRUE))), zero); F.append(x / C); W.append(w * C**2); I.append(row.SID); S.append(row.STACK)
    F, W, I, S = np.array(F), np.array(W), np.array(I, dtype=str), np.array(S, dtype=str); np.savez(cp, FLUX=F.astype(np.float32), W=W.astype(np.float32), IDS=I, STACKS=S, SEL=s.SID.to_numpy().astype(str), KEY=key); return F, W, I, S


def pca():
    """hr_build.pca; a pixel needs min(5, stacks of the group) covering stacks, so that a group of fewer than 5 stacks (WD_HE may be one) still gets a basis; groups of 5 or more are as in versions 6 and 7."""
    emp = B.S7.empca
    def empca(X, W, nvec, niter, **k):
        k.setdefault('min_cover', min(5, k.get('n_ref') or len(X))); return emp(X, W, nvec, niter, **k)
    B.S7.empca = empca; B.OUT = OUT; B.TAG = TAG; B.individuals = individuals; B.pca()
    gauge('empca', f'v8 (hr_all.py, DESI only): {len(glob.glob(os.path.join(OUT, "empca_desi_*.fits")))} group bases in {OUT}')


def gauge(name, note):
    """the shared step-6 and step-7 code leave the monitor's star gauges at 'waiting' for an SDSS half; v8 is DESI only, so mark them done."""
    import json; json.dump({'done': 1, 'total': 1, 'status': 'done', 'note': note, 'finished': time.time()}, open(os.path.expanduser(f'~/data/sdssdesi_status/step_star_{name}.json'), 'w'))


if __name__ == '__main__':
    {'select': select, 'store': store, 'lines': lines, 'stack': stack, 'pca': pca}[sys.argv[1]]()
