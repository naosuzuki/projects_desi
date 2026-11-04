#!/usr/bin/env python3
"""DESI stacks and sqrt(N) EM-PCA bases per 2D HR group (hr_groups.py; Nao 2026-11-03). Spectra: the v2 store on exdisk3 only, no coadd read; T0 as release v5 (DESI: C = 1 + T0). A stack is a 0.5 mag x 0.1 colour box inside one group, at least NMIN stars.
Outputs in OUT: stacks_desi_hr2d.fits, members_desi_hr2d.csv, tell_coeff_desi_hr2d.parquet, empca_desi_<group>.fits (format of step7_empca.py; NORM table, header STACKWT, GROUP).
STAR_HR=single (Nao 2026-11-03): only stars Gaia gives no sign of being multiple. The catalog already has non_single_star = 0, RUWE < 1.4, not VARIABLE; this adds ipd_frac_multi_peak = 0, ipd_frac_odd_win = 0, one Gaia source per TARGETID, rv_chisq_pvalue > 0.01 where it exists. Outputs in hr2d_single.
Usage: python3 hr_build.py stack | pca
"""
import os, sys, time, argparse
import numpy as np, pandas as pd, fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path[:0] = [HERE, os.path.join(os.path.dirname(HERE), 'programs_telluric')]
import hr_groups as H, step6_stack as S6, step7_empca as S7, star_hook, tell_hook, tellapply, starlib
SINGLE = os.environ.get('STAR_HR') == 'single'; OUT = '/Volumes/exdisk3/sdssdesi_stars/hr2d' + ('_single' if SINGLE else ''); GCSV = os.path.expanduser('~/github/projects_gaia/csvfiles/gaiadr3_desidr2pix_star{i}.csv'); TAG = 'desi_hr2d'; V2TAG = 'desi_0.5x0.1_gspphot'; V2 = tell_hook.V2; SEL = os.path.expanduser(f'~/data/sdssdesi_cache/stars/step6_selection_{V2TAG}.parquet'); NMIN = 30; SWAVE = starlib.SWAVE; NIND = 5; SN = 20.0; DROP = ('OTHER',)


def single_ids(h):
    """TARGETIDs of the catalog that pass the extra Gaia flags; prints what each flag removes."""
    c = ['targetid', 'source_id', 'ipd_frac_multi_peak', 'ipd_frac_odd_win', 'rv_chisq_pvalue']; g = pd.concat([pd.read_csv(GCSV.format(i=i), usecols=c) for i in range(1, 6)], ignore_index=True); n = g.groupby('targetid').source_id.nunique(); g = g.drop_duplicates('targetid').set_index('targetid').reindex(h.TARGETID.to_numpy()); ok = pd.Series(True, index=g.index); steps = (('one Gaia source', g.index.map(n).to_numpy() == 1), ('ipd_frac_multi_peak = 0', g.ipd_frac_multi_peak.to_numpy() == 0), ('ipd_frac_odd_win = 0', g.ipd_frac_odd_win.to_numpy() == 0), ('rv_chisq_pvalue > 0.01 or none', ~(g.rv_chisq_pvalue.to_numpy() <= 0.01)))
    for lab, m in steps: print(f'  {lab}: keeps {m.mean() * 100:.1f} % of the catalog', flush=True); ok &= m
    print(f'  all flags: {int(ok.sum()):,} of {len(ok):,}', flush=True); return set(ok.index[ok.to_numpy()])


class DirectT0:
    """the telluric factor of DESI, C = 1 + T0 at the star's velocity. tellapply fits a_k too but they do not enter C there, so there is nothing else to do."""
    def __init__(self, fam, tag, df):
        self.tag = tag; self.rows = []; self.b = tellapply.basis('desi'); self.zero = np.zeros(len(self.b.T)); self.checked = False

    def apply(self, row, spec, w):
        q = float(tell_hook.logz(row.get('Z_TRUE', row.Z))); wave = SWAVE * 10**q; C = self.b.factor_of(self.b.grid(wave), self.zero)
        if not self.checked:
            C1 = self.b.correct(wave, spec, w, np.where(spec > 0, spec, 1.0))[0]; assert np.allclose(C, C1), 'C differs from tellapply.correct'; self.checked = True
        self.rows.append((row.SID, 'desi', q, *self.zero, 0, 0.0, np.nan)); return spec / C, w * C**2

    def write(self, outdir):
        n = len(self.zero); t = pd.DataFrame(self.rows, columns=['SID', 'SURVEY', 'LOGZ'] + [f'A{k}' for k in range(n)] + ['NPIX', 'DCHI2', 'S']); t.to_parquet(os.path.join(outdir, f'tell_coeff_{self.tag}.parquet')); print(f'T0 applied to {len(t):,} stars', flush=True)


def stack():
    t0 = time.time(); os.makedirs(OUT, exist_ok=True); h = H.load(); p = H.ms_ridge(h); h['GROUP'] = H.assign(h, p)[0]; n_all = h.GROUP.value_counts()
    if SINGLE:
        h = h[h.TARGETID.isin(single_ids(h))]; print(pd.DataFrame(dict(catalog=n_all, single=h.GROUP.value_counts())).assign(frac=lambda t: (t.single / t.catalog).round(3)).to_string(), flush=True)
    keep = set(pd.read_csv(os.path.join(V2, f'members_{V2TAG}.csv'), usecols=['SID']).SID.astype(str))
    d = pd.read_parquet(SEL); d = d[d.SID.astype(str).isin(keep)].drop(columns=['CELL', 'STACK'], errors='ignore').merge(h, on='TARGETID'); have = {e.name[:-4] + '.fits' for e in os.scandir(star_hook.STORE_V2) if e.name.endswith('.npz') and '.tmp' not in e.name}; n0 = len(d); d = d[d.FILE.isin(have) & ~d.GROUP.isin(DROP)]
    d['CELL'] = np.floor(d.MG0 / 0.5).astype(int).astype(str) + '_' + np.floor(d.BPRP0 / 0.1).astype(int).astype(str) + '_' + d.GROUP.str.replace('_', ''); n = d.CELL.map(d.CELL.value_counts()); print(f'{n0:,} stars of v2 with Gaia colours, {len(d):,} in the store and a group, {int((n < NMIN).sum()):,} in boxes below {NMIN}', flush=True); d = d[n >= NMIN].copy(); d['STACK'] = d.CELL
    c = d.groupby('CELL').agg(N=('SID', 'size'), MG0_MED=('MG0', 'median'), COL0_MED=('BPRP0', 'median')).reset_index(); c['HOMOGENEOUS'] = True; c['CAUSES'] = ''; print(d.GROUP.value_counts().to_dict(), f'{len(c)} stacks', flush=True)
    d.to_parquet(os.path.join(OUT, f'selection_{TAG}.parquet')); S6.OUTDIR = OUT; tell_hook.StarTelluric = DirectT0; a = argparse.Namespace(survey='desi', tag=TAG, onepass=True, telluric=True, maxfiles=0, nsig=4.0, outdir=OUT)
    S6.run_stack(a, c, None, d, sorted(d.STACK.unique()), lambda x, m: S6.iter_store(x, m, star_hook.STORE_V2), t0)


def individuals(d):
    z = np.load(os.path.join(V2, f'individuals_{V2TAG}_{NIND}_{SN:g}.npz')); ids = z['IDS'].astype(str); t = pd.read_parquet(os.path.join(OUT, f'tell_coeff_{TAG}.parquet')).set_index('SID').reindex(ids); st = d.set_index(d.SID.astype(str)).STACK.reindex(ids).to_numpy(); k = np.where(t.SURVEY.notna().to_numpy() & pd.notna(st))[0]; F = z['FLUX'][k].astype(float); W = z['W'][k].astype(float); b = tellapply.basis('desi'); zero = np.zeros(len(b.T))
    for i, q in enumerate(t.LOGZ.to_numpy()[k]):
        C = b.factor_of(b.grid(SWAVE * 10**q), zero); F[i] /= C; W[i] *= C**2
    return F, W, ids[k], st[k].astype(str)


def pca():
    t0 = time.time(); S7.OUTDIR = OUT; fl, ivar, good, cells = S7.load_stacks(TAG); ns = fitsio.read(os.path.join(OUT, f'stacks_{TAG}.fits'), ext='NSPEC'); gf, mg, ab = S7.g_flux(fl, good); gid = cells.STACK.str.split('_').str[2].to_numpy(); d = pd.read_parquet(os.path.join(OUT, f'selection_{TAG}.parquet'), columns=['SID', 'STACK']); Fi, Wi, Ii, Si = individuals(d); gi, mi, ai = S7.g_flux(Fi, Wi > 0); print(f'{len(cells)} stacks, {len(Ii)} individuals, {time.time()-t0:.0f}s', flush=True); cmap = dict(zip(cells.STACK, gid))
    for g in sorted(set(gid)):
        sel = np.where((gid == g) & np.isfinite(gf))[0]
        if len(sel) < 3: print(f'{g}: {len(sel)} stacks, skipped', flush=True); continue
        X = fl[sel] / gf[sel, None]; W = ivar[sel] * gf[sel, None]**2 * good[sel] * (ns[sel] >= 0.5 * cells.N.to_numpy()[sel, None]); med = np.array([np.median(w[w > 0]) if (w > 0).any() else 1.0 for w in W]); n = ns[sel].max(axis=1).astype(float); W = W / med[:, None] * np.sqrt(n / np.median(n))[:, None]
        rows = [dict(ID=cells.STACK.iloc[i], KIND='stack', CLASS=g, NORM_GFLUX=gf[i], M_G_SYN=mg[i], M_G_SYN_AB=ab[i], M_G0_CELL=cells.MG0_MED.iloc[i], COL0_CELL=cells.COL0_MED.iloc[i], N=int(cells.N.iloc[i])) for i in sel]; j = [k for k, s in enumerate(Si) if cmap.get(s) == g and np.isfinite(gi[k])]
        if j:
            Xi = Fi[j] / gi[j, None]; Wj = Wi[j] * gi[j, None]**2; Wj = Wj * (W.sum() / max(Wj.sum(), 1e-30)); X = np.vstack([X, Xi]); W = np.vstack([W, Wj]); rows += [dict(ID=Ii[k], KIND='star', CLASS=g, NORM_GFLUX=gi[k], M_G_SYN=mi[k], M_G_SYN_AB=ai[k], M_G0_CELL=np.nan, COL0_CELL=np.nan, N=1) for k in j]
        X = np.nan_to_num(X); W = np.where(np.isfinite(W), W, 0.0); mean = np.nan_to_num(np.sum(W * X, axis=0) / np.where(W.sum(0) > 0, W.sum(0), np.nan)); Xc = X - mean; nv = min(10, len(X) - 1); P, C, ok = S7.empca(Xc, W, nv, 25, min_wfrac=0.0, n_ref=len(sel), trim_edge=3); v = S7.explained(np.where(ok[None, :], Xc, 0.0), np.where(ok[None, :], W, 0.0), P); mean = np.where(ok, mean, 0.0); f = os.path.join(OUT, f'empca_desi_{g}.fits')
        if os.path.exists(f): os.remove(f)
        with fitsio.FITS(f, 'rw') as F:
            F.write(mean.astype(np.float32), header=dict(COEFF0=S7.SRID_MIN * S7.DLOG, COEFF1=S7.DLOG, GROUP=g, NVEC=nv, NITER=25, TRIMEDGE=3, MINMEMB=0.5, NOBS=len(X), TAG=TAG, STACKWT='sqrtn'), extname='MEAN'); F.write(P.astype(np.float32), extname='EIGVEC'); F.write(C.astype(np.float32), extname='COEFF'); F.write(v.astype(np.float32), extname='VARIANCE'); F.write(ok.astype(np.int8), extname='PIXMASK'); F.write(pd.DataFrame(rows).to_records(index=False), extname='NORM')
        print(f'{g}: {len(sel)} stacks + {len(X) - len(sel)} stars, {int(ok.sum())} pixels; R2 gains ' + ' '.join(f'{x:.4f}' for x in np.diff(np.r_[0, v])[:5]) + f'; R2_{nv} {v[-1]:.4f}; {time.time()-t0:.0f}s', flush=True)


if __name__ == '__main__':
    {'stack': stack, 'pca': pca}[sys.argv[1]]()
