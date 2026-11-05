#!/usr/bin/env python3
"""Release v3 chaining with the link measured on the STACKS (Nao 2026-11-01, "option 1"), replacing the per-galaxy lookup of
23_chain_v3.py. Why: a stack normalized at 3000-3100 A is an ivar-weighted mean, and galaxies bright at 3050 relative to their
noise weigh more at every pixel, so the stack is bluer than its median galaxy: its own F(3050)/F(5450) is +0.07 dex (16-84 %
+0.03..+0.14) above the per-galaxy lookup (DESI v3 check). The link the stacks need is therefore measured on the stacks.
Per version (G or 5450), over the DESI and SDSS stacks together:
 1. every stack of the 3000-3100 slices (z >= 0.30) that covers both windows (>= COVER of its galaxies in each, F(5450) > 0)
    gives a direct ratio R_DIRECT with its error from the stack IVAR;
 2. per DB_EFF, log R_DIRECT in DM-mag bins of the version's axis (N-weighted medians, >= NBIN stacks) is fitted with the
    models of 18_gal_ratio_fit (line / softplus / sigmoid by BIC), held flat outside the fitted range: the stack-level lookup;
 3. each upper stack's link is the inverse-variance combination of its own log R_DIRECT (error from IVAR) and the lookup at
    its magnitude (error = the robust scatter of the stacks about the fit, per DB_EFF); stacks without a direct ratio (z > 0.8,
    thin coverage) take the lookup alone. LINK_SRC: 0 lookup only, 1 mostly direct (direct weight > 0.5), 2 mostly lookup.
Then as 23: renormalize in the own window, multiply by the link, write <stacks>_chain_stacks/ with cell_chain.csv, and the
stack-level lookup per version to ~/data/sdssdesi_products/gal_stack_link_v3<ver>.csv (DB_EFF, M, logR, SOURCE) and the per-survey
offset about the joint fit (applied to that survey's lookup) to gal_stack_link_v3<ver>_survey_offset.csv.
Usage: python3 27_chain_v3_stacks.py --version G|5450 --desi DIR --desi-tables GLOB --sdss DIR --sdss-tables GLOB
"""
import os, sys, glob, argparse, importlib
import numpy as np
import pandas as pd
import fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import v3lib
rf = importlib.import_module('18_gal_ratio_fit')

W5450, W3050 = (5400, 5500), (3000, 3100)
SWITCH_Z, DZ, COVER = 0.30, 0.05, 0.8
KEY = ['ZB', 'MB', 'CB', 'EB', 'DB']
DM, NBIN = 0.5, 5
MCOL = {'G': 'MABS_G_MED', '5450': 'M5450_MED'}
P = os.path.expanduser('~/data/sdssdesi_products')


def medians(tables, version):
    cols = ['ZB', 'CB', 'EB', 'DB', f'MB_{version}', f'ACC_{version}', 'NORM_OK', 'MABS_G_FULL', 'M_AB5450']
    import pyarrow.parquet as pq
    fs = [f for f in sorted(glob.glob(tables)) if set(cols) <= set(pq.read_schema(f).names)]
    t = pd.concat([pd.read_parquet(f, columns=cols) for f in fs], ignore_index=True)
    t = t[(t[f'ACC_{version}'] >= 0) & t.NORM_OK].rename(columns={f'MB_{version}': 'MB'})
    return t.groupby(KEY).agg(MABS_G_MED=('MABS_G_FULL', 'median'), M5450_MED=('M_AB5450', 'median')).reset_index()


def measure(fn, med):
    """per stack of one file: window means, coverages, errors, DB_EFF; returns (CELLS frame, own-window mask flag)."""
    F = fitsio.FITS(fn); f = F['FLUX'].read(); iv = F['IVAR'].read(); ns = F['NSPEC'].read(); w = F['WAVE'].read()
    c = F['CELLS'].read(); c = pd.DataFrame(c.astype(c.dtype.newbyteorder('='))); F.close()
    c = c.drop(columns=[x for x in ('MABS_G_MED', 'M5450_MED') if x in c]).merge(med, on=KEY, how='left')
    zb = int(c.ZB.iloc[0]); upper = (zb + 1) * DZ > SWITCH_Z + 1e-9
    s5 = (w >= W5450[0]) & (w <= W5450[1]); s3 = (w >= W3050[0]) & (w <= W3050[1]); own = s3 if upper else s5
    N = np.maximum(c.N_STACKED.to_numpy(), 1)
    with np.errstate(all='ignore'):
        def wmean(sel):
            g = ns[:, sel] > 0
            m = np.nanmean(np.where(g, f[:, sel], np.nan), axis=1)
            var = np.where(g & (iv[:, sel] > 0), 1 / np.where(iv[:, sel] > 0, iv[:, sel], 1), np.nan)
            n = np.isfinite(var).sum(1); e = np.sqrt(np.nansum(var, axis=1)) / np.maximum(n, 1)
            return m, np.where(n > 0, e, np.nan)
        m_own, _ = wmean(own); m5, e5 = wmean(s5); m3, e3 = wmean(s3)
    cov5 = ns[:, s5].min(axis=1) / N; cov3 = ns[:, s3].min(axis=1) / N
    c['UPPER'] = upper; c['M_OWN'] = m_own
    c['DB_EFF'] = v3lib.db_eff(c.Z_MED, c.DB, c.EW_MED, c.MABS_G_MED)
    okd = (cov5 >= COVER) & (cov3 >= COVER) & (m5 > 0) & (m3 > 0)
    c['R_DIRECT'] = np.where(okd, m3 / m5, np.nan)
    c['LOGR_DIRECT_ERR'] = np.where(okd, np.sqrt((e5 / m5)**2 + (e3 / m3)**2) / np.log(10), np.nan)
    return c, (f, iv, ns, w)


def fit_link(C, mcol):
    """per DB_EFF: binned N-weighted medians of log R_DIRECT vs M, best model by BIC (18_gal_ratio_fit.fit), robust scatter."""
    out, look, scat = {}, [], {}
    d = C[np.isfinite(C.R_DIRECT) & (C.R_DIRECT > 0) & np.isfinite(C[mcol])].copy(); d['lr'] = np.log10(d.R_DIRECT)
    for db, g in d.groupby('DB_EFF'):
        e = np.arange(np.floor(g[mcol].min() / DM) * DM, g[mcol].max() + DM, DM); b = np.digitize(g[mcol], e) - 1; rows = []
        for bi in range(len(e) - 1):
            v = g[b == bi]
            if len(v) < NBIN:
                continue
            o = np.argsort(v.lr.to_numpy()); cw = np.cumsum(v.N_STACKED.to_numpy()[o]); medv = v.lr.to_numpy()[o][np.searchsorted(cw, 0.5 * cw[-1])]
            mad = 1.4826 * np.median(np.abs(v.lr - medv)); rows.append((0.5 * (e[bi] + e[bi + 1]), medv, mad, len(v), max(mad, 0.02) / np.sqrt(len(v))))
        B = pd.DataFrame(rows, columns=['M', 'med', 'mad', 'n', 'err'])
        if len(B) < 3:
            continue
        name, fn, names, desc, p, Cov, chi2 = rf.fit(B)
        grid = np.arange(-25.0, -12.0 + 1e-9, 0.25); lo, hi = B.M.min(), B.M.max()
        lr = fn(np.clip(grid, lo, hi), *p)
        resid = g.lr - fn(np.clip(g[mcol], lo, hi), *p); scat[db] = float(1.4826 * np.median(np.abs(resid - np.median(resid))))
        d.loc[g.index, 'res'] = resid; out[db] = (grid, lr)
        look.append(pd.DataFrame(dict(DB_EFF=db, M=grid, logR=lr, SOURCE=np.where(grid < lo, 'held_bright', np.where(grid > hi, 'held_faint', 'fit')),
                                      MODEL=name, SCATTER=scat[db], NSTACKS=len(g))))
        print(f'  DB_EFF {db}: {len(g)} stacks, {len(B)} bins {lo:.2f}..{hi:.2f}, {name}, chi2 {chi2:.1f}, scatter {scat[db]:.3f} dex; '
              f'log R at -23/-21/-19: ' + ' '.join(f'{np.interp(m, grid, lr):+.3f}' for m in (-23, -21, -19)), flush=True)
    # per-survey offset of the stacks about the joint fit (the fit is DESI-dominated; SDSS stacks sit ~0.06 dex below it,
    # the DESI ivar-weighted stacks above the per-galaxy ratio): the lookup for a survey's stacks is the fit plus its offset
    off = {lab: float(np.median(v)) for lab, v in d.dropna(subset=['res']).groupby('SET').res}
    print('  survey offsets about the joint fit: ' + ', '.join(f'{k} {v:+.3f} dex' for k, v in off.items()), flush=True)
    return out, scat, pd.concat(look, ignore_index=True), off


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', required=True, choices=['G', '5450'])
    ap.add_argument('--desi', required=True); ap.add_argument('--desi-tables', required=True)
    ap.add_argument('--sdss', required=True); ap.add_argument('--sdss-tables', required=True)
    a = ap.parse_args(); mcol = MCOL[a.version]
    sets = [('DESI', a.desi, a.desi_tables), ('SDSS', a.sdss, a.sdss_tables)]
    meas = {}
    for lab, sdir, tabs in sets:
        med = medians(tabs, a.version)
        for fn in sorted(glob.glob(os.path.join(sdir, 'stacks_z*.fits'))):
            c, _ = measure(fn, med); c['SET'] = lab; meas[(lab, fn)] = c
    C = pd.concat(meas.values(), ignore_index=True); up = C[C.UPPER]
    print(f'version {a.version}: {len(C)} stacks, {len(up)} in the 3000-3100 slices, {int(np.isfinite(up.R_DIRECT).sum())} with a direct ratio', flush=True)
    link, scat, look, off = fit_link(up, mcol)
    look.to_csv(os.path.join(P, f'gal_stack_link_v3{a.version}.csv'), index=False)
    pd.DataFrame(dict(SET=list(off), OFFSET_DEX=list(off.values()))).to_csv(os.path.join(P, f'gal_stack_link_v3{a.version}_survey_offset.csv'), index=False)
    rows = []
    for (lab, fn), c in meas.items():
        F = fitsio.FITS(fn); f = F['FLUX'].read(); iv = F['IVAR'].read(); sc = F['SCATTER'].read(); ns = F['NSPEC'].read(); w = F['WAVE'].read(); F.close()
        upper = bool(c.UPPER.iloc[0])
        lrl = np.array([np.interp(m, *link[int(d)]) if (np.isfinite(m) and int(d) in link) else np.nan for m, d in zip(c[mcol], c.DB_EFF)]) + off.get(lab, 0.0)
        sl = np.array([scat.get(int(d), np.nan) for d in c.DB_EFF])
        lrd = np.log10(np.where(c.R_DIRECT > 0, c.R_DIRECT, np.nan)); ed = c.LOGR_DIRECT_ERR.to_numpy()
        wd = np.where(np.isfinite(lrd) & np.isfinite(ed) & (ed > 0), 1 / ed**2, 0.0); wl = np.where(np.isfinite(lrl) & np.isfinite(sl), 1 / sl**2, 0.0)
        with np.errstate(all='ignore'):
            lr = (np.nan_to_num(lrd) * wd + np.nan_to_num(lrl) * wl) / (wd + wl)
        lr = np.where(wd + wl > 0, lr, np.nan); fd = np.where(wd + wl > 0, wd / (wd + wl), 0.0)
        src = np.where(wd == 0, 0, np.where(fd > 0.5, 1, 2))
        renorm = 1.0 / c.M_OWN.to_numpy(); chain = 10**lr if upper else np.ones(len(c)); k = renorm * chain
        c = c.drop(columns=['UPPER', 'M_OWN'])
        c['NORM_WIN'] = '3000-3100' if upper else '5400-5500'; c['RENORM'] = renorm; c['R_LOOKUP'] = 10**lrl
        c['LINK_DIRECT_W'] = fd; c['LINK_SRC'] = src; c['CHAIN'] = chain
        out = sdir_out = (a.desi if lab == 'DESI' else a.sdss).rstrip('/') + '_chain_stacks'; os.makedirs(out, exist_ok=True)
        tmp = os.path.join(out, os.path.basename(fn)) + '.tmp'
        with fitsio.FITS(tmp, 'rw', clobber=True) as fo:
            fo.write((f * k[:, None]).astype(np.float32), extname='FLUX'); fo.write((iv / k[:, None]**2).astype(np.float32), extname='IVAR')
            fo.write((sc * k[:, None]).astype(np.float32), extname='SCATTER'); fo.write(ns, extname='NSPEC'); fo.write(w, extname='WAVE')
            fo.write(c.to_records(index=False), extname='CELLS')
        os.replace(tmp, os.path.join(out, os.path.basename(fn))); rows.append(c.assign(FILE=os.path.basename(fn)))
    C = pd.concat(rows, ignore_index=True)
    for lab, sdir, _ in sets:
        C[C.SET == lab].to_csv(os.path.join(sdir.rstrip('/') + '_chain_stacks', 'cell_chain.csv'), index=False)
    u = C[C.NORM_WIN == '3000-3100']
    for lab in ('DESI', 'SDSS'):
        s = u[u.SET == lab]; both = s[(s.R_DIRECT > 0) & np.isfinite(s.R_LOOKUP)]
        d = np.log10(both.R_DIRECT / both.R_LOOKUP); d2 = np.log10(both.R_DIRECT / both.CHAIN)
        print(f'{lab}: {len(s)} upper stacks, link source: lookup only {(s.LINK_SRC == 0).sum()}, mostly direct {(s.LINK_SRC == 1).sum()}, '
              f'mostly lookup {(s.LINK_SRC == 2).sum()}; direct vs stack lookup median {np.median(d):+.3f} (16-84 % {np.percentile(d, 16):+.3f}..'
              f'{np.percentile(d, 84):+.3f}); direct vs applied link {np.median(d2):+.3f} ({np.percentile(d2, 16):+.3f}..{np.percentile(d2, 84):+.3f}) dex', flush=True)


if __name__ == '__main__':
    main()
