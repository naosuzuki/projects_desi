#!/usr/bin/env python3
"""Galaxy normalization-window flux ratios (Nao 2026-10-31), the galaxy counterpart of programs_quasars/08_window_ratio.py.

Every spectrum is unit mean F_lambda in its slice's window (normwin): W6125 = 6000-6250 (z < 0.10), W5450 = 5400-5500
(0.10-0.25), W4700 = 4600-4800 (0.25-0.60), W3050 = 3000-3100 (0.60-1.45), W2150 = 2100-2200 (>= 1.45). A stack whose
members also cover a neighbouring window (NSPEC >= COVER x N_STACKED at every pixel of it) gives the ratio directly:
its mean flux there. The four links of the chain:
  R1 = F(5450)/F(6125)   switch at z = 0.10
  R2 = F(4700)/F(5450)   z = 0.25
  R3 = F(3050)/F(4700)   z = 0.60
  R4 = F(2150)/F(3050)   z = 1.45
Each is fit n-weighted as log R = a + b (Dn4000 - 1.5) + c log10(max(EW, 1)) (+ d z), with the cell medians.
Output: <out>/gal_window_ratio.csv (one row per stack and ratio), gal_window_ratio_fit.json, gal_window_ratio.png.
Usage: python3 14_window_ratio.py [--stacks DIR ...] [--out DIR]
"""
import os, sys, glob, json, argparse
import numpy as np
import pandas as pd
import fitsio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

WIN = {'W6125': (6000, 6250), 'W5450': (5400, 5500), 'W4700': (4600, 4800), 'W3050': (3000, 3100), 'W2150': (2100, 2200)}
RATIOS = {'R1': ('W5450', 'W6125', 0.10), 'R2': ('W4700', 'W5450', 0.25), 'R3': ('W3050', 'W4700', 0.60), 'R4': ('W2150', 'W3050', 1.45)}
DEF_STACKS = [os.path.expanduser('~/data/sdssdesi_passC_resp/stacks'), os.path.expanduser('~/data/sdssdesi_passG_sdss_v2/stacks_sdss')]


def window_of(zb):
    import normwin
    w = normwin.slice_window(zb * 0.05, (zb + 1) * 0.05)
    return next(k for k, v in WIN.items() if tuple(v) == tuple(w))


def measure(stacks, cover):
    rows = []; wide = []
    for d in stacks:
        src = 'SDSS' if 'passG' in d else 'DESI'
        for fn in sorted(glob.glob(os.path.join(d, 'stacks_z*.fits'))):
            F = fitsio.FITS(fn); f = F['FLUX'].read(); ns = F['NSPEC'].read(); w = F['WAVE'].read(); c = F['CELLS'].read(); F.close()
            for i in range(len(c)):
                zb = int(c['ZB'][i]); own = window_of(zb); N = int(c['N_STACKED'][i])
                cov = {}
                for k, (lo, hi) in WIN.items():
                    s = (w >= lo) & (w <= hi)
                    cov[k] = (ns[i, s].min() / max(N, 1), float(np.mean(f[i, s])))
                wide.append(dict(SRC=src, ZB=zb, MB=int(c['MB'][i]), CB=int(c['CB'][i]), EB=int(c['EB'][i]), DB=int(c['DB'][i]), N=N,
                                 Z=c['Z_MED'][i], DN=c['DN4000_MED'][i], EW=c['EW_MED'][i], own=own,
                                 **{k: v[1] for k, v in cov.items()}, **{'C' + k[1:]: v[0] for k, v in cov.items()}))
                for r, (num, den, zs) in RATIOS.items():
                    if own not in (num, den) or cov[num][0] < cover or cov[den][0] < cover:
                        continue
                    R = cov[num][1] / cov[den][1]          # one of the two is 1 by construction
                    if not R > 0:
                        continue
                    rows.append(dict(ratio=r, logR=np.log10(R), own=own, SRC=src, ZB=zb, MB=int(c['MB'][i]), CB=int(c['CB'][i]),
                                     EB=int(c['EB'][i]), DB=int(c['DB'][i]), N=N, Z=c['Z_MED'][i], MABS=c['MABS_MED'][i],
                                     COLOR=c['COLOR_MED'][i], EW=c['EW_MED'][i], DN=c['DN4000_MED'][i],
                                     COVER=min(cov[num][0], cov[den][0])))
    return pd.DataFrame(rows), pd.DataFrame(wide)


def wfit(X, y, wt):
    A = np.column_stack([np.ones(len(y))] + X); W = np.sqrt(wt)
    p, *_ = np.linalg.lstsq(A * W[:, None], y * W, rcond=None)
    res = y - A @ p
    return p, float(np.sqrt(np.sum(wt * res**2) / np.sum(wt)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stacks', nargs='+', default=DEF_STACKS); ap.add_argument('--cover', type=float, default=0.8)
    ap.add_argument('--out', default=os.path.expanduser('~/data/sdssdesi_passC_resp/qa'))
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    D, Wd = measure(a.stacks, a.cover); D.to_csv(os.path.join(a.out, 'gal_window_ratio.csv'), index=False)
    Wd.to_csv(os.path.join(a.out, 'gal_window_means.csv'), index=False)
    out = {}
    fig, axs = plt.subplots(4, 3, figsize=(14, 13))
    for k, (r, (num, den, zs)) in enumerate(RATIOS.items()):
        d = D[(D.ratio == r) & np.isfinite(D.DN) & np.isfinite(D.EW) & (D.DN > 0.8) & (D.DN < 2.6)]      # a few cells have unphysical Dn4000 medians
        ax = axs[k]
        if len(d) < 5:
            out[r] = dict(nstack=len(d)); [x.set_visible(False) for x in ax]; continue
        wt = d.N.to_numpy(float); y = d.logR.to_numpy()
        dn = d.DN.to_numpy() - 1.5; le = np.log10(np.clip(d.EW.to_numpy(), 1, None)); z = d.Z.to_numpy()
        p1, s1 = wfit([dn, le], y, wt); p2, s2 = wfit([dn, le, z], y, wt)
        m = float(np.sum(wt * y) / wt.sum())
        out[r] = dict(definition=f'F({WIN[num]})/F({WIN[den]})', z_switch=zs, nstack=len(d), ngal=int(wt.sum()),
                      nstack_DESI=int((d.SRC == 'DESI').sum()), nstack_SDSS=int((d.SRC == 'SDSS').sum()),
                      zrange=[float(z.min()), float(z.max())], logR_mean=m, R_mean=10**m,
                      scatter_dex=float(np.sqrt(np.sum(wt * (y - m)**2) / wt.sum())),
                      fit_dn_ew=dict(a=p1[0], b_dn=p1[1], c_logew=p1[2], scatter_dex=s1),
                      fit_dn_ew_z=dict(a=p2[0], b_dn=p2[1], c_logew=p2[2], d_z=p2[3], scatter_dex=s2),
                      by_DB=[dict(DB=int(g), n=int(x.N.sum()), logR=round(float(np.average(x.logR, weights=x.N)), 4))
                             for g, x in d.groupby('DB')],
                      by_MB=[dict(MB=int(g), n=int(x.N.sum()), logR=round(float(np.average(x.logR, weights=x.N)), 4))
                             for g, x in d.groupby('MB')])
        sz = 8 + 40 * wt / wt.max()
        for x_, xl, ax_ in ((d.DN, 'Dn4000 (cell median)', ax[0]), (d.MABS, 'M (cell median, band-switched)', ax[1]), (d.Z, 'z', ax[2])):
            sc = ax_.scatter(x_, y, c=le, s=sz, cmap='coolwarm_r', vmin=0, vmax=2, alpha=0.75)
            ax_.axhline(m, color='0.5', ls=':'); ax_.set_xlabel(xl); ax_.set_ylabel(f'log {r}')
        xx = np.linspace(d.DN.min(), d.DN.max(), 20)
        for e, ls in ((0, '-'), (1, '--'), (2, ':')):
            ax[0].plot(xx, p1[0] + p1[1] * (xx - 1.5) + p1[2] * e, 'k', ls=ls, lw=1)
        for (db, eb), g in d.groupby(['DB', 'EB']):                       # binned N-weighted medians (what 15_norm_factors uses)
            if g.N.sum() < 2000:
                continue
            o = np.argsort(g.logR.to_numpy()); cw = np.cumsum(g.N.to_numpy(float)[o]); mlr = g.logR.to_numpy()[o][np.searchsorted(cw, 0.5 * cw[-1])]
            ax[0].plot(np.average(g.DN, weights=g.N), mlr, 'D', ms=7, mfc=plt.get_cmap('coolwarm_r')(min(eb, 7) / 7 * 1.0), mec='k')
        ax[0].set_title(f'{r} = F({WIN[num][0]}-{WIN[num][1]})/F({WIN[den][0]}-{WIN[den][1]}), switch z={zs}  ({len(d)} stacks)', fontsize=9, loc='left')
    fig.colorbar(sc, ax=axs, label='log10 EW (Halpha z<0.46, [OII] above)', fraction=0.02)
    fig.savefig(os.path.join(a.out, 'gal_window_ratio.png'), dpi=100); plt.close(fig)
    json.dump(out, open(os.path.join(a.out, 'gal_window_ratio_fit.json'), 'w'), indent=1, default=float)
    for r, o in out.items():
        if 'R_mean' not in o:
            print(r, 'too few stacks', o); continue
        f1 = o['fit_dn_ew']; f2 = o['fit_dn_ew_z']
        print(f"{r} {o['definition']}: {o['nstack']} stacks ({o['nstack_DESI']} DESI, {o['nstack_SDSS']} SDSS), {o['ngal']:,} gal, "
              f"z {o['zrange'][0]:.2f}-{o['zrange'][1]:.2f}; mean R {o['R_mean']:.3f} (log {o['logR_mean']:+.3f}, scatter {o['scatter_dex']:.3f} dex)")
        print(f"     fit log R = {f1['a']:+.3f} {f1['b_dn']:+.3f}(Dn-1.5) {f1['c_logew']:+.3f} logEW   scatter {f1['scatter_dex']:.3f};"
              f"  + z: d_z {f2['d_z']:+.3f}, scatter {f2['scatter_dex']:.3f}")
        print('     by DB:', ' '.join(f"{b['DB']}:{10**b['logR']:.3f}(n{b['n']})" for b in o['by_DB']))
    print('outputs in', a.out)


if __name__ == '__main__':
    main()
