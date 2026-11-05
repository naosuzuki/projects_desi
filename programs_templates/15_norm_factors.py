#!/usr/bin/env python3
"""Galaxy normalization factors between the five windows, and their sanity checks (Nao 2026-10-31).
Input: <qa>/gal_window_means.csv from 14_window_ratio.py (per stack: mean flux in each window W*, coverage C*).

Normalization factor of a window = F(window) / F(5400-5500), the reference W5450, by spectral type (Dn4000 bin DB, and
EW bin EB). It is chained through the neighbouring links R1..R4 (14_window_ratio.py), each the N-weighted median over the
stacks of that type that cover both windows.
Sanity checks:
  S1 own window: each stack's mean in its own window must be 1.
  S2 closure: direct F(3050)/F(5450) from stacks covering both (z 0.6-0.8) vs the chain R2 x R3 from other stacks.
  S3 survey: DESI vs SDSS for the same cell (ZB, MB, CB, EB, DB) and the same ratio.
  S4 jump: the z = 0.6 jump of 13_norm_continuity (norm_jumps.csv) vs 1/R3 of the same type.
Output: <qa>/gal_norm_factors.png, gal_norm_factors.csv (factor per DB x window), printed checks.
"""
import os, json, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

WINS = ['W6125', 'W5450', 'W4700', 'W3050', 'W2150']
LAM = {'W6125': 6125, 'W5450': 5450, 'W4700': 4700, 'W3050': 3050, 'W2150': 2150}
LINKS = [('W5450', 'W6125'), ('W4700', 'W5450'), ('W3050', 'W4700'), ('W2150', 'W3050')]
COVER = 0.8
DBLAB = {0: 'Dn<1.15', 1: '1.15-1.40', 2: '1.40-1.65', 3: '1.65-1.90', 4: '>1.90'}


def wmed(x, w):
    o = np.argsort(x); c = np.cumsum(w[o]); return float(x[o][np.searchsorted(c, 0.5 * c[-1])])


def link(W, num, den, by):
    d = W[(W['C' + num[1:]] >= COVER) & (W['C' + den[1:]] >= COVER) & (W[num] > 0) & (W[den] > 0)].copy()
    d['logR'] = np.log10(d[num] / d[den])
    return d, {k: (wmed(g.logR.to_numpy(), g.N.to_numpy(float)), int(g.N.sum()), len(g)) for k, g in d.groupby(by)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--qa', default=os.path.expanduser('~/data/sdssdesi_passC_resp/qa'))
    a = ap.parse_args()
    W = pd.read_csv(os.path.join(a.qa, 'gal_window_means.csv'))

    # S1
    own = np.array([W.loc[i, W.own[i]] for i in W.index]); cov_own = np.array([W.loc[i, 'C' + W.own[i][1:]] for i in W.index])
    s = (cov_own >= COVER) & np.isfinite(own)
    dev = np.abs(own[s] - 1)
    print(f'S1 own window: {s.sum()} of {len(own)} stacks fully covered, mean in own window: median {np.median(own[s]):.4f}, '
          f'99th pct |dev| {np.percentile(dev, 99):.4f}, max |dev| {dev.max():.4f}')

    # links by DB and by (DB, EB)
    L = {}; Lde = {}; Ld = {}
    for num, den in LINKS:
        Ld[(num, den)], L[(num, den)] = link(W, num, den, 'DB')
        _, Lde[(num, den)] = link(W, num, den, ['DB', 'EB'])
    fac = {}
    for db in range(5):
        f = {'W5450': 0.0}
        f['W6125'] = -L[('W5450', 'W6125')][db][0] if db in L[('W5450', 'W6125')] else np.nan
        f['W4700'] = L[('W4700', 'W5450')][db][0] if db in L[('W4700', 'W5450')] else np.nan
        f['W3050'] = f['W4700'] + L[('W3050', 'W4700')][db][0] if db in L[('W3050', 'W4700')] else np.nan
        f['W2150'] = f['W3050'] + L[('W2150', 'W3050')][db][0] if db in L[('W2150', 'W3050')] else np.nan
        fac[db] = f
    F = pd.DataFrame([{'DB': db, 'Dn4000': DBLAB[db], **{w: 10**v for w, v in f.items()}} for db, f in fac.items()])[['DB', 'Dn4000'] + WINS]
    F.to_csv(os.path.join(a.qa, 'gal_norm_factors.csv'), index=False)
    print('\nnormalization factor F(window)/F(5400-5500) by Dn4000 bin (chained medians):')
    print(F.round(3).to_string(index=False))
    print('\nlinks (N-weighted median ratio; stacks):')
    for (num, den), v in L.items():
        print(f'  {num}/{den}: ' + '  '.join(f'DB{k} {10**m:.3f} ({ns})' for k, (m, n, ns) in sorted(v.items())))

    # S2 closure
    d = W[(W.C5450 >= COVER) & (W.C3050 >= COVER) & (W.W5450 > 0) & (W.W3050 > 0)].copy()
    d['direct'] = np.log10(d.W3050 / d.W5450)
    clo = []
    for (db, eb), g in d.groupby(['DB', 'EB']):
        r2 = Lde[('W4700', 'W5450')].get((db, eb)); r3 = Lde[('W3050', 'W4700')].get((db, eb))
        if r2 is None or r3 is None:
            continue
        clo.append(dict(DB=db, EB=eb, n=int(g.N.sum()), direct=wmed(g.direct.to_numpy(), g.N.to_numpy(float)), chain=r2[0] + r3[0]))
    clo = pd.DataFrame(clo)
    dd = clo.direct - clo.chain
    print(f'\nS2 closure F(3050)/F(5450): {len(clo)} (DB,EB) types, {len(d)} direct stacks; direct - chain: '
          f'median {np.median(dd):+.3f} dex, rms {np.sqrt(np.mean(dd**2)):.3f} dex ({100*(10**np.sqrt(np.mean(dd**2))-1):.1f} %)')

    # S3 DESI vs SDSS
    key = ['ZB', 'MB', 'CB', 'EB', 'DB']
    sv = []
    for num, den in LINKS:
        x = Ld[(num, den)]; p = x[x.SRC == 'DESI'].merge(x[x.SRC == 'SDSS'], on=key, suffixes=('_d', '_s'))
        for _, r in p.iterrows():
            sv.append(dict(link=f'{num}/{den}', d=r.logR_d, s=r.logR_s, n=min(r.N_d, r.N_s)))
    sv = pd.DataFrame(sv)
    for l, g in sv.groupby('link'):
        df_ = g.d - g.s
        print(f'S3 DESI-SDSS {l}: {len(g)} common cells, median {np.median(df_):+.4f} dex, rms {np.sqrt(np.mean(df_**2)):.4f} dex')

    # S4 jump at 0.6 vs 1/R3
    jp = os.path.join(a.qa, 'norm_jumps.csv'); s4 = None
    if os.path.exists(jp):
        J = pd.read_csv(jp); J = J[np.isclose(J.z_boundary, 0.6)]
        T = {'passive': [3, 4], 'green': [2], 'starforming': [0, 1]}
        x3 = Ld[('W3050', 'W4700')]; x3 = x3[(x3.ZB == 11) | (x3.ZB == 12)]
        s4 = []
        for t, g in J.groupby('type'):
            xs = x3[x3.DB.isin(T[t])]
            if len(xs):
                s4.append((t, float(np.median(g.k)), 10**-wmed(xs.logR.to_numpy(), xs.N.to_numpy(float))))
        for t, k, ir in s4:
            print(f'S4 jump at z=0.6, {t}: k {k:.2f} vs 1/R3 (z 0.55-0.65) {ir:.2f}')

    # figure
    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5)); cm = plt.get_cmap('coolwarm')
    a0 = ax[0, 0]
    for db in range(5):
        y = [10**fac[db][w] for w in WINS]; a0.plot([LAM[w] for w in WINS], y, 'o-', color=cm(db / 4), label=f'DB{db} Dn {DBLAB[db]}')
    for w in WINS:
        a0.axvspan(LAM[w] - 50, LAM[w] + 50, color='0.9', zorder=0)
    a0.set_xscale('log'); a0.set_yscale('log'); a0.axhline(1, color='0.5', lw=0.6)
    a0.set_xlabel('window centre, rest [A]'); a0.set_ylabel('F(window) / F(5400-5500)')
    a0.set_title('normalization factors by Dn4000 (chained window ratios)', fontsize=10); a0.legend(fontsize=8)
    for t_, (w, zr) in zip(WINS, [('6000-6250', 'z<0.10'), ('5400-5500', '0.10-0.25'), ('4600-4800', '0.25-0.60'), ('3000-3100', '0.60-1.45'), ('2100-2200', '>=1.45')]):
        a0.text(LAM[t_], 0.035, f'{w}\n{zr}', ha='center', fontsize=7)
    a0.set_ylim(0.03, 3)
    a1 = ax[0, 1]
    for i, (num, den) in enumerate(LINKS):
        x = Ld[(num, den)]
        a1.scatter(x.Z, 10**x.logR, c=x.DB, cmap='coolwarm', vmin=0, vmax=4, s=4 + 30 * x.N / x.N.max(), alpha=0.6)
    for z in (0.10, 0.25, 0.60, 1.45):
        a1.axvline(z, color='k', lw=0.8)
    a1.set_yscale('log'); a1.set_xlabel('z (cell median)'); a1.set_ylabel('ratio of the neighbouring windows')
    a1.set_title('per-stack ratios R1 5450/6125, R2 4700/5450, R3 3050/4700, R4 2150/3050 (colour DB)', fontsize=9)
    a2 = ax[1, 0]
    a2.scatter(10**clo.chain, 10**clo.direct, c=clo.DB, cmap='coolwarm', vmin=0, vmax=4, s=10 + 60 * clo.n / clo.n.max())
    lim = [0.05, 2.5]; a2.plot(lim, lim, 'k', lw=0.8); a2.set_xscale('log'); a2.set_yscale('log'); a2.set_xlim(lim); a2.set_ylim(lim)
    a2.set_xlabel('chain R2 x R3 (other stacks)'); a2.set_ylabel('direct F(3050)/F(5450)')
    a2.set_title(f'S2 closure: rms {np.sqrt(np.mean(dd**2)):.3f} dex over {len(clo)} (DB,EB) types', fontsize=10)
    a3 = ax[1, 1]
    for l, g in sv.groupby('link'):
        a3.scatter(10**g.s, 10**g.d, s=8, label=f'{l} ({len(g)})', alpha=0.7)
    lim = [0.1, 3]; a3.plot(lim, lim, 'k', lw=0.8); a3.set_xscale('log'); a3.set_yscale('log'); a3.set_xlim(lim); a3.set_ylim(lim)
    a3.set_xlabel('SDSS v2 stack'); a3.set_ylabel('DESI v1 stack'); a3.legend(fontsize=8); a3.set_title('S3 same cell, DESI vs SDSS', fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(a.qa, 'gal_norm_factors.png'), dpi=105); plt.close(fig)
    print('figure', os.path.join(a.qa, 'gal_norm_factors.png'))


if __name__ == '__main__':
    main()
