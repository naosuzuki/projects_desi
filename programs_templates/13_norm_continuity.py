#!/usr/bin/env python3
"""Normalization check of the galaxy cell stacks (Nao 2026-10-30): are the stacks continuous across the redshift slices?

Every spectrum, and so every stack, is scaled to unit mean F_lambda in ONE rest window per Delta z = 0.05 slice (normwin):
  z 0.00-0.10 6000-6250 | 0.10-0.25 5400-5500 | 0.25-0.60 4600-4800 | 0.60-1.45 3000-3100 | 1.45-1.60 2100-2200
and the M axis is band-switched (M_R z<0.35, M_G <0.72, M_U <1.45, M_NUV), the colour axis switches g-r -> r-z at z = 0.5.

Groups: one spectral type = (EB, DB) (line EW and Dn4000 bins mean the same at all z), split by MB; within a (ZB, MB, EB, DB)
group the cells (colour bins) are combined with weights N_STACKED (they share the slice window, so this is legitimate).

Figure 1 (stacks_by_M_<type>.png): one panel per MB, the group stacks of every z slice as stored (colour = z), so a jump
  between neighbouring slices is visible directly.
Figure 2 (norm_jumps.png): for each pair of neighbouring slices and each (MB, type) with both present, the ratio
  k = median(stack_{zb+1} / stack_{zb}) over their common rest pixels (weighted by both IVARs). k = 1 means the two
  slices join continuously; a window switch shifts k by the type's colour between the two windows. Vertical lines mark
  the window, M-band and colour switches.
Also prints the per-boundary median k and writes norm_jumps.csv.

Usage: python3 13_norm_continuity.py [--stacks ~/data/sdssdesi_passC_resp/stacks] [--out <stacks>/../qa]
"""
import os, sys, glob, argparse
import numpy as np
import pandas as pd
import fitsio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import normwin, gridcells as gc

TYPES = {'passive': dict(EB=[0, 1], DB=[3, 4]), 'green': dict(EB=[2, 3], DB=[2]), 'starforming': dict(EB=[4, 5, 6], DB=[0, 1])}
SWITCH_WIN = [0.10, 0.25, 0.60, 1.45]
SWITCH_M = [0.35, 0.72, 1.45]
SWITCH_C = [0.50]


def load(stacks):
    rows = []
    for fn in sorted(glob.glob(os.path.join(stacks, 'stacks_z*.fits'))):
        F = fitsio.FITS(fn); f = F['FLUX'].read(); iv = F['IVAR'].read(); c = F['CELLS'].read(); w = F['WAVE'].read(); F.close()
        for i in range(len(c)):
            rows.append(dict(ZB=int(c['ZB'][i]), MB=int(c['MB'][i]), EB=int(c['EB'][i]), DB=int(c['DB'][i]),
                             N=int(c['N_STACKED'][i]), f=f[i], iv=iv[i]))
    return w, rows


def groups(rows):
    """(type, ZB, MB) -> (flux, ivar, N): N-weighted mean of the colour cells, ivar summed (good pixels only)."""
    out = {}
    for t, sel in TYPES.items():
        for r in rows:
            if r['EB'] not in sel['EB'] or r['DB'] not in sel['DB']:
                continue
            k = (t, r['ZB'], r['MB'])
            g = r['iv'] > 0
            if k not in out:
                out[k] = [np.zeros_like(r['f'], float), np.zeros_like(r['f'], float), np.zeros_like(r['f'], float), 0]
            o = out[k]; o[0] += np.where(g, r['N'] * r['f'], 0); o[1] += np.where(g, r['N'], 0); o[2] += r['iv']; o[3] += r['N']
    return {k: (np.where(v[1] > 0, v[0] / np.where(v[1] > 0, v[1], 1), 0), v[2], v[3]) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stacks', default=os.path.expanduser('~/data/sdssdesi_passC_resp/stacks'))
    ap.add_argument('--out', default=None); ap.add_argument('--minn', type=int, default=200)
    a = ap.parse_args()
    out = a.out or os.path.join(os.path.dirname(a.stacks.rstrip('/')), 'qa'); os.makedirs(out, exist_ok=True)
    w, rows = load(a.stacks)
    G = {k: v for k, v in groups(rows).items() if v[2] >= a.minn}
    zc = lambda zb: (zb + 0.5) * gc.DZ
    win = {zb: normwin.slice_window(zb * gc.DZ, (zb + 1) * gc.DZ) for zb in range(32)}

    # Figure 1: stacks by M, per type
    cmap = plt.get_cmap('viridis')
    for t in TYPES:
        mbs = sorted({k[2] for k in G if k[0] == t})
        if not mbs:
            continue
        fig, axs = plt.subplots(len(mbs), 1, figsize=(11, 2.3 * len(mbs)), sharex=True, squeeze=False)
        for ax, mb in zip(axs[:, 0], mbs):
            for (tt, zb, m), (f, iv, n) in sorted(G.items(), key=lambda x: x[0][1]):
                if tt != t or m != mb:
                    continue
                g = iv > 0
                ax.plot(w[g], f[g], lw=0.5, color=cmap(zb / 31), label=f'z {zb*gc.DZ:.2f} N {n}')
                lo, hi = win[zb]; ax.axvspan(lo, hi, color=cmap(zb / 31), alpha=0.08)
            ax.set_xscale('log'); ax.set_ylim(0, 3.0); ax.set_ylabel(f'M {mb}..{mb+1}')
            ax.text(0.01, 0.85, f'{t}, MB {mb}', transform=ax.transAxes, fontsize=8)
        axs[-1, 0].set_xlabel('rest wavelength [A]')
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1.6)); fig.colorbar(sm, ax=axs[:, 0], label='z', fraction=0.02)
        fig.suptitle(f'{t} (EB {TYPES[t]["EB"]}, DB {TYPES[t]["DB"]}): stacks as stored, per-slice window shaded')
        fig.savefig(os.path.join(out, f'stacks_by_M_{t}.png'), dpi=110); plt.close(fig)

    # Figure 2: jump between neighbouring slices
    J = []
    for (t, zb, mb), (f1, i1, n1) in G.items():
        k2 = (t, zb + 1, mb)
        if k2 not in G:
            continue
        f2, i2, n2 = G[k2]
        g = (i1 > 0) & (i2 > 0) & (f1 > 0)
        if g.sum() < 200:
            continue
        r = f2[g] / f1[g]; wt = 1 / (1 / i1[g] / f1[g]**2 + 1 / i2[g] / f2[g]**2)
        o = np.argsort(r); c = np.cumsum(wt[o]); k = r[o][np.searchsorted(c, 0.5 * c[-1])]
        J.append(dict(type=t, MB=mb, ZB=zb, z_boundary=(zb + 1) * gc.DZ, k=k, N1=n1, N2=n2, npix=int(g.sum()),
                      win1=str(win[zb]), win2=str(win[zb + 1])))
    J = pd.DataFrame(J); J.to_csv(os.path.join(out, 'norm_jumps.csv'), index=False)
    fig, ax = plt.subplots(figsize=(11, 4.5))
    mk = {'passive': 'o', 'green': 's', 'starforming': '^'}
    for t, d in J.groupby('type'):
        sc = ax.scatter(d.z_boundary, d.k, c=d.MB, cmap='coolwarm', marker=mk[t], s=22, label=t, vmin=-24, vmax=-16)
    for z in SWITCH_WIN: ax.axvline(z, color='k', lw=1)
    for z in SWITCH_M: ax.axvline(z, color='r', ls='--', lw=1)
    for z in SWITCH_C: ax.axvline(z, color='g', ls=':', lw=1)
    ax.axhline(1, color='0.5'); ax.set_ylim(0.6, 1.6); ax.set_xlabel('slice boundary z'); ax.set_ylabel('k = stack(z+dz) / stack(z)')
    ax.legend(loc='upper left', fontsize=8); fig.colorbar(sc, label='MB')
    ax.set_title('jump between neighbouring z slices (black: window switch, red dashed: M band, green dotted: colour)')
    fig.tight_layout(); fig.savefig(os.path.join(out, 'norm_jumps.png'), dpi=110); plt.close(fig)
    s = J.groupby('z_boundary').k.agg(k_median='median', n='count', k_min='min', k_max='max')
    print(s.round(3).to_string()); print('outputs in', out)


if __name__ == '__main__':
    main()
