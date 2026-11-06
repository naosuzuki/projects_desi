#!/usr/bin/env python3
"""Figures and table of the 2D HR groups (hr_build.py; Nao 2026-11-03). Luminosity bins per group, e0 to e9 per group, and the two WD branches at the same colour with their H and He I equivalent widths.
STAR_HR=single reads hr2d_single and writes plots_star_hr2dsingle_*; STAR_HR=all reads hr2d_all (version 8, hr_all.py) and writes plots_star_hr2dall_*. Output: plots_luminosity/plots_star_hr2d_lum.png, plots_star_hr2d_empca_{wd,ms,other}.png, plots_star_hr2d_wdbranch.png; csv_files/plots_star_hr2d_summary.csv, plots_star_hr2d_wdbranch.csv
"""
import os, sys
import numpy as np, pandas as pd, fitsio
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path[:0] = [HERE, os.path.join(ROOT, 'programs_telluric')]
import plot_lum, tellapply
MODE = os.environ.get('STAR_HR', ''); SINGLE = MODE == 'single'; OUT = '/Volumes/exdisk3/sdssdesi_stars/hr2d' + {'single': '_single', 'all': '_all'}.get(MODE, ''); PRE = 'hr2d' + MODE; PNG = os.path.join(ROOT, 'plots_luminosity'); CSV = os.path.join(ROOT, 'csv_files'); WAVE = 10**(3.5502 + 1e-4 * np.arange(4447)); BLUE, GRAY, BAND, INK, MUTED, RED = '#2a78d6', '#a8a79f', '#ecebe6', '#0b0b0b', '#52514e', '#d6402a'
LAB = dict(WDHOT='WD hot', WDA='WD bright branch', WDB='WD faint branch', WDCOOL='WD cool', MSF='MS 0.55-0.75', MSG='MS 0.75-1.05', MSK='MS 1.05-1.9', MSM0='MS 1.9-2.7', MSM5='MS > 2.7', BINARY='binary sequence', SD='below MS', SGB='subgiants', WDHE='WD with He I', GIANT='giants', BLUE='blue, above WD')
SETS = dict(wd=['WDHOT', 'WDA', 'WDB', 'WDCOOL'], ms=['MSF', 'MSG', 'MSK', 'MSM0'], other=['MSM5', 'BINARY', 'SD', 'SGB'])
if MODE == 'all':                         # version 8 (hr_all.py): no BINARY, He I white dwarfs as WD_HE
    SETS = dict(wd=['WDHOT', 'WDHE', 'WDA', 'WDB'], ms=['MSF', 'MSG', 'MSK', 'MSM0'], other=['MSM5', 'SD', 'SGB', 'WDCOOL'])
LINES = {'He I 4471': (4471.5, 4460, 4485), 'He I 5876': (5877.2, 5865, 5890), 'He I 6678': (6679, 6670, 6690), 'H beta': (4862.7, 4840, 4885), 'H alpha': (6564.6, 6545, 6585), 'H delta': (4102.9, 4085, 4120)}


def style(q, size=6):
    q.tick_params(labelsize=size, colors=MUTED, length=2, pad=1.5, width=0.5); [q.spines[s].set_visible(False) for s in ('top', 'right')]; [q.spines[s].set_color(GRAY) for s in ('left', 'bottom')]; [q.spines[s].set_linewidth(0.5) for s in ('left', 'bottom')]


def load():
    F = fitsio.FITS(f'{OUT}/stacks_desi_hr2d.fits'); c = F['CELLS'].read(); fl = F['FLUX'].read().astype(float); ns = F['NSPEC'].read().astype(float); F.close(); s = np.char.strip(c['STACK'].astype(str)); g = np.array([x.split('_')[2] for x in s]); return c, fl, ns, s, g


def lum(c, fl, ns, g):
    cm = plt.get_cmap('rainbow_r'); norm = plt.Normalize(2, 15); names = [x for v in SETS.values() for x in v]; fig, ax = plt.subplots(3, 4, figsize=(7.3, 6.6), sharex=True)
    for q, n in zip(ax.ravel(), names):
        style(q); k = np.where(g == n)[0]; b = np.floor(c['MG0_MED'][k]).astype(int); lo, hi = np.inf, -np.inf
        for u in np.unique(b):
            i = k[b == u]; w = np.where((ns[i] >= 0.5 * c['N'][i, None]) & (fl[i] > 0), ns[i], 0.0); W = w.sum(0)
            with np.errstate(all='ignore'): y = plot_lum.smooth(np.where(W > 0, (w * np.nan_to_num(fl[i])).sum(0) / W, np.nan)); m = -2.5 * np.log10(y * WAVE**2 / plot_lum.C_AA) - 48.60
            q.plot(WAVE, m, color=cm(norm(u + 0.5)), lw=0.4); f = m[np.isfinite(m)]
            if len(f): lo, hi = min(lo, np.percentile(f, .2)), max(hi, np.percentile(f, 99.9))
        q.set_ylim(hi + .3, lo - .3); q.set_xlim(3550, 9900); q.locator_params(axis='x', nbins=3); q.set_title(f'{LAB[n]}, {int(c["N"][k].sum()):,}', fontsize=6.5, loc='left', color=INK, pad=2)
    for q in ax[:, 0]: q.set_ylabel('M_AB(λ)', fontsize=6.5, color=MUTED)
    for q in ax[-1]: q.set_xlabel('rest wavelength [Å]', fontsize=6, color=MUTED)
    fig.tight_layout(pad=.4, w_pad=.3, h_pad=.5, rect=(0, 0, .93, 1)); cax = fig.add_axes([.945, .12, .012, .76]); cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cm, norm=norm), cax=cax); cb.set_label('M_G0 of the bin', fontsize=6.5, color=MUTED); cb.ax.tick_params(labelsize=6, colors=MUTED, length=2); cb.outline.set_visible(False); fig.savefig(f'{PNG}/plots_star_{PRE}_lum.png', dpi=220); plt.close(fig)


def eigen(tag, names):
    band = [(n[0], n[-1]) for n in tellapply.basis('desi').bands]; fig, ax = plt.subplots(10, 4, figsize=(7.3, 8.6), sharex=True)
    for j, n in enumerate(names):
        p = f'{OUT}/empca_desi_{n}.fits'; M, E, V, P = (fitsio.read(p, ext=e).astype(float) for e in ('MEAN', 'EIGVEC', 'VARIANCE', 'PIXMASK')); A = np.where(P > 0, np.vstack([M, E[:9]]), np.nan); r2 = np.diff(np.r_[0, V]); ok = np.isfinite(A[0]); ax[0, j].set_title(LAB[n], fontsize=7, loc='left', color=INK, pad=3)
        for k in range(10):
            q = ax[k, j]; style(q, 5)
            for lo, hi in band: q.axvspan(lo, hi, color=BAND, lw=0, zorder=0)
            if k: q.axhline(0, color=GRAY, lw=.3, zorder=1)
            q.plot(WAVE, A[k], color=BLUE if k else INK, lw=.25, zorder=2); lo, hi = np.nanpercentile(A[k], [.2, 99.8]); q.set_ylim(lo - .12 * (hi - lo), hi + .25 * (hi - lo)); q.set_xlim(WAVE[ok][0], WAVE[ok][-1]); q.locator_params(axis='y', nbins=3); q.locator_params(axis='x', nbins=3); q.text(.985, .96, 'e0' if k == 0 else f'e{k}  +{r2[k - 1]:.4f}', transform=q.transAxes, fontsize=5, color=INK, ha='right', va='top')
        ax[-1, j].set_xlabel('rest wavelength [Å]', fontsize=6, color=MUTED)
    fig.tight_layout(pad=.3, w_pad=.2, h_pad=.15); fig.savefig(f'{PNG}/plots_star_{PRE}_empca_{tag}.png', dpi=220); plt.close(fig)


def ew(f, lam, lo, hi):
    i = (WAVE >= lo) & (WAVE <= hi); c = ((WAVE >= lo - 15) & (WAVE < lo)) | ((WAVE > hi) & (WAVE <= hi + 15)); C = np.nanmedian(f[c]); return np.nansum((1 - f[i] / C) * np.gradient(WAVE)[i])


def branches(c, fl, ns, s, g):
    bins = [(.1, .2), (.2, .3), (.3, .4), (.4, .5), (.5, .75)]; col = c['COL0_MED']; rows = []; fig, ax = plt.subplots(len(bins), 1, figsize=(7.3, 7.5), sharex=True); out = {}
    for q, (a, b) in zip(ax, bins):
        style(q)
        for n, cc in (('WDA', BLUE), ('WDB', RED)):
            i = np.where((g == n) & (col >= a) & (col < b))[0]
            if not len(i): continue
            w = np.where((ns[i] >= .5 * c['N'][i, None]) & (fl[i] > 0), ns[i], 0.0); W = w.sum(0); nm = np.nanmedian(fl[i][:, (WAVE > 5400) & (WAVE < 5500)], axis=1)[:, None]
            with np.errstate(all='ignore'): y = np.where(W > 0, (w * np.nan_to_num(fl[i] / nm)).sum(0) / W, np.nan)
            q.plot(WAVE, y, color=cc, lw=.45, label=LAB[n]); rows.append(dict(colour=f'{a}-{b}', group=n, stacks=len(i), stars=int(c['N'][i].sum()), **{k: ew(y, *v) for k, v in LINES.items()}))
        q.set_title(f'(BP-RP)0 {a}-{b}', fontsize=7, loc='left', color=INK, pad=2); q.set_xlim(3700, 7000)
        for k, v in LINES.items(): q.axvline(v[0], color=GRAY if k[0] == 'H' and k[1] == ' ' else '#7a3b1d', lw=.3, zorder=0)
    ax[0].legend(fontsize=7, frameon=False, ncol=2, loc='lower right'); ax[-1].set_xlabel('rest wavelength [Å]; brown: He I, gray: H', fontsize=7, color=MUTED); fig.supylabel('F_λ, scaled in 5400-5500 Å', fontsize=7, color=MUTED); fig.tight_layout(pad=.4); fig.savefig(f'{PNG}/plots_star_{PRE}_wdbranch.png', dpi=200); plt.close(fig); t = pd.DataFrame(rows); t.to_csv(f'{CSV}/plots_star_{PRE}_wdbranch.csv', index=False); return t


def main():
    pd.set_option('display.width', 180); pd.set_option('display.max_columns', 30); c, fl, ns, s, g = load(); lum(c, fl, ns, g)
    for k, v in SETS.items(): eigen(k, v)
    print(branches(c, fl, ns, s, g).round(2).to_string(index=False)); r = []
    for n in [x for x in LAB if os.path.exists(f'{OUT}/empca_desi_{x}.fits')]:
        p = f'{OUT}/empca_desi_{n}.fits'; V = fitsio.read(p, ext='VARIANCE'); k = np.char.strip(fitsio.read(p, ext='NORM')['KIND'].astype(str)); i = g == n; d = np.diff(np.r_[0, V]); N = ns[i].max(axis=1).astype(float); r.append(dict(group=n, stacks=int(i.sum()), stars=int(c['N'][i].sum()), individuals=int((k == 'star').sum()), N_min=int(N.min()), N_med=int(np.median(N)), N_max=int(N.max()), e1=d[0], e2=d[1], e3=d[2], r2_10=V[9], colmin=c['COL0_MED'][i].min(), colmax=c['COL0_MED'][i].max(), mgmin=c['MG0_MED'][i].min(), mgmax=c['MG0_MED'][i].max()))
    t = pd.DataFrame(r); t.to_csv(f'{CSV}/plots_star_{PRE}_summary.csv', index=False); print(t.round(3).to_string(index=False))


if __name__ == '__main__':
    main()
