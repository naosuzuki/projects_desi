#!/usr/bin/env python3
"""Telluric-corrected (v3t) vs uncorrected (v3s) galaxy stacks, in the telluric windows (Nao 2026-11-03).
For every stack cell present in both releases, r = F_v3t / F_v3s - 1 is moved to the observed frame at the stack's median
redshift and the median over stacks is taken per observed pixel, per survey (DESI, SDSS). A telluric correction fixed in the
observed frame shows up as a sharp feature there; its rms inside the O2 / H2O bands vs outside is printed.
Also a summary of the per-spectrum TELL_A0..6 / TELL_CHI2R columns over all v3t pass tables (median, 16-84 %).
Outputs: plots_luminosity/plots_gal_tell_compare_<V>.png, csv_files/plots_gal_tell_compare_<V>.csv,
         csv_files/gal_tell_amp_v3t.csv.
Usage: python3 28_tell_compare.py --version G [--no-tables]
"""
import os, glob, argparse
import numpy as np
import pandas as pd
import fitsio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PNG_DIR = os.path.join(ROOT, 'plots_luminosity'); CSV_DIR = os.path.join(ROOT, 'csv_files')
H = os.path.expanduser('~/data')
OBS = np.arange(3550.0, 10000.0, 1.0)
BANDS = {'B band': (6860, 6960), 'A band': (7590, 7700), 'H2O 7200': (7160, 7340), 'H2O 8200': (8100, 8400),
         'H2O 9000': (8900, 9700)}
KEY = ['ZB', 'MB', 'CB', 'EB', 'DB', 'SET']


def load(d):
    out = {}
    for fn in sorted(glob.glob(os.path.join(d, 'stacks_z*.fits'))):
        F = fitsio.FITS(fn); f = F['FLUX'].read(); ns = F['NSPEC'].read(); w = F['WAVE'].read(); c = F['CELLS'].read(); F.close()
        for i in range(len(c)):
            k = tuple(c[j][i] for j in KEY)
            N = max(int(c['N_STACKED'][i]), 1)
            good = np.isfinite(f[i]) & (ns[i] >= 0.5 * N) & (f[i] > 0)
            out[k] = (np.where(good, f[i], np.nan), float(c['Z_MED'][i]), str(c['SET'][i]))
    return out, w


def compare(dirs_s, dirs_t):
    S, T = {}, {}
    for d in dirs_s:
        o, w = load(d); S.update(o)
    for d in dirs_t:
        o, w = load(d); T.update(o)
    rows = {'DESI': [], 'SDSS': []}
    for k in S.keys() & T.keys():
        fs, z, lab = S[k]; ft = T[k][0]
        r = ft / fs - 1
        if np.isfinite(r).sum() < 500:
            continue
        rows[lab].append(np.interp(OBS, w * (1 + z), r, left=np.nan, right=np.nan))
    return {lab: np.array(v) for lab, v in rows.items()}


def band_stats(mo):
    inb = np.zeros(len(OBS), bool)
    s = {}
    for name, (a, b) in BANDS.items():
        m = (OBS >= a) & (OBS < b); inb |= m
        s[name] = np.nanstd(mo[m]), np.nanmedian(mo[m]), np.nanmin(mo[m]), np.nanmax(mo[m])
    s['outside'] = np.nanstd(mo[~inb]), np.nanmedian(mo[~inb]), np.nanmin(mo[~inb]), np.nanmax(mo[~inb])
    return s


def tables():
    cols = [f'TELL_A{k}' for k in range(7)] + ['TELL_CHI2R', 'TELL_NFIT']
    res = []
    for lab, pat in [('DESI', f'{H}/sdssdesi_passC_v3t/tables/*/*.parquet'), ('DR8', f'{H}/sdssdesi_passG_sdss_v3t/tables/dr8/*.parquet'),
                     ('DR17', f'{H}/sdssdesi_passG_sdss_v3t/tables/dr17/*.parquet')]:
        parts = []
        for fn in sorted(glob.glob(pat)):
            try:
                parts.append(pd.read_parquet(fn, columns=cols).astype('float32'))
            except Exception:
                pass
        d = pd.concat(parts, ignore_index=True)
        for c in cols:
            q = np.nanpercentile(d[c], [16, 50, 84])
            res.append(dict(SURVEY=lab, COL=c, N=int(np.isfinite(d[c]).sum()), P16=q[0], MED=q[1], P84=q[2],
                            RMS=float(np.nanstd(d[c])), FRAC_NONZERO=float((d[c] != 0).mean())))
        print(f'{lab}: {len(d)} spectra')
    return pd.DataFrame(res)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--version', default='G'); ap.add_argument('--no-tables', action='store_true')
    a = ap.parse_args(); V = a.version
    S = [f'{H}/sdssdesi_passC_v3/stacks_{V}_chain_stacks', f'{H}/sdssdesi_passG_sdss_v3/stacks_{V}_chain_stacks']
    T = [f'{H}/sdssdesi_passC_v3t/stacks_{V}_chain_stacks', f'{H}/sdssdesi_passG_sdss_v3t/stacks_{V}_chain_stacks']
    R = compare(S, T)
    fig, ax = plt.subplots(2, 1, figsize=(16, 7), sharex=True); out = []
    for k, lab in enumerate(['DESI', 'SDSS']):
        ro = R[lab]
        with np.errstate(all='ignore'):
            mo = np.nanmedian(ro, axis=0); no = np.isfinite(ro).sum(0)
        mo[no < 20] = np.nan
        out.append(pd.DataFrame(dict(SURVEY=lab, WAVE_OBS=OBS, MED=mo, N=no)))
        st = band_stats(mo)
        print(f'{lab} version {V}: {len(ro)} matched stacks; median(F_v3t/F_v3s - 1) in the observed frame:')
        for name, (sd, md, lo, hi) in st.items():
            print(f'  {name:9s} rms {sd:.5f}  median {md:+.5f}  min {lo:+.4f}  max {hi:+.4f}')
        ax[k].plot(OBS, mo, lw=0.6, color='k'); ax[k].set_ylim(-0.02, 0.02); ax[k].grid(alpha=0.3)
        ax[k].set_title(f'{lab}, version {V}: {len(ro)} stacks, median of F_v3t / F_v3s - 1 in the OBSERVED frame', fontsize=10)
        for name, (x0, x1) in BANDS.items():
            ax[k].axvspan(x0, x1, color='C3', alpha=0.08); ax[k].text(x0, 0.018, name, fontsize=7, va='top', color='C3')
        ax[k].axvline(9280, color='C0', lw=0.6, ls='--')
    ax[-1].set_xlabel('observed wavelength [A]'); ax[-1].set_xlim(OBS[0], OBS[-1])
    os.makedirs(PNG_DIR, exist_ok=True); os.makedirs(CSV_DIR, exist_ok=True)
    png = os.path.join(PNG_DIR, f'plots_gal_tell_compare_{V}.png'); fig.tight_layout(); fig.savefig(png, dpi=80); plt.close(fig)
    pd.concat(out).to_csv(os.path.join(CSV_DIR, f'plots_gal_tell_compare_{V}.csv'), index=False); print(png)
    if not a.no_tables:
        t = tables(); t.to_csv(os.path.join(CSV_DIR, 'gal_tell_amp_v3t.csv'), index=False)
        print(t.to_string(index=False, float_format=lambda x: f'{x:.4f}'))


if __name__ == '__main__':
    main()
