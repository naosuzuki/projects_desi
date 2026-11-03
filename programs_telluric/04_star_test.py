#!/usr/bin/env python3
"""Telluric step 2 on stars, before any restack: per class, what does the fit of tellapply gain on held-out pixels against the same eigenvectors moved off the atmosphere (--shift pixels), and does it
absorb stellar features? Sample: the stars of v2 in --nfiles random coadds / plates, from the stores, calibrated by step 6's machinery; model as tell_hook.
Feature test: r_TiO = mean f / model - 1 at rest 7055-7140 A minus 6620-6800 A, outside every telluric band: how much stronger a star's TiO is than its stack's. A rank correlation of a_k with it within a class
means the fit is following the star, not the atmosphere. r_CaT (Ca II triplet against 8400-8480 A, no band) does the same for the warmer classes.
Output: csv_files/plots_telluric_startest_<survey>.csv.  Usage: python3 04_star_test.py --survey desi|sdss [--nfiles 1500] [--shift 20]
"""
import os, sys, argparse, importlib.util, time
import numpy as np, pandas as pd, fitsio
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); STARS = os.path.join(ROOT, 'programs_stars'); sys.path[:0] = [HERE, STARS]
import tellapply, tell_hook, starlib, star_hook
from step7_empca import star_class
sp = importlib.util.spec_from_file_location('step6', os.path.join(STARS, 'step6_stack.py')); S6 = importlib.util.module_from_spec(sp); sp.loader.exec_module(S6)
def win(*r):
    return np.any([(SWAVE > lo) & (SWAVE < hi) for lo, hi in r], axis=0)


SWAVE = starlib.SWAVE; TIO = win((7055, 7140)); TIO_C = win((6620, 6800)); CAT = win((8490, 8506), (8532, 8552), (8652, 8672)); CAT_C = win((8400, 8480))


def rwin(spec, w, m, mask):
    """chi2-weighted mean of f / m - 1 in a window, without a common scale: differences of two windows measure the feature against its local continuum."""
    g = mask & (w > 0) & np.isfinite(m) & (m > 0)
    if g.sum() < 10: return np.nan
    return np.sum((w * m * (spec - m))[g]) / np.sum((w * m * m)[g])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--survey', default='desi', choices=['desi', 'sdss']); ap.add_argument('--nfiles', type=int, default=1500); ap.add_argument('--shift', type=int, default=20)
    a = ap.parse_args(); tag = f'{a.survey}_0.5x0.1_gspphot'; rng = np.random.default_rng(5); t0 = time.time(); keep = tell_hook.v2_sids(tag)
    if a.survey == 'desi':
        df = pd.read_parquet(os.path.join(S6.CACHE, f'step6_selection_{tag}.parquet')); df = df[df.SID.isin(keep)]; it = lambda d: S6.iter_store(d, 0, star_hook.STORE_V2)
    else:
        cells, mem = S6.members(tag, 100); ns = argparse.Namespace(survey='sdss', passdir=None, maxper=0, select_only=False, snmin=5.0, nmin=100, v2=True, from_store=False, onepass=True, store_only=False, telluric=True, tag=tag, outdir=tell_hook.OUT); df, _ = S6.select_stars(ns, cells, mem, t0); df = df[df.SID.isin(keep)]; it = lambda d: S6.iter_sdss(d, 0)
    model = tell_hook.models(tag, a.survey, df); files = df.FILE.unique(); df = df[df.FILE.isin(rng.choice(files, min(a.nfiles, len(files)), replace=False))]; print(f'{len(df):,} stars in {df.FILE.nunique():,} files, {time.time() - t0:.0f} s', flush=True)
    c = fitsio.read(os.path.join(tell_hook.V2, f'stacks_{tag}.fits'), ext='CELLS'); cls = {s.strip(): star_class(float(m), float(k)) for s, m, k in zip(c['STACK'].astype(str), c['MG0_MED'], c['COL0_MED'])}
    rows = []
    for row, f, iv, mk, rid0 in it(df):
        spec, w, med = S6.star_spectrum(row, f, iv, mk, rid0, raw=True); S = tell_hook.survey_of(a.survey, row); z = row.get('Z_TRUE', row.Z); z = np.clip(z if np.isfinite(z) else 0.0, -0.004, 0.004); q = round(np.log10(1 + z) / tell_hook.QUANT) * tell_hook.QUANT; wave = SWAVE * 10**q; m = model[row.STACK]
        b = tellapply.basis(S); gain, x, n = b.holdout(wave, spec, w, m); gs, _, _ = tellapply.basis(S, shift=a.shift).holdout(wave, spec, w, m); _, full, info = b.correct(wave, spec, w, m)
        rr = dict(TIO=rwin(spec, w, m, TIO) - rwin(spec, w, m, TIO_C), CAT=rwin(spec, w, m, CAT) - rwin(spec, w, m, CAT_C))
        rows.append((row.SID, S, cls.get(row.STACK, '?'), row.SNR_ALL, gain, gs, n, rr['TIO'], rr['CAT'], *full))
    k = len(rows[0]) - 9; t = pd.DataFrame(rows, columns=['SID', 'SURVEY', 'CLASS', 'SNR', 'GAIN', 'GAIN_SHIFT', 'NHELD', 'R_TIO', 'R_CAT'] + [f'A{j}' for j in range(k)]); out = []
    for (s, c_), g in t.groupby(['SURVEY', 'CLASS']):
        g_ = g[g.NHELD > 50]; b = tellapply.basis(s); d = dict(survey=s, cls=c_, n=len(g_), snr_med=round(g_.SNR.median(), 1), gain=round(g_.GAIN.mean(), 2), gain_err=round(g_.GAIN.std() / np.sqrt(max(len(g_), 1)), 2), gain_shift=round(g_.GAIN_SHIFT.mean(), 2), gain_shift_err=round(g_.GAIN_SHIFT.std() / np.sqrt(max(len(g_), 1)), 2))
        for j in b.sel:
            lab = b.labels[j].replace(' ', ''); ok = np.isfinite(g_.R_TIO); d[f'rms_{lab}'] = round(g_[f'A{j}'].std(), 4); d[f'mean_{lab}'] = round(g_[f'A{j}'].mean(), 4); d[f'rTiO_{lab}'] = round(spearmanr(g_[f'A{j}'][ok], g_.R_TIO[ok])[0], 3) if ok.sum() > 20 else np.nan
            okc = np.isfinite(g_.R_CAT); d[f'rCaT_{lab}'] = round(spearmanr(g_[f'A{j}'][okc], g_.R_CAT[okc])[0], 3) if okc.sum() > 20 else np.nan
        out.append(d)
    o = pd.DataFrame(out); p = os.path.join(ROOT, 'csv_files', f'plots_telluric_startest_{a.survey}.csv'); o.to_csv(p, index=False); print(o.T.to_string()); print(p, f'{time.time() - t0:.0f} s', flush=True)


if __name__ == '__main__':
    main()
