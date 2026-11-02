#!/usr/bin/env python3
"""Telluric step 2 for the stars (Nao 2026-11-01; programs_telluric/tellapply.py is the recipe all threads share). step6_stack.py --telluric calls apply() per star, on the rest grid at the star's own
observed wavelengths (log10(1 + z) rounded to 0.1 pixel, enough for the smooth vectors and it lets the grids be cached). Same stars as v2.
Model: the star's v2 stack with the mean telluric divided out. Telluric does not smear at z ~ 0, so the stack carries T0; on SDSS it carries a mix of the DR8 and DR17 T0 by the stack's share of members.
Then r = f / model - 1 - T0 is what tellapply fits, as for galaxies and quasars. A stack's own mean a_k cannot be seen against itself; it is small only if its stars span many observations.
Outputs in OUT: the stacks and EM-PCA bases of steps 6 and 7, and tell_coeff_<tag>.parquet (SID, SURVEY, LOGZ, A0-A6, NPIX, DCHI2, S); step 7 divides the v2 individuals by the same factors.
"""
import os, sys
import numpy as np, pandas as pd, fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path[:0] = [HERE, os.path.join(os.path.dirname(HERE), 'programs_telluric')]
import tellapply, starlib
SWAVE = starlib.SWAVE; V2 = os.path.expanduser('~/data/sdssdesi_stars/v2'); OUT = '/Volumes/exdisk3/sdssdesi_stars/v2_tell'; QUANT = 1e-5


def v2_sids(tag):
    return set(pd.read_csv(os.path.join(V2, f'members_{tag}.csv'), usecols=['SID']).SID.astype(str))


def survey_of(fam, row):
    return 'desi' if fam == 'desi' else tellapply.sdss_survey(row.MJD)


def models(tag, fam, shares):
    """stack name -> v2 stack / (1 + mean T0 of its members). shares: stack name -> fraction of DR8 members (SDSS)."""
    F = fitsio.FITS(os.path.join(V2, f'stacks_{tag}.fits')); fl = F['FLUX'].read().astype(float); names = np.char.strip(F['CELLS'].read()['STACK'].astype(str)); F.close()
    t0 = {s: tellapply.basis(s).grid(SWAVE)['C0'] - 1 for s in (['desi'] if fam == 'desi' else ['dr8', 'dr17'])}; out = {}
    for n, f in zip(names, fl):
        m = t0['desi'] if fam == 'desi' else shares.get(n, 0.0) * t0['dr8'] + (1 - shares.get(n, 0.0)) * t0['dr17']
        out[str(n)] = np.where(np.isfinite(f) & (f > 0), f / (1 + m), np.nan)
    return out


def dr8_shares(df):
    return (df.MJD < 55025).groupby(df.STACK).mean().to_dict() if 'MJD' in df else {}


class StarTelluric:
    def __init__(self, fam, tag, df):
        self.fam, self.tag = fam, tag; self.model = models(tag, fam, dr8_shares(df)); self.rows = []

    def apply(self, row, spec, w):
        S = survey_of(self.fam, row); z = row.get('Z_TRUE', row.Z); z = np.clip(z if np.isfinite(z) else 0.0, -0.004, 0.004); q = round(np.log10(1 + z) / QUANT) * QUANT
        C, a, i = tellapply.basis(S).correct(SWAVE * 10**q, spec, w, self.model[row.STACK]); self.rows.append((row.SID, S, q, *a, i['npix'], i['dchi2'], i['s']))
        return spec / C, w * C**2

    def write(self, outdir):
        n = len(tellapply.basis('desi' if self.fam == 'desi' else 'dr8').T); t = pd.DataFrame(self.rows, columns=['SID', 'SURVEY', 'LOGZ'] + [f'A{k}' for k in range(n)] + ['NPIX', 'DCHI2', 'S'])
        os.makedirs(outdir, exist_ok=True); t.to_parquet(os.path.join(outdir, f'tell_coeff_{self.tag}.parquet')); print(f'telluric: {len(t):,} stars, {int((t.NPIX > 0).sum()):,} fitted', flush=True)


def individuals(tag, v2dir, k, snmin):
    """the cached calibrated individuals of v2 divided by the factors step 6 fitted for the same stars."""
    z = np.load(os.path.join(v2dir, f'individuals_{tag}_{k}_{snmin:g}.npz')); F = z['FLUX'].astype(float); W = z['W'].astype(float); ids = z['IDS'].astype(str)
    t = pd.read_parquet(os.path.join(OUT, f'tell_coeff_{tag}.parquet')).set_index('SID').reindex(ids); cols = [c for c in t.columns if c[0] == 'A' and c[1:].isdigit()]; miss = 0
    for i, (s, q) in enumerate(zip(t.SURVEY, t.LOGZ)):
        if not isinstance(s, str): miss += 1; continue
        b = tellapply.basis(s); C = b.factor_of(b.grid(SWAVE * 10**q), t[cols].iloc[i].to_numpy(float)); F[i] /= C; W[i] *= C**2
    print(f'  individuals: {len(ids)} from v2, {miss} without a factor', flush=True); return F, W, ids.tolist(), z['STACKS'].astype(str).tolist()
