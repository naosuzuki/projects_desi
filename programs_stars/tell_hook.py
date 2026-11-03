#!/usr/bin/env python3
"""Telluric step 2 for the stars (Nao 2026-11-01; programs_telluric/tellapply.py is the recipe all threads share). step6_stack.py --telluric calls apply() per star, on the rest grid at the star's own
observed wavelengths (log10(1 + z) rounded to 0.1 pixel, enough for the smooth vectors and it lets the grids be cached). Same stars as v2.
Model: the star's v2 stack with the T0 of its members divided out: each member's T0 (DR8 or DR17 on SDSS) at its own velocity, -3 to +2 pixels, averaged. Dividing by T0 at rest moved the DR8 stack means of a_k by about 0.1 % in flux (08_smear_ak.py; SDSS restacked 2026-11-03). The DESI tell_coeff still holds a_k from that older model, unused since C = 1 + T0 there.
Then r = f / model - 1 - T0 is what tellapply fits, as for galaxies and quasars. A stack's own mean a_k cannot be seen against itself; it is small only if its stars span many observations.
OUT is v2_tell unless STAR_TELL_OUT says otherwise (stellar v5: /Volumes/exdisk3/sdssdesi_stars/v5). Outputs in OUT: the stacks and EM-PCA bases of steps 6 and 7, and tell_coeff_<tag>.parquet (SID, SURVEY, LOGZ, A0-A6, NPIX, DCHI2, S); step 7 divides the v2 individuals by the same factors.
"""
import os, sys
import numpy as np, pandas as pd, fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path[:0] = [HERE, os.path.join(os.path.dirname(HERE), 'programs_telluric')]
import tellapply, starlib
SWAVE = starlib.SWAVE; V2 = os.path.expanduser('~/data/sdssdesi_stars/v2'); OUT = os.environ.get('STAR_TELL_OUT', '/Volumes/exdisk3/sdssdesi_stars/v2_tell'); QUANT = 1e-5


def v2_sids(tag):
    """the stars of v2, less OUT/drop_<tag>.csv if it exists (stellar v5: DESI stars without BP or RP, and cells that fall below 100 without them)."""
    s = set(pd.read_csv(os.path.join(V2, f'members_{tag}.csv'), usecols=['SID']).SID.astype(str)); p = os.path.join(OUT, f'drop_{tag}.csv')
    return s - set(pd.read_csv(p, usecols=['SID']).SID.astype(str)) if os.path.exists(p) else s


def survey_of(fam, row):
    return 'desi' if fam == 'desi' else tellapply.sdss_survey(row.MJD)


def logz(z):
    """log10(1 + z) of stars, clipped as step 6 clips them and rounded to QUANT."""
    z = np.clip(np.nan_to_num(np.asarray(z, float), nan=0.0), -0.004, 0.004); return np.round(np.log10(1 + z) / QUANT) * QUANT


def smeared_t0(stack, survey, q):
    """stack name -> T0 its members carried: each star's own survey at its own velocity q, averaged. Every star counts once in a stack (step 6, w / med)."""
    g = pd.DataFrame(dict(STACK=stack, S=survey, Q=np.round(q, 6))).groupby(['STACK', 'S', 'Q']).size().rename('n').reset_index(); ll = np.log10(SWAVE); cache = {}; out = {}
    for name, x in g.groupby('STACK'):
        y = np.zeros(len(SWAVE))
        for s, u, n in zip(x.S, x.Q, x.n):
            if (s, u) not in cache: cache[(s, u)] = np.interp(ll + u, ll, tellapply.basis(s).Ffull[0], left=0.0, right=0.0)
            y += n * cache[(s, u)]
        out[str(name)] = y / x.n.sum()
    return out


def models(tag, fam, df):
    """stack name -> v2 stack / (1 + the T0 of its members, smeared_t0). df: the stars of the stacks."""
    F = fitsio.FITS(os.path.join(V2, f'stacks_{tag}.fits')); fl = F['FLUX'].read().astype(float); names = np.char.strip(F['CELLS'].read()['STACK'].astype(str)); F.close()
    S = np.full(len(df), 'desi') if fam == 'desi' else np.where(df.MJD < 55025, 'dr8', 'dr17'); T = smeared_t0(df.STACK.to_numpy(), S, logz(df['Z_TRUE'] if 'Z_TRUE' in df else df.Z))
    return {str(n): np.where(np.isfinite(f) & (f > 0), f / (1 + T[str(n)]), np.nan) for n, f in zip(names, fl) if str(n) in T}


class StarTelluric:
    def __init__(self, fam, tag, df):
        self.fam, self.tag = fam, tag; self.model = models(tag, fam, df); self.rows = []

    def apply(self, row, spec, w):
        S = survey_of(self.fam, row); q = float(logz(row.get('Z_TRUE', row.Z)))
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
    keep = t.SURVEY.map(lambda s: isinstance(s, str)).to_numpy(); print(f'  individuals: {len(ids)} from v2, {miss} without a factor, left out', flush=True); return F[keep], W[keep], ids[keep].tolist(), z['STACKS'].astype(str)[keep].tolist()
