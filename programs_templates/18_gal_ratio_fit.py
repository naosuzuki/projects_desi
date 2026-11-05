#!/usr/bin/env python3
"""Smooth fit of the galaxy normalization link log R = log F(3000-3100)/F(5400-5500) vs absolute magnitude, PER Dn4000 BIN
(DB 0-4; Nao 2026-10-31: do not mix the types), all overlap redshifts pooled, DESI only (the galaxy counterpart of
programs_quasars/11_ratio_fit.py). Per-galaxy ratios as 17_gal_ratio_vs_mag.py (COV >= 0.8 in both windows), from the
DESI normtab of 12_restack_response.py.
Two magnitude axes, each fit separately (one figure each):
  MABS_G     photometric rest g (whole galaxy), noise independent of the window fluxes  -> gal_ratio_fit_MABS_G.png
  M_AB(5450) from the spectrum, F_lambda = W5450/(NORM_SCALE x RESCALE) at 5450 A  -> gal_ratio_fit_MSPEC.png
             (at the denominator, not at 3050 as the quasars: W3050 can be negative and is kept)
Per 0.25-mag bin: median of the LINEAR ratio, negative W3050 kept (7-9 % of bright-program galaxies at z 0.2-0.5 sit at the
DESI blue edge; dropping them biased the medians high), then logged. Only bins holding >= FITFRAC (0.1 %) of the type are fit
(a few thousand dark-program galaxies at z ~ 0.25 make up the whole faint end of the red types alone; shown as open circles).
(0.5-mag bins, as docs/HANDOVER_normalization.md; every fitted bin holds >= 1,000 DESI galaxies.) Fit to those binned medians, weighted by their standard error (1.2533 sigma/sqrt n); models tried:
  line      logR = a + b (M + 21)
  softplus  two straight lines joined smoothly at M0 over width w (as the quasars):
            logR = a + s_bright (M - M0) + (s_faint - s_bright) w ln(1 + exp((M - M0)/w))
  sigmoid   logR = a + d tanh((M - M0)/w) + s (M - M0)
the one with the lowest BIC = chi2 + k ln(nbins) is kept (three starting M0 for the 5-parameter models).
Outputs (<out>): gal_ratio_fit_{MABS_G,MSPEC}.png, gal_ratio_fit.json ({axis: {DB: parameters, covariance, binned medians, residuals}}).
--uniform N: for the relation only, N galaxies per 0.25-mag bin drawn equally over the strata (--strata program_z: program x
z slice, default; z: z slices) of that bin (equal weight per
magnitude; the per-bin number of z slices is printed under the residuals). Outputs then carry _uniform<N>.
Lookup (MABS_G axis, Nao 2026-10-31): per DB, log R on a 0.25-mag grid LOOK_M from the fit inside the fitted range, held at
the end values outside it (SOURCE fit / held_faint / held_bright / held_plateau: DB2, DB3 held at the bright plateau). DB2-4 are fit only brighter than HOLD (the faint side of the
targeting gap holds blue interlopers, not faint red galaxies): gal_ratio_lookup<tag>.csv.
--sdss DIR adds the SDSS records (PROGRAM = DR8 / DR17 as strata); --sdss-mode fill (default) uses them only in the
(M bin, z slice) cells where DESI has < NFILL (50) galaxies, all = everywhere. Outputs then carry _sdss<mode>.
Usage: python3 18_gal_ratio_fit.py [--passc DIR] [--out DIR] [--zmax 9] [--uniform 2000] [--sdss DIR]
"""
import os, sys, glob, json, argparse
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib
pb = importlib.import_module('02_passB_measure')

COVER = 0.8; DM = 0.5; NMIN = 100; FITFRAC = 0.001; NFILL = 50; C_AA = 2.99792458e18
HOLD = {2: -20.0, 3: -20.0, 4: -20.5}
PLATEAU = {2: -21.0, 3: -21.0}           # DB2, DB3 held at the bright plateau (median of fitted bins brighter than this; Nao 2026-10-31)   # MABS_G: red types fit only brighter than this (DESI/SDSS targeting gap), held flat fainter
LOOK_M = (-25.0, -15.5)                 # lookup grid in MABS_G (covers the stack M bins)
LOOK_M_SPEC = (-25.0, -12.0)            # lookup grid in M_AB(5450) from the spectrum (fiber light: 0.3-0.8 mag fainter)
COLS = ['Z', 'DB', 'PROGRAM', 'NORM_SCALE', 'RESCALE', 'MABS_G', 'COV3050', 'COV5450', 'W3050', 'W5450']
DBLAB = {0: 'Dn4000 < 1.15', 1: 'Dn4000 1.15-1.40', 2: 'Dn4000 1.40-1.65', 3: 'Dn4000 1.65-1.90', 4: 'Dn4000 > 1.90'}
AXES = {'MABS_G': 'MABS_G (photometric, whole galaxy)', 'MSPEC': 'M_AB at 5450 A (from the spectrum)'}


def line(M, a, b):
    return a + b * (M + 21)


def softplus(M, a, s_faint, s_bright, M0, w):
    return a + s_bright * (M - M0) + (s_faint - s_bright) * w * np.logaddexp(0, (M - M0) / w)


def sigmoid(M, a, d, s, M0, w):
    return a + d * np.tanh((M - M0) / w) + s * (M - M0)


MODELS = {'line': (line, ['a', 'b'], 'line: a + b (M+21)'),
          'softplus': (softplus, ['a', 's_faint', 's_bright', 'M0', 'w'], 'softplus: a + s_bright(M-M0) + (s_faint-s_bright) w ln(1+exp((M-M0)/w))'),
          'sigmoid': (sigmoid, ['a', 'd', 's', 'M0', 'w'], 'sigmoid: a + d tanh((M-M0)/w) + s (M-M0)')}


def binned(m, R, nfit):
    """Median of the LINEAR ratio R (negative W3050 kept: dropping them biases faint bins high), logged afterwards.
    FIT = bin holds >= nfit galaxies (the main population; a few hundred odd objects must not set the faint end)."""
    e = np.arange(np.floor(np.quantile(m, 0.002) / DM) * DM, np.quantile(m, 0.998) + DM, DM)
    b = np.digitize(m, e) - 1; rows = []
    for bi in range(len(e) - 1):
        v = R[b == bi]
        if len(v) >= NMIN:
            q16, q50, q84 = np.percentile(v, [16, 50, 84])
            if q50 <= 0:
                continue
            rows.append((0.5 * (e[bi] + e[bi + 1]), np.log10(q50), np.log10(max(q16, 1e-3)), np.log10(max(q84, 1e-3)), len(v),
                         1.2533 * (q84 - q16) / 2 / np.sqrt(len(v)) / (q50 * np.log(10)), len(v) >= nfit))
    return pd.DataFrame(rows, columns=['M', 'med', 'p16', 'p84', 'n', 'err', 'FIT']), e


def uniform(s, mcol, e, nper, rng, strata='program_z'):
    """Uniform sampling across magnitude, for measuring the relation only (Nao 2026-10-31): nper galaxies per M bin, drawn
    equally from every stratum that bin has (strata 'z': Delta z = 0.05 slices; 'program_z': program x z slice, so the bright
    (BGS) and dark (LRG/ELG) programs weigh equally wherever both exist). A stratum short of its share gives all it has; the
    rest is topped up at random from the bin. Returns the sample and {bin centre: number of strata}."""
    b = np.digitize(s[mcol].to_numpy(), e) - 1; zs = np.floor(s.Z.to_numpy() / 0.05).astype(int)
    if strata == 'program_z':                        # PROGRAM: DESI bright / dark; SDSS rows carry their survey (DR8, DR17)
        zs = zs * 8 + pd.factorize(s.PROGRAM)[0]
    take, nz = [], {}
    for bi in range(len(e) - 1):
        ix = np.where(b == bi)[0]
        if not len(ix):
            continue
        sl = np.unique(zs[ix]); nz[round(0.5 * (e[bi] + e[bi + 1]), 4)] = len(sl); q = int(np.ceil(nper / len(sl))); got = []
        for z in sl:
            iz = ix[zs[ix] == z]; got.append(rng.choice(iz, min(q, len(iz)), replace=False))
        got = np.concatenate(got)
        if len(got) > nper:
            got = rng.choice(got, nper, replace=False)
        elif len(got) < nper:
            rest = np.setdiff1d(ix, got); got = np.concatenate([got, rng.choice(rest, min(nper - len(got), len(rest)), replace=False)])
        take.append(got)
    return s.iloc[np.concatenate(take)], nz


BOUNDS = {'softplus': ([-3, -1, -1, None, 0.3], [3, 1, 1, None, 4]), 'sigmoid': ([-3, -2, -0.5, None, 0.3], [3, 2, 0.5, None, 4])}


def fit(B):
    """Lowest-BIC model of MODELS on the binned medians: (name, fn, names, desc, p, C, chi2)."""
    best = None
    for mn, (fn, names, desc) in MODELS.items():
        if len(B) <= len(names) + 1:
            continue
        starts = [[B.med.mean(), 0.0]] if mn == 'line' else \
                 [[B.med.mean(), 0.05, -0.1, M0, 1.0] if mn == 'softplus' else [B.med.mean(), (B.med.iloc[-1] - B.med.iloc[0]) / 2, 0.0, M0, 1.0]
                  for M0 in np.percentile(B.M, [25, 50, 75])]
        for p0 in starts:
            kw = {}
            if mn in BOUNDS:                 # bounded so the parameters stay meaningful (M0 inside the fitted range, 0.3 <= w <= 4 mag)
                lo = [B.M.min() if v is None else v for v in BOUNDS[mn][0]]; hi = [B.M.max() if v is None else v for v in BOUNDS[mn][1]]
                kw['bounds'] = (lo, hi); p0 = list(np.clip(p0, np.array(lo) + 1e-6, np.array(hi) - 1e-6))
            try:
                p, C = curve_fit(fn, B.M, B.med, p0=p0, sigma=B.err, absolute_sigma=True, maxfev=50000, **kw)
            except RuntimeError:
                continue
            chi2 = float((((B.med - fn(B.M, *p)) / B.err)**2).sum()); bic = chi2 + len(p) * np.log(len(B))
            if best is None or bic < best[0]:
                best = (bic, mn, fn, names, desc, p, C, chi2)
    return best[1:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--passc', default=os.path.expanduser('~/data/sdssdesi_passC_v2_resp'))
    ap.add_argument('--out', default=os.path.expanduser('~/data/sdssdesi_passC_v2_resp/qa_2win')); ap.add_argument('--zmax', type=float, default=9.0)
    ap.add_argument('--sdss', default='', help='also use the SDSS records in <dir>/normtab (e.g. ~/data/sdssdesi_passG_sdss_v2_2win)')
    ap.add_argument('--sdss-mode', default='fill', choices=['fill', 'all'],
                    help="fill: SDSS only in (M bin, z slice) cells where DESI has < NFILL galaxies (Nao: just the region DESI misses)")
    ap.add_argument('--strata', default='program_z', choices=['z', 'program_z'])
    ap.add_argument('--uniform', type=int, default=0, help='galaxies drawn per M bin, equally over z slices (0 = all galaxies)')
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(20261031); tag = (f'_uniform{a.uniform}_{a.strata}' if a.uniform else '') + (f'_sdss{a.sdss_mode}' if a.sdss else '')
    t = pd.concat([pd.read_parquet(f, columns=COLS) for f in sorted(glob.glob(os.path.join(a.passc, 'normtab', 'norm_shard*.parquet')))],
                  ignore_index=True); t['SRC'] = 'DESI'
    if a.sdss:                                       # SDSS v2 per-galaxy records (16_norm_table_sdss.py); no program -> the survey
        for f in sorted(glob.glob(os.path.join(a.sdss, 'normtab', 'norm_dr*.parquet'))):
            x = pd.read_parquet(f, columns=[c for c in COLS if c != 'PROGRAM'] + ['SURVEY'])
            t = pd.concat([t, x.assign(PROGRAM=x.SURVEY, SRC=x.SURVEY).drop(columns='SURVEY')], ignore_index=True)
    t = t[(t.COV3050 >= COVER) & (t.COV5450 >= COVER) & (t.W5450 > 0) & (t.Z < a.zmax)].copy()
    t['R'] = t.W3050 / t.W5450
    with np.errstate(all='ignore'):
        t['logR'] = np.log10(t.R)
    fl = t.W5450.to_numpy(float) / (t.NORM_SCALE * t.RESCALE).to_numpy(float)
    t['MSPEC'] = -2.5 * np.log10(fl * 5450.0**2 / C_AA) - 48.60 - pb.COSMO.distmod(t.Z.to_numpy(float)).value
    zr = [round(float(t.Z.min()), 3), round(float(t.Z.max()), 3)]
    res = {}; look = {}
    for mcol, lab in AXES.items():
        res[mcol] = {}
        fig, ax = plt.subplots(2, 5, figsize=(27, 9.5), gridspec_kw=dict(height_ratios=[3, 1]), sharex='col')
        for db in range(5):
            s = t[(t.DB == db) & np.isfinite(t[mcol])]
            if a.sdss and a.sdss_mode == 'fill':  # SDSS only where DESI has no population: (M bin, z slice) cells with < NFILL DESI
                d_ = s[s.SRC == 'DESI']; ee = np.arange(np.floor(d_[mcol].quantile(0.002) / DM) * DM - 3, d_[mcol].quantile(0.998) + DM + 3, DM)
                key = lambda x: (np.digitize(x[mcol], ee) * 1000 + np.floor(x.Z / 0.05)).to_numpy()
                kd, nd = np.unique(key(d_), return_counts=True); full = set(kd[nd >= NFILL])
                sd = s[s.SRC != 'DESI']; sd = sd[~np.isin(key(sd), list(full))]
                s = pd.concat([d_, sd], ignore_index=True)
                if mcol == 'MABS_G':
                    print(f'DB{db}: SDSS fill {len(sd):,} galaxies ({sd.SRC.value_counts().to_dict()}), MABS_G 5/50/95 % '
                          f'{np.percentile(sd[mcol], [5, 50, 95]).round(2).tolist() if len(sd) else "-"}, z median {sd.Z.median():.2f}')
            B, e = binned(s[mcol].to_numpy(), s.R.to_numpy(), max(NMIN, FITFRAC * len(s)))
            nz = {}
            if a.sdss and mcol == 'MABS_G':  # each survey/program median minus the pooled one, same DB and M bin (n >= 200)
                mb = np.digitize(s[mcol], e) - 1; g = s.assign(mb=mb).groupby(['mb', 'PROGRAM']).R.agg(['median', 'size'])
                g = g[g['size'] >= 200]; pool = s.assign(mb=mb).groupby('mb').R.median()
                g['dev'] = np.log10(g['median'].clip(1e-3)) - np.log10(pool.reindex(g.index.get_level_values(0)).clip(1e-3).to_numpy())
                print(f'DB{db} survey - pooled (dex, median over M bins): ' + ', '.join(f'{k} {v.dev.median():+.3f} ({len(v)} bins, {int(v["size"].sum()):,})'
                                                                               for k, v in g.groupby(level=1)))
            if a.uniform:                    # the same bins, re-measured on a uniform draw; bins left out above stay out
                keep = set(B.M[B.FIT].round(4)); s, nz = uniform(s, mcol, e, a.uniform, rng, a.strata)
                B, _ = binned(s[mcol].to_numpy(), s.R.to_numpy(), NMIN); B['FIT'] = B.FIT & B.M.round(4).isin(keep)
            hold = HOLD.get(db) if mcol == 'MABS_G' else None
            if hold is not None:             # red types: fit only brighter than the targeting gap (faint side = faint blue interlopers)
                B['FIT'] = B.FIT & (B.M <= hold)
            mn, fn, names, desc, p, C, chi2 = fit(B[B.FIT].reset_index(drop=True))
            B['fit'] = fn(B.M, *p); B['resid'] = B.med - B.fit; F = B[B.FIT]
            r = res[mcol][db] = dict(Dn4000=DBLAB[db], model=desc, axis=lab, n_galaxies=int(len(s)), z_range=zr,
                                     params=dict(zip(names, map(float, p))), errors=dict(zip(names, map(float, np.sqrt(np.diag(C))))),
                                     cov=C.tolist(), chi2=chi2, dof=len(F) - len(p),
                                     rms_resid_dex=float(np.sqrt(np.average(F.resid**2, weights=F.n))),
                                     max_abs_resid_dex=float(F.resid.abs().max()), M_range=[float(F.M.min()), float(F.M.max())],
                                     logR_median=float(s.logR.median()), bins=B.round(5).to_dict(orient='list'))
            a0 = ax[0, db]
            a0.hexbin(s[mcol], s.logR, gridsize=(70, 60), extent=(e[0], e[-1], -1.6, 1.0), bins='log', cmap='Greys', mincnt=1)
            a0.fill_between(B.M, B.p16, B.p84, color='C0', alpha=0.15, label='16-84 %')
            a0.errorbar(F.M, F.med, F.err, fmt='o', ms=3, color='C0', label=('DESI+SDSS' if a.sdss else 'DESI') + f' median per {DM} mag ({len(s):,}' + (f', {a.uniform} per bin over {a.strata})' if a.uniform else ')'))
            X = B[~B.FIT]
            if len(X):
                a0.plot(X.M, X.med, 'o', ms=4, mfc='none', color='C0', label=f'bins < {FITFRAC:.1%} of the type, not fit')
            mm = np.linspace(F.M.min(), F.M.max(), 300); a0.plot(mm, fn(mm, *p), 'C3-', lw=2, label=f'fit ({mn})')
            if True:                         # lookup (both axes): the fit inside the fitted range, held at its end values outside it
                lm = LOOK_M if mcol == 'MABS_G' else LOOK_M_SPEC
                lo_, hi_ = float(F.M.min()), float(F.M.max()); g = np.round(np.arange(lm[0], lm[1] + 1e-9, DM), 3)
                v = fn(np.clip(g, lo_, hi_), *p); src = np.where(g > hi_, 'held_faint', np.where(g < lo_, 'held_bright', 'fit'))
                if db in PLATEAU and mcol == 'MABS_G':   # held at the bright plateau: median of the fitted bins brighter than PLATEAU[db], from the
                    pv = float(F.med[F.M <= PLATEAU[db]].median())     # faintest grid point where the fit is within 0.01 dex of it
                    ok_ = np.where((g <= hi_) & (np.abs(fn(np.clip(g, lo_, hi_), *p) - pv) <= 0.01))[0]
                    mc = g[ok_.max()] if len(ok_) else PLATEAU[db]
                    v = np.where(g > mc, pv, v); src = np.where(g > mc, 'held_plateau', src)
                bb = B.set_index(B.M.round(3))
                look.setdefault(mcol, []).append(pd.DataFrame(dict(DB=db, Dn4000=DBLAB[db], **{'MABS_G' if mcol == 'MABS_G' else 'M_AB5450': g}, logR=np.round(v, 4), R=np.round(10**v, 4), SOURCE=src,
                                              spread_16_84_dex=[round(float((bb.p84 - bb.p16).get(x, np.nan) / 2), 4) for x in g],
                                              n=[int(bb.n.get(x, 0)) for x in g], M_fit_lo=lo_, M_fit_hi=hi_, model=mn)))
                a0.plot(g, v, 'C1--', lw=1.5, label='lookup (held parts flat)')
                if hold is not None:
                    a0.axvline(hold, color='C1', lw=0.8, ls=':'); a0.text(hold, 0.9, f' hold {hold}', color='C1', fontsize=7)
            a0.set_title(f'DB{db}  {DBLAB[db]}\n{desc}\n' + ', '.join(f'{n}={v:.3f}' for n, v in zip(names, p)), fontsize=8)
            a0.set_ylim(-1.6, 1.0); a0.legend(fontsize=7, loc='upper left'); a0.invert_xaxis()
            sec = a0.secondary_yaxis('right', functions=(lambda y: -2.5 * y, lambda m: -m / 2.5))
            if db == 4:
                sec.set_ylabel('m(3050) - m(5450), f_lam mag')
            a1 = ax[1, db]
            a1.axhline(0, color='k', lw=0.5); a1.errorbar(F.M, F.resid, F.err, fmt='o', ms=3, color='C0'); a1.plot(X.M, X.resid, 'o', ms=4, mfc='none', color='C0')
            a1.set_ylim(-0.12, 0.12); a1.set_xlabel(lab)
            for x_, y_ in zip(B.M, B.resid):
                if round(x_, 4) in nz:
                    a1.text(x_, -0.115, str(nz[round(x_, 4)]), ha='center', fontsize=6, color='0.4')
            a1.set_title((f'number = {a.strata} strata per bin, ' if nz else '') + f'rms {r["rms_resid_dex"]:.4f} dex, max {r["max_abs_resid_dex"]:.4f} dex', fontsize=9)
        ax[0, 0].set_ylabel('log R = log F(3000-3100)/F(5400-5500)'); ax[1, 0].set_ylabel('median - fit (dex)')
        fig.suptitle(f'R = F(3000-3100)/F(5400-5500) vs {lab}, per Dn4000 bin, ' + ('DESI+SDSS' if a.sdss else 'DESI') + f' z {zr[0]}-{zr[1]}' + (f', uniform sampling {a.uniform} per 0.25 mag over {a.strata}' if a.uniform else ''), fontsize=13)
        fig.tight_layout(); fig.savefig(os.path.join(a.out, f'gal_ratio_fit_{mcol}{tag}.png'), dpi=80); plt.close(fig)
    json.dump(res, open(os.path.join(a.out, f'gal_ratio_fit{tag}.json'), 'w'), indent=1)
    for mcol, lk in look.items():                    # MABS_G: gal_ratio_lookup<tag>.csv (as before); M_AB(5450): ..._MSPEC.csv
        L = pd.concat(lk, ignore_index=True); fn_ = f'gal_ratio_lookup{tag}' + ('' if mcol == 'MABS_G' else '_MSPEC') + '.csv'
        L.to_csv(os.path.join(a.out, fn_), index=False); xc = 'MABS_G' if mcol == 'MABS_G' else 'M_AB5450'
        print(f'lookup: {fn_}'); print(L.pivot(index=xc, columns='DB', values='logR').iloc[::4].round(3).to_string())
    for mcol, d in res.items():
        for db, r in d.items():
            print(mcol, f'DB{db}', r['model'].split(':')[0], r['n_galaxies'], {n: round(v, 3) for n, v in r['params'].items()},
                  f'rms {r["rms_resid_dex"]:.4f} max {r["max_abs_resid_dex"]:.4f}', r['M_range'])


if __name__ == '__main__':
    main()
