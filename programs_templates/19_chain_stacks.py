#!/usr/bin/env python3
"""Chain the galaxy stacks to one normalization point (docs/HANDOVER_normalization.md step 4; Nao 2026-10-31).

Galaxies are normalized at two points (normwin v2): 5400-5500 A below z = 0.30, 3000-3100 A above. Every stack is put on
ONE scale, unit mean F_lambda in 5400-5500 A:
  1. renormalize the stack in its own window (stack mean over 3000-3100 or 5400-5500 = 1; for stacks built with the v1
     five-window normalization this is the stack-level approximation of the per-galaxy renormalization),
  2. stacks of the 3000-3100 slices: x R(DB, MABS_G) from the lookup (18_gal_ratio_fit.py, gal_ratio_lookup*.csv), R =
     F(3000-3100)/F(5400-5500), interpolated in the cell's median MABS_G (from the per-galaxy records, normtab); a stack
     with F(3050) = 1 then has F(5450) = 1/R, so multiplying by R gives F(5450) = 1.
Above z = 0.8 Dn4000 sits at the red end of the spectrum and scatters star-forming galaxies into DB2-4 ([OII] EW 50-80 A in
the z > 1.2 "passive" stacks): there a DB1-4 stack is re-typed by its [OII] EW (DB_EFF; boundaries halfway in log between
the per-DB EW of the z 0.46-0.8 stacks) and DB_EFF picks the lookup curve. DB2-4 stacks fainter than MABS_G = -19.5 (the
targeting gap; the faint side holds blue interlopers) are re-typed the same way at every z, by Halpha EW below z = 0.46.
Where the stack also covers 5400-5500 (NSPEC >= 0.8 N_STACKED), the direct ratio R_DIRECT is kept beside R_LOOKUP as a check.
Output <in>/<out>/stacks_z<zb>.fits: same HDUs as 09_merge_stacks.py (FLUX, IVAR, SCATTER scaled; NSPEC, WAVE), CELLS with
NORM_WIN, MABS_G_MED (MABS_G_SRC 0 = g, 1 = from U, 2 = from NUV, per-DB colours), R_LOOKUP, R_DIRECT, CHAIN (the factor applied after step 1), RENORM (step 1 factor).
Usage: python3 19_chain_stacks.py --stacks DIR --normtab GLOB [--lookup CSV] [--out stacks_chain]
"""
import os, glob, argparse
import numpy as np
import pandas as pd
import fitsio

W5450, W3050 = (5400, 5500), (3000, 3100)
SWITCH_Z, DZ, COVER = 0.30, 0.05, 0.8
KEY = ['ZB', 'MB', 'CB', 'EB', 'DB']
ZTYPE = 0.8; MFAINT = -19.5
LOOKUP = os.path.expanduser('~/data/sdssdesi_passC_v2_resp/qa_2win/gal_ratio_lookup_uniform2000_program_z_sdssfill.csv')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stacks', required=True); ap.add_argument('--normtab', required=True, help='glob of per-galaxy records')
    ap.add_argument('--lookup', default=LOOKUP); ap.add_argument('--out', default=None)
    a = ap.parse_args()
    out = a.out or os.path.join(os.path.dirname(a.stacks.rstrip('/')), os.path.basename(a.stacks.rstrip('/')) + '_chain')
    os.makedirs(out, exist_ok=True)
    nt = pd.concat([pd.read_parquet(f, columns=KEY + ['MABS_G', 'MABS_U', 'MABS_NUV']) for f in sorted(glob.glob(a.normtab))], ignore_index=True)
    # rest g leaves the spectrum above z ~ 0.8: MABS_G = MABS_U + <G-U>_DB (or MABS_NUV + <U-NUV>_DB + <G-U>_DB above z ~ 1.45),
    # the colours being the per-DB medians over the galaxies that have both
    gu = (nt.MABS_G - nt.MABS_U).groupby(nt.DB).median(); un = (nt.MABS_U - nt.MABS_NUV).groupby(nt.DB).median()
    est = np.where(nt.MABS_G.notna(), nt.MABS_G, np.where(nt.MABS_U.notna(), nt.MABS_U + nt.DB.map(gu),
                                                         nt.MABS_NUV + nt.DB.map(un) + nt.DB.map(gu)))
    src = np.where(nt.MABS_G.notna(), 0, np.where(nt.MABS_U.notna(), 1, 2))
    mg = nt.assign(MG=est, SRC=src).groupby(KEY).agg(MABS_G_MED=('MG', 'median'), MABS_G_SRC=('SRC', 'median')).reset_index()
    print('colours used to reach MABS_G: G-U ' + str(gu.round(3).to_dict()) + ', U-NUV ' + str(un.round(3).to_dict()))
    L = pd.read_csv(a.lookup); look = {db: (g.MABS_G.to_numpy(), g.logR.to_numpy()) for db, g in L.sort_values('MABS_G').groupby('DB')}
    # re-typing above z = ZTYPE (Dn4000 at the red end of the spectrum, unreliable): a DB2-4 stack is given the type its [OII] EW
    # says, with the per-DB median [OII] EW of the z 0.46-0.8 stacks (Dn4000 still reliable, EW already [OII]) as the reference
    # and below it, the faint red stacks (DB2-4 fainter than MFAINT, the DESI/SDSS targeting gap: blue interlopers) are re-typed
    # the same way by their own line: Halpha EW below z = 0.46 (reference: stacks z 0.05-0.46 brighter than MFAINT), [OII] above
    cal = []
    for fn in sorted(glob.glob(os.path.join(a.stacks, 'stacks_z*.fits'))):
        x = fitsio.read(fn, ext='CELLS', columns=KEY + ['Z_MED', 'EW_MED', 'N_STACKED']); x = pd.DataFrame(x.astype(x.dtype.newbyteorder('=')))
        cal.append(x.merge(mg, on=KEY, how='left'))
    cal = pd.concat(cal, ignore_index=True)
    edges = {}
    for line, sel in (('Halpha', (cal.Z_MED >= 0.05) & (cal.Z_MED < 0.46) & (cal.MABS_G_MED <= MFAINT)),
                      ('[OII]', (cal.Z_MED >= 0.46) & (cal.Z_MED < ZTYPE))):
        x = cal[sel]; lew = np.log10(x.EW_MED.clip(0.5))
        ref = {db: float(np.average(g, weights=x.N_STACKED[g.index])) for db, g in lew.groupby(x.DB)}
        edges[line] = np.array([(ref[d] + ref[d + 1]) / 2 for d in range(4)])     # log EW boundaries DB0|...|DB4 (EW falls with DB)
        print(f'{line} EW reference per DB: ' + ', '.join(f'DB{d} {10**v:.1f} A' for d, v in ref.items())
              + '; boundaries ' + ', '.join(f'{10**e:.1f}' for e in edges[line]))

    def db_eff(z, ew, db, mg_):
        hiz = z >= ZTYPE and db >= 1; faint = db >= 2 and np.isfinite(mg_) and mg_ > MFAINT
        if not ((hiz or faint) and np.isfinite(ew)):
            return db
        return min(db, int(np.sum(np.log10(max(ew, 0.5)) < edges['Halpha' if z < 0.46 else '[OII]'])))
    rows = []
    for fn in sorted(glob.glob(os.path.join(a.stacks, 'stacks_z*.fits'))):
        F = fitsio.FITS(fn); f = F['FLUX'].read(); iv = F['IVAR'].read(); sc = F['SCATTER'].read(); ns = F['NSPEC'].read()
        w = F['WAVE'].read(); c = pd.DataFrame(F['CELLS'].read().astype(F['CELLS'].read().dtype.newbyteorder('='))) ; F.close()
        zb = int(c.ZB.iloc[0]); upper = (zb + 1) * DZ > SWITCH_Z + 1e-9
        s5 = (w >= W5450[0]) & (w <= W5450[1]); s3 = (w >= W3050[0]) & (w <= W3050[1]); own = s3 if upper else s5
        c = c.merge(mg, on=KEY, how='left')
        N = np.maximum(c.N_STACKED.to_numpy(), 1)
        with np.errstate(all='ignore'):
            m_own = np.nanmean(np.where(ns[:, own] > 0, f[:, own], np.nan), axis=1)
            m5 = np.nanmean(np.where(ns[:, s5] > 0, f[:, s5], np.nan), axis=1); m3 = np.nanmean(np.where(ns[:, s3] > 0, f[:, s3], np.nan), axis=1)
        cov5 = ns[:, s5].min(axis=1) / N; cov3 = ns[:, s3].min(axis=1) / N
        renorm = 1.0 / m_own
        rdir = np.where((cov5 >= COVER) & (cov3 >= COVER) & (m5 > 0), m3 / m5, np.nan)
        c['DB_EFF'] = [db_eff(z, e, int(d), m) for z, e, d, m in zip(c.Z_MED, c.EW_MED, c.DB, c.MABS_G_MED)]
        rlook = np.array([10**np.interp(m, *look[int(d)]) if np.isfinite(m) else np.nan for m, d in zip(c.MABS_G_MED, c.DB_EFF)])
        chain = rlook if upper else np.ones(len(c))
        k = renorm * chain
        f2 = f * k[:, None]; iv2 = iv / k[:, None]**2; sc2 = sc * k[:, None]
        c['NORM_WIN'] = '3000-3100' if upper else '5400-5500'; c['RENORM'] = renorm; c['R_LOOKUP'] = rlook; c['R_DIRECT'] = rdir; c['CHAIN'] = chain
        with fitsio.FITS(os.path.join(out, os.path.basename(fn)) + '.tmp', 'rw', clobber=True) as fo:
            fo.write(f2.astype(np.float32), extname='FLUX'); fo.write(iv2.astype(np.float32), extname='IVAR')
            fo.write(sc2.astype(np.float32), extname='SCATTER'); fo.write(ns, extname='NSPEC'); fo.write(w, extname='WAVE')
            fo.write(c.to_records(index=False), extname='CELLS')
        os.replace(os.path.join(out, os.path.basename(fn)) + '.tmp', os.path.join(out, os.path.basename(fn)))
        rows.append(c.assign(FILE=os.path.basename(fn)))
    C = pd.concat(rows, ignore_index=True); C.to_csv(os.path.join(out, 'cell_chain.csv'), index=False)
    rt = C[C.DB_EFF != C.DB]
    print(f'  re-typed (z >= {ZTYPE} DB1-4, or DB2-4 fainter than {MFAINT}, by line EW): {len(rt)} stacks, {int(rt.N_STACKED.sum()):,} galaxies; '
          + ', '.join(f'DB{a_}->DB{b_} {n_}' for (a_, b_), n_ in rt.groupby(['DB', 'DB_EFF']).size().items()))
    up = C[C.NORM_WIN == '3000-3100']; both = up[(up.R_DIRECT > 0) & np.isfinite(up.R_LOOKUP)]
    d = np.log10(both.R_DIRECT / both.R_LOOKUP)
    print(f'{out}: {len(C)} stacks, {len(up)} chained from 3000-3100 (missing MABS_G {up.MABS_G_MED.isna().sum()}); '
          f'renorm median {np.nanmedian(C.RENORM):.4f} (16-84 % {np.nanpercentile(C.RENORM, 16):.3f}-{np.nanpercentile(C.RENORM, 84):.3f})')
    print(f'  check where 5400-5500 is covered too ({len(both)} stacks, z {both.Z_MED.min():.2f}-{both.Z_MED.max():.2f}): '
          f'log R_direct - log R_lookup median {np.median(d):+.3f}, 16-84 % {np.percentile(d, 16):+.3f}..{np.percentile(d, 84):+.3f} dex; by DB '
          + ', '.join(f'DB{k} {v.median():+.3f}' for k, v in d.groupby(both.DB)))


if __name__ == '__main__':
    main()
