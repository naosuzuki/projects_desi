"""Galaxy release v3 (Nao 2026-10-31): two versions of the stacking grid that differ only in the luminosity axis, built in ONE
re-read of the original spectra (DESI coadds, SDSS spPlates):
  version G     MB = floor(MABS_G), the photometric rest-g absolute magnitude (whole galaxy; K from the spectrum)
  version 5450  MB = floor(M_AB(5450)), the absolute AB magnitude of the galaxy's own spectrum at 5400-5500 A (fiber light)
Both: two-window normalization (normwin v2: 5400-5500 below z = 0.30, 3000-3100 above), the same z, color, EW and Dn4000 axes.

Where the band is not covered:
  MABS_G (rest g leaves the spectrum above z ~ 0.75): MABS_U + <G-U>, above z ~ 1.45 MABS_NUV + <U-NUV> + <G-U>, with the
          median colors of the galaxy's EFFECTIVE type DB_EFF (measured on the galaxies that have both, DESI v2)
  M_AB(5450) (above z ~ 0.78, or < 40 % of the window usable): from the flux at 3000-3100 A and the lookup ratio
          R = F(3000-3100)/F(5400-5500) on the M_AB(5450) axis (18_gal_ratio_fit.py, _MSPEC), solved self-consistently
DB_EFF: Dn4000 bin, re-typed by line EW above z = 0.8 (DB >= 1; Dn4000 at the red end of the spectrum) and for DB >= 2 galaxies
fainter than MABS_G = -19.5 (the targeting gap's blue interlopers): Halpha EW below z = 0.46, [OII] above, boundaries halfway in
log EW between the per-DB median EWs where Dn4000 is reliable (19_chain_stacks.py, DESI).
"""
import os
import numpy as np
import pandas as pd

C_AA = 2.99792458e18
GU = {0: -0.691, 1: -1.058, 2: -1.449, 3: -1.654, 4: -1.759}          # median MABS_G - MABS_U per DB (DESI v2)
UN = {0: -0.425, 1: -0.544, 2: -1.399, 3: -1.950, 4: -1.828}          # median MABS_U - MABS_NUV per DB
EW_EDGES_RETYPE = {'Halpha': np.log10([25.6, 9.6, 2.2, 0.9]), '[OII]': np.log10([28.3, 9.2, 3.4, 2.0])}
ZTYPE, MFAINT = 0.8, -19.5
LOOKUP_5450 = os.path.expanduser('~/data/sdssdesi_passC_v2_resp/qa_2win/gal_ratio_lookup_uniform2000_program_z_sdssfill_MSPEC.csv')
M5450_MIN, M5450_MAX = -24, -14     # version-5450 MB index range (fiber light runs 0.3-0.8 mag fainter than MABS_G)
W5450, W3050 = (5400, 5500), (3000, 3100)
MIN_GOOD = 0.4


def db_eff(z, db, ew, mabs_g):
    """Vectorized effective type (see module doc)."""
    z = np.asarray(z, float); db = np.asarray(db, int); ew = np.asarray(ew, float); mg = np.asarray(mabs_g, float)
    retype = ((z >= ZTYPE) & (db >= 1)) | ((db >= 2) & np.isfinite(mg) & (mg > MFAINT))
    retype &= np.isfinite(ew)
    le = np.log10(np.clip(np.nan_to_num(ew, nan=0.5), 0.5, None))
    by_ha = np.sum(le[:, None] < EW_EDGES_RETYPE['Halpha'][None, :], axis=1)
    by_oii = np.sum(le[:, None] < EW_EDGES_RETYPE['[OII]'][None, :], axis=1)
    t = np.where(z < 0.46, by_ha, by_oii)
    return np.where(retype, np.minimum(db, t), db)


def mabs_g_full(mg, mu, mnuv, dbe):
    """MABS_G where measured, else from MABS_U / MABS_NUV with the colors of the effective type; and the source (0 g, 1 U, 2 NUV)."""
    mg, mu, mnuv = (np.asarray(x, float) for x in (mg, mu, mnuv)); dbe = np.asarray(dbe, int)
    gu = np.array([GU[int(d)] for d in range(5)])[dbe]; un = np.array([UN[int(d)] for d in range(5)])[dbe]
    out = np.where(np.isfinite(mg), mg, np.where(np.isfinite(mu), mu + gu, mnuv + un + gu))
    src = np.where(np.isfinite(mg), 0, np.where(np.isfinite(mu), 1, 2))
    return out, src


def load_lookup_5450(path=LOOKUP_5450):
    L = pd.read_csv(path).sort_values('M_AB5450')
    return {int(db): (g.M_AB5450.to_numpy(), g.logR.to_numpy()) for db, g in L.groupby('DB')}


def mab(flam, lam, dm):
    """Absolute AB magnitude from the rest-frame F_lambda [erg/s/cm2/A] at lam and the distance modulus."""
    with np.errstate(all='ignore'):
        return np.where(flam > 0, -2.5 * np.log10(flam * lam**2 / C_AA) - 48.60 - dm, np.nan)


def m5450(f5450, f3050, dbe, dm, look, niter=6):
    """M_AB(5450) from the rest-frame flux in 5400-5500 A where available, else self-consistently from 3000-3100 A:
    F(5450) = F(3050) / R(DB_EFF, M_AB(5450)). Returns (M, source 0 = 5450 measured, 1 = from 3050, -1 = none)."""
    f5450, f3050, dm = (np.asarray(x, float) for x in (f5450, f3050, dm)); dbe = np.asarray(dbe, int)
    m = mab(f5450, 5450.0, dm); src = np.where(np.isfinite(m), 0, -1)
    need = ~np.isfinite(m) & (f3050 > 0)
    if need.any():
        m3 = mab(f3050[need], 3050.0, dm[need]); d = dbe[need]
        mm = m3 - 5 * np.log10(5450.0 / 3050.0) + 2.5 * (-0.5)            # start: R ~ 10^-0.5
        for _ in range(niter):
            lr = np.array([np.interp(x, *look[int(k)]) for x, k in zip(mm, d)])
            mm = mab(f3050[need] / 10**lr, 5450.0, dm[need])
        m[need] = mm; src[need] = 1
    return m, src


def window_mean(rwave, rflux, rivar, win, min_good=MIN_GOOD):
    """Mean F_lambda over the good pixels of `win` per row; NaN if < min_good of the window is usable."""
    sel = (rwave >= win[0]) & (rwave <= win[1]); f = np.atleast_2d(rflux)[:, sel]; g = np.atleast_2d(rivar)[:, sel] > 0
    n = g.sum(axis=1)
    with np.errstate(all='ignore'):
        mean = np.where(n > 0, np.where(g, f, 0.0).sum(axis=1) / np.where(n > 0, n, 1), np.nan)
    return np.where(n >= min_good * sel.sum(), mean, np.nan)


def window_mean_cov(rwave, rflux, rivar, win, min_good=MIN_GOOD):
    """(mean F_lambda over the good pixels of `win`, NaN if < min_good usable; good-pixel fraction) per row."""
    sel = (rwave >= win[0]) & (rwave <= win[1]); f = np.atleast_2d(rflux)[:, sel]; g = np.atleast_2d(rivar)[:, sel] > 0
    n = g.sum(axis=1); cov = n / sel.sum()
    with np.errstate(all='ignore'):
        mean = np.where(n > 0, np.where(g, f, 0.0).sum(axis=1) / np.where(n > 0, n, 1), np.nan)
    return np.where(cov >= min_good, mean, np.nan), cov


CELLS = {'G': os.path.expanduser('~/data/sdssdesi_products/galaxy_cell_table_v3G.csv'),
         '5450': os.path.expanduser('~/data/sdssdesi_products/galaxy_cell_table_v35450.csv')}
MRANGE = {'G': (-24, -16), '5450': (M5450_MIN, M5450_MAX)}
ACC_MIN_N = 200


class V3:
    """Both versions' cell maps and accumulators for one worker (08_passC_stack.py --v3, 11_passG_sdss.py --v3)."""

    def __init__(self, accpath, writer):
        import gridcells as gc
        self.gc = gc; self.look = load_lookup_5450(); self.accpath = accpath; self.cmap = {}; self.acc = {}; self.done = set()
        self.meta = {}
        for ver, p in CELLS.items():
            t = pd.read_csv(p); t = t[t.N >= ACC_MIN_N].reset_index(drop=True)
            codes = gc.cell_code({k: t[k].to_numpy() for k in gc.KEY}); self.cmap[ver] = {int(c): i for i, c in enumerate(codes)}
            nc = len(t); ap = accpath(ver)
            if os.path.exists(ap):
                z = np.load(ap, allow_pickle=False); self.acc[ver] = {k: z[k] for k in ('wf', 'w', 'wf2', 'n', 'ncell')}
                self.done |= set(z['done'].tolist())
            else:
                self.acc[ver] = dict(wf=np.zeros((nc, NREST_), np.float32), w=np.zeros((nc, NREST_), np.float32),
                                     wf2=np.zeros((nc, NREST_), np.float32), n=np.zeros((nc, NREST_), np.int32), ncell=np.zeros(nc, np.int64))
            self.meta[ver] = dict(cells=p, acc_min_n=ACC_MIN_N, nrest=NREST_, rid_min=RID_MIN_, dlog=1e-4, lines='vacuum', ncells=nc,
                                  writer=writer, version=ver, maxis=('MABS_G_FULL' if ver == 'G' else 'M_AB5450'),
                                  windows='5400-5500 (z<0.30) / 3000-3100', min_good=MIN_GOOD)

    def save(self, save_acc):
        for ver in self.acc:
            save_acc(self.accpath(ver), self.acc[ver], self.done, self.meta[ver])

    def process(self, df, rwave, rflux, rivar, z, dm, normwin):
        """df: measured galaxies (MABS_G/U/NUV, COLOR, EW, DN4000); rflux/rivar: calibrated rest-frame F_lambda (cgs).
        Adds the v3 columns to df, accumulates both versions; returns (n stacked G, n stacked 5450)."""
        gc = self.gc
        z = np.asarray(z, float); col = df.COLOR.to_numpy(float); ew = df.EW.to_numpy(float); dn = df.DN4000.to_numpy(float)
        db = np.searchsorted(gc.DN_EDGES, np.nan_to_num(dn, nan=-9.0)); eb = np.searchsorted(gc.EW_EDGES, np.nan_to_num(ew, nan=-99.0))
        zb = np.floor(np.clip(z, 0, None) / gc.DZ).astype(int); cb = np.floor(np.nan_to_num(col, nan=-99) / 0.25).astype(int)
        mg0, _ = mabs_g_full(df.MABS_G, df.MABS_U, df.MABS_NUV, db)
        dbe = db_eff(z, db, ew, mg0); mgf, mgs = mabs_g_full(df.MABS_G, df.MABS_U, df.MABS_NUV, dbe)
        f5, c5 = window_mean_cov(rwave, rflux, rivar, W5450); f3, c3 = window_mean_cov(rwave, rflux, rivar, W3050)
        m5, m5s = m5450(f5, f3, dbe, dm, self.look)
        s = np.full(len(df), np.nan); nok = np.zeros(len(df), bool)
        for b in np.unique(zb):
            r = zb == b
            if not (0 <= b < round(gc.ZMAX / gc.DZ)):
                continue
            win = normwin.slice_window(b * gc.DZ, (b + 1) * gc.DZ, 'GALAXY')
            s[r], nok[r] = normwin.normalize(rwave, rflux[r], rivar[r], win, min_good=MIN_GOOD)
        for k, v in dict(ZB=zb, CB=cb, EB=eb, DB=db, DB_EFF=dbe, MABS_G_FULL=mgf, MG_SRC=mgs, M_AB5450=m5, M5_SRC=m5s,
                         F5450=f5, COV5450=c5, F3050=f3, COV3050=c3, NORM_SCALE=s, NORM_OK=nok).items():
            df[k] = v
        nacc = {}
        base = np.isfinite(col) & (z > 0) & (z < gc.ZMAX) & (cb >= gc.C_MIN) & (cb <= gc.C_MAX)
        for ver, M in (('G', mgf), ('5450', m5)):
            lo, hi = MRANGE[ver]; mb = np.floor(np.nan_to_num(M, nan=99)).astype(int)
            ok = base & np.isfinite(M) & (mb >= lo) & (mb <= hi)
            code = gc.cell_code(dict(ZB=zb, MB=mb, CB=cb, EB=eb, DB=db))
            idx = np.array([self.cmap[ver].get(int(c), -1) if o else -1 for c, o in zip(code, ok)], dtype=np.int32)
            df[f'MB_{ver}'] = np.where(ok, mb, -99); df[f'ACC_{ver}'] = idx
            use = nok & (idx >= 0); nacc[ver] = int(use.sum())
            if use.any():
                fn = rflux[use] * s[use, None]; iv = rivar[use] / s[use, None]**2; ci = idx[use]; A = self.acc[ver]
                for c in np.unique(ci):
                    m = ci == c
                    A['wf'][c] += (iv[m] * fn[m]).sum(axis=0); A['w'][c] += iv[m].sum(axis=0)
                    A['wf2'][c] += (iv[m] * fn[m]**2).sum(axis=0); A['n'][c] += (iv[m] > 0).sum(axis=0); A['ncell'][c] += m.sum()
        return nacc['G'], nacc['5450']


NREST_, RID_MIN_ = None, None


def init_grid(nrest, rid_min):
    global NREST_, RID_MIN_
    NREST_, RID_MIN_ = nrest, rid_min
