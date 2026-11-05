#!/usr/bin/env python3
"""Galaxy normalization link from individual galaxies (docs/HANDOVER_normalization.md step 3), per-object medians, laid out as
programs_quasars/10_ratio_vs_mag.py.

Galaxies are normalized at two points (normwin v2, Nao 2026-10-31): 5400-5500 A below z = 0.30, 3000-3100 A above. The one
link is R = F(3000-3100)/F(5400-5500), measured on every galaxy covering both windows (COV >= 0.8 each; z ~ 0.2-0.8), from
its own spectrum (the ratio does not depend on where the spectrum was normalized).
Inputs: DESI <passC>/normtab/norm_shard*.parquet (12_restack_response.py), SDSS <passG>/normtab/norm_dr{8,17}.parquet
(16_norm_table_sdss.py). SNR_SDSSPIX is taken from the record, or joined from the pass tables for older records.
Magnitude axes (both shown):
  MABS_G   photometric rest g absolute magnitude (whole galaxy); its noise is independent of the window fluxes -> the lookup axis
  M_AB(5450) from the spectrum = -2.5 log10(F_lambda w^2/c) - 48.60 - DM(z), F_lambda = W5450/(NORM_SCALE x RESCALE): shares the
           noise of the denominator, so faint noisy spectra scatter to high R (kept to show that bias)
Cut: SNR_SDSSPIX >= --snr (default 2; 0 = no cut, no join). Bins 0.5 mag if every bin of the 1-99 % range holds >= NMIN, else 1.0.
Outputs (<out>): gal_ratio_vs_mag.png (one column per Dn4000 bin; rows: R vs MABS_G, R vs M_AB(5450), numbers per bin),
gal_ratio_vs_mag.csv (median, 16/84 %, n per DB x survey x z range x MABS_G bin), gal_ratio_lookup.csv (pooled over surveys
and z: DB x EB (-1 = all) x MABS_G bin), printed survey / z consistency.
Usage: python3 17_gal_ratio_vs_mag.py [--passc DIR] [--passg DIR] [--out DIR] [--snr 2]
"""
import os, glob, argparse, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib
pb = importlib.import_module('02_passB_measure')

NUM, DEN = '3050', '5450'
COVER = 0.8; NMIN = 100; C_AA = 2.99792458e18
SVCOL = {'DESI': 'C0', 'DR17': 'C3', 'DR8': 'C2'}; LS = ['-', '--', ':']
DBLAB = {0: 'Dn<1.15', 1: '1.15-1.40', 2: '1.40-1.65', 3: '1.65-1.90', 4: '>1.90'}


def attach_snr(t):
    """SNR_SDSSPIX for records written before it was a column: DESI from galaxy_cells (int64 TARGETID), SDSS from the plate tables."""
    d = t[t.SURVEY == 'DESI']
    if len(d):
        g = pd.read_parquet(os.path.expanduser('~/data/sdssdesi_products/galaxy_cells.parquet'), columns=['TARGETID', 'PROGRAM', 'SNR_SDSSPIX'])
        g = g.drop_duplicates(['TARGETID', 'PROGRAM'])
        d = d.merge(g, on=['TARGETID', 'PROGRAM'], how='left', validate='many_to_one')
    s = []
    for sv in ('dr8', 'dr17'):
        x = t[t.SURVEY == sv.upper()]
        if len(x):
            fs = glob.glob(os.path.expanduser(f'~/data/sdssdesi_passG_sdss_v2/tables/{sv}/*.parquet'))
            tb = pd.concat([pd.read_parquet(f, columns=['PLATE', 'MJD', 'FIBERID', 'SNR_SDSSPIX']) for f in fs], ignore_index=True)
            s.append(x.merge(tb.drop_duplicates(['PLATE', 'MJD', 'FIBERID']), on=['PLATE', 'MJD', 'FIBERID'], how='left', validate='many_to_one'))
    return pd.concat([d] + s, ignore_index=True)


def load(passc, passg, snr):
    parts = [pd.read_parquet(f) for f in sorted(glob.glob(os.path.join(passc, 'normtab', 'norm_shard*.parquet')))]
    parts += [pd.read_parquet(f) for f in sorted(glob.glob(os.path.join(passg, 'normtab', 'norm_dr*.parquet')))]
    for p in parts:                                  # SDSS rows have no TARGETID: nullable int, else the concat makes it float64
        if 'TARGETID' in p:
            p['TARGETID'] = p['TARGETID'].astype('Int64')
    t = pd.concat(parts, ignore_index=True)
    if snr <= 0:
        t['SNR_SDSSPIX'] = t.get('SNR_SDSSPIX', np.inf)
    elif 'SNR_SDSSPIX' not in t or t.SNR_SDSSPIX.isna().all():
        t = attach_snr(t.drop(columns='SNR_SDSSPIX', errors='ignore'))
    t = t[(t['COV' + NUM] >= COVER) & (t['COV' + DEN] >= COVER) & (t['W' + NUM] > 0) & (t['W' + DEN] > 0)].copy()
    t['logR'] = np.log10(t['W' + NUM] / t['W' + DEN])
    dm = pb.COSMO.distmod(np.clip(t.Z.to_numpy(float), 1e-4, None)).value
    fl = t['W' + DEN].to_numpy(float) / (t.NORM_SCALE * t.RESCALE).to_numpy(float)
    t['MSPEC'] = -2.5 * np.log10(fl * 5450.0**2 / C_AA) - 48.60 - dm
    return t


def choose_dm(x):
    """0.5 mag bins over the 1-99 % range if every bin holds >= NMIN, else 1.0 (as programs_quasars/10_ratio_vs_mag.py)."""
    lo, hi = np.percentile(x, [1, 99])
    for dm in (0.5, 1.0):
        e = np.arange(np.floor(lo / dm) * dm, hi + dm, dm); n, _ = np.histogram(x, e)
        if (n >= NMIN).all():
            return dm, e
    return 1.0, np.arange(np.floor(lo), hi + 1.0, 1.0)


def bstats(g, mcol, e, nmin=30):
    b = np.digitize(g[mcol], e) - 1; y = g.logR.to_numpy(); out = []
    for bi in range(len(e) - 1):
        v = y[b == bi]
        if len(v) >= nmin:
            q16, q50, q84 = np.percentile(v, [16, 50, 84])
            out.append(dict(M_lo=e[bi], M_hi=e[bi + 1], n=len(v), logR_med=q50, logR_16=q16, logR_84=q84))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--passc', default=os.path.expanduser('~/data/sdssdesi_passC_v2_2win'))
    ap.add_argument('--passg', default=os.path.expanduser('~/data/sdssdesi_passG_sdss_v2_2win'))
    ap.add_argument('--out', default=os.path.expanduser('~/data/sdssdesi_passC_v2_2win/qa')); ap.add_argument('--snr', type=float, default=2.0)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    t0 = load(a.passc, a.passg, a.snr)
    t = t0[(t0.SNR_SDSSPIX >= a.snr) & np.isfinite(t0.MABS_G)]
    print(f'R = F(3000-3100)/F(5400-5500): {len(t0):,} galaxies cover both windows, {len(t):,} with SNR_SDSSPIX >= {a.snr} and MABS_G; '
          f'z {t.Z.min():.2f}-{t.Z.max():.2f}; {t.SURVEY.value_counts().to_dict()}')
    dd = t.MSPEC - t.MABS_G
    print('M_AB(5450, spectrum) - MABS_G (fiber light lost): ' + ', '.join(f'{s} {g.median():+.2f}' for s, g in dd.groupby(t.SURVEY)))
    zq = np.round(np.percentile(t.Z, [33.3, 66.7]) / 0.05) * 0.05
    zr = [(np.floor(t.Z.min() / 0.05) * 0.05, zq[0]), (zq[0], zq[1]), (zq[1], np.ceil(t.Z.max() / 0.05) * 0.05)]

    rows, look = [], []
    fig, ax = plt.subplots(3, 5, figsize=(26, 14), gridspec_kw=dict(height_ratios=[2, 2, 1]))
    ylo, yhi = np.percentile(t.logR, [1, 99])
    for db in range(5):
        g0 = t[t.DB == db]
        if len(g0) < NMIN:
            [x.set_visible(False) for x in ax[:, db]]; continue
        for r, mcol in enumerate(('MABS_G', 'MSPEC')):
            dm, e = choose_dm(g0[mcol]); a0 = ax[r, db]
            for sv, col in SVCOL.items():
                for j, (z0, z1) in enumerate(zr):
                    g = g0[(g0.SURVEY == sv) & (g0.Z >= z0) & (g0.Z < z1)]
                    if len(g) < NMIN:
                        continue
                    st = bstats(g, mcol, e)
                    if not len(st):
                        continue
                    if mcol == 'MABS_G':
                        rows.append(st.assign(DB=db, SURVEY=sv, z_lo=z0, z_hi=z1))
                    xm = 0.5 * (st.M_lo + st.M_hi)
                    a0.plot(xm, st.logR_med, LS[j], color=col, marker='o', ms=3, label=f'{sv} z {z0:.2f}-{z1:.2f}')
                    if j == 1:
                        a0.fill_between(xm, st.logR_16, st.logR_84, color=col, alpha=0.08)
            a0.invert_xaxis(); a0.set_ylim(ylo, yhi); a0.legend(fontsize=7)
            a0.set_xlabel({'MABS_G': 'MABS_G (photometric, whole galaxy)', 'MSPEC': 'M_AB at 5450 A (from the spectrum: biased)'}[mcol])
            sec = a0.secondary_yaxis('right', functions=(lambda y: -2.5 * y, lambda m: -m / 2.5))
            if db == 4:
                sec.set_ylabel('m(3050) - m(5450), f_lam mag')
            if r == 0:
                a0.set_title(f'DB{db}  Dn4000 {DBLAB[db]}   ({dm} mag bins)', fontsize=11)
                for eb in [-1] + sorted(g0.EB.unique()):
                    ge = g0 if eb == -1 else g0[g0.EB == eb]
                    st = bstats(ge, 'MABS_G', e, NMIN)
                    if len(st):
                        look.append(st.assign(DB=db, EB=eb, R_med=10**st.logR_med,
                                              sem=[1.2533 * ge.logR[(ge.MABS_G >= l) & (ge.MABS_G < h)].std() / np.sqrt(n)
                                                   for l, h, n in zip(st.M_lo, st.M_hi, st.n)]))
                for sv, col in SVCOL.items():
                    gs = g0[g0.SURVEY == sv]
                    if len(gs):
                        ax[2, db].hist(gs.MABS_G, e, histtype='step', color=col, label=f'{sv} {len(gs):,}')
        ax[2, db].set_yscale('log'); ax[2, db].invert_xaxis(); ax[2, db].set_xlabel('MABS_G'); ax[2, db].legend(fontsize=8)
    for r in range(2):
        ax[r, 0].set_ylabel('log F(3000-3100)/F(5400-5500)  (median, 16-84 % band)')
    ax[2, 0].set_ylabel('galaxies per bin')
    fig.suptitle(f'R = F(3000-3100) / F(5400-5500) per galaxy, by Dn4000 bin (SNR_SDSSPIX >= {a.snr})', fontsize=14)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, 'gal_ratio_vs_mag.png'), dpi=80); plt.close(fig)
    B = pd.concat(rows, ignore_index=True); B.to_csv(os.path.join(a.out, 'gal_ratio_vs_mag.csv'), index=False)
    P = pd.concat(look, ignore_index=True); P.to_csv(os.path.join(a.out, 'gal_ratio_lookup.csv'), index=False)

    # consistency: each survey x z-range median minus the pooled (EB = -1) median at the same DB and MABS_G bin, n >= NMIN
    m = B[B.n >= NMIN].merge(P[P.EB == -1][['DB', 'M_lo', 'logR_med']].rename(columns={'logR_med': 'pool'}), on=['DB', 'M_lo'])
    m['dev'] = m.logR_med - m.pool
    print(f'consistency vs MABS_G: rms {np.sqrt(np.mean(m.dev**2)):.3f} dex over {len(m)} survey x z x DB x M cells')
    print('  by survey: ' + ', '.join(f'{s} {g.dev.median():+.3f} (rms {np.sqrt(np.mean(g.dev**2)):.3f}, {len(g)})' for s, g in m.groupby('SURVEY')))
    print('  by z range: ' + ', '.join(f'{z0:.2f}-{z1:.2f} {g.dev.median():+.3f} (rms {np.sqrt(np.mean(g.dev**2)):.3f})' for (z0, z1), g in m.groupby(['z_lo', 'z_hi'])))
    print('pooled median R by DB (all M): ' + '  '.join(f'DB{k} {10**v:.3f}' for k, v in t.groupby('DB').logR.median().items()))
    print('outputs in', a.out)


if __name__ == '__main__':
    main()
