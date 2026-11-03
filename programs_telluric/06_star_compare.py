#!/usr/bin/env python3
"""Telluric step 3 for the stars: the corrected stacks and bases (v2_tell) against v2, same stars, same pixels (Nao 2026-11-02). Each table in csv_files/plots_telluric_starcompare_<name>.csv
T0 here is what the members really carried: each star's T0 at its own observed wavelengths, averaged over the stack (stars are stacked in their rest frame, the 1 to 99 % range of velocities is -3 to +2 pixels).
sanity  corrected / v2 of each stack: 1 outside the bands, 1 / (1 + mean T0) inside, up to the DR8 a_k.
slope   T0 amplitude left in each stack: ln f regressed on T0 with the per-band nuisance curve of tellapply, 4 sigma clip. A featureless star should go from 1 to 0; WD and HOT are the test, the other classes have lines in the bands.
cross   independent of T0 fit: ln(SDSS / DESI) of the same cell, where the star cancels. rms in the bands and in control windows, noise removed; slope on T0(SDSS) - T0(DESI).
amp     rms amplitude of each telluric vector in the training set (stacks, stars) as held by the basis, against the same vector moved 20 band pixels off the atmosphere.
heldout coefficients fitted outside the bands and controls, rms of residual / mean inside each of them with 10 vectors (the first 5 are not comparable: EM-PCA does not sort, see space).
space   overlap of the first n vectors of the two bases (smallest cosine of the principal angles) and R2.
ak      DR8 a_k per class and against colour: does the fit take stellar features as atmosphere.
"""
import os, sys
import numpy as np, pandas as pd, fitsio

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path[:0] = [HERE, os.path.join(ROOT, 'programs_stars')]
import tellapply, tell_hook
TELL = tell_hook.OUT; V2 = tell_hook.V2; CSV = os.path.join(ROOT, 'csv_files'); WAVE = 10**(3.5502 + 1e-4 * np.arange(4447)); RCLIP = 0.3; SHIFT = 20
CLASSES = {'sdss': ['WD', 'AF', 'GK', 'M'], 'desi': ['WD', 'HOT', 'AF', 'GK', 'M']}; CONTROL = [(5950, 6200), (6330, 6540), (7010, 7150), (7720, 8110), (8340, 8900)]


def tag(s):
    return f'{s}_0.5x0.1_gspphot'


def stacks(d, s):
    F = fitsio.FITS(f'{d}/stacks_{tag(s)}.fits'); o = dict(f=F['FLUX'].read().astype(float), e=F['ERR'].read().astype(float), n=F['NSPEC'].read().astype(float), c=F['CELLS'].read()); F.close()
    o['name'] = np.char.strip(o['c']['STACK'].astype(str)); o['cell'] = np.char.strip(o['c']['CELL'].astype(str)); return o


def base(d, s, c, sfx=''):
    p = f'{d}/empca_{s}_0.5x0.1_gspphot_{c}{sfx}.fits'; b = {e: fitsio.read(p, ext=e).astype(float) for e in ('MEAN', 'EIGVEC', 'COEFF', 'VARIANCE', 'PIXMASK')}
    n = fitsio.read(p, ext='NORM'); b['id'] = np.char.strip(n['ID'].astype(str)); b['kind'] = np.char.strip(n['KIND'].astype(str)); b['norm'] = n['NORM_GFLUX'].astype(float); b['col'] = n['COL0_CELL'].astype(float); b['ok'] = b['PIXMASK'] > 0; return b


def hats(pix, nodes):
    """continuous piecewise-linear curve per band on the pixels pix: the nuisance of tellapply."""
    w = WAVE[pix]; D = []
    for n in nodes:
        m = (w >= n[0]) & (w <= n[-1])
        for j in range(len(n)):
            h = np.zeros(len(w)); h[m] = np.interp(w[m], n, np.eye(len(n))[j]); D.append(h)
    return np.array(D)


def lsq(D, y, W):
    A = (D * W) @ D.T; A[np.diag_indices_from(A)] += 1e-9 * max(np.trace(A) / len(A), 1e-30); return np.linalg.solve(A, (D * W) @ y)


def fit_slope(y, W, T, H):
    """(slope, kept pixels) of y on T with the nuisance H, two passes of a 4 sigma clip on the residual."""
    keep = (W > 0) & np.isfinite(y)
    if keep.sum() < 30: return np.nan, 0
    for _ in range(3):
        D = np.vstack([T, H])[:, keep]; p = lsq(D, y[keep], W[keep]); r = y - np.vstack([T, H]).T @ p; s = 1.4826 * np.nanmedian(np.abs(r[keep]))
        keep = (W > 0) & np.isfinite(y) & (np.abs(r) < 4 * s)
    return p[0], int(keep.sum())


def smeared_t0(s):
    """stack name -> mean T0 of its members at their own velocities, the factor step 6 divided them by when a_k = 0."""
    m = pd.read_csv(f'{TELL}/members_{tag(s)}.csv', usecols=['SID', 'STACK']); t = pd.read_parquet(f'{TELL}/tell_coeff_{tag(s)}.parquet', columns=['SID', 'SURVEY', 'LOGZ']); g = m.merge(t, on='SID')
    return tell_hook.smeared_t0(g.STACK.to_numpy(), g.SURVEY.to_numpy(), g.LOGZ.to_numpy())


def lnw(st, i, pix):
    f, e = st['f'][i, pix], st['e'][i, pix]; ok = (st['n'][i, pix] > 0) & np.isfinite(f) & np.isfinite(e) & (f > 0) & (e > 0)
    y = np.full(len(pix), np.nan); v = np.zeros(len(pix)); y[ok] = np.log(f[ok]); v[ok] = (f[ok] / e[ok])**2; return y, v


def classes_of(s):
    out = {}
    for c in CLASSES[s]:
        b = base(TELL, s, c)
        for i in b['id'][b['kind'] == 'stack']: out[i] = c
    return out


def sanity_slope(s, A, B, cl, Ts, B8, Bd):
    rows = []; bas = Bd if s == 'desi' else B8; nodes = bas.bands; pix = np.where(bas.ok)[0]; inb = np.zeros(4447, bool); inb[pix] = True
    span = np.zeros(4447, bool)
    for n in nodes: span |= (WAVE >= n[0]) & (WAVE <= n[-1])
    H = hats(pix, nodes)
    for i, name in enumerate(A['name']):
        c = cl.get(name)
        if c is None: continue
        T = Ts[name]; ya, wa = lnw(A, i, pix); yb, wb = lnw(B, i, pix); good = (A['n'][i] > 0) & (A['f'][i] > 0) & (B['f'][i] > 0)
        ra = A['f'][i][good] / B['f'][i][good]; o = ~span[good]; pred = 1 / (1 + T[good])
        row = dict(survey=s, cls=c, stack=name, n=np.nanmedian(A['n'][i][good]), out_med=np.nanmedian(ra[o] - 1), out_max=np.nanmax(np.abs(ra[o] - 1)), in_dev=np.sqrt(np.nanmean((ra[~o] / pred[~o] - 1)**2)), in_t0=np.sqrt(np.nanmean((pred[~o] - 1)**2)))
        row['slope_v2'], _ = fit_slope(yb, wb, T[pix], H); row['slope_tell'], row['npix'] = fit_slope(ya, wa, T[pix], H)
        for j, n in enumerate(nodes):
            m = (WAVE[pix] >= n[0]) & (WAVE[pix] <= n[-1]); Hb = hats(pix[m], [n]); row[f'v2_{n[0]:.0f}'], _ = fit_slope(yb[m], wb[m], T[pix][m], Hb); row[f'tell_{n[0]:.0f}'], _ = fit_slope(ya[m], wa[m], T[pix][m], Hb)
        rows.append(row)
    return pd.DataFrame(rows)


def to_cells(st):
    """N-weighted mean of the stacks of a cell: the splits differ between surveys, the cells do not."""
    out = {}
    for cell in np.unique(st['cell']):
        i = np.where(st['cell'] == cell)[0]; n = np.where(np.isfinite(st['f'][i]) & (st['f'][i] > 0) & np.isfinite(st['e'][i]) & (st['e'][i] > 0), st['n'][i], 0.0); N = n.sum(0)
        with np.errstate(all='ignore'):
            f = (n * np.nan_to_num(st['f'][i])).sum(0) / N; e = np.sqrt((n**2 * np.nan_to_num(st['e'][i])**2).sum(0)) / N
        out[str(cell)] = (f, e, N)
    return out


def cross(S, D, Sv, Dv, cl, Ts, Bd):
    """ln(SDSS / DESI) per cell, v2 and corrected, on DESI band pixels; returns rows and the stacks of residuals for the figure."""
    cs, cd, vs, vd = to_cells(S), to_cells(D), to_cells(Sv), to_cells(Dv); common = sorted(set(cs) & set(cd)); nb = Bd.bands; windows = [('band', nb), ('control', [np.linspace(lo, hi, int(np.ceil((hi - lo) / 300)) + 1) for lo, hi in CONTROL])]
    names, cells, dn, dc = S['name'], S['cell'], D['name'], D['cell']; rows = []; res = {'v2': [], 'tell': []}
    allb = np.zeros(4447, bool)
    for n in nb: allb |= (WAVE >= n[0]) & (WAVE <= n[-1])
    bp = np.where(allb & (WAVE < 9200))[0]; Hb = hats(bp, nb)
    for cell in common:
        (fa, ea, na), (fb, eb, nb2), (ga, ha, _), (gb, hb, _) = cs[cell], cd[cell], vs[cell], vd[cell]; dT = np.mean([Ts['sdss'][x] for x in names[cells == cell]], axis=0) - np.mean([Ts['desi'][x] for x in dn[dc == cell]], axis=0); row = dict(cell=cell, cls=cl.get(str(names[cells == cell][0])), nsdss=np.nanmedian(na), ndesi=np.nanmedian(nb2))
        for lab, (a, ae, b, be) in (('tell', (fa, ea, fb, eb)), ('v2', (ga, ha, gb, hb))):
            with np.errstate(all='ignore'):
                y = np.log(a / b); nv = (ae / a)**2 + (be / b)**2
            for wn, nodes in windows:
                pix = np.where(np.any([(WAVE >= n[0]) & (WAVE <= n[-1]) for n in nodes], axis=0) & (WAVE < 9200))[0]; ok = np.isfinite(y[pix]) & np.isfinite(nv[pix]) & (nv[pix] > 0) & (np.abs(y[pix] - np.nanmedian(y[pix])) < 1)
                if ok.sum() < 50: continue
                H = hats(pix, nodes); W = np.where(ok, 1 / np.where(ok, nv[pix], 1), 0.0); r = y[pix] - H.T @ lsq(H[:, ok], y[pix][ok], W[ok]); s = 1.4826 * np.median(np.abs(r[ok])); k = ok & (np.abs(r) < 4 * s)
                row[f'{wn}_rms_{lab}'] = np.sqrt(np.mean(r[k]**2)); row[f'{wn}_noise_{lab}'] = np.sqrt(np.mean(nv[pix][k]))
            ok = np.isfinite(y[bp]) & np.isfinite(nv[bp]) & (nv[bp] > 0); w = np.where(ok, 1 / np.where(ok, nv[bp], 1), 0.0); yy = np.where(ok, y[bp], np.nan); row[f'slope_{lab}'], _ = fit_slope(yy, w, dT[bp], Hb)
            r = yy - np.vstack([Hb]).T @ lsq(Hb[:, ok], yy[ok], w[ok]) if ok.sum() > 50 else yy * np.nan; res[lab].append(r)
        rows.append(row)
    return pd.DataFrame(rows), bp, {k: np.array(v) for k, v in res.items()}


def amp(s, c, a, b, B, Bs):
    """rms amplitude of each vector T in the data as the first 10 vectors hold them: x = coeff @ eigvec, T scaled by the mean spectrum because telluric multiplies."""
    rows = []
    for lab, q in (('v2', b), ('tell', a)):
        pix = np.where(q['ok'] & B.ok)[0]; H = hats(pix, B.bands); H = H[H.sum(1) > 0]; m = q['MEAN'][pix]; X = q['COEFF'] @ q['EIGVEC'][:, pix]; one = np.ones(len(pix))
        for name, T0, T1 in [('T0', B.Ffull[0], Bs.Ffull[0])] + [(B.labels[k], B.T[k], Bs.T[k]) for k in range(len(B.T))]:
            for ctl, T in (('real', T0), ('shifted', T1)):
                v = m * T[pix]; v = v - H.T @ lsq(H, v, one); nv = np.dot(v, v); z = X @ v / nv     # v is orthogonal to the nuisance H, so the projection needs no detrending of X
                for kind in ('stack', 'star'):
                    w = q['kind'] == kind
                    if w.any(): rows.append(dict(survey=s, cls=c, basis=lab, T=name, control=ctl, kind=kind, rms=np.std(z[w])))
    return rows


def heldout(s, c, a, b, st_a, st_b, ia, ib, B):
    """the band pixels and controls never enter the coefficients. rms of residual / mean there, per object, then the median; err gives the stack noise on the same scale."""
    rows = []; win = {f'{n[0]:.0f}-{n[-1]:.0f}': (WAVE >= n[0]) & (WAVE <= n[-1]) for n in B.bands}; band = np.any(list(win.values()), axis=0); ctrl = np.any([(WAVE >= lo) & (WAVE <= hi) for lo, hi in CONTROL], axis=0); win['all bands'] = band; win['control'] = ctrl
    for lab, q, st, ind in (('v2', b, st_b, ib), ('tell', a, st_a, ia)):
        fit = q['ok'] & ~band & ~ctrl; ids = {x: j for j, x in enumerate(st['name'])}; F, W, I = ind; idx = {x: j for j, x in enumerate(I)}; out = {}; E = q['EIGVEC'][:10]
        for i, kind, nm in zip(q['id'], q['kind'], q['norm']):
            if kind == 'stack':
                j = ids.get(i)
                if j is None: continue
                g = (st['n'][j] > 0) & np.isfinite(st['f'][j]) & (st['e'][j] > 0); x = np.where(g, st['f'][j], 0.0) / nm; w = np.where(g, (nm / np.where(g, st['e'][j], 1))**2, 0.0)
            else:
                j = idx.get(i)
                if j is None: continue
                x = np.nan_to_num(F[j]) / nm; w = W[j]
            r0 = x - q['MEAN']; u = fit & (w > 0); r = r0 - E.T @ lsq(E[:, u], r0[u], w[u]); m = np.abs(q['MEAN']).clip(1e-3)
            for wn, mask in win.items():
                k = mask & q['ok'] & (w > 0) & (np.abs(r / m) < RCLIP)
                if k.sum() > 20:
                    v = np.sum(w[k] * (r[k] / m[k])**2) / np.sum(w[k]); out.setdefault((kind, wn), []).append(np.sqrt(v))
        rows += [dict(survey=s, cls=c, basis=lab, kind=kk, window=wn, rms=np.median(v), n=len(v)) for (kk, wn), v in out.items()]
    return rows


def space(s, c, a, b):
    k = a['ok'] & b['ok']; row = dict(survey=s, cls=c)
    for n in (1, 2, 3, 5, 10):
        A = np.linalg.qr(a['EIGVEC'][:n, k].T)[0]; B = np.linalg.qr(b['EIGVEC'][:n, k].T)[0]; row[f'mincos{n}'] = np.linalg.svd(A.T @ B, compute_uv=False).min(); row[f'r2{n}_tell'] = a['VARIANCE'][n - 1]; row[f'r2{n}_v2'] = b['VARIANCE'][n - 1]
    return row


def ak(cl):
    m = pd.read_csv(f'{TELL}/members_{tag("sdss")}.csv', usecols=['SID', 'STACK']); t = pd.read_parquet(f'{TELL}/tell_coeff_{tag("sdss")}.parquet'); t = t[(t.SURVEY == 'dr8') & (t.NPIX > 0)].merge(m, on='SID')
    C = fitsio.read(f'{TELL}/stacks_{tag("sdss")}.fits', ext='CELLS'); col = dict(zip(np.char.strip(C['STACK'].astype(str)), C['COL0_MED'].astype(float))); t['cls'] = t.STACK.map(cl); t['col'] = t.STACK.map(col); B = tellapply.basis('dr8'); rows = []
    a = [f'A{k}' for k in range(len(B.T))]
    for c, g in t.groupby('cls'):
        for k, n in zip(a, B.labels): rows.append(dict(cls=c, group='all', T=n, n=len(g), med=g[k].median(), rms=1.4826 * (g[k] - g[k].median()).abs().median(), r_col=np.corrcoef(g[k].clip(*g[k].quantile([.01, .99])), g['col'])[0, 1]))
    M = t[t.cls == 'M'].copy(); M['bin'] = pd.qcut(M['col'], 4, labels=['M blue', 'M 2', 'M 3', 'M red'])
    for c, g in M.groupby('bin', observed=True):
        for k, n in zip(a, B.labels): rows.append(dict(cls='M', group=str(c), T=n, n=len(g), med=g[k].median(), rms=1.4826 * (g[k] - g[k].median()).abs().median()))
    return pd.DataFrame(rows)


def individuals(s):
    k = 5; sn = 20.0; z = np.load(os.path.join(V2, f'individuals_{tag(s)}_{k}_{sn:g}.npz')); v = (z['FLUX'].astype(float), z['W'].astype(float), z['IDS'].astype(str).tolist()); F, W, I, _ = tell_hook.individuals(tag(s), V2, k, sn); return (F, W, I), v


def boot(x, n=500):
    x = np.asarray(x, float); x = x[np.isfinite(x)]; g = np.random.default_rng(1); return np.median(x), np.std([np.median(g.choice(x, len(x))) for _ in range(n)])


def main():
    pd.set_option('display.width', 220); pd.set_option('display.max_columns', 40); pd.set_option('display.max_rows', 200)
    B8, B17, Bd = (tellapply.basis(x) for x in ('dr8', 'dr17', 'desi')); Bs = {'sdss': tellapply.Basis('dr8', shift=SHIFT), 'desi': tellapply.Basis('desi', shift=SHIFT)}; Ts = {s: smeared_t0(s) for s in CLASSES}; data = {s: (stacks(TELL, s), stacks(V2, s)) for s in CLASSES}; cl = {s: classes_of(s) for s in CLASSES}; ss = []; am = []; he = []; sp = []
    for x in Bs.values():
        i = np.where(x.ok)[0]; x.Ffull[0, i] = np.roll(x.Ffull[0, i], SHIFT)     # Basis moves only the eigenvectors
    for s in CLASSES:
        A, Bb = data[s]; print(s, 'same stacks', (A['name'] == Bb['name']).all(), flush=True); t = sanity_slope(s, A, Bb, cl[s], Ts[s], B8, Bd); ss.append(t); ia, ib = individuals(s); bas = Bd if s == 'desi' else B8
        for c in CLASSES[s]:
            a, b = base(TELL, s, c), base(V2, s, c); am += amp(s, c, a, b, bas, Bs[s]); he += heldout(s, c, a, b, A, Bb, ia, ib, bas); sp.append(space(s, c, a, b)); print(' ', s, c, flush=True)
    ss = pd.concat(ss); ss.to_csv(f'{CSV}/plots_telluric_starcompare_slope.csv', index=False); am = pd.DataFrame(am); am.to_csv(f'{CSV}/plots_telluric_starcompare_amp.csv', index=False); he = pd.DataFrame(he); he.to_csv(f'{CSV}/plots_telluric_starcompare_heldout.csv', index=False); sp = pd.DataFrame(sp); sp.to_csv(f'{CSV}/plots_telluric_starcompare_space.csv', index=False)
    cr, bp, res = cross(data['sdss'][0], data['desi'][0], data['sdss'][1], data['desi'][1], cl['sdss'], Ts, Bd); cr.to_csv(f'{CSV}/plots_telluric_starcompare_cross.csv', index=False); np.savez(os.path.join(TELL, 'starcompare_cross.npz'), wave=WAVE[bp], v2=res['v2'], tell=res['tell'], cell=cr.cell.values.astype(str), noise=cr.band_noise_tell.values)
    k = ak(cl['sdss']); k.to_csv(f'{CSV}/plots_telluric_starcompare_ak.csv', index=False)
    print('\nSANITY'); print(ss.groupby(['survey', 'cls'])[['out_med', 'out_max', 'in_dev', 'in_t0']].median().round(6))
    print('\nSLOPE on T0, median over stacks'); print(ss.groupby(['survey', 'cls'])[['slope_v2', 'slope_tell']].median().round(3).join(ss.groupby(['survey', 'cls']).size().rename('n')))
    print('\nSLOPE per band: T0 rms, v2 and corrected, median over stacks +- bootstrap')
    for (sv, c), g in ss.groupby(['survey', 'cls']):
        bas = Bd if sv == 'desi' else B8
        for n in bas.bands:
            i = (WAVE >= n[0]) & (WAVE <= n[-1]) & bas.ok; print(f' {sv} {c:3s} {n[0]:.0f}-{n[-1]:.0f} T0 rms {np.sqrt(np.mean(bas.Ffull[0][i]**2)):.5f}  v2 %6.2f+-%.2f  tell %6.2f+-%.2f' % (*boot(g[f'v2_{n[0]:.0f}']), *boot(g[f'tell_{n[0]:.0f}'])))
    print('\nCROSS, cells', len(cr)); e = lambda a, b: np.sqrt(np.maximum(a**2 - b**2, 0)); good = cr[cr.band_noise_tell < 0.02]; print(' with noise < 2 %', len(good))
    for g, lab in ((cr, 'all'), (good, 'bright')):
        for w in ('band', 'control'):
            print(f' {lab:6s} {w:7s} rms v2 %.5f+-%.5f tell %.5f+-%.5f noise %.5f; excess v2 %.5f tell %.5f' % (*boot(g[f'{w}_rms_v2']), *boot(g[f'{w}_rms_tell']), g[f'{w}_noise_tell'].median(), np.nanmedian(e(g[f'{w}_rms_v2'], g[f'{w}_noise_v2'])), np.nanmedian(e(g[f'{w}_rms_tell'], g[f'{w}_noise_tell']))))
        print(f' {lab:6s} slope on dT0: v2 %.3f+-%.3f tell %.3f+-%.3f' % (*boot(g.slope_v2), *boot(g.slope_tell)))
    print(cr.groupby('cls')[['band_rms_v2', 'band_rms_tell', 'band_noise_tell', 'control_rms_v2', 'control_rms_tell', 'slope_v2', 'slope_tell']].median().round(5).join(cr.groupby('cls').size().rename('n')))
    print('\nAMP, stacks and stars, rms'); p = am.pivot_table(index=['survey', 'cls', 'kind', 'T'], columns=['basis', 'control'], values='rms'); print(p.round(4))
    print('\nHELD OUT rms of residual / mean, 10 vectors'); print(he.pivot_table(index=['survey', 'cls', 'kind', 'window'], columns='basis', values='rms').round(5))
    print('\nSPACE'); print(sp.round(4)); print('\nAK'); print(k[k.group == 'all'].pivot(index='T', columns='cls', values='med').round(4)); print(k[k.group == 'all'].pivot(index='T', columns='cls', values='r_col').round(3)); print(k[(k.group != 'all')].pivot(index='T', columns='group', values='med').round(4))


if __name__ == '__main__':
    main()
