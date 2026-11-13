"""Round 8 (Nao 2026-11-12): a SALT-like multiplicative colour term for the quasar spectra.
   model:  F(lambda) = A * [e0(lambda) + sum_k C_k e_k(lambda)] * (lambda / 3025 A)^c
   The additive part is the non-BAL eigenbasis (v4_nonbal) already used for the specphot coefficients; the colour term is a
   power-law tilt in F_lambda, unity at the 3025 A pivot so A keeps the meaning of NORM.  Fitted on the DESI z-slice spectra of
   27_specphot_ml_spectra.py (rest 1350-3650 A, 4-pixel bins): for every c on a grid the coefficients are weighted linear least
   squares; chi2(c) is minimized and refined by a parabola through the three lowest grid points.
   Step 1 ("make a plot first"): the colour against the original coefficients, the colour-free coefficients and the brightness.
   --mode joint : colour and all coefficients fitted together (step 1; colour degenerate with the basis, see ms.tex)
   --mode hier  : colour from the mean spectrum alone, F = A e0 (lambda/3025)^c, then the coefficients on F / colour term
   --mode ortho : the components e1..e9 projected orthogonal to the tilt direction e0 ln(lambda/3025), then the joint fit
   --mode salt  : SALT-style retraining on the slice: divide every spectrum by its colour term, weighted PCA of the bright half
                  (mean + 9 components, components orthogonalized against the tilt), refit colour + coefficients, --niter times.
   Outputs: /Volumes/exdisk3/sdssdesi_ml/qso_color_z<zlo>-<zhi>_<tag>.npz (per quasar: COLOR, A, CC1..CC9, CHI2_C, CHI2_0, C01..C09),
            csv_files/qso_color_corr_z<zlo>-<zhi>_<tag>.csv (correlations), plots_luminosity/plots_qso_color_z<zlo>-<zhi>_<tag>.png."""
import argparse, os, sys, time, importlib.util
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); import qsolib as Q
_s = importlib.util.spec_from_file_location('eigenfit', os.path.join(HERE, '16_eigenfit.py')); EF = importlib.util.module_from_spec(_s); _s.loader.exec_module(EF)
REPO = os.path.dirname(HERE)
ap = argparse.ArgumentParser(); ap.add_argument('--zlo', type=float, default=1.65); ap.add_argument('--zhi', type=float, default=1.70); ap.add_argument('--tag', default='v1')
ap.add_argument('--inp', default='/Volumes/exdisk3/sdssdesi_ml'); ap.add_argument('--basis', default=os.path.join(REPO, 'data_eigenspectra', 'quasar_eigenspectra_v4_nonbal.fits'))
ap.add_argument('--nvec', type=int, default=10); ap.add_argument('--cmax', type=float, default=2.5); ap.add_argument('--cstep', type=float, default=0.05); ap.add_argument('--bin', type=int, default=4)
ap.add_argument('--lmin', type=float, default=1350.0); ap.add_argument('--lmax', type=float, default=3650.0); ap.add_argument('--threads', type=int, default=4); ap.add_argument('--plot-only', action='store_true'); ap.add_argument('--mode', default='joint', choices=['joint', 'hier', 'ortho', 'salt'])
ap.add_argument('--niter', type=int, default=4); ap.add_argument('--train-all', action='store_true', help='salt: train the PCA on all quasars of the slice, not the bright half (selection check)'); ap.add_argument('--plot-coefs', default='', help='production | own (default: production for joint, own otherwise)'); a = ap.parse_args()
os.environ.setdefault('OMP_NUM_THREADS', str(a.threads))
tagz = f'z{a.zlo}-{a.zhi}_{a.tag}'; sfx = ('' if a.mode == 'joint' else '_' + a.mode) + ('_all' if a.train_all else ''); outnpz = f'{a.inp}/qso_color_{tagz}{sfx}.npz'; pc = a.plot_coefs or ('production' if a.mode == 'joint' else 'own')

z = np.load(f'{a.inp}/qso_slice_{tagz}.npz'); w = z['WAVE'].astype(np.float64); F = z['FLUX'].astype(np.float64); IV = z['IVAR'].astype(np.float64)
M = z['MAB'].astype(np.float64); C0 = np.stack([z[f'C{k}'] for k in range(1, 10)], axis=1).astype(np.float64); n = len(F)
if not a.plot_only:
    # basis on the slice grid: same cut and 4-pixel binning as 27_specphot_ml_spectra.py
    E, ok = EF.load_basis(a.basis, a.nvec); sel = (Q.QWAVE >= a.lmin) & (Q.QWAVE <= a.lmax); i0 = np.argmax(sel); n0 = sel.sum() - sel.sum() % a.bin
    Eb = E[:, i0:i0 + n0].reshape(a.nvec, -1, a.bin).mean(axis=2); okb = ok[i0:i0 + n0].reshape(-1, a.bin).all(axis=1)
    assert Eb.shape[1] == len(w), (Eb.shape, len(w))
    IV = IV * okb[None, :]; good = (IV > 0).sum(axis=1) > 200
    print(f'{n:,} spectra, {good.sum():,} with > 200 usable pixels; basis {a.nvec} vectors on {len(w)} bins', flush=True)
    cgrid = np.arange(-a.cmax, a.cmax + 1e-9, a.cstep); x = np.log(w / 3025.0); t0 = time.time()

    def colour_grid(Eb, F, IV):
        """chi2(c) on the grid with linear coefficients at each c -> refined colour per spectrum."""
        chi = np.full((n, len(cgrid)), np.inf)
        for j, c in enumerate(cgrid):
            Ec = Eb * np.exp(c * x)[None, :]
            for s_ in range(0, n, 5000):
                sl = slice(s_, s_ + 5000); co = EF.wls(Ec, F[sl], IV[sl]); chi[sl, j] = np.sum(IV[sl] * (F[sl] - co @ Ec)**2, axis=1)
        jmin = np.argmin(chi, axis=1); jm = np.clip(jmin, 1, len(cgrid) - 2); ar = np.arange(n)
        y0, y1, y2 = chi[ar, jm - 1], chi[ar, jm], chi[ar, jm + 1]; den = y0 - 2 * y1 + y2
        off = np.where(den > 0, 0.5 * (y0 - y2) / np.where(den > 0, den, 1), 0.0)
        return np.clip(cgrid[jm] + off * a.cstep, -a.cmax, a.cmax)

    def refit(Eb, cbest, F, IV):
        """coefficients of the colour-multiplied basis at the refined colour of each spectrum, and chi2."""
        CC = np.zeros((n, Eb.shape[0])); CH = np.zeros(n)
        for s_ in range(0, n, 5000):
            sl = slice(s_, s_ + 5000); T = np.exp(np.nan_to_num(cbest[sl])[:, None] * x[None, :]); Ecs = Eb[None] * T[:, None, :]
            Mm = np.einsum('bkm,bm,blm->bkl', Ecs, IV[sl], Ecs, optimize=True); bb = np.einsum('bkm,bm->bk', Ecs, IV[sl] * F[sl], optimize=True)
            Mm += 1e-12 * np.trace(Mm, axis1=1, axis2=2)[:, None, None] * np.eye(Eb.shape[0])[None]
            co = np.linalg.solve(Mm, bb[..., None])[..., 0]; CC[sl] = co; CH[sl] = np.sum(IV[sl] * (F[sl] - np.einsum('bk,bkm->bm', co, Ecs))**2, axis=1)
        return CC, CH

    def ortho(Eb):
        """components 1.. made orthogonal to the tilt direction e0 * ln(lambda/3025) (and to e0) on the usable pixels."""
        t = Eb[0] * x * okb; t = t - okb * np.dot(t, Eb[0] * okb) / np.dot(Eb[0] * okb, Eb[0] * okb); t /= np.linalg.norm(t)
        E2 = Eb.copy()
        for k in range(1, Eb.shape[0]):
            E2[k] = E2[k] - np.dot(E2[k] * okb, t) * t
        return E2

    def wpca(F, IV, cbest, train, nvec):
        """weighted PCA of F / colour term over the training rows: mean (normalized to 1 at 3000-3050 A) + nvec-1 components."""
        T = np.exp(cbest[:, None] * x[None, :]); G = F / T; Wt = IV * T**2
        G, Wt = G[train], Wt[train]; Wt = np.where(np.isfinite(G), Wt, 0); G = np.nan_to_num(G)
        mean = np.sum(Wt * G, axis=0) / np.maximum(np.sum(Wt, axis=0), 1e-30); piv = (w > 3000) & (w < 3050); mean /= np.median(mean[piv])
        amp = np.sum(Wt * G * mean[None], axis=1) / np.maximum(np.sum(Wt * mean[None]**2, axis=1), 1e-30)
        R = G / amp[:, None] - mean[None]; R[Wt == 0] = 0                      # missing pixels -> mean (zero residual)
        sw = np.sqrt(np.sum(Wt, axis=1) / np.maximum(np.sum(Wt > 0, axis=1), 1))[:, None]  # per-spectrum weight ~ mean sqrt(ivar)
        _, _, Vt = np.linalg.svd(R * sw, full_matrices=False)
        return np.vstack([mean, Vt[:nvec - 1]]) * okb[None, :]

    C00 = np.zeros((n, a.nvec)); CH0 = np.zeros(n)
    for s_ in range(0, n, 5000):
        sl = slice(s_, s_ + 5000); co0 = EF.wls(Eb, F[sl], IV[sl]); C00[sl] = co0; CH0[sl] = np.sum(IV[sl] * (F[sl] - co0 @ Eb)**2, axis=1)
    Euse = Eb
    if a.mode == 'joint':
        cbest = colour_grid(Eb, F, IV)
    elif a.mode == 'hier':
        cbest = colour_grid(Eb[:1], F, IV)                                     # colour from the mean spectrum alone
    elif a.mode == 'ortho':
        Euse = ortho(Eb); cbest = colour_grid(Euse, F, IV)
    elif a.mode == 'salt':
        Euse = ortho(Eb); cbest = colour_grid(Euse, F, IV)
        train = good & ((M < np.median(M)) | a.train_all) & (CH0 / np.maximum((IV > 0).sum(1), 1) < 3)
        for it in range(a.niter):
            Euse = ortho(wpca(F, IV, np.nan_to_num(cbest), train, a.nvec)); cbest = colour_grid(Euse, F, IV)
            CCi, CHi = refit(Euse, cbest, F, IV); gi = good & (np.abs(cbest) < a.cmax - 1e-6)
            print(f'  salt iteration {it + 1}: colour 16-84 % {np.percentile(cbest[gi], 16):+.3f} .. {np.percentile(cbest[gi], 84):+.3f}, '
                  f'chi2/pix median {np.median(CHi[gi] / (IV[gi] > 0).sum(1)):.4f}, rho_S(c, c1) {spearmanr(cbest[gi], CCi[gi, 1] / CCi[gi, 0])[0]:+.3f}  {time.time() - t0:.0f}s', flush=True)
        np.savez(f'{a.inp}/qso_color_basis_{tagz}{sfx}.npz', WAVE=w, E=Euse, OK=okb)
    cbest[~good] = np.nan
    CC, CH = refit(Euse, cbest, F, IV)
    npix = (IV > 0).sum(axis=1)
    with np.errstate(all='ignore'):
        A = CC[:, 0]; CCn = CC[:, 1:] / A[:, None]; A0 = C00[:, 0]; C0n = C00[:, 1:] / A0[:, None]
    np.savez(outnpz, TARGETID=z['TARGETID'], Z=z['Z'], MAB=M, COLOR=cbest, A=A, A0=A0, CHI2_C=CH, CHI2_0=CH0, NPIX=npix,
             **{f'CC{k}': CCn[:, k - 1] for k in range(1, 10)}, **{f'C0{k}': C0n[:, k - 1] for k in range(1, 10)})
    print(f'written {outnpz}  ({time.time() - t0:.0f}s)', flush=True)

r = np.load(outnpz); c = r['COLOR']; CC = np.stack([r[f'CC{k}'] for k in range(1, 10)], 1); Cz = np.stack([r[f'C0{k}'] for k in range(1, 10)], 1)
g = np.isfinite(c) & np.isfinite(CC).all(1) & (np.abs(c) < a.cmax - 1e-6)
dchi = (r['CHI2_0'] - r['CHI2_C']) / np.maximum(r['NPIX'], 1)
print(f'colour: median {np.nanmedian(c[g]):+.3f}, 16-84 % {np.nanpercentile(c[g], 16):+.3f} .. {np.nanpercentile(c[g], 84):+.3f}; at the grid edge {(~g & np.isfinite(c)).sum():,}')
print(f'chi2/pix improvement from the colour term: median {np.median(dchi[g]):.4f}')
bright = M < np.median(M)
rows = []
for name, X in [('C (production)', C0), ('C (c=0 refit)', Cz), ('C (colour-free)', CC)]:
    for k in range(9):
        rows.append(dict(set=name, coef=f'c{k + 1}', pearson_colour=pearsonr(c[g], X[g, k])[0], spearman_colour=spearmanr(c[g], X[g, k])[0],
                         spearman_M=spearmanr(M[g], X[g, k])[0], spearman_M_bright=spearmanr(M[g & bright], X[g & bright, k])[0]))
rows.append(dict(set='colour', coef='c', pearson_colour=1, spearman_colour=1, spearman_M=spearmanr(M[g], c[g])[0], spearman_M_bright=spearmanr(M[g & bright], c[g & bright])[0]))
T = pd.DataFrame(rows); T.to_csv(f'{REPO}/csv_files/qso_color_corr_{tagz}{sfx}.csv', index=False); print(T.round(3).to_string())

def panel(ax, xx, yy, xl, yl, invert=False):
    ok = np.isfinite(xx) & np.isfinite(yy); xx, yy = xx[ok], yy[ok]
    xr = np.percentile(xx, [0.5, 99.5]); yr = np.percentile(yy, [0.5, 99.5])
    ax.hist2d(xx, yy, bins=120, range=[xr, yr], cmap='Greys', norm=matplotlib.colors.LogNorm())
    edges = np.linspace(*xr, 25); mid = 0.5 * (edges[1:] + edges[:-1]); med = np.full(len(mid), np.nan); lo = med.copy(); hi = med.copy()
    for i in range(len(mid)):
        s = (xx >= edges[i]) & (xx < edges[i + 1])
        if s.sum() > 50: lo[i], med[i], hi[i] = np.percentile(yy[s], [16, 50, 84])
    ax.plot(mid, med, 'r-', lw=1.5); ax.plot(mid, lo, 'r:', lw=1); ax.plot(mid, hi, 'r:', lw=1)
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.text(0.03, 0.95, f'$\\rho_S$ = {spearmanr(xx, yy)[0]:+.2f}', transform=ax.transAxes, va='top', fontsize=9)
    if invert: ax.invert_yaxis()

fig, axs = plt.subplots(3, 4, figsize=(16, 11.5))
ax = axs[0, 0]; ax.hist(c[g], bins=100, range=(-a.cmax, a.cmax), color='k', histtype='step'); ax.set_xlabel('colour $c$  [$F_\\lambda \\propto (\\lambda/3025)^c$]'); ax.set_ylabel('quasars')
ax.set_title(f'DESI z = {a.zlo}-{a.zhi}: {g.sum():,} quasars', fontsize=10)
panel(axs[0, 1], c[g], M[g], 'colour $c$', '$M_{\\rm AB}$(3025 Å)', invert=True)
panel(axs[0, 2], c[g], dchi[g], 'colour $c$', '$\\Delta\\chi^2$ / pixel (c = 0 minus fit)')
panel(axs[0, 3], C0[g, 0], CC[g, 0], '$c_1$ production', '$c_1$ colour-free')
Xp, lab = (C0, 'production') if pc == 'production' else (CC, 'colour-free')
for k in range(8):
    panel(axs[1 + k // 4, k % 4], Xp[g, k], c[g], f'$c_{k + 1}$ ({lab})', 'colour $c$')
fig.suptitle(f'Round 8 ({a.mode}): multiplicative colour term  $F = A\\,[e_0 + \\sum_k c_k e_k]\\,(\\lambda/3025)^c$ — colour against the {lab} coefficients', fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.97)); out = f'{REPO}/plots_luminosity/plots_qso_color_{tagz}{sfx}.png'; fig.savefig(out, dpi=110); print(out)

# ---- second figure: brightness (y, bright at the top) against the colour and the colour-free coefficients (x)
fig, axs = plt.subplots(2, 5, figsize=(19, 8))
panel(axs[0, 0], c[g], M[g], 'colour $c$', '$M_{\\rm AB}$(3025 Å)', invert=True)
for k in range(9):
    panel(axs.flat[k + 1], CC[g, k], M[g], f'$c_{k + 1}$ (colour-free)', '$M_{\\rm AB}$(3025 Å)', invert=True)
fig.suptitle(f'Round 8 ({a.mode}): brightness against the colour and the colour-free coefficients, DESI z = {a.zlo}-{a.zhi}, {g.sum():,} quasars', fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.97)); out = f'{REPO}/plots_luminosity/plots_qso_color_vs_mag_{tagz}{sfx}.png'; fig.savefig(out, dpi=110); print(out)
