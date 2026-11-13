"""Pixel-level brightness predictors on the rest-frame spectra of a redshift slice (27_specphot_ml_spectra.py output).
Input: the spectrum in units of its own 2950-3100 A median (shape only; the amplitude is the target). Variants of the
input, to separate shape information from the noise level (at fixed z the per-pixel noise measures the apparent flux):
  raw       flux and a second channel sqrt(ivar) * flux-scale (the noise level is visible);
  fluxonly  flux alone (the noise level is still readable from the pixel scatter);
  equalized flux alone after adding Gaussian noise so that every spectrum has the same S/N per pixel in the window,
            the S/N of the --snfloor quantile of the sample (the noise level no longer tells the brightness).
Models, 5-fold cross-validation with the same folds: ridge on the pixels; PCA(--npca) + ridge; PCA + HGB; PCA + MLP; a 1D
CNN (4 conv blocks, global pooling, 2 dense) in torch on CPU. Metrics on held-out quasars: R^2, rms, rms of the bright
half, 68 % |residual|. Rows appended to csv_files/qso_ml_pix_z<zlo>-<zhi>_<tag>.csv (read with the bench csv by 30_specphot_ml_plot.py) with method names 'pix: ...'."""
import argparse, os, sys, time
import numpy as np, pandas as pd
from sklearn.linear_model import Ridge
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(); ap.add_argument('--zlo', type=float, default=1.65); ap.add_argument('--zhi', type=float, default=1.70); ap.add_argument('--tag', default='v1')
ap.add_argument('--inp', default='/Volumes/exdisk3/sdssdesi_ml'); ap.add_argument('--nfold', type=int, default=5); ap.add_argument('--npca', type=int, default=50); ap.add_argument('--snfloor', type=float, default=0.1)
ap.add_argument('--epochs', type=int, default=30); ap.add_argument('--threads', type=int, default=6); ap.add_argument('--variants', default='raw,fluxonly,equalized'); ap.add_argument('--models', default='ridge,pca_ridge,pca_hgb,pca_mlp,cnn'); a = ap.parse_args()
z = np.load(f'{a.inp}/qso_slice_z{a.zlo}-{a.zhi}_{a.tag}.npz'); F = z['FLUX'].astype(np.float32); IV = z['IVAR'].astype(np.float32); w = z['WAVE']; M = z['MAB'].astype(np.float64)
n, npix = F.shape; print(f'{n:,} spectra x {npix} pixels, {w[0]:.0f}-{w[-1]:.0f} A; sigma(M) {M.std():.3f}', flush=True)
rng = np.random.default_rng(1); fold = rng.integers(0, a.nfold, n); bright = M < np.median(M)
F = np.clip(F, -5, 20); win = (w > 2950) & (w < 3100)
sn = np.array([np.median(np.sqrt(iv[win][iv[win] > 0])) if (iv[win] > 0).any() else 0 for iv in IV])     # S/N per pixel in the window (flux = 1 there)
sn_floor = np.quantile(sn[sn > 0], a.snfloor); print(f'S/N per pixel in the window: median {np.median(sn):.1f}, floor (q={a.snfloor}) {sn_floor:.1f}', flush=True)


def inputs(variant):
    if variant == 'raw':
        return np.stack([F, np.sqrt(IV) / 10.0], axis=1)                           # [n, 2, npix]
    if variant == 'fluxonly':
        return F[:, None, :]
    if variant == 'equalized':
        with np.errstate(all='ignore'):
            sig_have = np.where(IV > 0, 1.0 / np.sqrt(np.where(IV > 0, IV, 1.0)), 0.0)   # per-pixel sigma now
            sig_want = np.where(IV > 0, sig_have * (sn[:, None] / sn_floor), 0.0)         # scale every spectrum's noise up to the floor S/N
            add = np.sqrt(np.clip(sig_want**2 - sig_have**2, 0, None))
        g = np.random.default_rng(2).standard_normal(F.shape).astype(np.float32)
        return (F + add * g)[:, None, :]
    raise ValueError(variant)


class CNN:
    def __init__(self, epochs): self.epochs = epochs
    def fit(self, X, y):
        import torch, torch.nn as nn
        torch.manual_seed(0); torch.set_num_threads(a.threads); nch = X.shape[1]
        self.my, self.sy = y.mean(), y.std(); X = torch.tensor(X); Y = torch.tensor((y - self.my) / self.sy, dtype=torch.float32)[:, None]
        def block(i, o): return nn.Sequential(nn.Conv1d(i, o, 7, padding=3), nn.BatchNorm1d(o), nn.GELU(), nn.MaxPool1d(2))
        self.net = nn.Sequential(block(nch, 16), block(16, 32), block(32, 64), block(64, 96), nn.AdaptiveAvgPool1d(8), nn.Flatten(), nn.Linear(96 * 8, 128), nn.GELU(), nn.Dropout(0.2), nn.Linear(128, 1))
        nv = max(2000, len(X) // 10); idx = torch.randperm(len(X)); va, tr = idx[:nv], idx[nv:]
        opt = torch.optim.AdamW(self.net.parameters(), lr=2e-3, weight_decay=1e-4); sched = torch.optim.lr_scheduler.OneCycleLR(opt, 2e-3, total_steps=self.epochs * ((len(tr) + 255) // 256))
        best = (np.inf, None); bad = 0
        for ep in range(self.epochs):
            self.net.train(); perm = tr[torch.randperm(len(tr))]
            for i in range(0, len(perm), 256):
                b = perm[i:i + 256]; opt.zero_grad(); loss = nn.functional.huber_loss(self.net(X[b]), Y[b]); loss.backward(); opt.step(); sched.step()
            self.net.eval()
            with torch.no_grad():
                vl = float(np.mean([float(nn.functional.mse_loss(self.net(X[va[i:i + 1024]]), Y[va[i:i + 1024]])) for i in range(0, nv, 1024)]))
            if vl < best[0] - 1e-5:
                best = (vl, {k: v.clone() for k, v in self.net.state_dict().items()}); bad = 0
            else:
                bad += 1
                if bad >= 6:
                    break
        self.net.load_state_dict(best[1]); return self
    def predict(self, X):
        import torch
        self.net.eval(); X = torch.tensor(X)
        with torch.no_grad():
            return np.concatenate([self.net(X[i:i + 2048]).numpy()[:, 0] for i in range(0, len(X), 2048)]) * self.sy + self.my


def flat(X): return X.reshape(len(X), -1)
models = {
    'ridge': lambda: ('flat', make_pipeline(StandardScaler(), Ridge(alpha=1000.0))),
    'pca_ridge': lambda: ('flat', make_pipeline(StandardScaler(), PCA(a.npca, random_state=0), Ridge(alpha=1.0))),
    'pca_hgb': lambda: ('flat', make_pipeline(StandardScaler(), PCA(a.npca, random_state=0), HistGradientBoostingRegressor(max_iter=500, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=50, early_stopping=True, random_state=0))),
    'pca_mlp': lambda: ('flat', make_pipeline(StandardScaler(), PCA(a.npca, random_state=0), StandardScaler(), MLPRegressor(hidden_layer_sizes=(128, 128), early_stopping=True, max_iter=300, random_state=0))),
    'cnn': lambda: ('cnn', CNN(a.epochs)),
}
out = f'{REPO}/csv_files/qso_ml_pix_z{a.zlo}-{a.zhi}_{a.tag}.csv'; rows = []   # own file: 26_specphot_ml.py rewrites the bench csv from a start-of-run snapshot
for variant in a.variants.split(','):
    X = inputs(variant)
    for mname in a.models.split(','):
        t0 = time.time(); pred = np.empty_like(M)
        for f in range(a.nfold):
            kind, m = models[mname](); tr = fold != f; Xf = X if kind == 'cnn' else flat(X)
            m.fit(Xf[tr], M[tr]); pred[~tr] = m.predict(Xf[~tr])
        r = M - pred; r2 = 1 - r.var() / M.var(); rms = r.std(); rb = r[bright].std(); m68 = np.percentile(np.abs(r), 68)
        rows.append(dict(slice=f'{a.zlo}-{a.zhi}', method=f'pix: {mname}', features=variant, nfeat=X.shape[1] * npix, target='M', n=n, sigma=M.std(), r2=r2, rms=rms, rms_bright=rb, mad68=m68, seconds=time.time() - t0, error=''))
        print(f'{variant:10s} {mname:10s} R2 {r2:.3f} rms {rms:.3f} bright {rb:.3f} |r|68 {m68:.3f} {time.time() - t0:6.0f}s', flush=True)
        prev = pd.read_csv(out) if os.path.exists(out) else pd.DataFrame()
        prev = prev[~((prev.get('method', pd.Series(dtype=str)) == f'pix: {mname}') & (prev.get('features', pd.Series(dtype=str)) == variant))] if len(prev) else prev
        pd.concat([prev, pd.DataFrame(rows[-1:])], ignore_index=True).to_csv(out, index=False)
print('done')
