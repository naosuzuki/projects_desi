"""Benchmark of regressors that predict the brightness M_AB(3025) from the eigencoefficients (Nao 2026-11-12: 'try every
method you can think of and find the best ones'). Redshift slice of the non-BAL table (PCORR_MODE 0-2, NORM > 0, fit).
Feature sets:  D  the direction c/|c| of the coefficient vector only (round 2: the amplitude |c| carries the noise level);
               A2 c_1..c_9 plus the general-basis coefficients g_1..g_9 / g_0 (round 2: coefficients of both bases, nothing else);
               A  c_1..c_9 of the non-BAL basis (the shape; nonlinear methods can still read the noise level from the
                  amplitude of the coefficient vector, which grows toward the faint end);
               S  A + other shape measures: general-basis shape g_k/g_0 (k=1..9), BAL AI of C IV / Si IV / Mg II / Al III,
                  rest coverage, survey one-hot;
               N  noise proxies only: chi2/npix of the fit, clipped fraction, |c| = sqrt(sum c_k^2) (no shape: at fixed z
                  the noise level measures the apparent flux, which is what we try to predict - the control);
               B  S + N (everything the fit records);
               C  B + the photometric-correction record a_0, a_1, dm_g, dm_r, dm_i/z (variability and calibration).
Targets:  M     = M_AB(3025) of NORM (imaging-epoch brightness);
          M_raw = M - a_0, the uncorrected, spectroscopic-epoch brightness (same epoch as the shape).
All features are winsorized at the 0.2 / 99.8 % quantiles of the training fold before scaling (the coefficient tails
reach 100 sigma and wreck polynomial and linear extrapolation on held-out quasars); tree methods get the same input.
Methods (5-fold cross-validation, same folds for all; metrics on the held-out quasars: R^2, rms, rms of the bright half):
  OLS, ridge, quadratic and cubic polynomial ridge, kNN, random forest, extra trees, histogram gradient boosting (default
  and tuned), Nystroem RBF kernel ridge, SVR (RBF, trained on a 15k subsample per fold), Gaussian process (3k subsample),
  scikit-learn MLP, and a torch MLP (3 x 256, dropout, Adam, early stopping).
Outputs csv_files/qso_ml_bench_z<zlo>-<zhi>_<tag>.csv (one row per method x features x target) and a log of the run."""
import argparse, os, sys, time, warnings
import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.base import BaseEstimator, TransformerMixin


class Winsor(BaseEstimator, TransformerMixin):
    """clip every feature to its [q, 1-q] quantiles of the training fold (the coefficient tails reach 100 sigma and wreck polynomial and linear extrapolation)."""
    def __init__(self, q=0.002): self.q = q
    def fit(self, X, y=None): self.lo_, self.hi_ = np.quantile(X, self.q, axis=0), np.quantile(X, 1 - self.q, axis=0); return self
    def transform(self, X): return np.clip(X, self.lo_, self.hi_)
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.kernel_approximation import Nystroem
from sklearn.svm import SVR
from sklearn.gaussian_process import GaussianProcessRegressor, kernels as GK
from sklearn.neural_network import MLPRegressor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import qsolib
warnings.filterwarnings('ignore')

PROD = os.path.expanduser('~/data/sdssdesi_products/quasar'); REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(); ap.add_argument('--zlo', type=float, default=1.65); ap.add_argument('--zhi', type=float, default=1.70); ap.add_argument('--tag', default='v1')
ap.add_argument('--nfold', type=int, default=5); ap.add_argument('--only', default='', help='comma list of method names to run'); ap.add_argument('--threads', type=int, default=6)
ap.add_argument('--fsets', default='A,S,N,B,C', help='feature sets to run'); ap.add_argument('--save-pred', action='store_true', help='write the held-out predictions to csv_files/qso_ml_pred_<method>_<fset>_<target>.npz')
ap.add_argument('--train-zlo', type=float, default=None); ap.add_argument('--train-zhi', type=float, default=None, help='round 3: train on this wider redshift range (z becomes a feature); metrics on the held-out quasars inside --zlo..--zhi'); ap.add_argument('--targets', default='M,M_raw', help='M, M_raw, M_rel (round 6: apparent magnitude minus the median of its 0.05 redshift bin, cosmology-free), M_rel_raw'); ap.add_argument('--no-z', action='store_true', help='round 6: keep redshift out of the inputs even with --train-zlo'); ap.add_argument('--extra-npz', default='', help='round 8: qso_color_*.npz of 31_specphot_color.py merged on TARGETID (DESI slice only); adds feature sets X (colour-free c1..c9) and XC (+ colour)'); ap.add_argument('--out', default='', help='csv to write instead of csv_files/qso_ml_bench_z<slice>_<tag>.csv (a concurrent run must not share the file)'); ap.add_argument('--zbin', type=float, default=0.05, help='bin width for M_rel and for the per-bin rms table'); ap.add_argument('--round', default='1', help='label written to the csv; rows of the same round are replaced'); a = ap.parse_args()
os.environ.setdefault('OMP_NUM_THREADS', str(a.threads))

# ---- data
cols = ['TARGETID', 'SURVEY', 'Z', 'NORM', 'PCORR_MODE', 'PCORR_A0', 'PCORR_A1', 'PCORR_DM_G', 'PCORR_DM_R', 'PCORR_DM_I', 'PCORR_DM_Z', 'CHI2', 'NPIX', 'NCLIP', 'REST_LAMMIN', 'REST_LAMMAX',
        'FOREST_TRANS', 'BAL_AI_CIV', 'BAL_AI_SIIV', 'BAL_AI_MGII', 'BAL_AI_ALIII'] + [f'C{k}' for k in range(1, 10)] + [f'G{k}' for k in range(0, 10)]
d = pd.read_parquet(f'{PROD}/specphot/qso_specphot_coeff_{a.tag}.parquet', columns=cols)
tzlo, tzhi = (a.train_zlo, a.train_zhi) if a.train_zlo is not None else (a.zlo, a.zhi)
s = d[(d.Z >= tzlo) & (d.Z < tzhi) & (d.NORM > 0) & d.PCORR_MODE.isin([0, 1, 2])].copy()
s = s[np.isfinite(s[[f'C{k}' for k in range(1, 10)]]).all(axis=1) & (s.G0 > 0)].reset_index(drop=True)
if a.extra_npz:                                                          # round 8: colour-free coefficients and colour from the pixel refit
    x = np.load(a.extra_npz); xd = pd.DataFrame({'TARGETID': x['TARGETID'].astype('int64'), 'COLOR': x['COLOR'], **{f'X{k}': x[f'CC{k}'] for k in range(1, 10)}})
    xd = xd[np.isfinite(xd.drop(columns='TARGETID')).all(axis=1)]
    s = s[s.SURVEY == 'desi'].copy(); s['TARGETID'] = s.TARGETID.astype('int64'); s = s.merge(xd, on='TARGETID', how='inner').reset_index(drop=True)
    print(f'extra npz: {len(s):,} DESI quasars matched', flush=True)
K = -48.6 - 2.5 * np.log10(3025.0**2 / 2.998e18); s['M'] = -2.5 * (np.log10(s.NORM) + 0.4 * qsolib.passB.COSMO.distmod(s.Z.to_numpy()).value) + K
s['M_raw'] = s.M - s.PCORR_A0
s['m_app'] = -2.5 * np.log10(s.NORM) + K; s['ZBIN'] = np.floor(s.Z / a.zbin).astype(int)
s['M_rel'] = s.m_app - s.groupby('ZBIN').m_app.transform('median')          # round 6: relative brightness inside the z bin; no distance modulus enters
s['M_rel_raw'] = s.M_rel - s.PCORR_A0
s['DM_RED'] = s.PCORR_DM_Z.fillna(s.PCORR_DM_I)
FA = [f'C{k}' for k in range(1, 10)]
for k in range(1, 10):
    s[f'g{k}'] = s[f'G{k}'] / s.G0
s['CHI2N'] = s.CHI2 / s.NPIX; s['CLIPF'] = s.NCLIP / s.NPIX
for sv in ('desi', 'dr17', 'dr8'):
    s[f'S_{sv}'] = (s.SURVEY == sv).astype(float)
s['ABSC'] = np.sqrt((s[FA]**2).sum(axis=1))
FS = FA + [f'g{k}' for k in range(1, 10)] + ['BAL_AI_CIV', 'BAL_AI_SIIV', 'BAL_AI_MGII', 'BAL_AI_ALIII', 'REST_LAMMIN', 'REST_LAMMAX', 'S_desi', 'S_dr17', 'S_dr8']
FN = ['CHI2N', 'CLIPF', 'ABSC']
FB = FS + FN
for k in range(1, 10):
    s[f'u{k}'] = s[f'C{k}'] / s.ABSC                                    # direction of the coefficient vector (amplitude removed)
FD = [f'u{k}' for k in range(1, 10)]                                     # D: direction only - the amplitude carries the noise level
FA2 = FA + [f'g{k}' for k in range(1, 10)]                              # A2: coefficients of both bases (non-BAL c_k and general g_k/g_0)

FC = FB + ['PCORR_A0', 'PCORR_A1', 'PCORR_DM_G', 'PCORR_DM_R', 'DM_RED']
FSETS = {'A': FA, 'S': FS, 'N': FN, 'B': FB, 'C': FC, 'D': FD, 'A2': FA2}
if a.extra_npz:
    FSETS['X'] = [f'X{k}' for k in range(1, 10)]; FSETS['XC'] = FSETS['X'] + ['COLOR']; FSETS['AC'] = FA + ['COLOR']
if a.train_zlo is not None and not a.no_z:                               # round 3: redshift is an input, the slice is the evaluation window
    FSETS = {k: v + ['Z'] for k, v in FSETS.items()}
for c in sorted(set(FC + FD + FA2)):
    s[c] = s[c].astype(float).replace([np.inf, -np.inf], np.nan).fillna(0.0)
insl = ((s.Z >= a.zlo) & (s.Z < a.zhi)).to_numpy()
print(f'train {tzlo}-{tzhi}: {len(s):,} quasars, evaluation slice {a.zlo}-{a.zhi}: {insl.sum():,}; sigma(M) {s.M[insl].std():.3f}, sigma(M_raw) {s.M_raw[insl].std():.3f}', flush=True)
rng = np.random.default_rng(1); fold = rng.integers(0, a.nfold, len(s))
bright = (s.M < np.median(s.M[insl])).to_numpy()


class TorchMLP:
    def __init__(self, width=256, depth=3, dropout=0.1, lr=1e-3, epochs=60, batch=512, patience=8, seed=0):
        self.p = dict(width=width, depth=depth, dropout=dropout, lr=lr, epochs=epochs, batch=batch, patience=patience, seed=seed)
    def fit(self, X, y):
        import torch, torch.nn as nn
        p = self.p; torch.manual_seed(p['seed']); torch.set_num_threads(a.threads)
        self.sx = make_pipeline(Winsor(), StandardScaler()).fit(X); self.my, self.sy = y.mean(), y.std()
        X = torch.tensor(self.sx.transform(X), dtype=torch.float32); Y = torch.tensor((y - self.my) / self.sy, dtype=torch.float32)[:, None]
        nv = max(2000, len(X) // 10); idx = torch.randperm(len(X)); va, tr = idx[:nv], idx[nv:]
        layers = []; nin = X.shape[1]
        for _ in range(p['depth']):
            layers += [nn.Linear(nin, p['width']), nn.GELU(), nn.Dropout(p['dropout'])]; nin = p['width']
        self.net = nn.Sequential(*layers, nn.Linear(nin, 1)); opt = torch.optim.AdamW(self.net.parameters(), lr=p['lr'], weight_decay=1e-4)
        sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, p['epochs']); best = (np.inf, None); bad = 0
        for ep in range(p['epochs']):
            self.net.train(); perm = tr[torch.randperm(len(tr))]
            for i in range(0, len(perm), p['batch']):
                b = perm[i:i + p['batch']]; opt.zero_grad(); loss = nn.functional.huber_loss(self.net(X[b]), Y[b], delta=1.0); loss.backward(); opt.step()
            sched.step(); self.net.eval()
            with torch.no_grad():
                vl = float(nn.functional.mse_loss(self.net(X[va]), Y[va]))
            if vl < best[0] - 1e-5:
                best = (vl, {k: v.clone() for k, v in self.net.state_dict().items()}); bad = 0
            else:
                bad += 1
                if bad >= p['patience']:
                    break
        self.net.load_state_dict(best[1]); return self
    def predict(self, X):
        import torch
        self.net.eval()
        with torch.no_grad():
            return self.net(torch.tensor(self.sx.transform(X), dtype=torch.float32)).numpy()[:, 0] * self.sy + self.my


class Sub:
    """fit on a random subsample of the training fold (for O(n^2) methods)."""
    def __init__(self, model, n): self.model, self.n = model, n
    def fit(self, X, y):
        i = np.random.default_rng(0).choice(len(X), min(self.n, len(X)), replace=False); self.model.fit(X[i], y[i]); return self
    def predict(self, X): return self.model.predict(X)


def methods():
    sc = lambda: make_pipeline(Winsor(), StandardScaler())
    return {
        'OLS': lambda: make_pipeline(sc(), LinearRegression()),
        'ridge': lambda: make_pipeline(sc(), Ridge(alpha=1.0)),
        'poly2 ridge': lambda: make_pipeline(sc(), PolynomialFeatures(2), StandardScaler(), Ridge(alpha=10.0)),
        'poly3 ridge': lambda: make_pipeline(sc(), PolynomialFeatures(3), StandardScaler(), Ridge(alpha=100.0)),
        'kNN k=50': lambda: make_pipeline(sc(), KNeighborsRegressor(50, weights='distance', n_jobs=a.threads)),
        'random forest': lambda: RandomForestRegressor(300, min_samples_leaf=20, max_features=0.5, n_jobs=a.threads, random_state=0),
        'extra trees': lambda: ExtraTreesRegressor(300, min_samples_leaf=20, max_features=0.5, n_jobs=a.threads, random_state=0),
        'HGB default': lambda: HistGradientBoostingRegressor(random_state=0),
        'HGB tuned': lambda: HistGradientBoostingRegressor(max_iter=800, learning_rate=0.03, max_leaf_nodes=63, min_samples_leaf=50, l2_regularization=1.0, early_stopping=True, validation_fraction=0.1, random_state=0),
        'Nystroem RBF ridge': lambda: make_pipeline(sc(), Nystroem(gamma=0.1, n_components=1000, random_state=0), Ridge(alpha=1.0)),
        'SVR RBF (15k)': lambda: Sub(make_pipeline(sc(), SVR(C=3.0, epsilon=0.1, gamma='scale')), 15000),
        'GP Matern (3k)': lambda: Sub(make_pipeline(sc(), GaussianProcessRegressor(GK.ConstantKernel() * GK.Matern(nu=1.5, length_scale=np.ones(1)) + GK.WhiteKernel(), normalize_y=True, n_restarts_optimizer=0)), 3000),
        'torch MLP': lambda: TorchMLP(),
        'torch MLP wide': lambda: TorchMLP(width=512, depth=4, dropout=0.15, epochs=120, patience=12),
        'torch MLP ens5': lambda: Ens([TorchMLP(width=256, depth=3, epochs=80, patience=10, seed=i) for i in range(5)]),
        'HGB search': lambda: HGBSearch(),
        'torch MLP deep': lambda: TorchMLP(width=256, depth=6, dropout=0.05, epochs=120, patience=12),
        'torch MLP search': lambda: MLPSearch(),
        'torch MLP ens10 long': lambda: Ens([TorchMLP(width=384, depth=4, dropout=0.1, lr=7e-4, epochs=150, patience=15, seed=i) for i in range(10)]),
        'stack MLPsearch+HGBsearch': lambda: Ens([MLPSearch(), HGBSearch()]),
        'stack HGB+MLP+kNN': lambda: Ens([HistGradientBoostingRegressor(max_iter=800, learning_rate=0.03, max_leaf_nodes=63, min_samples_leaf=50, l2_regularization=1.0, early_stopping=True, random_state=0),
                                          TorchMLP(width=256, depth=3, epochs=80, patience=10), make_pipeline(sc(), KNeighborsRegressor(50, weights='distance', n_jobs=a.threads))]),
    }


class Ens:
    """mean of the predictions of several fitted models."""
    def __init__(self, models): self.models = models
    def fit(self, X, y):
        for m in self.models: m.fit(X, y)
        return self
    def predict(self, X): return np.mean([m.predict(X) for m in self.models], axis=0)


class HGBSearch:
    """small grid over learning rate / leaves / leaf size, chosen on a 15 % validation split of the training fold."""
    def fit(self, X, y):
        rng = np.random.default_rng(0); va = rng.random(len(X)) < 0.15; best = (np.inf, None)
        for lr in (0.02, 0.05, 0.1):
            for leaves in (15, 31, 63, 127):
                for leaf in (20, 50, 200):
                    m = HistGradientBoostingRegressor(max_iter=1500, learning_rate=lr, max_leaf_nodes=leaves, min_samples_leaf=leaf, l2_regularization=1.0, early_stopping=True, validation_fraction=0.1, random_state=0).fit(X[~va], y[~va])
                    err = np.mean((m.predict(X[va]) - y[va])**2)
                    if err < best[0]: best = (err, (lr, leaves, leaf))
        lr, leaves, leaf = best[1]; self.best = best[1]
        self.m = HistGradientBoostingRegressor(max_iter=1500, learning_rate=lr, max_leaf_nodes=leaves, min_samples_leaf=leaf, l2_regularization=1.0, early_stopping=True, validation_fraction=0.1, random_state=0).fit(X, y); return self
    def predict(self, X): return self.m.predict(X)


class MLPSearch:
    """grid over width / depth / dropout / learning rate on a 15 % validation split of the training fold, then a 3-seed ensemble refit."""
    def fit(self, X, y):
        rng = np.random.default_rng(0); va = rng.random(len(X)) < 0.15; best = (np.inf, None)
        for width in (256, 512):
            for depth in (3, 5):
                for drop in (0.05, 0.2):
                    for lr in (5e-4, 2e-3):
                        m = TorchMLP(width=width, depth=depth, dropout=drop, lr=lr, epochs=50, patience=8).fit(X[~va], y[~va])
                        err = np.mean((m.predict(X[va]) - y[va])**2)
                        if err < best[0]: best = (err, (width, depth, drop, lr))
        width, depth, drop, lr = best[1]; self.best = best[1]; print('   MLP search best', best[1], flush=True)
        self.m = Ens([TorchMLP(width=width, depth=depth, dropout=drop, lr=lr, epochs=120, patience=12, seed=i) for i in range(3)]).fit(X, y); return self
    def predict(self, X): return self.m.predict(X)


def cv(make, X, y):
    pred = np.empty_like(y)
    for f in range(a.nfold):
        tr = fold != f; m = make().fit(X[tr], y[tr]); pred[~tr] = m.predict(X[~tr])
    r = (y - pred)[insl]; yy = y[insl]; b = bright[insl]
    cv.last = (yy, pred[insl]); cv.full = pred; return 1 - r.var() / yy.var(), r.std(), r[b].std(), np.percentile(np.abs(r), 68)


rows = []; want = [w.strip() for w in a.only.split(',') if w.strip()]; fsets = a.fsets.split(','); targets = a.targets.split(',')
outcsv = a.out or f'{REPO}/csv_files/qso_ml_bench_z{a.zlo}-{a.zhi}_{a.tag}.csv'


def write_rows():
    """re-read at write time (not a start-of-run snapshot), replace this round's rows, append ours."""
    prev = pd.read_csv(outcsv) if os.path.exists(outcsv) else pd.DataFrame()
    if len(prev) and 'round' in prev:
        prev = prev[prev['round'].astype(str) != a.round]
    pd.concat([prev, pd.DataFrame(rows)], ignore_index=True).to_csv(outcsv, index=False)


for mname, make in methods().items():
    if want and mname not in want:
        continue
    for fname in fsets:
        feats = FSETS[fname]
        if mname in ('GP Matern (3k)',) and fname != 'A':
            continue                                               # isotropic GP on 9 features only
        for tname in targets:
            if a.round == '1' and tname == 'M_raw' and (fname not in ('A', 'S') or mname not in ('OLS', 'poly2 ridge', 'HGB tuned', 'torch MLP')):
                continue
            t0 = time.time(); X = s[feats].to_numpy(); y = s[tname].to_numpy()
            try:
                r2, rms, rms_b, mad68 = cv(make, X, y); err = ''
            except Exception as e:
                r2 = rms = rms_b = mad68 = np.nan; err = str(e)[:60]
            rows.append(dict(round=a.round, slice=f'{a.zlo}-{a.zhi}', train=f'{tzlo}-{tzhi}', method=mname, features=fname, nfeat=len(feats), target=tname, n=int(insl.sum()), ntrain=len(s), sigma=y[insl].std(), r2=r2, rms=rms, rms_bright=rms_b, mad68=mad68, seconds=time.time() - t0, error=err))
            if a.save_pred and not err and a.train_zlo is not None:        # held-out rms per redshift bin over the whole training range
                rr = pd.DataFrame(dict(ZBIN=s.ZBIN, r=y - cv.full, bright=bright)); g = rr.groupby('ZBIN')
                zb = pd.DataFrame(dict(zlo=g.r.count().index * a.zbin, n=g.r.count(), sigma=g.apply(lambda q: y[q.index].std()), rms=g.r.std(), rms_bright=g.apply(lambda q: q.r[q.bright].std()), mad68=g.r.apply(lambda q: np.percentile(np.abs(q), 68))))
                zb.to_csv(f"{REPO}/csv_files/qso_ml_zbins_{mname.replace(' ', '_').replace('+', '_')}_{fname}_{tname}_r{a.round}.csv", index=False); print(zb.round(3).to_string(), flush=True)
            if a.save_pred and not err:
                np.savez(f"{REPO}/csv_files/qso_ml_pred_{mname.replace(' ', '_').replace('+', '')}_{fname}_{tname}_r{a.round}.npz", M=cv.last[0], pred=cv.last[1])
            print(f'{mname:20s} {fname:2s} {tname:5s} R2 {r2:.3f} rms {rms:.3f} bright {rms_b:.3f} |r|68 {mad68:.3f}  {time.time() - t0:6.0f}s {err}', flush=True)
            write_rows()
print('done')
