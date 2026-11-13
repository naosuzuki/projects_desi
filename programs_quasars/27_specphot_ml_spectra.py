"""Rest-frame spectra of a redshift slice for the pixel-level brightness predictors (28_specphot_ml_cnn.py). DESI members
of the slice (non-BAL table, PCORR_MODE 0-2, fit) are re-read from the native-pixel store, corrected with the recorded
a_0, a_1 (same rule as 21_specphot_fit.correction: x = log10(lam/6000) clipped to the covered bands' effective-wavelength
range), dereddened (F99, 0.86 SFD), put on the quasar rest lattice as F_lam(1+z), cut to --lmin..--lmax A and binned by
--bin pixels (inverse-variance weighted). Stored: FLUX/IVAR [n, npix] float32 in units of the median flux in 2950-3100 A
(a measured, basis-free normalization; NORM of the fit is kept alongside), the target M_AB(3025) and the table columns.
Output /Volumes/exdisk3/sdssdesi_ml/qso_slice_z<zlo>-<zhi>_<tag>.npz (exdisk3 may be written; the internal disk is full)."""
import argparse, glob, importlib.util, os, struct, sys, time, zipfile
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import qsolib as Q
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'programs_templates')); import restframe as RF
spec = importlib.util.spec_from_file_location('sp', os.path.join(os.path.dirname(os.path.abspath(__file__)), '21_specphot_fit.py')); SP = importlib.util.module_from_spec(spec); spec.loader.exec_module(SP)
passB = Q.passB

ap = argparse.ArgumentParser(); ap.add_argument('--zlo', type=float, default=1.65); ap.add_argument('--zhi', type=float, default=1.70); ap.add_argument('--tag', default='v1')
ap.add_argument('--lmin', type=float, default=1350.0); ap.add_argument('--lmax', type=float, default=3650.0); ap.add_argument('--bin', type=int, default=4)
ap.add_argument('--out', default='/Volumes/exdisk3/sdssdesi_ml'); ap.add_argument('--nfiles', type=int, default=0); a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
cols = ['TARGETID', 'FILE', 'SURVEY', 'Z', 'EBV', 'PHOTSYS', 'NORM', 'PCORR_MODE', 'PCORR_A0', 'PCORR_A1', 'PCORR_NBAND', 'PCORR_DM_G', 'PCORR_DM_R', 'PCORR_DM_Z', 'MORPHTYPE', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'CHI2', 'NPIX'] + [f'C{k}' for k in range(1, 10)]
d = pd.read_parquet(f'{SP.OUT}/qso_specphot_coeff_{a.tag}.parquet', columns=cols)
s = d[(d.SURVEY == 'desi') & (d.Z >= a.zlo) & (d.Z < a.zhi) & (d.NORM > 0) & d.PCORR_MODE.isin([0, 1, 2])].copy()
K = -48.6 - 2.5 * np.log10(3025.0**2 / 2.998e18); s['MAB'] = -2.5 * (np.log10(s.NORM) + 0.4 * passB.COSMO.distmod(s.Z.to_numpy()).value) + K
files = sorted(s.FILE.unique()); files = files[:a.nfiles] if a.nfiles else files
print(f'{len(s):,} DESI quasars in {len(files)} coadd files', flush=True)
# rest sub-grid
rid = np.arange(Q.QRID_MIN, Q.QRID_MIN + Q.QNREST); sel = (Q.QWAVE >= a.lmin) & (Q.QWAVE <= a.lmax); i0 = np.argmax(sel); n0 = sel.sum() - sel.sum() % a.bin
wrest = Q.QWAVE[i0:i0 + n0].reshape(-1, a.bin); wout = 10**np.log10(wrest).mean(axis=1)
# effective wavelengths used by the correction's x-clipping
def xrange_of(names):
    xs = [np.log10(passB.band_range(passB.filt(nm))[0] * 0 + np.average(passB.filt(nm).wavelength, weights=passB.filt(nm).response) / SP.X0) for nm in names.values()]
    return min(xs), max(xs)
def npz_rows(p, rows):
    """TARGETID, W0, DW and the FLUX / IVAR rows of a native-store npz without reading the whole file: the members are stored
    uncompressed, so each .npy is memory-mapped at its offset inside the zip (the full read costs 2 s per 42 MB file; 15,630 files)."""
    zf = zipfile.ZipFile(p); out = {}
    with open(p, 'rb') as f:
        for zi in zf.infolist():
            f.seek(zi.header_offset); h = f.read(30); nlen, elen = struct.unpack('<HH', h[26:30]); start = zi.header_offset + 30 + nlen + elen
            f.seek(start); ver = np.lib.format.read_magic(f); shape, fortran, dtype = np.lib.format._read_array_header(f, ver); off = f.tell()
            name = zi.filename[:-4]
            if name in ('FLUX', 'IVAR'):
                mm = np.memmap(p, dtype=dtype, mode='r', offset=off, shape=shape); out[name] = np.array(mm[np.sort(rows)]) if rows is not None else mm
            else:
                out[name] = np.fromfile(f, dtype=dtype, count=int(np.prod(shape)) if shape else 1).reshape(shape) if shape else np.fromfile(f, dtype=dtype, count=1)[0]
    return out


byfile = {f: g for f, g in s.groupby('FILE')}; F_all, IV_all, meta = [], [], []; t0 = time.time()
for j, fn in enumerate(files):
    g = byfile[fn]; p = os.path.join(SP.NATIVE, os.path.basename(fn).replace('.fits', '') + '.npz')
    if not os.path.exists(p):
        continue
    hdr = npz_rows(p, None); tid = hdr['TARGETID']; pos = pd.Series(np.arange(len(tid)), index=tid); hit = pos.reindex(g.TARGETID.to_numpy()); ok = hit.notna().to_numpy(); g = g[ok]; rows = hit[ok].to_numpy(int)
    if len(g) == 0:
        continue
    order = np.argsort(rows); g = g.iloc[order]; rows = rows[order]; z = npz_rows(p, rows)
    w = float(hdr['W0']) + float(hdr['DW']) * np.arange(hdr['FLUX'].shape[1]); flux = z['FLUX'].astype(np.float64) * passB.FLUX_UNIT; ivar = z['IVAR'].astype(np.float64) / passB.FLUX_UNIT**2
    good = ivar > 0; ivar = np.where(good, ivar, 0.0)
    # photometric correction from the recorded a0, a1 (tilt clipped to the band range of the quasar's filter system)
    xl = np.log10(w / SP.X0); fac = np.ones_like(flux)
    for south in (True, False):
        m = (g.PHOTSYS.to_numpy() == 'S') == south
        if not m.any():
            continue
        xlo, xhi = xrange_of(SP.filters('desi', south)); xc = np.clip(xl, xlo, xhi)
        a0 = g.PCORR_A0.to_numpy(float)[m][:, None]; a1 = g.PCORR_A1.to_numpy(float)[m][:, None]
        fac[m] = 10.0**(-0.4 * (a0 + np.where(a1 != 0, a1 * xc[None, :], 0.0)))
    flux *= fac; ivar /= fac**2
    A1 = passB.extinction.fitzpatrick99(np.asarray(w, float), 1.0, 3.1); A = (3.1 * 0.86 * np.nan_to_num(g.EBV.to_numpy(float)))[:, None] * A1[None, :]
    flux *= 10**(0.4 * A); ivar *= 10**(-0.8 * A)
    zz = g.Z.to_numpy(float); rf, riv = RF.desi_rest(w, flux, ivar, (~good).astype(np.uint8), zz, Q.QRID_MIN, Q.QNREST)
    F = rf.astype(np.float64) * (1 + zz)[:, None]; IV = riv.astype(np.float64) / (1 + zz)[:, None]**2; F[~np.isfinite(F)] = 0; IV[~np.isfinite(IV)] = 0
    F = F[:, i0:i0 + n0].reshape(len(g), -1, a.bin); IV = IV[:, i0:i0 + n0].reshape(len(g), -1, a.bin)
    with np.errstate(all='ignore'):
        Fb = np.where(IV.sum(2) > 0, (F * IV).sum(2) / IV.sum(2), 0.0); IVb = IV.sum(2)
    win = (wout > 2950) & (wout < 3100); med = np.array([np.median(f[win][iv[win] > 0]) if (iv[win] > 0).sum() > 10 else np.nan for f, iv in zip(Fb, IVb)])
    F_all.append((Fb / med[:, None]).astype(np.float32)); IV_all.append((IVb * med[:, None]**2).astype(np.float32)); gm = g.copy(); gm['MED3025'] = med; meta.append(gm)
    if (j + 1) % 200 == 0:
        print(f'  {j + 1}/{len(files)} files, {sum(len(m) for m in meta):,} quasars, {(time.time() - t0) / 60:.1f} min', flush=True)
meta = pd.concat(meta, ignore_index=True); F = np.concatenate(F_all); IV = np.concatenate(IV_all)
keep = np.isfinite(meta.MED3025.to_numpy()) & (meta.MED3025.to_numpy() > 0); F, IV, meta = F[keep], IV[keep], meta[keep].reset_index(drop=True); F[~np.isfinite(F)] = 0
out = f'{a.out}/qso_slice_z{a.zlo}-{a.zhi}_{a.tag}.npz'
np.savez(out, WAVE=wout.astype(np.float32), FLUX=F, IVAR=IV, **{c: meta[c].to_numpy() for c in ['TARGETID', 'Z', 'NORM', 'MED3025', 'MAB', 'PCORR_A0', 'PCORR_MODE', 'CHI2', 'NPIX'] + [f'C{k}' for k in range(1, 10)]})
print(f'wrote {out}: {F.shape}, {(time.time() - t0) / 60:.1f} min', flush=True)
