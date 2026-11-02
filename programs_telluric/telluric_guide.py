"""The external atmosphere as a guide (Nao 2026-11-01): DECam standard atmosphere (MODTRAN, 5 A) from the LSST stack, on a given wavelength grid.
It says where telluric bands can be, and it is the shape the eigenspectra are compared to; the bases themselves come from the data."""
import numpy as np
ATM = '/opt/miniconda3/envs/lsst-scipipe-12.3.0/share/eups/Darwin/obs_decam_data/g9eeeb5e174/decam/transmission_atmosphere/2012-09-12T00:00:00.ecsv'
TAU_MIN = 0.008                                       # band edge: optical depth at airmass 1
O2 = [(6240, 6340), (6840, 7010), (7570, 7725)]       # gamma, B, A; everything else above TAU_MIN red of 5600 A is H2O


def tau(wave):
    t = np.loadtxt(ATM, comments='#', skiprows=16); order = np.argsort(t[:, 0])
    return -np.log(np.clip(np.interp(wave, t[order, 0], t[order, 1]), 1e-3, 1))


def smooth_continuum_removed(wave, step=150.0):
    """tau minus its slowly varying part (Rayleigh, aerosol, ozone): the lower envelope, linear between the minima of step-A bins."""
    g = np.arange(3500.0, 10500.0, 1.0); t = tau(g); xs, ys = [], []
    for lo in np.arange(g[0], g[-1], step):
        m = (g >= lo) & (g < lo + step); i = np.argmin(t[m]); xs.append(g[m][i]); ys.append(t[m][i])
    return tau(wave) - np.interp(wave, xs, ys)


def band_mask(wave, tau_min=TAU_MIN, grow=8.0):
    """True where the external model has a molecular band above tau_min, grown by +-grow A; red of 5600 A only."""
    t = smooth_continuum_removed(wave); m = (t > tau_min) & (wave > 5600); idx = np.where(m)[0]; out = np.zeros(len(wave), bool)
    for i in idx: out |= np.abs(wave - wave[i]) < grow
    return out


def species(wave):
    """+1 O2, +2 H2O for the band pixels."""
    b = band_mask(wave); s = np.where(b, 2, 0)
    for lo, hi in O2: s[(wave > lo) & (wave < hi) & b] = 1
    return s
