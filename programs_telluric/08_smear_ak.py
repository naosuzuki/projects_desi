#!/usr/bin/env python3
"""Why the DR8 stack means of a_k are not zero (Nao's step 3, 2026-11-02). tell_hook.models first divided a stack by 1 + T0 at rest, but the stack holds T0 smeared by the members' velocities, so each star's model kept T0(rest) - T0(members); fixed 2026-11-03.
This predicts a_k per stack from that old model: the difference fitted with the design of tellapply, no prior and no spectrum weights, so a rough estimate. Compared with the measured stack mean of the DR8 a_k (run on the fixed coefficients, r and slope should fall to 0). --tell DIR reads the coefficients and members from DIR. Output: csv_files/plots_telluric_starcompare_smear.csv, or stdout with --tell
"""
import os, sys, argparse, importlib.util
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path[:0] = [HERE, os.path.join(os.path.dirname(HERE), 'programs_stars')]
import tellapply
sp = importlib.util.spec_from_file_location('compare', os.path.join(HERE, '06_star_compare.py')); C = importlib.util.module_from_spec(sp); sp.loader.exec_module(C)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--tell', default=None); a_ = ap.parse_args(); d = a_.tell or C.TELL
    B8, B17 = tellapply.basis('dr8'), tellapply.basis('dr17'); tag = C.tag('sdss')
    g = pd.read_csv(f'{d}/members_{tag}.csv', usecols=['SID', 'STACK']).merge(pd.read_parquet(f'{d}/tell_coeff_{tag}.parquet'), on='SID'); share = (g.SURVEY == 'dr8').groupby(g.STACK).mean().to_dict(); Ts = C.tell_hook.smeared_t0(g.STACK.to_numpy(), g.SURVEY.to_numpy(), g.LOGZ.to_numpy())
    a = [f'A{k}' for k in range(len(B8.T))]; mean = g[(g.SURVEY == 'dr8') & (g.NPIX > 0)].groupby('STACK')[a].mean(); G = B8.grid(C.WAVE); p, D = G['p'], G['D']; ns = len(B8.sel); rows = []
    for name in mean.index:
        x = np.zeros(len(a)); x[B8.sel] = np.linalg.lstsq(D.T, (share[name] * B8.Ffull[0] + (1 - share[name]) * B17.Ffull[0] - Ts[name])[p], rcond=None)[0][:ns]; rows.append(x)
    P = pd.DataFrame(rows, index=mean.index, columns=a); r = mean - P
    t = pd.DataFrame([dict(T=n, measured_med=mean[k].median(), measured_rms=mean[k].std(), predicted_med=P[k].median(), predicted_rms=P[k].std(), r=np.corrcoef(P[k], mean[k])[0, 1], slope=np.polyfit(P[k], mean[k], 1)[0], left_med=r[k].median(), left_rms=r[k].std()) for k, n in zip(a, B8.labels)])
    if not a_.tell: t.to_csv(os.path.join(C.CSV, 'plots_telluric_starcompare_smear.csv'), index=False)
    print(t.round(4).to_string())


if __name__ == '__main__':
    main()
