#!/usr/bin/env python3
"""
02b_merge_shards.py — igmdust: merge sharded 02_paircount outputs.

02_paircount.py --nshard N --shard i splits foregrounds AND randoms into
N interleaved subsets while keeping one region indexing built from the
full samples, so the per-(color, region, bin) accumulations are directly
addable.  This script sums them, recomputes delta/sigma (jackknife) and
E(B-V) (diagonal solve) with the 02 functions, and writes the merged npz
under the unsharded name.

    python3 02b_merge_shards.py --name bgs_mqcsfd --nshard 6
reads  paircount_<name>_sN.i_qso.npz  (i = 0..N-1)
writes paircount_<name>_qso.npz
"""
import argparse
import importlib.util
import json
import os
import time

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "paircount02", os.path.join(_here, "02_paircount.py"))
pc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pc)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--name", required=True, help="e.g. bgs_mqcsfd")
    ap.add_argument("--nshard", type=int, required=True)
    ap.add_argument("--out", default=pc.OUT_DEFAULT)
    a = ap.parse_args()

    acc = {}
    edges = rp_mid = colors = None
    meta = {}
    for i in range(a.nshard):
        f = f"{a.out}/paircount_{a.name}_s{a.nshard}.{i}_qso.npz"
        d = np.load(f, allow_pickle=True)
        for k in ("jk_sw_d", "jk_swc_d", "jk_sw_r", "jk_swc_r",
                  "npair_data", "npair_ran"):
            acc[k] = acc.get(k, 0) + d[k]
        edges, rp_mid, colors = d["rp_edges"], d["rp_mid"], d["colors"]
        meta = json.loads(str(d["meta"]))
        print(f"  + shard {i}: {f}")

    dat = (acc["jk_sw_d"], acc["jk_swc_d"], None)
    rnd = (acc["jk_sw_r"], acc["jk_swc_r"], None)
    delta, sig = pc.jackknife(dat, rnd)
    ebv, ebv_err = pc.solve_ebv(delta, sig)

    print(f"\n== merged: {a.name} ({a.nshard} shards)")
    for i, rp in enumerate(rp_mid):
        print(f"  {rp*1e3:8.0f} kpc  E(B-V) = {ebv[i]:+9.2e} "
              f"+- {ebv_err[i]:8.1e}")

    meta.update(nshard=a.nshard, merged=time.strftime("%F %T"))
    fo = f"{a.out}/paircount_{a.name}_qso.npz"
    np.savez(fo, rp_edges=edges, rp_mid=rp_mid, delta=delta, sigma=sig,
             ebv=ebv, ebv_err=ebv_err,
             npair_data=acc["npair_data"], npair_ran=acc["npair_ran"],
             colors=colors,
             jk_sw_d=acc["jk_sw_d"], jk_swc_d=acc["jk_swc_d"],
             jk_sw_r=acc["jk_sw_r"], jk_swc_r=acc["jk_swc_r"],
             meta=json.dumps(meta))
    print(f"saved -> {fo}")


if __name__ == "__main__":
    main()
