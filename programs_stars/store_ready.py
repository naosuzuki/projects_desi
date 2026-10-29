#!/usr/bin/env python3
"""exit 0 when every coadd among the first N lines of star_priority_files.txt that holds step-6 stars is in star_hook's store."""
import os, sys
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import star_hook
n = int(sys.argv[1])
pr = [l.strip() for l in open(os.path.expanduser('~/data/sdssdesi_products/stars/star_priority_files.txt')) if l.strip()][:n]
need = set(pd.read_parquet(star_hook.SEL300, columns=['FILE']).FILE)
have = {e.name[:-4] + '.fits' for e in os.scandir(star_hook.STORE) if e.name.endswith('.npz') and '.tmp' not in e.name}
# coadds the readers cannot find on disk (logged MISSING by the quasar shards / fill_store; 44 coadds, 48 stars) never arrive
S = os.path.expanduser('~/data/sdssdesi_status'); gone = set()
for lg in [f'{S}/qso_stack_desi{i}.log' for i in range(3)] + [f'{S}/fill_store.log']:
    if os.path.exists(lg):
        gone |= {l.split()[1] for l in open(lg, errors='replace') if l.startswith('MISSING ')}
missing = [f for f in pr if f in need and f not in have and f not in gone]
print(f'{len(missing):,} of the first {n:,} priority coadds still missing')
sys.exit(0 if not missing else 1)
