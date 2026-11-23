"""Regenerate data/*_w2.csv from the raw COMSOL exports.   python3 make_w2_csv.py"""
from pathlib import Path
from comsol_io import read_export, write_csv
H = Path(__file__).parent / 'data'
for sub in ('raw', 'diagnostic'):
    for f in sorted((H / sub).glob('*_raw.txt')):
        xs, ws = read_export(f)
        out = (H if sub == 'raw' else H / sub) / f.name.replace('_raw.txt', '_w2.csv')
        write_csv(out, xs, ws); print(out.relative_to(H.parent), len(xs), 'rows')
