#!/usr/bin/env python3
"""Select non-failing candidates using the historical ranking rule."""
import argparse
import csv
import gzip
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--report', type=Path, required=True)
parser.add_argument('--target', type=int, default=5200)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
if args.target < 1:
    parser.error('--target must be positive')
opener = gzip.open if args.report.suffix == '.gz' else open
with opener(args.report, 'rt', newline='') as handle:
    rows = [r for r in csv.DictReader(handle) if r['classification'] in {'PASS', 'BORDERLINE'}]
rows.sort(key=lambda r: (r['classification'] != 'PASS', float(r['risk_score'])))
if len(rows) < args.target:
    raise SystemExit('Not enough non-failing candidates')
args.out.parent.mkdir(parents=True, exist_ok=True)
with args.out.open('w') as handle:
    for i, row in enumerate(rows[:args.target], 1):
        handle.write(f">oligo_{i:04d}\n{row['payload_sequence']}\n")
print(f'Wrote {args.target} sequences to {args.out}')
