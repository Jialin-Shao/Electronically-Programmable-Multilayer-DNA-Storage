#!/usr/bin/env python3
"""Select non-failing candidates using the historical ranking rule."""
import argparse
import csv
import gzip
import glob
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--report', nargs='+', required=True, help='Report files or quoted glob patterns')
parser.add_argument('--target', type=int, default=5200)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
if args.target < 1:
    parser.error('--target must be positive')
reports = []
for pattern in args.report:
    matches = sorted(glob.glob(pattern))
    if not matches:
        parser.error(f'No report files match: {pattern}')
    reports.extend(Path(name) for name in matches)
if len(reports) != len(set(reports)):
    parser.error('Report files must not be supplied more than once')
rows = []
fieldnames = None
for report in reports:
    opener = gzip.open if report.suffix == '.gz' else open
    with opener(report, 'rt', newline='') as handle:
        reader = csv.DictReader(handle)
        if fieldnames is None:
            fieldnames = reader.fieldnames
        elif reader.fieldnames != fieldnames:
            parser.error(f'Inconsistent report columns: {report}')
        rows.extend(r for r in reader if r['classification'] in {'PASS', 'BORDERLINE'})
rows.sort(key=lambda r: (r['classification'] != 'PASS', float(r['risk_score'])))
if len(rows) < args.target:
    raise SystemExit('Not enough non-failing candidates')
args.out.parent.mkdir(parents=True, exist_ok=True)
with args.out.open('w') as handle:
    for i, row in enumerate(rows[:args.target], 1):
        handle.write(f">oligo_{i:04d}\n{row['payload_sequence']}\n")
print(f'Wrote {args.target} sequences to {args.out}')
