#!/usr/bin/env python3
"""Verify archive integrity and the final library without third-party packages."""
import csv
import gzip
import hashlib
from pathlib import Path

root = Path(__file__).resolve().parents[1]
for line in (root / 'SHA256SUMS.txt').read_text().splitlines():
    expected, name = line.split('  ', 1)
    h = hashlib.sha256()
    with (root / name).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    assert h.hexdigest() == expected, f'Checksum mismatch: {name}'

def sequences(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('>')]

base = root / 'data/sequences'
raw = sequences(base / 'candidates_300000.fasta.gz')
survivors = sequences(base / 'survivors_10163_encoded.fasta')
final = sequences(base / 'final_5200_encoded.fasta')
full = sequences(base / 'final_5200_full_length.fasta')
assert len(raw) == len(set(raw)) == 300000
assert len(survivors) == len(set(survivors)) == 10163
assert len(final) == len(set(final)) == len(full) == 5200
assert set(final) <= set(survivors) <= set(raw)
assert all(len(s) == 120 and set(s) <= set('ACGT') for s in raw)
assert all(b == 'TCTGCAAGTAGCCAAGGGTAAGCAAGGATC' + a + 'CAAGGGTTCATAAGTGCTTCGTGCAGATTC' for a,b in zip(final,full))
with gzip.open(root / 'data/screening/screening_all_300000.csv.gz', 'rt') as f:
    rows = [r for r in csv.DictReader(f) if r['classification'] in {'PASS','BORDERLINE'}]
rows.sort(key=lambda r: (r['classification'] != 'PASS', float(r['risk_score'])))
assert [r['payload_sequence'] for r in rows[:5200]] == final
print('PASS: checksums, counts, lengths, uniqueness, membership, flanking regions and final ranking')
