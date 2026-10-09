# Data dictionary

- `candidate_id`: normalized ID in original candidate order (`candidate_000001` to `candidate_300000`).
- `oligo_id`: final library ID in increasing-risk selection order (`oligo_0001` to `oligo_5200`).
- `sequence_id`: original screening identifier; see candidate mapping table.
- `classification`: `PASS`, `BORDERLINE` or `FAIL`, as returned by the screening implementation.
- `risk_score`: dimensionless ranking heuristic; lower scores receive priority.
- `fail_reason`: semicolon-separated rejection codes; empty for survivors.
- `payload_sequence`: entire 120 nt encoded segment, including seed and RS redundancy.
- `full_sequence`: 180 nt strand including flanking regions.
- `MFE_structure`: full-strand dot-bracket structure.
- `MFE_deltaG`: predicted MFE in kcal/mol.
- Prefixes `F` and `R`: forward region and appended reverse-primer binding region.
- `3p_core`, `3p`: 5 nt core and 8 nt window at the primer 3′-corresponding end.
- `J_F`, `J_R`: 20 nt junction windows.
- `paired_count`, `pair_count`: MFE paired positions or pairs between indicated regions.
- `paired_fraction`, `gc_content`: fractions on a 0–1 scale.
- `mean_pair_probability`: mean of per-position row sums of the matrix returned by `fc.bpp()`, without explicit matrix symmetrization. See the screening-method description in README.
- `max_continuous_stem`, `max_stem_length`, `payload_internal_max_stem`: maximum consecutive MFE-paired run under the script's region-specific rules; partners are not checked for uninterrupted stacking.
- `homopolymer_max_length`: longest identical-base run in the encoded segment, nt.
- `terminal_primer_payload_complementarity`: longest common substring in the implemented comparisons between terminal 10 nt of the encoded segment and reverse complements of the selected 8 nt primer windows.
- `basic_filter_fail`: first detected basic-filter rejection code, otherwise empty.

`simulated_reads_55200.txt.gz` contains one simulated read per line (duplicates retained). `simulated_reads_unique_sorted.txt.gz` contains unique reads sorted by abundance for decoding. `simulated_reads_abundance_sorted.txt.gz` contains whitespace-separated count and sequence columns in that same order. The candidate mapping tables (`candidate_id_mapping.part001.tsv.gz` through `candidate_id_mapping.part003.tsv.gz`) use tab-separated fields. Each part contains its own header and 100,000 records; parts are ordered by candidate ID.

The complete screening report comprises six files named `screening_all_300000.part001.csv.gz` through `screening_all_300000.part006.csv.gz`. Each contains the same CSV header and 50,000 records, in original report order.
