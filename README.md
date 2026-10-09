# Electronically-Programmable-Multilayer-DNA-Storage

Mona Lisa DNA Fountain encoding and secondary-structure screening

Version: 1.0.0

This dataset accompanies the Mona Lisa DNA-storage experiment and provides the oligonucleotide sequences, secondary-structure screening results, analysis code and decoding validation data. A pool of 300,000 candidate encoded segments was screened to retain 10,163 sequences (3.39%). The 5,200 candidates with the lowest structural risk scores were selected from this set, yielding the final library (1.73% of the initial pool), from which the original image was successfully decoded.

## Dataset contents

| Path | Description |
| --- | --- |
| `data/input/mona_lisa.jpg` | Input image (97,530 bytes). |
| `data/sequences/candidates_300000.fasta.gz` | All 300,000 candidate encoded segments (120 nt), in generation order. |
| `data/sequences/survivors_10163_encoded.fasta` | The 10,163 retained encoded segments, in increasing risk-score order. |
| `data/sequences/final_5200_encoded.fasta` | Final 5,200 encoded segments (120 nt), used as decoder input. |
| `data/sequences/final_5200_full_length.fasta` | Corresponding full-length oligonucleotide designs (180 nt), including the fixed flanking regions. |
| `data/screening/screening_all_300000.part*.csv.gz` | Complete screening results with candidate identifiers. |
| `data/screening/screening_final_5200.csv` | Screening metrics for the final library. |
| `data/validation/` | Simulated reads, unique reads ordered by abundance and corresponding read counts. |
| `code/screen_secondary_structure.py` | Secondary-structure screening implementation. |
| `code/dna_fountain/` | DNA Fountain encoding and decoding source code. |
| `code/select_from_report.py` | Selection of the lowest-risk candidates from a screening report. |
| `code/verify_archive.py` | Verification of file checksums, sequence lengths, library membership and selection order. |
| `metadata/` | Identifier mappings, field definitions, software versions, validation results and file inventory. |
| `SHA256SUMS.txt` | SHA-256 checksums for all other files in the dataset. |

All sequences are provided in the 5′–3′ orientation. FASTA records use stable identifiers; the mapping tables link these identifiers to the original sequence headers. Detailed field definitions are provided in `metadata/data_dictionary.md`. Files ending in `.gz` are gzip-compressed.

## Oligonucleotide design

Each full-length strand comprises a 30 nt forward-primer region, a 120 nt encoded segment and a 30 nt reverse-primer binding region. The encoded segment contains a 16 nt seed, a 96 nt data payload and an 8 nt Reed–Solomon redundancy segment. In the screening code and report fields, `payload` denotes the entire 120 nt encoded segment.

| Region | Sequence (5′–3′) |
| --- | --- |
| Forward primer | `TCTGCAAGTAGCCAAGGGTAAGCAAGGATC` |
| Reverse primer | `GAATCTGCACGAAGCACTTATGAACCCTTG` |
| Reverse-primer binding region | `CAAGGGTTCATAAGTGCTTCGTGCAGATTC` |

The full-length sequences were constructed as forward primer + encoded segment + reverse complement of the reverse primer.

## Secondary-structure screening

Full-length strands were analysed using ViennaRNA 2.7.2 with the Mathews 2004 DNA thermodynamic parameter set, loaded through `RNA.params_load_DNA_Mathews2004()`. Other model settings used the library defaults; temperature and salt conditions were not explicitly specified in the script.

The forward 3′ window comprises the last 8 nt of the forward-primer region, and the reverse 3′ window comprises the first 8 nt of the appended reverse-primer binding region. The corresponding terminal cores comprise 5 nt. Each primer–encoded-segment junction window spans 10 nt on either side of the boundary.

Candidates were excluded if any of the following criteria were met:

- At least 3 paired positions in either 5 nt core in the minimum free-energy (MFE) structure, or a mean core pairing metric of at least 0.50.
- At least 4 consecutive paired positions in either 8 nt window.
- At least 3 pairs between either 3′ window and the opposite primer region, or at least 8 pairs between the two primer regions with involvement of a 3′ window.
- At least 6 consecutive paired positions in either junction window.
- MFE below −35 kcal/mol with at least 10 paired positions across the two primer regions.

Additional sequence filters required GC fractions of 0.40–0.60 in both the encoded segment and the full strand, encoded-segment homopolymers of at most 3 nt, and a terminal complementarity metric of at most 4 nt. The terminal comparison is defined in the data dictionary and screening code.

Retained sequences were ranked using a weighted score incorporating terminal-core and window pairing, consecutive paired runs, primer paired fractions, junction pairing and pairing between the flanking regions. Internal encoded-segment paired runs longer than 11 positions received a smaller penalty. The exact score and classification rules are specified in `code/screen_secondary_structure.py`. All 10,163 retained sequences were classified as `BORDERLINE` under these rules. Classification and risk scores are computational screening metrics and have not been calibrated against measured RPA efficiency.

The report fields containing `stem` denote consecutive paired positions; they do not require the pairing partners to form an uninterrupted stacked helix. The `mean_pair_probability` fields are calculated by averaging row sums of the matrix returned by `fc.bpp()`, without explicit matrix symmetrization. These definitions apply to the supplied screening results and code.

## Software and usage

The validation environment used Python 3.11.2 and the dependency versions listed in `requirements.txt` and `metadata/environment.json`. Building the Cython extensions requires a C compiler. The following commands use a POSIX shell and should be run from the dataset root directory.

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
(cd code/dna_fountain && python setup.py build_ext --inplace)
python code/verify_archive.py
mkdir -p results
```

### Generate candidate sequences

```bash
python code/dna_fountain/encode.py --file_in data/input/mona_lisa.jpg --size 24 --max_homopolymer 3 --gc 0.05 --rs 2 --delta 0.001 --c_dist 0.025 --alpha 0.20 --stop 300000 --out results/candidates.fasta
```

### Select the final library from the screening report

```bash
python code/select_from_report.py --report "data/screening/screening_all_300000.part*.csv.gz" --target 5200 --out results/selected_5200.fasta
```

This command reads all six report parts and reproduces selection without repeating structure prediction. Each part is an independently readable gzip-compressed CSV file with a header and 50,000 candidate records. Numerical part order preserves the original report order. Keep the glob pattern in quotation marks; the selection script expands it internally.

### Repeat secondary-structure screening

```bash
python -c "import gzip,shutil; shutil.copyfileobj(gzip.open('data/sequences/candidates_300000.fasta.gz','rb'),open('results/candidates.fasta','wb'))"
python code/screen_secondary_structure.py --input results/candidates.fasta --report results/screening.csv --selected-out results/selected.fasta --target 5200 --dna-param mathews2004 --workers 4
```

For a preliminary check, add `--limit 20` to screen the first 20 candidates. Regenerated reports use the normalized candidate identifiers; their correspondence to the original identifiers is provided in `metadata/candidate_id_mapping.part001.tsv.gz` through `metadata/candidate_id_mapping.part003.tsv.gz`.

### Decode the final library

```bash
python code/dna_fountain/decode.py --file_in data/sequences/final_5200_encoded.fasta --fasta --header_size 4 --rs 2 --delta 0.001 --c_dist 0.025 --chunk_num 4064 --size 24 --max_homopolymer 3 --gc 0.05 --max_hamming 0 --out results/decoded_padded.bin
python -c "from pathlib import Path; import hashlib; b=Path('results/decoded_padded.bin').read_bytes(); assert len(b)==97536; b=b[:97530]; assert b==Path('data/input/mona_lisa.jpg').read_bytes(); Path('results/decoded.jpg').write_bytes(b); print(hashlib.sha256(b).hexdigest())"
```

The decoded output contains 97,536 bytes, including 6 terminal padding bytes. Removing the padding recovers the original 97,530-byte image. The decoder also writes the diagnostic file `seen_barocdes.json` in its working directory.

### Decode simulated reads

Decompress `data/validation/simulated_reads_unique_sorted.txt.gz` and supply the resulting file to the decoder with the same parameters, omitting `--fasta`.

The simulated dataset contains 55,200 reads: the 5,200-sequence library with five records replaced by single-nucleotide variants, 30,000 additional single-nucleotide variants and 20,000 random sequences. Reads were shuffled and then sorted by abundance for decoding. The supplied files contain 54,938 unique sequences. The simulation generator and random seed are unavailable; the supplied reads enable repetition of the decoding analysis for this specific simulated dataset.

## Validation

The following checks were completed in the documented environment:

- The first 100 regenerated candidates matched the corresponding supplied sequences exactly.
- Repeated screening of the first 20 candidates reproduced all reported screening fields.
- Selection from the complete screening report reproduced the final 5,200-sequence FASTA exactly.
- Decoding both the final library and the supplied simulated reads recovered the original image after removal of the 6 padding bytes.
- Sequence counts, lengths, uniqueness, nested library membership and simulated-read abundance ordering were verified.

The recovered image has SHA-256 checksum `b7d510972c41453b710c268762d4b267129b3c3a210e21dcdb60af4d4a11c445`. Validation details are provided in `metadata/validation_summary.json` and `metadata/validation_checks.json`. Repeat encoding and folding checks were limited to the subsets specified above; complete regeneration of the 300,000-candidate dataset and installation in a fresh environment were not tested.

## Code attribution

DNA Fountain was developed by Yaniv Erlich and Dina Zielinski: https://github.com/TeamErlich/dna-fountain. The included Python 3 port credits Yihang Du, Wenrong Wu and Justin Brody. Its original documentation is provided in `code/dna_fountain/UPSTREAM_README.md`; the exact source revision is unavailable. The DNA Fountain code is distributed under GPLv3-or-later, with copyright notices and license text provided in `code/dna_fountain/COPYING` and the source files.

## Partitioned data files

The candidate mapping table is distributed as three independently readable gzip-compressed TSV files, each containing a header and 100,000 records. Parts 001, 002 and 003 cover `candidate_000001`–`candidate_100000`, `candidate_100001`–`candidate_200000` and `candidate_200001`–`candidate_300000`, respectively. To combine them, decompress and concatenate the parts in numerical order, retaining the header from the first part only.

The complete screening report is distributed as `screening_all_300000.part001.csv.gz` through `screening_all_300000.part006.csv.gz`. To reconstruct a single report, decompress and concatenate these files in numerical order, retaining only the first header. 
