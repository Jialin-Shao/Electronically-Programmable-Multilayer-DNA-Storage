#!/usr/bin/env python3
"""Screen fountain payloads for fixed-primer secondary-structure risk."""

from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import RNA


F_PRIMER = "TCTGCAAGTAGCCAAGGGTAAGCAAGGATC"
R_PRIMER = "GAATCTGCACGAAGCACTTATGAACCCTTG"
R_BINDING = "CAAGGGTTCATAAGTGCTTCGTGCAGATTC"


@dataclass(frozen=True)
class Record:
    seq_id: str
    payload: str
    header: str | None


def parse_records(path: Path) -> list[Record]:
    records: list[Record] = []
    header: str | None = None
    seq_parts: list[str] = []
    counter = 1

    def flush() -> None:
        nonlocal counter, header, seq_parts
        if not seq_parts:
            return
        payload = "".join(seq_parts).upper().replace("U", "T")
        seq_id = re.sub(r"[^A-Za-z0-9_.:-]+", "_", header[1:].strip()) if header else str(counter)
        records.append(Record(seq_id=seq_id, payload=payload, header=header))
        counter += 1
        seq_parts = []

    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line
            else:
                seq_parts.append(line)
        flush()
    return records


def max_homopolymer(seq: str) -> int:
    if not seq:
        return 0
    return max(len(m.group(0)) for m in re.finditer(r"([ACGT])\1*", seq))


def gc_fraction(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def longest_common_substring(a: str, b: str) -> int:
    best = 0
    for i in range(len(a)):
        for j in range(len(b)):
            k = 0
            while i + k < len(a) and j + k < len(b) and a[i + k] == b[j + k]:
                k += 1
            best = max(best, k)
    return best


def terminal_complementarity(payload: str) -> int:
    f_3p = F_PRIMER[-8:]
    r_3p = R_BINDING[:8]
    checks = [
        (payload[:10], revcomp(f_3p)),
        (payload[:10], revcomp(r_3p)),
        (payload[-10:], revcomp(f_3p)),
        (payload[-10:], revcomp(r_3p)),
    ]
    return max(longest_common_substring(a, b) for a, b in checks)


def parse_dot_bracket(structure: str) -> dict[int, int]:
    stack: list[int] = []
    pairs: dict[int, int] = {}
    for idx, ch in enumerate(structure):
        if ch == "(":
            stack.append(idx)
        elif ch == ")":
            if not stack:
                continue
            left = stack.pop()
            pairs[left] = idx
            pairs[idx] = left
    return pairs


def region_indices(start: int, end: int) -> set[int]:
    return set(range(start, end))


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def paired_count(region: set[int], pairs: dict[int, int]) -> int:
    return sum(1 for i in region if i in pairs)


def pair_count_between(a: set[int], b: set[int], pairs: dict[int, int]) -> int:
    return sum(1 for i in a if pairs.get(i) in b)


def max_paired_run(region: set[int], pairs: dict[int, int]) -> int:
    best = run = 0
    for i in sorted(region):
        if i in pairs:
            run += 1
            best = max(best, run)
        else:
            run = 0
    return best


def internal_payload_max_stem(payload_region: set[int], pairs: dict[int, int]) -> int:
    payload_pairs = {i: j for i, j in pairs.items() if i in payload_region and j in payload_region}
    return max_paired_run(payload_region, payload_pairs)


def base_pair_probabilities(fc: RNA.fold_compound, n: int) -> list[float]:
    fc.pf()
    matrix = fc.bpp()
    probs = [0.0] * n
    for i in range(1, n + 1):
        row = matrix[i]
        total = 0.0
        for j in range(1, n + 1):
            if i == j:
                continue
            total += row[j]
        probs[i - 1] = total
    return probs


def classify_and_score(metrics: dict[str, float | int | str]) -> tuple[str, str, float]:
    fail_reasons: list[str] = []
    if metrics["F_3p_core_paired_count"] >= 3:
        fail_reasons.append("FAIL_F_3P_CORE_MFE_PAIRED")
    if metrics["R_3p_core_paired_count"] >= 3:
        fail_reasons.append("FAIL_R_3P_CORE_MFE_PAIRED")
    if metrics["F_3p_core_mean_pair_probability"] >= 0.50:
        fail_reasons.append("FAIL_F_3P_CORE_HIGH_PAIR_PROBABILITY")
    if metrics["R_3p_core_mean_pair_probability"] >= 0.50:
        fail_reasons.append("FAIL_R_3P_CORE_HIGH_PAIR_PROBABILITY")
    if metrics["F_3p_max_continuous_stem"] >= 4:
        fail_reasons.append("FAIL_F_3P_CONTINUOUS_STEM")
    if metrics["R_3p_max_continuous_stem"] >= 4:
        fail_reasons.append("FAIL_R_3P_CONTINUOUS_STEM")
    if metrics["F_3p_to_R_pair_count"] >= 3:
        fail_reasons.append("FAIL_F_3P_TO_R_PRIMER_PAIRING")
    if metrics["R_3p_to_F_pair_count"] >= 3:
        fail_reasons.append("FAIL_R_3P_TO_F_PRIMER_PAIRING")
    if metrics["F_R_total_pair_count"] >= 8 and (
        metrics["F_3p_to_R_pair_count"] > 0 or metrics["R_3p_to_F_pair_count"] > 0
    ):
        fail_reasons.append("FAIL_PRIMER_PRIMER_PAIRING_WITH_3P")
    if metrics["J_F_max_stem_length"] >= 6:
        fail_reasons.append("FAIL_JUNCTION_F_STEM")
    if metrics["J_R_max_stem_length"] >= 6:
        fail_reasons.append("FAIL_JUNCTION_R_STEM")
    if metrics["MFE_deltaG"] < -35 and metrics["primer_paired_bases"] >= 10:
        fail_reasons.append("FAIL_EXTREME_MFE_WITH_PRIMER_PAIRING")
    if metrics["basic_filter_fail"]:
        fail_reasons.append(str(metrics["basic_filter_fail"]))

    score = (
        40 * float(metrics["F_3p_core_mean_pair_probability"])
        + 40 * float(metrics["R_3p_core_mean_pair_probability"])
        + 20 * float(metrics["F_3p_mean_pair_probability"])
        + 20 * float(metrics["R_3p_mean_pair_probability"])
        + 10 * int(metrics["F_3p_paired_count"])
        + 10 * int(metrics["R_3p_paired_count"])
        + 20 * int(metrics["F_3p_max_continuous_stem"])
        + 20 * int(metrics["R_3p_max_continuous_stem"])
        + 15 * float(metrics["F_primer_paired_fraction"])
        + 15 * float(metrics["R_primer_paired_fraction"])
        + 10 * int(metrics["J_F_max_stem_length"])
        + 10 * int(metrics["J_R_max_stem_length"])
        + 15 * int(metrics["F_R_total_pair_count"])
        + 3 * max(0, int(metrics["payload_internal_max_stem"]) - 11)
    )

    if fail_reasons:
        return "FAIL", ";".join(fail_reasons), score

    borderline = (
        metrics["F_3p_core_mean_pair_probability"] >= 0.35
        or metrics["R_3p_core_mean_pair_probability"] >= 0.35
        or metrics["F_primer_paired_fraction"] >= 0.40
        or metrics["R_primer_paired_fraction"] >= 0.40
        or metrics["J_F_max_stem_length"] >= 5
        or metrics["J_R_max_stem_length"] >= 5
        or 5 <= metrics["F_R_total_pair_count"] <= 7
        or score > 30
    )
    return ("BORDERLINE" if borderline else "PASS"), "", score


def screen_record(record: Record) -> dict[str, float | int | str]:
    payload = record.payload
    full_seq = F_PRIMER + payload + R_BINDING
    n = len(full_seq)
    f_start, f_end = 0, len(F_PRIMER)
    p_start, p_end = f_end, f_end + len(payload)
    r_start, r_end = p_end, n

    regions = {
        "F_all": region_indices(f_start, f_end),
        "F_3p": region_indices(f_end - 8, f_end),
        "F_3p_core": region_indices(f_end - 5, f_end),
        "R_all": region_indices(r_start, r_end),
        "R_3p": region_indices(r_start, r_start + 8),
        "R_3p_core": region_indices(r_start, r_start + 5),
        "payload": region_indices(p_start, p_end),
        "J_F": region_indices(max(0, f_end - 10), min(n, p_start + 10)),
        "J_R": region_indices(max(0, p_end - 10), min(n, r_start + 10)),
        "F_terminal10": region_indices(f_end - 10, f_end),
        "R_terminal10": region_indices(r_start, r_start + 10),
        "payload_left20": region_indices(p_start, min(p_end, p_start + 20)),
        "payload_right20": region_indices(max(p_start, p_end - 20), p_end),
    }

    fc = RNA.fold_compound(full_seq)
    structure, mfe = fc.mfe()
    probs = base_pair_probabilities(fc, n)
    pairs = parse_dot_bracket(structure)

    basic_fail = ""
    if gc_fraction(payload) < 0.40 or gc_fraction(payload) > 0.60:
        basic_fail = "FAIL_BASIC_PAYLOAD_GC"
    elif gc_fraction(full_seq) < 0.40 or gc_fraction(full_seq) > 0.60:
        basic_fail = "FAIL_BASIC_FULL_GC"
    elif max_homopolymer(payload) > 3:
        basic_fail = "FAIL_BASIC_HOMOPOLYMER"
    elif terminal_complementarity(payload) > 4:
        basic_fail = "FAIL_BASIC_TERMINAL_COMPLEMENTARITY"

    f_primer_paired = paired_count(regions["F_all"], pairs)
    r_primer_paired = paired_count(regions["R_all"], pairs)
    metrics: dict[str, float | int | str] = {
        "sequence_id": record.seq_id,
        "payload_sequence": payload,
        "full_sequence": full_seq,
        "MFE_structure": structure,
        "MFE_deltaG": float(mfe),
        "payload_gc_content": gc_fraction(payload),
        "full_sequence_gc_content": gc_fraction(full_seq),
        "homopolymer_max_length": max_homopolymer(payload),
        "terminal_primer_payload_complementarity": terminal_complementarity(payload),
        "basic_filter_fail": basic_fail,
        "F_3p_core_paired_count": paired_count(regions["F_3p_core"], pairs),
        "R_3p_core_paired_count": paired_count(regions["R_3p_core"], pairs),
        "F_3p_paired_count": paired_count(regions["F_3p"], pairs),
        "R_3p_paired_count": paired_count(regions["R_3p"], pairs),
        "F_3p_core_mean_pair_probability": mean([probs[i] for i in regions["F_3p_core"]]),
        "R_3p_core_mean_pair_probability": mean([probs[i] for i in regions["R_3p_core"]]),
        "F_3p_mean_pair_probability": mean([probs[i] for i in regions["F_3p"]]),
        "R_3p_mean_pair_probability": mean([probs[i] for i in regions["R_3p"]]),
        "F_primer_mean_pair_probability": mean([probs[i] for i in regions["F_all"]]),
        "R_primer_mean_pair_probability": mean([probs[i] for i in regions["R_all"]]),
        "J_F_mean_pair_probability": mean([probs[i] for i in regions["J_F"]]),
        "J_R_mean_pair_probability": mean([probs[i] for i in regions["J_R"]]),
        "F_3p_max_continuous_stem": max_paired_run(regions["F_3p"], pairs),
        "R_3p_max_continuous_stem": max_paired_run(regions["R_3p"], pairs),
        "F_primer_paired_fraction": f_primer_paired / len(regions["F_all"]),
        "R_primer_paired_fraction": r_primer_paired / len(regions["R_all"]),
        "J_F_max_stem_length": max_paired_run(regions["J_F"], pairs),
        "J_R_max_stem_length": max_paired_run(regions["J_R"], pairs),
        "J_F_paired_count": paired_count(regions["J_F"], pairs),
        "J_R_paired_count": paired_count(regions["J_R"], pairs),
        "F_R_total_pair_count": pair_count_between(regions["F_all"], regions["R_all"], pairs),
        "F_3p_to_R_pair_count": pair_count_between(regions["F_3p"], regions["R_all"], pairs),
        "R_3p_to_F_pair_count": pair_count_between(regions["R_3p"], regions["F_all"], pairs),
        "F_terminal_to_payload_pair_count": pair_count_between(regions["F_terminal10"], regions["payload_left20"], pairs),
        "R_terminal_to_payload_pair_count": pair_count_between(regions["R_terminal10"], regions["payload_right20"], pairs),
        "payload_internal_max_stem": internal_payload_max_stem(regions["payload"], pairs),
        "primer_paired_bases": f_primer_paired + r_primer_paired,
    }
    classification, fail_reason, risk_score = classify_and_score(metrics)
    metrics["classification"] = classification
    metrics["fail_reason"] = fail_reason
    metrics["risk_score"] = risk_score
    return metrics


def init_worker(dna_param: str) -> None:
    if dna_param == "mathews2004":
        RNA.params_load_DNA_Mathews2004()
    else:
        RNA.params_load_DNA_Mathews1999()


def write_selected(records: list[Record], rows: list[dict[str, float | int | str]], target: int, out: Path) -> None:
    by_id = {record.seq_id: record for record in records}
    eligible = [row for row in rows if row["classification"] in {"PASS", "BORDERLINE"}]
    eligible.sort(key=lambda row: (row["classification"] != "PASS", float(row["risk_score"])))
    selected = eligible[:target]
    with out.open("w") as handle:
        for idx, row in enumerate(selected, 1):
            record = by_id[str(row["sequence_id"])]
            handle.write(f">selected_{idx} source={record.seq_id} class={row['classification']} risk={float(row['risk_score']):.3f}\n")
            handle.write(f"{record.payload}\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path, help="Input .dna/.fasta payload file")
    parser.add_argument("--report", required=True, type=Path, help="CSV report path")
    parser.add_argument("--selected-out", type=Path, help="Optional selected payload output")
    parser.add_argument("--target", type=int, default=0, help="Number of selected payloads to write")
    parser.add_argument("--limit", type=int, default=0, help="Only screen the first N records")
    parser.add_argument("--dna-param", default="mathews2004", choices=["mathews2004", "mathews1999"])
    parser.add_argument("--workers", type=int, default=1, help="Parallel worker processes")
    args = parser.parse_args()

    init_worker(args.dna_param)

    records = parse_records(args.input)
    if args.limit:
        records = records[: args.limit]

    rows: list[dict[str, float | int | str]] = []
    if args.workers <= 1:
        iterator = map(screen_record, records)
    else:
        pool = mp.Pool(processes=args.workers, initializer=init_worker, initargs=(args.dna_param,))
        iterator = pool.imap(screen_record, records, chunksize=20)
    try:
        for idx, row in enumerate(iterator, 1):
            rows.append(row)
            if idx % 100 == 0:
                print(f"screened {idx}/{len(records)}", file=sys.stderr)
    finally:
        if args.workers > 1:
            pool.close()
            pool.join()

    fieldnames = [
        "sequence_id",
        "classification",
        "risk_score",
        "fail_reason",
        "MFE_deltaG",
        "MFE_structure",
        "payload_sequence",
        "full_sequence",
        "F_3p_core_paired_count",
        "R_3p_core_paired_count",
        "F_3p_core_mean_pair_probability",
        "R_3p_core_mean_pair_probability",
        "F_3p_paired_count",
        "R_3p_paired_count",
        "F_3p_mean_pair_probability",
        "R_3p_mean_pair_probability",
        "F_3p_max_continuous_stem",
        "R_3p_max_continuous_stem",
        "F_primer_paired_fraction",
        "R_primer_paired_fraction",
        "F_primer_mean_pair_probability",
        "R_primer_mean_pair_probability",
        "J_F_max_stem_length",
        "J_R_max_stem_length",
        "J_F_mean_pair_probability",
        "J_R_mean_pair_probability",
        "J_F_paired_count",
        "J_R_paired_count",
        "F_R_total_pair_count",
        "F_3p_to_R_pair_count",
        "R_3p_to_F_pair_count",
        "F_terminal_to_payload_pair_count",
        "R_terminal_to_payload_pair_count",
        "payload_internal_max_stem",
        "payload_gc_content",
        "full_sequence_gc_content",
        "homopolymer_max_length",
        "terminal_primer_payload_complementarity",
        "basic_filter_fail",
    ]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    if args.selected_out and args.target:
        write_selected(records, rows, args.target, args.selected_out)

    counts: dict[str, int] = {}
    for row in rows:
        counts[str(row["classification"])] = counts.get(str(row["classification"]), 0) + 1
    print("screened", len(rows), counts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
