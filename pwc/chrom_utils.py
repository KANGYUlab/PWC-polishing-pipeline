#!/usr/bin/env python3
"""Chromosome name utilities: normalization and cross-file compatibility checks."""

from __future__ import annotations

import os
import subprocess
import sys
from typing import Iterable


def read_fasta_chroms(fasta_path: str) -> list[str]:
    chroms: list[str] = []
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            chroms.append(line[1:].split()[0])
    if not chroms:
        raise ValueError(f"No FASTA sequences found in {fasta_path}")
    return chroms


def read_bed_chroms(bed_path: str) -> set[str]:
    chroms: set[str] = set()
    with open(bed_path) as fh:
        for line_num, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"{bed_path}:{line_num}: expected at least 3 BED columns")
            chroms.add(fields[0])
    return chroms


def read_bam_chroms(bam_path: str) -> list[str]:
    try:
        proc = subprocess.run(
            ["samtools", "idxstats", bam_path],
            check=True,
            capture_output=True,
            text=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        raise ValueError(f"Cannot read BAM chromosomes from {bam_path}: {exc}") from exc

    chroms: list[str] = []
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        chroms.append(line.split("\t", 1)[0])
    if not chroms:
        raise ValueError(f"No chromosomes reported by samtools idxstats for {bam_path}")
    return chroms


def chrom_aliases(chrom: str) -> set[str]:
    aliases = {chrom}
    if chrom.startswith("chr"):
        base = chrom[3:]
        aliases.add(base)
    else:
        aliases.add(f"chr{chrom}")

    mt_like = {"M", "MT", "chrM", "chrMT"}
    if chrom in mt_like:
        aliases.update(mt_like)
    return aliases


def build_chrom_mapper(query_chroms: Iterable[str], ref_chroms: Iterable[str]) -> dict[str, str]:
    ref_set = set(ref_chroms)
    mapper: dict[str, str] = {}
    for query in sorted(set(query_chroms)):
        if query in ref_set:
            mapper[query] = query
            continue
        matched = None
        for alias in chrom_aliases(query):
            if alias in ref_set:
                matched = alias
                break
        mapper[query] = matched if matched is not None else query
    return mapper


def normalize_bed_rows(rows: list[list[str]], chrom_mapper: dict[str, str]) -> tuple[list[list[str]], int]:
    normalized: list[list[str]] = []
    remapped = 0
    for row in rows:
        old_chrom = row[0]
        new_chrom = chrom_mapper.get(old_chrom, old_chrom)
        if new_chrom != old_chrom:
            remapped += 1
        normalized.append([new_chrom] + row[1:])
    return normalized, remapped


def chrom_overlap_stats(query_chroms: Iterable[str], ref_chroms: Iterable[str]) -> tuple[int, int, list[str]]:
    ref_set = set(ref_chroms)
    mapper = build_chrom_mapper(query_chroms, ref_chroms)
    query_unique = sorted(set(query_chroms))
    matched = [chrom for chrom in query_unique if mapper.get(chrom, chrom) in ref_set]
    unmatched = [chrom for chrom in query_unique if mapper.get(chrom, chrom) not in ref_set]
    return len(matched), len(query_unique), unmatched


def bam_is_indexed(bam_path: str) -> bool:
    for suffix in (".bai", ".csi"):
        if os.path.exists(bam_path + suffix):
            return True
    try:
        subprocess.run(
            ["samtools", "idxstats", bam_path],
            check=True,
            capture_output=True,
            text=True,
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def validate_bed_file(bed_path: str, label: str = "BED") -> int:
    row_count = 0
    with open(bed_path) as fh:
        for line_num, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"{label} {bed_path}:{line_num}: expected at least 3 columns")
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{label} {bed_path}:{line_num}: start/end must be integers ({fields[1]}, {fields[2]})"
                ) from exc
            if start < 0 or end < 0:
                raise ValueError(f"{label} {bed_path}:{line_num}: coordinates must be non-negative")
            if end < start:
                raise ValueError(f"{label} {bed_path}:{line_num}: end ({end}) < start ({start})")
            row_count += 1
    if row_count == 0:
        raise ValueError(f"{label} file has no data rows: {bed_path}")
    return row_count


def validate_pwc_inputs(
    ref_fasta: str,
    variants_bed: str,
    short_bam: str,
    hifi_bam: str,
    ont_bam: str,
) -> None:
    paths = {
        "reference FASTA": ref_fasta,
        "variants BED": variants_bed,
        "short-read BAM": short_bam,
        "HiFi BAM": hifi_bam,
        "ONT BAM": ont_bam,
    }
    for label, path in paths.items():
        if not os.path.isfile(path):
            raise ValueError(f"{label} not found: {path}")

    ref_chroms = read_fasta_chroms(ref_fasta)
    variant_rows = validate_bed_file(variants_bed, "variants BED")

    variant_chroms = read_bed_chroms(variants_bed)
    matched, total, unmatched = chrom_overlap_stats(variant_chroms, ref_chroms)
    if matched == 0:
        raise ValueError(
            "No variant BED chromosomes match the reference FASTA after chr-prefix normalization. "
            f"Reference example: {ref_chroms[:3]}. "
            f"Variant example: {sorted(variant_chroms)[:3]}."
        )
    if unmatched:
        print(
            f"Warning: {len(unmatched)} variant chromosome(s) have no reference match "
            f"(first few: {unmatched[:5]})",
            file=sys.stderr,
        )

    for label, bam_path in (
        ("short-read BAM", short_bam),
        ("HiFi BAM", hifi_bam),
        ("ONT BAM", ont_bam),
    ):
        if not bam_is_indexed(bam_path):
            raise ValueError(f"{label} is missing an index (.bai/.csi) or is unreadable: {bam_path}")
        bam_chroms = read_bam_chroms(bam_path)
        bam_matched, bam_total, bam_unmatched = chrom_overlap_stats(bam_chroms, ref_chroms)
        if bam_matched == 0:
            raise ValueError(
                f"{label} chromosomes do not match reference FASTA naming. "
                f"BAM example: {bam_chroms[:3]}. Reference example: {ref_chroms[:3]}."
            )
        if bam_unmatched:
            print(
                f"Warning: {label} has {len(bam_unmatched)} chromosome(s) not in reference "
                f"(first few: {bam_unmatched[:5]})",
                file=sys.stderr,
            )

    print(
        f"[validate] OK: reference={len(ref_chroms)} seqs, "
        f"variants={variant_rows} regions, "
        f"variant chrom overlap {matched}/{total}",
        file=sys.stderr,
    )
