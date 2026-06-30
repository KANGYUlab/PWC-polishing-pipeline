#!/usr/bin/env python3
"""
Split a variant BED into homo.bed and nonhomo.bed using homopolymer regions.

Homopolymer regions are produced by homo.py from the reference genome.
Variants overlapping those regions go to homo.bed; the rest go to nonhomo.bed.

Input BED tips:
  - SAS base-editing (BE) BED files can be used directly.
  - SNP / small INDEL calls from variant callers should be padded (e.g. --pad 10)
    so each site is expanded by 10 bp on both sides before splitting.
"""

import argparse
import os
import subprocess
import sys
import tempfile

from chrom_utils import build_chrom_mapper, normalize_bed_rows, read_bed_chroms, read_fasta_chroms


def read_bed3(path: str) -> list[list[str]]:
    rows: list[list[str]] = []
    with open(path) as fh:
        for line_num, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_num}: expected at least 3 BED columns")
            rows.append(fields)
    return rows


def write_bed3(path: str, rows: list[list[str]]) -> None:
    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(row[:3]) + "\n")


def pad_bed_rows(rows: list[list[str]], pad: int) -> list[list[str]]:
    padded: list[list[str]] = []
    for row in rows:
        chrom = row[0]
        start = max(0, int(row[1]) - pad)
        end = int(row[2]) + pad
        if end <= start:
            end = start + 1
        padded.append([chrom, str(start), str(end)] + row[3:])
    return padded


def homo_regions_to_bed3(homo_regions_path: str) -> str:
    """Convert homo.py TSV (>=3 columns) to a temporary 3-column BED file."""
    rows = read_bed3(homo_regions_path)
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False, prefix="homo_regions_")
    write_bed3(tmp.name, rows)
    tmp.close()
    return tmp.name


def split_with_bedtools(variant_bed: str, homo_bed: str, homo_out: str, nonhomo_out: str) -> None:
    with open(homo_out, "w") as homo_fh:
        subprocess.run(
            ["bedtools", "intersect", "-a", variant_bed, "-b", homo_bed, "-wa", "-u"],
            check=True,
            stdout=homo_fh,
        )
    with open(nonhomo_out, "w") as nonhomo_fh:
        subprocess.run(
            ["bedtools", "intersect", "-a", variant_bed, "-b", homo_bed, "-wa", "-v"],
            check=True,
            stdout=nonhomo_fh,
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split variant BED into homo.bed and nonhomo.bed using homopolymer regions from homo.py.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Examples:
  # SAS BE BED — use directly
  python3 pwc/split_variants.py homo_regions.tsv sas_be.bed -o 02_homo_nonhomo

  # Variant-caller SNP / small INDEL — pad 10 bp on each side first
  python3 pwc/split_variants.py homo_regions.tsv deepvariant.vcf.bed --pad 10 -o 02_homo_nonhomo
""",
    )
    parser.add_argument("homo_regions", help="Homopolymer regions from homo.py (BED/TSV, first 3 columns used).")
    parser.add_argument("variants", help="Variant BED to split (3+ columns).")
    parser.add_argument("-o", "--output-dir", default="02_homo_nonhomo", help="Output directory (default: 02_homo_nonhomo).")
    parser.add_argument(
        "--pad",
        type=int,
        default=0,
        help="Expand each variant interval by N bp on both sides before splitting (use 10 for SNP/small INDEL).",
    )
    parser.add_argument(
        "--ref-fasta",
        default=None,
        help="Reference FASTA for chromosome name normalization (chr1 <-> 1, MT <-> chrM).",
    )
    parser.add_argument(
        "--no-normalize-chrom",
        action="store_true",
        help="Disable automatic chromosome name normalization.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    variant_rows = read_bed3(args.variants)
    if not variant_rows:
        raise SystemExit(f"No variant records found in {args.variants}")

    if not args.no_normalize_chrom:
        if args.ref_fasta:
            ref_chroms = read_fasta_chroms(args.ref_fasta)
        else:
            ref_chroms = sorted(read_bed_chroms(args.homo_regions))
            if not ref_chroms:
                raise SystemExit(f"No homopolymer regions found in {args.homo_regions}")

        variant_chroms = {row[0] for row in variant_rows}
        chrom_mapper = build_chrom_mapper(variant_chroms, ref_chroms)
        variant_rows, remapped = normalize_bed_rows(variant_rows, chrom_mapper)
        if remapped:
            print(
                f"Normalized chromosome names on {remapped} variant row(s) "
                f"to match reference naming",
                file=sys.stderr,
            )
        unmatched = sorted({row[0] for row in variant_rows if row[0] not in set(ref_chroms)})
        if unmatched:
            print(
                f"Warning: {len(unmatched)} chromosome name(s) still unmatched after normalization: "
                f"{unmatched[:5]}",
                file=sys.stderr,
            )

    if args.pad > 0:
        variant_rows = pad_bed_rows(variant_rows, args.pad)

    padded_variant_bed = os.path.join(args.output_dir, "variants.padded.bed" if args.pad > 0 else "variants.input.bed")
    write_bed3(padded_variant_bed, variant_rows)

    homo_bed = homo_regions_to_bed3(args.homo_regions)
    homo_out = os.path.join(args.output_dir, "homo.bed")
    nonhomo_out = os.path.join(args.output_dir, "nonhomo.bed")

    try:
        split_with_bedtools(padded_variant_bed, homo_bed, homo_out, nonhomo_out)
    finally:
        os.unlink(homo_bed)

    homo_n = sum(1 for _ in open(homo_out))
    nonhomo_n = sum(1 for _ in open(nonhomo_out))
    if homo_n == 0 and nonhomo_n == 0:
        raise SystemExit("Error: split produced no regions; check chromosome naming and BED coordinates.")
    print(f"Wrote {homo_out} ({homo_n} regions) and {nonhomo_out} ({nonhomo_n} regions)", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        print(f"Error: bedtools command failed: {exc}", file=sys.stderr)
        sys.exit(1)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
