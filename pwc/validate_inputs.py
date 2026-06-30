#!/usr/bin/env python3
"""Validate PWC pipeline inputs before running run_pwc.sh."""

import argparse
import sys

from chrom_utils import validate_pwc_inputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate inputs for run_pwc.sh / polish.sh.")
    parser.add_argument("-r", "--ref", required=True, help="Reference FASTA")
    parser.add_argument("--variants", required=True, help="Variant BED")
    parser.add_argument("--short-bam", required=True, help="Short-read BAM")
    parser.add_argument("--hifi-bam", required=True, help="HiFi BAM")
    parser.add_argument("--ont-bam", required=True, help="ONT BAM")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    validate_pwc_inputs(
        ref_fasta=args.ref,
        variants_bed=args.variants,
        short_bam=args.short_bam,
        hifi_bam=args.hifi_bam,
        ont_bam=args.ont_bam,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
