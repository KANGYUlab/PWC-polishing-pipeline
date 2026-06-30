#!/usr/bin/env bash
#
# End-to-end PWC entry point: homopolymer calling, variant split, and polish.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PWC_DIR="${SCRIPT_DIR}/pwc"

usage() {
    cat <<'EOF'
Usage: run_pwc.sh [options]

Run the full PWC workflow:
  1) homo.py       — homopolymer regions from reference
  2) split_variants.py — split variants into homo.bed / nonhomo.bed
  3) polish.sh     — reference support + multi-platform consensus

Required:
  -r, --ref FASTA           Reference genome FASTA
  --variants BED            Variant BED to polish
  --short-bam BAM           Short-read BAM
  --hifi-bam BAM            HiFi BAM
  --ont-bam BAM             ONT BAM

Optional:
  -o, --outdir DIR          Output directory (default: pwc_out)
  -p, --threads N           Parallel jobs (default: nproc or 8)
  --pad N                   Pad each variant by N bp on both sides (default: 0)
                            Use 0 for SAS BE BED; use 10 for SNP/small INDEL from variant callers
  --no-normalize-chrom      Disable chr-prefix auto normalization for variant BED
  --skip-homo               Skip homo.py (reuse existing homo_regions.tsv in outdir)
  --skip-split              Skip split (reuse existing homo.bed / nonhomo.bed in outdir)
  --skip-validate           Skip preflight input validation
  -h, --help                Show this help

Examples:
  # SAS base-editing BED — no padding
  bash run_pwc.sh -r ref.fasta --variants sas_be.bed \
      --short-bam short.bam --hifi-bam hifi.bam --ont-bam ont.bam

  # Variant-caller SNP / small INDEL — pad 10 bp
  bash run_pwc.sh -r ref.fasta --variants snp_indel.bed --pad 10 \
      --short-bam short.bam --hifi-bam hifi.bam --ont-bam ont.bam
EOF
}

die() {
    echo "Error: $*" >&2
    exit 1
}

require_positive_int() {
    local name="$1"
    local value="$2"
    [[ "$value" =~ ^[0-9]+$ ]] || die "${name} must be a non-negative integer, got: ${value}"
}

require_positive_int_gt0() {
    local name="$1"
    local value="$2"
    require_positive_int "$name" "$value"
    [[ "$value" -gt 0 ]] || die "${name} must be > 0, got: ${value}"
}

REF_FASTA=""
VARIANTS_BED=""
BAM_SHORT=""
BAM_HIFI=""
BAM_ONT=""
OUTDIR="pwc_out"
N_JOBS="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 8)"
PAD=0
SKIP_HOMO=0
SKIP_SPLIT=0
SKIP_VALIDATE=0
NORMALIZE_CHROM=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--ref) REF_FASTA="$2"; shift 2 ;;
        --variants) VARIANTS_BED="$2"; shift 2 ;;
        --short-bam) BAM_SHORT="$2"; shift 2 ;;
        --hifi-bam) BAM_HIFI="$2"; shift 2 ;;
        --ont-bam) BAM_ONT="$2"; shift 2 ;;
        -o|--outdir) OUTDIR="$2"; shift 2 ;;
        -p|--threads) N_JOBS="$2"; shift 2 ;;
        --pad) PAD="$2"; shift 2 ;;
        --no-normalize-chrom) NORMALIZE_CHROM=0; shift ;;
        --skip-homo) SKIP_HOMO=1; shift ;;
        --skip-split) SKIP_SPLIT=1; shift ;;
        --skip-validate) SKIP_VALIDATE=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown option: $1. Run with -h for help." ;;
    esac
done

for var in REF_FASTA VARIANTS_BED BAM_SHORT BAM_HIFI BAM_ONT; do
  if [[ -z "${!var}" ]]; then
    usage >&2
    die "missing required argument: ${var}"
  fi
done

require_positive_int "pad" "$PAD"
require_positive_int_gt0 "threads" "$N_JOBS"

for tool in python3 samtools bedtools; do
  command -v "$tool" >/dev/null 2>&1 || die "required command not found: ${tool}"
done

for script in homo.py split_variants.py validate_inputs.py chrom_utils.py correct_read.py merge_ref_support.py; do
  [[ -f "${PWC_DIR}/${script}" ]] || die "missing pipeline script: ${PWC_DIR}/${script}"
done
[[ -f "${SCRIPT_DIR}/polish.sh" ]] || die "missing pipeline script: ${SCRIPT_DIR}/polish.sh"

if [[ "$SKIP_VALIDATE" -eq 0 ]]; then
    echo "[run_pwc] Preflight validation" >&2
    python3 "${PWC_DIR}/validate_inputs.py" \
        -r "$REF_FASTA" \
        --variants "$VARIANTS_BED" \
        --short-bam "$BAM_SHORT" \
        --hifi-bam "$BAM_HIFI" \
        --ont-bam "$BAM_ONT"
fi

mkdir -p "$OUTDIR"
HOMO_REGIONS="${OUTDIR}/01_homo_regions/homo_regions.tsv"
SPLIT_DIR="${OUTDIR}/02_homo_nonhomo"

SPLIT_ARGS=(--pad "$PAD" -o "$SPLIT_DIR" --ref-fasta "$REF_FASTA")
if [[ "$NORMALIZE_CHROM" -eq 0 ]]; then
    SPLIT_ARGS+=(--no-normalize-chrom)
fi

if [[ "$SKIP_HOMO" -eq 0 ]]; then
    echo "[run_pwc] Step 1/3: homopolymer regions (homo.py)" >&2
    mkdir -p "$(dirname "$HOMO_REGIONS")"
    python3 "${PWC_DIR}/homo.py" "$REF_FASTA" "$HOMO_REGIONS"
else
    [[ -f "$HOMO_REGIONS" ]] || die "--skip-homo set but homo_regions not found: ${HOMO_REGIONS}"
    echo "[run_pwc] Step 1/3: skipped (using $HOMO_REGIONS)" >&2
fi

if [[ "$SKIP_SPLIT" -eq 0 ]]; then
    echo "[run_pwc] Step 2/3: split variants (split_variants.py)" >&2
    if [[ "$PAD" -gt 0 ]]; then
        echo "[run_pwc] Padding variants by ${PAD} bp on each side (SNP/small INDEL mode)" >&2
    else
        echo "[run_pwc] No padding (SAS BE / pre-expanded BED mode)" >&2
    fi
    if [[ "$NORMALIZE_CHROM" -eq 1 ]]; then
        echo "[run_pwc] Auto-normalizing variant chromosome names to match reference" >&2
    fi
    python3 "${PWC_DIR}/split_variants.py" "$HOMO_REGIONS" "$VARIANTS_BED" "${SPLIT_ARGS[@]}"
else
    [[ -f "${SPLIT_DIR}/homo.bed" && -f "${SPLIT_DIR}/nonhomo.bed" ]] \
        || die "--skip-split set but homo.bed/nonhomo.bed not found in ${SPLIT_DIR}"
    echo "[run_pwc] Step 2/3: skipped (using $SPLIT_DIR)" >&2
fi

echo "[run_pwc] Step 3/3: polish (polish.sh)" >&2
bash "${SCRIPT_DIR}/polish.sh" \
    -r "$REF_FASTA" \
    --short-bam "$BAM_SHORT" \
    --hifi-bam "$BAM_HIFI" \
    --ont-bam "$BAM_ONT" \
    --homo-bed "${SPLIT_DIR}/homo.bed" \
    --nonhomo-bed "${SPLIT_DIR}/nonhomo.bed" \
    -o "$OUTDIR" \
    -p "$N_JOBS"

echo "[run_pwc] Complete. Results: ${OUTDIR}/" >&2
