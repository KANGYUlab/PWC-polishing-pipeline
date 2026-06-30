#!/usr/bin/env bash
#
# PWC multi-platform polishing pipeline.
#
# Prerequisites:
#   1. Run homo.py on the reference to get homopolymer regions.
#   2. Run split_variants.py to produce homo.bed and nonhomo.bed.
#   3. Run this script with aligned BAMs and the split BED files.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PWC_DIR="${SCRIPT_DIR}/pwc"
CORRECT_READ="${PWC_DIR}/correct_read.py"
MERGE_REF_SUPPORT="${PWC_DIR}/merge_ref_support.py"

usage() {
    cat <<'EOF'
Usage: polish.sh [options]

Run the PWC reference-support and multi-platform consensus pipeline.

Required:
  -r, --ref FASTA           Reference genome FASTA
  --short-bam BAM           Short-read BAM (indexed)
  --hifi-bam BAM            HiFi BAM (indexed)
  --ont-bam BAM             ONT BAM (indexed)
  --homo-bed BED            Homopolymer-overlapping variant regions (from split_variants.py)
  --nonhomo-bed BED         Non-homopolymer variant regions (from split_variants.py)

Optional:
  -o, --outdir DIR          Output directory (default: pwc_out)
  -p, --threads N           Parallel jobs (default: nproc or 8)
  -h, --help                Show this help

Upstream preparation (run once before polish.sh):

  # 1) Homopolymer regions from reference
  python3 pwc/homo.py reference.fasta homo_regions.tsv

  # 2) Split variants into homo.bed / nonhomo.bed
  #
  # SAS base-editing (BE) BED can be used directly:
  python3 pwc/split_variants.py homo_regions.tsv sas_be.bed -o 02_homo_nonhomo
  #
  # SNP / small INDEL from variant callers — pad 10 bp on each side:
  python3 pwc/split_variants.py homo_regions.tsv variants.bed --pad 10 -o 02_homo_nonhomo

  # 3) Polish
  bash polish.sh -r reference.fasta \
      --short-bam short.bam --hifi-bam hifi.bam --ont-bam ont.bam \
      --homo-bed 02_homo_nonhomo/homo.bed \
      --nonhomo-bed 02_homo_nonhomo/nonhomo.bed \
      -o pwc_out

Dependencies: python3, samtools, bedtools, bc, awk
EOF
}

REF_FASTA=""
BAM_SHORT=""
BAM_HIFI=""
BAM_ONT=""
HOMO_BED=""
NONHOMO_BED=""
OUTDIR="pwc_out"
N_JOBS="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 8)"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--ref) REF_FASTA="$2"; shift 2 ;;
        --short-bam) BAM_SHORT="$2"; shift 2 ;;
        --hifi-bam) BAM_HIFI="$2"; shift 2 ;;
        --ont-bam) BAM_ONT="$2"; shift 2 ;;
        --homo-bed) HOMO_BED="$2"; shift 2 ;;
        --nonhomo-bed) NONHOMO_BED="$2"; shift 2 ;;
        -o|--outdir) OUTDIR="$2"; shift 2 ;;
        -p|--threads) N_JOBS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

for var in REF_FASTA BAM_SHORT BAM_HIFI BAM_ONT HOMO_BED NONHOMO_BED; do
  if [[ -z "${!var}" ]]; then
    echo "Error: missing required argument for ${var}" >&2
    usage >&2
    exit 1
  fi
done

for f in "$REF_FASTA" "$BAM_SHORT" "$BAM_HIFI" "$BAM_ONT" "$HOMO_BED" "$NONHOMO_BED"; do
  if [[ ! -f "$f" ]]; then
    echo "Error: file not found: $f" >&2
    exit 1
  fi
done

for tool in python3 samtools bedtools bc awk; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Error: required command not found: $tool" >&2
    exit 1
  fi
done

if [[ ! -f "$CORRECT_READ" ]]; then
  echo "Error: correct_read.py not found at $CORRECT_READ" >&2
  exit 1
fi

HOMO_NONHOMO_DIR="${OUTDIR}/02_homo_nonhomo"
DIR_HOMO_ELEM="${OUTDIR}/03_homo_element"
DIR_HOMO_REVIO="${OUTDIR}/04_homo_revio"
DIR_NONHOMO_ELEM="${OUTDIR}/05_nonhomo_element"
DIR_NONHOMO_REVIO="${OUTDIR}/06_nonhomo_revio"
DIR_NONHOMO_ONT="${OUTDIR}/07_nonhomo_ont"

mkdir -p "$HOMO_NONHOMO_DIR" "$DIR_HOMO_ELEM" "$DIR_HOMO_REVIO" \
         "$DIR_NONHOMO_ELEM" "$DIR_NONHOMO_REVIO" "$DIR_NONHOMO_ONT"

HOMO_BED_LINK="${HOMO_NONHOMO_DIR}/homo.bed"
NONHOMO_BED_LINK="${HOMO_NONHOMO_DIR}/nonhomo.bed"
cp -f "$HOMO_BED" "$HOMO_BED_LINK"
cp -f "$NONHOMO_BED" "$NONHOMO_BED_LINK"

echo "[PWC] Step 1/3: reference support analysis" >&2

python3 "$CORRECT_READ" "$HOMO_BED_LINK" "$BAM_SHORT" "$REF_FASTA" "${DIR_HOMO_ELEM}/homo.elem.txt" -p "$N_JOBS"
python3 "$CORRECT_READ" "$HOMO_BED_LINK" "$BAM_HIFI" "$REF_FASTA" "${DIR_HOMO_REVIO}/homo.revio.txt" -p "$N_JOBS"
python3 "$CORRECT_READ" "$NONHOMO_BED_LINK" "$BAM_SHORT" "$REF_FASTA" "${DIR_NONHOMO_ELEM}/nonhomo.elem.txt" -p "$N_JOBS"
python3 "$CORRECT_READ" "$NONHOMO_BED_LINK" "$BAM_HIFI" "$REF_FASTA" "${DIR_NONHOMO_REVIO}/nonhomo.revio.txt" -p "$N_JOBS"
python3 "$CORRECT_READ" "$NONHOMO_BED_LINK" "$BAM_ONT" "$REF_FASTA" "${DIR_NONHOMO_ONT}/nonhomo.ont.txt" -p "$N_JOBS"

python3 "$MERGE_REF_SUPPORT" "$HOMO_BED_LINK" \
    "${DIR_HOMO_ELEM}/homo.elem.txt" "${DIR_HOMO_REVIO}/homo.revio.txt" \
    -o "${HOMO_NONHOMO_DIR}/homo_ref_support.txt"

python3 "$MERGE_REF_SUPPORT" "$NONHOMO_BED_LINK" \
    "${DIR_NONHOMO_ELEM}/nonhomo.elem.txt" \
    "${DIR_NONHOMO_REVIO}/nonhomo.revio.txt" \
    "${DIR_NONHOMO_ONT}/nonhomo.ont.txt" \
    -o "${HOMO_NONHOMO_DIR}/nonhomo_ref_support.txt"

awk '{
  c4 = $4 + 0; c5 = $5 + 0; c6 = $6 + 0; c7 = $7 + 0;
  ok4 = (c4 >= 10 && c5 > 0 && c4/c5 >= 0.55 && c5 <= 300);
  ok6 = (c6 >= 10 && c7 > 0 && c6/c7 >= 0.55);
  if (!(ok4 || ok6)) print
}' "${HOMO_NONHOMO_DIR}/homo_ref_support.txt" > "${HOMO_NONHOMO_DIR}/homo_unsupport.bed"

awk '{
  c4 = $4 + 0; c5 = $5 + 0; c6 = $6 + 0; c7 = $7 + 0;
  ok4 = (c4 >= 10 && c5 > 0 && c4/c5 >= 0.55 && c5 <= 300);
  ok6 = (c6 >= 10 && c7 > 0 && c6/c7 >= 0.55);
  if (NF >= 9) {
    c8 = $8 + 0; c9 = $9 + 0;
    ok8 = (c8 >= 10 && c9 > 0 && c8/c9 >= 0.55);
  } else {
    ok8 = 1;
  }
  if (!(ok4 || ok6 || ok8)) print
}' "${HOMO_NONHOMO_DIR}/nonhomo_ref_support.txt" > "${HOMO_NONHOMO_DIR}/nonhomo_unsupport.bed"

echo "[PWC] Step 2/3: homopolymer-region consensus (short-read + HiFi)" >&2

CONSENSUS_NORMAL="0.70"
QUALITY_NORMAL="70"
CONSENSUS_RELAXED="0.50"
QUALITY_RELAXED="50"
MIN_COVERAGE=10
SAMTOOLS_DEPTH=10

process_region_consensus() {
    local line="$1"
    [[ -z "$line" ]] && return 0

    local chrom start end
    read -r chrom start end <<<"$line"
    [[ -z "$chrom" ]] && return 0

    local region="${chrom}:$((start + 1))-${end}"
    local tmp_dir
    tmp_dir=$(mktemp -d)
    local tmp_bed="${tmp_dir}/region.bed"
    local pileup_file="${tmp_dir}/pileup.txt"
    local part_file="${CONSENSUS_PARTS_DIR}/${$}_${RANDOM}.txt"

    printf "%s\t%s\t%s\n" "$chrom" "$start" "$end" >"$tmp_bed"

    local ref_seq
    ref_seq=$(bedtools getfasta -fi "$REF_FASTA" -bed "$tmp_bed" -fo - 2>/dev/null | sed 1d | tr -d '\n' | tr '[:lower:]' '[:upper:]')
    if [[ -z "$ref_seq" ]]; then
        ref_seq="ERROR"
    fi

    local seq_len=${#ref_seq}
    local ag_freq="NA"
    local ct_freq="NA"
    local consensus_threshold="$CONSENSUS_NORMAL"
    local quality_threshold="$QUALITY_NORMAL"
    local mode="normal"

    if [[ "$ref_seq" == "ERROR" || $seq_len -eq 0 ]]; then
        printf "%s\t%s\t%s\tFAIL\t%s\t%s\t%s\n" "$chrom" "$start" "$end" "$ref_seq" "$ag_freq" "$ct_freq" >"${part_file}.fail"
        rm -rf "$tmp_dir"
        return 0
    fi

    local a_count g_count c_count t_count
    a_count=$(printf '%s' "$ref_seq" | tr -cd 'A' | wc -c | tr -d '[:space:]')
    g_count=$(printf '%s' "$ref_seq" | tr -cd 'G' | wc -c | tr -d '[:space:]')
    c_count=$(printf '%s' "$ref_seq" | tr -cd 'C' | wc -c | tr -d '[:space:]')
    t_count=$(printf '%s' "$ref_seq" | tr -cd 'T' | wc -c | tr -d '[:space:]')

    ag_freq=$(echo "scale=4; ($a_count + $g_count) / $seq_len" | bc -l)
    ct_freq=$(echo "scale=4; ($c_count + $t_count) / $seq_len" | bc -l)

    if [[ $(echo "$ag_freq >= 0.9" | bc -l) -eq 1 || $(echo "$ct_freq >= 0.9" | bc -l) -eq 1 ]]; then
        consensus_threshold="$CONSENSUS_RELAXED"
        quality_threshold="$QUALITY_RELAXED"
        mode="relaxed"
    fi

    if ! samtools consensus -a -c "$consensus_threshold" -r "$region" "$BAM_FILE" \
        -m simple -d "$SAMTOOLS_DEPTH" --show-del yes -f pileup >"$pileup_file" 2>/dev/null; then
        printf "%s\t%s\t%s\tFAIL\t%s\t%s\t%s\n" "$chrom" "$start" "$end" "$ref_seq" "$ag_freq" "$ct_freq" >"${part_file}.fail"
        rm -rf "$tmp_dir"
        return 0
    fi

    if ! awk -v q="$quality_threshold" -v d="$MIN_COVERAGE" 'NF>=6 {if($4 < d || $6 < q){exit 1}} END{if(NR==0) exit 1}' "$pileup_file"; then
        printf "%s\t%s\t%s\tFAIL\t%s\t%s\t%s\n" "$chrom" "$start" "$end" "$ref_seq" "$ag_freq" "$ct_freq" >"${part_file}.fail"
        rm -rf "$tmp_dir"
        return 0
    fi

    local consensus
    consensus=$(awk '{printf "%s", $5}' "$pileup_file")
    if [[ -z "$consensus" ]]; then
        printf "%s\t%s\t%s\tFAIL\t%s\t%s\t%s\n" "$chrom" "$start" "$end" "$ref_seq" "$ag_freq" "$ct_freq" >"${part_file}.fail"
        rm -rf "$tmp_dir"
        return 0
    fi

    local output_line
    output_line=$(printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$chrom" "$start" "$end" "$consensus" "$ref_seq" "$ag_freq" "$ct_freq")

    if [[ "$mode" == "normal" ]]; then
        printf "%s\n" "$output_line" >"${part_file}.normal"
    else
        printf "%s\n" "$output_line" >"${part_file}.relaxed"
    fi

    rm -rf "$tmp_dir"
}

export -f process_region_consensus

process_platform() {
    local platform="$1"
    local input_bed="$2"
    local bam_file="$3"
    local output_dir="$4"

    mkdir -p "$output_dir"

    local success_normal_file="$output_dir/consensus.success.normal.txt"
    local success_relaxed_file="$output_dir/consensus.success.relaxed.txt"
    local fail_file="$output_dir/consensus.fail.txt"
    local normal_ref_file="$output_dir/consensus.success.normal.ref.txt"
    local normal_nonref_file="$output_dir/consensus.success.normal.nonref.txt"
    local relaxed_ref_file="$output_dir/consensus.success.relaxed.ref.txt"
    local relaxed_nonref_file="$output_dir/consensus.success.relaxed.nonref.txt"
    local combined_bed="$output_dir/consensus.combined.bed"
    local parts_dir="$output_dir/consensus.parts"

    rm -rf "$parts_dir"
    mkdir -p "$parts_dir"

    : >"$success_normal_file"
    : >"$success_relaxed_file"
    : >"$fail_file"

    if [[ ! -s "$input_bed" ]]; then
        printf "[%s] No regions to process (input: %s)\n" "$platform" "$input_bed" >&2
        : >"$normal_ref_file"
        : >"$normal_nonref_file"
        : >"$relaxed_ref_file"
        : >"$relaxed_nonref_file"
        : >"$combined_bed"
        echo "$combined_bed"
        return 0
    fi

    export BAM_FILE="$bam_file"
    export REF_FASTA
    export SAMTOOLS_DEPTH
    export MIN_COVERAGE
    export CONSENSUS_NORMAL
    export QUALITY_NORMAL
    export CONSENSUS_RELAXED
    export QUALITY_RELAXED
    export CONSENSUS_PARTS_DIR="$parts_dir"

    printf "[%s] Processing regions from %s with %s\n" "$platform" "$input_bed" "$bam_file" >&2

    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" ]] && continue
        printf "%s\0" "$line"
    done <"$input_bed" | xargs -0 -P "$N_JOBS" -I {} bash -c 'process_region_consensus "$@"' _ {}

    if compgen -G "${parts_dir}/*.normal" >/dev/null; then
        cat "${parts_dir}"/*.normal >"$success_normal_file"
    fi
    if compgen -G "${parts_dir}/*.relaxed" >/dev/null; then
        cat "${parts_dir}"/*.relaxed >"$success_relaxed_file"
    fi
    if compgen -G "${parts_dir}/*.fail" >/dev/null; then
        cat "${parts_dir}"/*.fail >"$fail_file"
    fi
    rm -rf "$parts_dir"

    if [[ -s "$success_normal_file" ]]; then
        awk 'BEGIN{OFS="\t"} {a=$4; b=$5; gsub(/\*/, "", a); gsub(/\*/, "", b); if(a == b) print $0}' "$success_normal_file" >"$normal_ref_file"
        awk 'BEGIN{OFS="\t"} {a=$4; b=$5; gsub(/\*/, "", a); gsub(/\*/, "", b); if(a != b) print $0}' "$success_normal_file" >"$normal_nonref_file"
    else
        : >"$normal_ref_file"
        : >"$normal_nonref_file"
    fi

    if [[ -s "$success_relaxed_file" ]]; then
        awk 'BEGIN{OFS="\t"} {a=$4; b=$5; gsub(/\*/, "", a); gsub(/\*/, "", b); if(a == b) print $0}' "$success_relaxed_file" >"$relaxed_ref_file"
        awk 'BEGIN{OFS="\t"} {a=$4; b=$5; gsub(/\*/, "", a); gsub(/\*/, "", b); if(a != b) print $0}' "$success_relaxed_file" >"$relaxed_nonref_file"
    else
        : >"$relaxed_ref_file"
        : >"$relaxed_nonref_file"
    fi

    {
        if [[ -s "$fail_file" ]]; then
            awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$fail_file"
        fi
        if [[ -s "$normal_nonref_file" ]]; then
            awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$normal_nonref_file"
        fi
        if [[ -s "$relaxed_nonref_file" ]]; then
            awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$relaxed_nonref_file"
        fi
    } | sort -u >"$combined_bed" || : >"$combined_bed"

    printf "[%s] Results saved in %s\n" "$platform" "$output_dir" >&2
    echo "$combined_bed"
}

run_consensus_track() {
    local initial_bed="$1"
    local short_input_bed="$2"
    local short_skip_bed="$3"
    local elem_dir="$4"
    local revio_dir="$5"
    local ont_dir="$6"
    local bam_short="$7"
    local bam_hifi="$8"
    local bam_ont="$9"
    local include_ont="${10}"

    mkdir -p "$(dirname "$initial_bed")"
    : >"$short_input_bed"
    : >"$short_skip_bed"
    if [[ -s "$initial_bed" ]]; then
        awk -v keep="$short_input_bed" -v skip="$short_skip_bed" \
            'BEGIN{FS=OFS="\t"} {if (NF>=5 && $5+0> 300) print > skip; else print > keep}' "$initial_bed"
    fi

    local next_bed
    next_bed=$(process_platform "Short-read" "$short_input_bed" "$bam_short" "$elem_dir")

    local hifi_input_bed="$next_bed"
    if [[ -s "$short_skip_bed" ]]; then
        mkdir -p "$revio_dir"
        hifi_input_bed="${revio_dir}/consensus.input.bed"
        cat "$next_bed" "$short_skip_bed" | sort -u >"$hifi_input_bed"
    fi

    next_bed=$(process_platform "HiFi" "$hifi_input_bed" "$bam_hifi" "$revio_dir")

    if [[ "$include_ont" == "yes" ]]; then
        process_platform "ONT" "$next_bed" "$bam_ont" "$ont_dir" >/dev/null
    fi
}

run_consensus_track \
    "${HOMO_NONHOMO_DIR}/homo_unsupport.bed" \
    "${DIR_HOMO_ELEM}/element.short_consensus.bed" \
    "${DIR_HOMO_ELEM}/element.skip_short.bed" \
    "$DIR_HOMO_ELEM" "$DIR_HOMO_REVIO" "$DIR_NONHOMO_ONT" \
    "$BAM_SHORT" "$BAM_HIFI" "$BAM_ONT" "no"

echo "[PWC] Step 3/3: non-homopolymer-region consensus (short-read + HiFi + ONT)" >&2

run_consensus_track \
    "${HOMO_NONHOMO_DIR}/nonhomo_unsupport.bed" \
    "${DIR_NONHOMO_ELEM}/element.short_consensus.bed" \
    "${DIR_NONHOMO_ELEM}/element.skip_short.bed" \
    "$DIR_NONHOMO_ELEM" "$DIR_NONHOMO_REVIO" "$DIR_NONHOMO_ONT" \
    "$BAM_SHORT" "$BAM_HIFI" "$BAM_ONT" "yes"

echo "[PWC] Done. Outputs written under ${OUTDIR}/" >&2
echo "  Reference support: ${HOMO_NONHOMO_DIR}/" >&2
echo "  Homo consensus:      ${DIR_HOMO_ELEM}/ and ${DIR_HOMO_REVIO}/" >&2
echo "  Nonhomo consensus: ${DIR_NONHOMO_ELEM}/, ${DIR_NONHOMO_REVIO}/, ${DIR_NONHOMO_ONT}/" >&2
