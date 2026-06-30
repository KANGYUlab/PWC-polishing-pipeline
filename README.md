# PWC (Platform-integrated Window Consensus)

The PWC (Platform-integrated Window Consensus) is an innovative genome polishing method that integrates data from multiple sequencing platforms to achieve high-quality genome assembly validation and correction. Rather than a standalone software, PWC represents a comprehensive methodology that leverages the complementary strengths of different sequencing technologies.

## Table of contents

* [Overview](#overview)
* [Repository Layout](#repository-layout)
* [Quick Start](#quick-start)
* [Method Description](#method-description)
* [Platform Integration](#platform-integration)
* [Implementation](#implementation)
* [Output and Results](#output-and-results)

## Overview

PWC is designed to address the limitations of single-platform genome polishing approaches by integrating short-read, HiFi, and Oxford Nanopore Technologies (ONT) sequencing data. The method employs adaptive consensus thresholds based on sequence composition to identify and correct assembly errors more effectively than traditional single-platform methods.

## Repository Layout

```
PWC-polishing-pipeline/
├── README.md              # This file
├── consensus.md           # Method details for consensus stage
├── run_pwc.sh             # One-command entry (recommended)
├── polish.sh              # Multi-platform polish stage only
└── pwc/                   # Python utilities
    ├── homo.py            # Homopolymer regions from reference
    ├── split_variants.py  # Split variants → homo.bed / nonhomo.bed
    ├── validate_inputs.py # Preflight checks (called by run_pwc.sh)
    ├── chrom_utils.py     # Chromosome name normalization helpers
    ├── correct_read.py    # Reference-support pileup analysis
    └── merge_ref_support.py
```

Clone the repo, `cd` into it, and run `bash run_pwc.sh` from the repository root.

## Quick Start

### Dependencies

- Python 3 with [Biopython](https://biopython.org/) (`pip install biopython`)
- [samtools](http://www.htslib.org/) (≥1.17, for `samtools consensus`)
- [bedtools](https://bedtools.readthedocs.io/)
- `bc`, `awk`, `bash`

Align short-read, HiFi, and ONT data to the **same reference** and index all BAMs before running PWC.

### One-command workflow

```bash
# SAS base-editing (BE) BED — use directly, no padding
bash run_pwc.sh \
  -r reference.fasta \
  --variants sas_be.bed \
  --short-bam short.bam \
  --hifi-bam hifi.bam \
  --ont-bam ont.bam \
  -o pwc_out

# SNP / small INDEL from variant callers — expand each site by 10 bp
bash run_pwc.sh \
  -r reference.fasta \
  --variants snp_indel.bed \
  --pad 10 \
  --short-bam short.bam \
  --hifi-bam hifi.bam \
  --ont-bam ont.bam \
  -o pwc_out
```

### Step-by-step workflow

```bash
# 1. Homopolymer regions from reference
python3 pwc/homo.py reference.fasta homo_regions.tsv

# 2. Split variants into homo.bed / nonhomo.bed
#    Chromosome names are auto-normalized (chr1 <-> 1, MT <-> chrM) via --ref-fasta
#    SAS BE BED — direct input:
python3 pwc/split_variants.py homo_regions.tsv sas_be.bed --ref-fasta reference.fasta -o 02_homo_nonhomo
#    Variant-caller SNP / small INDEL — pad 10 bp on each side:
python3 pwc/split_variants.py homo_regions.tsv variants.bed --ref-fasta reference.fasta --pad 10 -o 02_homo_nonhomo

# 3. Multi-platform polish
bash polish.sh \
  -r reference.fasta \
  --short-bam short.bam \
  --hifi-bam hifi.bam \
  --ont-bam ont.bam \
  --homo-bed 02_homo_nonhomo/homo.bed \
  --nonhomo-bed 02_homo_nonhomo/nonhomo.bed \
  -o pwc_out
```

### Input BED notes

| Source | Recommendation |
|--------|----------------|
| SAS base-editing (BE) BED | Use as-is (`--pad 0`, default) |
| SNP / small INDEL from DeepVariant, GATK, etc. | Pad **10 bp** on each side (`--pad 10`) |

`run_pwc.sh` runs preflight validation (file existence, BAM index, BED format, chromosome overlap) and auto-normalizes variant BED chromosome names to match the reference (`chr1` ↔ `1`, `MT` ↔ `chrM`). Use `--no-normalize-chrom` to disable.

Homopolymer-overlapping sites are polished with short-read + HiFi only; other sites additionally use ONT.

### Output layout

```
pwc_out/
├── 01_homo_regions/homo_regions.tsv   # from homo.py (run_pwc.sh only)
├── 02_homo_nonhomo/                   # split BEDs + reference-support tables
├── 03_homo_element/                   # homo short-read consensus
├── 04_homo_revio/                     # homo HiFi consensus
├── 05_nonhomo_element/                # nonhomo short-read consensus
├── 06_nonhomo_revio/                  # nonhomo HiFi consensus
└── 07_nonhomo_ont/                    # nonhomo ONT consensus
```

See `consensus.md` for method details and `bash polish.sh -h` / `bash run_pwc.sh -h` for all options.

## Method Description

### Core Philosophy

PWC operates on the principle that **different sequencing platforms have complementary strengths**:
- **Short-read data**: High accuracy for base-level corrections
- **HiFi long-reads**: Improved contiguity and structural accuracy
- **ONT long-reads**: Enhanced detection of structural variants and complex regions

### Multi-Platform Integration Strategy

The method integrates these platforms sequentially rather than supporting them independently:

1. **Reference Support Analysis**: Identifies regions requiring consensus analysis
2. **Sequential Platform Processing**: Processes data through each platform in order
3. **Adaptive Threshold Selection**: Adjusts consensus criteria based on sequence composition
4. **Cross-Platform Validation**: Ensures consistency across all platforms

## Platform Integration

### Sequential Processing Approach

PWC processes genomic regions through three platforms in sequence:

```
Short-read Analysis → HiFi Analysis → ONT Analysis
```

**Key Integration Features:**
- **Sequential Refinement**: Each platform builds upon the previous platform's results
- **Shared Region Analysis**: Same genomic regions analyzed across all platforms
- **Consensus Propagation**: Results from one platform inform the next platform's analysis
- **Unified Output**: Single comprehensive result integrating all platform data

### Platform-Specific Roles

- **Short-read**: Initial consensus generation and error identification
- **HiFi**: Refinement of consensus using high-fidelity long-read data
- **ONT**: Final validation and detection of complex structural variations

## Implementation

### Reference Support Analysis

Before consensus generation, PWC performs reference support analysis to identify regions requiring attention:

```bash
# Regions classified as needing consensus analysis if they fail:
# Short-read: Coverage ≤2*mean_depth, ≥10 reads, support ratio ≥55%
# HiFi/ONT: ≥10 reads, support ratio ≥55%
```

### Adaptive Consensus Thresholds

PWC employs intelligent threshold selection based on sequence composition:

- **Standard thresholds**: Consensus ≥70%, Quality ≥70 (normal composition)
- **Relaxed thresholds**: Consensus ≥50%, Quality ≥50 (AG/CT ≥90% bias)

### Consensus Generation

```bash
# Multi-platform consensus generation
for platform in shortread hifi ont; do
    samtools consensus -a -c "$threshold" -r "$region" "$BAM_FILE" \
        -m simple -d 10 --show-del yes -f pileup > "$pileup_file"
done
```

## Output and Results

### Consensus Classification

PWC classifies results into four categories per platform:

1. **Normal Reference**: Standard thresholds, consensus matches reference
2. **Normal Non-reference**: Standard thresholds, consensus differs from reference
3. **Relaxed Reference**: Relaxed thresholds, consensus matches reference
4. **Relaxed Non-reference**: Relaxed thresholds, consensus differs from reference

### Output Structure

```
├── consensus.success.normal.txt      # Standard threshold results
├── consensus.success.relaxed.txt     # Relaxed threshold results  
├── consensus.fail.txt               # Failed consensus generation
├── consensus.success.normal.ref.txt  # Reference-matching sequences
├── consensus.success.normal.nonref.txt # Non-reference sequences
├── consensus.success.relaxed.ref.txt # Relaxed reference-matching
├── consensus.success.relaxed.nonref.txt # Relaxed non-reference
└── consensus.combined.bed           # Regions for next platform
```

### Quality Metrics

- **Consensus Threshold**: 70% (normal) / 50% (relaxed)
- **Quality Threshold**: 70 (normal) / 50 (relaxed)  
- **Composition Bias**: AG or CT frequency ≥90% triggers relaxed thresholds
- **Coverage Requirements**: Minimum 10 reads per position

## Implementation Requirements

### Dependencies

- **Samtools**: For consensus generation and BAM processing
- **Bedtools**: For genomic region manipulation
- **Bash scripting**: For pipeline orchestration
- **Parallel processing**: For efficient large-scale analysis

### Input Data

- **Short-read BAM files**: High-accuracy base-level data
- **HiFi BAM files**: Long-read high-fidelity data
- **ONT BAM files**: Long-read structural variant data
- **Reference genome**: Assembly to be polished

## Citation

If you use the PWC method in your research, please cite:

Yanan Chu, Zhuo Huang, Changjun Shao, Shuming Guo, Xinyao Yu, Jian Wang, Yabin Tian, Jing Chen, Ran Li, Yukun He, Jun Yu, Jie Huang, Zhancheng Gao, Yu Kang. Approaching an Error-Free Diploid Human Genome. bioRxiv. doi: https://doi.org/10.1101/2025.08.01.667781
