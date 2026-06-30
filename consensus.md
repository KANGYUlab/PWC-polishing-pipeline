# Multi-platform Consensus Analysis for Assembly Validation

## Overview
We implemented a multi-platform consensus analysis pipeline using short-read, HiFi, and ONT sequencing data. The pipeline employs adaptive consensus thresholds based on sequence composition to identify assembly errors.

## Method Description

### Reference Support Analysis
Prior to consensus generation, we performed reference support analysis using `correct_read.py` to identify regions with sufficient read support from multiple sequencing platforms. Regions were classified as having adequate support if they met the following criteria from any platform:
- Short-read: Coverage ≤2*mean_depth, ≥10 reads, support ratio ≥55%
- HiFi/ONT: ≥10 reads, support ratio ≥55%

**Only regions that failed to meet these reference support criteria were subjected to consensus analysis**, as regions with adequate support were considered validated and did not require further analysis.

### Consensus Generation Strategy
For each genomic region identified as potentially erroneous (i.e., regions lacking sufficient reference support), we generated consensus sequences from three sequencing platforms using the following approach:

1. **Reference Sequence Extraction**: Extract reference sequences using bedtools getfasta
2. **Sequence Composition Analysis**: Calculate AG/CT frequencies to determine consensus thresholds
3. **Adaptive Threshold Selection**: Apply relaxed thresholds (consensus ≥50%, quality ≥50) for regions with ≥90% AG or CT content; standard thresholds (consensus ≥70%, quality ≥70) for normal composition regions
4. **Quality Assessment**: Generate consensus using samtools consensus with quality filtering

### Platform-Specific Processing
The pipeline processes regions sequentially through three platforms:
- **Short-read Analysis**: Initial consensus generation from short-read data
- **HiFi Analysis**: Refined consensus using HiFi long-read data  
- **ONT Analysis**: Final consensus validation using ONT long-read data

### Consensus Classification
For each platform, consensus sequences were classified into four categories:
1. **Normal Reference**: Standard thresholds, consensus matches reference
2. **Normal Non-reference**: Standard thresholds, consensus differs from reference  
3. **Relaxed Reference**: Relaxed thresholds, consensus matches reference
4. **Relaxed Non-reference**: Relaxed thresholds, consensus differs from reference

### Implementation Details

#### Sequence Composition Analysis
```bash
# Calculate base frequencies
seq_len=${#ref_seq}
a_count=$(echo "$ref_seq" | tr -cd 'A' | wc -c)
t_count=$(echo "$ref_seq" | tr -cd 'T' | wc -c) 
g_count=$(echo "$ref_seq" | tr -cd 'G' | wc -c)
c_count=$(echo "$ref_seq" | tr -cd 'C' | wc -c)

# Calculate AG/CT ratios
ag_freq=$(echo "scale=4; ($a_count + $g_count) / $seq_len" | bc -l)
ct_freq=$(echo "scale=4; ($c_count + $t_count) / $seq_len" | bc -l)

# Apply relaxed thresholds for biased composition
if (( $(echo "$ag_freq >= 0.9" | bc -l) )) || (( $(echo "$ct_freq >= 0.9" | bc -l) )); then
    consensus_threshold=0.50
    quality_threshold=50
fi
```

#### Consensus Generation
```bash
# Generate consensus using samtools
samtools consensus -a -c "$consensus_threshold" -r "$region" "$BAM_FILE" \
    -m simple -d 10 --show-del yes -f pileup > "$pileup_file"

# Quality assessment and consensus extraction
if awk -v thresh="$quality_threshold" '{if($6<thresh){exit 1}}' "$pileup_file"; then
    consensus=$(awk '{printf "%s", $5}' "$pileup_file")
    # Output with metadata
    echo -e "${chrom}\t${start}\t${end}\t${consensus}\t${ref_seq}\t${ag_freq}\t${ct_freq}"
else
    echo -e "${chrom}\t${start}\t${end}\tFAIL\t${ref_seq}\t${ag_freq}\t${ct_freq}"
fi
```

#### Parallel Processing
```bash
# Process regions in parallel using xargs
N_JOBS=$(nproc 2>/dev/null || echo 64)
cat input.bed | while read -r line; do
    printf "%s\0" "$line"
done | xargs -0 -P "$N_JOBS" -I {} bash -c 'process_consensus_and_ref "$@"' _ {}
```

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

### Validation Criteria
1. **Reference Agreement**: Consensus sequences matching reference indicate correct assembly
2. **Platform Consistency**: Agreement across multiple platforms strengthens validation
3. **Quality Metrics**: High-quality consensus sequences provide reliable validation
4. **Composition Awareness**: Adaptive thresholds account for sequence bias

This multi-platform consensus approach provides robust assembly validation by leveraging the complementary strengths of different sequencing technologies while accounting for platform-specific biases and sequence composition effects.
