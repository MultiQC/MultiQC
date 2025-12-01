# Wittyer

Wittyer is a tool for benchmarking structural variant (SV) calls against a truth set, developed by Illumina. This MultiQC module parses the JSON statistics output from Wittyer to provide comprehensive performance metrics.

## Module Output

The module generates three main sections in the MultiQC report:

### 1. General Statistics

Added to the main MultiQC general statistics table:

- **Precision**: Event-level precision (%)
- **Recall**: Event-level recall (%)
- **F-score**: Event-level F-score (%)

### 2. Overall Statistics Table

A detailed table showing both event-level and base-level statistics:

**Event-level metrics:**
- True Positives (TP)
- False Negatives (FN)
- False Positives (FP)
- Recall, Precision, F-score

**Base-level metrics:**
- True Positives (TP) - number of bases
- False Negatives (FN) - number of bases
- False Positives (FP) - number of bases
- Recall, Precision, F-score

### 3. Variant Type Performance

Separate tables for each of the six main structural variant types:

1. **Deletions** - Deletion variants
2. **Insertions** - Insertion variants
3. **Duplications** - Duplication variants
4. **Inversions** - Inversion variants
5. **Copy Number Gains** - Copy number gain events
6. **Copy Number Losses** - Copy number loss events

**Each table shows:**
- True Positives (TP) - Number of correctly identified variants
- False Negatives (FN) - Number of missed true variants
- False Positives (FP) - Number of incorrectly called variants
- Recall (%) - Percentage of true variants detected
- Precision (%) - Percentage of called variants that are correct
- F-score (%) - Harmonic mean of precision and recall

**Benefits of separate tables:**
- Clear focus on each variant type's performance
- Easy comparison across samples for specific variant classes
- Better readability than a single large table
- Simplified analysis and troubleshooting

**Note**: Only tables with available data will be displayed. If a variant type is not present in your Wittyer output, its table will be automatically omitted.

## File Search Pattern

The module looks for JSON files containing the following keys:

- `EventPrecision`
- `EventRecall`
- `PerSampleStats`

## Example Usage

To generate Wittyer statistics for inclusion in MultiQC:

```bash
# Run Wittyer benchmarking
Wittyer \
  -t truth.vcf \
  -i query.vcf \
  -o output_dir \
  -b benchmark.bed

# The Wittyer.Stats.json file will be automatically detected by MultiQC
multiqc output_dir
```

## Interpreting Results

- **Precision**: Proportion of called variants that match the truth set (low false positives)
- **Recall**: Proportion of true variants that were successfully called (low false negatives)
- **F-score**: Harmonic mean of precision and recall, providing a balanced metric

High-quality variant calling should achieve:
- Precision > 90%
- Recall > 90%
- F-score > 90%

Lower metrics may indicate:
- Issues with variant calling parameters
- Technical problems with sequencing data
- Limitations of the calling algorithm for specific variant types

## More Information

For more details about Wittyer, see:

- [Wittyer GitHub Repository](https://github.com/Illumina/Wittyer)
- [Wittyer Publication](https://doi.org/10.1093/bioinformatics/btaa397)
