# bamUtil clipOverlap

**bamUtil clipOverlap** is a part of the bamUtil toolkit that clips overlapping ends of paired-end reads to avoid counting the same evidence twice during variant calling.

The clipOverlap utility detects cases where paired-end reads overlap and clips one of the reads at the overlapping region. This is particularly useful for variant calling pipelines where double-counting the same region could lead to artificially inflated confidence in a variant.

This module parses the output of the bamUtil clipOverlap tool when run with the `--stats` flag, typically found in stderr logs.

Example command:

```bash
/path/to/bamUtil/bin/bam clipOverlap --in input.bam --out output.bam --stats
```

The MultiQC module shows:

- The number of overlapping read pairs detected
- The average number of reference bases that overlapped in those pairs
- Statistics on forward vs reverse strand clipping when available
- Information about orientation-based additional clipping
- Variance of reference bases overlapped
- Missing overlapping mate information (if warnings were present)

**General Statistics Table Columns:**

- **Overlapping Pairs** - Number of overlapping pairs detected
- **Avg Overlap** - Average number of reference bases overlapped

**Hidden General Statistics Table Columns:**

- **Fwd Clipped** - Number of times the forward strand was clipped
- **Rev Clipped** - Number of times the reverse strand was clipped
- **Orientation Clip** - Number of times orientation causes additional clipping
- **Overlap Variance** - Variance of reference bases overlapped
- **Missing Mates** - Number of reads where expected overlapping mates were not found

When multiple samples are present, an additional **Clipping Distribution** bar plot is shown to visualize the differences in forward vs. reverse strand clipping across samples.

### Usage Note

The module searches for files containing "Overlap Statistics:" and may use the log filename as the sample name since the bamUtil output does not typically include a sample identifier in the log contents.
