---
name: Samtools
urls: ["http://www.htslib.org"]
summary: >
  Toolkit for interacting with BAM/CRAM files
---

Supported commands:

- `stats`
- `flagstats`
- `idxstats`
- `rmdup`
- `coverage`
- `markdup`

#### idxstats

The `samtools idxstats` prints its results to standard out (no consistent file name) and has no header lines
(no way to recognise from content of file). As such, `idxstats` result files must have the string `idxstat`
somewhere in the filename.

There are a few MultiQC config options that you can add to customise how the idxstats module works. A typical
configuration could look as follows:

```yaml
# Always include these chromosomes in the plot
samtools_idxstats_always:
  - X
  - Y

# Never include these chromosomes in the plot
samtools_idxstats_ignore:
  - MT

# Threshold where chromosomes are ignored in the plot.
# Should be a fraction, default is 0.001 (0.1% of total)
samtools_idxstats_fraction_cutoff: 0.001

# Name of the X and Y chromosomes.
# If not specified, MultiQC will search for any chromosome
# names that look like x, y, chrx or chry (case-insensitive search)
samtools_idxstats_xchr: myXchr
samtools_idxstats_ychr: myYchr
```

### File search patterns

```yaml
samtools/coverage:
  contents: "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\t\
    meanmapq"
  num_lines: 10
samtools/flagstat:
  contents: in total (QC-passed reads + QC-failed reads)
samtools/idxstats:
  fn: "*idxstat*"
samtools/markdup_json:
  contents: '"COMMAND": "samtools markdup'
  num_lines: 10
samtools/markdup_txt:
  contents: "COMMAND: samtools markdup"
  num_lines: 2
samtools/rmdup:
  contents: "[bam_rmdup"
samtools/stats:
  contents: This file was produced by samtools stats
```