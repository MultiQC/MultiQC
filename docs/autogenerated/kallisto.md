---
name: Kallisto
urls: ['http://pachterlab.github.io/kallisto/']
summary: >
  Quantifies abundances of transcripts (or more generally, of target sequences) from RNA-Seq data
---

**Note** - MultiQC parses the standard out from Kallisto, _not_ any of its output files
(`abundance.h5`, `abundance.tsv`, and `run_info.json`). As such, you must capture the
Kallisto stdout to a file when running to use the MultiQC module.

### File search patterns

```yaml
kallisto:
  contents: '[quant] finding pseudoalignments for the reads'
```