---
name: Prokka
urls: ["http://www.vicbioinformatics.com/software.prokka.shtml"]
summary: >
  Rapid annotation of prokaryotic genomes
---

The Prokka module accepts two configuration options:

- `prokka_table`: default `False`. Show a table in the report.
- `prokka_barplot`: default `True`. Show a barplot in the report.
- `prokka_fn_snames`: default `False`. Use filenames for sample names (see below).

Sample names are generated using the first line in the prokka reports:

```
organism: Helicobacter pylori Sample1
```

The module assumes that the first two words are the organism name and
the third is the sample name. So the above will give a sample name of
`Sample1`.

If you prefer, you can set `config.prokka_fn_snames` to `True` and MultiQC
will instead use the log filename as the sample name.

### File search patterns

```yaml
prokka:
  contents: "contigs:"
  num_lines: 2
```
