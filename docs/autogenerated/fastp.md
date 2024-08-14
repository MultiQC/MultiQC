---
name: fastp
urls: ["https://github.com/OpenGene/fastp"]
summary: >
  All-in-one FASTQ preprocessor (QC, adapters, trimming, filtering, splitting...)
---

Fastp goes through fastq files in a folder and perform a series of quality control and filtering.
Quality control and reporting are displayed both before and after filtering, allowing for a clear
depiction of the consequences of the filtering process. Notably, the latter can be conducted on a
variety of parameters including quality scores, length, as well as the presence of adapters, polyG,
or polyX tailing.

By default, the module generates the sample names based on the input FastQ file names in
the command line used by fastp. If you prefer, you can tell the module to use
the filenames as sample names instead. To do so, use the following config option:

```yaml
fastp:
  s_name_filenames: true
```

### File search patterns

```yaml
fastp:
  contents: '"before_filtering": {'
  fn: "*.json"
  num_lines: 50
```
