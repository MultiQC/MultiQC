---
name: VEP
urls: ["https://www.ensembl.org/info/docs/tools/vep/index.html"]
summary: >
  Determines the effect of variants on genes, transcripts and protein sequences, as well as regulatory regions
---

MultiQC parses the Ensembl VEP summary statistics stored in either HTML or plain text format.

Beside VEP's default naming convention, you can run VEP with one of the options below to use this module:

- `--stats_file [OUTPUT_FILENAME]_summary.html` _(VEP's default naming convention)_
- `--stats_file [SAMPLE_NAME].vep.html` _(without the `vep` or `summary` suffix, MultiQC will ignore the HTML files)_
- `--stats_file [SAMPLE_NAME]_vep.html`
- `--stats_text --stats_file [SAMPLE_NAME].vep.txt`
- `--stats_text --stats_file [SAMPLE_NAME]_vep.txt`

See the [VEP](https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#stats)
documentation for more information.

### File search patterns

```yaml
vep/vep_html:
  contents: VEP summary
  fn: "*.html"
  max_filesize: 1000000
  num_lines: 10
vep/vep_txt:
  contents: "[VEP run statistics]"
  max_filesize: 100000
  num_lines: 1
```