# MultiQC Configuration Reference

This document describes all configuration options available in MultiQC.

## Introduction

MultiQC configuration can be set in several ways:

1. **Command line parameters** - Command line flags are available for many options (run `multiqc --help` to see all available options)
2. **Configuration files** - MultiQC looks for configuration files in the following locations (in order of precedence):
   - `<current working directory>/multiqc_config.yaml`
   - `~/.multiqc_config.yaml`
   - `<installation_dir>/multiqc/utils/config_defaults.yaml`
3. **Environment variables** - MultiQC checks for environment variables that match configuration options prefixed with `MULTIQC_`, for example: `MULTIQC_TITLE="My Report"`

Configuration values are loaded in the following order of precedence (highest to lowest):

1. Command line parameters
2. Current working directory config file
3. User home directory config file
4. Environment variables
5. Default configuration values

The options below can be specified in your YAML configuration files.
For boolean options, use `true` or `false` (all lowercase) in your YAML files.

:::tip

If you'd rather build your config visually, the [Config Wizard](https://seqera.io/multiqc_config_wizard) renders every option below as a form field with the same descriptions and defaults, and validates as you type.

:::

## Report Meta

### Header text

#### `title`

**Type**: <code>str</code>

Title shown at the top of the report and used in the page title.

#### `subtitle`

**Type**: <code>str</code>

Subtitle shown under the report title. Plain text only.

#### `intro_text`

**Type**: <code>str</code>

Paragraph shown under the title. Useful for adding context about the analysis.

#### `report_comment`

**Type**: <code>str</code>

Free-text comment shown at the top of the report. HTML is allowed.

**Example**:

```yaml
report_comment: This report was generated from the RNA-seq pipeline on 2024-08-21.
```

#### `report_header_info`

**Type**: <code>List[Dict[str, str]]</code>

Extra key/value pairs shown in the report header, eg. contact name, run ID, pipeline version. Each list item is a single-key dictionary.

**Example**:

```yaml
report_header_info:
  - Contact E-mail: phil.ewels@seqera.io
  - Application Type: RNA-seq
  - Project Type: Application
  - Sequencing Platform: HiSeq 2500 High Output V4
```

### Report generation info

#### `show_analysis_paths`

**Type**: <code>bool</code> (default: `true`)

Show the absolute paths of analysed directories in the report header.

#### `show_analysis_time`

**Type**: <code>bool</code> (default: `true`)

Show the date and time the report was generated in the header.

## Report Appearance

### Template

#### `template`

**Type**: <code>str</code> (default: `"default"`)

Name of the report template. Built-in templates: default, original, simple, sections, gathered, geo, disco. Plugin packages can register additional templates via the `multiqc.templates.v1` entry point.

**Example**:

```yaml
template: default
```

#### `template_dark_mode`

**Type**: <code>bool</code> (default: `true`)

Enable the dark mode toggle in the report template.

#### `simple_output`

**Type**: <code>bool</code> (default: `false`)

Render a minimal HTML report without the toolbox or interactive widgets. Useful for very large reports.

### Logo

#### `custom_logo`

**Type**: <code>str</code>

Path to an image to show at the top of the report, replacing the MultiQC logo.

**Examples**:

```yaml
custom_logo: /path/to/logo.png
```

```yaml
custom_logo: ./assets/logo.svg
```

#### `custom_logo_dark`

**Type**: <code>str</code>

Path to an alternative logo for dark mode. Falls back to custom_logo if unset.

**Example**:

```yaml
custom_logo_dark: ./assets/logo_dark.svg
```

#### `custom_logo_url`

**Type**: <code>str</code>

URL the custom logo links to when clicked.

**Example**:

```yaml
custom_logo_url: https://www.scilifelab.se
```

#### `custom_logo_title`

**Type**: <code>str</code>

Tooltip text shown when hovering over the custom logo.

**Example**:

```yaml
custom_logo_title: Our institute name
```

#### `custom_logo_width`

**Type**: <code>int</code>

Logo width in pixels. Height scales proportionally.

**Example**:

```yaml
custom_logo_width: 200
```

### Branding

#### `custom_favicon`

**Type**: <code>str</code>

Path to a custom favicon image to show in the browser tab.

**Examples**:

```yaml
custom_favicon: /path/to/favicon.ico
```

```yaml
custom_favicon: ./assets/favicon.png
```

#### `custom_css_files`

**Type**: <code>List[str]</code>

Paths to additional CSS files to inline into the report. Useful for branding overrides.

**Example**:

```yaml
custom_css_files:
  - ./assets/custom.css
  - /path/to/branding.css
```

## Report Contents

### Custom content

#### `custom_content`

**Type**: <code>Dict[str, Any]</code>

Embed arbitrary plots, tables or text in the report. See the [Custom Content docs](https://docs.seqera.io/multiqc/custom_content) for the full structure.

**Example**:

```yaml
custom_content:
  data:
    my-section-id:
      data:
        sample1:
          col1: 100
        sample2:
          col1: 200
      id: my-section-id
      plot_type: table
      section_name: My Custom Section
  order:
    - my-section-id
    - my-other-section-id
```

#### `custom_content_modules`

**Type**: <code>List[str]</code>

Extra module IDs whose output should be parsed as custom content.

#### `custom_data`

**Type**: <code>Dict[str, Any]</code>

Inline custom content data keyed by section ID. Companion to custom_content for users who prefer splitting the metadata and the data across two top-level keys.

### Module ordering

#### `top_modules`

**Type**: <code>List[Union[str, Dict[str, <a href="#moduleoverride">ModuleOverride</a>]]]</code>

Module IDs to render before module_order. Useful for pinning a module to the top regardless of where it appears in module_order. Same shape as module_order entries.

**Example**:

```yaml
top_modules:
  - fastqc
  - cutadapt
```

#### `module_order`

**Type**: <code>List[Union[str, Dict[str, <a href="#moduleoverride">ModuleOverride</a>]]]</code>

Order in which modules appear in the report. Each entry is either a module ID, or a single-key dict mapping the ID to per-run overrides (eg. name, anchor, info, path_filters, path_filters_exclude, generalstats, custom_config).

<details><summary>Default value</summary>

```yaml
- custom_content
- ccs
- ngsderive
- purple
- conpair
- isoseq
- lima
- peddy
- percolator
- haplocheck
- somalier
- methylqa
- mosdepth
- phantompeakqualtools
- qualimap
- bamdst
- preseq
- hifiasm
- quast
- qorts
- rna_seqc
- rockhopper
- rsem
- rseqc
- busco
- checkm
- bustools
- goleft_indexcov
- gffcompare
- disambiguate
- supernova
- deeptools
- sargasso
- verifybamid
- mirtrace
- happy
- mirtop
- glimpse
- gopeaks
- homer
- hops
- macs2
- theta2
- snpeff
- gatk
- htseq
- bcftools
- featurecounts
- fgbio
- dragen
- dragen_fastqc
- dedup
- pbmarkdup
- damageprofiler
- mapdamage
- biobambam2
- jcvi
- mtnucratio
- picard
- vep
- bakta
- prokka
- checkm2
- qc3C
- nanoq
- nanostat
- samblaster
- samtools
- bamtools
- sambamba
- ngsbits
- pairtools
- sexdeterrmine
- seqera_cli
- eigenstratdatabasetools
- jellyfish
- vcftools
- longranger
- stacks
- varscan2
- snippy
- umicollapse
- umitools
- truvari
- megahit
- sincei
- ganon
- gtdbtk
- bbmap
- bismark
- biscuit
- diamond
- hicexplorer
- hicup
- hicpro
- salmon
- kallisto
- slamdunk
- star
- hisat2
- tophat
- bowtie2
- bowtie1
- hostile
- cellranger
- checkatlas
- snpsplit
- odgi
- vg
- pangolin
- nextclade
- freyja
- humid
- kat
- leehom
- librarian
- nonpareil
- adapterremoval
- bbduk
- clipandmerge
- cutadapt
- trim_galore
- flexbar
- sourmash
- kaiju
- kraken
- malt
- motus
- trimmomatic
- sickle
- skewer
- sortmerna
- ribodetector
- biobloomtools
- seqfu
- fastq_screen
- fastqe
- afterqc
- fastp
- fastqc
- sequali
- filtlong
- prinseqplusplus
- pychopper
- porechop
- pycoqc
- minionqc
- anglerfish
- multivcfanalyzer
- clusterflow
- checkqc
- bcl2fastq
- bclconvert
- interop
- ivar
- flash
- seqyclean
- optitype
- whatshap
- spaceranger
- xenome
- xengsort
- metaphlan
- sylphtax
- seqwho
- telseq
- ataqv
- mgikit
- mosaicatcher
```

</details>

**Example**:

```yaml
module_order:
  - fastqc
  - fastqc:
      name: FastQC (trimmed)
      path_filters:
        - "*_trimmed*"
  - fastqc:
      generalstats: false
      name: FastQC (raw)
  - cutadapt
```

#### `run_modules`

**Type**: <code>List[str]</code>

Module IDs to run. If set, only listed modules are processed (mirror of the --module CLI flag).

**Example**:

```yaml
run_modules:
  - fastqc
  - cutadapt
  - samtools
```

#### `exclude_modules`

**Type**: <code>List[str]</code>

Module IDs to skip (mirror of the --exclude CLI flag).

**Example**:

```yaml
exclude_modules:
  - fastqc
```

#### `remove_sections`

**Type**: <code>List[str]</code>

Module sections to hide. Use the section anchor as it appears in the URL.

**Example**:

```yaml
remove_sections:
  - fastqc_overrepresented_sequences
  - gatk-compare-overlap
```

#### `report_section_order`

**Type**: <code>Dict[str, Union[Literal["remove"], <a href="#sectionorderoverride">SectionOrderOverride</a>]]</code>

Reorder, group or hide report sections by ID. Values are either the literal string 'remove' (drops the section) or a dict with any combination of `order` (int), `before` (str) and `after` (str). See the [customisation docs](https://docs.seqera.io/multiqc/reports/customisation#order-of-module-and-module-subsection-output) for the full grammar.

**Example**:

```yaml
report_section_order:
  custom_content-my-section:
    before: fastqc
  fastqc:
    order: -10
```

### Section comments + indicators

#### `section_comments`

**Type**: <code>Dict[str, str]</code>

Markdown text shown under specific module sections. Keys are section anchors.

**Example**:

```yaml
section_comments:
  fastqc_overrepresented_sequences: "**This is** an important note about the overrepresented\
    \ sequences."
  samtools: Reviewed by *Phil* on 2024-08-21.
```

#### `section_status_checks`

**Type**: <code>Dict[str, Union[bool, Dict[str, bool]]]</code>

Enable or disable the green/yellow/red status indicators on report sections. Top-level keys are module IDs, values are either a bool or a dict mapping section ID to bool.

**Example**:

```yaml
section_status_checks:
  fastqc: true
  samtools:
    alignment_stats: false
```

## Output Options

### Report file

#### `force`

**Type**: <code>bool</code> (default: `false`)

Overwrite existing output files without prompting.

#### `output_fn_name`

**Type**: <code>str</code> (default: `"multiqc_report.html"`)

Filename for the generated HTML report. Defaults to multiqc_report.html.

#### `make_report`

**Type**: <code>bool</code> (default: `true`)

Generate the HTML report. Set to false to only produce data files.

### Data files

#### `make_data_dir`

**Type**: <code>bool</code> (default: `true`)

Write parsed data as files alongside the report.

#### `zip_data_dir`

**Type**: <code>bool</code> (default: `false`)

Compress the data directory into a single .zip file.

#### `data_dir_name`

**Type**: <code>str</code> (default: `"multiqc_data"`)

Name of the directory written alongside the report holding parsed data. Defaults to multiqc_data.

#### `data_format`

**Type**: <code>Literal["tsv", "csv", "json", "yaml"]</code> (default: `"tsv"`)

Format used when writing parsed data files.

#### `data_format_extensions`

**Type**: <code>Dict[str, str]</code> (default: `{"tsv":"txt","csv":"csv","json":"json","yaml":"yaml"}`)

Override the file extension used when writing each data format, eg. {tsv: txt} to write TSV as .txt.

**Example**:

```yaml
data_format_extensions:
  json: json
  tsv: txt
  yaml: yml
```

#### `parquet_format`

**Type**: <code>Literal["long", "wide"]</code> (default: `"long"`)

Parquet table layout. 'long' has rows of (sample_name, metric_name, val_raw, val_raw_type, val_str), easy to filter by metric. 'wide' uses one column per metric (prefixed with table name and namespace), easier for analytics but can hit column limits or mixed-type issues.

### Data dump

#### `data_dump_file`

**Type**: <code>bool</code> (default: `true`)

Write a single JSON file containing all parsed data, for re-running MultiQC later.

#### `data_dump_file_write_raw`

**Type**: <code>bool</code> (default: `true`)

Include raw values (before any normalisation or filtering) in the dumped JSON.

### Plot export

#### `export_plots`

**Type**: <code>bool</code> (default: `false`)

Save each plot as a static image (formats set by export_plot_formats).

#### `export_plot_formats`

**Type**: <code>List[Literal["png", "svg", "pdf"]]</code> (default: `["png","svg","pdf"]`)

Image formats to export when export_plots is on.

#### `export_plots_timeout`

**Type**: <code>int</code> (default: `60`)

Timeout for exporting each plot, in seconds.

#### `plots_dir_name`

**Type**: <code>str</code> (default: `"multiqc_plots"`)

Directory for exported plot images when export_plots is on. Defaults to multiqc_plots.

### PDF

#### `make_pdf`

**Type**: <code>bool</code> (default: `false`)

Also generate a PDF version of the report. Requires Pandoc to be installed.

#### `pandoc_template`

**Type**: <code>str</code>

Path to a Pandoc template used when exporting the report as PDF.

## Sample Names

### Prepend directory

#### `prepend_dirs`

**Type**: <code>bool</code> (default: `false`)

Prefix sample names with their parent directory. Useful when the same sample name occurs in multiple directories.

#### `prepend_dirs_depth`

**Type**: <code>int</code> (default: `0`)

How many parent directories to include. 0 means all the way to the root.

#### `prepend_dirs_sep`

**Type**: <code>str</code> (default: `" | "`)

String inserted between directory names and the sample name. Defaults to '|'.

**Examples**:

```yaml
prepend_dirs_sep: _
```

```yaml
prepend_dirs_sep: " - "
```

### Name cleaning

#### `fn_clean_sample_names`

**Type**: <code>bool</code> (default: `true`)

Apply the cleaning rules in fn_clean_exts and fn_clean_trim to sample names.

#### `extra_fn_clean_exts`

**Type**: <code>List[Union[str, <a href="#cleanpattern">CleanPattern</a>]]</code>

Extensions appended to the built-in list. Use to add custom suffixes without overriding defaults.

**Example**:

```yaml
extra_fn_clean_exts:
  - .mySuffix
  - module:
      - samtools
    pattern: _tmp
    type: remove
```

#### `extra_fn_clean_trim`

**Type**: <code>List[str]</code>

Strings appended to the built-in trim list, without overriding defaults.

**Example**:

```yaml
extra_fn_clean_trim:
  - sample_
  - _processed
```

#### `fn_clean_exts`

**Type**: <code>List[Union[str, <a href="#cleanpattern">CleanPattern</a>]]</code>

Extensions stripped from sample names, eg. .gz, .fastq. Replaces the built-in list.

<details><summary>Default value</summary>

```yaml
- .gz
- .fastq
- .fq
- .bam
- .cram
- .sam
- .sra
- .vcf
- .dat
- _tophat
- .pbmarkdup.log
- .log
- .stderr
- .out
- .spp
- .fa
- .fasta
- .png
- .jpg
- .jpeg
- .html
- Log.final
- ReadsPerGene
- .flagstat
- _star_aligned
- _fastqc
- .hicup
- .counts
- _counts
- .txt
- .tsv
- .csv
- .aligned
- Aligned
- .merge
- .deduplicated
- .dedup
- .clean
- .sorted
- .report
- "| stdin"
- .geneBodyCoverage
- .inner_distance_freq
- .junctionSaturation_plot.r
- .pos.DupRate.xls
- .GC.xls
- _slamdunk
- _bismark
- .conpair
- .concordance
- .contamination
- .BEST.results
- _peaks.xls
- .relatedness
- .cnt
- .aqhist
- .bhist
- .bincov
- .bqhist
- .covhist
- .covstats
- .ehist
- .gchist
- .idhist
- .ihist
- .indelhist
- .lhist
- .mhist
- .qahist
- .qchist
- .qhist
- .rpkm
- .selfSM
- .extendedFrags
- _SummaryStatistics
- .purple.purity
- .purple.qc
- .trim
- .bowtie2
- .mkD
- .highfreq
- .lowfreq
- .consensus
- .snpEff
- .snpeff
- .scaffolds
- .contigs
- .kraken2
- .ccurve
- .hisat2
- _duprate
- .markdup
- .read_distribution
- .junction_annotation
- .infer_experiment
- .biotype
- .ivar
- .mpileup
- .primer_trim
- .mapped
- .vep
- _vep
- ccs
- _NanoStats
- .cutadapt
- .qcML
- .mosdepth
- _gopeaks
- .readCounts
- .wgs_contig_mean_cov
- _overall_mean_cov
- _coverage_metrics
- .wgs_fine_hist
- .wgs_coverage_metrics
- .wgs_hist
- .vc_metrics
- .gvcf_metrics
- .ploidy_estimation_metrics
- _overall_mean_cov
- .fragment_length_hist
- .mapping_metrics
- .gc_metrics
- .trimmer_metrics
- .time_metrics
- .quant_metrics
- .quant.metrics
- .quant.transcript_coverage
- .scRNA_metrics
- .scRNA.metrics
- .scATAC_metrics
- .scATAC.metrics
- .fastqc_metrics
- .labels
- .bammetrics.metrics
- .filter_summary
- .cluster_report
- .error.spl
- .error.grp
- .vgstats
- _mapq_table
- _strand_table
- _isize_table
- _dup_report
- _cv_table
- _covdist_all
- _covdist_q40
- _CpGRetention
- _CpHRetentionByReadPos
- _totalBaseConversionRate
- _totalReadConversionRate
- .sylphmpa
- _qual
- _hifi_trimmer
- .hifi_trimmer
- _trimmer
```

</details>

**Example**:

```yaml
fn_clean_exts:
  - .gz
  - .fastq
  - .bam
  - pattern: _S\d+_L\d+
    type: regex
```

#### `fn_clean_trim`

**Type**: <code>List[str]</code>

Strings trimmed from the start or end of sample names. Replaces the built-in list.

<details><summary>Default value</summary>

```yaml
- .
- ":"
- _
- "-"
- .r
- _val
- .idxstats
- _trimmed
- .trimmed
- .csv
- .yaml
- .yml
- .json
- _mqc
- short_summary_
- _summary
- .summary
- .align
- .h5
- _matrix
- .stats
- .hist
- .phased
- .tar
- runs_
- .qc
```

</details>

**Example**:

```yaml
fn_clean_trim:
  - _R1
  - _R2
  - _001
```

#### `use_filename_as_sample_name`

**Type**: <code>Union[bool, List[str]]</code> (default: `false`)

Use the source filename as the sample name instead of any name parsed from the log. Set to true for all modules, or to a list of module IDs / patterns to apply selectively.

### Ignore samples

#### `sample_names_ignore`

**Type**: <code>List[str]</code>

Glob patterns. Matching samples are dropped from the report.

**Example**:

```yaml
sample_names_ignore:
  - "*_temp"
  - control_*
```

#### `sample_names_ignore_re`

**Type**: <code>List[str]</code>

Regex patterns. Matching samples are dropped from the report.

**Example**:

```yaml
sample_names_ignore_re:
  - ^test_.*
  - .*_neg_ctrl$
```

#### `sample_names_only_include`

**Type**: <code>List[str]</code>

Glob patterns. If set, only matching samples are kept.

**Example**:

```yaml
sample_names_only_include:
  - RNA_*
  - Sample_??
```

#### `sample_names_only_include_re`

**Type**: <code>List[str]</code>

Regex patterns. If set, only matching samples are kept.

**Example**:

```yaml
sample_names_only_include_re:
  - ^WGS_[0-9]+$
```

### Rename and replace

#### `sample_names_rename`

**Type**: <code>List[List[str]]</code>

Toolbox rename rows. Each entry is a list where the first element is the source sample name and each subsequent element is the rename for the corresponding button in `sample_names_rename_buttons` (so inner lists should have `1 + len(sample_names_rename_buttons)` elements).

**Example**:

```yaml
sample_names_rename:
  - - SMP001
    - Patient_A
  - - SMP002
    - Patient_B
  - - SMP003
    - Patient_C
```

#### `sample_names_rename_buttons`

**Type**: <code>List[str]</code>

Names of the toolbox buttons that switch between the rename groups defined in sample_names_rename.

**Example**:

```yaml
sample_names_rename_buttons:
  - Sample ID
  - Patient ID
  - Lane
```

#### `sample_names_replace`

**Type**: <code>Dict[str, str]</code>

Substring replacements applied to every sample name. Keys are matched, values are replacements.

**Example**:

```yaml
sample_names_replace:
  Sample_: S
  _001: ""
```

#### `sample_names_replace_complete`

**Type**: <code>bool</code> (default: `false`)

Replace the entire sample name when the key matches anywhere in it.

#### `sample_names_replace_exact`

**Type**: <code>bool</code> (default: `false`)

Only replace when the key matches the sample name exactly, not as a substring.

#### `sample_names_replace_regex`

**Type**: <code>bool</code> (default: `false`)

Treat keys in sample_names_replace as regex patterns.

## File Discovery

### Input source

#### `file_list`

**Type**: <code>bool</code> (default: `false`)

Treat the input path as a file containing a list of paths to scan, one per line.

#### `require_logs`

**Type**: <code>bool</code> (default: `false`)

Fail with an error if any module explicitly requested with `--module` has no log files found. Off by default, so missing inputs are skipped silently.

### Size limits

#### `log_filesize_limit`

**Type**: <code>int</code> (default: `50000000`)

Skip log files larger than this many bytes.

#### `filesearch_lines_limit`

**Type**: <code>int</code> (default: `1000`)

Stop reading a log file after this many lines.

### Skip patterns

#### `ignore_symlinks`

**Type**: <code>bool</code> (default: `false`)

Skip symlinked files and directories during the file search.

#### `ignore_images`

**Type**: <code>bool</code> (default: `true`)

Skip image files (PNG/JPEG/etc.) to avoid wasting time opening them.

#### `fn_ignore_dirs`

**Type**: <code>List[str]</code> (default: `["multiqc_data",".git","icarus_viewers","runs_per_reference","not_aligned","contigs_reports"]`)

Glob patterns for directory names to skip entirely during the file search.

**Example**:

```yaml
fn_ignore_dirs:
  - work
  - .nextflow
  - "*_logs"
```

#### `fn_ignore_paths`

**Type**: <code>List[str]</code> (default: `["*/work/??/??????????????????????????????","*/.snakemake","*/.singularity","*/__pycache__","*/site-packages/multiqc"]`)

Glob patterns for paths to skip during the file search.

**Example**:

```yaml
fn_ignore_paths:
  - "*/test_data/*"
  - "*/.snakemake/*"
```

#### `fn_ignore_files`

**Type**: <code>List[str]</code>

Glob patterns for file names to skip during the file search.

<details><summary>Default value</summary>

```yaml
- .DS_Store
- .py[cod]
- "*.bam"
- "*.bai"
- "*.sam"
- "*.fq.gz"
- "*.fastq.gz"
- "*.fq"
- "*.fastq"
- "*.fa"
- "*.gtf"
- "*.bed"
- "*.vcf"
- "*.tbi"
- "*.txt.gz"
- "*.pdf"
- "*.md5"
- "*.parquet"
- "*[!s][!u][!m][!_\\.m][!mva][!qer][!cpy].html"
- multiqc_data.json
- "*.gam"
- "*.gamp"
- "*.jar"
```

</details>

**Example**:

```yaml
fn_ignore_files:
  - "*.bai"
  - "*.bak"
  - "*.tmp"
```

#### `filesearch_file_shared`

**Type**: <code>List[str]</code>

Module IDs whose log files may be matched by multiple modules during the search.

### Search patterns

#### `sp`

**Type**: <code>Dict[str, Union[<a href="#searchpattern">SearchPattern</a>, List[<a href="#searchpattern">SearchPattern</a>]]]</code>

Override or add to the built-in module search patterns. Top-level keys are module IDs (eg. `fastqc`); values are a single `SearchPattern` dict or a list of them. See the [SearchPattern](#searchpattern) definition below for the accepted fields.

<details><summary>Default value</summary>

```yaml
multiqc_data:
  fn: "*multiqc.parquet"
adapterremoval:
  fn: "*.settings"
  contents: AdapterRemoval
  num_lines: 1
xenium/metrics:
  fn: metrics_summary.csv
  contents: num_cells_detected
  num_lines: 5
xenium/experiment:
  fn: experiment.xenium
  num_lines: 50
afterqc:
  fn: "*.json"
  contents: allow_mismatch_in_poly
  num_lines: 10000
anglerfish:
  fn: "*.json"
  contents: anglerfish_version
bakta:
  fn: "*.txt"
  contents: "Bakta:"
bamdst/coverage:
  contents: "## The file was created by bamdst"
  num_lines: 5
bamtools/stats:
  contents: "Stats for BAM file(s):"
  num_lines: 10
bases2fastq/run:
  fn: RunStats.json
  contents: SampleStats
  num_lines: 100
bases2fastq/project:
  fn: "*_RunStats.json"
  contents: SampleStats
  num_lines: 100
bases2fastq/manifest:
  fn: RunManifest.json
  contents: Settings
  num_lines: 100
bbduk:
  contents: Executing jgi.BBDuk
  num_lines: 2
bbmap/stats:
  contents:
    - "#File"
    - "#Total"
    - "#Matched"
    - "#Name\tReads\tReadsPct"
  num_lines: 10
bbmap/bbsplit:
  contents: "#name\t%unambiguousReads\tunambiguousMB\t%ambiguousReads"
  num_lines: 5
bbmap/aqhist:
  contents: "#Quality\tcount1\tfraction1\tcount2\tfraction2"
  num_lines: 10
bbmap/bhist:
  contents: "#Pos\tA\tC\tG\tT\tN"
  num_lines: 10
bbmap/bincov:
  contents: "#RefName\tCov\tPos\tRunningPos"
  num_lines: 10
bbmap/bqhist:
  contents: "#BaseNum\tcount_1\tmin_1\tmax_1\tmean_1\tQ1_1\tmed_1\tQ3_1\tLW_1\tRW_1\t\
    count_2\tmin_2\tmax_2\tmean_2\tQ1_2\tmed_2\tQ3_2\tLW_2\tRW_2"
  num_lines: 10
bbmap/covhist:
  contents: "#Coverage\tnumBases"
  num_lines: 10
bbmap/covstats:
  contents: "#ID\tAvg_fold"
  num_lines: 10
bbmap/ehist:
  contents: "#Errors\tCount"
  num_lines: 10
bbmap/gchist:
  contents:
    - "#Mean\t"
    - "#GC\tCount"
  num_lines: 10
bbmap/idhist:
  contents:
    - "#Mean_reads"
    - "#Identity\tReads\tBases"
  num_lines: 10
bbmap/ihist:
  contents:
    - "#Mean\t"
    - "#InsertSize\tCount"
  num_lines: 10
bbmap/indelhist:
  contents: "#Length\tDeletions\tInsertions"
  num_lines: 10
bbmap/lhist:
  contents: "#Length\tCount"
  num_lines: 10
bbmap/mhist:
  contents: "#BaseNum\tMatch1\tSub1\tDel1\tIns1\tN1\tOther1\tMatch2\tSub2\tDel2\t\
    Ins2\tN2\tOther2"
  num_lines: 10
bbmap/qahist:
  contents: "#Quality\tMatch\tSub\tIns\tDel"
  num_lines: 10
bbmap/qchist:
  contents_re: "#Quality\tcount1\tfraction1$"
  num_lines: 10
bbmap/qhist:
  contents: "#BaseNum\tRead1_linear\tRead1_log\tRead1_measured"
  num_lines: 10
bbmap/rpkm:
  contents:
    - "#File\t"
    - "#Reads\t"
    - "#Mapped\t"
    - "#RefSequences\t"
    - "#Name Length"
  num_lines: 10
bbmap/statsfile_machine:
  contents: Reads Used=
  num_lines: 10
bbmap/statsfile:
  contents:
    - "Reads Used:"
    - "Mapping:"
    - "Reads/sec:"
    - "kBases/sec:"
  num_lines: 10
bcftools/stats:
  contents: This file was produced by bcftools stats
bcl2fastq:
  fn: Stats.json
  contents: DemuxResults
  num_lines: 300
bclconvert/runinfo:
  fn: RunInfo.xml
bclconvert/demux:
  fn: Demultiplex_Stats.csv
bclconvert/quality_metrics:
  fn: Quality_Metrics.csv
bclconvert/adaptermetrics:
  fn: Adapter_Metrics.csv
bclconvert/unknown_barcodes:
  fn: Top_Unknown_Barcodes.csv
biobambam2/bamsormadup:
  contents: "# bamsormadup"
  num_lines: 2
biobloomtools:
  contents: "filter_id\thits\tmisses\tshared\trate_hit\trate_miss\trate_shared"
  num_lines: 2
biscuit/align_mapq:
  fn: "*_mapq_table.txt"
  contents: BISCUITqc Mapping Quality Table
  num_lines: 3
biscuit/align_strand:
  fn: "*_strand_table.txt"
  contents: BISCUITqc Strand Table
  num_lines: 3
biscuit/align_isize:
  fn: "*_isize_table.txt"
  contents: BISCUITqc Insert Size Table
  num_lines: 3
biscuit/dup_report:
  fn: "*_dup_report.txt"
  contents: BISCUITqc Read Duplication Table
  num_lines: 3
biscuit/qc_cv:
  fn: "*_cv_table.txt"
  contents: BISCUITqc Uniformity Table
  num_lines: 3
biscuit/covdist_all_base_botgc:
  fn: "*_covdist_all_base_botgc_table.txt"
biscuit/covdist_all_base:
  fn: "*_covdist_all_base_table.txt"
biscuit/covdist_all_base_topgc:
  fn: "*_covdist_all_base_topgc_table.txt"
biscuit/covdist_q40_base_botgc:
  fn: "*_covdist_q40_base_botgc_table.txt"
biscuit/covdist_q40_base:
  fn: "*_covdist_q40_base_table.txt"
biscuit/covdist_q40_base_topgc:
  fn: "*_covdist_q40_base_topgc_table.txt"
biscuit/covdist_all_cpg_botgc:
  fn: "*_covdist_all_cpg_botgc_table.txt"
biscuit/covdist_all_cpg:
  fn: "*_covdist_all_cpg_table.txt"
biscuit/covdist_all_cpg_topgc:
  fn: "*_covdist_all_cpg_topgc_table.txt"
biscuit/covdist_q40_cpg_botgc:
  fn: "*_covdist_q40_cpg_botgc_table.txt"
biscuit/covdist_q40_cpg:
  fn: "*_covdist_q40_cpg_table.txt"
biscuit/covdist_q40_cpg_topgc:
  fn: "*_covdist_q40_cpg_topgc_table.txt"
biscuit/cpg_retention_readpos:
  fn: "*_CpGRetentionByReadPos.txt"
biscuit/cph_retention_readpos:
  fn: "*_CpHRetentionByReadPos.txt"
biscuit/base_avg_retention_rate:
  fn: "*_totalBaseConversionRate.txt"
biscuit/read_avg_retention_rate:
  fn: "*_totalReadConversionRate.txt"
bismark/align:
  fn: "*_[SP]E_report.txt"
bismark/dedup:
  fn: "*.deduplication_report.txt"
bismark/meth_extract:
  fn: "*_splitting_report.txt"
bismark/m_bias:
  fn: "*M-bias.txt"
bismark/bam2nuc:
  fn: "*.nucleotide_stats.txt"
bowtie1:
  contents: "# reads processed:"
  exclude_fn:
    - bowtie.left_kept_reads.log
    - bowtie.left_kept_reads.m2g_um.log
    - bowtie.left_kept_reads.m2g_um_seg1.log
    - bowtie.left_kept_reads.m2g_um_seg2.log
    - bowtie.right_kept_reads.log
    - bowtie.right_kept_reads.m2g_um.log
    - bowtie.right_kept_reads.m2g_um_seg1.log
    - bowtie.right_kept_reads.m2g_um_seg2.log
  shared: true
bowtie2:
  contents: "reads; of these:"
  exclude_contents:
    - bisulfite
    - HiC-Pro
  shared: true
busco:
  fn: short_summary*
  contents: "BUSCO version is:"
  num_lines: 1
bustools:
  fn: "*inspect.json"
ccs/v4:
  contents: ZMWs generating CCS
  num_lines: 2
  max_filesize: 1024
ccs/v5:
  contents: '"id": "ccs_processing"'
  fn: "*.json"
checkatlas/summary:
  fn: "*.tsv"
  contents_re: ^AtlasFileType\tNbCells\tNbGenes
  num_lines: 1
checkatlas/adata:
  fn: "*.tsv"
  contents_re: ^atlas_obs\tobsm\tvar\tvarm\tuns
  num_lines: 1
checkatlas/qc:
  fn: "*.tsv"
  contents_re: cellrank_(total_counts|n_genes_by_counts|pct_counts_mt)
  num_lines: 1
checkatlas/cluster:
  fn: "*.tsv"
  contents_re: ^Clust_Sample\tobs
  num_lines: 1
checkatlas/annotation:
  fn: "*.tsv"
  contents_re: ^Annot_Sample\tReference\tobs
  num_lines: 1
checkatlas/dimred:
  fn: "*.tsv"
  contents_re: ^Dimred_Sample\tobsm
  num_lines: 1
cellranger/count_html:
  - fn: "*.html"
    contents: '"command":"Cell Ranger","subcommand":"count"'
    num_lines: 20
  - fn: "*.html"
    contents: '"command": "Cell Ranger", "subcommand": "count"'
    num_lines: 20
cellranger/vdj_html:
  - fn: "*.html"
    contents: '"command":"Cell Ranger","subcommand":"vdj"'
    num_lines: 20
  - fn: "*.html"
    contents: '"command": "Cell Ranger", "subcommand": "vdj"'
    num_lines: 20
cellranger_arc:
  - fn: "*.html"
    contents: Cell Ranger ARC
    num_lines: 250
cells2stats/run:
  fn: RunStats.json
  contents: '"AnalysisID": "c2s.'
  num_lines: 100
checkm:
  - contents_re: ".*Bin Id(?:\t| {3,})Marker lineage(?:\t| {3,})# genomes(?:\t| {3,})#\
      \ markers(?:\t| {3,})# marker sets.*"
    num_lines: 10
checkm2:
  contents: "Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used"
  num_lines: 10
checkqc:
  contents: instrument_and_reagent_type
  fn: "*.json"
custom_content:
  fn_re: .+_mqc\.(yaml|yml|json|txt|csv|tsv|log|out|png|jpg|jpeg|gif|webp|tiff|html|md)
clipandmerge:
  contents: ClipAndMerge (
  num_lines: 5
clusterflow/logs:
  fn: "*_clusterFlow.txt"
  shared: true
clusterflow/runfiles:
  fn: "*.run"
  contents: Cluster Flow Run File
  num_lines: 2
conpair/concordance:
  contents: markers (coverage per marker threshold
  num_lines: 3
conpair/contamination:
  contents: "Tumor sample contamination level: "
  num_lines: 3
cutadapt:
  - contents: This is cutadapt
    exclude_contents_re: "Trim Galore version: (?:[2-9]|\\d{2,})\\."
    num_lines: 100
  - fn: "*.json"
    contents: Cutadapt report
damageprofiler:
  fn: "*dmgprof.json"
deacon:
  fn: "*.json"
  contents: '"version": "deacon'
  num_lines: 30
dedup:
  fn: "*.json"
  contents: '"tool_name": "DeDup"'
  num_lines: 20
deeptools/bamPEFragmentSizeTable:
  contents: "\tFrag. Sampled\tFrag. Len. Min.\tFrag. Len. 1st. Qu.\tFrag. Len. Mean\t\
    Frag. Len. Median\tFrag. Len. 3rd Qu."
  num_lines: 1
deeptools/bamPEFragmentSizeDistribution:
  contents: "#bamPEFragmentSize"
  num_lines: 1
deeptools/estimateReadFiltering:
  contents: "Sample\tTotal Reads\tMapped Reads\tAlignments in blacklisted regions\t\
    Estimated mapped reads"
  num_lines: 1
deeptools/plotCorrelationData:
  contents: "#plotCorrelation --outFileCorMatrix"
  num_lines: 1
deeptools/plotCoverageStdout:
  contents: "sample\tmean\tstd\tmin\t25%\t50%\t75%\tmax"
  num_lines: 1
deeptools/plotCoverageOutRawCounts:
  contents: "#plotCoverage --outRawCounts"
  num_lines: 1
deeptools/plotEnrichment:
  contents: "file\tfeatureType\tpercent\tfeatureReadCount\ttotalReadCount"
  num_lines: 1
deeptools/plotFingerprintOutRawCounts:
  contents: "#plotFingerprint --outRawCounts"
  num_lines: 1
deeptools/plotFingerprintOutQualityMetrics:
  contents: "Sample\tAUC\tSynthetic AUC\tX-intercept\tSynthetic X-intercept\tElbow\
    \ Point\tSynthetic Elbow Point"
  num_lines: 1
deeptools/plotPCAData:
  contents: "#plotPCA --outFileNameData"
  num_lines: 1
deeptools/plotProfile:
  contents: bin labels
  num_lines: 1
diamond:
  fn: diamond.log
disambiguate:
  contents: unique species A pairs
  num_lines: 2
dragen/vc_metrics:
  fn: "*.vc_metrics.csv"
dragen/gvcf_metrics:
  fn: "*.gvcf_metrics.csv"
dragen/ploidy_estimation_metrics:
  fn: "*.ploidy_estimation_metrics.csv"
dragen/wgs_contig_mean_cov:
  fn_re: .*\.wgs_contig_mean_cov_?(tumor|normal)?\.csv
dragen/overall_mean_cov_metrics:
  fn_re: .*_overall_mean_cov.*\.csv
dragen/coverage_metrics:
  fn_re: .*_coverage_metrics.*\.csv
dragen/wgs_fine_hist:
  fn_re: .*\.wgs_fine_hist_?(tumor|normal)?\.csv
dragen/fragment_length_hist:
  fn: "*.fragment_length_hist.csv"
dragen/mapping_metrics:
  fn: "*.mapping_metrics.csv"
  contents: Number of unique reads (excl. duplicate marked reads)
  num_lines: 50
dragen/gc_metrics:
  fn: "*.gc_metrics.csv"
dragen/trimmer_metrics:
  fn: "*.trimmer_metrics.csv"
dragen/time_metrics:
  fn: "*.time_metrics.csv"
dragen/rna_quant_metrics:
  fn: "*.quant[._]metrics.csv"
dragen/rna_transcript_cov:
  fn: "*.quant.transcript_coverage.txt"
dragen/sc_rna_metrics:
  fn: "*.scRNA[._]metrics.csv"
dragen/sc_atac_metrics:
  fn: "*.scATAC[._]metrics.csv"
dragen_fastqc:
  fn: "*.fastqc_metrics.csv"
eigenstratdatabasetools:
  fn: "*_eigenstrat_coverage.json"
fastp:
  fn: "*.json"
  contents: '"before_filtering": {'
  num_lines: 50
fastq_screen:
  fn: "*_screen.txt"
fastqe:
  fn: "*fastqe*"
  contents: "Filename\tStatistic\tQualities"
  num_lines: 1
fastqc/data:
  fn: "*fastqc_data.txt"
fastqc/zip:
  fn: "*_fastqc.zip"
fastqc/theoretical_gc:
  fn: "*fastqc_theoretical_gc*"
featurecounts:
  fn: "*.summary"
  shared: true
fgbio/groupreadsbyumi:
  contents: fraction_gt_or_eq_family_size
  num_lines: 3
fgbio/errorratebyreadposition:
  contents: "read_number\tposition\tbases_total\terrors\terror_rate\ta_to_c_error_rate\t\
    a_to_g_error_rate\ta_to_t_error_rate\tc_to_a_error_rate\tc_to_g_error_rate\tc_to_t_error_rate"
  num_lines: 3
filtlong:
  contents: Scoring long reads
  contents_re: .*Filtering long reads.*
  num_lines: 5
flash/log:
  contents: "[FLASH]"
flash/hist:
  fn: "*flash*.hist"
flexbar:
  contents: Flexbar - flexible barcode and adapter removal
freyja:
  fn: "*.tsv"
  contents: "summarized\t["
  num_lines: 6
ganon:
  contents:
    - ganon-classify processed
  num_lines: 100
gatk/varianteval:
  contents: "#:GATKTable:TiTvVariantEvaluator"
gatk/base_recalibrator:
  - contents: "#:GATKTable:Arguments:Recalibration"
    num_lines: 3
  - contents: "#:SENTIEON_QCAL_TABLE:Arguments:Recalibration"
    num_lines: 3
gatk/analyze_saturation_mutagenesis:
  fn: "*.readCounts"
  contents: ">>Reads in disjoint pairs evaluated separately:"
  num_lines: 10
gffcompare:
  fn: "*.stats"
  contents: "# gffcompare"
  num_lines: 2
glimpse/err_spl:
  fn: "*.error.spl.txt.gz"
  num_lines: 1
glimpse/err_grp:
  fn: "*.error.grp.txt.gz"
  num_lines: 1
goleft_indexcov/roc:
  fn: "*-indexcov.roc"
goleft_indexcov/ped:
  fn: "*-indexcov.ped"
gopeaks:
  fn: "*_gopeaks.json"
gtdbtk:
  contents: "user_genome\tclassification\tclosest_genome_reference\tclosest_genome_reference_radius\t\
    closest_genome_taxonomy\tclosest_genome_ani"
  num_lines: 10
haplocheck:
  contents: "\"Sample\"\t\"Contamination Status\"\t\"Contamination Level\"\t\"Distance\"\
    \t\"Sample Coverage\""
  num_lines: 10
happy:
  fn: "*.summary.csv"
  contents: Type,Filter,TRUTH
htseq:
  - contents_re: ^feature\tcount$
    num_lines: 1
    shared: true
  - contents_re: ^\w+.*\t\d+$
    num_lines: 1
    shared: true
hicexplorer:
  contents: Min rest. site distance
  max_filesize: 4096
  num_lines: 26
hicup:
  fn: HiCUP_summary_report*
hicup/html:
  fn: "*HiCUP_summary_report*.html"
hicpro/mmapstat:
  fn: "*mapstat"
  contents: total_R
  num_lines: 10
hicpro/mpairstat:
  fn: "*pairstat"
  contents: Total_pairs_processed
  num_lines: 10
hicpro/mergestat:
  fn: "*.mergestat"
  contents: valid_interaction
  num_lines: 10
hicpro/mRSstat:
  fn: "*RSstat"
  contents: Valid_interaction_pairs
hicpro/assplit:
  fn: "*assplit.stat"
hicstuff/pipeline_stats:
  - fn: "*.txt"
    contents: "## hicstuff:"
    num_lines: 100
  - fn: "*.log"
    contents: "## hicstuff:"
    num_lines: 10
hicstuff/distancelaw:
  contents: "## distance_law"
  num_lines: 5
hifiasm:
  contents: "[M::ha_analyze_count]"
  num_lines: 1
hifi_trimmer:
  fn: "*.json"
  contents: '"total_reads_trimmed"'
  num_lines: 10
hisat2:
  contents: "HISAT2 summary stats:"
homer/findpeaks:
  contents: "# HOMER Peaks"
  num_lines: 3
homer/GCcontent:
  fn: tagGCcontent.txt
homer/genomeGCcontent:
  fn: genomeGCcontent.txt
homer/RestrictionDistribution:
  fn: petagRestrictionDistribution.*.txt
homer/LengthDistribution:
  fn: tagLengthDistribution.txt
homer/tagInfo:
  fn: tagInfo.txt
homer/FreqDistribution:
  fn: petag.FreqDistribution_1000.txt
hops:
  fn: heatmap_overview_Wevid.json
hostile:
  fn: "*.json"
  contents: '"reads_removed_proportion"'
  num_lines: 100
humid/stats:
  fn: stats.dat
  contents: "total: "
  num_lines: 1
humid/neighbours:
  fn: neigh.dat
  contents_re: "[0-9]+ [0-9]+"
  num_lines: 1
humid/counts:
  fn: counts.dat
  contents_re: "[0-9]+ [0-9]+"
  num_lines: 1
humid/clusters:
  fn: clusters.dat
  contents_re: "[0-9]+ [0-9]+"
  num_lines: 1
interop/summary:
  contents: Level,Yield,Projected Yield,Aligned,Error Rate,Intensity C1,%>=Q30
interop/index-summary:
  contents: Total Reads,PF Reads,% Read Identified (PF),CV,Min,Max
isoseq/refine-json:
  contents: '"num_reads_fl"'
  fn: "*.json"
isoseq/refine-csv:
  contents: id,strand,fivelen,threelen,polyAlen,insertlen,primer
  fn: "*.csv"
isoseq/cluster-csv:
  contents: cluster_id
  fn: "*cluster_report.csv"
  num_lines: 1
ivar/trim:
  contents: Number of references
  num_lines: 8
jcvi:
  contents: "     o    % GC    % of genome    Average size (bp)    Median size (bp)\
    \    Number    Total length (Mb)"
jellyfish:
  fn: "*_jf.hist"
kaiju:
  contents_re: file\tpercent\treads\ttaxon_id\ttaxon_name
  num_lines: 1
kallisto:
  contents: "[quant] finding pseudoalignments for the reads"
kat:
  fn: "*.dist_analysis.json"
kraken:
  contents_re: ^\s{0,2}(\d{1,3}\.\d{1,2})\t(\d+)\t(\d+)\t((\d+)\t(\d+)\t)?([URDKPCOFGS-]\d{0,2})\t(\d+)(\s+)[root|unclassified]
  num_lines: 2
librarian:
  fn: librarian_heatmap.txt
leehom:
  contents: Adapter dimers/chimeras
  num_lines: 100
lima/summary:
  contents: ZMWs above all thresholds
  num_lines: 2
  max_filesize: 1024
lima/counts:
  contents: "IdxFirst\tIdxCombined\tIdxFirstNamed\tIdxCombinedNamed\tCounts\tMeanScore"
  num_lines: 1
longranger/summary:
  fn: "*summary.csv"
  contents: longranger_version,instrument_ids,gems_detected,mean_dna_per_gem,bc_on_whitelist,bc_mean_qscore,n50_linked_reads_per_molecule
  num_lines: 2
longranger/invocation:
  fn: _invocation
  contents: call PHASER_SVCALLER_CS(
  max_filesize: 2048
macs2:
  fn: "*_peaks.xls"
malt:
  contents: MaltRun - Aligns sequences using MALT (MEGAN alignment tool)
  num_lines: 2
mapdamage:
  - fn: 3p*_freq.txt
  - fn: 5p*_freq.txt
  - fn: lgdistribution.txt
megahit:
  contents: " - MEGAHIT v"
  num_lines: 5
metaphlan:
  fn: "*.txt"
  contents: "#clade_name\tNCBI_tax_id\trelative_abundance\t"
methurator:
  fn: "*methurator_summary.yml"
methylqa:
  fn: "*.report"
  shared: true
mgikit/mgi_ambiguous_barcode:
  fn: "*.mgikit.ambiguous_barcode"
mgikit/mgi_sample_stats:
  fn: "*.mgikit.sample_stats"
mgikit/mgi_general_info:
  fn: "*.mgikit.general"
mgikit/mgi_sample_reads:
  fn: "*.mgikit.info"
mgikit/mgi_undetermined_barcode:
  fn: "*.mgikit.undetermined_barcode"
minionqc:
  fn: summary.yaml
  contents: total.gigabases
mirtop:
  fn: "*_mirtop_stats.log"
mirtrace/summary:
  fn: mirtrace-results.json
mirtrace/length:
  fn: mirtrace-stats-length.tsv
mirtrace/contaminationbasic:
  fn: mirtrace-stats-contamination_basic.tsv
mirtrace/mirnacomplexity:
  fn: mirtrace-stats-mirna-complexity.tsv
mtnucratio:
  fn: "*mtnuc.json"
mosdepth/summary:
  fn: "*.mosdepth.summary.txt"
mosdepth/global_dist:
  fn: "*.mosdepth.global.dist.txt"
mosdepth/region_dist:
  fn: "*.mosdepth.region.dist.txt"
motus:
  contents: Reads are aligned (by BWA) to marker gene sequences in the reference database
  num_lines: 2
multivcfanalyzer:
  fn: MultiVCFAnalyzer.json
nanostat:
  max_filesize: 4096
  contents_re: Metrics\s+dataset\s*
  num_lines: 1
nanostat/legacy:
  max_filesize: 4096
  contents_re: General summary:\s*
  num_lines: 1
nanoq:
  contents: Nanoq Read Summary
  num_lines: 3
nextclade:
  contents: seqName;clade;
  num_lines: 1
ngsbits/readqc:
  - fn: "*.qcML"
    contents: ReadQC
    num_lines: 20
  - fn: "*.qcML"
    contents: SeqPurge
    num_lines: 20
ngsbits/mappingqc:
  - fn: "*.qcML"
    contents: MappingQC
    num_lines: 20
ngsbits/samplegender:
  - fn: "*_ngsbits_sex.tsv"
ngsderive/strandedness:
  contents: "File\tTotalReads\tForwardPct\tReversePct\tPredicted"
  num_lines: 1
ngsderive/instrument:
  contents: "File\tInstrument\tConfidence\tBasis"
  num_lines: 1
ngsderive/readlen:
  contents: "File\tEvidence\tMajorityPctDetected\tConsensusReadLength"
  num_lines: 1
ngsderive/encoding:
  contents: "File\tEvidence\tProbableEncoding"
  num_lines: 1
ngsderive/junction_annotation:
  contents: "File\ttotal_junctions\ttotal_splice_events\tknown_junctions\tpartial_novel_junctions\t\
    complete_novel_junctions\tknown_spliced_reads\tpartial_novel_spliced_reads\tcomplete_novel_spliced_reads"
  num_lines: 1
nonpareil:
  - fn: "*.json"
    contents: LRstar
    num_lines: 50
    max_filesize: 1048576
optitype:
  contents: "\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective"
  num_lines: 1
pangolin:
  contents: pangolin_version
  num_lines: 1
odgi:
  - fn: "*.og.stats.yaml"
  - fn: "*.og.stats.yml"
  - fn: "*.odgi.stats.yaml"
  - fn: "*.odgi.stats.yml"
pairtools:
  contents:
    - "total_single_sided_mapped\t"
    - "cis\t"
    - "trans\t"
    - pair_types/
  num_lines: 20
peddy/summary_table:
  fn: "*.peddy.ped"
peddy/het_check:
  fn: "*.het_check.csv"
peddy/ped_check:
  fn: "*.ped_check.csv"
peddy/sex_check:
  fn: "*.sex_check.csv"
peddy/background_pca:
  fn: "*.background_pca.json"
percolator:
  fn: "*percolator_feature_weights.tsv"
seqera_cli/run_dump:
  fn: runs_*.tar.gz
seqera_cli/json:
  fn: workflow.json
sequali:
  fn: "*.json"
  contents: '"sequali_version"'
  num_lines: 10
somalier/somalier-ancestry:
  fn: "*.somalier-ancestry.tsv"
somalier/samples:
  fn: "*.samples.tsv"
  contents: "#family_id"
  num_lines: 5
somalier/pairs:
  fn: "*.pairs.tsv"
  contents: hom_concordance
  num_lines: 5
sourmash/compare:
  fn: "*.labels.txt"
sourmash/gather:
  contents: intersect_bp,f_orig_query,f_match,f_unique_to_query,f_unique_weighted,
  num_lines: 1
pbmarkdup:
  contents_re: LIBRARY +READS +UNIQUE MOLECULES +DUPLICATE READS
  num_lines: 5
phantompeakqualtools/out:
  fn: "*.spp.out"
picard/alignment_metrics:
  - contents: picard.analysis.AlignmentSummaryMetrics
  - contents: --algo AlignmentStat
picard/basedistributionbycycle:
  contents: BaseDistributionByCycleMetrics
picard/crosscheckfingerprints:
  contents: CrosscheckFingerprints
picard/gcbias:
  - contents: GcBiasDetailMetrics
  - contents: GcBiasSummaryMetrics
  - contents: --algo GCBias
picard/hsmetrics:
  - contents: HsMetrics
  - contents: --algo HsMetricAlgo
picard/insertsize:
  - contents: picard.analysis.InsertSizeMetrics
  - contents: --algo InsertSizeMetricAlgo
picard/markdups:
  - contents: picard.sam.MarkDuplicates
  - contents: picard.sam.DuplicationMetrics
  - contents: picard.sam.markduplicates.MarkDuplicates
  - contents: markduplicates.DuplicationMetrics
  - contents: MarkDuplicatesSpark
  - contents: markduplicates.GATKDuplicationMetrics
  - contents: --algo Dedup
picard/oxogmetrics:
  - contents: "# picard.analysis.CollectOxoGMetrics"
  - contents: "# CollectOxoGMetrics"
  - contents_re: "# CollectMultipleMetrics .*OxoGMetrics"
    shared: true
picard/pcr_metrics:
  - contents: "# picard.analysis.directed.CollectTargetedPcrMetrics"
  - contents_re: "# CollectMultipleMetrics .*TargetedPcrMetrics"
    shared: true
picard/quality_by_cycle:
  - contents: "# MeanQualityByCycle"
  - contents: --algo MeanQualityByCycle
  - contents_re: .*CollectMultipleMetrics.*MeanQualityByCycle
    shared: true
picard/quality_score_distribution:
  - contents: "# QualityScoreDistribution"
  - contents: --algo QualDistribution
  - contents_re: .*CollectMultipleMetrics.*QualityScoreDistribution
    shared: true
picard/quality_yield_metrics:
  - contents: "# CollectQualityYieldMetrics"
  - contents_re: .*CollectMultipleMetrics.*QualityYieldMetrics
    shared: true
picard/rnaseqmetrics:
  - contents: "# picard.analysis.Collectrnaseqmetrics"
  - contents: "# picard.analysis.CollectRnaSeqMetrics"
  - contents: "# CollectRnaSeqMetrics"
  - contents_re: "# CollectMultipleMetrics .*RnaSeqMetrics"
    shared: true
picard/rrbs_metrics:
  - contents: "# picard.analysis.CollectRrbsMetrics"
  - contents_re: "# CollectMultipleMetrics .*RrbsMetrics"
    shared: true
picard/sam_file_validation:
  fn: "*[Vv]alidate[Ss]am[Ff]ile*"
picard/variant_calling_metrics:
  contents_re: "## METRICS CLASS.*VariantCallingDetailMetrics"
picard/wgs_metrics:
  - contents: --algo WgsMetricsAlgo
  - contents_re: "## METRICS CLASS.*WgsMetrics"
    shared: true
picard/collectilluminabasecallingmetrics:
  contents: CollectIlluminaBasecallingMetrics
picard/collectilluminalanemetrics:
  contents: CollectIlluminaLaneMetrics
picard/extractilluminabarcodes:
  contents: ExtractIlluminaBarcodes
picard/markilluminaadapters:
  contents: MarkIlluminaAdapters
porechop:
  contents: Looking for known adapter sets
  num_lines: 10
preseq:
  - contents: EXPECTED_DISTINCT
    num_lines: 2
  - contents: distinct_reads
    num_lines: 2
preseq/real_counts:
  fn: "*preseq_real_counts*"
prinseqplusplus:
  - contents: reads removed by -
    num_lines: 2
prokka:
  contents: "contigs:"
  num_lines: 2
purple/qc:
  fn: "*.purple.qc"
purple/purity:
  fn: "*.purple.purity.tsv"
pycoqc:
  contents: '"pycoqc":'
  num_lines: 2
pychopper:
  contents: "Classification\tRescue"
  num_lines: 6
qc3C:
  fn: "*.qc3C.json"
qorts:
  contents: BENCHMARK_MinutesOnSamIteration
  num_lines: 100
qorts/log:
  fn: QC.*.log
  contents: Starting QoRTs
  num_lines: 2
qualimap/bamqc/genome_results:
  fn: genome_results.txt
qualimap/bamqc/coverage:
  fn: coverage_histogram.txt
qualimap/bamqc/insert_size:
  fn: insert_size_histogram.txt
qualimap/bamqc/genome_fraction:
  fn: genome_fraction_coverage.txt
qualimap/bamqc/gc_dist:
  fn: mapped_reads_gc-content_distribution.txt
qualimap/bamqc/html:
  fn: qualimapReport.html
  contents: "Qualimap report: BAM QC"
  num_lines: 10
qualimap/rnaseq/rnaseq_results:
  fn: rnaseq_qc_results.txt
qualimap/rnaseq/coverage:
  fn: coverage_profile_along_genes_(total).txt
qualimap/rnaseq/html:
  fn: qualimapReport.html
  contents: "Qualimap report: RNA Seq QC"
  num_lines: 10
quast:
  fn: report.tsv
  contents: "Assembly\t"
  num_lines: 2
rna_seqc/metrics_v1:
  fn: "*metrics.tsv"
  contents: "Sample\tNote\t"
rna_seqc/metrics_v2:
  fn: "*metrics.tsv"
  contents: High Quality Ambiguous Alignment Rate
rna_seqc/coverage:
  fn_re: meanCoverageNorm_(high|medium|low)\.txt
rna_seqc/correlation:
  fn_re: corrMatrix(Pearson|Spearman)\.txt
rna_seqc/html:
  fn: index.html
  contents: RNA-SeQC</a> v
  num_lines: 200
ribotish/qual:
  fn: "*_qual.txt"
  num_lines: 10
ribowaltz/psite_region:
  fn: "*ribowaltz*psite_region.tsv"
  contents_re: "sample[,\t]region[,\t]count[,\t]scaled_count"
  num_lines: 1
ribowaltz/frames:
  fn: "*ribowaltz*frames.tsv"
  contents_re: "sample[,\t]region[,\t]frame[,\t]count[,\t]scaled_count"
  num_lines: 1
ribowaltz/metaprofile:
  fn: "*ribowaltz*metaprofile_psite.tsv"
  contents_re: "sample[,\t]region[,\t]x[,\t]y"
  num_lines: 1
riker/alignment:
  fn: "*.alignment-metrics.txt"
  contents_re: ^sample\b.*\bcategory\b
  num_lines: 1
riker/basic_base_dist:
  fn: "*.base-distribution-by-cycle.txt"
  contents_re: ^sample\b.*\bfrac_a\b
  num_lines: 1
riker/basic_mean_quality:
  fn: "*.mean-quality-by-cycle.txt"
  contents_re: ^sample\b.*\bmean_quality\b
  num_lines: 1
riker/basic_quality_dist:
  fn: "*.quality-score-distribution.txt"
  contents_re: ^sample\b.*\bfrac_bases\b
  num_lines: 1
riker/gcbias_detail:
  fn: "*.gcbias-detail.txt"
  contents_re: ^sample\b.*\bnormalized_coverage\b
  num_lines: 1
riker/gcbias_summary:
  fn: "*.gcbias-summary.txt"
  contents_re: ^sample\b.*\bgc_0_19_normcov\b
  num_lines: 1
riker/hybcap_metrics:
  fn: "*.hybcap-metrics.txt"
  contents_re: ^sample\b.*\bbait_territory\b
  num_lines: 1
riker/isize_metrics:
  fn: "*.isize-metrics.txt"
  contents_re: ^sample\b.*\bpair_orientation\b
  num_lines: 1
riker/isize_histogram:
  fn: "*.isize-histogram.txt"
  contents_re: ^sample\b.*\bfr_count\b
  num_lines: 1
riker/wgs_metrics:
  fn: "*.wgs-metrics.txt"
  contents_re: ^sample\b.*\bgenome_territory\b
  num_lines: 1
riker/wgs_coverage:
  fn: "*.wgs-coverage.txt"
  contents_re: ^sample\b.*\bbases_at_or_above\b
  num_lines: 1
rockhopper:
  fn: summary.txt
  contents: Number of gene-pairs predicted to be part of the same operon
  max_filesize: 500000
ribodetector:
  contents: Writing output non-rRNA sequences into file
  num_lines: 20
rsem:
  fn: "*.cnt"
rseqc/bam_stat:
  contents: "Proper-paired reads map to different chrom:"
  max_filesize: 500000
rseqc/gene_body_coverage:
  fn: "*.geneBodyCoverage.txt"
rseqc/inner_distance:
  fn: "*.inner_distance_freq.txt"
rseqc/junction_annotation:
  contents: "Partial Novel Splicing Junctions:"
  max_filesize: 500000
rseqc/junction_saturation:
  fn: "*.junctionSaturation_plot.r"
rseqc/read_gc:
  fn: "*.GC.xls"
rseqc/read_distribution:
  contents: Group               Total_bases         Tag_count           Tags/Kb
  max_filesize: 500000
rseqc/read_duplication_pos:
  fn: "*.pos.DupRate.xls"
rseqc/infer_experiment:
  - fn: "*infer_experiment.txt"
  - contents: Fraction of reads explained by
    max_filesize: 500000
rseqc/tin:
  fn: "*.summary.txt"
  contents: TIN(median)
  num_lines: 1
salmon/meta:
  fn: meta_info.json
  contents: salmon_version
  num_lines: 10
  max_filesize: 50000
salmon/lfc:
  fn: lib_format_counts.json
salmon/fld:
  fn: flenDist.txt
sambamba/markdup:
  contents: finding positions of the duplicate reads in the file
  num_lines: 50
samblaster:
  contents: "samblaster: Version"
samtools/stats:
  contents: This file was produced by samtools stats
samtools/flagstat:
  contents: in total (QC-passed reads + QC-failed reads)
samtools/idxstats:
  fn: "*idxstat*"
samtools/rmdup:
  contents: "[bam_rmdup"
samtools/ampliconclip:
  contents:
    - "COMMAND:"
    - samtools ampliconclip
  num_lines: 11
samtools/coverage:
  contents: "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\t\
    meanmapq"
  num_lines: 10
samtools/markdup_txt:
  contents:
    - "^COMMAND:"
    - samtools markdup
  num_lines: 2
samtools/markdup_json:
  contents:
    - '"COMMAND":'
    - samtools markdup
  num_lines: 10
sargasso:
  fn: overall_filtering_summary.txt
seqfu/stats:
  contents: "File\t#Seq\tTotal bp\tAvg\tN50\tN75\tN90\tauN\tMin\tMax"
  num_lines: 1
seqkit/stats:
  contents_re: ^file\s+format\s+type\s+num_seqs\s+sum_len
  num_lines: 1
seqwho:
  contents: '  "Per Base Seq": ['
  num_lines: 10
seqyclean:
  fn: "*_SummaryStatistics.tsv"
sexdeterrmine:
  fn: sexdeterrmine.json
sickle:
  contents_re: "FastQ \\w*\\s?records kept: .*"
  num_lines: 2
sincei/scFilterStats:
  contents: "Cell_ID\tTotal_sampled\tFiltered\tBlacklisted\tLow_MAPQ\tMissing_Flags\t\
    Excluded_Flags"
  num_lines: 1
sincei/scCountQC:
  fn: "*.cells.tsv"
  contents: "Cell_ID\tbarcodes\tsample\tn_genes_by_counts\tlog1p_n_genes_by_counts\t\
    total_counts"
skewer:
  contents: "maximum error ratio allowed (-r):"
slamdunk/summary:
  contents: "# slamdunk summary"
  num_lines: 1
slamdunk/PCA:
  contents: "# slamdunk PCA"
  num_lines: 1
slamdunk/rates:
  contents: "# slamdunk rates"
  num_lines: 1
slamdunk/utrrates:
  contents: "# slamdunk utrrates"
  num_lines: 1
slamdunk/tcperreadpos:
  contents: "# slamdunk tcperreadpos"
  num_lines: 1
slamdunk/tcperutrpos:
  contents: "# slamdunk tcperutr"
  num_lines: 1
snippy/snippy:
  contents: snippy
  num_lines: 20
snippy/snippy-core:
  contents_re: ID\tLENGTH\tALIGNED\tUNALIGNED\tVARIANT\tHET\tMASKED\tLOWCOV
  num_lines: 1
snpeff:
  contents: SnpEff_version
  max_filesize: 5000000
snpsplit/old:
  contents: "Writing allele-flagged output file to:"
  num_lines: 2
snpsplit/new:
  fn: "*SNPsplit_report.yaml"
software_versions:
  fn_re: .+_mqc_versions\.(yaml|yml)
sompy:
  fn: "*.stats.csv"
  contents: ",sompyversion,sompycmd"
  num_lines: 2
sortmerna:
  contents: Minimal SW score based on E-value
spaceranger/count_html:
  - fn: "*.html"
    contents: '"command":"Space Ranger","subcommand":"count"'
    num_lines: 20
  - fn: "*.html"
    contents: '"command": "Space Ranger", "subcommand": "count"'
    num_lines: 20
stacks/gstacks:
  fn: gstacks.log.distribs
  contents: BEGIN effective_coverages_per_sample
stacks/populations:
  fn: populations.log.distribs
  contents: BEGIN missing_samples_per_loc_prefilters
stacks/sumstats:
  fn: "*.sumstats_summary.tsv"
  contents: "# Pop ID\tPrivate\tNum_Indv\tVar\tStdErr\tP\tVar"
  max_filesize: 1000000
star:
  fn: "*Log.final.out"
star/genecounts:
  fn: "*ReadsPerGene.out.tab"
supernova/report:
  fn: "*report*.txt"
  num_lines: 100
  contents: "- assembly checksum ="
supernova/summary:
  fn: summary.json
  num_lines: 120
  contents: '"lw_mean_mol_len":'
supernova/molecules:
  fn: histogram_molecules.json
  num_lines: 10
  contents: '"description": "molecules",'
supernova/kmers:
  fn: histogram_kmer_count.json
  num_lines: 10
  contents: '"description": "kmer_count",'
sylphtax:
  fn: "*.sylphmpa"
telseq:
  num_lines: 3
  contents: "ReadGroup\tLibrary\tSample\tTotal\tMapped\tDuplicates\tLENGTH_ESTIMATE"
theta2:
  fn: "*.BEST.results"
tophat:
  fn: "*align_summary.txt"
  shared: true
trim_galore:
  fn: "*_trimming_report.json"
trimmomatic:
  contents_re: ^Trimmomatic
truvari/bench:
  contents_re: .*truvari.* bench.*
  fn: log.txt
  num_lines: 10
umicollapse:
  num_lines: 100
  contents: "UMI collapsing finished in "
umitools/extract:
  contents: "# output generated by extract"
  num_lines: 100
umitools/dedup:
  contents: "# output generated by dedup"
  num_lines: 100
varscan2/mpileup2snp:
  contents: Only SNPs will be reported
  num_lines: 10
varscan2/mpileup2indel:
  contents: Only indels will be reported
  num_lines: 10
varscan2/mpileup2cns:
  contents: Only variants will be reported
  num_lines: 10
vcftools/relatedness2:
  fn: "*.relatedness2"
vcftools/tstv_by_count:
  fn: "*.TsTv.count"
vcftools/tstv_by_qual:
  fn: "*.TsTv.qual"
vcftools/tstv_summary:
  fn: "*.TsTv.summary"
vep/vep_html:
  fn: "*.html"
  contents: VEP summary
  num_lines: 10
  max_filesize: 1000000
vep/vep_txt:
  contents: "[VEP run statistics]"
  num_lines: 1
  max_filesize: 100000
verifybamid/selfsm:
  fn: "*.selfSM"
vg/stats:
  contents:
    - "Total perfect:"
    - "Total gapless (softclips allowed):"
    - "Total time:"
    - "Speed:"
  num_lines: 30
whatshap/stats:
  contents: "#sample\tchromosome\tfile_name\tvariants\tphased\tunphased\tsingletons"
  num_lines: 1
xenome:
  contents: "B\tG\tH\tM\tcount\tpercent\tclass"
  num_lines: 2
xengsort:
  contents: "# Xengsort classify"
  num_lines: 2
ataqv:
  fn: "*.json"
  contents: ataqv_version
  num_lines: 10
mosaicatcher:
  fn: "*.mosaicatcher_info_raw.txt"
```

</details>

**Example**:

```yaml
sp:
  fastqc/data:
    fn: fastqc_data.txt
  fastqc/zip:
    fn: "*_fastqc.zip"
```

## Plot Settings

### Rendering mode

#### `plots_force_flat`

**Type**: <code>bool</code> (default: `false`)

Render plots as static images instead of interactive Plotly. Useful for very large reports.

#### `plots_force_interactive`

**Type**: <code>bool</code> (default: `false`)

Force interactive plots even when MultiQC would normally fall back to flat images.

#### `plots_flat_numseries`

**Type**: <code>int</code> (default: `2000`)

If a plot has more than this many series, MultiQC switches it from interactive to flat image.

#### `plots_defer_loading_numseries`

**Type**: <code>int</code> (default: `100`)

Plots with more than this many series start collapsed. The user clicks a button to render them.

#### `num_datasets_plot_limit`

**Type**: <code>int</code> (default: `100`)

Deprecated. Use `plots_defer_loading_numseries` instead.

### Appearance

#### `plots_export_font_scale`

**Type**: <code>float</code> (default: `1.0`)

Multiplier applied to font sizes in exported plot images. Bump up for publication-quality output.

#### `plot_font_family`

**Type**: <code>str</code>

CSS font-family for plot text. Defaults to a system font stack.

#### `custom_plot_config`

**Type**: <code>Dict[str, Any]</code>

Override plot config options per plot. Top-level keys are plot IDs, values are option dicts.

**Example**:

```yaml
custom_plot_config:
  fastqc_per_base_sequence_quality_plot:
    title: "FastQC: Mean Quality Scores (custom)"
    yaxis:
      title: Phred score
```

#### `lineplot_number_of_points_to_hide_markers`

**Type**: <code>int</code> (default: `50`)

Hide individual data point markers in line plots once the total point count across samples exceeds this.

#### `barplot_legend_on_bottom`

**Type**: <code>bool</code> (default: `false`)

Place bar plot legends below the plot instead of to the side. Not recommended.

### Boxplot and violin

#### `boxplot_boxpoints`

**Type**: <code>Literal["outliers", "suspectedoutliers", "all", False]</code> (default: `"outliers"`)

How boxplot data points are drawn. Use false to hide individual points.

#### `box_min_threshold_outliers`

**Type**: <code>int</code> (default: `100`)

When a boxplot has more samples than this, only outlier points are drawn.

#### `box_min_threshold_no_points`

**Type**: <code>int</code> (default: `1000`)

When a boxplot has more samples than this, no individual points are drawn.

#### `violin_downsample_after`

**Type**: <code>int</code> (default: `2000`)

Start downsampling violin plot data once the sample count exceeds this. Keeps rendering snappy.

#### `violin_min_threshold_outliers`

**Type**: <code>int</code> (default: `100`)

When a violin plot has more samples than this, only outlier points are drawn.

#### `violin_min_threshold_no_points`

**Type**: <code>int</code> (default: `1000`)

When a violin plot has more samples than this, no individual points are drawn.

## Toolbox

### Highlighting

#### `highlight_patterns`

**Type**: <code>List[str]</code>

Substring (or regex) patterns. Matching samples are highlighted in plots and tables.

**Example**:

```yaml
highlight_patterns:
  - control
  - treated
```

#### `highlight_colors`

**Type**: <code>List[str]</code>

CSS colour for each entry in highlight_patterns, in the same order. Accepts hex (`#377eb8`), named colours (`red`), or any CSS colour function (`rgb(...)`, `hsl(...)`).

**Example**:

```yaml
highlight_colors:
  - "#377eb8"
  - "#e41a1c"
```

#### `highlight_regex`

**Type**: <code>bool</code> (default: `false`)

Treat highlight_patterns as regex instead of plain substring.

### Show/hide buttons

#### `show_hide_buttons`

**Type**: <code>List[str]</code>

Labels for the toolbox show/hide buttons. One per pattern set.

**Example**:

```yaml
show_hide_buttons:
  - Tumour samples
  - Normal samples
```

#### `show_hide_patterns`

**Type**: <code>List[Union[str, List[str]]]</code>

Patterns for each show/hide button. Each entry is a string or list of strings to match against sample names.

**Example**:

```yaml
show_hide_patterns:
  - - _T_
    - _tumour_
  - - _N_
    - _normal_
```

#### `show_hide_mode`

**Type**: <code>List[Literal["show", "hide", "show_re", "hide_re"]]</code>

Action for each show/hide button: 'show' (only show matches), 'hide' (hide matches), or their `_re` variants which signal regex patterns (set by the TSV loader).

**Example**:

```yaml
show_hide_mode:
  - show
  - show
```

#### `show_hide_regex`

**Type**: <code>List[Union[str, bool]]</code>

Whether each pattern set is treated as regex. List of bools aligned with show_hide_buttons.

**Example**:

```yaml
show_hide_regex:
  - false
  - false
```

## Table Settings

### General

#### `collapse_tables`

**Type**: <code>bool</code> (default: `true`)

Collapse module tables by default. Users click to expand.

#### `max_table_rows`

**Type**: <code>int</code> (default: `500`)

Tables larger than this many rows are rendered as a violin plot instead.

#### `max_configurable_table_columns`

**Type**: <code>int</code> (default: `200`)

Cap on the number of columns the user can toggle in the table-configure toolbox.

#### `decimalPoint_format`

**Type**: <code>str</code> (default: `"."`)

Decimal-point character used in formatted numbers. Defaults to `.`

**Example**:

```yaml
decimalPoint_format: ","
```

#### `thousandsSep_format`

**Type**: <code>str</code> (default: `" "`)

Thousands separator used in formatted numbers. Defaults to a single space, which is rendered as a small non-breaking space.

**Example**:

```yaml
thousandsSep_format: ","
```

### General Stats table

#### `general_stats_columns`

**Type**: <code>Dict[str, <a href="#generalstatsmoduleconfig">GeneralStatsModuleConfig</a>]</code>

Per-module overrides for General Stats columns. Top-level keys are module IDs.

**Example**:

```yaml
general_stats_columns:
  fastqc:
    columns:
      percent_duplicates:
        format: "{:,.1f}%"
        max: 100
        min: 0
        scale: RdYlGn-rev
        title: "% Dups"
```

#### `general_stats_helptext`

**Type**: <code>str</code>

Help text shown under the General Statistics heading at the top of the report.

#### `skip_generalstats`

**Type**: <code>bool</code> (default: `false`)

Hide the General Statistics table at the top of the report.

### Column overrides

#### `table_columns_name`

**Type**: <code>Dict[str, Union[str, Dict[str, str]]]</code>

Rename table columns. Top-level keys are module IDs, inner keys are column IDs, values are the new display name.

**Example**:

```yaml
table_columns_name:
  fastqc:
    percent_duplicates: "% Dups"
    percent_gc: "% GC"
```

#### `table_columns_placement`

**Type**: <code>Dict[str, Dict[str, float]]</code>

Reorder table columns. Top-level keys are module IDs, inner keys are column IDs, values are float sort weights (lower is further left).

**Example**:

```yaml
table_columns_placement:
  fastqc:
    percent_duplicates: 900
    percent_gc: 800
    total_sequences: 700
```

#### `table_columns_visible`

**Type**: <code>Dict[str, Union[bool, Dict[str, bool]]]</code>

Hide or show specific columns. Top-level keys are module IDs, values are either a bool (apply to all columns) or a dict mapping column ID to bool.

**Example**:

```yaml
table_columns_visible:
  fastqc: false
  samtools:
    error_rate: false
    raw_total_sequences: true
```

#### `custom_table_header_config`

**Type**: <code>Dict[str, Any]</code>

Override table column config. Same shape as custom_plot_config but for table headers.

**Example**:

```yaml
custom_table_header_config:
  general_stats_table:
    "% Dups":
      format: "{:,.1f}%"
      max: 100
      min: 0
```

### Conditional formatting

#### `table_cond_formatting_rules`

**Type**: <code>Dict[str, Dict[str, List[<a href="#condformattingrule">CondFormattingRule</a>]]]</code>

Conditional cell formatting. Nested dicts map table ID (or the literal 'all_columns') to colour ID to a list of rules. Each rule has exactly one operator: string operators (s_eq, s_ne, s_contains) compare case-insensitively; numeric operators (eq, ne, gt, lt, ge, le) cast both sides to float. See the customisation docs for the full grammar.

<details><summary>Default value</summary>

```yaml
all_columns:
  pass:
    - s_eq: pass
    - s_eq: "true"
    - s_eq: "yes"
    - s_eq: ok
  warn:
    - s_eq: warn
    - s_eq: unknown
  fail:
    - s_eq: fail
    - s_eq: "false"
    - s_eq: "no"
  male:
    - s_eq: male
    - s_eq: M
  female:
    - s_eq: female
    - s_eq: F
QCStatus:
  fail:
    - s_contains: fail
```

</details>

**Example**:

```yaml
table_cond_formatting_rules:
  all_columns:
    fail:
      - s_eq: fail
    pass:
      - s_eq: pass
      - s_eq: ok
    warn:
      - s_eq: warn
  mqc-generalstats-percent_duplicates:
    fail:
      - gt: 50
    warn:
      - gt: 20
```

#### `table_cond_formatting_colours`

**Type**: <code>List[Dict[str, str]]</code>

Background colours referenced by table_cond_formatting_rules. List of single-key dicts mapping a colour ID to a hex code.

<details><summary>Default value</summary>

```yaml
- blue: "#337ab7"
- lbue: "#5bc0de"
- pass: "#5cb85c"
- warn: "#f0ad4e"
- fail: "#d9534f"
- male: "#5bc0de"
- female: "#d9534f"
```

</details>

**Example**:

```yaml
table_cond_formatting_colours:
  - pass: "#5cb85c"
  - warn: "#f0ad4e"
  - fail: "#d9534f"
```

### Row merging

#### `table_sample_merge`

**Type**: <code>Dict[str, Union[str, <a href="#cleanpattern">CleanPattern</a>, List[Union[str, <a href="#cleanpattern">CleanPattern</a>]]]]</code>

Group samples by merging rows of supporting modules' tables, by collapsing samples that match a pattern. Keys are the merged group name; values are a clean-pattern entry (a string suffix, or a {type, pattern} dict) or a list of such entries.

**Examples**:

```yaml
table_sample_merge:
  R1: _1
  R2: _2
```

```yaml
table_sample_merge:
  R1:
    - _R1
    - pattern: "[_.-][rR]?1$"
      type: regex
  R2:
    - _R2
    - pattern: "[_.-][rR]?2$"
      type: regex
```

## Software Versions

### `software_versions`

**Type**: <code>Dict[str, Union[str, List[str], Dict[str, Union[str, List[str]]]]]</code>

Manually specify software versions for the Software Versions section. Top-level keys are group or software names. Values are a single version string, a list of version strings, or a dict mapping software name to a version string or list of version strings (when the group contains multiple tools).

**Examples**:

```yaml
software_versions:
  bwa: 0.7.17
  fastqc: 0.12.1
  samtools: "1.20"
```

```yaml
software_versions:
  quast:
    - 5.2.0
    - 5.1.0
```

```yaml
software_versions:
  samtools:
    htslib: "1.3"
    samtools: "1.11"
```

### `versions_table_group_header`

**Type**: <code>str</code> (default: `"Group"`)

Column header for the grouping column in the Software Versions table. Defaults to 'Group'.

### `disable_version_detection`

**Type**: <code>bool</code> (default: `false`)

Skip parsing software versions from module log files.

### `skip_versions_section`

**Type**: <code>bool</code> (default: `false`)

Hide the Software Versions section.

## Read & Base Counts

### Short reads

#### `read_count_multiplier`

**Type**: <code>float</code> (default: `1e-06`)

Multiplier applied to read counts before display. Default 0.000001 shows reads in millions.

**Example**:

```yaml
read_count_multiplier: 0.001
```

#### `read_count_prefix`

**Type**: <code>str</code> (default: `"M"`)

Suffix shown after formatted read counts, eg. 'M' for millions.

**Example**:

```yaml
read_count_prefix: K
```

#### `read_count_desc`

**Type**: <code>str</code> (default: `"millions"`)

Word used in plot/axis labels for read counts, eg. 'millions'.

**Examples**:

```yaml
read_count_desc: thousands
```

```yaml
read_count_desc: raw reads
```

### Long reads

#### `long_read_count_multiplier`

**Type**: <code>float</code> (default: `0.001`)

Multiplier for long-read counts. Default 0.001 shows counts in thousands.

**Example**:

```yaml
long_read_count_multiplier: 1.0e-06
```

#### `long_read_count_prefix`

**Type**: <code>str</code> (default: `"K"`)

Suffix shown after formatted long-read counts, eg. 'K' for thousands.

**Example**:

```yaml
long_read_count_prefix: M
```

#### `long_read_count_desc`

**Type**: <code>str</code> (default: `"thousands"`)

Word used in labels for long-read counts, eg. 'thousands'.

**Example**:

```yaml
long_read_count_desc: millions
```

### Bases

#### `base_count_multiplier`

**Type**: <code>float</code> (default: `1e-06`)

Multiplier for base counts. Default 0.000001 shows bases in megabases.

**Example**:

```yaml
base_count_multiplier: 0.001
```

#### `base_count_prefix`

**Type**: <code>str</code> (default: `"Mb"`)

Suffix shown after formatted base counts, eg. 'Mb' for megabases.

**Example**:

```yaml
base_count_prefix: Kb
```

#### `base_count_desc`

**Type**: <code>str</code> (default: `"millions"`)

Word used in labels for base counts, eg. 'megabases'.

**Example**:

```yaml
base_count_desc: kilobases
```

## AI Summary

### On/off

#### `ai_summary`

**Type**: <code>bool</code> (default: `false`)

Generate a short AI-written summary at the top of the report.

#### `ai_summary_full`

**Type**: <code>bool</code> (default: `false`)

Also generate a longer per-section AI summary. Requires ai_summary to be on.

#### `no_ai`

**Type**: <code>bool</code> (default: `false`)

Disable AI summaries entirely. Overrides ai_summary and ai_summary_full.

### Prompts

#### `ai_prompt_short`

**Type**: <code>str</code>

Custom prompt prepended to the short AI summary request. Use to steer tone, length, or focus.

**Example**:

```yaml
ai_prompt_short: Write the summary in one short paragraph aimed at a lab head, no
  jargon.
```

#### `ai_prompt_full`

**Type**: <code>str</code>

Custom prompt prepended to the full-section AI summary request.

**Example**:

```yaml
ai_prompt_full: Use bullet points and call out any sample that looks like an outlier.
```

### Privacy

#### `ai_anonymize_samples`

**Type**: <code>bool</code> (default: `false`)

Replace sample names with placeholders before sending data to the AI provider.

### Provider

#### `ai_provider`

**Type**: <code>Literal["seqera", "openai", "anthropic", "aws_bedrock", "custom"]</code> (default: `"seqera"`)

AI provider used for summaries. One of seqera, openai, anthropic, aws_bedrock, custom.

#### `ai_model`

**Type**: <code>str</code>

Model name. Provider-specific.

**Examples**:

```yaml
ai_model: gpt-4o
```

```yaml
ai_model: claude-sonnet-4-5.
```

#### `ai_custom_endpoint`

**Type**: <code>str</code>

Base URL for the 'custom' provider, eg. a self-hosted OpenAI-compatible API.

**Examples**:

```yaml
ai_custom_endpoint: http://localhost:11434/v1
```

```yaml
ai_custom_endpoint: https://api.example.com/v1
```

#### `ai_auth_type`

**Type**: <code>Literal["bearer", "api-key"]</code>

Authentication scheme used by the custom endpoint. 'bearer' sends an Authorization header, 'api-key' sends an api-key header.

#### `seqera_website`

**Type**: <code>str</code> (default: `"https://ai.seqera.io"`)

Base URL used for Seqera Platform links in the report.

#### `seqera_api_url`

**Type**: <code>str</code> (default: `"https://ai.seqera.io/v1/web"`)

Base URL for the Seqera Platform API. Defaults to the public instance.

### Tuning

#### `ai_retries`

**Type**: <code>int</code> (default: `3`)

Number of times to retry an AI request on transient errors.

#### `ai_extra_query_options`

**Type**: <code>Dict[str, Any]</code>

Extra request-body fields merged into the AI request payload (provider-specific).

**Example**:

```yaml
ai_extra_query_options:
  temperature: 0.3
  top_p: 0.9
```

#### `ai_custom_context_window`

**Type**: <code>int</code>

Override the model's context window in tokens. Set this if MultiQC's default for your model is wrong.

#### `ai_max_completion_tokens`

**Type**: <code>int</code>

Maximum completion tokens for OpenAI reasoning models.

#### `ai_reasoning_effort`

**Type**: <code>Literal["low", "medium", "high"]</code>

Reasoning effort for OpenAI reasoning models.

#### `ai_extended_thinking`

**Type**: <code>bool</code> (default: `false`)

Enable extended thinking on Anthropic Claude models that support it.

#### `ai_thinking_budget_tokens`

**Type**: <code>int</code>

Token budget for Anthropic extended thinking when enabled.

## MegaQC

### `megaqc_url`

**Type**: <code>str</code>

URL of a MegaQC instance to upload report data to after generation.

### `megaqc_access_token`

**Type**: <code>str</code>

Auth token for the MegaQC instance.

### `megaqc_timeout`

**Type**: <code>int</code> (default: `30`)

Upload timeout in seconds when posting to MegaQC.

### `megaqc_upload`

**Type**: <code>bool</code>

Upload report data to MegaQC after generation. Requires megaqc_url and megaqc_access_token.

## Performance & Debugging

### Profiling

#### `profile_runtime`

**Type**: <code>bool</code> (default: `false`)

Time each module and include the breakdown in the report.

#### `profile_memory`

**Type**: <code>bool</code> (default: `false`)

Track peak memory per module. Adds runtime overhead.

### Logging

#### `verbose`

**Type**: <code>bool</code> (default: `false`)

Print extra debug log messages to the terminal.

#### `no_ansi`

**Type**: <code>bool</code> (default: `false`)

Disable ANSI colour codes in terminal output.

#### `quiet`

**Type**: <code>bool</code> (default: `false`)

Suppress non-essential log messages.

### Linting

#### `strict`

**Type**: <code>bool</code> (default: `false`)

Treat module warnings as errors. Stricter than lint.

#### `lint`

**Type**: <code>bool</code> (default: `false`)

Deprecated. Run module linting and fail the build on issues. Used in MultiQC's own tests, rarely useful otherwise.

### Developer

#### `development`

**Type**: <code>bool</code> (default: `false`)

Enable developer-mode features such as live JS reloading. For internal use.

#### `report_readerrors`

**Type**: <code>bool</code> (default: `false`)

Surface file read errors in the log instead of silently skipping them.

#### `preserve_module_raw_data`

**Type**: <code>bool</code> (default: `false`)

Keep each module's raw parsed data in memory after report generation. Used by Python API consumers.

### Version check

#### `no_version_check`

**Type**: <code>bool</code> (default: `false`)

Skip the network check for newer MultiQC versions on startup.

#### `version_check_url`

**Type**: <code>str</code> (default: `"https://api.multiqc.info/version"`)

URL queried by MultiQC's own update check. Set to override the default endpoint.

## Special Types

### SearchPattern

Configuration for file search patterns used to find tool outputs.

The `SearchPattern` type is used in the `sp` configuration option to define patterns for finding and parsing tool output files.

Example:

```yaml
sp:
  fastqc:
    fn: "*_fastqc.zip"
  custom_tool:
    fn: "*.log"
    contents: "Started analysis"
```

Properties:

- **contents** (<code>Union[str, List[str]]</code>): File contents to match
- **contents_re** (<code>Union[str, List[str]]</code>): File contents regex pattern to match
- **exclude_contents** (<code>Union[str, List[str]]</code>): Exclude files containing this content
- **exclude_contents_re** (<code>Union[str, List[str]]</code>): Exclude files containing this regex content
- **exclude_fn** (<code>Union[str, List[str]]</code>): Exclude files matching this pattern
- **exclude_fn_re** (<code>Union[str, List[str]]</code>): Exclude files matching this regex pattern
- **fn** (<code>str</code>): Filename pattern to match
- **fn_re** (<code>str</code>): Filename regex pattern to match
- **max_filesize** (<code>int</code>): Maximum file size to process
- **num_lines** (<code>int</code>): Number of lines to search
- **shared** (<code>bool</code>): Allow file to be processed by multiple search patterns
- **skip** (<code>bool</code>): Skip this search pattern

### CleanPattern

Pattern for cleaning sample names.

The `CleanPattern` type is used in the `fn_clean_exts` and `extra_fn_clean_exts` configuration options to define patterns for cleaning sample names.

Example:

```yaml
fn_clean_exts:
  - type: truncate
    pattern: '_S\d+_L\d+'
  - type: regex
    pattern: '\d{4}-\d{2}-\d{2}'
```

Properties:

- **module** (<code>Union[str, List[str]]</code>): Module(s) to apply this pattern to
- **pattern** (<code>str</code>): Pattern to match
- **type** (<code>Literal["truncate", "remove", "regex", "regex_keep"]</code>): Type of pattern matching to use

### GeneralStatsModuleConfig

Per-module wrapper for General Stats column overrides.

The `GeneralStatsModuleConfig` type is the value of each module entry in the `general_stats_columns` configuration option. It has a single `columns` key mapping column IDs to `GeneralStatsColumnConfig` settings.

Example:

```yaml
general_stats_columns:
  fastqc:
    columns:
      percent_duplicates:
        title: "% Dups"
```

Properties:

- **columns** (<code>Dict[str, <a href="#generalstatscolumnconfig">GeneralStatsColumnConfig</a>]</code>): Columns to show in general stats table. Keys are column IDs.

### GeneralStatsColumnConfig

Configuration for columns in the general statistics table.

The `GeneralStatsColumnConfig` type is used in the `general_stats_columns` configuration option to customize the appearance and behavior of columns in the general statistics table.

Example:

```yaml
general_stats_columns:
  fastqc:
    columns:
      percent_duplicates:
        title: "% Dups"
        description: "Percentage of duplicate reads"
        scale: "RdYlGn-rev"
        max: 100
        min: 0
```

Properties:

- **ceiling** (<code>float</code>): Ceiling value
- **description** (<code>str</code>): Column description
- **floor** (<code>float</code>): Floor value
- **format** (<code>str</code>): Number format
- **hidden** (<code>bool</code>): Whether column is hidden by default
- **max** (<code>float</code>): Maximum value
- **min** (<code>float</code>): Minimum value
- **namespace** (<code>str</code>): Column namespace
- **placement** (<code>float</code>): Column placement order
- **scale** (<code>str</code>): Color scale
- **shared_key** (<code>str</code>): Shared key name
- **title** (<code>str</code>): Column title

### CondFormattingRule

One conditional-formatting comparison for a table cell.

Used in the `table_cond_formatting_rules` configuration option. Each rule is a dict with exactly one operator key paired with its comparison value. String operators (`s_eq`, `s_ne`, `s_contains`) compare case-insensitively; numeric operators (`eq`, `ne`, `gt`, `lt`, `ge`, `le`) cast both sides via `float()`.

Example:

```yaml
table_cond_formatting_rules:
  all_columns:
    pass:
      - s_eq: "pass"
    fail:
      - gt: 50
```

Properties:

- **eq** (<code>Union[float, int]</code>): Numeric equality
- **ge** (<code>Union[float, int]</code>): Greater than or equal to
- **gt** (<code>Union[float, int]</code>): Strictly greater than
- **le** (<code>Union[float, int]</code>): Less than or equal to
- **lt** (<code>Union[float, int]</code>): Strictly less than
- **ne** (<code>Union[float, int]</code>): Numeric inequality
- **s_contains** (<code>str</code>): Case-insensitive substring match
- **s_eq** (<code>str</code>): Case-insensitive string equality
- **s_ne** (<code>str</code>): Case-insensitive string inequality

### ModuleOverride

Per-module override values for `top_modules` and `module_order` entries.

Each entry in `top_modules` / `module_order` is either a module ID (string) or a single-key dict mapping the module ID to a `ModuleOverride` dict.

Example:

```yaml
module_order:
  - fastqc:
      name: "FastQC (trimmed)"
      anchor: "fastqc_trimmed"
      path_filters:
        - "*_trimmed*"
```

Properties:

- **anchor** (<code>str</code>): HTML/section anchor for this module run
- **comment** (<code>str</code>): Comment text rendered as markdown under the heading
- **custom_config** (<code>Dict[str, Any]</code>): Module-specific config values merged into config.<module_id>
- **doi** (<code>Union[str, List[str]]</code>): DOI or list of DOIs
- **extra** (<code>str</code>): Extra HTML appended after the intro
- **generalstats** (<code>bool</code>): Set to false to suppress this module's general-stats columns
- **href** (<code>Union[str, List[str]]</code>): Tool homepage URL, or list of URLs
- **info** (<code>str</code>): Intro text rendered as markdown under the section heading
- **name** (<code>str</code>): Display name for this module run
- **path_filters** (<code>Union[str, List[str]]</code>): Glob patterns restricting which files this module run sees
- **path_filters_exclude** (<code>Union[str, List[str]]</code>): Glob patterns excluding files from this module run

### SectionOrderOverride

Override dict accepted as a `report_section_order` value.

Each value in `report_section_order` is either the literal string `"remove"` (drops the section) or a `SectionOrderOverride` dict combining any of `order`, `before` and `after`.

Example:

```yaml
report_section_order:
  fastqc:
    order: -10
  custom_content-my-section:
    before: fastqc
  mod_section_2: remove
```

Properties:

- **after** (<code>str</code>): Section/module/anchor ID to position this entry after
- **before** (<code>str</code>): Section/module/anchor ID to position this entry before
- **order** (<code>int</code>): Explicit numeric order
