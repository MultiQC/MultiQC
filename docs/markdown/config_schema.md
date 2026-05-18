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

**Type**: `str`

Title shown at the top of the report and used in the page title.

#### `subtitle`

**Type**: `str`

Subtitle shown under the report title. Plain text only.

#### `intro_text`

**Type**: `str`

Paragraph shown under the title. Useful for adding context about the analysis.

#### `report_comment`

**Type**: `str`

Free-text comment shown at the top of the report. HTML is allowed.

**Example**:

```yaml
report_comment: This report was generated from the RNA-seq pipeline on 2024-08-21.
```

#### `report_header_info`

**Type**: `List[Dict[str, str]]`

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

**Type**: `bool` (default: `true`)

Show the absolute paths of analysed directories in the report header.

#### `show_analysis_time`

**Type**: `bool` (default: `true`)

Show the date and time the report was generated in the header.

## Report Appearance

### Template

#### `template`

**Type**: `str` (default: `"default"`)

Name of the report template. Built-in templates: default, original, simple, sections, gathered, geo, disco. Plugin packages can register additional templates via the `multiqc.templates.v1` entry point.

**Example**:

```yaml
template: default
```

#### `template_dark_mode`

**Type**: `bool` (default: `true`)

Enable the dark mode toggle in the report template.

#### `simple_output`

**Type**: `bool` (default: `false`)

Render a minimal HTML report without the toolbox or interactive widgets. Useful for very large reports.

### Logo

#### `custom_logo`

**Type**: `str`

Path to an image to show at the top of the report, replacing the MultiQC logo.

**Examples**:

```yaml
custom_logo: /path/to/logo.png
```

```yaml
custom_logo: ./assets/logo.svg
```

#### `custom_logo_dark`

**Type**: `str`

Path to an alternative logo for dark mode. Falls back to custom_logo if unset.

**Example**:

```yaml
custom_logo_dark: ./assets/logo_dark.svg
```

#### `custom_logo_url`

**Type**: `str`

URL the custom logo links to when clicked.

**Example**:

```yaml
custom_logo_url: https://www.scilifelab.se
```

#### `custom_logo_title`

**Type**: `str`

Tooltip text shown when hovering over the custom logo.

**Example**:

```yaml
custom_logo_title: Our institute name
```

#### `custom_logo_width`

**Type**: `int`

Logo width in pixels. Height scales proportionally.

**Example**:

```yaml
custom_logo_width: 200
```

### Branding

#### `custom_favicon`

**Type**: `str`

Path to a custom favicon image to show in the browser tab.

**Examples**:

```yaml
custom_favicon: /path/to/favicon.ico
```

```yaml
custom_favicon: ./assets/favicon.png
```

#### `custom_css_files`

**Type**: `List[str]`

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

**Type**: `Dict[str, Any]`

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

**Type**: `List[str]`

Extra module IDs whose output should be parsed as custom content.

#### `custom_data`

**Type**: `Dict[str, Any]`

Inline custom content data keyed by section ID. Companion to custom_content for users who prefer splitting the metadata and the data across two top-level keys.

### Module ordering

#### `top_modules`

**Type**: `List[Union[str, Dict[str, Dict[str, Any]]]]`

Module IDs to render before module_order. Useful for pinning a module to the top regardless of where it appears in module_order. Same shape as module_order entries.

**Example**:

```yaml
top_modules:
  - fastqc
  - cutadapt
```

#### `module_order`

**Type**: `List[Union[str, Dict[str, Dict[str, Any]]]]`

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

**Type**: `List[str]`

Module IDs to run. If set, only listed modules are processed (mirror of the --module CLI flag).

**Example**:

```yaml
run_modules:
  - fastqc
  - cutadapt
  - samtools
```

#### `exclude_modules`

**Type**: `List[str]`

Module IDs to skip (mirror of the --exclude CLI flag).

**Example**:

```yaml
exclude_modules:
  - fastqc
```

#### `remove_sections`

**Type**: `List[str]`

Module sections to hide. Use the section anchor as it appears in the URL.

**Example**:

```yaml
remove_sections:
  - fastqc_overrepresented_sequences
  - gatk-compare-overlap
```

#### `report_section_order`

**Type**: `Dict[str, Any]`

Reorder, group or hide report sections by ID. Values can be a position string ('before'/'after'), an explicit order number, or a dict of overrides. See the [customisation docs](https://docs.seqera.io/multiqc/reports/customisation#order-of-module-and-module-subsection-output) for the full grammar.

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

**Type**: `Dict[str, str]`

Markdown text shown under specific module sections. Keys are section anchors.

**Example**:

```yaml
section_comments:
  fastqc_overrepresented_sequences: "**This is** an important note about the overrepresented\
    \ sequences."
  samtools: Reviewed by *Phil* on 2024-08-21.
```

#### `section_status_checks`

**Type**: `Dict[str, Union[bool, Dict[str, bool]]]`

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

**Type**: `bool` (default: `false`)

Overwrite existing output files without prompting.

#### `output_fn_name`

**Type**: `str` (default: `"multiqc_report.html"`)

Filename for the generated HTML report. Defaults to multiqc_report.html.

#### `make_report`

**Type**: `bool` (default: `true`)

Generate the HTML report. Set to false to only produce data files.

### Data files

#### `make_data_dir`

**Type**: `bool` (default: `true`)

Write parsed data as files alongside the report.

#### `zip_data_dir`

**Type**: `bool` (default: `false`)

Compress the data directory into a single .zip file.

#### `data_dir_name`

**Type**: `str` (default: `"multiqc_data"`)

Name of the directory written alongside the report holding parsed data. Defaults to multiqc_data.

#### `data_format`

**Type**: `Literal["tsv", "csv", "json", "yaml"]` (default: `"tsv"`)

Format used when writing parsed data files.

#### `data_format_extensions`

**Type**: `Dict[str, str]` (default: `{"tsv":"txt","csv":"csv","json":"json","yaml":"yaml"}`)

Override the file extension used when writing each data format, eg. {tsv: txt} to write TSV as .txt.

**Example**:

```yaml
data_format_extensions:
  json: json
  tsv: txt
  yaml: yml
```

#### `parquet_format`

**Type**: `Literal["long", "wide"]` (default: `"long"`)

Parquet table layout. 'long' has rows of (sample_name, metric_name, val_raw, val_raw_type, val_str), easy to filter by metric. 'wide' uses one column per metric (prefixed with table name and namespace), easier for analytics but can hit column limits or mixed-type issues.

### Data dump

#### `data_dump_file`

**Type**: `bool` (default: `true`)

Write a single JSON file containing all parsed data, for re-running MultiQC later.

#### `data_dump_file_write_raw`

**Type**: `bool` (default: `true`)

Include raw values (before any normalisation or filtering) in the dumped JSON.

### Plot export

#### `export_plots`

**Type**: `bool` (default: `false`)

Save each plot as a static image (formats set by export_plot_formats).

#### `export_plot_formats`

**Type**: `List[Literal["png", "svg", "pdf"]]` (default: `["png","svg","pdf"]`)

Image formats to export when export_plots is on.

#### `export_plots_timeout`

**Type**: `int` (default: `60`)

Timeout for exporting each plot, in seconds.

#### `plots_dir_name`

**Type**: `str` (default: `"multiqc_plots"`)

Directory for exported plot images when export_plots is on. Defaults to multiqc_plots.

### PDF

#### `make_pdf`

**Type**: `bool` (default: `false`)

Also generate a PDF version of the report. Requires Pandoc to be installed.

#### `pandoc_template`

**Type**: `str`

Path to a Pandoc template used when exporting the report as PDF.

## Sample Names

### Prepend directory

#### `prepend_dirs`

**Type**: `bool` (default: `false`)

Prefix sample names with their parent directory. Useful when the same sample name occurs in multiple directories.

#### `prepend_dirs_depth`

**Type**: `int` (default: `0`)

How many parent directories to include. 0 means all the way to the root.

#### `prepend_dirs_sep`

**Type**: `str` (default: `" | "`)

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

**Type**: `bool` (default: `true`)

Apply the cleaning rules in fn_clean_exts and fn_clean_trim to sample names.

#### `extra_fn_clean_exts`

**Type**: `List[Union[str, CleanPattern]]`

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

**Type**: `List[str]`

Strings appended to the built-in trim list, without overriding defaults.

**Example**:

```yaml
extra_fn_clean_trim:
  - sample_
  - _processed
```

#### `fn_clean_exts`

**Type**: `List[Union[str, CleanPattern]]`

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

**Type**: `List[str]`

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

**Type**: `Union[bool, List[str]]` (default: `false`)

Use the source filename as the sample name instead of any name parsed from the log. Set to true for all modules, or to a list of module IDs / patterns to apply selectively.

### Ignore samples

#### `sample_names_ignore`

**Type**: `List[str]`

Glob patterns. Matching samples are dropped from the report.

**Example**:

```yaml
sample_names_ignore:
  - "*_temp"
  - control_*
```

#### `sample_names_ignore_re`

**Type**: `List[str]`

Regex patterns. Matching samples are dropped from the report.

**Example**:

```yaml
sample_names_ignore_re:
  - ^test_.*
  - .*_neg_ctrl$
```

#### `sample_names_only_include`

**Type**: `List[str]`

Glob patterns. If set, only matching samples are kept.

**Example**:

```yaml
sample_names_only_include:
  - RNA_*
  - Sample_??
```

#### `sample_names_only_include_re`

**Type**: `List[str]`

Regex patterns. If set, only matching samples are kept.

**Example**:

```yaml
sample_names_only_include_re:
  - ^WGS_[0-9]+$
```

### Rename and replace

#### `sample_names_rename`

**Type**: `List[List[str]]`

Toolbox rename pairs. Each entry is a [from, to] pair, grouped by the buttons in sample_names_rename_buttons.

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

**Type**: `List[str]`

Names of the toolbox buttons that switch between the rename groups defined in sample_names_rename.

**Example**:

```yaml
sample_names_rename_buttons:
  - Sample ID
  - Patient ID
  - Lane
```

#### `sample_names_replace`

**Type**: `Dict[str, str]`

Substring replacements applied to every sample name. Keys are matched, values are replacements.

**Example**:

```yaml
sample_names_replace:
  Sample_: S
  _001: ""
```

#### `sample_names_replace_complete`

**Type**: `bool` (default: `false`)

Replace the entire sample name when the key matches anywhere in it.

#### `sample_names_replace_exact`

**Type**: `bool` (default: `false`)

Only replace when the key matches the sample name exactly, not as a substring.

#### `sample_names_replace_regex`

**Type**: `bool` (default: `false`)

Treat keys in sample_names_replace as regex patterns.

## File Discovery

### Input source

#### `file_list`

**Type**: `bool` (default: `false`)

Treat the input path as a file containing a list of paths to scan, one per line.

#### `require_logs`

**Type**: `bool` (default: `false`)

Fail with an error if any module explicitly requested with `--module` has no log files found. Off by default, so missing inputs are skipped silently.

### Size limits

#### `log_filesize_limit`

**Type**: `int` (default: `50000000`)

Skip log files larger than this many bytes.

#### `filesearch_lines_limit`

**Type**: `int` (default: `1000`)

Stop reading a log file after this many lines.

### Skip patterns

#### `ignore_symlinks`

**Type**: `bool` (default: `false`)

Skip symlinked files and directories during the file search.

#### `ignore_images`

**Type**: `bool` (default: `true`)

Skip image files (PNG/JPEG/etc.) to avoid wasting time opening them.

#### `fn_ignore_dirs`

**Type**: `List[str]` (default: `["multiqc_data",".git","icarus_viewers","runs_per_reference","not_aligned","contigs_reports"]`)

Glob patterns for directory names to skip entirely during the file search.

**Example**:

```yaml
fn_ignore_dirs:
  - work
  - .nextflow
  - "*_logs"
```

#### `fn_ignore_paths`

**Type**: `List[str]` (default: `["*/work/??/??????????????????????????????","*/.snakemake","*/.singularity","*/__pycache__","*/site-packages/multiqc"]`)

Glob patterns for paths to skip during the file search.

**Example**:

```yaml
fn_ignore_paths:
  - "*/test_data/*"
  - "*/.snakemake/*"
```

#### `fn_ignore_files`

**Type**: `List[str]`

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

**Type**: `List[str]`

Module IDs whose log files may be matched by multiple modules during the search.

## Plot Settings

### Rendering mode

#### `plots_force_flat`

**Type**: `bool` (default: `false`)

Render plots as static images instead of interactive Plotly. Useful for very large reports.

#### `plots_force_interactive`

**Type**: `bool` (default: `false`)

Force interactive plots even when MultiQC would normally fall back to flat images.

#### `plots_flat_numseries`

**Type**: `int` (default: `2000`)

If a plot has more than this many series, MultiQC switches it from interactive to flat image.

#### `plots_defer_loading_numseries`

**Type**: `int` (default: `100`)

Plots with more than this many series start collapsed. The user clicks a button to render them.

#### `num_datasets_plot_limit`

**Type**: `int` (default: `100`)

Deprecated. Use `plots_defer_loading_numseries` instead.

### Appearance

#### `plots_export_font_scale`

**Type**: `float` (default: `1.0`)

Multiplier applied to font sizes in exported plot images. Bump up for publication-quality output.

#### `plot_font_family`

**Type**: `str`

CSS font-family for plot text. Defaults to a system font stack.

#### `custom_plot_config`

**Type**: `Dict[str, Any]`

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

**Type**: `int` (default: `50`)

Hide individual data point markers in line plots once the total point count across samples exceeds this.

#### `barplot_legend_on_bottom`

**Type**: `bool` (default: `false`)

Place bar plot legends below the plot instead of to the side. Not recommended.

### Boxplot and violin

#### `boxplot_boxpoints`

**Type**: `Literal["outliers", "suspectedoutliers", "all", False]` (default: `"outliers"`)

How boxplot data points are drawn. Use false to hide individual points.

#### `box_min_threshold_outliers`

**Type**: `int` (default: `100`)

When a boxplot has more samples than this, only outlier points are drawn.

#### `box_min_threshold_no_points`

**Type**: `int` (default: `1000`)

When a boxplot has more samples than this, no individual points are drawn.

#### `violin_downsample_after`

**Type**: `int` (default: `2000`)

Start downsampling violin plot data once the sample count exceeds this. Keeps rendering snappy.

#### `violin_min_threshold_outliers`

**Type**: `int` (default: `100`)

When a violin plot has more samples than this, only outlier points are drawn.

#### `violin_min_threshold_no_points`

**Type**: `int` (default: `1000`)

When a violin plot has more samples than this, no individual points are drawn.

## Toolbox

### Highlighting

#### `highlight_patterns`

**Type**: `List[str]`

Substring (or regex) patterns. Matching samples are highlighted in plots and tables.

**Example**:

```yaml
highlight_patterns:
  - control
  - treated
```

#### `highlight_colors`

**Type**: `List[str]`

Hex colour for each entry in highlight_patterns, in the same order.

**Example**:

```yaml
highlight_colors:
  - "#377eb8"
  - "#e41a1c"
```

#### `highlight_regex`

**Type**: `bool` (default: `false`)

Treat highlight_patterns as regex instead of plain substring.

### Show/hide buttons

#### `show_hide_buttons`

**Type**: `List[str]`

Labels for the toolbox show/hide buttons. One per pattern set.

**Example**:

```yaml
show_hide_buttons:
  - Tumour samples
  - Normal samples
```

#### `show_hide_patterns`

**Type**: `List[Union[str, List[str]]]`

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

**Type**: `List[str]`

Action for each show/hide button: 'show' (only show matches) or 'hide' (hide matches).

**Example**:

```yaml
show_hide_mode:
  - show
  - show
```

#### `show_hide_regex`

**Type**: `List[Union[str, bool]]`

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

**Type**: `bool` (default: `true`)

Collapse module tables by default. Users click to expand.

#### `max_table_rows`

**Type**: `int` (default: `500`)

Tables larger than this many rows are rendered as a violin plot instead.

#### `max_configurable_table_columns`

**Type**: `int` (default: `200`)

Cap on the number of columns the user can toggle in the table-configure toolbox.

#### `decimalPoint_format`

**Type**: `str` (default: `"."`)

Decimal-point character used in formatted numbers. Defaults to `.`

**Example**:

```yaml
decimalPoint_format: ","
```

#### `thousandsSep_format`

**Type**: `str` (default: `" "`)

Thousands separator used in formatted numbers. Defaults to a single space, which is rendered as a small non-breaking space.

**Example**:

```yaml
thousandsSep_format: ","
```

### General Stats table

#### `general_stats_columns`

**Type**: `Dict[str, GeneralStatsModuleConfig]`

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

**Type**: `str`

Help text shown under the General Statistics heading at the top of the report.

#### `skip_generalstats`

**Type**: `bool` (default: `false`)

Hide the General Statistics table at the top of the report.

### Column overrides

#### `table_columns_name`

**Type**: `Dict[str, Union[str, Dict[str, str]]]`

Rename table columns. Top-level keys are module IDs, inner keys are column IDs, values are the new display name.

**Example**:

```yaml
table_columns_name:
  fastqc:
    percent_duplicates: "% Dups"
    percent_gc: "% GC"
```

#### `table_columns_placement`

**Type**: `Dict[str, Dict[str, float]]`

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

**Type**: `Dict[str, Union[bool, Dict[str, bool]]]`

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

**Type**: `Dict[str, Any]`

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

**Type**: `Dict[str, Dict[str, List[Dict[str, Union[str, int, float]]]]]`

Conditional cell formatting. Nested dicts map table ID to column ID to a list of rules (eg. {s_eq: pass} matches an exact value). See the customisation docs for the full grammar.

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

**Type**: `List[Dict[str, str]]`

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

**Type**: `Dict[str, Union[str, Dict[str, Union[str, List[str]]], List[Union[str, Dict[str, Union[str, List[str]]]]]]]`

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

**Type**: `Dict[str, Any]`

Manually specify software versions for the Software Versions section. Top-level keys are tool names.

**Example**:

```yaml
software_versions:
  bwa: 0.7.17
  fastqc: 0.12.1
  samtools: "1.20"
```

### `versions_table_group_header`

**Type**: `str` (default: `"Group"`)

Column header for the grouping column in the Software Versions table. Defaults to 'Group'.

### `disable_version_detection`

**Type**: `bool` (default: `false`)

Skip parsing software versions from module log files.

### `skip_versions_section`

**Type**: `bool` (default: `false`)

Hide the Software Versions section.

## Read & Base Counts

### Short reads

#### `read_count_multiplier`

**Type**: `float` (default: `1e-06`)

Multiplier applied to read counts before display. Default 0.000001 shows reads in millions.

**Example**:

```yaml
read_count_multiplier: 0.001
```

#### `read_count_prefix`

**Type**: `str` (default: `"M"`)

Suffix shown after formatted read counts, eg. 'M' for millions.

**Example**:

```yaml
read_count_prefix: K
```

#### `read_count_desc`

**Type**: `str` (default: `"millions"`)

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

**Type**: `float` (default: `0.001`)

Multiplier for long-read counts. Default 0.001 shows counts in thousands.

**Example**:

```yaml
long_read_count_multiplier: 1.0e-06
```

#### `long_read_count_prefix`

**Type**: `str` (default: `"K"`)

Suffix shown after formatted long-read counts, eg. 'K' for thousands.

**Example**:

```yaml
long_read_count_prefix: M
```

#### `long_read_count_desc`

**Type**: `str` (default: `"thousands"`)

Word used in labels for long-read counts, eg. 'thousands'.

**Example**:

```yaml
long_read_count_desc: millions
```

### Bases

#### `base_count_multiplier`

**Type**: `float` (default: `1e-06`)

Multiplier for base counts. Default 0.000001 shows bases in megabases.

**Example**:

```yaml
base_count_multiplier: 0.001
```

#### `base_count_prefix`

**Type**: `str` (default: `"Mb"`)

Suffix shown after formatted base counts, eg. 'Mb' for megabases.

**Example**:

```yaml
base_count_prefix: Kb
```

#### `base_count_desc`

**Type**: `str` (default: `"millions"`)

Word used in labels for base counts, eg. 'megabases'.

**Example**:

```yaml
base_count_desc: kilobases
```

## AI Summary

### On/off

#### `ai_summary`

**Type**: `bool` (default: `false`)

Generate a short AI-written summary at the top of the report.

#### `ai_summary_full`

**Type**: `bool` (default: `false`)

Also generate a longer per-section AI summary. Requires ai_summary to be on.

#### `no_ai`

**Type**: `bool` (default: `false`)

Disable AI summaries entirely. Overrides ai_summary and ai_summary_full.

### Prompts

#### `ai_prompt_short`

**Type**: `str`

Custom prompt prepended to the short AI summary request. Use to steer tone, length, or focus.

**Example**:

```yaml
ai_prompt_short: Write the summary in one short paragraph aimed at a lab head, no
  jargon.
```

#### `ai_prompt_full`

**Type**: `str`

Custom prompt prepended to the full-section AI summary request.

**Example**:

```yaml
ai_prompt_full: Use bullet points and call out any sample that looks like an outlier.
```

### Privacy

#### `ai_anonymize_samples`

**Type**: `bool` (default: `false`)

Replace sample names with placeholders before sending data to the AI provider.

### Provider

#### `ai_provider`

**Type**: `Literal["seqera", "openai", "anthropic", "aws_bedrock", "custom"]` (default: `"seqera"`)

AI provider used for summaries. One of seqera, openai, anthropic, aws_bedrock, custom.

#### `ai_model`

**Type**: `str`

Model name. Provider-specific.

**Examples**:

```yaml
ai_model: gpt-4o
```

```yaml
ai_model: claude-sonnet-4-5.
```

#### `ai_custom_endpoint`

**Type**: `str`

Base URL for the 'custom' provider, eg. a self-hosted OpenAI-compatible API.

**Examples**:

```yaml
ai_custom_endpoint: http://localhost:11434/v1
```

```yaml
ai_custom_endpoint: https://api.example.com/v1
```

#### `ai_auth_type`

**Type**: `Literal["bearer", "api-key"]`

Authentication scheme used by the custom endpoint. 'bearer' sends an Authorization header, 'api-key' sends an api-key header.

#### `seqera_website`

**Type**: `str` (default: `"https://seqera.io"`)

Base URL used for Seqera Platform links in the report.

#### `seqera_api_url`

**Type**: `str` (default: `"https://intern.seqera.io"`)

Base URL for the Seqera Platform API. Defaults to the public instance.

### Tuning

#### `ai_retries`

**Type**: `int` (default: `3`)

Number of times to retry an AI request on transient errors.

#### `ai_extra_query_options`

**Type**: `Dict[str, Any]`

Extra request-body fields merged into the AI request payload (provider-specific).

**Example**:

```yaml
ai_extra_query_options:
  temperature: 0.3
  top_p: 0.9
```

#### `ai_custom_context_window`

**Type**: `int`

Override the model's context window in tokens. Set this if MultiQC's default for your model is wrong.

#### `ai_max_completion_tokens`

**Type**: `int`

Maximum completion tokens for OpenAI reasoning models.

#### `ai_reasoning_effort`

**Type**: `Literal["low", "medium", "high"]`

Reasoning effort for OpenAI reasoning models.

#### `ai_extended_thinking`

**Type**: `bool` (default: `false`)

Enable extended thinking on Anthropic Claude models that support it.

#### `ai_thinking_budget_tokens`

**Type**: `int`

Token budget for Anthropic extended thinking when enabled.

## MegaQC

### `megaqc_url`

**Type**: `str`

URL of a MegaQC instance to upload report data to after generation.

### `megaqc_access_token`

**Type**: `str`

Auth token for the MegaQC instance.

### `megaqc_timeout`

**Type**: `int` (default: `30`)

Upload timeout in seconds when posting to MegaQC.

### `megaqc_upload`

**Type**: `bool`

Upload report data to MegaQC after generation. Requires megaqc_url and megaqc_access_token.

## Performance & Debugging

### Profiling

#### `profile_runtime`

**Type**: `bool` (default: `false`)

Time each module and include the breakdown in the report.

#### `profile_memory`

**Type**: `bool` (default: `false`)

Track peak memory per module. Adds runtime overhead.

### Logging

#### `verbose`

**Type**: `bool` (default: `false`)

Print extra debug log messages to the terminal.

#### `no_ansi`

**Type**: `bool` (default: `false`)

Disable ANSI colour codes in terminal output.

#### `quiet`

**Type**: `bool` (default: `false`)

Suppress non-essential log messages.

### Linting

#### `strict`

**Type**: `bool` (default: `false`)

Treat module warnings as errors. Stricter than lint.

#### `lint`

**Type**: `bool` (default: `false`)

Deprecated. Run module linting and fail the build on issues. Used in MultiQC's own tests, rarely useful otherwise.

### Developer

#### `development`

**Type**: `bool` (default: `false`)

Enable developer-mode features such as live JS reloading. For internal use.

#### `report_readerrors`

**Type**: `bool` (default: `false`)

Surface file read errors in the log instead of silently skipping them.

#### `preserve_module_raw_data`

**Type**: `bool` (default: `false`)

Keep each module's raw parsed data in memory after report generation. Used by Python API consumers.

### Version check

#### `no_version_check`

**Type**: `bool` (default: `false`)

Skip the network check for newer MultiQC versions on startup.

#### `version_check_url`

**Type**: `str` (default: `"https://api.multiqc.info/version"`)

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

- **contents** (`Union[str, List[str]]`): File contents to match
- **contents_re** (`Union[str, List[str]]`): File contents regex pattern to match
- **exclude_contents** (`Union[str, List[str]]`): Exclude files containing this content
- **exclude_contents_re** (`Union[str, List[str]]`): Exclude files containing this regex content
- **exclude_fn** (`Union[str, List[str]]`): Exclude files matching this pattern
- **exclude_fn_re** (`Union[str, List[str]]`): Exclude files matching this regex pattern
- **fn** (`str`): Filename pattern to match
- **fn_re** (`str`): Filename regex pattern to match
- **max_filesize** (`int`): Maximum file size to process
- **num_lines** (`int`): Number of lines to search
- **shared** (`bool`): Allow file to be processed by multiple search patterns
- **skip** (`bool`): Skip this search pattern

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

- **module** (`Union[str, List[str]]`): Module(s) to apply this pattern to
- **pattern** (`str`): Pattern to match
- **type** (`Literal["truncate", "remove", "regex", "regex_keep"]`): Type of pattern matching to use

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

- **ceiling** (`float`): Ceiling value
- **description** (`str`): Column description
- **floor** (`float`): Floor value
- **format** (`str`): Number format
- **hidden** (`bool`): Whether column is hidden by default
- **max** (`float`): Maximum value
- **min** (`float`): Minimum value
- **namespace** (`str`): Column namespace
- **placement** (`float`): Column placement order
- **scale** (`str`): Color scale
- **shared_key** (`str`): Shared key name
- **title** (`str`): Column title
