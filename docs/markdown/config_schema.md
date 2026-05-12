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

## Report Appearance

### custom_css_files

**Type**: `Optional[List[str]]` (default: `[]`)

Paths to additional CSS files to inline into the report. Useful for branding overrides.

**Example**:

```yaml
custom_css_files:
  - ./assets/custom.css
  - ./assets/branding.css
```

### custom_logo

**Type**: `Optional[str]` (default: `None`)

Path to an image to show at the top of the report, replacing the MultiQC logo.

**Examples**:

```yaml
custom_logo: /path/to/logo.png
```

```yaml
custom_logo: ./assets/logo.svg
```

### custom_logo_title

**Type**: `Optional[str]` (default: `None`)

Tooltip text shown when hovering over the custom logo.

**Example**:

```yaml
custom_logo_title: Our institute name
```

### custom_logo_url

**Type**: `Optional[str]` (default: `None`)

URL the custom logo links to when clicked.

**Example**:

```yaml
custom_logo_url: https://www.scilifelab.se
```

### intro_text

**Type**: `Optional[str]` (default: `None`)

Paragraph shown under the title. Useful for adding context about the analysis.

### report_comment

**Type**: `Optional[str]` (default: `None`)

Free-text comment shown at the top of the report. HTML is allowed.

**Example**:

```yaml
report_comment: This report was generated from the RNA-seq pipeline on 2024-08-21.
```

### report_header_info

**Type**: `Optional[List[Dict[str, str]]]` (default: `None`)

Extra key/value pairs shown in the report header, eg. contact name, run ID, pipeline version. Each list item is a single-key dictionary.

**Example**:

```yaml
report_header_info:
  - Contact E-mail: phil.ewels@seqera.io
  - Application Type: RNA-seq
  - Project Type: Application
  - Sequencing Platform: HiSeq 2500 High Output V4
```

### show_analysis_paths

**Type**: `Optional[bool]` (default: `true`)

Show the absolute paths of analysed directories in the report header.

### show_analysis_time

**Type**: `Optional[bool]` (default: `true`)

Show the date and time the report was generated in the header.

### simple_output

**Type**: `Optional[bool]` (default: `false`)

Render a minimal HTML report without the toolbox or interactive widgets. Useful for very large reports.

### subtitle

**Type**: `Optional[str]` (default: `None`)

Subtitle shown under the report title. Plain text only.

### template

**Type**: `Optional[Literal["default", "original", "simple", "sections", "gathered", "geo", "disco"]]` (default: `"default"`)

Name of the report template.

### title

**Type**: `Optional[str]` (default: `None`)

Title shown at the top of the report and used in the page title.

## Output Options

### data_dir_name

**Type**: `Optional[str]` (default: `"multiqc_data"`)

Name of the directory written alongside the report holding parsed data. Defaults to multiqc_data.

### data_dump_file

**Type**: `Optional[bool]` (default: `true`)

Write a single JSON file containing all parsed data, for re-running MultiQC later.

### data_dump_file_write_raw

**Type**: `Optional[bool]` (default: `true`)

Include raw values (before any normalisation or filtering) in the dumped JSON.

### data_format

**Type**: `Optional[Literal["tsv", "csv", "json", "yaml"]]` (default: `"tsv"`)

Format used when writing parsed data files.

### export_plots

**Type**: `Optional[bool]` (default: `false`)

Save each plot as a static image (formats set by export_plot_formats).

### export_plots_timeout

**Type**: `Optional[int]` (default: `60`)

Timeout for exporting each plot, in seconds.

### force

**Type**: `Optional[bool]` (default: `false`)

Overwrite existing output files without prompting.

### make_data_dir

**Type**: `Optional[bool]` (default: `true`)

Write parsed data as files alongside the report.

### make_pdf

**Type**: `Optional[bool]` (default: `false`)

Also generate a PDF version of the report. Requires Pandoc to be installed.

### make_report

**Type**: `Optional[bool]` (default: `true`)

Generate the HTML report. Set to false to only produce data files.

### output_fn_name

**Type**: `Optional[str]` (default: `"multiqc_report.html"`)

Filename for the generated HTML report. Defaults to multiqc_report.html.

### plots_dir_name

**Type**: `Optional[str]` (default: `"multiqc_plots"`)

Directory for exported plot images when export_plots is on. Defaults to multiqc_plots.

### zip_data_dir

**Type**: `Optional[bool]` (default: `false`)

Compress the data directory into a single .zip file.

## MegaQC Integration

### megaqc_access_token

**Type**: `Optional[str]` (default: `None`)

Auth token for the MegaQC instance.

### megaqc_timeout

**Type**: `Optional[int]` (default: `30`)

Upload timeout in seconds when posting to MegaQC.

### megaqc_url

**Type**: `Optional[str]` (default: `None`)

URL of a MegaQC instance to upload report data to after generation.

## AI Summary

### ai_anonymize_samples

**Type**: `Optional[bool]` (default: `false`)

Replace sample names with placeholders before sending data to the AI provider.

### ai_auth_type

**Type**: `Optional[Literal["bearer", "api-key"]]` (default: `None`)

Authentication scheme used by the custom endpoint. 'bearer' sends an Authorization header, 'api-key' sends an api-key header.

### ai_custom_context_window

**Type**: `Optional[str]` (default: `None`)

Override the model's context window in tokens. Set this if MultiQC's default for your model is wrong.

### ai_custom_endpoint

**Type**: `Optional[str]` (default: `None`)

Base URL for the 'custom' provider, eg. a self-hosted OpenAI-compatible API.

**Examples**:

```yaml
ai_custom_endpoint: http://localhost:11434/v1
```

```yaml
ai_custom_endpoint: https://api.example.com/v1
```

### ai_extra_query_options

**Type**: `Optional[str]` (default: `None`)

Extra URL query parameters appended to AI requests. Format: key1=val1&key2=val2.

**Example**:

```yaml
ai_extra_query_options: temperature=0.3&top_p=0.9
```

### ai_model

**Type**: `Optional[str]` (default: `None`)

Model name, eg. gpt-4o or claude-sonnet-4-5. Provider-specific.

### ai_prompt_full

**Type**: `Optional[str]` (default: `None`)

Custom prompt prepended to the full-section AI summary request.

**Example**:

```yaml
ai_prompt_full: Use bullet points and call out any sample that looks like an outlier.
```

### ai_prompt_short

**Type**: `Optional[str]` (default: `None`)

Custom prompt prepended to the short AI summary request. Use to steer tone, length, or focus.

**Example**:

```yaml
ai_prompt_short: Write the summary in one short paragraph aimed at a lab head, no
  jargon.
```

### ai_provider

**Type**: `Optional[Literal["seqera", "openai", "anthropic", "aws_bedrock", "custom"]]` (default: `"seqera"`)

AI provider used for summaries. One of seqera, openai, anthropic, aws_bedrock, custom.

### ai_retries

**Type**: `Optional[int]` (default: `3`)

Number of times to retry an AI request on transient errors.

### ai_summary

**Type**: `Optional[bool]` (default: `false`)

Generate a short AI-written summary at the top of the report.

### ai_summary_full

**Type**: `Optional[bool]` (default: `false`)

Also generate a longer per-section AI summary. Requires ai_summary to be on.

### no_ai

**Type**: `Optional[bool]` (default: `false`)

Disable AI summaries entirely. Overrides ai_summary and ai_summary_full.

## Seqera Integration

### seqera_api_url

**Type**: `Optional[str]` (default: `"https://intern.seqera.io"`)

Base URL for the Seqera Platform API. Defaults to the public instance.

### seqera_website

**Type**: `Optional[str]` (default: `"https://seqera.io"`)

Base URL used for Seqera Platform links in the report.

## Plot Settings

### barplot_legend_on_bottom

**Type**: `Optional[bool]` (default: `false`)

Place bar plot legends below the plot instead of to the side. Not recommended.

### lineplot_number_of_points_to_hide_markers

**Type**: `Optional[int]` (default: `50`)

Hide individual data point markers in line plots once the total point count across samples exceeds this.

### plots_defer_loading_numseries

**Type**: `Optional[int]` (default: `100`)

Plots with more than this many series start collapsed. The user clicks a button to render them.

### plots_export_font_scale

**Type**: `Optional[float]` (default: `1.0`)

Multiplier applied to font sizes in exported plot images. Bump up for publication-quality output.

### plots_flat_numseries

**Type**: `Optional[int]` (default: `2000`)

If a plot has more than this many series, MultiQC switches it from interactive to flat image.

### plots_force_flat

**Type**: `Optional[bool]` (default: `false`)

Render plots as static images instead of interactive Plotly. Useful for very large reports.

### plots_force_interactive

**Type**: `Optional[bool]` (default: `false`)

Force interactive plots even when MultiQC would normally fall back to flat images.

### violin_downsample_after

**Type**: `Optional[int]` (default: `2000`)

Start downsampling violin plot data once the sample count exceeds this. Keeps rendering snappy.

### violin_min_threshold_no_points

**Type**: `Optional[int]` (default: `1000`)

When a violin plot has more samples than this, no individual points are drawn.

### violin_min_threshold_outliers

**Type**: `Optional[int]` (default: `100`)

When a violin plot has more samples than this, only outlier points are drawn.

## Table Settings

### collapse_tables

**Type**: `Optional[bool]` (default: `true`)

Collapse module tables by default. Users click to expand.

### decimalPoint_format

**Type**: `Optional[str]` (default: `None`)

Decimal-point character used in formatted numbers, eg. '.' (default) or ','.

**Example**:

```yaml
decimalPoint_format: ","
```

### max_configurable_table_columns

**Type**: `Optional[int]` (default: `200`)

Cap on the number of columns the user can toggle in the table-configure toolbox.

### max_table_rows

**Type**: `Optional[int]` (default: `500`)

Tables larger than this many rows are rendered as a violin plot instead.

### thousandsSep_format

**Type**: `Optional[str]` (default: `None`)

Thousands separator used in formatted numbers, eg. ',' (default), ' ', or '.

**Examples**:

```yaml
thousandsSep_format: " "
```

```yaml
thousandsSep_format: "'"
```

## Sample Names

### extra_fn_clean_exts

**Type**: `Optional[List[Union[str, CleanPattern]]]` (default: `None`)

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

### extra_fn_clean_trim

**Type**: `Optional[List[str]]` (default: `None`)

Strings appended to the built-in trim list, without overriding defaults.

**Example**:

```yaml
extra_fn_clean_trim:
  - sample_
  - _processed
```

### fn_clean_exts

**Type**: `Optional[List[Union[str, CleanPattern]]]` (default: `[".gz",".fastq",".fq",".bam",".cram",".sam",".sra",".vcf",".dat","_tophat",".pbmarkdup.log",".log",".stderr",".out",".spp",".fa",".fasta",".png",".jpg",".jpeg",".html","Log.final","ReadsPerGene",".flagstat","_star_aligned","_fastqc",".hicup",".counts","_counts",".txt",".tsv",".csv",".aligned","Aligned",".merge",".deduplicated",".dedup",".clean",".sorted",".report","| stdin",".geneBodyCoverage",".inner_distance_freq",".junctionSaturation_plot.r",".pos.DupRate.xls",".GC.xls","_slamdunk","_bismark",".conpair",".concordance",".contamination",".BEST.results","_peaks.xls",".relatedness",".cnt",".aqhist",".bhist",".bincov",".bqhist",".covhist",".covstats",".ehist",".gchist",".idhist",".ihist",".indelhist",".lhist",".mhist",".qahist",".qchist",".qhist",".rpkm",".selfSM",".extendedFrags","_SummaryStatistics",".purple.purity",".purple.qc",".trim",".bowtie2",".mkD",".highfreq",".lowfreq",".consensus",".snpEff",".snpeff",".scaffolds",".contigs",".kraken2",".ccurve",".hisat2","_duprate",".markdup",".read_distribution",".junction_annotation",".infer_experiment",".biotype",".ivar",".mpileup",".primer_trim",".mapped",".vep","_vep","ccs","_NanoStats",".cutadapt",".qcML",".mosdepth","_gopeaks",".readCounts",".wgs_contig_mean_cov","_overall_mean_cov","_coverage_metrics",".wgs_fine_hist",".wgs_coverage_metrics",".wgs_hist",".vc_metrics",".gvcf_metrics",".ploidy_estimation_metrics","_overall_mean_cov",".fragment_length_hist",".mapping_metrics",".gc_metrics",".trimmer_metrics",".time_metrics",".quant_metrics",".quant.metrics",".quant.transcript_coverage",".scRNA_metrics",".scRNA.metrics",".scATAC_metrics",".scATAC.metrics",".fastqc_metrics",".labels",".bammetrics.metrics",".filter_summary",".cluster_report",".error.spl",".error.grp",".vgstats","_mapq_table","_strand_table","_isize_table","_dup_report","_cv_table","_covdist_all","_covdist_q40","_CpGRetention","_CpHRetentionByReadPos","_totalBaseConversionRate","_totalReadConversionRate",".sylphmpa","_qual"]`)

Extensions stripped from sample names, eg. .gz, .fastq. Replaces the built-in list.

**Example**:

```yaml
fn_clean_exts:
  - .gz
  - .fastq
  - .bam
  - pattern: _S\d+_L\d+
    type: regex
```

### fn_clean_sample_names

**Type**: `Optional[bool]` (default: `true`)

Apply the cleaning rules in fn_clean_exts and fn_clean_trim to sample names.

### fn_clean_trim

**Type**: `Optional[List[str]]` (default: `[".",":","_","-",".r","_val",".idxstats","_trimmed",".trimmed",".csv",".yaml",".yml",".json","_mqc","short_summary_","_summary",".summary",".align",".h5","_matrix",".stats",".hist",".phased",".tar","runs_",".qc"]`)

Strings trimmed from the start or end of sample names. Replaces the built-in list.

**Example**:

```yaml
fn_clean_trim:
  - _R1
  - _R2
  - _001
```

### prepend_dirs

**Type**: `Optional[bool]` (default: `false`)

Prefix sample names with their parent directory. Useful when the same sample name occurs in multiple folders.

### prepend_dirs_depth

**Type**: `Optional[int]` (default: `0`)

How many parent directories to include. 0 means all the way to the root.

### prepend_dirs_sep

**Type**: `Optional[str]` (default: `" | "`)

String inserted between directory names and the sample name. Defaults to '|'.

**Examples**:

```yaml
prepend_dirs_sep: _
```

```yaml
prepend_dirs_sep: " - "
```

### sample_names_ignore

**Type**: `Optional[List[str]]` (default: `[]`)

Glob patterns. Matching samples are dropped from the report.

**Example**:

```yaml
sample_names_ignore:
  - "*_temp"
  - control_*
```

### sample_names_ignore_re

**Type**: `Optional[List[str]]` (default: `[]`)

Regex patterns. Matching samples are dropped from the report.

**Example**:

```yaml
sample_names_ignore_re:
  - ^test_.*
  - .*_neg_ctrl$
```

### sample_names_only_include

**Type**: `Optional[List[str]]` (default: `[]`)

Glob patterns. If set, only matching samples are kept.

**Example**:

```yaml
sample_names_only_include:
  - RNA_*
  - Sample_??
```

### sample_names_only_include_re

**Type**: `Optional[List[str]]` (default: `[]`)

Regex patterns. If set, only matching samples are kept.

**Example**:

```yaml
sample_names_only_include_re:
  - ^WGS_[0-9]+$
```

### sample_names_rename

**Type**: `Optional[List[List[str]]]` (default: `[]`)

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

### sample_names_rename_buttons

**Type**: `Optional[List[str]]` (default: `[]`)

Names of the toolbox buttons that switch between the rename groups defined in sample_names_rename.

**Example**:

```yaml
sample_names_rename_buttons:
  - Sample ID
  - Patient ID
  - Lane
```

### sample_names_replace

**Type**: `Optional[Dict[str, str]]` (default: `{}`)

Substring replacements applied to every sample name. Keys are matched, values are replacements.

**Example**:

```yaml
sample_names_replace:
  Sample_: S
  _001: ""
```

### sample_names_replace_complete

**Type**: `Optional[bool]` (default: `false`)

Replace the entire sample name when the key matches anywhere in it.

### sample_names_replace_exact

**Type**: `Optional[bool]` (default: `false`)

Only replace when the key matches the sample name exactly, not as a substring.

### sample_names_replace_regex

**Type**: `Optional[bool]` (default: `false`)

Treat keys in sample_names_replace as regex patterns.

### use_filename_as_sample_name

**Type**: `Optional[Union[bool, List[str]]]` (default: `false`)

Use the source filename as the sample name instead of any name parsed from the log. Set to true for all modules, or to a list of module IDs / patterns to apply selectively.

## Toolbox

### highlight_colors

**Type**: `Optional[List[str]]` (default: `[]`)

Hex colour for each entry in highlight_patterns, in the same order.

**Example**:

```yaml
highlight_colors:
  - "#377eb8"
  - "#e41a1c"
```

### highlight_patterns

**Type**: `Optional[List[str]]` (default: `[]`)

Substring (or regex) patterns. Matching samples are highlighted in plots and tables.

**Example**:

```yaml
highlight_patterns:
  - control
  - treated
```

### highlight_regex

**Type**: `Optional[bool]` (default: `false`)

Treat highlight_patterns as regex instead of plain substring.

### show_hide_buttons

**Type**: `Optional[List[str]]` (default: `[]`)

Labels for the toolbox show/hide buttons. One per pattern set.

**Example**:

```yaml
show_hide_buttons:
  - Tumour samples
  - Normal samples
```

### show_hide_mode

**Type**: `Optional[List[str]]` (default: `[]`)

Action for each show/hide button: 'show' (only show matches) or 'hide' (hide matches).

**Example**:

```yaml
show_hide_mode:
  - show
  - show
```

### show_hide_patterns

**Type**: `Optional[List[Union[str, List[str]]]]` (default: `[]`)

Patterns for each show/hide button. Each entry is a string or list of strings to match against sample names.

**Example**:

```yaml
show_hide_patterns:
  - - _T_
    - _tumour_
  - - _N_
    - _normal_
```

### show_hide_regex

**Type**: `Optional[List[Union[str, bool]]]` (default: `[]`)

Whether each pattern set is treated as regex. List of bools aligned with show_hide_buttons.

**Example**:

```yaml
show_hide_regex:
  - false
  - false
```

## Performance & Debugging

### development

**Type**: `Optional[bool]` (default: `false`)

Enable developer-mode features such as live JS reloading. Internal use.

### filesearch_lines_limit

**Type**: `Optional[int]` (default: `1000`)

Stop reading a log file after this many lines.

### lint

**Type**: `Optional[bool]` (default: `false`)

Run module linting and fail the build on issues. Used in MultiQC's own tests, rarely useful otherwise.

### log_filesize_limit

**Type**: `Optional[int]` (default: `50000000`)

Skip log files larger than this many bytes.

### no_ansi

**Type**: `Optional[bool]` (default: `false`)

Disable ANSI colour codes in terminal output.

### profile_memory

**Type**: `Optional[bool]` (default: `false`)

Track peak memory per module. Adds runtime overhead.

### profile_runtime

**Type**: `Optional[bool]` (default: `false`)

Time each module and include the breakdown in the report.

### quiet

**Type**: `Optional[bool]` (default: `false`)

Suppress non-essential log messages.

### report_readerrors

**Type**: `Optional[bool]` (default: `false`)

Surface file read errors in the log instead of silently skipping them.

### strict

**Type**: `Optional[bool]` (default: `false`)

Treat module warnings as errors. Stricter than lint.

### verbose

**Type**: `Optional[bool]` (default: `false`)

Print extra debug log messages to the terminal.

## File Discovery

### filesearch_file_shared

**Type**: `Optional[List[str]]` (default: `[]`)

Module IDs whose log files may be matched by multiple modules during the search.

### fn_ignore_dirs

**Type**: `Optional[List[str]]` (default: `["multiqc_data",".git","icarus_viewers","runs_per_reference","not_aligned","contigs_reports"]`)

Glob patterns for directory names to skip entirely during the file search.

**Example**:

```yaml
fn_ignore_dirs:
  - work
  - .nextflow
  - "*_logs"
```

### fn_ignore_paths

**Type**: `Optional[List[str]]` (default: `["*/work/??/??????????????????????????????","*/.snakemake","*/.singularity","*/__pycache__","*/site-packages/multiqc"]`)

Glob patterns for paths to skip during the file search.

**Example**:

```yaml
fn_ignore_paths:
  - "*/test_data/*"
  - "*/.snakemake/*"
```

### ignore_images

**Type**: `Optional[bool]` (default: `true`)

Skip image files (PNG/JPEG/etc.) to avoid wasting time opening them.

### ignore_symlinks

**Type**: `Optional[bool]` (default: `false`)

Skip symlinked files and directories during the file search.

### require_logs

**Type**: `Optional[bool]` (default: `false`)

Fail with an error if any module explicitly requested with --module has no log files found. Off by default, so missing inputs are skipped silently.

## Other

### ai_extended_thinking

**Type**: `Optional[bool]` (default: `false`)

Enable extended thinking on Anthropic Claude models that support it.

### ai_max_completion_tokens

**Type**: `Optional[int]` (default: `None`)

Maximum completion tokens for OpenAI reasoning models.

### ai_reasoning_effort

**Type**: `Optional[Literal["low", "medium", "high"]]` (default: `None`)

Reasoning effort for OpenAI reasoning models.

### ai_thinking_budget_tokens

**Type**: `Optional[int]` (default: `None`)

Token budget for Anthropic extended thinking when enabled.

### base_count_desc

**Type**: `Optional[str]` (default: `"millions"`)

Word used in labels for base counts, eg. 'gigabases'.

**Examples**:

```yaml
base_count_desc: megabases
```

```yaml
base_count_desc: kilobases
```

### base_count_multiplier

**Type**: `Optional[float]` (default: `1e-06`)

Multiplier for base counts. Default 0.000000001 shows bases in gigabases.

**Examples**:

```yaml
base_count_multiplier: 1.0e-06
```

```yaml
base_count_multiplier: 0.001
```

### base_count_prefix

**Type**: `Optional[str]` (default: `"Mb"`)

Suffix shown after formatted base counts, eg. 'Gb' for gigabases.

**Examples**:

```yaml
base_count_prefix: Mb
```

```yaml
base_count_prefix: Kb
```

### box_min_threshold_no_points

**Type**: `Optional[int]` (default: `1000`)

When a boxplot has more samples than this, no individual points are drawn.

### box_min_threshold_outliers

**Type**: `Optional[int]` (default: `100`)

When a boxplot has more samples than this, only outlier points are drawn.

### boxplot_boxpoints

**Type**: `Optional[Literal["outliers", "suspectedoutliers", "all", False]]` (default: `"outliers"`)

How boxplot data points are drawn. Use false to hide individual points.

### custom_content

**Type**: `Optional[Dict[str, Any]]` (default: `{"order":[]}`)

Embed arbitrary plots, tables or text in the report. See the Custom Content docs for the full structure.

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

### custom_logo_dark

**Type**: `Optional[str]` (default: `None`)

Path to an alternative logo for dark mode. Falls back to custom_logo if unset.

**Example**:

```yaml
custom_logo_dark: ./assets/logo_dark.svg
```

### custom_logo_width

**Type**: `Optional[int]` (default: `None`)

Logo width in pixels. Height scales proportionally.

**Example**:

```yaml
custom_logo_width: 200
```

### custom_plot_config

**Type**: `Optional[Dict[str, Any]]` (default: `{}`)

Override plot config options per plot. Top-level keys are plot IDs, values are option dicts.

**Example**:

```yaml
custom_plot_config:
  fastqc_per_base_sequence_quality_plot:
    title: "FastQC: Mean Quality Scores (custom)"
    yaxis:
      title: Phred score
```

### custom_table_header_config

**Type**: `Optional[Dict[str, Any]]` (default: `{}`)

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

### data_format_extensions

**Type**: `Optional[Dict[str, str]]` (default: `{"tsv":"txt","csv":"csv","json":"json","yaml":"yaml"}`)

Override the file extension used when writing each data format, eg. {tsv: txt} to write TSV as .txt.

**Example**:

```yaml
data_format_extensions:
  json: json
  tsv: txt
  yaml: yml
```

### disable_version_detection

**Type**: `Optional[bool]` (default: `false`)

Skip parsing software versions from module log files.

### export_plot_formats

**Type**: `Optional[List[Literal["png", "svg", "pdf"]]]` (default: `["png","svg","pdf"]`)

Image formats to export when export_plots is on.

### file_list

**Type**: `Optional[bool]` (default: `false`)

Treat the input path as a file containing a list of paths to scan, one per line.

### fn_ignore_files

**Type**: `Optional[List[str]]` (default: `[".DS_Store",".py[cod]","*.bam","*.bai","*.sam","*.fq.gz","*.fastq.gz","*.fq","*.fastq","*.fa","*.gtf","*.bed","*.vcf","*.tbi","*.txt.gz","*.pdf","*.md5","*.parquet","*[!s][!u][!m][!_\\.m][!mva][!qer][!cpy].html","multiqc_data.json","*.gam","*.gamp","*.jar"]`)

Glob patterns for file names to skip during the file search.

**Example**:

```yaml
fn_ignore_files:
  - "*.bai"
  - "*.bak"
  - "*.tmp"
```

### general_stats_columns

**Type**: `Dict[str, <class 'multiqc.utils.config_schema.GeneralStatsModuleConfig'>]` (default: `{}`)

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

### general_stats_helptext

**Type**: `Optional[str]` (default: `None`)

Help text shown under the General Statistics heading at the top of the report.

### long_read_count_desc

**Type**: `Optional[str]` (default: `"thousands"`)

Word used in labels for long-read counts, eg. 'thousands'.

**Examples**:

```yaml
long_read_count_desc: millions
```

```yaml
long_read_count_desc: reads
```

### long_read_count_multiplier

**Type**: `Optional[float]` (default: `0.001`)

Multiplier for long-read counts. Default 0.001 shows counts in thousands.

**Examples**:

```yaml
long_read_count_multiplier: 1.0e-06
```

```yaml
long_read_count_multiplier: 1
```

### long_read_count_prefix

**Type**: `Optional[str]` (default: `"K"`)

Suffix shown after formatted long-read counts, eg. 'K' for thousands.

**Examples**:

```yaml
long_read_count_prefix: M
```

```yaml
long_read_count_prefix: ""
```

### module_order

**Type**: `Optional[List[Union[str, Dict[str, Dict[str, Union[str, List[str]]]]]]]` (default: `["custom_content","ccs","ngsderive","purple","conpair","isoseq","lima","peddy","percolator","haplocheck","somalier","methylqa","mosdepth","phantompeakqualtools","qualimap","bamdst","preseq","hifiasm","quast","qorts","rna_seqc","rockhopper","rsem","rseqc","busco","checkm","bustools","goleft_indexcov","gffcompare","disambiguate","supernova","deeptools","sargasso","verifybamid","mirtrace","happy","mirtop","glimpse","gopeaks","homer","hops","macs2","theta2","snpeff","gatk","htseq","bcftools","featurecounts","fgbio","dragen","dragen_fastqc","dedup","pbmarkdup","damageprofiler","mapdamage","biobambam2","jcvi","mtnucratio","picard","vep","bakta","prokka","checkm2","qc3C","nanoq","nanostat","samblaster","samtools","bamtools","sambamba","ngsbits","pairtools","sexdeterrmine","seqera_cli","eigenstratdatabasetools","jellyfish","vcftools","longranger","stacks","varscan2","snippy","umicollapse","umitools","truvari","megahit","ganon","gtdbtk","bbmap","bismark","biscuit","diamond","hicexplorer","hicup","hicpro","salmon","kallisto","slamdunk","star","hisat2","tophat","bowtie2","bowtie1","hostile","cellranger","checkatlas","snpsplit","odgi","vg","pangolin","nextclade","freyja","humid","kat","leehom","librarian","nonpareil","adapterremoval","bbduk","clipandmerge","cutadapt","trim_galore","flexbar","sourmash","kaiju","kraken","malt","motus","trimmomatic","sickle","skewer","sortmerna","biobloomtools","seqfu","fastq_screen","fastqe","afterqc","fastp","fastqc","sequali","filtlong","prinseqplusplus","pychopper","porechop","pycoqc","minionqc","anglerfish","multivcfanalyzer","clusterflow","checkqc","bcl2fastq","bclconvert","interop","ivar","flash","seqyclean","optitype","whatshap","spaceranger","xenome","xengsort","metaphlan","sylphtax","seqwho","telseq","ataqv","mgikit","mosaicatcher"]`)

Order in which modules appear in the report. Each entry is either a module ID, or a single-key dict mapping the ID to per-run overrides (eg. name, path_filters).

**Example**:

```yaml
module_order:
  - fastqc
  - fastqc:
      name: FastQC (trimmed)
      path_filters:
        - "*_trimmed*"
  - cutadapt
```

### no_version_check

**Type**: `Optional[bool]` (default: `false`)

Skip the network check for newer MultiQC versions on startup.

### num_datasets_plot_limit

**Type**: `Optional[int]` (default: `100`)

Deprecated. Use plots_defer_loading_numseries instead.

### pandoc_template

**Type**: `Optional[str]` (default: `None`)

Path to a Pandoc template used when exporting the report as PDF.

### parquet_format

**Type**: `Optional[Literal["long", "wide"]]` (default: `"long"`)

Parquet table layout. 'long' has rows of (sample_name, metric_name, val_raw, val_raw_type, val_str), easy to filter by metric. 'wide' uses one column per metric (prefixed with table name and namespace), easier for analytics but can hit column limits or mixed-type issues.

### plot_font_family

**Type**: `Optional[str]` (default: `None`)

CSS font-family for plot text. Defaults to a system font stack.

### preserve_module_raw_data

**Type**: `Optional[bool]` (default: `false`)

Keep each module's raw parsed data in memory after report generation. Used by Python API consumers.

### read_count_desc

**Type**: `Optional[str]` (default: `"millions"`)

Word used in plot/axis labels for read counts, eg. 'millions'.

**Examples**:

```yaml
read_count_desc: thousands
```

```yaml
read_count_desc: raw reads
```

### read_count_multiplier

**Type**: `Optional[float]` (default: `1e-06`)

Multiplier applied to read counts before display. Default 0.000001 shows reads in millions.

**Examples**:

```yaml
read_count_multiplier: 0.001
```

```yaml
read_count_multiplier: 1
```

### read_count_prefix

**Type**: `Optional[str]` (default: `"M"`)

Suffix shown after formatted read counts, eg. 'M' for millions.

**Examples**:

```yaml
read_count_prefix: K
```

```yaml
read_count_prefix: ""
```

### remove_sections

**Type**: `Optional[List[str]]` (default: `[]`)

Module sections to hide. Use the section anchor as it appears in the URL.

**Example**:

```yaml
remove_sections:
  - fastqc_overrepresented_sequences
  - gatk-compare-overlap
```

### section_comments

**Type**: `Optional[Dict[str, str]]` (default: `{}`)

Markdown text shown under specific module sections. Keys are section anchors.

**Example**:

```yaml
section_comments:
  fastqc_overrepresented_sequences: "**This is** an important note about the overrepresented
    sequences."
  samtools: Reviewed by *Phil* on 2024-08-21.
```

### section_status_checks

**Type**: `Optional[Dict[str, Union[bool, Dict[str, bool]]]]` (default: `{}`)

Enable or disable the green/yellow/red status indicators on report sections. Top-level keys are module IDs, values are either a bool or a dict mapping section ID to bool.

**Example**:

```yaml
section_status_checks:
  fastqc: true
  samtools:
    alignment_stats: false
```

### skip_generalstats

**Type**: `Optional[bool]` (default: `false`)

Hide the General Statistics table at the top of the report.

### skip_versions_section

**Type**: `Optional[bool]` (default: `false`)

Hide the Software Versions section.

### software_versions

**Type**: `Optional[Dict[str, Any]]` (default: `{}`)

Manually specify software versions for the Software Versions section. Top-level keys are tool names.

**Example**:

```yaml
software_versions:
  bwa: 0.7.17
  fastqc: 0.12.1
  samtools: "1.20"
```

### sp

**Type**: `Optional[Dict[str, Union[SearchPattern, List[SearchPattern]]]]` (default: `None`)

Search patterns for finding tool outputs

### table_columns_name

**Type**: `Optional[Dict[str, Union[str, Dict[str, str]]]]` (default: `{}`)

Rename table columns. Top-level keys are module IDs, inner keys are column IDs, values are the new display name.

**Example**:

```yaml
table_columns_name:
  fastqc:
    percent_duplicates: "% Dups"
    percent_gc: "% GC"
```

### table_columns_placement

**Type**: `Optional[Dict[str, Dict[str, float]]]` (default: `{}`)

Reorder table columns. Top-level keys are module IDs, inner keys are column IDs, values are float sort weights (lower is further left).

**Example**:

```yaml
table_columns_placement:
  fastqc:
    percent_duplicates: 900
    percent_gc: 800
    total_sequences: 700
```

### table_columns_visible

**Type**: `Optional[Dict[str, Union[bool, Dict[str, bool]]]]` (default: `{}`)

Hide or show specific columns. Top-level keys are module IDs, values are either a bool (apply to all columns) or a dict mapping column ID to bool.

**Example**:

```yaml
table_columns_visible:
  fastqc: false
  samtools:
    error_rate: false
    raw_total_sequences: true
```

### table_cond_formatting_colours

**Type**: `Optional[List[Dict[str, str]]]` (default: `[{"blue":"#337ab7"},{"lbue":"#5bc0de"},{"pass":"#5cb85c"},{"warn":"#f0ad4e"},{"fail":"#d9534f"},{"male":"#5bc0de"},{"female":"#d9534f"}]`)

Background colours referenced by table_cond_formatting_rules. List of single-key dicts mapping a colour ID to a hex code.

**Example**:

```yaml
table_cond_formatting_colours:
  - pass: "#5cb85c"
  - warn: "#f0ad4e"
  - fail: "#d9534f"
```

### table_cond_formatting_rules

**Type**: `Optional[Dict[str, Dict[str, List[Dict[str, Union[str, int, float]]]]]]` (default: `{"all_columns":{"pass":[{"s_eq":"pass"},{"s_eq":"true"},{"s_eq":"yes"},{"s_eq":"ok"}],"warn":[{"s_eq":"warn"},{"s_eq":"unknown"}],"fail":[{"s_eq":"fail"},{"s_eq":"false"},{"s_eq":"no"}],"male":[{"s_eq":"male"},{"s_eq":"M"}],"female":[{"s_eq":"female"},{"s_eq":"F"}]},"QCStatus":{"fail":[{"s_contains":"fail"}]}}`)

Conditional cell formatting. Nested dicts map table ID to column ID to a list of rules (eg. {s_eq: pass} matches an exact value). See the customisation docs for the full grammar.

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

### table_sample_merge

**Type**: `Optional[Dict[str, List[Union[str, Dict[str, Union[str, List[str]]]]]]]` (default: `None`)

Group samples by merging rows of supporting modules' tables, by collapsing samples that match a pattern. Keys are the merged group name, values are clean-pattern entries (string or {type, pattern}).

**Example**:

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

### template_dark_mode

**Type**: `Optional[bool]` (default: `true`)

Enable the dark mode toggle in the report template.

### top_modules

**Type**: `Optional[List[Union[str, Dict[str, Dict[str, str]]]]]` (default: `[]`)

Module IDs to render before module_order. Useful for pinning a module to the top regardless of where it appears in module_order. Same shape as module_order entries.

**Example**:

```yaml
top_modules:
  - fastqc
  - cutadapt
```

### version_check_url

**Type**: `Optional[str]` (default: `"https://api.multiqc.info/version"`)

URL queried by MultiQC's own update check. Set to override the default endpoint.

### versions_table_group_header

**Type**: `Optional[str]` (default: `"Group"`)

Column header for the grouping column in the Software Versions table. Defaults to 'Group'.

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

- **contents** (`Optional[Union[str, List[str]]]`): File contents to match
- **contents_re** (`Optional[Union[str, List[str]]]`): File contents regex pattern to match
- **exclude_contents** (`Optional[Union[str, List[str]]]`): Exclude files containing this content
- **exclude_contents_re** (`Optional[Union[str, List[str]]]`): Exclude files containing this regex content
- **exclude_fn** (`Optional[Union[str, List[str]]]`): Exclude files matching this pattern
- **exclude_fn_re** (`Optional[Union[str, List[str]]]`): Exclude files matching this regex pattern
- **fn** (`Optional[str]`): Filename pattern to match
- **fn_re** (`Optional[str]`): Filename regex pattern to match
- **max_filesize** (`Optional[int]`): Maximum file size to process
- **num_lines** (`Optional[int]`): Number of lines to search
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

- **module** (`Optional[Union[str, List[str]]]`): Module(s) to apply this pattern to
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

- **ceiling** (`Optional[float]`): Ceiling value
- **description** (`Optional[str]`): Column description
- **floor** (`Optional[float]`): Floor value
- **format** (`Optional[str]`): Number format
- **hidden** (`Optional[bool]`): Whether column is hidden by default
- **max** (`Optional[float]`): Maximum value
- **min** (`Optional[float]`): Minimum value
- **namespace** (`Optional[str]`): Column namespace
- **placement** (`Optional[float]`): Column placement order
- **scale** (`Optional[str]`): Color scale
- **shared_key** (`Optional[str]`): Shared key name
- **title** (`Optional[str]`): Column title
