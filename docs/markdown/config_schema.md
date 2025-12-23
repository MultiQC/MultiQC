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

Custom CSS files to include

### custom_logo

**Type**: `Optional[str]` (default: `None`)

Path to custom logo image

### custom_logo_title

**Type**: `Optional[str]` (default: `None`)

Title for custom logo

### custom_logo_url

**Type**: `Optional[str]` (default: `None`)

URL for custom logo

### intro_text

**Type**: `Optional[str]` (default: `None`)

Report introduction text

### report_comment

**Type**: `Optional[str]` (default: `None`)

Report comment

### report_header_info

**Type**: `Optional[List[Dict[str, str]]]` (default: `None`)

Report header dictionary

### show_analysis_paths

**Type**: `Optional[bool]` (default: `true`)

Show analysis paths in the report

### show_analysis_time

**Type**: `Optional[bool]` (default: `true`)

Show analysis time in the report

### simple_output

**Type**: `Optional[bool]` (default: `false`)

Simple output

### subtitle

**Type**: `Optional[str]` (default: `None`)

Report subtitle

### template

**Type**: `Optional[str]` (default: `"default"`)

Report template to use

### title

**Type**: `Optional[str]` (default: `None`)

Report title

## Output Options

### data_dir_name

**Type**: `Optional[str]` (default: `"multiqc_data"`)

Data directory name

### data_dump_file

**Type**: `Optional[bool]` (default: `true`)

Write data to a file

### data_dump_file_write_raw

**Type**: `Optional[bool]` (default: `true`)

Write raw data to a file

### data_format

**Type**: `Optional[str]` (default: `"tsv"`)

Data format for output files

### export_plots

**Type**: `Optional[bool]` (default: `false`)

Export plots

### export_plots_timeout

**Type**: `Optional[int]` (default: `30`)

Export plots timeout

### force

**Type**: `Optional[bool]` (default: `false`)

Overwrite existing reports

### make_data_dir

**Type**: `Optional[bool]` (default: `true`)

Create data directory

### make_pdf

**Type**: `Optional[bool]` (default: `false`)

Make PDF

### make_report

**Type**: `Optional[bool]` (default: `true`)

Make report

### output_fn_name

**Type**: `Optional[str]` (default: `"multiqc_report.html"`)

Output filename

### plots_dir_name

**Type**: `Optional[str]` (default: `"multiqc_plots"`)

Plots directory name

### zip_data_dir

**Type**: `Optional[bool]` (default: `false`)

Zip data directory

## MegaQC Integration

### megaqc_access_token

**Type**: `Optional[str]` (default: `None`)

MegaQC access token

### megaqc_timeout

**Type**: `Optional[int]` (default: `30`)

MegaQC timeout

### megaqc_url

**Type**: `Optional[str]` (default: `None`)

MegaQC URL to upload to

## AI Summary

### ai_anonymize_samples

**Type**: `Optional[bool]` (default: `false`)

Anonymize samples

### ai_auth_type

**Type**: `Optional[str]` (default: `None`)

AI auth type

### ai_custom_context_window

**Type**: `Optional[str]` (default: `None`)

AI custom context window

### ai_custom_endpoint

**Type**: `Optional[str]` (default: `None`)

AI custom endpoint

### ai_extra_query_options

**Type**: `Optional[str]` (default: `None`)

AI extra query options

### ai_model

**Type**: `Optional[str]` (default: `None`)

AI model

### ai_prompt_full

**Type**: `Optional[str]` (default: `None`)

Prompt for full AI summary, put before the report details when sent to the provider

### ai_prompt_short

**Type**: `Optional[str]` (default: `None`)

Prompt for short AI summary, put before the report details when sent to the provider

### ai_provider

**Type**: `Optional[Literal["seqera", "openai", "anthropic", "custom"]]` (default: `"seqera"`)

AI provider

### ai_retries

**Type**: `Optional[int]` (default: `3`)

AI retries

### ai_summary

**Type**: `Optional[bool]` (default: `false`)

AI summary

### ai_summary_full

**Type**: `Optional[bool]` (default: `false`)

AI summary full

### no_ai

**Type**: `Optional[bool]` (default: `false`)

Disable AI

## Seqera Integration

### seqera_api_url

**Type**: `Optional[str]` (default: `"https://intern.seqera.io"`)

Seqera API URL

### seqera_website

**Type**: `Optional[str]` (default: `"https://seqera.io"`)

Seqera website

## Plot Settings

### barplot_legend_on_bottom

**Type**: `Optional[bool]` (default: `false`)

Place bar plot legend at the bottom (not recommended)

### lineplot_number_of_points_to_hide_markers

**Type**: `Optional[int]` (default: `50`)

Number of points to hide markers - sum of data points in all samples

### plots_defer_loading_numseries

**Type**: `Optional[int]` (default: `100`)

Number of series to defer loading - user will need to press button to render plot

### plots_export_font_scale

**Type**: `Optional[float]` (default: `1.0`)

Font scale for exported plots

### plots_flat_numseries

**Type**: `Optional[int]` (default: `2000`)

Number of series to show in flat plots

### plots_force_flat

**Type**: `Optional[bool]` (default: `false`)

Force static plot images

### plots_force_interactive

**Type**: `Optional[bool]` (default: `false`)

Force interactive plots

### violin_downsample_after

**Type**: `Optional[int]` (default: `2000`)

Downsample data for violin plot starting from this number os samples

### violin_min_threshold_no_points

**Type**: `Optional[int]` (default: `1000`)

For more than this number of samples, show no points

### violin_min_threshold_outliers

**Type**: `Optional[int]` (default: `100`)

For more than this number of samples, show only outliers

## Table Settings

### collapse_tables

**Type**: `Optional[bool]` (default: `true`)

Collapse tables

### decimalPoint_format

**Type**: `Optional[str]` (default: `None`)

Decimal point format

### max_configurable_table_columns

**Type**: `Optional[int]` (default: `200`)

Maximum number of columns to show in tables

### max_table_rows

**Type**: `Optional[int]` (default: `500`)

Maximum number of rows to show in tables

### thousandsSep_format

**Type**: `Optional[str]` (default: `None`)

Thousands separator format

## Sample Names

### extra_fn_clean_exts

**Type**: `Optional[List[Union[str, CleanPattern]]]` (default: `None`)

Additional extensions to clean from sample names

### extra_fn_clean_trim

**Type**: `Optional[List[str]]` (default: `None`)

Additional strings to trim from start/end of sample names

### fn_clean_exts

**Type**: `Optional[List[Union[str, CleanPattern]]]` (default: `[".gz",".fastq",".fq",".bam",".cram",".sam",".sra",".vcf",".dat","_tophat",".pbmarkdup.log",".log",".stderr",".out",".spp",".fa",".fasta",".png",".jpg",".jpeg",".html","Log.final","ReadsPerGene",".flagstat","_star_aligned","_fastqc",".hicup",".counts","_counts",".txt",".tsv",".csv",".aligned","Aligned",".merge",".deduplicated",".dedup",".clean",".sorted",".report","| stdin",".geneBodyCoverage",".inner_distance_freq",".junctionSaturation_plot.r",".pos.DupRate.xls",".GC.xls","_slamdunk","_bismark",".conpair",".concordance",".contamination",".BEST.results","_peaks.xls",".relatedness",".cnt",".aqhist",".bhist",".bincov",".bqhist",".covhist",".covstats",".ehist",".gchist",".idhist",".ihist",".indelhist",".lhist",".mhist",".qahist",".qchist",".qhist",".rpkm",".selfSM",".extendedFrags","_SummaryStatistics",".purple.purity",".purple.qc",".trim",".bowtie2",".mkD",".highfreq",".lowfreq",".consensus",".snpEff",".snpeff",".scaffolds",".contigs",".kraken2",".ccurve",".hisat2","_duprate",".markdup",".read_distribution",".junction_annotation",".infer_experiment",".biotype",".ivar",".mpileup",".primer_trim",".mapped",".vep","_vep","ccs","_NanoStats",".cutadapt",".qcML",".mosdepth","_gopeaks",".readCounts",".wgs_contig_mean_cov","_overall_mean_cov","_coverage_metrics",".wgs_fine_hist",".wgs_coverage_metrics",".wgs_hist",".vc_metrics",".gvcf_metrics",".ploidy_estimation_metrics","_overall_mean_cov",".fragment_length_hist",".mapping_metrics",".gc_metrics",".trimmer_metrics",".time_metrics",".quant_metrics",".quant.metrics",".quant.transcript_coverage",".scRNA_metrics",".scRNA.metrics",".scATAC_metrics",".scATAC.metrics",".fastqc_metrics",".labels",".bammetrics.metrics",".filter_summary",".cluster_report",".error.spl",".error.grp",".vgstats"]`)

Extensions to clean from sample names

### fn_clean_sample_names

**Type**: `Optional[bool]` (default: `true`)

Clean sample names

### fn_clean_trim

**Type**: `Optional[List[str]]` (default: `[".",":","_","-",".r","_val",".idxstats","_trimmed",".trimmed",".csv",".yaml",".yml",".json","_mqc","short_summary_","_summary",".summary",".align",".h5","_matrix",".stats",".hist",".phased",".tar","runs_"]`)

Strings to trim from start/end of sample names

### prepend_dirs

**Type**: `Optional[bool]` (default: `false`)

Prepend directories to sample names

### prepend_dirs_depth

**Type**: `Optional[int]` (default: `0`)

Depth to prepend directories

### prepend_dirs_sep

**Type**: `Optional[str]` (default: `" | "`)

Separator for prepended directories

### sample_names_ignore

**Type**: `Optional[List[str]]` (default: `[]`)

Sample names to ignore

### sample_names_ignore_re

**Type**: `Optional[List[str]]` (default: `[]`)

Sample names to ignore (regex)

### sample_names_only_include

**Type**: `Optional[List[str]]` (default: `[]`)

Sample names to include

### sample_names_only_include_re

**Type**: `Optional[List[str]]` (default: `[]`)

Sample names to include (regex)

### sample_names_rename

**Type**: `Optional[List[List[str]]]` (default: `[]`)

Sample names to rename

### sample_names_rename_buttons

**Type**: `Optional[List[str]]` (default: `[]`)

Sample names to rename

### sample_names_replace

**Type**: `Optional[Dict[str, str]]` (default: `{}`)

Sample names to replace

### sample_names_replace_complete

**Type**: `Optional[bool]` (default: `false`)

Sample names to replace (complete)

### sample_names_replace_exact

**Type**: `Optional[bool]` (default: `false`)

Sample names to replace (exact)

### sample_names_replace_regex

**Type**: `Optional[bool]` (default: `false`)

Sample names to replace (regex)

### use_filename_as_sample_name

**Type**: `Optional[bool]` (default: `false`)

Use filename as sample name

## Toolbox

### highlight_colors

**Type**: `Optional[List[str]]` (default: `[]`)

Colors to use for highlighting patterns

### highlight_patterns

**Type**: `Optional[List[str]]` (default: `[]`)

Patterns for highlighting samples

### highlight_regex

**Type**: `Optional[bool]` (default: `false`)

Whether to use regex mode for highlighting

### show_hide_buttons

**Type**: `Optional[List[str]]` (default: `[]`)

Show/hide buttons

### show_hide_mode

**Type**: `Optional[List[str]]` (default: `[]`)

Show/hide mode

### show_hide_patterns

**Type**: `Optional[List[str]]` (default: `[]`)

Show/hide patterns

### show_hide_regex

**Type**: `Optional[List[str]]` (default: `[]`)

Show/hide regex

## Performance & Debugging

### development

**Type**: `Optional[bool]` (default: `false`)

Development

### filesearch_lines_limit

**Type**: `Optional[int]` (default: `1000`)

Filesearch lines limit

### lint

**Type**: `Optional[bool]` (default: `false`)

Lint

### log_filesize_limit

**Type**: `Optional[int]` (default: `50000000`)

Log filesize limit

### no_ansi

**Type**: `Optional[bool]` (default: `false`)

Disable ANSI output

### profile_memory

**Type**: `Optional[bool]` (default: `false`)

Profile memory

### profile_runtime

**Type**: `Optional[bool]` (default: `false`)

Profile runtime

### quiet

**Type**: `Optional[bool]` (default: `false`)

Quiet output

### report_readerrors

**Type**: `Optional[int]` (default: `false`)

Report read errors

### strict

**Type**: `Optional[bool]` (default: `false`)

Strict

### verbose

**Type**: `Optional[bool]` (default: `false`)

Verbose output

## File Discovery

### filesearch_file_shared

**Type**: `Optional[List[str]]` (default: `[]`)

Filesearch file shared

### fn_ignore_dirs

**Type**: `Optional[List[str]]` (default: `["multiqc_data",".git","icarus_viewers","runs_per_reference","not_aligned","contigs_reports"]`)

Directories to ignore

### fn_ignore_paths

**Type**: `Optional[List[str]]` (default: `["*/work/??/??????????????????????????????","*/.snakemake","*/.singularity","*/__pycache__","*/site-packages/multiqc"]`)

Paths to ignore

### ignore_images

**Type**: `Optional[bool]` (default: `true`)

Ignore images

### ignore_symlinks

**Type**: `Optional[bool]` (default: `false`)

Ignore symlinks

### require_logs

**Type**: `Optional[bool]` (default: `false`)

Require logs for reports

## Other

### base_count_desc

**Type**: `Optional[str]` (default: `"millions"`)

Base count description

### base_count_multiplier

**Type**: `Optional[float]` (default: `1e-06`)

Base count multiplier

### base_count_prefix

**Type**: `Optional[str]` (default: `"Mb"`)

Base count prefix

### custom_content

**Type**: `Optional[Dict[str, Any]]` (default: `{"order":[]}`)

Custom content settings

### custom_plot_config

**Type**: `Optional[Dict[str, Any]]` (default: `{}`)

Custom plot config

### custom_table_header_config

**Type**: `Optional[Dict[str, Any]]` (default: `{}`)

Custom table header config

### data_format_extensions

**Type**: `Optional[Dict[str, str]]` (default: `{"tsv":"txt","csv":"csv","json":"json","yaml":"yaml"}`)

Data format extensions

### disable_version_detection

**Type**: `Optional[bool]` (default: `false`)

Disable version detection

### export_plot_formats

**Type**: `Optional[List[str]]` (default: `["png","svg","pdf"]`)

Export plot formats

### file_list

**Type**: `Optional[bool]` (default: `false`)

Create a file list

### general_stats_columns

**Type**: `Dict[str, <class 'multiqc.utils.config_schema.GeneralStatsModuleConfig'>]` (default: `{}`)

Configuration for general stats columns per module. Keys are module IDs.

### long_read_count_desc

**Type**: `Optional[str]` (default: `"thousands"`)

Long read count description

### long_read_count_multiplier

**Type**: `Optional[float]` (default: `0.001`)

Long read count multiplier

### long_read_count_prefix

**Type**: `Optional[str]` (default: `"K"`)

Long read count prefix

### no_version_check

**Type**: `Optional[bool]` (default: `false`)

No version check

### pandoc_template

**Type**: `Optional[str]` (default: `None`)

Pandoc template

### read_count_desc

**Type**: `Optional[str]` (default: `"millions"`)

Read count description

### read_count_multiplier

**Type**: `Optional[float]` (default: `1e-06`)

Read count multiplier

### read_count_prefix

**Type**: `Optional[str]` (default: `"M"`)

Read count prefix

### remove_sections

**Type**: `Optional[List[str]]` (default: `[]`)

Sections to remove

### section_comments

**Type**: `Optional[Dict[str, str]]` (default: `{}`)

Comments for sections

### skip_generalstats

**Type**: `Optional[bool]` (default: `false`)

Skip generalstats

### skip_versions_section

**Type**: `Optional[bool]` (default: `false`)

Skip versions section

### software_versions

**Type**: `Optional[Dict[str, Any]]` (default: `{}`)

Software versions

### sp

**Type**: `Optional[Dict[str, Union[SearchPattern, List[SearchPattern]]]]` (default: `None`)

Search patterns for finding tool outputs

### table_columns_name

**Type**: `Optional[Dict[str, Union[str, Dict[str, str]]]]` (default: `{}`)

Name of columns in tables

### table_columns_placement

**Type**: `Optional[Dict[str, Dict[str, float]]]` (default: `{}`)

Placement of columns in tables

### table_columns_visible

**Type**: `Optional[Dict[str, Union[bool, Dict[str, bool]]]]` (default: `{}`)

Which columns to show in tables

### table_cond_formatting_colours

**Type**: `Optional[List[Dict[str, str]]]` (default: `[{"blue":"#337ab7"},{"lbue":"#5bc0de"},{"pass":"#5cb85c"},{"warn":"#f0ad4e"},{"fail":"#d9534f"},{"male":"#5bc0de"},{"female":"#d9534f"}]`)

Colours to use for conditional formatting in tables

### table_cond_formatting_rules

**Type**: `Optional[Dict[str, Dict[str, List[Dict[str, Union[str, int, float]]]]]]` (default: `{"all_columns":{"pass":[{"s_eq":"pass"},{"s_eq":"true"},{"s_eq":"yes"},{"s_eq":"ok"}],"warn":[{"s_eq":"warn"},{"s_eq":"unknown"}],"fail":[{"s_eq":"fail"},{"s_eq":"false"},{"s_eq":"no"}],"male":[{"s_eq":"male"},{"s_eq":"M"}],"female":[{"s_eq":"female"},{"s_eq":"F"}]},"QCStatus":{"fail":[{"s_contains":"fail"}]}}`)

Rules for conditional formatting in tables

### version_check_url

**Type**: `Optional[str]` (default: `"https://api.multiqc.info/version"`)

Version check URL

### versions_table_group_header

**Type**: `Optional[str]` (default: `"Group"`)

Versions table group header

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
