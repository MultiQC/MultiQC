"""
JSON Schema for MultiQC config validation.
Generated from the config defaults and type hints.
"""

from contextlib import contextmanager
from contextvars import ContextVar
from typing import Any, Dict, Iterator, List, Literal, Optional, Union

from pydantic import BaseModel, ConfigDict, Field


class SearchPattern(BaseModel):
    """Search pattern configuration for finding tool outputs"""

    fn: Optional[str] = Field(None, description="Filename pattern to match")
    fn_re: Optional[str] = Field(None, description="Filename regex pattern to match")
    contents: Optional[Union[str, List[str]]] = Field(None, description="File contents to match")
    contents_re: Optional[Union[str, List[str]]] = Field(None, description="File contents regex pattern to match")
    num_lines: Optional[int] = Field(None, description="Number of lines to search")
    shared: bool = Field(False, description="Allow file to be processed by multiple search patterns")
    skip: bool = Field(False, description="Skip this search pattern")
    max_filesize: Optional[int] = Field(None, description="Maximum file size to process")
    exclude_fn: Optional[Union[str, List[str]]] = Field(None, description="Exclude files matching this pattern")
    exclude_fn_re: Optional[Union[str, List[str]]] = Field(
        None, description="Exclude files matching this regex pattern"
    )
    exclude_contents: Optional[Union[str, List[str]]] = Field(None, description="Exclude files containing this content")
    exclude_contents_re: Optional[Union[str, List[str]]] = Field(
        None, description="Exclude files containing this regex content"
    )


class CleanPattern(BaseModel):
    """Pattern for cleaning sample names"""

    type: Literal["truncate", "remove", "regex", "regex_keep"] = Field(
        "truncate", description="Type of pattern matching to use"
    )
    pattern: str = Field(..., description="Pattern to match")
    module: Optional[Union[str, List[str]]] = Field(None, description="Module(s) to apply this pattern to")


class GeneralStatsColumnConfig(BaseModel):
    """Configuration for a general stats column"""

    title: Optional[str] = Field(None, description="Column title")
    description: Optional[str] = Field(None, description="Column description")
    namespace: Optional[str] = Field(None, description="Column namespace")
    scale: Optional[str] = Field(None, description="Color scale")
    format: Optional[str] = Field(None, description="Number format")
    min: Optional[float] = Field(None, description="Minimum value")
    max: Optional[float] = Field(None, description="Maximum value")
    ceiling: Optional[float] = Field(None, description="Ceiling value")
    floor: Optional[float] = Field(None, description="Floor value")
    shared_key: Optional[str] = Field(None, description="Shared key name")
    hidden: Optional[bool] = Field(None, description="Whether column is hidden by default")
    placement: Optional[float] = Field(None, description="Column placement order")


class GeneralStatsModuleConfig(BaseModel):
    """Configuration for a module's general stats columns"""

    columns: Dict[str, GeneralStatsColumnConfig] = Field(
        default_factory=dict, description="Columns to show in general stats table. Keys are column IDs."
    )


AiProviderLiteral = Literal["seqera", "openai", "anthropic", "aws_bedrock", "custom"]

_current_section: ContextVar[Optional[str]] = ContextVar("_current_section", default=None)


@contextmanager
def section(name: str) -> Iterator[None]:
    """Group fields by section inside the MultiQCConfig class body.

    Used as a ``with`` block; every ``cfg()`` call inside inherits ``section=name``
    automatically. Nested blocks override outer ones. An explicit ``section=``
    on a ``cfg()`` call still wins. Python's ``with`` does not introduce a new
    scope, so annotations inside the block still land in the class's
    ``__annotations__`` in source order. Section display order in the docs and
    wizard follows the order in which each section name first appears below.
    """
    token = _current_section.set(name)
    try:
        yield
    finally:
        _current_section.reset(token)


def cfg(
    description: str,
    *,
    section: Optional[str] = None,
    advanced: bool = False,
    multiline: bool = False,
    default: Any = None,
    default_factory: Any = None,
    **kwargs: Any,
) -> Any:
    """Wrapper around ``Field`` that tags each option with a section.

    ``section`` is usually inherited from an enclosing ``with section(...)``
    block; pass ``section=`` explicitly to override or for a one-off field
    outside any block. ``advanced=True`` hides the field behind the wizard's
    "Show advanced options" toggle. ``multiline=True`` renders string fields
    as a textarea instead of a single-line input, use for prose fields like
    descriptions and prompts. All flags end up in the JSON schema under
    ``json_schema_extra`` and are read back by
    ``scripts/_config_schema_loader.py``. Any other ``Field`` kwargs
    (``examples``, validators, ``deprecated``) pass straight through via
    ``**kwargs`` - Pydantic surfaces ``deprecated`` as ``"deprecated": true``
    in the JSON schema, which the wizard's Validate YAML view picks up.
    """
    if section is None:
        section = _current_section.get()
    if section is None:
        raise ValueError('cfg() needs a section. Wrap the field in `with section("..."):` or pass section= explicitly.')
    extra: Dict[str, Any] = {"section": section}
    if advanced:
        extra["advanced"] = True
    if multiline:
        extra["multiline"] = True
    if default_factory is not None:
        return Field(default_factory=default_factory, description=description, json_schema_extra=extra, **kwargs)
    return Field(default, description=description, json_schema_extra=extra, **kwargs)


class MultiQCConfig(BaseModel):
    """Schema for MultiQC config validation"""

    with section("Report Meta"):
        title: Optional[str] = cfg("Title shown at the top of the report and used in the page title.")
        subtitle: Optional[str] = cfg("Subtitle shown under the report title. Plain text only.")
        intro_text: Optional[str] = cfg(
            "Paragraph shown under the title. Useful for adding context about the analysis.",
            multiline=True,
        )
        report_comment: Optional[str] = cfg(
            "Free-text comment shown at the top of the report. HTML is allowed.",
            multiline=True,
            examples=["This report was generated from the RNA-seq pipeline on 2024-08-21."],
        )
        report_header_info: Optional[List[Dict[str, str]]] = cfg(
            (
                "Extra key/value pairs shown in the report header, eg. contact name, run ID, pipeline version. "
                "Each list item is a single-key dictionary."
            ),
            examples=[
                [
                    {"Contact E-mail": "phil.ewels@seqera.io"},
                    {"Application Type": "RNA-seq"},
                    {"Project Type": "Application"},
                    {"Sequencing Platform": "HiSeq 2500 High Output V4"},
                ]
            ],
        )
        show_analysis_paths: Optional[bool] = cfg(
            "Show the absolute paths of analysed directories in the report header."
        )
        show_analysis_time: Optional[bool] = cfg("Show the date and time the report was generated in the header.")

    with section("Report Appearance"):
        template: Optional[Literal["default", "original", "simple", "sections", "gathered", "geo", "disco"]] = cfg(
            "Name of the report template.",
        )
        template_dark_mode: Optional[bool] = cfg(
            "Enable the dark mode toggle in the report template.",
            advanced=True,
        )
        custom_logo: Optional[str] = cfg(
            "Path to an image to show at the top of the report, replacing the MultiQC logo.",
            examples=["/path/to/logo.png", "./assets/logo.svg"],
        )
        custom_logo_dark: Optional[str] = cfg(
            "Path to an alternative logo for dark mode. Falls back to custom_logo if unset.",
            examples=["./assets/logo_dark.svg"],
        )
        custom_logo_url: Optional[str] = cfg(
            "URL the custom logo links to when clicked.",
            examples=["https://www.scilifelab.se"],
        )
        custom_logo_title: Optional[str] = cfg(
            "Tooltip text shown when hovering over the custom logo.",
            examples=["Our institute name"],
        )
        custom_logo_width: Optional[int] = cfg(
            "Logo width in pixels. Height scales proportionally.",
            examples=[200],
            gt=0,
        )
        custom_favicon: Optional[str] = cfg(
            "Path to a custom favicon image to show in the browser tab.",
            examples=["/path/to/favicon.ico", "./assets/favicon.png"],
        )
        custom_css_files: Optional[List[str]] = cfg(
            "Paths to additional CSS files to inline into the report. Useful for branding overrides.",
            examples=[["./assets/custom.css", "./assets/branding.css"]],
        )
        simple_output: Optional[bool] = cfg(
            "Render a minimal HTML report without the toolbox or interactive widgets. Useful for very large reports."
        )

    with section("Report Contents"):
        custom_content: Optional[Dict[str, Any]] = cfg(
            "Embed arbitrary plots, tables or text in the report. See the "
            "[Custom Content docs](https://docs.seqera.io/multiqc/custom_content) for the full structure.",
            examples=[
                {
                    "order": ["my-section-id", "my-other-section-id"],
                    "data": {
                        "my-section-id": {
                            "id": "my-section-id",
                            "section_name": "My Custom Section",
                            "plot_type": "table",
                            "data": {"sample1": {"col1": 100}, "sample2": {"col1": 200}},
                        }
                    },
                }
            ],
        )
        custom_content_modules: Optional[List[str]] = cfg(
            "Extra module IDs whose output should be parsed as custom content.",
            advanced=True,
        )
        custom_data: Optional[Dict[str, Any]] = cfg(
            "Inline custom content data keyed by section ID. Companion to custom_content for users who prefer "
            "splitting the metadata and the data across two top-level keys.",
        )
        top_modules: Optional[List[Union[str, Dict[str, Dict[str, str]]]]] = cfg(
            (
                "Module IDs to render before module_order. Useful for pinning a module to the top "
                "regardless of where it appears in module_order. Same shape as module_order entries."
            ),
            examples=[["fastqc", "cutadapt"]],
        )
        module_order: Optional[List[Union[str, Dict[str, Dict[str, Union[str, List[str]]]]]]] = cfg(
            (
                "Order in which modules appear in the report. Each entry is either a module ID, "
                "or a single-key dict mapping the ID to per-run overrides (eg. name, path_filters)."
            ),
            examples=[
                [
                    "fastqc",
                    {"fastqc": {"name": "FastQC (trimmed)", "path_filters": ["*_trimmed*"]}},
                    "cutadapt",
                ]
            ],
        )
        run_modules: Optional[List[str]] = cfg(
            "Module IDs to run. If set, only listed modules are processed (mirror of the --module CLI flag).",
            examples=[["fastqc", "cutadapt", "samtools"]],
        )
        exclude_modules: Optional[List[str]] = cfg(
            "Module IDs to skip (mirror of the --exclude CLI flag).",
            examples=[["fastqc"]],
        )
        remove_sections: Optional[List[str]] = cfg(
            "Module sections to hide. Use the section anchor as it appears in the URL.",
            examples=[["fastqc_overrepresented_sequences", "gatk-compare-overlap"]],
        )
        report_section_order: Optional[Dict[str, Any]] = cfg(
            (
                "Reorder, group or hide report sections by ID. Values can be a position string ('before'/'after'), "
                "an explicit order number, or a dict of overrides. See the [customisation docs](https://docs.seqera.io/multiqc/reports/customisation#order-of-module-and-module-subsection-output) for the full grammar."
            ),
            examples=[{"fastqc": {"order": -10}, "custom_content-my-section": {"before": "fastqc"}}],
        )
        section_comments: Optional[Dict[str, str]] = cfg(
            "Markdown text shown under specific module sections. Keys are section anchors.",
            examples=[
                {
                    "fastqc_overrepresented_sequences": "**This is** an important note about the overrepresented sequences.",
                    "samtools": "Reviewed by *Phil* on 2024-08-21.",
                }
            ],
        )
        section_status_checks: Optional[Dict[str, Union[bool, Dict[str, bool]]]] = cfg(
            (
                "Enable or disable the green/yellow/red status indicators on report sections. "
                "Top-level keys are module IDs, values are either a bool or a dict mapping section ID to bool."
            ),
            advanced=True,
            examples=[{"fastqc": True, "samtools": {"alignment_stats": False}}],
        )

    with section("Output Options"):
        force: Optional[bool] = cfg("Overwrite existing output files without prompting.")
        output_fn_name: Optional[str] = cfg("Filename for the generated HTML report. Defaults to multiqc_report.html.")
        data_dir_name: Optional[str] = cfg(
            "Name of the directory written alongside the report holding parsed data. Defaults to multiqc_data."
        )
        plots_dir_name: Optional[str] = cfg(
            "Directory for exported plot images when export_plots is on. Defaults to multiqc_plots.",
        )
        data_format: Optional[Literal["tsv", "csv", "json", "yaml"]] = cfg(
            "Format used when writing parsed data files.",
        )
        data_format_extensions: Optional[Dict[str, str]] = cfg(
            "Override the file extension used when writing each data format, eg. {tsv: txt} to write TSV as .txt.",
            advanced=True,
            examples=[{"tsv": "txt", "json": "json", "yaml": "yml"}],
        )
        parquet_format: Optional[Literal["long", "wide"]] = cfg(
            (
                "Parquet table layout. 'long' has rows of (sample_name, metric_name, val_raw, val_raw_type, val_str), "
                "easy to filter by metric. 'wide' uses one column per metric (prefixed with table name and namespace), "
                "easier for analytics but can hit column limits or mixed-type issues."
            ),
            advanced=True,
        )

        make_data_dir: Optional[bool] = cfg("Write parsed data as files alongside the report.")
        zip_data_dir: Optional[bool] = cfg("Compress the data directory into a single .zip file.")
        data_dump_file: Optional[bool] = cfg(
            "Write a single JSON file containing all parsed data, for re-running MultiQC later.",
            advanced=True,
        )
        data_dump_file_write_raw: Optional[bool] = cfg(
            "Include raw values (before any normalisation or filtering) in the dumped JSON.",
            advanced=True,
        )
        export_plots: Optional[bool] = cfg(
            "Save each plot as a static image (formats set by export_plot_formats).",
        )
        export_plot_formats: Optional[List[Literal["png", "svg", "pdf"]]] = cfg(
            "Image formats to export when export_plots is on.",
        )
        export_plots_timeout: Optional[int] = cfg("Timeout for exporting each plot, in seconds.")
        make_report: Optional[bool] = cfg(
            "Generate the HTML report. Set to false to only produce data files.",
        )
        make_pdf: Optional[bool] = cfg(
            "Also generate a PDF version of the report. Requires Pandoc to be installed.",
        )
        pandoc_template: Optional[str] = cfg(
            "Path to a Pandoc template used when exporting the report as PDF.",
            advanced=True,
        )

    with section("Sample Names"):
        prepend_dirs: Optional[bool] = cfg(
            "Prefix sample names with their parent directory. Useful when the same sample name occurs in multiple directories.",
        )
        prepend_dirs_depth: Optional[int] = cfg(
            "How many parent directories to include. 0 means all the way to the root.",
        )
        prepend_dirs_sep: Optional[str] = cfg(
            "String inserted between directory names and the sample name. Defaults to '|'.",
            advanced=True,
            examples=["_", " - "],
        )
        fn_clean_sample_names: Optional[bool] = cfg(
            "Apply the cleaning rules in fn_clean_exts and fn_clean_trim to sample names.",
        )
        extra_fn_clean_exts: Optional[List[Union[str, CleanPattern]]] = cfg(
            "Extensions appended to the built-in list. Use to add custom suffixes without overriding defaults.",
            examples=[[".mySuffix", {"type": "remove", "pattern": "_tmp", "module": ["samtools"]}]],
        )
        extra_fn_clean_trim: Optional[List[str]] = cfg(
            "Strings appended to the built-in trim list, without overriding defaults.",
            examples=[["sample_", "_processed"]],
        )
        fn_clean_exts: Optional[List[Union[str, CleanPattern]]] = cfg(
            "Extensions stripped from sample names, eg. .gz, .fastq. Replaces the built-in list.",
            examples=[[".gz", ".fastq", ".bam", {"type": "regex", "pattern": r"_S\d+_L\d+"}]],
            advanced=True,
        )
        fn_clean_trim: Optional[List[str]] = cfg(
            "Strings trimmed from the start or end of sample names. Replaces the built-in list.",
            examples=[["_R1", "_R2", "_001"]],
            advanced=True,
        )
        use_filename_as_sample_name: Optional[Union[bool, List[str]]] = cfg(
            (
                "Use the source filename as the sample name instead of any name parsed from the log. "
                "Set to true for all modules, or to a list of module IDs / patterns to apply selectively."
            ),
            default=False,
        )
        sample_names_ignore: Optional[List[str]] = cfg(
            "Glob patterns. Matching samples are dropped from the report.",
            examples=[["*_temp", "control_*"]],
        )
        sample_names_ignore_re: Optional[List[str]] = cfg(
            "Regex patterns. Matching samples are dropped from the report.",
            examples=[[r"^test_.*", r".*_neg_ctrl$"]],
        )
        sample_names_only_include: Optional[List[str]] = cfg(
            "Glob patterns. If set, only matching samples are kept.",
            advanced=True,
            examples=[["RNA_*", "Sample_??"]],
        )
        sample_names_only_include_re: Optional[List[str]] = cfg(
            "Regex patterns. If set, only matching samples are kept.",
            advanced=True,
            examples=[[r"^WGS_[0-9]+$"]],
        )
        sample_names_rename: Optional[List[List[str]]] = cfg(
            "Toolbox rename pairs. Each entry is a [from, to] pair, grouped by the buttons in sample_names_rename_buttons.",
            examples=[
                [
                    ["SMP001", "Patient_A"],
                    ["SMP002", "Patient_B"],
                    ["SMP003", "Patient_C"],
                ]
            ],
        )
        sample_names_rename_buttons: Optional[List[str]] = cfg(
            "Names of the toolbox buttons that switch between the rename groups defined in sample_names_rename.",
            examples=[["Sample ID", "Patient ID", "Lane"]],
        )
        sample_names_replace: Optional[Dict[str, str]] = cfg(
            "Substring replacements applied to every sample name. Keys are matched, values are replacements.",
            examples=[{"_001": "", "Sample_": "S"}],
        )
        sample_names_replace_complete: Optional[bool] = cfg(
            "Replace the entire sample name when the key matches anywhere in it.",
        )
        sample_names_replace_exact: Optional[bool] = cfg(
            "Only replace when the key matches the sample name exactly, not as a substring.",
        )
        sample_names_replace_regex: Optional[bool] = cfg(
            "Treat keys in sample_names_replace as regex patterns.",
        )

    with section("File Discovery"):
        file_list: Optional[bool] = cfg(
            "Treat the input path as a file containing a list of paths to scan, one per line.",
        )
        require_logs: Optional[bool] = cfg(
            (
                "Fail with an error if any module explicitly requested with `--module` has no log files found. "
                "Off by default, so missing inputs are skipped silently."
            ),
        )
        log_filesize_limit: Optional[int] = cfg("Skip log files larger than this many bytes.")
        filesearch_lines_limit: Optional[int] = cfg("Stop reading a log file after this many lines.")
        ignore_symlinks: Optional[bool] = cfg(
            "Skip symlinked files and directories during the file search.",
            advanced=True,
        )
        ignore_images: Optional[bool] = cfg(
            "Skip image files (PNG/JPEG/etc.) to avoid wasting time opening them.",
        )
        fn_ignore_dirs: Optional[List[str]] = cfg(
            "Glob patterns for directory names to skip entirely during the file search.",
            examples=[["work", ".nextflow", "*_logs"]],
        )
        fn_ignore_paths: Optional[List[str]] = cfg(
            "Glob patterns for paths to skip during the file search.",
            examples=[["*/test_data/*", "*/.snakemake/*"]],
        )
        fn_ignore_files: Optional[List[str]] = cfg(
            "Glob patterns for file names to skip during the file search.",
            examples=[["*.bai", "*.bak", "*.tmp"]],
        )
        filesearch_file_shared: Optional[List[str]] = cfg(
            "Module IDs whose log files may be matched by multiple modules during the search.",
        )

    with section("Plot Settings"):
        plots_force_flat: Optional[bool] = cfg(
            "Render plots as static images instead of interactive Plotly. Useful for very large reports.",
        )
        plots_force_interactive: Optional[bool] = cfg(
            "Force interactive plots even when MultiQC would normally fall back to flat images.",
        )
        plots_export_font_scale: Optional[float] = cfg(
            "Multiplier applied to font sizes in exported plot images. Bump up for publication-quality output.",
        )
        plots_flat_numseries: Optional[int] = cfg(
            "If a plot has more than this many series, MultiQC switches it from interactive to flat image.",
        )
        plots_defer_loading_numseries: Optional[int] = cfg(
            "Plots with more than this many series start collapsed. The user clicks a button to render them.",
        )
        num_datasets_plot_limit: Optional[int] = cfg(
            "Deprecated. Use `plots_defer_loading_numseries` instead.",
            advanced=True,
            deprecated="Use `plots_defer_loading_numseries` instead.",
        )
        plot_font_family: Optional[str] = cfg(
            "CSS font-family for plot text. Defaults to a system font stack.",
            advanced=True,
        )
        custom_plot_config: Optional[Dict[str, Any]] = cfg(
            "Override plot config options per plot. Top-level keys are plot IDs, values are option dicts.",
            examples=[
                {
                    "fastqc_per_base_sequence_quality_plot": {
                        "title": "FastQC: Mean Quality Scores (custom)",
                        "yaxis": {"title": "Phred score"},
                    }
                }
            ],
        )
        lineplot_number_of_points_to_hide_markers: Optional[int] = cfg(
            "Hide individual data point markers in line plots once the total point count across samples exceeds this.",
            advanced=True,
        )
        barplot_legend_on_bottom: Optional[bool] = cfg(
            "Place bar plot legends below the plot instead of to the side. Not recommended.",
            advanced=True,
        )
        boxplot_boxpoints: Optional[Literal["outliers", "suspectedoutliers", "all", False]] = cfg(
            "How boxplot data points are drawn. Use false to hide individual points.",
        )
        box_min_threshold_outliers: Optional[int] = cfg(
            "When a boxplot has more samples than this, only outlier points are drawn.",
        )
        box_min_threshold_no_points: Optional[int] = cfg(
            "When a boxplot has more samples than this, no individual points are drawn.",
        )
        violin_downsample_after: Optional[int] = cfg(
            "Start downsampling violin plot data once the sample count exceeds this. Keeps rendering snappy.",
        )
        violin_min_threshold_outliers: Optional[int] = cfg(
            "When a violin plot has more samples than this, only outlier points are drawn.",
        )
        violin_min_threshold_no_points: Optional[int] = cfg(
            "When a violin plot has more samples than this, no individual points are drawn.",
        )

    with section("Toolbox"):
        highlight_patterns: Optional[List[str]] = cfg(
            "Substring (or regex) patterns. Matching samples are highlighted in plots and tables.",
            examples=[["control", "treated"]],
        )
        highlight_colors: Optional[List[str]] = cfg(
            "Hex colour for each entry in highlight_patterns, in the same order.",
            examples=[["#377eb8", "#e41a1c"]],
        )
        highlight_regex: Optional[bool] = cfg(
            "Treat highlight_patterns as regex instead of plain substring.",
        )
        show_hide_buttons: Optional[List[str]] = cfg(
            "Labels for the toolbox show/hide buttons. One per pattern set.",
            examples=[["Tumour samples", "Normal samples"]],
        )
        show_hide_patterns: Optional[List[Union[str, List[str]]]] = cfg(
            "Patterns for each show/hide button. Each entry is a string or list of strings to match against sample names.",
            examples=[[["_T_", "_tumour_"], ["_N_", "_normal_"]]],
        )
        show_hide_mode: Optional[List[str]] = cfg(
            "Action for each show/hide button: 'show' (only show matches) or 'hide' (hide matches).",
            examples=[["show", "show"]],
        )
        show_hide_regex: Optional[List[Union[str, bool]]] = cfg(
            "Whether each pattern set is treated as regex. List of bools aligned with show_hide_buttons.",
            examples=[[False, False]],
        )

    with section("Table Settings"):
        collapse_tables: Optional[bool] = cfg(
            "Collapse module tables by default. Users click to expand.",
            advanced=True,
        )
        max_table_rows: Optional[int] = cfg(
            "Tables larger than this many rows are rendered as a violin plot instead.",
        )
        max_configurable_table_columns: Optional[int] = cfg(
            "Cap on the number of columns the user can toggle in the table-configure toolbox.",
        )
        decimalPoint_format: Optional[str] = cfg(
            "Decimal-point character used in formatted numbers, eg. `.` (default) or `,`.",
            examples=[","],
        )
        thousandsSep_format: Optional[str] = cfg(
            "Thousands separator used in formatted numbers, eg. `,` (default), ` ` (space), or `.`",
            examples=[" ", "'"],
        )
        general_stats_columns: Dict[str, GeneralStatsModuleConfig] = cfg(
            "Per-module overrides for General Stats columns. Top-level keys are module IDs.",
            default_factory=dict,
            examples=[
                {
                    "fastqc": {
                        "columns": {
                            "percent_duplicates": {
                                "title": "% Dups",
                                "min": 0,
                                "max": 100,
                                "scale": "RdYlGn-rev",
                                "format": "{:,.1f}%",
                            }
                        }
                    }
                }
            ],
        )
        general_stats_helptext: Optional[str] = cfg(
            "Help text shown under the General Statistics heading at the top of the report.",
            multiline=True,
        )
        skip_generalstats: Optional[bool] = cfg(
            "Hide the General Statistics table at the top of the report.",
        )
        table_columns_name: Optional[Dict[str, Union[str, Dict[str, str]]]] = cfg(
            (
                "Rename table columns. Top-level keys are module IDs, inner keys are column IDs, "
                "values are the new display name."
            ),
            examples=[{"fastqc": {"percent_duplicates": "% Dups", "percent_gc": "% GC"}}],
        )
        table_columns_placement: Optional[Dict[str, Dict[str, float]]] = cfg(
            (
                "Reorder table columns. Top-level keys are module IDs, inner keys are column IDs, "
                "values are float sort weights (lower is further left)."
            ),
            examples=[{"fastqc": {"percent_duplicates": 900, "percent_gc": 800, "total_sequences": 700}}],
        )
        table_columns_visible: Optional[Dict[str, Union[bool, Dict[str, bool]]]] = cfg(
            (
                "Hide or show specific columns. Top-level keys are module IDs, "
                "values are either a bool (apply to all columns) or a dict mapping column ID to bool."
            ),
            examples=[
                {
                    "fastqc": False,
                    "samtools": {"raw_total_sequences": True, "error_rate": False},
                }
            ],
        )
        table_cond_formatting_rules: Optional[Dict[str, Dict[str, List[Dict[str, Union[str, int, float]]]]]] = cfg(
            (
                "Conditional cell formatting. Nested dicts map table ID to column ID to a list of rules "
                "(eg. {s_eq: pass} matches an exact value). See the customisation docs for the full grammar."
            ),
            examples=[
                {
                    "all_columns": {
                        "pass": [{"s_eq": "pass"}, {"s_eq": "ok"}],
                        "warn": [{"s_eq": "warn"}],
                        "fail": [{"s_eq": "fail"}],
                    },
                    "mqc-generalstats-percent_duplicates": {
                        "fail": [{"gt": 50}],
                        "warn": [{"gt": 20}],
                    },
                }
            ],
        )
        table_cond_formatting_colours: Optional[List[Dict[str, str]]] = cfg(
            (
                "Background colours referenced by table_cond_formatting_rules. "
                "List of single-key dicts mapping a colour ID to a hex code."
            ),
            examples=[
                [
                    {"pass": "#5cb85c"},
                    {"warn": "#f0ad4e"},
                    {"fail": "#d9534f"},
                ]
            ],
        )
        custom_table_header_config: Optional[Dict[str, Any]] = cfg(
            "Override table column config. Same shape as custom_plot_config but for table headers.",
            examples=[
                {
                    "general_stats_table": {
                        "% Dups": {"min": 0, "max": 100, "format": "{:,.1f}%"},
                    }
                }
            ],
        )
        table_sample_merge: Optional[Dict[str, List[Union[str, Dict[str, Union[str, List[str]]]]]]] = cfg(
            (
                "Group samples by merging rows of supporting modules' tables, by collapsing samples that match a pattern. "
                "Keys are the merged group name, values are clean-pattern entries (string or {type, pattern})."
            ),
            advanced=True,
            examples=[
                {
                    "R1": ["_R1", {"type": "regex", "pattern": "[_.-][rR]?1$"}],
                    "R2": ["_R2", {"type": "regex", "pattern": "[_.-][rR]?2$"}],
                }
            ],
        )

    with section("Software Versions"):
        software_versions: Optional[Dict[str, Any]] = cfg(
            "Manually specify software versions for the Software Versions section. Top-level keys are tool names.",
            examples=[{"samtools": "1.20", "bwa": "0.7.17", "fastqc": "0.12.1"}],
        )
        versions_table_group_header: Optional[str] = cfg(
            "Column header for the grouping column in the Software Versions table. Defaults to 'Group'.",
        )
        disable_version_detection: Optional[bool] = cfg(
            "Skip parsing software versions from module log files.",
        )
        skip_versions_section: Optional[bool] = cfg("Hide the Software Versions section.")

    with section("Read & Base Counts"):
        read_count_multiplier: Optional[float] = cfg(
            "Multiplier applied to read counts before display. Default 0.000001 shows reads in millions.",
            examples=[0.001, 1],
        )
        read_count_prefix: Optional[str] = cfg(
            "Suffix shown after formatted read counts, eg. 'M' for millions.",
            examples=["K", ""],
        )
        read_count_desc: Optional[str] = cfg(
            "Word used in plot/axis labels for read counts, eg. 'millions'.",
            examples=["thousands", "raw reads"],
        )
        long_read_count_multiplier: Optional[float] = cfg(
            "Multiplier for long-read counts. Default 0.001 shows counts in thousands.",
            examples=[0.000001, 1],
        )
        long_read_count_prefix: Optional[str] = cfg(
            "Suffix shown after formatted long-read counts, eg. 'K' for thousands.",
            examples=["M", ""],
        )
        long_read_count_desc: Optional[str] = cfg(
            "Word used in labels for long-read counts, eg. 'thousands'.",
            examples=["millions", "reads"],
        )
        base_count_multiplier: Optional[float] = cfg(
            "Multiplier for base counts. Default 0.000000001 shows bases in gigabases.",
            examples=[0.000001, 0.001],
        )
        base_count_prefix: Optional[str] = cfg(
            "Suffix shown after formatted base counts, eg. 'Gb' for gigabases.",
            examples=["Mb", "Kb"],
        )
        base_count_desc: Optional[str] = cfg(
            "Word used in labels for base counts, eg. 'gigabases'.",
            examples=["megabases", "kilobases"],
        )

    with section("AI Summary"):
        ai_summary: Optional[bool] = cfg(
            "Generate a short AI-written summary at the top of the report.",
        )
        ai_summary_full: Optional[bool] = cfg(
            "Also generate a longer per-section AI summary. Requires ai_summary to be on.",
        )
        ai_provider: Optional[AiProviderLiteral] = cfg(
            "AI provider used for summaries. One of seqera, openai, anthropic, aws_bedrock, custom.",
        )
        ai_model: Optional[str] = cfg(
            "Model name. Provider-specific.",
            examples=["gpt-4o", "claude-sonnet-4-5."],
        )
        ai_custom_endpoint: Optional[str] = cfg(
            "Base URL for the 'custom' provider, eg. a self-hosted OpenAI-compatible API.",
            examples=["http://localhost:11434/v1", "https://api.example.com/v1"],
        )
        ai_auth_type: Optional[Literal["bearer", "api-key"]] = cfg(
            (
                "Authentication scheme used by the custom endpoint. "
                "'bearer' sends an Authorization header, 'api-key' sends an api-key header."
            ),
            advanced=True,
        )
        ai_retries: Optional[int] = cfg(
            "Number of times to retry an AI request on transient errors.",
            advanced=True,
        )
        ai_extra_query_options: Optional[Dict[str, Any]] = cfg(
            "Extra request-body fields merged into the AI request payload (provider-specific).",
            examples=[{"temperature": 0.3, "top_p": 0.9}],
            advanced=True,
        )
        ai_custom_context_window: Optional[int] = cfg(
            "Override the model's context window in tokens. Set this if MultiQC's default for your model is wrong.",
            advanced=True,
            gt=0,
        )
        ai_max_completion_tokens: Optional[int] = cfg(
            "Maximum completion tokens for OpenAI reasoning models.",
            advanced=True,
        )
        ai_reasoning_effort: Optional[Literal["low", "medium", "high"]] = cfg(
            "Reasoning effort for OpenAI reasoning models.",
            advanced=True,
        )
        ai_extended_thinking: Optional[bool] = cfg(
            "Enable extended thinking on Anthropic Claude models that support it.",
            advanced=True,
        )
        ai_thinking_budget_tokens: Optional[int] = cfg(
            "Token budget for Anthropic extended thinking when enabled.",
            advanced=True,
        )
        ai_prompt_short: Optional[str] = cfg(
            "Custom prompt prepended to the short AI summary request. Use to steer tone, length, or focus.",
            multiline=True,
            examples=["Write the summary in one short paragraph aimed at a lab head, no jargon."],
        )
        ai_prompt_full: Optional[str] = cfg(
            "Custom prompt prepended to the full-section AI summary request.",
            multiline=True,
            examples=["Use bullet points and call out any sample that looks like an outlier."],
        )
        no_ai: Optional[bool] = cfg(
            "Disable AI summaries entirely. Overrides ai_summary and ai_summary_full.",
        )
        ai_anonymize_samples: Optional[bool] = cfg(
            "Replace sample names with placeholders before sending data to the AI provider.",
        )

    with section("MegaQC Integration"):
        megaqc_url: Optional[str] = cfg(
            "URL of a MegaQC instance to upload report data to after generation.",
            advanced=True,
        )
        megaqc_access_token: Optional[str] = cfg(
            "Auth token for the MegaQC instance.",
            advanced=True,
        )
        megaqc_timeout: Optional[int] = cfg(
            "Upload timeout in seconds when posting to MegaQC.",
            advanced=True,
        )
        megaqc_upload: Optional[bool] = cfg(
            "Upload report data to MegaQC after generation. Requires megaqc_url and megaqc_access_token.",
            advanced=True,
        )

    with section("Seqera Integration"):
        seqera_website: Optional[str] = cfg(
            "Base URL used for Seqera Platform links in the report.",
            advanced=True,
        )
        seqera_api_url: Optional[str] = cfg(
            "Base URL for the Seqera Platform API. Defaults to the public instance.",
            advanced=True,
        )

    with section("Performance & Debugging"):
        profile_runtime: Optional[bool] = cfg(
            "Time each module and include the breakdown in the report.",
        )
        profile_memory: Optional[bool] = cfg(
            "Track peak memory per module. Adds runtime overhead.",
            advanced=True,
        )
        verbose: Optional[bool] = cfg("Print extra debug log messages to the terminal.")
        no_ansi: Optional[bool] = cfg("Disable ANSI colour codes in terminal output.")
        quiet: Optional[bool] = cfg("Suppress non-essential log messages.")
        lint: Optional[bool] = cfg(
            "Deprecated. Run module linting and fail the build on issues. Used in MultiQC's own tests, rarely useful otherwise.",
            deprecated="Use `--lint` command-line flag instead.",
        )
        strict: Optional[bool] = cfg("Treat module warnings as errors. Stricter than lint.", advanced=True)
        development: Optional[bool] = cfg(
            "Enable developer-mode features such as live JS reloading. Internal use.",
            advanced=True,
        )
        report_readerrors: Optional[bool] = cfg(
            "Surface file read errors in the log instead of silently skipping them.",
            advanced=True,
        )
        no_version_check: Optional[bool] = cfg(
            "Skip the network check for newer MultiQC versions on startup.",
        )
        version_check_url: Optional[str] = cfg(
            "URL queried by MultiQC's own update check. Set to override the default endpoint.",
            advanced=True,
        )
        preserve_module_raw_data: Optional[bool] = cfg(
            "Keep each module's raw parsed data in memory after report generation. Used by Python API consumers.",
            advanced=True,
        )

    # Search patterns. Excluded from the wizard via SKIP_PROPERTIES; rendered manually
    # by generate_config_docs.py under "Special Types". No section.
    sp: Optional[Dict[str, Union[SearchPattern, List[SearchPattern]]]] = Field(
        None, description="Search patterns for finding tool outputs"
    )

    model_config = ConfigDict(extra="allow")  # Allow additional fields that aren't in the schema


def config_to_schema() -> Dict[str, Any]:
    """Convert the config schema to a JSON Schema dict"""
    return MultiQCConfig.model_json_schema()
