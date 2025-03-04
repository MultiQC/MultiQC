"""
JSON Schema for MultiQC config validation.
Generated from the config defaults and type hints.
"""

from typing import Any, Dict, List, Literal, Optional, Union

from pydantic import BaseModel, Field

AiProviderLiteral = Literal["seqera", "openai", "anthropic", "custom"]


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


class MultiQCConfig(BaseModel):
    """Schema for MultiQC config validation"""

    title: Optional[str] = Field(None, description="Report title")
    subtitle: Optional[str] = Field(None, description="Report subtitle")
    intro_text: Optional[str] = Field(None, description="Report introduction text")
    report_comment: Optional[str] = Field(None, description="Report comment")
    report_header_info: Optional[List[Dict[str, str]]] = Field(None, description="Report header dictionary")
    show_analysis_paths: Optional[bool] = Field(None, description="Show analysis paths in the report")
    show_analysis_time: Optional[bool] = Field(None, description="Show analysis time in the report")
    custom_logo: Optional[str] = Field(None, description="Path to custom logo image")
    custom_logo_url: Optional[str] = Field(None, description="URL for custom logo")
    custom_logo_title: Optional[str] = Field(None, description="Title for custom logo")
    custom_css_files: Optional[List[str]] = Field(None, description="Custom CSS files to include")
    simple_output: Optional[bool] = Field(None, description="Simple output")
    template: Optional[str] = Field(None, description="Report template to use")
    profile_runtime: Optional[bool] = Field(None, description="Profile runtime")
    profile_memory: Optional[bool] = Field(None, description="Profile memory")
    pandoc_template: Optional[str] = Field(None, description="Pandoc template")
    read_count_multiplier: Optional[float] = Field(None, description="Read count multiplier")
    read_count_prefix: Optional[str] = Field(None, description="Read count prefix")
    read_count_desc: Optional[str] = Field(None, description="Read count description")
    long_read_count_multiplier: Optional[float] = Field(None, description="Long read count multiplier")
    long_read_count_prefix: Optional[str] = Field(None, description="Long read count prefix")
    long_read_count_desc: Optional[str] = Field(None, description="Long read count description")
    base_count_multiplier: Optional[float] = Field(None, description="Base count multiplier")
    base_count_prefix: Optional[str] = Field(None, description="Base count prefix")
    base_count_desc: Optional[str] = Field(None, description="Base count description")
    output_fn_name: Optional[str] = Field(None, description="Output filename")
    data_dir_name: Optional[str] = Field(None, description="Data directory name")
    plots_dir_name: Optional[str] = Field(None, description="Plots directory name")
    data_format: Optional[str] = Field(None, description="Data format for output files")
    force: Optional[bool] = Field(None, description="Overwrite existing reports")
    verbose: Optional[bool] = Field(None, description="Verbose output")
    no_ansi: Optional[bool] = Field(None, description="Disable ANSI output")
    quiet: Optional[bool] = Field(None, description="Quiet output")
    prepend_dirs: Optional[bool] = Field(None, description="Prepend directories to sample names")
    prepend_dirs_depth: Optional[int] = Field(None, description="Depth to prepend directories")
    prepend_dirs_sep: Optional[str] = Field(None, description="Separator for prepended directories")
    file_list: Optional[bool] = Field(None, description="Create a file list")
    require_logs: Optional[bool] = Field(None, description="Require logs for reports")
    version_check_url: Optional[str] = Field(None, description="Version check URL")

    make_data_dir: Optional[bool] = Field(None, description="Create data directory")
    zip_data_dir: Optional[bool] = Field(None, description="Zip data directory")
    data_dump_file: Optional[bool] = Field(None, description="Write data to a file")
    data_dump_file_write_raw: Optional[bool] = Field(None, description="Write raw data to a file")
    megaqc_url: Optional[bool] = Field(None, description="Upload to MegaQC")
    megaqc_access_token: Optional[str] = Field(None, description="MegaQC access token")
    megaqc_timeout: Optional[int] = Field(None, description="MegaQC timeout")
    export_plots: Optional[bool] = Field(None, description="Export plots")
    export_plots_timeout: Optional[int] = Field(None, description="Export plots timeout")
    make_report: Optional[bool] = Field(None, description="Make report")
    make_pdf: Optional[bool] = Field(None, description="Make PDF")

    ai_summary: Optional[bool] = Field(None, description="AI summary")
    ai_summary_full: Optional[bool] = Field(None, description="AI summary full")
    ai_provider: Optional[AiProviderLiteral] = Field(None, description="AI provider")
    ai_model: Optional[str] = Field(None, description="AI model")
    ai_custom_endpoint: Optional[str] = Field(None, description="AI custom endpoint")
    ai_extra_query_options: Optional[str] = Field(None, description="AI extra query options")
    ai_custom_context_window: Optional[str] = Field(None, description="AI custom context window")
    no_ai: Optional[bool] = Field(None, description="Disable AI")
    ai_anonymize_samples: Optional[bool] = Field(None, description="Anonymize samples")

    seqera_api_url: Optional[str] = Field(None, description="Seqera API URL")
    seqera_website: Optional[str] = Field(None, description="Seqera website")

    plots_force_flat: Optional[bool] = Field(None, description="Force static plot images")
    plots_force_interactive: Optional[bool] = Field(None, description="Force interactive plots")
    plots_export_font_scale: Optional[float] = Field(None, description="Font scale for exported plots")
    plots_flat_numseries: Optional[int] = Field(None, description="Number of series to show in flat plots")
    plots_defer_loading_numseries: Optional[int] = Field(
        None, description="Number of series to defer loading - user will need to press button to render plot"
    )
    lineplot_number_of_points_to_hide_markers: Optional[int] = Field(
        None, description="Number of points to hide markers - sum of data points in all samples"
    )
    barplot_legend_on_bottom: Optional[bool] = Field(
        None, description="Place bar plot legend at the bottom (not recommended)"
    )
    violin_downsample_after: Optional[int] = Field(
        None, description="Downsample data for violin plot starting from this number os samples"
    )
    violin_min_threshold_outliers: Optional[int] = Field(
        None, description="For more than this number of samples, show only outliers"
    )
    violin_min_threshold_no_points: Optional[int] = Field(
        None, description="For more than this number of samples, show no points"
    )

    collapse_tables: Optional[bool] = Field(None, description="Collapse tables")
    max_table_rows: Optional[int] = Field(None, description="Maximum number of rows to show in tables")
    max_configurable_table_columns: Optional[int] = Field(
        None, description="Maximum number of columns to show in tables"
    )
    table_columns_visible: Optional[Dict[str, Union[bool, Dict[str, bool]]]] = Field(
        None, description="Which columns to show in tables"
    )
    table_columns_placement: Optional[Dict[str, Dict[str, float]]] = Field(
        None, description="Placement of columns in tables"
    )
    table_columns_name: Optional[Dict[str, Union[str, Dict[str, str]]]] = Field(
        None, description="Name of columns in tables"
    )
    table_cond_formatting_colours: Optional[List[Dict[str, str]]] = Field(
        None, description="Colours to use for conditional formatting in tables"
    )
    table_cond_formatting_rules: Optional[Dict[str, Dict[str, List[Dict[str, Union[str, int, float]]]]]] = Field(
        None, description="Rules for conditional formatting in tables"
    )
    decimalPoint_format: Optional[str] = Field(None, description="Decimal point format")
    thousandsSep_format: Optional[str] = Field(None, description="Thousands separator format")
    remove_sections: Optional[List[str]] = Field(None, description="Sections to remove")
    section_comments: Optional[Dict[str, str]] = Field(None, description="Comments for sections")
    lint: Optional[bool] = Field(None, description="Lint")
    strict: Optional[bool] = Field(None, description="Strict")
    development: Optional[bool] = Field(None, description="Development")
    custom_plot_config: Optional[Dict[str, Any]] = Field(None, description="Custom plot config")
    custom_table_header_config: Optional[Dict[str, Any]] = Field(None, description="Custom table header config")
    software_versions: Optional[Dict[str, Any]] = Field(None, description="Software versions")
    ignore_symlinks: Optional[bool] = Field(None, description="Ignore symlinks")
    ignore_images: Optional[bool] = Field(None, description="Ignore images")
    fn_ignore_dirs: Optional[List[str]] = Field(None, description="Directories to ignore")
    fn_ignore_paths: Optional[List[str]] = Field(None, description="Paths to ignore")
    sample_names_ignore: Optional[List[str]] = Field(None, description="Sample names to ignore")
    sample_names_ignore_re: Optional[List[str]] = Field(None, description="Sample names to ignore (regex)")
    sample_names_rename_buttons: Optional[List[str]] = Field(None, description="Sample names to rename")
    sample_names_replace: Optional[Dict[str, str]] = Field(None, description="Sample names to replace")
    sample_names_replace_regex: Optional[bool] = Field(None, description="Sample names to replace (regex)")
    sample_names_replace_exact: Optional[bool] = Field(None, description="Sample names to replace (exact)")
    sample_names_replace_complete: Optional[bool] = Field(None, description="Sample names to replace (complete)")
    sample_names_rename: Optional[List[List[str]]] = Field(None, description="Sample names to rename")
    show_hide_buttons: Optional[List[str]] = Field(None, description="Show/hide buttons")
    show_hide_patterns: Optional[List[str]] = Field(None, description="Show/hide patterns")
    show_hide_regex: Optional[List[str]] = Field(None, description="Show/hide regex")
    show_hide_mode: Optional[List[str]] = Field(None, description="Show/hide mode")
    highlight_patterns: Optional[List[str]] = Field(None, description="Patterns for highlighting samples")
    highlight_colors: Optional[List[str]] = Field(None, description="Colors to use for highlighting patterns")
    highlight_regex: Optional[bool] = Field(None, description="Whether to use regex mode for highlighting")
    no_version_check: Optional[bool] = Field(None, description="No version check")
    log_filesize_limit: Optional[int] = Field(None, description="Log filesize limit")
    filesearch_lines_limit: Optional[int] = Field(None, description="Filesearch lines limit")
    report_readerrors: Optional[int] = Field(None, description="Report read errors")
    skip_generalstats: Optional[bool] = Field(None, description="Skip generalstats")
    skip_versions_section: Optional[bool] = Field(None, description="Skip versions section")
    disable_version_detection: Optional[bool] = Field(None, description="Disable version detection")
    versions_table_group_header: Optional[str] = Field(None, description="Versions table group header")
    data_format_extensions: Optional[Dict[str, str]] = Field(None, description="Data format extensions")
    export_plot_formats: Optional[List[str]] = Field(None, description="Export plot formats")
    filesearch_file_shared: Optional[List[str]] = Field(None, description="Filesearch file shared")
    custom_content: Optional[Dict[str, Any]] = Field(None, description="Custom content")
    fn_clean_sample_names: Optional[bool] = Field(None, description="Clean sample names")
    use_filename_as_sample_name: Optional[bool] = Field(None, description="Use filename as sample name")
    fn_clean_exts: Optional[List[Union[str, CleanPattern]]] = Field(
        None, description="Extensions to clean from sample names"
    )
    fn_clean_trim: Optional[List[str]] = Field(None, description="Strings to trim from start/end of sample names")
    extra_fn_clean_exts: Optional[List[Union[str, CleanPattern]]] = Field(
        None, description="Additional extensions to clean from sample names"
    )
    extra_fn_clean_trim: Optional[List[str]] = Field(
        None, description="Additional strings to trim from start/end of sample names"
    )

    # Search patterns
    sp: Optional[Dict[str, SearchPattern]] = Field(None, description="Search patterns for finding tool outputs")

    class Config:
        extra = "allow"  # Allow additional fields that aren't in the schema


def config_to_schema() -> Dict[str, Any]:
    """Convert the config schema to a JSON Schema dict"""
    return MultiQCConfig.model_json_schema()
