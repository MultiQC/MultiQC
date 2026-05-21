"""
Module for handling plot data storage in a single parquet file.
This replaces the individual plot parquet files with a single file that contains data from all plots.
It also stores all report metadata (modules, data sources, configs) to make reports fully reproducible.
"""

import json
import logging
import os
from re import Pattern
from typing import Any, Dict, List, Optional, Set

import polars as pl
from pydantic import ValidationError  # type: ignore

from multiqc import config, report
from multiqc.core import tmp_dir
from multiqc.types import Anchor, ColumnKey
from multiqc.utils.config_schema import MultiQCConfig

logger = logging.getLogger(__name__)

# Set to keep track of which anchors have been saved
_saved_anchors: Set[Anchor] = set()
# Keep track of metric column names
_metric_col_names: Set[ColumnKey] = set()
# Buffer for batched parquet writes to avoid O(n²) read-concat-write behavior
_pending_dataframes: List[pl.DataFrame] = []
# Buffer for wide table dataframes (need special merging by sample)
_pending_wide_tables: List[pl.DataFrame] = []


def wide_table_to_parquet(table_df: pl.DataFrame, metric_col_names: Set[ColumnKey]) -> None:
    """
    Buffer wide-format table data for later merging and writing.

    This function buffers table data instead of immediately writing to avoid
    O(n²) read-concat-write behavior. The actual merging by sample happens
    in flush_to_parquet().
    """
    global _pending_wide_tables
    table_df = fix_creation_date(table_df)
    _pending_wide_tables.append(table_df)


def fix_creation_date(df: pl.DataFrame) -> pl.DataFrame:
    """
    Fix for Iceberg. Iceberg never keeps an arbitrary zone offset in the data –
    a value that has a zone is normalised to UTC, and the zone itself is discarded.
    """
    return df.with_columns(pl.col("creation_date").dt.replace_time_zone(None))


def append_to_parquet(df: pl.DataFrame) -> None:
    """
    Buffer plot data for later writing to the parquet file.

    This function buffers data instead of immediately writing to avoid
    O(n²) read-concat-write behavior when many plots are saved.
    Call flush_to_parquet() to write all buffered data at once.
    """
    global _pending_dataframes
    df = fix_creation_date(df)
    _pending_dataframes.append(df)


def flush_to_parquet() -> None:
    """
    Write all buffered dataframes to the parquet file at once.

    This should be called at the end of report generation to efficiently
    write all accumulated plot data in a single operation.
    """
    global _pending_dataframes, _pending_wide_tables

    if not _pending_dataframes and not _pending_wide_tables:
        return

    # Read existing data from file (if any)
    existing_df = _read_or_create_df()

    # Start with existing non-table rows
    existing_other_rows = existing_df.filter(pl.col("type") != "table_row") if not existing_df.is_empty() else None

    # Process wide tables - merge all buffered wide tables by sample
    merged_wide_tables: Optional[pl.DataFrame] = None
    if _pending_wide_tables:
        # Get existing table rows
        existing_table_rows = existing_df.filter(pl.col("type") == "table_row") if not existing_df.is_empty() else None

        # Start with existing table rows or first pending table
        if existing_table_rows is not None and not existing_table_rows.is_empty():
            merged_wide_tables = existing_table_rows
        else:
            merged_wide_tables = None

        # Merge all pending wide tables
        for table_df in _pending_wide_tables:
            if merged_wide_tables is None:
                merged_wide_tables = table_df
            elif table_df.height > 0:
                # Merge by joining on sample and creation_date
                merged_wide_tables = merged_wide_tables.join(table_df, on=["sample", "creation_date"], how="outer")
                # Ensure all columns are present
                all_cols = merged_wide_tables.columns
                for col in table_df.columns:
                    if col not in all_cols:
                        all_cols.append(col)
                merged_wide_tables = merged_wide_tables.select([c for c in all_cols if c in merged_wide_tables.columns])

    # Build list of dataframes to concatenate
    all_dfs: List[pl.DataFrame] = []

    if existing_other_rows is not None and not existing_other_rows.is_empty():
        all_dfs.append(existing_other_rows)

    all_dfs.extend(_pending_dataframes)

    if merged_wide_tables is not None and not merged_wide_tables.is_empty():
        all_dfs.append(merged_wide_tables)

    # Write combined data
    if all_dfs:
        combined_df = pl.concat(all_dfs, how="diagonal")
        _write_parquet(combined_df)

    # Clear buffers
    _pending_dataframes = []
    _pending_wide_tables = []


def get_report_metadata(df: pl.DataFrame) -> Optional[Dict[str, Any]]:
    """
    Extract all report metadata from the parquet file.

    Args:
        file_path: Optional path to the parquet file. If None, uses the default path.

    Returns a dictionary with modules, data_sources, creation_date, and config.
    """
    if df.is_empty():
        return None

    try:
        # Read the metadata table from the parquet file
        metadata_df = df.filter(pl.col("anchor") == "run_metadata")

        # New method: get metadata from DataFrame
        result = {}

        # Read modules
        if "modules" in metadata_df.columns and not metadata_df.get_column("modules").is_empty():
            result["modules"] = json.loads(metadata_df.get_column("modules")[0])

        # Read data sources
        if "data_sources" in metadata_df.columns and not metadata_df.get_column("data_sources").is_empty():
            result["data_sources"] = json.loads(metadata_df.get_column("data_sources")[0])

        # Read creation date
        if "creation_date" in metadata_df.columns and not metadata_df.get_column("creation_date").is_empty():
            result["creation_date"] = metadata_df.get_column("creation_date")[0]

        # Read config
        if "config" in metadata_df.columns and not metadata_df.get_column("config").is_empty():
            result["config"] = json.loads(metadata_df.get_column("config")[0])

        # Read software versions
        if "software_versions" in metadata_df.columns and not metadata_df.get_column("software_versions").is_empty():
            result["software_versions"] = json.loads(metadata_df.get_column("software_versions")[0])

        return result
    except Exception as e:
        logger.error(f"Error extracting report metadata from parquet: {e}")
        return None


def save_report_metadata() -> None:
    """
    Save all report metadata to the parquet file.

    This includes modules, data sources, creation date, config, and plot data.
    """
    # Prepare metadata row
    modules_data: List[Dict[str, Any]] = []
    for mod in report.modules:
        module_dict: Dict[str, Any] = {
            "name": mod.name,
            "anchor": str(mod.anchor),
            "info": mod.info,
            "intro": mod.intro,
            "comment": mod.comment,
            "sections": [section.model_dump(mode="json", exclude_none=True) for section in mod.sections],
            "versions": {
                name: [version for _, version in versions_tuples] for name, versions_tuples in mod.versions.items()
            },
        }

        modules_data.append(module_dict)

    # Prepare data sources
    data_sources_dict: Dict[str, Dict[str, Dict[str, str]]] = {}
    for mod_id, source_dict in report.data_sources.items():
        data_sources_dict[mod_id] = {}
        for section, sources in source_dict.items():
            data_sources_dict[mod_id][section] = {}
            for sname, source in sources.items():
                data_sources_dict[mod_id][section][sname] = source

    def _clean_config_values(value: Any) -> Any:
        """
        Convert any patterns to strings, recursively.
        """
        if isinstance(value, set):
            value = list(value)

        if isinstance(value, Pattern):
            return str(value)
        elif isinstance(value, dict):
            return {k: _clean_config_values(v) for k, v in value.items()}
        elif isinstance(value, list):
            return [_clean_config_values(item) for item in value]

        return value

    # Create MultiqcConfig object from config and dump it to dict
    config_dict = {}
    try:
        # Extract fields from config that match the schema
        config_fields = {key: getattr(config, key) for key in MultiQCConfig.model_fields if hasattr(config, key)}
        # Convert any patterns to strings, recursively
        config_fields = _clean_config_values(config_fields)
        # Create and validate config object
        config_obj = MultiQCConfig(**config_fields)
        # Dump to dict
        config_dict = config_obj.model_dump(mode="json", exclude_none=True)
    except ValidationError as e:
        logger.error(f"Error validating config: {e}")

    # Create metadata DataFrame
    metadata_df = pl.DataFrame(
        {
            "type": ["run_metadata"],
            "anchor": ["run_metadata"],
            "creation_date": [report.creation_date],
            "config": [json.dumps(config_dict)],
            "data_sources": [json.dumps(data_sources_dict)],
            "multiqc_version": [config.version if hasattr(config, "version") else ""],
            "modules": [json.dumps(modules_data)],
            "software_versions": [json.dumps(dict(report.software_versions))],
        }
    )

    append_to_parquet(metadata_df)

    # Flush all buffered data to the parquet file
    flush_to_parquet()


def _write_parquet(df: pl.DataFrame) -> None:
    parquet_file = tmp_dir.parquet_file()
    # Ensure directory exists
    os.makedirs(parquet_file.parent, exist_ok=True)

    # Write to file
    try:
        df.write_parquet(parquet_file, compression="gzip")
    except AttributeError:  # 'builtins.PyDataFrame' object has no attribute 'write_parquet'
        # Pyodide polars doesn't support write_parquet, fall back to CSV
        csv_file = parquet_file.with_suffix(".csv")
        logger.debug(f"Parquet writing not supported in Pyodide, falling back to CSV: {csv_file}")
        df.write_csv(csv_file)
    except Exception as e:
        logger.error(f"Error writing parquet file: {e}")
        raise


def _read_or_create_df() -> pl.DataFrame:
    parquet_file = tmp_dir.parquet_file()
    csv_file = parquet_file.with_suffix(".csv")

    # Try to read parquet first, then fall back to CSV
    if parquet_file.exists():
        try:
            return pl.read_parquet(parquet_file)
        except Exception as e:
            logger.error(f"Error reading parquet file: {e}")
            if config.strict:
                raise e
    elif csv_file.exists():
        try:
            # Read CSV and convert creation_date back to datetime
            df = pl.read_csv(csv_file)
            if "creation_date" in df.columns:
                df = df.with_columns(pl.col("creation_date").str.to_datetime())
            return df
        except Exception as e:
            logger.error(f"Error reading CSV file: {e}")
            if config.strict:
                raise e
    else:
        # Create directory if needed
        os.makedirs(parquet_file.parent, exist_ok=True)

    return pl.DataFrame(
        {
            "anchor": [],
            "type": [],
            "creation_date": [],
            "plot_type": [],
            "plot_input_data": [],
            "sample": [],
        },
        schema_overrides={
            "anchor": pl.Utf8,
            "type": pl.Utf8,
            "creation_date": pl.Datetime(time_unit="us"),
            "plot_type": pl.Utf8,
            "plot_input_data": pl.Utf8,
            "sample": pl.Utf8,
        },
    )


# def _update_parquet(df: pl.DataFrame, anchor: Anchor) -> None:
#     """
#     Update the parquet file with new data.

#     This function handles both creating a new file and updating an existing one.
#     """
#     parquet_file = tmp_dir.parquet_file()

#     # Check if file already exists
#     if parquet_file.exists():
#         try:
#             # Read existing data
#             existing_df = pl.read_parquet(parquet_file)

#             # Remove any existing data for this anchor
#             if "anchor" in existing_df.columns:
#                 existing_df = existing_df.filter(pl.col("anchor") != str(anchor))

#             # Append new data
#             merged_df = pl.concat([existing_df, df], how="diagonal")
#         except Exception as e:
#             logger.error(f"Error reading existing parquet file: {e}")
#             # If there was an error, just use the new data
#             merged_df = df
#     else:
#         # Create new file with just this data
#         merged_df = df

#     # Ensure directory exists
#     os.makedirs(parquet_file.parent, exist_ok=True)

#     # Write to file
#     try:
#         merged_df.write_parquet(parquet_file, compression="gzip")
#     except Exception as e:
#         logger.error(f"Error writing parquet file: {e}")
#         raise


def reset():
    """
    Reset the module state.
    """
    global _saved_anchors, _metric_col_names, _pending_dataframes, _pending_wide_tables
    _saved_anchors = set()
    _metric_col_names = set()
    _pending_dataframes = []
    _pending_wide_tables = []


def parse_value(value: Any, value_type: str) -> Any:
    """
    Parse a string value back to its original type, handling special markers.
    """
    import math

    # Handle special NaN marker
    if isinstance(value, str) and value == "__NAN__MARKER__":
        return math.nan

    # Handle regular type conversion
    if value_type == "int":
        return int(float(value))  # Handle float strings that represent integers
    elif value_type == "float":
        return float(value)
    elif value_type == "bool":
        return value.lower() == "true"
    else:
        return value
