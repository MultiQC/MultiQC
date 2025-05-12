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

import pandas as pd
from pydantic import ValidationError  # type: ignore

from multiqc import config, report
from multiqc.core import tmp_dir
from multiqc.types import Anchor
from multiqc.utils.config_schema import MultiQCConfig

logger = logging.getLogger(__name__)

# Global cache of dataframes for each plot
_plot_dataframes: Dict[Anchor, pd.DataFrame] = {}
# Set to keep track of which anchors have been saved
_saved_anchors: Set[Anchor] = set()
# Cache of merged wide-format table data (sample-based tables)
_merged_wide_df: Optional[pd.DataFrame] = None

# Metadata keys
META_MODULES = "modules"
META_DATA_SOURCES = "data_sources"
META_CREATION_DATE = "creation_date"
META_CONFIG = "config"
META_MULTIQC_VERSION = "multiqc_version"


def merge_wide_tables(df: pd.DataFrame) -> None:
    """
    Merge wide-format table data with existing sample-based tables.

    This function extracts table rows from the dataframe and merges them with
    the existing global cache of wide-format table data. This ensures all tables
    that have the same samples get combined into a single row per sample.
    """
    global _merged_wide_df

    # Extract table rows
    if "type" not in df.columns:
        return

    wide_df = df[df["type"] == "table_row"].copy()
    if wide_df.empty:
        return

    # Check if we have necessary columns
    if "sample_name" not in wide_df.columns:
        return

    # Add table ID to column names to avoid collisions between tables
    table_id = str(df["anchor"].iloc[0])

    # If the dataframe already has a table_id column, use that instead
    if "table_id" in wide_df.columns and not wide_df["table_id"].isna().all():
        # Get the first non-null table_id
        table_ids = wide_df["table_id"].dropna().unique()
        if len(table_ids) > 0:
            table_id = str(table_ids[0])
            # Log a warning if there are multiple different table_ids
            if len(table_ids) > 1:
                logger.warning(f"Multiple table IDs found in dataframe: {table_ids}. Using {table_id}")

    col_prefix = f"tbl_{table_id}_"

    # Rename metric columns to include table ID
    for col in wide_df.columns:
        if col.startswith("col_"):
            wide_df.rename(columns={col: f"{col_prefix}{col}"}, inplace=True)

    # Merge with existing data if available
    if _merged_wide_df is None:
        _merged_wide_df = wide_df
    else:
        # Merge on sample_name, preserving all columns from both dataframes
        _merged_wide_df = pd.merge(_merged_wide_df, wide_df, on=["sample_name", "creation_date"], how="outer")

    # Remove the wide format rows from the original dataframe
    df.drop(df[df["type"] == "table_row"].index, inplace=True)


def save_merged_tables() -> None:
    """
    Save the merged wide-format table data to the parquet file.
    """
    global _merged_wide_df

    if _merged_wide_df is None or _merged_wide_df.empty:
        return

    parquet_file = tmp_dir.parquet_file()

    # Update existing file or create new one
    if parquet_file.exists():
        try:
            # Read existing data
            existing_df = pd.read_parquet(parquet_file)

            # Remove any existing merged table rows
            if "type" in existing_df.columns:
                existing_df = existing_df[existing_df["type"] != "merged_table_row"]

            # Update the type field in our merged table
            _merged_wide_df["type"] = "merged_table_row"

            # Append merged data
            merged_df = pd.concat([existing_df, _merged_wide_df], ignore_index=True)
        except Exception as e:
            logger.error(f"Error updating parquet file with merged tables: {e}")
            if config.strict:
                raise e
            return
    else:
        # Create directory if needed
        os.makedirs(parquet_file.parent, exist_ok=True)

        # Update the type field and prepare for writing
        _merged_wide_df["type"] = "merged_table_row"
        merged_df = _merged_wide_df

    # Fix for Iceberg. Iceberg never keeps an arbitrary zone offset in the data –
    # a value that has a zone is normalised to UTC, and the zone itself is discarded.
    merged_df["creation_date"] = (
        pd.to_datetime(merged_df["creation_date"], utc=True)
        .dt.floor("us")  # tz-aware (+02:00)
        .dt.tz_localize(None)  # …but drop the zone
        .astype("datetime64[us]")  # make it explicit
    )

    # Write to file
    try:
        merged_df.to_parquet(parquet_file, compression="gzip")
        logger.debug(f"Saved merged table data to parquet file {parquet_file}")
    except Exception as e:
        logger.error(f"Error writing merged table data to parquet file: {e}")
        if config.strict:
            raise e


def save_plot_data(anchor: Anchor, df: pd.DataFrame) -> None:
    """
    Save plot data to the parquet file.

    This function adds/updates data for a specific plot in the file.
    """
    # Update the global cache
    _plot_dataframes[anchor] = df
    _saved_anchors.add(anchor)

    # Merge wide-format table data if present
    merge_wide_tables(df)

    # Write to the file
    _update_parquet(df, anchor)


def get_report_metadata(df: pd.DataFrame) -> Optional[Dict[str, Any]]:
    """
    Extract all report metadata from the parquet file.

    Args:
        file_path: Optional path to the parquet file. If None, uses the default path.

    Returns a dictionary with modules, data_sources, creation_date, and config.
    """
    if df.empty:
        return None

    try:
        # Read the metadata table from the parquet file
        metadata_df = df[df["anchor"] == "run_metadata"]

        # New method: get metadata from DataFrame
        result = {}

        # Read modules
        if "modules" in metadata_df.columns and not metadata_df["modules"].empty:
            result["modules"] = json.loads(metadata_df["modules"].iloc[0])

        # Read data sources
        if "data_sources" in metadata_df.columns and not metadata_df["data_sources"].empty:
            result["data_sources"] = json.loads(metadata_df["data_sources"].iloc[0])

        # Read creation date
        if "creation_date" in metadata_df.columns and not metadata_df["creation_date"].empty:
            result["creation_date"] = metadata_df["creation_date"].iloc[0]

        # Read config
        if "config" in metadata_df.columns and not metadata_df["config"].empty:
            result["config"] = json.loads(metadata_df["config"].iloc[0])

        return result
    except Exception as e:
        logger.error(f"Error extracting report metadata from parquet: {e}")
        return None


def save_report_metadata() -> None:
    """
    Save all report metadata to the parquet file.

    This includes modules, data sources, creation date, config, and plot data.
    """
    parquet_file = tmp_dir.parquet_file()

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

    # Creation date
    # try:
    #     creation_date_str = report.creation_date.strftime("%Y-%m-%d, %H:%M %Z")
    # except UnicodeEncodeError:
    #     # Fall back to a format without timezone if we encounter encoding issues
    #     creation_date_str = report.creation_date.strftime("%Y-%m-%d, %H:%M")

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
    metadata_df = pd.DataFrame(
        {
            "type": ["run_metadata"],
            "anchor": ["run_metadata"],
            "creation_date": [report.creation_date],
            "config": [json.dumps(config_dict)],
            "data_sources": [json.dumps(data_sources_dict)],
            "multiqc_version": [config.version if hasattr(config, "version") else ""],
            "modules": [json.dumps(modules_data)],
        }
    )

    df = metadata_df

    # Update existing file or create new one
    if parquet_file.exists():
        try:
            existing_df = pd.read_parquet(parquet_file)

            # Remove any existing metadata
            existing_df = existing_df[existing_df["anchor"] != "run_metadata"]

            # Append new metadata
            merged_df = pd.concat([existing_df, metadata_df], ignore_index=True)
        except Exception as e:
            logger.error(f"Error updating parquet file with metadata: {e}")
            if config.strict:
                raise e
        else:
            df = merged_df
    else:
        # Create directory if needed
        os.makedirs(parquet_file.parent, exist_ok=True)

    # Fix for Iceberg. Iceberg never keeps an arbitrary zone offset in the data –
    # a value that has a zone is normalised to UTC, and the zone itself is discarded.
    df["creation_date"] = (
        pd.to_datetime(df["creation_date"], utc=True)
        .dt.floor("us")  # tz-aware (+02:00)
        .dt.tz_localize(None)  # …but drop the zone
        .astype("datetime64[us]")  # make it explicit
    )

    # Write to file
    df.to_parquet(parquet_file, compression="gzip")
    logger.debug(f"Saved parquet file {parquet_file}")

    # Save the merged table data
    save_merged_tables()


def _update_parquet(df: pd.DataFrame, anchor: Anchor) -> None:
    """
    Update the parquet file with new data.

    This function handles both creating a new file and updating an existing one.
    """
    parquet_file = tmp_dir.parquet_file()

    # Check if file already exists
    if parquet_file.exists():
        try:
            # Read existing data
            existing_df = pd.read_parquet(parquet_file)

            # Remove any existing data for this anchor
            if "anchor" in existing_df.columns:
                existing_df = existing_df[existing_df["anchor"] != str(anchor)]

            # Append new data
            merged_df = pd.concat([existing_df, df], ignore_index=True)
        except Exception as e:
            logger.error(f"Error reading existing parquet file: {e}")
            # If there was an error, just use the new data
            merged_df = df
    else:
        # Create new file with just this data
        merged_df = df

    # Ensure directory exists
    os.makedirs(parquet_file.parent, exist_ok=True)

    # Write to file
    try:
        merged_df.to_parquet(parquet_file, compression="gzip")
    except Exception as e:
        logger.error(f"Error writing parquet file: {e}")
        raise


def reset():
    """
    Reset the module state.
    """
    global _plot_dataframes, _saved_anchors, _merged_wide_df
    _plot_dataframes = {}
    _saved_anchors = set()
    _merged_wide_df = None


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
