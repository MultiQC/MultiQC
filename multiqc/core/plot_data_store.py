"""
Module for handling plot data storage in a single parquet file.
This replaces the individual plot parquet files with a single file that contains data from all plots.
It also stores all report metadata (modules, data sources, configs) to make reports fully reproducible.
"""

import json
import logging
import os
from pathlib import Path
from re import Pattern
from typing import Any, Dict, Optional, Sequence, Set, Union

import pandas as pd
import pyarrow.parquet as pq
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

# Metadata keys
META_MODULES = "modules"
META_DATA_SOURCES = "data_sources"
META_CREATION_DATE = "creation_date"
META_CONFIG = "config"
META_MULTIQC_VERSION = "multiqc_version"


def save_plot_data(anchor: Anchor, df: pd.DataFrame) -> None:
    """
    Save plot data to the parquet file.

    This function adds/updates data for a specific plot in the file.
    """
    # Update the global cache
    _plot_dataframes[anchor] = df
    _saved_anchors.add(anchor)

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
        if "timestamp" in metadata_df.columns and not metadata_df["timestamp"].empty:
            result["timestamp"] = metadata_df["timestamp"].iloc[0]

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
    modules_data = []
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

    # Update existing file or create new one
    if parquet_file.exists():
        try:
            existing_df = pd.read_parquet(parquet_file)

            # Remove any existing metadata
            existing_df = existing_df[existing_df["anchor"] != "run_metadata"]

            # Append new metadata
            merged_df = pd.concat([existing_df, metadata_df], ignore_index=True)

            # Write to file
            merged_df.to_parquet(parquet_file, compression="gzip")
        except Exception as e:
            logger.error(f"Error updating parquet file with metadata: {e}")
            if config.strict:
                raise e
            # If error, just write the metadata
            metadata_df.to_parquet(parquet_file, compression="gzip")
    else:
        # Create directory if needed
        os.makedirs(parquet_file.parent, exist_ok=True)

        # Write to file
        metadata_df.to_parquet(parquet_file, compression="gzip")

    logger.debug("Saved report metadata to parquet file")


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
    global _plot_dataframes, _saved_anchors
    _plot_dataframes = {}
    _saved_anchors = set()


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
