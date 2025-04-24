"""
Module for handling plot data storage in a single parquet file.
This replaces the individual plot parquet files with a single file that contains data from all plots.
It also stores all report metadata (modules, data sources, configs) to make reports fully reproducible.
"""

import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, Optional, Set, Union

import pandas as pd
import pyarrow.parquet as pq  # type: ignore
from cloudpathlib import CloudPath

from multiqc import config, report
from multiqc.core import tmp_dir
from multiqc.types import Anchor

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


def get_report_metadata(file_path: Optional[Union[str, Path]] = None) -> Optional[Dict[str, Any]]:
    """
    Extract all report metadata from the parquet file.

    Args:
        file_path: Optional path to the parquet file. If None, uses the default path.

    Returns a dictionary with modules, data_sources, creation_date, and config.
    """
    parquet_file: Union[CloudPath, Path]
    if file_path is None:
        parquet_file = tmp_dir.parquet_file()
    else:
        parquet_file = Path(file_path)

    if not parquet_file.exists():
        return None

    try:
        # Read the metadata table from the parquet file
        df = pd.read_parquet(parquet_file)
        metadata_df = df[df["anchor"] == "metadata"]

        if metadata_df.empty:
            # Try to read from metadata (legacy)
            table = pq.read_table(parquet_file, columns=[])
            if not table.schema.metadata:
                return None

            metadata = {k.decode("utf-8"): v.decode("utf-8") for k, v in table.schema.metadata.items()}

            result = {}

            # Extract metadata
            if META_MODULES in metadata:
                result["modules"] = json.loads(metadata[META_MODULES])

            if META_DATA_SOURCES in metadata:
                result["data_sources"] = json.loads(metadata[META_DATA_SOURCES])

            if META_CREATION_DATE in metadata:
                result["creation_date"] = metadata[META_CREATION_DATE]

            if META_CONFIG in metadata:
                result["config"] = json.loads(metadata[META_CONFIG])

            return result

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
    modules_data = []
    for mod in report.modules:
        module_dict: Dict[str, Any] = {
            "name": mod.name,
            "anchor": str(mod.anchor),
            "info": mod.info,
            "intro": mod.intro,
            "comment": mod.comment,
        }

        # Add sections
        module_dict["sections"] = [section.model_dump() for section in mod.sections]

        # Add software versions
        module_dict["versions"] = {
            name: [version for _, version in versions_tuples] for name, versions_tuples in mod.versions.items()
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
    try:
        creation_date_str = report.creation_date.strftime("%Y-%m-%d, %H:%M %Z")
    except UnicodeEncodeError:
        # Fall back to a format without timezone if we encounter encoding issues
        creation_date_str = report.creation_date.strftime("%Y-%m-%d, %H:%M")

    # Config data
    config_dict = {}
    for key in [
        "title",
        "subtitle",
        "intro_text",
        "report_comment",
        "report_header_info",
        "short_version",
        "version",
        "git_hash",
        "analysis_dir",
        "output_dir",
    ]:
        if hasattr(config, key):
            value = getattr(config, key)
            # Convert Path objects to strings
            if isinstance(value, Path):
                value = str(value)
            elif isinstance(value, list) and all(isinstance(x, Path) for x in value):
                value = [str(x) for x in value]
            config_dict[key] = value

    # Plot data
    plot_data_dict = {}
    for anchor, plot_data in report.plot_data.items():
        plot_data_dict[anchor] = plot_data

    # Create metadata DataFrame
    metadata_df = pd.DataFrame(
        {
            "anchor": ["metadata"],
            "type": ["metadata"],
            "modules": [json.dumps(modules_data)],
            "data_sources": [json.dumps(data_sources_dict)],
            "creation_date": [creation_date_str],
            "config": [json.dumps(config_dict)],
            "plot_data": [json.dumps(plot_data_dict)],
            "multiqc_version": [config.version if hasattr(config, "version") else ""],
        }
    )

    # Update existing file or create new one
    if parquet_file.exists():
        try:
            existing_df = pd.read_parquet(parquet_file)

            # Remove any existing metadata
            existing_df = existing_df[existing_df["anchor"] != "metadata"]

            # Append new metadata
            merged_df = pd.concat([existing_df, metadata_df], ignore_index=True)

            # Write to file
            merged_df.to_parquet(parquet_file, compression="gzip")
        except Exception as e:
            logger.error(f"Error updating parquet file with metadata: {e}")
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
        merged_df.to_parquet(str(parquet_file), compression="gzip")
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
    if value_type == "int":
        return int(value)
    elif value_type == "float":
        return float(value)
    elif value_type == "bool":
        return value.lower() == "true"
    else:
        return value
