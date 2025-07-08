"""
Special case MultiQC module to load multiqc.parquet
It allows rerunning MultiQC when original data is gone, as well as extend
existing reports with new data.
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union, cast

import packaging.version
import polars as pl

from multiqc import config, report
from multiqc.base_module import BaseMultiqcModule, Section
from multiqc.core import plot_data_store
from multiqc.core.software_versions import parse_version, sort_versions
from multiqc.plots.bargraph import BarPlot, BarPlotInputData
from multiqc.plots.box import BoxPlot, BoxPlotInputData
from multiqc.plots.heatmap import HeatmapNormalizedInputData, HeatmapPlot
from multiqc.plots.linegraph import LinePlot, LinePlotNormalizedInputData
from multiqc.plots.plot import NormalizedPlotInputData, Plot
from multiqc.plots.scatter import ScatterNormalizedInputData, ScatterPlot
from multiqc.plots.violin import ViolinPlot, ViolinPlotInputData
from multiqc.types import Anchor, PlotType

log = logging.getLogger(__name__)


def create_plot_input_data_only(plot_input_data_dict: Dict) -> NormalizedPlotInputData:
    """
    Create only the plot input data object from a dictionary, without creating the plot object.
    """

    # Process the JSON data to replace NaN markers with proper NaN values
    def replace_nan_markers(obj: Any) -> Any:
        import math

        if isinstance(obj, str) and obj == "__NAN__MARKER__":
            return math.nan
        elif isinstance(obj, dict):
            return {k: replace_nan_markers(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [replace_nan_markers(i) for i in obj]
        return obj

    # Ensure we're working with a dictionary and apply the replacement to each value
    if not isinstance(plot_input_data_dict, dict):
        raise TypeError(f"Expected a dictionary, got {type(plot_input_data_dict)}")

    processed_dict = {k: replace_nan_markers(v) for k, v in plot_input_data_dict.items()}
    plot_type = PlotType.from_str(processed_dict["plot_type"])

    plot_input: NormalizedPlotInputData

    if plot_type == PlotType.LINE:
        plot_input = LinePlotNormalizedInputData(**processed_dict)
    elif plot_type == PlotType.BAR:
        plot_input = BarPlotInputData(**processed_dict)
    elif plot_type == PlotType.BOX:
        plot_input = BoxPlotInputData(**processed_dict)
    elif plot_type == PlotType.HEATMAP:
        plot_input = HeatmapNormalizedInputData(**processed_dict)
    elif plot_type == PlotType.VIOLIN or plot_type == PlotType.TABLE:
        plot_input = ViolinPlotInputData(**processed_dict)
    elif plot_type == PlotType.SCATTER:
        plot_input = ScatterNormalizedInputData(**processed_dict)
    else:
        raise ValueError(f"Unknown plot type: {plot_type}")

    return plot_input


def create_plot_from_input_data(plot_input: NormalizedPlotInputData) -> Union[Plot, str, None]:
    """
    Create a plot object from plot input data.
    """
    if plot_input.plot_type == PlotType.LINE:
        from multiqc.plots.linegraph import LinePlot

        return LinePlot.from_inputs(cast(LinePlotNormalizedInputData, plot_input))
    elif plot_input.plot_type == PlotType.BAR:
        from multiqc.plots.bargraph import BarPlot

        return BarPlot.from_inputs(cast(BarPlotInputData, plot_input))
    elif plot_input.plot_type == PlotType.VIOLIN or plot_input.plot_type == PlotType.TABLE:
        from multiqc.plots.violin import ViolinPlot

        return ViolinPlot.from_inputs(cast(ViolinPlotInputData, plot_input))
    elif plot_input.plot_type == PlotType.BOX:
        from multiqc.plots.box import BoxPlot

        return BoxPlot.from_inputs(cast(BoxPlotInputData, plot_input))
    elif plot_input.plot_type == PlotType.SCATTER:
        from multiqc.plots.scatter import ScatterPlot

        return ScatterPlot.from_inputs(cast(ScatterNormalizedInputData, plot_input))
    elif plot_input.plot_type == PlotType.HEATMAP:
        from multiqc.plots.heatmap import HeatmapPlot

        return HeatmapPlot.from_inputs(cast(HeatmapNormalizedInputData, plot_input))
    else:
        log.warning(f"Unknown plot type {plot_input.plot_type}, cannot create plot object")
        return None


def load_plot_input(plot_input_data_dict: Dict) -> Tuple[NormalizedPlotInputData, Union[Plot, str, None]]:
    """
    Load plot input data and create plot object from a dictionary.
    This function combines create_plot_input_data_only and create_plot_from_input_data.
    """
    plot_input = create_plot_input_data_only(plot_input_data_dict)
    plot = create_plot_from_input_data(plot_input)
    return plot_input, plot


class LoadMultiqcData(BaseMultiqcModule):
    def __init__(self):
        super(LoadMultiqcData, self).__init__(
            name="MultiQC Data",
            anchor=Anchor("multiqc_data"),
            info="loads multiqc data",
        )

        # Dictionary to collect all software versions from all parquet files
        self.collected_software_versions: Dict[str, List[str]] = {}

        # First, try to find parquet file
        parquet_files = self.find_log_files("multiqc_data")
        if parquet_files:
            for f in parquet_files:
                self.load_parquet_file(Path(f["root"]) / f["fn"])

            # After loading all files, process and deduplicate software versions
            self._process_collected_software_versions()

    def load_parquet_file(self, path: Union[str, Path]):
        """
        Load a multiqc.parquet file containing all report data.
        """
        path = Path(path)
        assert path.suffix == ".parquet"
        log.debug(f"Loading report data from parquet file: {path}")

        try:
            # Read the entire parquet file
            df = pl.read_parquet(path)

            # Extract metadata from parquet
            metadata = plot_data_store.get_report_metadata(df)
            if metadata is None:
                log.error(f"Failed to extract metadata from parquet file: {path}")
                return

            # Load modules
            if "modules" in metadata:
                for mod_dict in metadata["modules"]:
                    # Extract module data first so we can use it for section defaults
                    anchor = mod_dict.get("anchor", "")
                    name = mod_dict.get("name", "")
                    info = mod_dict.get("info", "")
                    intro = mod_dict.get("intro", "")
                    comment = mod_dict.get("comment", "")

                    # Create sections, providing default values for missing required fields
                    sections = []
                    for section_data in mod_dict.pop("sections"):
                        # Ensure required fields have default values if missing
                        section_data.setdefault("id", section_data.get("anchor", ""))
                        section_data.setdefault("module", name)  # Use the module name
                        section_data.setdefault("module_anchor", anchor)  # Use the module anchor
                        section_data.setdefault("description", "")
                        sections.append(Section(**section_data))

                    # Convert versions to expected format
                    versions: Dict[str, List[Tuple[Optional[packaging.version.Version], str]]] = {}
                    if "versions" in mod_dict:
                        versions_data = mod_dict.pop("versions")
                        versions = {
                            name: [(None, version) for version in versions] for name, versions in versions_data.items()
                        }

                    # Remove already extracted data from mod_dict
                    mod_dict.pop("anchor", None)
                    mod_dict.pop("name", None)
                    mod_dict.pop("info", None)
                    mod_dict.pop("intro", None)
                    mod_dict.pop("comment", None)

                    # Special handling for Software Versions modules - skip them to avoid duplicates
                    if anchor == "multiqc_software_versions" and name == "Software Versions":
                        # Extract software versions from the HTML content
                        if sections:
                            log.debug("Extracting software versions from HTML content in parquet file")
                            import re

                            for section in sections:
                                if section.content:
                                    # Parse the HTML table to extract software versions
                                    # Look for table rows with pattern: <td>Software</td><td><samp>Version</samp></td>
                                    row_pattern = r"<tr><td>([^<]+)</td><td><samp>([^<]+)</samp></td></tr>"
                                    matches = re.findall(row_pattern, section.content)

                                    for software_name, version in matches:
                                        # Collect software versions for later processing
                                        if software_name not in self.collected_software_versions:
                                            self.collected_software_versions[software_name] = []
                                        self.collected_software_versions[software_name].append(version)
                                        log.debug(f"Collected software version: {software_name} = {version}")

                        log.debug(
                            "Skipping Software Versions module from parquet file - extracted data to global software_versions"
                        )
                        continue  # Skip adding this module to report.modules

                    # Create module
                    mod = BaseMultiqcModule(name=name, anchor=Anchor(anchor), info=info)
                    mod.sections = sections
                    mod.versions = versions
                    mod.intro = intro
                    mod.comment = comment

                    # Check for duplicate modules and merge if found (same logic as exec_modules.py)
                    existing_mod = None
                    for prev_mod in report.modules:
                        if prev_mod.anchor == mod.anchor:
                            existing_mod = prev_mod
                            break

                    if existing_mod:
                        # Merge the existing module into the new one
                        mod.merge(existing_mod)  # This only merges versions

                        # Debug: Log sections before merging
                        log.debug(f"Before merging - Existing module has {len(existing_mod.sections)} sections:")
                        for s in existing_mod.sections:
                            log.debug(f"  - Existing: {s.name} (anchor: {s.anchor})")
                        log.debug(f"Before merging - New module has {len(mod.sections)} sections:")
                        for s in mod.sections:
                            log.debug(f"  - New: {s.name} (anchor: {s.anchor})")

                        # Merge sections based on anchor - keep all unique sections from both modules
                        existing_sections = {s.anchor: s for s in existing_mod.sections}
                        new_sections = {s.anchor: s for s in mod.sections}

                        merged_sections = []
                        all_section_anchors = set(existing_sections.keys()) | set(new_sections.keys())

                        log.debug(f"All section anchors to process: {sorted(all_section_anchors)}")

                        for section_anchor in all_section_anchors:
                            if section_anchor in existing_sections and section_anchor in new_sections:
                                # Both modules have this section - merge content if different
                                log.debug(f"Merging content for section: {section_anchor}")
                                existing_section = existing_sections[section_anchor]
                                new_section = new_sections[section_anchor]

                                # Determine which section has actual data
                                existing_has_data = existing_section.plot_anchor is not None or (
                                    existing_section.content and existing_section.content.strip()
                                )
                                new_has_data = new_section.plot_anchor is not None or (
                                    new_section.content and new_section.content.strip()
                                )

                                if existing_has_data and not new_has_data:
                                    # Existing section has data, new section is empty - use existing as base
                                    log.debug(
                                        f"Using existing section as base (new section is empty): {section_anchor}"
                                    )
                                    merged_sections.append(existing_section)
                                elif new_has_data and not existing_has_data:
                                    # New section has data, existing section is empty - use new as base
                                    log.debug(
                                        f"Using new section as base (existing section is empty): {section_anchor}"
                                    )
                                    merged_sections.append(new_section)
                                elif existing_has_data and new_has_data:
                                    # Both sections have data - perform proper merging
                                    log.debug(f"Both sections have data, merging content: {section_anchor}")

                                    # Combine content if it's different
                                    merged_content = new_section.content
                                    if existing_section.content and existing_section.content != new_section.content:
                                        # If content is different, append existing content to new content
                                        if merged_content:
                                            merged_content = merged_content + "\n\n" + existing_section.content
                                        else:
                                            merged_content = existing_section.content

                                    # Preserve plot-related attributes from whichever section has them
                                    merged_plot_anchor = new_section.plot_anchor or existing_section.plot_anchor

                                    # Create merged section with combined content
                                    merged_section = Section(
                                        name=new_section.name,
                                        anchor=new_section.anchor,
                                        id=new_section.id,
                                        description=new_section.description,
                                        module=new_section.module,
                                        module_anchor=new_section.module_anchor,
                                        module_info=new_section.module_info,
                                        comment=new_section.comment,
                                        helptext=new_section.helptext,
                                        content_before_plot=new_section.content_before_plot,
                                        content=merged_content,
                                        print_section=new_section.print_section,
                                        plot_anchor=merged_plot_anchor,  # Use merged plot_anchor
                                        ai_summary=new_section.ai_summary,
                                    )
                                    merged_sections.append(merged_section)
                                else:
                                    # Both sections are empty - use new section as default
                                    log.debug(f"Both sections are empty, using new section: {section_anchor}")
                                    merged_sections.append(new_section)
                            elif section_anchor in existing_sections:
                                # Only in existing module
                                log.debug(f"Preserving section from existing module: {section_anchor}")
                                merged_sections.append(existing_sections[section_anchor])
                            else:
                                # Only in new module
                                log.debug(f"Adding section from new module: {section_anchor}")
                                merged_sections.append(new_sections[section_anchor])

                        mod.sections = merged_sections

                        # Debug: Log sections after merging
                        log.debug(f"After merging - Final module has {len(mod.sections)} sections:")
                        for s in mod.sections:
                            log.debug(f"  - Final: {s.name} (anchor: {s.anchor})")

                        log.debug(f'Updating module "{existing_mod.name}" with data from parquet')
                        report.modules.remove(existing_mod)
                    else:
                        log.debug(f"Loading module {mod.name} from parquet")

                    report.modules.append(mod)

            # Load global software versions data
            if "software_versions" in metadata:
                software_versions_data = metadata["software_versions"]
                log.debug("Loading global software versions data from parquet file")
                for group_name, group_versions in software_versions_data.items():
                    for software_name, versions_list in group_versions.items():
                        # Collect software versions for later processing
                        if software_name not in self.collected_software_versions:
                            self.collected_software_versions[software_name] = []
                        self.collected_software_versions[software_name].extend(versions_list)

            # Load data sources
            if "data_sources" in metadata:
                for mod_id, source_dict in metadata["data_sources"].items():
                    for section_name, sources in source_dict.items():
                        for sname, source in sources.items():
                            report.data_sources[mod_id][section_name][sname] = source

            # Set creation date
            if "creation_date" in metadata:
                try:
                    # Convert the datetime to a Python datetime object
                    if isinstance(metadata["creation_date"], datetime):
                        report.creation_date = metadata["creation_date"]
                    else:
                        creation_date_str = str(metadata["creation_date"])
                        # Use standard Python datetime parsing
                        report.creation_date = datetime.fromisoformat(creation_date_str.replace("Z", "+00:00"))
                except ValueError as e:
                    log.error(f"Could not parse creation date: {metadata['creation_date']}, error: {e}")

            if "config" in metadata:
                pass  # We do not load config, but keep the current one

            # Load plot input data from plot_input rows in the dataframe
            if "type" in df.columns and "plot_input_data" in df.columns:
                plot_input_rows = df.filter(pl.col("type") == "plot_input")
                for row in plot_input_rows.iter_rows(named=True):
                    anchor = Anchor(str(row["anchor"]))
                    plot_input_data = row["plot_input_data"]
                    plot_input_data_dict = json.loads(plot_input_data)
                    try:
                        # First, create only the plot input data (without the plot object)
                        plot_input = create_plot_input_data_only(plot_input_data_dict)
                        log.debug(f"Loading plot input data for {anchor}, type: {plot_input.__class__.__name__}")

                        # Check if we already have plot input data for this anchor (from previous parquet files)
                        if anchor in report.plot_input_data:
                            # Merge the existing data with the new data using the appropriate merge method
                            existing_plot_input = report.plot_input_data[anchor]
                            log.debug(
                                f"Found existing plot input data for {anchor}, merging {existing_plot_input.__class__.__name__} with {plot_input.__class__.__name__}"
                            )

                            # Use the merge class method from the plot input data class
                            merged_plot_input = existing_plot_input.__class__.merge(existing_plot_input, plot_input)
                            report.plot_input_data[anchor] = merged_plot_input
                            log.debug(f"Successfully merged plot input data for {anchor}")

                            # Create the plot object from the merged data (this ensures proper color assignment)
                            merged_plot: Union[Plot, str, None] = create_plot_from_input_data(merged_plot_input)
                            if merged_plot is not None:
                                report.plot_by_id[anchor] = merged_plot
                                log.debug(f"Updated plot object for {anchor} with merged data")

                        else:
                            # No existing data, just add the new data and create plot object
                            log.debug(f"No existing plot input data for {anchor}, adding new data")
                            report.plot_input_data[anchor] = plot_input

                            # Create the plot object from the input data
                            plot = create_plot_from_input_data(plot_input)
                            if plot is not None:
                                report.plot_by_id[anchor] = plot
                    except Exception as e:
                        log.error(f"Error loading plot input data {anchor}: {e}")
                        if config.strict:
                            raise e

        except Exception as e:
            log.error(f"Error loading data from parquet file: {e}")
            if config.strict:
                raise e

    def _process_collected_software_versions(self):
        """
        Process collected software versions: parse, deduplicate, sort, and add to report.software_versions
        """
        for software_name, version_strings in self.collected_software_versions.items():
            if version_strings:
                # Parse versions and create tuples of (Version object, version string)
                ver_and_verstr = []
                for version_str in version_strings:
                    # Split comma-separated versions (in case they're combined)
                    individual_versions = [v.strip() for v in version_str.split(",") if v.strip()]
                    for individual_version in individual_versions:
                        ver_and_verstr.append((parse_version(individual_version), individual_version))

                # Remove duplicates by converting to set and back to list
                ver_and_verstr = list(set(ver_and_verstr))

                # Sort versions using MultiQC's sorting logic
                sorted_versions = sort_versions(ver_and_verstr)

                # Add to global software_versions (using software name as group name)
                final_versions = [version_str for _, version_str in sorted_versions]
                report.software_versions[software_name][software_name].extend(final_versions)

                log.debug(f"Processed {len(final_versions)} versions for {software_name}: {', '.join(final_versions)}")
