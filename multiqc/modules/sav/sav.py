"""
MultiQC module for Illumina SAV (Sequencing Analysis Viewer) metrics.

This core module parses basic run information from XML files.
For advanced InterOp-based visualizations, install the multiqc_sav plugin.
"""

import logging
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.core import plugin_hooks
from multiqc.plots.table_object import ColumnDict

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Illumina SAV (Sequencing Analysis Viewer) provides quality metrics from Illumina sequencers.

    This module parses basic run information from `RunInfo.xml` and `RunParameters.xml` files
    found in Illumina sequencing run directories.

    For advanced visualizations including Q-score heatmaps, cluster density plots,
    intensity graphs, and lane summary tables, install the
    [multiqc_sav plugin](https://github.com/MultiQC/MultiQC_SAV):

    ```bash
    pip install multiqc_sav
    ```
    """

    def __init__(self):
        super().__init__(
            name="Illumina SAV",
            anchor="sav",
            href="https://support.illumina.com/sequencing/sequencing_software/sequencing_analysis_viewer_sav.html",
            info="Sequencing metrics from Illumina sequencers.",
        )

        # Store run data by sample/run name
        self.data_by_sample: Dict[str, Dict] = {}

        # Track directories we've processed to avoid duplicates
        processed_dirs: set = set()

        # Find and parse RunInfo.xml files
        for f in self.find_log_files("sav/runinfo"):
            run_dir = Path(f["root"])
            if run_dir in processed_dirs:
                continue

            parsed_data = self.parse_run_info(f)
            if parsed_data:
                # Use run ID or directory name as sample name
                s_name = self.clean_s_name(parsed_data.get("run_id") or run_dir.name, f)
                self.data_by_sample[s_name] = parsed_data
                self.add_data_source(f, s_name)
                processed_dirs.add(run_dir)

        # Find and parse RunParameters.xml files to augment data
        for f in self.find_log_files("sav/runparameters"):
            run_dir = Path(f["root"])
            parsed_params = self.parse_run_parameters(f)
            if parsed_params:
                # Find matching run by directory
                s_name = None
                for name, data in self.data_by_sample.items():
                    if data.get("_run_dir") == str(run_dir):
                        s_name = name
                        break

                if s_name:
                    # Merge parameters into existing data
                    self.data_by_sample[s_name].update(parsed_params)
                    self.add_data_source(f, s_name)
                else:
                    # Create new entry if RunInfo wasn't found
                    s_name = self.clean_s_name(run_dir.name, f)
                    parsed_params["_run_dir"] = str(run_dir)
                    self.data_by_sample[s_name] = parsed_params
                    self.add_data_source(f, s_name)
                    processed_dirs.add(run_dir)

        # Apply sample filtering
        self.data_by_sample = self.ignore_samples(self.data_by_sample)

        if len(self.data_by_sample) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.data_by_sample)} Illumina run(s)")

        # Extract software version if available
        for s_name, data in self.data_by_sample.items():
            version = data.get("rta_version") or data.get("application_version")
            self.add_software_version(version, s_name)

        # Add run info section first
        self.add_run_info_section()

        # Call plugin hook to allow extensions to add InterOp visualizations
        # Plugins can add sections and augment data_by_sample
        if "sav_extra" in plugin_hooks.hook_functions:
            log.debug("Calling sav_extra plugin hooks")
            for hook_fn in plugin_hooks.hook_functions["sav_extra"]:
                hook_fn(self)
        else:
            log.info("Install multiqc_sav for additional InterOp visualizations")

        # Setup and add general stats
        self.setup_general_stats_headers()
        self.general_stats_addcols(self.data_by_sample, self.genstat_headers)

        # Write parsed data to a file (at the end, after plugins may have augmented it)
        clean_data = {}
        for s_name, data in self.data_by_sample.items():
            clean_data[s_name] = {k: v for k, v in data.items() if not k.startswith("_")}
        self.write_data_file(clean_data, "multiqc_sav")

    def parse_run_info(self, f) -> Optional[Dict]:
        """Parse RunInfo.xml file to extract run configuration."""
        try:
            root = ET.fromstring(f["f"])
        except ET.ParseError as e:
            log.warning(f"Could not parse RunInfo.xml {f['fn']}: {e}")
            return None

        parsed_data: Dict = {"_run_dir": f["root"]}

        # Find Run element
        run_elem = root.find("Run")
        if run_elem is None:
            log.warning(f"No Run element found in {f['fn']}")
            return None

        # Extract run attributes
        parsed_data["run_number"] = run_elem.get("Number")
        parsed_data["run_id"] = run_elem.get("Id")

        # Flowcell
        flowcell = run_elem.find("Flowcell")
        if flowcell is not None and flowcell.text:
            parsed_data["flowcell"] = flowcell.text

        # Instrument
        instrument = run_elem.find("Instrument")
        if instrument is not None and instrument.text:
            parsed_data["instrument"] = instrument.text

        # Date
        date_elem = run_elem.find("Date")
        if date_elem is not None and date_elem.text:
            parsed_data["run_date_raw"] = date_elem.text
            parsed_data["run_date"] = self._parse_run_date(date_elem.text)

        # Read configuration
        reads_elem = run_elem.find("Reads")
        if reads_elem is not None:
            reads = []
            total_cycles = 0
            for read in reads_elem.findall("Read"):
                read_info = {
                    "number": read.get("Number"),
                    "cycles": int(read.get("NumCycles", 0)),
                    "is_indexed": read.get("IsIndexedRead", "N") == "Y",
                }
                reads.append(read_info)
                if not read_info["is_indexed"]:
                    total_cycles += read_info["cycles"]
            parsed_data["reads"] = reads
            parsed_data["total_cycles"] = total_cycles
            parsed_data["num_reads"] = len([r for r in reads if not r["is_indexed"]])
            parsed_data["num_indexes"] = len([r for r in reads if r["is_indexed"]])

        # FlowcellLayout
        layout = run_elem.find("FlowcellLayout")
        if layout is not None:
            parsed_data["lane_count"] = int(layout.get("LaneCount", 0))
            parsed_data["surface_count"] = int(layout.get("SurfaceCount", 0))
            parsed_data["swath_count"] = int(layout.get("SwathCount", 0))
            parsed_data["tile_count"] = int(layout.get("TileCount", 0))

        return parsed_data

    def parse_run_parameters(self, f) -> Optional[Dict]:
        """Parse RunParameters.xml file to extract additional run information."""
        try:
            root = ET.fromstring(f["f"])
        except ET.ParseError as e:
            log.warning(f"Could not parse RunParameters.xml {f['fn']}: {e}")
            return None

        parsed_data: Dict = {}

        # Try various element names used by different instruments
        rta_version = self._find_element_text(root, ["RTAVersion", "RtaVersion"])
        if rta_version:
            parsed_data["rta_version"] = rta_version

        app_version = self._find_element_text(root, ["ApplicationVersion", "Application Version"])
        if app_version:
            parsed_data["application_version"] = app_version

        exp_name = self._find_element_text(root, ["ExperimentName", "Experiment Name"])
        if exp_name:
            parsed_data["experiment_name"] = exp_name

        run_id = self._find_element_text(root, ["RunId", "RunID"])
        if run_id:
            parsed_data["run_id_params"] = run_id

        chemistry = self._find_element_text(root, ["Chemistry"])
        if chemistry:
            parsed_data["chemistry"] = chemistry

        # Planned read lengths (NovaSeq format)
        planned_read1 = self._find_element_text(root, ["PlannedRead1Cycles", "Read1NumberOfCycles", "Read1"])
        if planned_read1:
            try:
                parsed_data["planned_read1_cycles"] = int(planned_read1)
            except ValueError:
                pass

        planned_read2 = self._find_element_text(root, ["PlannedRead2Cycles", "Read2NumberOfCycles", "Read2"])
        if planned_read2:
            try:
                parsed_data["planned_read2_cycles"] = int(planned_read2)
            except ValueError:
                pass

        planned_index1 = self._find_element_text(root, ["PlannedIndex1ReadCycles", "IndexRead1NumberOfCycles"])
        if planned_index1:
            try:
                parsed_data["planned_index1_cycles"] = int(planned_index1)
            except ValueError:
                pass

        planned_index2 = self._find_element_text(root, ["PlannedIndex2ReadCycles", "IndexRead2NumberOfCycles"])
        if planned_index2:
            try:
                parsed_data["planned_index2_cycles"] = int(planned_index2)
            except ValueError:
                pass

        output_folder = self._find_element_text(root, ["OutputFolder", "OutputDirectory"])
        if output_folder:
            parsed_data["output_folder"] = output_folder

        return parsed_data if parsed_data else None

    def _find_element_text(self, root: ET.Element, names: list) -> Optional[str]:
        """Find first matching element from a list of possible names."""
        for name in names:
            # Try direct child
            elem = root.find(name)
            if elem is not None and elem.text:
                return elem.text.strip()
            # Try recursive search
            elem = root.find(f".//{name}")
            if elem is not None and elem.text:
                return elem.text.strip()
        return None

    def _parse_run_date(self, date_str: str) -> str:
        """Parse run date from various formats used by different Illumina instruments."""
        formats = [
            "%y%m%d",  # MiSeq/NextSeq500/HiSeq: 230415
            "%m/%d/%Y %I:%M:%S %p",  # NovaSeq6000: 4/15/2023 10:30:00 AM
            "%Y-%m-%dT%H:%M:%SZ",  # NextSeq2000: 2023-04-15T10:30:00Z
            "%Y-%m-%dT%H:%M:%S",  # Without Z
            "%Y%m%d",  # Some instruments: 20230415
        ]

        for fmt in formats:
            try:
                dt = datetime.strptime(date_str, fmt)
                return dt.strftime("%Y-%m-%d")
            except ValueError:
                continue

        return date_str

    def setup_general_stats_headers(self):
        """Configure headers for general stats table."""
        self.genstat_headers: Dict[str, ColumnDict] = {}

        self.genstat_headers["total_cycles"] = ColumnDict(
            {
                "title": "Cycles",
                "description": "Total sequencing cycles (excluding index reads)",
                "scale": "Blues",
                "format": "{:,.0f}",
                "hidden": False,
            }
        )

        self.genstat_headers["lane_count"] = ColumnDict(
            {
                "title": "Lanes",
                "description": "Number of lanes on the flowcell",
                "scale": "Greens",
                "format": "{:,.0f}",
                "hidden": True,
            }
        )

    def add_run_info_section(self):
        """Add a section showing run information."""
        runs_html = []
        for s_name, data in self.data_by_sample.items():
            # Build read configuration string
            read_info = ""
            if "reads" in data:
                for read in data["reads"]:
                    read_type = "Index" if read["is_indexed"] else "Read"
                    read_info += f"<li><b>{read_type} {read['number']}</b>: {read['cycles']} cycles</li>"

            instrument_info = f"""
                <li><b>Instrument:</b> {data.get('instrument', 'N/A')}</li>
                <li><b>Flowcell:</b> {data.get('flowcell', 'N/A')}</li>
                <li><b>Run Number:</b> {data.get('run_number', 'N/A')}</li>
                <li><b>Run Date:</b> {data.get('run_date', 'N/A')}</li>
            """

            experiment_info = ""
            if data.get("experiment_name"):
                experiment_info += f"<li><b>Experiment:</b> {data['experiment_name']}</li>"
            if data.get("chemistry"):
                experiment_info += f"<li><b>Chemistry:</b> {data['chemistry']}</li>"
            if data.get("rta_version"):
                experiment_info += f"<li><b>RTA Version:</b> {data['rta_version']}</li>"

            run_html = f"""
            <div class="card mb-3">
                <div class="card-header"><strong>{s_name}</strong></div>
                <div class="card-body">
                    <div class="row">
                        <div class="col-md-4">
                            <h5>Instrument</h5>
                            <ul class="list-unstyled">{instrument_info}</ul>
                        </div>
                        <div class="col-md-4">
                            <h5>Read Configuration</h5>
                            <ul>{read_info if read_info else '<li>N/A</li>'}</ul>
                        </div>
                        <div class="col-md-4">
                            <h5>Settings</h5>
                            <ul class="list-unstyled">{experiment_info if experiment_info else '<li>N/A</li>'}</ul>
                        </div>
                    </div>
                </div>
            </div>
            """
            runs_html.append(run_html)

        if runs_html:
            self.add_section(
                name="Run Information",
                anchor="sav-run-info",
                description="Basic run configuration from Illumina sequencer output.",
                content="".join(runs_html),
            )
