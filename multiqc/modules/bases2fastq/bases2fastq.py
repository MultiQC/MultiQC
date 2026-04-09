import copy
import json
import logging
import random
import re
import uuid
from collections import defaultdict
from itertools import chain
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from natsort import natsorted

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.modules.bases2fastq.plot_runs import (
    plot_base_quality_by_cycle,
    plot_base_quality_hist,
    plot_run_stats,
    tabulate_index_assignment_stats,
    tabulate_manifest_stats,
    tabulate_project_stats,
    tabulate_run_stats,
    tabulate_unassigned_index_stats,
)
from multiqc.modules.bases2fastq.plot_samples import (
    plot_adapter_content,
    plot_per_cycle_N_content,
    plot_per_read_gc_hist,
    sequence_content_plot,
    tabulate_sample_stats,
)
from multiqc.types import LoadedFileDict
from multiqc.utils import mqc_colour

log = logging.getLogger(__name__)

ELEMBIO_DOCS_URL = "https://docs.elembio.io/docs/bases2fastq/introduction/"

# Default minimum polony threshold - samples below this are skipped
DEFAULT_MIN_POLONIES = 1000


def _get_min_polonies() -> int:
    """
    Get the minimum polonies threshold from config or use default.

    Can be configured in multiqc_config.yaml:
        bases2fastq_config:
            min_polonies: 5000
    """
    cfg = getattr(config, "bases2fastq_config", {})
    if not isinstance(cfg, dict):
        return DEFAULT_MIN_POLONIES

    min_polonies = cfg.get("min_polonies", DEFAULT_MIN_POLONIES)
    try:
        min_polonies = int(min_polonies)
    except (ValueError, TypeError):
        log.warning(f"Invalid min_polonies value '{min_polonies}', using default {DEFAULT_MIN_POLONIES}")
        min_polonies = DEFAULT_MIN_POLONIES

    if min_polonies != DEFAULT_MIN_POLONIES:
        log.debug(f"Using custom min_polonies threshold: {min_polonies}")

    return min_polonies


class MultiqcModule(BaseMultiqcModule):
    """
    Bases2Fastq is Element Biosciences' secondary analysis software for demultiplexing
    sequencing data from AVITI systems and converting base calls into FASTQ files.

    Data Flow Overview
    ------------------
    The module handles three distinct data hierarchy levels:

    1. **Run Level**: Single sequencing run with all samples in one output
        - Directory: `<run_output>/`
        - Files: `RunStats.json`, `RunManifest.json`
        - Samples identified by: `{RunName}-{AnalysisID}__{SampleName}`

    2. **Project Level**: Demultiplexing by project, samples split into project subdirectories
        - Directory: `<run_output>/Samples/<ProjectName>/`
        - Files: Project-specific `RunStats.json`
        - Run-level `RunManifest.json` accessed via `../../RunManifest.json`
        - Samples identified by: `{RunName}-{AnalysisID}__{SampleName}`

    3. **Combined Level**: Both run and project data present (merged view)

    Parsing Flow
    ------------
    ```
    __init__()
        │
        ├─> _init_data_structures()     # Initialize empty dicts for all data levels
        │
        ├─> _parse_and_validate_data()  # Main parsing entry point
        │       │
        │       ├─> _parse_run_project_data("bases2fastq/run")     # Parse run-level RunStats.json
        │       │       └─> Populates: run_level_data, run_level_samples, run_level_samples_to_project
        │       │
        │       ├─> _parse_run_project_data("bases2fastq/project") # Parse project-level RunStats.json
        │       │       └─> Populates: project_level_data, project_level_samples, project_level_samples_to_project
        │       │
        │       └─> _determine_summary_path()  # Returns: "run_level" | "project_level" | "combined_level"
        │
        ├─> _select_data_by_summary_path()  # Route to appropriate data sources
        │       │
        │       ├─> _parse_run_manifest() or _parse_run_manifest_in_project()
        │       │       └─> Returns: manifest_data (lane settings, adapter info)
        │       │
        │       ├─> _parse_index_assignment() or _parse_index_assignment_in_project()
        │       │       └─> Returns: index_assignment_data (per-sample index stats)
        │       │
        │       └─> _parse_run_unassigned_sequences() (run_level only)
        │               └─> Returns: unassigned_sequences (unknown barcodes)
        │
        ├─> _setup_colors()             # Assign colors to runs/projects/samples
        │
        └─> _generate_plots()           # Create all report sections and plots
    ```

    Data Structures
    ---------------
    - `run_level_data`: Dict[run_name, run_stats] - Run-level QC metrics
    - `run_level_samples`: Dict[sample_id, sample_stats] - Sample metrics from run-level
    - `project_level_data`: Dict[project_name, project_stats] - Project-level QC metrics
    - `project_level_samples`: Dict[sample_id, sample_stats] - Sample metrics from project-level
    - `*_samples_to_project`: Dict[sample_id, project_name] - Maps samples to their projects

    Sample Naming Convention
    ------------------------
    Samples are uniquely identified as: `{RunName}-{AnalysisID[0:4]}__{SampleName}`
    This ensures uniqueness across multiple runs while keeping names readable.

    Files Parsed
    ------------
    - `RunStats.json`: Run/project QC metrics, sample statistics, lane data
    - `RunManifest.json`: Sample sheet info, index sequences, adapter settings

    Metrics Displayed
    -----------------
    - Polony counts and yields
    - Base quality distributions (histogram and by-cycle)
    - Index assignment statistics
    - Per-sample sequence content and GC distribution
    - Adapter content analysis
    - Unassigned/unknown barcode sequences (run-level only)
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Bases2Fastq",
            anchor="bases2fastq",
            href=ELEMBIO_DOCS_URL,
            info="Demultiplexes and converts Element AVITI base calls into FASTQ files",
            doi="10.1038/s41587-023-01750-7",
        )

        # Get configurable minimum polonies threshold
        self.min_polonies = _get_min_polonies()

        # Initialize data structures
        self._init_data_structures()

        # Parse and validate input data
        summary_path = self._parse_and_validate_data()

        # Select data based on summary path and parse additional sources
        run_data, sample_data, samples_to_projects, manifest_data, index_assignment_data, unassigned_sequences = (
            self._select_data_by_summary_path(summary_path)
        )

        # Set up color schemes for groups and samples
        self._setup_colors(sample_data, samples_to_projects, summary_path)

        # Generate all plots and sections
        self._generate_plots(
            summary_path,
            run_data,
            sample_data,
            samples_to_projects,
            manifest_data,
            index_assignment_data,
            unassigned_sequences,
        )

        # Write main data file at the very end after all sections are added
        self.write_data_file(sample_data, "bases2fastq")

    def _init_data_structures(self) -> None:
        """
        Initialize all data structures used by the module.

        Data structures are organized by hierarchy level:
        - Run level: Data from single-run Bases2Fastq output (no project splitting)
        - Project level: Data from project-split Bases2Fastq output
        - Combined: Merged data when both levels are present
        """
        # File cache to avoid reading the same JSON files multiple times
        # Key: resolved file path, Value: parsed JSON data
        self._file_cache: Dict[str, Any] = {}

        # === Run-level data structures ===
        # Populated from <run_output>/RunStats.json
        self.run_level_data: Dict[str, Any] = {}  # run_name -> full run stats
        self.run_level_samples: Dict[str, Any] = {}  # sample_id -> sample stats
        self.run_level_samples_to_project: Dict[str, str] = {}  # sample_id -> project name

        # === Project-level data structures ===
        # Populated from <run_output>/Samples/<project>/RunStats.json
        self.project_level_data: Dict[str, Any] = {}  # project_name -> project stats
        self.project_level_samples: Dict[str, Any] = {}  # sample_id -> sample stats
        self.project_level_samples_to_project: Dict[str, str] = {}  # sample_id -> project name

        # === Grouping structures for color assignment ===
        self.group_dict: Dict[str, Any] = {}  # group_name -> list of members
        self.group_lookup_dict: Dict[str, Any] = {}  # item -> group it belongs to
        self.project_lookup_dict: Dict[str, Any] = {}  # sample -> project mapping

    def _validate_path(self, file_path: Path, base_directory: Path) -> bool:
        """
        Validate that a file path doesn't escape outside the expected directory hierarchy.

        Args:
            file_path: Path to validate
            base_directory: The base directory that the path should stay within

        Returns:
            True if path is valid, False if it escapes the base directory
        """
        try:
            resolved_path = file_path.resolve()
            resolved_base = base_directory.resolve()
            # Check if the resolved path is within the base directory tree
            resolved_path.relative_to(resolved_base)
            return True
        except ValueError:
            # relative_to raises ValueError if path is not relative to base
            log.warning(
                f"Path {file_path} resolves outside expected directory {base_directory}. Skipping for security reasons."
            )
            return False

    def _read_json_file(self, file_path: Path, base_directory: Optional[Path] = None) -> Optional[Dict[str, Any]]:
        """
        Read and parse a JSON file with caching.

        Args:
            file_path: Path to the JSON file
            base_directory: Optional base directory to validate path against

        Returns:
            Parsed JSON data or None if reading failed
        """
        # Validate path doesn't escape expected directory if base is provided
        if base_directory is not None and not self._validate_path(file_path, base_directory):
            return None

        cache_key = str(file_path.resolve())

        if cache_key in self._file_cache:
            return self._file_cache[cache_key]

        if not file_path.exists():
            log.error(
                f"{file_path.name} does not exist at {file_path}.\n"
                f"Please visit Elembio online documentation for more information - "
                f"{ELEMBIO_DOCS_URL}"
            )
            return None

        try:
            with open(file_path) as _infile:
                data = json.load(_infile)
                self._file_cache[cache_key] = data
                return data
        except (json.JSONDecodeError, OSError) as e:
            log.error(f"Error reading {file_path}: {e}")
            return None

    def _parse_and_validate_data(self) -> str:
        """
        Parse input data and validate that samples were found.

        Returns:
            summary_path: The determined summary path ('run_level', 'project_level', or 'combined_level')
        """
        # Collect log files once per pattern (find_log_files returns a generator).
        # Stored as instance vars so downstream parsers can reuse them.
        self._run_level_log_files = list(self.find_log_files("bases2fastq/run"))
        self._project_level_log_files = list(self.find_log_files("bases2fastq/project"))
        run_level_log_files = self._run_level_log_files
        project_level_log_files = self._project_level_log_files

        if len(run_level_log_files) == 0 and len(project_level_log_files) == 0:
            error_msg = "No run- or project-level log files found within the Bases2Fastq results."
            log.error(error_msg)
            raise ModuleNoSamplesFound(error_msg)

        # Parse data from available sources
        if len(run_level_log_files) > 0:
            (self.run_level_data, self.run_level_samples, self.run_level_samples_to_project) = (
                self._parse_run_project_data("bases2fastq/run", log_files=run_level_log_files)
            )
        if len(project_level_log_files) > 0:
            (self.project_level_data, self.project_level_samples, self.project_level_samples_to_project) = (
                self._parse_run_project_data("bases2fastq/project", log_files=project_level_log_files)
            )

        # Count samples
        num_run_level_samples = len(self.run_level_samples)
        num_project_level_samples = len(self.project_level_samples)

        # Ensure at least some data was found
        if all(
            [
                len(self.run_level_data) == 0,
                num_run_level_samples == 0,
                len(self.project_level_data) == 0,
                num_project_level_samples == 0,
            ]
        ):
            error_msg = "No run-, project- or sample-level data found"
            log.error(error_msg)
            raise ModuleNoSamplesFound(error_msg)

        # Determine summary path
        summary_path = self._determine_summary_path()

        # Required call to confirm module is used (after confirming data was found)
        self.add_software_version(None)

        # Log what was found
        log.info(f"Found {len(self.run_level_data)} run(s) within the Bases2Fastq results.")
        log.info(f"Found {len(self.project_level_data)} project(s) within the Bases2Fastq results.")
        if summary_path == "run_level":
            log.info(f"Found {num_run_level_samples} sample(s) within the Bases2Fastq results.")
        else:
            log.info(f"Found {num_project_level_samples} sample(s) within the Bases2Fastq results.")

        # Warn if no data found
        if len(self.run_level_data) == 0 and len(self.project_level_data) == 0:
            log.warning("No run/project stats found!")
        if num_run_level_samples == 0 and num_project_level_samples == 0:
            log.warning("No sample stats found!")

        return summary_path

    def _determine_summary_path(self) -> str:
        """
        Determine which summary path to use based on available data.

        Returns:
            'run_level', 'project_level', or 'combined_level'
        """
        has_run_data = len(self.run_level_data) > 0
        has_project_data = len(self.project_level_data) > 0

        if has_run_data and not has_project_data:
            return "run_level"
        elif not has_run_data and has_project_data:
            return "project_level"
        elif has_run_data and has_project_data:
            return "combined_level"
        else:
            error_msg = "No run- or project-level data was retained. No report will be generated."
            log.error(error_msg)
            raise ModuleNoSamplesFound(error_msg)

    def _select_data_by_summary_path(
        self, summary_path: str
    ) -> Tuple[
        Dict[str, Any], Dict[str, Any], Dict[str, str], Dict[str, Any], Dict[str, Any], Dict[int, Dict[str, Any]]
    ]:
        """
        Select the appropriate data sources based on the summary path.

        Returns:
            Tuple of (run_data, sample_data, samples_to_projects, manifest_data,
                index_assignment_data, unassigned_sequences)
        """
        if summary_path == "run_level":
            manifest_log_files = list(self.find_log_files("bases2fastq/manifest"))
            return (
                self.run_level_data,
                self.run_level_samples,
                self.run_level_samples_to_project,
                self._parse_run_manifest("bases2fastq/manifest", log_files=manifest_log_files),
                self._parse_index_assignment("bases2fastq/manifest", log_files=manifest_log_files),
                self._parse_run_unassigned_sequences("bases2fastq/run", log_files=self._run_level_log_files),
            )
        elif summary_path == "project_level":
            return (
                self.project_level_data,
                self.project_level_samples,
                self.project_level_samples_to_project,
                self._parse_run_manifest_in_project("bases2fastq/project", log_files=self._project_level_log_files),
                self._parse_index_assignment_in_project("bases2fastq/project", log_files=self._project_level_log_files),
                {},  # No unassigned sequences for project level
            )
        elif summary_path == "combined_level":
            # Use run-level stats for the run table (more complete), but
            # project-level samples for per-sample plots (properly split by project).
            manifest_log_files = list(self.find_log_files("bases2fastq/manifest"))
            return (
                self.run_level_data,
                self.project_level_samples,
                self.project_level_samples_to_project,
                self._parse_run_manifest("bases2fastq/manifest", log_files=manifest_log_files),
                self._parse_index_assignment("bases2fastq/manifest", log_files=manifest_log_files),
                self._parse_run_unassigned_sequences("bases2fastq/run", log_files=self._run_level_log_files),
            )
        else:
            error_msg = "No run- or project-level data was retained. No report will be generated."
            log.error(error_msg)
            raise ModuleNoSamplesFound(error_msg)

    def _setup_colors(
        self, sample_data: Dict[str, Any], samples_to_projects: Dict[str, str], summary_path: str
    ) -> None:
        """Set up color schemes for groups and samples."""
        # Create run and project groups
        run_groups: Dict[str, List] = defaultdict(list)
        project_groups: Dict[str, List] = defaultdict(list)
        ind_sample_groups: Dict[str, List] = defaultdict(list)

        for sample in natsorted(sample_data.keys()):
            run_name, _ = sample.split("__", maxsplit=1)
            run_groups[run_name].append(sample)
            sample_project = samples_to_projects.get(sample, "DefaultProject")
            project_groups[sample_project].append(sample)
            ind_sample_groups[sample] = [sample]

        merged_groups = {**run_groups, **project_groups, **ind_sample_groups}

        # Build color palette
        self.color_getter = mqc_colour.mqc_colour_scale()
        self.palette = list(
            chain.from_iterable(
                self.color_getter.get_colours(hue)
                for hue in ["Set2", "Pastel1", "Accent", "Set1", "Set3", "Dark2", "Paired", "Pastel2"]
            )
        )

        # Add extra colors if needed
        if len(merged_groups) > len(self.palette):
            extra_colors = [
                f"#{random.randrange(0, 0xFFFFFF):06x}" for _ in range(len(self.palette), len(merged_groups))
            ]
            self.palette = self.palette + extra_colors

        # Assign colors to groups
        self.group_color = {
            group: color for group, color in zip(merged_groups.keys(), self.palette[: len(merged_groups)])
        }

        # Assign colors to samples
        self.sample_color: Dict[str, str] = {}
        for sample_name in natsorted(samples_to_projects.keys()):
            if summary_path == "project_level" or len(project_groups) == 1:
                sample_color = self.group_color[sample_name]
            else:
                sample_color = self.group_color[samples_to_projects[sample_name]]
            self.sample_color[sample_name] = sample_color

        # Copy group colors to run colors
        self.run_color = copy.deepcopy(self.group_color)
        self.palette = self.palette[len(merged_groups) :]

    def _generate_plots(
        self,
        summary_path: str,
        run_data: Dict[str, Any],
        sample_data: Dict[str, Any],
        samples_to_projects: Dict[str, str],
        manifest_data: Dict[str, Any],
        index_assignment_data: Dict[str, Any],
        unassigned_sequences: Dict[int, Dict[str, Any]],
    ) -> None:
        """Generate all plots and add sections to the report."""
        # QC metrics table
        qc_metrics_function = (
            tabulate_run_stats if summary_path in ["run_level", "combined_level"] else tabulate_project_stats
        )
        self.add_run_plots(data=run_data, plot_functions=[qc_metrics_function])

        # Manifest stats
        self.add_run_plots(data=manifest_data, plot_functions=[tabulate_manifest_stats])

        # Index assignment stats
        self.add_run_plots(data=index_assignment_data, plot_functions=[tabulate_index_assignment_stats])

        # Unassigned sequences (only for run_level and combined_level)
        if summary_path in ["run_level", "combined_level"]:
            self.add_run_plots(data=unassigned_sequences, plot_functions=[tabulate_unassigned_index_stats])

        # Run-level plots
        self.add_run_plots(
            data=run_data,
            plot_functions=[plot_run_stats, plot_base_quality_hist, plot_base_quality_by_cycle],
        )

        # Sample-level plots
        self.add_sample_plots(
            data=sample_data,
            group_lookup=samples_to_projects,
            project_lookup=samples_to_projects,
        )

    def get_uuid(self) -> str:
        return str(uuid.uuid4()).replace("-", "").lower()

    def _extract_run_analysis_name(
        self,
        data: Dict[str, Any],
        source_info: str = "RunStats.json",
    ) -> Optional[str]:
        """
        Extract and validate run_analysis_name from data dict.

        Args:
            data: Dictionary containing RunName and AnalysisID keys
            source_info: Description of the data source for error messages

        Returns:
            The run_analysis_name (RunName-AnalysisID[0:4]) or None if extraction failed
        """
        run_name = data.get("RunName")
        analysis_id = data.get("AnalysisID")

        if not run_name or not analysis_id:
            log.error(
                f"Error with {source_info}. Either RunName or AnalysisID is absent.\n"
                f"RunName: {run_name}, AnalysisID: {analysis_id}\n"
                f"Please visit Elembio online documentation for more information - "
                f"{ELEMBIO_DOCS_URL}"
            )
            return None

        return f"{run_name}-{analysis_id[0:4]}"

    def _parse_run_project_data(
        self, data_source: str, log_files: Optional[List[LoadedFileDict[Any]]] = None
    ) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, str]]:
        """
        Parse RunStats.json files to extract run/project and sample-level data.

        This is the primary parsing method that populates the core data structures.
        It handles both run-level and project-level RunStats.json files.

        Args:
            data_source: Search pattern key ("bases2fastq/run" or "bases2fastq/project")
            log_files: Optional pre-collected list of file dicts from find_log_files.
                When provided, used instead of calling find_log_files again.

        Returns:
            Tuple of:
            - runs_global_data: Dict[run_name, run_stats] - Run/project level metrics
            - runs_sample_data: Dict[sample_id, sample_stats] - Per-sample metrics
            - sample_to_project: Dict[sample_id, project_name] - Sample-to-project mapping

        Data Flow:
            RunStats.json -> parse -> filter samples by min_polonies -> populate dicts
        """
        runs_global_data: Dict[str, Any] = {}
        runs_sample_data: Dict[str, Any] = {}
        sample_to_project: Dict[str, str] = {}
        if data_source == "":
            return (runs_global_data, runs_sample_data, sample_to_project)

        files_to_process = log_files if log_files is not None else list(self.find_log_files(data_source))
        for f in files_to_process:
            data = json.loads(f["f"])

            # Copy incomind data and reset samples to include only desired
            data_to_return = copy.deepcopy(data)
            data_to_return["SampleStats"] = []

            # get run + analysis
            run_name = data.get("RunName")
            run_analysis_name = self._extract_run_analysis_name(data, source_info=f"RunStats.json ({f['fn']})")
            if run_analysis_name is None:
                continue
            run_analysis_name = self.clean_s_name(run_analysis_name, f)

            # skip run if in user provider ignore list
            if self.is_ignore_sample(run_analysis_name):
                log.info(f"Skipping <{run_analysis_name}> because it is present in ignore list.")
                continue

            # Check run is present in the final dictionaries
            if run_analysis_name not in runs_global_data:
                runs_global_data[run_analysis_name] = data_to_return

            project = self.clean_s_name(data.get("Project", "DefaultProject"), f)

            # map sample UUIDs to run_analysis_name
            for sample_data in data["SampleStats"]:
                sample_id = sample_data["SampleID"]
                sample_name = sample_data["SampleName"]
                sample_data["RunName"] = run_name
                run_analysis_sample_name = "__".join([run_analysis_name, sample_name])

                num_polonies = sample_data["NumPolonies"]
                if num_polonies < self.min_polonies:
                    log.warning(
                        f"Skipping {run_analysis_sample_name} because it has "
                        f"<{self.min_polonies} assigned reads [n={num_polonies}]."
                    )
                    continue

                # skip run if in user provider ignore list
                if self.is_ignore_sample(sample_id) or self.is_ignore_sample(run_analysis_sample_name):
                    log.info(
                        f"Skipping <{sample_id}> ({run_analysis_sample_name}) because it is present in ignore list."
                    )
                    continue

                # If sample passes all checks add it back
                runs_sample_data[run_analysis_sample_name] = sample_data
                sample_to_project[run_analysis_sample_name] = project

            self.add_data_source(f=f, s_name=run_analysis_name, module="bases2fastq")

        return (runs_global_data, runs_sample_data, sample_to_project)

    def _extract_manifest_lane_settings(
        self, run_manifest_data: Dict[str, Any], run_analysis_name: str
    ) -> Dict[str, Dict[str, Any]]:
        """
        Extract per-lane settings from a parsed RunManifest.json Settings section.

        Args:
            run_manifest_data: Parsed RunManifest.json (must contain "Settings" list)
            run_analysis_name: Run identifier for building run_lane keys

        Returns:
            Dict[run_lane, settings] where run_lane = "{run_analysis_name} | L{lane_id}"
            and settings contain Indexing, AdapterTrimType, R1/R2AdapterMinimumTrimmedLength
        """
        result: Dict[str, Dict[str, Any]] = {}
        if "Settings" not in run_manifest_data:
            return result
        for lane_data in run_manifest_data["Settings"]:
            lane_id = lane_data.get("Lane")
            if not lane_id:
                log.error("<Lane> not found in Settings section of RunManifest. Skipping lanes.")
                continue
            lane_name = f"L{lane_id}"
            run_lane = f"{run_analysis_name} | {lane_name}"
            result[run_lane] = {}

            indices = []
            indices_cycles = []
            mask_pattern = re.compile(r"^I\d+Mask$")
            matching_keys = [key for key in lane_data.keys() if mask_pattern.match(key)]
            for key in matching_keys:
                for mask_info in lane_data[key]:
                    if mask_info["Read"] not in indices:
                        indices.append(mask_info["Read"])
                    indices_cycles.append(str(len(mask_info["Cycles"])))
            indexing = f"{' + '.join(indices_cycles)}<br>{' + '.join(indices)}"
            result[run_lane]["Indexing"] = indexing
            result[run_lane]["AdapterTrimType"] = lane_data.get("AdapterTrimType", "N/A")
            result[run_lane]["R1AdapterMinimumTrimmedLength"] = lane_data.get("R1AdapterMinimumTrimmedLength", "N/A")
            result[run_lane]["R2AdapterMinimumTrimmedLength"] = lane_data.get("R2AdapterMinimumTrimmedLength", "N/A")
        return result

    def _parse_run_manifest(
        self, data_source: str, log_files: Optional[List[LoadedFileDict[Any]]] = None
    ) -> Dict[str, Any]:
        """
        Parse RunManifest.json for run-level analysis to extract lane and adapter settings.

        Data Flow:
            RunManifest.json (via data_source pattern)
            + RunStats.json (for run name) from same directory
            -> Extract per-lane: index masks, adapter settings, trim lengths

        Args:
            data_source: Search pattern key for RunManifest.json files
            log_files: Optional pre-collected list of file dicts from find_log_files.

        Returns:
            Dict[run_lane, settings] where run_lane = "{run_name} | L{lane_id}"
        """
        runs_manifest_data: Dict[str, Dict[str, Any]] = {}

        if data_source == "":
            return runs_manifest_data

        files_to_process = log_files if log_files is not None else list(self.find_log_files(data_source))
        for f in files_to_process:
            directory = f.get("root")
            if not directory:
                continue

            # Get RunName and RunID from RunStats.json
            run_stats_path = Path(directory) / "RunStats.json"
            run_stats = self._read_json_file(run_stats_path)
            if run_stats is None:
                continue

            run_analysis_name = self._extract_run_analysis_name(run_stats, source_info=str(run_stats_path))
            if run_analysis_name is None:
                continue

            run_manifest = json.loads(f["f"])
            if "Settings" not in run_manifest:
                log.warning(
                    f"<Settings> section not found in {directory}/RunManifest.json.\nSkipping RunManifest metrics."
                )
            else:
                runs_manifest_data.update(self._extract_manifest_lane_settings(run_manifest, run_analysis_name))

            self.add_data_source(f=f, s_name=run_analysis_name, module="bases2fastq")

        return runs_manifest_data

    def _parse_run_manifest_in_project(
        self, data_source: str, log_files: Optional[List[LoadedFileDict[Any]]] = None
    ) -> Dict[str, Any]:
        """
        Parse RunManifest.json for project-level analysis.

        Similar to _parse_run_manifest but navigates up from project directories
        to find the run-level RunManifest.json (via ../../RunManifest.json).

        Data Flow:
            Project RunStats.json (for run name)
            + ../../RunManifest.json (run-level manifest)
            -> Extract per-lane settings
        """
        project_manifest_data: Dict[str, Dict[str, Any]] = {}

        if data_source == "":
            return project_manifest_data

        files_to_process = log_files if log_files is not None else list(self.find_log_files(data_source))
        for f in files_to_process:
            directory = f.get("root")
            if not directory:
                continue

            # Resolve base_directory to the run output root (not the project subdirectory),
            # since RunManifest.json lives at the run root. Path validation in _read_json_file
            # will check the manifest path against this run root directory.
            base_directory = Path(directory).resolve()
            if not (base_directory / "RunManifest.json").exists():
                base_directory = base_directory.parent.parent
            run_manifest = base_directory / "RunManifest.json"
            project_stats = json.loads(f["f"])
            run_analysis_name = self._extract_run_analysis_name(
                project_stats, source_info=f"project RunStats.json ({f['fn']})"
            )
            if run_analysis_name is None:
                continue

            # skip run if in user provider ignore list
            if self.is_ignore_sample(run_analysis_name):
                log.info(f"Skipping <{run_analysis_name}> because it is present in ignore list.")
                continue

            run_manifest_data = self._read_json_file(run_manifest, base_directory=base_directory)
            if run_manifest_data is None:
                continue

            if "Settings" not in run_manifest_data:
                log.warning(f"<Settings> section not found in {run_manifest}.\nSkipping RunManifest metrics.")
            else:
                project_manifest_data.update(self._extract_manifest_lane_settings(run_manifest_data, run_analysis_name))
            data_source_info: LoadedFileDict[Any] = {
                "fn": str(run_manifest.name),
                "root": str(run_manifest.parent),
                "sp_key": data_source,
                "s_name": str(run_manifest.with_suffix("").name),
                "f": run_manifest_data,
            }
            self.add_data_source(f=data_source_info, s_name=run_analysis_name, module="bases2fastq")

        return project_manifest_data

    def _build_index_assignment_from_stats(
        self,
        stats_dict: Dict[str, Any],
        run_analysis_name: str,
        project: Optional[str] = None,
    ) -> Tuple[Dict[str, Dict[str, Any]], int]:
        """
        Build per-run index assignment dict from RunStats SampleStats/Occurrences.

        Returns:
            Tuple of (run_inner_dict, total_polonies). run_inner_dict is
            { merged_expected_sequence -> { SampleID, SamplePolonyCounts, PercentOfPolonies, Index1, Index2, ... } }
        """
        run_inner: Dict[str, Dict[str, Any]] = {}
        total_polonies = stats_dict.get("NumPoloniesBeforeTrimming", 0)
        if "SampleStats" not in stats_dict:
            return (run_inner, total_polonies)
        for sample_data in stats_dict["SampleStats"]:
            sample_name = sample_data.get("SampleName")
            sample_id = "__".join([run_analysis_name, sample_name]) if (run_analysis_name and sample_name) else None
            if "Occurrences" not in sample_data:
                log.error(f"Missing data needed to extract index assignment for sample {sample_id}. Skipping.")
                continue
            for occurrence in sample_data["Occurrences"]:
                sample_expected_seq = occurrence.get("ExpectedSequence")
                sample_counts = occurrence.get("NumPoloniesBeforeTrimming")
                if any(x is None for x in [sample_expected_seq, sample_counts, sample_id]):
                    log.error(f"Missing data needed to extract index assignment for sample {sample_id}. Skipping.")
                    continue
                if sample_expected_seq not in run_inner:
                    entry: Dict[str, Any] = {
                        "SampleID": sample_id,
                        "SamplePolonyCounts": 0,
                        "PercentOfPolonies": float("nan"),
                        "Index1": "",
                        "Index2": "",
                    }
                    if project is not None:
                        entry["Project"] = project
                    run_inner[sample_expected_seq] = entry
                run_inner[sample_expected_seq]["SamplePolonyCounts"] += sample_counts
        for entry in run_inner.values():
            if total_polonies > 0:
                entry["PercentOfPolonies"] = round(entry["SamplePolonyCounts"] / total_polonies * 100, 2)
        return (run_inner, total_polonies)

    def _merge_manifest_index_sequences(
        self,
        sample_to_index_assignment: Dict[str, Any],
        run_manifest_data: Dict[str, Any],
        run_analysis_name: str,
    ) -> None:
        """Merge Index1/Index2 from RunManifest Samples into sample_to_index_assignment (mutates)."""
        if "Samples" not in run_manifest_data or run_analysis_name not in sample_to_index_assignment:
            return
        run_data = sample_to_index_assignment[run_analysis_name]
        for sample_data in run_manifest_data["Samples"]:
            sample_name = sample_data.get("SampleName")
            if run_analysis_name is None or sample_name is None or "Indexes" not in sample_data:
                continue
            sample_id = "__".join([run_analysis_name, sample_name])
            for index_data in sample_data["Indexes"]:
                index_1 = index_data.get("Index1", "")
                index_2 = index_data.get("Index2", "")
                merged_indices = f"{index_1}{index_2}"
                if merged_indices not in run_data:
                    log.error(f"Index assignment information not found for sample {sample_id}. Skipping.")
                    continue
                if sample_id != run_data[merged_indices]["SampleID"]:
                    log.error(
                        f"RunManifest SampleID <{sample_id}> does not match "
                        f"RunStats SampleID {run_data[merged_indices]['SampleID']}. Skipping."
                    )
                    continue
                run_data[merged_indices]["Index1"] = index_1
                run_data[merged_indices]["Index2"] = index_2

    def _parse_run_unassigned_sequences(
        self, data_source: str, log_files: Optional[List[LoadedFileDict[Any]]] = None
    ) -> Dict[int, Dict[str, Any]]:
        """
        Parse unassigned/unknown barcode sequences from run-level data.

        Only available for run-level analysis. Extracts sequences that could not
        be assigned to any sample, useful for troubleshooting index issues.

        Data Flow:
            RunStats.json -> Lanes -> UnassignedSequences
            -> Extract: sequence, count, percentage of total polonies
        """
        run_unassigned_sequences: Dict[int, Dict[str, Any]] = {}
        if data_source == "":
            return run_unassigned_sequences

        files_to_process = log_files if log_files is not None else list(self.find_log_files(data_source))
        for f in files_to_process:
            data = json.loads(f["f"])

            # Get RunName and AnalysisID
            run_analysis_name = self._extract_run_analysis_name(data, source_info=f"RunStats.json ({f['fn']})")
            if run_analysis_name is None:
                continue
            run_analysis_name = self.clean_s_name(run_analysis_name, f)

            # skip run if in user provider ignore list
            if self.is_ignore_sample(run_analysis_name):
                log.info(f"Skipping <{run_analysis_name}> because it is present in ignore list.")
                continue

            # Get total polonies and build unassigned indices dictionary
            total_polonies = data.get("NumPoloniesBeforeTrimming", 0)
            if "Lanes" not in data:
                log.error(
                    f"Missing lane information in RunStats.json for run {run_analysis_name}."
                    f"Skipping building unassigned indices table."
                )
                continue
            index_number = 1
            for lane in data["Lanes"]:
                lane_id = lane.get("Lane")
                if lane_id:
                    lane_id = f"L{lane_id}"
                for sequence in lane.get("UnassignedSequences", []):
                    run_unassigned_sequences[index_number] = {
                        "Run Name": run_analysis_name,
                        "Lane": lane_id,
                        "I1": sequence["I1"],
                        "I2": sequence["I2"],
                        "Number of Polonies": sequence["Count"],
                        "% Polonies": float("nan"),
                    }
                    if total_polonies > 0:
                        run_unassigned_sequences[index_number]["% Polonies"] = round(
                            sequence["Count"] / total_polonies, 2
                        )
                    index_number += 1

        return run_unassigned_sequences

    def _parse_index_assignment(
        self, manifest_data_source: str, log_files: Optional[List[LoadedFileDict[Any]]] = None
    ) -> Dict[str, Any]:
        """
        Parse index assignment statistics for run-level analysis.

        Combines data from RunStats.json (polony counts) and RunManifest.json
        (index sequences) to show how well each sample's index performed.

        Data Flow:
            RunStats.json -> SampleStats -> per-sample polony counts
            + RunManifest.json -> Samples -> index sequences (Index1, Index2)
            -> Combined index assignment table
        """
        sample_to_index_assignment: Dict[str, Dict[str, Dict[str, Any]]] = {}

        if manifest_data_source == "":
            return sample_to_index_assignment

        files_to_process = log_files if log_files is not None else list(self.find_log_files(manifest_data_source))
        for f in files_to_process:
            directory = f.get("root")
            if not directory:
                continue

            # Get RunName and RunID from RunStats.json
            run_stats_path = Path(directory) / "RunStats.json"
            run_stats = self._read_json_file(run_stats_path)
            if run_stats is None:
                continue

            run_analysis_name = self._extract_run_analysis_name(run_stats, source_info=str(run_stats_path))
            if run_analysis_name is None:
                continue

            # skip run if in user provider ignore list
            if self.is_ignore_sample(run_analysis_name):
                log.info(f"Skipping <{run_analysis_name}> because it is present in ignore list.")
                continue

            if "SampleStats" not in run_stats:
                log.error(
                    f"Error, missing SampleStats in RunStats.json. Skipping index assignment metrics.\n"
                    f"Available keys: {list(run_stats.keys())}\n"
                    f"Please visit Elembio online documentation for more information - "
                    f"{ELEMBIO_DOCS_URL}"
                )
                continue

            run_inner, _ = self._build_index_assignment_from_stats(run_stats, run_analysis_name)
            sample_to_index_assignment[run_analysis_name] = run_inner

            run_manifest = json.loads(f["f"])
            if "Samples" not in run_manifest:
                log.warning(
                    f"<Samples> section not found in {directory}/RunManifest.json.\n"
                    f"Skipping RunManifest sample index assignment metrics."
                )
            elif len(sample_to_index_assignment) == 0:
                log.warning("Index assignment data missing. Skipping creation of index assignment metrics.")
            else:
                self._merge_manifest_index_sequences(sample_to_index_assignment, run_manifest, run_analysis_name)

        return sample_to_index_assignment

    def _parse_index_assignment_in_project(
        self, data_source: str, log_files: Optional[List[LoadedFileDict[Any]]] = None
    ) -> Dict[str, Any]:
        """
        Parse index assignment statistics for project-level analysis.

        Similar to _parse_index_assignment but works with project-split output,
        navigating up to find the run-level RunManifest.json.

        Data Flow:
            Project RunStats.json -> SampleStats -> polony counts
            + ../../RunManifest.json -> Samples -> index sequences
            -> Combined index assignment table
        """
        sample_to_index_assignment: Dict[str, Dict[str, Dict[str, Any]]] = {}

        if data_source == "":
            return sample_to_index_assignment

        files_to_process = log_files if log_files is not None else list(self.find_log_files(data_source))
        for f in files_to_process:
            directory = f.get("root")
            if not directory:
                continue

            # Get RunManifest.json from run output root (two levels up from project directory)
            base_directory = Path(directory).parent.parent
            run_manifest = base_directory / "RunManifest.json"

            project_stats = json.loads(f["f"])
            project = self.clean_s_name(project_stats.get("Project", "DefaultProject"), f)

            run_analysis_name = self._extract_run_analysis_name(
                project_stats, source_info=f"project RunStats.json ({f['fn']})"
            )
            if run_analysis_name is None:
                continue

            # skip run if in user provider ignore list
            if self.is_ignore_sample(run_analysis_name):
                log.info(f"Skipping <{run_analysis_name}> because it is present in ignore list.")
                continue

            if "SampleStats" not in project_stats:
                log.error(
                    f"Error, missing SampleStats in RunStats.json. Skipping index assignment metrics.\n"
                    f"Available keys: {list(project_stats.keys())}\n"
                    f"Please visit Elembio online documentation for more information - "
                    f"{ELEMBIO_DOCS_URL}"
                )
                continue

            run_inner, _ = self._build_index_assignment_from_stats(project_stats, run_analysis_name, project=project)
            sample_to_index_assignment[run_analysis_name] = run_inner

            run_manifest_data = self._read_json_file(run_manifest, base_directory=base_directory)
            if run_manifest_data is None:
                continue

            if "Samples" not in run_manifest_data:
                log.warning(
                    f"<Samples> section not found in {run_manifest}.\n"
                    f"Skipping RunManifest sample index assignment metrics."
                )
            elif len(sample_to_index_assignment) == 0:
                log.warning("Index assignment data missing. Skipping creation of index assignment metrics.")
            else:
                self._merge_manifest_index_sequences(sample_to_index_assignment, run_manifest_data, run_analysis_name)

        return sample_to_index_assignment

    def add_run_plots(self, data: Dict[Any, Any], plot_functions: List[Callable]) -> None:
        if not data:
            return
        for func in plot_functions:
            plot_html, plot_name, anchor, description, helptext, plot_data = func(data, self.run_color)
            if plot_html is not None:
                self.add_section(
                    name=plot_name, plot=plot_html, anchor=anchor, description=description, helptext=helptext
                )
                self.write_data_file(plot_data, f"base2fastq:{plot_name}")

    def add_sample_plots(
        self, data: Dict[str, Any], group_lookup: Dict[str, str], project_lookup: Dict[str, str]
    ) -> None:
        if not data:
            return
        plot_functions: List[Callable] = [
            tabulate_sample_stats,
            sequence_content_plot,
            plot_per_cycle_N_content,
            plot_adapter_content,
            plot_per_read_gc_hist,
        ]
        for func in plot_functions:
            plot_html, plot_name, anchor, description, helptext, plot_data = func(
                data, group_lookup, project_lookup, self.sample_color
            )
            if plot_html is not None:
                self.add_section(
                    name=plot_name, plot=plot_html, anchor=anchor, description=description, helptext=helptext
                )
                self.write_data_file(plot_data, f"base2fastq:{plot_name}")
