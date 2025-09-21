from collections import defaultdict
import copy
import csv
import json
import logging
import random
from typing import Any, Dict, List
import uuid

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.utils import mqc_colour

from multiqc.modules.bases2fastq.plot_runs import (
    plot_run_stats,
    tabulate_run_stats,
    tabulate_project_stats,
    plot_base_quality_hist,
    plot_base_quality_by_cycle,
    plot_lane_cycle_stats,
)
from multiqc.modules.bases2fastq.plot_samples import (
    tabulate_sample_stats,
    sequence_content_plot,
    plot_per_cycle_N_content,
    plot_adapter_content,
    plot_per_read_gc_hist,
    plot_sample_read_length,
)

log = logging.getLogger(__name__)


MIN_POLONIES = 10000


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Bases2Fastq",
            anchor="bases2fastq",
            href="https://docs.elembio.io/docs/bases2fastq/introduction/",
            info="Demultiplexes and converts Element AVITI base calls into FASTQ files",
            doi="10.1038/s41587-023-01750-7",
        )

        # Initialize run, project and sample level structures
        self.run_level_data = {}
        self.run_level_samples = {}
        self.run_level_samples_to_project = {}
        self.project_level_data = {}
        self.project_level_samples = {}
        self.project_level_samples_to_project = {}
        num_run_level_samples = 0
        num_project_level_samples = 0

        # Initialize run and project groups
        self.group_dict = dict()
        self.group_lookup_dict = dict()
        self.project_lookup_dict = dict()


        self.b2f_sample_data = dict()
        self.b2f_run_data = dict()
        self.b2f_run_project_data = dict()
        self.b2f_run_project_sample_data = dict()
        self.missing_runs = set()
        self.sample_id_to_run = dict()

        # Define if call is project- or run-level
        run_level_log_files = len(list(self.find_log_files("bases2fastq/run")))
        project_level_log_files = len(list(self.find_log_files("bases2fastq/project")))
        
        if run_level_log_files == 0 and project_level_log_files == 0:
            error_msg = "No run- or project-level log files found within the Bases2Fastq results."
            log.error(error_msg)
            raise ModuleNoSamplesFound(error_msg)
        
        # Parse data
        if run_level_log_files > 0:
            (
                self.run_level_data, self.run_level_samples, self.run_level_samples_to_project
            ) = self._parse_run_project_data("bases2fastq/run")
        if project_level_log_files > 0:
            (
                self.project_level_data, self.project_level_samples, self.project_level_samples_to_project
            ) = self._parse_run_project_data("bases2fastq/project")

        # Get run- and project-level samples
        for data in self.run_level_samples.values():
            num_run_level_samples += len(data.keys())
        for data in self.project_level_samples.values():
            num_project_level_samples += len(data.keys())

        # Ensure run/sample data found
        if all([
            len(self.run_level_data) == 0,
            num_run_level_samples == 0,
            len(self.project_level_data),
            num_project_level_samples == 0,
        ]):
            error_msg = "No run-, project- or sample-level data found"
            log.error(error_msg)
            raise ModuleNoSamplesFound(error_msg)
        
        # Log runs, projects and samples found
        log.info(f"Found {len(self.run_level_data)} run(s) within the Bases2Fastq results.")
        log.info(f"Found {len(self.project_level_data)} project(s) within the Bases2Fastq results.")
        log.info(f"Found {num_project_level_samples} sample(s) within the Bases2Fastq results.")

        # Superfluous function call to confirm that it is used in this module
        self.add_software_version(None)

        # Warn user if run-level/project-level or sample-level metrics were not found
        if len(self.run_level_data) == 0 and len(self.project_level_data) == 0:
            log.warning("No run/project stats found!")
        if num_project_level_samples == 0:
            log.warning("No sample stats found!")
        
        # Choose path to take, if project use only project-level data, otherwise use run- and project-level
        summary_path = ""
        if len(self.run_level_data) > 0 and len(self.project_level_data) == 0:
            summary_path = "run_level"
        if len(self.run_level_data) == 0 and len(self.project_level_data) > 0:
            summary_path = "project_level"
        elif len(self.run_level_data) > 0 and len(self.project_level_data) > 0:
            summary_path = "combined_level"
        
        # Define data to use
        run_data = {}
        sample_data = {}
        samples_to_projects = {}
        if summary_path == "run_level":
            run_data = self.run_level_data
            sample_data = self.project_level_samples
            samples_to_projects = self.run_level_samples_to_project
        elif summary_path == "project_level":
            run_data = self.project_level_data
            sample_data = self.project_level_samples
            samples_to_projects = self.project_level_samples_to_project
        elif summary_path == "combined_level":
            run_data = self.run_level_data
            sample_data = self.project_level_samples
            samples_to_projects = self.project_level_samples_to_project
        else:
            error_msg = "No run- or project-level data was retained. No report will be generated."
            log.error(error_msg)
            return

        # Create run and project groups
        run_groups = defaultdict(list)
        project_groups = defaultdict(list)
        sample_to_run_group = {}
        for sample in sample_data.keys():
            (_run_name, _) = sample.split("__")
            run_groups[_run_name].append(sample)
            sample_to_run_group[sample] = _run_name
            sample_project = samples_to_projects[sample]
            project_groups[sample_project].append(sample)
        merged_groups = dict(run_groups) | dict(project_groups)

        # Assign color for each group
        self.color_getter = mqc_colour.mqc_colour_scale()
        self.palette = sum(
            [
                self.color_getter.get_colours(hue)
                for hue in ["Set2", "Pastel1", "Accent", "Set1", "Set3", "Dark2", "Paired", "Pastel2"]
            ],
            [],
        )
        if len(merged_groups) > len(self.palette):
            hex_range = 2**24
            extra_colors = [hex(random.randrange(0, hex_range)) for _ in range(len(merged_groups), len(self.palette))]
            self.palette = self.palette + extra_colors
        self.group_color = {g: c for g, c in zip(merged_groups.keys(), self.palette[: len(merged_groups)])}
        self.sample_color = dict()
        for s_name in samples_to_projects.keys():
            self.sample_color.update({s_name: self.group_color[samples_to_projects[s_name]]})
        self.run_color = copy.deepcopy(self.group_color)  # Make sure that run colors and group colors match
        self.palette = self.palette[len(merged_groups) :]


        # Plot metrics
        qc_metrics_function = (
            tabulate_run_stats if summary_path in ["run_level", "combined_level"] else tabulate_project_stats
        )
        self.add_run_plots(
            data=run_data,
            plot_functions=[
                qc_metrics_function,
                plot_lane_cycle_stats,
                plot_run_stats,
                plot_base_quality_hist,
                plot_base_quality_by_cycle
            ]
        )
        self.add_sample_plots(
            data=sample_data, group_lookup=samples_to_projects, project_lookup=samples_to_projects
        )

    def get_uuid(self):
        return str(uuid.uuid4()).replace("-", "").lower()

    def _parse_run_project_data(self, data_source: str) -> List[Dict[str, Any]]:
        runs_global_data = {}
        runs_sample_data = {}
        sample_to_project = {}
        if data_source == "":
            return [runs_global_data, runs_sample_data, sample_to_project]

        for f in self.find_log_files(data_source):
            data = json.loads(f["f"])

            # Copy incomind data and reset samples to include only desired
            data_to_return = copy.deepcopy(data)
            data_to_return["SampleStats"] = []

            # get run + analysis
            run_name = data.get("RunName", None)
            analysis_id = data.get("AnalysisID", None)[0:4]

            if not run_name or not analysis_id:
                log.error(
                    "Error with RunStats.json. Either RunName or AnalysisID is absent.\n"
                    "Please visit Elembio online documentation for more information - "
                    "https://docs.elembio.io/docs/bases2fastq/introduction/"
                )
                continue
        
            run_analysis_name = "-".join([run_name, analysis_id])
            run_analysis_name = self.clean_s_name(run_analysis_name, f)

            # skip run if in user provider ignore list
            if self.is_ignore_sample(run_analysis_name):
                log.info(
                    f"Skipping <{run_analysis_name}> because it is present in ignore list."
                )
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
                if num_polonies < 1000:
                    log.warning(
                        f"Skipping {run_analysis_sample_name} because it has"
                        f" <{MIN_POLONIES} assigned reads [n={num_polonies}]."
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

        return [runs_global_data, runs_sample_data, sample_to_project]

    def add_run_plots(self, data, plot_functions):
        for func in plot_functions:
            plot_html, plot_name, anchor, description, helptext, plot_data = func(data, self.run_color)
            self.add_section(name=plot_name, plot=plot_html, anchor=anchor, description=description, helptext=helptext)
            self.write_data_file(plot_data, f"base2fastq:{plot_name}")

    def add_sample_plots(self, data, group_lookup, project_lookup):
        plot_functions = [
            tabulate_sample_stats,
            plot_sample_read_length,
            sequence_content_plot,
            plot_per_cycle_N_content,
            plot_adapter_content,
            plot_per_read_gc_hist,
        ]
        for func in plot_functions:
            plot_html, plot_name, anchor, description, helptext, plot_data = func(
                data, group_lookup, project_lookup, self.sample_color
            )
            self.add_section(name=plot_name, plot=plot_html, anchor=anchor, description=description, helptext=helptext)
            self.write_data_file(plot_data, f"base2fastq:{plot_name}")
