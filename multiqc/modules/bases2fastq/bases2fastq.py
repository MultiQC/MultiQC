import copy
import csv
import json
import logging
import random
import uuid

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.utils import mqc_colour

from multiqc.modules.bases2fastq.plot_runs import (
    plot_run_stats,
    tabulate_run_stats,
    plot_base_quality_hist,
    plot_base_quality_by_cycle,
)
from multiqc.modules.bases2fastq.plot_project_runs import tabulate_project_run_stats
from multiqc.modules.bases2fastq.plot_samples import (
    tabulate_sample_stats,
    sequence_content_plot,
    plot_per_cycle_N_content,
    plot_adapter_content,
    plot_per_read_gc_hist,
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

        self.b2f_sample_data = dict()
        self.b2f_run_data = dict()
        self.b2f_run_project_data = dict()
        self.missing_runs = set()
        self.sample_id_to_run = dict()

        # Group by run name
        self.group_dict = dict()
        self.group_lookup_dict = dict()
        self.project_lookup_dict = dict()

        # bases2fastq/run
        num_runs = 0
        num_samples = 0
        for f in self.find_log_files("bases2fastq/run"):
            data = json.loads(f["f"])

            # get run + analysis
            run_name = data.get("RunName", None)
            analysis_id = data.get("AnalysisID", None)[0:4]

            if not run_name or not analysis_id:
                log.error("Error with RunStats.json. Either RunName or AnalysisID is absent.")
                log.error(
                    "Please visit Elembio online documentation for more information - https://docs.elembio.io/docs/bases2fastq/introduction/"
                )
                continue

            run_analysis_name = "-".join([run_name, analysis_id])
            run_analysis_name = self.clean_s_name(run_analysis_name, f)

            # map sample UUIDs to run_analysis_name
            for sample_data in data["SampleStats"]:
                sample_id = sample_data["SampleID"]
                sample_name = sample_data["SampleName"]
                sample_data["RunName"] = run_name

                run_analysis_sample_name = "__".join([run_analysis_name, sample_name])

                num_polonies = sample_data["NumPolonies"]
                if num_polonies < MIN_POLONIES:
                    log.warning(
                        f"Skipping {run_analysis_sample_name} because it has <{MIN_POLONIES} assigned reads [n={num_polonies}]."
                    )
                    continue

                # skip run if in user provider ignore list
                if self.is_ignore_sample(sample_id):
                    continue
                if self.is_ignore_sample(run_analysis_sample_name):
                    continue

                self.sample_id_to_run[sample_id] = run_analysis_name
                self.b2f_sample_data[run_analysis_sample_name] = sample_data
                num_samples += 1

            # skip run if in user provider ignore list
            if self.is_ignore_sample(run_analysis_name):
                continue

            num_runs += 1
            self.b2f_run_data[run_analysis_name] = data
            self.add_data_source(f=f, s_name=run_analysis_name, module="bases2fastq")

        # Checking if run lengths configurations are the same for all samples.
        self.run_r1r2_lens = []
        for s in self.b2f_run_data.keys():
            read_lens = str(len(self.b2f_run_data[s]["Reads"][0]["Cycles"]))
            if len(self.b2f_run_data[s]["Reads"]) > 1:
                read_lens += "+" + str(len(self.b2f_run_data[s]["Reads"][1]["Cycles"]))
            self.run_r1r2_lens.append(read_lens)

        run_r1r2_lens_dict = {}
        for nn, rl in enumerate(self.run_r1r2_lens):
            if not run_r1r2_lens_dict.get(rl):
                run_r1r2_lens_dict[rl] = []
            run_r1r2_lens_dict[rl].append(list(self.b2f_run_data.keys())[nn])

        #
        # bases2fastq/project
        #
        num_projects = 0
        for f in self.find_log_files("bases2fastq/project"):
            data = json.loads(f["f"])
            samples = data["Samples"]

            # get run + analysis
            run_name = data.get("RunName", None)
            analysis_id = data.get("AnalysisID", None)[0:4]

            run_analysis_name = "-".join([run_name, analysis_id])
            run_analysis_name = self.clean_s_name(run_analysis_name, f)

            if not run_name or not analysis_id:
                log.error(f"Error with {f['root']}.  Either RunName or AnalysisID is absent.")
                log.error("Please visit Elembio online documentation for more information -")
                continue

            project = self.clean_s_name(data.get("Project", "DefaultProject"), f)

            run_analysis_project_name = "__".join([run_name, project, analysis_id])
            run_analysis_project_name = self.clean_s_name(run_analysis_project_name, f)

            # skip project if in user provider ignore list
            if self.is_ignore_sample(run_analysis_project_name):
                continue

            for sample_name in samples:
                run_analysis_sample_name = self.clean_s_name("__".join([run_analysis_name, sample_name]), f)
                self.project_lookup_dict[run_analysis_sample_name] = project
            num_projects += 1

            # remove samples
            del data["Samples"]

            self.b2f_run_project_data[run_analysis_project_name] = data
            self.add_data_source(f=f, s_name=project, module="bases2fastq")

        # if all RunStats.json too large, none will be found.  Guide customer and Exit at this point.
        if len(self.sample_id_to_run) != 0:
            log.info(f"Found {num_runs} total RunStats.json")

        # ensure run/sample data found
        if num_projects == 0 and num_samples == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {num_samples} samples and {num_projects} projects within the bases2fastq results")

        # Superfluous function call to confirm that it is used in this module
        self.add_software_version(None)

        # process groups / projects
        for s_name in self.b2f_sample_data.keys():
            s_group = self.b2f_sample_data[s_name]["RunName"]

            if not self.group_dict.get(s_group):
                self.group_dict.update({s_group: []})

            self.group_dict[s_group].append(s_name)
            self.group_lookup_dict.update({s_name: s_group})

        # Assign project
        for s_name in self.b2f_sample_data.keys():
            if self.project_lookup_dict.get(s_name):
                s_group = self.project_lookup_dict[s_name]
                if not self.group_dict.get(s_group):
                    self.group_dict.update({s_group: []})
                self.group_dict[s_group].append(s_name)
                self.group_lookup_dict.update({s_name: s_group})

        # Assign color for each group
        self.color_getter = mqc_colour.mqc_colour_scale()
        self.palette = sum(
            [
                self.color_getter.get_colours(hue)
                for hue in ["Set2", "Pastel1", "Accent", "Set1", "Set3", "Dark2", "Paired", "Pastel2"]
            ],
            [],
        )
        if len(self.group_dict) > len(self.palette):
            hex_range = 2**24
            extra_colors = [hex(random.randrange(0, hex_range)) for _ in range(len(self.group_dict), len(self.palette))]
            self.palette = self.palette + extra_colors
        self.group_color = {g: c for g, c in zip(self.group_dict.keys(), self.palette[: len(self.group_dict)])}
        self.sample_color = dict()
        for s_name in self.b2f_sample_data.keys():
            self.sample_color.update({s_name: self.group_color[self.group_lookup_dict[s_name]]})
        self.run_color = copy.deepcopy(self.group_color)  # Make sure that run colors and group colors match
        self.palette = self.palette[len(self.group_dict) :]

        # Read custom group info
        self.group_info_exist = False
        for f in self.find_log_files("bases2fastq/group"):
            if self.group_info_exist:
                log.warning(
                    "More than one group assignment files are found. Please only keep "
                    "one assignment file in the analysis folder. Bases2Fastq stats will "
                    "not be plotted"
                )
            for row in csv.DictReader(f["f"]):
                s_group = row["Group"]
                s_name = row["Sample Name"]
                if self.group_dict.get(s_group) is None:
                    self.group_dict[s_group] = []
                self.group_dict[s_group].append(s_name)
                self.group_lookup_dict[s_name] = s_group
        for group in self.group_dict.keys():
            if group not in self.run_color:
                if len(self.palette) > 0:
                    self.group_color[group] = self.palette.pop(0)
                else:
                    hex_range = 2**24
                    extra_color = hex(random.randrange(0, hex_range))
                    self.group_color[group] = extra_color
        self.sample_color = dict()
        for s_name in self.b2f_sample_data.keys():
            self.sample_color.update({s_name: self.group_color[self.group_lookup_dict[s_name]]})

        # sort run
        data_keys = list(self.b2f_run_data.keys())
        data_keys.sort()
        sorted_data = {s_name: self.b2f_run_data[s_name] for s_name in data_keys}
        self.b2f_run_data = sorted_data
        # sort projects
        data_keys = list(self.b2f_run_project_data.keys())
        data_keys.sort()
        sorted_data = {s_name: self.b2f_run_project_data[s_name] for s_name in data_keys}
        self.b2f_run_project_data = sorted_data
        # sort samples
        data_keys = list(self.b2f_sample_data.keys())
        sorted_keys = sorted(data_keys, key=lambda x: (self.group_lookup_dict[x], x))
        sorted_data = {s_name: self.b2f_sample_data[s_name] for s_name in sorted_keys}
        self.b2f_sample_data = sorted_data

        if len(self.b2f_run_data) == 0:
            log.warning("No run stats file found!")
        if len(self.b2f_sample_data) == 0:
            log.warning("No sample stats file found!")

        # Add sections
        self.add_run_plots()
        if num_projects > 0:
            self.add_project_run_plots()
        self.add_sample_plots()

    def get_uuid(self):
        return str(uuid.uuid4()).replace("-", "").lower()

    def add_run_plots(self):
        plot_functions = [tabulate_run_stats, plot_run_stats, plot_base_quality_hist, plot_base_quality_by_cycle]
        for func in plot_functions:
            plot_html, plot_name, anchor, description, helptext, plot_data = func(self.b2f_run_data, self.run_color)
            self.add_section(name=plot_name, plot=plot_html, anchor=anchor, description=description, helptext=helptext)
            self.write_data_file(plot_data, f"base2fastq:{plot_name}")

    def add_project_run_plots(self):
        plot_functions = [tabulate_project_run_stats]
        for func in plot_functions:
            plot_html, plot_name, anchor, description, helptext, plot_data = func(
                self.b2f_run_project_data, self.run_color
            )
            self.add_section(name=plot_name, plot=plot_html, anchor=anchor, description=description, helptext=helptext)
            self.write_data_file(plot_data, f"base2fastq_projects:{plot_name}")

    def add_sample_plots(self):
        plot_functions = [
            tabulate_sample_stats,
            sequence_content_plot,
            plot_per_cycle_N_content,
            plot_adapter_content,
            plot_per_read_gc_hist,
        ]
        for func in plot_functions:
            plot_html, plot_name, anchor, description, helptext, plot_data = func(
                self.b2f_sample_data, self.group_lookup_dict, self.project_lookup_dict, self.sample_color
            )
            self.add_section(name=plot_name, plot=plot_html, anchor=anchor, description=description, helptext=helptext)
            self.write_data_file(plot_data, f"base2fastq:{plot_name}")
