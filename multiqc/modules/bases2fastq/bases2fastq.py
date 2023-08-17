import copy
import json
import logging
import os
import uuid
from io import StringIO

import numpy as np
import pandas as pd
import seaborn as sns

from multiqc.modules.base_module import BaseMultiqcModule

from .plot_runs import *
from .plot_samples import *

log = logging.getLogger(__name__)


### Change documentations after development
class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="bases2fastq",
            anchor="bases2fastq",
            href="https://www.elementbiosciences.com/resources",
            info="is used to call sequences from element AVITI sequencing images",
            doi="10.1038/s41587-023-01750-7",
        )

        self.minimum_polonies = 10000

        self.b2f_data = dict()
        self.b2f_run_data = dict()
        self.missing_runs = set()

        self.sample_id_to_run = dict()

        # Read overall stats json as dictionaries
        run_prefix_len = None
        for f in self.find_log_files("bases2fastq/run"):
            data_dict = json.loads(f["f"])

            # get run + analysis
            run_name = data_dict.get("RunName", None)
            analysis_id = data_dict.get("AnalysisID", None)[0:4]

            if not run_name or not analysis_id:
                log.error("Error with RunStats.json.  Either RunName or AnalysisID is absent.")
                log.error("Please visit Elembio docs for more information - https://docs.elembio.io/docs/bases2fastq/")
                raise UserWarning

            run_analysis_name = "__".join([run_name, analysis_id])

            # map sample UUIDs to run_analysis_name
            for sample in data_dict["SampleStats"]:
                self.sample_id_to_run[sample["SampleID"]] = run_analysis_name

            # sample stats not needed at run level - save on memory
            del data_dict["SampleStats"]

            self.b2f_run_data[run_analysis_name] = data_dict
            self.add_data_source(f=f, s_name=run_analysis_name, module="bases2fastq")
        # if all RunStats.json too large, none will be found.  Guide customer and Exit at this point.
        if len(self.sample_id_to_run) == 0:
            log.error("No run-stats were found.  Either file-size above limit or RunStats.json does not exist.")
            log.error("Please visit Elembio docs for more information - https://docs.elembio.io/docs/bases2fastq/")
            raise UserWarning

        log.info(f"Found {len(self.b2f_run_data)} total RunStats.json")
        run_r1r2_lens = [
            str(len(self.b2f_run_data[s]["Reads"][0]["Cycles"]))
            + "+"
            + str(len(self.b2f_run_data[s]["Reads"][1]["Cycles"]))
            for s in self.b2f_run_data.keys()
        ]
        run_r1r2_lens_set = set(run_r1r2_lens)
        run_r1r2_lens_dict = {}
        for nn, rl in enumerate(run_r1r2_lens):
            if not run_r1r2_lens_dict.get(rl):
                run_r1r2_lens_dict[rl] = []
            run_r1r2_lens_dict[rl].append(list(self.b2f_run_data.keys())[nn])

        if len(run_r1r2_lens_set) > 1:
            log.warning(
                f"More than one read length configurations are found in the dataset:{','.join(run_r1r2_lens_set)}"
            )
            log.warning(
                f"Runnning MultiQC with different read length configurations may cause unusual plotting behavior. If possible, please split runs with different read length configurations into different folder and run MultiQC on each"
            )
            for rl in run_r1r2_lens_dict.keys():
                log.warning(f"These runs have {rl} read length configuration:{','.join(run_r1r2_lens_dict[rl])}")

        # Read project info and make it into a lookup dictionary of {sample:project}:
        project_num = 0
        projectLookupDict = {}
        for f in self.find_log_files("bases2fastq/project"):
            data_dict = json.loads(f["f"])
            samples = data_dict["Samples"]

            # get run + analysis
            run_name = data_dict.get("RunName", None)
            analysis_id = data_dict.get("AnalysisID", None)[0:4]

            if not run_name or not analysis_id:
                log.error(f"Error with {f['root']}.  Either RunName or AnalysisID is absent.")
                log.error("Please visit Elembio docs for more information - https://docs.elembio.io/docs/bases2fastq/")
                raise UserWarning

            run_analysis_name = "__".join([run_name, analysis_id])

            project = data_dict.get("Project", "DefaultProject")

            # run stats no longer needed - save on memory
            del data_dict

            for sample_name in samples:
                run_analysis_sample_name = "__".join([run_analysis_name, sample_name])
                projectLookupDict[run_analysis_sample_name] = project
            project_num += 1
            self.add_data_source(f=f, s_name=project, module="bases2fastq")

        log.info(f"Found {project_num} Projects within bases2fastq results")

        total_sample = 0
        # Read per sample stats json as dictionaries
        for f in self.find_log_files("bases2fastq/persample"):
            total_sample += 1
            data_dict = json.loads(f["f"])

            # ensure sample UUID is known to list of runs
            sample_id = data_dict["SampleID"]
            if sample_id not in self.sample_id_to_run:
                log.warning(f"{data_dict['RunName']} RunStats.json is missing for sample, {f['root']}")
                continue

            run_analysis_name = self.sample_id_to_run[sample_id]

            sample_name = data_dict["SampleName"]
            run_analysis_sample_name = "__".join([run_analysis_name, sample_name])

            num_polonies = data_dict["NumPolonies"]
            if num_polonies < self.minimum_polonies:
                log.warning(
                    f"Skipping {run_analysis_sample_name} because it has <{self.minimum_polonies} assigned reads [n={num_polonies}]."
                )
                continue

            self.b2f_data[run_analysis_sample_name] = data_dict
            self.b2f_data[run_analysis_sample_name]["RunName"] = run_analysis_name

            self.add_data_source(f=f, s_name=run_analysis_sample_name, module="bases2fastq")
        log.info(
            f"Found {total_sample} samples within bases2fastq results, and {len(self.b2f_data)} samples have run information and enough polonies"
        )

        # Group by run name
        self.groupDict = dict()
        self.groupLookupDict = dict()
        for s_name in self.b2f_data.keys():
            s_group = self.b2f_data[s_name]["RunName"]

            if not self.groupDict.get(s_group):
                self.groupDict.update({s_group: []})

            self.groupDict[s_group].append(s_name)
            self.groupLookupDict.update({s_name: s_group})

        # Assign color for each run
        n_colors = len(self.groupDict.keys())
        palette = [
            "rgba({r},{g},{b},0.5)".format(r=rgb[0] * 255, g=rgb[1] * 255, b=rgb[2] * 255)
            for rgb in sns.color_palette("bright", n_colors)
        ]
        groupColor = {g: c for g, c in zip(self.groupDict.keys(), palette[: len(self.groupDict.keys())])}
        self.runColor = copy.deepcopy(groupColor)  #
        self.sampleColor = dict()
        for s_name in self.b2f_data.keys():
            self.sampleColor.update({s_name: groupColor[self.groupLookupDict[s_name]]})

        for s_name in self.b2f_data.keys():
            self.sampleColor.update({s_name: groupColor[self.groupLookupDict[s_name]]})

        # Assign color for each project
        for s_name in self.b2f_data.keys():
            if projectLookupDict.get(s_name):
                s_group = projectLookupDict[s_name]
                if not self.groupDict.get(s_group):
                    self.groupDict.update({s_group: []})
                self.groupDict[s_group].append(s_name)
                self.groupLookupDict.update({s_name: s_group})

        # Read custom group info
        self.group_info_exist = False
        for f in self.find_log_files("bases2fastq/group"):
            if self.group_info_exist:
                log.warning(
                    "More than one group assignment files are found. Please only keep one assignment file in the analysis folder. Bases2fastq stats will not be plotted"
                )
            groupInfo = pd.read_csv(StringIO(f["f"]))
            for nn in groupInfo.index:
                s_group = groupInfo.loc[nn, "Group"]
                # if groupInfo.loc[nn, "Run Name"] in groupInfo.loc[nn, "Sample Name"]:
                s_name = groupInfo.loc[nn, "Sample Name"]
                # else:
                # s_name = groupInfo.loc[nn, "Run Name"] + "_" + groupInfo.loc[nn, "Sample Name"]
                if self.groupDict.get(s_group) is None:
                    self.groupDict.update({s_group: []})
                self.groupDict[s_group].append(s_name)
                self.groupLookupDict.update({s_name: s_group})

        n_colors = len(self.groupDict.keys())
        palette = [
            "rgba({r},{g},{b},0.5)".format(r=rgb[0] * 255, g=rgb[1] * 255, b=rgb[2] * 255)
            for rgb in sns.color_palette("Paired", n_colors)
        ]
        groupColor = {g: c for g, c in zip(self.groupDict.keys(), palette[: len(self.groupDict.keys())])}
        self.sampleColor = dict()
        for s_name in self.b2f_data.keys():
            self.sampleColor.update({s_name: groupColor[self.groupLookupDict[s_name]]})

        # Sort samples alphabetically
        data_keys = list(self.b2f_run_data.keys())
        data_keys.sort()
        sorted_data = {s_name: self.b2f_run_data[s_name] for s_name in data_keys}
        self.b2f_run_data = sorted_data
        data_keys = list(self.b2f_data.keys())
        sorted_keys = sorted(data_keys, key=lambda x: (self.groupLookupDict[x], x))
        sorted_data = {s_name: self.b2f_data[s_name] for s_name in sorted_keys}
        self.b2f_data = sorted_data

        if len(self.b2f_run_data) == 0:
            log.warning("No run stats file found!")
        if len(self.b2f_data) == 0:
            log.warning("No sample stats file found!")

        # Add sections
        self.add_run_plots()
        self.add_sample_plots()

        # Add css and js
        self.css = {
            "assets/css/multiqc_fastqc.css": os.path.join(
                os.path.dirname(__file__), "..", "fastqc", "assets", "css", "multiqc_fastqc.css"
            )
        }
        self.js = {
            "assets/js/multiqc_dragen_fastqc.js": os.path.join(
                os.path.dirname(__file__), "assets", "js", "multiqc_bases2fastq.js"
            )
        }
        self.intro += '<script type="application/json" class="fastqc_passfails">["fastqc", {"per_base_sequence_content": {"TEST": "pass"}}]</script>'

    def get_uuid(self):
        return str(uuid.uuid4()).replace("-", "").lower()

    def add_run_plots(self):
        plot_functions = [tabulate_run_stats, plot_run_stats, plot_base_quality_hist, plot_base_quality_by_cycle]
        for func in plot_functions:
            plotHtml, plot_name, anchor, description, helptext, plot_data = func(self.b2f_run_data, self.runColor)
            self.add_section(name=plot_name, plot=plotHtml, anchor=anchor, description=description, helptext=helptext)
            self.write_data_file(plot_data, f"base2fastq:{plot_name}")

    def add_sample_plots(self):
        plot_functions = [
            tabulate_sample_stats,
            plot_adapter_content,
            sequence_content_plot,
            plot_per_cycle_N_content,
            plot_per_read_gc_hist,
        ]
        for func in plot_functions:
            plotHtml, plot_name, anchor, description, helptext, plot_data = func(
                self.b2f_data, self.groupLookupDict, self.sampleColor
            )
            self.add_section(name=plot_name, plot=plotHtml, anchor=anchor, description=description, helptext=helptext)
            self.write_data_file(plot_data, f"base2fastq:{plot_name}")
