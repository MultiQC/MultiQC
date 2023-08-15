from multiqc.modules.base_module import BaseMultiqcModule
import os
from io import StringIO
import re
import json
import numpy as np
import pandas as pd
from .plot_runs import *
from .plot_samples import *
import warnings
import seaborn as sns
import copy
import uuid

### Change documentations after development
class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="base2fastq",
            anchor="base2fastq",
            href="https://www.elementbiosciences.com/resources",
            info="is used to call sequences from element AVITI sequencing images",
            doi="10.1038/s41587-023-01750-7",
        )
        self.b2f_data = dict()
        self.b2f_run_data = dict()

        # Read overall stats json as dictionaries
        for f in self.find_log_files("base2fastq/run"):
            data_dict = json.loads(f["f"])
            # sample stats not needed at run level - save on memory
            del data_dict["SampleStats"]
            dir_name = os.path.basename(f["root"])
            run_name = data_dict.get("RunName","UNKNOWN")
            analysis_id = data_dict.get("AnalysisID",self.get_uuid())[0:8]

            sample_name = "__".join([run_name,analysis_id])
            s_name = self.clean_s_name(sample_name, f)
            self.b2f_run_data[s_name] = data_dict

        # Read project info and make it into a lookup dictionary of {sample:project}:
        projectLookupDict = {}
        for f in self.find_log_files("base2fastq/project"):
            data_dict = json.loads(f["f"])
            samples = data_dict["Samples"]
            runName = data_dict["RunName"]
            # run stats no longer needed - save on memory
            del data_dict
            
            project = f["s_name"].replace("_RunStats", "")
            for s in samples:
                projectLookupDict[runName + "_" + s] = project

        # Read per sample stats json as dictionaries
        s_names = []
        for f in self.find_log_files("base2fastq/persample"):
            data_dict = json.loads(f["f"])
            s_name = data_dict["RunName"] + "_" + data_dict["SampleName"]
            if len(data_dict["Reads"]) == 0:
                warnings.warn("Skipping {s} because it does not have any assigned reads.".format(s=s_name))
                continue
            self.b2f_data[s_name] = data_dict
            s_names.append(s_name)

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
        for f in self.find_log_files("base2fastq/group"):
            if self.group_info_exist:
                warnings.warn(
                    "More than one group assignment files are found. Please only keep one assignment file in the analysis folder. Base2fastq stats will not be plotted"
                )
            groupInfo = pd.read_csv(StringIO(f["f"]))
            for nn in groupInfo.index:
                s_group = groupInfo.loc[nn, "Group"]
                if groupInfo.loc[nn, "Run Name"] in groupInfo.loc[nn, "Sample Name"]:
                    s_name = groupInfo.loc[nn, "Sample Name"]
                else:
                    s_name = groupInfo.loc[nn, "Run Name"] + "_" + groupInfo.loc[nn, "Sample Name"]
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
            warnings.warn("No run stats file found!")
        if len(self.b2f_data) == 0:
            warnings.warn("No sample stats file found!")

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
                os.path.dirname(__file__), "assets", "js", "multiqc_base2fastq.js"
            )
        }
        self.intro += '<script type="application/json" class="fastqc_passfails">["fastqc", {"per_base_sequence_content": {"TEST": "pass"}}]</script>'

    def get_uuid(self):
        return str(uuid.uuid4()).replace("-", "").lower()

    def add_run_plots(self):
        plot_functions = [tabulate_run_stats, plot_run_stats, plot_base_quality_hist, plot_base_quality_by_cycle]
        for func in plot_functions:
            plotHtml, plot_name, anchor, description, helptext = func(self.b2f_run_data, self.runColor)
            self.add_section(name=plot_name, plot=plotHtml, anchor=anchor, description=description, helptext=helptext)

    def add_sample_plots(self):
        plot_functions = [
            tabulate_sample_stats,
            plot_adapter_content,
            sequence_content_plot,
            plot_per_cycle_N_content,
            plot_per_read_gc_hist,
        ]
        for func in plot_functions:
            plotHtml, plot_name, anchor, description, helptext = func(
                self.b2f_data, self.groupLookupDict, self.sampleColor
            )
            self.add_section(name=plot_name, plot=plotHtml, anchor=anchor, description=description, helptext=helptext)
