from multiqc.modules.base_module import BaseMultiqcModule
import os
from io import StringIO
import re
import json
import numpy as np
import pandas as pd
from .plot_runs import *
from .plot_samples import *
import seaborn as sns
import copy
import uuid
import logging

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

        root_to_analysis_id = dict()
        
        # Read overall stats json as dictionaries
        run_index = None
        for f in self.find_log_files("bases2fastq/run"):
            data_dict = json.loads(f["f"])
            
            # sample stats not needed at run level - save on memory
            del data_dict["SampleStats"]

            run_name = data_dict.get("RunName","UNKNOWN")
            analysis_id = data_dict.get("AnalysisID","123")[0:3]
            run_analysis_name = "__".join([run_name,analysis_id])

            run_root_arr = f['root'].rstrip('/').split("/")
            run_index = len(run_root_arr)-1
            run_root = run_root_arr[run_index]

            root_to_analysis_id[run_root] = analysis_id
        
            self.b2f_run_data[run_analysis_name] = data_dict

        # if all RunStats.json too large, none will be found.  Guide customer and Exit at this point.
        if not run_index:
            log.error("No run-stats were found.  Either file-size above limit or RunStats.json does not exist.")
            log.error("Please visit Elembio docs for more information - https://docs.elembio.io/docs/bases2fastq/")
            raise UserWarning

        # Read project info and make it into a lookup dictionary of {sample:project}:
        projectLookupDict = {}
        for f in self.find_log_files("bases2fastq/project"):
            data_dict = json.loads(f["f"])
            samples = data_dict["Samples"]

            run_root_arr = f['root'].rstrip('/').split("/")
            run_root = run_root_arr[run_index]

            run_name = data_dict.get("RunName","UNKNOWN")
            analysis_id = root_to_analysis_id[run_root]
            run_analysis_name = "__".join([run_name,analysis_id])
            
            if run_name not in self.b2f_run_data and run_name not in self.missing_runs:
                log.warning(f"{run_name} is missing from run_data - file size is too large, modify config or remove run")
                self.missing_runs.add(run_name)
                continue

            project = data_dict.get("Project","DefaultProject")
            
            # run stats no longer needed - save on memory
            del data_dict
            
            for sample_name in samples:
                run_analysis_sample_name = "__".join([run_analysis_name,sample_name])
                projectLookupDict[run_analysis_sample_name] = project

        # Read per sample stats json as dictionaries
        for f in self.find_log_files("bases2fastq/persample"):
            data_dict = json.loads(f["f"])

            run_root_arr = f['root'].rstrip('/').split("/")
            run_root = run_root_arr[run_index]

            run_name = data_dict.get("RunName","UNKNOWN")
            analysis_id = root_to_analysis_id[run_root]
            run_analysis_name = "__".join([run_name,analysis_id])

            if run_name not in self.b2f_run_data and run_name not in self.missing_runs:
                log.warning(f"{run_name} is missing from run_data - file size is too large, modify config or remove run")
                self.missing_runs.add(run_name)
                continue

            sample_name = data_dict["SampleName"]
            run_analysis_sample_name = "__".join([run_analysis_name,sample_name])

            # todo - check for read length, ensure all samples are the same otherwise exit
            
            num_polonies = data_dict["NumPolonies"]
            if num_polonies < self.minimum_polonies:
                log.warning(f"Skipping {run_analysis_sample_name} because it has {num_polonies} < {self.minimum_polonies} assigned reads.")
                continue

            self.b2f_data[run_analysis_sample_name] = data_dict
            self.b2f_data[run_analysis_sample_name]["RunName"] = run_analysis_name 

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
                #if groupInfo.loc[nn, "Run Name"] in groupInfo.loc[nn, "Sample Name"]:
                s_name = groupInfo.loc[nn, "Sample Name"]
                #else:
                    #s_name = groupInfo.loc[nn, "Run Name"] + "_" + groupInfo.loc[nn, "Sample Name"]
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
