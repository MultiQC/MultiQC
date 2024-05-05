#!/usr/bin/env python

"""MultiQC module to parse output from HTStream"""

from __future__ import print_function
from collections import OrderedDict
import logging
import json
import os

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.modules.base_module import ModuleNoSamplesFound

# Import modules
from .apps import __hts_import
from .apps import htstream_utils

#################################################

# Logger Initialization
log = logging.getLogger(__name__)

# Config Initialization
hconfig = {}

# get supported apps
supported_apps = __hts_import.supported_apps


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        self.sample_statistics = {}

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="HTStream",
            anchor="htstream",
            href="https://s4hts.github.io/HTStream/",
            info=" quality control and processing pipeline for High Throughput Sequencing data.",
            doi="",
        )

        # Initialize ordered dictionary (key: samples, values: their respective json files)
        self.htstream_data = OrderedDict()
        self.overview_stats = {}
        self.report_sections = {}
        self.app_order = []
        self.repeated_apps = []
        self.add_software_version(None)

        # Import js and css functions.
        self.js = {"assets/js/htstream.js": os.path.join(os.path.dirname(__file__), "assets", "js", "htstream.js")}
        self.css = {"assets/css/htstream.css": os.path.join(os.path.dirname(__file__), "assets", "css", "htstream.css")}

        # iterates through files found by "find_log_files" (located in base_module.py, re patterns found in search_patterns.yml)
        for file in self.find_log_files("htstream"):
            # clean sample name
            s_name = self.clean_s_name(file["s_name"], file["root"])  # sample name

            # do not parse excluded samples
            if self.is_ignore_sample(s_name):
                continue

            # exclude sample name if necessary (other implementation of this function wasn't working)
            if s_name in self.htstream_data.keys():
                log.debug("Duplicate sample name found! Overwriting: {}".format(file["s_name"]))

            # parse json file
            file_data, repeated_apps = htstream_utils.parse_json(
                file["s_name"], file["f"]
            )  # parse stats file. Should return json directory of apps and their stats

            # Add data to MultiQC and data structures
            self.add_data_source(file)
            self.htstream_data[s_name] = file_data
            self.repeated_apps = repeated_apps

        self.htstream_data = self.ignore_samples(self.htstream_data)

        # make sure samples are being processed
        if len(self.htstream_data) == 0:
            log.warning("Found " + str(len(self.htstream_data)) + " reports.")
            raise ModuleNoSamplesFound

        log.info("Found " + str(len(self.htstream_data)) + " reports.")

        # report number of files found
        log.info("Found " + str(len(self.htstream_data)) + " reports.")

        # parse config file
        hconfig = getattr(config, "htstream_config", {})

        # if 'sample_colors' parameter is specified, apply colors
        if "sample_colors" in hconfig.keys():
            try:
                color_file = hconfig["sample_colors"]
                hconfig["sample_colors"] = {}

                # Read in specified sample colors
                with open(color_file, "r") as f:
                    lines = f.readlines()

                    for line in lines:
                        if line.strip() != "":
                            split_line = line.split("\t")
                            hconfig["sample_colors"][split_line[0]] = split_line[1].strip()

            except:
                log.warning(
                    "Sample Coloring file could not be parsed. Check for proper path, TSV format, and HTML color codes."
                )
                raise

        # number of smaples
        hconfig["htstream_number_of_samples"] = len(self.htstream_data.keys())

        # intro json holding HTStream specific parameters, FASTQC does similar things but with actual data
        self.intro += '<div id="htstream_config" style="display: none;">{}</div>'.format(json.dumps(hconfig))

        # parse json containing stats on each sample and generate reports
        self.process_data(self.htstream_data)
        self.generate_reports()

    #########################
    # Iterate through files and generate report
    def process_data(self, json):
        # checks that order is consistent within stats files
        for key in json.keys():
            temp = list(json[key].keys())

            if self.app_order == []:
                self.app_order = temp

            elif self.app_order == temp:
                continue

            else:
                log.error("Inconsistent order of HTStream applications.")

        # sort list of samples
        sample_keys = list(sorted(json.keys()))

        # initialize some more useful variables
        self.overview_stats = {"Pipeline Input": {}, "details": {"read_reducer": [], "bp_reducer": []}}
        supported_apps_keys = list(supported_apps.keys())

        pipeline_input = True

        # process data
        for i in range(len(self.app_order)):
            # clean up app name
            app = self.app_order[i]
            program = app.split("hts_")[-1].split("_")[0]

            # Check if app is supported
            if program not in supported_apps_keys:
                log.warning(
                    "hts_"
                    + program
                    + " is currently not supported by MultiQC: HTStrean. Apps currently supported: "
                    + htstream_utils.key_print(self.programs)
                )
                continue

            # creat app specific dictionary, each entry will be a sample
            stats_dict = OrderedDict()

            # get input reads
            for key in sample_keys:
                # get app
                stats_dict[key] = json[key][app]

                # populate info for first app in pipeline
                if pipeline_input and len(self.app_order) > 1:
                    self.overview_stats["Pipeline Input"][key] = {
                        "Input_Reads": stats_dict[key]["Fragment"]["in"],
                        "Input_Bps": stats_dict[key]["Fragment"]["basepairs_in"],
                    }

            pipeline_input = False

            # if data exists for app, execute app specific stats processing
            if len(stats_dict.keys()) != 0:
                app_name = app
                index = app_name.split("_")[-1]
                app = supported_apps[program]()

                # add app to list of read or bp reducers
                if app.type == "both":
                    self.overview_stats["details"]["read_reducer"].append(app_name)
                    self.overview_stats["details"]["bp_reducer"].append(app_name)

                else:
                    self.overview_stats["details"][app.type].append(app_name)

                # dictionary of subsections
                description = app.info
                section_dict = app.execute(stats_dict, index)
                self.overview_stats[app_name] = section_dict["Overview"]

                try:
                    notes = stats_dict[list(stats_dict.keys())[1]]["Program_details"]["options"]["notes"]
                except:
                    notes = ""
                    raise

                if program != "Stats":
                    # if dictionary is not empty
                    if len(section_dict.keys()) != 0:
                        figure = section_dict["Figure"]
                        html = section_dict["Content"]

                        # add report section to dictionary, to be added later
                        self.report_sections[app_name] = {
                            "name": app_name,
                            "description": description,
                            "figure": figure,
                            "html": html,
                            "comment": notes,
                        }

                else:
                    self.report_sections[app_name + " (Summary Table)"] = {
                        "name": app_name + " (Table)",
                        "description": "Basic statistics about the reads in a dataset.",
                        "figure": section_dict["Table"],
                        "html": "",
                        "comment": notes,
                    }

                    read_types = []
                    if "PE" in section_dict.keys():
                        read_types.append("PE")

                    if "SE" in section_dict.keys():
                        read_types.append("SE")

                    for rtype in read_types:
                        self.report_sections[app_name + " (" + rtype + " Read Length)"] = {
                            "name": app_name + " (" + rtype + " Read Length)",
                            "description": "Distribution of read lengths for each sample.",
                            "figure": section_dict[rtype]["Read_Length"],
                            "html": "",
                            "comment": "",
                        }

                        self.report_sections[app_name + " (" + rtype + " PE Base by Cycle)"] = {
                            "name": app_name + " (" + rtype + " Base by Cycle)",
                            "description": """Provides a measure of the uniformity of a distribution. 
                                                                                                    The higher the average deviation from 25% is, 
                                                                                                    the more unequal the base pair composition. N's are excluded from this calculation.""",
                            "figure": section_dict[rtype]["Base_by_Cycle"],
                            "html": "",
                            "comment": "",
                        }

                        self.report_sections[app_name + " (" + rtype + " PE Quality by Cycle)"] = {
                            "name": app_name + " (" + rtype + " Quality by Cycle)",
                            "description": """Mean quality score for each position along the read. 
                                                                                                    Sample is colored red if less than 60% of bps have mean score of at least Q30, 
                                                                                                    orange if between 60% and 80%, and green otherwise.""",
                            "figure": section_dict[rtype]["Quality_by_Cycle"],
                            "html": "",
                            "comment": "",
                        }

    ############################
    # add sections
    def generate_reports(self):
        # add pipeline overview section if appropriate
        if self.overview_stats != {} and len(self.app_order) > 1:
            # try adding overview section
            try:
                app = supported_apps["OverviewStats"]()

                description = "Plots reduction of reads and basepairs across the preprocessing pipeline."
                plot = app.execute(self.overview_stats, self.app_order)

                self.add_section(name="Processing Overview", description=description, plot=plot)

            except:
                log.warning("Report Section for Processing Overview Failed.")
                raise

            # try adding cols to general table
            try:
                # collect data for general table in proper format
                last_app = list(self.overview_stats.keys())[-1]
                self.gen_report_stats = {}
                for k, v in self.overview_stats["Pipeline Input"].items():
                    self.gen_report_stats[k] = {
                        "Input Reads": self.overview_stats["Pipeline Input"][k]["Input_Reads"] / 1000000,
                        "Output Reads": self.overview_stats[last_app][k]["Output_Reads"] / 1000000,
                    }

                # create header dict for table
                headers = OrderedDict()
                headers["Input Reads"] = {
                    "title": "M Input",
                    "description": "Reads input to HTStream (millions).",
                    "format": "{:,.2f}",
                    "scale": "Blues",
                }
                headers["Output Reads"] = {
                    "title": "M Output",
                    "description": "Reads remaining after HTStream preprocessing (millions).",
                    "format": "{:,.2f}",
                    "scale": "Greens",
                }

                # add data to general table
                self.general_stats_addcols(self.gen_report_stats, headers)

            except:
                log.warning("Adding Columns to General Table Failed.")
                raise

        # add app sections
        for section, content in self.report_sections.items():
            if "hts_" + section.split("_")[1] not in self.repeated_apps:
                content["name"] = section[:-2]

            try:
                self.add_section(
                    name=content["name"],
                    description=content["description"],
                    plot=content["figure"],
                    comment=content["comment"],
                    content=content["html"],
                )

            except:
                msg = "Report Section for " + section + " Failed."
                log.warning(msg)
                raise

        self.write_data_file(self.htstream_data, "multiqc_htstream")
