# !/usr/bin/env python

""" MultiQC module to parse output from isoseq3 refine """

import json
import pandas as pd
import logging
import re

from collections import OrderedDict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table

# Initialise the logger

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialse the parent object
        super(MultiqcModule, self).__init__(
            name="isoseq3 refine",
            anchor="isoseq3_refine",
            href="https://github.com/PacificBiosciences/IsoSeq",
            info="is for poly(A) tails detection/trimming and concatemer removal.",
            doi="",
        )

        self.isoseq3_refine_data = dict()
        self.isoseq3_refine_data_json = dict()
        self.isoseq3_refine_data_csv = dict()
        self.isoseq3_refine_data_csv_primer = dict()
        self.isoseq3_refine_data_csv_strand = dict()
        self.parse_log()

        self.isoseq3_refine_data_json = self.ignore_samples(self.isoseq3_refine_data_json)
        self.isoseq3_refine_data_csv = self.ignore_samples(self.isoseq3_refine_data_csv)

        # If we found no data
        if not self.isoseq3_refine_data:
            raise UserWarning

        log.info("Found {} reports".format(len(self.isoseq3_refine_data)))

        self.write_data_file(self.isoseq3_refine_data, "multiqc_isoseq3_refine_report")
        self.add_general_stats()
        self.add_sections()

    def parse_log(self):
        for f in self.find_log_files("isoseq3_refine/json", filehandles=True):
            data = json.load(f["f"])
            filename = f["s_name"]
            self.isoseq3_refine_data[filename] = data
            self.isoseq3_refine_data_json[filename] = data
            self.add_data_source(f)

        for f in self.find_log_files("isoseq3_refine/csv", filehandles=True):
            data = pd.read_csv(f["f"])
            filename = f["s_name"]
            data_csv = {}

            data_min = data[["fivelen", "threelen", "polyAlen", "insertlen"]].min().to_frame().transpose()
            data_means = data[["fivelen", "threelen", "polyAlen", "insertlen"]].mean().to_frame().transpose()
            data_std = data[["fivelen", "threelen", "polyAlen", "insertlen"]].std().to_frame().transpose()
            data_max = data[["fivelen", "threelen", "polyAlen", "insertlen"]].max().to_frame().transpose()
            data_csv_strand = data["primer"].value_counts().to_dict()
            data_csv_primer = data["strand"].value_counts().to_dict()

            data_csv["min_fivelen"] = data_min["fivelen"][0]
            data_csv["mean_fivelen"] = data_means["fivelen"][0]
            data_csv["std_fivelen"] = data_std["fivelen"][0]
            data_csv["max_fivelen"] = data_max["fivelen"][0]
            data_csv["min_threelen"] = data_min["threelen"][0]
            data_csv["mean_threelen"] = data_means["threelen"][0]
            data_csv["std_threelen"] = data_std["threelen"][0]
            data_csv["max_threelen"] = data_max["threelen"][0]
            data_csv["min_polyAlen"] = data_min["polyAlen"][0]
            data_csv["mean_polyAlen"] = data_means["polyAlen"][0]
            data_csv["std_polyAlen"] = data_std["polyAlen"][0]
            data_csv["max_polyAlen"] = data_max["polyAlen"][0]
            data_csv["min_insertlen"] = data_min["insertlen"][0]
            data_csv["mean_insertlen"] = data_means["insertlen"][0]
            data_csv["std_insertlen"] = data_std["insertlen"][0]
            data_csv["max_insertlen"] = data_max["insertlen"][0]

            self.isoseq3_refine_data[filename] = data_csv
            self.isoseq3_refine_data_csv[filename] = data_csv
            self.isoseq3_refine_data_csv_strand[filename] = data_csv_strand
            self.isoseq3_refine_data_csv_primer[filename] = data_csv_primer
            self.add_data_source(f)

    def add_general_stats(self):
        gstats_data = {}
        for s_name, attrs in self.isoseq3_refine_data_json.items():
            gstats_data[s_name] = {}
            gstats_data[s_name]["num_reads_fl"] = attrs["num_reads_fl"]
            gstats_data[s_name]["num_reads_flnc"] = attrs["num_reads_flnc"]
            gstats_data[s_name]["num_reads_flnc_polya"] = attrs["num_reads_flnc_polya"]

        headers = OrderedDict()
        headers["num_reads_fl"] = {
            "title": "Number of full length reads",
            "description": "Number of CCS where both primers have been detected",
            "suffix": " fl",
            "format": "{:,}",
            "scale": "spectral",
        }
        headers["num_reads_flnc"] = {
            "title": "Number of non chimeric full length reads",
            "description": "Number of non chimeric CCS where both primers have been detected",
            "suffix": " flnc",
            "format": "{:,}",
            "scale": "RdYlGn",
        }
        headers["num_reads_flnc_polya"] = {
            "title": "Number of poly(A) free non chimeric full length reads",
            "description": "Number of non chimeric CCS where both primers have been detected and the poly(A) tail has been removed",
            "suffix": " flnc_polya",
            "format": "{:,}",
            "scale": "spectral",
        }
        self.general_stats_addcols(gstats_data, headers)

    def add_sections(self):
        self.add_section_fivelen()
        self.add_section_insertlen()
        self.add_section_polyalen()
        self.add_section_threelen()

    def add_section_fivelen(self):
        plot_data = dict()
        for s_name, data in self.isoseq3_refine_data_csv.items():
            plot_data[s_name] = dict()
            plot_data[s_name]["min_fivelen"] = data["min_fivelen"]
            plot_data[s_name]["mean_fivelen"] = data["mean_fivelen"]
            plot_data[s_name]["std_fivelen"] = data["std_fivelen"]
            plot_data[s_name]["max_fivelen"] = data["max_fivelen"]

        headers = OrderedDict()
        headers["min_fivelen"] = {
            "title": "Min 5' primer length",
            "description": "The minimum 5' primer length in base pair",
            "scale": "spectral",
        }
        headers["mean_fivelen"] = {
            "title": "Mean 5' primer length",
            "description": "The mean 5' primer length in base pair",
            "scale": "RdYlGn",
        }
        headers["std_fivelen"] = {
            "title": "Std of 5' primer length",
            "description": "The standard deviation of 5' primer length in base pair",
            "scale": "spectral",
        }
        headers["max_fivelen"] = {
            "title": "Max 5' primer length",
            "description": "The maximum 5' primer length in base pair",
            "scale": "RdYlGn",
        }

        config = {
            "id": "isoseq3_refine_5p_tab",
            "title": "isoseq3 refine: 5 prime primer",
        }

        self.add_section(
            name="5' primer",
            anchor="fivep-stats",
            description="Descriptive stats about 5' primer length.",
            helptext="""
                    The .report.csv file contains information about 5' prime and 3' primers length, insert length,\\ 
                    poly(A) length, and couple of primers detected for each CCS.
                    The table presents min, max, mean, standard deviation for each parameters.
                    """,
            plot=table.plot(plot_data, headers, config),
        )

    def add_section_threelen(self):
        plot_data = dict()
        for s_name, data in self.isoseq3_refine_data_csv.items():
            plot_data[s_name] = dict()
            plot_data[s_name]["min_threelen"] = data["min_threelen"]
            plot_data[s_name]["mean_threelen"] = data["mean_threelen"]
            plot_data[s_name]["std_threelen"] = data["std_threelen"]
            plot_data[s_name]["max_threelen"] = data["max_threelen"]

        headers = OrderedDict()
        headers["min_threelen"] = {
            "title": "Min 3' primer length",
            "description": "The minimum 3' primer length in base pair",
            "scale": "spectral",
        }
        headers["mean_threelen"] = {
            "title": "Mean 3' primer length",
            "description": "The mean 3' primer length in base pair",
            "scale": "RdYlGn",
        }
        headers["std_threelen"] = {
            "title": "Std of 3' primer length",
            "description": "The standard deviation of 3' primer length in base pair",
            "scale": "spectral",
        }
        headers["max_threelen"] = {
            "title": "Max 3' primer length",
            "description": "The maximum 3' primer length in base pair",
            "scale": "RdYlGn",
        }

        config = {
            "id": "isoseq3_refine_3p_tab",
            "title": "isoseq3 refine: 3 prime primer",
        }

        self.add_section(
            name="3' primer",
            anchor="threep-stats",
            description="Descriptive stats about 3' primer length.",
            helptext="""
                    The .report.csv file contains information about 5' prime and 3' primers length, insert length,\\ 
                    poly(A) length, and couple of primers detected for each CCS.
                    The table presents min, max, mean, standard deviation for each parameters.
                    """,
            plot=table.plot(plot_data, headers, config),
        )

    def add_section_polyalen(self):
        plot_data = dict()
        for s_name, data in self.isoseq3_refine_data_csv.items():
            plot_data[s_name] = dict()
            plot_data[s_name]["min_polyAlen"] = data["min_polyAlen"]
            plot_data[s_name]["mean_polyAlen"] = data["mean_polyAlen"]
            plot_data[s_name]["std_polyAlen"] = data["std_polyAlen"]
            plot_data[s_name]["max_polyAlen"] = data["max_polyAlen"]

            headers = OrderedDict()
            headers["min_polyAlen"] = {
                "title": "Min polyA tail length",
                "description": "The minimum polyA tail length in base pair",
                "scale": "spectral",
            }
            headers["mean_polyAlen"] = {
                "title": "Mean polyA tail length",
                "description": "The mean polyA tail length in base pair",
                "scale": "RdYlGn",
            }
            headers["std_polyAlen"] = {
                "title": "Std of polyA tail length",
                "description": "The standard deviation of polyA tail length in base pair",
                "scale": "spectral",
            }
            headers["max_polyAlen"] = {
                "title": "Max polyA tail length",
                "description": "The maximum polyA tail length in base pair",
                "scale": "RdYlGn",
            }

            config = {
                "id": "isoseq3_refine_polya_tab",
                "title": "isoseq3 refine: polyA tail",
            }

        self.add_section(
            name="polyA tails",
            anchor="polya-stats",
            description="Descriptive stats about polyA tail length.",
            helptext="""
                    The .report.csv file contains information about 5' prime and 3' primers length, insert length,\\ 
                    poly(A) length, and couple of primers detected for each CCS.
                    The table presents min, max, mean, standard deviation for each parameters.
                    """,
            plot=table.plot(plot_data, headers, config),
        )

    def add_section_insertlen(self):
        plot_data = dict()
        for s_name, data in self.isoseq3_refine_data_csv.items():
            plot_data[s_name] = dict()
            plot_data[s_name]["min_insertlen"] = data["min_insertlen"]
            plot_data[s_name]["mean_insertlen"] = data["mean_insertlen"]
            plot_data[s_name]["std_insertlen"] = data["std_insertlen"]
            plot_data[s_name]["max_insertlen"] = data["max_insertlen"]

            headers = OrderedDict()
            headers["min_insertlen"] = {
                "title": "Min insert length",
                "description": "The minimum insert length in base pair",
                "scale": "spectral",
            }
            headers["mean_insertlen"] = {
                "title": "Mean insert length",
                "description": "The mean insert length in base pair",
                "scale": "RdYlGn",
            }
            headers["std_insertlen"] = {
                "title": "Std of insert length",
                "description": "The standard deviation of insert length in base pair",
                "scale": "spectral",
            }
            headers["max_insertlen"] = {
                "title": "Max insert length",
                "description": "The maximum insert length in base pair",
                "scale": "RdYlGn",
            }

            config = {
                "id": "isoseq3_refine_insert_tab",
                "title": "isoseq3 refine: insert",
            }

        self.add_section(
            name="Insert length",
            anchor="insert-stats",
            description="Descriptives stats about insert length.",
            helptext="""
                    The .report.csv file contains information about 5' prime and 3' primers length, insert length,\\ 
                    poly(A) length, and couple of primers detected for each CCS.
                    The table presents min, max, mean, standard deviation for each parameters.
                    """,
            plot=table.plot(plot_data, headers, config),
        )
