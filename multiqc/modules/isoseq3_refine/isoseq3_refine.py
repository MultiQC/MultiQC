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

        # self.add_section(
        #     name="Full length CCS report",
        #     anchor="fl-report",
        #     description="Detail the CCS selected and trimmed by isoseq3 refine",
        #     helptext="""
        #              _fl / Full Length_
        #              * Both primers have been detected and trimmed from CCS
        #              _flnc / Full Length Non Chimeric_
        #              * Both primers have been detected and trimmed from CCS
        #              * No chimeric structure has been detected (ex: primer inside the CCS)
        #              _flnc_polya / Full Length Non Chimeric no poly(A) tail_
        #              * Both primers have been detected and trimmed from CCS
        #              * No chimeric structure has been detected (ex: primer inside the CCS)
        #              * The poly(A) tails has been detected and trimmed
        #             """,
        #     plot=table.plot(self.isoseq3_refine_data_json),
        # )

        self.add_section(
            name="Filtered and trimmed CCS descriptive stats",
            anchor="ccs-stats",
            description="Descriptives stats about 5'/3' primerlen, insert and polyA length.",
            helptext="""
            The .report.csv file contains information about 5' prime and 3' primers length, insert length, poly(A) length, and couple of primers detected for each CCS.  
            The table presents min, max, mean, standard deviation for each parameters.
            """,
            plot=table.plot(self.isoseq3_refine_data_csv),
        )

        self.add_section(
            name="Strand counting",
            anchor="strand-counting",
            description="Count of number of CCS on forward strand and reverse strand",
            helptext="Simple counting of + / - strand",
            plot=table.plot(self.isoseq3_refine_data_csv_primer),
        )

        self.add_section(
            name="Primers pairs counting",
            anchor="primer-counting",
            description="Count the number of CCS for each couples of primers detected",
            helptext="Number of CCS per couple of primers detected",
            plot=table.plot(self.isoseq3_refine_data_csv_strand),
        )

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
            data_csv_strand = data['primer'].value_counts().to_dict()
            data_csv_primer = data['strand'].value_counts().to_dict()

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
            'format': '{:,}',
        }
        headers["num_reads_flnc"] = {
            "title": "Number of non chimeric full length reads",
            "description": "Number of non chimeric CCS where both primers have been detected",
            "suffix": " flnc",
            'format': '{:,}',
        }
        headers["num_reads_flnc_polya"] = {
            "title": "Number of poly(A) free non chimeric full length reads",
            "description": "Number of non chimeric CCS where both primers have been detected and the poly(A) tail has been removed",
            "suffix": " flnc_polya",
            'format': '{:,}',
        }
        self.general_stats_addcols(gstats_data, headers)
