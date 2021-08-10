#!/usr/bin/env python

""" MultiQC module to parse output from Nextclade """

from __future__ import print_function
from collections import OrderedDict
import logging
from multiqc.modules import nextclade
import jinja2
import csv

from multiqc import config
from multiqc.utils import mqc_colour
from multiqc.plots import bargraph, table
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Nextclade",
            anchor="nextclade",
            href="https://github.com/nextstrain/nextclade",
            info="does viral genome alignment, clade assignment, mutation calling, and quality checks",
        )

        # Parse logs
        self.nextclade_data = dict()
        for f in self.find_log_files("nextclade", filehandles=True):
            self.parse_nextclade_log(f)
            self.add_data_source(f)

        # Remove ignored samples
        self.nextclade_data = self.ignore_samples(self.nextclade_data)

        # Check if we found any samples
        if len(self.nextclade_data) == 0:
            raise UserWarning

        log.info("Found {} samples".format(len(self.nextclade_data)))
        self.write_data_file(self.nextclade_data, "multiqc_nextclade")

        # Add data to general statistics table
        self.nextclade_general_stats_table()

        # Create the Nextclade table
        self.add_section(
            name="Run table", anchor="nextclade-run", description="", helptext="", plot=self.nextclade_table()
        )

    def parse_nextclade_log(self, f):

        for sample in csv.DictReader(f["f"], delimiter=";"):
            try:
                s_name = sample.pop("seqName")
                s_name = self.clean_s_name(s_name, f)
                if s_name in self.nextclade_data:
                    log.error(f"Found duplicate sample '{s_name}' - overwriting")

                # Clean column names: make lower case and
                # replace '.' with '_'
                for key in list(sample.keys()):
                    new_key = key.lower()
                    if "." in new_key:
                        new_key = "_".join(new_key.split("."))
                    sample[new_key] = sample.pop(key)

                for key in list(sample.keys()):
                    new_key = key.lower()
                    sample[new_key] = sample.pop(key)

                self.nextclade_data[s_name] = sample
            except KeyError:
                log.error(f"File '{f['fn']}' could not be parsed - no sequence name found")

    def nextclade_general_stats_table(self):
        """Takes the parsed sample data and adds the sequences' clade to the general stats table"""

        headers = OrderedDict()
        headers["clade"] = {
            "title": "Clade",
            "description": "Clade",
        }
        self.general_stats_addcols(self.nextclade_data, headers)

    def nextclade_table(self):
        """Creates the table of data for the samples"""

        headers = OrderedDict()
        headers["clade"] = {
            "title": "Clade",
            "description": "The inferred clade from the input sequence and reference tree",
        }
        headers["qc_overallscore"] = {
            "title": "QC overall score",
            "description": "Summarizing score of all QC test",
            "hidden": True,
        }

        headers["qc_overallstatus"] = {
            "title": "QC overall status",
            "description": "Summarizing status of all QC tests",
        }
        headers["totalgaps"] = {
            "title": "Total gaps",
            "description": "Total number of detected nucleotide gaps",
            "min": 0,
            "scale": "RdBl-rev",
        }

        headers["totalinsertions"] = {
            "title": "Total insertions",
            "description": "Total number of detected nucleotide insertions",
            "min": 0,
            "scale": "RdYlGr-rev",
        }

        headers["totalmissing"] = {
            "title": "Total missing",
            "description": "Total number of detected missing nucleotides",
            "min": 0,
            "scale": "RdYlBl-rev",
        }

        headers["totalmutations"] = {
            "title": "Total mutations",
            "description": "Total number of detected nucleotide mutations",
            "min": 0,
            "scale": "RdBl-rev",
        }

        headers["totalnonacgtns"] = {
            "title": "Total non ACGTNs",
            "description": "Total number of detected ambiguous nucleotides",
            "min": 0,
            "scale": "RdYlGr-rev",
        }
        headers["totalpcrprimerchanges"] = {
            "title": "Total PCR primer changes",
            "description": "Total number of nucleotide mutations detected in PCR primer regions",
            "min": 0,
            "scale": "RdYlBl-rev",
        }
        headers["totalaminoacidsubstitutions"] = {
            "title": "Total amino acid substitutions",
            "description": "Total number of detected aminoacid substitutions",
            "min": 0,
            "scale": "RdBl-rev",
        }
        headers["totalaminoaciddeletions"] = {
            "title": "Total amino acid deletions",
            "description": "Total number of detected aminoacid substitutions",
            "min": 0,
            "scale": "RdYlGr-rev",
        }
        headers["alignmentscore"] = {
            "title": "Alignment score",
            "description": "Indicates to what degree the input sequence and the reference sequence correspond",
            "min": 0,
            "scale": "RdYlBl-rev",
        }
        headers["alignmentstart"] = {
            "title": "Alignment start",
            "description": "Beginning of the sequenced region",
            "hidden": True,
        }
        headers["alignmentend"] = {
            "title": "Alignment end",
            "description": "End of the sequenced region",
            "hidden": True,
        }
        headers["qc_missingdata_missingdatathreshold"] = {
            "title": "QC missing data: missing data threshold",
            "description": "The threshold used for the 'Missing data' QC test",
            "scale": "GrBl",
            "hidden": True,
        }
        headers["qc_missingdata_score"] = {
            "title": "QC missingdata: score",
            "description": "Score from the 'Missing data' QC test",
            "scale": "RdYlGr-rev",
            "hidden": True,
        }
        headers["qc_missingdata_status"] = {
            "title": "QC missingdata: status",
            "description": "Status for 'Missing data' QC test",
        }
        headers["qc_missingdata_totalmissing"] = {
            "title": "QC missingdata: total missing",
            "description": "Total number of missing nucleotides used in 'Missing data' QC test",
            "scale": "RdYlBl-rev",
            "hidden": True,
        }
        headers["qc_mixedsites_mixedsitesthreshold"] = {
            "title": "QC mixedsites: mixed sites threshold",
            "description": "Threshold used for 'Mixed sites' QC test",
            "scale": "GrBl",
            "hidden": True,
        }
        headers["qc_mixedsites_score"] = {
            "title": "QC mixed sites: score",
            "description": "Score from the 'Missing data' QC testScore from the 'Missing data' QC test",
            "scale": "RdYlGr-rev",
            "hidden": True,
        }
        headers["qc_mixedsites_status"] = {
            "title": "QC mixedsites: status",
            "description": "Status for 'Missing data' QC test",
        }
        headers["qc_mixedsites_totalmixedsites"] = {
            "title": "QC mixedsites: total mixed sites",
            "description": "Total number of ambiguous nucleotides used for 'Mixed sites' QC test",
            "scale": "RdYlBl-rev",
            "hidden": True,
        }
        headers["qc_privatemutations_cutoff"] = {
            "title": "QC private mutations: cutoff",
            "description": "Cutoff parameter used for 'Private mutations' QC test",
            "scale": "GrBl",
            "hidden": True,
        }
        headers["qc_privatemutations_excess"] = {
            "title": "QC private mutations: excess",
            "description": "Excess parameter used for 'Private mutations' QC test",
            "scale": "OrRd",
            "hidden": True,
        }
        headers["qc_privatemutations_score"] = {
            "title": "QC private mutations: score",
            "description": "Score for 'Private mutations' QC rule",
            "scale": "RdYlGr-rev",
            "hidden": True,
        }
        headers["qc_privatemutations_status"] = {
            "title": "QC private mutations: status",
            "description": "Status for 'Private mutations' QC rule",
        }
        headers["qc_privatemutations_total"] = {
            "title": "QC private mutations: total",
            "description": "Total number of private mutations used for 'Private mutations' QC rule",
            "scale": "RdYlBl-rev",
            "hidden": True,
        }
        headers["qc_snpclusters_clusteredsnps"] = {
            "title": "QC SNP clusters: clustered SNPs",
            "description": "Clustered SNP detected for 'SNP clusters' QC test",
            "hidden": True,
        }
        headers["qc_snpclusters_score"] = {
            "title": "QC SNP clusters: score",
            "description": "Score for 'SNP clusters' QC test",
            "scale": "RdYlGr-rev",
            "hidden": True,
        }
        headers["qc_snpclusters_status"] = {
            "title": "QC SNP clusters: status",
            "description": "Status for 'SNP clusters' QC test",
        }
        headers["qc_snpclusters_totalsnps"] = {
            "title": "QC SNP clusters: total SNPs",
            "description": "Total number of SNPs for 'SNP clusters' QC test",
            "scale": "RdYlBl-rev",
            "hidden": True,
        }
        headers["errors"] = {"title": "errors", "description": "List of errors that occurred during processing"}

        # Main table config
        table_config = {
            "namespace": "Nextclade",
            "id": "nextclade_run_table",
            "table_title": "Nextclade Run details",
        }

        return table.plot(self.nextclade_data, headers, table_config)
