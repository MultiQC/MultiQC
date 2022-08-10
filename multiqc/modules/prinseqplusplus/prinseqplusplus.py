from multiqc.modules.base_module import BaseMultiqcModule
import logging
import re
from collections import OrderedDict
from multiqc import config
from multiqc.plots import beeswarm

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="PRINSEQ++",
            anchor="prinseqplusplus",
            href="https://github.com/Adrian-Cantu/PRINSEQ-plus-plus",
            info="PRINSEQ++ is a C++ implementation of the prinseq-lite.pl program. It can be used to filter, reformat or trim genomic and metagenomic sequence data.",
            doi="10.7287/peerj.preprints.27553v1",
        )

        # Find and load reports
        self.prinseqplusplus_data = dict()

        # Find all files for prinseqplusplus
        for f in self.find_log_files("prinseqplusplus", filehandles=True):
            self.parse_logs(f)

        self.prinseqplusplus_data = self.ignore_samples(self.prinseqplusplus_data)

        if len(self.prinseqplusplus_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.prinseqplusplus_data)))

        # Write data to file
        self.write_data_file(self.prinseqplusplus_data, "prinseqplusplus")

        self.prinseqplusplus_general_stats()
        self.prinseqplusplus_beeswarm_plot()

    def parse_logs(self, logfile):
        """Parsing Logs."""
        s_name = logfile["fn"]
        s_name = self.clean_s_name(s_name, logfile)

        if self.prinseqplusplus_data.get(s_name) is not None:
            log.warn("Duplicate sample name found! Overwriting: {}".format(s_name))

        self.prinseqplusplus_data[s_name] = {}
        file_content = logfile["f"]
        for l in file_content:
            ## Find line after loading reads, and remove suffixes for sample name

            if "-min_len" in l:
                self.prinseqplusplus_data[s_name]["min_len"] = {}
                self.prinseqplusplus_data[s_name]["min_len"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-max_len" in l:
                self.prinseqplusplus_data[s_name]["max_len"] = {}
                self.prinseqplusplus_data[s_name]["max_len"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-min_cg" in l:  ## This seems to be a typo in the log file of PRINSEQ++, should be gc
                self.prinseqplusplus_data[s_name]["min_gc"] = {}
                self.prinseqplusplus_data[s_name]["min_gc"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-max_cg" in l:
                self.prinseqplusplus_data[s_name]["max_gc"] = {}
                self.prinseqplusplus_data[s_name]["max_gc"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-min_qual_score" in l:
                self.prinseqplusplus_data[s_name]["min_qual_score"] = {}
                self.prinseqplusplus_data[s_name]["min_qual_score"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-min_qual_mean" in l:
                self.prinseqplusplus_data[s_name]["min_qual_mean"] = {}
                self.prinseqplusplus_data[s_name]["min_qual_mean"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-ns_max_n" in l:
                self.prinseqplusplus_data[s_name]["ns_max_n"] = {}
                self.prinseqplusplus_data[s_name]["ns_max_n"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-noiupac" in l:
                self.prinseqplusplus_data[s_name]["noiupac"] = {}
                self.prinseqplusplus_data[s_name]["noiupac"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-derep" in l:
                self.prinseqplusplus_data[s_name]["derep"] = {}
                self.prinseqplusplus_data[s_name]["derep"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-lc_entropy" in l:
                self.prinseqplusplus_data[s_name]["lc_entropy"] = {}
                self.prinseqplusplus_data[s_name]["lc_entropy"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-lc_dust" in l:
                self.prinseqplusplus_data[s_name]["lc_dust"] = {}
                self.prinseqplusplus_data[s_name]["lc_dust"] = int(l.lstrip().rstrip().split(" ")[0])

    def prinseqplusplus_general_stats(self):
        """PRINSEQ++ General Stats Table"""
        headers = OrderedDict()
        headers["min_len"] = {
            "title": "Min. Length Filtered ({})".format(config.read_count_prefix),
            "description": "Reads falling below min length filter ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        headers["max_len"] = {
            "title": "Max. Length Filtered ({})".format(config.read_count_prefix),
            "description": "Reads execeed max length filter ({})".format(config.read_count_prefix),
            "scale": "Purples",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        headers["min_gc"] = {
            "title": "Min. GC Filtered ({})".format(config.read_count_prefix),
            "description": "Reads execeed min GC filter ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        headers["max_gc"] = {
            "title": "Max. GC Filtered ({})".format(config.read_count_prefix),
            "description": "Reads exceeding min GC filter ({})".format(config.read_count_prefix),
            "scale": "Purples",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        headers["min_qual_score"] = {
            "title": "Min. Qual. Score Filtered ({})".format(config.read_count_prefix),
            "description": "Reads falling below min quality score ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        headers["min_qual_mean"] = {
            "title": "Min. Qual. Mean Filtered ({})".format(config.read_count_prefix),
            "description": "Reads falling below min mean quality score ({})".format(config.read_count_prefix),
            "scale": "Purples",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        headers["ns_max_n"] = {
            "title": "Max N Filtered ({})".format(config.read_count_prefix),
            "description": "Reads with more than max Ns ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        headers["noiupac"] = {
            "title": "IUPAC Read Filtered ({})".format(config.read_count_prefix),
            "description": "Reads with bases other than ACTGN ({})".format(config.read_count_prefix),
            "scale": "Purples",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        headers["derep"] = {
            "title": "Duplicated Reads Filtered ({})".format(config.read_count_prefix),
            "description": "Duplicated reads removed ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        headers["lc_entropy"] = {
            "title": "Entropy Filtered Reads ({})".format(config.read_count_prefix),
            "description": "Entropy filtered reads removed ({})".format(config.read_count_prefix),
            "scale": "Purples",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }
        headers["lc_dust"] = {
            "title": "DUST Filtered Reads ({})".format(config.read_count_prefix),
            "description": "DUST_score filtered reads removed ({})".format(config.read_count_prefix),
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x * config.read_count_multiplier,
        }

        self.general_stats_addcols(self.prinseqplusplus_data, headers)

    def prinseqplusplus_beeswarm_plot(self):
        """Beeswarm plot of all possible filtering results"""
        headers = OrderedDict()

        reads = {
            "min": 0,
            "modify": lambda x: float(x) * config.read_count_multiplier,
            "suffix": "{} reads".format(config.read_count_prefix),
            "decimalPlaces": 0,
            "shared_key": "read_count",
        }
        headers["min_len"] = dict(reads, title="min_len")
        headers["max_len"] = dict(reads, title="max_len")
        headers["min_gc"] = dict(reads, title="min_gc")
        headers["max_gc"] = dict(reads, title="max_gc")
        headers["lc_entropy"] = dict(reads, title="lc_entropy")
        headers["min_qual_score"] = dict(reads, title="min_qual_score")
        headers["min_qual_mean"] = dict(reads, title="min_qual_mean")
        headers["ns_max_n"] = dict(reads, title="ns_max_n")
        headers["noiupac"] = dict(reads, title="noiupac")
        headers["derep"] = dict(reads, title="derep")
        headers["lc_entropy"] = dict(reads, title="lc_entropy")
        headers["lc_dust"] = dict(reads, title="lc_dust")

        self.add_section(
            name="Filtered Reads",
            anchor="prinseqplusplus",
            description="Shows the number of reads removed by the various PRINSEQ++ filter options",
            plot=beeswarm.plot(self.prinseqplusplus_data, headers, {"id": "prinseplusplus"}),
        )
