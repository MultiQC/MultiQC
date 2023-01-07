import logging
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
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

    def parse_logs(self, f):
        """Parsing Logs."""
        s_name = f["s_name"]
        if self.prinseqplusplus_data.get(s_name) is not None:
            log.warn("Duplicate sample name found! Overwriting: {}".format(s_name))

        self.prinseqplusplus_data[s_name] = {}
        self.add_data_source(f, s_name=s_name)

        for l in f["f"]:
            ## Find line after loading reads, and remove suffixes for sample name

            if "-min_len" in l:
                self.prinseqplusplus_data[s_name]["min_len"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-max_len" in l:
                self.prinseqplusplus_data[s_name]["max_len"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-min_cg" in l:  ## This seems to be a typo in the log file of PRINSEQ++, should be gc
                self.prinseqplusplus_data[s_name]["min_gc"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-max_cg" in l:
                self.prinseqplusplus_data[s_name]["max_gc"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-min_qual_score" in l:
                self.prinseqplusplus_data[s_name]["min_qual_score"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-min_qual_mean" in l:
                self.prinseqplusplus_data[s_name]["min_qual_mean"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-ns_max_n" in l:
                self.prinseqplusplus_data[s_name]["ns_max_n"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-noiupac" in l:
                self.prinseqplusplus_data[s_name]["noiupac"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-derep" in l:
                self.prinseqplusplus_data[s_name]["derep"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-lc_entropy" in l:
                self.prinseqplusplus_data[s_name]["lc_entropy"] = int(l.lstrip().rstrip().split(" ")[0])
            elif "-lc_dust" in l:
                self.prinseqplusplus_data[s_name]["lc_dust"] = int(l.lstrip().rstrip().split(" ")[0])

    def prinseqplusplus_general_stats(self):
        """PRINSEQ++ General Stats Table"""
        data = {}
        for s_name, d in self.prinseqplusplus_data.items():
            data[s_name] = {"prinseqplusplus_total": sum(d.values())}

        self.general_stats_addcols(
            data,
            {
                "prinseqplusplus_total": {
                    "title": "Filtered Reads ({})".format(config.read_count_prefix),
                    "description": "Sum of filtered reads ({})".format(config.read_count_desc),
                    "scale": "Oranges",
                    "shared_key": "read_count",
                    "modify": lambda x: x * config.read_count_multiplier,
                }
            },
        )

    def prinseqplusplus_beeswarm_plot(self):
        """Beeswarm plot of all possible filtering results"""
        # This would be nicer as a stacked-bar plot, but as we don't have
        # the total read count it doesn't really make sense.

        reads = {
            "min": 0,
            "modify": lambda x: float(x) * config.read_count_multiplier,
            "suffix": "{} reads".format(config.read_count_prefix),
            "decimalPlaces": 0,
            "shared_key": "read_count",
        }
        headers = OrderedDict()
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
            anchor="prinseqplusplus-filtered-reads",
            description="Shows the number of reads removed by the various PRINSEQ++ filter options",
            plot=beeswarm.plot(
                self.prinseqplusplus_data,
                headers,
                {
                    "id": "prinseplusplus-filtered-reads-beeswarm",
                },
            ),
        )
