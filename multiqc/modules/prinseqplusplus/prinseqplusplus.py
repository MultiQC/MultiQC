import logging

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import violin

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    This module requires that PRINSEQ++ has been run with the flag `-VERBOSE 1`.

    It uses the log file name as the sample name.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="PRINSEQ++",
            anchor="prinseqplusplus",
            href="https://github.com/Adrian-Cantu/PRINSEQ-plus-plus",
            info="C++ implementation of the prinseq-lite.pl program. Filters, reformats, and trims genomic and "
            "metagenomic reads.",
            doi="10.7287/peerj.preprints.27553v1",
        )

        # Find and load reports
        self.prinseqplusplus_data = dict()

        # Find all files for prinseqplusplus
        for f in self.find_log_files("prinseqplusplus", filehandles=True):
            self.parse_logs(f)

        self.prinseqplusplus_data = self.ignore_samples(self.prinseqplusplus_data)

        if len(self.prinseqplusplus_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.prinseqplusplus_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write data to file
        self.write_data_file(self.prinseqplusplus_data, "prinseqplusplus")

        self.prinseqplusplus_general_stats()
        self.prinseqplusplus_violin_plot()

    def parse_logs(self, f):
        """Parsing Logs."""
        s_name = f["s_name"]
        if self.prinseqplusplus_data.get(s_name) is not None:
            log.warn(f"Duplicate sample name found! Overwriting: {s_name}")

        self.prinseqplusplus_data[s_name] = {}
        self.add_data_source(f, s_name=s_name)

        for line in f["f"]:
            # Find line after loading reads, and remove suffixes for sample name
            if "-min_len" in line:
                self.prinseqplusplus_data[s_name]["min_len"] = int(line.lstrip().rstrip().split(" ")[0])
            elif "-max_len" in line:
                self.prinseqplusplus_data[s_name]["max_len"] = int(line.lstrip().rstrip().split(" ")[0])
            elif "-min_cg" in line:  ## This seems to be a typo in the log file of PRINSEQ++, should be gc
                self.prinseqplusplus_data[s_name]["min_gc"] = int(line.lstrip().rstrip().split(" ")[0])
            elif "-max_cg" in line:
                self.prinseqplusplus_data[s_name]["max_gc"] = int(line.lstrip().rstrip().split(" ")[0])
            elif "-min_qual_score" in line:
                self.prinseqplusplus_data[s_name]["min_qual_score"] = int(line.lstrip().rstrip().split(" ")[0])
            elif "-min_qual_mean" in line:
                self.prinseqplusplus_data[s_name]["min_qual_mean"] = int(line.lstrip().rstrip().split(" ")[0])
            elif "-ns_max_n" in line:
                self.prinseqplusplus_data[s_name]["ns_max_n"] = int(line.lstrip().rstrip().split(" ")[0])
            elif "-noiupac" in line:
                self.prinseqplusplus_data[s_name]["noiupac"] = int(line.lstrip().rstrip().split(" ")[0])
            elif "-derep" in line:
                self.prinseqplusplus_data[s_name]["derep"] = int(line.lstrip().rstrip().split(" ")[0])
            elif "-lc_entropy" in line:
                self.prinseqplusplus_data[s_name]["lc_entropy"] = int(line.lstrip().rstrip().split(" ")[0])
            elif "-lc_dust" in line:
                self.prinseqplusplus_data[s_name]["lc_dust"] = int(line.lstrip().rstrip().split(" ")[0])

    def prinseqplusplus_general_stats(self):
        """PRINSEQ++ General Stats Table"""
        data = {}
        for s_name, d in self.prinseqplusplus_data.items():
            data[s_name] = {"prinseqplusplus_total": sum(d.values())}

        self.general_stats_addcols(
            data,
            {
                "prinseqplusplus_total": {
                    "title": f"Filtered Reads ({config.read_count_prefix})",
                    "description": f"Sum of filtered reads ({config.read_count_desc})",
                    "scale": "Oranges",
                    "shared_key": "read_count",
                    "modify": lambda x: x * config.read_count_multiplier,
                }
            },
        )

    def prinseqplusplus_violin_plot(self):
        """Violin plot of all possible filtering results"""
        # This would be nicer as a stacked-bar plot, but as we don't have
        # the total read count it doesn't really make sense.

        reads = {
            "min": 0,
            "modify": lambda x: float(x) * config.read_count_multiplier,
            "suffix": config.read_count_prefix,
            "shared_key": "read_count",
        }
        headers = {
            "min_len": dict(reads, title="min_len"),
            "max_len": dict(reads, title="max_len"),
            "min_gc": dict(reads, title="min_gc"),
            "max_gc": dict(reads, title="max_gc"),
            "lc_entropy": dict(reads, title="lc_entropy"),
            "min_qual_score": dict(reads, title="min_qual_score"),
            "min_qual_mean": dict(reads, title="min_qual_mean"),
            "ns_max_n": dict(reads, title="ns_max_n"),
            "noiupac": dict(reads, title="noiupac"),
            "derep": dict(reads, title="derep"),
            "lc_dust": dict(reads, title="lc_dust"),
        }

        self.add_section(
            name="Filtered Reads",
            anchor="prinseqplusplus-filtered-reads",
            description="Shows the number of reads removed by the various PRINSEQ++ filter options",
            plot=violin.plot(
                self.prinseqplusplus_data,
                headers,
                {
                    "id": "prinseplusplus-filtered-reads-violin",
                    "title": "PRINSEQ++: Filtered Reads",
                },
            ),
        )
