import json
import logging
from pathlib import Path

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module recognizes files with the `*_gopeaks.json` suffix (which is the default behavior), and will report
    the number of peaks called per sample via the general table and the bar plot.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="GoPeaks",
            anchor="gopeaks",
            href="https://github.com/maxsonBraunLab/gopeaks",
            info="Calls peaks in CUT&TAG/CUT&RUN datasets.",
            extra="""
            Gopeaks uses a binomial distribution to model the read counts in sliding windows across 
            the genome and calculate peak regions that are enriched over the background.""",
            doi="10.1186/s13059-022-02707-w",
        )

        # data vars ---------------------------------------------------------------------

        # Find and load any gopeaks reports
        self.gopeaks_data = dict()
        for f in self.find_log_files("gopeaks", filehandles=True):
            parsed = self.parse_gopeaks_log(f)

            if parsed is not None:
                if f["s_name"] in self.gopeaks_data:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {f['s_name']}")

                if "gopeaks_version" in parsed:
                    self.add_software_version(parsed["gopeaks_version"], f["s_name"])

                self.gopeaks_data[f["s_name"]] = parsed

                # add gopeaks data to multiqc_source.txt
                self.add_data_source(f)

        # filter away samples if MultiQC user does not want them
        self.gopeaks_data = self.ignore_samples(self.gopeaks_data)

        if len(self.gopeaks_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.gopeaks_data)} samples")

        self.write_data_file(self.gopeaks_data, "multiqc_gopeaks")

        # Add sample log info to basic stats table
        self.gopeaks_general_stats_table()

        # Add sample log info to bargraph
        self.gopeaks_bargraph()

    # parsing functions -------------------------------------------------------------

    def parse_gopeaks_log(self, f):
        """
        Read gopeaks json log file and extract number of peaks.
        """
        with open(Path(f["root"]) / Path(f["fn"])) as f:
            sample_log = json.load(f)

        return sample_log

    def gopeaks_general_stats_table(self):
        """
        Put peak counts to the general table.
        """

        headers = {
            "peak_counts": {
                "title": "Peak Counts",
                "description": "Number of peaks per sample",
                "min": 0,
                "scale": "YlGnBu",
                "format": "{:,.0f}",
            }
        }
        self.general_stats_addcols(self.gopeaks_data, headers)

    def gopeaks_bargraph(self):
        """
        Put peak counts to a bargraph.
        """

        cats = {"peak_counts": {"name": "Peak Counts"}}
        config = {
            "id": "GoPeaksBarGraph",
            "title": "GoPeaks: Number of Peaks by Sample",
            "ylab": "Sample",
            "logswitch": True,
            "cpswitch": False,
        }

        self.add_section(
            name="GoPeaks",
            anchor="gopeaks_bargraph",
            description="Number of peaks called by GoPeaks.",
            plot=bargraph.plot(self.gopeaks_data, cats, config),
        )
