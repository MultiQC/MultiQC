import logging
import re

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)

VERSION_REGEX = r"Flexbar - flexible barcode and adapter removal, version ([\d\.]+)"


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Flexbar",
            anchor="flexbar",
            href="https://github.com/seqan/flexbar",
            info="Barcode and adapter removal tool.",
            extra="""
            Flexbar efficiently preprocesses high-throughput sequencing data. It demultiplexes 
            barcoded runs and removes adapter sequences. Moreover, trimming and filtering features are provided.
            Flexbar increases read mapping rates and improves genome as well as transcriptome assemblies.
            """,
            doi="10.1093/bioinformatics/btx330",
        )

        # Parse logs
        self.flexbar_data = dict()
        for f in self.find_log_files("flexbar", filehandles=True):
            self.parse_flexbar(f)
            self.add_data_source(f)

        # Filter to strip out ignored sample names
        self.flexbar_data = self.ignore_samples(self.flexbar_data)

        if len(self.flexbar_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.flexbar_data)} logs")
        self.write_data_file(self.flexbar_data, "multiqc_flexbar")

        # Add drop rate to the general stats table
        headers = {}
        headers["removed_bases_pct"] = {
            "title": "% bp Trimmed",
            "description": "% Total Base Pairs removed",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "YlOrRd",
        }
        self.general_stats_addcols(self.flexbar_data, headers)

        # Make barplot
        self.flexbar_barplot()

    def parse_flexbar(self, f):
        def _save_data(parsed_data):
            if len(parsed_data) > 0:
                # Calculate removed_bases
                if "processed_bases" in parsed_data and "remaining_bases" in parsed_data:
                    parsed_data["removed_bases"] = parsed_data["processed_bases"] - parsed_data["remaining_bases"]
                    try:
                        parsed_data["removed_bases_pct"] = (
                            float(parsed_data["removed_bases"]) / float(parsed_data["processed_bases"])
                        ) * 100.0
                    except ZeroDivisionError:
                        parsed_data["removed_bases_pct"] = 0
                if s_name in self.flexbar_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.flexbar_data[s_name] = parsed_data

        regexes = {
            "output_filename": r"Read file:\s+(.+)$",
            "processed_reads": r"Processed reads\s+(\d+)",
            "skipped_due_to_uncalled_bases": r"skipped due to uncalled bases\s+(\d+)",
            "short_prior_to_adapter_removal": r"short prior to adapter removal\s+(\d+)",
            "finally_skipped_short_reads": r"finally skipped short reads\s+(\d+)",
            "discarded_reads_overall": r"Discarded reads overall\s+(\d+)",
            "remaining_reads": r"Remaining reads\s+(\d+)",
            "processed_bases": r"Processed bases:?\s+(\d+)",
            "remaining_bases": r"Remaining bases:?\s+(\d+)",
        }
        s_name = f["s_name"]
        parsed_data = dict()
        version = None
        for line in f["f"]:
            # The version appears in the before the sample name in the log so we
            # assign it to a variable for now.
            version_match = re.search(VERSION_REGEX, line)
            if version_match:
                version = version_match.group(1)

            for k, r in regexes.items():
                match = re.search(r, line)
                if match:
                    if k == "output_filename":
                        s_name = self.clean_s_name(match.group(1), f)
                        if version is not None:
                            self.add_software_version(version, s_name)
                    else:
                        parsed_data[k] = int(match.group(1))

            # End of log output. Save and reset in case of more logs.
            if "Flexbar completed" in line:
                _save_data(parsed_data)
                s_name = f["s_name"]
                parsed_data = dict()

        # Pick up any partial logs
        _save_data(parsed_data)

    def flexbar_barplot(self):
        # Specify the order of the different possible categories
        keys = {
            "remaining_reads": {"color": "#437bb1", "name": "Remaining reads"},
            "skipped_due_to_uncalled_bases": {"color": "#e63491", "name": "Skipped due to uncalled bases"},
            "short_prior_to_adapter_removal": {"color": "#b1084c", "name": "Short prior to adapter removal"},
            "finally_skipped_short_reads": {"color": "#7f0000", "name": "Finally skipped short reads"},
        }

        # Config for the plot
        pconfig = {
            "id": "flexbar_plot",
            "title": "Flexbar: Processed Reads",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
            "hide_empty": False,
        }

        self.add_section(plot=bargraph.plot(self.flexbar_data, keys, pconfig))
