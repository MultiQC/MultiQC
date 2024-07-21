import logging
import re

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph

log = logging.getLogger(__name__)

VERSION_REGEX = r"\[M::main\] Version: ([\d\.r\-]+)"


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="HiFiasm",
            anchor="hifiasm",
            href="https://github.com/chhylp123/hifiasm",
            info="Haplotype-resolved assembler for accurate Hifi reads",
            doi="10.1038/s41592-020-01056-5",
        )

        # To store the mod data
        self.hifiasm_data = dict()
        self.parse_hifiasm_log_files()
        self.hifiasm_data = self.ignore_samples(self.hifiasm_data)

        # If we found no data
        if not self.hifiasm_data:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.hifiasm_data)} reports")

        self.write_data_file(self.hifiasm_data, "multiqc_hifiasm_report")
        self.add_sections()

    def parse_hifiasm_log_files(self):
        for f in self.find_log_files("hifiasm", filehandles=True):
            self.add_data_source(f)
            if f["s_name"] in self.hifiasm_data:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            data = self.extract_kmer_graph(f["f"])
            if data:
                self.hifiasm_data[f["s_name"]] = data

            version = self.extract_version(f["f"])
            if version is not None:
                self.add_software_version(version, f["s_name"])

    def add_sections(self):
        # Plot configuration
        config = {
            "id": "hifiasm-kmer-graph",
            "title": "HiFiasm: kmer graph",
            "ylab": "Count of kmer occurrence",
            "xlab": "Kmer occurrence",
            "logswitch": True,
            "logswitch_active": True,
        }

        self.add_section(
            name="HiFiasm kmer graph",
            anchor="hifiasm-kmer-section",
            description="Kmer counts in the input data",
            helptext="""
                The kmer distribution graph for the input data. For homozygous
                samples, there should be one peak around read coverage. For
                heterozygous samples, there should be two peaks, see the
                [HiFiasm documentation](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html#hifiasm-log-interpretation)
                for details.
                """,
            plot=linegraph.plot(self.hifiasm_data, config),
        )

    def extract_version(self, fin):
        """Extract the Hifiasm version from file contents"""
        for line in fin:
            if not line.startswith("[M::main]"):
                continue

            version_match = re.search(VERSION_REGEX, line)
            if version_match:
                return version_match.group(1)
        return None

    def extract_kmer_graph(self, fin):
        """Extract the kmer graph from file in"""
        data = dict()

        found_histogram = False

        for line in fin:
            if line.startswith("[M::ha_hist_line]"):
                found_histogram = True
                spline = line.strip().split()
                # Occurrence of kmer
                occurrence = spline[1][:-1]
                # Special case
                if occurrence == "rest":
                    continue
                # Count of the occurrence, checking for lines with no asterisk before count.
                if "*" in spline[2]:
                    count = int(spline[3])
                else:
                    count = int(spline[2])
                data[int(occurrence)] = count
            # If we are no longer in the histogram
            elif found_histogram:
                return data
