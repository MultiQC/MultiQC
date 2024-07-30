import logging
import re

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

log = logging.getLogger(__name__)

VERSION_REGEX = r"diamond v([\d\.]+)"


class MultiqcModule(BaseMultiqcModule):
    """
    The module takes summary statistics from the `diamond.log` file (`--log` option). It parses and reports
    the number of sequences aligned and displays them in the General Stats table.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="DIAMOND",
            anchor="diamond",
            href="https://github.com/bbuchfink/diamond",
            info="Sequence aligner for protein and translated DNA searches, a drop-in replacement for the NCBI BLAST",
            extra="""
            Key features are:
            - Pairwise alignment of proteins and translated DNA at 100x-10,000x speed of BLAST.
            - Frameshift alignments for long read analysis.
            - Low resource requirements and suitable for running on standard desktops or laptops.
            - Various output formats, including BLAST pairwise, tabular and XML, as well as taxonomic classification.
            """,
            doi="10.1038/s41592-021-01101-x",
        )

        # Find and load any DIAMOND reports
        self.diamond_data = dict()

        for f in self.find_log_files("diamond", filehandles=True):
            self.parse_logs(f)

        # Filter to strip out ignored sample names
        self.diamond_data = self.ignore_samples(self.diamond_data)

        if len(self.diamond_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.diamond_data)} reports")

        # Write parsed report data to file
        self.write_data_file(self.diamond_data, "diamond")
        self.diamond_general_stats()

    def parse_logs(self, f):
        """Parsing logs"""
        s_name = self.clean_s_name(f["root"], f)
        for line in f["f"]:
            # Try to get the sample name from --out, otherwise --query, fallback to directory name
            if "diamond blastx" in line and "--out" in line:
                s_name = line.split("--out ")[1].split(" ")[0]
                s_name = self.clean_s_name(s_name, f)
            elif "diamond blastx" in line and "--query" in line:
                s_name = line.split("--query ")[1].split(" ")[0]
                s_name = self.clean_s_name(s_name, f)

            # Get version
            version_match = re.search(VERSION_REGEX, line)
            if version_match:
                self.add_software_version(version_match.group(1), s_name)

            if "queries aligned" in line:
                self.add_data_source(f)
                if s_name in self.diamond_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.diamond_data[s_name] = {"queries_aligned": int(line.split(" ")[0])}

    def diamond_general_stats(self):
        """Diamond General Stats Table"""
        headers = {
            "queries_aligned": {
                "title": "Queries aligned",
                "description": "number of queries aligned",
                "scale": "YlGn",
                "format": "{:,.0f}",
            }
        }
        self.general_stats_addcols(self.diamond_data, headers)
