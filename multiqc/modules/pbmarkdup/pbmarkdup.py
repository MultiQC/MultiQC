import logging
import re

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module adds the **% Unique Molecules** and **%Duplicate Reads** (hidden) to the General Statistics
    table.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="pbmarkdup",
            anchor="pbmarkdup",
            href="https://github.com/PacificBiosciences/pbmarkdup",
            info="Takes one or multiple sequencing chips of an amplified libray as HiFi reads and marks or "
            "removes duplicates.",
            # Can't find a DOI // doi=
        )

        self.pbmarkdup = dict()

        for logfile in self.find_log_files("pbmarkdup", filehandles=True):
            data = self.parse_logfile(logfile)
            if data:
                # Clean the file name
                s_name = self.clean_s_name(logfile["s_name"], logfile)
                self.pbmarkdup[s_name] = data
                self.add_data_source(logfile)

                # Superfluous function call to confirm that it is used in this module
                # Replace None with actual version if it is available
                self.add_software_version(None, s_name)

        # Filter to strip out ignored sample names
        self.pbmarkdup = self.ignore_samples(self.pbmarkdup)

        # Raise ModuleNoSamplesFound if we did not find any data
        if not self.pbmarkdup:
            raise ModuleNoSamplesFound

        # Write the parsed data to file
        self.write_data_file(self.pbmarkdup, "multiqc_pbmarkdup")

        self.pbmarkdup_add_general_stats()

    def parse_logfile(self, logfile):
        """
        Parse the standard output from pbmarkdup

        This is currently quite ugly, since pbmarkdup does not have an easily
        parsable output format. See this github issue:
        https://github.com/PacificBiosciences/pbbioconda/issues/365

        For now, we just assume that the data fields are always in the same
        order, without any missing values.
        """

        file_content = logfile["f"]

        # The number of spaces in the header can vary, based on the length of
        # the name of the library. We therefore need a regex to match the
        # header.
        header_pattern = "LIBRARY +READS +UNIQUE MOLECULES +DUPLICATE READS"
        pattern = re.compile(header_pattern)

        # The file header
        file_header = next(file_content).strip()

        # Log an error if the header doesn't match the expected pattern
        if not re.match(pattern, file_header):
            fname = logfile["fn"]
            log.error(f"Can't parse file '{fname}', unknown header: '{file_header}'")
            return False

        data = dict()

        # Each parsable line is either for a library, or for 'TOTAL', the sum
        # of all parsed libraries
        for line in file_content:
            # Lines which start with --- denote separators in the file, and do
            # not need to be parsed
            if line.startswith("-"):
                continue

            # Not very nice, we assume that all fields are always present
            lib_name, reads, unique_mol_count, unique_mol_perc, duplicate_count, duplicate_perc = line.split()

            # We are only interested in the counts, not the percentages
            data[lib_name] = {
                "READS": int(reads),
                "UNIQUE MOLECULES": int(unique_mol_count),
                "DUPLICATE READS": int(duplicate_count),
            }

        return data

    def pbmarkdup_add_general_stats(self):
        """Add pbmarkdup duplicates to the general stats table"""

        general_stats_headers = {
            "unique_molecules": {
                "title": "% Unique Molecules",
                "description": "Percentage of unique molecules",
                "suffix": "%",
                "min": 0,
                "max": 100,
                "modify": lambda x: x * 100,
                "scale": "RdYlGn",
            },
            "duplicate_reads": {
                "title": "% Duplicate Reads",
                "description": "Percentage of duplicate reads",
                "suffix": "%",
                "min": 0,
                "max": 100,
                "modify": lambda x: x * 100,
                "scale": "RdYlGn-rev",
                "hidden": True,
            },
        }

        general = dict()

        for sample, sample_stats in self.pbmarkdup.items():
            # We are only interested in the total counts per sample
            total = sample_stats["TOTAL"]

            # Calculate the percentages unique and duplicated
            perc_unique = total["UNIQUE MOLECULES"] / total["READS"]
            perc_duplicate = total["DUPLICATE READS"] / total["READS"]

            general[sample] = {"unique_molecules": perc_unique, "duplicate_reads": perc_duplicate}

        self.general_stats_addcols(general, general_stats_headers)
