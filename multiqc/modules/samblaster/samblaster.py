import logging
import os
import re

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)

VERSION_REGEX = r"Version\ (\d{1}.\d+.\d+)"


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Samblaster",
            anchor="samblaster",
            href="https://github.com/GregoryFaust/samblaster",
            info="Marks duplicates and extracts discordant and split reads from sam files.",
            doi="10.1093/bioinformatics/btu314",
        )

        self.samblaster_data = dict()
        for f in self.find_log_files("samblaster", filehandles=True):
            self.parse_samblaster(f)

        # Filter to strip out ignored sample names
        self.samblaster_data = self.ignore_samples(self.samblaster_data)

        if len(self.samblaster_data) == 0:
            raise ModuleNoSamplesFound

        headers = {
            "pct_dups": {
                "title": "% Dups",
                "description": "Percent Duplication",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "OrRd",
            }
        }

        self.general_stats_addcols(self.samblaster_data, headers)

        # Write parsed report data to a file
        self.write_data_file(self.samblaster_data, "multiqc_samblaster")

        log.info(f"Found {len(self.samblaster_data)} reports")

        self.add_barplot()

    def add_barplot(self):
        """Generate the Samblaster bar plot."""
        cats = {
            "n_nondups": {"name": "Non-duplicates"},
            "n_dups": {"name": "Duplicates"},
        }

        pconfig = {
            "id": "samblaster_duplicates",
            "title": "Samblaster: Number of duplicate reads",
            "ylab": "Number of reads",
        }
        self.add_section(plot=bargraph.plot(self.samblaster_data, cats, pconfig))

    def parse_samblaster(self, f):
        """Go through log file looking for samblaster output.
        Grab the name from the RG tag of the preceding bwa command"""

        # Should capture the following:
        # samblaster: Marked     1134898 of   43791982 (2.592%) total read ids as duplicates using 753336k \
        # memory in 1M1S(60.884S) CPU seconds and 3M53S(233S) wall time.
        dups_regex = (
            r"samblaster: (Removed|Marked)\s+(\d+)\s+of\s+(\d+) \((\d+.\d+)%\)\s*(total)?\s*read ids as duplicates"
        )

        input_file_regex = r"samblaster: Opening (\S+) for read."
        rgtag_name_regex = r"\\\\tID:(\S*?)\\\\t"
        data = {}
        s_name = None
        version = None
        fh = f["f"]
        for line in fh:
            match = re.search(VERSION_REGEX, line)
            if match is not None:
                version = match.group(1)
            # try to find name from RG-tag. If bwa mem is used upstream samblaster with pipes, then the bwa mem command
            # including the read group will be written in the log
            match = re.search(rgtag_name_regex, line)
            if match:
                s_name = self.clean_s_name(match.group(1), f)

            # try to find name from the input file name, if used
            match = re.search(input_file_regex, line)
            if match:
                basefn = os.path.basename(match.group(1))
                fname, ext = os.path.splitext(basefn)
                # if it's stdin, then try bwa RG-tag instead
                if fname != "stdin":
                    s_name = self.clean_s_name(fname, f)

            match = re.search(dups_regex, line)
            if match:
                data["n_dups"] = int(match.group(2))
                data["n_tot"] = int(match.group(3))
                data["n_nondups"] = data["n_tot"] - data["n_dups"]
                data["pct_dups"] = float(match.group(4))

        if s_name is None:
            s_name = f["s_name"]

        if version is not None:
            self.add_software_version(version, s_name)

        if len(data) > 0:
            if s_name in self.samblaster_data:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
            self.add_data_source(f, s_name)
            self.samblaster_data[s_name] = data
