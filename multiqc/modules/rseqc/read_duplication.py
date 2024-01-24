""" MultiQC submodule to parse output from RSeQC read_duplication.py
http://rseqc.sourceforge.net/#read-duplication-py """

import logging

from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find RSeQC read_duplication reports and parse their data"""

    # Set up vars
    self.read_dups = dict()

    # Go through files and parse data
    for f in self.find_log_files("rseqc/read_duplication_pos"):
        if f["f"].startswith("Occurrence	UniqReadNumber"):
            if f["s_name"] in self.read_dups:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            self.add_data_source(f, section="read_duplication")
            self.read_dups[f["s_name"]] = dict()
            for line in f["f"].splitlines():
                s = line.split()
                try:
                    if int(s[0]) <= 500:
                        self.read_dups[f["s_name"]][int(s[0])] = int(s[1])
                except Exception:
                    pass

    # Filter to strip out ignored sample names
    self.read_dups = self.ignore_samples(self.read_dups)

    if len(self.read_dups) == 0:
        return 0

    # Write data to file
    self.write_data_file(self.read_dups, "rseqc_read_dups")

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    self.add_software_version(None)

    # Add line graph to section
    pconfig = {
        "smooth_points": 200,
        "id": "rseqc_read_dups_plot",
        "title": "RSeQC: Read Duplication",
        "ylab": "Number of Reads (log10)",
        "xlab": "Occurrence of read",
        "yLog": True,
        "tt_label": "<strong>{point.x} occurrences</strong>: {point.y} reads",
    }

    self.add_section(
        name="Read Duplication",
        anchor="rseqc-read_dups",
        description='<a href="http://rseqc.sourceforge.net/#read-duplication-py" target="_blank">read_duplication.py</a>'
        " calculates how many alignment positions have a certain number of exact duplicates."
        " Note - plot truncated at 500 occurrences and binned.</p>",
        plot=linegraph.plot(self.read_dups, pconfig),
    )

    # Return number of samples found
    return len(self.read_dups)
