""" MultiQC submodule to parse output from RSeQC geneBody_coverage.py
http://rseqc.sourceforge.net/#genebody-coverage-py """

import logging

from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find RSeQC gene_body_coverage reports and parse their data"""

    # Set up vars
    self.gene_body_cov_hist_counts = dict()
    self.gene_body_cov_hist_percent = dict()

    # TODO - Do separate parsing step to find skewness values
    # and add these to the general stats table?

    # Go through files and parse data
    for f in self.find_log_files("rseqc/gene_body_coverage"):
        # geneBodyCoverage.py
        # RSeQC >= v2.4
        # NB: Capitilisation
        if f["f"].startswith("Percentile"):
            keys = []
            nrows = 0
            for line in f["f"].splitlines():
                s = line.split("\t")
                # Check that this is the right file type (detection by capitilisation is pretty weak!)
                if len(s) < 3:
                    break
                if len(keys) == 0:
                    keys = s[1:]
                else:
                    nrows += 1
                    s_name = self.clean_s_name(s[0], f)
                    if s_name in self.gene_body_cov_hist_counts:
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    self.add_data_source(f, s_name, section="gene_body_coverage")
                    self.gene_body_cov_hist_counts[s_name] = dict()
                    for k, var in enumerate(s[1:]):
                        self.gene_body_cov_hist_counts[s_name][int(keys[k])] = float(var)
            if nrows == 0:
                log.warning(f"Empty geneBodyCoverage file found: {f['fn']}")

        # geneBodyCoverage2.py
        #   AND
        # geneBodyCoverage.py
        # RSeQC < v2.4
        if f["f"].startswith("Total reads") or f["f"].startswith("percentile"):
            if f["s_name"] in self.gene_body_cov_hist_counts:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            self.add_data_source(f, section="gene_body_coverage")
            self.gene_body_cov_hist_counts[f["s_name"]] = dict()
            nrows = 0
            for line in f["f"].splitlines():
                s = line.split("\t")
                # Check that this is the right file type (detection by capitilisation is pretty weak!)
                if len(s) > 2:
                    break
                try:
                    nrows += 1
                    self.gene_body_cov_hist_counts[f["s_name"]][int(s[0])] = float(s[1])
                except (IndexError, ValueError):
                    # Header lines don't use tabs, so get IndexError.
                    # Other unexpected weird stuff will not play nicely with int() and float()
                    pass
            if nrows == 0:
                del self.gene_body_cov_hist_counts[f["s_name"]]
                log.warning(f"Empty geneBodyCoverage file found: {f['fn']}")

    # Filter to strip out ignored sample names
    self.gene_body_cov_hist_counts = self.ignore_samples(self.gene_body_cov_hist_counts)

    if len(self.gene_body_cov_hist_counts) == 0:
        return 0

    # Write data to file
    self.write_data_file(self.gene_body_cov_hist_counts, "rseqc_gene_body_cov")

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    self.add_software_version(None)

    # Make a normalised coverage for plotting using the formula (cov - min_cov) / (max_cov - min_cov)
    for s_name in self.gene_body_cov_hist_counts:
        self.gene_body_cov_hist_percent[s_name] = dict()
        # min_cov and max_cov are required to compute the normalized coverage
        min_cov = min(self.gene_body_cov_hist_counts[s_name].values())
        max_cov = max(self.gene_body_cov_hist_counts[s_name].values())
        for k, v in self.gene_body_cov_hist_counts[s_name].items():
            self.gene_body_cov_hist_percent[s_name][k] = (v - min_cov) / (max_cov - min_cov) * 100.0

    # Add line graph to section
    pconfig = {
        "id": "rseqc_gene_body_coverage_plot",
        "title": "RSeQC: Gene Body Coverage",
        "ylab": "Coverage",
        "xlab": "Gene Body Percentile (5' -> 3')",
        "xmin": 0,
        "xmax": 100,
        "tt_label": "<strong>{point.x}% from 5'</strong>: {point.y:.2f}",
        "data_labels": [
            {"name": "Percentages", "ylab": "Percentage Coverage"},
            {"name": "Counts", "ylab": "Coverage"},
        ],
    }
    self.add_section(
        name="Gene Body Coverage",
        anchor="rseqc-gene_body_coverage",
        description='<a href="http://rseqc.sourceforge.net/#genebody-coverage-py" target="_blank">Gene Body Coverage</a>'
        " calculates read coverage over gene bodies."
        " This is used to check if reads coverage is uniform and"
        " if there is any 5' or 3' bias.",
        plot=linegraph.plot([self.gene_body_cov_hist_percent, self.gene_body_cov_hist_counts], pconfig),
    )

    # Return number of samples found
    return len(self.gene_body_cov_hist_counts)
