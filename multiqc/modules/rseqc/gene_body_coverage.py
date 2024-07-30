"""MultiQC submodule to parse output from RSeQC geneBody_coverage.py
http://rseqc.sourceforge.net/#genebody-coverage-py"""

import logging
from typing import Dict, List

from multiqc import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


def parse_reports(module: BaseMultiqcModule) -> int:
    """Find RSeQC gene_body_coverage reports and parse their data"""

    # Set up vars
    gene_body_cov_hist_counts: Dict = dict()
    gene_body_cov_hist_percent: Dict = dict()

    # TODO - Do separate parsing step to find skewness values
    # and add these to the general stats table?

    # Go through files and parse data
    for f in module.find_log_files("rseqc/gene_body_coverage"):
        # geneBodyCoverage.py
        # RSeQC >= v2.4
        # NB: Capitalisation
        if f["f"].startswith("Percentile"):
            keys: List = []
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
                    s_name = module.clean_s_name(s[0], f)
                    if s_name in gene_body_cov_hist_counts:
                        log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                    module.add_data_source(f, s_name, section="gene_body_coverage")
                    gene_body_cov_hist_counts[s_name] = dict()
                    for k, var in enumerate(s[1:]):
                        gene_body_cov_hist_counts[s_name][int(keys[k])] = float(var)
            if nrows == 0:
                log.warning(f"Empty geneBodyCoverage file found: {f['fn']}")

        # geneBodyCoverage2.py
        #   AND
        # geneBodyCoverage.py
        # RSeQC < v2.4
        if f["f"].startswith("Total reads") or f["f"].startswith("percentile"):
            if f["s_name"] in gene_body_cov_hist_counts:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            module.add_data_source(f, section="gene_body_coverage")
            gene_body_cov_hist_counts[f["s_name"]] = dict()
            nrows = 0
            for line in f["f"].splitlines():
                s = line.split("\t")
                # Check that this is the right file type (detection by capitilisation is pretty weak!)
                if len(s) > 2:
                    break
                try:
                    nrows += 1
                    gene_body_cov_hist_counts[f["s_name"]][int(s[0])] = float(s[1])
                except (IndexError, ValueError):
                    # Header lines don't use tabs, so get IndexError.
                    # Other unexpected weird stuff will not play nicely with int() and float()
                    pass
            if nrows == 0:
                del gene_body_cov_hist_counts[f["s_name"]]
                log.warning(f"Empty geneBodyCoverage file found: {f['fn']}")

    # Filter to strip out ignored sample names
    gene_body_cov_hist_counts = module.ignore_samples(gene_body_cov_hist_counts)

    if len(gene_body_cov_hist_counts) == 0:
        return 0

    # Write data to file
    module.write_data_file(gene_body_cov_hist_counts, "rseqc_gene_body_cov")

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Make a normalised coverage for plotting using the formula (cov - min_cov) / (max_cov - min_cov)
    for s_name in gene_body_cov_hist_counts:
        gene_body_cov_hist_percent[s_name] = dict()
        # min_cov and max_cov are required to compute the normalized coverage
        min_cov = min(gene_body_cov_hist_counts[s_name].values())
        max_cov = max(gene_body_cov_hist_counts[s_name].values())
        for k, v in gene_body_cov_hist_counts[s_name].items():
            gene_body_cov_hist_percent[s_name][k] = (v - min_cov) / (max_cov - min_cov) * 100.0

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
    module.add_section(
        name="Gene Body Coverage",
        anchor="rseqc-gene_body_coverage",
        description='<a href="http://rseqc.sourceforge.net/#genebody-coverage-py" target="_blank">Gene Body Coverage</a>'
        " calculates read coverage over gene bodies."
        " This is used to check if reads coverage is uniform and"
        " if there is any 5' or 3' bias.",
        plot=linegraph.plot([gene_body_cov_hist_percent, gene_body_cov_hist_counts], pconfig),
    )

    # Return number of samples found
    return len(gene_body_cov_hist_counts)
