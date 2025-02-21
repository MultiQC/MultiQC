"""MultiQC submodule to parse output from RSeQC infer_experiment.py
http://rseqc.sourceforge.net/#infer-experiment-py"""

import logging
import re
from typing import Dict

from multiqc import BaseMultiqcModule
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


def parse_reports(module: BaseMultiqcModule) -> int:
    """Find RSeQC infer_experiment reports and parse their data"""

    infer_exp: Dict = dict()
    regexes = {
        "pe_sense": r"\"1\+\+,1--,2\+-,2-\+\": (\d\.\d+)",
        "pe_antisense": r"\"1\+-,1-\+,2\+\+,2--\": (\d\.\d+)",
        "se_sense": r"\"\+\+,--\": (\d\.\d+)",
        "se_antisense": r"\+-,-\+\": (\d\.\d+)",
        "failed": r"Fraction of reads failed to determine: (\d\.\d+)",
    }

    # Go through files and parse data using regexes
    for f in module.find_log_files("rseqc/infer_experiment"):
        d = dict()
        for k, r in regexes.items():
            r_search = re.search(r, f["f"], re.MULTILINE)
            if r_search:
                d[k] = float(r_search.group(1))

        if len(d) > 0:
            if f["s_name"] in infer_exp:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            module.add_data_source(f, section="infer_experiment")
            infer_exp[f["s_name"]] = d

    # Filter to strip out ignored sample names
    infer_exp = module.ignore_samples(infer_exp)

    if len(infer_exp) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write to file
    module.write_data_file(infer_exp, "multiqc_rseqc_infer_experiment")

    # Merge PE and SE for plot
    pdata = dict()
    for s_name, vals in infer_exp.items():
        pdata[s_name] = dict()
        for k, v in vals.items():
            v *= 100.0  # Multiply to get percentage
            if k[:2] == "pe" or k[:2] == "se":
                k = k[3:]
            pdata[s_name][k] = v + pdata[s_name].get(k, 0)

    # Plot bar graph of groups
    keys = {
        "sense": {"name": "Sense"},
        "antisense": {"name": "Antisense"},
        "failed": {"name": "Undetermined"},
    }
    # Config for the plot
    pconfig = {
        "id": "rseqc_infer_experiment_plot",
        "title": "RSeQC: Infer experiment",
        "ylab": "% Tags",
        "ymin": 0,
        "ymax": 100,
        "ysuffix": "%",
        "cpswitch": False,
    }

    module.add_section(
        name="Infer experiment",
        anchor="rseqc-infer_experiment",
        description='<a href="http://rseqc.sourceforge.net/#infer-experiment-py" target="_blank">Infer experiment</a>'
        " counts the percentage of reads and read pairs that match the strandedness of overlapping transcripts."
        " It can be used to infer whether RNA-seq library preps are stranded (sense or antisense).",
        plot=bargraph.plot(pdata, keys, pconfig),
    )

    # Return number of samples found
    return len(infer_exp)
