""" MultiQC submodule to parse output from Picard RrbsSummaryMetrics """

import logging
from collections import OrderedDict

from multiqc.modules.picard import util
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(module):
    """Find Picard RrbsSummaryMetrics reports and parse their data"""

    data_by_sample = dict()

    # Go through logs and find Metrics
    for f in module.find_log_files("picard/rrbs_metrics", filehandles=True):
        s_name = None
        for line in f["f"]:
            maybe_s_name = util.extract_sample_name(
                module,
                line,
                f,
                picard_tool="CollectRrbsMetrics",
            )
            if maybe_s_name:
                s_name = maybe_s_name
                continue

            if s_name is None:
                continue

            if util.is_line_right_before_table(line, picard_class="RrbsSummaryMetrics"):
                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
                module.add_data_source(f, s_name, section="RnaSeqMetrics")
                data_by_sample[s_name] = dict()

                keys = f["f"].readline().strip("\n").split("\t")
                vals = f["f"].readline().strip("\n").split("\t")
                if len(vals) != len(keys):
                    continue

                for k, v in zip(keys, vals):
                    if not v:
                        v = "NA"
                    else:
                        try:
                            v = float(v)
                        except ValueError:
                            pass
                    data_by_sample[s_name][k] = v

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)
    if len(data_by_sample) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed data to a file
    module.write_data_file(data_by_sample, "multiqc_picard_RrbsSummaryMetrics")

    # Add to general stats table
    headers = {
        "PCT_CPG_BASES_CONVERTED": {
            "title": "CpG Methylated",
            "description": "Percentage of times a CpG cytosine was converted",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "format": "{:,.0f}",
            "scale": "RdYlGn-rev",
            "modify": lambda x: 100 - util.multiply_hundred(x),
        },
        "PCT_NON_CPG_BASES_CONVERTED": {
            "title": "Non-CpG Methylated",
            "description": "Percentage of times a non-CpG cytosine was converted",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "format": "{:,.0f}",
            "scale": "RdYlGn",
            "modify": lambda x: 100 - util.multiply_hundred(x),
        },
        "MEDIAN_CPG_COVERAGE": {
            "title": "Median CpG Cov",
            "description": "Median coverage of CpG sites",
            "min": 0,
        },
    }
    module.general_stats_addcols(data_by_sample, headers, namespace="RrbsSummaryMetrics")

    # Make the bar plot of converted bases
    pdata_cpg = dict()
    pdata_noncpg = dict()
    for s_name in data_by_sample.keys():
        pdata_cpg[s_name] = dict()
        pdata_cpg[s_name]["converted"] = data_by_sample[s_name]["CPG_BASES_CONVERTED"]
        pdata_cpg[s_name]["not_converted"] = (
            data_by_sample[s_name]["CPG_BASES_SEEN"] - data_by_sample[s_name]["CPG_BASES_CONVERTED"]
        )
        pdata_noncpg[s_name] = dict()
        pdata_noncpg[s_name]["converted"] = data_by_sample[s_name]["NON_CPG_BASES"]
        pdata_noncpg[s_name]["not_converted"] = (
            data_by_sample[s_name]["NON_CPG_BASES"] - data_by_sample[s_name]["NON_CPG_CONVERTED_BASES"]
        )

    keys = dict()
    keys["not_converted"] = {"name": "Unconverted Bases (Methylated)"}
    keys["converted"] = {"name": "Converted Bases (Unmethylated)"}

    # Config for the plot
    pconfig = {
        "id": "picard_rrbs_converted_bases_plot",
        "title": "Picard: Converted Bases",
        "ylab": "# CpG Bases",
        "cpswitch_counts_label": "Number of Bases",
        "data_labels": [{"name": "CpG", "ylab": "# CpG Bases"}, {"name": "Non-CpG", "ylab": "# Non-CpG Bases"}],
    }

    module.add_section(
        name="RRBS Converted Bases",
        anchor="picard-rrbssummary-convertedbases",
        plot=bargraph.plot([pdata_cpg, pdata_noncpg], [keys, keys], pconfig),
    )

    # Make the bar plot of processed reads
    pdata = dict()
    for s_name in data_by_sample:
        pdata[s_name] = dict()
        pdata[s_name]["with_no_cpg"] = data_by_sample[s_name]["READS_WITH_NO_CPG"]
        pdata[s_name]["ignored_short"] = data_by_sample[s_name]["READS_IGNORED_SHORT"]
        pdata[s_name]["ignored_mismatches"] = data_by_sample[s_name]["READS_IGNORED_MISMATCHES"]
        pdata[s_name]["not_ignored"] = (
            data_by_sample[s_name]["READS_ALIGNED"]
            - pdata[s_name]["with_no_cpg"]
            - pdata[s_name]["ignored_short"]
            - pdata[s_name]["ignored_mismatches"]
        )

    keys = OrderedDict()
    keys["not_ignored"] = {"name": "Utilised reads"}
    keys["with_no_cpg"] = {"name": "Ignored (no CpG sites)"}
    keys["ignored_short"] = {"name": "Ignored (too short)"}
    keys["ignored_mismatches"] = {"name": "Ignored (exceeded mismatch threshold)"}

    # Config for the plot
    pconfig = {
        "id": "picard_rrbs_ignored_reads_plot",
        "title": "Picard: RRBS Read Counts",
        "ylab": "# Reads",
        "cpswitch_counts_label": "Number of Reads",
    }

    module.add_section(
        name="RRBS Read Counts", anchor="picard-rrbssummary-readcounts", plot=bargraph.plot(pdata, keys, pconfig)
    )

    # Return the number of detected samples to the parent module
    return len(data_by_sample)
