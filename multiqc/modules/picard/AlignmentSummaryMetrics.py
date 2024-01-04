""" MultiQC submodule to parse output from Picard AlignmentSummaryMetrics """

import logging
from collections import OrderedDict

from multiqc.modules.picard import util
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(module):
    """Find Picard AlignmentSummaryMetrics reports and parse their data"""

    data_by_sample = dict()

    # Go through logs and find Metrics
    for f in module.find_log_files("picard/alignment_metrics", filehandles=True):
        # Sample name from input file name by default.
        s_name = f["s_name"]
        keys = None

        for line in f["f"]:
            maybe_s_name = util.extract_sample_name(
                module,
                line,
                f,
                picard_tool="CollectAlignmentSummaryMetrics",
                sentieon_algo="AlignmentStat",
            )
            if maybe_s_name:
                s_name = maybe_s_name
                keys = None

            if s_name is None:
                continue

            if util.is_line_right_before_table(
                line, picard_class="AlignmentSummaryMetrics", sentieon_algo="AlignmentStat"
            ):
                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: " f"{s_name}")
                data_by_sample[s_name] = dict()
                module.add_data_source(f, s_name, section="AlignmentSummaryMetrics")
                keys = f["f"].readline().strip("\n").split("\t")

            elif keys:
                vals = line.strip("\n").split("\t")
                if len(vals) != len(keys):
                    keys = None
                    continue

                # Ignore the FIRST_OF_PAIR / SECOND_OF_PAIR data to simplify things
                if vals[0] == "PAIR" or vals[0] == "UNPAIRED":
                    for k, v in zip(keys, vals):
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
    module.write_data_file(data_by_sample, f"multiqc_{module.anchor}_AlignmentSummaryMetrics")

    # Add to general stats table
    module.general_stats_headers["PCT_PF_READS_ALIGNED"] = {
        "title": "% Aligned",
        "description": "Percent of aligned reads",
        "max": 100,
        "min": 0,
        "suffix": "%",
        "format": "{:,.0f}",
        "scale": "RdYlGn",
        "modify": lambda x: util.multiply_hundred(x),
    }
    for s_name in data_by_sample:
        if s_name not in module.general_stats_data:
            module.general_stats_data[s_name] = dict()
        module.general_stats_data[s_name].update(data_by_sample[s_name])

    # Make the bar plot of alignment read count + # aligned bases
    pdata = dict()
    for s_name in data_by_sample.keys():
        pdata[s_name] = dict()
        # Picard reports both reads for PE data. Divide it by two as most people will
        # expect # clusters
        if data_by_sample[s_name]["CATEGORY"] == "PAIR":
            pdata[s_name]["total_reads"] = data_by_sample[s_name]["TOTAL_READS"] / 2
            pdata[s_name]["aligned_reads"] = data_by_sample[s_name]["PF_READS_ALIGNED"] / 2
        else:
            pdata[s_name]["total_reads"] = data_by_sample[s_name]["TOTAL_READS"]
            pdata[s_name]["aligned_reads"] = data_by_sample[s_name]["PF_READS_ALIGNED"]
        pdata[s_name]["unaligned_reads"] = pdata[s_name]["total_reads"] - pdata[s_name]["aligned_reads"]

    keys = [OrderedDict(), OrderedDict()]
    keys[0]["aligned_reads"] = {"name": "Aligned Reads"}
    keys[0]["unaligned_reads"] = {"name": "Unaligned Reads"}
    keys[1]["PF_ALIGNED_BASES"] = {"name": "Aligned Bases"}

    # Config for the plot
    pconfig = {
        "id": f"{module.anchor}_alignment_summary",
        "title": f"{module.name}: Alignment Summary",
        "ylab": "# Reads",
        "data_labels": [
            {
                "name": "Aligned Reads",
                "ylab": "# Reads",
                "cpswitch_counts_label": "Number of Reads",
            },
            {
                "name": "Aligned Bases",
                "ylab": "# Bases",
                "cpswitch_counts_label": "Number of Bases",
            },
        ],
    }

    # The different data sets we want to plot
    module.add_section(
        name="Alignment Summary",
        anchor=f"{module.anchor}-alignmentsummary",
        description=f"Please note that {module.name}'s read counts are divided by two "
        f"for paired-end data. Total bases (including unaligned) is not "
        f"provided.",
        plot=bargraph.plot([pdata, data_by_sample], keys, pconfig),
    )

    # Make a bar plot of mean read length
    keys = {"MEAN_READ_LENGTH": {"name": "Mean Read Length"}}
    pconfig = {
        "id": f"{module.anchor}_alignment_readlength_plot",
        "title": f"{module.name}: Mean Read Length",
        "ylab": "Base pairs",
        "cpswitch": False,
    }

    # The different data sets we want to plot
    module.add_section(
        name="Mean read length",
        anchor=f"{module.anchor}_alignment_readlength",
        description="The mean read length of the set of reads examined.",
        plot=bargraph.plot(data_by_sample, keys, pconfig),
    )

    # Return the number of detected samples to the parent module
    return len(data_by_sample)
