"""MultiQC submodule to parse output from Picard TargetedPcrMetrics"""

import logging
from typing import Dict

from multiqc.modules.picard import util
from multiqc.plots import bargraph
from multiqc import config

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(module):
    """Find Picard TargetedPcrMetrics reports and parse their data"""

    data_by_sample: Dict = dict()
    histogram_by_sample: Dict = dict()

    picard_config = getattr(config, "picard_config", {})
    skip_histo = picard_config.get("targeted_pcr_skip_histogram", False)

    # Go through logs and find Metrics
    for f in module.find_log_files("picard/pcr_metrics", filehandles=True):
        # Sample name from input file name by default.
        s_name = f["s_name"]
        in_hist = False

        for line in f["f"]:
            maybe_s_name = util.extract_sample_name(
                module,
                line,
                f,
                picard_tool="TargetedPcrMetrics",
            )
            if maybe_s_name:
                s_name = maybe_s_name

            if s_name is None:
                continue

            # Catch the histogram values
            if in_hist and not skip_histo:
                try:
                    sections = line.split("\t")
                    cov = int(sections[0])
                    count = int(sections[1])
                    histogram_by_sample[s_name][cov] = count
                except ValueError:
                    # Reset in case we have more in this log file
                    s_name = None
                    in_hist = False

            if util.is_line_right_before_table(line, picard_class="TargetedPcrMetrics"):
                keys = f["f"].readline().strip("\n").split("\t")
                vals = f["f"].readline().strip("\n").split("\t")
                if len(vals) != len(keys):
                    continue

                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")

                module.add_data_source(f, s_name, section="TargetedPcrMetrics")
                data_by_sample[s_name] = dict()
                histogram_by_sample[s_name] = dict()

                for k, v in zip(keys, vals):
                    try:
                        # Multiply percentages by 100
                        if k.startswith("PCT_"):
                            v = float(v) * 100.0
                    except ValueError:
                        pass
                    data_by_sample[s_name][k] = v

            elif line.startswith("## HISTOGRAM"):
                keys = f["f"].readline().strip("\n").split("\t")
                assert len(keys) >= 2, (keys, f)
                in_hist = True
                histogram_by_sample[s_name] = dict()

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)
    histogram_by_sample = module.ignore_samples(histogram_by_sample)
    if len(data_by_sample) == 0:
        return set()

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed data to a file
    module.write_data_file(data_by_sample, "multiqc_picard_pcrmetrics")

    # Add to general stats table
    headers = {
        "PCT_AMPLIFIED_BASES": {
            "title": "Amplified Bases",
            "description": "The fraction of aligned bases that mapped to or near an amplicon.",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "BrBG",
        },
        "MEDIAN_TARGET_COVERAGE": {
            "title": "Median Target Coverage",
            "description": "The median coverage of reads that mapped to target regions of an experiment.",
            "min": 0,
            "suffix": "X",
            "scale": "GnBu",
        },
    }
    module.general_stats_addcols(data_by_sample, headers, namespace="TargetedPcrMetrics")

    # Bar plot of ignored bases
    keys = {
        "ON_AMPLICON_BASES": {"name": "On-amplicon bases"},
        "NEAR_AMPLICON_BASES": {"name": "Near-amplicon bases"},
        "OFF_AMPLICON_BASES": {"name": "Off-amplicon bases", "color": "#f28f43"},
    }

    # Config for the plot
    pconfig = {
        "id": "picard_pcr_metrics_bases",
        "title": "Picard: PCR Amplicon Bases",
        "ylab": "# Bases",
        "cpswitch_counts_label": "# Bases",
        "hide_empty": False,
    }

    module.add_section(
        name="PCR Amplicon Bases",
        anchor="picard-pcrmetrics-bases",
        description="Metrics about reads obtained from targeted PCR experiments.",
        helptext="""
        This plot shows the number of bases aligned on or near to amplified regions of the genome.

        * `ON_AMPLICON_BASES`: The number of `PF_BASES_ALIGNED` that mapped to an amplified region of the genome.
        * `NEAR_AMPLICON_BASES`: The number of `PF_BASES_ALIGNED` that mapped to within a fixed interval of an amplified region, but not on a baited region.
        * `OFF_AMPLICON_BASES`: The number of `PF_BASES_ALIGNED` that mapped neither on or near an amplicon.

        For more information see the [Picard documentation](https://broadinstitute.github.io/picard/picard-metric-definitions.html#TargetedPcrMetrics).""",
        plot=bargraph.plot(data_by_sample, keys, pconfig),
    )

    # Return the number of detected samples to the parent module
    return data_by_sample.keys()
