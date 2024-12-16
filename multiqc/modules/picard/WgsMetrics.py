"""MultiQC submodule to parse output from Picard WgsMetrics"""

import logging
from typing import Dict

from multiqc import config
from multiqc.modules.picard import util
from multiqc.plots import bargraph, linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(module):
    """Find Picard WgsMetrics reports and parse their data"""

    # Set up vars
    data_by_sample: Dict = dict()
    histogram_by_sample: Dict = dict()

    picard_config = getattr(config, "picard_config", {})
    skip_histo = picard_config.get("wgsmetrics_skip_histogram", False)

    # Go through logs and find Metrics
    for f in module.find_log_files("picard/wgs_metrics", filehandles=True):
        # Sample name from input file name by default
        s_name = f["s_name"]
        in_hist = False

        for line in f["f"]:
            maybe_s_name = util.extract_sample_name(
                module,
                line,
                f,
                picard_tool="CollectWgsMetrics",
            )
            if maybe_s_name:
                s_name = maybe_s_name

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

            if s_name is None:
                continue

            if util.is_line_right_before_table(line, picard_class="WgsMetrics"):
                keys = f["f"].readline().strip("\n").split("\t")
                vals = f["f"].readline().strip("\n").split("\t")
                if len(vals) != len(keys):
                    continue

                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")

                module.add_data_source(f, s_name, section="WgsMetrics")
                data_by_sample[s_name] = dict()

                for k, v in zip(keys, vals):
                    try:
                        v = float(v)
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
    if len(data_by_sample) == 0:
        return set()

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed data to a file
    module.write_data_file(data_by_sample, "multiqc_picard_wgsmetrics")

    # Add to general stats table
    headers = dict()
    headers["MEDIAN_COVERAGE"] = {
        "title": "Median Coverage",
        "description": "The median coverage in bases of the genome territory, after all filters are applied.",
        "min": 0,
        "suffix": "X",
        "scale": "GnBu",
    }
    headers["MEAN_COVERAGE"] = {
        "title": "Mean Coverage",
        "description": "The mean coverage in bases of the genome territory, after all filters are applied.",
        "min": 0,
        "suffix": "X",
        "scale": "GnBu",
        "hidden": True,
    }
    headers["SD_COVERAGE"] = {
        "title": "SD Coverage",
        "description": "The standard deviation coverage in bases of the genome territory, after all filters are applied.",
        "min": 0,
        "suffix": "X",
        "scale": "GnBu",
        "hidden": True,
    }
    # user configurable coverage level
    picard_config = getattr(config, "picard_config", {})
    covs = picard_config.get("general_stats_target_coverage", [])
    if isinstance(covs, list) and len(covs) > 0:
        covs = [str(i) for i in covs]
        log.debug(f"Custom Picard coverage thresholds: {', '.join([i for i in covs])}")
    else:
        covs = ["30"]
    for c in covs:
        headers[f"PCT_{c}X"] = {
            "rid": f"picard_target_bases_{c}X",
            "title": f"Bases ≥ {c}X",
            "description": f"Percent of target bases with coverage ≥ {c}X",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "format": "{:,.0f}",
            "scale": "RdYlGn",
            "modify": lambda x: util.multiply_hundred(x),
        }
    module.general_stats_addcols(data_by_sample, headers, namespace="WgsMetrics")

    # Section with histogram plot
    if histogram_by_sample and not skip_histo:
        # Figure out where to cut histogram tail
        max_cov = picard_config.get("wgsmetrics_histogram_max_cov")
        if max_cov is None:
            max_cov = 10
            for s_name, hist in histogram_by_sample.items():
                total = float(sum(hist.values()))
                running_total = 0
                for k, v in hist.items():
                    running_total += v
                    if running_total > total * 0.99:
                        max_cov = max(k, max_cov)
                        break

        # Cut histogram tail and make a normalised percentage version of the data plus dropoff
        data: Dict = {}
        data_percent: Dict = {}
        maxval = 0
        for s_name, hist in histogram_by_sample.items():
            data[s_name] = dict()
            data_percent[s_name] = dict()
            total = float(sum(hist.values()))
            cumulative = 0
            for k, v in hist.items():
                if k <= max_cov:
                    cumulative += v
                    data[s_name][k] = v
                    maxval = max(maxval, v)
                    pct_bases_cumulative = (cumulative / total) * 100
                    pct_bases_with_greater_cov = 100 - pct_bases_cumulative
                    pct_bases_current = (v / total) * 100
                    pct_bases_with_greater_or_equal_cov = pct_bases_with_greater_cov + pct_bases_current
                    data_percent[s_name][k] = pct_bases_with_greater_or_equal_cov
                else:
                    break

        # Plot the histogram data and add section
        pconfig = {
            "id": "picard_wgs_metrics_histogram",
            "title": "Picard: WGS Coverage",
            "ylab": "Percentage of Bases",
            "xlab": "Fold Coverage",
            "x_decimals": False,
            "tt_label": "<b>{point.x} X</b>: {point.y:.1f}",
            "ymin": 0,
            "ymax": 100,
            "smooth_points": picard_config.get("wgsmetrics_histogram_smooth", 1000),
            "data_labels": [
                {"name": "Percentage Drop-Off", "ylab": "Percentage of Bases", "ymax": 100},
                {"name": "Counts Histogram", "ylab": "Coverage", "ymax": maxval},
            ],
        }
        module.add_section(
            name="WGS Coverage",
            anchor="picard-wgsmetrics-cov",
            description="The number of bases in the genome territory for each fold coverage. "
            + "Note that final 1% of data is hidden to prevent very long tails.",
            plot=linegraph.plot([data_percent, data], pconfig),
        )

    # Bar plot of ignored bases
    pdata: Dict = dict()
    for s_name, data in data_by_sample.items():
        pdata[s_name] = dict()
        pdata[s_name]["PCT_EXC_MAPQ"] = data["PCT_EXC_MAPQ"] * 100.0
        pdata[s_name]["PCT_EXC_DUPE"] = data["PCT_EXC_DUPE"] * 100.0
        pdata[s_name]["PCT_EXC_UNPAIRED"] = data["PCT_EXC_UNPAIRED"] * 100.0
        pdata[s_name]["PCT_EXC_BASEQ"] = data["PCT_EXC_BASEQ"] * 100.0
        pdata[s_name]["PCT_EXC_OVERLAP"] = data["PCT_EXC_OVERLAP"] * 100.0
        pdata[s_name]["PCT_EXC_CAPPED"] = data["PCT_EXC_CAPPED"] * 100.0

    keys = dict()
    keys["PCT_EXC_MAPQ"] = {"name": "Low mapping quality"}
    keys["PCT_EXC_DUPE"] = {"name": "Duplicates reads"}
    keys["PCT_EXC_UNPAIRED"] = {"name": "No mapped mate pair"}
    keys["PCT_EXC_BASEQ"] = {"name": "Low base quality"}
    keys["PCT_EXC_OVERLAP"] = {"name": "Overlapping insert"}
    keys["PCT_EXC_CAPPED"] = {"name": "Over capped coverage"}

    # Config for the plot
    pconfig = {
        "id": "picard_wgs_metrics_bases",
        "title": "Picard: WGS Filtered Bases",
        "cpswitch": False,
        "ylab": "% Bases",
        "ymax": 100,
    }

    module.add_section(
        name="WGS Filtered Bases",
        anchor="picard-wgsmetrics-bases",
        description="For more information about the filtered categories, see the "
        + '<a href="http://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectWgsMetrics.WgsMetrics" target="_blank">Picard documentation</a>.',
        plot=bargraph.plot(pdata, keys, pconfig),
    )

    # Return the number of detected samples to the parent module
    return data_by_sample.keys()
