"""MultiQC submodule to parse output from Picard InsertSizeMetrics"""

import logging
from typing import Dict

from multiqc import config
from multiqc.modules.picard import util
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(module):
    """Find Picard InsertSizeMetrics reports and parse their data"""

    data_by_sample: Dict = dict()
    histogram_by_sample: Dict = dict()
    samplestats_by_sample: Dict = dict()

    # Go through logs and find Metrics
    for f in module.find_log_files("picard/insertsize", filehandles=True):
        # Sample name from input file name by default
        s_name = f["s_name"]
        in_hist = False

        for line in f["f"]:
            maybe_s_name = util.extract_sample_name(
                module,
                line,
                f,
                picard_tool="CollectInsertSizeMetrics",
                sentieon_algo="InsertSizeMetricAlgo",
            )
            if maybe_s_name:
                s_name = maybe_s_name

            if s_name is None:
                continue

            # Catch the histogram values
            if in_hist:
                try:
                    sections = line.split("\t")
                    ins = int(sections[0])
                    tot_count = sum([int(x) for x in sections[1:]])
                    histogram_by_sample[s_name][ins] = tot_count
                    samplestats_by_sample[s_name]["total_count"] += tot_count
                except ValueError:
                    # Reset in case we have more in this log file
                    s_name = None
                    in_hist = False

            if util.is_line_right_before_table(
                line, picard_class="InsertSizeMetrics", sentieon_algo="InsertSizeMetricAlgo"
            ):
                keys = f["f"].readline().strip("\n").split("\t")
                vals = f["f"].readline().strip("\n").split("\t")
                if len(vals) != len(keys):
                    continue

                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")

                module.add_data_source(f, s_name, section="InsertSizeMetrics")
                samplestats_by_sample[s_name] = {"total_count": 0, "meansum": 0, "total_pairs": 0}
                orientation_idx = keys.index("PAIR_ORIENTATION")

                while len(vals) == len(keys):
                    pair_orientation = vals[orientation_idx]
                    rowkey = f"{s_name}_{pair_orientation}"
                    data_by_sample[rowkey] = dict()
                    data_by_sample[rowkey]["SAMPLE_NAME"] = s_name
                    for i, k in enumerate(keys):
                        try:
                            data_by_sample[rowkey][k] = float(vals[i])
                        except ValueError:
                            try:
                                unfixed = vals[i]
                                fixed = unfixed.replace(",", ".")
                                data_by_sample[rowkey][k] = float(fixed)
                                log.debug(f"Switching commas for points in '{f['fn']}': {unfixed} -> {fixed}")
                            except ValueError:
                                data_by_sample[rowkey][k] = vals[i]
                        except IndexError:
                            pass  # missing data
                    # Add to mean sums
                    rp = data_by_sample[rowkey]["READ_PAIRS"]
                    mis = data_by_sample[rowkey]["MEAN_INSERT_SIZE"]
                    samplestats_by_sample[s_name]["meansum"] += rp * mis
                    samplestats_by_sample[s_name]["total_pairs"] += rp

                    vals = f["f"].readline().strip("\n").split("\t")

            if line.startswith("## HISTOGRAM"):
                keys = f["f"].readline().strip("\n").split("\t")
                assert len(keys) >= 2, (keys, f)
                in_hist = True
                histogram_by_sample[s_name] = dict()

    # Calculate summed mean values for all read orientations
    for s_name, v in samplestats_by_sample.items():
        try:
            samplestats_by_sample[s_name]["summed_mean"] = v["meansum"] / v["total_pairs"]
        except ZeroDivisionError:
            # The mean of zero elements is zero
            samplestats_by_sample[s_name]["summed_mean"] = 0

    # Calculate summed median values for all read orientations
    for s_name in histogram_by_sample:
        j = 0
        for idx, c in histogram_by_sample[s_name].items():
            j += c
            if j > (samplestats_by_sample[s_name]["total_count"] / 2):
                samplestats_by_sample[s_name]["summed_median"] = idx
                break

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)
    if len(data_by_sample) == 0:
        return set()

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed data to a file
    module.write_data_file(data_by_sample, f"multiqc_{module.anchor}_insertSize")

    # Do we have median insert sizes?
    missing_medians = False
    for v in samplestats_by_sample.values():
        if "summed_median" not in v:
            missing_medians = True

    # Add to general stats table
    headers = {
        "summed_median": {
            "title": "Insert Size",
            "description": "Median Insert Size, all read orientations (bp)",
            "min": 0,
            "suffix": " bp",
            "format": "{:,.0f}",
            "scale": "GnBu",
        },
        "summed_mean": {
            "title": "Mean Insert Size",
            "description": "Mean Insert Size, all read orientations (bp)",
            "min": 0,
            "suffix": " bp",
            "format": "{:,.0f}",
            "scale": "GnBu",
            "hidden": False if missing_medians else True,
        },
    }
    module.general_stats_addcols(samplestats_by_sample, headers, namespace="InsertSizeMetrics")

    # Section with histogram plot
    if len(histogram_by_sample) > 0:
        # Make a normalised percentage version of the data
        data_percent: Dict = {}
        for s_name, data in histogram_by_sample.items():
            data_percent[s_name] = dict()
            total = float(sum(data.values()))
            for k, v in data.items():
                data_percent[s_name][k] = (v / total) * 100

        # Allow customisation of how smooth the plot is
        insertsize_smooth_points = getattr(config, "picard_config", {}).get("insertsize_smooth_points")
        if insertsize_smooth_points is not None:
            log.debug(f"Custom Picard insert size smoothing: {insertsize_smooth_points}")
            insertsize_smooth_points = int(insertsize_smooth_points)
        else:
            insertsize_smooth_points = 500

        # Plot the data and add section
        pconfig = {
            "smooth_points": insertsize_smooth_points,
            "smooth_points_sumcounts": True,
            "id": f"{module.anchor}_insert_size",
            "title": f"{module.name}: Insert Size",
            "ylab": "Count",
            "xlab": "Insert Size (bp)",
            "x_decimals": False,
            "tt_label": "<b>{point.x} bp</b>: {point.y:.0f}",
            "ymin": 0,
            "data_labels": [
                {"name": "Counts", "ylab": "Coverage"},
                {"name": "Percentages", "ylab": "Percentage of Counts"},
            ],
        }
        insertsize_xmax = getattr(config, "picard_config", {}).get("insertsize_xmax")
        if insertsize_xmax is not None:
            pconfig["xmax"] = int(insertsize_xmax)

        module.add_section(
            name="Insert Size",
            anchor=f"{module.anchor}-insertsize",
            description="Plot shows the number of reads at a given insert size. Reads "
            "with different orientations are summed.",
            plot=linegraph.plot([histogram_by_sample, data_percent], pconfig),
        )

    # Return the number of detected samples to the parent module
    return data_by_sample.keys()
