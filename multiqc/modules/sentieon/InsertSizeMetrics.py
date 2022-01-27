#!/usr/bin/env python

""" MultiQC submodule to parse output from Sentieon InsertSizeMetrics (based
 on the Picard module of the same name """

from collections import OrderedDict
import logging
import os

from multiqc import config
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find Sentieon InsertSizeMetrics reports and parse their data"""

    # Set up vars
    self.sentieon_insertSize_data = dict()
    self.sentieon_insertSize_histogram = dict()
    self.sentieon_insertSize_samplestats = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files("sentieon/insertsize", filehandles=True):
        s_name = None
        in_hist = False
        for line in f["f"]:

            # Catch the histogram values
            if s_name is not None and in_hist is True:
                try:
                    sections = line.split("\t")
                    ins = int(sections[0])
                    tot_count = sum([int(x) for x in sections[1:]])
                    self.sentieon_insertSize_histogram[s_name][ins] = tot_count
                    self.sentieon_insertSize_samplestats[s_name]["total_count"] += tot_count
                except ValueError:
                    # Reset in case we have more in this log file
                    s_name = None
                    in_hist = False

            # New log starting
            if s_name is None and "InsertSizeMetricAlgo" in line:
                # Pull sample name from filename
                s_name = os.path.basename(f["s_name"])
                s_name = self.clean_s_name(s_name, f)

            if s_name is not None:
                if "InsertSizeMetricAlgo" in line and "#SentieonCommandLine" in line:
                    if s_name in self.sentieon_insertSize_data:
                        log.debug(
                            "Duplicate sample name found in {}!\
                             Overwriting: {}".format(
                                f["fn"], s_name
                            )
                        )
                    self.add_data_source(f, s_name, section="InsertSizeMetrics")
                    keys = f["f"].readline().strip("\n").split("\t")
                    vals = f["f"].readline().strip("\n").split("\t")
                    self.sentieon_insertSize_samplestats[s_name] = {"total_count": 0, "meansum": 0, "total_pairs": 0}
                    orientation_idx = keys.index("PAIR_ORIENTATION")
                    while len(vals) == len(keys):
                        pair_orientation = vals[orientation_idx]
                        rowkey = "{}_{}".format(s_name, pair_orientation)
                        self.sentieon_insertSize_data[rowkey] = OrderedDict()
                        self.sentieon_insertSize_data[rowkey]["SAMPLE_NAME"] = s_name
                        for i, k in enumerate(keys):
                            try:
                                self.sentieon_insertSize_data[rowkey][k] = float(vals[i])
                            except ValueError:
                                try:
                                    self.sentieon_insertSize_data[rowkey][k] = float(vals[i].replace(",", "."))
                                    log.debug(
                                        "Switching commas for points in '{}':\
                                             {} - {}".format(
                                            f["fn"], vals[i], vals[i].replace(",", ".")
                                        )
                                    )
                                except ValueError:
                                    self.sentieon_insertSize_data[rowkey][k] = vals[i]
                            except IndexError:
                                pass  # missing data
                        # Add to mean sums
                        rp = self.sentieon_insertSize_data[rowkey]["READ_PAIRS"]
                        mis = self.sentieon_insertSize_data[rowkey]["MEAN_INSERT_SIZE"]
                        self.sentieon_insertSize_samplestats[s_name]["meansum"] += rp * mis
                        self.sentieon_insertSize_samplestats[s_name]["total_pairs"] += rp

                        vals = f["f"].readline().strip("\n").split("\t")

                    # Skip lines on to histogram (variable used in later loops)
                    line = f["f"].readline().strip("\n")
                    line = f["f"].readline().strip("\n")

                    self.sentieon_insertSize_histogram[s_name] = OrderedDict()
                    in_hist = True

        for key in list(self.sentieon_insertSize_data.keys()):
            if len(self.sentieon_insertSize_data[key]) == 0:
                self.sentieon_insertSize_data.pop(key, None)
        for s_name in list(self.sentieon_insertSize_histogram.keys()):
            if len(self.sentieon_insertSize_histogram[s_name]) == 0:
                self.sentieon_insertSize_histogram.pop(s_name, None)
                log.debug("Ignoring '{}' histogram as no data parsed".format(s_name))

    # Calculate summed mean values for all read orientations
    for s_name, v in self.sentieon_insertSize_samplestats.items():
        try:
            self.sentieon_insertSize_samplestats[s_name]["summed_mean"] = v["meansum"] / v["total_pairs"]
        except ZeroDivisionError:
            # The mean of zero elements is zero
            self.sentieon_insertSize_samplestats[s_name]["summed_mean"] = 0

    # Calculate summed median values for all read orientations
    for s_name in self.sentieon_insertSize_histogram:
        j = 0
        for idx, c in self.sentieon_insertSize_histogram[s_name].items():
            j += c
            if j > (self.sentieon_insertSize_samplestats[s_name]["total_count"] / 2):
                self.sentieon_insertSize_samplestats[s_name]["summed_median"] = idx
                break

    # Filter to strip out ignored sample names
    self.sentieon_insertSize_data = self.ignore_samples(self.sentieon_insertSize_data)

    if len(self.sentieon_insertSize_data) > 0:

        # Write parsed data to a file
        self.write_data_file(self.sentieon_insertSize_data, "multiqc_sentieon_insertSize")

        # Do we have median insert sizes?
        missing_medians = False
        for v in self.sentieon_insertSize_samplestats.values():
            if "summed_median" not in v:
                missing_medians = True

        # Add to general stats table
        self.general_stats_headers["summed_median"] = {
            "title": "Insert Size",
            "description": "Median Insert Size, all read orientations (bp)",
            "min": 0,
            "suffix": " bp",
            "format": "{:,.0f}",
            "scale": "GnBu",
        }
        self.general_stats_headers["summed_mean"] = {
            "title": "Mean Insert Size",
            "description": "Mean Insert Size, all read orientations (bp)",
            "min": 0,
            "suffix": " bp",
            "format": "{:,.0f}",
            "scale": "GnBu",
            "hidden": False if missing_medians else True,
        }
        for s_name in self.sentieon_insertSize_samplestats:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update(self.sentieon_insertSize_samplestats[s_name])

        # Section with histogram plot
        if len(self.sentieon_insertSize_histogram) > 0:
            # Make a normalised percentage version of the data
            data_percent = {}
            for s_name, data in self.sentieon_insertSize_histogram.items():
                data_percent[s_name] = OrderedDict()
                total = float(sum(data.values()))
                for k, v in data.items():
                    data_percent[s_name][k] = (v / total) * 100

            # Allow customisation of how smooth the the plot is
            try:
                insertsize_smooth_points = int(config.sentieon_config["insertsize_smooth_points"])
                log.debug(
                    "Custom Sentieon insert size smoothing:\
                     {}".format(
                        insertsize_smooth_points
                    )
                )
            except (AttributeError, KeyError, ValueError):
                insertsize_smooth_points = 500

            # Plot the data and add section
            pconfig = {
                "smooth_points": insertsize_smooth_points,
                "smooth_points_sumcounts": [True, False],
                "id": "sentieon_insert_size",
                "title": "Sentieon: Insert Size",
                "ylab": "Count",
                "xlab": "Insert Size (bp)",
                "xDecimals": False,
                "tt_label": "<b>{point.x} bp</b>: {point.y:.0f}",
                "ymin": 0,
                "data_labels": [
                    {"name": "Counts", "ylab": "Coverage"},
                    {"name": "Percentages", "ylab": "Percentage of Counts"},
                ],
            }
            self.add_section(
                name="Insert Size",
                anchor="sentieon-insertsize",
                description="Plot shows the number of reads at a given insert\
                     size. Reads with different orientations are summed.",
                plot=linegraph.plot([self.sentieon_insertSize_histogram, data_percent], pconfig),
            )

    # Return the number of detected samples to the parent module
    return len(self.sentieon_insertSize_data)
