""" MultiQC submodule to parse output from Picard InsertSizeMetrics """

import logging
from collections import OrderedDict

from multiqc import config
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find Picard InsertSizeMetrics reports and parse their data"""

    # Set up vars
    self.picard_insertSize_data = dict()
    self.picard_insertSize_histogram = dict()
    self.picard_insertSize_samplestats = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files(f"{self.anchor}/insertsize", filehandles=True):
        s_name = None
        in_hist = False
        data = dict()
        histogram = dict()
        samplestats = dict()
        for line in f["f"]:
            maybe_s_name = self.extract_sample_name(line, f)
            if maybe_s_name:
                # Starts information for a new sample
                s_name = maybe_s_name

            # Catch the histogram values
            if s_name is not None and in_hist is True:
                try:
                    sections = line.split("\t")
                    ins = int(sections[0])
                    tot_count = sum([int(x) for x in sections[1:]])
                    histogram[s_name][ins] = tot_count
                    samplestats[s_name]["total_count"] += tot_count
                except ValueError:
                    # Reset in case we have more in this log file
                    s_name = None
                    in_hist = False

            if s_name is not None and self.is_line_right_before_table(line):
                if s_name in data or s_name in self.picard_insertSize_data:
                    log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f["fn"], s_name))
                keys = f["f"].readline().strip("\n").split("\t")
                vals = f["f"].readline().strip("\n").split("\t")
                samplestats[s_name] = {"total_count": 0, "meansum": 0, "total_pairs": 0}
                orientation_idx = keys.index("PAIR_ORIENTATION")
                while len(vals) == len(keys):
                    pair_orientation = vals[orientation_idx]
                    rowkey = "{}_{}".format(s_name, pair_orientation)
                    data[rowkey] = OrderedDict()
                    data[rowkey]["SAMPLE_NAME"] = s_name
                    for i, k in enumerate(keys):
                        try:
                            data[rowkey][k] = float(vals[i])
                        except ValueError:
                            try:
                                data[rowkey][k] = float(vals[i].replace(",", "."))
                                log.debug(
                                    "Switching commas for points in '{}': {} - {}".format(
                                        f["fn"], vals[i], vals[i].replace(",", ".")
                                    )
                                )
                            except ValueError:
                                data[rowkey][k] = vals[i]
                        except IndexError:
                            pass  # missing data
                    # Add to mean sums
                    rp = data[rowkey]["READ_PAIRS"]
                    mis = data[rowkey]["MEAN_INSERT_SIZE"]
                    samplestats[s_name]["meansum"] += rp * mis
                    samplestats[s_name]["total_pairs"] += rp

                    vals = f["f"].readline().strip("\n").split("\t")

                # Skip lines on to histogram
                f["f"].readline().strip("\n")
                f["f"].readline().strip("\n")

                histogram[s_name] = dict()
                in_hist = True

        for key in list(data.keys()):
            if len(data[key]) == 0:
                data.pop(key, None)
        for s_name in list(histogram.keys()):
            if len(histogram[s_name]) == 0:
                histogram.pop(s_name, None)
                log.debug("Ignoring '{}' histogram as no data parsed".format(s_name))

        # When there is only one sample, using the file name to extract the sample name.
        if max(len(data), len(histogram), len(samplestats)) == 1:
            data = {f["s_name"]: list(data.values())[0]} if data else dict()
            histogram = {f["s_name"]: list(histogram.values())[0]} if histogram else dict()
            samplestats = {f["s_name"]: list(samplestats.values())[0]} if samplestats else dict()

        self.picard_insertSize_data.update(data)
        self.picard_insertSize_histogram.update(histogram)
        self.picard_insertSize_samplestats.update(samplestats)

        for s_name in data:
            self.add_data_source(f, s_name, section="InsertSizeMetrics")

    # Calculate summed mean values for all read orientations
    for s_name, v in self.picard_insertSize_samplestats.items():
        try:
            self.picard_insertSize_samplestats[s_name]["summed_mean"] = v["meansum"] / v["total_pairs"]
        except ZeroDivisionError:
            # The mean of zero elements is zero
            self.picard_insertSize_samplestats[s_name]["summed_mean"] = 0

    # Calculate summed median values for all read orientations
    for s_name in self.picard_insertSize_histogram:
        j = 0
        for idx, c in self.picard_insertSize_histogram[s_name].items():
            j += c
            if j > (self.picard_insertSize_samplestats[s_name]["total_count"] / 2):
                self.picard_insertSize_samplestats[s_name]["summed_median"] = idx
                break

    # Filter to strip out ignored sample names
    self.picard_insertSize_data = self.ignore_samples(self.picard_insertSize_data)

    if len(self.picard_insertSize_data) > 0:
        # Write parsed data to a file
        self.write_data_file(self.picard_insertSize_data, f"multiqc_{self.anchor}_insertSize")

        # Do we have median insert sizes?
        missing_medians = False
        for v in self.picard_insertSize_samplestats.values():
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
        for s_name in self.picard_insertSize_samplestats:
            if s_name not in self.general_stats_data:
                self.general_stats_data[s_name] = dict()
            self.general_stats_data[s_name].update(self.picard_insertSize_samplestats[s_name])

        # Section with histogram plot
        if len(self.picard_insertSize_histogram) > 0:
            # Make a normalised percentage version of the data
            data_percent = {}
            for s_name, data in self.picard_insertSize_histogram.items():
                data_percent[s_name] = OrderedDict()
                total = float(sum(data.values()))
                for k, v in data.items():
                    data_percent[s_name][k] = (v / total) * 100

            # Allow customisation of how smooth the plot is
            try:
                insertsize_smooth_points = int(config.picard_config["insertsize_smooth_points"])
                log.debug("Custom Picard insert size smoothing: {}".format(insertsize_smooth_points))
            except (AttributeError, KeyError, ValueError):
                insertsize_smooth_points = 500

            # Plot the data and add section
            pconfig = {
                "smooth_points": insertsize_smooth_points,
                "smooth_points_sumcounts": [True, False],
                "id": f"{self.anchor}_insert_size",
                "title": f"{self.name}: Insert Size",
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
            try:
                pconfig["xmax"] = config.picard_config["insertsize_xmax"]
            except (AttributeError, KeyError):
                pass

            self.add_section(
                name="Insert Size",
                anchor=f"{self.anchor}-insertsize",
                description="Plot shows the number of reads at a given insert size. Reads with different orientations are summed.",
                plot=linegraph.plot([self.picard_insertSize_histogram, data_percent], pconfig),
            )

    # Return the number of detected samples to the parent module
    return len(self.picard_insertSize_data)
