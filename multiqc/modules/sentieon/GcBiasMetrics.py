""" MultiQC submodule to parse output from Sentieon GcBiasMetrics
 (based on the Picard module of the same name """

import logging
import os

from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find Sentieon GcBiasMetrics reports and parse their data"""

    # Set up vars
    self.sentieon_GCbias_data = dict()
    self.sentieon_GCbiasSummary_data = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files("sentieon/gcbias", filehandles=True):
        s_name = None
        gc_col = None
        cov_col = None
        for l in f["f"]:
            # New log starting
            if s_name is None and "--algo GCBias" in l:
                # Pull sample name from filename
                s_name = os.path.basename(f["s_name"])
                s_name = self.clean_s_name(s_name, f)

            if s_name is not None:
                if gc_col is not None and cov_col is not None:
                    try:
                        # Note that GC isn't always the first column.
                        s = l.strip("\n").split("\t")
                        self.sentieon_GCbias_data[s_name][int(s[gc_col])] = float(s[cov_col])
                    except IndexError:
                        s_name = None
                        gc_col = None
                        cov_col = None

                if "#SentieonCommandLine" in l and "--algo GCBias" in l and "Summary" not in s_name:
                    if s_name in self.sentieon_GCbias_data:
                        log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f["fn"], s_name))
                    self.add_data_source(f, s_name, section="GcBiasDetailMetrics")
                    self.sentieon_GCbias_data[s_name] = dict()
                    # Get header - find columns with the data we want
                    line = f["f"].readline()
                    s = line.strip("\n").split("\t")
                    try:
                        gc_col = s.index("GC")
                    except ValueError:
                        pass
                    try:
                        cov_col = s.index("NORMALIZED_COVERAGE")
                    except ValueError:
                        pass

                if "#SentieonCommandLine" in l and "--algo GCBias" in l and "Summary" in s_name:
                    if s_name in self.sentieon_GCbiasSummary_data:
                        log.debug(
                            "Duplicate sample name found in {}!\
                             Overwriting: {}".format(
                                f["fn"], s_name
                            )
                        )
                    self.add_data_source(f, s_name, section="GcBiasSummaryMetrics")
                    self.sentieon_GCbiasSummary_data[s_name] = dict()

                    keys = f["f"].readline().rstrip("\n").split("\t")
                    vals = f["f"].readline().rstrip("\n").split("\t")
                    for i, k in enumerate(keys):
                        try:
                            self.sentieon_GCbiasSummary_data[s_name][k] = float(vals[i])
                        except ValueError:
                            self.sentieon_GCbiasSummary_data[s_name][k] = vals[i]

        for s_name in list(self.sentieon_GCbias_data.keys()):
            if len(self.sentieon_GCbias_data[s_name]) == 0:
                self.sentieon_GCbias_data.pop(s_name, None)
                log.debug("Removing {} as no data parsed".format(s_name))

        for s_name in list(self.sentieon_GCbiasSummary_data.keys()):
            if len(self.sentieon_GCbiasSummary_data[s_name]) == 0:
                self.sentieon_GCbiasSummary_data.pop(s_name, None)
                log.debug("Removing {} as no data parsed".format(s_name))

    # Filter to strip out ignored sample names
    self.sentieon_GCbias_data = self.ignore_samples(self.sentieon_GCbias_data)

    if len(self.sentieon_GCbias_data) > 0:
        # Plot the graph

        pconfig = {
            "id": "sentieon_gcbias_plot",
            "title": "Sentieon: GC Coverage Bias",
            "ylab": "Normalized Coverage",
            "xlab": "% GC",
            "xmin": 0,
            "xmax": 100,
            "xDecimals": False,
            "ymin": 0,
            "yCeiling": 10,
            "tt_label": "<b>{point.x} %GC</b>: {point.y:.2f}",
            "yPlotLines": [
                {"value": 1, "color": "#999999", "width": 2, "dashStyle": "LongDash"},
            ],
        }
        self.add_section(
            name="GC Coverage Bias",
            anchor="sentieon-gcbias",
            description="This plot shows bias in coverage across regions of "
            "the genome with varying GC content. A perfect "
            "library would be a flat line at <code>y = 1</code>.",
            plot=linegraph.plot(self.sentieon_GCbias_data, pconfig),
        )

    if len(self.sentieon_GCbiasSummary_data) > 0:
        # Write parsed summary data to a file
        self.write_data_file(self.sentieon_GCbiasSummary_data, "multiqc_sentieon_gcbias")

    # Return the number of detected samples to the parent module
    return len(self.sentieon_GCbias_data)
