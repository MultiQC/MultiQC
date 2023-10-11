""" MultiQC submodule to parse output from Picard GcBiasMetrics """

import logging
import os
import re

from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find Picard GcBiasMetrics reports and parse their data"""

    # Set up vars
    self.picard_gc_bias_data = dict()
    self.picard_gc_bias_summary_data = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files(f"{self.anchor}/gcbias", filehandles=True):
        s_name = None
        gc_col = None
        cov_col = None
        for l in f["f"]:
            if s_name is None and "--algo GCBias" in l:  # Sentieon
                s_name = os.path.basename(f["s_name"])
                s_name = self.clean_s_name(s_name, f)

            if "GcBiasMetrics" in l and "INPUT" in l:
                s_name = None
                # Pull sample name from input
                fn_search = re.search(r"INPUT(?:=|\s+)(\[?[^\s]+\]?)", l, flags=re.IGNORECASE)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1).strip("[]"))
                    s_name = self.clean_s_name(s_name, f)

            if s_name is not None:
                if gc_col is not None and cov_col is not None:
                    try:
                        # Note that GC isn't always the first column.
                        s = l.strip("\n").split("\t")
                        self.picard_gc_bias_data[s_name][int(s[gc_col])] = float(s[cov_col])
                    except IndexError:
                        s_name = None
                        gc_col = None
                        cov_col = None

                if ("GcBiasDetailMetrics" in l and "## METRICS CLASS" in l) or (
                    "#SentieonCommandLine" in l and "--algo GCBias" in l and "Summary" not in s_name
                ):
                    if s_name in self.picard_gc_bias_data:
                        log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f["fn"], s_name))
                    self.add_data_source(f, s_name, section="GcBiasDetailMetrics")
                    self.picard_gc_bias_data[s_name] = dict()
                    # Get header - find columns with the data we want
                    l = f["f"].readline()
                    s = l.strip("\n").split("\t")
                    try:
                        gc_col = s.index("GC")
                    except ValueError:
                        pass
                    try:
                        cov_col = s.index("NORMALIZED_COVERAGE")
                    except ValueError:
                        pass

                if any(kw in l for kw in ("ACCUMULATION_LEVEL", "GC_DROPOUT")):
                    # Summary metrics
                    # if "GcBiasSummaryMetrics" in l and "## METRICS CLASS" in l:
                    if s_name in self.picard_gc_bias_summary_data:
                        log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
                    self.add_data_source(f, s_name, section="GcBiasSummaryMetrics")
                    self.picard_gc_bias_summary_data[s_name] = dict()

                    keys = f["f"].readline().rstrip("\n").split("\t")
                    vals = f["f"].readline().rstrip("\n").split("\t")
                    for i, k in enumerate(keys):
                        try:
                            self.picard_gc_bias_summary_data[s_name][k] = float(vals[i])
                        except ValueError:
                            self.picard_gc_bias_summary_data[s_name][k] = vals[i]

            for s_name in list(self.picard_gc_bias_data.keys()):
                if len(self.picard_gc_bias_data[s_name]) == 0:
                    self.picard_gc_bias_data.pop(s_name, None)
                    log.debug("Removing {} as no data parsed".format(s_name))

            for s_name in list(self.picard_gc_bias_summary_data.keys()):
                if len(self.picard_gc_bias_summary_data[s_name]) == 0:
                    self.picard_gc_bias_summary_data.pop(s_name, None)
                    log.debug("Removing {} as no data parsed".format(s_name))

            for s_name in set(self.picard_gc_bias_data.keys()) | set(self.picard_gc_bias_summary_data.keys()):
                # Superfluous function call to confirm that it is used in this module
                # Replace None with actual version if it is available
                self.add_software_version(None, s_name)

    # Filter to strip out ignored sample names
    self.picard_gc_bias_data = self.ignore_samples(self.picard_gc_bias_data)

    if len(self.picard_gc_bias_data) > 0:
        # Plot the graph

        pconfig = {
            "id": f"{self.anchor}_gcbias_plot",
            "title": f"{self.name}: GC Coverage Bias",
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
            anchor=f"{self.anchor}-gcbias",
            description="This plot shows bias in coverage across regions of the genome with varying GC content."
            " A perfect library would be a flat line at <code>y = 1</code>.",
            plot=linegraph.plot(self.picard_gc_bias_data, pconfig),
        )

    if len(self.picard_gc_bias_summary_data) > 0:
        # Write parsed summary data to a file
        self.write_data_file(self.picard_gc_bias_summary_data, f"multiqc_{self.anchor}_gcbias")

    # Return the number of detected samples to the parent module
    return len(self.picard_gc_bias_data)
