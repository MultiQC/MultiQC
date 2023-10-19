""" MultiQC submodule to parse output from Picard GcBiasMetrics """

import logging

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
        s_name = f["s_name"]
        gc_col = None
        cov_col = None
        data = dict()
        summary_data = dict()
        for l in f["f"]:
            maybe_s_name = self.extract_sample_name(l, f, picard_tool="CollectGcBiasMetrics", sentieon_algo="GCBias")
            if maybe_s_name:
                # Starts information for a new sample
                s_name = maybe_s_name

            if s_name is not None:
                if gc_col is not None and cov_col is not None:
                    try:
                        # Note that GC isn't always the first column.
                        s = l.strip("\n").split("\t")
                        data[s_name][int(s[gc_col])] = float(s[cov_col])
                    except IndexError:
                        s_name = None
                        gc_col = None
                        cov_col = None

                elif self.is_line_right_before_table(
                    l, picard_class=["GcBiasDetailMetrics", "GcBiasSummaryMetrics"], sentieon_algo="GCBias"
                ):
                    # Get header - find columns with the data we want
                    l = f["f"].readline()
                    keys = l.strip("\n").split("\t")

                    if "GC" in keys and "NORMALIZED_COVERAGE" in keys:
                        # Detail metrics: one line per GC percentage
                        if s_name in self.picard_gc_bias_data or s_name in data:
                            log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f["fn"], s_name))
                        data[s_name] = dict()
                        gc_col = keys.index("GC")
                        cov_col = keys.index("NORMALIZED_COVERAGE")

                    elif "ACCUMULATION_LEVEL" in keys and "GC_DROPOUT" in keys:
                        # Summary metrics - just one line below the header
                        if s_name in self.picard_gc_bias_summary_data or s_name in summary_data:
                            log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
                        summary_data[s_name] = dict()
                        vals = f["f"].readline().rstrip("\n").split("\t")
                        assert len(keys) == len(vals), (keys, vals, f)
                        for k, v in zip(keys, vals):
                            try:
                                summary_data[s_name][k] = float(v)
                            except ValueError:
                                summary_data[s_name][k] = v

        # When there is only one sample, using the file name to extract the sample name.
        # if len(data) <= 1 and len(summary_data):
        #     data = {f["s_name"]: list(data.values())[0]}
        #     summary_data = {f["s_name"]: list(summary_data.values())[0]}

        self.picard_gc_bias_data.update(data)
        self.picard_gc_bias_summary_data.update(summary_data)

        for s_name in set(data.keys()) | set(summary_data.keys()):
            self.add_data_source(f, s_name, section="GcBiasMetrics")

    for s_name in list(self.picard_gc_bias_data.keys()):
        if len(self.picard_gc_bias_data[s_name]) == 0:
            self.picard_gc_bias_data.pop(s_name, None)
            log.debug("Removing {} as no data parsed".format(s_name))

    for s_name in list(self.picard_gc_bias_summary_data.keys()):
        if len(self.picard_gc_bias_summary_data[s_name]) == 0:
            self.picard_gc_bias_summary_data.pop(s_name, None)
            log.debug("Removing {} as no data parsed".format(s_name))

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
