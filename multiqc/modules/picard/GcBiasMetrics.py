"""MultiQC submodule to parse output from Picard GcBiasMetrics"""

import logging
from typing import Dict

from multiqc.modules.picard import util
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(module):
    """
    Find Picard GcBiasMetrics reports and parse their data. There are two types of
    GC bias files:
    * detail: one set of metrics per GC percentage,
    * summary: one set of metrics per input or sample.
    In Picard, they are generated separately with GcBiasDetailMetrics and
    GcBiasSummaryMetrics tools correspondingly. In Sentieon, they are produced with
    the same command, and the summary is only written when the `--summary` option
    is provided.
    """

    data_by_sample: Dict[str, Dict] = dict()
    summary_data_by_sample: Dict[str, Dict] = dict()

    # Go through logs and find Metrics
    for f in module.find_log_files("picard/gcbias", filehandles=True):
        # Sample name from input file name by default.
        s_name = f["s_name"]
        gc_col = None
        cov_col = None

        for line in f["f"]:
            maybe_s_name = util.extract_sample_name(
                module,
                line,
                f,
                picard_tool="CollectGcBiasMetrics",
                sentieon_algo="GCBias",
            )
            if maybe_s_name:
                s_name = maybe_s_name

            if s_name is None:
                continue

            if util.is_line_right_before_table(
                line, picard_class=["GcBiasDetailMetrics", "GcBiasSummaryMetrics"], sentieon_algo="GCBias"
            ):
                # Get header - find columns with the data we want
                line = f["f"].readline()
                keys = line.strip("\n").split("\t")

                if "GC" in keys and "NORMALIZED_COVERAGE" in keys:
                    # Detail metrics: one line per GC percentage
                    if s_name in data_by_sample:
                        log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
                    data_by_sample[s_name] = dict()
                    gc_col = keys.index("GC")
                    cov_col = keys.index("NORMALIZED_COVERAGE")

                elif "ACCUMULATION_LEVEL" in keys and "GC_DROPOUT" in keys:
                    # Summary metrics - just one line below the header
                    if s_name in summary_data_by_sample:
                        log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
                    summary_data_by_sample[s_name] = dict()
                    vals = f["f"].readline().rstrip("\n").split("\t")
                    if len(keys) != len(vals):
                        s_name = None
                        continue

                    for k, v in zip(keys, vals):
                        try:
                            summary_data_by_sample[s_name][k] = float(v)
                        except ValueError:
                            summary_data_by_sample[s_name][k] = v

            elif gc_col is not None and cov_col is not None:
                try:
                    # Note that GC isn't always the first column.
                    s = line.strip("\n").split("\t")
                    data_by_sample[s_name][int(s[gc_col])] = float(s[cov_col])
                except IndexError:
                    gc_col = None
                    cov_col = None

        for s_name in set(data_by_sample.keys()) | set(summary_data_by_sample.keys()):
            module.add_data_source(f, s_name, section="GcBiasMetrics")

    for s_name in list(data_by_sample.keys()):
        if len(data_by_sample[s_name]) == 0:
            data_by_sample.pop(s_name, None)
            log.debug(f"Removing {s_name} as no data parsed")

    for s_name in list(summary_data_by_sample.keys()):
        if len(summary_data_by_sample[s_name]) == 0:
            summary_data_by_sample.pop(s_name, None)
            log.debug(f"Removing {s_name} as no data parsed")

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)
    summary_data_by_sample = module.ignore_samples(summary_data_by_sample)

    samples = data_by_sample.keys() | summary_data_by_sample.keys()
    if not samples:
        return set()

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    if data_by_sample:
        # Plot the graph
        pconfig = {
            "id": f"{module.anchor}_gcbias_plot",
            "title": f"{module.name}: GC Coverage Bias",
            "ylab": "Normalized Coverage",
            "xlab": "% GC",
            "xmin": 0,
            "xmax": 100,
            "ymin": 0,
            "y_clipmax": 10,
            "tt_label": "<b>{point.x} %GC</b>: {point.y:.2f}",
            "y_lines": [
                {"value": 1, "color": "#999999", "width": 2, "dash": "longdash"},
            ],
        }
        module.add_section(
            name="GC Coverage Bias",
            anchor=f"{module.anchor}-gcbias",
            description="This plot shows bias in coverage across regions of the genome with varying GC content."
            " A perfect library would be a flat line at <code>y = 1</code>.",
            plot=linegraph.plot(data_by_sample, pconfig),
        )

    if summary_data_by_sample:
        # Write parsed summary data to a file
        module.write_data_file(summary_data_by_sample, f"multiqc_{module.anchor}_gcbias")

    # Return the number of detected samples to the parent module
    return samples
