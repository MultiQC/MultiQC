"""MultiQC submodule to parse output from Picard BaseDistributionByCycleMetrics"""

import logging
from collections import defaultdict
from typing import Dict, List

from multiqc.modules.picard import util
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(module):
    """Find Picard BaseDistributionByCycleMetrics reports and parse their data"""

    data_by_sample = dict()
    samplestats_by_sample = dict()

    # Go through logs and find Metrics
    for f in module.find_log_files("picard/basedistributionbycycle", filehandles=True):
        # Sample name from input file name by default.
        s_name = f["s_name"]

        # A file can be concatenated from multiple samples, so we need to keep track of
        # the current sample name and header.
        keys = None
        data_by_read_end: Dict[str, Dict] = defaultdict(dict)
        max_cycle_r1 = 0

        def _finalize_sample(data_by_read_end, s_name: str):
            """
            Populate `data_by_sample` and `samplestats_by_sample`
            """
            # set up the set of s_names
            if 2 in set(data_by_read_end):
                s_names = {1: f"{s_name}_R1", 2: f"{s_name}_R2"}
            else:
                s_names = {1: s_name}

            for read_end in s_names:
                data_by_cycle = data_by_read_end[read_end]
                s_name = s_names[read_end]
                data_by_sample[s_name] = data_by_cycle
                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
                module.add_data_source(f, s_name, section="BaseDistributionByCycle")
                sample_stats = {
                    "sum_pct_a": 0,
                    "sum_pct_c": 0,
                    "sum_pct_g": 0,
                    "sum_pct_t": 0,
                    "sum_pct_n": 0,
                    "cycle_count": 0,
                }
                samplestats_by_sample[s_name] = sample_stats
                for c, row in data_by_cycle.items():
                    pct_a, pct_c, pct_g, pct_t, pct_n = row
                    sample_stats["sum_pct_a"] += pct_a
                    sample_stats["sum_pct_c"] += pct_c
                    sample_stats["sum_pct_g"] += pct_g
                    sample_stats["sum_pct_t"] += pct_t
                    sample_stats["sum_pct_n"] += pct_n
                sample_stats["cycle_count"] += len(data_by_cycle.keys())

        for line in f["f"]:
            maybe_s_name = util.extract_sample_name(module, line, f, picard_tool="CollectBaseDistributionByCycle")
            if maybe_s_name:
                # Starts information for a new sample
                s_name = maybe_s_name
                keys = None

            if s_name is None:
                continue

            if util.is_line_right_before_table(line, picard_class="BaseDistributionByCycleMetrics"):
                keys = f["f"].readline().strip("\n").split("\t")
                assert keys == ["READ_END", "CYCLE", "PCT_A", "PCT_C", "PCT_G", "PCT_T", "PCT_N"]

            elif keys:
                raw_vals = line.strip("\n").split("\t")
                if len(raw_vals) != len(keys):
                    # Finalize previous sample
                    if data_by_read_end:
                        _finalize_sample(data_by_read_end, s_name)
                    # Reset for next sample
                    max_cycle_r1 = 0
                    s_name = None
                    data_by_read_end = defaultdict(dict)
                    continue

                vals = []
                for v in raw_vals:
                    try:
                        v = float(v)
                    except ValueError:
                        pass
                    vals.append(v)

                read_end, cycle, pct_a, pct_c, pct_g, pct_t, pct_n = vals
                try:
                    cycle = int(cycle)
                except ValueError:
                    continue
                if read_end == 1.0:
                    max_cycle_r1 = max(max_cycle_r1, cycle)
                else:
                    cycle -= max_cycle_r1
                data_by_read_end[read_end][cycle] = (pct_a, pct_c, pct_g, pct_t, pct_n)

        if data_by_read_end and s_name:
            _finalize_sample(data_by_read_end, s_name)

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)
    samplestats_by_sample = module.ignore_samples(samplestats_by_sample)
    if len(data_by_sample) == 0:
        return set()

    # Calculate summed mean values for all read orientations
    for s_name, v in samplestats_by_sample.items():
        v["mean_pct_a"] = v["sum_pct_a"] / v["cycle_count"]
        v["mean_pct_c"] = v["sum_pct_c"] / v["cycle_count"]
        v["mean_pct_g"] = v["sum_pct_g"] / v["cycle_count"]
        v["mean_pct_t"] = v["sum_pct_t"] / v["cycle_count"]

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed data to a file
    module.write_data_file(samplestats_by_sample, "multiqc_picard_baseContent")

    # Plot the data and add section
    pconfig = {
        "id": "picard_base_distribution_by_cycle",
        "title": "Picard: Base Distribution",
        "ylab": "%",
        "xlab": "Cycle #",
        "x_decimals": False,
        "tt_label": "<b>cycle {point.x}</b>: {point.y:.2f} %",
        "ymax": 100,
        "ymin": 0,
        "data_labels": [
            {"name": "% Adenine", "ylab": "% Adenine"},
            {"name": "% Cytosine", "ylab": "% Cytosine"},
            {"name": "% Guanine", "ylab": "% Guanine"},
            {"name": "% Thymine", "ylab": "% Thymine"},
            {"name": "% Undetermined", "ylab": "% Undetermined"},
        ],
    }

    # build list of line graphs
    linegraph_data: List[Dict] = [{}, {}, {}, {}, {}]
    for s_name, cycles in data_by_sample.items():
        for lg, index in zip(linegraph_data, range(5)):
            lg[s_name] = {cycle: tup[index] for cycle, tup in cycles.items()}

    module.add_section(
        name="Base Distribution",
        anchor="picard-base-distribution-by-cycle",
        description="Plot shows the distribution of bases by cycle.",
        plot=linegraph.plot(linegraph_data, pconfig),
    )

    # Return the number of detected samples to the parent module
    return data_by_sample.keys()
