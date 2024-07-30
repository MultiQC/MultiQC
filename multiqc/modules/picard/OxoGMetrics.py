"""MultiQC submodule to parse output from Picard OxoGMetrics"""

import logging
from collections import defaultdict
from typing import Dict

from multiqc.modules.picard import util

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(module):
    """Find Picard OxoGMetrics reports and parse their data"""

    # Set up vars
    data_by_sample: Dict = dict()

    # Go through logs and find Metrics
    for f in module.find_log_files("picard/oxogmetrics", filehandles=True):
        # Sample name from input file name by default.
        s_name = f["s_name"]
        keys = None
        context_col = None

        for line in f["f"]:
            maybe_s_name = util.extract_sample_name(
                module,
                line,
                f,
                picard_tool=["CollectOxoGMetrics", "ConvertSequencingArtifactToOxoG"],
                picard_opt=["INPUT", "INPUT_BASE"],
            )
            if maybe_s_name:
                s_name = maybe_s_name
                keys = None
                context_col = None

            if s_name is None:
                continue

            if util.is_line_right_before_table(line, picard_class="CollectOxoGMetrics"):
                keys = f["f"].readline().strip("\n").split("\t")
                context_col = keys.index("CONTEXT")
                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: " f"{s_name}")
                data_by_sample[s_name] = defaultdict()
                module.add_data_source(f, s_name, section="OxoGMetrics")

            elif keys:
                vals = line.strip("\n").split("\t")
                if len(vals) != len(keys) or context_col is None:
                    s_name = None
                    keys = None
                    context_col = None
                    continue

                context = vals[context_col]
                data_by_sample[s_name][context] = dict()
                for i, k in enumerate(keys):
                    k = k.strip()
                    try:
                        val = float(vals[i])
                    except ValueError:
                        val = vals[i].strip()
                    data_by_sample[s_name][context][k] = val

    data_by_sample = module.ignore_samples(data_by_sample)
    if len(data_by_sample) == 0:
        return set()

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed data to a file
    # Collapse into 2D structure with sample_context keys
    print_data = {f"{s}_{c}": v for s in data_by_sample.keys() for c, v in data_by_sample[s].items()}
    module.write_data_file(print_data, "multiqc_picard_OxoGMetrics")

    # Add to general stats table
    general_stats_data: Dict = dict()
    for s_name in data_by_sample:
        general_stats_data[s_name] = dict()
        try:
            val = data_by_sample[s_name]["CCG"]["OXIDATION_ERROR_RATE"]
            general_stats_data[s_name]["CCG_OXIDATION_ERROR_RATE"] = val
        except KeyError:
            log.warning(f"Couldn't find picard CCG oxidation error rate for {s_name}")
    headers = {
        "CCG_OXIDATION_ERROR_RATE": {
            "title": "CCG Oxidation",
            "description": "CCG-CAG Oxidation Error Rate",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "format": "{:,.0f}",
            "scale": "RdYlGn-rev",
            "modify": lambda x: util.multiply_hundred(x),
        }
    }
    module.general_stats_addcols(general_stats_data, headers, namespace="OxoGMetrics")

    # Return the number of detected samples to the parent module
    return data_by_sample.keys()
