#!/usr/bin/env python
""" MultiQC submodule to parse output from truvari bench """
import json
import os
import logging
import random
from collections import OrderedDict

import matplotlib.colors as mcolors

from multiqc.plots import table, scatter

# Initialise the logger
log = logging.getLogger(__name__)


class BenchSummary:
    def parse_bench_stats(self):
        """Find truvari bench logs and parse their data"""
        data = {}
        versions = set()  # Not used at the moment
        for f in self.find_log_files("truvari/bench"):
            collect_stats = False
            stats = ""
            for line in f["f"].splitlines():
                if "Stats:" in line:
                    collect_stats = True
                    stats = "{\n"
                    continue

                if collect_stats:
                    stats += line
                    if line.startswith("}"):
                        collect_stats = False

                if "Truvari version:" in line:
                    versions.add(line.strip().split(":")[-1])

            if stats:
                # Use output directory as sample name
                f["s_name"] = os.path.basename(f["root"])
                f["s_name"] = self.clean_s_name(f["s_name"], f, root=os.path.dirname(f["root"]))

                # Load stats
                stats = json.loads(str(stats))

                # Convert entries with NaN to 0.0
                stats = {k: v if v != "NaN" else 0.0 for k, v in stats.items()}

                if f["s_name"] in data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))

                self.add_data_source(f, section="bench")
                data[f["s_name"]] = stats

        # Filter to strip out ignored sample names
        data = self.ignore_samples(data)

        # Return if no samples
        if len(data) == 0:
            return len(data)

        # Write parsed report data to a file (restructure first)
        self.write_data_file(data, "multiqc_truvari_bench")

        # General Stats Table
        bench_headers = dict()
        bench_headers["precision"] = {
            "title": "Precision",
            "description": "Precision for SV calls. Definition: TP-call / (TP-call + FP)",
            "min": -0.001,
            "max": 1.0,
            "format": "{:.1%}",
            "placement": 100,
            "scale": "RdYlGn",
        }
        bench_headers["recall"] = {
            "title": "Recall",
            "description": "Recall for SV calls. Definition: TP-base / (TP-base + FN)",
            "min": -0.001,
            "max": 1.0,
            "format": "{:.1%}",
            "placement": 101,
            "scale": "RdYlGn",
        }
        bench_headers["f1"] = {
            "title": "F1",
            "description": "F1 score for SV calls. Definition: 2 * ((recall * precision) / (recall + precision))",
            "min": -0.001,
            "max": 1.0,
            "format": "{:.1%}",
            "placement": 102,
            "scale": "RdYlGn",
        }
        self.general_stats_addcols(data, bench_headers)

        # Make table
        # Descriptions taken from: https://github.com/ACEnglish/truvari/wiki/bench
        keys = OrderedDict()
        keys["TP-base"] = {
            "title": "TP-base",
            "description": "Number of matching calls from the base vcf",
            "scale": "BuPu",
            "format": "{:,.0f}",
        }
        keys["TP-call"] = {
            "title": "TP-call",
            "description": "Number of matching calls from the comp vcf",
            "scale": "PuBu",
            "format": "{:,.0f}"
        }
        keys["FP"] = {
            "title": "FP",
            "description": "Number of non-matching calls from the comp vcf",
            "scale": "OrRd",
            "format": "{:,.0f}"
        }
        keys["FN"] = {
            "title": "FN",
            "description": "Number of non-matching calls from the base vcf",
            "scale": "BuGn",
            "format": "{:,.0f}"
        }

        # Reuse info from bench_headers but hide by default
        keys["precision"] = dict(bench_headers["precision"], hidden=True)
        keys["recall"] = dict(bench_headers["recall"], hidden=True)
        keys["f1"] = dict(bench_headers["f1"], hidden=True)

        keys["base cnt"] = {
            "title": "nBase",
            "description": "Number of calls in the base vcf",
            "scale": "Oranges",
            "format": "{:,.0f}"
        }
        keys["call cnt"] = {
            "title": "nCall",
            "description": "Number of calls in the comp vcf",
            "scale": "Blues",
            "format": "{:,.0f}"
        }
        keys["TP-call_TP-gt"] = {
            "title": "TP-call_TP-gt",
            "description": "TP-call with genotype match",
            "hidden": True,
            "scale": "Greens",
            "format": "{:,.0f}"
        }
        keys["TP-call_FP-gt"] = {
            "title": "TP-call_FP-gt",
            "description": "TP-call without genotype match",
            "hidden": True,
            "scale": "YlGn",
            "format": "{:,.0f}"
        }
        keys["TP-base_TP-gt"] = {
            "title": "TP-base_TP-gt",
            "description": "TP-base with genotype match",
            "hidden": True,
            "scale": "RdPu",
            "format": "{:,.0f}"
        }
        keys["TP-base_FP-gt"] = {
            "title": "TP-call_FP-gt",
            "description": "TP-base without genotype match",
            "hidden": True,
            "scale": "Purples",
            "format": "{:,.0f}"
        }
        keys["gt_concordance"] = {
            "title": "gt Concordance",
            "description": "Genotype concordance. Definition: TP-call_TP-gt / (TP-call_TP-gt + TP-call_FP-gt)",
            "hidden": True,
            "scale": "GnBu",
            "format": "{:.1%}",
        }

        self.add_section(
            name="Truvari bench",
            anchor="truvari-bench",
            description="This module parses the output from <code>truvari bench</code>.",
            plot=table.plot(data, keys, {"id": "truvari-bench-summary"}),
        )
        # Generate shuffled list of colors to label samples
        colors = list(mcolors.CSS4_COLORS)
        random.shuffle(colors)

        # Check that not too many samples for scatter plot
        if len(data) > len(colors):
            log.debug("Too many samples ({}) to generate scatter plot.".format(len(data)))
            return len(data)

        # Make scatter plot
        scatter_data = {}
        for (sample, sample_data), color in zip(data.items(), colors):
            scatter_data[sample] = {
                "x": sample_data["precision"] * 100,
                "y": sample_data["recall"] * 100,
                "color": mcolors.to_hex(color),
            }

        scatter_config = {
            "marker_size": 5,
            "height": 560,  # increase height slightly to fit title.
            "ymax": 100,
            "ymin": 0,
            "xmax": 100,
            "xmin": 0,
            "xlab": "Precision (%)",
            "ylab": "Recall (%)",
            "square": True,
            "title": "Truvari bench: Precision-Recall",
            "tt_label": "Precision: {point.x:.1f}%, Recall: {point.y:.1f}%",
        }
        self.add_section(
            name="Truvari bench - Precision-Recall",
            anchor="truvari-bench-pre-rec",
            description="This module parses the output from <code>truvari bench</code>.",
            plot=scatter.plot(scatter_data, scatter_config),
        )

        # Return the number of logs that were found
        return len(data)
