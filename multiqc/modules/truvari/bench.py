#!/usr/bin/env python
""" MultiQC submodule to parse output from truvari bench """
import json
import logging
import os
from collections import OrderedDict

from multiqc.plots import scatter, table
from multiqc.utils.mqc_colour import mqc_colour_scale

# Initialise the logger
log = logging.getLogger(__name__)


class BenchSummary:
    def parse_bench_stats(self):
        """Find truvari bench logs and parse their data"""
        data = {}
        for f in self.find_log_files("truvari/bench"):
            collect_stats = False
            stats = ""
            version = None
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
                    version = line.strip().split(":")[-1]

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

                if version is not None:
                    self.add_software_version(version, f["s_name"])

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
            "description": "Precision of the SV calls. Definition: TP-call / (TP-call + FP)",
            "suffix": "%",
            "modify": lambda x: x * 100,
            "placement": 100,
            "scale": "RdYlGn",
        }
        bench_headers["recall"] = {
            "title": "Recall",
            "description": "Recall of the SV calls. Definition: TP-base / (TP-base + FN)",
            "suffix": "%",
            "modify": lambda x: x * 100,
            "placement": 101,
            "scale": "RdYlGn",
        }
        bench_headers["f1"] = {
            "title": "F1",
            "description": "F1 score of the SV calls. Definition: 2 * ((recall * precision) / (recall + precision))",
            "suffix": "%",
            "modify": lambda x: x * 100,
            "placement": 102,
            "scale": "RdYlGn",
        }
        self.general_stats_addcols(data, bench_headers)

        # Make table
        # Descriptions taken from: https://github.com/ACEnglish/truvari/wiki/bench
        keys = OrderedDict()
        keys["TP-base"] = {
            "title": "TP (base)",
            "description": "Number of the base VCF calls matching the comp VCF",
            "scale": "BuPu",
            "format": "{:,.d}",
            "min": 0,
        }
        keys["TP-call"] = {
            "title": "TP (comp)",
            "description": "Number of the comp VCF calls matching the base VCF",
            "scale": "PuBu",
            "format": "{:,.d}",
            "min": 0,
        }
        keys["FP"] = {
            "title": "FP",
            "description": "Number of the comp VCF calls non-matching the base VCF",
            "scale": "OrRd",
            "format": "{:,.d}",
            "min": 0,
        }
        keys["FN"] = {
            "title": "FN",
            "description": "Number of the base VCF calls non-matching the comp VCF",
            "scale": "BuGn",
            "format": "{:,.d}",
            "min": 0,
        }

        # Reuse info from bench_headers but hide by default
        keys["precision"] = dict(bench_headers["precision"], hidden=True)
        keys["recall"] = dict(bench_headers["recall"], hidden=True)
        keys["f1"] = dict(bench_headers["f1"], hidden=True)

        keys["base cnt"] = {
            "title": "Base calls",
            "description": "Number of calls in the base VCF",
            "scale": "Oranges",
            "format": "{:,.d}",
            "min": 0,
        }
        keys["call cnt"] = {
            "title": "Comp calls",
            "description": "Number of calls in the comp VCF",
            "scale": "Blues",
            "format": "{:,.d}",
            "min": 0,
        }
        keys["TP-call_TP-gt"] = {
            "title": "TP-call with GT match",
            "description": "TP-call with genotype match",
            "hidden": True,
            "scale": "Greens",
            "format": "{:,.d}",
            "min": 0,
        }
        keys["TP-call_FP-gt"] = {
            "title": "TP-call w/o GT match",
            "description": "TP-call without genotype match",
            "hidden": True,
            "scale": "YlGn",
            "format": "{:,.d}",
            "min": 0,
        }
        keys["TP-base_TP-gt"] = {
            "title": "TP-base with GT match",
            "description": "TP-base with genotype match",
            "hidden": True,
            "scale": "RdPu",
            "format": "{:,.d}",
            "min": 0,
        }
        keys["TP-base_FP-gt"] = {
            "title": "TP-base w/o GT match",
            "description": "TP-base without genotype match",
            "hidden": True,
            "scale": "Purples",
            "format": "{:,.d}",
            "min": 0,
        }
        keys["gt_concordance"] = {
            "title": "GT concordance",
            "description": "Genotype concordance. Definition: TP-call with GT / (TP-call with GT + TP-call w/o GT)",
            "hidden": True,
            "scale": "GnBu",
            "suffix": "%",
            "modify": lambda x: x * 100,
        }

        self.add_section(
            name="Truvari bench",
            anchor="truvari-bench",
            description="Concordance statistics parsed from the output from <code>truvari bench</code>.",
            plot=table.plot(data, keys, {"id": "truvari-bench-summary"}),
        )
        # Generate color scale to label samples. The "plot_defaults"
        # scale contains 10 colors so if there are more samples than
        # that, we will reuse colors.
        color_scale = mqc_colour_scale("plot_defaults")

        # Make scatter plot
        scatter_data = {}
        for i, (sample, sample_data) in enumerate(data.items()):
            scatter_data[sample] = {
                "x": sample_data["precision"] * 100.0,
                "y": sample_data["recall"] * 100.0,
                "color": color_scale.get_colour(i, lighten=0.9),
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
            "title": "Truvari bench: precision-recall",
            "tt_label": "Precision: {point.x:>4.1f}%<br/> Recall: {point.y:>6.1f}%",
        }
        self.add_section(
            name="Truvari bench: precision vs. recall",
            anchor="truvari-bench-pre-rec",
            description="Precision vs. recall for each sample. Parsed from the <code>truvari bench</code> output.",
            plot=scatter.plot(scatter_data, scatter_config),
        )

        # Return the number of logs that were found
        return len(data)
