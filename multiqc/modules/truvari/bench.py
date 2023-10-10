#!/usr/bin/env python
""" MultiQC submodule to parse output from truvari bench """
import json
import logging
import os
import re
from collections import OrderedDict

from multiqc.plots import bargraph, scatter, table
from multiqc.utils.mqc_colour import mqc_colour_scale

# Initialise the logger
log = logging.getLogger(__name__)

VERSION_REGEX = r"(\d+\.\d+\.[\d\.\-\w]+)"


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

                # Get version from log
                # Log lines look like:
                #   2022-12-01 18:25:54,715 [INFO] Truvari v4.0.0.dev0+detached
                # or
                #   2022-08-25 00:31:16,781 [INFO] Truvari version: 3.5.0-dev
                if "Truvari" in line:
                    match = re.search(VERSION_REGEX, line)
                    if match:
                        version = match.group(1)

            if stats:
                # Use output directory as sample name
                f["s_name"] = os.path.basename(f["root"])
                f["s_name"] = self.clean_s_name(f["s_name"], f, root=os.path.dirname(f["root"]))

                # Load stats
                try:
                    stats = json.loads(str(stats))
                except json.decoder.JSONDecodeError as e:
                    log.debug(e)
                    log.warning("Could not parse stats from file: {}".format(f["fn"]))

                if f["s_name"] in data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(f["s_name"]))

                # Some stats were renamed in truvari 4.0.0 (commit 6e37058)
                # This renames them back to the old names for backwards compatibility
                if all(key in stats for key in ["TP-call_TP-gt", "TP-call_FP-gt", "TP-call", "call cnt"]):
                    stats["TP-comp_TP-gt"] = stats["TP-call_TP-gt"]
                    stats["TP-comp_FP-gt"] = stats["TP-call_FP-gt"]
                    stats["TP-comp"] = stats["TP-call"]
                    stats["comp cnt"] = stats["call cnt"]

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
            "description": "Precision of the SV calls. Definition: TP-comp / (TP-comp + FP)",
            "suffix": "%",
            "modify": lambda x: x * 100,
            "placement": 100,
            "max": 100,
            "min": 0,
            "scale": "PRGn",
        }
        bench_headers["recall"] = {
            "title": "Recall",
            "description": "Recall of the SV calls. Definition: TP-base / (TP-base + FN)",
            "suffix": "%",
            "modify": lambda x: x * 100,
            "placement": 101,
            "max": 100,
            "min": 0,
            "scale": "BrBG",
        }
        bench_headers["f1"] = {
            "title": "F1",
            "description": "F1 score of the SV calls. Definition: 2 * ((Recall * Precision) / (Recall + Precision))",
            "suffix": "%",
            "modify": lambda x: x * 100,
            "placement": 102,
            "max": 100,
            "min": 0,
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
        keys["TP-comp"] = {
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
        keys["comp cnt"] = {
            "title": "Comp calls",
            "description": "Number of calls in the comp VCF",
            "scale": "Blues",
            "format": "{:,.d}",
            "min": 0,
        }
        keys["TP-comp_TP-gt"] = {
            "title": "TP (comp) with GT match",
            "description": "TP (comp) with genotype match",
            "hidden": True,
            "scale": "Greens",
            "format": "{:,.d}",
            "min": 0,
        }
        keys["TP-comp_FP-gt"] = {
            "title": "TP (comp) w/o GT match",
            "description": "TP (comp) without genotype match",
            "hidden": True,
            "scale": "YlGn",
            "format": "{:,.d}",
            "min": 0,
        }
        keys["TP-base_TP-gt"] = {
            "title": "TP (base) with GT match",
            "description": "TP (base) with genotype match",
            "hidden": True,
            "scale": "RdPu",
            "format": "{:,.d}",
            "min": 0,
        }
        keys["TP-base_FP-gt"] = {
            "title": "TP (base) w/o GT match",
            "description": "TP (base) without genotype match",
            "hidden": True,
            "scale": "Purples",
            "format": "{:,.d}",
            "min": 0,
        }
        keys["gt_concordance"] = {
            "title": "GT concordance",
            "description": "Genotype concordance. Definition: TP-comp with GT / (TP-comp with GT + TP-comp w/o GT)",
            "hidden": True,
            "scale": "GnBu",
            "suffix": "%",
            "modify": lambda x: x * 100,
        }

        self.add_section(
            name="Statistics",
            anchor="truvari-bench",
            description="Performance statistics parsed from the output from `truvari bench`.",
            helptext="""
            Performance metrics from comparison of two VCFs, one truth set ("base") and one to be evaluated ("call").

            - **TP (base)**:	Number of matching calls from the base VCF
            - **TP (comp)**: Number of matching calls from the comp VCF
            - **FP**: Number of non-matching calls from the comp VCF
            - **FN**: Number of non-matching calls from the base VCF
            - **Precision**: Precision metric. Definition: TP-comp / (TP-comp + FP)
            - **Recall**: Recall metric. Definition: TP-base / (TP-base + FN)
            - **F1**: F1 score metric. Definition: 2 * ((Recall * Precision) / (Recall + Precision))
            - **Base calls**: Number of calls in the base vcf. Should be same as TP-base + FN
            - **Comp calls**: Number of calls in the comp vcf. Should be same as TP-comp + FP
            - **TP-comp_TP-gt**: TP (comp) with genotype match
            - **TP-comp_FP-gt**: TP (comp) without genotype match
            - **TP-base_TP-gt**: TP (base) with genotype match
            - **TP-base_FP-gt**: TP (base) without genotype match
            - **GT concordance**: Genotype concordance. Definition: TP-comp_TP-gt / (TP-comp_TP-gt + TP-comp_FP-gt)

            For more information, see the [truvari benwiki](https://github.com/acenglish/truvari/wiki/bench)
            """,
            plot=table.plot(data, keys, {"id": "truvari-bench-summary"}),
        )

        # Make bar graph
        bar_data = [{}, {}]
        for sample, sample_data in data.items():
            bar_data[0][sample] = {
                "TP": sample_data["TP-comp"],
                "FP": sample_data["FP"],
            }
            bar_data[1][sample] = {
                "TP": sample_data["TP-base"],
                "FN": sample_data["FN"],
            }

        bar_categories = [
            {
                "TP": {"name": "TP", "color": "#648FFF"},
                "FP": {"name": "FP", "color": "#DC267F"},
            },
            {
                "TP": {"name": "TP", "color": "#648FFF"},
                "FN": {"name": "FN", "color": "#FFB000"},
            },
        ]

        bar_config = {
            "id": "truvari-bench-classifications-plot",
            "title": "Truvari: Classifications",
            "ylab": "Number of calls",
            "tt_suffix": " calls",
            "data_labels": ["Comp", "Base"],
        }

        self.add_section(
            name="Classifications",
            anchor="truvari-bench-classifications",
            description="Bargraph of SV call classifications parsed from the output from `truvari bench`",
            helptext="""
            Bargraph of SV call classifications from the perspectives of the *comp* and *base* ("truth") VCFs.

            For *comp* calls, the classifications are:

            - **TP**: Comp call matches a base call
            - **FP**: Comp call does not match a base call

            For *base* calls, the classifications are:

            - **TP**: Base call matches a comp call
            - **FN**: Base call does not match a comp call
            """,
            plot=bargraph.plot(bar_data, bar_categories, pconfig=bar_config),
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
            "title": "Truvari: Precision vs. Recall",
            "tt_label": "Precision: {point.x:>4.1f}%<br/> Recall: {point.y:>6.1f}%",
        }
        self.add_section(
            name="Precision vs. Recall",
            anchor="truvari-bench-pre-rec",
            description="Precision vs. Recall for each sample. Parsed from the <code>truvari bench</code> output.",
            helptext="""
            Scatter plot of precision vs. recall for each sample. The precision and recall values are
            calculated as follows:

            - **Precision**: TP-comp / (TP-comp + FP)
            - **Recall**: TP-base / (TP-base + FN)
            """,
            plot=scatter.plot(scatter_data, scatter_config),
        )

        # Return the number of logs that were found
        return len(data)
