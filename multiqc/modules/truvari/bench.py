""" MultiQC submodule to parse output from truvari bench """

import json
import logging
import os
import re

from multiqc.plots import bargraph, scatter

# Initialise the logger
log = logging.getLogger(__name__)

VERSION_REGEX = r"(\d+\.\d+\.[\d\.\-\w]+)"


class BenchSummary:
    def parse_bench_stats(self):
        """Find truvari bench logs and parse their data"""
        data = {}
        for f in self.find_log_files("truvari/bench"):
            inside_json_block = False
            json_text = ""
            version = None
            for line in f["f"].splitlines():
                if "Stats:" in line:
                    inside_json_block = True
                    json_text = "{\n"
                    continue

                if inside_json_block:
                    json_text += line
                    if line.startswith("}"):
                        inside_json_block = False

                # Get version from log
                # Log lines look like:
                #   2022-12-01 18:25:54,715 [INFO] Truvari v4.0.0.dev0+detached
                # or
                #   2022-08-25 00:31:16,781 [INFO] Truvari version: 3.5.0-dev
                if "Truvari" in line:
                    match = re.search(VERSION_REGEX, line)
                    if match:
                        version = match.group(1)

            if not json_text:
                log.warning(f"Could not find the 'Stats' JSON block in file: {f['fn']}")
                continue

            # Load stats
            try:
                stats = json.loads(str(json_text))
            except json.decoder.JSONDecodeError as e:
                log.debug(e)
                log.warning(f"Could not parse the 'Stats' JSON block in file: {f['fn']}")
                continue

            # Use output directory as sample name
            f["s_name"] = os.path.basename(f["root"])
            f["s_name"] = self.clean_s_name(f["s_name"], f, root=os.path.dirname(f["root"]))
            if f["s_name"] in data:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")

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
        bench_headers["gt_concordance"] = {
            "title": "GT concordance",
            "description": "Genotype concordance. Definition: TP-comp with GT / (TP-comp with GT + TP-comp w/o GT)",
            "scale": "GnBu",
            "suffix": "%",
            "placement": 103,
            "max": 100,
            "min": 0,
            "hidden": True,
            "modify": lambda x: x * 100,
        }

        self.general_stats_addcols(data, bench_headers)

        # Make bar graph
        bar_data = [{}, {}, {}, {}]
        for sample, sample_data in data.items():
            # Comp
            bar_data[0][sample] = {
                "TP": sample_data["TP-comp"],
                "FP": sample_data["FP"],
            }
            # Base
            bar_data[1][sample] = {
                "TP": sample_data["TP-base"],
                "FN": sample_data["FN"],
            }
            # TP (comp) GT
            bar_data[2][sample] = {
                "Match": sample_data["TP-comp_TP-gt"],
                "No match": sample_data["TP-comp_FP-gt"],
            }
            # TP (base) GT
            bar_data[3][sample] = {
                "Match": sample_data["TP-base_TP-gt"],
                "No match": sample_data["TP-base_FP-gt"],
            }

        bar_categories = [
            # Comp
            {
                "TP": {"name": "TP", "color": "#648FFF"},
                "FP": {"name": "FP", "color": "#DC267F"},
            },
            # Base
            {
                "TP": {"name": "TP", "color": "#648FFF"},
                "FN": {"name": "FN", "color": "#FFB000"},
            },
            # TP (comp) GT
            {
                "Match": {"name": "Match", "color": "#785EF0"},
                "No match": {"name": "No match", "color": "#FE6100"},
            },
            # TP (base) GT
            {
                "Match": {"name": "Match", "color": "#785EF0"},
                "No match": {"name": "No match", "color": "#FE6100"},
            },
        ]

        bar_config = {
            "id": "truvari-bench-classifications-plot",
            "title": "Truvari: Classifications",
            "ylab": "Counts",
            "tt_suffix": " calls",
            "data_labels": [
                {"name": "Comp", "ylab": "Counts"},
                {"name": "Base", "ylab": "Counts"},
                {"name": "TP (comp) GT", "ylab": "Counts"},
                {"name": "TP (base) GT", "ylab": "Counts"},
            ],
        }

        self.add_section(
            name="Classifications",
            anchor="truvari-bench-classifications",
            description="Bargraph of SV call classifications parsed from the output from `truvari bench`",
            helptext="""
            Bargraph of SV call classifications from the perspectives of the *comp* and
            *base* ("truth") VCFs. Four different groups of calls are shown:

            #### Comp
            Compares TP and FP calls in the *comp* VCF relative the *base*. The
            classifications are:

            - **TP**: Comp call matches a base call
            - **FP**: Comp call does not match a base call

            #### Base
            Compares TP and FN calls in the *base* VCF relative the *comp*. The
            classifications are:

            - **TP**: Base call matches a comp call
            - **FN**: Base call does not match a comp call

            #### TP (comp) GT
            Compares TP calls in the *comp* VCF relative to the *base* VCF for
            matching genotypes. The classifications are:

            - **Match**: TP call in comp has matching genotype in base
            - **No match**: TP call in comp does not have matching genotype in base

            #### TP (base) GT
            Compares TP calls in the *base* VCF relative to the *comp* VCF for
            matching genotypes. The classifications are:

            - **Match**: TP call in base has matching genotype in comp
            - **No match**: TP call in base does not have matching genotype in comp

            For more information, see the [truvari bench wiki](https://github.com/acenglish/truvari/wiki/bench)
            """,
            plot=bargraph.plot(bar_data, bar_categories, pconfig=bar_config),
        )

        # Make scatter plot
        scatter_data = {}
        for i, (sample, sample_data) in enumerate(data.items()):
            scatter_data[sample] = {
                "x": sample_data["precision"] * 100.0,
                "y": sample_data["recall"] * 100.0,
            }

        scatter_config = {
            "id": "truvari-bench-pre-rec-plot",
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
            "tt_label": "{point.x:.1f}% precision<br/>{point.y:.1f}% recall",
        }
        self.add_section(
            name="Precision vs. Recall",
            anchor="truvari-bench-pre-rec",
            description="Precision vs. Recall for each sample. Parsed from the output of `truvari bench`",
            helptext="""
            Scatter plot of precision vs. recall comparing SV calls between two VCFs,
            one truth set ("base") and one to be evaluated ("comp"). The precision and
            recall values are calculated as follows:

            - **Precision**: TP-comp / (TP-comp + FP)
            - **Recall**: TP-base / (TP-base + FN)

            The TP, FP, and FN values are intrun defined as follows:

            - **TP (base)**: Number of matching calls from the base ('truth') VCF
            - **TP (comp)**: Number of matching calls from the comp VCF
            - **FP**: Number of non-matching calls from the comp VCF
            - **FN**: Number of non-matching calls from the base ('truth') VCF

            For more information, see the [truvari bench wiki](https://github.com/acenglish/truvari/wiki/bench)
            """,
            plot=scatter.plot(scatter_data, scatter_config),
        )

        # Return the number of logs that were found
        return len(data)
