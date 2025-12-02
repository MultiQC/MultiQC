"""MultiQC module to parse output from Wittyer"""

import json
import logging
from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Wittyer is a tool for benchmarking structural variant (SV) calls against a truth set.
    This module parses the JSON stats output from Wittyer to provide event-level statistics different variant types.
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Wittyer",
            anchor="wittyer",
            target="Wittyer",
            href="https://github.com/Illumina/Wittyer",
            info="A tool for benchmarking structural variant calls against a truth set.",
            doi="10.1093/bioinformatics/btaa397",
        )

        # Find and parse JSON files
        self.wittyer_data = dict()

        for f in self.find_log_files("wittyer", filehandles=True):
            try:
                data = json.load(f["f"])
                self.parse_wittyer_json(data, f)
            except json.JSONDecodeError:
                log.debug(f"Could not parse Wittyer JSON: {f['fn']}")
                continue

        # Filter to strip out ignored sample names
        self.wittyer_data = self.ignore_samples(self.wittyer_data)

        if len(self.wittyer_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.wittyer_data)} Wittyer reports")

        # Write parsed data to a file
        self.write_data_file(self.wittyer_data, "multiqc_wittyer")

        # Add sections to the report
        self.add_variant_type_section()

    def parse_wittyer_json(self, data: Dict, f: Dict) -> None:
        """Parse Wittyer JSON output"""

        # Get sample name from PerSampleStats
        if "PerSampleStats" not in data or len(data["PerSampleStats"]) == 0:
            return

        sample_stats = data["PerSampleStats"][0]
        query_sample = sample_stats.get("QuerySampleName", f["s_name"])
        truth_sample = sample_stats.get("TruthSampleName", "unknown")

        # Use query sample name as the main identifier
        s_name = self.clean_s_name(query_sample, f)

        if s_name in self.wittyer_data:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")

        # Store the parsed data
        self.wittyer_data[s_name] = {
            "command": data.get("Command", ""),
            "truth_sample": truth_sample,
            "event_precision": data.get("EventPrecision", 0),
            "event_recall": data.get("EventRecall", 0),
            "event_fscore": data.get("EventFscore", 0),
            "overall_stats": sample_stats.get("OverallStats", []),
            "detailed_stats": sample_stats.get("DetailedStats", []),
        }

    def add_variant_type_section(self) -> None:
        """Add separate table sections for each main variant type"""

        # Define the main variant types to show (6 most common/important)
        main_variant_types = [
            ("Deletion", "Deletions"),
            ("Insertion", "Insertions"),
            ("Duplication", "Duplications"),
            ("Inversion", "Inversions"),
            ("CopyNumberGain", "Copy Number Gains"),
            ("CopyNumberLoss", "Copy Number Losses"),
        ]

        # Organize data by variant type
        variant_data = {}
        for variant_type, _ in main_variant_types:
            variant_data[variant_type] = {}

        # Parse data from all samples
        for s_name, data in self.wittyer_data.items():
            for variant_stat in data["detailed_stats"]:
                variant_type = variant_stat.get("VariantType", "Unknown")

                # Only process main variant types
                if variant_type not in variant_data:
                    continue

                overall_stats = variant_stat.get("OverallStats", [])

                # Get event-level stats
                for stat in overall_stats:
                    if stat.get("StatsType") == "Event":
                        # Handle NaN values
                        recall = stat.get("Recall", 0)
                        precision = stat.get("Precision", 0)
                        fscore = stat.get("Fscore", 0)

                        if recall == "NaN" or recall is None:
                            recall = 0
                        else:
                            recall = float(recall) * 100

                        if precision == "NaN" or precision is None:
                            precision = 0
                        else:
                            precision = float(precision) * 100

                        if fscore == "NaN" or fscore is None:
                            fscore = 0
                        else:
                            fscore = float(fscore) * 100

                        variant_data[variant_type][s_name] = {
                            "precision": precision,
                            "recall": recall,
                            "fscore": fscore,
                            "ttp": stat.get("TruthTpCount", 0),
                            "tfn": stat.get("TruthFnCount", 0),
                            "total_truth": stat.get("TruthTotalCount", 0),
                            "qtp": stat.get("QueryTpCount", 0),
                            "qfp": stat.get("QueryFpCount", 0),
                            "total_query": stat.get("QueryTotalCount", 0),
                        }
                        break

        # Create a separate table for each variant type
        for variant_type, display_name in main_variant_types:
            if variant_type not in variant_data or not variant_data[variant_type]:
                # Skip if no data for this variant type
                continue

            # Prepare table data for this variant type
            table_data = variant_data[variant_type]

            # Define headers
            table_headers = {
                "precision": {
                    "title": "Precision (%)",
                    "description": f"Precision for {display_name.lower()}",
                    "scale": "Blues",
                    "format": "{:.2f}",
                },
                "recall": {
                    "title": "Recall (%)",
                    "description": f"Recall for {display_name.lower()}",
                    "scale": "Blues",
                    "format": "{:.2f}",
                },
                "fscore": {
                    "title": "F1-score (%)",
                    "description": f"F-score for {display_name.lower()}",
                    "scale": "Blues",
                    "format": "{:.2f}",
                },
                "ttp": {
                    "title": "True Positives - Truth",
                    "description": f"Number of correctly identified {display_name.lower()}",
                    "scale": "Greens",
                    "format": "{:,.0f}",
                },
                "tfn": {
                    "title": "False Negatives - Truth",
                    "description": f"Number of missed {display_name.lower()}",
                    "scale": "Reds",
                    "format": "{:,.0f}",
                },
                "total_truth": {
                    "title": "Total Truth Events",
                    "description": f"Total number of {display_name.lower()} in truth set",
                    "scale": "Greys",
                    "format": "{:,.0f}",
                },
                "qtp": {
                    "title": "True Positives - Query",
                    "description": f"Number of correctly called {display_name.lower()}",
                    "scale": "Greens",
                    "format": "{:,.0f}",
                },
                "qfp": {
                    "title": "False Positives - Query",
                    "description": f"Number of incorrectly called {display_name.lower()}",
                    "scale": "Reds",
                    "format": "{:,.0f}",
                },
                "total_query": {
                    "title": "Total Query Events",
                    "description": f"Total number of {display_name.lower()} called in query set",
                    "scale": "Greys",
                    "format": "{:,.0f}",
                },
            }

            # Add section for this variant type
            self.add_section(
                name=f"{display_name}",
                anchor=f"wittyer-{variant_type.lower()}",
                description=f"Event-level benchmarking metrics for {display_name.lower()}",
                plot=table.plot(
                    table_data,
                    table_headers,
                    pconfig={"id": f"wittyer-{variant_type.lower()}-plot", "title": f"wittyer: {variant_type.lower()}"},
                ),
            )
