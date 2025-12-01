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
    This module parses the JSON stats output from Wittyer to provide:
    - Overall precision, recall, and F-score metrics
    - Per-variant-type statistics (deletions, insertions, duplications, etc.)
    - Event and base-level statistics
    - Size-binned performance metrics
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
        self.add_overall_stats_section()
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
            "detailed_stats": sample_stats.get("DetailedStats", [])
        }

    def add_overall_stats_section(self) -> None:
        """Add section showing overall performance metrics"""
        
        # Prepare data for general stats table
        general_stats_data = {}
        general_stats_headers = {}
        
        for s_name, data in self.wittyer_data.items():
            general_stats_data[s_name] = {
                "event_precision": data["event_precision"] * 100,
                "event_recall": data["event_recall"] * 100,
                "event_fscore": data["event_fscore"] * 100,
            }
        
        general_stats_headers["event_precision"] = {
            "title": "Precision",
            "description": "Event-level precision (%)",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "RdYlGn",
            "format": "{:,.2f}",
        }
        general_stats_headers["event_recall"] = {
            "title": "Recall",
            "description": "Event-level recall (%)",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "RdYlGn",
            "format": "{:,.2f}",
        }
        general_stats_headers["event_fscore"] = {
            "title": "F-score",
            "description": "Event-level F-score (%)",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "RdYlGn",
            "format": "{:,.2f}",
        }
        
        self.general_stats_addcols(general_stats_data, general_stats_headers)
        
        # Create detailed stats table
        table_data = {}
        for s_name, data in self.wittyer_data.items():
            table_data[s_name] = {}
            for stat in data["overall_stats"]:
                stats_type = stat.get("StatsType", "")
                prefix = f"{stats_type}_"
                
                table_data[s_name][f"{prefix}tp"] = stat.get("TruthTpCount", 0)
                table_data[s_name][f"{prefix}fn"] = stat.get("TruthFnCount", 0)
                table_data[s_name][f"{prefix}fp"] = stat.get("QueryFpCount", 0)
                table_data[s_name][f"{prefix}recall"] = stat.get("Recall", 0) * 100
                table_data[s_name][f"{prefix}precision"] = stat.get("Precision", 0) * 100
                table_data[s_name][f"{prefix}fscore"] = stat.get("Fscore", 0) * 100
        
        table_headers = {
            "Event_tp": {
                "title": "Event TP",
                "description": "True Positive events",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "Event_fn": {
                "title": "Event FN",
                "description": "False Negative events",
                "scale": "Reds",
                "format": "{:,.0f}",
            },
            "Event_fp": {
                "title": "Event FP",
                "description": "False Positive events",
                "scale": "Reds",
                "format": "{:,.0f}",
            },
            "Event_recall": {
                "title": "Event Recall",
                "description": "Event recall (%)",
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.2f}",
                "max": 100,
            },
            "Event_precision": {
                "title": "Event Precision",
                "description": "Event precision (%)",
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.2f}",
                "max": 100,
            },
            "Event_fscore": {
                "title": "Event F-score",
                "description": "Event F-score (%)",
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.2f}",
                "max": 100,
            },
            "Base_tp": {
                "title": "Base TP",
                "description": "True Positive bases",
                "scale": "Greens",
                "format": "{:,.0f}",
            },
            "Base_fn": {
                "title": "Base FN",
                "description": "False Negative bases",
                "scale": "Reds",
                "format": "{:,.0f}",
            },
            "Base_fp": {
                "title": "Base FP",
                "description": "False Positive bases",
                "scale": "Reds",
                "format": "{:,.0f}",
            },
            "Base_recall": {
                "title": "Base Recall",
                "description": "Base-level recall (%)",
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.2f}",
                "max": 100,
            },
            "Base_precision": {
                "title": "Base Precision",
                "description": "Base-level precision (%)",
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.2f}",
                "max": 100,
            },
            "Base_fscore": {
                "title": "Base F-score",
                "description": "Base-level F-score (%)",
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.2f}",
                "max": 100,
            },
        }
        
        self.add_section(
            name="Overall Statistics",
            anchor="wittyer-overall",
            description="Overall event and base-level statistics for structural variant benchmarking",
            plot=table.plot(table_data, table_headers),
        )

    def add_variant_type_section(self) -> None:
        """Add section showing per-variant-type statistics"""
        
        # Prepare data for table - showing all metrics by variant type
        table_data = {}
        
        for s_name, data in self.wittyer_data.items():
            table_data[s_name] = {}
            
            for variant_stat in data["detailed_stats"]:
                variant_type = variant_stat.get("VariantType", "Unknown")
                overall_stats = variant_stat.get("OverallStats", [])
                
                # Get event-level stats
                for stat in overall_stats:
                    if stat.get("StatsType") == "Event":
                        # Create prefixed keys for this variant type
                        prefix = variant_type
                        
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
                        
                        table_data[s_name][f"{prefix}_tp"] = stat.get("TruthTpCount", 0)
                        table_data[s_name][f"{prefix}_fn"] = stat.get("TruthFnCount", 0)
                        table_data[s_name][f"{prefix}_fp"] = stat.get("QueryFpCount", 0)
                        table_data[s_name][f"{prefix}_recall"] = recall
                        table_data[s_name][f"{prefix}_precision"] = precision
                        table_data[s_name][f"{prefix}_fscore"] = fscore
                        break
        
        # Define headers for all variant types
        variant_types = [
            "Deletion",
            "Insertion", 
            "Duplication",
            "Inversion",
            "CopyNumberGain",
            "CopyNumberLoss",
            "CopyNumberReference",
            "CopyNumberTandemRepeat",
            "IntraChromosomeBreakend",
            "TranslocationBreakend"
        ]
        
        table_headers = {}
        
        for variant_type in variant_types:
            # Add headers for this variant type
            table_headers[f"{variant_type}_tp"] = {
                "title": f"{variant_type} TP",
                "description": f"True Positive {variant_type} events",
                "scale": "Greens",
                "format": "{:,.0f}",
                "hidden": True,  # Hide by default, can be shown
            }
            table_headers[f"{variant_type}_fn"] = {
                "title": f"{variant_type} FN",
                "description": f"False Negative {variant_type} events",
                "scale": "Reds",
                "format": "{:,.0f}",
                "hidden": True,
            }
            table_headers[f"{variant_type}_fp"] = {
                "title": f"{variant_type} FP",
                "description": f"False Positive {variant_type} events",
                "scale": "Reds",
                "format": "{:,.0f}",
                "hidden": True,
            }
            table_headers[f"{variant_type}_recall"] = {
                "title": f"{variant_type} Recall",
                "description": f"{variant_type} recall (%)",
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.2f}",
                "max": 100,
            }
            table_headers[f"{variant_type}_precision"] = {
                "title": f"{variant_type} Precision",
                "description": f"{variant_type} precision (%)",
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.2f}",
                "max": 100,
            }
            table_headers[f"{variant_type}_fscore"] = {
                "title": f"{variant_type} F-score",
                "description": f"{variant_type} F-score (%)",
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.2f}",
                "max": 100,
            }
        
        # Only create the table if we have data
        if table_data:
            self.add_section(
                name="Variant Type Performance",
                anchor="wittyer-variant-types",
                description="Event-level metrics by structural variant type. TP/FN/FP columns are hidden by default but can be shown using 'Configure Columns'.",
                plot=table.plot(table_data, table_headers),
            )
