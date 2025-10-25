import json
import logging
from pathlib import Path
from typing import Dict, Tuple

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph
from multiqc.plots.table_object import ColumnDict, TableConfig

log = logging.getLogger(__name__)


# Define gene categories for coloring based on Xenium naming conventions
GENE_CATS = {
    "Pre-designed": {"color": "rgba(31, 119, 180, 0.8)"},  # Standard gene names - blue with transparency
    "Custom": {"color": "rgba(255, 127, 14, 0.8)"},  # Orange with transparency
    "Negative Control Probe": {"color": "rgba(214, 39, 40, 0.8)"},  # Red with transparency
    "Negative Control Codeword": {"color": "rgba(255, 153, 0, 0.8)"},  # Yellow/Orange with transparency
    "Genomic Control Probe": {"color": "rgba(227, 119, 194, 0.8)"},  # Pink with transparency
    "Unassigned Codeword": {"color": "rgba(127, 127, 127, 0.8)"},  # Gray with transparency
    "Deprecated Codeword": {"color": "rgba(188, 189, 34, 0.8)"},  # Olive with transparency
}


def categorize_feature(feature_name) -> Tuple[str, str]:
    """Categorize a feature based on its name
    Splits the feature name into category and feature id"""
    # Check prefixes directly instead of using regex for better performance
    category = ""
    feature_id = feature_name.split("_")[1] if "_" in feature_name else feature_name
    if feature_name.startswith("Custom_"):
        category = "Custom"
    elif feature_name.startswith("NegControlProbe_"):
        category = "Negative Control Probe"
    elif feature_name.startswith("NegControlCodeword_"):
        category = "Negative Control Codeword"
    elif feature_name.startswith("GenomicControlProbe_"):
        category = "Genomic Control Probe"
    elif feature_name.startswith("UnassignedCodeword_"):
        category = "Unassigned Codeword"
    else:
        category = "Pre-designed"  # Default category for standard gene names
    return category, feature_id


class MultiqcModule(BaseMultiqcModule):
    """
    Xenium is a spatial transcriptomics platform from 10x Genomics that provides subcellular resolution.

    :::note
    This module provides basic quality metrics from the Xenium pipeline. For advanced visualizations
    including transcript quality distributions, cell area analysis, and field-of-view quality plots,
    install the `multiqc-xenium-advanced` plugin.
    :::

    The MultiQC module is tested with outputs from xenium-3.x, older versions of xenium output are
    not supported and may even cause MultiQC to crash (see https://github.com/MultiQC/MultiQC/issues/3344).
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Xenium",
            anchor="xenium",
            href="https://www.10xgenomics.com/platforms/xenium",
            info="Spatial transcriptomics platform from 10x Genomics that provides subcellular resolution.",
            # doi=,
        )

        data_by_sample = {}
        for f in self.find_log_files("xenium/metrics"):
            parsed_data = self.parse_xenium_metrics(f)
            if parsed_data:
                # Use parent directory name as sample name
                parent_dir = Path(f["root"]).name if f["root"] else f["s_name"]
                data_by_sample[parent_dir] = parsed_data
                self.add_data_source(f, parent_dir)

        # Parse experiment.xenium files for additional metrics
        for f in self.find_log_files("xenium/experiment"):
            parsed_experiment_data = self.parse_experiment_json(f)
            if parsed_experiment_data:
                # Use parent directory name as sample name
                parent_dir = Path(f["root"]).name if f["root"] else f["s_name"]
                if parent_dir in data_by_sample:
                    # Merge with existing metrics data
                    data_by_sample[parent_dir].update(parsed_experiment_data)
                else:
                    # Create new entry if metrics file wasn't found
                    data_by_sample[parent_dir] = parsed_experiment_data
                self.add_data_source(f, parent_dir)

        data_by_sample = self.ignore_samples(data_by_sample)

        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(data_by_sample)} Xenium reports")

        # Check for QC issues and add warnings
        self.check_qc_warnings(data_by_sample)

        # Add software version info (Xenium files don't contain version info)
        for s_name in data_by_sample.keys():
            self.add_software_version(None, s_name)

        # Call plugin hook to allow extensions to add data from parquet/H5 files
        # This must happen BEFORE writing data files so plugins can augment data_by_sample
        from multiqc.core import plugin_hooks

        # Check if xenium_extend hook is available
        if "xenium_extend" in plugin_hooks.hook_functions:
            log.debug("Calling xenium_extend plugin hooks")
            for hook_fn in plugin_hooks.hook_functions["xenium_extend"]:
                hook_fn(self, data_by_sample)
        else:
            log.info("Run 'pip install multiqc-xenium-extra' for additional visualizations")

        # Write parsed data to a file
        self.write_data_file(data_by_sample, "multiqc_xenium")

        # Add key metrics to general stats
        self.xenium_general_stats_table(data_by_sample)

    def parse_xenium_metrics(self, f) -> Dict:
        """Parse Xenium metrics_summary.csv file"""
        lines = f["f"].splitlines()
        if len(lines) < 2:
            return {}

        # Get header and data row
        header = lines[0].split(",")
        data_row = lines[1].split(",")

        if len(header) != len(data_row):
            log.warning(f"Header and data row lengths don't match in {f['fn']}")
            return {}

        # Create mapping of headers to values
        metrics = dict(zip(header, data_row))

        # Convert numeric fields
        numeric_fields = [
            "region_area",
            "total_cell_area",
            "total_high_quality_decoded_transcripts",
            "fraction_transcripts_decoded_q20",
            "fraction_predesigned_transcripts_decoded_q20",
            "fraction_custom_transcripts_decoded_q20",
            "nuclear_transcripts_per_100um2",
            "decoded_transcripts_per_100um2",
            "adjusted_negative_control_probe_rate",
            "adjusted_negative_control_codeword_rate",
            "adjusted_genomic_control_probe_rate",
            "negative_control_probe_counts_per_control_per_cell",
            "genomic_control_probe_counts_per_control_per_cell",
            "estimated_number_of_false_positive_transcripts_per_cell",
            "estimated_number_of_false_positive_transcripts_per_cell_including_genomic_counts",
            "num_cells_detected",
            "fraction_empty_cells",
            "cells_per_100um2",
            "fraction_transcripts_assigned",
            "median_genes_per_cell",
            "median_predesigned_genes_per_cell",
            "median_custom_genes_per_cell",
            "median_transcripts_per_cell",
            "median_predesigned_transcripts_per_cell",
            "median_custom_transcripts_per_cell",
            "thickness_transcripts_high_quality",
            "segmented_cell_stain_frac",
            "segmented_cell_boundary_frac",
            "segmented_cell_interior_frac",
            "segmented_cell_nuc_expansion_frac",
            "segmented_cell_boundary_count",
            "segmented_cell_interior_count",
            "segmented_cell_nuc_expansion_count",
            "fraction_of_ambiguous_nucleus_mask_pixels",
            "fraction_of_ambiguous_cell_mask_pixels",
            "fraction_of_nucleus_polygons_removed",
            "fraction_of_cell_polygons_removed",
            "fraction_of_nuclei_without_cell",
            "number_of_cell_non_simple_polygons",
            "number_of_cell_multi_polygons",
            "number_of_nucleus_non_simple_polygons",
            "number_of_nucleus_multi_polygons",
            "segmented_cell_imported_frac",
            "segmented_cell_imported_count",
        ]

        parsed_metrics = {}
        for field in numeric_fields:
            if field in metrics and metrics[field] != "":
                try:
                    parsed_metrics[field] = float(metrics[field])
                except ValueError:
                    log.warning(f"Could not convert {field}='{metrics[field]}' to float")

        # Keep string fields
        string_fields = [
            "run_name",
            "cassette_name",
            "region_name",
            "panel_name",
            "panel_design_id",
            "predesigned_panel_id",
            "stain_definition",
        ]
        for field in string_fields:
            if field in metrics:
                parsed_metrics[field] = metrics[field]

        return parsed_metrics

    def parse_experiment_json(self, f) -> Dict:
        """Parse Xenium experiment.xenium JSON file for additional metrics"""
        try:
            experiment_data = json.loads(f["f"])

            # Extract key metrics that aren't in the CSV
            parsed_experiment = {}

            # Total transcript count - this is the main missing metric from notebooks
            if "num_transcripts" in experiment_data:
                parsed_experiment["num_transcripts"] = experiment_data["num_transcripts"]

            # High quality transcript count
            if "num_transcripts_high_quality" in experiment_data:
                parsed_experiment["num_transcripts_high_quality"] = experiment_data["num_transcripts_high_quality"]

            # Software version info
            if "analysis_sw_version" in experiment_data:
                parsed_experiment["analysis_sw_version"] = experiment_data["analysis_sw_version"]

            # Panel information
            if "panel_num_targets_predesigned" in experiment_data:
                parsed_experiment["panel_num_targets_predesigned"] = experiment_data["panel_num_targets_predesigned"]
            if "panel_num_targets_custom" in experiment_data:
                parsed_experiment["panel_num_targets_custom"] = experiment_data["panel_num_targets_custom"]

            return parsed_experiment

        except (json.JSONDecodeError, KeyError) as e:
            log.warning(f"Could not parse experiment.xenium file {f['fn']}: {e}")
            return {}

    def check_qc_warnings(self, data_by_sample):
        """Check for common QC issues and add warnings to samples"""
        for s_name, data in data_by_sample.items():
            # Check for low transcript assignment rate
            if data.get("fraction_transcripts_assigned", 1.0) < 0.7:
                log.warning(
                    f"Sample '{s_name}' has low transcript assignment rate: {data['fraction_transcripts_assigned']:.3f} (< 0.7). Cell segmentation likely needs refinement."
                )

    def xenium_general_stats_table(self, data_by_sample):
        """Add key Xenium metrics to the general statistics table"""
        headers = {}

        # Columns are PREPENDED (added to the left), so add in REVERSE order
        # Target order: Total Transcripts, Cells, Transcripts Assigned, Genes/Cell, Q20+ Transcripts, Median Cell, Median Nucleus, Nucleus/Cell
        # Add in reverse: Nucleus/Cell, Median Nucleus, Median Cell, Q20+ Transcripts, Genes/Cell, Transcripts Assigned, Cells, Total Transcripts

        headers["median_transcripts_per_cell"] = ColumnDict(
            {
                "title": "Transcripts/Cell",
                "description": "Median transcripts per cell",
                "min": 0,
                "scale": "Greens",
                "format": "{:,.0f}",
                "hidden": True,
            }
        )

        headers["adjusted_negative_control_probe_rate"] = ColumnDict(
            {
                "title": "Neg Ctrl Rate",
                "description": "Adjusted negative control probe rate",
                "max": 0.1,
                "min": 0,
                "scale": "OrRd",
                "format": "{:,.3f}",
                "hidden": True,
            }
        )

        self.general_stats_addcols(data_by_sample, headers)

    def xenium_segmentation_plot(self, data_by_sample):
        """Create stacked bar chart showing cell segmentation methods"""
        keys = {
            "segmented_cell_boundary_frac": {"name": "Boundary", "color": "#c72eba"},
            "segmented_cell_interior_frac": {"name": "Interior", "color": "#bbbf34"},
            "segmented_cell_nuc_expansion_frac": {"name": "Nuclear Expansion", "color": "#426cf5"},
        }

        plot_config = {
            "id": "xenium_segmentation",
            "title": "Xenium: Cell Segmentation Method",
            "ylab": "Fraction",
            "stacking": "normal",
            "ymax": 1.0,
            "cpswitch": False,
        }

        return bargraph.plot(data_by_sample, keys, plot_config)
