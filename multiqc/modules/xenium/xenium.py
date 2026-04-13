import json
import logging
from pathlib import Path
from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.core import plugin_hooks
from multiqc.plots.table_object import ColumnDict

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Xenium is a spatial transcriptomics platform from 10x Genomics that provides subcellular resolution.

    :::note
    This module provides basic quality metrics from the Xenium pipeline (total transcripts, cells detected,
    transcript assignment rates, and median genes per cell).

    For advanced visualizations including:

    - Transcript quality distributions by gene category
    - Cell and nucleus area distributions
    - Field-of-view quality plots
    - Segmentation method breakdown
    - Transcripts per gene distributions

    Install the [multiqc-xenium-extra](https://pypi.org/project/multiqc-xenium-extra/) plugin:

    ```bash
    pip install multiqc multiqc-xenium-extra
    ```

    The plugin automatically adjusts the log filesize limit to parse large Xenium files (`.parquet` and `.h5`),
    so you don't need to manually configure `log_filesize_limit` in your MultiQC config when using the plugin.
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

        self.data_by_sample = self.ignore_samples(data_by_sample)

        if len(self.data_by_sample) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.data_by_sample)} Xenium reports")

        # Check for QC issues and add warnings
        self.check_qc_warnings()

        # Add software version info from experiment.xenium if available
        for s_name, data in self.data_by_sample.items():
            version = data.get("analysis_sw_version")
            self.add_software_version(version, s_name)

        # Configure initial headers for general stats table
        self.setup_general_stats_headers()

        # Call plugin hook to allow extensions to add data from parquet/H5 files
        # This must happen BEFORE writing data files so plugins can augment data_by_sample
        if "xenium_extra" in plugin_hooks.hook_functions:
            log.debug("Calling xenium_extra plugin hooks")
            for hook_fn in plugin_hooks.hook_functions["xenium_extra"]:
                hook_fn(self)
        else:
            log.info("Run 'pip install multiqc-xenium-extra' for additional visualizations")

        # Write parsed data to a file
        self.write_data_file(self.data_by_sample, "multiqc_xenium")

        # Add key metrics to general stats
        self.xenium_general_stats_table()

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

    def check_qc_warnings(self):
        """Check for common QC issues and add warnings to samples"""
        for s_name, data in self.data_by_sample.items():
            # Check for low transcript assignment rate
            if data.get("fraction_transcripts_assigned", 1.0) < 0.7:
                log.warning(
                    f"Sample '{s_name}' has low transcript assignment rate: {data['fraction_transcripts_assigned']:.3f} (< 0.7). Cell segmentation likely needs refinement."
                )

    def setup_general_stats_headers(self):
        self.genstat_headers = {}

        # Add basic metrics from metrics_summary.csv (always available)
        self.genstat_headers["num_transcripts"] = ColumnDict(
            {
                "title": "Total Transcripts",
                "description": "Total number of transcripts detected",
                "scale": "YlOrRd",
                "format": "{:,.0f}",
                "hidden": False,
            }
        )

        self.genstat_headers["num_cells_detected"] = ColumnDict(
            {
                "title": "Cells",
                "description": "Number of cells detected",
                "scale": "Blues",
                "format": "{:,.0f}",
                "hidden": False,
            }
        )

        self.genstat_headers["fraction_transcripts_assigned"] = ColumnDict(
            {
                "title": "Transcripts Assigned",
                "description": "Fraction of transcripts assigned to cells",
                "suffix": "%",
                "scale": "RdYlGn",
                "modify": lambda x: x * 100.0,
                "max": 100.0,
                "hidden": False,
            }
        )

        self.genstat_headers["median_genes_per_cell"] = ColumnDict(
            {
                "title": "Genes/Cell",
                "description": "Median number of genes per cell",
                "scale": "Purples",
                "format": "{:,.0f}",
                "hidden": False,
            }
        )

        self.genstat_headers["median_transcripts_per_cell"] = ColumnDict(
            {
                "title": "Transcripts/Cell",
                "description": "Median transcripts per cell",
                "min": 0,
                "scale": "Greens",
                "format": "{:,.0f}",
                "hidden": True,
            }
        )

        self.genstat_headers["adjusted_negative_control_probe_rate"] = ColumnDict(
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

    def xenium_general_stats_table(self):
        """Add key Xenium metrics to the general statistics table"""

        self.general_stats_addcols(self.data_by_sample, self.genstat_headers)
