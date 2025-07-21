import logging
from typing import Dict, Any

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Xenium",
            anchor="xenium",
            href="https://www.10xgenomics.com/products/xenium-in-situ",
            info="is a spatial transcriptomics platform from 10x Genomics that provides subcellular resolution.",
            doi="10.1038/s41587-022-01483-z",
        )

        data_by_sample = {}
        for f in self.find_log_files("xenium"):
            parsed_data = self.parse_xenium_metrics(f)
            if parsed_data:
                data_by_sample[f["s_name"]] = parsed_data
                self.add_data_source(f, f["s_name"])

        data_by_sample = self.ignore_samples(data_by_sample)

        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(data_by_sample)} Xenium reports")

        # Write parsed data to a file
        self.write_data_file(data_by_sample, "multiqc_xenium")

        # Add key metrics to general stats
        self.xenium_general_stats_table(data_by_sample)

        # Create plots
        self.add_section(
            name="Cell Detection",
            anchor="xenium-cells",
            description="Number of cells detected and cells per 100μm²",
            plot=self.xenium_cells_plot(data_by_sample),
        )

        self.add_section(
            name="Transcript Assignment",
            anchor="xenium-transcripts",
            description="Fraction of transcripts assigned to cells and transcripts per cell",
            plot=self.xenium_transcripts_plot(data_by_sample),
        )

        self.add_section(
            name="Segmentation Methods",
            anchor="xenium-segmentation",
            description="Distribution of cell segmentation methods used",
            plot=self.xenium_segmentation_plot(data_by_sample),
        )

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

    def xenium_general_stats_table(self, data_by_sample):
        """Add key Xenium metrics to the general statistics table"""
        headers: Dict[str, Dict[str, Any]] = {
            "num_cells_detected": {
                "title": "Cells",
                "description": "Number of cells detected",
                "scale": "Blues",
                "format": "{:,.0f}",
            },
            "fraction_transcripts_assigned": {
                "title": "% Transcripts Assigned",
                "description": "Fraction of transcripts assigned to cells",
                "suffix": "%",
                "scale": "RdYlGn",
                "modify": lambda x: x * 100.0,
            },
            "median_genes_per_cell": {
                "title": "Genes/Cell",
                "description": "Median number of genes per cell",
                "scale": "Purples",
                "format": "{:,.0f}",
            },
            "fraction_transcripts_decoded_q20": {
                "title": "% Q20+ Transcripts",
                "description": "Fraction of transcripts decoded with Q20+",
                "suffix": "%",
                "scale": "Greens",
                "modify": lambda x: x * 100.0,
            },
        }
        self.general_stats_addcols(data_by_sample, headers)

    def xenium_cells_plot(self, data_by_sample):
        """Create bar plot for cell detection metrics"""
        plot_data = {}
        for s_name, data in data_by_sample.items():
            plot_data[s_name] = {
                "num_cells_detected": data.get("num_cells_detected", 0),
                "cells_per_100um2": data.get("cells_per_100um2", 0),
            }

        keys = {
            "num_cells_detected": {"name": "Total Cells", "color": "#1f77b4"},
            "cells_per_100um2": {"name": "Cells per 100μm²", "color": "#ff7f0e"},
        }

        config = {
            "id": "xenium_cells",
            "title": "Xenium: Cell Detection",
            "ylab": "Count / Density",
            "cpswitch_counts_label": "Counts",
        }

        return bargraph.plot(plot_data, keys, config)

    def xenium_transcripts_plot(self, data_by_sample):
        """Create bar plot for transcript metrics"""
        plot_data = {}
        for s_name, data in data_by_sample.items():
            plot_data[s_name] = {
                "fraction_transcripts_assigned": data.get("fraction_transcripts_assigned", 0) * 100,
                "median_transcripts_per_cell": data.get("median_transcripts_per_cell", 0),
            }

        keys = {
            "fraction_transcripts_assigned": {"name": "% Transcripts Assigned", "color": "#2ca02c"},
            "median_transcripts_per_cell": {"name": "Median Transcripts/Cell", "color": "#d62728"},
        }

        config = {
            "id": "xenium_transcripts",
            "title": "Xenium: Transcript Assignment",
            "ylab": "Percentage / Count",
            "cpswitch_counts_label": "Values",
        }

        return bargraph.plot(plot_data, keys, config)

    def xenium_segmentation_plot(self, data_by_sample):
        """Create stacked bar plot for segmentation methods"""
        plot_data = {}
        for s_name, data in data_by_sample.items():
            plot_data[s_name] = {
                "segmented_cell_boundary_frac": data.get("segmented_cell_boundary_frac", 0),
                "segmented_cell_interior_frac": data.get("segmented_cell_interior_frac", 0),
                "segmented_cell_nuc_expansion_frac": data.get("segmented_cell_nuc_expansion_frac", 0),
            }

        keys = {
            "segmented_cell_boundary_frac": {"name": "Boundary", "color": "#1f77b4"},
            "segmented_cell_interior_frac": {"name": "Interior", "color": "#ff7f0e"},
            "segmented_cell_nuc_expansion_frac": {"name": "Nuclear Expansion", "color": "#2ca02c"},
        }

        config = {
            "id": "xenium_segmentation",
            "title": "Xenium: Cell Segmentation Methods",
            "ylab": "Fraction",
            "stacking": "normal",
            "ymax": 1.0,
            "cpswitch": False,
        }

        return bargraph.plot(plot_data, keys, config)
