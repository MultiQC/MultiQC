import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Optional

import polars as pl

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph, scatter

log = logging.getLogger(__name__)


# Define gene categories for coloring based on Xenium naming conventions
GENE_CATS = {
    "pre-designed": {"color": "#1f77b4"},  # Standard gene names
    "custom": {"color": "#ff7f0e"},
    "negative_control_probe": {"color": "#d62728"},
    "negative_control_codeword": {"color": "#ff9900"},
    "genomic_control": {"color": "#e377c2"},
    "unassigned": {"color": "#7f7f7f"},
}


def categorize_feature(feature_name):
    """Categorize a feature based on its name"""
    # Check prefixes directly instead of using regex for better performance
    if feature_name.startswith("Custom_"):
        return "custom"
    elif feature_name.startswith("NegControlProbe_"):
        return "negative_control_probe"
    elif feature_name.startswith("NegControlCodeword_"):
        return "negative_control_codeword"
    elif feature_name.startswith("GenomicControlProbe_"):
        return "genomic_control"
    elif feature_name.startswith("UnassignedCodeword_"):
        return "unassigned"
    else:
        return "pre-designed"  # Default category for standard gene names


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
        for f in self.find_log_files("xenium/metrics"):
            parsed_data = self.parse_xenium_metrics(f)
            if parsed_data:
                # Use parent directory name as sample name
                parent_dir = Path(f["root"]).name if f["root"] else f["s_name"]
                data_by_sample[parent_dir] = parsed_data
                self.add_data_source(f, parent_dir)

        # Parse transcript quality data
        transcript_data_by_sample = {}
        for transcript_f in self.find_log_files("xenium/transcripts", filecontents=False, filehandles=False):
            parsed_transcript_data = self.parse_transcripts_parquet(transcript_f)
            if parsed_transcript_data:
                # Use parent directory name as sample name
                parent_dir = Path(transcript_f["root"]).name if transcript_f["root"] else transcript_f["s_name"]
                transcript_data_by_sample[parent_dir] = parsed_transcript_data
                self.add_data_source(transcript_f, parent_dir)

        data_by_sample = self.ignore_samples(data_by_sample)

        if len(data_by_sample) == 0 and len(transcript_data_by_sample) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(data_by_sample)} Xenium reports")

        # Add software version info (Xenium files don't contain version info)
        for s_name in data_by_sample.keys():
            self.add_software_version(None, s_name)

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

        # Add transcript quality section if transcript data is available
        if transcript_data_by_sample:
            self.add_section(
                name="Transcript Quality Distribution",
                anchor="xenium-transcript-quality",
                description="Distribution of transcript quality values by codeword category across samples",
                plot=self.xenium_transcript_quality_plot(transcript_data_by_sample),
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

    def parse_transcripts_parquet(self, f) -> Optional[Dict]:
        """Parse Xenium transcripts.parquet file to extract quality distribution by codeword"""
        # Read the parquet file content - sample more for better scatter plots
        file_path = Path(f["root"]) / f["fn"]
        df = pl.read_parquet(file_path)

        # Check if required columns exist
        required_cols = ["qv", "feature_name"]
        if not all(col in df.columns for col in required_cols):
            log.warning(f"Missing required columns in {f['fn']}: {required_cols}")
            return None

        # Group by feature_name and calculate both quality distribution and transcript counts
        quality_dist = {}
        transcript_counts = {}

        grouped = df.group_by("feature_name").agg(
            [
                pl.col("qv").value_counts().alias("qv_counts"),
                pl.col("qv").mean().alias("mean_qv"),
                pl.len().alias("transcript_count"),
            ]
        )

        for row in grouped.iter_rows(named=True):
            feature = str(row["feature_name"])
            qv_counts_df = row["qv_counts"]
            transcript_counts[feature] = {"count": row["transcript_count"], "mean_quality": row["mean_qv"]}

            # Convert to dictionary format for backward compatibility
            qv_dict = {}
            for qv_row in qv_counts_df:
                qv_dict[qv_row["qv"]] = qv_row["count"]
            quality_dist[feature] = qv_dict

        return {
            "quality_distribution": quality_dist,
            "transcript_counts": transcript_counts,
            "total_transcripts": df.height,
            "unique_features": len(quality_dist),
        }

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

    def xenium_transcript_quality_plot(self, transcript_data_by_sample):
        """Create adaptive transcript quality plots based on sample count"""
        if not transcript_data_by_sample:
            return None

        num_samples = len(transcript_data_by_sample)

        if num_samples == 1:
            # Single sample: scatter plot of transcript count vs mean quality
            return self._create_single_sample_scatter(transcript_data_by_sample)

        else:
            # Many samples: violin plots with tabs for categories
            return self._create_multi_sample(transcript_data_by_sample)

    def _create_single_sample_scatter(self, transcript_data_by_sample):
        """Create scatter plot - handles both single and multiple samples"""
        # Prepare scatter data - create individual points for each gene from all samples
        plot_data: Dict[str, Any] = {}

        for sample_name, sample_data in transcript_data_by_sample.items():
            if "transcript_counts" not in sample_data:
                continue

            for feature, counts_data in sample_data["transcript_counts"].items():
                category = categorize_feature(feature)

                if category not in plot_data:
                    plot_data[category] = []

                # Each point is a separate data point
                # For multiple samples, include sample name in the hover text
                if len(transcript_data_by_sample) > 1:
                    point_name = f"{feature} ({sample_name})"
                else:
                    point_name = feature

                plot_data[category].append(
                    {
                        "x": counts_data["count"],
                        "y": counts_data["mean_quality"],
                        "name": point_name,  # Use gene name (+ sample) for hover text
                    }
                )

        # Filter out empty categories and add colors to each point
        final_plot_data = {}
        for category, points in plot_data.items():
            if points:  # Only include categories with data
                # Add color to each point in the category
                for point in points:
                    point["color"] = GENE_CATS[category]["color"]
                final_plot_data[category] = points

        # Adjust title based on number of samples
        if len(transcript_data_by_sample) == 1:
            title = "Xenium: Gene-Specific Transcript Quality"
        else:
            title = f"Xenium: Gene-Specific Transcript Quality ({len(transcript_data_by_sample)} samples)"

        config = {
            "id": "xenium_transcript_quality_combined",
            "title": title,
            "xlab": "Total transcripts per gene",
            "ylab": "Mean calibrated quality of gene transcripts",
            "xlog": True,
            "marker_size": 5,
            "series_label": "transcripts",
        }

        return scatter.plot(final_plot_data, config)

    def _create_multi_sample(self, transcript_data_by_sample):
        """Create multi-dataset line plot with categories as datasets"""
        # Prepare data for each category as a separate dataset
        # First, collect all categories that have data across samples
        all_categories = set()
        for sample_data in transcript_data_by_sample.values():
            if "transcript_counts" in sample_data:
                for feature in sample_data["transcript_counts"].keys():
                    category = categorize_feature(feature)
                    all_categories.add(category)

        # Create a dataset for "all transcripts" first (combining all categories)
        datasets = {"All transcripts": {}}
        for cat in all_categories:
            datasets[cat] = {}

        for sample_name, sample_data in transcript_data_by_sample.items():
            if "transcript_counts" not in sample_data:
                continue

            datasets["All transcripts"][sample_name] = {}
            for cat in all_categories:
                datasets[cat][sample_name] = {}

            for feature, counts_data in sample_data["transcript_counts"].items():
                count = counts_data["count"]
                quality = counts_data["mean_quality"]
                datasets["All transcripts"][sample_name][count] = quality
                for cat in all_categories:
                    if categorize_feature(feature) == cat:
                        datasets[cat][sample_name][count] = quality

        if not datasets:
            return None

        config = {
            "id": "xenium_transcript_quality_multi",
            "title": "Xenium: Transcript Quality by Category",
            "xlab": "Total transcripts per gene",
            "ylab": "Mean calibrated quality of gene transcripts",
            "xlog": True,
            "data_labels": [{"name": cat} for cat in datasets.keys()],
        }

        return linegraph.plot(list(datasets.values()), config)
