import json
import logging
from pathlib import Path
from typing import Any, Dict, Optional

import polars as pl

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph, scatter, box

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

        # Parse transcript quality data
        transcript_data_by_sample = {}
        for transcript_f in self.find_log_files("xenium/transcripts", filecontents=False, filehandles=False):
            parsed_transcript_data = self.parse_transcripts_parquet(transcript_f)
            if parsed_transcript_data:
                # Use parent directory name as sample name
                parent_dir = Path(transcript_f["root"]).name if transcript_f["root"] else transcript_f["s_name"]
                transcript_data_by_sample[parent_dir] = parsed_transcript_data
                self.add_data_source(transcript_f, parent_dir)

        # Parse cells.parquet files for cell-level metrics
        cells_data_by_sample = {}
        for cells_f in self.find_log_files("xenium/cells", filecontents=False, filehandles=False):
            parsed_cells_data = self.parse_cells_parquet(cells_f)
            if parsed_cells_data:
                # Use parent directory name as sample name
                parent_dir = Path(cells_f["root"]).name if cells_f["root"] else cells_f["s_name"]
                cells_data_by_sample[parent_dir] = parsed_cells_data
                self.add_data_source(cells_f, parent_dir)

        data_by_sample = self.ignore_samples(data_by_sample)

        if len(data_by_sample) == 0 and len(transcript_data_by_sample) == 0 and len(cells_data_by_sample) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(data_by_sample)} Xenium reports")

        # Check for QC issues and add warnings
        self.check_qc_warnings(data_by_sample)

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
            helptext="""
            This plot shows fundamental cell detection metrics from Xenium spatial transcriptomics:
            
            * **Total Cells**: Total number of cells detected via segmentation algorithms
            * **Cells per 100μm²**: Cell density normalized by tissue area
            
            **What to look for:**
            * Consistent cell counts across similar tissue types
            * Cell densities appropriate for the tissue type (e.g., brain ~1000-2000 cells/100μm², liver ~500-1000)
            * Large variations might indicate segmentation issues or genuine biological differences
            
            **Troubleshooting:**
            * Very low cell counts: Check segmentation parameters, tissue quality, or imaging issues
            * Very high cell counts: May indicate over-segmentation or dense cell types
            """,
            plot=self.xenium_cells_plot(data_by_sample),
        )

        self.add_section(
            name="Segmentation Methods",
            anchor="xenium-segmentation",
            description="Distribution of cell segmentation methods used",
            helptext="""
            This stacked bar chart shows the fraction of cells segmented by each method:
            
            * **Boundary**: Cells segmented using boundary staining (e.g., ATP1A1/E-cadherin/CD45)
            * **Interior**: Cells segmented using interior staining (e.g., 18S RNA)  
            * **Nuclear Expansion**: Cells segmented by expanding from nucleus boundaries
            
            **What to look for:**
            * **Boundary segmentation** typically provides the most accurate cell boundaries
            * **High nuclear expansion fraction** may indicate poor membrane staining
            * Consistent ratios across samples of the same tissue type
            
            **Interpretation:**
            * >80% boundary segmentation: Excellent membrane staining and segmentation
            * >50% nuclear expansion: Consider optimizing membrane staining protocols
            * Large sample-to-sample variation: Check staining consistency
            """,
            plot=self.xenium_segmentation_plot(data_by_sample),
        )

        # Add transcript quality section if transcript data is available
        if transcript_data_by_sample:
            self.add_section(
                name="Transcript Quality Distribution",
                anchor="xenium-transcript-quality",
                description="Distribution of transcript quality values by codeword category across samples",
                helptext="""
                This plot shows transcript quality (QV score) vs. transcript count by gene category:
                
                * **Pre-designed genes**: Standard genes from Xenium panels (blue)
                * **Custom genes**: User-added custom targets (orange)  
                * **Negative controls**: Control probes for background estimation (red/yellow)
                * **Genomic controls**: Genomic DNA controls (pink)
                
                **Quality Score (QV) Interpretation:**
                * QV ≥20: High-quality transcripts (≥99% accuracy)
                * QV 10-20: Medium quality (90-99% accuracy)
                * QV <10: Low-quality transcripts (<90% accuracy)
                
                **What to look for:**
                * Pre-designed genes should cluster at high QV (>20) and reasonable counts
                * Negative controls should have low counts and variable quality
                * Outlier genes with very low quality or unexpectedly high/low counts
                
                **Single vs. Multiple samples:**
                * Single sample: Scatter plot showing individual genes
                * Multiple samples: Line plot showing category trends
                """,
                plot=self.xenium_transcript_quality_plot(transcript_data_by_sample),
            )

            # Add Field of View quality section if FoV data is available
            fov_plot = self.xenium_fov_quality_plot(transcript_data_by_sample)
            if fov_plot is not None:
                self.add_section(
                    name="Field of View Quality",
                    anchor="xenium-fov-quality",
                    description="Transcript quality distribution by Field of View (FoV)",
                    helptext="""
                    This plot shows transcript quality distributions across different imaging fields (FoVs - Fields of View):
                    
                    **What is a Field of View?**
                    * Each FoV represents one microscope imaging area/tile
                    * Large tissue sections are imaged as multiple overlapping FoVs
                    * FoVs are systematically captured in a grid pattern across the tissue
                    
                    **Plot types:**
                    * **Single sample**: Box plots showing quality distributions per FoV (median, quartiles, outliers)
                    * **Multiple samples**: Box plots with aggregated quality distributions per FoV across all samples
                    
                    **Box plot interpretation:**
                    * **Box boundaries**: 25th and 75th percentiles (Q1 and Q3)
                    * **Center line**: Median quality score
                    * **Whiskers**: Extend to 1.5 × IQR or the most extreme data point
                    * **Points**: Individual transcript quality scores (outliers or small datasets)
                    
                    **What to look for:**
                    * **Consistent quality** across FoVs (similar median QV values around 30-40)
                    * **Tight distributions** (narrow boxes indicate consistent quality within FoVs)
                    * **No systematic patterns**: Random variation is normal, systematic gradients are not
                    * **Outlier FoVs**: Any FoV with notably poor median quality (<20 QV)
                    
                    **Quality thresholds:**
                    * QV >30: Excellent imaging quality
                    * QV 20-30: Good quality  
                    * QV <20: Poor quality, investigate issues
                    
                    **Troubleshooting:**
                    * Specific low-quality FoVs: Focus/illumination issues, debris, tissue damage
                    * Edge effects: FoVs at tissue edges often have lower quality
                    * Systematic gradients: Temperature, timing, or optical alignment issues
                    """,
                    plot=fov_plot,
                )

        # Add cell metrics section if cells data is available
        if cells_data_by_sample:
            self.add_section(
                name="Cell Area Distribution",
                anchor="xenium-cell-area",
                description="Distribution of cell and nucleus areas",
                helptext="""
                This plot shows cell morphology metrics derived from segmentation:
                
                **Metrics displayed:**
                * **Mean Cell Area**: Average cell size in μm² 
                * **Mean Nucleus Area**: Average nucleus size in μm²
                * **Nucleus/Cell Area Ratio**: Proportion of cell occupied by nucleus
                
                **Typical ranges (tissue-dependent):**
                * **Cell area**: 50-200 μm² for most tissues
                * **Nucleus area**: 20-80 μm² 
                * **Nucleus/cell ratio**: 0.1-0.4 (10-40%)
                
                **What to look for:**
                * **Consistent metrics** across samples of same tissue type
                * **Biologically reasonable values** for your tissue
                * **Nucleus/cell ratio** should be <0.5 (nucleus shouldn't dominate cell)
                
                **Troubleshooting:**
                * Very large cells: Over-segmentation issues, cell doublets
                * Very small cells: Under-segmentation, debris detection
                * High nucleus/cell ratio: Nuclear expansion artifacts, poor membrane staining
                * Missing metrics: Insufficient cells for statistical calculation
                """,
                plot=self.xenium_cell_area_plot(cells_data_by_sample),
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

        result = {
            "quality_distribution": quality_dist,
            "transcript_counts": transcript_counts,
            "total_transcripts": df.height,
            "unique_features": len(quality_dist),
        }

        # Add FoV quality analysis if fov_name column is present
        if "fov_name" in df.columns:
            fov_quality_stats = {}
            fov_quality_distributions = {}

            # Group by FoV and calculate quality stats and distributions
            fov_grouped = df.group_by("fov_name").agg(
                [
                    pl.col("qv").mean().alias("mean_qv"),
                    pl.col("qv").median().alias("median_qv"),
                    pl.col("qv").std().alias("std_qv"),
                    pl.col("qv").alias("qv_values"),  # Keep all QV values for distributions
                    pl.len().alias("transcript_count"),
                ]
            )

            for row in fov_grouped.iter_rows(named=True):
                fov_name = str(row["fov_name"])
                fov_quality_stats[fov_name] = {
                    "mean_quality": row["mean_qv"],
                    "median_quality": row["median_qv"],
                    "std_quality": row["std_qv"],
                    "transcript_count": row["transcript_count"],
                }

                # Store quality values for violin plots (limit to reasonable sample size)
                qv_values = row["qv_values"]

                # Check if it's a polars Series or already a list
                if hasattr(qv_values, "to_list"):
                    # It's a polars Series
                    if len(qv_values) > 1000:
                        sampled_qv = qv_values.sample(1000, seed=42)
                        fov_quality_distributions[fov_name] = sampled_qv.to_list()
                    else:
                        fov_quality_distributions[fov_name] = qv_values.to_list()
                else:
                    # It's already a Python list
                    if len(qv_values) > 1000:
                        import random

                        random.seed(42)
                        sampled_qv = random.sample(qv_values, 1000)
                        fov_quality_distributions[fov_name] = sampled_qv
                    else:
                        fov_quality_distributions[fov_name] = qv_values

            result["fov_quality_stats"] = fov_quality_stats
            result["fov_quality_distributions"] = fov_quality_distributions

        return result

    def parse_cells_parquet(self, f) -> Optional[Dict]:
        """Parse Xenium cells.parquet file to extract cell-level metrics"""
        file_path = Path(f["root"]) / f["fn"]

        try:
            # Read cells parquet file
            df = pl.read_parquet(file_path)

            # Check for required columns
            required_cols = ["cell_area", "nucleus_area", "total_counts"]
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                log.warning(f"Missing columns in {f['fn']}: {missing_cols}")
                return None

            # Calculate summary statistics
            cell_stats = {}

            # Cell area distribution stats
            cell_area_stats = df["cell_area"].drop_nulls()
            if cell_area_stats.len() > 0:
                cell_stats.update(
                    {
                        "cell_area_mean": cell_area_stats.mean(),
                        "cell_area_median": cell_area_stats.median(),
                        "cell_area_std": cell_area_stats.std(),
                        "cell_area_min": cell_area_stats.min(),
                        "cell_area_max": cell_area_stats.max(),
                    }
                )

            # Nucleus area distribution stats
            nucleus_area_stats = df["nucleus_area"].drop_nulls()
            if nucleus_area_stats.len() > 0:
                cell_stats.update(
                    {
                        "nucleus_area_mean": nucleus_area_stats.mean(),
                        "nucleus_area_median": nucleus_area_stats.median(),
                        "nucleus_area_std": nucleus_area_stats.std(),
                    }
                )

                # Nucleus to cell area ratio (only for non-null values)
                valid_ratio_df = df.filter(
                    (pl.col("cell_area").is_not_null())
                    & (pl.col("nucleus_area").is_not_null())
                    & (pl.col("cell_area") > 0)
                )
                if valid_ratio_df.height > 0:
                    ratio = valid_ratio_df["nucleus_area"] / valid_ratio_df["cell_area"]
                    cell_stats.update(
                        {
                            "nucleus_to_cell_area_ratio_mean": ratio.mean(),
                            "nucleus_to_cell_area_ratio_median": ratio.median(),
                        }
                    )

            # Total cell count
            cell_stats["total_cells"] = df.height

            # Add nucleus RNA fraction if nucleus_count is available
            if "nucleus_count" in df.columns:
                # Filter out cells with zero total counts to avoid division by zero
                valid_cells = df.filter(pl.col("total_counts") > 0)
                if valid_cells.height > 0:
                    nucleus_rna_fraction = valid_cells["nucleus_count"] / valid_cells["total_counts"]
                    cell_stats.update(
                        {
                            "nucleus_rna_fraction_mean": nucleus_rna_fraction.mean(),
                            "nucleus_rna_fraction_median": nucleus_rna_fraction.median(),
                        }
                    )

            return cell_stats

        except Exception as e:
            log.warning(f"Could not parse cells.parquet file {f['fn']}: {e}")
            return None

    def check_qc_warnings(self, data_by_sample):
        """Check for quality control issues and log warnings"""
        low_assignment_threshold = 0.7  # 70% threshold as mentioned in notebooks

        for s_name, data in data_by_sample.items():
            if "fraction_transcripts_assigned" in data:
                assignment_rate = data["fraction_transcripts_assigned"]
                if assignment_rate < low_assignment_threshold:
                    log.warning(
                        f"Sample '{s_name}' has low transcript assignment rate: "
                        f"{assignment_rate:.3f} (< {low_assignment_threshold}). "
                        f"Cell segmentation likely needs refinement."
                    )

    def xenium_general_stats_table(self, data_by_sample):
        """Add key Xenium metrics to the general statistics table"""
        headers: Dict[str, Dict[str, Any]] = {
            "num_transcripts": {
                "title": "Total Transcripts",
                "description": "Total number of transcripts detected",
                "scale": "YlOrRd",
                "format": "{:,.0f}",
            },
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
        datasets: Dict[str, Dict[str, Dict]] = {"All transcripts": {}}
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

    def xenium_cell_area_plot(self, cells_data_by_sample):
        """Create bar plot for cell area metrics"""
        plot_data = {}
        for s_name, data in cells_data_by_sample.items():
            plot_data[s_name] = {}

            # Only add metrics that exist and are not None/NaN
            if (
                "cell_area_mean" in data
                and data["cell_area_mean"] is not None
                and str(data["cell_area_mean"]).lower() != "nan"
            ):
                try:
                    plot_data[s_name]["cell_area_mean"] = float(data["cell_area_mean"])
                except (ValueError, TypeError):
                    pass

            if (
                "nucleus_area_mean" in data
                and data["nucleus_area_mean"] is not None
                and str(data["nucleus_area_mean"]).lower() != "nan"
            ):
                try:
                    plot_data[s_name]["nucleus_area_mean"] = float(data["nucleus_area_mean"])
                except (ValueError, TypeError):
                    pass

            if (
                "nucleus_to_cell_area_ratio_mean" in data
                and data["nucleus_to_cell_area_ratio_mean"] is not None
                and str(data["nucleus_to_cell_area_ratio_mean"]).lower() != "nan"
            ):
                try:
                    plot_data[s_name]["nucleus_to_cell_ratio"] = float(data["nucleus_to_cell_area_ratio_mean"])
                except (ValueError, TypeError):
                    pass

        # Check if we have any data to plot
        has_data = any(bool(sample_data) for sample_data in plot_data.values())
        if not has_data:
            return None

        keys = {
            "cell_area_mean": {"name": "Mean Cell Area (μm²)", "color": "#1f77b4"},
            "nucleus_area_mean": {"name": "Mean Nucleus Area (μm²)", "color": "#ff7f0e"},
            "nucleus_to_cell_ratio": {"name": "Nucleus/Cell Area Ratio", "color": "#2ca02c"},
        }

        config = {
            "id": "xenium_cell_area",
            "title": "Xenium: Cell Area Metrics",
            "ylab": "Area (μm²) / Ratio",
            "cpswitch_counts_label": "Values",
        }

        return bargraph.plot(plot_data, keys, config)

    def xenium_fov_quality_plot(self, transcript_data_by_sample):
        """Create adaptive FoV quality plot - violin plots for single sample, summary for multiple"""
        fov_data_found = False
        samples_with_fov = []

        # Check which samples have FoV data
        for s_name, data in transcript_data_by_sample.items():
            if "fov_quality_distributions" in data and "fov_quality_stats" in data:
                fov_data_found = True
                samples_with_fov.append(s_name)

        if not fov_data_found:
            return None

        num_samples = len(samples_with_fov)

        if num_samples == 1:
            # Single sample: Create violin plot showing quality distributions per FoV
            return self._create_single_sample_fov_box(transcript_data_by_sample[samples_with_fov[0]])
        else:
            # Multiple samples: Create bar plot showing mean quality per FoV across samples
            return self._create_multi_sample_fov_summary(transcript_data_by_sample, samples_with_fov)

    def _create_single_sample_fov_box(self, sample_data):
        """Create box plot showing quality distributions for single sample FoVs"""
        if "fov_quality_distributions" not in sample_data:
            return None

        plot_data = {}

        # Use the raw quality distributions for proper box plots
        for fov_name, qv_values in sample_data["fov_quality_distributions"].items():
            if qv_values and len(qv_values) > 0:
                # Box plot expects the raw data points
                plot_data[fov_name] = qv_values

        if not plot_data:
            return None

        config = {
            "id": "xenium_fov_quality_single",
            "title": "Xenium: transcript quality distribution by field of view",
            "xlab": "Field of view",
            "series_label": "fields of view",
            "ylab": "Quality value (QV)",
            "sort_by_median": True,  # Use the new core box plot sorting feature
            "sort_switch_sorted_active": True,  # Start with sorted view active
        }

        return box.plot(plot_data, config)

    def _create_multi_sample_fov_summary(self, transcript_data_by_sample, samples_with_fov):
        """Create box plot showing quality distributions for each FoV aggregated across all samples"""
        fov_quality_data = {}

        # Aggregate quality distributions for each FoV across all samples
        for s_name in samples_with_fov:
            data = transcript_data_by_sample[s_name]
            if "fov_quality_distributions" in data:
                fov_distributions = data["fov_quality_distributions"]
                for fov_name, quality_values in fov_distributions.items():
                    if fov_name not in fov_quality_data:
                        fov_quality_data[fov_name] = []
                    # Add all quality values from this sample's FoV to the aggregated distribution
                    fov_quality_data[fov_name].extend(quality_values)

        if not fov_quality_data:
            return None

        config = {
            "id": "xenium_fov_quality_multi",
            "title": "Xenium: Transcript quality distribution by field of view (averaged across samples)",
            "xlab": "Quality Score (QV)",
            "ylab": "Field of View",
            "series_label": "fields of view",
            "sort_by_median": True,  # Use the new core box plot sorting feature
            "sort_switch_sorted_active": True,  # Start with sorted view active
        }

        return box.plot(fov_quality_data, config)
