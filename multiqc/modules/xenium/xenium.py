import json
import logging
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import polars as pl

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, box, linegraph, scatter, violin

# Try importing scipy, fallback gracefully if not available
try:
    import scipy
    import scipy.stats

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

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

    NOTE: parsing huge files is not an intended MultiQC usage. By default, MultiQC will ignore the `*.parquet` files
    as they are gigabyte-sized. To enable parsing those, make sure to have this line in your config:

    ```
    log_filesize_limit: 5000000000 # 5GB
    ```
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
        transcript_data_by_sample = self.ignore_samples(transcript_data_by_sample)
        cells_data_by_sample = self.ignore_samples(cells_data_by_sample)

        if len(data_by_sample) == 0 and len(transcript_data_by_sample) == 0 and len(cells_data_by_sample) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(data_by_sample)} Xenium reports")

        # Check for QC issues and add warnings
        self.check_qc_warnings(data_by_sample)

        # Add software version info (Xenium files don't contain version info)
        for s_name in data_by_sample.keys():
            self.add_software_version(None, s_name)

        # Merge cell area metrics into main data for general stats
        for sample_name, cell_data in cells_data_by_sample.items():
            if sample_name in data_by_sample:
                # Add cell area metrics to existing sample data
                if "cell_area_mean" in cell_data:
                    data_by_sample[sample_name]["cell_area_mean"] = cell_data["cell_area_mean"]
                if "nucleus_area_mean" in cell_data:
                    data_by_sample[sample_name]["nucleus_area_mean"] = cell_data["nucleus_area_mean"]
                if "nucleus_to_cell_area_ratio_mean" in cell_data:
                    data_by_sample[sample_name]["nucleus_to_cell_area_ratio_mean"] = cell_data[
                        "nucleus_to_cell_area_ratio_mean"
                    ]
            elif cell_data:
                # Create new sample entry if only cell data exists
                data_by_sample[sample_name] = {}
                if "cell_area_mean" in cell_data:
                    data_by_sample[sample_name]["cell_area_mean"] = cell_data["cell_area_mean"]
                if "nucleus_area_mean" in cell_data:
                    data_by_sample[sample_name]["nucleus_area_mean"] = cell_data["nucleus_area_mean"]
                if "nucleus_to_cell_area_ratio_mean" in cell_data:
                    data_by_sample[sample_name]["nucleus_to_cell_area_ratio_mean"] = cell_data[
                        "nucleus_to_cell_area_ratio_mean"
                    ]

        # Write parsed data to a file
        self.write_data_file(data_by_sample, "multiqc_xenium")

        # Add key metrics to general stats
        self.xenium_general_stats_table(data_by_sample)

        # Create plots - Cell detection metrics are already in general stats table

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

            # Add violin plot for transcript quality by codeword category if available
            violin_plot = self.xenium_transcript_category_violin_plot(transcript_data_by_sample)
            if violin_plot is not None:
                self.add_section(
                    name="Transcript Quality by Codeword Category",
                    anchor="xenium-transcript-category-violin",
                    description="Distribution of transcript quality values by codeword category (averaged across samples)",
                    helptext="""
                    This violin plot shows transcript quality (QV score) distributions for each codeword category:
                    
                    **Codeword Categories:**
                    * **predesigned_gene**: Standard genes from Xenium panels  
                    * **custom_gene**: User-added custom targets
                    * **negative_control_probe**: Control probes for background estimation
                    * **negative_control_codeword**: Control codewords for background estimation
                    * **genomic_control_probe**: Genomic DNA controls
                    * **unassigned_codeword**: Unassigned transcripts
                    
                    **Quality Score (QV) Interpretation:**
                    * QV ≥20: High-quality transcripts (≥99% accuracy)
                    * QV 10-20: Medium quality (90-99% accuracy) 
                    * QV <10: Low-quality transcripts (<90% accuracy)
                    
                    **What to look for:**
                    * **predesigned_gene** should have the highest QV distributions (centered around 30-40)
                    * **custom_gene** typically has slightly lower but still good QV distributions
                    * **Negative controls** should have variable quality and lower overall counts
                    * **Large differences** between categories may indicate panel design issues
                    
                    **Multiple samples:**
                    * For multiple samples, data is averaged across all samples to show overall category patterns
                    * The violin shape shows the density of QV values at each quality level
                    * Wider sections indicate more transcripts at that quality level
                    """,
                    plot=violin_plot,
                )

            # Add molecules per gene distribution if available
            molecules_plot = self.xenium_molecules_per_gene_plot(transcript_data_by_sample)
            if molecules_plot is not None:
                self.add_section(
                    name="Molecules per Gene Distribution",
                    anchor="xenium-molecules-per-gene",
                    description="Distribution of transcript counts per gene",
                    helptext="""
                    This histogram shows the distribution of transcript counts per gene across all samples:
                    
                    **What it shows:**
                    * **X-axis**: Number of transcripts per gene (log scale)
                    * **Y-axis**: Number of genes with that transcript count
                    * **Two categories**: Genes vs. non-genes (controls, unassigned, etc.)
                    
                    **Interpretation:**
                    * **Most genes** should have moderate transcript counts (hundreds to thousands)
                    * **Controls and non-genes** typically have lower counts
                    * **Very high counts** may indicate highly expressed genes or technical artifacts
                    * **Very low counts** may indicate poorly detected genes
                    
                    **What to look for:**
                    * **Smooth distribution** for genes with a peak in the hundreds-thousands range
                    * **Lower counts** for non-gene features (controls)  
                    * **No extreme outliers** unless biologically expected
                    * **Consistent patterns** across similar tissue types
                    
                    **Quality indicators:**
                    * Peak gene expression around 100-10,000 transcripts per gene is typical
                    * Clear separation between gene and non-gene distributions
                    * Absence of unusual spikes or gaps in the distribution
                    """,
                    plot=molecules_plot,
                )

        # Add cell area distribution section if cells data is available
        if cells_data_by_sample and SCIPY_AVAILABLE:
            area_plot = self.xenium_cell_area_distribution_plot(cells_data_by_sample)
            if area_plot:
                self.add_section(
                    name="Cell Area Distribution",
                    anchor="xenium-cell-area-distribution",
                    description="Distribution of cell areas across samples",
                    helptext="""
                    This plot shows the distribution of cell areas in the sample(s):
                    
                    **Single sample**: Density plot with vertical lines showing mean and median cell area
                    **Multiple samples**: Violin plots showing the distribution for each sample
                    
                    **Typical cell area ranges (tissue-dependent):**
                    * **Most tissues**: 50-200 μm²
                    * **Large cells** (e.g., neurons): 200-500 μm²
                    * **Small cells** (e.g., lymphocytes): 20-80 μm²
                    
                    **What to look for:**
                    * **Consistent distributions** across samples of the same tissue type
                    * **Biologically reasonable values** for your tissue
                    * **Outliers**: Very large or small cells may indicate segmentation issues
                    
                    **Troubleshooting:**
                    * Bimodal distributions: May indicate mixed cell types or segmentation artifacts
                    * Very large cells: Over-segmentation, cell doublets, or debris
                    * Very small cells: Under-segmentation, nuclear fragments
                    """,
                    plot=area_plot,
                )

            # Add nucleus RNA fraction distribution plot (scipy required)
            nucleus_plot = self.xenium_nucleus_rna_fraction_plot(cells_data_by_sample)
            if nucleus_plot:
                self.add_section(
                    name="Distribution of fractions of molecules in nucleus per cell",
                    anchor="xenium-nucleus-rna-fraction",
                    description="Distribution of nucleus RNA molecule fractions across cells",
                    helptext="""
                    This plot shows the distribution of the fraction of RNA molecules located in the nucleus versus cytoplasm for each cell:
                    
                    **Single sample**: Density plot showing the distribution of nucleus RNA fractions
                    **Multiple samples**: Box plots comparing distributions across samples
                    
                    **Biological interpretation:**
                    * **Low values (0.0-0.2)**: Most RNA is cytoplasmic (expected for mature mRNAs)
                    * **High values (>0.5)**: High nuclear retention (may indicate processing issues)
                    * **Peak around 0.0-0.1**: Normal for most cell types with efficient RNA export
                    
                    **What to look for:**
                    * **Consistent distributions** across samples of the same tissue type
                    * **Biologically reasonable values** for your cell types
                    * **Sample differences**: May reflect cell type composition or processing efficiency
                    
                    **Troubleshooting:**
                    * Very high nuclear fractions: Check for nuclear segmentation issues
                    * Bimodal distributions: May indicate different cell types or states
                    * Outliers: Individual cells with unusual RNA localization patterns
                    """,
                    plot=nucleus_plot,
                )

            # Add nucleus-to-cell area ratio distribution plot
            ratio_plot = self.xenium_nucleus_cell_area_ratio_plot(cells_data_by_sample)
            if ratio_plot:
                self.add_section(
                    name="Nucleus to cell area distribution",
                    anchor="xenium-nucleus-cell-area-ratio",
                    description="Distribution of nucleus-to-cell area ratios across cells",
                    helptext="""
                    This plot shows the distribution of the ratio between nucleus area and total cell area for each cell:
                    
                    **Single sample**: Density plot showing the distribution of nucleus-to-cell area ratios
                    **Multiple samples**: Box plots comparing distributions across samples
                    
                    **Biological interpretation:**
                    * **Typical range**: 0.2-0.6 for most cell types
                    * **Low values (<0.2)**: Small nucleus relative to cell (may indicate active/mature cells)
                    * **High values (>0.6)**: Large nucleus relative to cell (may indicate dividing or stressed cells)
                    * **Peak around 0.3-0.5**: Normal for most healthy cell types
                    
                    **What to look for:**
                    * **Consistent distributions** across samples of the same tissue type
                    * **Biologically reasonable values** for your cell types
                    * **Sample differences**: May reflect different cell states or tissue composition
                    
                    **Quality assessment:**
                    * Very low ratios: May indicate over-segmented cells or debris
                    * Very high ratios: May indicate under-segmented cells or nuclear fragments
                    * Bimodal distributions: May indicate different cell types or segmentation artifacts
                    
                    **Troubleshooting:**
                    * Unusual distributions may suggest issues with nuclear or cell segmentation parameters
                    * Consider tissue-specific expected ranges when evaluating results
                    """,
                    plot=ratio_plot,
                )

            # Add combined cell distribution plot (transcripts and genes per cell)
            combined_plot = self.xenium_cell_distributions_combined_plot(cells_data_by_sample)
            if combined_plot:
                self.add_section(
                    name="Cell Distribution Analysis",
                    anchor="xenium-cell-distributions",
                    description="Distribution of transcripts and detected genes per cell",
                    helptext="""
                    This plot shows two key cell-level distributions with separate tabs/datasets:
                    
                    **Tab 1: Transcripts per cell** - Shows the distribution of total transcript counts per cell
                    **Tab 2: Detected genes per cell** - Shows the distribution of unique genes detected per cell
                    
                    **Plot types:**
                    * **Single sample**: Density plots showing the distribution shapes
                    * **Multiple samples**: Box plots comparing distributions across samples
                    
                    **Transcripts per cell interpretation:**
                    * **Typical range**: 100-5000 transcripts per cell for most tissues
                    * **High transcript counts**: Metabolically active cells or large cell types
                    * **Low transcript counts**: Less active cells, technical dropouts, or small cell fragments
                    * **Quality thresholds**: <50 may indicate poor segmentation, >10,000 may indicate doublets
                    
                    **Detected genes per cell interpretation:**
                    * **Typical range**: 50-2000 genes per cell depending on cell type and panel size
                    * **High gene counts**: Metabolically active cells or cells with high expression diversity
                    * **Low gene counts**: Specialized cells, inactive cells, or technical dropouts
                    * **Quality thresholds**: <20 may indicate poor cells or debris
                    
                    **What to look for:**
                    * **Unimodal distributions**: Expected for homogeneous cell populations
                    * **Multimodal distributions**: May indicate different cell types or technical artifacts
                    * **Sample consistency**: Similar distributions expected for replicate samples
                    * **Positive correlation**: Generally expect transcripts and genes per cell to correlate
                    
                    **Panel considerations:**
                    * **Pre-designed panels**: Gene counts limited by panel design (typically 100-1000 genes)
                    * **Custom panels**: Consider gene selection bias when interpreting results
                    * **Detection efficiency**: Some genes may be harder to detect than others
                    
                    **Quality assessment:**
                    * **Transcripts**: Very low (<50) or very high (>10,000) may indicate segmentation issues
                    * **Genes**: Very low (<20) may indicate poor cells, counts near panel size may indicate artifacts
                    * **Shoulder distributions**: May indicate presence of different cell types
                    
                    **Troubleshooting:**
                    * Unusual distributions may suggest issues with transcript detection or cell segmentation
                    * Consider cell type and tissue context when evaluating expected ranges
                    * Low gene detection may suggest transcript assignment issues
                    """,
                    plot=combined_plot,
                )

        # Add Field of View quality section at the end if FoV data is available
        if transcript_data_by_sample:
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
        file_path = Path(f["root"]) / f["fn"]

        # Use lazy loading to avoid reading entire file into memory
        try:
            df_lazy = pl.scan_parquet(file_path)

            # Check if required columns exist by scanning schema (avoid performance warning)
            schema = df_lazy.collect_schema()
            required_cols = ["qv", "feature_name"]
            if not all(col in schema for col in required_cols):
                log.warning(f"Missing required columns in {f['fn']}: {required_cols}")
                return None

            # Get total row count efficiently without loading full data
            total_transcripts = df_lazy.select(pl.len()).collect().item()

            # Group by feature_name and calculate statistics efficiently
            grouped = (
                df_lazy.group_by("feature_name")
                .agg(
                    [
                        pl.col("qv").value_counts().alias("qv_counts"),
                        pl.col("qv").mean().alias("mean_qv"),
                        pl.len().alias("transcript_count"),
                    ]
                )
                .collect(streaming=True)  # Use streaming to reduce memory usage
            )

            # Process results
            quality_dist = {}
            transcript_counts = {}

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
                "total_transcripts": total_transcripts,
                "unique_features": len(quality_dist),
            }
        except Exception as e:
            log.warning(f"Lazy processing failed for {file_path}: {e}. Falling back to direct read.")
            # Fallback to direct read for problematic files
            try:
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

                # Handle optional columns with fallback
                if "codeword_category" in df.columns:
                    try:
                        category_grouped = df.group_by("codeword_category").agg(pl.col("qv").alias("qv_values"))
                        category_quality_distributions = {}
                        for row in category_grouped.iter_rows(named=True):
                            category = str(row["codeword_category"])
                            qv_values = row["qv_values"]
                            if len(qv_values) > 2000:
                                sampled_qv = qv_values.sample(2000, seed=42)
                                category_quality_distributions[category] = sampled_qv.to_list()
                            else:
                                category_quality_distributions[category] = qv_values.to_list()
                        result["category_quality_distributions"] = category_quality_distributions
                    except Exception:
                        log.warning("Could not process codeword category quality data")

                if "is_gene" in df.columns:
                    try:
                        gene_grouped = df.group_by("feature_name").agg(
                            [pl.len().alias("molecule_count"), pl.col("is_gene").first().alias("is_gene")]
                        )
                        molecules_per_gene = {}
                        for row in gene_grouped.iter_rows(named=True):
                            feature_name = str(row["feature_name"])
                            molecules_per_gene[feature_name] = {
                                "count": row["molecule_count"],
                                "is_gene": row["is_gene"],
                            }
                        result["molecules_per_gene"] = molecules_per_gene
                    except Exception:
                        log.warning("Could not process molecules per gene data")

                if "fov_name" in df.columns:
                    try:
                        fov_grouped = df.group_by("fov_name").agg(
                            [
                                pl.col("qv").mean().alias("mean_qv"),
                                pl.col("qv").median().alias("median_qv"),
                                pl.col("qv").std().alias("std_qv"),
                                pl.len().alias("transcript_count"),
                            ]
                        )
                        fov_quality_stats = {}
                        for row in fov_grouped.iter_rows(named=True):
                            fov_name = str(row["fov_name"])
                            fov_quality_stats[fov_name] = {
                                "mean_quality": row["mean_qv"],
                                "median_quality": row["median_qv"],
                                "std_quality": row["std_qv"],
                                "transcript_count": row["transcript_count"],
                            }
                        result["fov_quality_stats"] = fov_quality_stats

                        # Sample FoV distributions
                        fov_dist_grouped = df.group_by("fov_name").agg(pl.col("qv").alias("qv_values"))
                        fov_quality_distributions = {}
                        for row in fov_dist_grouped.iter_rows(named=True):
                            fov_name = str(row["fov_name"])
                            qv_values = row["qv_values"]
                            if len(qv_values) > 1000:
                                sampled_qv = qv_values.sample(1000, seed=42)
                                fov_quality_distributions[fov_name] = sampled_qv.to_list()
                            else:
                                fov_quality_distributions[fov_name] = qv_values.to_list()
                        result["fov_quality_distributions"] = fov_quality_distributions
                    except Exception:
                        log.warning("Could not process FoV quality data")

                return result

            except Exception as fallback_error:
                log.error(f"Both lazy and direct processing failed for {file_path}: {fallback_error}")
                return None

        # Add codeword category quality analysis if codeword_category column is present
        if "codeword_category" in schema:
            try:
                # Group by codeword_category and sample QV values efficiently for violin plots
                category_grouped = (
                    df_lazy.group_by("codeword_category")
                    .agg(pl.col("qv").sample(n=2000, seed=42, with_replacement=True).alias("qv_sample"))
                    .collect(streaming=True)
                )

                category_quality_distributions = {}
                for row in category_grouped.iter_rows(named=True):
                    category = str(row["codeword_category"])
                    qv_values = row["qv_sample"]
                    category_quality_distributions[category] = (
                        qv_values.to_list() if hasattr(qv_values, "to_list") else qv_values
                    )

                result["category_quality_distributions"] = category_quality_distributions
            except Exception as e:
                log.warning(f"Could not process codeword category quality data: {e}")

        # Add molecules per gene analysis if feature_name and is_gene columns are present
        if "is_gene" in schema:
            try:
                gene_grouped = (
                    df_lazy.group_by("feature_name")
                    .agg([pl.len().alias("molecule_count"), pl.col("is_gene").first().alias("is_gene")])
                    .collect(streaming=True)
                )

                molecules_per_gene = {}
                for row in gene_grouped.iter_rows(named=True):
                    feature_name = str(row["feature_name"])
                    molecule_count = row["molecule_count"]
                    is_gene = row["is_gene"]
                    molecules_per_gene[feature_name] = {"count": molecule_count, "is_gene": is_gene}

                result["molecules_per_gene"] = molecules_per_gene
            except Exception as e:
                log.warning(f"Could not process molecules per gene data: {e}")

        # Add FoV quality analysis if fov_name column is present
        if "fov_name" in schema:
            try:
                # Group by FoV and calculate quality stats efficiently
                fov_stats_grouped = (
                    df_lazy.group_by("fov_name")
                    .agg(
                        [
                            pl.col("qv").mean().alias("mean_qv"),
                            pl.col("qv").median().alias("median_qv"),
                            pl.col("qv").std().alias("std_qv"),
                            pl.len().alias("transcript_count"),
                        ]
                    )
                    .collect(streaming=True)
                )

                fov_quality_stats = {}
                for row in fov_stats_grouped.iter_rows(named=True):
                    fov_name = str(row["fov_name"])
                    fov_quality_stats[fov_name] = {
                        "mean_quality": row["mean_qv"],
                        "median_quality": row["median_qv"],
                        "std_quality": row["std_qv"],
                        "transcript_count": row["transcript_count"],
                    }

                # Separate query for sampled QV distributions to avoid memory issues
                fov_dist_grouped = (
                    df_lazy.group_by("fov_name")
                    .agg(pl.col("qv").sample(n=1000, seed=42, with_replacement=True).alias("qv_sample"))
                    .collect(streaming=True)
                )

                fov_quality_distributions = {}
                for row in fov_dist_grouped.iter_rows(named=True):
                    fov_name = str(row["fov_name"])
                    qv_values = row["qv_sample"]
                    fov_quality_distributions[fov_name] = (
                        qv_values.to_list() if hasattr(qv_values, "to_list") else qv_values
                    )

                result["fov_quality_stats"] = fov_quality_stats
                result["fov_quality_distributions"] = fov_quality_distributions
            except Exception as e:
                log.warning(f"Could not process FoV quality data: {e}")

        return result

    def parse_cells_parquet(self, f) -> Optional[Dict]:
        """Parse Xenium cells.parquet file to extract cell-level metrics"""
        file_path = Path(f["root"]) / f["fn"]

        # Due to persistent Polars panic issues with certain parquet files,
        # we'll use direct reading approach for cells files
        log.info(f"Processing cells parquet file with direct read: {file_path}")
        try:
            df = pl.read_parquet(file_path)

            # Check for required columns
            required_cols = ["cell_area", "nucleus_area", "total_counts", "transcript_counts"]
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                log.warning(f"Missing columns in {f['fn']}: {missing_cols}")
                return None

            # Calculate summary statistics
            cell_stats = {"total_cells": df.height}

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
                # Sample cell area values for distribution plots
                if cell_area_stats.len() > 10000:
                    cell_stats["cell_area_values"] = cell_area_stats.sample(10000, seed=42).to_list()
                else:
                    cell_stats["cell_area_values"] = cell_area_stats.to_list()

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
                    # Sample ratio values for distribution plots
                    if ratio.len() > 10000:
                        cell_stats["nucleus_to_cell_area_ratio_values"] = ratio.sample(10000, seed=42).to_list()
                    else:
                        cell_stats["nucleus_to_cell_area_ratio_values"] = ratio.to_list()

            # Store transcript counts per cell (total_counts) for distribution plots
            total_counts_stats = df["total_counts"].drop_nulls()
            if total_counts_stats.len() > 0:
                if total_counts_stats.len() > 10000:
                    cell_stats["transcript_counts_values"] = total_counts_stats.sample(10000, seed=42).to_list()
                else:
                    cell_stats["transcript_counts_values"] = total_counts_stats.to_list()

            # Store detected genes per cell (transcript_counts) for distribution plots
            detected_genes_stats = df["transcript_counts"].drop_nulls()
            if detected_genes_stats.len() > 0:
                if detected_genes_stats.len() > 10000:
                    cell_stats["detected_genes_values"] = detected_genes_stats.sample(10000, seed=42).to_list()
                else:
                    cell_stats["detected_genes_values"] = detected_genes_stats.to_list()

            # Add nucleus RNA fraction if nucleus_count is available
            if "nucleus_count" in df.columns:
                valid_cells = df.filter(pl.col("total_counts") > 0)
                if valid_cells.height > 0:
                    nucleus_rna_fraction = valid_cells["nucleus_count"] / valid_cells["total_counts"]
                    cell_stats.update(
                        {
                            "nucleus_rna_fraction_mean": nucleus_rna_fraction.mean(),
                            "nucleus_rna_fraction_median": nucleus_rna_fraction.median(),
                        }
                    )
                    # Sample nucleus fraction values for distribution plots
                    if nucleus_rna_fraction.len() > 10000:
                        cell_stats["nucleus_rna_fraction_values"] = nucleus_rna_fraction.sample(
                            10000, seed=42
                        ).to_list()
                    else:
                        cell_stats["nucleus_rna_fraction_values"] = nucleus_rna_fraction.to_list()

            return cell_stats

        except Exception as e:
            log.error(f"Error processing cells parquet file {file_path}: {e}")
            return None

            # Sample distribution data for plots using separate efficient queries
            # Cell area distribution sample
            try:
                cell_area_sample = (
                    df_lazy.filter(pl.col("cell_area").is_not_null())
                    .select(pl.col("cell_area").sample(n=10000, seed=42, with_replacement=True))
                    .collect()
                )
                if cell_area_sample.height > 0:
                    cell_stats["cell_area_values"] = cell_area_sample["cell_area"].to_list()
            except Exception as e:
                log.warning(f"Could not sample cell area values: {e}")

            # Transcript counts distribution sample
            try:
                transcript_sample = (
                    df_lazy.filter(pl.col("total_counts").is_not_null())
                    .select(pl.col("total_counts").sample(n=10000, seed=42, with_replacement=True))
                    .collect()
                )
                if transcript_sample.height > 0:
                    cell_stats["transcript_counts_values"] = transcript_sample["total_counts"].to_list()
            except Exception as e:
                log.warning(f"Could not sample transcript counts: {e}")

            # Detected genes distribution sample
            try:
                genes_sample = (
                    df_lazy.filter(pl.col("transcript_counts").is_not_null())
                    .select(pl.col("transcript_counts").sample(n=10000, seed=42, with_replacement=True))
                    .collect()
                )
                if genes_sample.height > 0:
                    cell_stats["detected_genes_values"] = genes_sample["transcript_counts"].to_list()
            except Exception as e:
                log.warning(f"Could not sample detected genes: {e}")

            # Nucleus RNA fraction if nucleus_count is available
            if "nucleus_count" in schema:
                try:
                    nucleus_stats = (
                        df_lazy.filter(pl.col("total_counts") > 0)
                        .with_columns((pl.col("nucleus_count") / pl.col("total_counts")).alias("nucleus_fraction"))
                        .select(
                            [
                                pl.col("nucleus_fraction").mean().alias("nucleus_rna_fraction_mean"),
                                pl.col("nucleus_fraction").median().alias("nucleus_rna_fraction_median"),
                                pl.col("nucleus_fraction")
                                .sample(n=10000, seed=42, with_replacement=True)
                                .alias("nucleus_fraction_sample"),
                            ]
                        )
                        .collect()
                    )

                    if nucleus_stats.height > 0:
                        for stat, value in (
                            nucleus_stats.select(["nucleus_rna_fraction_mean", "nucleus_rna_fraction_median"])
                            .to_dicts()[0]
                            .items()
                        ):
                            if value is not None:
                                cell_stats[stat] = value
                        cell_stats["nucleus_rna_fraction_values"] = nucleus_stats["nucleus_fraction_sample"].to_list()
                except Exception as e:
                    log.warning(f"Could not process nucleus RNA fraction: {e}")

            # Nucleus to cell area ratio distribution sample
            try:
                ratio_sample = (
                    df_lazy.filter(
                        (pl.col("cell_area").is_not_null())
                        & (pl.col("nucleus_area").is_not_null())
                        & (pl.col("cell_area") > 0)
                    )
                    .with_columns((pl.col("nucleus_area") / pl.col("cell_area")).alias("ratio"))
                    .select(pl.col("ratio").sample(n=10000, seed=42, with_replacement=True))
                    .collect()
                )
                if ratio_sample.height > 0:
                    cell_stats["nucleus_to_cell_area_ratio_values"] = ratio_sample["ratio"].to_list()
            except Exception as e:
                log.warning(f"Could not sample nucleus to cell area ratio: {e}")

            return cell_stats

        except Exception as e:
            log.warning(f"Lazy processing failed for cells file {file_path}: {e}. Falling back to direct read.")
            # Fallback to direct read for problematic files
            try:
                df = pl.read_parquet(file_path)

                # Check for required columns
                required_cols = ["cell_area", "nucleus_area", "total_counts", "transcript_counts"]
                missing_cols = [col for col in required_cols if col not in df.columns]
                if missing_cols:
                    log.warning(f"Missing columns in {f['fn']}: {missing_cols}")
                    return None

                # Calculate summary statistics
                cell_stats = {"total_cells": df.height}

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
                    # Sample cell area values
                    if cell_area_stats.len() > 10000:
                        cell_stats["cell_area_values"] = cell_area_stats.sample(10000, seed=42).to_list()
                    else:
                        cell_stats["cell_area_values"] = cell_area_stats.to_list()

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

                    # Nucleus to cell area ratio
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
                        # Sample ratio values
                        if ratio.len() > 10000:
                            cell_stats["nucleus_to_cell_area_ratio_values"] = ratio.sample(10000, seed=42).to_list()
                        else:
                            cell_stats["nucleus_to_cell_area_ratio_values"] = ratio.to_list()

                # Store transcript counts per cell for distribution plots
                total_counts_stats = df["total_counts"].drop_nulls()
                if total_counts_stats.len() > 0:
                    if total_counts_stats.len() > 10000:
                        cell_stats["transcript_counts_values"] = total_counts_stats.sample(10000, seed=42).to_list()
                    else:
                        cell_stats["transcript_counts_values"] = total_counts_stats.to_list()

                # Store detected genes per cell for distribution plots
                detected_genes_stats = df["transcript_counts"].drop_nulls()
                if detected_genes_stats.len() > 0:
                    if detected_genes_stats.len() > 10000:
                        cell_stats["detected_genes_values"] = detected_genes_stats.sample(10000, seed=42).to_list()
                    else:
                        cell_stats["detected_genes_values"] = detected_genes_stats.to_list()

                # Add nucleus RNA fraction if nucleus_count is available
                if "nucleus_count" in df.columns:
                    valid_cells = df.filter(pl.col("total_counts") > 0)
                    if valid_cells.height > 0:
                        nucleus_rna_fraction = valid_cells["nucleus_count"] / valid_cells["total_counts"]
                        cell_stats.update(
                            {
                                "nucleus_rna_fraction_mean": nucleus_rna_fraction.mean(),
                                "nucleus_rna_fraction_median": nucleus_rna_fraction.median(),
                            }
                        )
                        # Sample nucleus fraction values
                        if nucleus_rna_fraction.len() > 10000:
                            cell_stats["nucleus_rna_fraction_values"] = nucleus_rna_fraction.sample(
                                10000, seed=42
                            ).to_list()
                        else:
                            cell_stats["nucleus_rna_fraction_values"] = nucleus_rna_fraction.to_list()

                return cell_stats

            except Exception as fallback_error:
                log.error(f"Both lazy and direct processing failed for cells file {file_path}: {fallback_error}")
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
            "cell_area_mean": {
                "title": "Cell Area",
                "description": "Mean cell area",
                "suffix": " μm²",
                "scale": "Blues",
                "format": "{:,.1f}",
            },
            "nucleus_area_mean": {
                "title": "Nucleus Area",
                "description": "Mean nucleus area",
                "suffix": " μm²",
                "scale": "Oranges",
                "format": "{:,.1f}",
            },
            "nucleus_to_cell_area_ratio_mean": {
                "title": "Nucleus/Cell Ratio",
                "description": "Mean nucleus to cell area ratio",
                "scale": "Greens",
                "format": "{:.3f}",
                "max": 1.0,
            },
        }
        self.general_stats_addcols(data_by_sample, headers)

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
                category, feature_id = categorize_feature(feature)

                if category not in plot_data:
                    plot_data[category] = []

                # Each point is a separate data point
                # For multiple samples, include sample name in the hover text
                if len(transcript_data_by_sample) > 1:
                    point_name = f"{feature_id} ({sample_name})"
                else:
                    point_name = feature_id

                plot_data[category].append(
                    {
                        "x": counts_data["count"],
                        "y": counts_data["mean_quality"],
                        "name": point_name,  # Use gene name (+ sample) for hover text
                        "group": category,
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

        # Define desired category order for legend
        category_order = [
            "Pre-designed",
            "Custom",
            "Genomic Control Probe",
            "Negative Control Probe",
            "Negative Control Codeword",
            "Unassigned Codeword",
            "Deprecated Codeword",
        ]

        config = {
            "id": "xenium_transcript_quality_combined",
            "title": title,
            "xlab": "Total transcripts per gene",
            "ylab": "Mean calibrated quality of gene transcripts",
            "marker_size": 5,
            "series_label": "transcripts",
            "xlog": True,
            "showlegend": True,
            "groups": category_order,
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
                    category, _ = categorize_feature(feature)
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
                    if categorize_feature(feature)[0] == cat:
                        datasets[cat][sample_name][count] = quality

        if not datasets:
            return None

        config = {
            "id": "xenium_transcript_quality_multi",
            "title": "Xenium: Transcript Quality by Category",
            "xlab": "Total transcripts per gene",
            "ylab": "Mean calibrated quality of gene transcripts",
            "data_labels": [{"name": cat} for cat in datasets.keys()],
            "xlog": True,
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

    def xenium_cell_area_distribution_plot(self, cells_data_by_sample):
        """Create cell area distribution plot - line plot for single sample, violin plots for multiple"""
        # Check which samples have cell area data
        samples_with_areas = []
        for s_name, data in cells_data_by_sample.items():
            if "cell_area_values" in data:
                samples_with_areas.append(s_name)

        if not samples_with_areas:
            return None

        num_samples = len(samples_with_areas)

        if num_samples == 1:
            # Single sample: Create line plot (density) with vertical lines for mean/median
            return self._create_single_sample_area_density(cells_data_by_sample[samples_with_areas[0]])
        else:
            # Multiple samples: Create violin plots
            return self._create_multi_sample_area_violins(cells_data_by_sample, samples_with_areas)

    def _create_single_sample_area_density(self, cell_data):
        """Create density plot for single sample with mean/median lines"""
        if not SCIPY_AVAILABLE:
            log.warning("scipy not available, skipping density plots. Install scipy for enhanced plotting.")
            return None

        import numpy as np

        if SCIPY_AVAILABLE:
            from scipy.stats import gaussian_kde

        cell_areas = cell_data["cell_area_values"]
        if not cell_areas or len(cell_areas) < 10:
            return None

        # Create density estimation
        kde = gaussian_kde(cell_areas)

        # Create x values for density curve
        min_area = min(cell_areas)
        max_area = max(cell_areas)
        x_range = max_area - min_area
        x_vals = np.linspace(max(0, min_area - 0.1 * x_range), max_area + 0.1 * x_range, 200)
        density_vals = kde(x_vals)

        # Prepare data for linegraph
        datasets = {}
        density_data = {}
        for x, y in zip(x_vals, density_vals):
            density_data[float(x)] = float(y)
        datasets["Density"] = density_data

        config = {
            "id": "xenium_cell_area_distribution",
            "title": "Xenium: Cell Area Distribution",
            "xlab": "Cell area",
            "ylab": "Density",
            "xsuffix": " μm²",
        }

        # Add vertical lines for mean and median
        if "cell_area_mean" in cell_data and "cell_area_median" in cell_data:
            config["x_lines"] = [
                {
                    "value": float(cell_data["cell_area_mean"]),
                    "color": "red",
                    "dash": "dash",
                    "width": 2,
                    "label": f"Mean ({cell_data['cell_area_mean']:.1f} μm²)",
                },
                {
                    "value": float(cell_data["cell_area_median"]),
                    "color": "green",
                    "dash": "dash",
                    "width": 2,
                    "label": f"Median ({cell_data['cell_area_median']:.1f} μm²)",
                },
            ]

        return linegraph.plot(datasets, config)

    def _create_multi_sample_area_violins(self, cells_data_by_sample, samples_with_areas):
        """Create box plots for multiple samples - one box per sample"""
        from multiqc.plots import box

        # For box plots, we provide the raw data points grouped by sample
        data = {}

        for s_name in samples_with_areas:
            cell_data = cells_data_by_sample[s_name]
            if "cell_area_values" in cell_data:
                # Store all cell area values for this sample
                cell_areas = cell_data["cell_area_values"]
                if cell_areas:
                    # Box plots expect raw data points as a list
                    data[s_name] = [float(area) for area in cell_areas]

        if not data:
            return None

        config = {
            "id": "xenium_cell_area_distribution",
            "title": "Xenium: Cell Area Distribution",
            "ylab": "Cell area (μm²)",
            "xlab": "Sample",
            "boxpoints": False,
        }

        return box.plot(data, config)

    def xenium_nucleus_rna_fraction_plot(self, cells_data_by_sample):
        """Create nucleus RNA fraction distribution plot - density for single sample, box plots for multiple"""
        # Check which samples have nucleus RNA fraction data
        samples_with_nucleus_data = []
        for s_name, data in cells_data_by_sample.items():
            if "nucleus_rna_fraction_values" in data and data["nucleus_rna_fraction_values"]:
                samples_with_nucleus_data.append(s_name)

        if not samples_with_nucleus_data:
            return None

        num_samples = len(samples_with_nucleus_data)

        if num_samples == 1:
            # Single sample: Create density plot
            return self._create_single_sample_nucleus_density(cells_data_by_sample[samples_with_nucleus_data[0]])
        else:
            # Multiple samples: Create box plots
            return self._create_multi_sample_nucleus_boxes(cells_data_by_sample, samples_with_nucleus_data)

    def _create_single_sample_nucleus_density(self, cell_data):
        """Create density plot for single sample nucleus RNA fractions"""
        if not SCIPY_AVAILABLE:
            log.warning("scipy not available, skipping nucleus density plots. Install scipy for enhanced plotting.")
            return None

        import numpy as np
        from scipy import stats

        from multiqc.plots import linegraph

        nucleus_fractions = cell_data["nucleus_rna_fraction_values"]
        if not nucleus_fractions:
            return None

        # Use a more appropriate range for nucleus RNA fractions (0 to 1)
        x_range = np.linspace(0, 1, 200)

        # Calculate kernel density estimation
        try:
            kde = stats.gaussian_kde(nucleus_fractions)
            density = kde(x_range)
        except Exception:
            # Fallback to histogram if KDE fails
            hist, bin_edges = np.histogram(nucleus_fractions, bins=50, range=(0, 1), density=True)
            x_range = (bin_edges[:-1] + bin_edges[1:]) / 2
            density = hist

        # Create the density plot data
        data = {}
        data["Nucleus RNA Fraction Density"] = {str(x): y for x, y in zip(x_range, density)}

        # Note: Could add statistical lines (mean/median) in future if desired

        config = {
            "id": "xenium_nucleus_rna_fraction_single",
            "title": "Distribution of fractions of molecules in nucleus per cell",
            "xlab": "Fraction of molecules in nucleus per cell",
            "ylab": "Density",
            "data_labels": [
                {"name": "Density", "ylab": "Density"},
            ],
        }

        # Add vertical lines for mean and median
        plot = linegraph.plot(data, config)

        return plot

    def _create_multi_sample_nucleus_boxes(self, cells_data_by_sample, samples_with_nucleus_data):
        """Create box plots for multiple samples - one box per sample"""
        from multiqc.plots import box

        # For box plots, we provide the raw data points grouped by sample
        data = {}

        for s_name in samples_with_nucleus_data:
            cell_data = cells_data_by_sample[s_name]
            if "nucleus_rna_fraction_values" in cell_data:
                nucleus_fractions = cell_data["nucleus_rna_fraction_values"]
                if nucleus_fractions:
                    # Box plots expect raw data points as a list
                    data[s_name] = [float(fraction) for fraction in nucleus_fractions]

        if not data:
            return None

        config = {
            "id": "xenium_nucleus_rna_fraction_multi",
            "title": "Distribution of fractions of molecules in nucleus per cell",
            "ylab": "Fraction of molecules in nucleus per cell",
            "xlab": "Sample",
            "boxpoints": False,
        }

        return box.plot(data, config)

    def xenium_nucleus_cell_area_ratio_plot(self, cells_data_by_sample):
        """Create nucleus-to-cell area ratio distribution plot - density for single sample, box plots for multiple"""
        # Check which samples have nucleus-to-cell area ratio data
        samples_with_ratio_data = []
        for s_name, data in cells_data_by_sample.items():
            if "nucleus_to_cell_area_ratio_values" in data and data["nucleus_to_cell_area_ratio_values"]:
                samples_with_ratio_data.append(s_name)

        if not samples_with_ratio_data:
            return None

        num_samples = len(samples_with_ratio_data)

        if num_samples == 1:
            # Single sample: Create density plot
            return self._create_single_sample_ratio_density(cells_data_by_sample[samples_with_ratio_data[0]])
        else:
            # Multiple samples: Create box plots
            return self._create_multi_sample_ratio_boxes(cells_data_by_sample, samples_with_ratio_data)

    def _create_single_sample_ratio_density(self, cell_data):
        """Create density plot for single sample nucleus-to-cell area ratios"""
        if not SCIPY_AVAILABLE:
            log.warning("scipy not available, skipping plots. Install scipy for enhanced plotting.")
            return None

        import numpy as np
        from scipy import stats

        from multiqc.plots import linegraph

        ratio_values = cell_data["nucleus_to_cell_area_ratio_values"]
        if not ratio_values:
            return None

        # Use a reasonable range for nucleus-to-cell area ratios (0 to 1.0)
        x_range = np.linspace(0, 1.0, 200)

        # Calculate kernel density estimation
        try:
            kde = stats.gaussian_kde(ratio_values)
            density = kde(x_range)
        except Exception:
            # Fallback to histogram if KDE fails
            hist, bin_edges = np.histogram(ratio_values, bins=50, range=(0, 1.0), density=True)
            x_range = (bin_edges[:-1] + bin_edges[1:]) / 2
            density = hist

        # Create the density plot data
        data = {}
        data["Nucleus-to-Cell Area Ratio Density"] = {str(x): y for x, y in zip(x_range, density)}

        config = {
            "id": "xenium_nucleus_cell_area_ratio_single",
            "title": "Nucleus to cell area distribution",
            "xlab": "Nucleus-to-cell area ratio",
            "ylab": "Density",
            "data_labels": [
                {"name": "Density", "ylab": "Density"},
            ],
        }

        plot = linegraph.plot(data, config)

        return plot

    def _create_multi_sample_ratio_boxes(self, cells_data_by_sample, samples_with_ratio_data):
        """Create box plots for multiple samples - one box per sample"""
        from multiqc.plots import box

        # For box plots, we provide the raw data points grouped by sample
        data = {}

        for s_name in samples_with_ratio_data:
            cell_data = cells_data_by_sample[s_name]
            if "nucleus_to_cell_area_ratio_values" in cell_data:
                ratio_values = cell_data["nucleus_to_cell_area_ratio_values"]
                if ratio_values:
                    # Box plots expect raw data points as a list
                    data[s_name] = [float(ratio) for ratio in ratio_values]

        if not data:
            return None

        config = {
            "id": "xenium_nucleus_cell_area_ratio_multi",
            "title": "Nucleus to cell area distribution",
            "ylab": "Nucleus-to-cell area ratio",
            "xlab": "Sample",
            "boxpoints": False,
        }

        return box.plot(data, config)

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
            "boxpoints": False,  # Do not show individual data points
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
            "boxpoints": False,  # Do not show individual data points
        }

        return box.plot(fov_quality_data, config)

    def xenium_transcript_category_violin_plot(self, transcript_data_by_sample):
        """Create violin plot showing transcript quality distribution by codeword category"""
        # Check if any sample has category quality distributions
        samples_with_categories = []
        for s_name, data in transcript_data_by_sample.items():
            if "category_quality_distributions" in data:
                samples_with_categories.append(s_name)

        if not samples_with_categories:
            return None

        # Aggregate quality distributions across all samples for each category
        category_quality_data = {}

        for s_name in samples_with_categories:
            data = transcript_data_by_sample[s_name]
            if "category_quality_distributions" in data:
                category_distributions = data["category_quality_distributions"]
                for category, quality_values in category_distributions.items():
                    if category not in category_quality_data:
                        category_quality_data[category] = []
                    # Add all quality values from this sample's category to the aggregated distribution
                    category_quality_data[category].extend(quality_values)

        if not category_quality_data:
            return None

        # Create headers for the violin plot - each category is a "metric"
        headers = {}
        for category in category_quality_data.keys():
            headers[category] = {
                "title": category,
                "description": f"Quality score distribution for {category}",
                "suffix": " QV",
                "color": GENE_CATS.get(category, {}).get("color", "#888888"),
                "min": 0,
                "max": 50,  # QV scores typically range 0-50
            }

        # Create data dict - for violin plots we need sample -> category -> values
        # For multiple samples, we want to show category distributions per sample
        # but group all samples under each category tab
        data = {}

        # For each sample, create entries with category-wise data
        for s_name in samples_with_categories:
            sample_data = transcript_data_by_sample[s_name]
            if "category_quality_distributions" in sample_data:
                data[s_name] = sample_data["category_quality_distributions"]

        # If we only have one sample but want to show all categories together,
        # we can flatten the structure, otherwise keep sample-wise structure
        if len(samples_with_categories) == 1:
            # For single sample, use the aggregated data
            sample_name = samples_with_categories[0]
            data = {sample_name: category_quality_data}

        config = {
            "id": "xenium_transcript_category_violin",
            "title": "Xenium: Transcript Quality by Codeword Category (averaged across samples)",
            "col1_header": "Category",
            "series_label": "transcripts",
        }

        return violin.plot(data, headers, config)

    def xenium_molecules_per_gene_plot(self, transcript_data_by_sample):
        """Create histogram plot showing distribution of molecules per gene with separate lines per sample"""
        # Check if any sample has molecules per gene data
        samples_with_molecules = []
        for s_name, data in transcript_data_by_sample.items():
            if "molecules_per_gene" in data:
                samples_with_molecules.append(s_name)

        if not samples_with_molecules:
            return None

        import numpy as np

        # Calculate noise threshold based on aggregated data (for vertical line)
        aggregated_gene_counts = {}
        for s_name in samples_with_molecules:
            data = transcript_data_by_sample[s_name]
            molecules_data = data["molecules_per_gene"]
            for gene_name, gene_info in molecules_data.items():
                if gene_name not in aggregated_gene_counts:
                    aggregated_gene_counts[gene_name] = {"count": 0, "is_gene": gene_info["is_gene"]}
                aggregated_gene_counts[gene_name]["count"] += gene_info["count"]

        n_mols_threshold = self.calculate_noise_threshold(aggregated_gene_counts)

        # Determine global bins based on all samples' data
        all_gene_counts = []
        all_non_gene_counts = []

        for s_name in samples_with_molecules:
            data = transcript_data_by_sample[s_name]
            molecules_data = data["molecules_per_gene"]

            for gene_name, gene_info in molecules_data.items():
                count = gene_info["count"]
                if count > 0:
                    if gene_info["is_gene"]:
                        all_gene_counts.append(count)
                    else:
                        all_non_gene_counts.append(count)

        # Create consistent bins for all samples
        all_counts = all_gene_counts + all_non_gene_counts
        if not all_counts:
            return None

        min_count = max(1, min(all_counts))
        max_count = max(all_counts)
        bins = np.logspace(np.log10(min_count), np.log10(max_count), 50)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        # Check if single sample - handle differently
        num_samples = len(samples_with_molecules)

        if num_samples == 1:
            # Single sample: Put both Gene and Non-gene lines on the same plot
            return self._create_single_sample_molecules_plot(
                transcript_data_by_sample[samples_with_molecules[0]], bins, bin_centers, n_mols_threshold
            )
        else:
            # Multiple samples: Use tabs for Genes vs Non-genes
            return self._create_multi_sample_molecules_plot(
                transcript_data_by_sample, samples_with_molecules, bins, bin_centers, n_mols_threshold
            )

    def _create_single_sample_molecules_plot(self, sample_data, bins, bin_centers, n_mols_threshold):
        """Create single plot with both Gene and Non-gene lines for single sample"""
        import numpy as np

        molecules_data = sample_data["molecules_per_gene"]

        # Separate counts by gene type
        gene_counts = []
        non_gene_counts = []

        for gene_name, gene_info in molecules_data.items():
            count = gene_info["count"]
            if count > 0:
                if gene_info["is_gene"]:
                    gene_counts.append(count)
                else:
                    non_gene_counts.append(count)

        # Create plot data with both lines
        plot_data = {}

        if gene_counts:
            gene_hist, _ = np.histogram(gene_counts, bins=bins)
            gene_line_data = {}
            for i, count in enumerate(gene_hist):
                gene_line_data[float(bin_centers[i])] = int(count)
            plot_data["Genes"] = gene_line_data

        if non_gene_counts:
            non_gene_hist, _ = np.histogram(non_gene_counts, bins=bins)
            non_gene_line_data = {}
            for i, count in enumerate(non_gene_hist):
                non_gene_line_data[float(bin_centers[i])] = int(count)
            plot_data["Non-genes"] = non_gene_line_data

        if not plot_data:
            return None

        config = {
            "id": "xenium_molecules_per_gene",
            "title": "Xenium: Distribution of Molecules per Gene",
            "xlab": "Number of molecules per gene",
            "ylab": "Number of features",
        }

        # Add vertical line for noise threshold if calculated
        if n_mols_threshold is not None and n_mols_threshold > 0:
            config["x_lines"] = [
                {
                    "value": n_mols_threshold,
                    "color": "grey",
                    "dash": "dash",
                    "width": 1,
                    "label": f"Noise threshold ({n_mols_threshold:.0f})",
                }
            ]

        return linegraph.plot(plot_data, config)

    def _create_multi_sample_molecules_plot(
        self, transcript_data_by_sample, samples_with_molecules, bins, bin_centers, n_mols_threshold
    ):
        """Create tabbed plot with separate lines per sample for multiple samples"""
        import numpy as np

        # Create separate datasets for Genes and Non-genes
        genes_dataset = {}
        non_genes_dataset = {}

        for s_name in samples_with_molecules:
            data = transcript_data_by_sample[s_name]
            molecules_data = data["molecules_per_gene"]

            # Separate this sample's counts by gene type
            sample_gene_counts = []
            sample_non_gene_counts = []

            for gene_name, gene_info in molecules_data.items():
                count = gene_info["count"]
                if count > 0:
                    if gene_info["is_gene"]:
                        sample_gene_counts.append(count)
                    else:
                        sample_non_gene_counts.append(count)

            # Create histograms for this sample
            if sample_gene_counts:
                gene_hist, _ = np.histogram(sample_gene_counts, bins=bins)
                gene_data = {}
                for i, count in enumerate(gene_hist):
                    gene_data[float(bin_centers[i])] = int(count)
                genes_dataset[s_name] = gene_data

            if sample_non_gene_counts:
                non_gene_hist, _ = np.histogram(sample_non_gene_counts, bins=bins)
                non_gene_data = {}
                for i, count in enumerate(non_gene_hist):
                    non_gene_data[float(bin_centers[i])] = int(count)
                non_genes_dataset[s_name] = non_gene_data

        # Create datasets list for multiple tabs
        datasets = []
        data_labels = []

        if genes_dataset:
            datasets.append(genes_dataset)
            data_labels.append({"name": "Genes", "ylab": "Number of genes"})

        if non_genes_dataset:
            datasets.append(non_genes_dataset)
            data_labels.append({"name": "Non-genes", "ylab": "Number of non-gene features"})

        if not datasets:
            return None

        config = {
            "id": "xenium_molecules_per_gene",
            "title": "Xenium: Distribution of Molecules per Gene",
            "xlab": "Number of molecules per gene",
            "ylab": "Number of genes",
            "data_labels": data_labels,
        }

        # Add vertical line for noise threshold if calculated
        if n_mols_threshold is not None and n_mols_threshold > 0:
            config["x_lines"] = [
                {
                    "value": n_mols_threshold,
                    "color": "grey",
                    "dash": "dash",
                    "width": 1,
                    "label": f"Noise threshold ({n_mols_threshold:.0f})",
                }
            ]

        return linegraph.plot(datasets, config)

    def calculate_noise_threshold(self, gene_molecule_counts, quantile=0.99):
        """
        Calculate noise threshold based on negative control molecules.

        Args:
            gene_molecule_counts: Dict of {gene_name: {"count": int, "is_gene": bool}}
            quantile: Quantile for threshold calculation (default 0.99)

        Returns:
            Float threshold value or None if insufficient data
        """
        import numpy as np

        # Extract counts for negative control features (non-genes starting with "NegControl")
        neg_control_counts = []
        neg_control_found = 0
        for gene_name, gene_info in gene_molecule_counts.items():
            if not gene_info["is_gene"] and gene_name.startswith("NegControl"):
                neg_control_found += 1
                if gene_info["count"] > 0:
                    neg_control_counts.append(gene_info["count"])

        if len(neg_control_counts) < 3:  # Need at least 3 data points for meaningful statistics
            return None

        try:
            # Calculate threshold using log-space statistics (similar to notebook)
            log_counts = np.log10(neg_control_counts)

            # Use median absolute deviation as robust estimate of standard deviation
            median_log = np.median(log_counts)
            mad = np.median(np.abs(log_counts - median_log))
            # Convert MAD to standard deviation equivalent (normal distribution scaling factor)
            std_log = mad * 1.4826

            # Calculate upper bound using quantile
            if SCIPY_AVAILABLE:
                from scipy.stats import norm

            z_score = norm.ppf(quantile)
            threshold_log = median_log + z_score * std_log

            threshold = 10**threshold_log
            return threshold

        except (ImportError, ValueError):
            # Fallback to simple percentile if scipy not available
            return np.percentile(neg_control_counts, quantile * 100)

        except Exception:
            # Return None if calculation fails
            return None

    def xenium_cell_distributions_combined_plot(self, cells_data_by_sample):
        """Create combined plot for transcripts and detected genes per cell distributions"""
        # Check if we have data for either transcripts or genes
        samples_with_transcripts = {}
        samples_with_genes = {}

        for s_name, data in cells_data_by_sample.items():
            if data and "transcript_counts_values" in data and data["transcript_counts_values"]:
                samples_with_transcripts[s_name] = data["transcript_counts_values"]
            if data and "detected_genes_values" in data and data["detected_genes_values"]:
                samples_with_genes[s_name] = data["detected_genes_values"]

        # If neither dataset is available, return None
        if not samples_with_transcripts and not samples_with_genes:
            return None

        num_samples = max(len(samples_with_transcripts), len(samples_with_genes))

        if num_samples == 1:
            # Single sample: Create combined density plots
            return self._create_single_sample_combined_density(samples_with_transcripts, samples_with_genes)
        else:
            # Multiple samples: Create combined box plots
            return self._create_multi_sample_combined_boxes(samples_with_transcripts, samples_with_genes)

    def _create_single_sample_combined_density(self, samples_with_transcripts, samples_with_genes):
        """Create single sample combined density plots for transcripts and genes per cell"""
        plot_data = []
        data_labels = []

        # Handle transcripts per cell data
        if samples_with_transcripts:
            s_name, transcript_values = next(iter(samples_with_transcripts.items()))
            try:
                import numpy as np

                if SCIPY_AVAILABLE:
                    from scipy.stats import gaussian_kde

                transcript_values = np.array(transcript_values)
                kde = gaussian_kde(transcript_values)
                x_min, x_max = transcript_values.min(), transcript_values.max()
                x_range = np.linspace(x_min, x_max, 1000)
                density = kde(x_range)

                # Add to plot data with dataset identifier
                transcripts_data = {}
                for x, y in zip(x_range, density):
                    transcripts_data[float(x)] = float(y)
                plot_data.append({s_name: transcripts_data})
                data_labels.append({"name": "Transcripts per cell", "xlab": "Number of transcripts per cell"})

            except ImportError:
                # Fallback to histogram if scipy not available
                import numpy as np

                bins = min(50, len(transcript_values) // 20)
                hist, bin_edges = np.histogram(transcript_values, bins=bins)
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

                transcripts_data = {}
                for x, y in zip(bin_centers, hist):
                    transcripts_data[float(x)] = float(y)
                plot_data.append({s_name: transcripts_data})
                data_labels.append({"name": "Transcripts per cell", "xlab": "Number of transcripts per cell"})

        # Handle detected genes per cell data
        if samples_with_genes:
            s_name, gene_values = next(iter(samples_with_genes.items()))
            try:
                import numpy as np

                if SCIPY_AVAILABLE:
                    from scipy.stats import gaussian_kde

                gene_values = np.array(gene_values)
                kde = gaussian_kde(gene_values)
                x_min, x_max = gene_values.min(), gene_values.max()
                x_range = np.linspace(x_min, x_max, 1000)
                density = kde(x_range)

                # Add to plot data with dataset identifier
                genes_data = {}
                for x, y in zip(x_range, density):
                    genes_data[float(x)] = float(y)
                plot_data.append({s_name: genes_data})
                data_labels.append({"name": "Detected genes per cell", "xlab": "Number of detected genes per cell"})

            except ImportError:
                # Fallback to histogram if scipy not available
                import numpy as np

                bins = min(50, len(gene_values) // 20)
                hist, bin_edges = np.histogram(gene_values, bins=bins)
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

                genes_data = {}
                for x, y in zip(bin_centers, hist):
                    genes_data[float(x)] = float(y)
                plot_data.append({s_name: genes_data})
                data_labels.append({"name": "Detected genes per cell", "xlab": "Number of detected genes per cell"})

        config = {
            "id": "xenium_cell_distributions_combined",
            "title": "Xenium: Cell Distribution Analysis",
            "ylab": "Density",
            "smooth_points": 100,
            "data_labels": data_labels,
        }

        return linegraph.plot(plot_data, config)

    def _create_multi_sample_combined_boxes(self, samples_with_transcripts, samples_with_genes):
        """Create multi-sample combined box plots for transcripts and genes per cell"""
        from multiqc.plots import box

        plot_data = []
        data_labels = []

        # Add transcripts per cell data
        if samples_with_transcripts:
            transcripts_data = {}
            for s_name, transcript_values in samples_with_transcripts.items():
                transcripts_data[s_name] = transcript_values
            plot_data.append(transcripts_data)
            data_labels.append({"name": "Transcripts per cell", "ylab": "Number of transcripts per cell"})

        # Add detected genes per cell data
        if samples_with_genes:
            genes_data = {}
            for s_name, gene_values in samples_with_genes.items():
                genes_data[s_name] = gene_values
            plot_data.append(genes_data)
            data_labels.append({"name": "Detected genes per cell", "ylab": "Number of detected genes per cell"})

        config = {
            "id": "xenium_cell_distributions_combined",
            "title": "Xenium: Cell Distribution Analysis",
            "boxpoints": False,
            "data_labels": data_labels,
        }

        return box.plot(plot_data, config)

    def xenium_transcripts_per_cell_plot(self, cells_data_by_sample):
        """Create transcripts per cell distribution plot"""
        # Filter samples with transcript count data
        samples_with_transcripts = {}
        for s_name, data in cells_data_by_sample.items():
            if data and "transcript_counts_values" in data and data["transcript_counts_values"]:
                samples_with_transcripts[s_name] = data["transcript_counts_values"]

        if not samples_with_transcripts:
            return None

        num_samples = len(samples_with_transcripts)

        if num_samples == 1:
            # Single sample: Create density plot
            return self._create_single_sample_transcripts_density(samples_with_transcripts)
        else:
            # Multiple samples: Create box plots
            return self._create_multi_sample_transcripts_boxes(samples_with_transcripts)

    def _create_single_sample_transcripts_density(self, samples_with_transcripts):
        """Create single sample transcripts per cell density plot"""
        s_name, transcript_values = next(iter(samples_with_transcripts.items()))

        # Create kernel density estimation
        try:
            import numpy as np

            if SCIPY_AVAILABLE:
                from scipy.stats import gaussian_kde

            transcript_values = np.array(transcript_values)
            kde = gaussian_kde(transcript_values)

            # Create x range for density plot
            x_min, x_max = transcript_values.min(), transcript_values.max()
            x_range = np.linspace(x_min, x_max, 1000)
            density = kde(x_range)

            # Create plot data
            plot_data = {s_name: {}}
            for x, y in zip(x_range, density):
                plot_data[s_name][float(x)] = float(y)

            config = {
                "id": "xenium_transcripts_per_cell_single",
                "title": "Xenium: Distribution of Transcripts per Cell",
                "xlab": "Number of transcripts per cell",
                "ylab": "Density",
                "smooth_points": 100,
            }

            return linegraph.plot(plot_data, config)

        except ImportError:
            # Fallback to histogram if scipy not available
            import numpy as np

            bins = min(50, len(transcript_values) // 20)
            hist, bin_edges = np.histogram(transcript_values, bins=bins)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            plot_data = {s_name: {}}
            for x, y in zip(bin_centers, hist):
                plot_data[s_name][float(x)] = float(y)

            config = {
                "id": "xenium_transcripts_per_cell_single",
                "title": "Xenium: Distribution of Transcripts per Cell",
                "xlab": "Number of transcripts per cell",
                "ylab": "Number of cells",
            }

            return linegraph.plot(plot_data, config)

    def _create_multi_sample_transcripts_boxes(self, samples_with_transcripts):
        """Create multi-sample transcripts per cell box plots"""
        from multiqc.plots import box

        # Prepare data for box plot
        plot_data = {}
        for s_name, transcript_values in samples_with_transcripts.items():
            plot_data[s_name] = transcript_values

        config = {
            "id": "xenium_transcripts_per_cell_multi",
            "title": "Xenium: Distribution of Transcripts per Cell",
            "ylab": "Number of transcripts per cell",
            "boxpoints": False,
        }

        return box.plot(plot_data, config)

    def xenium_detected_genes_per_cell_plot(self, cells_data_by_sample):
        """Create detected genes per cell distribution plot"""
        # Filter samples with detected genes data
        samples_with_genes = {}
        for s_name, data in cells_data_by_sample.items():
            if data and "detected_genes_values" in data and data["detected_genes_values"]:
                samples_with_genes[s_name] = data["detected_genes_values"]

        if not samples_with_genes:
            return None

        num_samples = len(samples_with_genes)

        if num_samples == 1:
            # Single sample: Create density plot
            return self._create_single_sample_genes_density(samples_with_genes)
        else:
            # Multiple samples: Create box plots
            return self._create_multi_sample_genes_boxes(samples_with_genes)

    def _create_single_sample_genes_density(self, samples_with_genes):
        """Create single sample detected genes per cell density plot"""
        s_name, gene_values = next(iter(samples_with_genes.items()))

        # Create kernel density estimation
        try:
            import numpy as np

            if SCIPY_AVAILABLE:
                from scipy.stats import gaussian_kde

            gene_values = np.array(gene_values)
            kde = gaussian_kde(gene_values)

            # Create x range for density plot
            x_min, x_max = gene_values.min(), gene_values.max()
            x_range = np.linspace(x_min, x_max, 1000)
            density = kde(x_range)

            # Create plot data
            plot_data = {s_name: {}}
            for x, y in zip(x_range, density):
                plot_data[s_name][float(x)] = float(y)

            config = {
                "id": "xenium_detected_genes_per_cell_single",
                "title": "Xenium: Distribution of Detected Genes per Cell",
                "xlab": "Number of detected genes per cell",
                "ylab": "Density",
                "smooth_points": 100,
            }

            return linegraph.plot(plot_data, config)

        except ImportError:
            # Fallback to histogram if scipy not available
            import numpy as np

            bins = min(50, len(gene_values) // 20)
            hist, bin_edges = np.histogram(gene_values, bins=bins)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            plot_data = {s_name: {}}
            for x, y in zip(bin_centers, hist):
                plot_data[s_name][float(x)] = float(y)

            config = {
                "id": "xenium_detected_genes_per_cell_single",
                "title": "Xenium: Distribution of Detected Genes per Cell",
                "xlab": "Number of detected genes per cell",
                "ylab": "Number of cells",
            }

            return linegraph.plot(plot_data, config)

    def _create_multi_sample_genes_boxes(self, samples_with_genes):
        """Create multi-sample detected genes per cell box plots"""
        from multiqc.plots import box

        # Prepare data for box plot
        plot_data = {}
        for s_name, gene_values in samples_with_genes.items():
            plot_data[s_name] = gene_values

        config = {
            "id": "xenium_detected_genes_per_cell_multi",
            "title": "Xenium: Distribution of Detected Genes per Cell",
            "ylab": "Number of detected genes per cell",
            "boxpoints": False,
        }

        return box.plot(plot_data, config)
