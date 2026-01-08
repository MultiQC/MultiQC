import ast
import logging
from typing import Dict

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, heatmap, linegraph
from multiqc.plots.table_object import ColumnDict

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    Ribo-TISH is a tool for identifying translated ORFs from Ribo-seq data.
    This module parses the `*_qual.txt` output files to visualize reading frame
    quality metrics across different read lengths.

    The module creates one of two visualizations:
    1. A stacked bar chart showing the proportion of reads in each reading frame
       (Frame 0, 1, 2) for read lengths 25-34nt
    2. A heatmap showing the percentage distribution of read lengths within each sample
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Ribo-TISH",
            anchor="ribotish",
            href="https://github.com/zhpn1024/ribotish",
            info="Identifies translated ORFs from Ribo-seq data and reports reading frame quality metrics.",
            doi="10.1186/s13059-017-1316-1",
        )

        # Store parsed data
        self.ribotish_data: Dict[str, Dict] = {}
        self.frame_proportions: Dict[str, Dict] = {}

        # Parse all *_qual.txt files
        for f in self.find_log_files("ribotish/qual"):
            parsed_data = self.parse_ribotish_qual(f)
            if parsed_data:
                sample_name = f["s_name"]
                self.ribotish_data[sample_name] = parsed_data
                self.add_data_source(f)

        # Filter to strip out ignored sample names
        self.ribotish_data = self.ignore_samples(self.ribotish_data)

        if not self.ribotish_data:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.ribotish_data)} Ribo-TISH reports")

        # Calculate frame proportions for all samples
        self.calculate_frame_proportions()

        # Call add_software_version even though version is not available
        self.add_software_version(None)

        # Add sections with plots
        self.add_frame_proportion_bargraph()
        self.add_read_length_distribution()
        self.add_general_stats()

        # Write data file at the end
        self.write_data_file(self.ribotish_data, "multiqc_ribotish")

    def parse_ribotish_qual(self, f) -> Dict:
        """
        Parse Ribo-TISH *_qual.txt file.

        The relevant information is on line 4 (0-indexed line 3) in Python dict format:
        {25: [frame0_count, frame1_count, frame2_count], 26: [...], ...}

        See: https://github.com/zhpn1024/ribotish#offset-parameter-file

        Returns a dict with read lengths as keys and [f0, f1, f2] counts as values.
        """
        try:
            lines = f["f"].splitlines()
            if len(lines) < 4:
                log.warning(f"File {f['fn']} has fewer than 4 lines, skipping")
                return {}

            input_string = lines[3].strip()

            # Remove optional $ prefix if present
            if input_string.startswith("$"):
                input_string = input_string[1:].strip()

            # Parse Python dict literal safely with ast.literal_eval
            ribo_seq_counts = ast.literal_eval(input_string)

            # Validate the parsed data structure
            if not isinstance(ribo_seq_counts, dict):
                log.warning(f"File {f['fn']}: expected dict, got {type(ribo_seq_counts).__name__}")
                return {}

            # Validate that all values are lists with exactly 3 elements
            for read_length, counts in ribo_seq_counts.items():
                if not isinstance(counts, list):
                    log.warning(
                        f"File {f['fn']}: expected list for read length {read_length}, got {type(counts).__name__}"
                    )
                    return {}
                if len(counts) != 3:
                    log.warning(
                        f"File {f['fn']}: expected 3 frame counts for read length {read_length}, got {len(counts)}"
                    )
                    return {}
                # Validate that all counts are numeric
                if not all(isinstance(c, (int, float)) for c in counts):
                    log.warning(f"File {f['fn']}: non-numeric counts found for read length {read_length}")
                    return {}

            log.debug(f"Successfully parsed {f['fn']} with {len(ribo_seq_counts)} read lengths")
            return ribo_seq_counts

        except (SyntaxError, ValueError, IndexError) as e:
            log.warning(f"Error parsing {f['fn']}: {e}")
            return {}

    def calculate_frame_proportions(self):
        """Calculate frame proportions for each sample and read length."""
        for sample_name, ribo_seq_counts in self.ribotish_data.items():
            self.frame_proportions[sample_name] = {}

            for length, counts in ribo_seq_counts.items():
                total = sum(counts)
                if total > 0:
                    self.frame_proportions[sample_name][length] = {
                        "f0_prop": counts[0] / float(total),
                        "f1_prop": counts[1] / float(total),
                        "f2_prop": counts[2] / float(total),
                        "total": total,
                    }
                else:
                    self.frame_proportions[sample_name][length] = {
                        "f0_prop": 0.0,
                        "f1_prop": 0.0,
                        "f2_prop": 0.0,
                        "total": 0.0,
                    }

    def add_frame_proportion_bargraph(self):
        """
        Create a stacked bar graph showing frame proportions for all read lengths.
        All read lengths are shown side-by-side in a single plot for visual comparison.
        Each sample-length combination gets its own bar.
        """
        # First, collect all read lengths across all samples
        all_lengths_set: set[int] = set()
        for sample_data in self.frame_proportions.values():
            all_lengths_set.update(sample_data.keys())
        all_lengths = sorted(all_lengths_set)

        plot_data = {}
        sample_groups: Dict[str, list] = {}

        for length in all_lengths:
            group_name = f"{length}nt"
            length_group = []
            for sample_name in sorted(self.frame_proportions.keys()):
                if length in self.frame_proportions[sample_name]:
                    props = self.frame_proportions[sample_name][length]
                    sample_key = f"{sample_name}_{length}nt"
                    plot_data[sample_key] = {
                        "Frame 0": props["f0_prop"] * 100.0,
                        "Frame 1": props["f1_prop"] * 100.0,
                        "Frame 2": props["f2_prop"] * 100.0,
                    }
                    length_group.append([sample_key, sample_name])
            if length_group:
                sample_groups[group_name] = length_group

        pconfig = {
            "id": "ribotish_frame_proportions",
            "title": "Ribo-TISH: Reading Frame Proportions by Read Length",
            "ylab": "Proportion of Reads (%)",
            "stacking": "normal",
            "hide_zero_cats": False,
            "ymax": 100,
            "use_legend": True,
            "cpswitch": False,
            "sample_groups": sample_groups,
            "x_lines": [
                {
                    "color": "#0000ff",
                    "dash": "dash",
                    "value": 33.33,
                    "width": 2,
                    "label": "Random distribution (1/3)",
                },
                {
                    "color": "#0000ff",
                    "dash": "dash",
                    "value": 66.67,
                    "width": 2,
                    "label": "Random distribution (2/3)",
                },
            ],
        }

        # Define categories (let MultiQC assign default colors)
        cats = {
            "Frame 0": {"name": "Frame 0"},
            "Frame 1": {"name": "Frame 1"},
            "Frame 2": {"name": "Frame 2"},
        }

        plot_html = bargraph.plot(plot_data, cats, pconfig)

        self.add_section(
            name="Reading Frame Proportions",
            anchor="ribotish_frame_proportions_section",
            description="Proportion of reads in each reading frame (Frame 0, 1, 2) for different read lengths (25-34nt). "
            "Frame assignment is based on P-site positions as determined by Ribo-TISH. "
            "Some degree of frame preference (enrichment in Frame 0) is typically expected in Ribo-seq data.",
            helptext="""
            This plot shows the distribution of reads across the three reading frames for each read length,
            based on P-site positions.

            * **Frame 0**: The primary reading frame
            * **Frame 1**: Offset by 1 nucleotide from Frame 0
            * **Frame 2**: Offset by 2 nucleotides from Frame 0

            Ribo-seq data typically shows enrichment in Frame 0, particularly for read lengths 28-30nt,
            though the degree of frame preference can vary depending on the experimental protocol
            (e.g., RNase I vs. MNase treatment).
            Read lengths are shown as separate bars for each sample.
            """,
            plot=plot_html,
        )

    def add_read_length_distribution(self):
        """
        Create plots showing the percentage distribution of read lengths.
        Provides both line graph and heatmap views with a switcher.
        """
        # Collect all read lengths
        all_lengths_set: set[int] = set()
        for sample_data in self.frame_proportions.values():
            all_lengths_set.update(sample_data.keys())
        all_lengths = sorted(all_lengths_set)

        # Calculate sample totals once for efficiency
        sample_totals = {}
        for sample_name in self.frame_proportions.keys():
            sample_totals[sample_name] = sum(props["total"] for props in self.frame_proportions[sample_name].values())

        if len(sample_totals) <= 30:
            # Prepare data for line graph
            line_data: dict[str, dict[int, float]] = {}
            for sample_name in sorted(self.frame_proportions.keys()):
                sample_total = sample_totals[sample_name]
                line_data[sample_name] = {}
                for length in all_lengths:
                    if length in self.frame_proportions[sample_name]:
                        count = self.frame_proportions[sample_name][length]["total"]
                        percentage = (count / sample_total * 100.0) if sample_total > 0 else 0
                        line_data[sample_name][length] = percentage
                    else:
                        line_data[sample_name][length] = 0

            # Create line graph plot
            line_pconfig = {
                "id": "ribotish_read_length_line",
                "title": "Ribo-TISH: Read Length Distribution",
                "xlab": "Read Length (nt)",
                "ylab": "% of Total Reads",
                "smooth_points": 50,
                "smooth_points_sumcounts": False,
                "tt_label": "{point.x}nt: {point.y:.1f}%",
            }
            line_plot = linegraph.plot(line_data, line_pconfig)

            # Combine plots with data_labels for switching
            self.add_section(
                name="Read Length Distribution",
                anchor="ribotish_read_length_dist_section",
                description="Percentage of reads at each read length for each sample. "
                "Ribo-seq data typically shows enrichment around 28-30nt, representing ribosome-protected fragments.",
                helptext="""
                This plot shows what percentage of total reads each read length represents for each sample.

                * Each line represents a different sample
                * Peaks indicate the most common read lengths
                * Multiple samples can be easily compared

                The expected read length distribution can vary depending on the experimental protocol and organism.
                """,
                plot=line_plot,
            )

        else:
            # Prepare data for heatmap
            samples = sorted(self.frame_proportions.keys())
            heatmap_data = []
            for sample_name in samples:
                sample_total = sample_totals[sample_name]
                row = []
                for length in all_lengths:
                    if length in self.frame_proportions[sample_name]:
                        count = self.frame_proportions[sample_name][length]["total"]
                        percentage = (count / sample_total * 100.0) if sample_total > 0 else 0
                        row.append(percentage)
                    else:
                        row.append(0)
                heatmap_data.append(row)

            # Create heatmap plot
            heatmap_pconfig = {
                "id": "ribotish_read_length_heatmap",
                "title": "Ribo-TISH: Read Length Distribution (Heatmap)",
                "xlab": "Read Length (nt)",
                "ylab": "Sample",
                "zlab": "% of Total Reads",
                "square": False,
                "tt_decimals": 1,
                "legend": True,
                "xcats_samples": False,
                "ycats_samples": False,
                "cluster_rows": False,
                "cluster_cols": False,
                "display_values": False,
                "colstops": [[0, "#ffffff"], [0.5, "#4575b4"], [1, "#313695"]],
            }
            xcats = [f"{length}nt" for length in all_lengths]
            ycats = samples
            heatmap_plot = heatmap.plot(heatmap_data, xcats=xcats, ycats=ycats, pconfig=heatmap_pconfig)

            # Add heatmap as a second section for alternative view
            self.add_section(
                name="Read Length Distribution (Heatmap)",
                anchor="ribotish_read_length_heatmap_section",
                description="Alternative heatmap view of read length distribution. "
                "Useful for comparing many samples at once.",
                helptext="""
                This heatmap shows what percentage of total reads each read length represents for each sample.

                * Rows are samples, columns are read lengths
                * Darker blue indicates higher percentage of reads
                * Lighter colors indicate fewer reads

                This view is particularly useful when comparing many samples simultaneously.
                """,
                plot=heatmap_plot,
            )

    def add_general_stats(self):
        """Add key metrics to the general statistics table."""
        headers = {}

        # Weighted average Frame 0 proportion (most important metric)
        headers["weighted_f0_prop"] = ColumnDict(
            {
                "title": "Weighted F0 %",
                "description": "Frame 0 proportion weighted by read counts across all lengths - primary quality metric",
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.1f}",
                "min": 33,
                "max": 100,
            }
        )

        # Percentage of reads in optimal range
        headers["optimal_range_pct"] = ColumnDict(
            {
                "title": "% in 28-30nt",
                "description": "Percentage of reads in optimal ribosome footprint range (28-30nt)",
                "suffix": "%",
                "scale": "Blues",
                "format": "{:,.1f}",
                "max": 100,
            }
        )

        # Best Frame 0 proportion and length
        headers["best_f0_prop"] = ColumnDict(
            {
                "title": "Best F0 %",
                "description": "Highest Frame 0 proportion achieved at any read length",
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.1f}",
                "min": 33,
                "max": 100,
            }
        )

        headers["best_f0_length"] = ColumnDict(
            {
                "title": "Best F0 Length (nt)",
                "description": "Read length with highest Frame 0 proportion",
                "scale": "Greens",
                "format": "{:,.0f}",
            }
        )

        headers["peak_length"] = ColumnDict(
            {
                "title": "Peak Length (nt)",
                "description": "Read length with highest total read count",
                "scale": "Blues",
                "format": "{:,.0f}",
            }
        )

        headers["length_range"] = ColumnDict(
            {
                "title": "Length Range (nt)",
                "description": "Range of read lengths detected (min-max)",
                "scale": "Greys",
            }
        )

        # Total reads (hidden by default)
        headers["total_reads"] = ColumnDict(
            {
                "title": "Total Reads",
                "description": f"Total number of reads across all lengths ({config.read_count_desc})",
                "scale": "Purples",
                "hidden": True,
                "shared_key": "read_count",
            }
        )

        # Calculate statistics for general stats
        stats_data = {}
        for sample_name, length_props in self.frame_proportions.items():
            # Skip samples with no data
            if not length_props:
                continue

            # Calculate total reads and weighted Frame 0
            total_reads = 0
            weighted_f0_sum = 0
            for length, props in length_props.items():
                count = props["total"]
                total_reads += count
                weighted_f0_sum += props["f0_prop"] * count

            weighted_f0_prop = (weighted_f0_sum / total_reads * 100.0) if total_reads > 0 else 0

            # Calculate percentage of reads in optimal range (28-30nt)
            optimal_lengths = [28, 29, 30]
            optimal_reads = 0
            for length in optimal_lengths:
                if length in length_props:
                    optimal_reads += length_props[length]["total"]

            optimal_range_pct = (optimal_reads / total_reads * 100.0) if total_reads > 0 else 0

            # Find length with best frame 0 proportion
            best_f0_length = max(length_props.keys(), key=lambda k: length_props[k]["f0_prop"])
            best_f0_prop = length_props[best_f0_length]["f0_prop"] * 100.0

            # Find length with highest total count
            peak_length = max(length_props.keys(), key=lambda k: length_props[k]["total"])

            # Calculate read length range
            min_length = min(length_props.keys())
            max_length = max(length_props.keys())
            length_range = f"{min_length}-{max_length}"

            stats_data[sample_name] = {
                "weighted_f0_prop": weighted_f0_prop,
                "optimal_range_pct": optimal_range_pct,
                "best_f0_prop": best_f0_prop,
                "best_f0_length": best_f0_length,
                "peak_length": peak_length,
                "length_range": length_range,
                "total_reads": total_reads,
            }

        self.general_stats_addcols(stats_data, headers)
