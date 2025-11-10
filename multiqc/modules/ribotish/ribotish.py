import ast
import logging
import re
from typing import Dict, List

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, heatmap
from multiqc.plots.table_object import ColumnDict

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    RiboTish is a tool for identifying translated ORFs from Ribo-seq data.
    This module parses the `*_qual.txt` output files to visualize reading frame
    quality metrics across different read lengths.

    The module creates two visualizations:
    1. A stacked bar chart showing the proportion of reads in each reading frame
       (Frame 0, 1, 2) for read lengths 25-34nt
    2. A heatmap showing the percentage distribution of read lengths within each sample
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="RiboTish",
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

        # Calculate frame proportions for all samples
        self.calculate_frame_proportions()

        # Call add_software_version even though version is not available
        self.add_software_version(None)

        # Add sections with plots
        self.add_frame_proportion_bargraph()
        self.add_read_length_heatmap()
        self.add_general_stats()

        # Write data file at the end
        self.write_data_file(self.ribotish_data, "multiqc_ribotish")

    def parse_ribotish_qual(self, f) -> Dict:
        """
        Parse RiboTish *_qual.txt file.

        The relevant information is on line 4 (0-indexed line 3) in the format:
        $ {25: [frame0_count, frame1_count, frame2_count], 26: [...], ...}

        Returns a dict with read lengths as keys and [f0, f1, f2] counts as values.
        """
        try:
            lines = f["f"].splitlines()
            if len(lines) < 4:
                log.warning(f"File {f['fn']} has fewer than 4 lines, skipping")
                return {}

            input_string = lines[3].strip()

            # Clean up the input string
            input_string = input_string.strip()
            if input_string.startswith("$"):
                input_string = input_string[1:].strip()

            # Convert to valid Python dictionary format
            input_string = input_string.replace("{", "").replace("}", "")
            input_string = re.sub(r"(\d+):\s*\[", r'"\1": [', input_string)
            input_string = re.sub(r"\s+", " ", input_string)
            input_string = "{" + input_string.strip() + "}"

            # Parse the dictionary
            input_string = input_string.replace("'", '"')
            input_string = re.sub(r"(\d+),\s*(?=\d)", r"\1, ", input_string)
            ribo_seq_counts = ast.literal_eval(input_string)

            # Convert keys to integers for consistent handling
            return {int(k): v for k, v in ribo_seq_counts.items()}

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
                        "f0_prop": counts[0] / total,
                        "f1_prop": counts[1] / total,
                        "f2_prop": counts[2] / total,
                        "total": total,
                    }
                else:
                    self.frame_proportions[sample_name][length] = {
                        "f0_prop": 0,
                        "f1_prop": 0,
                        "f2_prop": 0,
                        "total": 0,
                    }

    def add_frame_proportion_bargraph(self):
        """
        Create a stacked bar graph showing frame proportions for all read lengths.
        Each sample gets multiple bars (one per read length).
        """
        # Prepare data for bargraph
        # Format: {sample_name: {category: value}}
        # We'll create categories like "25nt_Frame0", "25nt_Frame1", etc.

        # First, collect all read lengths across all samples
        all_lengths = set()
        for sample_data in self.frame_proportions.values():
            all_lengths.update(sample_data.keys())
        all_lengths = sorted(all_lengths)

        # Create separate plot data for each read length
        # Each dataset must have unique sample keys to avoid merging
        plot_data = []
        data_labels = []
        for length in all_lengths:
            length_data = {}
            for sample_name, length_props in self.frame_proportions.items():
                if length in length_props:
                    props = length_props[length]
                    # Add length suffix to make sample keys unique across datasets
                    # This prevents MultiQC from merging datasets
                    unique_sample_key = f"{sample_name}||{length}nt"
                    length_data[unique_sample_key] = {
                        "Frame 0": props["f0_prop"] * 100,  # Convert to percentage
                        "Frame 1": props["f1_prop"] * 100,
                        "Frame 2": props["f2_prop"] * 100,
                    }
            plot_data.append(length_data)
            data_labels.append({"name": f"{length}nt", "ylab": "Proportion of Reads (%)"})
            log.debug(f"Length {length}nt: {len(length_data)} samples")

        # Create configuration
        pconfig = {
            "id": "ribotish_frame_proportions",
            "title": "RiboTish: Reading Frame Proportions",
            "ylab": "Proportion of Reads (%)",
            "cpswitch_counts_label": "Counts",
            "stacking": "normal",
            "hide_zero_cats": False,
            "ymax": 100,
            "use_legend": True,
            "data_labels": data_labels,
        }

        # Define category colors (Frame 0, 1, 2) matching viridis-like colors
        # Need to provide the same categories for each dataset
        cats = {
            "Frame 0": {"name": "Frame 0", "color": "#440154"},
            "Frame 1": {"name": "Frame 1", "color": "#31688e"},
            "Frame 2": {"name": "Frame 2", "color": "#35b779"},
        }
        cats_list = [cats] * len(all_lengths)  # Same categories for each read length

        plot_html = bargraph.plot(plot_data, cats_list, pconfig)

        self.add_section(
            name="Reading Frame Proportions",
            anchor="ribotish_frame_proportions",
            description="Proportion of reads in each reading frame (Frame 0, 1, 2) for different read lengths (25-34nt). "
            "High Frame 0 enrichment (>70%) indicates good quality Ribo-seq data.",
            helptext="""
            This plot shows the distribution of reads across the three reading frames for each read length.

            * **Frame 0** (purple/dark): The correct reading frame - high values indicate good data quality
            * **Frame 1** (blue): Off by 1 nucleotide
            * **Frame 2** (green): Off by 2 nucleotides

            For high-quality Ribo-seq data, Frame 0 should typically be >70% for read lengths 28-30nt.
            Read lengths are shown as separate bars for each sample.
            """,
            plot=plot_html,
        )

    def add_read_length_heatmap(self):
        """
        Create a heatmap showing the percentage distribution of read lengths.
        Rows are samples, columns are read lengths.
        """
        # Collect all read lengths
        all_lengths = set()
        for sample_data in self.frame_proportions.values():
            all_lengths.update(sample_data.keys())
        all_lengths = sorted(all_lengths)

        # Prepare data for heatmap
        samples = sorted(self.frame_proportions.keys())
        heatmap_data = []

        for sample_name in samples:
            sample_total = sum(props["total"] for props in self.frame_proportions[sample_name].values())
            row = []
            for length in all_lengths:
                if length in self.frame_proportions[sample_name]:
                    count = self.frame_proportions[sample_name][length]["total"]
                    percentage = (count / sample_total * 100) if sample_total > 0 else 0
                    row.append(percentage)
                else:
                    row.append(0)
            heatmap_data.append(row)

        # Create heatmap configuration
        pconfig = {
            "id": "ribotish_read_length_dist",
            "title": "RiboTish: Read Length Distribution",
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
        }

        xcats = [f"{length}nt" for length in all_lengths]
        ycats = samples

        plot_html = heatmap.plot(heatmap_data, xcats=xcats, ycats=ycats, pconfig=pconfig)

        self.add_section(
            name="Read Length Distribution",
            anchor="ribotish_read_length_dist",
            description="Percentage of reads at each read length for each sample. "
            "Good Ribo-seq data typically shows enrichment at 28-30nt.",
            helptext="""
            This heatmap shows what percentage of total reads each read length represents for each sample.

            * **Darker blue** indicates higher percentage of reads
            * **Lighter colors** indicate fewer reads

            Typical Ribo-seq data shows strong enrichment around 28-30nt, which represents ribosome-protected fragments.
            """,
            plot=plot_html,
        )

    def add_general_stats(self):
        """Add key metrics to the general statistics table."""
        headers = {}

        # Add columns for best frame 0 proportion
        headers["best_f0_length"] = ColumnDict(
            {
                "title": "Best F0 Length",
                "description": "Read length with highest Frame 0 proportion",
                "suffix": "nt",
                "scale": "Greens",
                "format": "{:,.0f}",
            }
        )

        headers["best_f0_prop"] = ColumnDict(
            {
                "title": "Best F0 %",
                "description": "Highest Frame 0 proportion across all read lengths",
                "suffix": "%",
                "scale": "RdYlGn",
                "format": "{:,.1f}",
                "min": 33,
                "max": 100,
            }
        )

        headers["peak_length"] = ColumnDict(
            {
                "title": "Peak Length",
                "description": "Read length with highest total read count",
                "suffix": "nt",
                "scale": "Blues",
                "format": "{:,.0f}",
            }
        )

        headers["total_reads"] = ColumnDict(
            {
                "title": "Total Reads",
                "description": "Total number of reads across all lengths",
                "scale": "Purples",
                "format": "{:,.0f}",
                "hidden": True,
            }
        )

        # Calculate statistics for general stats
        stats_data = {}
        for sample_name, length_props in self.frame_proportions.items():
            # Find length with best frame 0 proportion
            best_f0_length = max(length_props.keys(), key=lambda k: length_props[k]["f0_prop"])
            best_f0_prop = length_props[best_f0_length]["f0_prop"] * 100

            # Find length with highest total count
            peak_length = max(length_props.keys(), key=lambda k: length_props[k]["total"])

            # Calculate total reads
            total_reads = sum(props["total"] for props in length_props.values())

            stats_data[sample_name] = {
                "best_f0_length": best_f0_length,
                "best_f0_prop": best_f0_prop,
                "peak_length": peak_length,
                "total_reads": total_reads,
            }

        self.general_stats_addcols(stats_data, headers)
