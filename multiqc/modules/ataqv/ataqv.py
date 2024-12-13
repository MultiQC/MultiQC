import json
import logging
from typing import Dict, Any, Optional

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph, scatter
from multiqc.plots.linegraph import LinePlotConfig

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="ATAQV",
            anchor="ataqv",
            href="https://github.com/ParkerLab/ataqv/",
            info="Toolkit for quality control and visualization of ATAC-seq data",
            doi="10.1093/bioinformatics/btx865",
        )

        # Find and parse JSON files
        self.ataqv_data = dict()

        for f in self.find_log_files("ataqv"):
            parsed_data = self.parse_ataqv_json(f)
            if parsed_data is not None:
                self.ataqv_data.update(parsed_data)

        # Filter to strip out ignored sample names
        self.ataqv_data = self.ignore_samples(self.ataqv_data)
        if len(self.ataqv_data) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.ataqv_data)} reports")

        self.general_stats_table()

        try:
            self.add_fragment_length_distribution_plot()
        except Exception as e:
            log.warning(f"Error adding fragment length distribution plot: {e}")

        try:
            self.add_peak_percentile_plot()
        except Exception as e:
            log.warning(f"Error adding peak percentile plot: {e}")

        try:
            self.add_mapq_distribution_plot()
        except Exception as e:
            log.warning(f"Error adding MAPQ distribution plot: {e}")

        try:
            self.add_chromosome_distribution_plot()
        except Exception as e:
            log.warning(f"Error adding chromosome distribution plot: {e}")

        try:
            self.add_fld_distance_plot()
        except Exception as e:
            log.warning(f"Error adding FLD distance plot: {e}")

        # Write data to file
        self.write_data_file(self.ataqv_data, "multiqc_ataqv")

    def parse_ataqv_json(self, f) -> Optional[Dict]:
        """Parse ataqv JSON report"""
        try:
            content = json.loads(f["f"])
            if not isinstance(content, list) or len(content) == 0:
                return None

            parsed_data = dict()

            # Get the metrics from the first record
            metrics = content[0]["metrics"]
            s_name = self.clean_s_name(metrics["name"], f)

            if s_name in parsed_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")

            self.add_data_source(f, s_name)

            self.add_software_version(
                version=content[0]["ataqv_version"],
                sample=s_name,
                software_name="ataqv",
            )

            # Extract key metrics
            parsed_data[s_name] = {
                "tss_enrichment": metrics["tss_enrichment"],
                "hqaa_percent": (metrics["hqaa"] / metrics["paired_reads"] * 100) if metrics["paired_reads"] > 0 else 0,
                "duplicate_fraction": metrics["duplicate_fraction_in_peaks"],
                "total_reads": metrics["paired_reads"],
                "hqaa": metrics["hqaa"],
                "mean_mapq": metrics["mean_mapq"],
                "median_mapq": metrics["median_mapq"],
                "hqaa_mononucleosomal_count": metrics["hqaa_mononucleosomal_count"],
                "hqaa_tf_count": metrics["hqaa_tf_count"],
                "max_fraction_reads_from_single_autosome": metrics["max_fraction_reads_from_single_autosome"],
                "duplicate_fraction_not_in_peaks": metrics["duplicate_fraction_not_in_peaks"],
                "peak_percentiles": {
                    "hqaa": metrics["peak_percentiles"]["cumulative_fraction_of_hqaa"],
                    "territory": metrics["peak_percentiles"]["cumulative_fraction_of_territory"],
                },
                "mapq_distribution": {int(entry[0]): int(entry[1]) for entry in metrics["mapq_counts"]},
                "chromosome_counts": {entry[0]: int(entry[1]) for entry in metrics["chromosome_counts"]},
                "fragment_length_counts": {int(entry[0]): int(entry[1]) for entry in metrics["fragment_length_counts"]},
            }

            return parsed_data

        except Exception as e:
            log.error(f"Error parsing ataqv JSON: {e}")
            return None

    def general_stats_table(self):
        """Add key metrics to the General Stats table"""

        self.general_stats_addcols(
            self.ataqv_data,
            {
                "tss_enrichment": {
                    "title": "TSS Enrichment",
                    "description": "Transcription start site enrichment score",
                    "min": 0,
                    "scale": "YlOrRd",
                    "format": "{:,.2f}",
                },
                "hqaa_percent": {
                    "title": "HQAA in peaks",
                    "description": "Percentage of high-quality autosomal alignments overlapping peaks",
                    "max": 100,
                    "min": 0,
                    "suffix": "%",
                    "scale": "YlGn",
                    "format": "{:,.1f}",
                },
                "duplicate_fraction": {
                    "title": "Dup Rate",
                    "description": "Duplicate read fraction in peaks",
                    "max": 1,
                    "min": 0,
                    "scale": "RdYlBu-rev",
                    "format": "{:,.3f}",
                },
                # Hidden by default - additional metrics
                "total_reads": {
                    "title": "Total Reads",
                    "description": "Total number of paired reads",
                    "min": 0,
                    "scale": "Blues",
                    "shared_key": "read_count",
                    "hidden": True,
                },
                "hqaa": {
                    "title": "HQAA",
                    "description": "Number of high-quality autosomal alignments",
                    "min": 0,
                    "scale": "Purples",
                    "shared_key": "read_count",
                    "hidden": True,
                },
                "mean_mapq": {
                    "title": "Mean MAPQ",
                    "description": "Mean mapping quality score",
                    "min": 0,
                    "max": 60,
                    "scale": "YlOrRd",
                    "format": "{:,.1f}",
                    "hidden": True,
                },
                "median_mapq": {
                    "title": "Median MAPQ",
                    "description": "Median mapping quality score",
                    "min": 0,
                    "max": 60,
                    "scale": "YlOrRd",
                    "format": "{:,.1f}",
                    "hidden": True,
                },
                "hqaa_mononucleosomal_count": {
                    "title": "Mono-nuc",
                    "description": "Number of high-quality autosomal alignments in mononucleosomal length range",
                    "min": 0,
                    "scale": "Greens",
                    "shared_key": "read_count",
                    "hidden": True,
                },
                "hqaa_tf_count": {
                    "title": "TF Footprint",
                    "description": "Number of high-quality autosomal alignments in transcription factor footprint length range",
                    "min": 0,
                    "scale": "Oranges",
                    "shared_key": "read_count",
                    "hidden": True,
                },
                "max_fraction_reads_from_single_autosome": {
                    "title": "Max Chr%",
                    "description": "Maximum fraction of reads from a single autosome",
                    "min": 0,
                    "max": 1,
                    "scale": "RdYlBu-rev",
                    "format": "{:,.3f}",
                    "hidden": True,
                },
                "duplicate_fraction_not_in_peaks": {
                    "title": "Dup Rate (non-peak)",
                    "description": "Duplicate read fraction outside of peaks",
                    "max": 1,
                    "min": 0,
                    "scale": "RdYlBu-rev",
                    "format": "{:,.3f}",
                    "hidden": True,
                },
            },
        )

    def add_fragment_length_distribution_plot(self):
        """Add fragment length distribution plot"""

        counts = {}
        for s_name, d in self.ataqv_data.items():
            if "fragment_length_counts" in d:
                counts[s_name] = d["fragment_length_counts"]

        percentages = {}
        for s_name, d in self.ataqv_data.items():
            if "fragment_length_counts" in d:
                percentages[s_name] = {
                    k: v / sum(d["fragment_length_counts"].values()) * 100
                    for k, v in d["fragment_length_counts"].items()
                }

        if len(counts) > 0:
            self.add_section(
                name="Fragment Length Distribution",
                anchor="ataqv-fragment-length-dist",
                description="Distribution of fragment lengths from ATAC-seq reads.",
                helptext="""
                Fragment length distributions can help assess ATAC-seq data quality:
                * A strong nucleosome-free peak should be visible ~50-100bp
                * Mono-nucleosome peak should be visible ~180-247bp
                * Di-nucleosome peak may be visible ~350-430bp
                """,
                plot=linegraph.plot(
                    [counts, percentages],
                    pconfig=LinePlotConfig(
                        id="ataqv_fld_plot",
                        title="ataqv: Fragment Length Distribution",
                        xlab="Fragment Length (bp)",
                        ylab="Read Count",
                        logswitch=True,
                        logswitch_active=True,
                        data_labels=[
                            {"name": "Read Count"},
                            {"name": "Percentage", "ysuffix": "%"},
                        ],
                    ),
                ),
            )

    def add_peak_percentile_plot(self):
        """Plot cumulative distribution of reads in peaks"""

        data_hqaa = {}
        data_territory = {}
        for s_name, d in self.ataqv_data.items():
            if "peak_percentiles" in d:
                data_hqaa[s_name] = d["peak_percentiles"]["hqaa"]
                data_territory[s_name] = d["peak_percentiles"]["territory"]

        if len(data_hqaa) > 0 or len(data_territory) > 0:
            self.add_section(
                name="Peak percentiles",
                anchor="ataqv-peak-percentiles",
                description="Cumulative distribution of reads in peaks and of genomic territory.",
                helptext="""
                This plot shows the cumulative fraction of high-quality autosomal alignments, 
                and the cumulative fraction of the genome territory covered by those peaks.

                * A steeper curve indicates better enrichment
                * More vertical early rise suggests stronger peak signals
                """,
                plot=linegraph.plot(
                    [data_hqaa, data_territory],
                    pconfig=LinePlotConfig(
                        id="ataqv_peak_percentiles",
                        title="ataqv: Peak Enrichment Distribution",
                        xlab="Peak percentile",
                        data_labels=[
                            {"name": "HQAA", "ylab": "Cumulative fraction of HQAA"},
                            {"name": "Territory", "ylab": "Cumulative fraction of territory"},
                        ],
                    ),
                ),
            )

    def add_mapq_distribution_plot(self):
        """Plot distribution of mapping quality scores"""

        plot_data = {}
        for s_name, d in self.ataqv_data.items():
            if "mapq_distribution" in d:
                plot_data[s_name] = d["mapq_distribution"]

        if len(plot_data) > 0:
            self.add_section(
                name="Mapping Quality Distribution",
                anchor="ataqv-mapq-dist",
                description="Distribution of read mapping quality scores.",
                helptext="""
                The mapping quality (MAPQ) distribution shows the confidence of read alignments:
                * Higher MAPQ scores indicate more unique, confident alignments
                * MAPQ=0 often indicates multi-mapped reads
                * Most high-quality alignments should have high MAPQ scores
                """,
                plot=linegraph.plot(
                    plot_data,
                    pconfig=LinePlotConfig(
                        id="ataqv_mapq_dist",
                        title="ataqv: Mapping Quality Scores",
                        xlab="MAPQ Score",
                        ylab="Number of Reads",
                        logswitch=True,
                        logswitch_active=True,
                    ),
                ),
            )

    def add_chromosome_distribution_plot(self):
        """Plot distribution of reads across chromosomes"""
        counts_data = {}
        percentages_data = {}
        for s_name, d in self.ataqv_data.items():
            if "chromosome_counts" not in d:
                continue

            # Sort chromosomes naturally (1,2,3...10,11 instead of 1,10,11,2,3...)
            chrom_counts = d["chromosome_counts"]
            sorted_chroms = sorted(
                chrom_counts.keys(),
                key=lambda x: (
                    # Sort chrM to the end
                    x == "chrM",
                    # Sort chrX, chrY after numbers
                    x in ("chrX", "chrY"),
                    # Extract number from chromosome name and convert to int for numerical sorting
                    int("".join(filter(str.isdigit, x))) if any(c.isdigit() for c in x) else float("inf"),
                    # Use original string as final tie-breaker
                    x,
                ),
            )
            counts_data[s_name] = [(chrom, chrom_counts[chrom]) for chrom in sorted_chroms]

            # Calculate percentages
            if (total := sum(chrom_counts.values())) > 0:
                percentages_data[s_name] = [(chrom, (chrom_counts[chrom] / total) * 100) for chrom in sorted_chroms]

        if len(counts_data) > 0:
            self.add_section(
                name="Chromosome Distribution",
                anchor="ataqv-chrom-dist",
                description="Distribution of reads across chromosomes.",
                helptext="""
                The proportion of reads mapping to each chromosome:
                * Mitochondrial reads (chrM) often indicate poor nuclear enrichment
                * Large deviations in chromosomal proportions may indicate copy number variations
                * Expected proportions should roughly match chromosome sizes
                """,
                plot=linegraph.plot(
                    [counts_data, percentages_data],
                    pconfig=LinePlotConfig(
                        id="ataqv_chrom_dist",
                        title="ataqv: Read Distribution by Chromosome",
                        xlab="Chromosome",
                        ylab="Reads",
                        ymin=0,
                        categories=True,  # Treat x values as categories
                        data_labels=[
                            {"name": "Reads"},
                            {"name": "Percentages", "ysuffix": "%"},
                        ],
                    ),
                ),
            )

    def add_fld_distance_plot(self):
        """Plot distance from reference fragment length distribution vs HQAA percentage"""

        # Calculate distances from reference FLD for each sample
        distances_data = {}
        for s_name, d in self.ataqv_data.items():
            if "fragment_length_counts" not in d and "fragment_length_distance" not in d:
                continue

            # Calculate HQAA percentage
            hqaa_pct = d["hqaa_percent"]
            fld_distance = d["fragment_length_distance"]

            # Here you would calculate the distance from reference FLD
            # This is a placeholder - implement actual distance calculation
            fld_distance = 0.0  # Replace with actual calculation

            distances_data[s_name] = {"x": fld_distance, "y": hqaa_pct}

        if len(distances_data) > 0:
            self.add_section(
                name="Distance from Reference FLD",
                anchor="ataqv-fld-distance",
                description="Comparison of fragment length distributions to a reference profile.",
                helptext="""
                This plot shows how each sample's fragment length distribution compares to a reference:
                * X-axis shows the distance from the reference distribution
                * Y-axis shows the percentage of high-quality autosomal alignments
                * Samples clustering together have similar fragment length profiles
                """,
                plot=scatter.plot(
                    distances_data,
                    pconfig={
                        "id": "ataqv_fld_distance",
                        "title": "ataqv: Distance from Reference FLD",
                        "xlab": "Distance from reference fragment length distribution",
                        "ylab": "High-quality autosomal (%)",
                        "ymin": 0,
                        "ymax": 100,
                    },
                ),
            )
