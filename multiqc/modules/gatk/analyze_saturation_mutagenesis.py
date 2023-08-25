# -*- coding: utf-8 -*-

""" MultiQC submodule to parse output from GATK tool AnalyzeSaturationMutagenesis """

import logging
from collections import OrderedDict

from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class AnalyzeSaturationMutagenesisMixin:
    """Mixin class to add GATK AnalyzeSaturationMutagenesis reports to MultiQC."""

    def parse_gatk_analyze_saturation_mutagenesis(self):
        """Find GATK AnalyzeSaturationMutagenesis logs and parse their data"""

        self.gatk_analyze_saturation_mutagenesis = {}

        for file_handle in self.find_log_files("gatk/analyze_saturation_mutagenesis", filehandles=True):
            parsed_data = self.parse_read_counts_file(file_handle["f"])
            if len(parsed_data) > 1:
                if file_handle["s_name"] in self.gatk_analyze_saturation_mutagenesis:
                    log.debug(f"Duplicate sample name found! Overwriting: {file_handle['s_name']}")
                self.add_data_source(file_handle, section="analyze_saturation_mutagenesis")
                self.gatk_analyze_saturation_mutagenesis[file_handle["s_name"]] = parsed_data

        # Filter to strip out ignored sample names
        self.gatk_analyze_saturation_mutagenesis = self.ignore_samples(self.gatk_analyze_saturation_mutagenesis)

        n_reports_found = len(self.gatk_analyze_saturation_mutagenesis)
        if n_reports_found == 0:
            return 0

        log.info("Found %d AnalyzeSaturationMutagenesis reports", n_reports_found)

        # Write parsed report data to a file (restructure first)
        self.write_data_file(self.gatk_analyze_saturation_mutagenesis, "multiqc_gatk_analyze_saturation_mutagenesis")

        self.gatk_analyze_saturation_mutagenesis_table(self.gatk_analyze_saturation_mutagenesis)
        # Add plots
        self.gatk_analyze_saturation_mutagenesis_plot_reads(self.gatk_analyze_saturation_mutagenesis)
        self.gatk_analyze_saturation_mutagenesis_plot_base_calls(self.gatk_analyze_saturation_mutagenesis)

        return n_reports_found

    def parse_read_counts_file(self, file_handle):
        """Parse a readCounts output file from GATK AnalyzeSaturationMutagenesis
        These files are tab delimited, with some hierarchical structuring.

        This does not parse percentages, since those are readily calculated
        and MultiQC automatically generates them for plotting."""

        keys = [
            "disjoint_pairs",
            "overlapping_pairs",
            "total_reads",
            "unmapped_reads",
            "lowq_reads",
            "evaluable_reads",
            "wt_reads_disjoint",
            "wt_reads_overlapping",
            "called_variants_disjoint",
            "called_variants_overlapping",
            "mate_ignored_disjoint",
            "inconsistent_overlapping",
            "low_quality_variation_disjoint",
            "low_quality_variation_overlapping",
            "insufficient_flank_disjoint",
            "insufficient_flank_overlapping",
            "total_base_calls",
            "evaluated_base_calls",
            "unevaluated_base_calls",
            "wt_total_reads",
            "variants_total_reads",
            "filtered_total_reads",
        ]

        # Initialize dictionary with all keys, since some may be missing
        data = {key: 0 for key in keys}

        # For keeping track of whether we are in disjoint or overlapping reads
        # Disjoint are first, so set to True initially
        disjoint_reads = True

        for line in file_handle:
            fields = line.split("\t")

            # First, note whether we are in disjoint or overlapping reads block

            if fields[0] == ">>Reads in disjoint pairs evaluated separately:":
                data["disjoint_pairs"] = int(fields[1])
                disjoint_reads = True

            elif fields[0] == ">>Reads in overlapping pairs evaluated together:":
                data["overlapping_pairs"] = int(fields[1])
                disjoint_reads = False

            # Proceed with the rest of the data
            elif fields[0] == "Total Reads:":
                data["total_reads"] = int(fields[1])

            elif fields[0] == ">Unmapped Reads:":
                data["unmapped_reads"] = int(fields[1])

            elif fields[0] == ">LowQ Reads:":
                data["lowq_reads"] = int(fields[1])

            elif fields[0] == ">Evaluable Reads:":
                data["evaluable_reads"] = int(fields[1])

            elif fields[0] == ">>>Wild type:":
                if disjoint_reads:
                    data["wt_reads_disjoint"] = int(fields[1])
                else:
                    data["wt_reads_overlapping"] = int(fields[1])

            elif fields[0] == ">>>Called variants:":
                if disjoint_reads:
                    data["called_variants_disjoint"] = int(fields[1])
                else:
                    data["called_variants_overlapping"] = int(fields[1])

            elif fields[0] == ">>>Mate ignored:":
                data["mate_ignored_disjoint"] = int(fields[1])

            elif fields[0] == ">>>Inconsistent pair:":
                data["inconsistent_overlapping"] = int(fields[1])

            elif fields[0] == ">>>Low quality variation:":
                if disjoint_reads:
                    data["low_quality_variation_disjoint"] = int(fields[1])
                else:
                    data["low_quality_variation_overlapping"] = int(fields[1])

            elif fields[0] == ">>>Insufficient flank:":
                if disjoint_reads:
                    data["insufficient_flank_disjoint"] = int(fields[1])
                else:
                    data["insufficient_flank_overlapping"] = int(fields[1])

            elif fields[0] == "Total base calls:":
                data["total_base_calls"] = int(fields[1])

            elif fields[0] == ">Base calls evaluated for variants:":
                data["evaluated_base_calls"] = int(fields[1])

            elif fields[0] == ">Base calls unevaluated:":
                data["unevaluated_base_calls"] = int(fields[1])

        # Create some summary fields

        data["wt_total_reads"] = data["wt_reads_disjoint"] + data["wt_reads_overlapping"]
        data["variants_total_reads"] = data["called_variants_disjoint"] + data["called_variants_overlapping"]
        data["filtered_reads"] = (
            data["lowq_reads"]
            + data["mate_ignored_disjoint"]
            + data["inconsistent_overlapping"]
            + data["low_quality_variation_disjoint"]
            + data["low_quality_variation_overlapping"]
            + data["insufficient_flank_disjoint"]
            + data["insufficient_flank_overlapping"]
        )

        return data

    def gatk_analyze_saturation_mutagenesis_plot_reads(self, data):
        """Make the plot for GATK AnalyzeSaturationMutagenesis read counts and add the section."""
        cats = OrderedDict()

        cats["filtered_reads"] = {"name": "Filtered reads", "color": "#3182bd"}
        cats["wt_total_reads"] = {"name": "WT reads", "color": "#9ecae1"}
        cats["variants_total_reads"] = {"name": "Variant reads", "color": "#deebf7"}
        config = {
            "id": "gatk_ASM_reads_plot",
            "title": "GATK AnalyzeSaturationMutagenesis: Read counts",
            "ylab": "Number of reads",
            "cpswitch_counts_label": "Counts",
        }

        self.add_section(
            name="Read counts",
            anchor="gatk-asm-read-counts",
            description="Read counts and read fate. Filtered reads include unmapped, low quality, and other pathologies.",
            helptext="""Reads can be filtered by GATK AnalyzeSaturationMutagenesis for a number of reasons, including low quality, insufficient flank, and other pathologies.
            This plot shows the number of reads that were mapped to WT, called as variants, and the number that were filtered.""",
            plot=bargraph.plot(data, cats, config),
        )

    def gatk_analyze_saturation_mutagenesis_plot_base_calls(self, data):
        """Make the plot for GATK AnalyzeSaturationMutagenesis base calls and add the section."""
        cats = OrderedDict()

        cats["evaluated_base_calls"] = {"name": "Base calls evaluated for variants", "648FFF": "#3182bd"}
        cats["unevaluated_base_calls"] = {"name": "Base calls not evaluated for variants", "color": "#deebf7"}
        config = {
            "id": "gatk_ASM_base_calls_plot",
            "title": "GATK AnalyzeSaturationMutagenesis: Base calls",
            "ylab": "Number of bases",
            "cpswitch_counts_label": "Counts",
        }

        self.add_section(
            name="Base calls",
            anchor="gatk-asm-bases",
            description="Base calls evaluated for variants and base calls not evaluated for variants.",
            helptext="""Bases can be filtered by GATK AnalyzeSaturationMutagenesis for a number of reasons, including low quality, insufficient flank, and other pathologies.
             This plot shows the number of base calls that were evaluated for variants and the number of base calls that were not evaluated for variants.""",
            plot=bargraph.plot(data, cats, config),
        )

    def gatk_analyze_saturation_mutagenesis_table(self, data):
        """Make the table for GATK AnalyzeSaturationMutagenesis and add the section."""
        asm_headers = {}

        asm_headers["total_reads"] = {
            "title": "Total reads (M)",
            "description": "Total reads in sample (millions)",
            "min": 0,
            "scale": "Greys",
            "shared_key": "read_count",
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["wt_total_reads"] = {
            "title": "WT reads (M)",
            "description": "Total evaluated reads mapped to WT (millions)",
            "min": 0,
            "scale": "Blues",
            "shared_key": "read_count",
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["variants_total_reads"] = {
            "title": "Variant reads (M)",
            "description": "Reads with a variant called (millions)",
            "min": 0,
            "scale": "Greens",
            "shared_key": "read_count",
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["filtered_reads"] = {
            "title": "Filtered reads (M)",
            "description": "Reads filtered from sample (millions)",
            "min": 0,
            "scale": "Reds-rev",
            "shared_key": "read_count",
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["unmapped_reads"] = {
            "title": "Unmapped (M)",
            "description": "Unmapped reads in sample (millions)",
            "min": 0,
            "scale": "Purples-rev",
            "shared_key": "read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["lowq_reads"] = {
            "title": "LowQ (M)",
            "description": "Low quality reads in sample (millions)",
            "min": 0,
            "scale": "Blues-rev",
            "shared_key": "read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["evaluable_reads"] = {
            "title": "Evaluable (M)",
            "description": "Evaluable reads in sample (millions)",
            "min": 0,
            "scale": "Greys",
            "shared_key": "read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["disjoint_pairs"] = {
            "title": "Disjoint (M)",
            "description": "Reads from disjoint (non-overlapping) paired-end reads in sample (millions)",
            "min": 0,
            "scale": "Oranges",
            "shared_key": "read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["wt_reads_disjoint"] = {
            "title": "WT (disjoint) (M)",
            "description": "WT reads called from disjoint pairs (millions)",
            "min": 0,
            "scale": "Purples",
            "shared_key": "disjoint_read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["called_variants_disjoint"] = {
            "title": "Called variants (disjoint) (M)",
            "description": "Reads with variants called from disjoint pairs (millions)",
            "min": 0,
            "scale": "Blues",
            "shared_key": "disjoint_read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["mate_ignored_disjoint"] = {
            "title": "Mateless (M)",
            "description": "Reads with ignored mates called from disjoint pairs (millions)",
            "min": 0,
            "scale": "Greens-rev",
            "shared_key": "disjoint_read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["low_quality_variation_disjoint"] = {
            "title": "LowQ (disjoint) (M)",
            "description": "Reads with low quality variation from disjoint pairs (millions)",
            "min": 0,
            "scale": "Reds-rev",
            "shared_key": "disjoint_read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["insufficient_flank_disjoint"] = {
            "title": "No flank (disjoint) (M)",
            "description": "Reads with insufficient flank from disjoint pairs (millions)",
            "min": 0,
            "scale": "Blues-rev",
            "shared_key": "disjoint_read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }

        asm_headers["overlapping_pairs"] = {
            "title": "Overlapping (M)",
            "description": "Reads from overlapping paired-end reads in sample (millions)",
            "min": 0,
            "scale": "Oranges",
            "shared_key": "read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["wt_reads_overlapping"] = {
            "title": "WT (overlapping) (M)",
            "description": "WT reads called from overlapping pairs (millions)",
            "min": 0,
            "scale": "Purples",
            "shared_key": "overlapping_read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["called_variants_overlapping"] = {
            "title": "Called variants (overlapping) (M)",
            "description": "Reads with variants called from overlapping pairs (millions)",
            "min": 0,
            "scale": "Blues",
            "shared_key": "overlapping_read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["inconsistent_overlapping"] = {
            "title": "Inconsistent (overlapping) (M)",
            "description": "Reads with inconsistent pairs from overlapping pairs (millions)",
            "min": 0,
            "scale": "Greens-rev",
            "shared_key": "overlapping_read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["low_quality_variation_overlapping"] = {
            "title": "Low quality reads (overlapping) (M)",
            "description": "Reads with low quality variation from overlapping pairs (millions)",
            "min": 0,
            "scale": "Reds-rev",
            "shared_key": "overlapping_read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }
        asm_headers["insufficient_flank_overlapping"] = {
            "title": "No flank (overlapping) (M)",
            "description": "Reads with insufficient flank from overlapping pairs (millions)",
            "min": 0,
            "scale": "Blues-rev",
            "shared_key": "overlapping_read_count",
            "hidden": True,
            "modify": lambda x: x / 1000000,
            "namespace": "GATK",
        }

        asm_headers["total_base_calls"] = {
            "title": "Total bases (B)",
            "description": "Total base calls in sample (billions)",
            "min": 0,
            "scale": "Greys",
            "shared_key": "base_calls",
            "hidden": True,
            "modify": lambda x: x / 1000000000,
            "namespace": "GATK",
        }
        asm_headers["evaluated_base_calls"] = {
            "title": "Evaluated bases (B)",
            "description": "Evaluated base calls in sample (billions)",
            "min": 0,
            "scale": "Blues",
            "hidden": True,
            "shared_key": "base_calls",
            "modify": lambda x: x / 1000000000,
            "namespace": "GATK",
        }
        asm_headers["unevaluated_base_calls"] = {
            "title": "Unevaluated bases (B)",
            "description": "Unevaluated base calls in sample (billions)",
            "min": 0,
            "scale": "Reds-rev",
            "hidden": True,
            "shared_key": "base_calls",
            "modify": lambda x: x / 1000000000,
            "namespace": "GATK",
        }
        config = {"id": "gatk_asm_stats", "namespace": "GATK", "table_title": "GATK ASM counts"}

        self.add_section(
            name="GATK ASM counts",
            anchor="gatk-asm-stats",
            description="Per-sample read count and base call fates from GATK AnalyzeSaturationMutagenesis.",
            helptext="""
            This table shows the distribution of calls (for reads or for bases) across all samples.
            Reads are categorized as WT, a variant, or filtered. Bases can be either evaluated or unevaluated, corresponding to the reads they come from.

            Reads are filtered for the following reasons:

            * Unmapped: the map quality is below the minimum MapQ (default = 4)'
            * Low quality reads: Reads are trimmed for quality > minQ (default 30) before calling variants. If the final size is less than the minimum length (default 15), or if the remaining segment does not cover the ORF, the read is filtered.

            Paired reads are also split into overlapping and disjoint sets, with further filters for both.

            * Inconsistent: If overlapping reads disagree on the called variant, the read is filtered.
            * Ignored mate: If the pairs are disjoint, the first read of the pair is used for analysis, and the second is ignored.
            * Low quality variation: If the variant includes ambiguous bases (not A, C, G, or T, or -), the read is filtered.
            * Insufficient flank: If the variant does not include a certain number of WT bases (default 2) flanking the variant, the read is filtered.
            """,
            plot=table.plot(data, asm_headers, config),
        )
