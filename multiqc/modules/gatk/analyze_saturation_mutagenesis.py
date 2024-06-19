"""MultiQC submodule to parse output from GATK tool AnalyzeSaturationMutagenesis"""

import logging
from collections import defaultdict

from multiqc import config
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

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

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

        label_to_key = {
            "Total Reads:": "total_reads",
            ">Unmapped Reads:": "unmapped_reads",
            ">LowQ Reads:": "lowq_reads",
            ">Evaluable Reads:": "evaluable_reads",
            ">>Reads in disjoint pairs evaluated separately:": "disjoint_pairs",
            ">>Reads in overlapping pairs evaluated together:": "overlapping_pairs",
            ">>>Wild type:": "wt_reads",
            ">>>Called variants:": "called_variants",
            ">>>Mate ignored:": "mate_ignored",
            ">>>Inconsistent pair:": "inconsistent",
            ">>>Low quality variation:": "low_quality_variation",
            ">>>Insufficient flank:": "insufficient_flank",
            "Total base calls:": "total_base_calls",
            ">Base calls evaluated for variants:": "evaluated_base_calls",
            ">Base calls unevaluated:": "unevaluated_base_calls",
        }

        data = defaultdict(int)

        # For keeping track of whether we are in disjoint or overlapping reads
        # Disjoint are first, so set to True initially
        disjoint_reads = True

        for line in file_handle:
            fields = line.split("\t")
            label = fields[0]
            key = label_to_key[label]
            if key == "disjoint_pairs":
                disjoint_reads = True
            elif key == "overlapping_pairs":
                disjoint_reads = False
            elif label.startswith(">>>"):
                key = key + "_" + ("disjoint" if disjoint_reads else "overlapping")
            data[key] = int(fields[1])

        # Create some summary fields
        data["wt_total_reads"] = data["wt_reads_disjoint"] + data["wt_reads_overlapping"]
        data["variants_total_reads"] = data["called_variants_disjoint"] + data["called_variants_overlapping"]
        data["filtered_reads"] = sum(
            [
                data["lowq_reads"],
                data["mate_ignored_disjoint"],
                data["inconsistent_overlapping"],
                data["low_quality_variation_disjoint"],
                data["low_quality_variation_overlapping"],
                data["insufficient_flank_disjoint"],
                data["insufficient_flank_overlapping"],
            ]
        )
        return data

    def gatk_analyze_saturation_mutagenesis_plot_reads(self, data):
        """Make the plot for GATK AnalyzeSaturationMutagenesis read counts and add the section."""
        cats = {
            "filtered_reads": {"name": "Filtered reads"},
            "wt_total_reads": {"name": "WT reads"},
            "variants_total_reads": {"name": "Variant reads"},
        }

        pconfig = {
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
            plot=bargraph.plot(data, cats, pconfig),
        )

    def gatk_analyze_saturation_mutagenesis_plot_base_calls(self, data):
        """Make the plot for GATK AnalyzeSaturationMutagenesis base calls and add the section."""
        cats = {
            "evaluated_base_calls": {"name": "Base calls evaluated for variants"},
            "unevaluated_base_calls": {"name": "Base calls not evaluated for variants"},
        }

        pconfig = {
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
            plot=bargraph.plot(data, cats, pconfig),
        )

    def gatk_analyze_saturation_mutagenesis_table(self, data):
        """Make the table for GATK AnalyzeSaturationMutagenesis and add the section."""
        asm_headers = {
            "total_reads": {
                "title": f"Total reads ({config.read_count_prefix})",
                "description": f"Total reads in sample ({config.read_count_desc})",
                "min": 0,
                "scale": "Greys",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "wt_total_reads": {
                "title": f"WT reads ({config.read_count_prefix})",
                "description": f"Total evaluated reads mapped to WT ({config.read_count_desc})",
                "min": 0,
                "scale": "Blues",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "variants_total_reads": {
                "title": f"Variant reads ({config.read_count_prefix})",
                "description": f"Reads with a variant called ({config.read_count_desc})",
                "min": 0,
                "scale": "Greens",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "filtered_reads": {
                "title": f"Filtered reads ({config.read_count_prefix})",
                "description": f"Reads filtered from sample ({config.read_count_desc})",
                "min": 0,
                "scale": "Reds-rev",
                "shared_key": "read_count",
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "unmapped_reads": {
                "title": f"Unmapped ({config.read_count_prefix})",
                "description": f"Unmapped reads in sample ({config.read_count_desc})",
                "min": 0,
                "scale": "Purples-rev",
                "shared_key": "read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "lowq_reads": {
                "title": f"LowQ ({config.read_count_prefix})",
                "description": f"Low quality reads in sample ({config.read_count_desc})",
                "min": 0,
                "scale": "Blues-rev",
                "shared_key": "read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "evaluable_reads": {
                "title": f"Evaluable ({config.read_count_prefix})",
                "description": f"Evaluable reads in sample ({config.read_count_desc})",
                "min": 0,
                "scale": "Greys",
                "shared_key": "read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "disjoint_pairs": {
                "title": f"Disjoint ({config.read_count_prefix})",
                "description": f"Reads from disjoint (non-overlapping) paired-end reads in sample ({config.read_count_desc})",
                "min": 0,
                "scale": "Oranges",
                "shared_key": "read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "wt_reads_disjoint": {
                "title": f"WT (disjoint) ({config.read_count_prefix})",
                "description": f"WT reads called from disjoint pairs ({config.read_count_desc})",
                "min": 0,
                "scale": "Purples",
                "shared_key": "disjoint_read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "called_variants_disjoint": {
                "title": f"Called variants (disjoint) ({config.read_count_prefix})",
                "description": f"Reads with variants called from disjoint pairs ({config.read_count_desc})",
                "min": 0,
                "scale": "Blues",
                "shared_key": "disjoint_read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "mate_ignored_disjoint": {
                "title": f"Mateless ({config.read_count_prefix})",
                "description": f"Reads with ignored mates called from disjoint pairs ({config.read_count_desc})",
                "min": 0,
                "scale": "Greens-rev",
                "shared_key": "disjoint_read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "low_quality_variation_disjoint": {
                "title": f"LowQ (disjoint) ({config.read_count_prefix})",
                "description": f"Reads with low quality variation from disjoint pairs ({config.read_count_desc})",
                "min": 0,
                "scale": "Reds-rev",
                "shared_key": "disjoint_read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "insufficient_flank_disjoint": {
                "title": f"No flank (disjoint) ({config.read_count_prefix})",
                "description": f"Reads with insufficient flank from disjoint pairs ({config.read_count_desc})",
                "min": 0,
                "scale": "Blues-rev",
                "shared_key": "disjoint_read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "overlapping_pairs": {
                "title": f"Overlapping ({config.read_count_prefix})",
                "description": f"Reads from overlapping paired-end reads in sample ({config.read_count_desc})",
                "min": 0,
                "scale": "Oranges",
                "shared_key": "overlapping_read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "wt_reads_overlapping": {
                "title": f"WT (overlapping) ({config.read_count_prefix})",
                "description": f"WT reads called from overlapping pairs ({config.read_count_desc})",
                "min": 0,
                "scale": "Purples",
                "shared_key": "overlapping_read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "called_variants_overlapping": {
                "title": f"Called variants (overlapping) ({config.read_count_prefix})",
                "description": f"Reads with variants called from overlapping pairs ({config.read_count_desc})",
                "min": 0,
                "scale": "Blues",
                "shared_key": "overlapping_read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "inconsistent_overlapping": {
                "title": f"Inconsistent (overlapping) ({config.read_count_prefix})",
                "description": f"Reads with inconsistent pairs from overlapping pairs ({config.read_count_desc})",
                "min": 0,
                "scale": "Greens-rev",
                "shared_key": "overlapping_read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "low_quality_variation_overlapping": {
                "title": f"Low quality reads (overlapping) ({config.read_count_prefix})",
                "description": f"Reads with low quality variation from overlapping pairs ({config.read_count_desc})",
                "min": 0,
                "scale": "Reds-rev",
                "shared_key": "overlapping_read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "insufficient_flank_overlapping": {
                "title": f"No flank (overlapping) ({config.read_count_prefix})",
                "description": f"Reads with insufficient flank from overlapping pairs ({config.read_count_desc})",
                "min": 0,
                "scale": "Blues-rev",
                "shared_key": "overlapping_read_count",
                "hidden": True,
                "modify": lambda x: x * config.read_count_multiplier,
                "namespace": "GATK",
            },
            "total_base_calls": {
                "title": f"Total bases ({config.base_count_prefix})",
                "description": f"Total base calls in sample ({config.base_count_desc})",
                "min": 0,
                "scale": "Greys",
                "shared_key": "base_calls",
                "hidden": True,
                "modify": lambda x: x * config.base_count_multiplier,
                "namespace": "GATK",
            },
            "evaluated_base_calls": {
                "title": f"Evaluated bases ({config.base_count_prefix})",
                "description": f"Evaluated base calls in sample ({config.base_count_desc})",
                "min": 0,
                "scale": "Blues",
                "hidden": True,
                "shared_key": "base_calls",
                "modify": lambda x: x * config.base_count_multiplier,
                "namespace": "GATK",
            },
            "unevaluated_base_calls": {
                "title": f"Unevaluated bases ({config.base_count_prefix})",
                "description": f"Unevaluated base calls in sample ({config.base_count_desc})",
                "min": 0,
                "scale": "Reds-rev",
                "hidden": True,
                "shared_key": "base_calls",
                "modify": lambda x: x * config.base_count_multiplier,
                "namespace": "GATK",
            },
        }

        # Add module specific prefix to all keys to be safe
        prefix = "gatk_ask_"
        asm_headers = {f"{prefix}{k}": v for k, v in asm_headers.items()}
        data = {sn: {f"{prefix}{k}": v for k, v in d.items()} for sn, d in data.items()}

        pconfig = {"id": f"{prefix}stats", "namespace": "GATK", "title": "GATK ASM counts"}

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
            plot=table.plot(data, asm_headers, pconfig),
        )
