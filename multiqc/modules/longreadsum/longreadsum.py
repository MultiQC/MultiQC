import copy
import logging
import os
import re
import json

from typing import Any, Union

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots.table_object import ColumnMeta
from multiqc.plots.table_object import ColumnDict

# from multiqc.base_module import SampleName, ColumnKey
from multiqc.plots import table
from multiqc.plots.table_object import TableConfig

log = logging.getLogger(__name__)

# https://docs.seqera.io/multiqc/reports
# https://github.com/MultiQC/MultiQC/blob/main/.github/CONTRIBUTING.md
# https://docs.seqera.io/multiqc/development/modules

SampleName = str
ColumnKey = str


class MultiqcModule(BaseMultiqcModule):
    """
    The module parses data in the `summary.txt` LongReadSum output files.
    """

    def __init__(self) -> None:
        super(MultiqcModule, self).__init__(
            name="LongReadSum",
            anchor="longreadsum",
            href="https://github.com/WGLab/LongReadSum.git",
            info="Quality control for long-read sequencing data.",
            extra="""
             LongReadSum: A fast and flexible quality control and signal
             summarization tool for long-read sequencing data.
            """,
            doi="10.1016/j.csbj.2025.01.019",
        )

        # Get data by sample
        data_by_sample: dict[SampleName | str, dict[ColumnKey | str, int | float | str | bool]] = dict()
        for f in self.find_log_files("longreadsum/summary"):
            # sample_name = f["s_name"]
            sample_name = SampleName(f["s_name"])
            log.debug("Processing file for sample: %s", sample_name)
            if sample_name in data_by_sample:
                log.debug("Duplicate sample name found! Overwriting: %s", sample_name)

            data_by_sample[sample_name] = parse_file(f["f"])

            # Add version information
            if ColumnKey("longreadsum_version") in data_by_sample[sample_name]:
                version = data_by_sample[sample_name][ColumnKey("longreadsum_version")]
                self.add_software_version(str(version), sample=sample_name)
                log.debug("Found LongReadSum version %s for sample %s", version, sample_name)
                # Remove the version from the data, not needed in the table
                del data_by_sample[sample_name][ColumnKey("longreadsum_version")]

        # Remove samples to ignore
        data_by_sample = self.ignore_samples(data_by_sample)

        # If no data found, raise an error
        if len(data_by_sample) == 0:
            log.warning("No LongReadSum summary files found.")
            raise ModuleNoSamplesFound(
                "No LongReadSum summary files found. "
                "Please ensure the output directory is correct and contains the expected files."
            )

        # If the filetype is sequence_summary.txt, it will contain keys "all",
        # "passed", and "failed" with their respective statistics.
        for sample, sample_data in data_by_sample.items():
            if sample_data.get(ColumnKey("File Type")) == "sequencing_summary":
                original_data = copy.deepcopy(sample_data)
                sample_data.clear()
                for key, value in original_data.items():
                    if key == "all" and isinstance(value, dict):
                        for sub_key, sub_value in value.items():
                            readable_sub_key = get_basic_stat_label(sub_key)
                            sample_data[readable_sub_key] = sub_value
                    else:
                        sample_data[key] = value

        # Basic statistics table
        log.debug("Creating basic statistics table.")
        # self.general_stats_addcols(data_by_sample)
        self.create_stats_table(data_by_sample)

        # If sequencing_summary file type, add table for passed and failed reads
        if any(sample_data["File Type"] == "sequencing_summary" for sample_data in data_by_sample.values()):
            self.add_passed_failed_reads_table(data_by_sample)

        # Aligment statistics tables
        if any("mapped" in sample_data for sample_data in data_by_sample.values()):
            log.debug("Creating read alignment statistics table.")
            self.add_read_alignment_stats(data_by_sample)

        if any("alignments" in sample_data for sample_data in data_by_sample.values()):
            log.debug("Creating read alignment type statistics table.")
            self.add_read_alignment_type_table(data_by_sample)

        if any("base_alignment" in sample_data for sample_data in data_by_sample.values()):
            log.debug("Creating base alignment type statistics table.")
            self.add_base_alignment_type_table(data_by_sample)

        # Base modification statistics
        if any("base_modifications" in sample_data for sample_data in data_by_sample.values()):
            log.debug("Creating base modification statistics table.")
            self.add_base_modification_table(data_by_sample)

        # TIN summary statistics
        if any("tin_data" in sample_data for sample_data in data_by_sample.values()):
            log.debug("Creating TIN summary statistics table.")
            self.add_tin_summary_table(data_by_sample)

        log.debug("Data processing complete. Adding sections to MultiQC report.")

    def create_stats_table(self, data: dict[SampleName | str, dict[ColumnKey | str, int | float | str | bool]]) -> None:
        """
        Create a basic statistics table from the parsed data.
        """
        headers: dict[str, dict[str, Any]] = {
            "File Type": {
                "title": "File Type",
                "description": "Type of the file (e.g., fasta, fastq, bam, sequencing_summary)",
            },
            "GC Content (%)": {
                "title": "GC Content",
                "description": "Percentage of G and C bases in the reads",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "GnBu",
                "shared_key": "gc_content",
            },
            "Total Reads": {
                "title": "Total Reads",
                "description": "Total number of reads in the file",
                "min": 0,
                "scale": "Oranges",
            },
            "Total Bases": {
                "title": "Total Bases",
                "description": "Total number of bases in the reads",
                "min": 0,
                "scale": "Purples",
            },
            "Mean Read Length": {
                "title": "Mean Read Length",
                "description": "Average length of the reads",
                "min": 0,
                "suffix": "bp",
                "scale": "GnBu",
            },
            "N50 Read Length": {
                "title": "N50 Read Length",
                "description": "N50 length of the reads",
                "min": 0,
                "suffix": "bp",
                "scale": "GnBu",
            },
            "Median Read Length": {
                "title": "Median Read Length",
                "description": "Median length of the reads",
                "min": 0,
                "suffix": "bp",
                "scale": "GnBu",
            },
            "Longest Read Length": {
                "title": "Longest Read Length",
                "description": "Length of the longest read",
                "min": 0,
                "suffix": "bp",
                "scale": "GnBu",
            },
        }
        self.general_stats_addcols(data, headers)

    def add_tin_summary_table(
        self, data: dict[SampleName | str, dict[ColumnKey | str, int | float | str | bool]]
    ) -> None:
        """
        Add a TIN summary statistics table.
        """
        tin_data = {}
        for sample, sample_data in data.items():
            tin_field_raw: int | float | str | bool | dict[str, Any] = sample_data.get("tin_data", {})
            if isinstance(tin_field_raw, dict):
                tin_field: dict[str, Any] = tin_field_raw
                for filepath in tin_field.keys():
                    filename = os.path.basename(filepath)

                    # Extract TIN statistics from the sample data
                    tin_stats_raw: int | float | str | bool | dict[str, Any] = tin_field.get(filepath, {})
                    if isinstance(tin_stats_raw, dict):
                        tin_stats: dict[str, Any] = tin_stats_raw

                        # Create a readable key for the sample
                        readable_sample = f"{sample} ({filename})"
                        tin_data[readable_sample] = {
                            "Total Transcripts": tin_stats.get("total_transcripts", 0),
                            "TIN Mean": tin_stats.get("mean", 0.0),
                            "TIN Median": tin_stats.get("median", 0.0),
                            "TIN StdDev": tin_stats.get("stddev", 0.0),
                        }

        # Prepare the headers for the TIN summary statistics table
        headers: dict[str | ColumnKey, ColumnDict] | None = {
            "Total Transcripts": {
                "title": "Total Transcripts",
                "description": "Total number of transcripts",
                "min": 0,
                "scale": "Blues",
            },
            "TIN Mean": {
                "title": "TIN Mean",
                "description": "Mean Transcript Integrity Number",
                "min": 0,
                "max": 100,
                "scale": "YlGn",
            },
            "TIN Median": {
                "title": "TIN Median",
                "description": "Median Transcript Integrity Number",
                "min": 0,
                "max": 100,
                "scale": "YlGn",
            },
            "TIN StdDev": {
                "title": "TIN StdDev",
                "description": "Standard Deviation of Transcript Integrity Number",
                "min": 0,
                "max": 100,
                "scale": "YlGn",
            },
        }

        # Add the TIN summary statistics table
        self.add_section(
            name="TIN Summary Statistics",
            anchor="tin-summary-stats",
            description="Transcript Integrity Number (TIN) summary statistics.",
            helptext="TIN statistics include total transcripts, mean, median, and standard deviation of TIN values.",
            plot=table.plot(
                tin_data,
                pconfig=TableConfig(id="longreadsum-tin-summary", title="TIN Summary Statistics"),
                headers=headers,
            ),
        )

    def add_passed_failed_reads_table(
        self, data: dict[SampleName | str, dict[ColumnKey | str, int | float | str | bool]]
    ) -> None:
        """
        Add a table for passed and failed reads statistics from sequencing_summary.txt files.
        """
        log.debug("Creating passed and failed reads table.")
        passed_data_by_sample = {}
        failed_data_by_sample = {}

        for sample, sample_data in data.items():
            if "File Type" in sample_data and sample_data["File Type"] == "sequencing_summary":
                passed_data: Any = sample_data.get("passed", {})
                failed_data: Any = sample_data.get("failed", {})

                # Remove the "all" key
                sample_data.pop("all", None)
                # Add passed and failed data to new dictionaries
                passed_data_by_sample[sample] = {k: v for k, v in passed_data.items()}
                failed_data_by_sample[sample] = {k: v for k, v in failed_data.items()}

        # Update the keys for each sample in passed and failed data
        for sample, passed_data in passed_data_by_sample.items():
            passed_data_fmt = {}
            for key, value in passed_data.items():
                readable_key = get_basic_stat_label(key)
                # Ensure values are always integers or floats
                if isinstance(value, str) and value.isdigit():
                    value = int(value)
                elif isinstance(value, str) and re.match(r"^\d+(\.\d+)?$", value):
                    value = float(value)
                passed_data_fmt[readable_key] = value
            passed_data_by_sample[sample] = passed_data_fmt

        for sample, failed_data in failed_data_by_sample.items():
            failed_data_fmt = {}
            for key, value in failed_data.items():
                readable_key = get_basic_stat_label(key)
                # Ensure values are always integers or floats
                if isinstance(value, str) and value.isdigit():
                    value = int(value)
                elif isinstance(value, str) and re.match(r"^\d+(\.\d+)?$", value):
                    value = float(value)
                failed_data_fmt[readable_key] = value
            failed_data_by_sample[sample] = failed_data_fmt

        headers: dict[str | ColumnKey, ColumnDict] | None = {
            "GC Content (%)": {
                "title": "GC Content",
                "description": "Percentage of G and C bases in the reads",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "GnBu",
                "shared_key": "gc_content",
            },
            "Total Reads": {
                "title": "Total Reads",
                "description": "Total number of reads in the file",
                "min": 0,
                "scale": "Oranges",
            },
            "Total Bases": {
                "title": "Total Bases",
                "description": "Total number of bases in the reads",
                "min": 0,
                "scale": "Purples",
            },
            "Mean Read Length": {
                "title": "Mean Read Length",
                "description": "Average length of the reads",
                "min": 0,
                "suffix": "bp",
                "scale": "GnBu",
            },
            "N50 Read Length": {
                "title": "N50 Read Length",
                "description": "N50 length of the reads",
                "min": 0,
                "suffix": "bp",
                "scale": "GnBu",
            },
            "Median Read Length": {
                "title": "Median Read Length",
                "description": "Median length of the reads",
                "min": 0,
                "suffix": "bp",
                "scale": "GnBu",
            },
            "Longest Read Length": {
                "title": "Longest Read Length",
                "description": "Length of the longest read",
                "min": 0,
                "suffix": "bp",
                "scale": "GnBu",
            },
        }

        # Add the passed reads table
        self.add_section(
            name="Passed Reads Statistics",
            anchor="passed-reads-stats",
            description="Statistics for reads that passed quality control.",
            helptext="Passed reads statistics include counts of reads that passed quality control.",
            plot=table.plot(
                passed_data_by_sample,
                headers=headers,
                pconfig=TableConfig(id="longreadsum-passed-reads", title="Passed Reads Statistics"),
            ),
        )

        # Add the failed reads table
        self.add_section(
            name="Failed Reads Statistics",
            anchor="failed-reads-stats",
            description="Statistics for reads that failed quality control.",
            helptext="Failed reads statistics include counts of reads that failed quality control.",
            plot=table.plot(
                failed_data_by_sample,
                headers=headers,
                pconfig=TableConfig(id="longreadsum-failed-reads", title="Failed Reads Statistics"),
            ),
        )

    def add_base_modification_table(
        self, data: dict[SampleName | str, dict[ColumnKey | str, int | float | str | bool]]
    ) -> None:
        """
        Add a base modification statistics table.
        """

        # Prepare data for the table
        base_modification_data = {
            sample: sample_data.get("base_modifications", 0) for sample, sample_data in data.items()
        }

        # There is a key "base_mod_counts" that contains a dictionary of counts
        # for each modification type. Flatten this into the main dictionary
        for sample, mod_data in base_modification_data.items():
            if isinstance(mod_data, dict):
                base_mod_counts = mod_data.get("base_mod_counts")
                if isinstance(base_mod_counts, dict):
                    for mod_type, count in base_mod_counts.items():
                        # Convert the modification type to a readable format
                        readable_mod_type = f"Modified Base Count ({mod_type})"
                        mod_data[readable_mod_type] = count
                    # Remove the original base_mod_counts key
                    del mod_data["base_mod_counts"]

                base_mod_counts_forward = mod_data.get("base_mod_counts_forward")
                if isinstance(base_mod_counts_forward, dict):
                    for mod_type, count in base_mod_counts_forward.items():
                        readable_mod_type = f"Modified Base Count ({mod_type}, Forward Strand)"
                        mod_data[readable_mod_type] = count
                    del mod_data["base_mod_counts_forward"]

                base_mod_counts_reverse = mod_data.get("base_mod_counts_reverse")
                if isinstance(base_mod_counts_reverse, dict):
                    for mod_type, count in base_mod_counts_reverse.items():
                        readable_mod_type = f"Modified Base Count ({mod_type}, Reverse Strand)"
                        mod_data[readable_mod_type] = count
                    del mod_data["base_mod_counts_reverse"]

        # Update the keys for each sample
        base_modification_data_fmt: dict[SampleName | str, dict[ColumnKey | str, int | float | str | bool]] = {}
        for sample, mod_data in base_modification_data.items():
            base_modification_data_fmt[sample] = {}
            if isinstance(mod_data, dict):
                for key, value in mod_data.items():
                    readable_key = get_base_mod_label(key)
                    # Ensure values are always integers or floats
                    if isinstance(value, str) and value.isdigit():
                        value = int(value)
                    elif isinstance(value, str) and re.match(r"^\d+(\.\d+)?$", value):
                        value = float(value)
                    base_modification_data_fmt[sample][readable_key] = value

        # Create a table configuration
        pconfig = TableConfig(
            id="longreadsum-base-modifications",
            title="Base Modifications",
        )

        # Define headers for the base modification statistics table
        headers: dict[str | ColumnKey, ColumnDict] | None = {
            "Total Predictions": {
                "title": "Total Predictions",
                "description": "Total number of base modification predictions",
                "min": 0,
                "scale": "Blues",
            },
            "Filter Threshold": {
                "title": "Filter Threshold",
                "description": "Filter threshold for base modification predictions",
                "min": 0,
                "scale": "Greens",
            },
            "Modified Base Count": {
                "title": "Modified Base Count",
                "description": "Total count of modified bases across all strands",
                "min": 0,
                "scale": "Reds",
            },
            "Modified Base Count (Forward Strand)": {
                "title": "Modified Base Count (Forward Strand)",
                "description": "Count of modified bases on the forward strand",
                "min": 0,
                "scale": "Purples",
            },
            "Modified Base Count (Reverse Strand)": {
                "title": "Modified Base Count (Reverse Strand)",
                "description": "Count of modified bases on the reverse strand",
                "min": 0,
                "scale": "Oranges",
            },
            "Modified Base Count (CpG Sites, Forward Strand)": {
                "title": "Modified Base Count (CpG Sites, Forward Strand)",
                "description": "Count of modified bases at CpG sites on the forward strand",
                "min": 0,
                "scale": "YlGn",
            },
            "Modified Base Count (CpG Sites, Reverse Strand)": {
                "title": "Modified Base Count (CpG Sites, Reverse Strand)",
                "description": "Count of modified bases at CpG sites on the reverse strand",
                "min": 0,
                "scale": "YlOrRd",
            },
        }

        # Add the base modification statistics table
        self.add_section(
            name="Base Modification Statistics",
            anchor="base-modification-stats",
            description="Statistics related to base modifications.",
            helptext="Base modification statistics include counts of modified bases, "
            "filter thresholds, and modified base counts for forward and reverse strands.",
            plot=table.plot(base_modification_data_fmt, pconfig=pconfig, headers=headers),
        )

    def add_base_alignment_type_table(
        self, data: dict[SampleName | str, dict[ColumnKey | str, int | float | str | bool]]
    ) -> None:
        """
        Add a base alignment type statistics table.
        """
        mapped_dict = {sample: sample_data.get("base_alignment", 0) for sample, sample_data in data.items()}

        # Update the keys for each sample
        mapped_dict_fmt: dict[SampleName | str, dict[ColumnKey | str, int | float | str | bool]] = {}
        for sample, aligned_data in mapped_dict.items():
            mapped_dict_fmt[sample] = {}
            if isinstance(aligned_data, dict):
                for key, value in aligned_data.items():
                    readable_key = get_base_alignment_label(key)
                    mapped_dict_fmt[sample][readable_key] = value

        # Define headers for the base alignment statistics table
        headers: dict[str | ColumnKey, ColumnDict] | None = {
            "Matched Bases": {
                "title": "Matched Bases",
                "description": "Number of bases that matched the reference",
                "min": 0,
                "scale": "Blues",
            },
            "Mismatched Bases": {
                "title": "Mismatched Bases",
                "description": "Number of bases that mismatched the reference",
                "min": 0,
                "scale": "Reds",
            },
            "Insertions": {
                "title": "Insertions",
                "description": "Number of inserted bases in the alignment",
                "min": 0,
                "scale": "Greens",
            },
            "Deletions": {
                "title": "Deletions",
                "description": "Number of deleted bases in the alignment",
                "min": 0,
                "scale": "Purples",
            },
            "Clipped Bases": {
                "title": "Clipped Bases",
                "description": "Number of clipped bases in the alignment",
                "min": 0,
                "scale": "Oranges",
            },
        }

        # Add a base alignment statistics table if alignment data is present
        self.add_section(
            name="Base Alignment Types",
            anchor="base-alignment-stats",
            description="Statistics related to base alignment.",
            helptext="Base alignment statistics include counts of matched, mismatched, inserted, deleted, and clipped bases.",
            plot=table.plot(
                mapped_dict_fmt,
                headers=headers,
                pconfig=TableConfig(id="longreadsum-base-alignment-types", title="Base Alignment Types"),
            ),
        )

    def add_read_alignment_type_table(
        self, data: dict[SampleName | str, dict[ColumnKey | str, int | float | str | bool]]
    ) -> None:
        """
        Add a read alignment type statistics table.
        """
        mapped_dict = {sample: sample_data.get("alignments", 0) for sample, sample_data in data.items()}

        # Update the keys for each sample
        mapped_dict_fmt: dict[SampleName | str, dict[ColumnKey | str, int | float | str | bool]] = {}
        for sample, aligned_data in mapped_dict.items():
            mapped_dict_fmt[sample] = {}
            if isinstance(aligned_data, dict):
                for key, value in aligned_data.items():
                    readable_key = get_read_alignment_label(key)
                    mapped_dict_fmt[sample][readable_key] = value

        headers: dict[str | ColumnKey, ColumnDict] | None = {
            "Primary Alignments": {
                "title": "Primary Alignments",
                "description": "Number of primary alignments for each read",
                "min": 0,
                "scale": "Blues",
            },
            "Secondary Alignments": {
                "title": "Secondary Alignments",
                "description": "Number of secondary alignments for each read",
                "min": 0,
                "scale": "Greens",
            },
            "Supplementary Alignments": {
                "title": "Supplementary Alignments",
                "description": "Number of supplementary alignments for each read",
                "min": 0,
                "scale": "Reds",
            },
            "Reads with Secondary Alignments": {
                "title": "Reads with Secondary Alignments",
                "description": "Number of reads that have secondary alignments",
                "min": 0,
                "scale": "Purples",
            },
            "Reads with Supplementary Alignments": {
                "title": "Reads with Supplementary Alignments",
                "description": "Number of reads that have supplementary alignments",
                "min": 0,
                "scale": "Oranges",
            },
            "Reads with Both Secondary and Supplementary Alignments": {
                "title": "Reads with Both Secondary and Supplementary Alignments",
                "description": "Number of reads that have both secondary and supplementary alignments",
                "min": 0,
                "scale": "Greys",
            },
            "Forward Alignments": {
                "title": "Forward Alignments",
                "description": "Number of forward alignments for each read",
                "min": 0,
                "scale": "YlGn",
            },
            "Reverse Alignments": {
                "title": "Reverse Alignments",
                "description": "Number of reverse alignments for each read",
                "min": 0,
                "scale": "YlOrRd",
            },
        }

        # Add a read alignment statistics table if alignment data is present
        self.add_section(
            name="Read Alignment Types",
            anchor="read-alignment-stats",
            description="Statistics related to read alignment types.",
            helptext="Read alignment types include primary, secondary, and supplementary alignments.",
            plot=table.plot(
                mapped_dict_fmt,
                headers=headers,
                pconfig=TableConfig(id="longreadsum-read-alignment-types", title="Read Alignment Types"),
            ),
        )

    def add_read_alignment_stats(
        self, data: dict[SampleName | str, dict[ColumnKey | str, int | float | str | bool]]
    ) -> None:
        """
        Add a read alignment statistics table.
        """
        mapped_dict = {sample: sample_data.get("mapped", 0) for sample, sample_data in data.items()}

        # Update the keys for each sample
        mapped_dict_fmt: dict[SampleName | str, dict[ColumnKey | str, int | float | str | bool]] = {}
        for sample, aligned_data in mapped_dict.items():
            mapped_dict_fmt[sample] = {}
            if isinstance(aligned_data, dict):
                for key, value in aligned_data.items():
                    readable_key = get_basic_stat_label(key)
                    mapped_dict_fmt[sample][readable_key] = value

        # Define headers for the read alignment statistics table
        headers: dict[str | ColumnKey, ColumnDict] | None = {
            "GC Content (%)": {
                "title": "GC Content",
                "description": "Percentage of G and C bases in the reads",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "GnBu",
                "shared_key": "gc_content",
            },
            "Total Reads": {
                "title": "Total Reads",
                "description": "Total number of reads in the file",
                "min": 0,
                "scale": "Oranges",
            },
            "Total Bases": {
                "title": "Total Bases",
                "description": "Total number of bases in the reads",
                "min": 0,
                "scale": "Purples",
            },
            "Mean Read Length": {
                "title": "Mean Read Length",
                "description": "Average length of the reads",
                "min": 0,
                "suffix": "bp",
                "scale": "GnBu",
            },
            "N50 Read Length": {
                "title": "N50 Read Length",
                "description": "N50 length of the reads",
                "min": 0,
                "suffix": "bp",
                "scale": "GnBu",
            },
            "Median Read Length": {
                "title": "Median Read Length",
                "description": "Median length of the reads",
                "min": 0,
                "suffix": "bp",
                "scale": "GnBu",
            },
            "Longest Read Length": {
                "title": "Longest Read Length",
                "description": "Length of the longest read",
                "min": 0,
                "suffix": "bp",
                "scale": "GnBu",
            },
        }

        # Add a read alignment statistics table if alignment data is present
        self.add_section(
            name="Read Alignment Statistics",
            anchor="read-alignment-stats",
            description="Basic statistics for aligned reads.",
            helptext="Read alignment statistics include counts of primary, secondary, and supplementary alignments.",
            plot=table.plot(
                mapped_dict_fmt,
                headers=headers,
                pconfig=TableConfig(id="longreadsum-read-alignment", title="Read Alignment Statistics"),
            ),
        )


def get_basic_stat_label(key: ColumnKey) -> ColumnKey:
    """
    Get the readable label for a basic statistic key.
    """
    basic_stats_labels = {
        "filetype": "File Type",
        "gc_percent": "GC Content (%)",
        "total_num_reads": "Total Reads",
        "total_num_bases": "Total Bases",
        "mean_read_length": "Mean Read Length",
        "n50_read_length": "N50 Read Length",
        "median_read_length": "Median Read Length",
        "longest_read_length": "Longest Read Length",
    }
    label = basic_stats_labels.get(key, key)
    return ColumnKey(label) if isinstance(label, str) else label


def get_read_alignment_label(key: ColumnKey) -> ColumnKey:
    """
    Get the readable label for a read alignment key.
    """
    read_alignment_labels = {
        "primary": "Primary Alignments",
        "secondary": "Secondary Alignments",
        "supplementary": "Supplementary Alignments",
        "reads_with_secondary": "Reads with Secondary Alignments",
        "reads_with_supplementary": "Reads with Supplementary Alignments",
        "reads_with_both": "Reads with Both Secondary and Supplementary Alignments",
        "forward": "Forward Alignments",
        "reverse": "Reverse Alignments",
    }
    label = read_alignment_labels.get(key, key)
    return ColumnKey(label) if isinstance(label, str) else label


def get_base_alignment_label(key: ColumnKey) -> ColumnKey:
    """
    Get the readable label for a base alignment key.
    """
    base_alignment_labels = {
        "matched": "Matched Bases",
        "mismatched": "Mismatched Bases",
        "insertions": "Insertions",
        "deletions": "Deletions",
        "clipped": "Clipped Bases",
    }
    label = base_alignment_labels.get(key, key)
    return ColumnKey(label) if isinstance(label, str) else label


def get_base_mod_label(key: ColumnKey) -> ColumnKey:
    """
    Get the readable label for a base modification key.
    """
    base_mod_labels = {
        "unfiltered_modifications": "Total Predictions",
        "filter_threshold": "Filter Threshold",
        "sample_modified_base_count": "Modified Base Count",
        "sample_modified_base_count_forward": "Modified Base Count (Forward Strand)",
        "sample_modified_base_count_reverse": "Modified Base Count (Reverse Strand)",
        "cpg_forward": "Modified Base Count (CpG Sites, Forward Strand)",
        "cpg_reverse": "Modified Base Count (CpG Sites, Reverse Strand)",
    }
    # return base_mod_labels.get(key, key)
    label = base_mod_labels.get(key, key)
    return ColumnKey(label) if isinstance(label, str) else label


# def parse_file(raw_text: str) -> Dict[ColumnKey, Union[float, int]]:
def parse_file(raw_text: str) -> dict[ColumnKey | str, int | float | str | bool]:
    """
    Parse the summary.json content from a raw text string and return a dictionary of values.
    """
    data: dict[ColumnKey | str, int | float | str | bool] = {}
    json_data = json.loads(raw_text)
    for key, value in json_data.items():
        # Convert keys to readable format
        readable_key = get_basic_stat_label(key)

        # Convert values to appropriate types
        if isinstance(value, str) and value.isdigit():
            value = int(value)
        elif isinstance(value, str) and re.match(r"^\d+(\.\d+)?$", value):
            value = float(value)
        # Add to data dictionary
        if readable_key == "File Type":
            # Ensure file type is always a string
            data[readable_key] = str(value)
        else:
            data[readable_key] = value

    return data
