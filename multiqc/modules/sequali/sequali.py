import json
import logging
import textwrap
from collections import defaultdict
from typing import Dict, Any

import multiqc
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph, table

log = logging.getLogger(__name__)


PHRED_SCORE_EXPLANATION = textwrap.dedent(
    """
    As Phred scores are logarithmic, the means are calculated by 
    calculating the probability for each base and then averaging that 
    over the total number of bases. The probability is then converted
    back into a Phred score. Tools that average Phred scores naively
    are prone to overestimate the average quality by orders of 
    magnitude. As such Sequali might give a different plot here than 
    other QC tools.
"""
)

DUPLICATION_EXPLANATION = textwrap.dedent(
    """
    [The methodology to estimate duplication uses fingerprinting 
    with subsampling based on the fingerprints themselves](
    https://www.usenix.org/system/files/conference/atc13/atc13-xie.pdf). 
    This mitigates biases that might occur in estimates that only 
    look at the first reads.
    
    Sequali fingerprints by combining an 8 bp fragment at an offset 
    of 64 bp from the beginning with an 8 bp fragment offset at 64
    bp from the end. The offsets were chosen to limit the chance 
    of adapter sequences contaminating the fingerprint.
    """
)


def avg_x_label(x_label: str):
    if "-" not in x_label:
        return int(x_label)
    min_x, max_x = x_label.split("-")
    min_x = int(min_x)
    max_x = int(max_x)
    return min_x + ((max_x - min_x) // 2)


def prune_sample_dict(sample_dict: Dict[str, Any]):
    """
    Function to remove unused keys from the parsed data. This prevents loading
    them into memory long-term wich reduces memory usage.
    """
    keys_to_delete = [
        # This is a per sample data gathering of all the base qualities
        # per position for a stacked bar plot to clearly see the distribution.
        # Too hard to aggregate for MultiQC.
        "per_position_quality_distribution",
        "per_position_quality_distribution_read2",
        "per_tile_quality",
        # Per tile quality is not analysed by multiqc and uses a lot of space
        "per_tile_quality_read2"
        # Nanopore metrics for pore data do not have modules yet
        "nanopore_metrics",
    ]
    for key in keys_to_delete:
        if key in sample_dict:
            del sample_dict[key]

    for key in ("per_position_mean_quality_and_spread", "per_position_mean_quality_and_spread_read2"):
        if key in sample_dict:
            # Only the mean is used for the per position mean quality and spread data,
            # so remove the rest of the data
            percentiles_list = sample_dict[key]["percentiles"]
            new_percentiles_list = [
                (percentile, values) for percentile, values in percentiles_list if percentile == "mean"
            ]
            sample_dict[key]["percentiles"] = new_percentiles_list


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Sequali",
            anchor="sequali",
            href="https://github.com/rhpvorderman/sequali",
            info="Sequencing quality control for both long-read and short-read data",
            extra="""
            Features adapter search, overrepresented sequence  analysis and duplication analysis and supports
            FASTQ and uBAM inputs.
            """,
            doi="10.5281/zenodo.10822485",
        )

        versions = set()
        data = {}
        min_lengths = set()
        max_lengths = set()
        self.paired_end = False
        for sample_file in self.find_log_files("sequali", filehandles=True):
            sample_name = sample_file["s_name"]
            if self.is_ignore_sample(sample_name):
                continue
            self.add_data_source(sample_file)
            filename = sample_file["fn"]
            filehandle = sample_file["f"]
            try:
                sample_dict = json.load(filehandle)
            except json.JSONDecodeError:
                log.error(f"Could not decode JSON data in {filename}")
                continue
            # Save memory by pruning the sample dict's unused keys.
            prune_sample_dict(sample_dict)
            try:
                meta = sample_dict["meta"]
                filename1 = meta["filename"]
                filename2 = meta.get("filename_read2")
                if filename2:
                    if not self.paired_end:
                        self.paired_end = True
                    sample_name = self.clean_s_name([filename1, filename2])
            except KeyError:
                log.error("JSON file is not a proper Sequali report")
                continue
            data[sample_name] = sample_dict
        if len(data) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(data)} reports")

        for sample_name, sample_dict in data.items():
            sequali_version = sample_dict["meta"]["sequali_version"]
            versions.add(sequali_version)
            min_lengths.add(sample_dict["summary"]["minimum_length"])
            max_lengths.add(sample_dict["summary"]["maximum_length"])
            self.add_software_version(sequali_version, sample_name)

        summary_data = {sample_name: sample_dict["summary"] for sample_name, sample_dict in data.items()}
        self.write_data_file(summary_data, "multiqc_sequali")

        if len(versions) != 1:
            log.warning(f"Multiple Sequali versions found: {','.join(versions)}")
        max_length = sorted(max_lengths, reverse=True)[0]
        self.max_length = max_length
        min_length = sorted(min_lengths)[0]
        self.lengths_differ = True
        self.use_xlog = False
        if len(min_lengths) == 1 and len(max_lengths) == 1 and min_length == max_length:
            self.lengths_differ = False
        if max_length >= 1000:
            self.use_xlog = True

        self.sequali_general_stats(data)
        self.read_count_plot(data)
        self.per_position_quality_plot(data)
        self.per_sequence_quality_plot(data)
        self.per_position_gc_content_plot(data)
        self.per_sequence_gc_content_plot(data)
        self.sequence_length_distribution_plot(data)
        self.sequence_duplication_levels_plot(data)
        self.top_overrepresented_sequences_table(data)
        self.adapter_content_plot(data)
        self.adapter_content_from_overlap(data)

    def sequali_general_stats(self, data):
        general_stats = dict()
        for sample_name, sample_dict in data.items():
            summary = sample_dict["summary"]
            stats_entry = {
                "sequali_gc_percentage": (
                    100 * summary["total_gc_bases"] / max(summary["total_bases"] - summary["total_n_bases"], 1)
                ),
                "sequali_mean_sequence_length": summary["mean_length"],
                "sequali_total_reads": summary["total_reads"],
                "sequali_duplication_percentage": (1 - sample_dict["duplication_fractions"]["remaining_fraction"])
                * 100,
            }
            general_stats[sample_name] = stats_entry
        headers = {
            "sequali_gc_percentage": {
                "title": "GC %",
                "description": "Average GC content %",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "format": "{:.2f}",
                # GC percentage is not "wrong" in any case. We want a color
                # scale that clearly distuingishes GC rich vs AT-rich. The
                # middle 50% should be white, as it is neutral. The
                # color scale should not include red or green which are coupled
                # with the 'bad' or 'good' value judgement respectively.
                # This way we samples with a anomalous GC percentage will
                # clearly stand out from the others, which is desirable.
                # The PuOr scale fits these criteria.
                "scale": "PuOr",
            },
            "sequali_mean_sequence_length": {
                "title": "Mean length",
                "description": "Geometric mean of all lengths",
                "min": 0,
                "suffix": " bp",
                "format": "{:,.1f}",
                # Longer reads is considered better for long-read technologies.
                # Use green to signal that.
                "scale": "Greens",
            },
            "sequali_total_reads": {
                "title": "Total reads",
                "description": f"Total Sequences ({multiqc.config.read_count_desc})",
                "shared_key": "read_count",
                # More data is better, generally. Unfortunately we cannot use
                # greens here also, as it is in the adjacent column.
                # So we resort to another unicolor scale that is not red.
                # Between blues and purples, using Deep Purple to signal for
                # good quality has a distinct advantage:
                # https://www.youtube.com/watch?v=Q2FzZSBD5LE
                "scale": "Purples",
            },
            "sequali_duplication_percentage": {
                "title": "% est. dups.",
                "description": "Estimated duplication percentage",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "format": "{:.2f}",
                # The more the worse. Use Reds to signal.
                "scale": "Reds",
            },
        }
        self.general_stats_addcols(general_stats, headers)

    def read_count_plot(self, data):
        """Stacked bar plot showing counts of reads"""
        plot_data = {}
        for sample_name, sample_dict in data.items():
            total_reads = sample_dict["summary"]["total_reads"]
            remaining_percentage = sample_dict["duplication_fractions"]["remaining_fraction"]
            unique_reads = round(remaining_percentage * total_reads)
            duplicated_reads = total_reads - unique_reads
            plot_data[sample_name] = {
                "Unique Reads": unique_reads,
                "Duplicate Reads": duplicated_reads,
            }
        plot = bargraph.plot(
            plot_data,
            cats=["Unique Reads", "Duplicate Reads"],
            pconfig={
                "id": "sequali_sequence_counts_plot",
                "title": "Sequali: Sequence Counts",
                "ylab": "Number of reads",
                "cpswitch_counts_label": "Number of reads",
                "hide_empty": False,
            },
        )

        self.add_section(
            name="Sequence Counts",
            anchor="sequali_sequence_counts",
            description="Sequence counts for each sample.  Duplicate read counts are an estimate.",
            helptext=textwrap.dedent(
                """
                This plots shows the total number of reads broken down into 
                unique and duplicate reads. 
                """
                + DUPLICATION_EXPLANATION
            ),
            plot=plot,
        )

    def per_position_quality_plot(self, data):
        all_x_labels = set()

        for read in (1, 2):
            if read == 2 and not self.paired_end:
                continue
            if read == 1:
                key = "per_position_mean_quality_and_spread"
                id_suffix = "_read1" if self.paired_end else ""
                title_suffix = ": Read 1" if self.paired_end else ""
            else:
                id_suffix = "_read2"
                title_suffix = ": Read 2"
                key = "per_position_mean_quality_and_spread_read2"
            plot_data = {}
            for sample_name, sample_dict in data.items():
                per_pos_qual = sample_dict.get(key)
                if not per_pos_qual:
                    continue
                x_labels = per_pos_qual["x_labels"]
                all_x_labels.update(set())
                percentiles = dict(per_pos_qual["percentiles"])
                mean_series = percentiles["mean"]

                plot_data[sample_name] = {avg_x_label(x_label): m for x_label, m in zip(x_labels, mean_series)}

            plot_config = {
                "id": "sequali_per_position_quality_plot" + id_suffix,
                "title": "Sequali: Per Position Mean Quality" + title_suffix,
                "ylab": "Phred Score",
                "xlab": "Position (bp)",
                "xlog": self.use_xlog,
                "ymin": 0,
                "xmin": 0,
            }

            self.add_section(
                name="Sequence Quality Per Position" + title_suffix,
                anchor="sequali_sequence_quality_per_position" + id_suffix,
                description="The mean quality value across each base position.",
                helptext=textwrap.dedent(
                    """
                Only mean scores are plotted. The means are approximated as Sequali
                stores 12 phred categories per position: 0-3, 4-7, etc up to 44 and 
                higher. It does not store all 94 discrete phred score counts for 
                each position. For context, Illumina FASTQ files only 
                utilize four different phred scores.
                """
                )
                + PHRED_SCORE_EXPLANATION,
                # Note: no need to publicly shame those other QC tools. Bug reports
                # have been submitted.
                plot=linegraph.plot(plot_data, plot_config),
            )

    def per_sequence_quality_plot(self, data):
        for read in (1, 2):
            if read == 2 and not self.paired_end:
                continue
            if read == 1:
                key = "per_sequence_quality_scores"
                id_suffix = "_read1" if self.paired_end else ""
                title_suffix = ": Read 1" if self.paired_end else ""
            else:
                id_suffix = "_read2"
                title_suffix = ": Read 2"
                key = "per_sequence_quality_scores_read2"
            plot_data = {}
            for sample_name, sample_dict in data.items():
                qual_dict = sample_dict.get(key)
                if qual_dict is None:
                    continue
                x_labels = qual_dict["x_labels"]
                average_quality_counts = qual_dict["average_quality_counts"]
                plot_data[sample_name] = {int(phred): count for phred, count in zip(x_labels, average_quality_counts)}

            plot_config = {
                "id": "sequali_per_sequence_quality_scores_plot" + id_suffix,
                "title": "Sequali: Per Sequence Average Quality Scores" + title_suffix,
                "ylab": "Number of Sequences",
                "xlab": "Mean Sequence Quality (Phred Score)",
                "ymin": 0,
                "xmin": 0,
            }

            self.add_section(
                name="Per Sequence Average Quality Scores" + title_suffix,
                anchor="sequali_per_sequence_quality_scores" + id_suffix,
                description="The number of reads with average quality scores.",
                helptext=textwrap.dedent(
                    """
                    Shows the quality score profile on a read level. As Illumina
                    FASTQ files only utilize four different phred scores, the plot
                    may look a bit erratic at times. Due to the logarithmic nature
                    of Phred scores, lower Phred scores have a more significant
                    impact on the average quality as than higher phred scores.
                """
                )
                + PHRED_SCORE_EXPLANATION,
                plot=linegraph.plot(plot_data, plot_config),
            )

    def per_position_gc_content_plot(self, data):
        for read in (1, 2):
            if read == 2 and not self.paired_end:
                continue
            if read == 1:
                key = "per_position_base_content"
                id_suffix = "_read1" if self.paired_end else ""
                title_suffix = ": Read 1" if self.paired_end else ""
            else:
                id_suffix = "_read2"
                title_suffix = ": Read 2"
                key = "per_position_base_content_read2"
            plot_data = {}
            for sample_name, sample_dict in data.items():
                per_position_base_content = sample_dict.get(key)
                if per_position_base_content is None:
                    continue
                x_labels = per_position_base_content["x_labels"]
                x_positions = [avg_x_label(x_label) for x_label in x_labels]
                c_fractions = per_position_base_content["C"]
                g_fractions = per_position_base_content["G"]
                sample_gc_percentages = (
                    (c_fraction + g_fraction) * 100 for c_fraction, g_fraction in zip(c_fractions, g_fractions)
                )
                plot_data[sample_name] = dict(zip(x_positions, sample_gc_percentages))

            plot_config = {
                "id": "sequali_per_position_gc_content_plot" + id_suffix,
                "title": "Sequali: Per Position GC Content" + title_suffix,
                "xlab": "Position (bp)",
                "ylab": "% GC",
                "ymin": 0,
                "ymax": 100,
                "xlog": self.use_xlog,
            }

            self.add_section(
                name="Per Position GC Content" + title_suffix,
                anchor="per_posistion_gc_content" + id_suffix,
                description="The GC content percentage at each position for each sample.",
                plot=linegraph.plot(plot_data, plot_config),
            )

    def per_sequence_gc_content_plot(self, data):
        for read in (1, 2):
            if read == 2 and not self.paired_end:
                continue
            if read == 1:
                key = "per_sequence_gc_content"
                id_suffix = "_read1" if self.paired_end else ""
                title_suffix = ": Read 1" if self.paired_end else ""
            else:
                id_suffix = "_read2"
                title_suffix = ": Read 2"
                key = "per_sequence_gc_content_read2"
            plot_data = {}
            for sample_name, sample_dict in data.items():
                gc_dict = sample_dict.get(key)
                if gc_dict is None:
                    continue
                # Take the smoothened results here. Less resolution, but also less
                # confusing at first glance.
                x_labels = gc_dict["smoothened_x_labels"]
                gc_content_counts = gc_dict["smoothened_gc_content_counts"]
                total = max(sum(gc_content_counts), 1)
                plot_data[sample_name] = {
                    avg_x_label(x_label): count * 100 / total for x_label, count in zip(x_labels, gc_content_counts)
                }

            plot_config = {
                "id": "sequali_per_sequence_gc_content_plot" + id_suffix,
                "title": "Sequali: Per Sequence GC Content" + title_suffix,
                "xlab": "% GC",
                "ylab": "Percentage",
                "ymin": 0,
                "xmin": 0,
                "xmax": 100,
            }

            self.add_section(
                name="Per Sequence GC Content" + title_suffix,
                anchor="sequali_per_sequence_gc_content" + id_suffix,
                description="The GC content distribution of the sequences for each sample.",
                plot=linegraph.plot(plot_data, plot_config),
            )

    def sequence_length_distribution_plot(self, data):
        if not self.lengths_differ:
            self.add_section(
                name="Sequence Length Distribution",
                anchor="sequali_sequence_length_distribution",
                description=f"All the sequences are the same length at {self.max_length} bp.",
            )
            return
        for read in (1, 2):
            if read == 2 and not self.paired_end:
                continue
            if read == 1:
                key = "sequence_length_distribution"
                id_suffix = "_read1" if self.paired_end else ""
                title_suffix = ": Read 1" if self.paired_end else ""
            else:
                id_suffix = "_read2"
                title_suffix = ": Read 2"
                key = "sequence_length_distribution_read2"
            plot_data = dict()
            for sample_name, sample_dict in data.items():
                seqlength_dict = sample_dict.get(key)
                if seqlength_dict is None:
                    continue
                x_labels = seqlength_dict["length_ranges"]
                counts = seqlength_dict["counts"]
                plot_data[sample_name] = {avg_x_label(x_label): count for x_label, count in zip(x_labels, counts)}

            plot_config = {
                "id": "sequali_sequence_length_distribution_plot" + id_suffix,
                "title": "Sequali: Sequence Length Distribution" + title_suffix,
                "ylab": "Read Count",
                "ymin": 0,
                "xlab": "Sequence Length (bp)",
                "xlog": self.use_xlog,
            }

            self.add_section(
                name="Sequence Length Distribution" + title_suffix,
                anchor="sequali_sequence_length_distribution" + id_suffix,
                description="The distribution of read lengths found.",
                plot=linegraph.plot(plot_data, plot_config),
            )

    def sequence_duplication_levels_plot(self, data):
        plot_data = {}
        for sample_name, sample_dict in data.items():
            estimated_fractions = sample_dict["duplication_fractions"]["estimated_duplication_fractions"]
            plot_data[sample_name] = {label: fraction * 100 for label, fraction in estimated_fractions.items()}

        plot_config = {
            "id": "sequali_sequence_duplication_levels_plot",
            "title": "Sequali: Sequence Duplication Levels",
            "categories": True,
            "ylab": "% of Library",
            "xlab": "Sequence Duplication Level",
            "ymin": 0,
        }

        self.add_section(
            name="Sequence Duplication Levels",
            anchor="sequali_sequence_duplication_levels",
            description="The relative level of duplication found for every sequence.",
            helptext=DUPLICATION_EXPLANATION,
            plot=linegraph.plot(plot_data, plot_config),
        )

    def top_overrepresented_sequences_table(self, data):
        for read in (1, 2):
            if read == 2 and not self.paired_end:
                continue
            if read == 1:
                key = "overrepresented_sequences"
                id_suffix = "_read1" if self.paired_end else ""
                title_suffix = ": Read 1" if self.paired_end else ""
            else:
                id_suffix = "_read2"
                title_suffix = ": Read 2"
                key = "overrepresented_sequences_read2"
            sequence_matches = {}
            sequence_counts = defaultdict(lambda: 0)
            sequence_fractions = defaultdict(lambda: 0.0)
            for sample_data, sample_dict in data.items():
                overrepr_dict = sample_dict.get(key)
                if overrepr_dict is None:
                    continue
                overrepresented_sequences = overrepr_dict["overrepresented_sequences"]
                for entry in overrepresented_sequences:
                    sequence = entry["sequence"]
                    best_match = entry["best_match"]
                    fraction = entry["fraction"]
                    sequence_matches[sequence] = best_match
                    sequence_counts[sequence] += 1
                    sequence_fractions[sequence] += fraction
            # When sequences occur in an equal number of samples, use the cumulative
            # fraction dictionary to get the most present.
            most_present = sorted(
                sequence_counts.items(), key=lambda x: (x[1], sequence_fractions[x[0]]), reverse=True
            )[:20]
            total = len(data)
            table_data = {}
            for sequence, count in most_present:
                table_data[sequence] = {
                    "best_match": sequence_matches[sequence],
                    "libraries_affected": 100 * count / total,
                }

            table_headers = {
                "best_match": {
                    "title": "Best Match",
                    "description": "Best match as found by kmer analysis",
                },
                "libraries_affected": {
                    "title": "Libraries Affected (%)",
                    "description": "The percentage of libraries where this sequence is overrepresented.",
                    "max": 100,
                    "min": 0,
                    "suffix": "%",
                    "format": "{:.2f}",
                    # The more, the worse; use Red scaling.
                    "scale": "Reds",
                },
            }

            table_config = {
                "id": "sequali_top_overrepresented_sequences_table" + id_suffix,
                "title": "Sequali: top overrepresented sequences" + title_suffix,
                "defaultsort": [{"column": "libraries_affected", "order": "desc"}],
            }
            self.add_section(
                name="Top overrepresented sequences" + title_suffix,
                anchor="sequali_top_overrepresented_sequences" + id_suffix,
                description="The top 20 overrepresented sequences in all libraries",
                plot=table.plot(table_data, table_headers, table_config),
            )

    def adapter_content_plot(self, data):
        plot_data = {}

        for sample_name, sample_dict in data.items():
            adapter_content = sample_dict.get("adapter_content")
            if adapter_content is None:
                continue
            x_labels = adapter_content["x_labels"]
            x_points = [avg_x_label(x_label) for x_label in x_labels]
            adapters_and_series = adapter_content["adapter_content"]
            for adapter, series in adapters_and_series:
                if max(series) < 0.001:
                    # Probably false positives for long reads.
                    continue
                series_name = f"{sample_name}-{adapter}"
                plot_data[series_name] = dict(zip(x_points, series))

        if not plot_data:
            return

        plot_config = {
            "id": "sequali_adapter_content_plot",
            "title": "Sequali: Adapter Content",
            "ylab": "% of Sequences",
            "xlab": "Position (bp)",
            "ymax": 100,
            "ymin": 0,
            "xlog": self.use_xlog,
        }

        self.add_section(
            name="Adapter Content",
            anchor="sequali_adapter_content",
            description="The cumulative percentage count of the found adapter sequences",
            helptext="""
            Note that only samples with >= 0.1% adapter contamination are shown. 
            There may be several adapters detected per sample. For long read data there
            maybe more adapters per sample, this is a result of the false positive detection
            rate increasing with longer read length.
            """,
            plot=linegraph.plot(plot_data, plot_config),
        )

    def adapter_content_from_overlap(self, data):
        adapter_table_data = {}
        common_adapters_read1 = defaultdict(lambda: 0)
        common_adapters_read2 = defaultdict(lambda: 0)
        total_samples = 0
        for sample_name, sample_dict in data.items():
            content = sample_dict.get("adapter_content_from_overlap")
            if content is None:
                continue
            total_samples += 1
            number_of_adapters_read1 = content["number_of_adapters_read1"]
            number_of_adapters_read2 = content["number_of_adapters_read2"]
            most_common_adapter_read1 = content["longest_adapter_read1"]
            most_common_adapter_read2 = content["longest_adapter_read2"]
            most_common_adapter_read1_match = content["longest_adapter_read1_match"]
            most_common_adapter_read2_match = content["longest_adapter_read2_match"]
            common_adapters_read1[(most_common_adapter_read1, most_common_adapter_read1_match)] += 1
            common_adapters_read2[(most_common_adapter_read2, most_common_adapter_read2_match)] += 1
            total_reads = max(content["total_reads"], 1)
            adapter_table_data[sample_name] = {
                "adapter_content_read1": number_of_adapters_read1 / total_reads,
                "adapter_content_read2": number_of_adapters_read2 / total_reads,
                "most_common_adapter_read1": most_common_adapter_read1,
                "most_common_adapter_read2": most_common_adapter_read2,
            }
        if total_samples < 1:
            return

        def count_dict_to_table_data(count_dict):
            table_data = {}
            for (adapter_sequence, adapter_match), adapter_count in count_dict.items():
                table_data[adapter_sequence] = {
                    "Adapter best match": adapter_match,
                    "Percentage of libraries with this adapter": adapter_count / total_samples,
                }
            return table_data

        headers = {
            "Adapter best match": {},
            "Percentage of libraries with this adapter": {"min": 0, "max": 1, "format": "{:.2%}"},
        }
        plot_config = {
            "id": "sequali_most_common_adapters_read_1",
            "title": "Sequali: Most common adapters for read 1",
            "defaultsort": [{"column": "Percentage of libraries with this adapter", "order": "desc"}],
        }

        self.add_section(
            name="Common Adapters: read 1",
            anchor="sequali_common_adapters_read1",
            plot=table.plot(count_dict_to_table_data(common_adapters_read1), pconfig=plot_config, headers=headers),
        )
        plot_config = {
            "id": "sequali_most_common_adapters_read_2",
            "title": "Sequali: Most common adapters for read 2",
            "defaultsort": [{"column": "Percentage of libraries with this adapter", "order": "desc"}],
        }

        self.add_section(
            name="Common Adapters: read 2",
            anchor="sequali_common_adapters_read2",
            plot=table.plot(count_dict_to_table_data(common_adapters_read2), pconfig=plot_config, headers=headers),
        )
        plot_config = {
            "id": "sequali_adapter_content_from_overlap_table",
            "title": "Sequali: Adapter Content calculated from overlap",
            "defaultsort": [
                {"column": "adapter_content_read1", "order": "desc"},
                {"column": "adapter_content_read2", "order": "desc"},
            ],
        }

        headers = {
            "adapter_content_read1": {
                "title": "Read 1: adapter content",
                "description": "The amount of reads that contain an adapter as based on "
                "the overlap between read 1 and read 2.",
                "scale": "Reds",
                "format": "{:.2%}",
                "max": 1,
                "min": 0,
            },
            "most_common_adapter_read1": {
                "title": "Most common adapter for read 1.",
                "description": "The most common adapter for read 1 as based on the overlap "
                "between read 1 and read 2.",
            },
            "adapter_content_read2": {
                "title": "Read 2: adapter content",
                "description": "The amount of reads that contain an adapter as based on "
                "the overlap between read 1 and read 2.",
                "scale": "Reds",
                "format": "{:.2%}",
                "max": 1,
                "min": 0,
            },
            "most_common_adapter_read2": {
                "title": "Most common adapter for read 2.",
                "description": "The most common adapter for read 1 as based on the overlap "
                "between read 1 and read 2.",
            },
        }
        self.add_section(
            name="Adapter Content",
            anchor="sequali_adapter_content_from_overlap",
            description="The cumulative percentage count of the found adapter sequences",
            plot=table.plot(adapter_table_data, headers=headers, pconfig=plot_config),
        )
