""" MultiQC submodule to parse output from Samtools stats """

import logging
import re

from multiqc import config
from multiqc.plots import bargraph, violin

# Initialise the logger
log = logging.getLogger(__name__)

# Regex to grab version number from samtools stats contents
VERSION_REGEX = r"# This file was produced by samtools stats \(([\d\.]+)"
HTSLIB_REGEX = r"\+htslib-([\d\.]+)"


class StatsReportMixin:
    """Mixin class, loaded by main samtools MultiqcModule class."""

    def parse_samtools_stats(self):
        """Find Samtools stats logs and parse their data"""

        self.samtools_stats = dict()
        for f in self.find_log_files("samtools/stats"):
            parsed_data = dict()
            for line in f["f"].splitlines():
                # Get version number from file contents
                if line.startswith("# This file was produced by samtools stats"):
                    # Look for Samtools version
                    version_match = re.search(VERSION_REGEX, line)
                    if version_match is None:
                        continue

                    # Add Samtools version
                    samtools_version = version_match.group(1)
                    self.add_software_version(samtools_version, f["s_name"])

                    # Look for HTSlib version
                    htslib_version_match = re.search(HTSLIB_REGEX, line)
                    if htslib_version_match is None:
                        continue

                    # Add HTSlib version if different from Samtools version
                    htslib_version = htslib_version_match.group(1)
                    if htslib_version != samtools_version:
                        self.add_software_version(htslib_version, f["s_name"], "HTSlib")

                if not line.startswith("SN"):
                    continue
                sections = line.split("\t")
                field = sections[1].strip()[:-1]
                field = field.replace(" ", "_")
                value = float(sections[2].strip())
                parsed_data[field] = value

            if len(parsed_data) > 0:
                # Work out some percentages
                if "raw_total_sequences" in parsed_data:
                    for k in list(parsed_data.keys()):
                        if (
                            k.startswith("reads_")
                            and k != "raw_total_sequences"
                            and parsed_data["raw_total_sequences"] > 0
                        ):
                            parsed_data[f"{k}_percent"] = (parsed_data[k] / parsed_data["raw_total_sequences"]) * 100

                if f["s_name"] in self.samtools_stats:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                self.add_data_source(f, section="stats")
                self.samtools_stats[f["s_name"]] = parsed_data

        # Filter to strip out ignored sample names
        self.samtools_stats = self.ignore_samples(self.samtools_stats)

        if len(self.samtools_stats) == 0:
            return 0

        # Write parsed report data to a file
        self.write_data_file(self.samtools_stats, "multiqc_samtools_stats")

        # General Stats Table
        stats_headers = {
            "error_rate": {
                "title": "Error rate",
                "description": "Error rate: mismatches (NM) / bases mapped (CIGAR)",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "OrRd",
                "format": "{:,.2f}",
                "modify": lambda x: x * 100.0,
            },
            "non-primary_alignments": {
                "title": "Non-primary",
                "description": f"Non-primary alignments ({config.read_count_desc})",
                "scale": "PuBu",
                "shared_key": "read_count",
            },
            "reads_mapped": {
                "title": "Reads mapped",
                "description": f"Reads mapped in the bam file ({config.read_count_desc})",
                "shared_key": "read_count",
            },
            "reads_mapped_percent": {
                "title": "% Mapped",
                "description": "% Mapped reads",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "reads_properly_paired_percent": {
                "title": "% Proper pairs",
                "description": "% Properly paired reads",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
                "hidden": True
                if (max([x["reads_mapped_and_paired"] for x in self.samtools_stats.values()]) == 0)
                else False,
            },
            "reads_MQ0_percent": {
                "title": "% MapQ 0 reads",
                "description": "% of reads that are ambiguously placed (MapQ=0)",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "OrRd",
                "hidden": True,
            },
            "raw_total_sequences": {
                "title": "Total seqs",
                "description": f"Total sequences in the bam file ({config.read_count_desc})",
                "shared_key": "read_count",
            },
        }
        self.general_stats_addcols(self.samtools_stats, stats_headers, namespace="stats")

        # Make bargraph plot of mapped/unmapped reads
        self.alignment_section(self.samtools_stats)

        # Make dot plot of counts
        keys = {}
        reads = {
            "min": 0,
            "modify": lambda x: float(x) * config.read_count_multiplier,
            "suffix": config.read_count_prefix,
            "tt_decimals": 2,
            "shared_key": "read_count",
        }
        bases = {
            "min": 0,
            "modify": lambda x: float(x) * config.base_count_multiplier,
            "suffix": config.base_count_prefix,
            "tt_decimals": 2,
            "shared_key": "base_count",
        }
        keys["raw_total_sequences"] = dict(reads, **{"title": "Total sequences"})
        keys["reads_mapped_and_paired"] = dict(
            reads,
            **{"title": "Mapped &amp; paired", "description": "Paired-end technology bit set + both mates mapped"},
        )
        keys["reads_properly_paired"] = dict(
            reads, **{"title": "Properly paired", "description": "Proper-pair bit set"}
        )
        keys["reads_duplicated"] = dict(
            reads, **{"title": "Duplicated", "description": "PCR or optical duplicate bit set"}
        )
        keys["reads_QC_failed"] = dict(reads, **{"title": "QC Failed"})
        keys["reads_MQ0"] = dict(reads, **{"title": "Reads MQ0", "description": "Reads mapped and MQ=0"})
        keys["bases_mapped_(cigar)"] = dict(
            bases, **{"title": "Mapped bases (CIGAR)", "description": "Mapped bases (CIGAR)"}
        )
        keys["bases_trimmed"] = dict(bases, **{"title": "Bases Trimmed"})
        keys["bases_duplicated"] = dict(bases, **{"title": "Duplicated bases"})
        keys["pairs_on_different_chromosomes"] = dict(
            reads, **{"title": "Diff chromosomes", "description": "Pairs on different chromosomes"}
        )
        keys["pairs_with_other_orientation"] = dict(
            reads, **{"title": "Other orientation", "description": "Pairs with other orientation"}
        )
        keys["inward_oriented_pairs"] = dict(reads, **{"title": "Inward pairs", "description": "Inward oriented pairs"})
        keys["outward_oriented_pairs"] = dict(
            reads, **{"title": "Outward pairs", "description": "Outward oriented pairs"}
        )

        self.add_section(
            name="Alignment stats",
            anchor="samtools-stats",
            description="This module parses the output from <code>samtools stats</code>. All numbers in millions.",
            plot=violin.plot(
                self.samtools_stats,
                keys,
                {
                    "id": "samtools-stats-dp",
                    "title": "Samtools stats: Alignment Stats",
                },
            ),
        )

        # Return the number of logs that were found
        return len(self.samtools_stats)

    def alignment_section(self, samples_data):
        bedgraph_data = {}
        for sample_id, data in samples_data.items():
            # Breaking up the mapped reads count into MQ0 and >MQ1 counts
            data["reads_mapped_MQ1"] = data["reads_mapped"] - data["reads_MQ0"]
            # Asserting the bar plot keys sum up to the total
            expected_total = data["raw_total_sequences"]
            read_sum = (
                data["reads_mapped_MQ1"] + data["reads_MQ0"] + data["reads_unmapped"] + data["filtered_sequences"]
            )
            if read_sum == expected_total:
                bedgraph_data[sample_id] = data
            else:
                log.warning(
                    "sum of mapped/unmapped/filtered reads not matching total, "
                    "skipping samtools plot for: {}".format(sample_id)
                )
        self.add_section(
            name="Percent mapped",
            anchor="samtools-stats-alignment",
            description="Alignment metrics from <code>samtools stats</code>; mapped vs. unmapped reads vs. reads mapped with MQ0.",
            helptext="""
            For a set of samples that have come from the same multiplexed library,
            similar numbers of reads for each sample are expected. Large differences in numbers might
            indicate issues during the library preparation process. Whilst large differences in read
            numbers may be controlled for in downstream processings (e.g. read count normalisation),
            you may wish to consider whether the read depths achieved have fallen below recommended
            levels depending on the applications.

            Low alignment rates could indicate contamination of samples (e.g. adapter sequences),
            low sequencing quality or other artefacts. These can be further investigated in the
            sequence level QC (e.g. from FastQC).

            Reads mapped with MQ0 often indicate that the reads are ambiguously mapped to multiple
            locations in the reference sequence. This can be due to repetitive regions in the genome,
            the presence of alternative contigs in the reference, or due to reads that are too short
            to be uniquely mapped. These reads are often filtered out in downstream analyses.
            """,
            plot=alignment_chart(bedgraph_data),
        )


def alignment_chart(data):
    keys = {
        "reads_mapped_MQ1": {"color": "#437bb1", "name": "Mapped (with MQ>0)"},
        "reads_MQ0": {"color": "#FF9933", "name": "MQ0"},
        "reads_unmapped": {"color": "#b1084c", "name": "Unmapped"},
    }

    # Config for the plot
    plot_conf = {
        "id": "samtools_alignment_plot",
        "title": "Samtools stats: Alignment Scores",
        "ylab": "# Reads",
        "cpswitch_counts_label": "Number of Reads",
    }
    return bargraph.plot(data, keys, plot_conf)
