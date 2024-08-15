import logging
import os
import re
from copy import deepcopy
from typing import Dict, Optional

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    This module parses summary statistics from the `Log.final.out` log files.
    Sample names are taken either from the filename prefix (`sampleNameLog.final.out`)
    when set with `--outFileNamePrefix` in STAR. If there is no filename prefix,
    the sample name is set as the name of the directory containing the file.

    In addition to this summary log file, the module parses `ReadsPerGene.out.tab`
    files generated with `--quantMode GeneCounts`, if found.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="STAR",
            anchor="star",
            href="https://github.com/alexdobin/STAR",
            info="Universal RNA-seq aligner.",
            doi="10.1093/bioinformatics/bts635",
        )

        # Find and load any STAR reports
        data_by_sample: Dict[str, Dict[str, float]] = dict()
        for f in self.find_log_files("star"):
            parsed_data = parse_star_report(f["f"])
            if parsed_data is not None:
                s_name = f["s_name"]
                if s_name == "" or s_name == "Log.final.out":
                    s_name = self.clean_s_name(os.path.basename(f["root"]), f, root=os.path.dirname(f["root"]))
                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(f, section="SummaryLog")
                data_by_sample[s_name] = parsed_data

        # Find and load any STAR gene count tables
        genecounts_unstranded: Dict[str, Dict[str, int]] = dict()
        genecounts_first_strand: Dict[str, Dict[str, int]] = dict()
        genecounts_second_strand: Dict[str, Dict[str, int]] = dict()
        for f in self.find_log_files("star/genecounts", filehandles=True):
            gc_parsed_data = parse_star_genecounts_report(f)
            if gc_parsed_data is not None:
                s_name = f["s_name"]
                if s_name == "" or s_name == "ReadsPerGene.out.tab":
                    s_name = self.clean_s_name(os.path.basename(f["root"]), f, root=os.path.dirname(f["root"]))
                if s_name in genecounts_unstranded:
                    log.debug(f"Duplicate ReadsPerGene sample name found! Overwriting: {s_name}")
                self.add_data_source(f, section="ReadsPerGene")
                genecounts_unstranded[s_name] = gc_parsed_data["unstranded"]
                genecounts_first_strand[s_name] = gc_parsed_data["first_strand"]
                genecounts_second_strand[s_name] = gc_parsed_data["second_strand"]

        # Filter to strip out ignored sample names
        data_by_sample = self.ignore_samples(data_by_sample)
        genecounts_unstranded = self.ignore_samples(genecounts_unstranded)
        genecounts_first_strand = self.ignore_samples(genecounts_first_strand)
        genecounts_second_strand = self.ignore_samples(genecounts_second_strand)
        if len(data_by_sample) == 0 and len(genecounts_unstranded) == 0:
            raise ModuleNoSamplesFound

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        if len(data_by_sample) > 0:
            if len(genecounts_unstranded) > 0:
                log.info(
                    "Found {} reports and {} gene count files".format(len(data_by_sample), len(genecounts_unstranded))
                )
            else:
                log.info(f"Found {len(data_by_sample)} reports")
        else:
            log.info(f"Found {len(genecounts_unstranded)} gene count files")

        if len(data_by_sample) > 0:
            # Write parsed report data to a file
            self.write_data_file(data_by_sample, "multiqc_star")

            # Basic Stats Table
            self.stats_tables(data_by_sample)

            # Alignment bar plot
            self.add_section(
                name="Alignment Scores", anchor="star_alignments", plot=star_alignment_chart(data_by_sample)
            )

        if len(genecounts_unstranded) > 0:
            self.add_section(
                name="Gene Counts",
                anchor="star_geneCounts",
                description="Statistics from results generated using <code>--quantMode GeneCounts</code>. "
                + "The three tabs show counts for unstranded RNA-seq, counts for the 1st read strand "
                + "aligned with RNA and counts for the 2nd read strand aligned with RNA.",
                plot=star_genecounts_chart(genecounts_unstranded, genecounts_first_strand, genecounts_second_strand),
            )

    def stats_tables(self, data_by_sample):
        """Take the parsed stats from the STAR report and add them to the
        basic stats table at the top of the report"""

        headers: Dict[str, Dict] = {
            "total_reads": {
                "title": "Total Reads",
                "description": "Number of input reads",
                "scale": "Blues",
                "shared_key": "read_count",
                "hidden": True,
            },
            "mapped": {
                "title": "Aligned",
                "description": "Mapped reads",
                "scale": "PuRd",
                "shared_key": "read_count",
                "hidden": True,
            },
            "mapped_percent": {
                "title": "Aligned",
                "description": "% Mapped reads",
                "suffix": "%",
                "scale": "PuRd",
            },
            "uniquely_mapped": {
                "title": "Uniq aligned",
                "description": "Uniquely mapped reads",
                "scale": "YlGn",
                "shared_key": "read_count",
                "hidden": True,
            },
            "uniquely_mapped_percent": {
                "title": "Uniq aligned",
                "description": "% Uniquely mapped reads",
                "suffix": "%",
                "scale": "YlGn",
            },
            "multimapped": {
                "title": "Multimapped",
                "description": "Multiple mapped reads",
                "scale": "PuRd",
                "shared_key": "read_count",
                "hidden": True,
            },
        }
        self.general_stats_addcols(data_by_sample, headers)

        all_headers: Dict[str, Dict] = deepcopy(headers)
        all_headers["total_reads"]["hidden"] = False

        all_headers.update(
            {
                "avg_input_read_length": {
                    "title": "Avg. read len",
                    "description": "Average input read length",
                    "suffix": "bp",
                    "scale": "Blues",
                    "hidden": True,
                },
                "avg_mapped_read_length": {
                    "title": "Avg. mapped len",
                    "description": "Average mapped length",
                    "suffix": "bp",
                    "scale": "Blues",
                },
                "num_splices": {
                    "title": "Splices",
                    "description": "Number of splices: Total",
                    "scale": "Blues",
                    "shared_key": "read_count",
                    "hidden": True,
                },
                "num_annotated_splices": {
                    "title": "Annotated splices",
                    "description": "Number of splices: Annotated (sjdb)",
                    "scale": "Blues",
                    "shared_key": "read_count",
                },
                "num_GTAG_splices": {
                    "title": "GT/AG splices",
                    "description": "Number of splices: GT/AG",
                    "scale": "Blues",
                    "shared_key": "read_count",
                    "hidden": True,
                },
                "num_GCAG_splices": {
                    "title": "GC/AG splices",
                    "description": "Number of splices: GC/AG",
                    "scale": "Blues",
                    "shared_key": "read_count",
                    "hidden": True,
                },
                "num_ATAC_splices": {
                    "title": "AT/AC splices",
                    "description": "Number of splices: AT/AC",
                    "scale": "Blues",
                    "shared_key": "read_count",
                    "hidden": True,
                },
                "num_noncanonical_splices": {
                    "title": "Non-canonical splices",
                    "description": "Number of splices: Non-canonical",
                    "scale": "Blues",
                    "shared_key": "read_count",
                    "hidden": True,
                },
                "mismatch_rate": {
                    "title": "Mismatch rate",
                    "description": "Mismatch rate per base",
                    "suffix": "%",
                    "scale": "Blues",
                },
                "deletion_rate": {
                    "title": "Del rate",
                    "description": "Deletion rate per base",
                    "suffix": "%",
                    "scale": "Blues",
                },
                "deletion_length": {
                    "title": "Del len",
                    "description": "Deletion average length",
                    "suffix": "bp",
                    "scale": "Blues",
                },
                "insertion_rate": {
                    "title": "Ins rate",
                    "description": "Insertion rate per base",
                    "suffix": "%",
                    "scale": "Blues",
                },
                "insertion_length": {
                    "title": "Ins len",
                    "description": "Insertion average length",
                    "suffix": "bp",
                    "scale": "Blues",
                },
            }
        )

        self.add_section(
            name="Summary Statistics",
            anchor="star_summary",
            description="Summary statistics from the STAR alignment",
            plot=table.plot(
                data_by_sample,
                all_headers,
                pconfig={"id": "star_summary_table", "title": "STAR: Summary Statistics", "namespace": "STAR"},
            ),
        )


def parse_star_genecounts_report(f) -> Optional[Dict[str, Dict[str, int]]]:
    """Parse a STAR gene counts output file"""
    # Three numeric columns: unstranded, stranded/first-strand, stranded/second-strand
    keys = ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"]
    unstranded: Dict[str, int] = {"N_genes": 0}
    first_strand: Dict[str, int] = {"N_genes": 0}
    second_strand: Dict[str, int] = {"N_genes": 0}
    num_errors = 0
    num_genes = 0
    for line in f["f"]:
        s = line.split("\t")
        try:
            for i in [1, 2, 3]:
                s[i] = float(s[i])
            if s[0] in keys:
                unstranded[s[0]] = s[1]
                first_strand[s[0]] = s[2]
                second_strand[s[0]] = s[3]
            else:
                unstranded["N_genes"] += s[1]
                first_strand["N_genes"] += s[2]
                second_strand["N_genes"] += s[3]
                num_genes += 1
        except IndexError:
            # Tolerate a few errors in case there is something random added at the top of the file
            num_errors += 1
            if num_errors > 10 and num_genes == 0:
                log.warning(f"Error parsing {f['fn']}")
                return None
    if num_genes > 0:
        return {"unstranded": unstranded, "first_strand": first_strand, "second_strand": second_strand}
    else:
        return None


def star_genecounts_chart(genecounts_unstranded, genecounts_first_strand, genecounts_second_strand):
    """Make a plot for the ReadsPerGene output"""

    # Specify the order of the different possible categories
    keys = {
        "N_genes": {"color": "#2f7ed8", "name": "Overlapping Genes"},
        "N_noFeature": {"color": "#0d233a", "name": "No Feature"},
        "N_ambiguous": {"color": "#492970", "name": "Ambiguous Features"},
        "N_multimapping": {"color": "#f28f43", "name": "Multimapping"},
        "N_unmapped": {"color": "#7f0000", "name": "Unmapped"},
    }

    # Config for the plot
    pconfig = {
        "id": "star_gene_counts",
        "title": "STAR: Gene Counts",
        "ylab": "# Reads",
        "cpswitch_counts_label": "Number of Reads",
        "data_labels": ["Unstranded", "Same Stranded", "Reverse Stranded"],
    }
    datasets = [
        genecounts_unstranded,
        genecounts_first_strand,
        genecounts_second_strand,
    ]
    return bargraph.plot(datasets, [keys, keys, keys], pconfig)


def star_alignment_chart(data_by_sample):
    """Make the plot showing alignment rates"""

    # Specify the order of the different possible categories
    keys = {
        "uniquely_mapped": {"color": "#437bb1", "name": "Uniquely mapped"},
        "multimapped": {"color": "#7cb5ec", "name": "Mapped to multiple loci"},
        "multimapped_toomany": {"color": "#f7a35c", "name": "Mapped to too many loci"},
        "unmapped_mismatches": {"color": "#e63491", "name": "Unmapped: too many mismatches"},
        "unmapped_tooshort": {"color": "#b1084c", "name": "Unmapped: too short"},
        "unmapped_other": {"color": "#7f0000", "name": "Unmapped: other"},
    }

    # Config for the plot
    pconfig = {
        "id": "star_alignment_plot",
        "title": "STAR: Alignment Scores",
        "ylab": "# Reads",
        "cpswitch_counts_label": "Number of Reads",
    }

    return bargraph.plot(data_by_sample, keys, pconfig)


def parse_star_report(contents: str) -> Optional[Dict[str, float]]:
    """Parse the final STAR log file."""

    regexes = {
        "total_reads": r"Number of input reads \|\s+(\d+)",
        "avg_input_read_length": r"Average input read length \|\s+([\d\.]+)",
        "uniquely_mapped": r"Uniquely mapped reads number \|\s+(\d+)",
        "uniquely_mapped_percent": r"Uniquely mapped reads % \|\s+([\d\.]+)",
        "avg_mapped_read_length": r"Average mapped length \|\s+([\d\.]+)",
        "num_splices": r"Number of splices: Total \|\s+(\d+)",
        "num_annotated_splices": r"Number of splices: Annotated \(sjdb\) \|\s+(\d+)",
        "num_GTAG_splices": r"Number of splices: GT/AG \|\s+(\d+)",
        "num_GCAG_splices": r"Number of splices: GC/AG \|\s+(\d+)",
        "num_ATAC_splices": r"Number of splices: AT/AC \|\s+(\d+)",
        "num_noncanonical_splices": r"Number of splices: Non-canonical \|\s+(\d+)",
        "mismatch_rate": r"Mismatch rate per base, % \|\s+([\d\.]+)",
        "deletion_rate": r"Deletion rate per base \|\s+([\d\.]+)",
        "deletion_length": r"Deletion average length \|\s+([\d\.]+)",
        "insertion_rate": r"Insertion rate per base \|\s+([\d\.]+)",
        "insertion_length": r"Insertion average length \|\s+([\d\.]+)",
        "multimapped": r"Number of reads mapped to multiple loci \|\s+(\d+)",
        "multimapped_percent": r"% of reads mapped to multiple loci \|\s+([\d\.]+)",
        "multimapped_toomany": r"Number of reads mapped to too many loci \|\s+(\d+)",
        "multimapped_toomany_percent": r"% of reads mapped to too many loci \|\s+([\d\.]+)",
        "unmapped_mismatches_percent": r"% of reads unmapped: too many mismatches \|\s+([\d\.]+)",
        "unmapped_tooshort_percent": r"% of reads unmapped: too short \|\s+([\d\.]+)",
        "unmapped_other_percent": r"% of reads unmapped: other \|\s+([\d\.]+)",
    }
    parsed_data: Dict[str, float] = dict()
    for k, r in regexes.items():
        r_search = re.search(r, contents, re.MULTILINE)
        if r_search:
            parsed_data[k] = float(r_search.group(1))
    # Figure out the numbers for unmapped as for some reason only the percentages are given
    try:
        total_mapped = parsed_data["uniquely_mapped"] + parsed_data["multimapped"] + parsed_data["multimapped_toomany"]
        unmapped_count = parsed_data["total_reads"] - total_mapped
        total_unmapped_percent = (
            parsed_data["unmapped_mismatches_percent"]
            + parsed_data["unmapped_tooshort_percent"]
            + parsed_data["unmapped_other_percent"]
        )
        try:
            parsed_data["unmapped_mismatches"] = int(
                round(unmapped_count * (parsed_data["unmapped_mismatches_percent"] / total_unmapped_percent), 0)
            )
            parsed_data["unmapped_tooshort"] = int(
                round(unmapped_count * (parsed_data["unmapped_tooshort_percent"] / total_unmapped_percent), 0)
            )
            parsed_data["unmapped_other"] = int(
                round(unmapped_count * (parsed_data["unmapped_other_percent"] / total_unmapped_percent), 0)
            )
        except ZeroDivisionError:
            parsed_data["unmapped_mismatches"] = 0
            parsed_data["unmapped_tooshort"] = 0
            parsed_data["unmapped_other"] = 0
    except KeyError:
        pass

    if len(parsed_data) == 0:
        return None

    parsed_data["mapped_percent"] = parsed_data["uniquely_mapped_percent"] + parsed_data["multimapped_percent"]
    parsed_data["mapped"] = parsed_data["uniquely_mapped"] + parsed_data["multimapped"]

    return parsed_data
