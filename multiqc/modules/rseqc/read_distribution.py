""" MultiQC submodule to parse output from RSeQC read_distribution.py
http://rseqc.sourceforge.net/#read-distribution-py """

import logging
import re

from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(self):
    """Find RSeQC read_distribution reports and parse their data"""

    # Set up vars
    self.read_dist = dict()
    first_regexes = {
        "total_reads": r"Total Reads\s+(\d+)\s*",
        "total_tags": r"Total Tags\s+(\d+)\s*",
        "total_assigned_tags": r"Total Assigned Tags\s+(\d+)\s*",
    }
    second_regexes = {
        "cds_exons": r"CDS_Exons\s+(\d+)\s+(\d+)\s+([\d\.]+)\s*",
        "5_utr_exons": r"5'UTR_Exons\s+(\d+)\s+(\d+)\s+([\d\.]+)\s*",
        "3_utr_exons": r"3'UTR_Exons\s+(\d+)\s+(\d+)\s+([\d\.]+)\s*",
        "introns": r"Introns\s+(\d+)\s+(\d+)\s+([\d\.]+)\s*",
        "tss_up_1kb": r"TSS_up_1kb\s+(\d+)\s+(\d+)\s+([\d\.]+)\s*",
        "tss_up_5kb": r"TSS_up_5kb\s+(\d+)\s+(\d+)\s+([\d\.]+)\s*",
        "tss_up_10kb": r"TSS_up_10kb\s+(\d+)\s+(\d+)\s+([\d\.]+)\s*",
        "tes_down_1kb": r"TES_down_1kb\s+(\d+)\s+(\d+)\s+([\d\.]+)\s*",
        "tes_down_5kb": r"TES_down_5kb\s+(\d+)\s+(\d+)\s+([\d\.]+)\s*",
        "tes_down_10kb": r"TES_down_10kb\s+(\d+)\s+(\d+)\s+([\d\.]+)\s*",
    }

    # Go through files and parse data using regexes
    for f in self.find_log_files("rseqc/read_distribution"):
        d = dict()
        for k, r in first_regexes.items():
            r_search = re.search(r, f["f"], re.MULTILINE)
            if r_search:
                d[k] = int(r_search.group(1))
        for k, r in second_regexes.items():
            r_search = re.search(r, f["f"], re.MULTILINE)
            if r_search:
                d[f"{k}_total_bases"] = int(r_search.group(1))
                d[f"{k}_tag_count"] = int(r_search.group(2))
                d[f"{k}_tags_kb"] = float(r_search.group(3))

        d["other_intergenic_tag_count"] = d["total_tags"] - d["total_assigned_tags"]

        # Calculate some percentages for parsed file
        if "total_tags" in d:
            t = float(d["total_tags"])
            pcts = dict()
            for k in d:
                if k.endswith("_tag_count"):
                    pk = f"{k[:-10]}_tag_pct"
                    try:
                        pcts[pk] = (float(d[k]) / t) * 100.0
                    except ZeroDivisionError:
                        pcts[pk] = 0
            d.update(pcts)

        if len(d) > 0:
            if f["s_name"] in self.read_dist:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            self.add_data_source(f, section="read_distribution")
            self.read_dist[f["s_name"]] = d

    # Filter to strip out ignored sample names
    self.read_dist = self.ignore_samples(self.read_dist)

    if len(self.read_dist) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    self.add_software_version(None)

    # Fix double counting of TSS_up and TES_down
    prefixes = ("tss_up", "tes_down")
    sizes = ("10kb", "5kb", "1kb")
    suffixes = ("total_bases", "tag_count")
    combinations = [
        (prefix, big, small, suffix)
        for prefix in prefixes
        for suffix in suffixes
        for big, small in zip(sizes, sizes[1:])
    ]

    for sample_name, sample in self.read_dist.items():
        for prefix, big, small, suffix in combinations:
            try:
                big_val = sample.pop(f"{prefix}_{big}_{suffix}")
                small_val = sample[f"{prefix}_{small}_{suffix}"]
                sample[f"{prefix}_{small}-{big}_{suffix}"] = big_val - small_val
            except KeyError as e:
                log.debug(f"Field {e} is missing from {sample_name}")

    # Write to file
    self.write_data_file(self.read_dist, "multiqc_rseqc_read_distribution")

    # Plot bar graph of groups
    keys = {
        "cds_exons_tag_count": {"name": "CDS_Exons"},
        "5_utr_exons_tag_count": {"name": "5'UTR_Exons"},
        "3_utr_exons_tag_count": {"name": "3'UTR_Exons"},
        "introns_tag_count": {"name": "Introns"},
        "tss_up_1kb_tag_count": {"name": "TSS_up_1kb"},
        "tss_up_1kb-5kb_tag_count": {"name": "TSS_up_1kb-5kb"},
        "tss_up_5kb-10kb_tag_count": {"name": "TSS_up_5kb-10kb"},
        "tes_down_1kb_tag_count": {"name": "TES_down_1kb"},
        "tes_down_1kb-5kb_tag_count": {"name": "TES_down_1kb-5kb"},
        "tes_down_5kb-10kb_tag_count": {"name": "TES_down_5kb-10kb"},
        "other_intergenic_tag_count": {"name": "Other_intergenic"},
    }

    # Config for the plot
    pconfig = {
        "id": "rseqc_read_distribution_plot",
        "title": "RSeQC: Read Distribution",
        "ylab": "# Tags",
        "cpswitch_counts_label": "Number of Tags",
        "cpswitch_c_active": False,
    }

    if max([d["total_tags"] for d in self.read_dist.values()]) > 0:
        self.add_section(
            name="Read Distribution",
            anchor="rseqc-read_distribution",
            description='<a href="http://rseqc.sourceforge.net/#read-distribution-py" target="_blank">Read Distribution</a>'
            " calculates how mapped reads are distributed over genome features.",
            plot=bargraph.plot(self.read_dist, keys, pconfig),
        )
    else:
        self.add_section(
            name="Read Distribution",
            anchor="rseqc-read_distribution",
            description='<a href="http://rseqc.sourceforge.net/#read-distribution-py" target="_blank">Read Distribution</a>'
            " calculates how mapped reads are distributed over genome features.",
            content='<div class="alert alert-info">All samples had zero tags.</div>',
        )

    # Return number of samples found
    return len(self.read_dist)
