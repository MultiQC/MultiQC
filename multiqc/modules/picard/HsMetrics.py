""" MultiQC submodule to parse output from Picard HsMetrics """

import logging
import os
import re
from collections import OrderedDict, defaultdict

from multiqc import config
from multiqc.plots import linegraph, table

# Initialise the logger
log = logging.getLogger(__name__)

FIELD_DESCRIPTIONS = {
    "AT_DROPOUT": "A measure of how undercovered <= 50% GC regions are relative to the mean.",
    "BAIT_DESIGN_EFFICIENCY": "Target territory / bait territory. 1 == perfectly efficient, 0.5 = half of baited bases are not target.",
    "BAIT_SET": "The name of the bait set used in the hybrid selection.",
    "BAIT_TERRITORY": "The number of bases which have one or more baits on top of them.",
    "FOLD_80_BASE_PENALTY": 'The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to the mean coverage level in those targets.',
    "FOLD_ENRICHMENT": "The fold by which the baited region has been amplified above genomic background.",
    "GC_DROPOUT": "A measure of how undercovered >= 50% GC regions are relative to the mean.",
    "GENOME_SIZE": "The number of bases in the reference genome used for alignment.",
    "HET_SNP_Q": "The Phred Scaled Q Score of the theoretical HET SNP sensitivity.",
    "HET_SNP_SENSITIVITY": "The theoretical HET SNP sensitivity.",
    "HS_LIBRARY_SIZE": "The estimated number of unique molecules in the selected part of the library.",
    "HS_PENALTY_100X": "The 'hybrid selection penalty' incurred to get 80% of target bases to 100X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 100X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 100 * HS_PENALTY_100X.",
    "HS_PENALTY_10X": "The 'hybrid selection penalty' incurred to get 80% of target bases to 10X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 10X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 10 * HS_PENALTY_10X.",
    "HS_PENALTY_20X": "The 'hybrid selection penalty' incurred to get 80% of target bases to 20X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 20X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 20 * HS_PENALTY_20X.",
    "HS_PENALTY_30X": "The 'hybrid selection penalty' incurred to get 80% of target bases to 30X. This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 30X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 30 * HS_PENALTY_30X.",
    "HS_PENALTY_40X": "The 'hybrid selection penalty' incurred to get 80% of target bases to 40X.  This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 40X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 40 * HS_PENALTY_40X.",
    "HS_PENALTY_50X": "The 'hybrid selection penalty' incurred to get 80% of target bases to 50X.  This metric should be interpreted as: if I have a design with 10 megabases of target, and want to get 50X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 50 * HS_PENALTY_50X.",
    "MAX_TARGET_COVERAGE": "The maximum coverage of reads that mapped to target regions of an experiment.",
    "MEAN_BAIT_COVERAGE": "The mean coverage of all baits in the experiment.",
    "MEAN_TARGET_COVERAGE": "The mean coverage of targets.",
    "MEDIAN_TARGET_COVERAGE": "The median coverage of targets.",
    "NEAR_BAIT_BASES": "The number of PF aligned bases that mapped to within a fixed interval of a baited region, but not on a baited region.",
    "OFF_BAIT_BASES": "The number of PF aligned bases that mapped to neither on or near a bait.",
    "OLD_80_BASE_PENALTY": 'The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to the mean coverage level in those targets.',
    "ON_BAIT_BASES": "The number of PF aligned bases that mapped to a baited region of the genome.",
    "ON_BAIT_VS_SELECTED": "The percentage of on+near bait bases that are on as opposed to near.",
    "ON_TARGET_BASES": "The number of PF aligned bases that mapped to a targeted region of the genome.",
    "PCT_EXC_BASEQ": "The fraction of aligned bases that were filtered out because they were of low base quality.",
    "PCT_EXC_DUPE": "The fraction of aligned bases that were filtered out because they were in reads marked as duplicates.",
    "PCT_EXC_MAPQ": "The fraction of aligned bases that were filtered out because they were in reads with low mapping quality.",
    "PCT_EXC_OFF_TARGET": "The fraction of aligned bases that were filtered out because they did not align over a target base.",
    "PCT_EXC_OVERLAP": "The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads.",
    "PCT_OFF_BAIT": "The percentage of aligned PF bases that mapped neither on or near a bait.",
    "PCT_PF_READS": "PF reads / total reads. The percent of reads passing filter.",
    "PCT_PF_UQ_READS_ALIGNED": "PF Reads Aligned / PF Reads.",
    "PCT_PF_UQ_READS": "PF Unique Reads / Total Reads.",
    "PCT_SELECTED_BASES": "On+Near Bait Bases / PF Bases Aligned.",
    "PCT_TARGET_BASES_100X": "The fraction of all target bases achieving 100X or greater coverage.",
    "PCT_TARGET_BASES_10X": "The fraction of all target bases achieving 10X or greater coverage.",
    "PCT_TARGET_BASES_1X": "The fraction of all target bases achieving 1X or greater coverage.",
    "PCT_TARGET_BASES_20X": "The fraction of all target bases achieving 20X or greater coverage.",
    "PCT_TARGET_BASES_2X": "The fraction of all target bases achieving 2X or greater coverage.",
    "PCT_TARGET_BASES_30X": "The fraction of all target bases achieving 30X or greater coverage.",
    "PCT_TARGET_BASES_40X": "The fraction of all target bases achieving 40X or greater coverage.",
    "PCT_TARGET_BASES_50X": "The fraction of all target bases achieving 50X or greater coverage.",
    "PCT_USABLE_BASES_ON_BAIT": "The number of aligned, de-duped, on-bait bases out of the PF bases available.",
    "PCT_USABLE_BASES_ON_TARGET": "The number of aligned, de-duped, on-target bases out of the PF bases available.",
    "PF_BASES_ALIGNED": "The number of PF unique bases that are aligned with mapping score > 0 to the reference genome.",
    "PF_READS": "The number of reads that pass the vendor's filter.",
    "PF_UNIQUE_READS": "The number of PF reads that are not marked as duplicates.",
    "PF_UQ_BASES_ALIGNED": "The number of bases in the PF aligned reads that are mapped to a reference base. Accounts for clipping and gaps.",
    "PF_UQ_READS_ALIGNED": "The number of PF unique reads that are aligned with mapping score > 0 to the reference genome.",
    "TARGET_TERRITORY": "The unique number of target bases in the experiment where target is usually exons etc.",
    "TOTAL_READS": "The total number of reads in the SAM or BAM file examine.",
    "ZERO_CVG_TARGETS_PCT": "The fraction of targets that did not reach coverage=1 over any base.",
}


def parse_reports(self):
    """Find Picard HsMetrics reports and parse their data"""

    # Set up vars
    self.picard_HsMetrics_data = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files("picard/hsmetrics", filehandles=True):
        parsed_data = dict()
        s_name = None
        keys = None
        commadecimal = None
        for l in f["f"]:
            # New log starting
            if "CalculateHsMetrics" in l or "CollectHsMetrics" in l and "INPUT" in l:
                s_name = None
                keys = None

                # Pull sample name from input
                fn_search = re.search(r"INPUT(?:=|\s+)(\[?[^\s]+\]?)", l, flags=re.IGNORECASE)
                if fn_search:
                    s_name = os.path.basename(fn_search.group(1).strip("[]"))
                    s_name = self.clean_s_name(s_name, f)
                    parsed_data[s_name] = dict()

            if s_name is not None:
                if "HsMetrics" in l and "## METRICS CLASS" in l:
                    keys = f["f"].readline().strip("\n").split("\t")
                elif keys:
                    vals = l.strip("\n").split("\t")
                    if len(vals) == len(keys):
                        j = "NA"
                        if keys[0] == "BAIT_SET":
                            j = vals[0]
                        parsed_data[s_name][j] = dict()
                        # Check that we're not using commas for decimal places
                        if commadecimal is None:
                            commadecimal = False
                            for i, k in enumerate(keys):
                                if "PCT" in k or "BAIT" in k or "MEAN" in k:
                                    if "," in vals[i]:
                                        commadecimal = True
                                        break
                        for i, k in enumerate(keys):
                            try:
                                if commadecimal:
                                    vals[i] = vals[i].replace(".", "")
                                    vals[i] = vals[i].replace(",", ".")
                                parsed_data[s_name][j][k] = float(vals[i])
                            except ValueError:
                                parsed_data[s_name][j][k] = vals[i]
                    else:
                        s_name = None
                        keys = None

        # Remove empty dictionaries
        for s_name in list(parsed_data.keys()):
            for j in parsed_data[s_name].keys():
                if len(parsed_data[s_name][j]) == 0:
                    parsed_data[s_name].pop(j, None)
            if len(parsed_data[s_name]) == 0:
                parsed_data.pop(s_name, None)

        # Manipulate sample names if multiple baits found
        for s_name in parsed_data.keys():
            for j in parsed_data[s_name].keys():
                this_s_name = s_name
                if len(parsed_data[s_name]) > 1:
                    this_s_name = "{}: {}".format(s_name, j)
                if this_s_name in self.picard_HsMetrics_data:
                    log.debug("Duplicate sample name found in {}! Overwriting: {}".format(f["fn"], this_s_name))
                self.add_data_source(f, this_s_name, section="HsMetrics")
                self.picard_HsMetrics_data[this_s_name] = parsed_data[s_name][j]

    # Filter to strip out ignored sample names
    self.picard_HsMetrics_data = self.ignore_samples(self.picard_HsMetrics_data)

    if len(self.picard_HsMetrics_data) > 0:
        # Write parsed data to a file
        self.write_data_file(self.picard_HsMetrics_data, "multiqc_picard_HsMetrics")

        # Swap question marks with -1
        data = self.picard_HsMetrics_data
        for s_name in data:
            if data[s_name]["FOLD_ENRICHMENT"] == "?":
                data[s_name]["FOLD_ENRICHMENT"] = -1

        # Add to general stats table
        general_stats_table(self, data)

        # Add report section
        self.add_section(
            name="HSMetrics",
            anchor="picard_hsmetrics",
            plot=table.plot(
                data,
                _get_table_headers(),
                {
                    "id": "picard_hsmetrics_table",
                    "namespace": "HsMetrics",
                    "scale": "RdYlGn",
                    "min": 0,
                },
            ),
        )
        tbases = _add_target_bases(data)
        self.add_section(
            name=tbases["name"], anchor=tbases["anchor"], description=tbases["description"], plot=tbases["plot"]
        )
        hs_pen_plot = hs_penalty_plot(data)
        if hs_pen_plot is not None:
            self.add_section(
                name="HS Penalty",
                anchor="picard_hsmetrics_hs_penalty",
                description='The "hybrid selection penalty" incurred to get 80% of target bases to a given coverage.',
                helptext="""
                    Can be used with the following formula:

                    ```
                    required_aligned_bases = bait_size_bp * desired_coverage * hs_penalty
                    ```
                """,
                plot=hs_pen_plot,
            )

    # Return the number of detected samples to the parent module
    return len(self.picard_HsMetrics_data)


def general_stats_table(self, data):
    """
    Generate table header configs for the General Stats table,
    add config and data to the base module.
    """
    # Look for a user config of which table columns we should use
    picard_config = getattr(config, "picard_config", {})
    HsMetrics_genstats_table_cols = picard_config.get("HsMetrics_genstats_table_cols", [])
    HsMetrics_genstats_table_cols_hidden = picard_config.get("HsMetrics_genstats_table_cols_hidden", [])

    # Custom general stats columns
    if len(HsMetrics_genstats_table_cols) or len(HsMetrics_genstats_table_cols_hidden):
        for k, v in _generate_table_header_config(
            HsMetrics_genstats_table_cols, HsMetrics_genstats_table_cols_hidden
        ).items():
            self.general_stats_headers[k] = v

    # Default General Stats headers
    else:
        self.general_stats_headers["FOLD_ENRICHMENT"] = {
            "title": "Fold Enrichment",
            "min": 0,
            "format": "{:,.0f}",
            "scale": "Blues",
            "suffix": " X",
        }
        try:
            covs = picard_config["general_stats_target_coverage"]
            assert type(covs) == list
            assert len(covs) > 0
            covs = [str(i) for i in covs]
            log.debug("Custom Picard coverage thresholds: {}".format(", ".join([i for i in covs])))
        except (KeyError, AttributeError, TypeError, AssertionError):
            covs = ["30"]
        for c in covs:
            self.general_stats_headers["PCT_TARGET_BASES_{}X".format(c)] = {
                "id": "picard_target_bases_{}X".format(c),
                "title": "% Target Bases {}X".format(c),
                "description": "Percent of target bases with coverage &ge; {}X".format(c),
                "max": 100,
                "min": 0,
                "suffix": "%",
                "format": "{:,.0f}",
                "scale": "RdYlGn",
                "modify": lambda x: self.multiply_hundred(x),
            }

    # Add data to general stats table
    for s_name in data:
        if s_name not in self.general_stats_data:
            self.general_stats_data[s_name] = dict()
        self.general_stats_data[s_name].update(data[s_name])


def _get_table_headers():
    # Look for a user config of which table columns we should use
    picard_config = getattr(config, "picard_config", {})
    HsMetrics_table_cols = picard_config.get("HsMetrics_table_cols")
    HsMetrics_table_cols_hidden = picard_config.get("HsMetrics_table_cols_hidden")

    # Default table columns
    if not HsMetrics_table_cols:
        HsMetrics_table_cols = [
            "AT_DROPOUT",
            "BAIT_DESIGN_EFFICIENCY",
            "BAIT_TERRITORY",
            "FOLD_80_BASE_PENALTY",
            "FOLD_ENRICHMENT",
            "GC_DROPOUT",
            "HET_SNP_Q",
            "HET_SNP_SENSITIVITY",
            "NEAR_BAIT_BASES",
            "OFF_BAIT_BASES",
            "ON_BAIT_BASES",
            "ON_TARGET_BASES",
            "PCT_USABLE_BASES_ON_BAIT",
            "PCT_USABLE_BASES_ON_TARGET",
            "PF_BASES_ALIGNED",
            "PF_READS",
            "PCT_SELECTED_BASES",
            "PF_UNIQUE_READS",
            "PF_UQ_BASES_ALIGNED",
            "PF_UQ_READS_ALIGNED",
            "TOTAL_READS",
            "MAX_TARGET_COVERAGE",
            "MEAN_BAIT_COVERAGE",
            "MEAN_TARGET_COVERAGE",
            "MEDIAN_TARGET_COVERAGE",
            "ON_BAIT_VS_SELECTED",
            "TARGET_TERRITORY",
            "ZERO_CVG_TARGETS_PCT",
        ]
    if not HsMetrics_table_cols_hidden:
        HsMetrics_table_cols_hidden = ["BAIT_TERRITORY", "TOTAL_READS", "TARGET_TERRITORY", "AT_DROPOUT", "GC_DROPOUT"]

    return _generate_table_header_config(HsMetrics_table_cols, HsMetrics_table_cols_hidden)


def _generate_table_header_config(table_cols=[], hidden_table_cols=[]):
    """
    Automatically generate some nice table header configs based on what we know about
    the different types of Picard data fields.
    """
    title_cleanup = [
        ("CVG", "coverage"),
        ("UQ", "unique"),
        ("ON_", "On-"),
        ("OFF_", "Off-"),
        ("NEAR_", "Near-"),
        ("_", " "),
        ("PCT", ""),
    ]

    # Warn if we see anything unexpected
    for c in table_cols + hidden_table_cols:
        if c not in FIELD_DESCRIPTIONS and c[:17] != "PCT_TARGET_BASES_":
            log.error(f"Field '{c}' not found in expected Picard fields. Please check your config.")

    headers = OrderedDict()
    for h in table_cols + hidden_table_cols:
        # Set up the configuration for each column
        if h not in headers:
            # Generate a nice string for the column title
            h_title = h
            for s, r in title_cleanup:
                h_title = h_title.replace(s, r)

            headers[h] = {
                "title": h_title.strip().lower().capitalize(),
                "description": FIELD_DESCRIPTIONS[h] if h in FIELD_DESCRIPTIONS else None,
            }
            if h.find("PCT") > -1:
                headers[h]["title"] = "% {}".format(headers[h]["title"])
                headers[h]["modify"] = lambda x: x * 100.0
                headers[h]["max"] = 100
                headers[h]["suffix"] = "%"

            elif h.find("READS") > -1:
                headers[h]["title"] = "{} {}".format(config.read_count_prefix, headers[h]["title"])
                headers[h]["modify"] = lambda x: x * config.read_count_multiplier
                headers[h]["shared_key"] = "read_count"

            elif h.find("BASES") > -1:
                headers[h]["title"] = "{} {}".format(config.base_count_prefix, headers[h]["title"])
                headers[h]["modify"] = lambda x: x * config.base_count_multiplier
                headers[h]["shared_key"] = "base_count"

            # Manual capitilisation for some strings
            headers[h]["title"] = headers[h]["title"].replace("Pf", "PF").replace("snp", "SNP")

            if h in hidden_table_cols:
                headers[h]["hidden"] = True

    return headers


def _add_target_bases(data):
    data_clean = defaultdict(dict)
    for s in data:
        for h in data[s]:
            if h.startswith("PCT_TARGET"):
                data_clean[s][int(h.replace("PCT_TARGET_BASES_", "")[:-1])] = data[s][h] * 100.0

    pconfig = {
        "id": "picard_percentage_target_bases",
        "title": "Picard: Percentage of target bases",
        "xlab": "Fold Coverage",
        "ylab": "Pct of bases",
        "ymax": 100,
        "ymin": 0,
        "xmin": 0,
        "tt_label": "<b>{point.x}X</b>: {point.y:.2f}%",
    }
    return {
        "name": "Target Region Coverage",
        "anchor": "picard_hsmetrics_target_bases",
        "description": "The percentage of all target bases with at least <code>x</code> fold coverage.",
        "plot": linegraph.plot(data_clean, pconfig),
    }


def hs_penalty_plot(data):
    data_clean = defaultdict(dict)
    any_non_zero = False
    for s in data:
        for h in data[s]:
            if h.startswith("HS_PENALTY"):
                data_clean[s][int(h.lstrip("HS_PENALTY_").rstrip("X"))] = data[s][h]
                if data[s][h] > 0:
                    any_non_zero = True

    pconfig = {
        "id": "picard_hybrid_selection_penalty",
        "title": "Picard: Hybrid Selection Penalty",
        "xlab": "Fold Coverage",
        "ylab": "Penalty",
        "ymin": 0,
        "xmin": 0,
        "xDecimals": False,
        "tt_label": "<b>{point.x}X</b>: {point.y:.2f}%",
    }

    if any_non_zero:
        return linegraph.plot(data_clean, pconfig)
