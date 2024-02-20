""" MultiQC submodule to parse output from Picard HsMetrics """

import logging
from collections import defaultdict

from multiqc import config
from multiqc.modules.picard import util
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


def parse_reports(module):
    """Find Picard HsMetrics reports and parse their data"""

    data_by_bait_by_sample = dict()

    # Go through logs and find Metrics
    for f in module.find_log_files("picard/hsmetrics", filehandles=True):
        s_name = f["s_name"]
        keys = None
        commadecimal = None

        for line in f["f"]:
            maybe_s_name = util.extract_sample_name(
                module,
                line,
                f,
                picard_tool="CollectHsMetrics",
                sentieon_algo="HsMetricAlgo",
            )
            if maybe_s_name:
                s_name = maybe_s_name
                keys = None

            if s_name is None:
                continue

            if util.is_line_right_before_table(line, picard_class="HsMetrics", sentieon_algo="HsMetricAlgo"):
                keys = f["f"].readline().strip("\n").split("\t")
                if s_name in data_by_bait_by_sample:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
                data_by_bait_by_sample[s_name] = dict()

            elif keys:
                vals = line.strip("\n").split("\t")
                if len(vals) != len(keys):
                    keys = None
                    continue

                bait = "NA"
                if keys[0] == "BAIT_SET":
                    bait = vals[0]
                data_by_bait_by_sample[s_name][bait] = dict()
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
                        data_by_bait_by_sample[s_name][bait][k] = float(vals[i])
                    except ValueError:
                        data_by_bait_by_sample[s_name][bait][k] = vals[i]

    # Remove empty dictionaries
    for s_name in data_by_bait_by_sample:
        for bait in data_by_bait_by_sample[s_name]:
            if len(data_by_bait_by_sample[s_name][bait]) == 0:
                data_by_bait_by_sample[s_name].pop(bait, None)
        if len(data_by_bait_by_sample[s_name]) == 0:
            data_by_bait_by_sample.pop(s_name, None)

    data_by_sample = dict()
    # Manipulate sample names if multiple baits found
    for s_name in data_by_bait_by_sample:
        for bait in data_by_bait_by_sample[s_name]:
            s_bait_name = s_name
            # If there are multiple baits, append the bait name to the sample name
            if len(data_by_bait_by_sample[s_name]) > 1:
                s_bait_name = f"{s_name}: {bait}"
            if s_bait_name in data_by_sample:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_bait_name}")
            data_by_sample[s_bait_name] = data_by_bait_by_sample[s_name][bait]
            module.add_data_source(f, s_bait_name, section="HsMetrics")

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)
    if len(data_by_sample) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed data to a file
    module.write_data_file(data_by_sample, f"multiqc_{module.anchor}_HsMetrics")

    # Swap question marks with -1
    for s_name in data_by_sample:
        if data_by_sample[s_name]["FOLD_ENRICHMENT"] == "?":
            data_by_sample[s_name]["FOLD_ENRICHMENT"] = -1

    # Add to general stats table
    _general_stats_table(module, data_by_sample)

    # Add report section
    module.add_section(
        name="HSMetrics",
        anchor=f"{module.anchor}_hsmetrics",
        plot=table.plot(
            data_by_sample,
            _get_table_headers(),
            {
                "id": f"{module.anchor}_hsmetrics_table",
                "namespace": "HsMetrics",
                "scale": "RdYlGn",
                "min": 0,
                "title": "Picard HsMetrics",
            },
        ),
    )
    tbases = _add_target_bases(module, data_by_sample)
    module.add_section(
        name=tbases["name"],
        anchor=tbases["anchor"],
        description=tbases["description"],
        plot=tbases["plot"],
    )
    hs_pen_plot = hs_penalty_plot(module, data_by_sample)
    if hs_pen_plot is not None:
        module.add_section(
            name="HS Penalty",
            anchor=f"{module.anchor}_hsmetrics_hs_penalty",
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
    return len(data_by_sample)


def _general_stats_table(module, data):
    """
    Generate table header configs for the General Stats table,
    add config and data to the base module.
    """
    # Look for a user config of which table columns we should use
    picard_config = getattr(config, "picard_config", {})
    genstats_table_cols = picard_config.get("HsMetrics_genstats_table_cols", [])
    genstats_table_cols_hidden = picard_config.get("HsMetrics_genstats_table_cols_hidden", [])

    headers = {}
    # Custom general stats columns
    if len(genstats_table_cols) or len(genstats_table_cols_hidden):
        for k, v in _generate_table_header_config(genstats_table_cols, genstats_table_cols_hidden).items():
            headers[k] = v

    # Default General Stats headers
    else:
        headers["FOLD_ENRICHMENT"] = {
            "title": "Fold Enrichment",
            "min": 0,
            "format": "{:,.0f}",
            "scale": "Blues",
            "suffix": " X",
        }
        headers["MEDIAN_TARGET_COVERAGE"] = {
            "title": "Median Target Coverage",
            "description": "The median coverage of reads that mapped to target regions of an experiment.",
            "min": 0,
            "suffix": "X",
            "scale": "GnBu",
        }
        try:
            covs = picard_config["general_stats_target_coverage"]
            assert isinstance(covs, list)
            assert len(covs) > 0
            covs = [str(i) for i in covs]
            log.debug(f"Custom picad coverage thresholds: {', '.join([i for i in covs])}")
        except (KeyError, AttributeError, TypeError, AssertionError):
            covs = ["30"]
        for c in covs:
            headers[f"PCT_TARGET_BASES_{c}X"] = {
                "id": f"{module.anchor}_target_bases_{c}X",
                "title": f"Target Bases ≥ {c}X",
                "description": f"Percent of target bases with coverage ≥ {c}X",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "format": "{:,.0f}",
                "scale": "RdYlGn",
                "modify": lambda x: util.multiply_hundred(x),
            }
    module.general_stats_addcols(data, headers, namespace="HsMetrics")


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
        HsMetrics_table_cols_hidden = [
            "BAIT_TERRITORY",
            "TOTAL_READS",
            "TARGET_TERRITORY",
            "AT_DROPOUT",
            "GC_DROPOUT",
        ]

    return _generate_table_header_config(HsMetrics_table_cols, HsMetrics_table_cols_hidden)


def _generate_table_header_config(table_cols, hidden_table_cols):
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

    headers = dict()
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
                headers[h]["title"] = headers[h]["title"]
                headers[h]["modify"] = lambda x: x * 100.0
                headers[h]["max"] = 100
                headers[h]["suffix"] = "%"

            elif h.find("READS") > -1:
                headers[h]["title"] = f"{config.read_count_prefix} {headers[h]['title']}"
                headers[h]["modify"] = lambda x: x * config.read_count_multiplier
                headers[h]["shared_key"] = "read_count"

            elif h.find("BASES") > -1:
                headers[h]["title"] = f"{config.base_count_prefix} {headers[h]['title']}"
                headers[h]["modify"] = lambda x: x * config.base_count_multiplier
                headers[h]["shared_key"] = "base_count"

            # Manual capitilisation for some strings
            headers[h]["title"] = headers[h]["title"].replace("Pf", "PF").replace("snp", "SNP")

            if h in hidden_table_cols:
                headers[h]["hidden"] = True

    return headers


def _add_target_bases(self, data):
    data_clean = defaultdict(dict)
    max_non_zero_cov = 0
    for s in data:
        for h in data[s]:
            if h.startswith("PCT_TARGET"):
                cov = int(h.replace("PCT_TARGET_BASES_", "")[:-1])
                bases_pct = data[s][h]
                data_clean[s][cov] = bases_pct * 100.0
                if bases_pct > 0 and cov > max_non_zero_cov:
                    max_non_zero_cov = cov

    pconfig = {
        "id": f"{self.anchor}_percentage_target_bases",
        "title": f"{self.name}: Percentage of target bases",
        "xlab": "Fold Coverage",
        "ylab": "Pct of bases",
        "ymax": 100,
        "ymin": 0,
        "xmin": 0,
        "xmax": max_non_zero_cov,
        "tt_label": "<b>{point.x}X</b>: {point.y:.2f}%",
    }
    return {
        "name": "Target Region Coverage",
        "anchor": f"{self.anchor}_hsmetrics_target_bases",
        "description": "The percentage of all target bases with at least <code>x</code> fold coverage.",
        "plot": linegraph.plot(data_clean, pconfig),
    }


def hs_penalty_plot(self, data):
    data_clean = defaultdict(dict)
    any_non_zero = False
    for s in data:
        for h in data[s]:
            if h.startswith("HS_PENALTY"):
                data_clean[s][int(h.lstrip("HS_PENALTY_").rstrip("X"))] = data[s][h]
                if data[s][h] > 0:
                    any_non_zero = True

    pconfig = {
        "id": f"{self.anchor}_hybrid_selection_penalty",
        "title": f"{self.name}: Hybrid Selection Penalty",
        "xlab": "Fold Coverage",
        "ylab": "Penalty",
        "ymin": 0,
        "xmin": 0,
        "xDecimals": False,
        "tt_label": "<b>{point.x}X</b>: {point.y:.2f}%",
    }

    if any_non_zero:
        return linegraph.plot(data_clean, pconfig)
