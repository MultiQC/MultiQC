"""MultiQC submodule to parse output from GLIMPSE_concordance"""

import gzip
import logging
import os
from typing import Dict, Union

from multiqc import BaseMultiqcModule
from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)

EXPECTED_COLUMNS = [
    "#GCsS",
    "id",
    "sample_name",
    "#val_gt_RR",
    "#val_gt_RA",
    "#val_gt_AA",
    "#filtered_gp",
    "RR_hom_matches",
    "RA_het_matches",
    "AA_hom_matches",
    "RR_hom_mismatches",
    "RA_het_mismatches",
    "AA_hom_mismatches",
    "RR_hom_mismatches_rate_percent",
    "RA_het_mismatches_rate_percent",
    "AA_hom_mimatches",
    "non_reference_discordanc_rate_percent",
    "best_gt_rsquared",
    "imputed_ds_rsquared",
]


def parse_glimpse_err_spl(module: BaseMultiqcModule) -> int:
    """Find Glimpse concordance by samples logs and parse their data"""

    glimpse_err_spl: Dict = dict()
    allfiles: Dict = dict()
    for f in module.find_log_files("glimpse/err_spl", filecontents=False, filehandles=False):
        with gzip.open(os.path.join(f["root"], f["fn"])) as f_gz:
            lines = [line.decode() for line in f_gz.readlines()]

        parsed_data = parse_err_spl_report(lines)
        if len(parsed_data) > 1:
            if f["s_name"] in allfiles:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            module.add_data_source(f, section="err_spl")
            # Filter to strip out ignored sample names
            allfiles[f["s_name"]] = f

            shared_keys = set(glimpse_err_spl.keys()) & set(parsed_data.keys())
            if len(shared_keys) > 0:
                log.debug(f"Duplicate sample name found! Overwriting: {shared_keys}")
            glimpse_err_spl.update(module.ignore_samples(parsed_data))

    n_reports_found = len(glimpse_err_spl)
    if n_reports_found == 0:
        return 0

    log.info(f"Found {n_reports_found} report(s) by samples")

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed report data to a file (restructure first)
    module.write_data_file(glimpse_err_spl, "multiqc_glimpse_err_spl")

    headers = {
        "val_gt_RR": {
            "title": "Genotype Reference-Reference",
            "description": "Number of genotypes classified as Reference-Reference",
            "min": 0,
            "scale": "BuGn",
        },
        "val_gt_RA": {
            "title": "Genotype Reference-Alternate",
            "description": "Number of genotypes classified as Reference-Alternate",
            "min": 0,
            "scale": "BuGn",
        },
        "val_gt_AA": {
            "title": "Genotype Alternate-Alternate",
            "description": "Number of genotypes classified as Alternate-Alternate",
            "min": 0,
            "scale": "BuGn",
        },
        "filtered_gp": {
            "title": "Genotype filtered",
            "description": "Number of genotypes filtered",
            "min": 0,
            "scale": "OrRd",
        },
        "RR_hom_matches": {
            "title": "Reference-Reference homozygous matches",
            "description": "Number of Reference-Reference hom matches",
            "min": 0,
            "scale": "PuBuGn",
            "hidden": True,
        },
        "RA_het_matches": {
            "title": "Reference-Alternate heterozygous matches",
            "description": "Number of Reference-Alternate het matches",
            "min": 0,
            "scale": "PuBuGn",
            "hidden": True,
        },
        "AA_hom_matches": {
            "title": "Alternate-Alternate homozygous matches",
            "description": "Number of Alternate-Alternate hom matches",
            "min": 0,
            "scale": "PuBuGn",
            "hidden": True,
        },
        "RR_hom_mismatches": {
            "title": "Reference-Reference homozygous mismatches",
            "description": "Number of Reference-Reference hom mismatches",
            "min": 0,
            "scale": "PuRd",
            "hidden": True,
        },
        "RA_het_mismatches": {
            "title": "Reference-Alternate heterozygous mismatches",
            "description": "Number of Reference-Alternate het mismatches",
            "min": 0,
            "scale": "PuRd",
            "hidden": True,
        },
        "AA_hom_mismatches": {
            "title": "Alternate-Alternate homozygous mismatches",
            "description": "Number of Alternate-Alternate hom mismatches",
            "min": 0,
            "scale": "PuRd",
            "hidden": True,
        },
        "RR_hom_mismatches_rate_percent": {
            "title": "Reference-Reference homozygous mismatches rate",
            "description": "Rate of Reference-Reference hom mismatches",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "YlOrRd",
        },
        "RA_het_mismatches_rate_percent": {
            "title": "Reference-Alternate heterozygous mismatches rate",
            "description": "Rate of Reference-Alternate het mismatches",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "YlOrRd",
        },
        "AA_hom_mismatches_rate_percent": {
            "title": "Alternate-Alternate homozygous mismatches rate",
            "description": "Rate of Alternate-Alternate hom mismatches",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "YlOrRd",
        },
        "non_reference_discordance_rate_percent": {
            "title": "Non-reference discordance rate",
            "description": "Rate of non-reference discordance",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "YlOrRd",
        },
        "best_gt_rsquared": {
            "title": "Best genotype r-squared",
            "description": "Best genotype r-squared",
            "min": 0,
            "max": 1,
            "scale": "YlGn",
        },
        "imputed_ds_rsquared": {
            "title": "Imputed dosage r-squared",
            "description": "Imputed dosage r-squared",
            "min": 0,
            "max": 1,
            "scale": "YlGn",
        },
    }

    # Keep only items from all variants (SNPs + indels)
    data = {sample: val for sample, vtype in glimpse_err_spl.items() for key, val in vtype.items() if key == "GCsV"}

    # Make a table summarising the stats across all samples
    summary_table(module, data, headers)

    # Add to the general stats table
    general_headers = ["best_gt_rsquared", "imputed_ds_rsquared"]
    module.general_stats_addcols(filter_err_spl_data(data, general_headers), headers)

    # Return the number of logs that were found
    return n_reports_found


def summary_table(module, data, headers):
    module.add_section(
        name="Genotype concordance by samples",
        anchor="glimpse-err-spl-table-section",
        description=("Stats parsed from <code>GLIMPSE2_concordance</code> output, and summarized across all samples."),
        plot=table.plot(
            data,
            headers,
            pconfig={
                "id": "glimpse-err-spl-table",
                "title": "Glimpse concordance by sample",
            },
        ),
    )


def filter_err_spl_data(data, keys):
    """Filter out unwanted keys from the parsed data"""
    return {sample: {k: v for k, v in datas.items() if k in keys} for sample, datas in data.items()}


def parse_err_spl_report(lines) -> Dict[str, Dict[str, Union[int, float]]]:
    """
    Example:
    #Genotype concordance by sample (SNPs)
    #GCsS id sample_name #val_gt_RR #val_gt_RA #val_gt_AA #filtered_gp RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches non_reference_discordanc_rate_percent best_gt_rsquared imputed_ds_rsquared
    GCsS 0 NA12878 851 3 0 0 602 0 0 3 1 0 0.496 100.000 0.000 100.000 0.000008 0.000008
    #Genotype concordance by sample (indels)
    #GCsI id sample_name #val_gt_RR #val_gt_RA #val_gt_AA #filtered_gp RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches non_reference_discordanc_rate_percent best_gt_rsquared imputed_ds_rsquared
    GCsI 0 NA12878 0 0 0 0 0 0 0 0 0 0 0.000 0.000 0.000 0.000 0.000000 0.000000
    #Genotype concordance by sample (Variants: SNPs + indels)
    #GCsV id sample_name #val_gt_RR #val_gt_RA #val_gt_AA #filtered_gp RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches non_reference_discordanc_rate_percent best_gt_rsquared imputed_ds_rsquared
    GCsV 0 NA12878 851 3 0 0 602 0 0 3 1 0 0.496 100.000 0.000 100.000 0.000008 0.000008

    Returns a dictionary with the contig name (rname) as the key and the rest of the fields as a dictionary
    """
    parsed_data = {}
    expected_header = "#Genotype concordance by sample (SNPs)\n"
    if lines[0] != expected_header:
        logging.warning(f"Expected header for GLIMPSE2_concordance: {expected_header}, got: {lines[0]}.")
        return {}

    for line in lines[1:]:
        fields = line.strip().split(" ")
        if fields[0][0] == "#":  # Skip comments
            continue
        if len(fields) != len(EXPECTED_COLUMNS):
            logging.warning(f"Skipping line with {len(fields)} fields, expected {len(EXPECTED_COLUMNS)}: {line}")
        (
            variants,
            idn,
            sample_name,
            val_gt_RR,
            val_gt_RA,
            val_gt_AA,
            filtered_gp,
            RR_hom_matches,
            RA_het_matches,
            AA_hom_matches,
            RR_hom_mismatches,
            RA_het_mismatches,
            AA_hom_mismatches,
            RR_hom_mismatches_rate_percent,
            RA_het_mismatches_rate_percent,
            AA_hom_mismatches_rate_percent,
            non_reference_discordance_rate_percent,
            best_gt_rsquared,
            imputed_ds_rsquared,
        ) = fields
        if sample_name not in parsed_data:
            parsed_data[sample_name] = {}
        parsed_data[sample_name][variants] = dict(
            idn=int(idn),
            val_gt_RR=int(val_gt_RR),
            val_gt_RA=int(val_gt_RA),
            val_gt_AA=int(val_gt_AA),
            filtered_gp=int(filtered_gp),
            RR_hom_matches=int(RR_hom_matches),
            RA_het_matches=int(RA_het_matches),
            AA_hom_matches=int(AA_hom_matches),
            RR_hom_mismatches=int(RR_hom_mismatches),
            RA_het_mismatches=int(RA_het_mismatches),
            AA_hom_mismatches=int(AA_hom_mismatches),
            RR_hom_mismatches_rate_percent=float(RR_hom_mismatches_rate_percent),
            RA_het_mismatches_rate_percent=float(RA_het_mismatches_rate_percent),
            AA_hom_mismatches_rate_percent=float(AA_hom_mismatches_rate_percent),
            non_reference_discordance_rate_percent=float(non_reference_discordance_rate_percent),
            best_gt_rsquared=float(best_gt_rsquared),
            imputed_ds_rsquared=float(imputed_ds_rsquared),
        )

    return parsed_data
