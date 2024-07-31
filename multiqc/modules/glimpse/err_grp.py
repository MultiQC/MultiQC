"""MultiQC submodule to parse output from GLIMPSE_concordance"""

import gzip
import logging
import os
from typing import Dict, Union

from multiqc import BaseMultiqcModule
from multiqc.plots import linegraph

log = logging.getLogger(__name__)

EXPECTED_COLUMNS = [
    "#GCsS",
    "id",
    "n_genotypes",
    "mean_AF",
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
    "best_gt_rsquared",
    "imputed_ds_rsquared",
]


def parse_glimpse_err_grp(module: BaseMultiqcModule) -> int:
    """Find Glimpse concordance errors by allele frequency bin groups logs and parse their data"""

    metrics_by_idn_by_var_by_sample: Dict[str, Dict[str, Dict[str, Dict[str, Union[int, float, str]]]]] = dict()
    for f in module.find_log_files("glimpse/err_grp", filecontents=False, filehandles=False):
        with gzip.open(os.path.join(f["root"], f["fn"])) as f_gz:
            lines = [line.decode().rstrip() for line in f_gz.readlines()]

        file_data = parse_err_grp_report(lines)
        if not file_data:
            continue

        module.add_data_source(f, section="err_grp")

        if f["s_name"] in metrics_by_idn_by_var_by_sample:
            log.debug(f"Duplicate sample name '{f['s_name']}' found from {f['fn']}. Overwriting")
        metrics_by_idn_by_var_by_sample[f["s_name"]] = file_data

    metrics_by_idn_by_var_by_sample = module.ignore_samples(metrics_by_idn_by_var_by_sample)
    n_reports_found = len(metrics_by_idn_by_var_by_sample)
    if n_reports_found == 0:
        return 0
    log.info(f"Found {n_reports_found} report(s) by samples")

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed report data to a file (restructure first)
    module.write_data_file(metrics_by_idn_by_var_by_sample, "multiqc_glimpse_err_grp")

    GCsSAF = "GCsSAF"
    GCsIAF = "GCsIAF"
    GCsVAF = "GCsVAF"
    vtypes = [GCsSAF, GCsIAF, GCsVAF]
    samples = metrics_by_idn_by_var_by_sample.keys()
    data_best_gt_rsquared: Dict[str, Dict[str, Dict]] = {v: {s: {} for s in samples} for v in vtypes}
    data_imputed_ds_rsquared: Dict[str, Dict[str, Dict]] = {v: {s: {} for s in samples} for v in vtypes}
    for sname, metrics_by_idn_by_var in metrics_by_idn_by_var_by_sample.items():
        for var_type, metrics_by_idn in metrics_by_idn_by_var.items():
            for idn, metrics in metrics_by_idn.items():
                data_best_gt_rsquared[var_type][sname].update({metrics["mean_af"]: metrics["best_gt_rsquared"]})
                data_imputed_ds_rsquared[var_type][sname].update({metrics["mean_af"]: metrics["imputed_ds_rsquared"]})

    # Make a table summarising the stats across all samples
    accuracy_plot(
        module,
        [
            data_best_gt_rsquared["GCsSAF"],
            data_imputed_ds_rsquared["GCsSAF"],
            data_best_gt_rsquared["GCsIAF"],
            data_imputed_ds_rsquared["GCsIAF"],
            data_best_gt_rsquared["GCsVAF"],
            data_imputed_ds_rsquared["GCsVAF"],
        ],
    )

    # Return the number of logs that were found
    return n_reports_found


def accuracy_plot(module, data):
    module.add_section(
        name="Genotype concordance by allele frequency bin",
        anchor="glimpse-err-grp-plot-section",
        description=(
            "Stats parsed from <code>GLIMPSE2_concordance</code> output, and summarized across allele frequency bins."
        ),
        plot=linegraph.plot(
            data,
            pconfig={
                "data_labels": [
                    {
                        "name": "Best genotype r-squared (SNPs)",
                        "ylab": "Best genotype r-squared (SNPs)",
                    },
                    {
                        "name": "Imputed dosage r-squared (SNPs)",
                        "ylab": "Imputed dosage r-squared (SNPs)",
                    },
                    {
                        "name": "Best genotype r-squared (indels)",
                        "ylab": "Best genotype r-squared (indels)",
                    },
                    {
                        "name": "Imputed dosage r-squared (indels)",
                        "ylab": "Imputed dosage r-squared (indels)",
                    },
                    {
                        "name": "Best genotype r-squared (SNPs + indels)",
                        "ylab": "Best genotype r-squared (SNPs + indels)",
                    },
                    {
                        "name": "Imputed dosage r-squared (SNPs + indels)",
                        "ylab": "Imputed dosage r-squared (SNPs + indels)",
                    },
                ],
                "id": "glimpse-err-grp-plot",
                "xlab": "Minor allele frequency",
                "logswitch": True,  # Show the 'Log10' switch?
                "logswitch_active": True,  # Initial display with 'Log10' active?
                "logswitch_label": "Log10",  # Label for 'Log10' button
                "xmin": 0,
                "xmax": 0.5,
                "ymin": 0,
                "ymax": 1.1,
                "title": "Glimpse concordance by allele frequency bins",
            },
        ),
    )


def parse_err_grp_report(lines) -> Dict[str, Dict[str, Dict[str, Union[int, float, str]]]]:
    """
    Example:
    #Genotype concordance by allele frequency bin (SNPs)
    #GCsSAF id n_genotypes mean_AF #val_gt_RR #val_gt_RA #val_gt_AA filtered_gp RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches_rate_percent
    GCsSAF 0 2838 0.001590 2830 8 0 0 2829 8 0 1 0 0 0.035 0.000 0.000 0.888575 0.882978
    GCsSAF 1 367 0.022783 335 31 1 0 335 30 0 0 1 1 0.000 3.226 100.000 0.938664 0.956999
    GCsSAF 2 130 0.070898 112 17 1 0 112 17 0 0 0 1 0.000 0.000 100.000 0.948173 0.918263
    GCsSAF 3 100 0.142570 63 36 1 0 63 36 1 0 0 0 0.000 0.000 0.000 1.000000 1.000000
    GCsSAF 4 98 0.362634 51 21 26 0 51 21 26 0 0 0 0.000 0.000 0.000 1.000000 0.999847
    #Genotype concordance by allele frequency bin (indels)

    Returns a dictionary with the contig name (rname) as the key and the rest of the fields as a dictionary
    """
    parsed_data: Dict[str, Dict[str, Dict[str, Union[int, float, str]]]] = {}
    expected_header = "#Genotype concordance by allele frequency bin (SNPs)"
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
            n_genotypes,
            mean_af,
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
            best_gt_rsquared,
            imputed_ds_rsquared,
        ) = fields
        if variants not in parsed_data:
            parsed_data[variants] = {}
        parsed_data[variants][idn] = dict(
            variants=str(variants),
            idn=int(idn),
            n_genotypes=int(n_genotypes),
            mean_af=float(mean_af),
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
            best_gt_rsquared=float(best_gt_rsquared),
            imputed_ds_rsquared=float(imputed_ds_rsquared),
        )

    return parsed_data
