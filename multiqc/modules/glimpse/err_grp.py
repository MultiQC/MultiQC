""" MultiQC submodule to parse output from GLIMPSE_concordance """

import logging
from typing import Dict, Union

from multiqc.plots import linegraph

# Initialise the loggerq
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


class ErrGrpReportMixin:
    """Mixin class, loaded by main Glimpse MultiqcModule class."""
    # parsing functions -------------------------------------------------------------

    def parse_glimpse_err_grp(self):
        """Find Glimpse concordance errors by allele frequency bin groups logs and parse their data"""

        self.glimpse_err_grp = dict()
        for f in self.find_log_files("glimpse/err_grp", filehandles=True):
            parsed_data = parse_err_grp_report([line.rstrip() for line in f["f"]])
            if len(parsed_data) > 1:
                if f["s_name"] in self.glimpse_err_grp:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                self.add_data_source(f, section="err_grp")
                # Filter to strip out ignored sample names
                self.glimpse_err_grp[f["s_name"]] = self.ignore_samples(parsed_data)

        n_reports_found = len(self.glimpse_err_grp)
        if n_reports_found == 0:
            return 0

        log.info(f"Found {n_reports_found} report(s) by allele frequency bin.")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file (restructure first)
        self.write_data_file(self.glimpse_err_grp, "multiqc_glimpse_err_grp")

        data_GCsVAF = {
            sname: data
            for sname, dataf in self.glimpse_err_grp.items()
            for vtype, data in dataf.items()
            if vtype == "GCsVAF"
        }

        data_best_gt_rsquared = {}
        data_imputed_ds_rsquared = {}
        for sname, dataf in data_GCsVAF.items():
            data_best_gt_rsquared[sname] = {}
            data_imputed_ds_rsquared[sname] = {}
            for idn, data in dataf.items():
                data_best_gt_rsquared[sname].update({data["mean_af"]: data["best_gt_rsquared"]})
                data_imputed_ds_rsquared[sname].update({data["mean_af"]: data["imputed_ds_rsquared"]})

        # Make a table summarising the stats across all samples
        self.accuracy_plot(
            [
                data_best_gt_rsquared,
                data_imputed_ds_rsquared,
            ]
        )

        # Return the number of logs that were found
        return n_reports_found

    def accuracy_plot(self, data):
        self.add_section(
            name="Genotype concordance accuracy by allele frequency bins",
            anchor="glimpse-err-grp-plot-section",
            description=(
                "Stats parsed from <code>GLIMPSE2_concordance</code> output, and summarized across allele frequency bins."
            ),
            plot=linegraph.plot(
                data,
                pconfig={
                    "data_labels": [
                        {
                            "name": "Best genotype r-squared",
                            "ylab": "Best genotype r-squared",
                            "xlab": "Minor allele frequency",
                            "xlog": True,
                            "xmin": 0,
                            "xmax": 0.5,
                            "ylog": True,
                            "ymin": 0,
                            "ymax": 1.1,
                        },
                        {
                            "name": "Imputed dosage r-squared",
                            "ylab": "Imputed dosage r-squared",
                            "xlab": "Minor allele frequency",
                            "xlog": True,
                            "xmin": 0,
                            "xmax": 0.5,
                            "ylog": True,
                            "ymin": 0,
                            "ymax": 1.1,
                        },
                    ],
                    "id": "glimpse-err-spl-table",
                    "title": "Glimpse concordance: accuracy by allele frequency bins",
                },
            ),
        )


def parse_err_grp_report(lines) -> Dict[str, Dict[str, Union[int, float]]]:
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
    parsed_data = {}
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
