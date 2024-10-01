"""MultiQC submodule to parse output from Picard VariantCallingMetrics"""

import logging
from typing import Dict

from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(module):
    """Find Picard VariantCallingMetrics reports and process their data"""

    # get data
    data_by_sample = collect_data(module)

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)
    if len(data_by_sample) == 0:
        return set()

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    derive_data(data_by_sample)

    # Write parsed data to a file
    module.write_data_file(data_by_sample, "multiqc_picard_variantCalling")

    headers = {
        "DBSNP_TITV": {
            "title": "TiTV ratio (known)",
            "description": "The Transition/Transversion ratio of the passing bi-allelic SNP calls made at SNP-database sites.",
            "min": 0,
            "scale": "Blues",
            "shared_key": "titv_ratio",
        },
        "NOVEL_TITV": {
            "title": "TiTV ratio (novel)",
            "description": "The Transition/Transversion ratio of the passing bi-allelic SNP calls made at non-SNP-database sites.",
            "min": 0,
            "scale": "Blues",
            "shared_key": "titv_ratio",
        },
        "DBSNP_INS_DEL_RATIO": {
            "title": "InDel ratio (known)",
            "description": "The Insertion / Deletion ratio of the passing bi-allelic SNP calls made at SNP-database sites.",
            "min": 0,
            "scale": "Greens",
            "shared_key": "indel_ratio",
            "hidden": True,
        },
        "NOVEL_INS_DEL_RATIO": {
            "title": "InDel ratio (novel)",
            "description": "The Insertion / Deletion ratio of the passing bi-allelic SNP calls made at non-SNP-database sites.",
            "min": 0,
            "scale": "Greens",
            "shared_key": "indel_ratio",
            "hidden": True,
        },
        "total_called_variants_known": {
            "title": "Called Variants (known)",
            "description": "Total counts of variants in SNP-database sites.",
            "shared_key": "variant_count",
            "min": 0,
            "format": "{0:,.0f}",
            "hidden": True,
        },
        "total_called_variants_novel": {
            "title": "Called Variants (novel)",
            "description": "Total counts of variants in non-SNP-database sites.",
            "shared_key": "variant_count",
            "min": 0,
            "format": "{0:,.0f}",
            "hidden": True,
        },
    }
    module.general_stats_addcols(data_by_sample, headers, namespace="VariantCallingMetrics")

    # Variant Calling Metrics Table
    module.add_section(
        name="Variant Calling Metrics",
        anchor=f"{module.id}_variantcallingmetrics",
        plot=table.plot(
            data=data_by_sample,
            headers=get_table_headers(),
            pconfig={
                "id": f"{module.id}_variantcallingmetrics_table",
                "namespace": "VariantCallingMetrics",
                "title": "Variant Calling Metrics",
            },
        ),
    )

    # Variant Counts Bargraph
    module.add_section(
        name="Variant Types",
        anchor="picard-variants-types",
        description="Variants that have been called, looking at variant types. Optionally filtered on label.",
        helptext="""
        Only passing variants are shown (i.e. non-filtered).\n
        SNPs are bi-allelic.\n
        Complex InDels are both an insertion and a deletion.
        """,
        plot=compare_variant_type_plot(data_by_sample),
    )

    # Variant Counts Table
    module.add_section(
        name="Variant Labels",
        anchor="picard-variants-labels",
        description="Variants that have been called, comparing with known variant sites.",
        helptext="""
        Only passing variants are shown (i.e. non-filtered).\n
        Variants contain bi-allelic SNPs, multi-allelic SNPs, simple and complex inserts and deletions.
        """,
        plot=compare_variants_label_plot(data_by_sample),
    )

    return data_by_sample.keys()


def collect_data(module):
    """Find Picard VariantCallingMetrics reports and parse their data"""

    data: Dict = dict()
    for f in module.find_log_files("picard/variant_calling_metrics", filehandles=True):
        s_name = None
        for header, value in table_in(f["f"], pre_header_string="## METRICS CLASS"):
            if header == "SAMPLE_ALIAS":
                s_name = value
                if s_name in data:
                    log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
                data[s_name] = dict()
            else:
                data[s_name][header] = value
    return data


def table_in(filehandle, pre_header_string):
    """Generator that assumes a table starts the line after a given string"""

    in_histogram = False
    next_is_header = False
    headers = list()
    for line in stripped(filehandle):
        if not in_histogram and line.startswith(pre_header_string):
            in_histogram = True
            next_is_header = True
        elif in_histogram and next_is_header:
            next_is_header = False
            headers = line.split("\t")
        elif in_histogram:
            values = line.split("\t")
            if values != [""]:
                for couple in zip(headers, values):
                    yield couple


def derive_data(data):
    """Based on the data derive additional data"""

    for s_name, values in data.items():
        # setup holding variable

        # Sum all variants that have been called
        total_called_variants = 0
        for value_name in [
            "TOTAL_SNPS",
            "TOTAL_COMPLEX_INDELS",
            "TOTAL_MULTIALLELIC_SNPS",
            "TOTAL_INDELS",
        ]:
            total_called_variants = total_called_variants + int(values[value_name])
        values["total_called_variants"] = total_called_variants

        # Sum all variants that have been called and are known
        total_called_variants_known = 0
        for value_name in [
            "NUM_IN_DB_SNP",
            "NUM_IN_DB_SNP_COMPLEX_INDELS",
            "NUM_IN_DB_SNP_MULTIALLELIC",
        ]:
            total_called_variants_known = total_called_variants_known + int(values[value_name])
        total_called_variants_known = (
            total_called_variants_known + int(values["TOTAL_INDELS"]) - int(values["NOVEL_INDELS"])
        )
        values["total_called_variants_known"] = total_called_variants_known

        # Extrapolate the total novel variants
        values["total_called_variants_novel"] = total_called_variants - total_called_variants_known


def stripped(iterator):
    """Generator to strip string of whitespace"""
    for item in iterator:
        yield item.strip()


def get_table_headers():
    """return metrics table header dictionary"""
    headers = {
        "HET_HOMVAR_RATIO": {
            "title": "Het/Hom Ratio",
            "description": "(count of hets)/(count of homozygous non-ref) for this sample",
            "hidden": False,
        },
        "PCT_GQ0_VARIANTS": {
            "title": "%GQ0",
            "description": "The percentage of variants in a particular sample that have a GQ score of 0.",
            "hidden": True,
        },
        "TOTAL_GQ0_VARIANTS": {
            "title": "Total GQ0",
            "description": "The total number of variants in a particular sample that have a GQ score of 0.",
            "hidden": True,
        },
        "TOTAL_HET_DEPTH": {
            "title": "Total Het Depth",
            "description": "total number of reads (from AD field) for passing bi-allelic SNP hets for this sample",
            "hidden": True,
        },
        "TOTAL_SNPS": {
            "title": "Total SNPs",
            "description": "	The number of passing bi-allelic SNPs calls (i.e. non-reference genotypes) that were examined",
            "hidden": False,
        },
        "NUM_IN_DB_SNP": {
            "title": "dbSNP SNP count",
            "description": "The number of passing bi-allelic SNPs found in dbSNP",
            "hidden": True,
        },
        "NOVEL_SNPS": {
            "title": "Novel SNPs",
            "description": "The number of passing bi-allelic SNPS called that were not found in dbSNP",
            "hidden": True,
        },
        "FILTERED_SNPS": {
            "title": "Filtered SNPs",
            "description": "The number of SNPs that are filtered",
            "hidden": True,
        },
        "PCT_DBSNP": {
            "title": "Fraction in dbSNP (SNPs)",
            "description": "The fraction of passing bi-allelic SNPs in dbSNP",
            "hidden": False,
        },
        "DBSNP_TITV": {
            "title": "Ti/Tv ratio (dbSNP)",
            "description": "The Transition/Transversion ratio of the passing bi-allelic SNP calls made at dbSNP sites",
            "hidden": False,
        },
        "NOVEL_TITV": {
            "title": "Ti/Tv ratio (Novel)",
            "description": "The Transition/Transversion ratio of the passing bi-allelic SNP calls made at non-dbSNP sites",
            "hidden": False,
        },
        "TOTAL_INDELS": {
            "title": "Indels (Total)",
            "description": "The number of passing indel calls that were examined",
            "hidden": True,
        },
        "NOVEL_INDELS": {
            "title": "Indels (Novel)",
            "description": "The number of passing indels called that were not found in dbSNP",
            "hidden": True,
        },
        "FILTERED_INDELS": {
            "title": "Indels (Filtered)",
            "description": "The number of indels that are filtered",
            "hidden": True,
        },
        "PCT_DBSNP_INDELS": {
            "title": "% in dbSNP (Indels)",
            "description": "The fraction of passing indels in dbSNP",
            "hidden": True,
        },
        "NUM_IN_DB_SNP_INDELS": {
            "title": "Indels (dbSNP)",
            "description": "The number of passing indels found in dbSNP",
            "hidden": True,
        },
        "DBSNP_INS_DEL_RATIO": {
            "title": "Indel Ratio (dbSNP)",
            "description": "The Insertion/Deletion ratio of the indel calls made at dbSNP sites",
            "hidden": False,
        },
        "NOVEL_INS_DEL_RATIO": {
            "title": "Indel Ratio (Novel)",
            "description": "The Insertion/Deletion ratio of the indel calls made at non-dbSNP sites",
            "hidden": False,
        },
        "TOTAL_MULTIALLELIC_SNPS": {
            "title": "Multiallelic SNPs (Total)",
            "description": "The number of passing multi-allelic SNP calls that were examined",
            "hidden": True,
        },
        "NUM_IN_DB_SNP_MULTIALLELIC": {
            "title": "Multiallelic SNPs (dbSNP)",
            "description": "The number of passing multi-allelic SNPs found in dbSNP",
            "hidden": True,
        },
        "TOTAL_COMPLEX_INDELS": {
            "title": "Complex Indels (Total)",
            "description": "The number of passing complex indel calls that were examined",
            "hidden": True,
        },
        "NUM_IN_DB_SNP_COMPLEX_INDELS": {
            "title": "Complex Indels (dbSNP)",
            "description": "The number of passing complex indels found in dbSNP",
            "hidden": True,
        },
        "SNP_REFERENCE_BIAS": {
            "title": "SNP Ref. Bias",
            "description": "The rate at which reference bases are observed at ref/alt heterozygous SNP sites",
            "hidden": True,
        },
        "NUM_SINGLETONS": {
            "title": "Singletons",
            "description": "For summary metrics, the number of variants that appear in only one sample. For detail metrics, the number of variants that appear only in the current sample.",
            "hidden": True,
        },
        # derived statistics
        "total_called_variants": {"title": "Total Called", "description": "Total called variants", "hidden": True},
        "total_called_variants_known": {
            "title": "Total Called (known)",
            "description": "Total called variants with membership in dbSNP",
            "hidden": True,
        },
        "total_called_variants_novel": {
            "title": "Total Called (novel)",
            "description": "Total called variants at non-dbSNP sites",
            "hidden": True,
        },
    }
    return headers


def compare_variant_type_plot(data):
    """Return HTML for the Variant Counts barplot"""
    keys = dict()
    keys["snps"] = {"name": "SNPs", "color": "#7cb5ec"}
    keys["indels"] = {"name": "InDels", "color": "#90ed7d"}
    keys["multiallelic_snps"] = {"name": "multi-allelic SNP", "color": "orange"}
    keys["complex_indels"] = {"name": "Complex InDels", "color": "#8085e9"}

    total_variants = dict()
    known_variants = dict()
    novel_variants = dict()
    for s_name, values in data.items():
        total_variants[s_name] = {
            "snps": values["TOTAL_SNPS"],
            "indels": values["TOTAL_INDELS"],
            "multiallelic_snps": values["TOTAL_MULTIALLELIC_SNPS"],
            "complex_indels": values["TOTAL_COMPLEX_INDELS"],
        }

        known_variants[s_name] = {
            "snps": values["NUM_IN_DB_SNP"],
            "indels": int(values["TOTAL_INDELS"]) - int(values["NOVEL_INDELS"]),
            "multiallelic_snps": values["NUM_IN_DB_SNP_MULTIALLELIC"],
            "complex_indels": values["NUM_IN_DB_SNP_COMPLEX_INDELS"],
        }

        novel_variants[s_name] = {
            "snps": values["NOVEL_SNPS"],
            "indels": values["NOVEL_INDELS"],
            "multiallelic_snps": int(values["TOTAL_MULTIALLELIC_SNPS"]) - int(values["NUM_IN_DB_SNP_MULTIALLELIC"]),
            "complex_indels": int(values["TOTAL_COMPLEX_INDELS"]) - int(values["NUM_IN_DB_SNP_COMPLEX_INDELS"]),
        }

    plot_conf = {
        "id": "picard_variantCallingMetrics_variant_type",
        "title": "Picard: Variants Called",
        "ylab": "Counts of Variants",
        "hide_empty": False,
        "data_labels": [{"name": "Total"}, {"name": "Known"}, {"name": "Novel"}],
    }
    return bargraph.plot(
        data=[total_variants, known_variants, novel_variants],
        cats=[keys, keys, keys],
        pconfig=plot_conf,
    )


def compare_variants_label_plot(data):
    """Return HTML for the Compare variants plot"""
    keys = dict()
    keys["total_called_variants_known"] = {"name": "Known Variants"}
    keys["total_called_variants_novel"] = {"name": "Novel Variants"}

    pconfig = {
        "id": "picard_variantCallingMetrics_variant_label",
        "title": "Picard: Variants Called",
        "ylab": "Counts of Variants",
    }

    return bargraph.plot(data, cats=keys, pconfig=pconfig)


def merge_two_dicts(x, y):
    z = x.copy()  # start with x's keys and values
    z.update(y)  # modifies z with y's keys and values & returns None
    return z
