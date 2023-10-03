""" MultiQC submodule to parse output from Picard VariantCallingMetrics """

import logging
from collections import OrderedDict

from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


def parse_reports(parent_module):
    """Find Picard VariantCallingMetrics reports and process their data"""

    # get data
    data = collect_data(parent_module)

    # Filter to strip out ignored sample names
    data = parent_module.ignore_samples(data)

    # Reference data in parent module
    parent_module.picard_variantCalling_data = data

    if len(data) > 0:
        derive_data(data)

        # Write parsed data to a file
        parent_module.write_data_file(data, "multiqc_picard_variantCalling")

        parent_module.general_stats_headers["DBSNP_TITV"] = {
            "title": "TiTV ratio (known)",
            "description": "The Transition/Transversion ratio of the passing bi-allelic SNP calls made at SNP-database sites.",
            "min": 0,
            "scale": "Blues",
            "shared_key": "titv_ratio",
        }

        parent_module.general_stats_headers["NOVEL_TITV"] = {
            "title": "TiTV ratio (novel)",
            "description": "The Transition/Transversion ratio of the passing bi-allelic SNP calls made at non-SNP-database sites.",
            "min": 0,
            "scale": "Blues",
            "shared_key": "titv_ratio",
        }

        parent_module.general_stats_headers["DBSNP_INS_DEL_RATIO"] = {
            "title": "InDel ratio (known)",
            "description": "The Insertion / Deletion ratio of the passing bi-allelic SNP calls made at SNP-database sites.",
            "min": 0,
            "scale": "Greens",
            "shared_key": "indel_ratio",
            "hidden": True,
        }

        parent_module.general_stats_headers["NOVEL_INS_DEL_RATIO"] = {
            "title": "InDel ratio (novel)",
            "description": "The Insertion / Deletion ratio of the passing bi-allelic SNP calls made at non-SNP-database sites.",
            "min": 0,
            "scale": "Greens",
            "shared_key": "indel_ratio",
            "hidden": True,
        }

        parent_module.general_stats_headers["total_called_variants_known"] = {
            "title": "Called Variants (known)",
            "description": "Total counts of variants in SNP-database sites.",
            "shared_key": "variant_count",
            "min": 0,
            "format": "{0:,.0f}",
            "hidden": True,
        }

        parent_module.general_stats_headers["total_called_variants_novel"] = {
            "title": "Called Variants (novel)",
            "description": "Total counts of variants in non-SNP-database sites.",
            "shared_key": "variant_count",
            "min": 0,
            "format": "{0:,.0f}",
            "hidden": True,
        }

        for s_name in data:
            if s_name not in parent_module.general_stats_data:
                parent_module.general_stats_data[s_name] = dict()
            parent_module.general_stats_data[s_name].update(data[s_name])

        # Variant Counts Bargraph
        parent_module.add_section(
            name="Variant Types",
            anchor="picard-variants-types",
            description="Variants that have been called, looking at variant types. Optionally filtered on label.",
            helptext="""
            Only passing variants are shown (i.e. non-filtered).\n
            SNPs are bi-allelic.\n
            Complex InDels are both an insertion and a deletion.
            """,
            plot=compare_variant_type_plot(data),
        )

        # Variant Counts Table
        parent_module.add_section(
            name="Variant Labels",
            anchor="picard-variants-labels",
            description="Variants that have been called, comparing with known variant sites.",
            helptext="""
            Only passing variants are shown (i.e. non-filtered).\n
            Variants contain bi-allelic SNPs, multi-allelic SNPs, simple and complex inserts and deletions.
            """,
            plot=compare_variants_label_plot(data),
        )

    return len(data)


def collect_data(parent_module):
    """Find Picard VariantCallingMetrics reports and parse their data"""

    data = dict()
    for file_meta in parent_module.find_log_files("picard/variant_calling_metrics", filehandles=True):
        s_name = None
        for header, value in table_in(file_meta["f"], pre_header_string="## METRICS CLASS"):
            if header == "SAMPLE_ALIAS":
                s_name = value
                if s_name in data:
                    log.debug("Duplicate sample name found in {}! Overwriting: {}".format(file_meta["fn"], s_name))
                data[s_name] = OrderedDict()
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
        for value_name in ["TOTAL_SNPS", "TOTAL_COMPLEX_INDELS", "TOTAL_MULTIALLELIC_SNPS", "TOTAL_INDELS"]:
            total_called_variants = total_called_variants + int(values[value_name])
        values["total_called_variants"] = total_called_variants

        # Sum all variants that have been called and are known
        total_called_variants_known = 0
        for value_name in ["NUM_IN_DB_SNP", "NUM_IN_DB_SNP_COMPLEX_INDELS", "NUM_IN_DB_SNP_MULTIALLELIC"]:
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


def compare_variant_type_plot(data):
    """Return HTML for the Variant Counts barplot"""
    keys = OrderedDict()
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
        "hide_zero_cats": False,
        "data_labels": [{"name": "Total"}, {"name": "Known"}, {"name": "Novel"}],
    }
    return bargraph.plot(
        data=[total_variants, known_variants, novel_variants], cats=[keys, keys, keys], pconfig=plot_conf
    )


def compare_variants_label_plot(data):
    """Return HTML for the Compare variants plot"""
    keys = OrderedDict()

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
